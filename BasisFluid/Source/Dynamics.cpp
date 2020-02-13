
#include "Application.h"

#include "glm/ext.hpp"

using namespace std;

// add forces from particles (buoyancy) using a grid. Clumps particles together in the grid, then integrates against the grid to obtain basis coefficient.
void Application::AddParticleForcesToBasisFlows()
{

    // reset forces to zero
    _forceField->populateWithFunction([](float, float) {return vec2(0); });


    // add particle forces in forces grid
    //_partPos->refreshCpuData();
    vec2* particlesPointer = _partPos->getCpuDataPointer();
    //_partAges->refreshCpuData();
    float* agesPointer = _partAges->getCpuDataPointer();


    for (unsigned int i = 0; i < _partPos->_nbElements; i++) { //  everywhere use const to store end of loop, so compiler can optimize the check?
        uvec2 index = _forceField->pointToClosestIndex(particlesPointer[i]);
        index = glm::clamp(index, uvec2(0), uvec2(_nbParticlesPerCell->_nbElementsX - 1,
            _nbParticlesPerCell->_nbElementsY - 1));
        // TODO: speed up addVectorCpuData
        _forceField->addVectorCpuData(index.x, index.y, _dt * powf(_buoyancyDecayRatioWithAge, agesPointer[i]) * _buoyancyPerParticle*vec2(0, 1));
    }




    //_basisFlowParams->refreshCpuData();
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    // set flag for whether or not to use basis for forces
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        if (AllBitsSet(b.bitFlags, INTERIOR)) {
            b.bitFlags = SetBits(b.bitFlags, BASIS_FLAGS::FORCE_PROJECTION);
        }
        else {
            b.bitFlags = UnsetBits(b.bitFlags, BASIS_FLAGS::FORCE_PROJECTION);
        }
    }
    _basisFlowParams->_sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;
    //_basisFlowParams->dirtyData();

    //_basisFlowParams->refreshCpuData();
    //BasisFlow* basisFlowParamsPointer = basisFlowParams->getCpuDataPointer();



    // project forces into basis space  
    scalar_inversion_storage* vecBPointer = _vecB->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        if (AllBitsSet(basisFlowParamsPointer[iBasis].bitFlags, BASIS_FLAGS::FORCE_PROJECTION)) {
            vecBPointer[iBasis] = IntegrateBasisGrid(basisFlowParamsPointer[iBasis], _forceField.get());
        }
        else {
            vecBPointer[iBasis] = 0;
        }
    }
    _vecB->_sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    //_vecB->dirtyData();


    // inverse to obtain base weights
    InverseBBMatrix(_vecXForces.get(), _vecB.get(), BASIS_FLAGS::FORCE_PROJECTION);




    // add force weights to current basis weights
    //_vecXForces->refreshCpuData();
    scalar_inversion_storage* vecXForcesPointer = _vecXForces->getCpuDataPointer();
    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].coeff += _dt * float(vecXForcesPointer[i]);
    }
    _vecXForces->_sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    //vecXForces->dirtyData();

}




float Application::ComputeNewCenterProportions(vec2& newCenter, BasisFlow& bi, BasisFlow& bj, vec2& interBasisDist)
{
    float tempCoeff;
    tempCoeff = glm::max<float>(0.f, (1.f - abs(newCenter.x - bj.center.x) / interBasisDist.x)) *
        glm::max<float>(0.f, (1.f - abs(newCenter.y - bj.center.y) / interBasisDist.y));
    bj.newCoeff += bi.coeff * tempCoeff;
    return tempCoeff;
}






void Application::ComputeBasisAdvection()
{
    const unsigned int nbBasisFlows = _basisFlowParams->_nbElements;


    //basisFlowParams->refreshCpuData();
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    //intersectingBasesIds->refreshCpuData();
    //accelBasisCentersIds->refreshCpuData();



    for (uint i = 0; i < nbBasisFlows; i++) {
        BasisFlow& bi = basisFlowParamsPointer[i];
        bi.newCoeff = 0;
        //bi.coeffBoundary = 0;
    }

    unsigned int nbBasisTransportTooFar = 0;
    float minTotalTempCoeffs = 99999999999.f;
    uint minTotalTempCoeffsId = -1;
    float maxTotalTempCoeffs = 0.f;
    uint maxTotalTempCoeffsId = -1;


    //        float maxCoeffNorm = 0;

    for (uint i = 0; i < nbBasisFlows; i++) {

        //BasisFlow bi = basisFlowParams->getCpuData(i);
        BasisFlow& bi = basisFlowParamsPointer[i];
        vec2 avgDisplacement(0);

        if (!AllBitsSet(bi.bitFlags, INTERIOR) && !AllBitsSet(bi.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) { continue; }

        // compute displacement (I is transported by J)
        vector<CoeffTDecompressedIntersectionInfo>& intersectionInfos = _coeffsTDecompressedIntersections[i];
        for (CoeffTDecompressedIntersectionInfo& inter : intersectionInfos) {
            avgDisplacement += inter.coeff * (basisFlowParamsPointer[inter.j].coeff + _obstacleBoundaryFactorTransferOnly * basisFlowParamsPointer[inter.j].coeffBoundary);
        }


        // compute new center
        vec2 newCenter = bi.center + _dt * avgDisplacement;

        float totalTempCoeffs = 0;

        vec2 freqI(1 << bi.freqLvl.x, 1 << bi.freqLvl.y);

        vec2 interBasisDist = 0.5f * _densityMultiplierBasisHalfSize / freqI;


        if (
            abs(newCenter.x - bi.center.x) > bi.supportHalfSize().x*_densityMultiplierBasisHalfSize*0.99 ||
            abs(newCenter.y - bi.center.y) > bi.supportHalfSize().y*_densityMultiplierBasisHalfSize*0.99
            ) {
            // new center not located within immediate neighbours, must look further using basis acceleration grid
            // TODO: test this is okay by comparing a simulation where the basis does not go away from neighbors to the same sim forcing to use this block instead.
            nbBasisTransportTooFar++;
            uvec2 minId(9999);
            uvec2 maxId(0);
            for (int iX = -1; iX <= 1; iX += 2) {
                for (int iY = -1; iY <= 1; iY += 2) {
                    {
                        vec2 corner = newCenter + vec2(iX, iY)*bi.supportHalfSize()*_densityMultiplierBasisHalfSize;
                        uvec2 cellId;
                        cellId.x = glm::clamp<int>(int(floor((corner.x - _domainLeft) / (_domainRight - _domainLeft)*_accelBasisRes)), 0, _accelBasisRes - 1);
                        cellId.y = glm::clamp<int>(int(floor((corner.y - _domainBottom) / (_domainTop - _domainBottom)*_accelBasisRes)), 0, _accelBasisRes - 1);
                        minId = glm::min<uint>(cellId, minId);
                        maxId = glm::max<uint>(cellId, maxId);
                    }
                }
            }


            for (uint iX = minId.x; iX <= maxId.x; iX++) {
                for (uint iY = minId.y; iY <= maxId.y; iY++) {
                    {
                        for (unsigned int basisId : *(_accelBasisCentersIds->getCpuData_noRefresh(iX, iY)))
                        {
                            BasisFlow& bj = basisFlowParamsPointer[basisId];
                            if (bj.freqLvl == bi.freqLvl) {
                                totalTempCoeffs += ComputeNewCenterProportions(newCenter, bi, bj, interBasisDist);
                            }
                        }
                    }
                }
            }

        }
        else
        {

            // new center within immediate neighbours
            vector<unsigned int>* localIntersectingBasesIdsTransport = _intersectingBasesIdsTransport->getCpuData_noRefresh(i);
            for (auto itJ = localIntersectingBasesIdsTransport->begin(); itJ != localIntersectingBasesIdsTransport->end(); ++itJ) {
                BasisFlow& bj = basisFlowParamsPointer[*itJ];
                totalTempCoeffs += ComputeNewCenterProportions(newCenter, bi, bj, interBasisDist);
            }

        }

        // TODO: don't compute if not printing
        if (totalTempCoeffs < minTotalTempCoeffs) {
            minTotalTempCoeffsId = i;
            minTotalTempCoeffs = totalTempCoeffs;
        }
        if (totalTempCoeffs > maxTotalTempCoeffs) {
            maxTotalTempCoeffsId = i;
            maxTotalTempCoeffs = totalTempCoeffs;
        }

    }


    // add transport cofficients and set new coefficients
    // TODO: Should we also ass the boundar coefficient here?
    //     - if not, we are using the boundary flow to push the flow, but the boundary flow is not added to what is pushed. Might be okay, not sure.
    for (uint i = 0; i < nbBasisFlows; i++) {
        BasisFlow& bi = basisFlowParamsPointer[i];

        if (!AllBitsSet(bi.bitFlags, INTERIOR) && !AllBitsSet(bi.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) { continue; }

        bi.coeff = bi.newCoeff;
        bi.newCoeff = 0;
    }



    _basisFlowParams->_sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;

    basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();





    _basisFlowParams->_sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;

    basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();



    //
    // transfer energy from .coeff to .newCoeff
    //

    // normalize coefficients
    float transferCoeffs[_nbExplicitTransferFreqs];
    transferCoeffs[0] = _explicitTransfer_10;
    transferCoeffs[1] = _explicitTransfer_01;
    transferCoeffs[2] = _explicitTransfer_11;
    transferCoeffs[3] = _explicitTransfer_m10; // reverse cascade
    transferCoeffs[4] = _explicitTransfer_0m1;
    transferCoeffs[5] = _explicitTransfer_m1m1;

    float sum = 0;
    for (int iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++) {
        sum += transferCoeffs[iRelFreq];
    }
    for (int iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++) {
        transferCoeffs[iRelFreq] /= sum;
    }


    float maxAlpha = 0;
    float minAlpha = 999999999999.f;

    for (uint iSubstep = 0; iSubstep < _substepsDeformation; iSubstep++)
    {
        for (uint i = 0; i < nbBasisFlows; i++) {
            BasisFlow& bi = basisFlowParamsPointer[i];
            bi.newCoeff = 0;
        }

        for (uint i = 0; i < nbBasisFlows; i++) {

            BasisFlow& bi = basisFlowParamsPointer[i];

            if (!AllBitsSet(bi.bitFlags, INTERIOR) && !AllBitsSet(bi.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) { continue; }

            float alpha = _dt * _explicitTransferSpeed * powf(WavenumberBasis(bi), -_explicitTransferExponent) / _substepsDeformation; // not exact substepping, but that's what we do by multiplying by dt anyways, so close enough I guess.
            maxAlpha = glm::max<float>(maxAlpha, alpha); // TODO: remove this if not printing
            minAlpha = glm::min<float>(minAlpha, alpha);


            for (uint iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++)
            {
                vector<CoeffBBDecompressedIntersectionInfo>* localIntersectingBasesIdsDeformation = _intersectingBasesIdsDeformation[iRelFreq]->getCpuData_noRefresh(i);

                for (CoeffBBDecompressedIntersectionInfo& inter : *localIntersectingBasesIdsDeformation) {
                    BasisFlow& bj = basisFlowParamsPointer[inter.j];

                    float alphaBiCoeff;
                    alphaBiCoeff = alpha * (bi.coeff + _obstacleBoundaryFactorTransferOnly * bi.coeffBoundary);

                    bj.newCoeff += alphaBiCoeff * transferCoeffs[iRelFreq] * inter.coeff / _coeffBBExplicitTransferSum_abs[i].coeffs[iRelFreq];

                    bi.newCoeff -= alpha * bi.coeff * transferCoeffs[iRelFreq] * abs(inter.coeff) / _coeffBBExplicitTransferSum_abs[i].coeffs[iRelFreq];
                }
            }

        }

        // set .newCoeff as .coeff
        for (uint i = 0; i < nbBasisFlows; i++) {
            BasisFlow& bi = basisFlowParamsPointer[i];

            if (!AllBitsSet(bi.bitFlags, INTERIOR) && !AllBitsSet(bi.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) { continue; }

            bi.coeff += _factorDeformation * bi.newCoeff;
        }

    }





}





