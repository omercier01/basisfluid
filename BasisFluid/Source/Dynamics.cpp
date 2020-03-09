
#include "Application.h"

#include "glm/ext.hpp"

using namespace std;

// add forces from particles (buoyancy) using a grid. Clumps particles together in the grid, then integrates against the grid to obtain basis coefficient.
void Application::AddParticleForcesToBasisFlows()
{

    // reset forces to zero
    _forceField->populateWithFunction([](float, float) {return vec2(0); });


    // add particle forces in forces grid
    vec2* particlesPointer = _partPos->getCpuDataPointer();
    float* agesPointer = _partAges->getCpuDataPointer();


    float cellSizeX = (_forceField->_boundXMax-_forceField->_boundXMin)/_forceField->_nbCellsX;
    float cellSizeY = (_forceField->_boundYMax-_forceField->_boundYMin)/_forceField->_nbCellsY;
    for (unsigned int i = 0; i < _partPos->_nbElements; i++) {
        uvec2 indexCell = _forceField->pointToCellIndex(particlesPointer[i]);
        vec2 lowerCellCornerPosition = _forceField->indexToPosition(indexCell);
        vec2 increment = _dt * powf(_buoyancyDecayRatioWithAge, agesPointer[i]) * _buoyancyPerParticle*vec2(0, 1);
        for(int iX = 0; iX <= 1; iX++) {
            for(int iY = 0; iY <= 1; iY++) {
                _forceField->addVectorCpuData(
                    indexCell.x+iX, indexCell.y+iY,
                    increment * (1-abs(lowerCellCornerPosition.x-particlesPointer[i].x)/cellSizeX) * ((1-abs(lowerCellCornerPosition.y-particlesPointer[i].y)/cellSizeY)));
            }
        }
    }

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

    // project forces into basis space  
    double* vecBPointer = _vecB->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        if (AllBitsSet(basisFlowParamsPointer[iBasis].bitFlags, BASIS_FLAGS::FORCE_PROJECTION)) {
            vecBPointer[iBasis] = IntegrateBasisGrid(basisFlowParamsPointer[iBasis], _forceField.get());
        }
        else {
            vecBPointer[iBasis] = 0;
        }
    }


    // inverse to obtain base weights
    InverseBBMatrix(_vecXForces.get(), _vecB.get(), BASIS_FLAGS::FORCE_PROJECTION);


    // add force weights to current basis weights
    double* vecXForcesPointer = _vecXForces->getCpuDataPointer();
    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].coeff += _dt * float(vecXForcesPointer[i]);
    }

}



void Application::ComputeNewCenterProportions(vec2& newCenter, BasisFlow& bi, BasisFlow& bj, vec2& interBasisDist)
{
    float tempCoeff = glm::max<float>(0.f, (1.f - abs(newCenter.x - bj.center.x) / interBasisDist.x)) *
        glm::max<float>(0.f, (1.f - abs(newCenter.y - bj.center.y) / interBasisDist.y));
    bj.newCoeff += bi.coeff * tempCoeff;
}



void Application::ComputeBasisAdvection()
{
    const unsigned int nbBasisFlows = _basisFlowParams->_nbElements;

    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    for (uint i = 0; i < nbBasisFlows; i++) {
        BasisFlow& bi = basisFlowParamsPointer[i];
        bi.newCoeff = 0;
    }

    for (uint i = 0; i < nbBasisFlows; i++) {

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

        vec2 freqI(1 << bi.freqLvl.x, 1 << bi.freqLvl.y);

        vec2 interBasisDist = 0.5f * 0.5f / freqI;


        if (
            abs(newCenter.x - bi.center.x) > bi.supportHalfSize().x*0.5f*0.99 ||
            abs(newCenter.y - bi.center.y) > bi.supportHalfSize().y*0.5f*0.99
            ) {
            // new center not located within immediate neighbours, must look further using basis acceleration grid
            uvec2 minId(9999);
            uvec2 maxId(0);
            for (int iX = -1; iX <= 1; iX += 2) {
                for (int iY = -1; iY <= 1; iY += 2) {
                    vec2 corner = newCenter + vec2(iX, iY)*bi.supportHalfSize()*0.5f;
                    uvec2 cellId;
                    cellId.x = glm::clamp<int>(int(floor((corner.x - _domainLeft) / (_domainRight - _domainLeft)*_accelBasisRes)), 0, _accelBasisRes - 1);
                    cellId.y = glm::clamp<int>(int(floor((corner.y - _domainBottom) / (_domainTop - _domainBottom)*_accelBasisRes)), 0, _accelBasisRes - 1);
                    minId = glm::min<uint>(cellId, minId);
                    maxId = glm::max<uint>(cellId, maxId);
                }
            }

            for (uint iX = minId.x; iX <= maxId.x; iX++) {
                for (uint iY = minId.y; iY <= maxId.y; iY++) {
                    {
                        for (unsigned int basisId : *(_accelBasisCentersIds->getCpuData(iX, iY)))
                        {
                            BasisFlow& bj = basisFlowParamsPointer[basisId];
                            if (bj.freqLvl == bi.freqLvl) {
                                ComputeNewCenterProportions(newCenter, bi, bj, interBasisDist);
                            }
                        }
                    }
                }
            }

        }
        else
        {

            // new center within immediate neighbours
            vector<unsigned int>* localIntersectingBasesIdsTransport = _intersectingBasesIdsTransport->getCpuData(i);
            for (auto itJ = localIntersectingBasesIdsTransport->begin(); itJ != localIntersectingBasesIdsTransport->end(); ++itJ) {
                BasisFlow& bj = basisFlowParamsPointer[*itJ];
                ComputeNewCenterProportions(newCenter, bi, bj, interBasisDist);
            }

        }

    }


    // add transport cofficients and set new coefficients
    for (uint i = 0; i < nbBasisFlows; i++) {
        BasisFlow& bi = basisFlowParamsPointer[i];

        if (!AllBitsSet(bi.bitFlags, INTERIOR) && !AllBitsSet(bi.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) { continue; }

        bi.coeff = bi.newCoeff;
        bi.newCoeff = 0;
    }

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

            for (uint iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++)
            {
                vector<CoeffBBDecompressedIntersectionInfo>* localIntersectingBasesIdsDeformation = _intersectingBasesIdsDeformation[iRelFreq]->getCpuData(i);

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

            bi.coeff += bi.newCoeff;
        }

    }



}





