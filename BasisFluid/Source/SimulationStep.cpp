
#include "Application.h"
#include "Obstacles.h"

void Application::SimulationStep()
{
    if (_stepSimulation) {
        ++_frameCount;
    }

    if (_stepSimulation) {
        // update dynamic obstacles
        for (Obstacle* obs : _obstacles)
        {
            if (obs->dynamic)
            {
                obs->prevPhi = obs->phi;
                obs->prevGradPhi = obs->gradPhi;
                obs->updatePhi();

                _obstacleDisplayNeedsUpdating = true;
                _basisStretchedUpdateRequired = true;
            }
        }
    }

    // set prevBitFlags
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        b.prevBitFlags = b.bitFlags;
        b.coeffBoundary = 0;
    }



    // compute stretch ratio, "valid" bool, and stretchedCorners
    if (_basisStretchedUpdateRequired)
    {

        // Reset flags
        BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();
        for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
            basisFlowParamsPointer[i].bitFlags = 0;
            //basisFlowParamsPointer[i].coeffBoundary = 0;
        }

        ComputeStretches();
        _basisStretchedUpdateRequired = false;

        _basisFlowParams->_sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;

        // save stretch flags
        for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
            basisFlowParamsPointer[i].stretchBitFlags = basisFlowParamsPointer[i].bitFlags;
        }
    }

    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].bitFlags = basisFlowParamsPointer[i].stretchBitFlags;
    }




    if (_stepSimulation && _seedParticles)
    {
        SetParticlesInAccelGrid();
    }



    // project dynamic obstacle boundary on bases (reusing forces grid)
    // TODO create a reference to the forces grid instead
    float boundaryPhiBandDecrease = 99999;
    float obstacleBoundaryFactor = 1.5f;
    float obstacleBoundaryFactorTransferOnly = 1.5f;


    // clear forces
    _forceField->populateWithFunction([=](float /*x*/, float /*y*/) { return vec2(0); });

    for (Obstacle* obs : _obstacles)
    {
        if (obs->dynamic)
        {
            _forceField->addFunction([=](float x, float y) {
                float phi = obs->phi(vec2(x, y));
                return
                    -(obstacleBoundaryFactor * (phi - obs->prevPhi(vec2(x, y))) / _dt * obs->gradPhi(vec2(x, y))) * glm::max<float>(1.f - abs(phi) / boundaryPhiBandDecrease, 0.f)
                    ;
        });
    }
}
    _forceField->_vectors._sourceStorageType = DataBuffer2D<vec2>::StorageType::CPU;
    //_forceField->mVectors.dirtyData();


    //BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();
    basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    scalar_inversion_storage* vecBPointer = _vecB->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis)
    {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        vecBPointer[iBasis] = IntegrateBasisGrid(b, _forceField.get());
    }

    //inverseBBMatrix(vecX,vecB,inversionPrecision);
    InverseBBMatrix(_vecXBoundaryForces.get(), _vecB.get(), 0, BASIS_FLAGS::DYNAMIC_BOUNDARY_PROJECTION);

    scalar_inversion_storage* vecXBoundaryForcesPointer = _vecXBoundaryForces->getCpuDataPointer();
    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].coeffBoundary = float(vecXBoundaryForcesPointer[i]);
    }



    ComputeBasisAdvection();


    AddParticleForcesToBasisFlows();


    ComputeParticleAdvection();


    if (_seedParticles)
    {
        SeedParticles();
    }


}