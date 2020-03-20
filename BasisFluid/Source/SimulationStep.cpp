
#include "Application.h"
#include "Obstacles.h"

void Application::SimulationStep()
{
    _velocityGridNeedsUpdating = true;

    // update dynamic obstacles
    for (Obstacle* obs : _obstacles)
    {
        if (obs->dynamic)
        {
            obs->prevPhi = obs->phi;
            obs->prevGradPhi = obs->gradPhi;
            if (_moveObstacles) {
                obs->updatePhi();
            }

            _obstacleDisplayNeedsUpdating = true;
            _basisStretchedUpdateRequired = true;
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
        }

        ComputeStretches();
        _basisStretchedUpdateRequired = false;

        // save stretch flags
        for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
            basisFlowParamsPointer[i].stretchBitFlags = basisFlowParamsPointer[i].bitFlags;
        }
    }

    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].bitFlags = basisFlowParamsPointer[i].stretchBitFlags;
    }

    SetParticlesInAccelGrid();

    // project dynamic obstacle boundary on bases (reusing forces grid)
    _forceField->populateWithFunction([=](float /*x*/, float /*y*/) { return vec2(0); });
    for (Obstacle* obs : _obstacles)
    {
        if (obs->dynamic)
        {
            _forceField->addFunction([=](float x, float y) {
                float phi = obs->phi(vec2(x, y));
                return -(_obstacleBoundaryFactor * (phi - obs->prevPhi(vec2(x, y))) / _dt * obs->gradPhi(vec2(x, y))) * glm::max<float>(1.f - abs(phi) / _boundarySDFBandDecrease, 0.f);
            });
        }
    }

    basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    double* vecBPointer = _vecB->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis)
    {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        vecBPointer[iBasis] = IntegrateBasisGrid(b, _forceField.get());
    }

    InverseBBMatrix(_vecXBoundaryForces.get(), _vecB.get(), BASIS_FLAGS::DYNAMIC_BOUNDARY_PROJECTION);

    double* vecXBoundaryForcesPointer = _vecXBoundaryForces->getCpuDataPointer();
    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        basisFlowParamsPointer[i].coeffBoundary = float(vecXBoundaryForcesPointer[i]);
    }

    ComputeBasisAdvection();

    if (_useForcesFromParticles) {
        AddParticleForcesToBasisFlows();
    }

    ComputeParticleAdvection();

    if (_seedParticles)
    {
        SeedParticles();
    }

}