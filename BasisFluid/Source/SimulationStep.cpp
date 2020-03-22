
#include "Application.h"
#include "Obstacles.h"

void Application::SimulationStep()
{
    _velocityGridNeedsUpdating = true;

    // update dynamic obstacles
    for (Obstacle* obs : _obstacles) {
        if (obs->dynamic) {
            obs->prevPhi = obs->phi;
            obs->prevGradPhi = obs->gradPhi;
            if (_moveObstacles) {
                obs->updatePhi();
            }

            _obstacleDisplayNeedsUpdating = true;
            _basisStretchedUpdateRequired = true;
        }
    }

    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    // set prevBitFlags
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        b.coeffBoundary = 0;
    }

    // compute stretch ratio, "valid" bool, and stretchedCorners
    if (_basisStretchedUpdateRequired)
    {
        // Reset flags
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

    ProjectDynamicObstacleBoundaryMotion();

    // Advect bases
    ComputeBasisAdvection();

    // project particle buoyancy to basis flows
    if (_useForcesFromParticles) {
        AddParticleForcesToBasisFlows();
    }

    ComputeParticleAdvection();

    if (_seedParticles) {
        SeedParticles();
    }

}