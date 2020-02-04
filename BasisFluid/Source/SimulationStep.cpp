
#include "Application.h"
#include "Obstacles.h"

void Application::SimulationStep()
{
    if(_stepSimulation) {
        ++_frameCount;
    }

    if( _stepSimulation ) {
        // update dynamic obstacles
        for(Obstacle* obs : _obstacles)
        {
            if(obs->dynamic)
            {
                obs->prevPhi = obs->phi;
                obs->prevGradPhi = obs->gradPhi;
                obs->updatePhi();
                
                _obstacleDisplayNeedsUpdating = true;
                _basisStretchedUpdateRequired = true;
            }
        }
    }
}