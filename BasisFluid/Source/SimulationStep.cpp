
#include "Application.h"

void Application::SimulationStep()
{
    if(_stepSimulation) {
        ++_frameCount;
    }
}