
#include "Application.h"

#include <iostream>

using namespace std;

bool Application::Run() {

    if (!app->Init()) {
        system("pause");
        return 0;
    }

    glfwSwapInterval(0);
    while (!app->_readyToQuit) {
        glfwPollEvents();
        if (_stepSimulation) {
            SimulationStep();
        }
        Draw();
    }

    return 0;
}