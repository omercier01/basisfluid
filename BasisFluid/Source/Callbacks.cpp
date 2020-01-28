
#include "Application.h"
#include <GLFW\glfw3.h>
#include <iostream>

using namespace std;

void Application::CallbackWindowClosed(GLFWwindow* /*pGlfwWindow*/)
{
    app->readyToQuit = true;
}

inline string TrueFalseMessage(bool b) {
    return b  ?"TRUE" : "FALSE";
}

void Application::CallbackKey(GLFWwindow* /*pGlfwWindow*/, int key, int /*scancode*/, int action, int /*mods*/)
{
    if(action != GLFW_PRESS) return;

    switch(key) {
    case 'S':
        app->seedParticles = !app->seedParticles;
        cout << "Seed particles = " << TrueFalseMessage(app->seedParticles) << endl;
        break;
    case 'V':
        app->showVelocityGrid = !app->showVelocityGrid;
        cout << "Show velocity = " << TrueFalseMessage(app->showVelocityGrid) << endl;
        break;
    case 'F':
        app->useForcesFromParticles = !app->useForcesFromParticles;
        cout << "Use particle buoyancy = " << TrueFalseMessage(app->useForcesFromParticles) << endl;
        break;
    case 'P':
        app->drawParticles = !app->drawParticles;
        cout << "Draw particles = " << TrueFalseMessage(app->drawParticles) << endl;
        break;
    case 'G':
        app->drawGrid = !app->drawGrid;
        cout << "Draw grid = " << TrueFalseMessage(app->drawGrid) << endl;
        break;
    case 'O':
        app->drawObstacles = !app->drawObstacles;
        cout << "Draw obstacles = " << TrueFalseMessage(app->drawObstacles) << endl;
        break;
    //case '0':
    //    resetSimulation();
    //    break;
    case GLFW_KEY_SPACE:
        app->stepSimulation = !app->stepSimulation;
        cout << "Step simulation = " << TrueFalseMessage(app->stepSimulation) << endl;
        break;
    }
}


void Application::CallbackMouseScroll(GLFWwindow* /*pGlfwWindow*/, double /*xOffset*/, double yOffset)
{
    static float factor = 0.1f;
    if(yOffset > 0) {
        factor *= 1.2f;
    } else {
        factor /= 1.2f;
    }
    cout << "Arrow display factor =" << factor << "." << endl;
    //pipelineArrows->in_arrowLengthFactor.receive(factor);
}
