
#include "Application.h"

#include <GLFW\glfw3.h>

#include <iostream>
#include <string>

using namespace std;

void Application::CallbackWindowClosed(GLFWwindow* /*pGlfwWindow*/)
{
    app->_readyToQuit = true;
}

inline string TrueFalseMessage(bool b) {
    return (b ? "TRUE" : "FALSE");
}

void Application::CallbackKey(GLFWwindow* /*pGlfwWindow*/, int key, int /*scancode*/, int action, int /*mods*/)
{
    if(action != GLFW_PRESS) return;

    switch(key) {
    case 'S':
        app->_seedParticles = !app->_seedParticles;
        cout << "Seed particles = " << TrueFalseMessage(app->_seedParticles) << endl;
        break;
    case 'V':
        app->_showVelocity = !app->_showVelocity;
        cout << "Show velocity = " << TrueFalseMessage(app->_showVelocity) << endl;
        break;
    case 'F':
        app->_useForcesFromParticles = !app->_useForcesFromParticles;
        cout << "Use particle buoyancy = " << TrueFalseMessage(app->_useForcesFromParticles) << endl;
        break;
    case 'P':
        app->_drawParticles = !app->_drawParticles;
        cout << "Draw particles = " << TrueFalseMessage(app->_drawParticles) << endl;
        break;
    case 'G':
        app->_drawGrid = !app->_drawGrid;
        cout << "Draw grid = " << TrueFalseMessage(app->_drawGrid) << endl;
        break;
    case 'O':
        app->_drawObstacles = !app->_drawObstacles;
        cout << "Draw obstacles = " << TrueFalseMessage(app->_drawObstacles) << endl;
        break;
    case 'M':
        app->_moveObstacles = !app->_moveObstacles;
        cout << "Move obstacles = " << TrueFalseMessage(app->_moveObstacles) << endl;
        break;
    //case '0':
    //    resetSimulation();
    //    break;
    case GLFW_KEY_SPACE:
        app->_stepSimulation = !app->_stepSimulation;
        cout << "Step simulation = " << TrueFalseMessage(app->_stepSimulation) << endl;
        break;
    }
}


void Application::CallbackMouseScroll(GLFWwindow* /*pGlfwWindow*/, double /*xOffset*/, double yOffset)
{
    if(yOffset > 0) {
        app->_velocityArrowFactor *= 1.2f;
    } else {
        app->_velocityArrowFactor /= 1.2f;
    }
    cout << "Velocity arrow length = " << app->_velocityArrowFactor << "." << endl;
    //pipelineArrows->in_arrowLengthFactor.receive(factor);
}


void Application::DebugCallback(GLenum source, GLenum /*type*/, GLuint /*id*/, GLenum severity, string message)
{
    if(severity == GL_DEBUG_SEVERITY_HIGH &&
       source != GL_DEBUG_SOURCE_APPLICATION){
        cout << message << endl;
    }
}