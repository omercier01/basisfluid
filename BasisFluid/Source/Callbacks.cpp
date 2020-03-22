
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
        // toggle particle seeding
        app->_seedParticles = !app->_seedParticles;
        cout << "Seed particles = " << TrueFalseMessage(app->_seedParticles) << endl;
        break;
    case 'V':
        // toggle draw velocity grid
        app->_showVelocity = !app->_showVelocity;
        cout << "Show velocity = " << TrueFalseMessage(app->_showVelocity) << endl;
        break;
    case 'F':
        // toggle use particle buoyancy
        app->_useForcesFromParticles = !app->_useForcesFromParticles;
        cout << "Use particle buoyancy = " << TrueFalseMessage(app->_useForcesFromParticles) << endl;
        break;
    case 'P':
        // toggle draw particles
        app->_drawParticles = !app->_drawParticles;
        cout << "Draw particles = " << TrueFalseMessage(app->_drawParticles) << endl;
        break;
    case 'M':
        // toggle move obstacles
        app->_moveObstacles = !app->_moveObstacles;
        cout << "Move obstacles = " << TrueFalseMessage(app->_moveObstacles) << endl;
        break;
    case GLFW_KEY_SPACE:
        // toggle step simulation
        app->_stepSimulation = !app->_stepSimulation;
        cout << "Step simulation = " << TrueFalseMessage(app->_stepSimulation) << endl;
        break;
    }
}


void Application::CallbackMouseScroll(GLFWwindow* /*pGlfwWindow*/, double /*xOffset*/, double yOffset)
{
    // stretch arrow length for velocity draw
    if(yOffset > 0) {
        app->_velocityArrowFactor *= 1.2f;
    } else {
        app->_velocityArrowFactor /= 1.2f;
    }
    cout << "Velocity arrow length = " << app->_velocityArrowFactor << "." << endl;
}


void Application::DebugCallback(GLenum source, GLenum /*type*/, GLuint /*id*/, GLenum severity, string message)
{
    if(severity == GL_DEBUG_SEVERITY_HIGH &&
       source != GL_DEBUG_SOURCE_APPLICATION){
        cout << message << endl;
    }
}