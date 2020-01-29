
#include "Application.h"

#include "glm/glm.hpp"
#include <GL/glew.h>

#include <algorithm>
#include <iostream>
#include <time.h>

using namespace std;
using namespace glm;

bool Application::Init() {

    // Init GLFW
    if( !glfwInit() )
    {
        cerr << "Failed to initialize GLFW.\n";
        return false;
    } else {
        clog << "Initialized GLFW." << endl;
    }

    // Init Window
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    _glfwWindow = glfwCreateWindow(1000, 1000, "Local Bases for Model-reduced Smoke Simulations", 0, 0);
    glfwMakeContextCurrent(_glfwWindow);
    if( _glfwWindow == NULL )
    {
        cerr << "Failed to open GLFW window.\n";
        glfwTerminate();
        return false;
    }

    glfwSetWindowCloseCallback( _glfwWindow, Application::CallbackWindowClosed );
    glfwSetKeyCallback( _glfwWindow, Application::CallbackKey );
    glfwSetScrollCallback( _glfwWindow, Application::CallbackMouseScroll );

    // Init GLEW
    glewExperimental = GL_TRUE; // Needed for core profile
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        cerr << "Failed to initialize GLEW\n" << glewGetErrorString(err) << endl;
        return false;
    } else {
        clog << "Initialized GLEW." << endl;
    }

    // OpenGL
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   
    // Random
    srand((unsigned int)(time(0)));

    //
    // Init basis flows
    //  

    if(!Init_DataBuffers()) {
        cerr << "Error initializing data buffers." << endl;
        return false;
    }
    if(!Init_Obstacles()) {
        cerr << "Error initializing obstacles." << endl;
        return false;
    }
    if(!Init_ComputeShaders()) {
        cerr << "Error initializing compute shaders." << endl;
        return false;
    }
    if(!Init_BasisFlows()) {
        cerr << "Error initializing basis flows." << endl;
        return false;
    }

    return true;
}