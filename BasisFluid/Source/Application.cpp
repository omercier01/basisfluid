
#include "Application.h"

#include <iostream>
#include "glm/glm.hpp"

using namespace std;

bool Application::Run() {

    if(!app->Init()) {
        system("pause");
        return 0;
    }

    while(!app->readyToQuit) {
        glfwPollEvents();
    }
}


bool Application::Init() {

    // Initialise GLFW
    if( !glfwInit() )
    {
        cerr << "Failed to initialize GLFW.\n";
        return false;
    } else {
        clog << "Initialized GLFW." << endl;
    }

    // Initialize GLFW
    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* glfwWindow = 
        glfwCreateWindow(1000, 1000, "Local Bases for Model-reduced Smoke Simulations", 0, 0);
    glfwMakeContextCurrent(glfwWindow);
    if( glfwWindow == NULL )
    {
        cerr << "Failed to open GLFW window.\n";
        glfwTerminate();
        return false;
    }

    glfwSetWindowCloseCallback( glfwWindow, Application::CallbackWindowClosed );
    glfwSetKeyCallback( glfwWindow, Application::CallbackKey );
    glfwSetScrollCallback( glfwWindow, Application::CallbackMouseScroll );

    // Initialize GLEW
    //GLEWContext* glewContext = new GLEWContext();
    glewExperimental = GL_TRUE; // Needed for core profile
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        cerr << "Failed to initialize GLEW\n" << glewGetErrorString(err) << endl;
        return false;
    } else {
        clog << "Initialized GLEW." << endl;
    }

    return true;
}