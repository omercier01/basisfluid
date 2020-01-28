
#ifndef APPLICATION_H
#define APPLICATION_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <string>

// Global class to manage program execution
class Application {
public:
    bool Run();
    bool Init();

public:
    static void CallbackWindowClosed(GLFWwindow* pGlfwWindow);
    static void CallbackKey(GLFWwindow* pGlfwWindow, int key, int scancode, int action, int mods);
    static void CallbackMouseScroll(GLFWwindow* pGlfwWindow, double xOffset, double yOffset);

public:
    bool readyToQuit = false;
    bool seedParticles = true;
    bool showVelocityGrid = false;
    bool useForcesFromParticles = true;
    bool drawParticles = true;
    bool drawGrid = true;
    bool drawObstacles = true;
    bool stepSimulation = true;
};

extern Application* app;

#endif /*APPLICATION_H*/


