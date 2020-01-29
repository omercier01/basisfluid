
#ifndef APPLICATION_H
#define APPLICATION_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <string>

// Global class to manage program execution
class Application {
public:
    bool Run();
    bool Init();
    bool Init_DataBuffers();
    bool Init_Obstacles();
    bool Init_BasisFlows();
    
    bool Init_ComputeShaders();
    bool Init_CSIntegrateAvgBasis();
    bool Init_CSIntegrateBasisBasis();
    bool Init_CSIntegrateBasisGradBasis();
    bool Init_CSIntegrateBasisGrid_onePerDispatch();
    bool Init_CSIntegrateBasisGrid_onePerInvocation();

    void SimulationStep();
    void Draw();

public:
    static void CallbackWindowClosed(GLFWwindow* pGlfwWindow);
    static void CallbackKey(GLFWwindow* pGlfwWindow, int key, int scancode, int action, int mods);
    static void CallbackMouseScroll(GLFWwindow* pGlfwWindow, double xOffset, double yOffset);

public:
    GLFWwindow* _glfwWindow ;

    const glm::ivec2 _explicitTransferFreqs[6] = {{1,0}, {0,1}, {1,1}, {-1,0}, {0,-1}, {-1,-1}};

    bool _readyToQuit = false;
    bool _seedParticles = true;
    bool _showVelocityGrid = false;
    bool _useForcesFromParticles = true;
    bool _drawParticles = true;
    bool _drawGrid = true;
    bool _drawObstacles = true;
    bool _stepSimulation = true;
    float _velocityArrowFactor = 0.1f;
};

extern Application* app;

#endif /*APPLICATION_H*/


