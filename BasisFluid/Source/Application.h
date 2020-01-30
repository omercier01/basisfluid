
#ifndef APPLICATION_H
#define APPLICATION_H

#include "VectorField2D.h"
#include "GridData2D.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <string>
#include <memory>
#include <vector>

#define BASE_GRID_SIZE 64
#define ACCEL_GRID_SIZE 64
#define INTEGRAL_GRID_SIZE 32

#define INVERSION_OPENMP 1
#define USE_DECOMPRESSED_COEFFICIENTS   1
#define INVERSION_INNER_DOUBLE_PRECISION 1
#define INVERSION_STORAGE_DOUBLE_PRECISION 1
#define DEF_COEFF_COMPUTE_GPU 1 // compute coefficients (A, BB, T) and forces on GPU
#define EXPLICIT_ENERGY_TRANSFER 1
#define EXPLICIT_TRANSPORT_ROTATION 1
#define USE_PRECISE_BASIS_EVAL 0
#define INTEGRATE_BASIS_ONE_BASIS_PER_DISPATCH   0
#define INTEGRATE_BASIS_ONE_BASIS_PER_INVOCATION 1
#define SAFETY_ASSERTS 0   // useful for debugging, otherwise turn off for performance


#if INVERSION_INNER_DOUBLE_PRECISION
    typedef double scalar_inversion_inner;
#else
    typedef float scalar_inversion_inner;
#endif

#if INVERSION_STORAGE_DOUBLE_PRECISION
    typedef double scalar_inversion_storage;
#else
    typedef float scalar_inversion_storage;
#endif

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

    const glm::vec2 _domainCenter = { 0,0 };
    const glm::vec2 _domainHalfSize = { 1.f,1.f };
    const float _domainLeft = _domainCenter.x - _domainHalfSize.x;
    const float _domainRight = _domainCenter.x + _domainHalfSize.x;
    const float _domainBottom = _domainCenter.y - _domainHalfSize.y;
    const float _domainTop = _domainCenter.y + _domainHalfSize.y;

    const int nbCells = BASE_GRID_SIZE - 1; // 64 corner data points, so 63 cells
    const int _nbCellsXTotal = nbCells;
    const int _nbCellsXBasis = nbCells;
    const int _nbCellsYTotal = nbCells;
    const int _nbCellsYBasis = nbCells;

    uint obstacleDisplayRes = BASE_GRID_SIZE;
    uint integralGridRes = INTEGRAL_GRID_SIZE - 1;//32-1; // e.g. 32 grid points, 31 cells
    uint accelBasisRes = ACCEL_GRID_SIZE;//32; // raw data, so 32 cells
    uint accelParticlesRes = accelBasisRes;//32; // stored at cell center, so 32 grid points == 32 cells
    uint forcesGridRes = BASE_GRID_SIZE - 1;//32-1; // 32 grid points, 31 cells

    const int _minFreqLvl = 1;
    const int _maxFreqLvl = 1;
    const int _minAnisoLvl = 0;
    const int _maxAnisoLvl = 1; //  MAXIMUM 2, OTHER BASES ARE NOT DEFINED
    //const int _maxFreqLvlFileToUse =  maxFreqLvl;
    //const int _maxAnisoLvlFileToUse = maxAnisoLvl;
    //const float _diffusionRatioLastFrequency = 0.1f*float(1<<maxFreqLvl);
    //const float _lengthLvl0 = 1.0f; // WARNING: not used everywhere, better to set it to 1 for now.


public:

    GLFWwindow* _glfwWindow;

    static const unsigned int _nbExplicitTransferFreqs = 6;
    const glm::ivec2 _explicitTransferFreqs[_nbExplicitTransferFreqs] = { {1,0}, {0,1}, {1,1}, {-1,0}, {0,-1}, {-1,-1} };

    // simulation buffers
    std::unique_ptr<VectorField2D> _velocityField = nullptr;
    std::unique_ptr<VectorField2D> _prevVelocityField = nullptr;
    std::unique_ptr<DataBuffer2D<unsigned int>> _nbParticlesPerCell = nullptr;
    std::unique_ptr<VectorField2D> _forceField = nullptr;
    std::unique_ptr<VectorField2D>* _basisFlowTemplates = nullptr;

    // particles buffers
    std::unique_ptr<DataBuffer1D<vec2>> _partPos = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _partVecs = nullptr;
    std::unique_ptr<DataBuffer1D<float>> _partAges = nullptr;
    //std::unique_ptr<GridData2D<std::vector<unsigned int>*>> _accelParticles = nullptr;
    std::unique_ptr<GridData2D<std::vector<unsigned int>>> _accelParticles = nullptr;

    // draw buffers
    std::unique_ptr<DataBuffer1D<vec2>> _bufferGridPoints = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _bufferArrows = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _obstacleLines = nullptr;

    // force projection buffers
    std::unique_ptr<DataBuffer1D<scalar_inversion_storage>> _vecX = nullptr;
    std::unique_ptr<DataBuffer1D<scalar_inversion_storage>> _vecTemp = nullptr;
    std::unique_ptr<DataBuffer1D<scalar_inversion_storage>> _vecXForces = nullptr;
    std::unique_ptr<DataBuffer1D<scalar_inversion_storage>> _vecXBoundaryForces = nullptr;
    std::unique_ptr<DataBuffer1D<scalar_inversion_storage>> _vecB = nullptr;

    // integration buffers
    //DataBuffer2D<std::vector<unsigned int>*> _accelBasisCentersIds = nullptr;
    std::unique_ptr<DataBuffer2D<std::vector<unsigned int>>> _accelBasisCentersIds = nullptr;
    std::unique_ptr<DataBuffer2D<vec4>> _integrationGridGpu = nullptr;
    std::unique_ptr<DataBuffer1D<vec4>> _integrationTransferBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<float>> _integrationMultipleTransferBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _integrationBasisCentersBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<vector<unsigned int>>> _intersectingBasesIds = nullptr;
    std::unique_ptr<DataBuffer1D<vector<unsigned int>>> _intersectingBasesSignificantBBIds = nullptr;
    std::unique_ptr<DataBuffer1D<vector<unsigned int>>> _intersectingBasesIdsTransport = nullptr;
    std::unique_ptr<DataBuffer1D<std::vector<CoeffBBDecompressedIntersectionInfo>>>
        _intersectingBasesIdsDeformation[_nbExplicitTransferFreqs] = nullptr;


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


