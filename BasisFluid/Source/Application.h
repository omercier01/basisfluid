
#ifndef APPLICATION_H
#define APPLICATION_H

#include "VectorField2D.h"
#include "GridData2D.h"
#include "BasisFlows.h"

#define GLM_FORCE_RADIANS
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <string>
#include <memory>
#include <vector>

class Obstacle;
class ObstacleShaderPipeline;
class ParticleShaderPipeline;
class VelocityArrowShaderPipeline;

enum class ObstacleType { None, Circle, Bar };

// Global class to manage program execution
class Application {

public:

    //================================================
    // Parameters

    // Frequency levels
    // Note: Pick maximum basis flows sizes (minimum frequency levels) that are not too large
    // compared to domain features, otherwise basis stretch can sometimes fail
    const int _minFreqLvl = 1;
    const int _maxFreqLvl = 2;
    const int _minAnisoLvl = 0;
    const int _maxAnisoLvl = 1; // MAXIMUM 2, OTHER BASES ARE NOT DEFINED

    // Simulation parameters

    // nb iterations to project forces
    const uint _maxNbItMatBBInversion = 10;

    // nb iterations when computing basis stretches
    const uint _nbStretchLoops = 2;

    // nb iterations to compute stretched coordinated (paper Section 6.1)
    const unsigned int _nbNewtonInversionIterations = 3;

    // time step
    const float _dt = 0.0325f;

    // buoyancy per particle
    const float _buoyancyPerParticle = 0.1f;
    
    // buoyancy reduction exponent with time (usually in [0,1])
    const float _buoyancyDecayRatioWithAge = 1.f;

    // simulation viualization viewpoint
    const glm::mat4 _viewProjMat = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 };

    // nb steps to advect particles
    const uint _substepsParticles = 1;

    // particle seeding region
    const float _seedCenterX = 0.f;
    const float _seedCenterY = -0.75f;
    const float _seedRadius = 0.1f;

    // nb substeps when evolving bases
    const uint _substepsDeformation = 1;

    // lambda and epsilon for energy cascade (paper Section 5.3)
    /*const float _explicitTransferSpeed = 0.1f;
    const float _explicitTransferExponent = -1.66f;
    const float _factorDeformation = 0.5f;*/
    const float _explicitTransferSpeed = 0.05f;
    const float _explicitTransferExponent = -1.66f;

    // multiplicator for projected flow around moving obstacles
    const float _obstacleBoundaryFactorTransferOnly = 1.5f;


    // the region of allowed basis corner movement has width _stretchBandRatio times the basis support half size, 
    // i.e. the basis can be stretched or sqquished by at most this ratio, otherwise it is discarded
    const float _stretchBandRatio = 0.5f; 

    // particle life time, in nb of frames
    const unsigned int _particleLifeTime = 300;

    // nb of particles seeded each frame
    const unsigned int _nbParticlesPerSeedGroupPerDimension = 200;

    // transfer ratio to each neighboring frequency layers (paper Figure 6)
    // Note: later normalized so total of ratios is 1
    const float _explicitTransfer_10 = 1.f;
    const float _explicitTransfer_01 = 1.f;
    const float _explicitTransfer_11 = 1.f;
    const float _explicitTransfer_m10 = 0.f;
    const float _explicitTransfer_0m1 = 0.f;
    const float _explicitTransfer_m1m1 = 0.f;

    // obstacle parameters

    //ObstacleType _obstacleType = ObstacleType::None;
    ObstacleType _obstacleType = ObstacleType::Circle;
    //ObstacleType _obstacleType = ObstacleType::Bar;

    const float _obstacleCircleRadius = 0.2f;
    const float _obstacleCircleMotionSpeed = 1.f;
    const float _obstacleCircleMotionRadius = 0.25f;

    const float _obstacleBarWidth = 0.1f;
    const float _obstacleBarHeight = 0.2f;
    const float _obstacleBarRotationSpeed = 0.789f;
    const float _obstacleBarMotionSpeed = 0.75f;
    const float _obstacleBarMotionAmplitude = 0.25f;

    // simulation domain
    const glm::vec2 _domainCenter = { 0,0 };
    const glm::vec2 _domainHalfSize = { 1.f,1.f };

    // simulation grid sizes

    // velocity grid for visualization
    const unsigned int _nbCellsVelocity = 64 - 1;

    // grid to store basis templates
    const unsigned int _nbCellsBasisTemplates = 64 - 1;

    // obstacle marching squares
    const unsigned int _obstacleDisplayRes = 256;

    // grid to integrate basis when computing coefficients
    const unsigned int _integralGridRes = 32 - 1; 

    // acceleration structure for basis centers
    const unsigned int _accelBasisRes = 32;

    // particle acceleration structure
    const unsigned int _accelParticlesRes = 64;

    // grid to compute forces
    const unsigned int _forcesGridRes = 64 - 1;

    // Width of support of first basis frequency level.
    // Note: not tested with _lengthLvl0 != 1
    const float _lengthLvl0 = 1.0f;

    // the vector growth factor using mouse scroll
    float _velocityArrowFactor = 0.1f;

    // initial application controls
    bool _seedParticles = true;
    bool _showVelocity = false;
    bool _useForcesFromParticles = true;
    bool _drawParticles = true;
    bool _moveObstacles = true;
    bool _stepSimulation = true;

    //================================================

    // Size of grid to snap basis centers when saving coefficient dictionary to avoid float errors
    const float _coeffSnapSize = _lengthLvl0 / float(1 << _maxFreqLvl) / 32.0f;

    const float _domainLeft = _domainCenter.x - _domainHalfSize.x;
    const float _domainRight = _domainCenter.x + _domainHalfSize.x;
    const float _domainBottom = _domainCenter.y - _domainHalfSize.y;
    const float _domainTop = _domainCenter.y + _domainHalfSize.y;

public:
    Application();
    ~Application();

    bool Run();
    bool Init();
    bool Init_DataBuffers();
    bool Init_Obstacles();
    bool Init_BasisFlows();

    bool Init_Shaders();

    void SimulationStep();
    void Draw();
    void ComputeVelocityGridForDisplay();

    float MatBBCoeff(int i, int j);
    float MatBBCoeff(const BasisFlow& b1, const BasisFlow& b2);

    glm::vec2 MatTCoeff(int i, int j);
    glm::vec2 MatTCoeff(BasisFlow bTransported, BasisFlow bTransporting);

    void InverseBBMatrix(
        DataBuffer1D<double>* vecX,
        DataBuffer1D<double>* vecB,
        unsigned int basisBitMask);
    void InverseBBMatrixMain(
        unsigned int iRow, double* vecX, double* vecB,
        BasisFlow* basisDataPointer, unsigned int basisBitMask);

    void SaveCoeffsBB(std::string filename);
    void LoadCoeffsBB(std::string filename);
    void SaveCoeffsT(std::string filename);
    void LoadCoeffsT(std::string filename);

    float IntegrateBasisBasis(BasisFlow b1, BasisFlow b2);
    float IntegrateBasisGrid(BasisFlow& b, VectorField2D* velField);
    glm::vec2 TranslatedBasisEval(const glm::vec2 p, const glm::ivec2 freqLvl, const glm::vec2 center);
    glm::vec2 AverageBasisOnSupport(BasisFlow bVec, BasisFlow bSupport);

    BasisFlow ComputeStretch(BasisFlow b, bool staticObstaclesOnly = false, bool noStretch = false);
    void ComputeStretches();
    vec2 QuadCoord(vec2 p, BasisFlow const& b);
    mat2 QuadCoordInvDeriv(vec2 uv, BasisFlow const& b);
    vec2 VecObstacle_stretch(vec2 p, BasisFlow const& b);

    void SetParticlesInAccelGrid();
    void ComputeParticleAdvection();
    void SeedParticles();

    void AddParticleForcesToBasisFlows();
    void ComputeNewCenterProportions(vec2& newCenter, BasisFlow& bi, BasisFlow& bj, vec2& interBasisDist);
    void ComputeBasisAdvection();

    inline float WavenumberBasis(BasisFlow& b)
    {
        return powf(2.f, 0.5f*(b.freqLvl.x + b.freqLvl.y));
    }

public:
    static void CallbackWindowClosed(GLFWwindow* pGlfwWindow);
    static void CallbackKey(GLFWwindow* pGlfwWindow, int key, int scancode, int action, int mods);
    static void CallbackMouseScroll(GLFWwindow* pGlfwWindow, double xOffset, double yOffset);
    static void DebugCallback(GLenum source, GLenum /*type*/, GLuint /*id*/, GLenum severity, std::string message);

public:

    GLFWwindow * _glfwWindow;

    static const unsigned int _nbExplicitTransferFreqs = 6;
    const glm::ivec2 _explicitTransferFreqs[_nbExplicitTransferFreqs] = { {1,0}, {0,1}, {1,1}, {-1,0}, {0,-1}, {-1,-1} };

    // simulation buffers
    std::unique_ptr<VectorField2D> _velocityField = nullptr;
    std::unique_ptr<VectorField2D> _prevVelocityField = nullptr;
    std::unique_ptr<VectorField2D> _forceField = nullptr;
    std::unique_ptr<VectorField2D>* _basisFlowTemplates = nullptr;
    std::unique_ptr<DataBuffer1D<BasisFlow>> _basisFlowParams = nullptr;
    std::vector<ivec2> _freqLvls;

    // particles buffers
    std::unique_ptr<DataBuffer1D<vec2>> _partPos = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _partVecs = nullptr;
    std::unique_ptr<DataBuffer1D<float>> _partAges = nullptr;
    std::unique_ptr<GridData2D<std::vector<unsigned int>*>> _accelParticles = nullptr;

    // draw buffers
    std::unique_ptr<DataBuffer1D<vec2>> _bufferGridPoints = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _bufferArrows = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _obstacleLines = nullptr;

    // force projection buffers
    std::unique_ptr<DataBuffer1D<double>> _vecX = nullptr;
    std::unique_ptr<DataBuffer1D<double>> _vecTemp = nullptr;
    std::unique_ptr<DataBuffer1D<double>> _vecXForces = nullptr;
    std::unique_ptr<DataBuffer1D<double>> _vecXBoundaryForces = nullptr;
    std::unique_ptr<DataBuffer1D<double>> _vecB = nullptr;

    // integration buffers
    std::unique_ptr<DataBuffer2D<std::vector<unsigned int>*>> _accelBasisCentersIds = nullptr;
    std::unique_ptr<DataBuffer1D<vec4>> _integrationTransferBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<float>> _integrationMultipleTransferBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<vec2>> _integrationBasisCentersBufferGpu = nullptr;
    std::unique_ptr<DataBuffer1D<std::vector<unsigned int>*>> _intersectingBasesIds = nullptr;
    std::unique_ptr<DataBuffer1D<std::vector<unsigned int>*>> _intersectingBasesSignificantBBIds = nullptr;
    std::unique_ptr<DataBuffer1D<std::vector<unsigned int>*>> _intersectingBasesIdsTransport = nullptr;
    std::unique_ptr<DataBuffer1D<std::vector<CoeffBBDecompressedIntersectionInfo>*>>
        _intersectingBasesIdsDeformation[_nbExplicitTransferFreqs];

    struct ExplicitTransferCoeffs {
        float coeffs[_nbExplicitTransferFreqs];
    };
    std::vector<ExplicitTransferCoeffs> _coeffBBExplicitTransferSum_abs;
    std::vector<ExplicitTransferCoeffs> _coeffBBExplicitTransferSum_sqr;

    std::vector<Obstacle*> _obstacles;

    std::vector<std::vector<unsigned int>> _orthogonalBasisGroupIds;
    std::vector<std::vector<unsigned int>> _sameBasisTemplateGroupIds;

    MapTypeBB _coeffsBB;
    MapTypeT _coeffsT;
    std::vector<std::vector<CoeffBBDecompressedIntersectionInfo>> _coeffsBBDecompressedIntersections; // Does not include the basis itself.
    std::vector<std::vector<CoeffTDecompressedIntersectionInfo>>  _coeffsTDecompressedIntersections;

    // shader pipelines
    std::unique_ptr<ObstacleShaderPipeline> _pipelineObstacle;
    std::unique_ptr<ParticleShaderPipeline> _pipelineParticle;
    std::unique_ptr<VelocityArrowShaderPipeline> _pipelineVelocityArrow;

    // Status
    bool _readyToQuit = false;
    bool _newBBCoeffComputed = false;
    bool _newACoeffComputed = false;
    bool _newTCoeffComputed = false;
    bool _newRCoeffComputed = false;
    bool _obstacleDisplayNeedsUpdating = true;
    bool _basisStretchedUpdateRequired = true;
    unsigned int _particleCircularSeedId = 0;
    bool _particleSeedBufferLooped = false;

};

extern Application* app;

#endif /*APPLICATION_H*/


