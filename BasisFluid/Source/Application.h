
// Some comments refer to the paper "Local Bases for Model-reduced Smoke Simulations" for more details.

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

    // transfer ratio to each neighboring frequency layers (paper Figure 6). "m" means "minus".
    // coefficient 10 sends energy down the cascade, e.g. from frequency (2,2) to frequency (4,2).
    // coefficient m10 sends energy back up , e.g. from frequency (4,2) to frequency (2,2).
    // These ratios are later normalized to have unit sum.
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

    // Main loop
    bool Run();

    // Initialization calls, return false if failed
    bool Init();
    bool Init_DataBuffers();
    bool Init_Obstacles();
    bool Init_BasisFlows();
    bool Init_Shaders();

    // Main simulation loop
    void SimulationStep();

    // Render loop
    void Draw();

    // Render velocity grid
    void ComputeVelocityGridForDisplay();

    // Computes the B^T.B coefficient, using cached value if coefficient has been previously computed
    // i,j: indices of coefficient to compute
    // b1,b2: basis info of coefficient to compute
    float MatBBCoeff(int i, int j);
    float MatBBCoeff(const BasisFlow& b1, const BasisFlow& b2);

    // Computes the T coefficient, using cached value if coefficient has been previously computed
    // iTransported, bTransported: index or basis info of basis being transported
    // iTransporting, bTransporting: index or basis info of basis acting on the other
    glm::vec2 MatTCoeff(int iTransported, int iTransporting);
    glm::vec2 MatTCoeff(BasisFlow bTransported, BasisFlow bTransporting);

    // Solves B^T.B.vecX = vecB for vecX, only using the basis flows that have all bits basisBitMask
    // turned on. This is used to project forces onto the basis (where boundary basis flows are ignored)
    // or to project a moving obstacle's motion onto the boundary bases (in which cases nly boundary
    // basis flows are used).
    void InverseBBMatrix(
        DataBuffer1D<double>* vecX,
        DataBuffer1D<double>* vecB,
        unsigned int basisBitMask);
    void InverseBBMatrixMain(
        unsigned int iRow, double* vecX, double* vecB,
        BasisFlow* basisDataPointer, unsigned int basisBitMask);

    // Saves/loads the coefficient dictionaries to/from file
    void SaveCoeffsBB(std::string filename);
    void LoadCoeffsBB(std::string filename);
    void SaveCoeffsT(std::string filename);
    void LoadCoeffsT(std::string filename);

    // Evaluates a basis at a given point from its basis template (i.e. scaling and translating the
    // basis template to the right frequency level and center)
    // p: point to evaluate
    // freqLevel: Frequency of the basis
    // center: center of the basis
    glm::vec2 TranslatedBasisEval(const glm::vec2 p, const glm::ivec2 freqLvl, const glm::vec2 center);

    // Computes \int(b1.b2), see Equation 1.
    float IntegrateBasisBasis(BasisFlow b1, BasisFlow b2);

    // Computes \int(b.vecField), where velField is a vector field defined on the simulation domain.
    float IntegrateBasisGrid(BasisFlow& b, VectorField2D* velField);

    // Computes \int_{S}(bVec)/\int_{S}, where S is bSupport's support. See Equation 18.
    glm::vec2 AverageBasisOnSupport(BasisFlow bVec, BasisFlow bSupport);

    // Computes the stretch of a given basis flow, see Section 6.1
    BasisFlow ComputeStretch(BasisFlow b);
    
    // Compute stretches for all basis flows
    void ComputeStretches();

    // Evaluated a stretched basis flow at point p.
    // p: evaluatiom point
    // b: stretched basis
    vec2 VecObstacle_stretch(vec2 p, BasisFlow const& b);

    // Uses Newton iterations to inverse the bilinear coefficients for stretched coordinates,
    // going from stretched world space (p) to unstretched UV space. Based on
    // http://stackoverflow.com/questions/808441/inverse-bilinear-interpolation .
    vec2 QuadCoord(vec2 p, BasisFlow const& b);

    // Computes the Jacobian of the deformation form unstretched UV space to stretched world space.
    mat2 QuadCoordInvDeriv(vec2 uv, BasisFlow const& b);

    // Stores all particles in an acceleration grid for easy retrieval
    void SetParticlesInAccelGrid();

    // Advects all particles
    void ComputeParticleAdvection();

    // Adds new particles
    void SeedParticles();

    // Projects all buoyancy forces from particles onto the basis flows
    void AddParticleForcesToBasisFlows();

    // Computes blinear weights for basis advection. See Equation 20.
    // newCenter: position where basis bi is being moved
    // bi: basis being moved
    // bj: one basis on a corner of the cell where newCenter lands, which will receive part of
    // bi's contribution
    void ComputeNewCenterProportions(vec2& newCenter, BasisFlow& bi, BasisFlow& bj, vec2& interBasisDist);
    
    // Compute all basis flows advection
    void ComputeBasisAdvection();

    // Computed wavenumber from the basis flow's frequency. See Section 5.3 .
    inline float WavenumberBasis(BasisFlow& b)
    {
        return powf(2.f, 0.5f*(b.freqLvl.x + b.freqLvl.y));
    }

public:
    // GLFW callbacks
    static void CallbackWindowClosed(GLFWwindow* pGlfwWindow);
    static void CallbackKey(GLFWwindow* pGlfwWindow, int key, int scancode, int action, int mods);
    static void CallbackMouseScroll(GLFWwindow* pGlfwWindow, double xOffset, double yOffset);
    static void DebugCallback(GLenum source, GLenum /*type*/, GLuint /*id*/, GLenum severity, std::string message);

public:

    GLFWwindow * _glfwWindow;

    // Relative neighboring basis frequecies used during energy transfers
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
    
    // stores the B^T.B for basis flows for all of their neighbors, to more easily compute the energy transfer
    // of Equation 24.
    std::unique_ptr<DataBuffer1D<std::vector<CoeffBBDecompressedIntersectionInfo>*>>
        _intersectingBasesIdsDeformation[_nbExplicitTransferFreqs];

    struct ExplicitTransferCoeffs {
        float coeffs[_nbExplicitTransferFreqs];
    };
    // total contribution of beighboring basis, see denominator of Equation 24
    std::vector<ExplicitTransferCoeffs> _coeffBBExplicitTransferSum_abs;

    // list of all obstacles, inclusing simulation domain walls
    std::vector<Obstacle*> _obstacles;

    // list of all orthogonal groups of basis flows. Used when inverting the B^T.B matrix
    // with the multicolor scheme, see Section 5.1 .
    std::vector<std::vector<unsigned int>> _orthogonalBasisGroupIds;

    // coefficient dictionaries
    MapTypeBB _coeffsBB;
    MapTypeT _coeffsT;

    // Instead of evaluating the coefficient dictionaries, we store, for each basis, the list of it's
    // neighboring coefficients and corresponding interaction coefficients. See last paragraph of
    // section 5.4 . The decompressed BB coefficients do not include the basis itself in its list of
    // neihbors, since we do not need it when inverting the B^T.B matrix. The B^T.B coefficient of a 
    // basis with itself is also already stored in the basis's normSquared field. For the
    // decopressed T coefficient, the basis itself is inclused in the list, since we need to account
    // for self advection when computing basis advection.
    std::vector<std::vector<CoeffBBDecompressedIntersectionInfo>> _coeffsBBDecompressedIntersections;
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


