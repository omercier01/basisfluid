#ifndef BASISFLOWS_H
#define BASISFLOWS_H

#include <glm\glm.hpp>

#include <tuple>
#include <unordered_map>
#include <functional>


// coefficients for BB matrix
// keys are freqLvlX1, freqLvlY1, freqLvlX2, freqLvlY2, centerDiffX, centerDiffY
typedef std::tuple<int, int, int, int, float, float> KeyTypeBB;
struct KeyTypeBB_hash// : public std::unary_function<KeyTypeBB, std::size_t>
{
    std::size_t operator()(const KeyTypeBB& k) const
    {
        std::hash<float> hasher_float;
        std::hash<int> hasher_uint;
        return hasher_uint (std::get<0>(k)) ^ hasher_uint (std::get<1>(k)) ^
            hasher_uint (std::get<2>(k)) ^ hasher_uint (std::get<3>(k)) ^
            hasher_float(std::get<4>(k)) ^ hasher_float(std::get<5>(k));
    }
};
typedef std::unordered_map<KeyTypeBB,float,KeyTypeBB_hash> MapTypeBB;

// coefficients for T matrix
// keys are freqLvlX1, freqLvlY1, freqLvlX2, freqLvlY2, centerDiffX, centerDiffY
typedef std::tuple<int, int, int, int, float, float> KeyTypeT;
struct KeyTypeT_hash// : public std::unary_function<KeyTypeT, std::size_t>
{
    std::size_t operator()(const KeyTypeT& k) const
    {
        std::hash<float> hasher_float;
        std::hash<int> hasher_uint;
    return hasher_uint (std::get<0>(k)) ^ hasher_uint (std::get<1>(k)) ^
            hasher_uint (std::get<2>(k)) ^ hasher_uint (std::get<3>(k)) ^
            hasher_float(std::get<4>(k)) ^ hasher_float(std::get<5>(k));
    }
};
typedef std::unordered_map<KeyTypeT,glm::vec2,KeyTypeT_hash> MapTypeT;



// stored in a vector indexed by _k_, so _k_ does not need to be stored here.
struct CoeffADecompressedIntersectionInfo {
    int i;
    int j;
    float coeff;

    CoeffADecompressedIntersectionInfo(int i, int j, float coeff) {
        this->i = i;
        this->j = j;
        this->coeff = coeff;
    }
};


// stored in a vector indexed by _i_, so _i_ does not need to be stored here.
struct CoeffTDecompressedIntersectionInfo {
    int j;
    glm::vec2 coeff;

    CoeffTDecompressedIntersectionInfo(int j, glm::vec2 coeff) {
        this->j = j;
        this->coeff = coeff;
    }
};

// stored in a vector indexed by _i_, so _i_ does not need to be stored here.
struct CoeffBBDecompressedIntersectionInfo {
    int j;
    float coeff;

    CoeffBBDecompressedIntersectionInfo(int j, float coeff) {
        this->j = j;
        this->coeff = coeff;
    }
};


struct BasisSupport {
    float left;
    float right;
    float bottom;
    float top;

    // supports are rectangular. They are used for computations, while the stretched corners are used to compute the effective basis velocity field (but if the basis is not stretched, then the basis support can be used too).
    BasisSupport(
        float left,
        float right,
        float bottom,
        float top
    ) {
        this->left = left;
        this->right = right;
        this->bottom = bottom;
        this->top = top;
    }
};

bool IntersectionInteriorEmpty(const BasisSupport &sup1, const BasisSupport &sup2);



enum BASIS_FLAGS {
    INTERIOR = 1 << 0,
    FORCE_PROJECTION = 1 << 1,  // used to restrict force projection to coarser frequencies
    DYNAMIC_BOUNDARY_PROJECTION = 1 << 2
};


struct BasisFlow {

    glm::ivec2 freqLvl;       // leave at this position (for pipelineBasisContour)
    glm::vec2 center;         // leave at this position
    //uint valid;     // replaced with bitflags. // LEAVE HERE. false if the basis is too compressed or is inside the obstacle.
    unsigned int bitFlags; // LEAVE HERE for pipelines.
    unsigned int prevBitFlags;
    unsigned int tempBitFlags;
    unsigned int stretchBitFlags;
    float coeff;
    float newCoeff;
    float coeffBoundary; // stores the coeficient not parallel to the boundary, necessary to represent object movement flow.
    unsigned int orthoGroup;
    glm::vec2 stretchedCornerLB; // left-bottom
    glm::vec2 stretchedCornerLT; // left-top
    glm::vec2 stretchedCornerRB; // right-bottom
    glm::vec2 stretchedCornerRT; // right-top
    float normSquared; // matBBCoeff(this, this)
    bool stretched; // true if the basis is close to an obstacle and we must compute its stretched vector field when using it.

    BasisFlow(glm::ivec2 freq, glm::vec2 center, unsigned int orthogonalityGroup = -1) {
        coeff = 0;
        newCoeff = 0;
        coeffBoundary = 0;
        orthoGroup = orthogonalityGroup;
        bitFlags = 0;
        this->freqLvl = freq;
        this->center = center;
        this->normSquared = 0; // computed during basis setup.
        this->stretched = false;
    }

    BasisFlow() {
        coeff = 0;
        newCoeff = 0;
        coeffBoundary = 0;
        orthoGroup = -1;
        bitFlags = 0;
        this->freqLvl = glm::vec2(0);
        this->center = glm::vec2(0);
        this->normSquared = 0;
    }

    BasisSupport getSupport() const;
    glm::vec2 supportHalfSize() const;
    glm::vec2 normalizedPositionInSupport(glm::vec2 p);
    bool pointIsInSupport(glm::vec2 p);
    bool EmptyIntersectionWithBasis(const BasisFlow& b) const;
};

glm::dvec2 flowBasisHat(glm::dvec2 p, int log2Aniso);
glm::dmat2 flowBasisHatGrad(glm::dvec2 p, int log2Aniso);

#endif // BASISFLOWS_H
