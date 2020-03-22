#ifndef BASISFLOWS_H
#define BASISFLOWS_H

#include <glm\glm.hpp>

#include <tuple>
#include <unordered_map>
#include <functional>


// Dictionary structure for BB coefficients
// keys are freqLvlX1, freqLvlY1, freqLvlX2, freqLvlY2, centerDiffX, centerDiffY
typedef std::tuple<int, int, int, int, float, float> KeyTypeBB;
struct KeyTypeBB_hash
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


// Dictionary structure for T coefficients
// keys are freqLvlX1, freqLvlY1, freqLvlX2, freqLvlY2, centerDiffX, centerDiffY
typedef std::tuple<int, int, int, int, float, float> KeyTypeT;
struct KeyTypeT_hash
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


// Struture to store neighboring T coefficients
struct CoeffTDecompressedIntersectionInfo {
    int j;
    glm::vec2 coeff;

    CoeffTDecompressedIntersectionInfo(int j, glm::vec2 coeff) {
        this->j = j;
        this->coeff = coeff;
    }
};

// Struture to store neighboring BB coefficients
struct CoeffBBDecompressedIntersectionInfo {
    int j;
    float coeff;

    CoeffBBDecompressedIntersectionInfo(int j, float coeff) {
        this->j = j;
        this->coeff = coeff;
    }
};


// Rectangular basis flow support (outside of which the basis is zero)
struct BasisSupport {
    float left;
    float right;
    float bottom;
    float top;

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


// True is supports sup1 and sup2 do not intersect
bool IntersectionInteriorEmpty(const BasisSupport &sup1, const BasisSupport &sup2);


enum BASIS_FLAGS {
    // Basis flow is in the interior od the domain. Could not be the case if a dunamic obstacle covers
    // the basis flow during the simualtion.
    INTERIOR = 1 << 0, 

    // Basis is used for force projection
    FORCE_PROJECTION = 1 << 1,

    // Basis flow is used for dynamic obstacle projection
    DYNAMIC_BOUNDARY_PROJECTION = 1 << 2
};


struct BasisFlow {

    glm::ivec2 freqLvl; // Log2 of basis frequency. Corresponds to \alpha and \beta in Section 4.1 .
    glm::vec2 center; // basis flow center. Corresponds to c_x and c_y in Section 4.1 .
    unsigned int bitFlags; // bitfield of BASIS_FLAGS
    unsigned int stretchBitFlags; // bitfield of BASIS_FLAGS for the stretched basis
    float coeff; // Basis flow coefficient, corresponds to \tilde{u} in Section 3
    float newCoeff; //  used when computing new coefficient
    float coeffBoundary; // basis coefficient or projected dynamic obstacle motion
    glm::vec2 stretchedCornerLB; // left-bottom stretched support corner
    glm::vec2 stretchedCornerLT; // left-top stretched support corner
    glm::vec2 stretchedCornerRB; // right-bottom stretched support corner
    glm::vec2 stretchedCornerRT; // right-top stretched support corner
    float normSquared; // diagonal term in B^T.B matrix
    bool stretched; // true if the basis is near an obstacle and must be evaluated with a stretch

    BasisFlow(glm::ivec2 freq, glm::vec2 center) {
        coeff = 0;
        newCoeff = 0;
        coeffBoundary = 0;
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
        bitFlags = 0;
        this->freqLvl = glm::vec2(0);
        this->center = glm::vec2(0);
        this->normSquared = 0;
    }

    // Returns  basis flow's support
    BasisSupport getSupport() const;

    // Distance from basis flow center to edge of support, in each direction
    glm::vec2 supportHalfSize() const;
};


// Basis template of Equation 14.
glm::dvec2 flowBasisHat(glm::dvec2 p, int log2Aniso);


#endif // BASISFLOWS_H
