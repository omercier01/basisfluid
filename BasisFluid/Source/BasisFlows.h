#ifndef BASISFLOWS_H
#define BASISFLOWS_H

#include "Application.h"


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
#if SIM3D
    float back;
    float front;
#endif

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

bool intersectionInteriorEmpty(BasisSupport &sup1, BasisSupport &sup2);


#if SIM3D
enum class AXIS { X, Y, Z };
#endif


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
    float coeffLostInTransport;
    unsigned int orthoGroup;
    glm::vec2 stretchedCornerLB; // left-bottom
    glm::vec2 stretchedCornerLT; // left-top
    glm::vec2 stretchedCornerRB; // right-bottom
    glm::vec2 stretchedCornerRT; // right-top
//    glm::VECD stretchedCenter; // new center of stretched basis
    float normSquared; // matBBCoeff(this, this)
    bool stretched; // true if the basis is close to an obstacle and we must compute its stretched vector field when using it.

    BasisFlow(glm::ivec2 freq, glm::vec2 center, unsigned int orthogonalityGroup = -1) {
        coeff = 0;
        newCoeff = 0;
        coeffBoundary = 0;
        coeffLostInTransport = 0;
        orthoGroup = orthogonalityGroup;
        bitFlags = 0;
        this->freqLvl = freq;
        this->center = center;
        this->normSquared = 0; //1; // computed during basis setup.
        this->stretched = false;
    }

    BasisFlow() {
        coeff = 0;
        newCoeff = 0;
        coeffBoundary = 0;
        coeffLostInTransport = 0;
        orthoGroup = -1;
        bitFlags = 0;
        this->freqLvl = glm::vec2(0);
        this->center = glm::vec2(0);
#if SIM3D
        this->axis = AXIS(-1);
#endif
        this->normSquared = 0; //1;
    }

    BasisSupport getSupport() const;
    glm::vec2 supportHalfSize() const;
    glm::vec2 normalizedPositionInSupport(glm::vec2 p);
    bool pointIsInSupport(glm::vec2 p);
    bool emptyIntersectionWithBasis(BasisFlow& b) const;
};

glm::dvec2 flowBasisHat(glm::dvec2 p, int log2Aniso);
glm::dmat2 flowBasisHatGrad(glm::dvec2 p, int log2Aniso);


// TEST: to see the basis domains more clearly
glm::vec2 squareFlow(glm::vec2 p, BasisFlow b);


void InverseBBMatrix(
    DataBuffer1D<scalar_inversion_storage>* vecX,
    DataBuffer1D<scalar_inversion_storage>* vecB,
    float tol, unsigned int basisBitMask, unsigned int minFreq);

float IntegrateBasisGrid(BasisFlow& b, VectorField2D* velField);
float IntegrateBasisBasis(BasisFlow b1, BasisFlow b2);
glm::dvec2 TranslatedBasisEvalPrecise(const glm::dvec2 p, const glm::ivec2 freqLvl, const glm::dvec2 center);
glm::mat2 TranslatedBasisGradEval(const glm::dvec2 p, const glm::ivec2 freqLvl, const glm::vec2 center);
glm::dmat2 TranslatedBasisGradEvalPrecise(const glm::dvec2 p, const glm::ivec2 freqLvl, const glm::dvec2 center);

#endif // BASISFLOWS_H
