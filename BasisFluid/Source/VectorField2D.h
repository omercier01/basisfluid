// vector field defined by vectors on a regular grid.

#ifndef VECTORFIELD2D_H
#define VECTORFIELD2D_H

#include "DataBuffer2D.h"
#include "DataBuffer1D.h"

#include "glm/glm.hpp"

#include <functional>


class VectorField2D
{
public:
    enum class BoundaryCondition{ LINEAR, FLAT, ZERO, MIRROR, PERIODIC };
    enum class GridNodeLocation{ CENTER, CORNER };
    enum class InterpolationMethod{ LINEAR };
public:
    DataBuffer2D<glm::vec2> _vectors;
    float _boundXMin, _boundYMin, _boundXMax, _boundYMax;
    unsigned int _nbCellsX, _nbCellsY;
    BoundaryCondition _boundaryCondition;
    GridNodeLocation _gridNodeLocation;

public:
    VectorField2D(float boundXMin, float boundXMax, float boundYMin,
                  float boundYMax, unsigned int nbCellsX, unsigned int nbCellsY,
                  GridNodeLocation gridNodeLocation = GridNodeLocation::CENTER,
                  BoundaryCondition boundaryCondition = BoundaryCondition::LINEAR);

    void populateWithFunction(std::function<glm::vec2 (float x, float y)> function);
    void addFunction(std::function<glm::vec2 (float x, float y)> function);
    vec2 interp(glm::vec2 pos);

    void addVectorCpuData(unsigned int i, unsigned int j, glm::vec2 data);
    void setVectorCpuData(unsigned int i, unsigned int j, glm::vec2 data);
    glm::vec2 getVectorCpuData(unsigned int i, unsigned int j);

    void createVectorCpuStorage();
    void createVectorTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                      GLenum externalFormat, GLenum sizedExternalFormat,
                                      unsigned int nbMipmapLevels);
    void createVectorTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                      GLenum externalFormat, GLenum sizedExternalFormat);
    void vectorTexture2DStorageAdjustBoundaryCondition();

    void setVectorsFromImage(Metadata2DImage2D metadataImage2D, unsigned int level);

    unsigned int nbElementsX();
    unsigned int nbElementsY();
    
    void setBounds(float inBoundXMin, float inBoundXMax,
                   float inBoundYMin, float inBoundYMax);

    glm::uvec2 pointToClosestIndex(vec2 point);
    glm::uvec2 pointToFlooredIndex(vec2 point);
    vec2 indexToPosition(uvec2 index);
};

#endif // VECTORFIELD2D_H
