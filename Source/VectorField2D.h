// vector field defined by vectors on a regular grid. Data is stored on cell corners.

#ifndef VECTORFIELD2D_H
#define VECTORFIELD2D_H

#include "DataBuffer2D.h"
#include "DataBuffer1D.h"

#include "glm/glm.hpp"

#include <functional>

class VectorField2D
{
public:
    DataBuffer2D<glm::vec2> _vectors;
    float _boundXMin, _boundYMin, _boundXMax, _boundYMax;
    unsigned int _nbCellsX, _nbCellsY;

public:
    VectorField2D(float boundXMin, float boundXMax, float boundYMin,
        float boundYMax, unsigned int nbCellsX, unsigned int nbCellsY);

    void populateWithFunction(std::function<glm::vec2(float x, float y)> function);
    void addFunction(std::function<glm::vec2(float x, float y)> function);
    glm::vec2 interp(glm::vec2 pos);

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

    glm::uvec2 pointToCellIndex(glm::vec2 point);
    glm::uvec2 pointToClosestNodeIndex(glm::vec2 point);
    glm::vec2 indexToPosition(glm::uvec2 index);

    Metadata1DCpu GenerateGridNodeLocations();
};

#endif // VECTORFIELD2D_H
