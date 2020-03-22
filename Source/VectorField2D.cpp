#include "VectorField2D.h"


using namespace glm;


VectorField2D::VectorField2D(
    float boundXMin, float boundXMax, float boundYMin, float boundYMax,
    unsigned int nbCellsX, unsigned int nbCellsY) :
    _nbCellsX(nbCellsX),
    _nbCellsY(nbCellsY),
    _vectors(1, 1) {

    _boundXMin = boundXMin;
    _boundXMax = boundXMax;
    _boundYMin = boundYMin;
    _boundYMax = boundYMax;

    // resizes, but does not necessarily allocate the data.
    _vectors.resize(nbCellsX + 1, nbCellsY + 1);
}


void VectorField2D::populateWithFunction(std::function<vec2(float x, float y)> function)
{
    vec2* _vectorsPointer = _vectors.getCpuDataPointer();
    unsigned int nxVec = _vectors._nbElementsX;

    for (unsigned int i = 0; i < _nbCellsX + 1; i++)
        for (unsigned int j = 0; j < _nbCellsY + 1; j++) {
            float x = _boundXMin + float(i) / _nbCellsX * (_boundXMax - _boundXMin);
            float y = _boundYMin + float(j) / _nbCellsY * (_boundYMax - _boundYMin);
            _vectorsPointer[nxVec*j + i] = function(x, y);
        }
}


void VectorField2D::addFunction(std::function<vec2(float, float)> function)
{
    vec2* _vectorsPointer = _vectors.getCpuDataPointer();
    unsigned int nxVec = _vectors._nbElementsX;

    for (unsigned int i = 0; i < _nbCellsX + 1; i++)
        for (unsigned int j = 0; j < _nbCellsY + 1; j++) {
            float x = _boundXMin + float(i) / _nbCellsX * (_boundXMax - _boundXMin);
            float y = _boundYMin + float(j) / _nbCellsY * (_boundYMax - _boundYMin);
            _vectorsPointer[nxVec*j + i] += function(x, y);
        }
}


vec2 VectorField2D::interp(vec2 pos)
{
    float x = pos.x;
    float y = pos.y;

    int indexXLeft, indexXRight;
    int indexYLeft, indexYRight;
    float weightX;
    float weightY;

    float normalizedX = (x - _boundXMin) / (_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (y - _boundYMin) / (_boundYMax - _boundYMin)*_nbCellsY;

    if (normalizedX < 0) {
        return vec2(0, 0);
    }
    else if (normalizedX >= _nbCellsX) {
        return vec2(0, 0);
    }
    else {
        indexXLeft = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
        indexXRight = indexXLeft + 1;
        weightX = normalizedX - indexXLeft;
    }
    if (normalizedY < 0) {
        return vec2(0, 0);
    }
    else if (normalizedY >= _nbCellsY) {
        return vec2(0, 0);
    }
    else {
        indexYLeft = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
        indexYRight = indexYLeft + 1;
        weightY = normalizedY - indexYLeft;
    }

    vec2* _vectorsPointer = _vectors.getCpuDataPointer();

    const unsigned int nx = _vectors._nbElementsX;

    return (1 - weightX)*(
        (1 - weightY)*_vectorsPointer[nx*indexYLeft + indexXLeft] +
        (weightY)*_vectorsPointer[nx*indexYRight + indexXLeft]
        )
        +
        (weightX)*(
        (1 - weightY)*_vectorsPointer[nx*indexYLeft + indexXRight] +
            (weightY)*_vectorsPointer[nx*indexYRight + indexXRight]
            );
}


void VectorField2D::addVectorCpuData(unsigned int i, unsigned int j, vec2 data)
{
    _vectors.addCpuData(i, j, data);
}


void VectorField2D::setVectorCpuData(unsigned int i, unsigned int j, vec2 data)
{
    _vectors.setCpuData(i, j, data);
}


glm::vec2 VectorField2D::getVectorCpuData(unsigned int i, unsigned int j)
{
    return _vectors.getCpuData(i, j);
}


void VectorField2D::createVectorCpuStorage()
{
    _vectors.createCpuStorage();
}


void VectorField2D::createVectorTexture2DStorage(
    GLenum internalFormat,
    GLenum sizedInternalFormat,
    GLenum externalFormat,
    GLenum sizedExternalFormat,
    unsigned int nbMipmapLevels)
{
    _vectors.createTexture2DStorage(internalFormat, sizedInternalFormat,
        externalFormat, sizedExternalFormat,
        nbMipmapLevels);
    vectorTexture2DStorageAdjustBoundaryCondition();
}


void VectorField2D::createVectorTexture2DStorage(
    GLenum internalFormat,
    GLenum sizedInternalFormat,
    GLenum externalFormat,
    GLenum sizedExternalFormat)
{
    _vectors.createTexture2DStorage(internalFormat, sizedInternalFormat,
        externalFormat, sizedExternalFormat);
    vectorTexture2DStorageAdjustBoundaryCondition();
}


void VectorField2D::vectorTexture2DStorageAdjustBoundaryCondition()
{
    float borderColor[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTextureParameterfv(_vectors._glidTexture2D, GL_TEXTURE_BORDER_COLOR, borderColor);
}


void VectorField2D::setVectorsFromImage(Metadata2DImage2D metadataImage2D, unsigned int level)
{
    _vectors.setFromImage(metadataImage2D, level);
}


unsigned int VectorField2D::nbElementsX()
{
    return _nbCellsX + 1;
}


unsigned int VectorField2D::nbElementsY()
{
    return _nbCellsY + 1;
}


void VectorField2D::setBounds(float in_boundXMin, float in_boundXMax, float in_boundYMin, float in_boundYMax)
{
    _boundXMin = in_boundXMin;
    _boundXMax = in_boundXMax;
    _boundYMin = in_boundYMin;
    _boundYMax = in_boundYMax;
}


uvec2 VectorField2D::pointToClosestNodeIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - _boundXMin) / (_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (point.y - _boundYMin) / (_boundYMax - _boundYMin)*_nbCellsY;

    if (normalizedX < 0) {
        index.x = 0;
    }
    else if (normalizedX >= _nbCellsX) {
        index.x = _nbCellsX;
    }
    else {
        index.x = clamp<int>(int(floor(normalizedX + 0.5)), 0, _nbCellsX);
    }
    if (normalizedY < 0) {
        index.y = 0;
    }
    else if (normalizedY >= _nbCellsY) {
        index.y = _nbCellsY;
    }
    else {
        index.y = clamp<int>(int(floor(normalizedY + 0.5)), 0, _nbCellsY);
    }

    return index;
}


uvec2 VectorField2D::pointToCellIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - _boundXMin) / (_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (point.y - _boundYMin) / (_boundYMax - _boundYMin)*_nbCellsY;

    index.x = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX);
    index.y = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY);

    return index;
}


vec2 VectorField2D::indexToPosition(uvec2 index)
{
    vec2 pos;

    pos.x = _boundXMin + (index.x)*(_boundXMax - _boundXMin) / _nbCellsX;
    pos.y = _boundYMin + (index.y)*(_boundYMax - _boundYMin) / _nbCellsY;

    return pos;
}


// copies and returns a linear buffer containing the grid positions.
Metadata1DCpu VectorField2D::GenerateGridNodeLocations() {

    Metadata1DCpu metadata;
    vec2* data = nullptr;
    float cellWidthX = (_boundXMax - _boundXMin) / _nbCellsX;
    float cellWidthY = (_boundYMax - _boundYMin) / _nbCellsY;

    data = new vec2[(_nbCellsX + 1)*(_nbCellsY + 1)];
    for (unsigned int i = 0; i <= _nbCellsX; i++) {
        for (unsigned int j = 0; j <= _nbCellsY; j++) {
            data[(_nbCellsY + 1)*j + i] = vec2(_boundXMin + i * cellWidthX,
                _boundYMin + j * cellWidthY);
        }
    }

    metadata.dataPointer = data;
    return metadata;
}

