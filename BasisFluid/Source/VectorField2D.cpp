#include "VectorField2D.h"

using namespace glm;

VectorField2D::VectorField2D(
    float boundXMin, float boundXMax, float boundYMin, float boundYMax,
    unsigned int nbCellsX, unsigned int nbCellsY,
    VectorField2D::GridNodeLocation gridNodeLocation,
    VectorField2D::BoundaryCondition boundaryCondition) :
    _nbCellsX(nbCellsX),
    _nbCellsY(nbCellsY),
    _vectors(1, 1),
    _gridNodeLocation(gridNodeLocation),
    _boundaryCondition(boundaryCondition) {

    _boundXMin = boundXMin;
    _boundXMax = boundXMax;
    _boundYMin = boundYMin;
    _boundYMax = boundYMax;

    // resizes, but does not necessarily allocate the data.
    switch (gridNodeLocation) {
    case GridNodeLocation::CENTER:
        _vectors.resize(nbCellsX, nbCellsY);
        break;
    case GridNodeLocation::CORNER:
        _vectors.resize(nbCellsX + 1, nbCellsY + 1);
        break;
    }
}

void VectorField2D::populateWithFunction(std::function<vec2(float x, float y)> function)
{
    vec2* _vectorsPointer = _vectors.getCpuDataPointer();
    unsigned int nxVec = _vectors._nbElementsX;

    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        for (unsigned int i = 0; i < _nbCellsX; i++)
            for (unsigned int j = 0; j < _nbCellsY; j++) {
                float x = _boundXMin + (i + 0.5f) / _nbCellsX * (_boundXMax - _boundXMin);
                float y = _boundYMin + (j + 0.5f) / _nbCellsY * (_boundYMax - _boundYMin);
                //            _vectors(i,j) = function(x, y);
                            // TODO: setting dtaa one by one like this is probably sloooow!
                            //_vectors.setCpuData(i, j, function(x, y));
                _vectorsPointer[nxVec*j + i] = function(x, y);
            }
        break;
    case GridNodeLocation::CORNER:
        for (unsigned int i = 0; i < _nbCellsX + 1; i++)
            for (unsigned int j = 0; j < _nbCellsY + 1; j++) {
                float x = _boundXMin + float(i) / _nbCellsX * (_boundXMax - _boundXMin);
                float y = _boundYMin + float(j) / _nbCellsY * (_boundYMax - _boundYMin);
                //            _vectors(i,j) = function(x, y);
                            //_vectors.setCpuData(i, j, function(x, y));
                _vectorsPointer[nxVec*j + i] = function(x, y);
            }
        break;
    }

    _vectors._sourceStorageType = DataBuffer2D<vec2>::StorageType::CPU;
    //_vectors.dirtyData();

}

void VectorField2D::addFunction(std::function<vec2(float, float)> function)
{

    vec2* _vectorsPointer = _vectors.getCpuDataPointer();
    unsigned int nxVec = _vectors._nbElementsX;

    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        for (unsigned int i = 0; i < _nbCellsX; i++)
            for (unsigned int j = 0; j < _nbCellsY; j++) {
                float x = _boundXMin + (i + 0.5f) / _nbCellsX * (_boundXMax - _boundXMin);
                float y = _boundYMin + (j + 0.5f) / _nbCellsY * (_boundYMax - _boundYMin);
                // TODO: setting data one by one like this is probably sloooow!
                //_vectors.addCpuData(i, j, function(x, y));
                _vectorsPointer[nxVec*j + i] += function(x, y);
            }
        break;
    case GridNodeLocation::CORNER:
        for (unsigned int i = 0; i < _nbCellsX + 1; i++)
            for (unsigned int j = 0; j < _nbCellsY + 1; j++) {
                float x = _boundXMin + float(i) / _nbCellsX * (_boundXMax - _boundXMin);
                float y = _boundYMin + float(j) / _nbCellsY * (_boundYMax - _boundYMin);
                //_vectors.addCpuData(i, j, function(x, y));
                _vectorsPointer[nxVec*j + i] += function(x, y);
            }
        break;
    }

    _vectors._sourceStorageType = DataBuffer2D<vec2>::StorageType::CPU;
    //_vectors.dirtyData();

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


    switch (_boundaryCondition) {
    case BoundaryCondition::LINEAR:
        switch (_gridNodeLocation) {
        case GridNodeLocation::CENTER:
            if (normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = normalizedX - 0.5f;
            }
            else if (normalizedX >= _nbCellsX - 0.5f) {
                indexXLeft = _nbCellsX - 2;
                indexXRight = _nbCellsX - 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, _nbCellsX - 2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if (normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = normalizedY - 0.5f;
            }
            else if (normalizedY >= _nbCellsY - 0.5f) {
                indexYLeft = _nbCellsY - 2;
                indexYRight = _nbCellsY - 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, _nbCellsY - 2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER:
            if (normalizedX < 0) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = normalizedX;
            }
            else if (normalizedX >= _nbCellsX) {
                indexXLeft = _nbCellsX - 1;
                indexXRight = _nbCellsX;
                weightX = normalizedX - indexXLeft;
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - indexXLeft;
            }
            if (normalizedY < 0) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = normalizedY;
            }
            else if (normalizedY >= _nbCellsY) {
                indexYLeft = _nbCellsY - 1;
                indexYRight = _nbCellsY;
                weightY = normalizedY - indexYLeft;
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - indexYLeft;
            }
            break;
        }
        break;
    case BoundaryCondition::FLAT:
        switch (_gridNodeLocation) {
        case GridNodeLocation::CENTER:
            if (normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = 0;
            }
            else if (normalizedX >= _nbCellsX - 0.5f) {
                indexXLeft = _nbCellsX - 2;
                indexXRight = _nbCellsX - 1;
                weightX = 1;
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, _nbCellsX - 2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if (normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = 0;
            }
            else if (normalizedY >= _nbCellsY - 0.5f) {
                indexYLeft = _nbCellsY - 2;
                indexYRight = _nbCellsY - 1;
                weightY = 1;
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, _nbCellsY - 2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER:
            if (normalizedX < 0) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = 0;
            }
            else if (normalizedX >= _nbCellsX) {
                indexXLeft = _nbCellsX - 1;
                indexXRight = _nbCellsX;
                weightX = 1;
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - indexXLeft;
            }
            if (normalizedY < 0) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = 0;
            }
            else if (normalizedY >= _nbCellsY) {
                indexYLeft = _nbCellsY - 1;
                indexYRight = _nbCellsY;
                weightY = 1;
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - indexYLeft;
            }
            break;
        }
        break;
    case BoundaryCondition::ZERO:
        switch (_gridNodeLocation) {
        case GridNodeLocation::CENTER:
            if (normalizedX < 0.5) {
                return vec2(0, 0);
            }
            else if (normalizedX >= _nbCellsX - 0.5f) {
                return vec2(0, 0);
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, _nbCellsX - 2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if (normalizedY < 0.5) {
                return vec2(0, 0);
            }
            else if (normalizedY >= _nbCellsY - 0.5f) {
                return vec2(0, 0);
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, _nbCellsY - 2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER:
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
            break;
        }
        break;
    case BoundaryCondition::MIRROR:
        switch (_gridNodeLocation) {
        case GridNodeLocation::CENTER:
            if (normalizedX < 0) {
                normalizedX = -normalizedX;
            }
            if (normalizedX >= _nbCellsX - 0.5f) {
                float xInDomain01 = normalizedX / _nbCellsX;
                unsigned int temp = int(floor(xInDomain01)) % 2;
                if (temp == 0) {
                    normalizedX = _nbCellsX * fract(xInDomain01);
                }
                else {
                    normalizedX = _nbCellsX * (1 - fract(xInDomain01));
                }
            }
            if (normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = indexXLeft;
                weightX = 0;
            }
            else if (normalizedX >= _nbCellsX - 0.5f) {
                indexXLeft = _nbCellsX - 1;
                indexXRight = indexXLeft;
                weightX = 0;
            }
            else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, _nbCellsX - 2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }

            if (normalizedY < 0) {
                normalizedY = -normalizedY;
            }
            if (normalizedY >= _nbCellsY - 0.5f) {
                float yInDomain01 = normalizedY / _nbCellsY;
                unsigned int temp = int(floor(yInDomain01)) % 2;
                if (temp == 0) {
                    normalizedY = _nbCellsY * fract(yInDomain01);
                }
                else {
                    normalizedY = _nbCellsY * (1 - fract(yInDomain01));
                }
            }
            if (normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = indexYLeft;
                weightY = normalizedY - indexYLeft;
            }
            else if (normalizedY >= _nbCellsY - 0.5f) {
                indexYLeft = _nbCellsY - 1;
                indexYRight = indexYLeft;
                weightY = normalizedY - indexYLeft;
            }
            else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, _nbCellsY - 2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER:
            if (normalizedX < 0) {
                normalizedX = -normalizedX;
            }
            if (normalizedX >= _nbCellsX) {
                float xInDomain01 = normalizedX / _nbCellsX;
                unsigned int temp = int(floor(xInDomain01)) % 2;
                if (temp == 0) {
                    normalizedX = _nbCellsX * fract(xInDomain01);
                }
                else {
                    normalizedX = _nbCellsX * (1 - fract(xInDomain01));
                }
            }
            indexXLeft = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - indexXLeft;

            if (normalizedY < 0) {
                normalizedY = -normalizedY;
            }
            if (normalizedY >= _nbCellsY) {
                float yInDomain01 = normalizedY / _nbCellsY;
                unsigned int temp = int(floor(yInDomain01)) % 2;
                if (temp == 0) {
                    normalizedY = _nbCellsY * fract(yInDomain01);
                }
                else {
                    normalizedY = _nbCellsY * (1 - fract(yInDomain01));
                }
            }
            indexYLeft = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - indexYLeft;
            break;
        }
        break;
    case BoundaryCondition::PERIODIC:
        switch (_gridNodeLocation) {
        case GridNodeLocation::CENTER:
            normalizedX -= _nbCellsX * floor(normalizedX / _nbCellsX);
            indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, _nbCellsX - 2);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - (indexXLeft + 0.5f);
            normalizedY -= _nbCellsY * floor(normalizedY / _nbCellsY);
            indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, _nbCellsY - 2);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - (indexYLeft + 0.5f);
            break;
        case GridNodeLocation::CORNER:
            normalizedX -= _nbCellsX * floor(normalizedX / _nbCellsX);
            indexXLeft = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - indexXLeft;
            normalizedY -= _nbCellsY * floor(normalizedY / _nbCellsY);
            indexYLeft = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - indexYLeft;
            break;
        }
        break;
    }

    //_vectors.refreshCpuData();
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
    //LINEAR, FLAT, ZERO, MIRROR, PERIODIC
    float borderColor[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    switch (_boundaryCondition) {
    case VectorField2D::BoundaryCondition::ZERO:
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glTextureParameterfv(_vectors._glidTexture2D, GL_TEXTURE_BORDER_COLOR, borderColor);
        break;
    case VectorField2D::BoundaryCondition::LINEAR:
        // not supported by OpenGL, substiture with FLAT boundary condition
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        break;
    case VectorField2D::BoundaryCondition::FLAT:
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        break;
    case VectorField2D::BoundaryCondition::MIRROR:
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
        break;
    case VectorField2D::BoundaryCondition::PERIODIC:
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTextureParameteri(_vectors._glidTexture2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        break;
    default:
        // TODO: should not happen
        break;
    }
}


void VectorField2D::setVectorsFromImage(Metadata2DImage2D metadataImage2D, unsigned int level)
{
    _vectors.setFromImage(metadataImage2D, level);
}

unsigned int VectorField2D::nbElementsX()
{
    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        return _nbCellsX;
        break;
    case GridNodeLocation::CORNER:
        return _nbCellsX + 1;
        break;
    default:
        // TODO: error, unknows grid location type.
        return 0;
    }
}

unsigned int VectorField2D::nbElementsY()
{
    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        return _nbCellsY;
        break;
    case GridNodeLocation::CORNER:
        return _nbCellsY + 1;
        break;
    default:
        // TODO: error, unknows grid location type.
        return 0;
    }
}


void VectorField2D::setBounds(float in_boundXMin, float in_boundXMax, float in_boundYMin, float in_boundYMax)
{
    _boundXMin = in_boundXMin;
    _boundXMax = in_boundXMax;
    _boundYMin = in_boundYMin;
    _boundYMax = in_boundYMax;
}


uvec2 VectorField2D::pointToClosestIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - _boundXMin) / (_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (point.y - _boundYMin) / (_boundYMax - _boundYMin)*_nbCellsY;


    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        if (normalizedX < 0) {
            index.x = 0;
        }
        else if (normalizedX >= _nbCellsX) {
            index.x = _nbCellsX - 1;
        }
        else {
            index.x = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX - 1);
        }
        if (normalizedY < 0) {
            index.y = 0;
        }
        else if (normalizedY >= _nbCellsY) {
            index.y = _nbCellsY - 1;
        }
        else {
            index.y = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY - 1);
        }
        break;
    case GridNodeLocation::CORNER:
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
        break;
    }

    return index;
}

uvec2 VectorField2D::pointToFlooredIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - _boundXMin) / (_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (point.y - _boundYMin) / (_boundYMax - _boundYMin)*_nbCellsY;


    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        if (normalizedX < 0) {
            index.x = 0;
        }
        else if (normalizedX >= _nbCellsX) {
            index.x = _nbCellsX - 1;
        }
        else {
            index.x = clamp<int>(int(floor(normalizedX - 0.5)), 0, _nbCellsX - 1);
        }
        if (normalizedY < 0) {
            index.y = 0;
        }
        else if (normalizedY >= _nbCellsY) {
            index.y = _nbCellsY - 1;
        }
        else {
            index.y = clamp<int>(int(floor(normalizedY - 0.5)), 0, _nbCellsY - 1);
        }
        break;
    case GridNodeLocation::CORNER:
        if (normalizedX < 0) {
            index.x = 0;
        }
        else if (normalizedX >= _nbCellsX) {
            index.x = _nbCellsX;
        }
        else {
            index.x = clamp<int>(int(floor(normalizedX)), 0, _nbCellsX);
        }
        if (normalizedY < 0) {
            index.y = 0;
        }
        else if (normalizedY >= _nbCellsY) {
            index.y = _nbCellsY;
        }
        else {
            index.y = clamp<int>(int(floor(normalizedY)), 0, _nbCellsY);
        }
        break;
    }

    return index;
}



vec2 VectorField2D::indexToPosition(uvec2 index)
{
    vec2 pos;

    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        pos.x = _boundXMin + (index.x + 0.5f)*(_boundXMax - _boundXMin) / _nbCellsX;
        pos.y = _boundYMin + (index.y + 0.5f)*(_boundYMax - _boundYMin) / _nbCellsY;
        break;
    case GridNodeLocation::CORNER:
        pos.x = _boundXMin + (index.x)*(_boundXMax - _boundXMin) / _nbCellsX;
        pos.y = _boundYMin + (index.y)*(_boundYMax - _boundYMin) / _nbCellsY;
        break;
    }

    return pos;
}




// copies and returns a linear buffer containing the grid positions.
Metadata1DCpu VectorField2D::GenerateGridNodeLocations() {

    //_vectors._dataCpu.refresh();

    Metadata1DCpu metadata;
    vec2* data = nullptr;
    float cellWidthX = (_boundXMax - _boundXMin) / _nbCellsX;
    float cellWidthY = (_boundYMax - _boundYMin) / _nbCellsY;

    switch (_gridNodeLocation) {
    case GridNodeLocation::CENTER:
        data = new vec2[_nbCellsX*_nbCellsY];
        for (unsigned int i = 0; i < _nbCellsX; i++) {
            for (unsigned int j = 0; j < _nbCellsY; j++) {
                data[_nbCellsY*j + i] = vec2(_boundXMin + (i + 0.5f)*cellWidthX,
                    _boundYMin + (j + 0.5f)*cellWidthY);
            }
        }
        break;
    case GridNodeLocation::CORNER:
        data = new vec2[(_nbCellsX + 1)*(_nbCellsY + 1)];
        for (unsigned int i = 0; i <= _nbCellsX; i++) {
            for (unsigned int j = 0; j <= _nbCellsY; j++) {
                data[(_nbCellsY + 1)*j + i] = vec2(_boundXMin + i * cellWidthX,
                    _boundYMin + j * cellWidthY);
            }
        }
        break;
    }

    metadata.dataPointer = data;
    return metadata;
}





//
////==============================================================================
//// PLUGS ACTIONS
////==============================================================================
//
//// copies and returns a linear buffer containing the grid positions.
//Metadata1DCpu VectorField2D::out_gridNodeLocationsCpu_getData() {
//
//    _vectors.dataCpu.refresh();
//
//    Metadata1DCpu metadata;
//    vec2* data;
//    float cellWidthX = (_boundXMax - _boundXMin) / _nbCellsX;
//    float cellWidthY = (_boundYMax - _boundYMin) / _nbCellsY;
//
//    switch (_gridNodeLocation) {
//    case GridNodeLocation::CENTER:
//        data = new vec2[_nbCellsX*_nbCellsY];
//        for (unsigned int i = 0; i < _nbCellsX; i++) {
//            for (unsigned int j = 0; j < _nbCellsY; j++) {
//                data[_nbCellsY*j + i] = vec2(_boundXMin + (i + 0.5f)*cellWidthX,
//                    _boundYMin + (j + 0.5f)*cellWidthY);
//            }
//        }
//        break;
//    case GridNodeLocation::CORNER:
//        data = new vec2[(_nbCellsX + 1)*(_nbCellsY + 1)];
//        for (unsigned int i = 0; i <= _nbCellsX; i++) {
//            for (unsigned int j = 0; j <= _nbCellsY; j++) {
//                data[(_nbCellsY + 1)*j + i] = vec2(_boundXMin + i * cellWidthX,
//                    _boundYMin + j * cellWidthY);
//            }
//        }
//        break;
//    }
//
//    metadata.dataPointer = data;
//    return metadata;
//}
//void VectorField2D::out_gridNodeLocationsCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//GLuint VectorField2D::out_gridNodeLocationsBuffer_getData() {
//    //TODO: generate grid node postions.
//    return 0;
//}
//void VectorField2D::out_gridNodeLocationsBuffer_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//// creates and return a linear copy of the vector data.
//Metadata1DCpu VectorField2D::out_vectorValuesCpu_getData() {
//
//    Metadata1DCpu metadata;
//    vec2* data;
//
//    switch (_gridNodeLocation) {
//    case GridNodeLocation::CENTER:
//        data = new vec2[_nbCellsX*_nbCellsY];
//        break;
//    case GridNodeLocation::CORNER:
//        data = new vec2[(_nbCellsX + 1)*(_nbCellsY + 1)];
//        break;
//    }
//
//    //    _vectors.sourceStorageType = DataBuffer2D<vec2>::StorageType::TEXTURE2D;
//    _vectors.dataCpu.dirtyItselfAndDependents();
//    memcpy(data, _vectors.dataCpu.get(), _vectors.dataSizeInBytes());
//
//    metadata.dataPointer = data;
//    return metadata;
//}
//void VectorField2D::out_vectorValuesCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//Metadata2DCpu VectorField2D::out_divergenceCpu_getData() {
//    Metadata2DCpu metadata;
//    float* data;
//
//    float hx = (_boundXMax - _boundXMin) / _nbCellsX;
//    float hy = (_boundYMax - _boundYMin) / _nbCellsY;
//    vec2 vecHx = vec2(hx, 0);
//    vec2 vecHy = vec2(0, hy);
//
//    switch (_gridNodeLocation) {
//    case GridNodeLocation::CENTER:
//        data = new float[_nbCellsX*_nbCellsY];
//
//        for (unsigned int i = 0; i < _nbCellsX; ++i) {
//            for (unsigned int j = 0; j < _nbCellsY; ++j) {
//                vec2 c = vec2(_boundXMin + (i + 0.5)*hx, _boundYMin + (j + 0.5)*hy);
//                data[_nbCellsY*j + i] = 0.5f / hx * (interp(c + vecHx).x - interp(c - vecHx).x) +
//                    0.5f / hy * (interp(c + vecHy).y - interp(c - vecHy).y);
//                //            data[_nbCellsY*j + i] *= 100;
//                //            data[_nbCellsY*j + i] = 0.5;
//            }
//        }
//
//        break;
//    case GridNodeLocation::CORNER:
//        data = new float[(_nbCellsX + 1)*(_nbCellsY + 1)];
//
//        for (unsigned int i = 0; i < _nbCellsX + 1; ++i) {
//            for (unsigned int j = 0; j < _nbCellsY + 1; ++j) {
//                vec2 c = vec2(_boundXMin + i * hx, _boundYMin + j * hy);
//                data[(_nbCellsY + 1)*j + i] = 0.5f / hx * (interp(c + vecHx).x - interp(c - vecHx).x) +
//                    0.5f / hy * (interp(c + vecHy).y - interp(c - vecHy).y);
//                //            data[(_nbCellsY+1)*j + i] = 1;
//            }
//        }
//
//        break;
//    }
//
//    //    memcpy(data, _vectors.dataCpu.get(), _vectors.dataSizeInBytes());
//
//
//    metadata.dataPointer = data;
//    return metadata;
//}
//void VectorField2D::out_divergenceCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//Metadata2DTexture2D VectorField2D::out_vectorsMetadataTexture2D_getData()
//{
//    return _vectors.out_metadataTexture2D.getData();
//}
//void VectorField2D::out_vectorsMetadataTexture2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//Metadata2DImage2D VectorField2D::out_vectorsMetadataImage2D_getData(unsigned int level)
//{
//    return _vectors.out_metadataImage2D.getData(level);
//}
//void VectorField2D::out_vectorsMetadataImage2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//void VectorField2D::in_vectorsMetadataImage2D_dirtyDependents()
//{
//    out_int_vectorsMetadataImage2D.dirtyConnections();
//}
//void VectorField2D::in_vectorsMetadataImage2D_receive(Metadata2DImage2D data)
//{
//    out_int_vectorsMetadataImage2D.pushUpdate(data);
//}
//void VectorField2D::in_vectorsMetadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//Metadata2DImage2D VectorField2D::out_int_vectorsMetadataImage2D_getData()
//{
//    return in_vectorsMetadataImage2D.pullUpdate();
//}
//void VectorField2D::out_int_vectorsMetadataImage2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//void VectorField2D::in_int_vectorsMetadataImage2D_dirtyDependents()
//{
//    out_int_vectorsMetadataImage2D.dirtyConnections();
//}
//void VectorField2D::in_int_vectorsMetadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
////<<goglu>> {
////    Beacon VectorField2D_dirtyblesDefs;
////} <<goglu>>;
//
//GOGLU_BEGIN(
//    component_definitions VectorField2D
//)GOGLU_END
//#include "../../gogluGeneratedCode/src/dataStructures/VectorField2D.cpp_d/snippet1.goglu"
