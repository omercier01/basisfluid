
#include "GridData2D.h"

using namespace glm;

template <class T>
GridData2D<T>::GridData2D(float parBoundXMin, float parBoundXMax, float parBoundYMin,
                          float parBoundYMax, unsigned int nbCellsX,
                          unsigned int nbCellsY,
                          GridNodeLocation gridNodeLocation) :
    _nbCellsX(nbCellsX),
    _nbCellsY(nbCellsY),
    _data(1, 1),
    _gridNodeLocation(gridNodeLocation)
{
    _boundXMin = parBoundXMin;
    _boundXMax = parBoundXMax;
    _boundYMin = parBoundYMin;
    _boundYMax = parBoundYMax;

    // resizes, but does not necessarily allocate the data.
    switch(gridNodeLocation) {
    case GridNodeLocation::CENTER:
        _data.resize(nbCellsX, nbCellsY);
        break;
    case GridNodeLocation::CORNER:
        _data.resize(nbCellsX+1, nbCellsY+1);
        break;
    }

}


template <class T>
void GridData2D<T>::setCpuData(unsigned int i, unsigned int j, T data)
{
    _data.setCpuData(i, j, data);
}

template <class T>
T GridData2D<T>::getCpuData(unsigned int i, unsigned int j)
{
    return _data.getCpuData(i, j);
}

template <class T>
T GridData2D<T>::getCpuData_noRefresh(unsigned int i, unsigned int j)
{
    return _data.getCpuData_noRefresh(i, j);
}

template <class T>
void GridData2D<T>::createCpuStorage()
{
    _data.createCpuStorage();
}

template <class T>
unsigned int GridData2D<T>::nbElementsX()
{
    switch(_gridNodeLocation) {
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

template <class T>
unsigned int GridData2D<T>::nbElementsY()
{
    switch(_gridNodeLocation) {
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

template <class T>
uvec2 GridData2D<T>::pointToClosestIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - _boundXMin)/(_boundXMax - _boundXMin)*_nbCellsX;
    float normalizedY = (point.y - _boundYMin)/(_boundYMax - _boundYMin)*_nbCellsY;


    switch(_gridNodeLocation) {
    case GridNodeLocation::CENTER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= _nbCellsX) {
            index.x = _nbCellsX-1;
        } else {
            index.x = clamp<unsigned int>((unsigned int)(std::floor(normalizedX)), 0, _nbCellsX-1);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= _nbCellsY) {
            index.y = _nbCellsY-1;
        } else {
            index.y = clamp<unsigned int>((unsigned int)(std::floor(normalizedY)), 0, _nbCellsY-1);
        }
        break;
    case GridNodeLocation::CORNER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= _nbCellsX) {
            index.x = _nbCellsX;
        } else {
            index.x = clamp<int>(int(floor(normalizedX + 0.5)), 0, _nbCellsX);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= _nbCellsY) {
            index.y = _nbCellsY;
        } else {
            index.y = clamp<int>(int(floor(normalizedY + 0.5)), 0, _nbCellsY);
        }
        break;
    }

    return index;
}


template <class T>
vec2 GridData2D<T>::indexToPosition(uvec2 index)
{
    vec2 pos;

    switch(_gridNodeLocation) {
    case GridNodeLocation::CENTER :
        pos.x = _boundXMin + (index.x + 0.5f)*(_boundXMax-_boundXMin)/_nbCellsX;
        pos.y = _boundYMin + (index.y + 0.5f)*(_boundYMax-_boundYMin)/_nbCellsY;
        break;
    case GridNodeLocation::CORNER :
        pos.x = _boundXMin + (index.x)*(_boundXMax-_boundXMin)/_nbCellsX;
        pos.y = _boundYMin + (index.y)*(_boundYMax-_boundYMin)/_nbCellsY;
        break;
    }

    return pos;
}


