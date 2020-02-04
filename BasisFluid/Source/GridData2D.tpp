
#include "GridData2D.h"

using namespace glm;

template <class T>
GridData2D<T>::GridData2D(float parBoundXMin, float parBoundXMax, float parBoundYMin,
                          float parBoundYMax, unsigned int nbCellsX,
                          unsigned int nbCellsY,
                          GridNodeLocation gridNodeLocation) :
    mNbCellsX(nbCellsX),
    mNbCellsY(nbCellsY),
    mData(1, 1),
    mGridNodeLocation(gridNodeLocation)
{
    mBoundXMin = parBoundXMin;
    mBoundXMax = parBoundXMax;
    mBoundYMin = parBoundYMin;
    mBoundYMax = parBoundYMax;

    // resizes, but does not necessarily allocate the data.
    switch(gridNodeLocation) {
    case GridNodeLocation::CENTER:
        mData.resize(nbCellsX, nbCellsY);
        break;
    case GridNodeLocation::CORNER:
        mData.resize(nbCellsX+1, nbCellsY+1);
        break;
    }

}


template <class T>
void GridData2D<T>::setCpuData(unsigned int i, unsigned int j, T data)
{
    mData.setCpuData(i, j, data);
}

template <class T>
T GridData2D<T>::getCpuData(unsigned int i, unsigned int j)
{
    return mData.getCpuData(i, j);
}

template <class T>
T GridData2D<T>::getCpuData_noRefresh(unsigned int i, unsigned int j)
{
    return mData.getCpuData_noRefresh(i, j);
}

template <class T>
void GridData2D<T>::createCpuStorage()
{
    mData.createCpuStorage();
}

template <class T>
unsigned int GridData2D<T>::nbElementsX()
{
    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        return mNbCellsX;
        break;
    case GridNodeLocation::CORNER:
        return mNbCellsX + 1;
        break;
    default:
        // TODO: error, unknows grid location type.
        return 0;
    }
}

template <class T>
unsigned int GridData2D<T>::nbElementsY()
{
    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        return mNbCellsY;
        break;
    case GridNodeLocation::CORNER:
        return mNbCellsY + 1;
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

    float normalizedX = (point.x - mBoundXMin)/(mBoundXMax - mBoundXMin)*mNbCellsX;
    float normalizedY = (point.y - mBoundYMin)/(mBoundYMax - mBoundYMin)*mNbCellsY;


    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX-1;
        } else {
            index.x = clamp<unsigned int>((unsigned int)(std::floor(normalizedX)), 0, mNbCellsX-1);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY-1;
        } else {
            index.y = clamp<unsigned int>((unsigned int)(std::floor(normalizedY)), 0, mNbCellsY-1);
        }
        break;
    case GridNodeLocation::CORNER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX;
        } else {
            index.x = clamp<int>(int(floor(normalizedX + 0.5)), 0, mNbCellsX);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY;
        } else {
            index.y = clamp<int>(int(floor(normalizedY + 0.5)), 0, mNbCellsY);
        }
        break;
    }

    return index;
}


template <class T>
vec2 GridData2D<T>::indexToPosition(uvec2 index)
{
    vec2 pos;

    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER :
        pos.x = boundXMin.get() + (index.x + 0.5f)*(boundXMax.get()-boundXMin.get())/mNbCellsX;
        pos.y = boundYMin.get() + (index.y + 0.5f)*(boundYMax.get()-boundYMin.get())/mNbCellsY;
        break;
    case GridNodeLocation::CORNER :
        pos.x = boundXMin.get() + (index.x)*(boundXMax.get()-boundXMin.get())/mNbCellsX;
        pos.y = boundYMin.get() + (index.y)*(boundYMax.get()-boundYMin.get())/mNbCellsY;
        break;
    }

    return pos;
}


