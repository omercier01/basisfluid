// grid with world space dimensions storing any type of data

#ifndef GRIDDATA2D_H
#define GRIDDATA2D_H

//#include "glm/glm.hpp"
#include "../DataBuffer/DataBuffer2D.h"



namespace goglu {

template <class T>
class GridData2D
{
public:
    enum class GridNodeLocation{ CENTER, CORNER };
public:
//private:
    DataBuffer2D<T> mData;
//    DataBuffer2D<glm::vec2> mNodePositions;
    float mBoundXMin, mBoundXMax, mBoundYMin, mBoundYMax;
    unsigned int mNbCellsX, mNbCellsY;
    GridNodeLocation mGridNodeLocation;
public:
    GridData2D(float boundXMin, float boundXMax, float boundYMin,
             float boundYMax, unsigned int nbCellsX, unsigned int nbCellsY,
             GridNodeLocation gridNodeLocation = GridNodeLocation::CENTER);

    void setCpuData(unsigned int i, unsigned int j, T data);
    T getCpuData(unsigned int i, unsigned int j);
    T getCpuData_noRefresh(unsigned int i, unsigned int j);
    
    void createCpuStorage();
    unsigned int nbElementsX();
    unsigned int nbElementsY();
    
    glm::uvec2 pointToClosestIndex(vec2 point);
    glm::vec2 indexToPosition(uvec2 index);
};

} // namespace goglu.

// include definitions because the class is templated.
#include "GridData2D.cpp"

#endif // GRIDDATA2D_H
