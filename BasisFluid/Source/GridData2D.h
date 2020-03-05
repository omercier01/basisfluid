
// grid with world space dimensions storing any type of data. Data is at cell centers.

#ifndef GRIDDATA2D_H
#define GRIDDATA2D_H

#include "DataBuffer2D.h"


template <class T>
class GridData2D
{
public:
    DataBuffer2D<T> _data;
    float _boundXMin, _boundXMax, _boundYMin, _boundYMax;
    unsigned int _nbCellsX, _nbCellsY;
public:
    GridData2D(float boundXMin, float boundXMax, float boundYMin,
        float boundYMax, unsigned int nbCellsX, unsigned int nbCellsY);

    void setCpuData(unsigned int i, unsigned int j, T data);
    T getCpuData(unsigned int i, unsigned int j);
    T getCpuData_noRefresh(unsigned int i, unsigned int j);

    void createCpuStorage();
    unsigned int nbElementsX();
    unsigned int nbElementsY();

    glm::uvec2 pointToClosestIndex(glm::vec2 point);
    glm::vec2 indexToPosition(glm::uvec2 index);
};

// include definitions because the class is templated.
#include "GridData2D.tpp"

#endif // GRIDDATA2D_H
