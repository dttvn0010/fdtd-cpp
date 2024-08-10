#include "grid_items.h"

Detector::Detector(const std::string& name)
{
    this->name = name;
}

void Detector::detect_E()
{

}

void Detector::detect_H()
{

}

// ===========================================================================

LineDetector::LineDetector(const std::string& name) : Detector(name)
{

}


void LineDetector::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    int x0 = x.start.is_null ? 0 : x.start.value;
    int x1 = x.stop.is_null ? grid.Nx : x.stop.value;

    int y0 = y.start.is_null ? 0 : y.start.value;
    int y1 = y.stop.is_null ? grid.Nx : y.stop.value;

    int z0 = z.start.is_null ? 0 : z.start.value;
    int z1 = z.stop.is_null ? grid.Nx : z.stop.value;

    int m = x1-x0;
    if(y1-y0 > m) m = y1-y0;
    if(z1-z0 > m) m = z1-z0;

    if(m == x1-x0)
    {
        for(int i = x0; i < x1; i ++) this->x.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->x.append(x0);
    }

    if(m == y1-y0)
    {
        for(int i = y0; i < y1; i ++) this->y.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->y.append(y0);
    }

    if(m == z1-z0)
    {
        for(int i = z0; i < z1; i ++) this->x.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->z.append(z0);
    }
}

void LineDetector::detect_E()
{
    E.append(grid->E.get(x, y, z));
}

void LineDetector::detect_H()
{
    H.append(grid->H.get(x, y, z));
}