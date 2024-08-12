#include "grid_items.h"

// Detector constructor
// Input:
//  + name: name of the detector
Detector::Detector(const std::string& name)
{
    this->name = name;
}

// Collect electric field at the location of the detector
// This is a virtual function, the real implemenation in inside each children class
void Detector::detect_E()
{

}

// Collect magnetic field at the location of the detector
// This is a virtual function, the real implemenation in inside each children class
void Detector::detect_H()
{

}

// ===========================================================================
// LineDetector constructor
// Input:
//  + name: name of the detector
LineDetector::LineDetector(const std::string& name) : Detector(name)
{

}

// LineDetector - register grid: add detector to the  grid
// Input:
//   + grid: the system grid
//   + x: x-range of the
//   + y: y-range of the
//   + z: z-range of the
void LineDetector::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;

    // resolving the start/end of location in x axis
    int x0 = x.start.is_null ? 0 : x.start.value;
    int x1 = x.stop.is_null ? grid.Nx : x.stop.value;

    // resolving the start/end of location in y axis
    int y0 = y.start.is_null ? 0 : y.start.value;
    int y1 = y.stop.is_null ? grid.Ny : y.stop.value;

    // resolving the start/end of location in z axis
    int z0 = z.start.is_null ? 0 : z.start.value;
    int z1 = z.stop.is_null ? grid.Nz : z.stop.value;

    // m = max(x1-x0, y1-y0, z1-z0)
    int m = x1-x0;
    if(y1-y0 > m) m = y1-y0;
    if(z1-z0 > m) m = z1-z0;

    // Generate list of x coordinates
    if(m == x1-x0)
    {
        for(int i = x0; i < x1; i ++) this->x.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->x.append(x0);
    }

    // Generate list of y coordinates
    if(m == y1-y0)
    {
        for(int i = y0; i < y1; i ++) this->y.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->y.append(y0);
    }

    // Generate list of z coordinates
    if(m == z1-z0)
    {
        for(int i = z0; i < z1; i ++) this->x.append(i);
    }
    else 
    {
        for(int i = 0; i < m; i++) this->z.append(z0);
    }
}

// Collect electric field at the location of the line detector
void LineDetector::detect_E()
{
    E.append(grid->E.get(x, y, z));
}

// Collect magnetic field at the location of the line detector
void LineDetector::detect_H()
{
    H.append(grid->H.get(x, y, z));
}