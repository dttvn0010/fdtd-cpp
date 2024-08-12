#include <math.h>
#include "grid_items.h"

const double light_speed = 299792458.0;

/* Caculate the rotar vector of E (∇xE)
  -Input: 
   + E : Electric field 
 - Output: 
   + The rotar vector of E
*/
Array4<double> curl_E(const Array4<double>& E)
{
    Array4<double> curl = Array4<double>::zeros(E.shape());
    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);

    // curl.x = ∂Ez/∂y - ∂Ey/∂z
    curl.get(all, head_1, all, 0) += E.get(all, tail_1, all, 2) - E.get(all, head_1, all, 2);
    curl.get(all, all, head_1, 0) -= E.get(all, all, tail_1, 1) - E.get(all, all, head_1, 1); 

    // curl.y = ∂Ex/∂z - ∂Ez/∂x
    curl.get(all, all, head_1, 1) += E.get(all, all, tail_1, 0) - E.get(all, all, head_1, 0);
    curl.get(head_1, all, all, 1) -= E.get(tail_1, all, all, 2) - E.get(head_1, all, all, 2);

    // curl.z = ∂Ey/∂x - ∂Ex/∂y
    curl.get(head_1, all, all, 2) += E.get(tail_1, all, all, 1) - E.get(head_1, all, all, 1);
    curl.get(all, head_1, all, 2) -= E.get(all, tail_1, all, 0) - E.get(all, head_1, all, 0);

    return curl;
}

/* Caculate the rotar vector of H (∇xH)
  -Input: 
   + E : Magnetic field 
 - Output: 
   + The rotar vector of H
*/
Array4<double> curl_H(const Array4<double>& H)
{
    Array4<double> curl = Array4<double>::zeros(H.shape());
    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);

    // curl.x = ∂Hz/∂y - ∂Hy/∂z
    curl.get(all, tail_1, all, 0) += H.get(all, tail_1, all, 2) - H.get(all, head_1, all, 2);
    curl.get(all, all, tail_1, 0) -= H.get(all, all, tail_1, 1) - H.get(all, all, head_1, 1); 

    // curl.y = ∂Hx/∂z - ∂Hz/∂x
    curl.get(all, all, tail_1, 1) += H.get(all, all, tail_1, 0) - H.get(all, all, head_1, 0);
    curl.get(tail_1, all, all, 1) -= H.get(tail_1, all, all, 2) - H.get(head_1, all, all, 2);

    // curl.z = ∂Hy/∂x - ∂Hx/∂y
    curl.get(tail_1, all, all, 2) += H.get(tail_1, all, all, 1) - H.get(head_1, all, all, 1);
    curl.get(all, tail_1, all, 2) -= H.get(all, tail_1, all, 0) - H.get(all, head_1, all, 0);

    return curl;
}

/*
Grid initialization
- Input:
  + shape : the xyz dimension of the grid
  + grid_spacing: distance between 2 grid nodes
  + permittivity: permittivity of the medium
  + permeability: permeability of the medium
  + courant_number: courant number, defaults to the inverse of square root of the number of dimension
*/
Grid::Grid(
    Tuple<double, double, double> shape,
    double grid_spacing,
    double permittivity,
    double permeability,
    Optional<double> courant_number_
)
{
    this->grid_spacing = grid_spacing;
    
    // Caculate the number of nodes for each dimension
    Nx = shape.item1 > 0 ? int(shape.item1 / grid_spacing + 0.5) : 1;
    Ny = shape.item2 > 0 ? int(shape.item2 / grid_spacing + 0.5) : 1;
    Nz = shape.item3 > 0 ? int(shape.item3 / grid_spacing + 0.5) : 1;
    
    // number of dimension
    D = int(Nx > 1) + int(Ny > 1) + int(Nz > 1);

    // max courant number set to square root of D
    double max_courant_number = pow(D, -0.5);
    if(courant_number_.is_null)
    {
        courant_number = 0.99 * max_courant_number;
    }
    else if(courant_number_.value > max_courant_number)
    {
        panic("courant_number %f too high for a %d simulation", courant_number_.value, D);
    }
    else 
    {
        courant_number = courant_number_.value;
    }

    // timestep of the simulation
    time_step = courant_number * grid_spacing / light_speed;
    
    // Tensors to store the electric & magnetic field
    E = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    H = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));

    // Tensor to store the inverse permittivity of the medium
    inverse_permittivity = Array4<double>::full(
        Tuple<int, int, int, int>(Nx, Ny, Nz, 3),
        1.0/permittivity
    );

    // Tensor to store the inverse permeability of the medium
    inverse_permeability = Array4<double>::full(
        Tuple<int, int, int, int>(Nx, Ny, Nz, 3),
        1.0/permeability
    );

    // Number of simulation step
    time_steps_passed = 0;
}

// The xyz-shape of the grid
Tuple<int, int, int> Grid::shape()
{
    return Tuple<int, int, int>(Nx, Ny, Nz);
}

// Amount of time that has passed since the start of the simulation
double Grid::time_passed()
{
    return time_steps_passed * time_step;
}

// Run the grid simulation
// Input:
//   + total_steps: number of step to run
void Grid::run(int total_steps)
{
    for(int step = 0; step < total_steps; step++)
    {
        update_E();
        update_H();
        time_steps_passed += 1;
        print("Step:", time_steps_passed);
    }
}

// Update the electric field for each simulation step
void Grid::update_E()
{
    // update electric field derivatives for each boundary
    for(auto& boundary : boundaries)
    {
        boundary->update_phi_E();
    }

    // Caculate ∇xH
    auto curl = curl_H(H);
    
    // ∂E/∂t = ∇xH/ε
    E +=  inverse_permittivity * courant_number * curl;

    // update electric field for each obstacle
    for(auto& obstacle : obstacles)
    {
        obstacle->update_E(curl);
    }

    // update electric field for each boundary
    for(auto& boundary : boundaries)
    {
        boundary->update_E();
    }

    // update electric field for each source
    for(auto& src : sources)
    {
        src->update_E();
    }

    // collect electric field at each detector
    for(auto& det : detectors)
    {
        det->detect_E();
    }
}

// Update the magnetic field for each simulation step
void Grid::update_H()
{
    // update magnetic field derivatives for each boundary
    for(auto & boundary : boundaries)
    {
        boundary->update_phi_H();
    }
    
    // Caculate ∇xE
    auto curl = curl_E(E);
    
    // ∂H/∂t = -∇xE/µ
    H -=  inverse_permeability * courant_number * curl;

    // update magnetic field for each obstacle
    for(auto& obstacle : obstacles)
    {
        obstacle->update_H(curl);
    }

    // update magnetic field for each boundary
    for(auto& boundary : boundaries)
    {
        boundary->update_H();
    }

    // update magnetic field for each source
    for(auto& src : sources)
    {
        src->update_H();
    }

    // collect magnetic field at each detector
    for(auto& det : detectors)
    {
        det->detect_H();
    }
}


// reset the grid
void Grid::reset()
{
    H *= 0.0;
    E *= 0.0;
    time_steps_passed = 0;
}

// Add a boundary to the grid
// Input:
//   + x: x-range of the boundary
//   + y: y-range of the boundary
//   + z: z-range of the boundary
void Grid::add_boundary(Boundary* boundary, Slice x, Slice y, Slice z)
{
    boundary->register_grid(*this, x, y, z);
    boundaries.append(boundary);
}


// Add a source to the grid
// Input:
//   + x: x-range of the source
//   + y: y-range of the source
//   + z: z-range of the source
void Grid::add_source(Source* source, Slice x, Slice y, Slice z)
{
    source->register_grid(*this, x, y, z);
    sources.append(source);
}

// Add a detector to the grid
// Input:
//   + x: x-range of the detector
//   + y: y-range of the detector
//   + z: z-range of the detector
void Grid::add_detector(Detector* detector, Slice x, Slice y, Slice z)
{
    detector->register_grid(*this, x, y, z);
    detectors.append(std::move(detector));
}

// Add an obstacle to the grid
// Input:
//   + x: x-range of the obstacle
//   + y: y-range of the obstacle
//   + z: z-range of the obstacle
void Grid::add_obstacle(Obstacle* obstacle, Slice x, Slice y, Slice z)
{
    obstacle->register_grid(*this, x, y, z);
    obstacles.append(std::move(obstacle));
}

// Get the courant number of the grid
double Grid::get_courant_number()
{
    return courant_number;
}

// Convert time to number of simulation step
int Grid::handle_time(double t)
{
    return int(t/time_step + 0.5);
}

// Get number of simulation step that has been executed
int Grid::get_time_steps_passed()
{
    return time_steps_passed;
}

// Grid destructor
Grid::~Grid()
{
    // Free all obstacles in memory
    for(auto& obstacle : obstacles)
    {
        delete obstacle;
    }

    // Free all boundaries in memory
    for(auto& boundary : boundaries)
    {
        delete boundary;
    }

    // Free all sources in memory
    for(auto& src : sources)
    {
        delete src;
    }

    // Free all detectors in memory
    for(auto& det : detectors)
    {
        delete det;
    }
}


// Register a grid item (boundary/source/detector/obstacle) to the grid
// This is a virtual function, the real implemenation in inside each children class
void GridItem::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    
}
