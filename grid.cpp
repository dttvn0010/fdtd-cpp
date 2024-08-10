#include <math.h>
#include "grid_items.h"

const double light_speed = 299792458.0;

Array4<double> curl_E(const Array4<double>& E)
{
    Array4<double> curl = Array4<double>::zeros(E.shape());
    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);

    // curl[:, :-1, :, 0] += E[:, 1:, :, 2] - E[:, :-1, :, 2]
    curl.get(all, head_1, all, 0) += E.get(all, tail_1, all, 2) - E.get(all, head_1, all, 2);

    // curl[:, :, :-1, 0] -= E[:, :, 1:, 1] - E[:, :, :-1, 1]
    curl.get(all, all, head_1, 0) -= E.get(all, all, tail_1, 1) - E.get(all, all, head_1, 1); 

    // curl[:, :, :-1, 1] += E[:, :, 1:, 0] - E[:, :, :-1, 0]
    curl.get(all, all, head_1, 1) += E.get(all, all, tail_1, 0) - E.get(all, all, head_1, 0);

    // curl[:-1, :, :, 1] -= E[1:, :, :, 2] - E[:-1, :, :, 2]
    curl.get(head_1, all, all, 1) -= E.get(tail_1, all, all, 2) - E.get(head_1, all, all, 2);

    // curl[:-1, :, :, 2] += E[1:, :, :, 1] - E[:-1, :, :, 1]
    curl.get(head_1, all, all, 2) += E.get(tail_1, all, all, 1) - E.get(head_1, all, all, 1);

    // curl[:, :-1, :, 2] -= E[:, 1:, :, 0] - E[:, :-1, :, 0]
    curl.get(all, head_1, all, 2) -= E.get(all, tail_1, all, 0) - E.get(all, head_1, all, 0);

    return curl;
}

Array4<double> curl_H(const Array4<double>& H)
{
    Array4<double> curl = Array4<double>::zeros(H.shape());
    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);

    //curl[:, 1:, :, 0] += H[:, 1:, :, 2] - H[:, :-1, :, 2]
    curl.get(all, tail_1, all, 0) += H.get(all, tail_1, all, 2) - H.get(all, head_1, all, 2);

    //curl[:, :, 1:, 0] -= H[:, :, 1:, 1] - H[:, :, :-1, 1]
    curl.get(all, all, tail_1, 0) -= H.get(all, all, tail_1, 1) - H.get(all, all, head_1, 1); 

    //curl[:, :, 1:, 1] += H[:, :, 1:, 0] - H[:, :, :-1, 0]
    curl.get(all, all, tail_1, 1) += H.get(all, all, tail_1, 0) - H.get(all, all, head_1, 0);

    //curl[1:, :, :, 1] -= H[1:, :, :, 2] - H[:-1, :, :, 2]
    curl.get(tail_1, all, all, 1) -= H.get(tail_1, all, all, 2) - H.get(head_1, all, all, 2);

    //curl[1:, :, :, 2] += H[1:, :, :, 1] - H[:-1, :, :, 1]
    curl.get(tail_1, all, all, 2) += H.get(tail_1, all, all, 1) - H.get(head_1, all, all, 1);

    //curl[:, 1:, :, 2] -= H[:, 1:, :, 0] - H[:, :-1, :, 0]
    curl.get(all, tail_1, all, 2) -= H.get(all, tail_1, all, 0) - H.get(all, head_1, all, 0);

    return curl;
}

Grid::Grid(
    Tuple<double, double, double> shape,
    double grid_spacing,
    double permittivity,
    double permeability,
    Optional<double> courant_number_
)
{
    this->grid_spacing = grid_spacing;
    Nx = shape.item1 > 0 ? int(shape.item1 / grid_spacing + 0.5) : 1;
    Ny = shape.item2 > 0 ? int(shape.item2 / grid_spacing + 0.5) : 1;
    Nz = shape.item3 > 0 ? int(shape.item3 / grid_spacing + 0.5) : 1;
    D = int(Nx > 1) + int(Ny > 1) + int(Nz > 1);

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

    time_step = courant_number * grid_spacing / light_speed;
    E = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    H = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));

    inverse_permittivity = Array4<double>::full(
        Tuple<int, int, int, int>(Nx, Ny, Nz, 3),
        1.0/permittivity
    );

    inverse_permeability = Array4<double>::full(
        Tuple<int, int, int, int>(Nx, Ny, Nz, 3),
        1.0/permeability
    );

    time_steps_passed = 0;
}

Tuple<int, int, int> Grid::shape()
{
    return Tuple<int, int, int>(Nx, Ny, Nz);
}

double Grid::time_passed()
{
    return time_steps_passed * time_step;
}

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

void Grid::update_E()
{

    for(auto& boundary : boundaries)
    {
        boundary->update_phi_E();
    }

    auto curl = curl_H(H);
    E +=  inverse_permittivity * courant_number * curl;

    for(auto& obstacle : obstacles)
    {
        obstacle->update_E(curl);
    }

    for(auto& boundary : boundaries)
    {
        boundary->update_E();
    }

    int idx = 0;
    for(auto& src : sources)
    {
        src->update_E();
        idx += 1;
    }

    for(auto& det : detectors)
    {
        det->detect_E();
    }
}

void Grid::update_H()
{
    for(auto & boundary : boundaries)
    {
        boundary->update_phi_H();
    }
    auto curl = curl_E(E);
    H -=  inverse_permeability * courant_number * curl;

    for(auto& obstacle : obstacles)
    {
        obstacle->update_H(curl);
    }

    for(auto& boundary : boundaries)
    {
        boundary->update_H();
    }

    for(auto& src : sources)
    {
        src->update_H();
    }

    for(auto& det : detectors)
    {
        det->detect_H();
    }
}

void Grid::reset()
{
    H *= 0.0;
    E *= 0.0;
    time_steps_passed = 0;
}

void Grid::add_boundary(Boundary* boundary, Slice x, Slice y, Slice z)
{
    boundary->register_grid(*this, x, y, z);
    boundaries.append(std::move(boundary));
}

void Grid::add_source(Source* source, Slice x, Slice y, Slice z)
{
    source->register_grid(*this, x, y, z);
    sources.append(source);
}

void Grid::add_detector(Detector* detector, Slice x, Slice y, Slice z)
{
    detector->register_grid(*this, x, y, z);
    detectors.append(std::move(detector));
}

void Grid::add_obstacle(Obstacle* obstacle, Slice x, Slice y, Slice z)
{
    obstacle->register_grid(*this, x, y, z);
    obstacles.append(std::move(obstacle));
}

double Grid::get_courant_number()
{
    return courant_number;
}

int Grid::handle_time(double t)
{
    return int(t/time_step + 0.5);
}

int Grid::get_time_steps_passed()
{
    return time_steps_passed;
}

void GridItem::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    
}

Grid::~Grid()
{
    for(auto& obstacle : obstacles)
    {
        delete obstacle;
    }

    for(auto& boundary : boundaries)
    {
        delete boundary;
    }

    for(auto& src : sources)
    {
        delete src;
    }

    for(auto& det : detectors)
    {
        delete det;
    }
}