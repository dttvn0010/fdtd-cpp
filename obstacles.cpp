#include "np.h"
#include "grid_items.h"

// Obstacle constructor
// Input:
//  + name: name of the obstacle
Obstacle::Obstacle(const std::string& name)
{
    this->name = name;
}

// Update electric field at the obstacle
// This is a virtual function, the real implemenation in inside each children class
void Obstacle::update_E(const Array4<double>& curl_E)
{
}

// Update magnetic field at the obstacle
// This is a virtual function, the real implemenation in inside each children class
void Obstacle::update_H(const Array4<double>& curl_E)
{

}

// Add the obstacle to the grid
// This is a virtual function, the real implemenation in inside each children class
void Obstacle::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{

}

// ===========================================================================================================================================
//AnisotropicObstacle constructor
// Input:
//   + name: name of the obstacle
//   + permittivity: relative permittivity at the obstacle
AnisotropicObstacle::AnisotropicObstacle(const std::string& name,  double permittivity):
    Obstacle(name)
{
    this->permittivity = permittivity;
}

// AnisotropicObstacle - register grid: add obstacle to the  grid
// Input:
//   + grid: the system grid
//   + x: x-range of the
//   + y: y-range of the
//   + z: z-range of the
void AnisotropicObstacle::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    this->x = x;
    this->y = y;
    this->z = z;
    
    // Calculate obstacle' size 3 dimensions
    Nx = x.stop.value - x.start.value;
    Ny = y.stop.value - y.start.value;
    Nz = z.stop.value - z.start.value;
    
    // Calculate inverse of permittivity
    auto inverse_permittivity = Array4<double>::ones(Tuple<int, int, int, int>(Nx, Ny, Nz, 3)) / permittivity;

    auto all = Slice::all();

    // set the permittivity values of the object at its border to be equal to the grid permittivity
    if(Nx > 1)
    {
        inverse_permittivity.get(-1, all, all, 0).copy_from(grid.inverse_permittivity.get(-1, y, z, 0));
    }

    if(Ny > 1)
    {
        inverse_permittivity.get(all, -1, all, 1).copy_from(grid.inverse_permittivity.get(x, -1, z, 1));
    }

    if(Nz > 1)
    {
        inverse_permittivity.get(all, all, -1, 2).copy_from(grid.inverse_permittivity.get(x,y,-1,2));
    }

    grid.inverse_permittivity.get(x, y, z, all).set_all(0);

    int N = Nx * Ny * Nz;
    auto eye = Array3<double>::zeros(Tuple<int, int, int> (N, 3, 3));

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < 3; j++) eye.at(i, j, j) = 1.0;
    }

    // reshape inverse permittivity for easier latter calculations
    this->inverse_permittivity = np::reshape<double>(
        inverse_permittivity, Tuple<int, int, int>(N, 1, 3)
    ) * eye;
}

// Batch matrix multiply two arrays ~ Einstein sum: ijk,ikl->ijl
Array3<double> bmm(const Array3<double>& arr1, const Array3<double>& arr2)
{
    auto shape1 = arr1.shape();
    auto shape2 = arr2.shape();
    assert(shape1.item1 == shape2.item1 && shape1.item3 == shape2.item2);

    auto res = Array3<double>::zeros(Tuple<int, int, int> (shape1.item1, shape1.item2, shape2.item3));
    for(int i = 0 ; i < shape1.item1; i++)
    {
        for(int j = 0; j < shape1.item2; j++)
        {
            for(int l = 0; l < shape2.item3; l++)
            {
                double s = 0.0;
                for(int k = 0; k < shape1.item3; k++)
                {
                    s += arr1.at(i, j, k) * arr2.at(i, k, l);
                }
                res.at(i,j,l) = s;
            }
        }
    }
    return res;
}

// Update electric field at the obstacle
void AnisotropicObstacle::update_E(const Array4<double>& curl_H)
{
    auto all = Slice::all();
    auto shape3 = Tuple<int,int,int>(Nx*Ny*Nz, 3, 1); 
    auto shape4 = Tuple<int,int,int,int>(Nx, Ny, Nz, 3);

    // ∂E/∂t = ∇xH/ε
    grid->E.get(x,y,z,all) += np::reshape<double>(
        bmm(inverse_permittivity, np::reshape<double>(curl_H.get(x,y,z,all), shape3)),
        shape4
    );

}

// Update magnetic field at the obstacle
void AnisotropicObstacle::update_H(const Array4<double>& curl_E)
{

}