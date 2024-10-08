#include "array.h"
#include "grid_items.h"

#define WAVELENGTH 1550e-9
#define SPEED_LIGHT 299792458.0

int main()
{
    // Initialize a new grid for simulation
    Grid grid(
        Tuple<double, double, double>(2.5e-5, 1.5e-5, 0),
        0.1 * WAVELENGTH,
        1.0,
        1.0
    );

    auto all = Slice::all();

    // Add x-low boundary
    grid.add_boundary(
        new PML("pml_xlow"),
        Slice(0,10),
        all,
        all
    );

    // Add x-high boundary
    grid.add_boundary(
        new PML("pml_xhigh"),
        Slice::tail(-10),
        all,
        all
    );

    // Add y-low boundary
    grid.add_boundary(
        new PML("pml_ylow"),
        all,
        Slice(0, 10),
        all
    );

    // Add y-high boundary
    grid.add_boundary(
        new PML("pml_yhigh"),
        all,
        Slice::tail(-10),
        all
    );

    // Add z boundary
    grid.add_boundary(
        new PeriodicBoundary("zbounds"),
        all,
        all,
        Slice(0, 0)
    );

    // Add a line source
    grid.add_source(
        new LineSource("linesource", WAVELENGTH / SPEED_LIGHT),
        Slice(50,55),
        Slice(70, 75),
        Slice(0, 0)
    );

    // Add a point source
    grid.add_source(
        new PointSource("pointsource", WAVELENGTH / SPEED_LIGHT),
        Slice(100, 100),
        Slice(60, 60),
        Slice(0, 0)
    );

    // Add a detector
    grid.add_detector(
        new LineDetector("detector"),
        Slice(77,77),
        all,
        Slice(0, 0)
    );

    // Add an anisotropic obstacle
    grid.add_obstacle(new AnisotropicObstacle("object", 2.5),
        Slice(11, 32),
        Slice(30,84),
        Slice(0, 1)
    );

    // Run the simulation
    grid.run(50);

    // Collect electric & magnetic strength
    auto Ez = grid.E.get(all, all, 0, 2);
    auto Hx = grid.H.get(all, all, 0, 0);
    auto Hy = grid.H.get(all, all, 0, 1);

    print("Ez=", Ez);
    print("Hx=", Hx);
    print("Hy=", Hy);

    return  0;
}