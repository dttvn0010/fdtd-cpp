#include "array.h"
#include "array4.h"
#include "base.h"
#include <string>
#ifndef __GRID_ITEMS_H__

class Boundary;
class Source;
class Detector;
class Obstacle;

// The Grid Object : a 3D xyz coordinates for simulation
class Grid: public Object
{
    // distance between 2 grid nodes
    double grid_spacing;
    
    // number of dimensions of the grid
    int D = 0;
    
    // courant number, defaults to the inverse of square root of the number of dimension
    double courant_number;
    
    // timestep of the simulation
    double time_step;
    
    // time passes since start of the simulation
    
    int time_steps_passed;
    
    // electric sources in the grid
    List<Source*> sources;
    
    // boundaries of the grid
    List<Boundary*> boundaries;
    
    // detectors in the grid
    List<Detector*> detectors;
    
    // obstacles in the grid
    List<Obstacle*> obstacles;
    
    // Update electric field at each node in the grid after each simulation step
    void update_E();
    
    // Update magnetic field at each node in the grid after each simulation step
    void update_H();
public:    
    // size of the grid (in nodes) in xyz dimension
    int Nx = 0;
    int Ny = 0;
    int Nz = 0;

    // Tensor to store electric field at each node of the grid
    Array4<double> E;
    
    // Tensor to store magnetic field at each node of the grid
    Array4<double> H;
    
    // Tensor to store the inverse of permittivity at each node of the grid
    Array4<double> inverse_permittivity;
    
    // Tensor to store the inverse of permeability at each node of the grid
    Array4<double> inverse_permeability;

    // Grid constructor
    Grid(
        Tuple<double, double, double> shape,
        double grid_spacing=155e-9,
        double permittivity=1.0,
        double permeability=1.0,
        Optional<double> courant_number_=Optional<double>::none()
    );
    
    // Add a boundary to the grid
    void add_boundary(Boundary* boundary, Slice x, Slice y, Slice z);
    
    // Add a source to the grid
    void add_source(Source* source, Slice x, Slice y, Slice z);
    
    // Add a detector to the grid
    void add_detector(Detector* detector, Slice x, Slice y, Slice z);
    
    // Add an obstacle to the grid
    void add_obstacle(Obstacle* obstacle, Slice x, Slice y, Slice z);
    
    // the shape of the grid in xyz dimension
    Tuple<int, int, int> shape();
    
    // time passed since the start of the simulation
    double time_passed();
    
    // Run the grid simulation
    void run(int total_steps);
    
    // Reset the grid
    void reset();
    
    // Get the courant number of the grid
    double get_courant_number();
    
    // Convert time in seconds to the number of simulation step
    int handle_time(double t);
    
    // time passed since the start of the simulation
    int get_time_steps_passed();
    ~Grid();
};

// GridItem: Abstract Item in the Grid
class GridItem: public Object
{
protected:
    // The parent grid
    Ptr<Grid> grid;
    
    // Coordinates of the item
    Slice x, y, z;
    
    // Name of the item
    std::string name;
public:

    // Add the item to the grid
    virtual void register_grid(Grid& grid, Slice x, Slice y, Slice z);
};

// Boundary : used for representing the boundaries of the grid
class Boundary: public GridItem
{
public:
    // Boundary constructor
    Boundary(const std::string& name);
    
    // update electric field derivative at the boundary
    virtual void update_phi_E();
    
    // update magnetic field derivative at the boundary
    virtual void update_phi_H();
    
    // update electric field at the boundary
    virtual void update_E();
    
    // update magnetic field at the boundary
    virtual void update_H();
    
    // Destructor
    virtual ~Boundary() = default;
};

// PML: Perfectly matched layer boundary
class PML: public Boundary
{
    // Layer parameters
    double k;
    double a;
    int thickness;
    
    // type: PMLXlow / PMLXHigh/ PMLYlow / PMLYHigh/ PMLZlow / PMLZHigh
    std::string type;
    
    // Shape of the PML in xyz dimension
    Tuple<int, int, int> shape = Tuple<int, int, int>(0,0,0);
    
    // Tensor to store electric & magnetic field
    Array4<double> phi_E, phi_H, psi_Ex, psi_Ey, psi_Ez, Hx, Hy, Hz, psi_Hx, psi_Hy, psi_Hz, bE, cE, bH, cH, sigmaE, sigmaH;

    // Location of the PML
    Tuple<Slice, Slice, Slice, Slice> loc = Tuple<Slice, Slice, Slice, Slice>(Slice::all(), Slice::all(), Slice::all(), Slice::all());
    Tuple<Slice, Slice, Slice, int> locx = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);
    Tuple<Slice, Slice, Slice, int> locy = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);
    Tuple<Slice, Slice, Slice, int> locz = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);

    // Calculate the PML' parameters from thickness
    void _calculate_parameters(int thickness=10);
    
    // Calculate the PML' location
    void _set_locations();
    
    // Calculate the PML' shape
    void _set_shape();
    
    // Calculate the PML' electric conductivity
    void _set_sigmaE();
    
    // Calculate the PML' magnetic conductivity
    void _set_sigmaH();
public:
    // Constructor
    PML(const std::string& name="", double a=1e-8);
    
    // Add the PML to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    
    // Calculate the conductivity
    Array<double> get_sigma(const Array<double>& vect);
    
    // Update electric field
    void update_E() override;
    
    // Update magnetic field
    void update_H() override;
    
    // Update electric field derivative
    void update_phi_E() override;
    
    // Update magnetic field derivative
    void update_phi_H() override;
};

// Periodic Boundary
class PeriodicBoundary: public Boundary
{
    // type : PeriodicBoundaryX/ PeriodicBoundaryY/ PeriodicBoundaryZ
    std::string type;
public:
    // Constructor
    PeriodicBoundary(const std::string& name="");
    
    // Add the item to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    
    // Update electric field
    void update_E() override;
    
    // Update magnetic field
    void update_H() override;
    
    // Update electric field derivative
    void update_phi_E() override;
    
    // Update magnetic field derivative
    void update_phi_H() override;
};

// Detector: collect electric & magnetic field at certain locations in the grid
class Detector : public GridItem
{
public:
    // Constructor
    Detector(const std::string& name);
    
    // Collect electric field
    virtual void detect_E();
    
    // Collect magnetic field
    virtual void detect_H();
    
    // Destructor
    virtual ~Detector() = default;
};

// LineDetector: detect electric & magnetic field in a line
class LineDetector : public Detector
{
    // List of coordinates of the locations to be detected
    List<int> x;
    List<int> y;
    List<int> z;
    
    // Tensor to store detected electric field values
    List<Array2<double>> E;
    
    // Tensor to store detected magnetic field values
    List<Array2<double>> H;
public:
    // Constructor
    LineDetector(const std::string& name);
    
    // Add the detector to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    
    // Collect electric field
    void detect_E() override;
    
    // Collect magnetic field
    void detect_H() override;
};


// Electric source
class Source : public GridItem
{
public:
    // Constructor
    Source(const std::string& name);
    
    // Update electric field
    virtual void update_E();
    
    // Update magnetic field
    virtual void update_H();
    
    // Destructor
    virtual ~Source() = default;
};

// Electric Source as a line
class LineSource: public Source
{
    // List of coordinates of the source
    List<int> x;
    List<int> y;
    List<int> z;
    
    // period of the source (in seconds)
    double period;
    
    // period of the source in time step
    int period_n;
    
    // Frequency of the source
    double frequency;
    
    // Amplitude of the source
    double amplitude;
    
    // Phase offset of the source
    double phase_shift;
    
    // Hanning pulse or continuos waveform
    bool pulse;
    
    // Cycle for Hanning pulse
    int cycle;
    
    // timestep used for Hanning pulse
    double hanning_dt;
    
    // Gaussian profile
    Array<double> profile;
public:

    // Constructor
    LineSource(const std::string& name, double period, double amplitude=1.0, double phase_shift=0.0, bool pulse=false, int cycle=5, double hanning_dt=10.0);
    
    // Add the line source to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    
    // Update electric field
    void update_E() override;
    
    // Update magnetic field
    void update_H() override;
};

// Electric Source as a point
class PointSource: public Source
{
    // Coordinates of the source
    int x, y, z;
    
    // period of the source (in seconds)
    double period;
    
    // period of the source in time step
    int period_n;
    
    // Frequency of the source
    double frequency;
    
    // Amplitude of the source
    double amplitude;
    
    // Phase offset of the source
    double phase_shift;
    
    // Hanning pulse or continuos waveform
    bool pulse;
    
    // Cycle for Hanning pulse
    int cycle;
    
    // timestep used for Hanning pulse
    double hanning_dt;
public:
    // Constructor
    PointSource(const std::string& name, double period, double amplitude=1.0, double phase_shift=0.0, bool pulse=false, int cycle=5, double hanning_dt=10.0);
    
    // Add the source to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    
    // Update electric field
    void update_E() override;
    
    // Update magnetic field
    void update_H() override;

};

// An obstacle in the grid
class Obstacle : public GridItem
{   
public:
    // Constructor
    Obstacle(const std::string& name);
    
    // Update electric field
    virtual void update_E(const Array4<double>& curl_H);
    
    // Update magnetic field
    virtual void update_H(const Array4<double>& curl_E);
    
    // Add the obstacle to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z);
    
    // Destructor
    virtual ~Obstacle() = default;
};

// Anisotropic Obstacle
class AnisotropicObstacle: public Obstacle
{
    // Relative permittivity
    double permittivity;
    
    // Inverse of permittivity
    Array3<double> inverse_permittivity;
    
    // Size of the obstacle in grid_spacing
    int Nx, Ny, Nz;
public:
    // Constructor
    AnisotropicObstacle(const std::string& name,  double permittivity);
    
    // Update electric field
    void update_E(const Array4<double>& curl_H) override;
    
    // Update magnetic field
    void update_H(const Array4<double>& curl_E) override;
    
    // Add the obstacle to the grid
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
};
#endif