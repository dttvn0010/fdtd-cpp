#include "array.h"
#include "array4.h"
#include "base.h"
#include <string>
#ifndef __GRID_ITEMS_H__

class Boundary;
class Source;
class Detector;
class Obstacle;

class Grid: public Object
{
    double grid_spacing;
    int D = 0;
    double courant_number;
    double time_step;
    int time_steps_passed;
    List<Source*> sources;
    List<Boundary*> boundaries;
    List<Detector*> detectors;
    List<Obstacle*> obstacles;
    void update_E();
    void update_H();
public:    
    int Nx = 0;
    int Ny = 0;
    int Nz = 0;

    Array4<double> E;
    Array4<double> H;
    Array4<double> inverse_permittivity;
    Array4<double> inverse_permeability;

    Grid(
        Tuple<double, double, double> shape,
        double grid_spacing=155e-9,
        double permittivity=1.0,
        double permeability=1.0,
        Optional<double> courant_number_=Optional<double>::none()
    );
    void add_boundary(Boundary* boundary, Slice x, Slice y, Slice z);
    void add_source(Source* source, Slice x, Slice y, Slice z);
    void add_detector(Detector* detector, Slice x, Slice y, Slice z);
    void add_obstacle(Obstacle* obstacle, Slice x, Slice y, Slice z);
    Tuple<int, int, int> shape();
    double time_passed();
    void run(int total_steps);
    void reset();
    double get_courant_number();
    int handle_time(double t);
    int get_time_steps_passed();
    ~Grid();
};


class GridItem: public Object
{
protected:
    Ptr<Grid> grid;
    Slice x, y, z;
    std::string name;
public:
    virtual void register_grid(Grid& grid, Slice x, Slice y, Slice z);
};

class Boundary: public GridItem
{
public:
    Boundary(const std::string& name);
    virtual void update_phi_E();
    virtual void update_phi_H();
    virtual void update_E();
    virtual void update_H();
    virtual ~Boundary() = default;
};

class PML: public Boundary
{
    double k, a;
    int thickness;
    std::string type;
    Tuple<int, int, int> shape = Tuple<int, int, int>(0,0,0);
    Array4<double> phi_E, phi_H, psi_Ex, psi_Ey, psi_Ez, Hx, Hy, Hz, psi_Hx, psi_Hy, psi_Hz, bE, cE, bH, cH, sigmaE, sigmaH;
    Tuple<Slice, Slice, Slice, Slice> loc = Tuple<Slice, Slice, Slice, Slice>(Slice::all(), Slice::all(), Slice::all(), Slice::all());
    Tuple<Slice, Slice, Slice, int> locx = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);
    Tuple<Slice, Slice, Slice, int> locy = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);
    Tuple<Slice, Slice, Slice, int> locz = Tuple<Slice, Slice, Slice, int>(Slice::all(), Slice::all(), Slice::all(), 0);

    void _calculate_parameters(int thickness=10);
    void _set_locations();
    void _set_shape();
    void _set_sigmaE();
    void _set_sigmaH();
public:
    PML(const std::string& name="", double a=1e-8);
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    Array<double> get_sigma(const Array<double>& vect);
    void update_E() override;
    void update_H() override;
    void update_phi_E() override;
    void update_phi_H() override;
};

class PeriodicBoundary: public Boundary
{
    std::string type;
public:
    PeriodicBoundary(const std::string& name="");
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    void update_E() override;
    void update_H() override;
    void update_phi_E() override;
    void update_phi_H() override;
};

class Detector : public GridItem
{
public:
    Detector(const std::string& name);
    virtual void detect_E();
    virtual void detect_H();
    virtual ~Detector() = default;
};

class LineDetector : public Detector
{
    List<int> x;
    List<int> y;
    List<int> z;
    List<Array2<double>> E;
    List<Array2<double>> H;
public:
    LineDetector(const std::string& name);
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    void detect_E() override;
    void detect_H() override;
};

class Source : public GridItem
{
public:
    Source(const std::string& name);
    virtual void update_E();
    virtual void update_H();
    virtual ~Source() = default;
};

class LineSource: public Source
{
    List<int> x;
    List<int> y;
    List<int> z;
    double period;
    int period_n;
    double frequency;
    double amplitude;
    double phase_shift;
    bool pulse;
    int cycle;
    double hanning_dt;
    Array<double> profile;
public:
    LineSource(const std::string& name, double period, double amplitude=1.0, double phase_shift=0.0, bool pulse=false, int cycle=5, double hanning_dt=10.0);
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    void update_E() override;
    void update_H() override;
};

class PointSource: public Source
{
    int x, y, z;
    double period;
    int period_n;
    double frequency;
    double amplitude;
    double phase_shift;
    bool pulse;
    int cycle;
    double hanning_dt;
public:
    PointSource(const std::string& name, double period, double amplitude=1.0, double phase_shift=0.0, bool pulse=false, int cycle=5, double hanning_dt=10.0);
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
    void update_E() override;
    void update_H() override;

};
class Obstacle : public GridItem
{   
public:
    Obstacle(const std::string& name);
    virtual void update_E(const Array4<double>& curl_H);
    virtual void update_H(const Array4<double>& curl_E);
    void register_grid(Grid& grid, Slice x, Slice y, Slice z);
    virtual ~Obstacle() = default;
};


class AnisotropicObstacle: public Obstacle
{
    double permittivity;
    Array3<double> inverse_permittivity;
    int Nx, Ny, Nz;
public:
    AnisotropicObstacle(const std::string& name,  double permittivity);
    void update_E(const Array4<double>& curl_H) override;
    void update_H(const Array4<double>& curl_E) override;
    void register_grid(Grid& grid, Slice x, Slice y, Slice z) override;
};
#endif