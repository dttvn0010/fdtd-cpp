#include "array.h"
#include "grid_items.h"
#include <math.h>
#include <string>

// Source constructor
// Input:
//  + name: name of the detector
Source::Source(const std::string& name)
{
    this->name = name;
}

// Update electric field at the source
// This is a virtual function, the real implemenation in inside each children class
void Source::update_E()
{

}

// Update magnetic field at the source
// This is a virtual function, the real implemenation in inside each children class
void Source::update_H()
{

}

// ======================================================================================================================================
// LineSource constructor
// Input:
//  + name: name of the detector
//  + period: The period of the source (in seconds)
//  + amplitude: The electric field amplitude (in simulation units)
//  + phase_shift: The phase offset of the source.
//  + pulse: Hanning pulse or continuous waveform
//  + cycle: cycles for Hanning pulse.
//  + hanning_dt: timestep used for Hanning pulse
LineSource::LineSource(const std::string& name, double period, double amplitude, double phase_shift, bool pulse, int cycle, double hanning_dt) : 
    Source(name)
{
    this->period = period;
    this->amplitude = amplitude;
    this->phase_shift = phase_shift;
    this->pulse = pulse;
    this->cycle = cycle;
    this->frequency = 1.0 / period;
    this->hanning_dt = hanning_dt;
}

// LineSource - register grid: add line source to the  grid
// Input:
//   + grid: the system grid
//   + x: x-range of the
//   + y: y-range of the
//   + z: z-range of the
void LineSource::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    period_n = grid.handle_time(period);

    // resolving the start/end of location in x axis
    int x0 = x.start.is_null ? 0 : x.start.value;
    int x1 = x.stop.is_null ? 0 : x.stop.value;
    
    // resolving the start/end of location in y axis
    int y0 = y.start.is_null ? 0 : y.start.value;
    int y1 = y.stop.is_null ? 0 : y.stop.value;
    
    // resolving the start/end of location in z axis
    int z0 = z.start.is_null ? 0 : z.start.value;
    int z1 = z.stop.is_null ? 0 : z.stop.value;
    
    // m = max(x1-x0, y1-y0, z1-z0)
    int m = x1-x0;
    if(m < y1-y0) m = y1-y0;
    if(m < z1-z0) m = z1-z0;

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

    int L = this->x.size();
    auto vect = (Array<int>::from_list(this->x) - this->x[L/2]).as_type<double>().pow(2) +
        (Array<int>::from_list(this->y) - this->y[L/2]).as_type<double>().pow(2) +
        (Array<int>::from_list(this->z) - this->z[L/2]).as_type<double>().pow(2);

    profile = (-vect * vect / (2 * pow(0.5 * vect.reduce_max(), 2))).exp();
    profile /= profile.reduce_sum<double>();
    profile *= amplitude;
}


// Hanning pulse function
double hanning(double f, double t, double n)
{
    return 0.5 * (1 - cos(f*t/n)) * sin(f*t);
}

// Update electric field at the line source
void LineSource::update_E()
{
    int q = grid->get_time_steps_passed();
    Array<double> vect = Array<double>::zeros(profile.size());
    if(pulse)
    {
        int t1 = int(2 * M_PI / (frequency * hanning_dt / cycle));
        if(q < t1)
        {
            vect = profile * hanning(frequency, q * hanning_dt, cycle);
        }
    }
    else 
    {
        vect = profile * sin(2 * M_PI * q / period_n + phase_shift);
    }
    for(int i = 0; i < x.size(); i++)
    {
        grid->E.at(x[i], y[i], z[i], 2) += vect[i];
    }
}

// Update magnetic field at the line source
void LineSource::update_H()
{

}

// ======================================================================================================================================
// PointSource constructor
// Input:
//  + name: name of the detector
//  + period: The period of the source (in seconds)
//  + amplitude: The electric field amplitude (in simulation units)
//  + phase_shift: The phase offset of the source.
//  + pulse: Hanning pulse or continuous waveform
//  + cycle: cycles for Hanning pulse.
//  + hanning_dt: timestep used for Hanning pulse
PointSource::PointSource(const std::string& name, double period, double amplitude, double phase_shift, bool pulse, int cycle, double hanning_dt) : 
    Source(name)
{
    this->period = period;
    this->amplitude = amplitude;
    this->phase_shift = phase_shift;
    this->pulse = pulse;
    this->cycle = cycle;
    this->frequency = 1.0 / period;
    this->hanning_dt = hanning_dt;
}

// PointSource - register grid: add point source to the  grid
// Input:
//   + grid: the system grid
//   + x: x-range of the
//   + y: y-range of the
//   + z: z-range of the
void PointSource::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    this->x = x.start.value;
    this->y = y.start.value;
    this->z = z.start.value;
    period_n = grid.handle_time(period);
}

// Update electric field at the line source
void PointSource::update_E()
{
    int q = grid->get_time_steps_passed();
    double src = 0;
    if(pulse)
    {
        int t1 = int(2 * M_PI / (frequency * hanning_dt / cycle));
        if(q < t1)
        {
            src = amplitude * hanning(frequency, q * hanning_dt, cycle);
        }
    }
    else
    {
        src = amplitude * sin(2 * M_PI * q / period_n + phase_shift);
    }
    grid->E.at(x, y, z, 2) += src;
}

// Update magnetic field at the line source
void PointSource::update_H()
{
    
}