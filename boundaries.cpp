#include "array.h"
#include "array3.h"
#include "array4.h"
#include "base.h"
#include "grid_items.h"
#include <math.h>
#include <string>

Boundary::Boundary(const std::string& name)
{
    this->name = name;
}

void Boundary::update_phi_E()
{

}

void Boundary::update_phi_H()
{
    
}

void Boundary::update_E()
{
    
}

void Boundary::update_H()
{
    
}
// ================================================================================================
PML::PML(const std::string& name, double a) : Boundary(name)
{
    k = 1.0;
    thickness = 0;
    this->a = a;
}

void PML::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    this->x = x;
    this->y = y;
    this->z = z;
    if(
        (x.start.is_null || x.start.value == 0) &&
        !x.stop.is_null && x.stop.value > 0
    )
    {
        type = "PMLXlow";
        _calculate_parameters(x.stop.value);
    }
    else if(!x.start.is_null && x.stop.is_null && x.start.value < 0)
    {
        type = "PMLXhigh";
        _calculate_parameters(-x.start.value);
    }
    else if(
        (y.start.is_null || y.start.value == 0) &&
        !y.stop.is_null && y.stop.value > 0
    )
    {
        type = "PMLYlow";
        _calculate_parameters(y.stop.value);
    }
    else if(!y.start.is_null && y.stop.is_null && y.start.value < 0)
    {
        type = "PMLYhigh";
        _calculate_parameters(-y.start.value);
    }
    else if(
        (z.start.is_null || z.start.value == 0) &&
        !z.stop.is_null && z.stop.value > 0
    )
    {
        type = "PMLZlow";
        _calculate_parameters(z.stop.value);
    }
    else if(!z.start.is_null && z.stop.is_null && z.start.value < 0)
    {
        type = "PMLZhigh";
        _calculate_parameters(-z.start.value);
    }
}

Array<double> PML::get_sigma(const Array<double>& vect)
{
    return vect.pow(3) * 40.0 / pow(thickness + 1, 4);
}

void PML::_set_locations()
{
    auto all = Slice::all();
    auto head = Slice::head(thickness);
    auto tail = Slice::tail(-thickness);

    if(type == "PMLXlow")
    {
        loc = Tuple<Slice, Slice, Slice, Slice>(head, all, all, all);
        locx = Tuple<Slice, Slice, Slice, int>(head, all, all, 0);
        locy = Tuple<Slice, Slice, Slice, int>(head, all, all, 1);
        locz = Tuple<Slice, Slice, Slice, int>(head, all, all, 2);
    }
    else if(type == "PMLXhigh")
    {   
        loc = Tuple<Slice, Slice, Slice, Slice>(tail, all, all, all);
        locx = Tuple<Slice, Slice, Slice, int>(tail, all, all, 0);
        locy = Tuple<Slice, Slice, Slice, int>(tail, all, all, 1);
        locz = Tuple<Slice, Slice, Slice, int>(tail, all, all, 2);
    }
    else if(type == "PMLYlow")
    {
        loc = Tuple<Slice, Slice, Slice, Slice>(all, head, all, all);
        locx = Tuple<Slice, Slice, Slice, int>(all, head, all, 0);
        locy = Tuple<Slice, Slice, Slice, int>(all, head, all, 1);
        locz = Tuple<Slice, Slice, Slice, int>(all, head, all, 2);
    }
    else if(type == "PMLYhigh")
    {   
        loc = Tuple<Slice, Slice, Slice, Slice>(all, tail, all, all);
        locx = Tuple<Slice, Slice, Slice, int>(all, tail, all, 0);
        locy = Tuple<Slice, Slice, Slice, int>(all, tail, all, 1);
        locz = Tuple<Slice, Slice, Slice, int>(all, tail, all, 2);
    }
    else if(type == "PMLZlow")
    {
        loc = Tuple<Slice, Slice, Slice, Slice>(all, all, head, all);
        locx = Tuple<Slice, Slice, Slice, int>(all, all, head, 0);
        locy = Tuple<Slice, Slice, Slice, int>(all,all, head, 1);
        locz = Tuple<Slice, Slice, Slice, int>(all,all, head,  2);
    }
    else if(type == "PMLZhigh")
    {   
        loc = Tuple<Slice, Slice, Slice, Slice>(all,all, tail, all);
        locx = Tuple<Slice, Slice, Slice, int>(all,all, tail, 0);
        locy = Tuple<Slice, Slice, Slice, int>(all,all, tail, 1);
        locz = Tuple<Slice, Slice, Slice, int>(all,all, tail, 2);
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }
}

void PML::_set_shape()
{
    if(type == "PMLXlow" || type == "PMLXhigh")
    {
        shape = Tuple<int, int, int>(thickness, grid->Ny, grid->Nz);
    }
    else if(type == "PMLYlow" || type == "PMLYhigh")
    {
        shape = Tuple<int, int, int>(grid->Nx, thickness, grid->Nz);
    }
    else if(type == "PMLZlow" || type == "PMLZhigh")
    {
        shape = Tuple<int, int, int>(grid->Nx, grid->Ny, thickness);
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }
}

Array<double> arange(double start, double stop, double step)
{
    List<double> lst;
    if(step > 0)
    {
        while(start < stop)
        {
            lst.append(start);
            start += step;
        }
    }
    else if(step < 0) 
    {
        while(start > stop)
        {
            lst.append(start);
            start += step;
        }
    }
    else 
    {
        panic("Step 0 is now allowed");
    }
    return Array<double>::from_list(lst);
}

void PML::_set_sigmaE()
{
    auto all = Slice::all();
    
    auto vect = (type == "PMLXlow"  || type == "PMLYlow" || type == "PMLZlow")?
            arange(thickness - 0.5, -0.5, -1.0) :
            arange(0.5, thickness + 0.5, 1.0);

    auto sigma = get_sigma(vect);

    if(type == "PMLXlow" || type == "PMLXhigh")
    {   
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(thickness, grid->Ny, grid->Nz));
        for(int i1 = 0; i1 < grid->Ny; i1++)
        {
            for(int i2 = 0; i2 < grid->Nz; i2++)
            {
                sigma3.get(all, i1, i2).copy_from(sigma);
            }
        } 
        sigmaE = Array4<double>::zeros(Tuple<int, int, int, int>(thickness, grid->Ny, grid->Nz, 3));
        sigmaE.get(all, all, all, 0).copy_from(sigma3);
    }
    else if(type == "PMLYlow" || type == "PMLYhigh")
    {
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(grid->Nx, thickness, grid->Nz));

        for(int i1 = 0; i1 < grid->Nx; i1++)
        {
            for(int i2 = 0; i2 < grid->Nz; i2++)
            {
                sigma3.get(i1, all, i2).copy_from(sigma);
            }
        }

        sigmaE = Array4<double>::zeros(Tuple<int, int, int, int>(grid->Nx, thickness, grid->Nz, 3));
        sigmaE.get(all, all, all, 1).copy_from(sigma3);
        
    }
    else if(type == "PMLZlow" || type == "PMLZhigh")
    {
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(grid->Nx, grid->Ny, thickness));

        for(int i1 = 0; i1 < grid->Nx; i1++)
        {
            for(int i2 = 0; i2 < grid->Ny; i2++)
            {
                sigma3.get(i1, i2, all).copy_from(sigma);
            }
        }

        sigmaE = Array4<double>::zeros(Tuple<int, int, int, int>(grid->Nx, grid->Ny, thickness, 3));
        sigmaE.get(all, all, all, 2).copy_from(sigma3);
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }
}

void PML::_set_sigmaH()
{
    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    
    auto vect = (type == "PMLXlow"  || type == "PMLYlow" || type == "PMLZlow")?
            arange(thickness - 1.0, 0.0, -1.0) :
            arange(1.0, thickness, 1.0);

    auto sigma = get_sigma(vect);

    if(type == "PMLXlow" || type == "PMLXhigh")
    {   
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(thickness-1, grid->Ny, grid->Nz));
        for(int i1 = 0; i1 < grid->Ny; i1++)
        {
            for(int i2 = 0; i2 < grid->Nz; i2++)
            {
                sigma3.get(all, i1, i2).copy_from(sigma);
            }
        } 
        sigmaH = Array4<double>::zeros(Tuple<int, int, int, int>(thickness, grid->Ny, grid->Nz, 3));
        sigmaH.get(head_1, all, all, 0).copy_from(sigma3);
    }
    else if(type == "PMLYlow" || type == "PMLYhigh")
    {
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(grid->Nx, thickness-1, grid->Nz));

        for(int i1 = 0; i1 < grid->Nx; i1++)
        {
            for(int i2 = 0; i2 < grid->Nz; i2++)
            {
                sigma3.get(i1, all, i2).copy_from(sigma);
            }
        }

        sigmaH = Array4<double>::zeros(Tuple<int, int, int, int>(grid->Nx, thickness, grid->Nz, 3));
        sigmaH.get(all, head_1, all, 1).copy_from(sigma3);
        
    }
    else if(type == "PMLZlow" || type == "PMLZhigh")
    {
        auto sigma3 = Array3<double>::zeros(Tuple<int, int, int>(grid->Nx, grid->Ny, thickness-1));

        for(int i1 = 0; i1 < grid->Nx; i1++)
        {
            for(int i2 = 0; i2 < grid->Ny; i2++)
            {
                sigma3.get(i1, i2, all).copy_from(sigma);
            }
        }

        sigmaH = Array4<double>::zeros(Tuple<int, int, int, int>(grid->Nx, grid->Ny, thickness, 3));
        sigmaH.get(all, all, head_1, 2).copy_from(sigma3);
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }

    
}
void PML::_calculate_parameters(int thickness)
{
    this->thickness = thickness;
    _set_locations();
    _set_shape();
    _set_sigmaE();
    _set_sigmaH();
    int Nx = shape.item1, Ny = shape.item2, Nz = shape.item3;
    phi_E = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    phi_H = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Ex = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Ey = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Ez = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Hx = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Hy = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));
    psi_Hz = Array4<double>::zeros(Tuple<int, int, int, int>(Nx, Ny, Nz, 3));

    bE = (-(sigmaE/k + a) * grid->get_courant_number()).exp();
    cE = (bE - 1.0) * sigmaE / (sigmaE * k + a * k * k);

    bH = (-(sigmaH/k + a) * grid->get_courant_number()).exp();
    cH = (bH - 1.0) * sigmaH / (sigmaH * k + a * k * k);
}

void PML::update_E()
{
    grid->E.get(loc) += grid->inverse_permittivity.get(loc) * grid->get_courant_number() * phi_E;
}

void PML::update_H()
{
    grid->H.get(loc) += grid->inverse_permeability.get(loc) * grid->get_courant_number() * phi_H;
}

void PML::update_phi_E()
{
    psi_Ex *= bE;
    psi_Ey *= bE;
    psi_Ez *= bE;
    
    auto Hx = grid->H.get(locx);
    auto Hy = grid->H.get(locy);
    auto Hz = grid->H.get(locz);

    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);

    //psi_Ex[:, tail_1, :, 1] += (Hz[:, tail_1, :] - Hz[:, head_1, :]) * c[:, tail_1, :, 1]
    psi_Ex.get(all, tail_1, all, 1) += (Hz.get(all, head_1, all) - Hz.get(all, tail_1, all)) * cE.get(all, tail_1, all, 1);
    
    //psi_Ex[:, :, tail_1, 2] += (Hy[:, :, tail_1] - Hy[:, :, head_1]) * c[:, :, tail_1, 2]
    psi_Ex.get(all, all, tail_1, 2) += (Hy.get(all, all, tail_1) - Hy.get(all, all, head_1)) * cE.get(all, all, tail_1, 2);

    //psi_Ey[:, :, tail_1, 2] += (Hx[:, :, tail_1] - Hx[:, :, head_1]) * c[:, :, tail_1, 2]
    psi_Ey.get(all, all, tail_1, 2) += (Hx.get(all, all, tail_1) - Hx.get(all, all, head_1)) * cE.get(all, all, tail_1, 2);

    //psi_Ey[tail_1, :, :, 0] += (Hz[tail_1, :, :] - Hz[head_1, :, :]) * c[tail_1, :, :, 0]
    psi_Ey.get(tail_1, all, all, 0) += (Hz.get(tail_1, all, all) - Hz.get(head_1, all, all)) * cE.get(tail_1, all, all, 0);

    //psi_Ez[tail_1, :, :, 0] += (Hy[tail_1, :, :] - Hy[head_1, :, :]) * c[tail_1, :, :, 0]
    psi_Ex.get(tail_1, all, all, 0) += (Hy.get(tail_1, all, all) - Hy.get(head_1, all, all)) * cE.get(tail_1, all, all, 0);

    //psi_Ez[:, tail_1, :, 1] += (Hx[:, tail_1, :] - Hx[:, head_1, :]) * c[:, tail_1, :, 1]
    psi_Ez.get(all, tail_1, all, 1) += (Hx.get(all, tail_1, all) - Hx.get(all, head_1, all)) * cE.get(all, tail_1, all, 1);

    //phi_E[..., 0] = psi_Ex[..., 1] - psi_Ex[..., 2]
    phi_E.get(all, all, all, 0) = psi_Ex.get(all, all, all, 1) - psi_Ex.get(all, all, all, 2);

    //phi_E[..., 1] = psi_Ey[..., 2] - psi_Ey[..., 0]
    phi_E.get(all, all, all, 1) = psi_Ey.get(all, all, all, 2) - psi_Ey.get(all, all, all, 0);

    //phi_E[..., 2] = psi_Ez[..., 0] - psi_Ez[..., 1]
    phi_E.get(all, all, all, 2) = psi_Ez.get(all, all, all, 0) - psi_Ez.get(all, all, all, 1);
}

void PML::update_phi_H()
{
    psi_Hx *= bH;
    psi_Hy *= bH;
    psi_Hz *= bH;

    auto all = Slice::all();
    auto head_1 = Slice::head(-1);
    auto tail_1 = Slice::tail(1);


    auto Ex = grid->E.get(locx);
    auto Ey = grid->E.get(locy);
    auto Ez = grid->E.get(locz);

    //psi_Hx[:, head_1, :, 1] += (Ez[:, tail_1, :] - Ez[:, head_1, :]) * c[:, head_1, :, 1]
    psi_Hx.get(all, head_1, all, 1) += (Ez.get(all, tail_1, all) - Ez.get(all, head_1, all)) * cH.get(all, head_1, all, 1);

    //psi_Hx[:, :, head_1, 2] += (Ey[:, :, tail_1] - Ey[:, :, head_1]) * c[:, :, head_1, 2]
    psi_Hx.get(all, all, head_1, 2) +=  (Ey.get(all, all, tail_1) - Ey.get(all, all, head_1)) * cH.get(all, all, head_1, 2);    

    //psi_Hy[:, :, head_1, 2] += (Ex[:, :, tail_1] - Ex[:, :, head_1]) * c[:, :, head_1, 2]
    psi_Hy.get(all, all, head_1, 2) += (Ex.get(all, all, tail_1) - Ex.get(all, all, head_1)) * cH.get(all, all, head_1, 2);
    
    //psi_Hy[head_1, :, :, 0] += (Ez[tail_1, :, :] - Ez[head_1, :, :]) * c[head_1, :, :, 0]
    psi_Hy.get(head_1, all, all, 0) += (Ez.get(tail_1, all, all) - Ex.get(head_1, all, all)) * cH.get(head_1, all, all, 0);

    //psi_Hz[head_1, :, :, 0] += (Ey[tail_1, :, :] - Ey[head_1, :, :]) * c[head_1, :, :, 0]
    psi_Hz.get(head_1, all, all, 0) += (Ey.get(tail_1, all, all) - Ey.get(head_1, all, all)) * cH.get(head_1, all, all, 0);

    //psi_Hz[:, head_1, :, 1] += (Ex[:, tail_1, :] - Ex[:, head_1, :]) * c[:, head_1, :, 1]
    psi_Hz.get(all, head_1, all, 1) += (Ex.get(all, tail_1, all) - Ex.get(all, head_1, all)) * cH.get(all, head_1, all, 1);

    //phi_H[..., 0] = psi_Hx[..., 1] - psi_Hx[..., 2]
    phi_H.get(all, all, all, 0) = psi_Hx.get(all, all, all, 1) - psi_Hx.get(all, all, all, 2);
    
    //phi_H[..., 1] = psi_Hy[..., 2] - psi_Hy[..., 0]
    phi_H.get(all, all, all, 1) = psi_Hy.get(all, all, all, 2) - psi_Hy.get(all, all, all, 0);

    //phi_H[..., 2] = psi_Hz[..., 0] - psi_Hz[..., 1]
    phi_H.get(all, all, all, 2) = psi_Hz.get(all, all, all, 0) - psi_Hz.get(all, all, all, 1);
}

PeriodicBoundary::PeriodicBoundary(const std::string& name) : Boundary(name)
{

}

void PeriodicBoundary::register_grid(Grid& grid, Slice x, Slice y, Slice z)
{
    this->grid = &grid;
    this->x = x;
    this->y = y;
    this->z = z;
    if(!x.start.is_null)
    {
        type = "PeriodicBoundaryX";
    }
    else if(!y.start.is_null)
    {
        type = "PeriodicBoundaryY";
    }
    else if(!z.start.is_null)
    {
        type = "PeriodicBoundaryZ";
    }
    else 
    {
        panic("A periodic boundary should be placed at the boundary of the grid using a single index");
    }
}

void PeriodicBoundary::update_E()
{
    auto all = Slice::all();
    if(type == "PeriodicBoundaryX")
    {
        grid->E[0].copy_from(grid->E[-1]);
    }
    else if (type == "PeriodicBoundaryY")
    {
        grid->E.get(all, 0).copy_from(grid->E.get(all, -1));
    }
    else if (type == "PeriodicBoundaryZ")
    {
        grid->E.get(all, all, 0).copy_from(grid->E.get(all, all, -1));
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }
}

void PeriodicBoundary::update_H()
{
    auto all = Slice::all();
    if(type == "PeriodicBoundaryX")
    {
        grid->H[-1].copy_from(grid->H[0]);
    }
    else if (type == "PeriodicBoundaryY")
    {
        grid->H.get(all, -1).copy_from(grid->E.get(all, 0));
    }
    else if (type == "PeriodicBoundaryZ")
    {
        grid->H.get(all, all, -1).copy_from(grid->H.get(all, all, 0));
    }
    else 
    {
        panic("Unknown type :% s", type.c_str());
    }
}

void PeriodicBoundary::update_phi_E()
{

}

void PeriodicBoundary::update_phi_H()
{
    
}