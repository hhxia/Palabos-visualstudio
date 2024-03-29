/* This file is part of the Palabos library.
 * Copyright (C) 2009 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas
 */

#ifndef ADVECTION_DIFFUSION_UNITS_H
#define ADVECTION_DIFFUSION_UNITS_H

#include "core/globalDefs.h"
#include "core/globalDefs.h"
#include <string>
#include <fstream>

namespace plb {

/// A useful class for the conversion between dimensionless and lattice units.
template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
class RayleighBenardFlowParam {
public:
    /// Constructor
    /** \param Re_  Reynolds number
     * \param Ra_  Raylegh number
     *  \param Pr_  Prandtl number
     *  \param coldTemperature_  minimum temperature
     *  \param hotTemperature_  maximum temperature
     *  \param deltaT_ time discretization number
     *  \param N_  resolution (a lattice of size 1 has N_+1 cells)
     *  \param lx_ x-length in dimensionless units (e.g. 1)
     *  \param ly_ y-length in dimensionless units (e.g. 1)
     *  \param lz_ z-length in dimensionless units (e.g. 1)
     */
    RayleighBenardFlowParam(T Ra_, T Pr_, T uMax_, T coldTemperature_, 
                            T hotTemperature_, T resolution_,
                            T lx_, T ly_, T lz_=T() )
        : uMax(uMax_), Ra(Ra_), Pr(Pr_), coldTemperature(coldTemperature_), hotTemperature(hotTemperature_),
          resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    { }
    /// Reynolds number
    T getRe() const      { return sqrt(getRa()/getPr()); }
    /// Rayleigh number
    T getRa() const      { return Ra; }
    /// Prandlt number
    T getPr() const      { return Pr; }
    /// delta temperature number
    T getColdTemperature() const   { return coldTemperature; }
    /// delta temperature number
    T getHotTemperature() const    { return hotTemperature; }
    /// delta temperature number
    T getDeltaTemperature() const { return (hotTemperature-coldTemperature); }
    /// delta temperature number
    T getAverageTemperature() const      { return (hotTemperature+coldTemperature)/(T)2; }
    /// resolution (a lattice of size 1 has getN()+1 cells)
    T getResolution() const { return resolution; }
    /// x-length in dimensionless units
    T getLx() const      { return lx; }
    /// y-length in dimensionless units
    T getLy() const      { return ly; }
    /// z-length in dimensionless units
    T getLz() const      { return lz; }
    /// lattice spacing in dimensionless units
    T getDeltaX() const  { return (T)1/resolution; }
    /// time step in dimensionless units
    T getDeltaT() const  { return getLatticeU() / (T)resolution; }
    /// conversion from dimensionless to lattice units for space coordinate
    plint nCell(T l) const { return (plint)(l/getDeltaX()+(T)0.5); }
    /// conversion from dimensionless to lattice units for time coordinate
    plint nStep(T t) const { return (plint)(t/getDeltaT()+(T)0.5); }
    /// number of lattice cells in x-direction
    plint getNx() const    { return nCell(lx)+1; }
    /// number of lattice cells in y-direction
    plint getNy() const    { return nCell(ly)+1; }
    /// number of lattice cells in z-direction
    plint getNz() const    { return nCell(lz)+1; }
    /// velocity in lattice units (proportional to Mach number)
    T getLatticeU() const       { return uMax  ; }
    /// viscosity in lattice units
    T getLatticeNu() const      { return getDeltaT()/(getDeltaX()*getDeltaX()*getRe()); }
    /// thermal conductivity in lattice units
    T getLatticeKappa() const   { return getLatticeNu() / getPr() ; }
    /// viscosity in lattice units
    T getLatticeGravity() const { return getDeltaT() * getDeltaT() / getDeltaX(); }
    /// relaxation time
    T getSolventTau() const   { return nsDescriptor<T>::invCs2*getLatticeNu()+(T)0.5; }
    /// relaxation frequency
    T getSolventOmega() const { return (T)1 / getSolventTau(); }
    /// relaxation time
    T getTemperatureTau() const    { return adDescriptor<T>::invCs2*getLatticeKappa()+(T)0.5; }
    /// relaxation frequency
    T getTemperatureOmega() const  { return (T)1 / getTemperatureTau(); }
private:
    T uMax, Ra, Pr, coldTemperature, hotTemperature, resolution, lx, ly, lz;
};

template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
void writeLogFile(RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> const& parameters,
                  std::string const& title)
{
    std::string fullName = global::directories().getLogOutDir() + "olbLog.dat";
    std::ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    ofile << "Raynleigh number:          Ra=" << parameters.getRa() << "\n";
    ofile << "Prandlt number:            Pr=" << parameters.getPr() << "\n";
    ofile << "Kinematic viscosity:       Nu=" << parameters.getLatticeNu() << "\n";
    ofile << "Thermal conductivity:   Kappa=" << parameters.getLatticeKappa() << "\n";
    ofile << "Lattice resolution:         N=" << parameters.getResolution() << "\n";
    ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
    ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
    ofile << "Solvent omega:        omega_S=" << parameters.getSolventOmega() << "\n";
    ofile << "Temperature omega:    omega_T=" << parameters.getTemperatureOmega() << "\n";
    ofile << "Caracteristic vel:       uMax=" << parameters.getLatticeU() << "\n";
}

}  // namespace plb

#endif
