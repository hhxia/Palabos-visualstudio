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

#ifndef CARREAU_UNITS_H
#define CARREAU_UNITS_H

#include "core/globalDefs.h"
#include "io/parallelIO.h"
#include <string>
#include <fstream>

/* Main author: Orestis Malaspinas
 */

namespace plb {

/// Numeric parameters for isothermal, incompressible flow.
template<typename T>
class CarreauFlowParam {
public:
    /// Constructor
    /** \param latticeU_  Characteristic velocity in lattice units (proportional to Mach number).
     *  \param Re_ Reynolds number.
     *  \param N_  Resolution (a lattice of size 1 has N_+1 cells).
     *  \param lx_ x-length in dimensionless units (e.g. 1).
     *  \param ly_ y-length in dimensionless units (e.g. 1).
     *  \param lz_ z-length in dimensionless units (e.g. 1).
     */
    CarreauFlowParam(T latticeU_, T Re_, T Cu_, T nuInf_, T n_,
                          plint resolution_, T lx_, T ly_, T lz_=T() )
        : latticeU(latticeU_), Re(Re_), Cu(Cu_), nuInf(nuInf_), n(n_),
          resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    { }
    /// velocity in lattice units (proportional to Mach number)
    T getLatticeU() const { return latticeU; }
    /// Reynolds number
    T getRe() const      { return Re; }
    /// Carreau number (lambda*u/L)
    T getCu() const      { return Cu; }
    /// The exponent of the power-law
    T getExponent() const      { return n; }
    /// resolution
    plint getResolution() const { return resolution; }
    /// x-length in dimensionless units
    T getLx() const      { return lx; }
    /// y-length in dimensionless units
    T getLy() const      { return ly; }
    /// z-length in dimensionless units
    T getLz() const      { return lz; }
    /// lattice spacing in dimensionless units
    T getDeltaX() const  { return (T)1/(T)getResolution(); }
    /// time step in dimensionless units
    T getDeltaT() const  { return getDeltaX()*getLatticeU(); }
    /// conversion from dimensionless to lattice units for space coordinate
    plint nCell(T l) const { return (int)(l/getDeltaX()+(T)0.5); }
    /// conversion from dimensionless to lattice units for time coordinate
    plint nStep(T t) const { return (int)(t/getDeltaT()+(T)0.5); }
    /// number of lattice cells in x-direction
    plint getNx(bool offLattice=false) const { return nCell(lx)+1+(int)offLattice; }
    /// number of lattice cells in y-direction
    plint getNy(bool offLattice=false) const { return nCell(ly)+1+(int)offLattice; }
    /// number of lattice cells in z-direction
    plint getNz(bool offLattice=false) const { return nCell(lz)+1+(int)offLattice; }
    /// solvent viscosity at zero shear rate in lattice units
    T getLatticeNu0() const { return getLatticeU()*getResolution()/getRe(); }
	/// solvent viscosity at infinite shear rate in lattice units
	T getLatticeNuInf() const { return nuInf * getDeltaT() / (getDeltaX() * getDeltaX()) ; }
    /// lambda paramter in lattice units
    T getLatticeLambda() const { return getResolution()/getLatticeU()*getCu(); }
    /// solvent relaxation time
    T getTau0() const       { return (T)3*getLatticeNu0()+(T)0.5; }
    /// solvent relaxation frequency
    T getOmega0() const     { return (T)1 / getTau0(); }
private:
    T latticeU, Re, Cu, nuInf, n;
    plint resolution;
    T lx, ly, lz;
};

template<typename T>
void writeLogFile(CarreauFlowParam<T> const& parameters,
                  std::string const& title)
{
    std::string fullName = global::directories().getLogOutDir() + "olbLog.dat";
    plb_ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
    ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    ofile << "Carreau number:            Cu=" << parameters.getCu() << "\n";
    ofile << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
    ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
    ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
    ofile << "Exponent:                  n=" << parameters.getExponent() << "\n";
    ofile << "Zero viscosity:            nu0=" << parameters.getLatticeNu0() << "\n";
	ofile << "Inf viscosity:             nuInf=" << parameters.getLatticeNuInf() << "\n";
    ofile << "Zero Omega:                omega0=" << parameters.getOmega0() << "\n";
    ofile << "Lambda:                    lambda=" << parameters.getLatticeLambda() << "\n";
}

}  // namespace plb

#endif
