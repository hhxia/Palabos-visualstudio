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

#ifndef BOUSSINESQ_THERMAL_PROCESSOR_3D_H
#define BOUSSINESQ_THERMAL_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

/**
* Multiphysics class for coupling between Navier-Stokes and advection-diffusion
* equations using the Boussinesq approximation.
*/
template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
class BoussinesqThermalProcessor3D :
    public BoxProcessingFunctional3D_LL<T,FluidDescriptor,TemperatureDescriptor>
{
public:
    
    BoussinesqThermalProcessor3D(T gravity_, T T0_, T deltaTemp_,
                                 Array<T,FluidDescriptor<T>::d> dir_);
    
    virtual void process( Box3D domain,
                          BlockLattice3D<T,FluidDescriptor>& fluid,
                          BlockLattice3D<T,TemperatureDescriptor>& temperature );
    virtual BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>* clone() const;
    
private:
    T gravity, T0, deltaTemp;
    Array<T,FluidDescriptor<T>::d> dir;
};

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_3D_H
