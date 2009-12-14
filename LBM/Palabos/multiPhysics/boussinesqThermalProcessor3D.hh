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

#ifndef BOUSSINESQ_THERMAL_PROCESSOR_3D_HH
#define BOUSSINESQ_THERMAL_PROCESSOR_3D_HH

#include "multiPhysics/boussinesqThermalProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::
        BoussinesqThermalProcessor3D(T gravity_, T T0_, T deltaTemp_, Array<T,FluidDescriptor<T>::d> dir_)
    :  gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
       dir(dir_)
{
    // We normalize the direction of the force vector.
    T normDir = sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
void BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,FluidDescriptor>& fluid,
        BlockLattice3D<T,TemperatureDescriptor>& temperature )
{
    typedef FluidDescriptor<T> D;
    enum {
        velOffset   = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt,
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };
    
    Array<T,D::d> gravOverDeltaTemp (
            gravity*dir[0]/deltaTemp,
            gravity*dir[1]/deltaTemp,
            gravity*dir[2]/deltaTemp );
           
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) 
            {
                // Velocity coupling
                T *u = temperature.get(iX,iY,iZ).getExternal(velOffset);
                Array<T,FluidDescriptor<T>::d> vel;
                fluid.get(iX,iY,iZ).computeVelocity(vel);
                vel.to_cArray(u);
                
                // Computation of the Boussinesq force
                T *force = fluid.get(iX,iY,iZ).getExternal(forceOffset);
                // Temperature is the order-0 moment of the advection-diffusion lattice.
                //   You can compute it with the method computeDensity().
                T localTemperature = temperature.get(iX,iY,iZ).computeDensity();
                T rho = fluid.get(iX,iY,iZ).computeDensity();
                const T diffT = rho*(localTemperature - T0);
                for (pluint iD = 0; iD < D::d; ++iD) 
                {
                    force[iD] = gravOverDeltaTemp[iD] * diffT;
                }
            }
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>*
    BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>::clone() const
{
    return new BoussinesqThermalProcessor3D<T,FluidDescriptor,TemperatureDescriptor>(*this);
}

}  // namespace plb

#endif
