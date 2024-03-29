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

/* Main author: Orestis Malaspinas.
 */

/** \file
 * Neumann boundary conditions for temperature -- generic implementation.
 */
#ifndef ADIABATIC_BOUNDARY_PROCESSOR_2D_HH
#define ADIABATIC_BOUNDARY_PROCESSOR_2D_HH

#include "complexDynamics/adiabaticBoundaryProcessor2D.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation> 
void FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint iX_prev = iX + ( (direction==0) ? (-orientation) : 0 );
            plint iY_prev = iY + ( (direction==1) ? (-orientation) : 0 );
            
            plint iX_prev_2 = iX + ( (direction==0) ? (-2*orientation) : 0 );
            plint iY_prev_2 = iY + ( (direction==1) ? (-2*orientation) : 0 );
            
            T temperature_1 = lattice.get(iX_prev,iY_prev).computeDensity();
            T temperature_2 = lattice.get(iX_prev_2,iY_prev_2).computeDensity();
            
            T temperature = fd::boundaryZeroGradient(temperature_1, temperature_2);
            lattice.get(iX,iY).defineDensity(temperature);
        }
    }
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation> 
FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>*
    FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>::clone() const
{
    return new FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation> 
void FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation> 
BlockDomain::DomainT FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>::appliesTo() const 
{
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // NEUMANN_CONDITION_2D_HH
