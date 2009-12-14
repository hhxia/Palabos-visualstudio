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

/** \file
 * Dirichlet boundary condition which imposes equilibrium (but computes
 * density properly from velocity, or vice versa)
 */
#ifndef EQUILIBRIUM_BOUNDARY_DYNAMICS_HH
#define EQUILIBRIUM_BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "core/cell.h"

namespace plb {


/* *************** Class EquilibriumVelocityBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
EquilibriumVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::
    EquilibriumVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
EquilibriumVelocityBoundaryDynamics<T,Descriptor,direction,orientation>*
    EquilibriumVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new EquilibriumVelocityBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void EquilibriumVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::
    completePopulations(Cell<T,Descriptor>& cell) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    this -> computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}


/* *************** Class EquilibriumDensityBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
EquilibriumDensityBoundaryDynamics<T,Descriptor,direction,orientation>::
    EquilibriumDensityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
EquilibriumDensityBoundaryDynamics<T,Descriptor,direction,orientation>*
    EquilibriumDensityBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new EquilibriumDensityBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void EquilibriumDensityBoundaryDynamics<T,Descriptor,direction,orientation>::
    completePopulations(Cell<T,Descriptor>& cell) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    this -> computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->getBaseDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

}  // namespace plb

#endif  // EQUILIBRIUM_BOUNDARY_DYNAMICS_HH
