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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_HH
#define ADVECTION_DIFFUSION_DYNAMICS_HH

#include "core/latticeStatistics.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "latticeBoltzmann/advectionDiffusionMomentTemplates.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"


namespace plb {
    
/* *************** Class AdvectionDiffusionDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionDynamics<T,Descriptor>::AdvectionDiffusionDynamics(T omega_)
  : BasicBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionDynamics<T,Descriptor>::regularize (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{
    // jAdvDiff is the first order moment of 
    typedef Descriptor<T> D;
    Array<T,Descriptor<T>::d> jEq;
    
    advectionDiffusionMomentTemplates<T,Descriptor>::get_jEq(cell, rhoBar, jEq);
    
    advectionDiffusionDynamicsTemplates<T,Descriptor>::regularize(cell,rhoBar,j,jEq);
}

template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}


/* *************** Class AdvectionDiffusionRLBdynamics *************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T,Descriptor>::AdvectionDiffusionRLBdynamics (T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionRLBdynamics<T,Descriptor>* AdvectionDiffusionRLBdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionRLBdynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionRLBdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics<T>& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq, jNeq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq_jNeq(cell, rhoBar, jEq, jNeq);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_rlb_collision(cell, rhoBar, jEq, jNeq, this->getOmega() );
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionRLBdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);
}


/* *************** Class AdvectionDiffusionBGKdynamics *************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T,Descriptor>::AdvectionDiffusionBGKdynamics (
        T omega_ )
    : AdvectionDiffusionDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
AdvectionDiffusionBGKdynamics<T,Descriptor>* AdvectionDiffusionBGKdynamics<T,Descriptor>::clone() const {
    return new AdvectionDiffusionBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AdvectionDiffusionBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics<T>& statistics ) 
{
    T rhoBar;
    Array<T,Descriptor<T>::d> jEq;
    advectionDiffusionMomentTemplates<T,Descriptor>::get_rhoBar_jEq(cell, rhoBar, jEq);
    
    T uSqr = advectionDiffusionDynamicsTemplates<T,Descriptor>::
            no_corr_bgk_collision(cell, rhoBar, jEq, this->getOmega());
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

/** \param j The parameter j is defined as j = j_advDiff = rho_advDiff*u_fluid
 */
template<typename T, template<typename U> class Descriptor>
T AdvectionDiffusionBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{ 
    return advectionDiffusionDynamicsTemplates<T,Descriptor>::
            bgk_ma1_equilibrium(iPop, rhoBar, j);

}


} // namespace plb

#endif
