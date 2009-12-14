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
#ifndef EXTERNAL_FORCE_DYNAMICS_HH
#define EXTERNAL_FORCE_DYNAMICS_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class ExternalForceDynamics *********************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
ExternalForceDynamics<T,Descriptor>::ExternalForceDynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
void ExternalForceDynamics<T,Descriptor>::computeVelocity( Cell<T,Descriptor> const& cell, 
                                  Array<T,Descriptor<T>::d>& u ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> force, j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar, j);
    force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        u[iD] = (j[iD] + force[iD] /(T)2) * invRho;
    }
    
}


/* *************** Class GuoExternalForceBGKdynamics ********************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
GuoExternalForceBGKdynamics<T,Descriptor>::GuoExternalForceBGKdynamics(T omega_ )
    : ExternalForceDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
GuoExternalForceBGKdynamics<T,Descriptor>* GuoExternalForceBGKdynamics<T,Descriptor>::clone() const {
    return new GuoExternalForceBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void GuoExternalForceBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    T rhoBar = this->computeRhoBar(cell);
    Array<T,Descriptor<T>::d> u,j;
    this->computeVelocity(cell,u);
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (plint iD = 0; iD < Descriptor<T>::d; ++iD)
    {
        j[iD] = rho * u[iD];
    }
    
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    externalForceTemplates<T,Descriptor>::addGuoForce(cell, u, this->getOmega(), (T)1);
    
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T GuoExternalForceBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


}

#endif
