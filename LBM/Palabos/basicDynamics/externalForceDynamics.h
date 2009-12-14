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
 * can be instantiated -- header file.
 */
#ifndef EXTERNAL_FORCE_DYNAMICS_H
#define EXTERNAL_FORCE_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/** Implementation of the computation of the velocity in the 
  * presence of an external force
*/
template<typename T, template<typename U> class Descriptor>
class ExternalForceDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ExternalForceDynamics(T omega_);

/* *************** Velocity computation ************************* */

    /// Implementation of velocity computation
    virtual void computeVelocity( Cell<T,Descriptor> const& cell, 
                                  Array<T,Descriptor<T>::d>& u ) const;
};


/// Implementation of O(Ma^2) BGK dynamics with an external force (Guo style)
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceBGKdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

}

#endif
