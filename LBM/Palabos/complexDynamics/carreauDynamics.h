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

#ifndef CARREAU_DYNAMICS_H
#define CARREAU_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"

namespace plb {

	
/// This class recomputes omega for a generalied newtonian fluid with a carreau constitutive equation.
/** The constitutive equation is nu=nu0*(1+(lambda*|gamma|)^2)^((n-1)/2).
 * Note that in order to be numerically efficient there is no resolution of the implicit equation
 * rather we try to get to the fixed point by auto-replacement in the solution.(??? need to be rewritten)
*/
template<typename T, template<typename U> class Descriptor, int N>
class CarreauDynamics : public OmegaFromPiDynamics<T,Descriptor> {
public:
    CarreauDynamics(Dynamics<T,Descriptor>* baseDynamics_);
    virtual T getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const;
    CarreauDynamics<T,Descriptor,N>* clone() const;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template<typename T, template<typename U> class Descriptor, int N>
class BGKCarreauDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    BGKCarreauDynamics();

    /// Clone the object on its dynamic type.
    virtual BGKCarreauDynamics<T,Descriptor,N>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template<typename T, template<typename U> class Descriptor, int N>
class RegularizedBGKCarreauDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    RegularizedBGKCarreauDynamics();

    /// Clone the object on its dynamic type.
    virtual RegularizedBGKCarreauDynamics<T,Descriptor,N>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

} // namespace plb

#endif  // VARIABLE_OMEGA_DYNAMICS_H
