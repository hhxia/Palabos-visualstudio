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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_DYNAMICS_HH
#define SMAGORINSKY_DYNAMICS_HH

#include "complexDynamics/smagorinskyDynamics.h"
#include "core/util.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct SmagoOperations {
    static T computePrefactor(T omega0, T cSmago) {
        return (T)0.5 * util::sqr(cSmago*omega0*Descriptor<T>::invCs2);
    }
    static T recomputePrefactor(T oldOmega0, T newOmega0, T oldPrefactor) {
        return oldPrefactor * util::sqr(newOmega0/oldOmega0);
    }
    static T computeOmega(T omega0, T preFactor, T rhoBar, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq)
    {
        T PiNeqNormSqr = SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq);
        T PiNeqNorm    = sqrt(PiNeqNormSqr);
        T alpha        = preFactor * Descriptor<T>::invRho(rhoBar);
        T linearTerm   = alpha*PiNeqNorm;
        T squareTerm   = (T)2*alpha*alpha*PiNeqNormSqr;
        // In the following formula, the square-root appearing in the explicit form of
        //   omega is developed to second-order.
        return omega0*(1-linearTerm+squareTerm);
    }
};

template<typename T, template<typename U> class Descriptor>
SmagorinskyDynamics<T,Descriptor>::SmagorinskyDynamics(Dynamics<T,Descriptor>* baseDynamics_, T omega0_, T cSmago_)
    : OmegaFromPiDynamics<T,Descriptor>(baseDynamics_),
      omega0(omega0_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0_, cSmago_))
{ }

template<typename T, template<typename U> class Descriptor>
T SmagorinskyDynamics<T,Descriptor>::getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const
{
    return SmagoOperations<T,Descriptor>::computeOmega(omega0, preFactor, rhoBar, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyDynamics<T,Descriptor>* SmagorinskyDynamics<T,Descriptor>::clone() const {
    return new SmagorinskyDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyDynamics<T,Descriptor>::setOmega(T omega0_)
{
    preFactor = SmagoOperations<T,Descriptor>::recomputePrefactor (
            this->getOmega(), omega0_, preFactor);
}

/* *************** Class SmagorinskyBGKdynamics ************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyBGKdynamics<T,Descriptor>::SmagorinskyBGKdynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0_,cSmago_))
{ }

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::setOmega(T omega0_)
{
    preFactor = SmagoOperations<T,Descriptor>::recomputePrefactor (
            this->getOmega(), omega0_, preFactor);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyBGKdynamics<T,Descriptor>* SmagorinskyBGKdynamics<T,Descriptor>::clone() const {
    return new SmagorinskyBGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyBGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, omega);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyBGKdynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}


/* *************** Class SmagorinskyRegularizedDynamics ************************************ */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T,Descriptor>::SmagorinskyRegularizedDynamics (
        T omega0_, T cSmago_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega0_),
      omega0(omega0_),
      preFactor(SmagoOperations<T,Descriptor>::computePrefactor(omega0_,cSmago_))
{ }

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::setOmega(T omega0_)
{
    preFactor = SmagoOperations<T,Descriptor>::recomputePrefactor (
            this->getOmega(), omega0_, preFactor);
}

template<typename T, template<typename U> class Descriptor>
SmagorinskyRegularizedDynamics<T,Descriptor>* SmagorinskyRegularizedDynamics<T,Descriptor>::clone() const
{
    return new SmagorinskyRegularizedDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void SmagorinskyRegularizedDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
    T omega = SmagoOperations<T,Descriptor>::computeOmega (
            omega0, preFactor, rhoBar, PiNeq );
    IsoThermalBulkDynamics<T,Descriptor>::setOmega(omega);
    T uSqr = dynamicsTemplates<T,Descriptor>::rlb_collision (
                 cell, rhoBar, j, PiNeq, omega );
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T SmagorinskyRegularizedDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar ) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

} // namespace plb

#endif  // SMAGORINSKY_DYNAMICS_HH
