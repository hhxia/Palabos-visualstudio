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

#ifndef SMAGORINSKY_DYNAMICS_H
#define SMAGORINSKY_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "complexDynamics/variableOmegaDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class SmagorinskyDynamics : public OmegaFromPiDynamics<T,Descriptor> {
public:
    SmagorinskyDynamics(Dynamics<T,Descriptor>* baseDynamics_, T omega0_, T cSmago_);
    /// Clone the object on its dynamic type.
    SmagorinskyDynamics<T,Descriptor>* clone() const;
    /// Modify the value of omega, using the Smagorinsky algorithm based on omega0.
    virtual T getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
private:
    T omega0;    //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T preFactor; //< A factor depending on the Smagorinky constant, used to compute the effective viscosity.
};

template<typename T, template<typename U> class Descriptor>
class SmagorinskyBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    SmagorinskyBGKdynamics(T omega0_, T cSmago_);
    /// Clone the object on its dynamic type.
    virtual SmagorinskyBGKdynamics<T,Descriptor>* clone() const;
    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
private:
    T omega0;    //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T preFactor; //< A factor depending on the Smagorinky constant, used to compute the effective viscosity.
};

template<typename T, template<typename U> class Descriptor>
class SmagorinskyRegularizedDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
    SmagorinskyRegularizedDynamics(T omega0_, T cSmago_);

    /// Clone the object on its dynamic type.
    virtual SmagorinskyRegularizedDynamics<T,Descriptor>* clone() const;
    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
    /// With this method, you can modify the constant value omega0 (not the actual value of omega,
    ///  which is computed during run-time from omega0 and the local strain-rate).
    virtual void setOmega(T omega_);
private:
    T omega0;    //< "Laminar" relaxation parameter, used when the strain-rate is zero.
    T preFactor; //< A factor depending on the Smagorinky constant, used to compute the effective viscosity.
};

} // namespace plb

#endif  // SMAGORINSKY_DYNAMICS_H
