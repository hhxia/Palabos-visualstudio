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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

/// Implementation of the entropic collision step
template<typename T, template<typename U> class Descriptor>
class MRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    MRTdynamics(T omega_);
    MRTdynamics(T omega_, T lambda_);

    /// Clone the object on its dynamic type.
    virtual MRTdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

/* *************** Configurable parameters *************************** */

    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;
    /// Get local relaxation parameter for shear viscosity
    virtual T getOmega() const;
    /// Set local relaxation parameter for shear viscosity
    virtual void setOmega(T omega_);
    /// Get local relaxation parameter for bulk viscosity
    T getLambda() const;
    /// Set local relaxation parameter for bulk viscosity
    void setLambda(T lambda_);
private:
    /// Precompute relaxation time matrix
    void precompute_invM_S();
private:
    T invM_S[Descriptor<T>::q][Descriptor<T>::q]; // relaxation time matrix.
    T omega; // the shear viscosity relaxation time
    T lambda;// the bulk viscosity relaxation time
};

}

#endif
