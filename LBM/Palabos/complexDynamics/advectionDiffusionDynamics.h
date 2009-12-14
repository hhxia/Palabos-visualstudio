/* This file is part of the Palabos library.
 * Copyright (C) 2009 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
e* published by the Free Software Foundation, either version 3 of the
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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_H
#define ADVECTION_DIFFUSION_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {
    
/// Common base iso-thermal (or athermal) bulk dynamics
template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionDynamics : public BasicBulkDynamics<T,Descriptor> {
public:
    AdvectionDiffusionDynamics(T omega_);

/* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Additional moments, intended for internal use ***** */

    /// Returns 0, as a default value for isothermal flow.
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;

/* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    virtual plint numDecomposedVariables(plint order) const { return 0; }

    /// Decompose from population representation into moment representation.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const { }

    /// Recompose from moment representation to population representation.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const { }

    /// Change the space and time scales of the variables in moment representation.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const { }
    
};

/// Regularized Advection-Diffusion dynamics
/** It uses the regularized approximation that can be found in
 *   the thesis of J. Latt (2007).
 */
template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionRLBdynamics : public AdvectionDiffusionDynamics <T,Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionRLBdynamics(T omega_);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionRLBdynamics<T,Descriptor>* clone() const;
    /// Collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics );
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// BGK Advection-Diffusion dynamics
/** This approach contains a slight error in the diffusion
 *  term.
 */
template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionBGKdynamics : public AdvectionDiffusionDynamics <T,Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionBGKdynamics(T omega_);
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionBGKdynamics<T,Descriptor>* clone() const;
    /// Collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics );
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

/* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    virtual plint numDecomposedVariables(plint order) const { return 0; }

    /// Decompose from population representation into moment representation.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const { }

    /// Recompose from moment representation to population representation.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const { }

    /// Change the space and time scales of the variables in moment representation.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const { }
    
};

} // namespace plb

#endif
