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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef THERMAL_DYNAMICS_H
#define THERMAL_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/// Common base iso-thermal (or athermal) bulk dynamics
template<typename T, template<typename U> class Descriptor>
class ThermalBulkDynamics : public BasicBulkDynamics<T,Descriptor> {
public:
    ThermalBulkDynamics(T omega_);

/* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ************** */

    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;

    /// Compute the deviatoric stress tensor ("off-equilibrium part of Pi")
    virtual void computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;

    /// Compute the heat flux in lattice units
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

/* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    /** In the present implementation, the decomposition is carried out up to order-1 in the
     *    Chapman-Enskog expansion. Example: Take the D2Q9 lattice. A decomposition means:
     *    - At order 0: Decompose into rho, u, and fNeq (1+2+9=12 variables)
     *    - At order 1: Decompose into rho, u, and PiNeq (1+2+3=6 variables)
     *    - At higher order: Decompose according to order 1.
     */
    virtual plint numDecomposedVariables(plint order) const;

    /// Decompose from population representation into moment representation.
    /**   \sa numDecomposedVariables()
     */
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Recompose from moment representation to population representation.
    /**   \sa numDecomposedVariables()
     *    This process is also known as "regularization step", and this function is therefore
     *    equivalent to regularize(), although one or the other function may be more useful
     *    in a specific context, due to the form of the parameters.
     */
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Change the space and time scales of the variables in moment representation.
    /**   \sa numDecomposedVariables()
     *    \param xDxInv Inverse of the factor by which space scale is multiplied.
     *    \param xDt Factor by which time scale is multiplied.
     */
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    

/* *************** Additional moments, intended for internal use ***** */

    /// Returns 0, as a default value for isothermal flow.
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;

private:
    virtual void decomposeOrder0(Cell<T,Descriptor> const& cell, std::vector<T>& rawData) const;
    virtual void decomposeOrder1(Cell<T,Descriptor> const& cell, std::vector<T>& rawData) const;
    virtual void recomposeOrder0(Cell<T,Descriptor>& cell, std::vector<T> const& rawData) const;
    virtual void recomposeOrder1(Cell<T,Descriptor>& cell, std::vector<T> const& rawData) const;
    virtual void rescaleOrder0(std::vector<T>& rawData, T xDxInv, T xDt) const;
    virtual void rescaleOrder1(std::vector<T>& rawData, T xDxInv, T xDt) const;
};

}

#endif
