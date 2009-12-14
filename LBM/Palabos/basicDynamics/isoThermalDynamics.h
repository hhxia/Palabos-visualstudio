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
#ifndef ISO_THERMAL_DYNAMICS_H
#define ISO_THERMAL_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/// Common base iso-thermal (or athermal) bulk dynamics
template<typename T, template<typename U> class Descriptor>
class IsoThermalBulkDynamics : public BasicBulkDynamics<T,Descriptor> {
public:
    IsoThermalBulkDynamics(T omega_);

/* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Re-compute particle populations from the leading moments
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ************** */

    /// Returns 1, as a default value for isothermal flow.
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;

    /// Compute the deviatoric stress tensor ("off-equilibrium part of Pi")
    virtual void computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;

    /// Returns 0, as a default value for isothermal flow.
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

/// Implementation of O(Ma^2) BGK dynamics
template<typename T, template<typename U> class Descriptor>
class BGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    BGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual BGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of O(Ma^2) BGK dynamics, density and momentum taken from external scalars
template<typename T, template<typename U> class Descriptor>
class ExternalMomentBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ExternalMomentBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual ExternalMomentBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of incompressible BGK dynamics
template<typename T, template<typename U> class Descriptor>
class IncBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    IncBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual IncBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
template<typename T, template<typename U> class Descriptor>
class ConstRhoBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ConstRhoBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual ConstRhoBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Generic implementation of the Regularized BGK dynamics
/** This implementation is valid for isothermal models only.
 *  This model is substantially more stable than plain BGK, and has roughly
 *  the same efficiency. However, it cuts out the modes at higher Knudsen
 *  numbers and therefore cannot be used in the regime of rarefied gases.
 */
template<typename T, template<typename U> class Descriptor>
class RLBdynamics : public BulkCompositeDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    RLBdynamics(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object on its dynamic type.
    virtual RLBdynamics<T,Descriptor>* clone() const;

/* *************** Completion algorithm and base dynamics ************ */

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
};

/// Implementation of O(Ma^2) BGK dynamics with constant average density
/** Semantically, this class is equivalent to RLBdynamics< . , . , BGKdynamics<.,.> >,
 *  but the implementation is more efficient.
 */
template<typename T, template<typename U> class Descriptor>
class RegularizedBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    RegularizedBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual RegularizedBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of O(Ma^2) BGK dynamics, density and momentum taken from external scalars
template<typename T, template<typename U> class Descriptor>
class ExternalMomentRegularizedBGKdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ExternalMomentRegularizedBGKdynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual ExternalMomentRegularizedBGKdynamics<T,Descriptor>* clone() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
};

/// Implementation of O(Ma^2) BGK dynamics with adjustable speed of sound
template<typename T, template<typename U> class Descriptor>
class ChopardDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ChopardDynamics(T vs2_, T omega_);

    /// Clone the object on its dynamic type.
    virtual ChopardDynamics<T,Descriptor>* clone() const;

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
    /// Set local speed of sound
    void setVs2(T vs2_);
    /// Get local speed of sound
    T    getVs2() const;

private:
/* *************** Static implementation methods********************** */

    /// Implementation of collision operator
    static T chopardBgkCollision (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T vs2, T omega);
    /// Implementation of equilibrium
    static T chopardEquilibrium (
        plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T vs2);
private:
    T vs2;    ///< speed of sound
};

}

#endif
