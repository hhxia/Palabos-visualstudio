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
#ifndef BOUNCE_BACK_MODELS_H
#define BOUNCE_BACK_MODELS_H

#include "core/dynamics.h"
#include <vector>

namespace plb {

/// Implementation of "full-way bounce-back" dynamics which computes momentum exchange.
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics.
 *
 * Once they have been instantiated, MomentumExchangeBounceBack objects
 * do _not_ compute the momentum exchange right away. First, they need
 * to be initialized through a function call to one of the flavors of
 * initializeMomentumExchange().
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Descriptor>
class MomentumExchangeBounceBack : public Dynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */

    /** You may fix a fictitious density value on bounce-back nodes via the constructor.
     *  \param forceIds_ Contains identifiers to access the reductive variables
     *  in the BlockStatistics objects, and to update the value of the total momentum
     *  exchange on the obstacle. The value of the force-ids must be determined
     *  previously by the user through a call to the method subscribeSum() of
     *  the BlockStatistics object in the used lattice. 
     */
    MomentumExchangeBounceBack(Array<plint, Descriptor<T>::d> forceIds_, T rho_=T() );

    /// Clone the object on its dynamic type.
    virtual MomentumExchangeBounceBack<T,Descriptor>* clone() const;

/* *************** Collision, Equilibrium, and Non-equilibrium ******* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);

    /// Yields 0
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    /// Does nothing
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ************** */

    /// Yields fictitious density
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Yields 0
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Yields 0
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

    /// Does nothing
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Access to Dynamics variables, e.g. omega ********** */

    /// Yields 0
    virtual T getOmega() const;

    /// Does nothing
    virtual void setOmega(T omega_);

/* *************** Switch between population and moment representation ****** */

    /// Yields Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Nothing happens here.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;

/* *************** Additional moments, intended for internal use ***** */

    /// Yields fictitious density
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Yields fictitious density and 0
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;
    /// Yields 0
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
public:
    /// Define the directions which point from the current cell into a fluid node.
    void setFluidDirections(std::vector<plint> const& fluidDirections_);
    /// Get the directions which point from the current cell into a fluid node.
    std::vector<plint> const& getFluidDirections() const;
private:
    std::vector<plint> fluidDirections;
    Array<plint, Descriptor<T>::d> forceIds;
    T rho;
};

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_H
