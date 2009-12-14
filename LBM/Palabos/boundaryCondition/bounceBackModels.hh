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
 * can be instantiated -- generic implementation.
 */
#ifndef BOUNCE_BACK_MODELS_HH
#define BOUNCE_BACK_MODELS_HH

#include "boundaryCondition/bounceBackModels.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include "core/array.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class MomentumExchangeBounceBack ************************* */

template<typename T, template<typename U> class Descriptor>
MomentumExchangeBounceBack<T,Descriptor>::MomentumExchangeBounceBack (
        Array<plint, Descriptor<T>::d> forceIds_, T rho_)
    : fluidDirections(),
      forceIds(forceIds_),
      rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
MomentumExchangeBounceBack<T,Descriptor>* MomentumExchangeBounceBack<T,Descriptor>::clone() const {
    return new MomentumExchangeBounceBack<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    enum { numDim=Descriptor<T>::d };

    // Before collision, the amount of exchanged momentum is computed.
    Array<T,numDim> momentum;
    momentum.resetToZero();
    // Sum over all populations which are incoming from the fluid.
    for (pluint iFluid=0; iFluid < fluidDirections.size(); ++iFluid) {
        plint iPop = fluidDirections[iFluid];
        for (plint iD=0; iD<numDim; ++iD) {
            // The momentum contribution is multiplied by two:
            //   One contribution for the momentum loss into the obstacle,
            //   and one contribution for the momentum gain in the subsequent
            //   after-bounce-back streaming step.
            momentum[iD] += 2.*Descriptor<T>::c[iPop][iD]*cell[iPop];
        }
    }
    // Add the momentum exchange for this cell to the total balance due to the
    //   obstacle.
    if (cell.takesStatistics()) {
        for (plint iD=0; iD<numDim; ++iD) {
            // Add a negative sign, to get the force acting on the obstacle,
            //   and not on the fluid.
            statistics.gatherSum(forceIds[iD], -momentum[iD]);
        }
    }

    // Finally, do the bounce-back operation, which replaces the usual collision.
    for (plint iPop=1; iPop <= Descriptor<T>::q/2; ++iPop) {
        std::swap(cell[iPop], cell[iPop+Descriptor<T>::q/2]);
    }
}

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                               T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return rho;
}

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::getOmega() const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T MomentumExchangeBounceBack<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint MomentumExchangeBounceBack<T,Descriptor>::numDecomposedVariables(plint order) const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    cell.serialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    PLB_PRECONDITION( rawData.size() == numDecomposedVariables(order) );
    cell.unSerialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }

template<typename T, template<typename U> class Descriptor>
void MomentumExchangeBounceBack<T,Descriptor>::setFluidDirections (
        std::vector<plint> const& fluidDirections_ )
{
    fluidDirections = fluidDirections_;
}

template<typename T, template<typename U> class Descriptor>
std::vector<plint> const&  MomentumExchangeBounceBack<T,Descriptor>::getFluidDirections() const
{
    return fluidDirections;
}

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_HH
