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
#ifndef THERMAL_DYNAMICS_HH
#define THERMAL_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "basicDynamics/thermalDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "core/latticeStatistics.h"

namespace plb {

/* *************** Class ThermalBulkDynamics ************************************ */

template<typename T, template<typename U> class Descriptor>
ThermalBulkDynamics<T,Descriptor>::ThermalBulkDynamics(T omega_)
  : BasicBulkDynamics<T,Descriptor>(omega_)
{ }

// TODO needs to be extended from iso-thermal to thermal case
template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::regularize (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{
    typedef Descriptor<T> L;
    cell[0] = this->computeEquilibrium(0, rhoBar, j, jSqr);
                + offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop=1; iPop<=L::q/2; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        cell[iPop+L::q/2] = this->computeEquilibrium(iPop+L::q/2, rhoBar, j, jSqr);
        T fNeq = offEquilibriumTemplatesImpl<T,L>::fromPiToFneq(iPop, PiNeq);
        cell[iPop] += fNeq;
        cell[iPop+L::q/2] += fNeq;
    }
}

template<typename T, template<typename U> class Descriptor>
T ThermalBulkDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    computeRhoBarJ(rhoBar,j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    return momentTemplates<T,Descriptor>::compute_theta(cell, rhoBar, jSqr);
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::computeDeviatoricStress (
    Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    computeRhoBarJ(rhoBar,j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T thetaBar = momentTemplates<T,Descriptor>::compute_theta(cell, rhoBar, jSqr)-(T)1;
    momentTemplates<T,Descriptor>::compute_thermal_PiNeq(cell, rhoBar, thetaBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& q ) const
{
    // TODO needs to be implemented
    PLB_PRECONDITION( false );
}

template<typename T, template<typename U> class Descriptor>
T ThermalBulkDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::get_eBar(cell);
}

// TODO needs to be extended from iso-thermal to thermal case
template<typename T, template<typename U> class Descriptor>
plint ThermalBulkDynamics<T,Descriptor>::numDecomposedVariables(plint order) const {
    plint numVariables =
                         // Order 0: density + velocity + fNeq
        ( order == 0 ) ? ( 1 + Descriptor<T>::d + Descriptor<T>::q )
                         // Order >=1: density + velocity + PiNeq
                       : ( 1 + Descriptor<T>::d + SymmetricTensor<T,Descriptor>::n );

    numVariables += Descriptor<T>::ExternalField::numScalars;
    return numVariables;
}

// TODO needs to be extended from iso-thermal to thermal case
template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));

    if (order==0) {
        decomposeOrder0(cell, rawData);
    }
    else {
        decomposeOrder1(cell, rawData);
    }
}

// TODO needs to be extended from iso-thermal to thermal case
template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    PLB_PRECONDITION( rawData.size() == numDecomposedVariables(order) );

    if (order==0) {
        recomposeOrder0(cell, rawData);
    }
    else {
        recomposeOrder1(cell, rawData);
    }
}

// TODO needs to be extended from iso-thermal to thermal case
template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{
    PLB_PRECONDITION( rawData.size()==numDecomposedVariables(order) );

    if (order==0) {
        rescaleOrder0(rawData, xDxInv, xDt);
    }
    else {
        rescaleOrder1(rawData, xDxInv, xDt);
    }
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::decomposeOrder0 (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d>& j;
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    momentTemplates<T,Descriptor>::compute_rhoBar_j(cell, rhoBar, j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        rawData[1+Descriptor<T>::d+iPop] =
            cell[iPop] - this->computeEquilibrium(iPop, rhoBar, j, jSqr);
    }

    T* externalField = cell.getExternal(0);
    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset+iExt] = *(externalField+iExt);
    }
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::decomposeOrder1 (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d>& j;
    Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq;
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);

    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    PiNeq.to_cArray(&rawData[1+Descriptor<T>::d]);

    T* externalField = cell.getExternal(0);
    int offset = 1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset+iExt] = *(externalField+iExt);
    }
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::recomposeOrder0 (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData ) const
{
    T rhoBar = rawData[0];
    Array<T,Descriptor<T>::d>& j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr)
                      + rawData[1+Descriptor<T>::d+iPop];
    }

    T* externalField = cell.getExternal(0);
    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *(externalField+iExt) = rawData[offset+iExt];
    }
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::recomposeOrder1 (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData ) const
{
    typedef Descriptor<T> L;

    T rhoBar;
    Array<T,Descriptor<T>::d>& j;
    Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq;

    rhoBar = rawData[0];
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    PiNeq.from_cArray(&rawData[1+Descriptor<T>::d]);

    cell[0] = this->computeEquilibrium(0, rhoBar, j, jSqr)
                + offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop=1; iPop<=L::q/2; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        cell[iPop+L::q/2] = this->computeEquilibrium(iPop+L::q/2, rhoBar, j, jSqr);
        T fNeq = offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, PiNeq);
        cell[iPop] += fNeq;
        cell[iPop+L::q/2] += fNeq;
    }

    T* externalField = cell.getExternal(0);
    int offset = 1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *(externalField+iExt) = rawData[offset+iExt];
    }
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::rescaleOrder0 (
        std::vector<T>& rawData, T xDxInv, T xDt ) const
{
    // Don't change rho (rawData[0]), because it is invariant

    // Change velocity, according to its units dx/dt
    T velScale = xDt * xDxInv;
    for (plint iVel=0; iVel<Descriptor<T>::d; ++iVel) {
        rawData[1+iVel] *= velScale;
    }

    // Change off-equilibrium, according to its units 1/dt
    T fNeqScale = xDt;
    for (plint iFneq=0; iFneq<Descriptor<T>::q; ++iFneq) {
        rawData[1+Descriptor<T>::d+iFneq] *= fNeqScale;
    }

    // Don't change external fields; their scaling must be taken care of
    //   in specialized versions of this class.
}

template<typename T, template<typename U> class Descriptor>
void ThermalBulkDynamics<T,Descriptor>::rescaleOrder1 (
        std::vector<T>& rawData, T xDxInv, T xDt ) const
{
    // Don't change rho (rawData[0]), because it is invariant

    // Change velocity, according to its units dx/dt
    T velScale = xDt * xDxInv;
    for (plint iVel=0; iVel<Descriptor<T>::d; ++iVel) {
        rawData[1+iVel] *= velScale;
    }

    // Change off-equilibrium stress, according to its units 1/dt
    T PiNeqScale = xDt;
    for (plint iPi=0; iPi<SymmetricTensor<T,Descriptor>::n; ++iPi) {
        rawData[1+Descriptor<T>::d+iPi] *= PiNeqScale;
    }

    // Don't change external fields; their scaling must be taken care of
    //   in specialized versions of this class.
}


}

#endif
