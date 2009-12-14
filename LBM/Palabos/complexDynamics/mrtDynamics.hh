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
 * MRT dynamics -- generic implementation.
 */
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include "latticeBoltzmann/mrtTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class MRTdynamics *********************************************** */

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_),
      lambda(omega_)
{
    precompute_invM_S();
}

/** \param omega_ relaxation parameter related to the dynamic viscosity
 *  \param lambda_ relaxation parameter related to the bulk viscosity
 */
template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>::MRTdynamics(T omega_, T lambda_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_),
      lambda(lambda_)
{
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
MRTdynamics<T,Descriptor>* MRTdynamics<T,Descriptor>::clone() const {
    return new MRTdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    typedef mrtTemplates<T,Descriptor> mrtTemp;

    T rho;
    Array<T,Descriptor<T>::d> u;
    momentTemplates<T,Descriptor>::compute_rho_uLb(cell, rho, u);

    T uSqr = mrtTemp::mrtCollision(cell, rho, u, invM_S);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rho), uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::getParameter(plint whichParameter) const 
{
    switch (whichParameter) {
        case dynamicParams::omega_shear : return getOmega();
        case dynamicParams::omega_bulk  : return getLambda();
    };
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
        case dynamicParams::omega_shear : setOmega(value);
        case dynamicParams::omega_bulk  : setLambda(value);
    };
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::getOmega() const 
{
    return BasicBulkDynamics<T,Descriptor>::getOmega();
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::setOmega(T omega_) 
{
    BasicBulkDynamics<T,Descriptor>::setOmega(omega_);
    precompute_invM_S();
}

template<typename T, template<typename U> class Descriptor>
T MRTdynamics<T,Descriptor>::getLambda() const 
{
    return lambda;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::setLambda(T lambda_) 
{
    lambda = lambda_;
}

template<typename T, template<typename U> class Descriptor>
void MRTdynamics<T,Descriptor>::precompute_invM_S() {
    Array<T,Descriptor<T>::q> rt;  // relaxation times vector.
    for (plint iPop  = 0; iPop < Descriptor<T>::q; ++iPop)
    {
        rt[iPop] = Descriptor<T>::S[iPop];
    }
    for (plint iPop  = 0; iPop < Descriptor<T>::shearIndexes; ++iPop)
    {
        rt[Descriptor<T>::shearViscIndexes[iPop]] = this->getOmega();
    }
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
    {
        for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop)
        {
            invM_S[iPop][jPop] = T();
            for (plint kPop = 0; kPop < Descriptor<T>::q; ++kPop)
            {
                if (kPop == jPop)
                {
                    invM_S[iPop][jPop] += Descriptor<T>::invM[iPop][kPop] * 
                            rt[kPop];
                }
            }
        }
    }
}

}

#endif

