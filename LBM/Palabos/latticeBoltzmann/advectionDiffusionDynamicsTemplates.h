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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_H
#define ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/offEquilibriumAdvectionDiffusionTemplates.h"

namespace plb {

template<typename T, class Descriptor> struct advectionDiffusionDynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct advectionDiffusionDynamicsTemplates {

static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& jEq) 
{
    return advectionDiffusionDynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma1_equilibrium(iPop, rhoBar, jEq);
}

static void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& jAdvDiff,
                       Array<T,Descriptor<T>::d> const& jEq)
{
    advectionDiffusionDynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::regularize(cell.getRawPopulations(), rhoBar, jAdvDiff, jEq);
}

static T no_corr_bgk_collision(
        Cell<T,Descriptor>& cell, T rhoBar, 
        Array<T,Descriptor<T>::d> const& jEq, T omega)
{
    return advectionDiffusionDynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::no_corr_bgk_collision(cell.getRawPopulations(), rhoBar, jEq, omega);
}

static T no_corr_rlb_collision(Cell<T,Descriptor>& cell, T rhoBar,
                             Array<T,Descriptor<T>::d> const& jEq,
                             Array<T,Descriptor<T>::d> const& jNeq,
                             T omega )
{
    return advectionDiffusionDynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::no_corr_rlb_collision(cell.getRawPopulations(), rhoBar, jEq, jNeq, omega);
}

};  // struct advectionDiffusionDynamicsTemplates


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct advectionDiffusionDynamicsTemplatesImpl 
{

static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T,Descriptor::d> const& jEq) 
{
    T c_j = Descriptor::c[iPop][0]*jEq[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*jEq[iD];
    }
    return Descriptor::t[iPop] * (rhoBar + Descriptor::invCs2 * c_j);
}

/// Regularization
static void regularize( Array<T,Descriptor::q>& f, T rhoBar,
                        Array<T,Descriptor::d> const& jAdvDiff,
                        Array<T,Descriptor::d> const& jEq )
{
    // Off-equilibrium j
    Array<T,Descriptor::d> jNeq;
    for (int iD=0; iD<Descriptor::d; ++iD) 
    {
        jNeq[iD] = jAdvDiff[iD] - jEq[iD];
    }

    // Regularize each population
    for (plint iPop=0; iPop<Descriptor::q; ++iPop) 
    {
        T fEq = bgk_ma1_equilibrium(iPop, rhoBar, jEq);
        T fNeq = (T)Descriptor::c[iPop][0] * jNeq[0];
        for (plint iD = 0; iD < Descriptor::d; ++iD)
        {
            fNeq += (T)Descriptor::c[iPop][iD] * jNeq[iD];
        }
        fNeq *= Descriptor::t[iPop] * Descriptor::invCs2;
        f[iPop] = fEq + fNeq;
    }
}

static T no_corr_bgk_collision(
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq, 
        T omega) 
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(jEq);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * advectionDiffusionDynamicsTemplatesImpl<T,Descriptor>::bgk_ma1_equilibrium (
                                iPop, rhoBar, jEq);
    }
    return jSqr*invRho*invRho;
}

static T no_corr_rlb_collision (
    Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq,
    Array<T,Descriptor::d> const& jNeq,T omega )
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(jEq);
    
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) 
    {
        f[iPop] = advectionDiffusionDynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma1_equilibrium(iPop, rhoBar, jEq);
        T fNeq = ((T)1-omega) *
                 offEquilibriumAdvectionDiffusionTemplatesImpl<T,Descriptor>::fromJtoFneq(iPop, jNeq);
        f[iPop] += fNeq;
    }
    return jSqr*invRho*invRho;
}

};  // struct advectionDiffusionDynamicsTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates2D.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates3D.h"

#endif
