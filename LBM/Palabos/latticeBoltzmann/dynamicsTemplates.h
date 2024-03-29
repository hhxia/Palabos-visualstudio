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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef DYNAMICS_TEMPLATES_H
#define DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

template<typename T, class Descriptor> struct dynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct dynamicsTemplates {

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

static T bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision(cell.getRawPopulations(), rhoBar, j, omega);
}

static T bgk_inc_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_inc_collision(cell.getRawPopulations(), rhoBar, j, omega);
}

static T rlb_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                       Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T omega )
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::rlb_collision(cell.getRawPopulations(), rhoBar, j, PiNeq, omega);
}

static T bgk_ma2_constRho_collision(Cell<T,Descriptor>& cell,
                                T rhoBar, Array<T,Descriptor<T>::d> const& j, T ratioRho, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
               ::bgk_ma2_constRho_collision(cell.getRawPopulations(), rhoBar, j, ratioRho, omega);
}

};  // struct dynamicsTemplates


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct dynamicsTemplatesImpl {

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/(T)2 * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr*invRho*invRho;
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    // In incompressible BGK, the Ma^2 term is preceeded by 1 instead of 1/rho.
    T invRho = (T)1;
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr;
}

static T rlb_collision (
    Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j,
    Array<T,SymmetricTensorImpl<T,Descriptor::d>::n> const& PiNeq, T omega )
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    f[0] = dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) *
                        offEquilibriumTemplatesImpl<T,Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop=1; iPop<=Descriptor::q/2; ++iPop) {
        f[iPop] = dynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
        f[iPop+Descriptor::q/2] = dynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma2_equilibrium(iPop+Descriptor::q/2, rhoBar, invRho, j, jSqr);
        T fNeq = ((T)1-omega) *
                 offEquilibriumTemplatesImpl<T,Descriptor>::fromPiToFneq(iPop, PiNeq);
        f[iPop] += fNeq;
        f[iPop+Descriptor::q/2] += fNeq;
    }
    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T ratioRho, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        T feq = dynamicsTemplatesImpl<T,Descriptor>::
                    bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + Descriptor::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}

};  // struct dynamicsTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/dynamicsTemplates2D.h"
#include "latticeBoltzmann/dynamicsTemplates3D.h"

#endif
