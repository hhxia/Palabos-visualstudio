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
 * Specialized helper functions for advanced techniques around LB
 * implementations. They implement the physics of the first-order terms
 * of the Chapman-Enskog expansion and are useful whenever a transition 
 * from hydrodynamical variables (rho, u) to kinetic variables (f) si to
 * be implemented. Additionally, they are used for the implementation of
 * the stable RLB dynamics.
 *
 * This file is all about efficiency. The generic
 * template code is specialized for commonly used Lattices, so that a
 * maximum performance can be taken out of each case.
 */
#ifndef ADVECTION_DIFFUSION_OFF_EQUILIBRIUM_TEMPLATES_H
#define ADVECTION_DIFFUSION_OFF_EQUILIBRIUM_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"

namespace plb {

template<typename T, class Descriptor> struct offEquilibriumAdvectionDiffusionTemplatesImpl;

/// General first-order functions
template<typename T, template<typename U> class Descriptor>
struct offEquilibriumAdvectionDiffusionTemplates {

/// Compute off-equilibrium part of the f's from the current j.
/** Implements the following formula (with Einstein index contraction):
 * /f[ f_i^{neq} = t_i / (c_s^4) *
 *                 (c_{ia} j_a /f]
 * By Pi we mean the tensor computed from the off-equilibrium functions:
 * /f[ j_a = \sum c_{ia} f_i^{neq}
 *         = \sum c_{ia} f_i - \rho u_a  Id /f]
 * where u_a is an external vectorial field (the velocity for example)
 */
static T fromJtoFneq(plint iPop, Array<T,Descriptor<T>::d> const& jNeq) {
    return offEquilibriumAdvectionDiffusionTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::fromJtoFneq(iPop, jNeq);
}

};  // struct offEquilibriumAdvectionDiffusionTemplates

template<typename T, class Descriptor>
struct offEquilibriumAdvectionDiffusionTemplatesImpl {

static T fromJtoFneq (
    plint iPop, Array<T,Descriptor::d> const& jNeq )
{
    T fNeq = (T)Descriptor::c[iPop][0]*jNeq[0];
    for (plint iD=1; iD<Descriptor::d; ++iD) 
    {
        fNeq += (T)Descriptor::c[iPop][iD]*jNeq[iD];
    }
    fNeq *= Descriptor::t[iPop] * Descriptor::invCs2;
            
    return fNeq;
}

};  // struct offEquilibriumAdvectionDiffusionTemplates

}  // namespace plb

#endif
