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
 * Helper functions for the implementation of force terms in LB dynamics.
 * This file is all about efficiency. The generic template code is specialized
 * for commonly used Lattices, so that a maximum performance can be taken out of
 * each case.
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_H
#define EXTERNAL_FORCE_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct externalForceTemplates {

/// Add a force term after BGK collision
static void addGuoForce(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u, T omega, T amplitude) {
    static const int forceBeginsAt = Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = cell.getExternal(forceBeginsAt);
    for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) 
    {
        T c_u = VectorTemplate<T,Descriptor>::scalarProduct(Descriptor<T>::c[iPop],u);
        c_u *= Descriptor<T>::invCs2 * Descriptor<T>::invCs2;
        T forceTerm = T();
        for (int iD=0; iD < Descriptor<T>::d; ++iD) {
            forceTerm +=
                (   ((T)Descriptor<T>::c[iPop][iD]-u[iD]) * Descriptor<T>::invCs2
                     + c_u * (T)Descriptor<T>::c[iPop][iD]
                )
                * force[iD];
        }
        forceTerm *= Descriptor<T>::t[iPop];
        forceTerm *= 1-omega/(T)2;
        forceTerm *= amplitude;
        cell[iPop] += forceTerm;
    }
}

};  // struct externalForceTemplates

}  // namespace plb

#include "latticeBoltzmann/externalForceTemplates2D.h"
#include "latticeBoltzmann/externalForceTemplates3D.h"

#endif
