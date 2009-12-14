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
 * Helper functions for domain initialization -- header file.
 */
#include "core/globalDefs.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/dynamics.h"

#ifndef ATOMIC_LATTICE_INITIALIZER_3D_H
#define ATOMIC_LATTICE_INITIALIZER_3D_H

namespace plb {

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& boolMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, bool whichFlag );

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& boolMask,
                     Dynamics<T,Descriptor>* dynamics, bool whichFlag );


template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& intMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, int whichFlag );

template<typename T, template<typename U> class Descriptor>
void defineDynamics( BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& intMask,
                     Dynamics<T,Descriptor>* dynamics, int whichFlag );

}  // namespace plb

#endif  // ATOMIC_LATTICE_INITIALIZER_3D_H
