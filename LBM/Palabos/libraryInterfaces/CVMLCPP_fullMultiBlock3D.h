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
 * Automatic creation of a 3D non-sparse multi-block -- header file.
 */

#ifdef PLB_USE_CVMLCPP

#ifndef FULL_MULTI_BLOCK_3D_H
#define FULL_MULTI_BLOCK_3D_H

#include "core/dynamics.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "libraryInterfaces/CVMLCPP_treeReader3D.h"
#include <cvmlcpp/base/Matrix>
#include <string>

namespace plb {

template<typename T>
ScalarField3D<T>* boolMaskFromSTL(std::string stlFile, plint N);

template<typename T, template <typename U> class Descriptor>
void dynamicsFromBoolMask (
        ScalarField3D<T>& boolMask,
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Dynamics<T,Descriptor>* dynamics, bool flag );

}  // namespace plb

#endif  // FULL_MULTI_BLOCK_3D_H

#endif  // PLB_USE_CVMLCPP
