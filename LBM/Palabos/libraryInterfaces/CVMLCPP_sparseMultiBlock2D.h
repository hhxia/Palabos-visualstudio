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

/* Main author: Daniel Lagrava
 */

/** \file
 * Automatic creation of a 2D sparse mmulti-block -- header file.
 */

#ifdef PLB_USE_CVMLCPP

#ifndef SPARSE_MULTI_BLOCK_2D_H
#define SPARSE_MULTI_BLOCK_2D_H

#include "core/dynamics.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "libraryInterfaces/CVMLCPP_treeReader2D.h"
#include <string>

namespace plb {

/// Use an octree-based tree-reader to create a multi-block and instantiate wall nodes,
///   and thus to set up a simulation representing the desired geometry.
template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>* generateSparseMultiBlock (
        TreeReader2D<int>& treeReader,
        Dynamics<T,Descriptor>* backgroundDynamics );

template<typename T, template<typename U> class Descriptor>
void wallsFromOctree( MultiBlockLattice2D<T,Descriptor>& multiBlock, 
                      TreeReader2D<int>& treeReader );

}  // namespace plb

#endif  // SPARSE_MULTI_BLOCK_2D_H

#endif  // PLB_USE_CVMLCPP
