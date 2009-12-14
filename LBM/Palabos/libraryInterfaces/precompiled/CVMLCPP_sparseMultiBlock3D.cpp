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
 * Automatic creation of a 3D sparse mmulti-block -- template instatiation.
 */
#ifdef PLB_USE_CVMLCPP

#include "libraryInterfaces/CVMLCPP_sparseMultiBlock3D.h"
#include "libraryInterfaces/CVMLCPP_sparseMultiBlock3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

template MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>*
    generateSparseMultiBlock<double,descriptors::D3Q19Descriptor> (
        TreeReader3D<int>& treeReader,
        Dynamics<double,descriptors::D3Q19Descriptor>* backgroundDynamics );

template void wallsFromOctree<double,descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& multiBlock, 
                      TreeReader3D<int>& treeReader );
}  // namespace plb

#endif  // PLB_USE_CVMLCPP
