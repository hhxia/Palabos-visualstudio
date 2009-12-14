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
 * Automatic creation of a 3D sparse multi-block -- generic code.
 */

#ifdef PLB_USE_CVMLCPP

#ifndef SPARSE_MULTI_BLOCK_3D_HH
#define SPARSE_MULTI_BLOCK_3D_HH

#include "parallelism/mpiManager.h"
#include "libraryInterfaces/CVMLCPP_sparseMultiBlock3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiLatticeInitializer3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"

namespace plb {

/// Use an octree-based tree-reader to create a multi-block and instantiate wall nodes,
///   and thus to set up a simulation representing the desired geometry.
template<typename T, template<typename U> class Descriptor>
MultiBlockLattice3D<T,Descriptor>* generateSparseMultiBlock (
        TreeReader3D<int>& treeReader,
        Dynamics<T,Descriptor>* backgroundDynamics )
{
    // The envelope width in the multi-block is chosen to be equal to the width of the
    //   neighborhood relation in the lattice.
    plint envelopeWidth = Descriptor<T>::vicinity;

    MultiBlockDistribution3D geometry (
            treeReader.getNx(), treeReader.getNy(), treeReader.getNz() );
    std::vector<Box3D> subBlocks = treeReader.getDomains();

    // Naive load-balancing strategy: assume all blocks have roughly the same size,
    //   and distribute them as evenly as possible over the processors.
    plint numProc      = defaultMultiBlockPolicy3D().getNumProcesses();
    plint blockPerProc = subBlocks.size() / numProc;

    plint iBlock = 0;
    for (int iProc=0; iProc<numProc; ++iProc) {
        plint numBlocks = blockPerProc;
        if (iProc < subBlocks.size() % numProc) {
            ++numBlocks;
        }
        for (int iLocalBlock=0; iLocalBlock<numBlocks; ++iLocalBlock) {
            geometry.addBlock(subBlocks[iBlock++], envelopeWidth, iProc);
        }
    }

    // Create the multi-block with default parameters, which are appropriate in
    //   serial as well as in parallel.
    MultiBlockLattice3D<T,Descriptor>* multiBlock =
        new MultiBlockLattice3D<T,Descriptor> (
            MultiBlockManagement3D( geometry,
                                    defaultMultiBlockPolicy3D().getThreadAttribution() ),
            defaultMultiBlockPolicy3D().getBlockCommunicator<T>(),
            defaultMultiBlockPolicy3D().getCombinedStatistics<T>(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T,Descriptor>(),
            backgroundDynamics );
    return multiBlock;
}


template<typename T, template<typename U> class Descriptor>
void dynamicsFromOctree( MultiBlockLattice3D<T,Descriptor>& multiBlock, 
                         TreeReader3D<int>& treeReader,
                         Dynamics<T,Descriptor>* dynamics,
                         bool flag )
{
    plint envelopeWidth = Descriptor<T>::vicinity;

    // The bool-mask is of type T, because during scalar-field - lattice couplings, the 
    //   data types of the two blocks must be identical.
    MultiScalarField3D<T> boolMask(multiBlock);

    // Get a full domain representation from the octree for all blocks which are local
    //   to the current processor.
    std::vector<plint> const& localBlocks = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint iLocal=0; iLocal<localBlocks.size(); ++iLocal) {
        plint iBlock = localBlocks[iLocal];
        ScalarField3D<T>& atomicMask =
            dynamic_cast<ScalarField3D<T>& > ( boolMask.getComponent(iBlock) );
        treeReader.createBooleanMask(iBlock, atomicMask, envelopeWidth);
    }

    // The treeReader has attributed values to the bulk of the atomic-blocks only. The method
    //   duplicateOverlaps() is invoked to copy them into the envelope.
    boolMask.getBlockCommunicator().duplicateOverlaps(boolMask);
    // Nodes are instantiated, based on the prescription found in the bool-mask.
    defineDynamics(multiBlock, boolMask, dynamics, flag);
}

}  // namespace plb

#endif  // SPARSE_MULTI_BLOCK_3D_HH

#endif  // PLB_USE_CVMLCPP
