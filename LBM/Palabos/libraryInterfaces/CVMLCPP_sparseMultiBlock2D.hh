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
 * Automatic creation of a 2D sparse mmulti-block -- generic code.
 */

#ifdef PLB_USE_CVMLCPP

#ifndef SPARSE_MULTI_BLOCK_2D_HH
#define SPARSE_MULTI_BLOCK_2D_HH

#include "parallelism/mpiManager.h"
#include "libraryInterfaces/CVMLCPP_sparseMultiBlock2D.h"
#include "multiBlock/multiDataCouplingWrapper2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/util.h"

namespace plb {

/// Define dynamics for wall nodes, which are located with help of a boolean mask.
/** Note that the boolean mask is of type MultiScalarField2D<T> instead of MultiScalarField2D<bool>,
 *  because coupling lattices of different data types is not possible.
 */
template<typename T, template<typename U> class Descriptor>
class WallFromMaskFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor> {
public:
    WallFromMaskFunctional2D(Dynamics<T,Descriptor>* wallDynamics_)
        : wallDynamics(wallDynamics_)
    { }
    WallFromMaskFunctional2D(WallFromMaskFunctional2D<T,Descriptor> const& rhs)
        : wallDynamics(rhs.wallDynamics->clone())
    { }
    WallFromMaskFunctional2D<T,Descriptor>& operator= (
            WallFromMaskFunctional2D<T,Descriptor> const& rhs )
    {
        delete wallDynamics; wallDynamics = rhs.wallDynamics->clone();
    }

    ~WallFromMaskFunctional2D() {
        delete wallDynamics;
    }
    virtual void process (
            Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                          ScalarField2D<T>& mask )
    {
        Dot2D offset = computeRelativeDisplacement(lattice, mask);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                bool isLiquid = (bool) util::roundToInt (
                        mask.get(iX+offset.x, iY+offset.y) );
                if (!isLiquid) {
                    lattice.attributeDynamics(iX,iY, wallDynamics->clone());
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Dynamics needs to be instantiated everywhere, including envelope.
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
        isWritten[1] = false;
    }
    virtual WallFromMaskFunctional2D<T,Descriptor>* clone() const 
    {
        return new WallFromMaskFunctional2D<T,Descriptor>(*this);
    }
private:
    Dynamics<T,Descriptor>* wallDynamics;
};


/// Use an octree-based tree-reader to create a multi-block and instantiate wall nodes,
///   and thus to set up a simulation representing the desired geometry.
template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>* generateSparseMultiBlock (
        TreeReader2D<int>& treeReader,
        Dynamics<T,Descriptor>* backgroundDynamics )
{
    // The envelope width in the multi-block is chosen to be equal to the width of the
    //   neighborhood relation in the lattice.
    plint envelopeWidth = Descriptor<T>::vicinity;


    MultiBlockDistribution2D geometry (
            treeReader.getNx(), treeReader.getNy() );
    std::vector<Box2D> subBlocks = treeReader.getDomains();

    // Naive load-balancing strategy: assume all blocks have roughly the same size,
    //   and distribute them as evenly as possible over the processors.
    plint numProc      = global::mpi().getSize();
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
    MultiBlockLattice2D<T,Descriptor>* multiBlock =
        new MultiBlockLattice2D<T,Descriptor> (
            MultiBlockManagement2D( geometry,
                                    defaultMultiBlockPolicy2D().getThreadAttribution() ),
            defaultMultiBlockPolicy2D().getBlockCommunicator<T>(),
            defaultMultiBlockPolicy2D().getCombinedStatistics<T>(),
            defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(),
            backgroundDynamics );
    return multiBlock;
}


template<typename T, template<typename U> class Descriptor>
void wallsFromOctree( MultiBlockLattice2D<T,Descriptor>& multiBlock, 
                      TreeReader2D<int>& treeReader )
{
    plint envelopeWidth = Descriptor<T>::vicinity;

    // The bool-mask is of type T, because during scalar-field - lattice couplings, the 
    //   data types of the two blocks must be identical.
    MultiScalarField2D<T> boolMask(multiBlock);

    // Get a full domain representation from the octree for all blocks which are local
    //   to the current processor.
    std::vector<plint> const& localBlocks = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint iLocal=0; iLocal<localBlocks.size(); ++iLocal) {
        plint iBlock = localBlocks[iLocal];
        ScalarField2D<T>& atomicMask =
            dynamic_cast<ScalarField2D<T>& > ( boolMask.getComponent(iBlock) );
        treeReader.createBooleanMask(iBlock, atomicMask, envelopeWidth);
    }

    // The treeReader has attributed values to the bulk of the atomic-blocks only. The method
    //   duplicateOverlaps() is invoked to copy them into the envelope.
    boolMask.getBlockCommunicator().duplicateOverlaps(boolMask);

    // Finally, bounce-back nodes are instantiated, based on the prescription found in the bool-mask.
    applyProcessingFunctional ( new WallFromMaskFunctional2D<T,Descriptor>(new BounceBack<T,Descriptor>),
                                multiBlock.getBoundingBox(),
                                multiBlock, boolMask );

}

}  // namespace plb

#endif  // SPARSE_MULTI_BLOCK_2D_HH

#endif  // PLB_USE_CVMLCPP
