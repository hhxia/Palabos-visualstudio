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
 * A 2D multiblock lattice -- generic implementation.
 */
#ifndef MULTI_BLOCK_LATTICE_2D_HH
#define MULTI_BLOCK_LATTICE_2D_HH

#include "multiBlock/multiBlockLattice2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace plb {


////////////////////// Class MultiBlockLattice2D /////////////////////////

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::MultiBlockLattice2D (
        MultiBlockManagement2D const& multiBlockManagement_,
        BlockCommunicator2D<T>* blockCommunicator_,
        CombinedStatistics<T>* combinedStatistics_,
        MultiCellAccess2D<T,Descriptor>* multiCellAccess_,
        Dynamics<T,Descriptor>* backgroundDynamics )
    : MultiBlock2D<T>(multiBlockManagement_, blockCommunicator_, combinedStatistics_ ),
      multiCellAccess(multiCellAccess_)
{
    allocateBlocks(backgroundDynamics);
    eliminateStatisticsInEnvelope();
    this->evaluateStatistics(); // Reset statistics to default.
}

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::MultiBlockLattice2D (
        plint nx, plint ny,
        Dynamics<T,Descriptor>* backgroundDynamics )
    : MultiBlock2D<T>(nx,ny),
      multiCellAccess(defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>())
{
    allocateBlocks(backgroundDynamics);
    eliminateStatisticsInEnvelope();
    this->evaluateStatistics(); // Reset statistics to default.
}

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::~MultiBlockLattice2D() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        delete blockLattices[iBlock];
    }
    delete multiCellAccess;
}

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::MultiBlockLattice2D(MultiBlockLattice2D<T,Descriptor> const& rhs)
    : BlockLatticeBase2D<T,Descriptor>(rhs),
      MultiBlock2D<T>(rhs),
      multiCellAccess(rhs.multiCellAccess->clone())
{
    for (plint iBlock=0; iBlock<getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            blockLattices.push_back (
                    new BlockLattice2D<T,Descriptor>(*rhs.blockLattices[iBlock]) );
        }
        else {
            blockLattices.push_back( 0 );
        }
    }
}

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::MultiBlockLattice2D(MultiBlock2D<T> const& rhs)
    : MultiBlock2D<T>(rhs),
      multiCellAccess(defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>())
{
    allocateBlocks(new NoDynamics<T,Descriptor>());
    eliminateStatisticsInEnvelope();
}

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>::MultiBlockLattice2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop)
    : MultiBlock2D<T>(rhs, subDomain, crop),
      multiCellAccess(defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>())
{
    allocateBlocks(new NoDynamics<T,Descriptor>());
    eliminateStatisticsInEnvelope();
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>& MultiBlockLattice2D<T,Descriptor>::get(plint iX, plint iY) {
    return multiCellAccess -> getDistributedCell(iX,iY, this->getMultiBlockManagement(), blockLattices);
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor> const& MultiBlockLattice2D<T,Descriptor>::get(plint iX, plint iY) const {
    return multiCellAccess -> getDistributedCell(iX,iY, this->getMultiBlockManagement(), blockLattices);
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::specifyStatisticsStatus (Box2D domain, bool status) {
    Box2D inters;
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        if (intersect(domain, params.getBulk(), inters ) ) {
            inters = params.toLocal(inters);
            blockLattices[iBlock] -> specifyStatisticsStatus(inters, status);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::collide(Box2D domain) {
    Box2D inters;
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        if (intersect(domain, params.getEnvelope(), inters ) ) {
            inters = params.toLocal(inters);
            blockLattices[iBlock] -> collide(inters);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::collide() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        blockLattices[iBlock] -> collide();
    }
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::stream(Box2D domain) {
    Box2D inters;
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        if (intersect(domain, params.getNonPeriodicEnvelope(), inters ) ) {
            inters = params.toLocal(inters);
            blockLattices[iBlock] -> stream(inters);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
Box2D MultiBlockLattice2D<T,Descriptor>::extendPeriodic(Box2D const& box, plint envelopeWidth) const
{
    Box2D boundingBox(this->getBoundingBox());
    Box2D periodicBox(box);
    bool periodicX = this->periodicity().get(0);
    bool periodicY = this->periodicity().get(1);
    if (periodicX) {
        if (periodicBox.x0 == boundingBox.x0) {
            periodicBox.x0 -= envelopeWidth;
        }
        if (periodicBox.x1 == boundingBox.x1) {
            periodicBox.x1 += envelopeWidth;
        }
    }
    if (periodicY) {
        if (periodicBox.y0 == boundingBox.y0) {
            periodicBox.y0 -= envelopeWidth;
        }
        if (periodicBox.y1 == boundingBox.y1) {
            periodicBox.y1 += envelopeWidth;
        }
    }
    return periodicBox;
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::stream() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        // Stream must be applied to full domain, including currently active envelopes.
        Box2D domain = extendPeriodic(params.getNonPeriodicEnvelope(), params.getEnvelopeWidth());
        blockLattices[iBlock] -> stream( params.toLocal(domain) );
    }

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::collideAndStream(Box2D domain) {
    Box2D inters;
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        if (intersect(domain, params.getNonPeriodicEnvelope(), inters ) ) {
            inters = params.toLocal(inters);
            blockLattices[iBlock] -> collideAndStream(inters);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::collideAndStream() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        BlockParameters2D const& params = getParameters(iBlock);
        // CollideAndStream must be applied to full domain,
        //   including currently active envelopes.
        Box2D domain = extendPeriodic(params.getNonPeriodicEnvelope(), params.getEnvelopeWidth());
        blockLattices[iBlock] -> collideAndStream( params.toLocal(domain) );
    }
    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::incrementTime() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock <  this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        blockLattices[iBlock] -> incrementTime();
    }
    this->getTimeCounter().incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::allocateBlocks(Dynamics<T,Descriptor>* backgroundDynamics )
{
    for (plint iBlock=0; iBlock<getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            Box2D const& envelope = this->getMultiBlockManagement().getEnvelope(iBlock);
            BlockLattice2D<T,Descriptor>* newLattice =
                new BlockLattice2D<T,Descriptor> (
                        envelope.getNx(), envelope.getNy(),
                        backgroundDynamics->clone() );
            newLattice -> setLocation(Dot2D(envelope.x0, envelope.y0));
            blockLattices.push_back(newLattice);
        }
        else {
            blockLattices.push_back( 0 );
        }
    }
    delete backgroundDynamics;
}

template<typename T, template<typename U> class Descriptor>
void MultiBlockLattice2D<T,Descriptor>::eliminateStatisticsInEnvelope() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        plint envelopeWidth = getParameters(iBlock).getEnvelopeWidth();
        BlockLattice2D<T,Descriptor>& block = *blockLattices[iBlock];
        plint maxX = block.getNx()-1;
        plint maxY = block.getNy()-1;

        block.specifyStatisticsStatus(Box2D(0, maxX, 0, envelopeWidth-1), false);
        block.specifyStatisticsStatus(Box2D(0, maxX, maxY-envelopeWidth+1, maxY), false);
        block.specifyStatisticsStatus(Box2D(0, envelopeWidth-1,         0, maxY), false);
        block.specifyStatisticsStatus(Box2D(maxX-envelopeWidth+1, maxX, 0, maxY), false);
    }
}

template<typename T, template<typename U> class Descriptor>
std::vector<BlockLattice2D<T,Descriptor>*>& MultiBlockLattice2D<T,Descriptor>::getBlockLattices() {
    return blockLattices;
}

template<typename T, template<typename U> class Descriptor>
std::vector<BlockLattice2D<T,Descriptor>*> const& MultiBlockLattice2D<T,Descriptor>::getBlockLattices() const {
    return blockLattices;
}

template<typename T, template<typename U> class Descriptor>
MultiBlockDistribution2D const& MultiBlockLattice2D<T,Descriptor>::getMultiBlockDistribution() const {
    return this->getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T, template<typename U> class Descriptor>
BlockParameters2D const& MultiBlockLattice2D<T,Descriptor>::getParameters(plint iParam) const {
    return getMultiBlockDistribution().getBlockParameters(iParam);
}

template<typename T, template<typename U> class Descriptor>
Overlap2D const& MultiBlockLattice2D<T,Descriptor>::getNormalOverlap(plint iOverlap) const {
    return getMultiBlockDistribution().getNormalOverlap(iOverlap);
}

template<typename T, template<typename U> class Descriptor>
PeriodicOverlap2D const& MultiBlockLattice2D<T,Descriptor>::getPeriodicOverlap(plint iOverlap) const {
    return getMultiBlockDistribution().getPeriodicOverlap(iOverlap);
}

template<typename T, template<typename U> class Descriptor>
AtomicBlock2D<T>& MultiBlockLattice2D<T,Descriptor>::getComponent(plint iBlock) {
    PLB_PRECONDITION( iBlock<(plint)getBlockLattices().size() );
    return *getBlockLattices()[iBlock];
}

template<typename T, template<typename U> class Descriptor>
AtomicBlock2D<T> const& MultiBlockLattice2D<T,Descriptor>::getComponent(plint iBlock) const {
    PLB_PRECONDITION( iBlock<(plint)getBlockLattices().size() );
    return *getBlockLattices()[iBlock];
}

template<typename T, template<typename U> class Descriptor>
plint MultiBlockLattice2D<T,Descriptor>::sizeOfCell() const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
identifiers::BlockId MultiBlockLattice2D<T,Descriptor>::getBlockId() const {
    return identifiers::getLatticeId<T,Descriptor>();
}

}  // namespace plb

#endif
