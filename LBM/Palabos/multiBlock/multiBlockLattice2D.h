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
 * A 2D multiblock lattice -- header file.
 */
#ifndef MULTI_BLOCK_LATTICE_2D_H
#define MULTI_BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"
#include "core/blockLatticeBase2D.h"
#include "core/blockStatistics.h"
#include "core/cell.h"
#include "core/dynamics.h"
#include <vector>


namespace plb {

template<typename T, template<typename U> class Descriptor> class BlockLattice2D;

template<typename T, template<typename U> class Descriptor>
struct MultiCellAccess2D {
    virtual ~MultiCellAccess2D() { }
    virtual Cell<T,Descriptor>& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<BlockLattice2D<T,Descriptor>*>& lattices ) =0;
    virtual Cell<T,Descriptor> const& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<BlockLattice2D<T,Descriptor>*> const& lattices ) const =0;
    virtual void broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                               MultiBlockManagement2D const& multiBlockManagement) const=0;
    virtual MultiCellAccess2D<T,Descriptor>* clone() const =0;
    static MultiCellAccess2D<T,Descriptor>* generateCellAccess (
            BlockCommunicator2D<T> const& communicator );
};

/// A complex LatticeBase, itself decomposed into smaller components.
/** This extensible class can be used for example for cache-optimized
 * lattices, irregular domains (no memory allocation in areas exterior to
 * the domain) and parallel lattices. The actual behavior of the lattice
 * is parametrizable by a multiBlockHandler instance, which is given to
 * the constructor.
 *
 * The MultiBlockLattice does not itself possess LatticeProcessors. The Lattice-
 * Processors are delegated to the respective LatticeBases.
 */
template<typename T, template<typename U> class Descriptor>
class MultiBlockLattice2D : public BlockLatticeBase2D<T,Descriptor>, public MultiBlock2D<T> {
public:
    MultiBlockLattice2D(MultiBlockManagement2D const& multiBlockManagement_,
                        BlockCommunicator2D<T>* blockCommunicator_,
                        CombinedStatistics<T>* combinedStatistics_,
                        MultiCellAccess2D<T,Descriptor>* multiCellAccess_,
                        Dynamics<T,Descriptor>* backgroundDynamics );
    MultiBlockLattice2D(plint nx, plint ny, Dynamics<T,Descriptor>* backgroundDynamics);
    ~MultiBlockLattice2D();
    MultiBlockLattice2D(MultiBlockLattice2D<T,Descriptor> const& rhs);
    MultiBlockLattice2D(MultiBlock2D<T> const& rhs);
    MultiBlockLattice2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop=true);

    virtual Cell<T,Descriptor>& get(plint iX, plint iY);
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY) const;
    virtual void specifyStatisticsStatus(Box2D domain, bool status);
    virtual void collide(Box2D domain);
    virtual void collide();
    virtual void stream(Box2D domain);
    virtual void stream();
    virtual void collideAndStream(Box2D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    virtual AtomicBlock2D<T>& getComponent(plint iBlock);
    virtual AtomicBlock2D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
public:
    MultiBlockDistribution2D const& getMultiBlockDistribution() const;
    std::vector<BlockLattice2D<T,Descriptor>*>& getBlockLattices();
    std::vector<BlockLattice2D<T,Descriptor>*> const& getBlockLattices() const;
private:
    MultiBlockLattice2D<T,Descriptor>& operator=(MultiBlockLattice2D<T,Descriptor> const& rhs);
    void allocateBlocks(Dynamics<T,Descriptor>* backgroundDynamics);
    void eliminateStatisticsInEnvelope();
    Box2D extendPeriodic(Box2D const& box, plint envelopeWidth) const;
private:
    BlockParameters2D const& getParameters(plint iParam) const;
    Overlap2D const& getNormalOverlap(plint iOverlap) const;
    PeriodicOverlap2D const& getPeriodicOverlap(plint iOverlap) const;
private:
    MultiCellAccess2D<T,Descriptor>* multiCellAccess;
    std::vector<BlockLattice2D<T,Descriptor>*> blockLattices;
};

}  // namespace plb

#endif
