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
 * A 3D multiblock lattice -- header file.
 */
#ifndef MULTI_BLOCK_LATTICE_3D_H
#define MULTI_BLOCK_LATTICE_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "core/blockLatticeBase3D.h"
#include "core/blockStatistics.h"
#include "core/cell.h"
#include "core/dynamics.h"
#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor> class BlockLattice3D;


template<typename T, template<typename U> class Descriptor>
struct MultiCellAccess3D {
    virtual ~MultiCellAccess3D() { }
    virtual Cell<T,Descriptor>& getDistributedCell (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<BlockLattice3D<T,Descriptor>*>& lattices ) =0;
    virtual Cell<T,Descriptor> const& getDistributedCell (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<BlockLattice3D<T,Descriptor>*> const& lattices ) const =0;
    virtual void broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                               MultiBlockManagement3D const& multiBlockManagement) const=0;
    virtual MultiCellAccess3D<T,Descriptor>* clone() const =0;
};

/// A complex LatticeBase, itself decomposed into smaller components.
/** This extensible class can be used for example for cache-optimized
 * lattices, irregular domains (no memory allocation in areas exterior to
 * the domain) and parallel lattices. The actual behavior of the lattice
 * is parametrizable by a multiBlockHandler instance, which is given to
 * the constructor.
 *
 * The MultiBlockLattice does not itself possess PostProcessors. The Post-
 * Processors are delegated to the respective LatticeBases.
 */
template<typename T, template<typename U> class Descriptor>
class MultiBlockLattice3D : public BlockLatticeBase3D<T,Descriptor>, public MultiBlock3D<T> {
public:
    MultiBlockLattice3D(MultiBlockManagement3D const& multiBlockManagement,
                        BlockCommunicator3D<T>* blockCommunicator_,
                        CombinedStatistics<T>* combinedStatistics_,
                        MultiCellAccess3D<T,Descriptor>* multiCellAccess_,
                        Dynamics<T,Descriptor>* backgroundDynamics);
    MultiBlockLattice3D(plint nx, plint ny, plint nz, Dynamics<T,Descriptor>* backgroundDynamics);
    ~MultiBlockLattice3D();
    MultiBlockLattice3D(MultiBlockLattice3D<T,Descriptor> const& rhs);
    MultiBlockLattice3D(MultiBlock3D<T> const& rhs);
    MultiBlockLattice3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop=true);

    virtual Cell<T,Descriptor>& get(plint iX, plint iY, plint iZ);
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY, plint iZ) const;
    virtual void specifyStatisticsStatus(Box3D domain, bool status);
    virtual void collide(Box3D domain);
    virtual void collide();
    virtual void stream(Box3D domain);
    virtual void stream();
    virtual void collideAndStream(Box3D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    virtual AtomicBlock3D<T>& getComponent(plint iBlock);
    virtual AtomicBlock3D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
public:
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    std::vector<BlockLattice3D<T,Descriptor>*>& getBlockLattices();
    std::vector<BlockLattice3D<T,Descriptor>*> const& getBlockLattices() const;
private:
    MultiBlockLattice3D<T,Descriptor>& operator=(MultiBlockLattice3D<T,Descriptor> const& rhs);
    void allocateBlocks(Dynamics<T,Descriptor>* backgroundDynamics);
    void eliminateStatisticsInEnvelope();
    Box3D extendPeriodic(Box3D const& box, plint envelopeWidth) const;
private:
    BlockParameters3D const& getParameters(plint iParam) const;
    Overlap3D const& getNormalOverlap(plint iOverlap) const;
    PeriodicOverlap3D const& getPeriodicOverlap(plint iOverlap) const;
private:
    MultiCellAccess3D<T,Descriptor>* multiCellAccess;
    std::vector<BlockLattice3D<T,Descriptor>*> blockLattices;
};

}  // namespace plb

#endif
