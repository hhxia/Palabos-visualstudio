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
 * The dynamics of a 2D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_2D_H
#define BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/cell.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/identifiers.h"
#include <vector>


/// All OpenLB code is contained in this namespace.
namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class BlockLattice2D;


template<typename T, template<typename U> class Descriptor>
class BlockLatticeDataTransfer2D : public BlockDataTransfer2D<T> {
public:
    BlockLatticeDataTransfer2D(BlockLattice2D<T,Descriptor>& lattice_);
    virtual plint sizeOfCell() const;
    virtual void send(Box2D domain, T* buffer) const;
    virtual void receive(Box2D domain, T const* buffer);
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from);
private:
    BlockLattice2D<T,Descriptor>& lattice;
};


/// A regular lattice for highly efficient 2D LB dynamics.
/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Descriptor>
class BlockLattice2D : public BlockLatticeBase2D<T,Descriptor>, public AtomicBlock2D<T> {
public:
    /// Construction of an nx_ by ny_ lattice
    BlockLattice2D(plint nx_, plint ny_, Dynamics<T,Descriptor>* backgroundDynamics);
    /// Destruction of the lattice
    ~BlockLattice2D();
    /// Copy construction
    BlockLattice2D(BlockLattice2D<T,Descriptor> const& rhs);
    /// Copy assignment
    BlockLattice2D& operator=(BlockLattice2D<T,Descriptor> const& rhs);
    /// Swap the content of two BlockLattices
    void swap(BlockLattice2D& rhs);
public:
    /// Get bounding box of the lattice
    virtual Box2D getBoundingBox() const;
    /// Read/write access to lattice cells
    virtual Cell<T,Descriptor>& get(plint iX, plint iY) {
        PLB_PRECONDITION(iX<nx);
        PLB_PRECONDITION(iY<ny);
        return grid[iX][iY];
    }
    /// Read only access to lattice cells
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY) const {
        PLB_PRECONDITION(iX<nx);
        PLB_PRECONDITION(iY<ny);
        return grid[iX][iY];
    }
    /// Initialize the lattice cells to get ready for simulation
    virtual void initialize();
    /// Specify wheter statistics measurements are done on given rect. domain
    virtual void specifyStatisticsStatus(Box2D domain, bool status);
    /// Apply collision step to a rectangular domain
    virtual void collide(Box2D domain);
    /// Apply collision step to the whole domain
    virtual void collide();
    /// Apply streaming step to a rectangular domain
    virtual void stream(Box2D domain);
    /// Apply streaming step to the whole domain
    virtual void stream();
    /// Apply first collision, then streaming step to a rectangular domain
    virtual void collideAndStream(Box2D domain);
    /// Apply first collision, then streaming step to the whole domain
    virtual void collideAndStream();
    /// Increment time counter
    virtual void incrementTime();
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual BlockLatticeDataTransfer2D<T,Descriptor>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockLatticeDataTransfer2D<T,Descriptor> const& getDataTransfer() const;
public:
    /// Attribute dynamics to a cell.
    void attributeDynamics(plint iX, plint iY, Dynamics<T,Descriptor>* dynamics);
    /// Get a reference to the background dynamics
    Dynamics<T,Descriptor>& getBackgroundDynamics();
    /// Get a const reference to the background dynamics
    Dynamics<T,Descriptor> const& getBackgroundDynamics() const;
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(Box2D domain);
    /// Apply streaming step to boundary cells
    void boundaryStream(Box2D bound, Box2D domain);
    /// Apply collision and streaming step to bulk (non-boundary) cells
    void bulkCollideAndStream(Box2D domain);
private:
    /// Helper method for memory allocation
    void allocateMemory();
    /// Helper method for memory de-allocation
    void releaseMemory();
    void implementPeriodicity();
private:
    void periodicDomain(Box2D domain);
private:
    plint                    nx, ny;
    Dynamics<T,Descriptor>* backgroundDynamics;
    Cell<T,Descriptor>     *rawData;
    Cell<T,Descriptor>    **grid;
    BlockLatticeDataTransfer2D<T,Descriptor> dataTransfer;
public:
    static CachePolicy2D& cachePolicy();
};

}  // namespace plb

#endif
