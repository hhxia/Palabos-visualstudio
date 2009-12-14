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
 * The dynamics of a 3D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_3D_H
#define BLOCK_LATTICE_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/cell.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockLatticeBase3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "core/identifiers.h"
#include <vector>

/// All OpenLB code is contained in this namespace.
namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class BlockLattice3D;


template<typename T, template<typename U> class Descriptor>
class BlockLatticeDataTransfer3D : public BlockDataTransfer3D<T> {
public:
    BlockLatticeDataTransfer3D(BlockLattice3D<T,Descriptor>& lattice_);
    virtual plint sizeOfCell() const;
    virtual void send(Box3D domain, T* buffer) const;
    virtual void receive(Box3D domain, T const* buffer);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D<T> const& from);
private:
    BlockLattice3D<T,Descriptor>& lattice;
};


/// A regular lattice for highly efficient 3D LB dynamics.
/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Descriptor>
class BlockLattice3D : public BlockLatticeBase3D<T,Descriptor>, public AtomicBlock3D<T> {
public:
    /// Construction of an nx_ by ny_ by nz_ lattice
    BlockLattice3D(plint nx_, plint ny_, plint nz_, Dynamics<T,Descriptor>* backgroundDynamics_);
    /// Destruction of the lattice
    ~BlockLattice3D();
    /// Copy construction
    BlockLattice3D(BlockLattice3D<T,Descriptor> const& rhs);
    /// Copy assignment
    BlockLattice3D& operator=(BlockLattice3D<T,Descriptor> const& rhs);
    /// Swap the content of two BlockLattices
    void swap(BlockLattice3D& rhs);
public:
    /// Get bounding box of the lattice
    virtual Box3D getBoundingBox() const;
    /// Read/write access to lattice cells
    virtual Cell<T,Descriptor>& get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX<nx);
        PLB_PRECONDITION(iY<ny);
        PLB_PRECONDITION(iZ<nz);
        return grid[iX][iY][iZ];
    }
    /// Read only access to lattice cells
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX<nx);
        PLB_PRECONDITION(iY<ny);
        PLB_PRECONDITION(iZ<nz);
        return grid[iX][iY][iZ];
    }
    /// Initialize the lattice cells to get ready for simulation
    virtual void initialize();
    /// Specify wheter statistics measurements are done on a rect. domain
    virtual void specifyStatisticsStatus (
        Box3D domain, bool status );
    /// Apply collision step to a 3D sub-box
    virtual void collide(Box3D domain);
    /// Apply collision step to the whole domain
    virtual void collide();
    /// Apply streaming step to a 3D sub-box
    virtual void stream(Box3D domain);
    /// Apply streaming step to the whole domain
    virtual void stream();
    /// Apply first collision, then streaming step to a 3D sub-box
    virtual void collideAndStream(Box3D domain);
    /// Apply first collision, then streaming step to the whole domain
    virtual void collideAndStream();
    /// Increment time counter
    virtual void incrementTime();
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual BlockLatticeDataTransfer3D<T,Descriptor>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockLatticeDataTransfer3D<T,Descriptor> const& getDataTransfer() const;
public:
    /// Attribute dynamics to a cell.
    void attributeDynamics(plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics);
    /// Get a reference to the background dynamics
    Dynamics<T,Descriptor>& getBackgroundDynamics();
    /// Get a const reference to the background dynamics
    Dynamics<T,Descriptor> const& getBackgroundDynamics() const;
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(Box3D domain);
    /// Apply streaming step to boundary cells
    void boundaryStream(Box3D bound, Box3D domain);
    /// Apply collision and streaming step to bulk (non-boundary) cells
    void bulkCollideAndStream(Box3D domain);
private:
    /// Helper method for memory allocation
    void allocateMemory();
    /// Helper method for memory de-allocation
    void releaseMemory();
    void implementPeriodicity();
private:
    void periodicDomain(Box3D domain);
private:
    plint                    nx, ny, nz;
    Dynamics<T,Descriptor>* backgroundDynamics;
    Cell<T,Descriptor>     *rawData;
    Cell<T,Descriptor>   ***grid;
    BlockLatticeDataTransfer3D<T,Descriptor> dataTransfer;
public:
    static CachePolicy3D& cachePolicy();
};

}  // namespace plb

#endif
