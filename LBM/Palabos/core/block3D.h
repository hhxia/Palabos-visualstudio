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


#ifndef BLOCK_3D_H
#define BLOCK_3D_H

#include "core/globalDefs.h"
#include "core/identifiers.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/serializer.h"
#include "core/blockStatistics.h"
#include "core/periodicity3D.h"

namespace plb {

template<typename T>
class Block3D {
public:
    Block3D();
    Block3D(Block3D<T> const& rhs);
    void swap(Block3D<T>& rhs);
    virtual ~Block3D() { }
    virtual identifiers::BlockId getBlockId() const =0;
    virtual Box3D getBoundingBox() const =0;
    virtual void initialize() =0;
    virtual void executeDataProcessor(DataProcessorGenerator3D<T> const& generator) =0;
    virtual void executeDataProcessor(ReductiveDataProcessorGenerator3D<T>& generator) =0;
    virtual void addInternalProcessor(DataProcessorGenerator3D<T> const& generator, plint level=0) =0;
    virtual void executeInternalProcessors() =0;
    virtual void executeInternalProcessors(plint level) =0;
    virtual DataSerializer<T>* getBlockSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering ) const =0;
    virtual DataUnSerializer<T>* getBlockUnSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering ) =0;
    virtual StatSubscriber<T>& internalStatSubscription() =0;
    virtual void evaluateStatistics() =0;
    plint getNx() const;
    plint getNy() const;
    plint getNz() const;
    BlockStatistics<T>& getInternalStatistics();
    BlockStatistics<T> const& getInternalStatistics() const;
    PeriodicitySwitch3D<T> const& periodicity() const;
    PeriodicitySwitch3D<T>& periodicity();
public:
    virtual void signalPeriodicity();
private:
    BlockStatistics<T> internalStatistics;
    PeriodicitySwitch3D<T> periodicitySwitch;
};

template <typename T>
void copySerializedBlock( Block3D<T> const& from, Block3D<T>& to,
                          IndexOrdering::OrderingT ordering = IndexOrdering::forward );

/// Some end-user implementations of the Block3D have a static cache-policy class,
///   which can be access to fine-tune the performance on a given platform.
class CachePolicy3D {
public:
    CachePolicy3D(plint blockSize_) : blockSize(blockSize_)
    { }
    void setBlockSize(plint blockSize_) {
        blockSize = blockSize_;
    }
    plint getBlockSize() const {
        return blockSize;
    }
private:
    plint blockSize;
};

} // namespace plb

#endif
