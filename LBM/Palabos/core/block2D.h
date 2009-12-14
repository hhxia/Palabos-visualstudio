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


#ifndef BLOCK_2D_H
#define BLOCK_2D_H

#include "core/globalDefs.h"
#include "core/identifiers.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/serializer.h"
#include "core/blockStatistics.h"
#include "core/periodicity2D.h"

namespace plb {

template<typename T>
class Block2D {
public:
    Block2D();
    Block2D(Block2D<T> const& rhs);
    void swap(Block2D<T>& rhs);
    virtual ~Block2D() { }
    virtual identifiers::BlockId getBlockId() const =0;
    virtual Box2D getBoundingBox() const =0;
    virtual void initialize() =0;
    virtual void executeDataProcessor(DataProcessorGenerator2D<T> const& generator) =0;
    virtual void executeDataProcessor(ReductiveDataProcessorGenerator2D<T>& generator) =0;
    virtual void addInternalProcessor(DataProcessorGenerator2D<T> const& generator, plint level=0) =0;
    virtual void executeInternalProcessors() =0;
    virtual void executeInternalProcessors(plint level) =0;
    virtual DataSerializer<T>* getBlockSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) const =0;
    virtual DataUnSerializer<T>* getBlockUnSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) =0;
    virtual StatSubscriber<T>& internalStatSubscription() =0;
    virtual void evaluateStatistics() =0;
    plint getNx() const;
    plint getNy() const;
    BlockStatistics<T>& getInternalStatistics();
    BlockStatistics<T> const& getInternalStatistics() const;
    PeriodicitySwitch2D<T> const& periodicity() const;
    PeriodicitySwitch2D<T>& periodicity();
public:
    virtual void signalPeriodicity();
private:
    BlockStatistics<T> internalStatistics;
    PeriodicitySwitch2D<T> periodicitySwitch;
};

template <typename T>
void copySerializedBlock( Block2D<T> const& from, Block2D<T>& to,
                          IndexOrdering::OrderingT ordering = IndexOrdering::forward );

/// Some end-user implementations of the Block2D have a static cache-policy class,
///   which can be access to fine-tune the performance on a given platform.
class CachePolicy2D {
public:
    CachePolicy2D(plint blockSize_) : blockSize(blockSize_)
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
