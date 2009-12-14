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


#ifndef ATOMIC_BLOCK_2D_H
#define ATOMIC_BLOCK_2D_H

#include "core/globalDefs.h"
#include "core/identifiers.h"
#include "core/block2D.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"

namespace plb {

// Forward declaration
template<typename T> class AtomicBlock2D;
    
/// Handles statistics subscriptions for the BlockLattice2D
template<typename T>
class StatSubscriber2D : public StatSubscriber<T> {
public:
    StatSubscriber2D(AtomicBlock2D<T>& block_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();
private:
    AtomicBlock2D<T>& block;
};

template<typename T>
struct BlockDataTransfer2D {
    virtual ~BlockDataTransfer2D() { }
    virtual plint sizeOfCell() const =0;
    virtual void send(Box2D domain, T* buffer) const =0;
    virtual void receive(Box2D domain, T const* buffer) =0;
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from) =0;
};

template<typename T>
class AtomicBlock2D : virtual public Block2D<T> {
public:
    AtomicBlock2D();
    ~AtomicBlock2D();
    AtomicBlock2D(AtomicBlock2D<T> const& rhs);
    void swap(AtomicBlock2D<T>& rhs);
public:
    /// Get access to data transfer between blocks
    virtual BlockDataTransfer2D<T>& getDataTransfer() =0;
    /// Get access to data transfer between blocks (const version)
    virtual BlockDataTransfer2D<T> const& getDataTransfer() const =0;
public:
    /// Data processors which act only on current block (no partners)
    virtual void executeDataProcessor(DataProcessorGenerator2D<T> const& generator);
    /// Data processors which act only on current block (no partners), and which return one or several values
    virtual void executeDataProcessor(ReductiveDataProcessorGenerator2D<T>& generator);
    /// Data processors which act only on current block (no partners)
    virtual void addInternalProcessor(DataProcessorGenerator2D<T> const& generator, plint level=0);
    /// Execute all internal dataProcessors at positive or zero level
    virtual void executeInternalProcessors();
    /// Execute all internal dataProcessors at a given level
    virtual void executeInternalProcessors(plint level);
    /// Get an object through which the atomic block can be serialized
    virtual DataSerializer<T>* getBlockSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) const;
    /// Get an object through which the atomic block can be un-serialized
    virtual DataUnSerializer<T>* getBlockUnSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering );
    /// Evaluate internal statistics
    virtual void evaluateStatistics();
    /// Get object to subscribe new internal statistics
    virtual StatSubscriber<T>& internalStatSubscription();
    void setLocation(Dot2D const& location_);
    Dot2D getLocation() const;
    /// Add a dataProcessor, which is executed after each iteration
    void integrateDataProcessor(DataProcessor2D<T>* processor, plint level);
private:
    typedef std::vector<std::vector<DataProcessor2D<T>*> > DataProcessorVector;
private:
    /// Common implementation for explicit/automatic processors
    void integrateDataProcessor (
                     DataProcessor2D<T>* processor, plint level, DataProcessorVector& processors );
    /// Common implementation for explicit/automatic processors
    void executeInternalProcessors(plint level, DataProcessorVector& processors);
    /// Copy processors from one vector to another
    void copyDataProcessors(DataProcessorVector const& from, DataProcessorVector& to);
    /// Release memory for lattice processors
    void clearDataProcessors();
    /// Release memory for a given species of lattice processors
    void clearDataProcessors(DataProcessorVector& processors);
private:
    DataProcessorVector explicitInternalProcessors;
    DataProcessorVector automaticInternalProcessors;
    StatSubscriber2D<T> statisticsSubscriber;
    plint gridLevel;
    Dot2D location;
};

template<typename T, typename U>
Dot2D computeRelativeDisplacement(AtomicBlock2D<T> const& block1, AtomicBlock2D<U> const& block2);

} // namespace plb

#endif  // ATOMIC_BLOCK_2D
