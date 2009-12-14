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
 * The 3D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_3D_H
#define MULTI_BLOCK_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/combinedStatistics.h"
#include "core/block3D.h"

namespace plb {

template<typename T> class AtomicBlock3D;
template<typename T> class MultiBlock3D;

/// Handles statistics subscriptions for the MultiBlockLattice3D
template<typename T>
class MultiStatSubscriber3D : public StatSubscriber<T> {
public:
    MultiStatSubscriber3D(MultiBlock3D<T>& multiBlock_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();
private:
    MultiBlock3D<T>& multiBlock;
};

template<typename T> class MultiScalarField3D;
template<typename T, int nDim> class MultiTensorField3D;

template<typename T>
class MultiBlock3D : virtual public Block3D<T> {
public:
    MultiBlock3D ( MultiBlockManagement3D const& multiBlockManagement_,
                   BlockCommunicator3D<T>* blockCommunicator_,
                   CombinedStatistics<T>* combinedStatistics_ );
    MultiBlock3D(MultiBlock3D<T> const& rhs);
    MultiBlock3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop);
    MultiBlock3D(plint nx, plint ny, plint nz);
    void swap(MultiBlock3D<T>& rhs);
    ~MultiBlock3D();
    virtual Box3D getBoundingBox() const;
    virtual void initialize();
    MultiBlockManagement3D const& getMultiBlockManagement() const;
    BlockCommunicator3D<T> const& getBlockCommunicator() const;
    CombinedStatistics<T> const& getCombinedStatistics() const;
    std::vector<plint> const& getRelevantBlocks() const;
    plint getNumRelevantBlocks() const;
    virtual DataSerializer<T>* getBlockSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering ) const;
    virtual DataUnSerializer<T>* getBlockUnSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering );
    virtual StatSubscriber<T>& internalStatSubscription();
    virtual void evaluateStatistics();
    virtual AtomicBlock3D<T>& getComponent(plint iBlock) =0;
    virtual AtomicBlock3D<T> const& getComponent(plint iBlock) const =0;
    virtual plint sizeOfCell() const =0;
public:
    void toggleInternalStatistics(bool statisticsOn_);
    bool isInternalStatisticsOn() const;
    BlockParameters3D const& getParameters(plint iBlock) const {
        return getMultiBlockManagement().getMultiBlockDistribution().getBlockParameters(iBlock);
    }
public:
    virtual void executeDataProcessor(DataProcessorGenerator3D<T> const& generator);
    virtual void executeDataProcessor(ReductiveDataProcessorGenerator3D<T>& generator);
    virtual void addInternalProcessor(DataProcessorGenerator3D<T> const& generator, plint level=0);
    virtual void executeInternalProcessors();
    virtual void executeInternalProcessors(plint level);
    void subscribeProcessor(plint level, std::vector<MultiBlock3D<T>*> modifiedBlocks,
                            bool includesEnvelope);
public:
    virtual void signalPeriodicity();
private:
    void reduceStatistics();
    void addModifiedBlocks(plint level, std::vector<MultiBlock3D<T>*> modifiedBlocks,
                           std::vector<std::vector<MultiBlock3D<T>*> >& multiBlockCollection,
                           bool includesEnvelope);
    void duplicateOverlapsInModifiedMultiBlocks(plint level);
    void duplicateOverlapsInModifiedMultiBlocks(std::vector<MultiBlock3D<T>*>& multiBlocks);
private:
    MultiBlockManagement3D multiBlockManagement;
    BlockCommunicator3D<T>* blockCommunicator;
    CombinedStatistics<T>* combinedStatistics;
    MultiStatSubscriber3D<T> statSubscriber;
    /// List of MultiBlocks which are modified by the manual processors and require an update of their envelope.
    std::vector<std::vector<MultiBlock3D<T>*> > multiBlocksChangedByManualProcessors;
    /// List of MultiBlocks which are modified by the automatic processors and require an update of their envelope.
    std::vector<std::vector<MultiBlock3D<T>*> > multiBlocksChangedByAutomaticProcessors;
    plint maxProcessorLevel;
    bool statisticsOn;
};

} // namespace plb

#endif
