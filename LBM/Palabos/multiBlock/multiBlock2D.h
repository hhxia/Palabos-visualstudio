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
 * The 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_2D_H
#define MULTI_BLOCK_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/combinedStatistics.h"
#include "core/block2D.h"

namespace plb {

template<typename T> class AtomicBlock2D;
template<typename T> class MultiBlock2D;

/// Handles statistics subscriptions for the MultiBlockLattice2D
template<typename T>
class MultiStatSubscriber2D : public StatSubscriber<T> {
public:
    MultiStatSubscriber2D(MultiBlock2D<T>& multiBlock_);
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum();
private:
    MultiBlock2D<T>& multiBlock;
};

template<typename T> class MultiScalarField2D;
template<typename T, int nDim> class MultiTensorField2D;

template<typename T>
class MultiBlock2D : virtual public Block2D<T> {
public:
    MultiBlock2D ( MultiBlockManagement2D const& multiBlockManagement_,
                   BlockCommunicator2D<T>* blockCommunicator_,
                   CombinedStatistics<T>* combinedStatistics_ );
    MultiBlock2D(MultiBlock2D<T> const& rhs);
    MultiBlock2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop);
    MultiBlock2D(plint nx, plint ny);
    void swap(MultiBlock2D<T>& rhs);
    ~MultiBlock2D();
    virtual Box2D getBoundingBox() const;
    virtual void initialize();
    MultiBlockManagement2D const& getMultiBlockManagement() const;
    BlockCommunicator2D<T> const& getBlockCommunicator() const;
    CombinedStatistics<T> const& getCombinedStatistics() const;
    std::vector<plint> const& getRelevantBlocks() const;
    plint getNumRelevantBlocks() const;
    virtual DataSerializer<T>* getBlockSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) const;
    virtual DataUnSerializer<T>* getBlockUnSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering );
    virtual StatSubscriber<T>& internalStatSubscription();
    virtual void evaluateStatistics();
    virtual AtomicBlock2D<T>& getComponent(plint iBlock) =0;
    virtual AtomicBlock2D<T> const& getComponent(plint iBlock) const =0;
    virtual plint sizeOfCell() const =0;
public:
    void toggleInternalStatistics(bool statisticsOn_);
    bool isInternalStatisticsOn() const;
    BlockParameters2D const& getParameters(plint iBlock) const {
        return getMultiBlockManagement().getMultiBlockDistribution().getBlockParameters(iBlock);
    }
public:
    virtual void executeDataProcessor(DataProcessorGenerator2D<T> const& generator);
    virtual void executeDataProcessor(ReductiveDataProcessorGenerator2D<T>& generator);
    virtual void addInternalProcessor(DataProcessorGenerator2D<T> const& generator, plint level=0);
    virtual void executeInternalProcessors();
    virtual void executeInternalProcessors(plint level);
    void subscribeProcessor(plint level, std::vector<MultiBlock2D<T>*> modifiedBlocks,
                            bool includesEnvelope);
public:
    virtual void signalPeriodicity();
private:
    void reduceStatistics();
    void addModifiedBlocks(plint level, std::vector<MultiBlock2D<T>*> modifiedBlocks,
                           std::vector<std::vector<MultiBlock2D<T>*> >& multiBlockCollection,
                           bool includesEnvelope);
    void duplicateOverlapsInModifiedMultiBlocks(plint level);
    void duplicateOverlapsInModifiedMultiBlocks(std::vector<MultiBlock2D<T>*>& multiBlocks);
private:
    MultiBlockManagement2D multiBlockManagement;
    BlockCommunicator2D<T>* blockCommunicator;
    CombinedStatistics<T>* combinedStatistics;
    MultiStatSubscriber2D<T> statSubscriber;
    /// List of MultiBlocks which are modified by the manual processors and require an update of their envelope.
    std::vector<std::vector<MultiBlock2D<T>*> > multiBlocksChangedByManualProcessors;
    /// List of MultiBlocks which are modified by the automatic processors and require an update of their envelope.
    std::vector<std::vector<MultiBlock2D<T>*> > multiBlocksChangedByAutomaticProcessors;
    plint maxProcessorLevel;
    bool statisticsOn;
};

} // namespace plb

#endif
