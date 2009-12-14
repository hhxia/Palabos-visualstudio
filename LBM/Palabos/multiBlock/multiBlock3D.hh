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
 * The 3D multiblock -- generic implementation.
 */
#ifndef MULTI_BLOCK_3D_HH
#define MULTI_BLOCK_3D_HH

#include "multiBlock/multiBlock3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include <cmath>
#include <algorithm>

using namespace std;

namespace plb {

////////////////////// Class MultiStatSubscriber3D /////////////////////////

template<typename T>
MultiStatSubscriber3D<T>::MultiStatSubscriber3D(MultiBlock3D<T>& multiBlock_)
    : multiBlock(multiBlock_)
{ }

template<typename T>
plint MultiStatSubscriber3D<T>::subscribeAverage() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeAverage();
    }
    return multiBlock.getInternalStatistics().subscribeAverage();
}

template<typename T>
plint MultiStatSubscriber3D<T>::subscribeSum() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeSum();
    }
    return multiBlock.getInternalStatistics().subscribeSum();
}

template<typename T>
plint MultiStatSubscriber3D<T>::subscribeMax() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeMax();
    }
    return multiBlock.getInternalStatistics().subscribeMax();
}

template<typename T>
plint MultiStatSubscriber3D<T>::subscribeIntSum() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeIntSum();
    }
    return multiBlock.getInternalStatistics().subscribeIntSum();
}


////////////////////// Class MultiBlock3D /////////////////////////

template<typename T>
MultiBlock3D<T>::MultiBlock3D( MultiBlockManagement3D const& multiBlockManagement_,
                               BlockCommunicator3D<T>* blockCommunicator_,
                               CombinedStatistics<T>* combinedStatistics_ )
    : multiBlockManagement(multiBlockManagement_),
      blockCommunicator(blockCommunicator_),
      combinedStatistics(combinedStatistics_),
      statSubscriber(*this),
      maxProcessorLevel(-1),
      statisticsOn(true)
{ }

template<typename T>
MultiBlock3D<T>::MultiBlock3D(plint nx, plint ny, plint nz)
    : multiBlockManagement( defaultMultiBlockPolicy3D().getMultiBlockManagement(nx,ny,nz) ),
      blockCommunicator( defaultMultiBlockPolicy3D().getBlockCommunicator<T>() ),
      combinedStatistics( defaultMultiBlockPolicy3D().getCombinedStatistics<T>() ),
      statSubscriber(*this),
      maxProcessorLevel(-1),
      statisticsOn(true)
{ }

template<typename T>
MultiBlock3D<T>::MultiBlock3D(MultiBlock3D<T> const& rhs)
  : Block3D<T>(rhs),
    multiBlockManagement(rhs.multiBlockManagement),
    blockCommunicator(rhs.blockCommunicator -> clone()),
    combinedStatistics(rhs.combinedStatistics -> clone()),
    statSubscriber(*this),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn)
{ }

template<typename T>
MultiBlock3D<T>::MultiBlock3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop)
  : Block3D<T>(rhs),
    multiBlockManagement(extractMultiBlockManagement(rhs.multiBlockManagement, subDomain, crop)),
    blockCommunicator(rhs.blockCommunicator -> clone()),
    combinedStatistics(rhs.combinedStatistics -> clone()),
    statSubscriber(*this),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn)
{ }

template<typename T>
void MultiBlock3D<T>::swap(MultiBlock3D<T>& rhs) {
    Block3D<T>::swap(rhs);
    multiBlockManagement.swap(rhs.multiBlockManagement);
    std::swap(blockCommunicator, rhs.blockCommunicator);
    std::swap(combinedStatistics, rhs.combinedStatistics);
    std::swap(maxProcessorLevel, rhs.maxProcessorLevel);
    std::swap(statisticsOn, rhs.statisticsOn);
}

template<typename T>
MultiBlock3D<T>::~MultiBlock3D() {
    delete blockCommunicator;
    delete combinedStatistics;
}

template<typename T>
Box3D MultiBlock3D<T>::getBoundingBox() const {
    return multiBlockManagement.getBoundingBox();
}

template<typename T>
void MultiBlock3D<T>::initialize() {
    // Invoke duplicateOverlaps(), which fills the envelope of
    //   each sub-block. This needs to be done in the first place,
    //   because the method initialize() of each sub-block may want
    //   to access the envelope. An additional duplicateOverlaps()
    //   is invoked in the end to copy the result of initialize() to
    //   all processors
    this->executeInternalProcessors();
}

template<typename T>
MultiBlockManagement3D const& MultiBlock3D<T>::getMultiBlockManagement() const {
    return multiBlockManagement;
}

template<typename T>
BlockCommunicator3D<T> const& MultiBlock3D<T>::getBlockCommunicator() const {
    return *blockCommunicator;
}

template<typename T>
CombinedStatistics<T> const& MultiBlock3D<T>::getCombinedStatistics() const {
    return *combinedStatistics;
}

template<typename T>
std::vector<plint> const& MultiBlock3D<T>::getRelevantBlocks() const {
    return getMultiBlockManagement().getRelevantIndexes().getBlocks();
}

template<typename T>
plint MultiBlock3D<T>::getNumRelevantBlocks() const {
    return getRelevantBlocks().size();
}

template<typename T>
DataSerializer<T>* MultiBlock3D<T>::getBlockSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering ) const
{
    return new MultiBlockSerializer3D<T>(*this, domain, ordering);
}

template<typename T>
DataUnSerializer<T>* MultiBlock3D<T>::getBlockUnSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering )
{
    return new MultiBlockUnSerializer3D<T>(*this, domain, ordering);
}

template<typename T>
void MultiBlock3D<T>::evaluateStatistics() {
    std::vector<plint> const& relevantBlocks
        = getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock < relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        getComponent(iBlock).evaluateStatistics();
    }
    if (isInternalStatisticsOn()) reduceStatistics();
}


template<typename T>
StatSubscriber<T>& MultiBlock3D<T>::internalStatSubscription() {
    return statSubscriber;
}

template<typename T>
void MultiBlock3D<T>::toggleInternalStatistics(bool statisticsOn_) {
    statisticsOn = statisticsOn_;
}

template<typename T>
bool MultiBlock3D<T>::isInternalStatisticsOn() const {
    return statisticsOn;
}

template<typename T>
void MultiBlock3D<T>::executeDataProcessor(DataProcessorGenerator3D<T> const& generator) {
    std::vector<MultiBlock3D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void MultiBlock3D<T>::executeDataProcessor (
        ReductiveDataProcessorGenerator3D<T>& generator)
{
    std::vector<MultiBlock3D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void MultiBlock3D<T>::addInternalProcessor(DataProcessorGenerator3D<T> const& generator, plint level) {
    std::vector<MultiBlock3D<T>*> objects;
    objects.push_back(this);
    plb::addInternalProcessor(generator, objects, level);
}

template<typename T>
void MultiBlock3D<T>::executeInternalProcessors() {
    // Execute all automatic internal processors.
    for (plint iLevel=0; iLevel<=maxProcessorLevel; ++iLevel) {
        executeInternalProcessors(iLevel);
    }
    // Duplicate boundaries at least once in case there is no automatic processor.
    if (maxProcessorLevel==-1) {
        this->getBlockCommunicator().duplicateOverlaps(*this);
    }
}

template<typename T>
void MultiBlock3D<T>::executeInternalProcessors(plint level) {
    std::vector<plint> const& relevantBlocks
        = this -> getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock < relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        this -> getComponent(iBlock).executeInternalProcessors(level);
    }
    // At level 0, the expected behavior is to update overlaps in current MultiBlock only.
    if (level==0) {
        this->getBlockCommunicator().duplicateOverlaps(*this);
    }
    // At the other levels, all affected MultiBlocks get their overlaps updated.
    else {
        duplicateOverlapsInModifiedMultiBlocks(level);
    }
}

template<typename T>
void MultiBlock3D<T>::subscribeProcessor(plint level, std::vector<MultiBlock3D<T>*> modifiedBlocks,
                                         bool includesEnvelope)
{
    maxProcessorLevel = max(level, maxProcessorLevel);

    // Don't add any blocks if level=0, because this is the rule.
    if (level>0) {
        addModifiedBlocks(level, modifiedBlocks,
                          multiBlocksChangedByAutomaticProcessors, includesEnvelope);
    }
    else if (level<0) {
        addModifiedBlocks(-level, modifiedBlocks, 
                          multiBlocksChangedByManualProcessors, includesEnvelope);
    }
}

template<typename T>
void MultiBlock3D<T>::reduceStatistics() {
    std::vector<plint> const& relevantBlocks
        = getMultiBlockManagement().getRelevantIndexes().getBlocks();
    std::vector<BlockStatistics<T> const*> individualStatistics;
    // Prepare a vector containing the BlockStatistics of all components
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        individualStatistics.push_back(&getComponent(iBlock).getInternalStatistics());
    }
    // Execute reduction operation on all individual statistics and store result into
    //   statistics of current MultiBlock.
    combinedStatistics -> combine(individualStatistics, this->getInternalStatistics());
    // Copy result to each individual statistics
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        (getComponent(iBlock).getInternalStatistics()) = (this->getInternalStatistics());
    }
}

template<typename T>
void MultiBlock3D<T>::addModifiedBlocks ( plint level,
                                          std::vector<MultiBlock3D<T>*> modifiedBlocks,
                                          std::vector<std::vector<MultiBlock3D<T>*> >& multiBlockCollection,
                                          bool includesEnvelope )
{
    // Resize vector which collects modified blocks (resize needs to be
    //   done even when no block is added, to avoid memory violations
    //   during read access to the vector).
    if ((pluint)level >= multiBlockCollection.size()) {
        multiBlockCollection.resize(level+1);
    }
    // Unless envelope is already included in the domain of application of the data
    //   processor, subscribe modified blocks for an update of the envelope.
    if (!includesEnvelope) {
        // Add new blocks to existing ones.
        multiBlockCollection[level].insert (
                multiBlockCollection[level].end(),
                modifiedBlocks.begin(), modifiedBlocks.end() );
        // Make them unique.
        std::sort(multiBlockCollection[level].begin(), multiBlockCollection[level].end());
        multiBlockCollection[level].erase (
                std::unique( multiBlockCollection[level].begin(), multiBlockCollection[level].end() ),
                multiBlockCollection[level].end() );
    }
}

template<typename T>
void MultiBlock3D<T>::duplicateOverlapsInModifiedMultiBlocks(plint level) {
    if (level>=0) {
        duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByAutomaticProcessors[level]);
    }
    else {
        duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByManualProcessors[level]);
    }
}

template<typename T>
void MultiBlock3D<T>::duplicateOverlapsInModifiedMultiBlocks(std::vector<MultiBlock3D<T>*>& multiBlocks)
{
    for (pluint iBlock=0; iBlock<multiBlocks.size(); ++iBlock) {
        multiBlocks[iBlock]->getBlockCommunicator().
            duplicateOverlaps(*multiBlocks[iBlock]);
    }
}

template<typename T>
void MultiBlock3D<T>::signalPeriodicity() {
    getBlockCommunicator().signalPeriodicity();
}

}  // namespace plb

#endif
