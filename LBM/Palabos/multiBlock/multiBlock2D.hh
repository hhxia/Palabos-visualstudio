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
 * The 2D multiblock -- generic implementation.
 */
#ifndef MULTI_BLOCK_2D_HH
#define MULTI_BLOCK_2D_HH

#include "multiBlock/multiBlock2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include <cmath>
#include <algorithm>
using namespace std;

namespace plb {

////////////////////// Class MultiStatSubscriber2D /////////////////////////

template<typename T>
MultiStatSubscriber2D<T>::MultiStatSubscriber2D(MultiBlock2D<T>& multiBlock_)
    : multiBlock(multiBlock_)
{ }

template<typename T>
plint MultiStatSubscriber2D<T>::subscribeAverage() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeAverage();
    }
    return multiBlock.getInternalStatistics().subscribeAverage();
}

template<typename T>
plint MultiStatSubscriber2D<T>::subscribeSum() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeSum();
    }
    return multiBlock.getInternalStatistics().subscribeSum();
}

template<typename T>
plint MultiStatSubscriber2D<T>::subscribeMax() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeMax();
    }
    return multiBlock.getInternalStatistics().subscribeMax();
}

template<typename T>
plint MultiStatSubscriber2D<T>::subscribeIntSum() {
    std::vector<plint> const& relevantBlocks
        = multiBlock.getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock<relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        multiBlock.getComponent(iBlock).internalStatSubscription().subscribeIntSum();
    }
    return multiBlock.getInternalStatistics().subscribeIntSum();
}


////////////////////// Class MultiBlock2D /////////////////////////

template<typename T>
MultiBlock2D<T>::MultiBlock2D( MultiBlockManagement2D const& multiBlockManagement_,
                               BlockCommunicator2D<T>* blockCommunicator_,
                               CombinedStatistics<T>* combinedStatistics_ )
    : multiBlockManagement(multiBlockManagement_),
      blockCommunicator(blockCommunicator_),
      combinedStatistics(combinedStatistics_),
      statSubscriber(*this),
      maxProcessorLevel(-1),
      statisticsOn(true)
{ }

template<typename T>
MultiBlock2D<T>::MultiBlock2D(plint nx, plint ny)
    : multiBlockManagement( defaultMultiBlockPolicy2D().getMultiBlockManagement(nx,ny) ),
      blockCommunicator( defaultMultiBlockPolicy2D().getBlockCommunicator<T>() ),
      combinedStatistics( defaultMultiBlockPolicy2D().getCombinedStatistics<T>() ),
      statSubscriber(*this),
      maxProcessorLevel(-1),
      statisticsOn(true)
{ }

template<typename T>
MultiBlock2D<T>::MultiBlock2D(MultiBlock2D<T> const& rhs)
  : Block2D<T>(rhs),
    multiBlockManagement(rhs.multiBlockManagement),
    blockCommunicator(rhs.blockCommunicator -> clone()),
    combinedStatistics(rhs.combinedStatistics -> clone()),
    statSubscriber(*this),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn)
{ }

template<typename T>
MultiBlock2D<T>::MultiBlock2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop)
  : Block2D<T>(rhs),
    multiBlockManagement(extractMultiBlockManagement(rhs.multiBlockManagement, subDomain, crop)),
    blockCommunicator(rhs.blockCommunicator -> clone()),
    combinedStatistics(rhs.combinedStatistics -> clone()),
    statSubscriber(*this),
    maxProcessorLevel(rhs.maxProcessorLevel),
    statisticsOn(rhs.statisticsOn)
{ }

template<typename T>
void MultiBlock2D<T>::swap(MultiBlock2D<T>& rhs) {
    Block2D<T>::swap(rhs);
    multiBlockManagement.swap(rhs.multiBlockManagement);
    std::swap(blockCommunicator, rhs.blockCommunicator);
    std::swap(combinedStatistics, rhs.combinedStatistics);
    std::swap(maxProcessorLevel, rhs.maxProcessorLevel);
    std::swap(statisticsOn, rhs.statisticsOn);
}

template<typename T>
MultiBlock2D<T>::~MultiBlock2D() {
    delete blockCommunicator;
    delete combinedStatistics;
}

template<typename T>
Box2D MultiBlock2D<T>::getBoundingBox() const {
    return multiBlockManagement.getBoundingBox();
}

template<typename T>
void MultiBlock2D<T>::initialize() {
    this->executeInternalProcessors();
}

template<typename T>
MultiBlockManagement2D const& MultiBlock2D<T>::getMultiBlockManagement() const {
    return multiBlockManagement;
}

template<typename T>
BlockCommunicator2D<T> const& MultiBlock2D<T>::getBlockCommunicator() const {
    return *blockCommunicator;
}

template<typename T>
CombinedStatistics<T> const& MultiBlock2D<T>::getCombinedStatistics() const {
    return *combinedStatistics;
}

template<typename T>
std::vector<plint> const& MultiBlock2D<T>::getRelevantBlocks() const {
    return getMultiBlockManagement().getRelevantIndexes().getBlocks();
}

template<typename T>
plint MultiBlock2D<T>::getNumRelevantBlocks() const {
    return getRelevantBlocks().size();
}

template<typename T>
DataSerializer<T>* MultiBlock2D<T>::getBlockSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) const
{
    return new MultiBlockSerializer2D<T>(*this, domain, ordering);
}

template<typename T>
DataUnSerializer<T>* MultiBlock2D<T>::getBlockUnSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering )
{
    return new MultiBlockUnSerializer2D<T>(*this, domain, ordering);
}

template<typename T>
void MultiBlock2D<T>::evaluateStatistics() {
    std::vector<plint> const& relevantBlocks
        = getMultiBlockManagement().getRelevantIndexes().getBlocks();
    for (pluint rBlock=0; rBlock < relevantBlocks.size(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        getComponent(iBlock).evaluateStatistics();
    }
    if (isInternalStatisticsOn()) reduceStatistics();
}


template<typename T>
StatSubscriber<T>& MultiBlock2D<T>::internalStatSubscription() {
    return statSubscriber;
}

template<typename T>
void MultiBlock2D<T>::toggleInternalStatistics(bool statisticsOn_) {
    statisticsOn = statisticsOn_;
}

template<typename T>
bool MultiBlock2D<T>::isInternalStatisticsOn() const {
    return statisticsOn;
}

template<typename T>
void MultiBlock2D<T>::executeDataProcessor(DataProcessorGenerator2D<T> const& generator) {
    std::vector<MultiBlock2D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void MultiBlock2D<T>::executeDataProcessor (
        ReductiveDataProcessorGenerator2D<T>& generator)
{
    std::vector<MultiBlock2D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void MultiBlock2D<T>::addInternalProcessor(DataProcessorGenerator2D<T> const& generator, plint level) {
    std::vector<MultiBlock2D<T>*> objects;
    objects.push_back(this);
    plb::addInternalProcessor(generator, objects, level);
}

template<typename T>
void MultiBlock2D<T>::executeInternalProcessors() {
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
void MultiBlock2D<T>::executeInternalProcessors(plint level) {
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
void MultiBlock2D<T>::subscribeProcessor(plint level, std::vector<MultiBlock2D<T>*> modifiedBlocks,
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
void MultiBlock2D<T>::reduceStatistics() {
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
void MultiBlock2D<T>::addModifiedBlocks ( plint level,
                                          std::vector<MultiBlock2D<T>*> modifiedBlocks,
                                          std::vector<std::vector<MultiBlock2D<T>*> >& multiBlockCollection,
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
void MultiBlock2D<T>::duplicateOverlapsInModifiedMultiBlocks(plint level) {
    if (level>=0) {
        duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByAutomaticProcessors[level]);
    }
    else {
        duplicateOverlapsInModifiedMultiBlocks(multiBlocksChangedByManualProcessors[level]);
    }
}

template<typename T>
void MultiBlock2D<T>::duplicateOverlapsInModifiedMultiBlocks(std::vector<MultiBlock2D<T>*>& multiBlocks)
{
    for (pluint iBlock=0; iBlock<multiBlocks.size(); ++iBlock) {
        multiBlocks[iBlock]->getBlockCommunicator().
            duplicateOverlaps(*multiBlocks[iBlock]);
    }
}

template<typename T>
void MultiBlock2D<T>::signalPeriodicity() {
    getBlockCommunicator().signalPeriodicity();
}

}  // namespace plb

#endif
