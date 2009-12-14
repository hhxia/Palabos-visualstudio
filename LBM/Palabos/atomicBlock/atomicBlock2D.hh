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
 * Atomic block -- generic implementation.
 */
#ifndef ATOMIC_BLOCK_2D_HH
#define ATOMIC_BLOCK_2D_HH

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/atomicBlockOperations2D.h"
#include "atomicBlock/atomicBlockSerializer2D.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

////////////////////// Class StatSubscriber2D /////////////////////////

template<typename T>
StatSubscriber2D<T>::StatSubscriber2D(AtomicBlock2D<T>& block_)
    : block(block_)
{ }

template<typename T>
plint StatSubscriber2D<T>::subscribeAverage() {
    return block.getInternalStatistics().subscribeAverage();
}

template<typename T>
plint StatSubscriber2D<T>::subscribeSum() {
    return block.getInternalStatistics().subscribeSum();
}

template<typename T>
plint StatSubscriber2D<T>::subscribeMax() {
    return block.getInternalStatistics().subscribeMax();
}

template<typename T>
plint StatSubscriber2D<T>::subscribeIntSum() {
    return block.getInternalStatistics().subscribeIntSum();
}


////////////////////// Class AtomicBlock2D /////////////////////////

template<typename T>
AtomicBlock2D<T>::AtomicBlock2D()
    : statisticsSubscriber(*this)
{ }

template<typename T>
AtomicBlock2D<T>::AtomicBlock2D(AtomicBlock2D<T> const& rhs)
    : Block2D<T>(rhs),
      statisticsSubscriber(*this),
      location(rhs.location)
{
    copyDataProcessors(rhs.explicitInternalProcessors, explicitInternalProcessors);
    copyDataProcessors(rhs.automaticInternalProcessors, automaticInternalProcessors);
}

template<typename T>
AtomicBlock2D<T>::~AtomicBlock2D()
{
    clearDataProcessors();
}

template<typename T>
void AtomicBlock2D<T>::swap(AtomicBlock2D<T>& rhs) {
    Block2D<T>::swap(rhs);
    explicitInternalProcessors.swap(rhs.explicitInternalProcessors);
    automaticInternalProcessors.swap(rhs.automaticInternalProcessors);
    std::swap(location, rhs.location);
}

template<typename T>
void AtomicBlock2D<T>::integrateDataProcessor (
    DataProcessor2D<T>* processor, plint level )
{
    // Negative level numbers account for explicit internal BlockProcessors
    if (level<0) {
        integrateDataProcessor(processor, -level-1, explicitInternalProcessors);
    }
    // Positive-or-zero level numbers account for automatic internal BlockProcessors
    else {
        integrateDataProcessor(processor, level, automaticInternalProcessors);
    }
}

template<typename T>
void AtomicBlock2D<T>::integrateDataProcessor (
    DataProcessor2D<T>* processor, plint level, DataProcessorVector& processors )
{
    if (level >= (plint)processors.size()) {
        processors.resize(level+1);
    }
    processors[level].push_back(processor);
}

template<typename T>
void AtomicBlock2D<T>::executeDataProcessor(DataProcessorGenerator2D<T> const& generator) {
    std::vector<AtomicBlock2D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void AtomicBlock2D<T>::executeDataProcessor(ReductiveDataProcessorGenerator2D<T>& generator) {
    std::vector<AtomicBlock2D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void AtomicBlock2D<T>::addInternalProcessor(DataProcessorGenerator2D<T> const& generator, plint level) {
    std::vector<AtomicBlock2D<T>*> objects;
    objects.push_back(this);
    plb::addInternalProcessor(generator, objects, level);
}

template<typename T>
void AtomicBlock2D<T>::copyDataProcessors(DataProcessorVector const& from, DataProcessorVector& to) {
    clearDataProcessors(to);
    to.resize(from.size());
    for (pluint iLevel=0; iLevel<from.size(); ++iLevel) {
        to[iLevel].resize(from[iLevel].size());
        for (pluint iProc=0; iProc<from[iLevel].size(); ++iProc) {
            to[iLevel][iProc] = from[iLevel][iProc]->clone();
        }
    }
}

template<typename T>
void AtomicBlock2D<T>::clearDataProcessors() {
    clearDataProcessors(explicitInternalProcessors);
    clearDataProcessors(automaticInternalProcessors);
}

template<typename T>
void AtomicBlock2D<T>::clearDataProcessors(DataProcessorVector& processors) {
    for (pluint iLevel=0; iLevel<processors.size(); ++iLevel) {
        for (pluint iProc=0; iProc<processors[iLevel].size(); ++iProc) {
            delete processors[iLevel][iProc];
        }
    }
    processors.clear();
}

template<typename T>
void AtomicBlock2D<T>::executeInternalProcessors() {
    for (pluint iLevel=0; iLevel<automaticInternalProcessors.size(); ++iLevel) {
        executeInternalProcessors(iLevel, automaticInternalProcessors);
    }
}

template<typename T>
void AtomicBlock2D<T>::executeInternalProcessors(plint level)
{
    // Negative level numbers account for explicit internal BlockProcessors
    if (level<0) {
        executeInternalProcessors(-level-1, explicitInternalProcessors);
    }
    // Positive-or-zero level numbers account for automatic internal BlockProcessors
    else {
        executeInternalProcessors(level, automaticInternalProcessors);
    }
}

template<typename T>
void AtomicBlock2D<T>::executeInternalProcessors(plint level, DataProcessorVector& processors)
{
    if (level<(plint)processors.size()) {
        for (pluint iProc=0; iProc<processors[level].size(); ++iProc) {
            processors[level][iProc] -> process();
        }
    }
}

template<typename T>
DataSerializer<T>* AtomicBlock2D<T>::getBlockSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering ) const
{
    return new AtomicBlockSerializer2D<T>(*this, domain, ordering);
}

template<typename T>
DataUnSerializer<T>* AtomicBlock2D<T>::getBlockUnSerializer (
            Box2D const& domain, IndexOrdering::OrderingT ordering )
{
    return new AtomicBlockUnSerializer2D<T>(*this, domain, ordering);
}

template<typename T>
void AtomicBlock2D<T>::evaluateStatistics() {
    // Copy running statistics to public statistics.
    this->getInternalStatistics().evaluate();
}

template<typename T>
StatSubscriber<T>& AtomicBlock2D<T>::internalStatSubscription() {
    return statisticsSubscriber;
}

template<typename T>
void AtomicBlock2D<T>::setLocation(Dot2D const& location_) {
    location = location_;
}

template<typename T>
Dot2D AtomicBlock2D<T>::getLocation() const {
    return location;
}

template<typename T, typename U>
Dot2D computeRelativeDisplacement(AtomicBlock2D<T> const& block1, AtomicBlock2D<U> const& block2) {
    return Dot2D(block1.getLocation().x-block2.getLocation().x,
                 block1.getLocation().y-block2.getLocation().y);
}


}  // namespace plb

#endif
