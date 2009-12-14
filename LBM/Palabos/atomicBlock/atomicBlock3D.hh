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
#ifndef ATOMIC_BLOCK_3D_HH
#define ATOMIC_BLOCK_3D_HH

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/atomicBlockSerializer3D.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

////////////////////// Class StatSubscriber3D /////////////////////////

template<typename T>
StatSubscriber3D<T>::StatSubscriber3D(AtomicBlock3D<T>& block_)
    : block(block_)
{ }

template<typename T>
plint StatSubscriber3D<T>::subscribeAverage() {
    return block.getInternalStatistics().subscribeAverage();
}

template<typename T>
plint StatSubscriber3D<T>::subscribeSum() {
    return block.getInternalStatistics().subscribeSum();
}

template<typename T>
plint StatSubscriber3D<T>::subscribeMax() {
    return block.getInternalStatistics().subscribeMax();
}

template<typename T>
plint StatSubscriber3D<T>::subscribeIntSum() {
    return block.getInternalStatistics().subscribeIntSum();
}


////////////////////// Class AtomicBlock3D /////////////////////////

template<typename T>
AtomicBlock3D<T>::AtomicBlock3D()
    : statisticsSubscriber(*this)
{ }

template<typename T>
AtomicBlock3D<T>::AtomicBlock3D(AtomicBlock3D<T> const& rhs)
    : Block3D<T>(rhs),
      statisticsSubscriber(*this),
      location(rhs.location)
{
    copyDataProcessors(rhs.explicitInternalProcessors, explicitInternalProcessors);
    copyDataProcessors(rhs.automaticInternalProcessors, automaticInternalProcessors);
}

template<typename T>
AtomicBlock3D<T>::~AtomicBlock3D()
{
    clearDataProcessors();
}

template<typename T>
void AtomicBlock3D<T>::swap(AtomicBlock3D<T>& rhs) {
    Block3D<T>::swap(rhs);
    explicitInternalProcessors.swap(rhs.explicitInternalProcessors);
    automaticInternalProcessors.swap(rhs.automaticInternalProcessors);
    std::swap(location, rhs.location);
}

template<typename T>
void AtomicBlock3D<T>::integrateDataProcessor (
    DataProcessor3D<T>* processor, plint level )
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
void AtomicBlock3D<T>::integrateDataProcessor (
    DataProcessor3D<T>* processor, plint level, DataProcessorVector& processors )
{
    if (level >= (plint)processors.size()) {
        processors.resize(level+1);
    }
    processors[level].push_back(processor);
}

template<typename T>
void AtomicBlock3D<T>::executeDataProcessor(DataProcessorGenerator3D<T> const& generator) {
    std::vector<AtomicBlock3D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void AtomicBlock3D<T>::executeDataProcessor(ReductiveDataProcessorGenerator3D<T>& generator) {
    std::vector<AtomicBlock3D<T>*> objects;
    objects.push_back(this);
    plb::executeDataProcessor(generator, objects);
}

template<typename T>
void AtomicBlock3D<T>::addInternalProcessor(DataProcessorGenerator3D<T> const& generator, plint level) {
    std::vector<AtomicBlock3D<T>*> objects;
    objects.push_back(this);
    plb::addInternalProcessor(generator, objects, level);
}

template<typename T>
void AtomicBlock3D<T>::copyDataProcessors(DataProcessorVector const& from, DataProcessorVector& to) {
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
void AtomicBlock3D<T>::clearDataProcessors() {
    clearDataProcessors(explicitInternalProcessors);
    clearDataProcessors(automaticInternalProcessors);
}

template<typename T>
void AtomicBlock3D<T>::clearDataProcessors(DataProcessorVector& processors) {
    for (pluint iLevel=0; iLevel<processors.size(); ++iLevel) {
        for (pluint iProc=0; iProc<processors[iLevel].size(); ++iProc) {
            delete processors[iLevel][iProc];
        }
    }
    processors.clear();
}

template<typename T>
void AtomicBlock3D<T>::executeInternalProcessors() {
    for (pluint iLevel=0; iLevel<automaticInternalProcessors.size(); ++iLevel) {
        executeInternalProcessors(iLevel, automaticInternalProcessors);
    }
}

template<typename T>
void AtomicBlock3D<T>::executeInternalProcessors(plint level)
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
void AtomicBlock3D<T>::executeInternalProcessors(plint level, DataProcessorVector& processors)
{
    if (level<(plint)processors.size()) {
        for (pluint iProc=0; iProc<processors[level].size(); ++iProc) {
            processors[level][iProc] -> process();
        }
    }
}

template<typename T>
DataSerializer<T>* AtomicBlock3D<T>::getBlockSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering ) const
{
    return new AtomicBlockSerializer3D<T>(*this, domain, ordering);
}

template<typename T>
DataUnSerializer<T>* AtomicBlock3D<T>::getBlockUnSerializer (
            Box3D const& domain, IndexOrdering::OrderingT ordering )
{
    return new AtomicBlockUnSerializer3D<T>(*this, domain, ordering);
}

template<typename T>
void AtomicBlock3D<T>::evaluateStatistics() {
    // Copy running statistics to public statistics.
    this->getInternalStatistics().evaluate();
}

template<typename T>
StatSubscriber<T>& AtomicBlock3D<T>::internalStatSubscription() {
    return statisticsSubscriber;
}

template<typename T>
void AtomicBlock3D<T>::setLocation(Dot3D const& location_) {
    location = location_;
}

template<typename T>
Dot3D AtomicBlock3D<T>::getLocation() const {
    return location;
}

template<typename T, typename U>
Dot3D computeRelativeDisplacement(AtomicBlock3D<T> const& block1, AtomicBlock3D<U> const& block2) {
    return Dot3D(block1.getLocation().x-block2.getLocation().x,
                 block1.getLocation().y-block2.getLocation().y,
                 block1.getLocation().z-block2.getLocation().z);
}

}  // namespace plb

#endif
