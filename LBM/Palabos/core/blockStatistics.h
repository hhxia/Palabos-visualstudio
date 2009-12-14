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
 * The BlockStatistics class -- header file.
 */
#ifndef BLOCK_STATISTICS_H
#define BLOCK_STATISTICS_H

#include "core/globalDefs.h"
#include <vector>
#include <algorithm>

namespace plb {

// Forward declaration

template<typename T> class BlockStatistics;

/// A counter for keeping track of discrete time evolution
class TimeCounter {
public:
    TimeCounter() : latticeTime(0) { }
    TimeCounter(plint iniTime) : latticeTime(iniTime) { }
    /// Increment the value of time step in this lattice
    void incrementTime() { ++latticeTime; };
    /// Reset the value of time step in this lattice
    void resetTime(pluint value=0) { latticeTime=value; } ;
    /// Get the value of time step in this lattice
    pluint getTime() const { return latticeTime; };
private:
    plint latticeTime;
};

/// A polymorphic type to handle subscription to a BlockStatistics class
template<typename T>
class StatSubscriber {
public:
    virtual ~StatSubscriber() { }
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage() =0;
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum() =0;
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax() =0;
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum() =0;
};

/// Store instances of observables, and compute their statistics.
/** This class is not intended to be inherited from
 */
template<typename T>
class BlockStatistics {
public:
    BlockStatistics();
    void swap(BlockStatistics<T>& rhs);
    /// evaluate() must be called after each lattice iteration.
    void evaluate();
    /// Attribute a value to the public statistics, and reset running statistics to default.
    void evaluate(std::vector<T> const& average, std::vector<T> const& sum,
                  std::vector<T> const& max, std::vector<plint> const& intSum, pluint numCells_);
    /// Contribute the values of the current cell to the statistics of an "average observable"
    void gatherAverage(plint whichAverage, T value);
    /// Contribute the values of the current cell to the statistics of a "sum observable"
    void gatherSum(plint whichSum, T value);
    /// Contribute the values of the current cell to the statistics of a "max observable"
    void gatherMax(plint whichMax, T value);
    /// Contribute the values of the current cell to the statistics of an integer "sum observable"
    void gatherIntSum(plint whichSum, plint value);
    /// Call this function once all statistics for a cell have been added
    void incrementStats();
    /// Return number of cells for which statistics have been added so far
    pluint const& getNumCells() const { return numCells; }

    /// Get the public value for any "average observable"
    T getAverage(plint whichAverage) const;
    /// Get the public value for any "sum observable"
    T getSum(plint whichSum) const;
    /// Get the public value for any "max observable"
    T getMax(plint whichMax) const;
    /// Get the public value for any integer "sum observable"
    plint getIntSum(plint whichSum) const;

    /// Get a handle to the vector with all "average observables"
    std::vector<T>& getAverageVect() { return averageVect; }
    /// Get a handle to the vector with all "sum observables"
    std::vector<T>& getSumVect() { return sumVect; }
    /// Get a handle to the vector with all "max observables"
    std::vector<T>& getMaxVect() { return maxVect; }
    /// Get a handle to the vector with all integer "sum observables"
    std::vector<plint>& getIntSumVect() { return intSumVect; }

    /// Subscribe a new observable for which the average value is computed.
    plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    plint subscribeIntSum();
private:
    /// Variables to store running statistics of type T
    std::vector<T> tmpAv, tmpSum, tmpMax;
    /// Variables to store summed integer observables
    std::vector<plint> tmpIntSum;
    /// Running value for number of cells over which statistics has been computed
    pluint tmpNumCells;
    /// Variables containing the public result of type T
    std::vector<T> averageVect, sumVect, maxVect;
    /// Variables containing the public result for the summed integer observables
    std::vector<plint> intSumVect;
    /// Public result for number of cells over which statistics has been computed
    pluint numCells;
};

}  // namespace plb

#endif
