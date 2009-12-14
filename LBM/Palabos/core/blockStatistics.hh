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
 * The BlockStatistics class -- generic implementation.
 */
#ifndef BLOCK_STATISTICS_HH
#define BLOCK_STATISTICS_HH

#include "core/util.h"
#include "core/plbDebug.h"
#include <cmath>
#include <numeric>
#include <limits>


namespace plb {

////////////////////// Class BlockStatistics /////////////////

template<typename T>
BlockStatistics<T>::BlockStatistics()
  : tmpNumCells(0)
{ }

template<typename T>
void BlockStatistics<T>::swap(BlockStatistics<T>& rhs)
{
    tmpAv.swap    (rhs.tmpAv);
    tmpSum.swap   (rhs.tmpSum);
    tmpMax.swap   (rhs.tmpMax);
    tmpIntSum.swap(rhs.tmpIntSum);
    std::swap(tmpNumCells, rhs.tmpNumCells);

    averageVect.swap(rhs.averageVect);
    sumVect.swap    (rhs.sumVect);
    maxVect.swap    (rhs.maxVect);
    intSumVect.swap (rhs.intSumVect);
    std::swap(numCells, rhs.numCells);
}

/** In the reset function, running statistics are copied to public statistics,
 *  and running statistics are reset to zero.
 */
template<typename T>
void BlockStatistics<T>::evaluate() {
    // First step: copy running statistics to public statistics

    // Avoid division by zero while evaluating average: if no cell has
    //   been accounted for so far, simply reset to zero.
    if (tmpNumCells == 0) {   
        for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
            averageVect[iVect] = T();
        }
    }
    else {
        for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
            averageVect[iVect] = tmpAv[iVect] / (T)tmpNumCells;
        }
    }
    for (pluint iVect=0; iVect<sumVect.size(); ++iVect) {
        sumVect[iVect]     = tmpSum[iVect];
    }
    for (pluint iVect=0; iVect<maxVect.size(); ++iVect) {
        maxVect[iVect]     = tmpMax[iVect];
    }
    for (pluint iVect=0; iVect<intSumVect.size(); ++iVect) {
        intSumVect[iVect]  = tmpIntSum[iVect];
    }
    numCells = tmpNumCells;

    // Second step: reset the running statistics, in order to be ready
    //   for next lattice iteration
    for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
        tmpAv[iVect]     = T();
    }
    for (pluint iVect=0; iVect<sumVect.size(); ++iVect) {
        tmpSum[iVect]    = T();
    }
    for (pluint iVect=0; iVect<maxVect.size(); ++iVect) {
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        tmpMax[iVect]    = -std::numeric_limits<T>::max();
    }
    for (pluint iVect=0; iVect<intSumVect.size(); ++iVect) {
        tmpIntSum[iVect] = 0;
    }

    tmpNumCells = 0;
}

template<typename T>
void BlockStatistics<T>::evaluate (
        std::vector<T> const& average, std::vector<T> const& sum,
        std::vector<T> const& max, std::vector<plint> const& intSum, pluint numCells_ )
{
    PLB_PRECONDITION ( averageVect.size() == average.size() );
    PLB_PRECONDITION ( sumVect.size()     == sum.size() );
    PLB_PRECONDITION ( maxVect.size()     == max.size() );

    for (pluint iAverage=0; iAverage<averageVect.size(); ++iAverage) {
        averageVect[iAverage] = average[iAverage];
        tmpAv[iAverage] = T();
    }
    for (pluint iSum=0; iSum<sumVect.size(); ++iSum) {
        sumVect[iSum] = sum[iSum];
        tmpSum[iSum]  = T();
    }
    for (pluint iMax=0; iMax<maxVect.size(); ++iMax) {
        maxVect[iMax] = max[iMax];
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        tmpMax[iMax]  = -std::numeric_limits<T>::max();
    }
    for (pluint iSum=0; iSum<intSumVect.size(); ++iSum) {
        intSumVect[iSum] = intSum[iSum];
        tmpIntSum[iSum]  = 0;
    }
    numCells    = numCells_;
    tmpNumCells = 0;
}

/** \return Identifier for this observable, to be used for gatherAverage() and getAverage().
 */
template<typename T>
plint BlockStatistics<T>::subscribeAverage() {
    plint newSize = tmpAv.size()+1;
    tmpAv.resize(newSize);
    averageVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpAv[newIndex] = T();
    averageVect[newIndex] = T();
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherSum() and getSum().
 */
template<typename T>
plint BlockStatistics<T>::subscribeSum() {
    plint newSize = tmpSum.size()+1;
    tmpSum.resize(newSize);
    sumVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpSum[newIndex] = T();
    sumVect[newIndex] = T();
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherMax() and getMax().
 */
template<typename T>
plint BlockStatistics<T>::subscribeMax() {
    plint newSize = tmpMax.size()+1;
    tmpMax.resize(newSize);
    maxVect.resize(newSize);
    plint newIndex = newSize-1;
    // Use -max() instead of min(), because min<float> yields a positive value close to zero.
    tmpMax[newIndex] = -std::numeric_limits<T>::max();
    // Use -max() instead of min(), because min<float> yields a positive value close to zero.
    maxVect[newIndex] = -std::numeric_limits<T>::max();
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherIntSum() and getIntSum().
 */
template<typename T>
plint BlockStatistics<T>::subscribeIntSum() {
    plint newSize = tmpIntSum.size()+1;
    tmpIntSum.resize(newSize);
    intSumVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpIntSum[newIndex] = 0;
    intSumVect[newIndex] = 0;
    return newIndex;
}

template<typename T>
void BlockStatistics<T>::gatherAverage(plint whichAverage, T value) {
    PLB_PRECONDITION( whichAverage < (plint) tmpAv.size() );
    tmpAv[whichAverage] += value;
}

template<typename T>
void BlockStatistics<T>::gatherSum(plint whichSum, T value) {
    PLB_PRECONDITION( whichSum < (plint) tmpSum.size() );
    tmpSum[whichSum] += value;
}

template<typename T>
void BlockStatistics<T>::gatherMax(plint whichMax, T value) {
    PLB_PRECONDITION( whichMax < (plint) tmpMax.size() );
    if (value > tmpMax[whichMax]) {
        tmpMax[whichMax] = value;
    }
}

template<typename T>
void BlockStatistics<T>::gatherIntSum(plint whichSum, plint value) {
    PLB_PRECONDITION( whichSum < (plint) tmpIntSum.size() );
    tmpIntSum[whichSum] += value;
}

template<typename T>
void BlockStatistics<T>::incrementStats() {
    ++tmpNumCells;
}

template<typename T>
T BlockStatistics<T>::getAverage(plint whichAverage) const {
    PLB_PRECONDITION( whichAverage < (plint) tmpAv.size() );
    return averageVect[whichAverage];
}

template<typename T>
T BlockStatistics<T>::getSum(plint whichSum) const {
    PLB_PRECONDITION( whichSum < (plint) tmpSum.size() );
    return sumVect[whichSum];
}

template<typename T>
T BlockStatistics<T>::getMax(plint whichMax) const {
    PLB_PRECONDITION( whichMax < (plint) tmpMax.size() );
    return maxVect[whichMax];
}

template<typename T>
plint BlockStatistics<T>::getIntSum(plint whichSum) const {
    PLB_PRECONDITION( whichSum < (plint) tmpIntSum.size() );
    return intSumVect[whichSum];
}


}  // namespace plb

#endif
