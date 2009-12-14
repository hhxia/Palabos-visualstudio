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
 * The CombinedStatistics class -- generic implementation.
 */
#ifndef COMBINED_STATISTICS_HH
#define COMBINED_STATISTICS_HH

#include "multiBlock/combinedStatistics.h"
#include <cmath>
#include <numeric>
#include <limits>

namespace plb {

template<typename T>
CombinedStatistics<T>::~CombinedStatistics()
{ }

template<typename T>
void CombinedStatistics<T>::computeLocalAverage (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights ) const
{
    // For each "average observable"
    for (pluint iAverage=0; iAverage<averageObservables.size(); ++iAverage) {
        averageObservables[iAverage] = T();
        sumWeights[iAverage] = T();
        // Compute local average
        for (pluint iStat=0; iStat<individualStatistics.size(); ++iStat) {
            T newElement = individualStatistics[iStat]->getAverage(iAverage);
            T newWeight  = individualStatistics[iStat]->getNumCells();
            averageObservables[iAverage] += newWeight * newElement;
            sumWeights[iAverage] += newWeight;
        }
        // Avoid division by zero
        if (std::fabs(sumWeights[iAverage]) > (T)0.5) {
            averageObservables[iAverage] /= sumWeights[iAverage];
        }
    }
}

template<typename T>
void CombinedStatistics<T>::computeLocalSum (
        std::vector<BlockStatistics<T> const*> const& individualStatistics,
        std::vector<T>& sumObservables ) const
{
    // For each "sum observable"
    for (pluint iSum=0; iSum<sumObservables.size(); ++iSum) {
        sumObservables[iSum] = T();
        // Compute local sum
        for (pluint iStat=0; iStat<individualStatistics.size(); ++iStat) {
            sumObservables[iSum] += individualStatistics[iStat]->getSum(iSum);
        }
    }
}

template<typename T>
void CombinedStatistics<T>::computeLocalMax (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<T>& maxObservables ) const
{
    // For each "max observable"
    for (pluint iMax=0; iMax<maxObservables.size(); ++iMax) {
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        maxObservables[iMax] = -std::numeric_limits<T>::max();
        // Compute local max
        for (pluint iStat=0; iStat<individualStatistics.size(); ++iStat) {
            T newMax = individualStatistics[iStat]->getMax(iMax);
            if (newMax > maxObservables[iMax]) {
                maxObservables[iMax] = newMax;
            }
        }
    }
}

template<typename T>
void CombinedStatistics<T>::computeLocalIntSum (
        std::vector<BlockStatistics<T> const*> const& individualStatistics,
        std::vector<plint>& intSumObservables ) const
{
    // For each integer "sum observable"
    for (pluint iSum=0; iSum<intSumObservables.size(); ++iSum) {
        intSumObservables[iSum] = 0;
        // Compute local sum
        for (pluint iStat=0; iStat<individualStatistics.size(); ++iStat) {
            intSumObservables[iSum] += individualStatistics[iStat]->getIntSum(iSum);
        }
    }
}


template<typename T>
void CombinedStatistics<T>::combine (
            std::vector<BlockStatistics<T> const*>& individualStatistics,
            BlockStatistics<T>& result ) const
{
    // Local averages
    std::vector<T> averageObservables(result.getAverageVect().size());
    std::vector<T> sumWeights(result.getAverageVect().size());
    computeLocalAverage(individualStatistics, averageObservables, sumWeights);

    // Local sums
    std::vector<T> sumObservables(result.getSumVect().size());
    computeLocalSum(individualStatistics, sumObservables);

    // Local maxima
    std::vector<T> maxObservables(result.getMaxVect().size());
    computeLocalMax(individualStatistics, maxObservables);

    // Local integer sums
    std::vector<plint> intSumObservables(result.getIntSumVect().size());
    computeLocalIntSum(individualStatistics, intSumObservables);

    // Compute global, cross-core statistics
    this->reduceStatistics (
            averageObservables, sumWeights,
            sumObservables,
            maxObservables,
            intSumObservables );

    // Update public statistics in resulting block
    result.evaluate (
        averageObservables, sumObservables, maxObservables, intSumObservables, 0 );
}


template<typename T>
SerialCombinedStatistics<T>* SerialCombinedStatistics<T>::clone() const {
    return new SerialCombinedStatistics<T>(*this);
}

template<typename T>
void SerialCombinedStatistics<T>::reduceStatistics (
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights,
            std::vector<T>& sumObservables,
            std::vector<T>& maxObservables,
            std::vector<plint>& intSumObservables ) const
{
    // Do nothing in serial case
};


}  // namespace plb

#endif
