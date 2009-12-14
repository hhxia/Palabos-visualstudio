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
#ifndef PARALLEL_STATISTICS_HH
#define PARALLEL_STATISTICS_HH

#include "parallelism/parallelStatistics.h"
#include <cmath>

namespace plb {

#ifdef PLB_MPI_PARALLEL

template<typename T>
ParallelCombinedStatistics<T>* ParallelCombinedStatistics<T>::clone() const
{
    return new ParallelCombinedStatistics<T>(*this);
}

template<typename T>
void ParallelCombinedStatistics<T>::reduceStatistics (
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights,
            std::vector<T>& sumObservables,
            std::vector<T>& maxObservables,
            std::vector<plint>& intSumObservables ) const
{
    // Averages
    for (pluint iAverage=0; iAverage<averageObservables.size(); ++iAverage) {
        T globalAverage, globalWeight;
        global::mpi().reduce(averageObservables[iAverage]*sumWeights[iAverage], globalAverage, MPI_SUM);
        global::mpi().reduce(sumWeights[iAverage], globalWeight, MPI_SUM);
        if (global::mpi().isMainProcessor() && fabs(globalWeight) > (T)0.5) {
            globalAverage /= globalWeight;
        }
        global::mpi().bCast(&globalAverage, 1);
        averageObservables[iAverage] = globalAverage;
    }

    // Sum
    for (pluint iSum=0; iSum<sumObservables.size(); ++iSum) {
        T globalSum;
        global::mpi().reduce(sumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        sumObservables[iSum] = globalSum;
    }

    // Max
    for (pluint iMax=0; iMax<maxObservables.size(); ++iMax) {
        T globalMax;
        global::mpi().reduce(maxObservables[iMax], globalMax, MPI_MAX);
        global::mpi().bCast(&globalMax, 1);
        maxObservables[iMax] = globalMax;
    }

    // Integer sum
    for (pluint iSum=0; iSum<intSumObservables.size(); ++iSum) {
        plint globalSum;
        global::mpi().reduce(intSumObservables[iSum], globalSum, MPI_SUM);
        global::mpi().bCast(&globalSum, 1);
        intSumObservables[iSum] = globalSum;
    }
}

#endif

}  // namespace plb

#endif
