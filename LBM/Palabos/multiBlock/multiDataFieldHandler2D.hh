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
 * Handler for 2D multi data field structure -- generic implementation.
 */

#ifndef MULTI_DATA_FIELD_HANDLER_2D_HH
#define MULTI_DATA_FIELD_HANDLER_2D_HH

#include "parallelism/mpiManager.h"
#include "multiBlock/multiDataFieldHandler2D.h"
#include "parallelism/parallelDynamics.h"
#include "core/util.h"
#include <algorithm>
#include <numeric>


namespace plb {

////////////////////// Class SerialMultiDataFieldHandler2D /////////////////////

template<typename T>
SerialMultiDataFieldHandler2D<T>::SerialMultiDataFieldHandler2D (
        MultiBlockDistribution2D const& multiBlockDistribution_ )
    : multiBlockDistribution(multiBlockDistribution_)
{ }

template<typename T>
plint SerialMultiDataFieldHandler2D<T>::getNx() const {
    return multiBlockDistribution.getBoundingBox().getNx();
}
template<typename T>
plint SerialMultiDataFieldHandler2D<T>::getNy() const {
    return multiBlockDistribution.getBoundingBox().getNy();
}

template<typename T>
MultiBlockDistribution2D const& SerialMultiDataFieldHandler2D<T>::getMultiBlockDistribution() const {
    return multiBlockDistribution;
}

template<typename T>
bool SerialMultiDataFieldHandler2D<T>::getLocalEnvelope(plint iBlock, plint& lx, plint& ly) const {
    BlockParameters2D const& parameters = multiBlockDistribution.getBlockParameters(iBlock);
    lx = parameters.getEnvelopeLx();
    ly = parameters.getEnvelopeLy();
    return true;
}

template<typename T>
T SerialMultiDataFieldHandler2D<T>::reduceSum (T localSum) const {
    return localSum;
}

template<typename T>
T SerialMultiDataFieldHandler2D<T>::reduceAverage (T localAverage, T localWeight) const {
    return localAverage;
}

template<typename T>
T SerialMultiDataFieldHandler2D<T>::reduceMin (T localMin) const {
    return localMin;
}

template<typename T>
T SerialMultiDataFieldHandler2D<T>::reduceMax (T localMax) const {
    return localMax;
}

template<typename T>
void SerialMultiDataFieldHandler2D<T>::broadCastScalar(T& scalar, plint fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T>
void SerialMultiDataFieldHandler2D<T>::broadCastVector(T* vect, plint size, plint fromBlock) const {
    // Nothing to do in the serial case
}


#ifdef PLB_MPI_PARALLEL

////////////////////// Class ParallelMultiDataFieldHandler2D /////////////////////

template<typename T>
ParallelMultiDataFieldHandler2D<T>::ParallelMultiDataFieldHandler2D (
        MultiBlockDistribution2D const& multiBlockDistribution_ )
    : multiBlockDistribution(multiBlockDistribution_)
{ }

template<typename T>
plint ParallelMultiDataFieldHandler2D<T>::getNx() const {
    return multiBlockDistribution.getBoundingBox().getNx();
}
template<typename T>
plint ParallelMultiDataFieldHandler2D<T>::getNy() const {
    return multiBlockDistribution.getBoundingBox().getNy();
}

template<typename T>
MultiBlockDistribution2D const& ParallelMultiDataFieldHandler2D<T>::getMultiBlockDistribution() const {
    return multiBlockDistribution;
}

template<typename T>
bool ParallelMultiDataFieldHandler2D<T>::getLocalEnvelope(plint iBlock, plint& lx, plint& ly) const {
    BlockParameters2D const& parameters = multiBlockDistribution.getBlockParameters(iBlock);
    if ( parameters.getProcId() == global::mpi().getRank() ) {
        lx = parameters.getEnvelopeLx();
        ly = parameters.getEnvelopeLy();
        return true;
    }
    else {
        lx = ly = 0;
        return false;
    }
}

template<typename T>
T ParallelMultiDataFieldHandler2D<T>::reduceSum(T localSum) const {
    T globalSum;
    global::mpi().reduce(localSum, globalSum, MPI_SUM);
    global::mpi().bCast(&globalSum, 1);
    return globalSum;
}

template<typename T>
T ParallelMultiDataFieldHandler2D<T>::reduceAverage(T localAverage, T localWeight) const {
    T sumAverage, sumWeights;
    global::mpi().reduce(localAverage*localWeight, sumAverage, MPI_SUM);
    global::mpi().reduce(localWeight, sumWeights, MPI_SUM);
    if (global::mpi().isMainProcessor() && sumWeights>1.e-12) {
        sumAverage /= sumWeights;
    }
    global::mpi().bCast(&sumAverage, 1);
    return sumAverage;
}

template<typename T>
T ParallelMultiDataFieldHandler2D<T>::reduceMin(T localMin) const {
    T globalMin;
    global::mpi().reduce(localMin, globalMin, MPI_MIN);
    global::mpi().bCast(&globalMin, 1);
    return globalMin;
}

template<typename T>
T ParallelMultiDataFieldHandler2D<T>::reduceMax(T localMax) const {
    T globalMax;
    global::mpi().reduce(localMax, globalMax, MPI_MAX);
    global::mpi().bCast(&globalMax, 1);
    return globalMax;
}

template<typename T>
void ParallelMultiDataFieldHandler2D<T>::broadCastScalar(T& scalar, plint fromBlock) const {
    plint fromProc = multiBlockDistribution.getBlockParameters(fromBlock).getProcId();
    global::mpi().bCast(&scalar, 1, fromProc);
}

template<typename T>
void ParallelMultiDataFieldHandler2D<T>::broadCastVector(T* vect, plint size, plint fromBlock) const {
    plint fromProc = multiBlockDistribution.getBlockParameters(fromBlock).getProcId();
    global::mpi().bCast(vect, size, fromProc);
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif
