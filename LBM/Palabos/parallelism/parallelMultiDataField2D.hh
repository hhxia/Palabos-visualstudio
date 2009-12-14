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
 * Parallel access to elements of a scalar/tensor field -- generic implementation.
 */
#ifndef PARALLEL_MULTI_DATA_FIELD_2D_HH
#define PARALLEL_MULTI_DATA_FIELD_2D_HH

#include "parallelMultiDataField2D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

/* *************** Class ParallelScalarAccess2D ************************ */

template<typename T>
ParallelScalarAccess2D<T>::ParallelScalarAccess2D()
    : locatedBlock(0)
{ }

template<typename T>
T& ParallelScalarAccess2D<T>::getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*>& fields )
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]] -> get(foundX[0], foundY[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}

template<typename T>
T const& ParallelScalarAccess2D<T>::getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*> const& fields ) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]] -> get(foundX[0], foundY[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;
}


template<typename T>
ParallelScalarAccess2D<T>* ParallelScalarAccess2D<T>::clone() const {
    return new ParallelScalarAccess2D(*this);
}


/* *************** Class ParallelTensorAccess2D ************************ */

template<typename T, int nDim>
ParallelTensorAccess2D<T,nDim>::ParallelTensorAccess2D()
    : locatedBlock(0)
{ }

template<typename T, int nDim>
Array<T,nDim>& ParallelTensorAccess2D<T,nDim>::getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*>& fields )
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        Array<T,nDim> const& foundTensor = fields[foundId[0]] -> get(foundX[0], foundY[0]);
        for (int iD=0; iD<nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template<typename T, int nDim>
Array<T,nDim> const& ParallelTensorAccess2D<T,nDim>::getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*> const& fields ) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    if (hasBulkCell) {
        Array<T,nDim> const& foundTensor = fields[foundId[0]] -> get(foundX[0], foundY[0]);
        for (int iD=0; iD<nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;

}

template<typename T, int nDim>
ParallelTensorAccess2D<T,nDim>* ParallelTensorAccess2D<T,nDim>::clone() const {
    return new ParallelTensorAccess2D(*this);
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_2D_HH


