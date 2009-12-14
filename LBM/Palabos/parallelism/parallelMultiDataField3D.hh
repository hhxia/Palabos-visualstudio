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
#ifndef PARALLEL_MULTI_DATA_FIELD_3D_HH
#define PARALLEL_MULTI_DATA_FIELD_3D_HH

#include "parallelMultiDataField3D.h"

/* *************** Class ParallelScalarAccess3D ************************ */

#ifdef PLB_MPI_PARALLEL

namespace plb {

template<typename T>
ParallelScalarAccess3D<T>::ParallelScalarAccess3D()
    : locatedBlock(0)
{ }

template<typename T>
T& ParallelScalarAccess3D<T>::getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<ScalarField3D<T>*>& fields )
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations (
            iX,iY,iZ, foundId, foundX, foundY, foundZ );
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]] -> get(foundX[0], foundY[0], foundZ[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;

}

template<typename T>
T const& ParallelScalarAccess3D<T>::getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<ScalarField3D<T>*> const& fields ) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations (
            iX,iY,iZ, foundId, foundX, foundY, foundZ );
    if (hasBulkCell) {
        distributedScalar = fields[foundId[0]] -> get(foundX[0], foundY[0], foundZ[0]);
    }
    global::mpi().bCastThroughMaster(&distributedScalar, 1, hasBulkCell);
    return distributedScalar;

}

template<typename T>
ParallelScalarAccess3D<T>* ParallelScalarAccess3D<T>::clone() const {
    return new ParallelScalarAccess3D(*this);
}


/* *************** Class ParallelTensorAccess3D ************************ */

template<typename T, int nDim>
ParallelTensorAccess3D<T,nDim>::ParallelTensorAccess3D()
    : locatedBlock(0)
{ }

template<typename T, int nDim>
Array<T,nDim>& ParallelTensorAccess3D<T,nDim>::getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<TensorField3D<T,nDim>*>& fields )
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations (
            iX,iY,iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        Array<T,nDim> const& foundTensor = fields[foundId[0]] -> get(foundX[0], foundY[0], foundZ[0]);
        for (int iD=0; iD<nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template<typename T, int nDim>
Array<T,nDim> const& ParallelTensorAccess3D<T,nDim>::getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<TensorField3D<T,nDim>*> const& fields ) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY, foundZ;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations (
            iX,iY,iZ, foundId, foundX, foundY, foundZ);
    if (hasBulkCell) {
        Array<T,nDim> const& foundTensor = fields[foundId[0]] -> get(foundX[0], foundY[0], foundZ[0]);
        for (int iD=0; iD<nDim; ++iD) {
            distributedTensor[iD] = foundTensor[iD];
        }
    }
    global::mpi().bCastThroughMaster(&distributedTensor[0], nDim, hasBulkCell);
    return distributedTensor;
}

template<typename T, int nDim>
ParallelTensorAccess3D<T,nDim>* ParallelTensorAccess3D<T,nDim>::clone() const {
    return new ParallelTensorAccess3D(*this);
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_3D_HH
