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
 * Parallel access to elements of a scalar/tensor field -- header file.
 */
#ifndef PARALLEL_MULTI_DATA_FIELD_2D_H
#define PARALLEL_MULTI_DATA_FIELD_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiDataField2D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

template<typename T>
class ParallelScalarAccess2D : public MultiScalarAccess2D<T> {
public:
    ParallelScalarAccess2D();
    virtual T& getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*>& fields );
    virtual T const& getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*> const& fields ) const;
    virtual ParallelScalarAccess2D<T>* clone() const;
private:
    mutable plint locatedBlock;
    mutable T distributedScalar;
};


template<typename T, int nDim>
class ParallelTensorAccess2D : public MultiTensorAccess2D<T,nDim> {
public:
    ParallelTensorAccess2D();
    virtual Array<T,nDim>& getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*>& fields );
    virtual Array<T,nDim> const& getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*> const& fields ) const;
    virtual ParallelTensorAccess2D<T,nDim>* clone() const;
private:
    mutable plint locatedBlock;
    mutable Array<T,nDim> distributedTensor;
};

}  // namespace plb
 
#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_DATA_FIELD_2D_H
