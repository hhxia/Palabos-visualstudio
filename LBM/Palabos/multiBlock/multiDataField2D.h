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
 * Scalar, vector and tensor fields for 2D data fields -- header file.
 */

#ifndef MULTI_DATA_FIELD_2D_H
#define MULTI_DATA_FIELD_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/dataFieldBase2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiDataFieldHandler2D.h"
#include "multiBlock/multiBlock2D.h"
#include <vector>

namespace plb {

template<typename T> class MultiScalarField2D;

template<typename T>
struct MultiScalarAccess2D {
    virtual ~MultiScalarAccess2D() { }
    virtual T& getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*>& fields ) =0;
    virtual T const& getDistributedScalar (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<ScalarField2D<T>*> const& fields ) const =0;
    virtual MultiScalarAccess2D<T>* clone() const=0;
};
    
template<typename T>
class MultiScalarField2D : public ScalarFieldBase2D<T>, public MultiBlock2D<T> {
public:
    MultiScalarField2D(MultiBlockManagement2D const& multiBlockManagement_,
                       BlockCommunicator2D<T>* blockCommunicator_,
                       CombinedStatistics<T>* combinedStatistics_,
                       MultiScalarAccess2D<T>* multiScalarAccess_);
    MultiScalarField2D(plint nx, plint ny);
    ~MultiScalarField2D();
    MultiScalarField2D(MultiScalarField2D<T> const& rhs);
    MultiScalarField2D(MultiBlock2D<T> const& rhs);
    MultiScalarField2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop=true);
    MultiScalarField2D<T>& operator=(MultiScalarField2D<T> const& rhs);
    void swap(MultiScalarField2D<T>& rhs);
public: 
    virtual void reset();
    virtual T& get(plint iX, plint iY);
    virtual T const& get(plint iX, plint iY) const;
public:
    virtual AtomicBlock2D<T>& getComponent(plint iBlock);
    virtual AtomicBlock2D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
private:
    void allocateFields();
    void deAllocateFields();
    BlockParameters2D const& getParameters(plint iParam) const;
    MultiBlockDistribution2D const& getMultiBlockDistribution() const;
private:
    std::vector<ScalarField2D<T>*> fields;
    MultiScalarAccess2D<T>* multiScalarAccess;
};


template<typename T, int nDim> class MultiTensorField2D;

template<typename T, int nDim>
struct MultiTensorAccess2D {
    virtual ~MultiTensorAccess2D() { }
    virtual Array<T,nDim>& getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*>& fields ) =0;
    virtual Array<T,nDim> const& getDistributedTensor (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<TensorField2D<T,nDim>*> const& fields ) const =0;
    virtual MultiTensorAccess2D<T,nDim>* clone() const=0;
};
 
template<typename T, int nDim>
class MultiTensorField2D : public TensorFieldBase2D<T,nDim>, public MultiBlock2D<T> {
public:
    MultiTensorField2D(MultiBlockManagement2D const& multiBlockManagement_,
                       BlockCommunicator2D<T>* blockCommunicator_,
                       CombinedStatistics<T>* combinedStatistics_,
                       MultiTensorAccess2D<T,nDim>* multiTensorAccess_);
    MultiTensorField2D(plint nx, plint ny);
    ~MultiTensorField2D();
    MultiTensorField2D(MultiTensorField2D<T,nDim> const& rhs);
    MultiTensorField2D(MultiBlock2D<T> const& rhs);
    MultiTensorField2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop=true);
    MultiTensorField2D<T,nDim>& operator=(MultiTensorField2D<T,nDim> const& rhs);
    void swap(MultiTensorField2D<T,nDim>& rhs);
public:
    virtual void reset();
    virtual Array<T,nDim>& get(plint iX, plint iY);
    virtual Array<T,nDim> const& get(plint iX, plint iY) const;
public:
    virtual AtomicBlock2D<T>& getComponent(plint iBlock);
    virtual AtomicBlock2D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
private:
    void allocateFields();
    void deAllocateFields();
    BlockParameters2D const& getParameters(plint iParam) const;
    MultiBlockDistribution2D const& getMultiBlockDistribution() const;
private:
    std::vector<TensorField2D<T,nDim>*> fields;
    MultiTensorAccess2D<T,nDim>* multiTensorAccess;
};

}  // namespace plb

#endif  // MULTI_DATA_FIELD_2D_H
