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
 * Scalar, vector and tensor fields for 3D data fields -- header file.
 */

#ifndef MULTI_DATA_FIELD_3D_H
#define MULTI_DATA_FIELD_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/dataFieldBase2D.h"
#include "core/dataFieldBase3D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiDataFieldHandler3D.h"
#include "multiBlock/multiBlock3D.h"
#include <vector>

namespace plb {

template<typename T> class MultiScalarField3D;

template<typename T>
struct MultiScalarAccess3D {
    virtual ~MultiScalarAccess3D() { }
    virtual T& getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<ScalarField3D<T>*>& fields ) =0;
    virtual T const& getDistributedScalar (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<ScalarField3D<T>*> const& fields ) const =0;
    virtual MultiScalarAccess3D<T>* clone() const=0;
};
 
template<typename T>
class MultiScalarField3D : public ScalarFieldBase3D<T>, public MultiBlock3D<T> {
public:
    MultiScalarField3D(MultiBlockManagement3D const& multiBlockManagement_,
                       BlockCommunicator3D<T>* blockCommunicator_,
                       CombinedStatistics<T>* combinedStatistics_,
                       MultiScalarAccess3D<T>* multiScalarAccess_);
    MultiScalarField3D(plint nx, plint ny, plint nz);
    ~MultiScalarField3D();
    MultiScalarField3D(MultiScalarField3D<T> const& rhs);
    MultiScalarField3D(MultiBlock3D<T> const& rhs);
    MultiScalarField3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop=true);
    MultiScalarField3D<T>& operator=(MultiScalarField3D<T> const& rhs);
    void swap(MultiScalarField3D<T>& rhs);
public: 
    virtual void reset();
    virtual T& get(plint iX, plint iY, plint iZ);
    virtual T const& get(plint iX, plint iY, plint iZ) const;
public:
    virtual AtomicBlock3D<T>& getComponent(plint iBlock);
    virtual AtomicBlock3D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
private:
    void allocateFields();
    void deAllocateFields();
    BlockParameters3D const& getParameters(plint iParam) const;
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
private:
    std::vector<ScalarField3D<T>*> fields;
    MultiScalarAccess3D<T>* multiScalarAccess;
};


template<typename T, int nDim> class MultiTensorField3D;

template<typename T, int nDim>
struct MultiTensorAccess3D {
    virtual ~MultiTensorAccess3D() { }
    virtual Array<T,nDim>& getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<TensorField3D<T,nDim>*>& fields ) =0;
    virtual Array<T,nDim> const& getDistributedTensor (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::vector<TensorField3D<T,nDim>*> const& fields ) const =0;
    virtual MultiTensorAccess3D<T,nDim>* clone() const=0;
};
 

template<typename T, int nDim>
class MultiTensorField3D : public TensorFieldBase3D<T,nDim>, public MultiBlock3D<T> {
public:
    MultiTensorField3D(MultiBlockManagement3D const& multiBlockManagement_,
                       BlockCommunicator3D<T>* blockCommunicator_,
                       CombinedStatistics<T>* combinedStatistics_,
                       MultiTensorAccess3D<T,nDim>* multiTensorAccess_);
    MultiTensorField3D(plint nx, plint ny, plint nz);
    ~MultiTensorField3D();
    MultiTensorField3D(MultiTensorField3D<T,nDim> const& rhs);
    MultiTensorField3D(MultiBlock3D<T> const& rhs);
    MultiTensorField3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop=true);
    MultiTensorField3D<T,nDim>& operator=(MultiTensorField3D<T,nDim> const& rhs);
    void swap(MultiTensorField3D<T,nDim>& rhs);
public:
    virtual void reset();
    virtual Array<T,nDim>& get(plint iX, plint iY, plint iZ);
    virtual Array<T,nDim> const& get(plint iX, plint iY, plint iZ) const;
public:
    virtual AtomicBlock3D<T>& getComponent(plint iBlock);
    virtual AtomicBlock3D<T> const& getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual identifiers::BlockId getBlockId() const;
private:
    void allocateFields();
    void deAllocateFields();
    BlockParameters3D const& getParameters(plint iParam) const;
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
private:
    std::vector<TensorField3D<T,nDim>*> fields;
    MultiTensorAccess3D<T,nDim>* multiTensorAccess;
};

}  // namespace plb

#endif  // MULTI_DATA_FIELD_3D_H
