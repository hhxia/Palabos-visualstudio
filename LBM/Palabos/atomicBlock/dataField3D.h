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
 * Serial implementation of scalar, vector and tensor fields for 3D data analysis.
 * -- header file
 */

#ifndef DATA_FIELD_3D_H
#define DATA_FIELD_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "atomicBlock/dataField2D.h"
#include "core/dataFieldBase2D.h"
#include "core/dataFieldBase3D.h"
#include "atomicBlock/atomicBlock3D.h"

namespace plb {

template<typename T> class ScalarField3D;

template<typename T>
class ScalarFieldDataTransfer3D : public BlockDataTransfer3D<T> {
public:
    ScalarFieldDataTransfer3D(ScalarField3D<T>& field_);
    virtual plint sizeOfCell() const;
    virtual void send(Box3D domain, T* buffer) const;
    virtual void receive(Box3D domain, T const* buffer);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D<T> const& from);
private:
    ScalarField3D<T>& field;
};


template<typename T>
class ScalarField3D : public ScalarFieldBase3D<T>, public AtomicBlock3D<T> {
public:
    ScalarField3D(plint nx_, plint ny_, plint nz_);
    ~ScalarField3D();
    ScalarField3D(ScalarField3D<T> const& rhs);
    ScalarField3D<T>& operator=(ScalarField3D<T> const& rhs);
    void swap(ScalarField3D<T>& rhs);
public:
    virtual void initialize() { }
    virtual void reset();
    virtual Box3D getBoundingBox() const;
    virtual pluint getSize() const { return (pluint)nx*(pluint)ny*(pluint)nz; }
    virtual T& get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        PLB_PRECONDITION(iZ>=0 && iZ<nz);
        return field[iX][iY][iZ];
    }
    virtual T const& get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        PLB_PRECONDITION(iZ>=0 && iZ<nz);
        return field[iX][iY][iZ];
    }
    T& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        return rawData[ind];
    }
    T const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        return rawData[ind];
    }
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual ScalarFieldDataTransfer3D<T>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual ScalarFieldDataTransfer3D<T> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    plint nx;
    plint ny;
    plint nz;
    T   *rawData;
    T   ***field;
    ScalarFieldDataTransfer3D<T> dataTransfer;
};

template<typename T, int nDim> class TensorField3D;

template<typename T, int nDim>
class TensorFieldDataTransfer3D : public BlockDataTransfer3D<T> {
public:
    TensorFieldDataTransfer3D(TensorField3D<T,nDim>& field_);
    virtual plint sizeOfCell() const;
    virtual void send(Box3D domain, T* buffer) const;
    virtual void receive(Box3D domain, T const* buffer);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D<T> const& from);
private:
    TensorField3D<T,nDim>& field;
};


template<typename T, int nDim>
class TensorField3D : public TensorFieldBase3D<T,nDim>, public AtomicBlock3D<T> {
public:
    TensorField3D(plint nx_, plint ny_, plint nz_);
    ~TensorField3D();
    TensorField3D(TensorField3D<T,nDim> const& rhs);
    TensorField3D<T,nDim>& operator=(TensorField3D<T,nDim> const& rhs);
    void swap(TensorField3D<T,nDim>& rhs);
public:
    virtual void initialize() { }
    virtual void reset();
    virtual Box3D getBoundingBox() const;
    virtual Array<T,nDim>& get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        PLB_PRECONDITION(iZ>=0 && iZ<nz);
        return field[iX][iY][iZ];
    }
    virtual Array<T,nDim> const& get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        PLB_PRECONDITION(iZ>=0 && iZ<nz);
        return field[iX][iY][iZ];
    }
    virtual Array<T,nDim>& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        return rawData[ind];
    }
    virtual Array<T,nDim> const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        return rawData[ind];
    }
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual TensorFieldDataTransfer3D<T,nDim>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual TensorFieldDataTransfer3D<T,nDim> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    plint nx;
    plint ny;
    plint nz;
    Array<T,nDim> *rawData;
    Array<T,nDim> ***field;
    TensorFieldDataTransfer3D<T,nDim> dataTransfer;
};

}  // namespace plb


#endif
