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
 * Serial implementation of scalar, vector and tensor fields for 2D data analysis.
 * -- header file
 */

#ifndef DATA_FIELD_2D_H
#define DATA_FIELD_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/dataFieldBase2D.h"
#include "atomicBlock/atomicBlock2D.h"

namespace plb {

template<typename T> class ScalarField2D;

template<typename T>
class ScalarFieldDataTransfer2D : public BlockDataTransfer2D<T> {
public:
    ScalarFieldDataTransfer2D(ScalarField2D<T>& field_);
    virtual plint sizeOfCell() const;
    virtual void send(Box2D domain, T* buffer) const;
    virtual void receive(Box2D domain, T const* buffer);
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from);
private:
    ScalarField2D<T>& field;
};

template<typename T>
class ScalarField2D : public ScalarFieldBase2D<T>, public AtomicBlock2D<T> {
public:
    ScalarField2D(plint nx_, plint ny_);
    ~ScalarField2D();
    ScalarField2D(ScalarField2D<T> const& rhs);
    ScalarField2D<T>& operator=(ScalarField2D<T> const& rhs);
    void swap(ScalarField2D<T>& rhs);
public:
    virtual void initialize() { }
    virtual void reset();
    virtual Box2D getBoundingBox() const;
    virtual pluint getSize() const { return (pluint)nx*(pluint)ny; }
    virtual T& get(plint iX, plint iY) {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        return field[iX][iY];
    }
    virtual T const& get(plint iX, plint iY) const {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        return field[iX][iY];
    }
    T& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny);
        return rawData[ind];
    }
    T const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny);
        return rawData[ind];
    }
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual ScalarFieldDataTransfer2D<T>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual ScalarFieldDataTransfer2D<T> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    plint nx;
    plint ny;
    T   *rawData;
    T   **field;
    ScalarFieldDataTransfer2D<T> dataTransfer;
};



template<typename T, int nDim> class TensorField2D;

template<typename T, int nDim>
class TensorFieldDataTransfer2D : public BlockDataTransfer2D<T> {
public:
    TensorFieldDataTransfer2D(TensorField2D<T,nDim>& field_);
    virtual plint sizeOfCell() const;
    virtual void send(Box2D domain, T* buffer) const;
    virtual void receive(Box2D domain, T const* buffer);
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from);
private:
    TensorField2D<T,nDim>& field;
};

template<typename T, int nDim>
class TensorField2D : public TensorFieldBase2D<T,nDim>, public AtomicBlock2D<T> {
public:
    TensorField2D(plint nx_, plint ny_);
    ~TensorField2D();
    TensorField2D(TensorField2D<T,nDim> const& rhs);
    TensorField2D<T,nDim>& operator=(TensorField2D<T,nDim> const& rhs);
    void swap(TensorField2D<T,nDim>& rhs);
public:
    virtual void initialize() { }
    virtual void reset();
    virtual Box2D getBoundingBox() const;
    virtual Array<T,nDim>& get(plint iX, plint iY) {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        return field[iX][iY];
    }
    virtual Array<T,nDim> const& get(plint iX, plint iY) const {
        PLB_PRECONDITION(iX>=0 && iX<nx);
        PLB_PRECONDITION(iY>=0 && iY<ny);
        return field[iX][iY];
    }
    Array<T,nDim>& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny);
        return rawData[ind];
    }
    Array<T,nDim> const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<nx*ny);
        return rawData[ind];
    }
    /// This ID is used to restore the full identity of a Block
    virtual identifiers::BlockId getBlockId() const;
    /// Get access to data transfer between blocks
    virtual TensorFieldDataTransfer2D<T,nDim>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual TensorFieldDataTransfer2D<T,nDim> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    plint nx;
    plint ny;
    Array<T,nDim> *rawData;
    Array<T,nDim> **field;
    TensorFieldDataTransfer2D<T,nDim> dataTransfer;
};

}


#endif
