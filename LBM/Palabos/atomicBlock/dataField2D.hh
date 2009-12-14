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
 * Scalar, vector and tensor fields for 2D data analysis -- generic implementation.
 */

#ifndef DATA_FIELD_2D_HH
#define DATA_FIELD_2D_HH

#include "atomicBlock/dataField2D.h"
#include "core/block2D.hh"
#include "atomicBlock/atomicBlock2D.hh"
#include <algorithm>
#include <typeinfo>

namespace plb {

/////// Class ScalarField2D //////////////////////////////////

template<typename T>
ScalarField2D<T>::ScalarField2D(plint nx_, plint ny_)
    : nx(nx_), ny(ny_),
      dataTransfer(*this)
{
    allocateMemory();
}

template<typename T>
ScalarField2D<T>::~ScalarField2D() {
    releaseMemory();
}

template<typename T>
ScalarField2D<T>::ScalarField2D(ScalarField2D<T> const& rhs) 
    : AtomicBlock2D<T>(rhs),
      nx(rhs.nx), ny(rhs.ny),
      dataTransfer(*this)
{
    allocateMemory();
    for (pluint iData=0; iData<getSize(); ++iData) {
        (*this)[iData] = rhs[iData];
    }
}

template<typename T>
ScalarField2D<T>& ScalarField2D<T>::operator=(ScalarField2D<T> const& rhs) {
    ScalarField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
void ScalarField2D<T>::swap(ScalarField2D<T>& rhs) {
    AtomicBlock2D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}


template<typename T>
void ScalarField2D<T>::reset() {
    for (plint index=0; index<nx*ny; ++index) {
        (*this)[index] = T();
    }
}

template<typename T>
Box2D ScalarField2D<T>::getBoundingBox() const {
    return Box2D(0, nx-1, 0, ny-1);
}

template<typename T>
identifiers::BlockId ScalarField2D<T>::getBlockId() const {
    return identifiers::getScalarId<T>();
}

template<typename T>
ScalarFieldDataTransfer2D<T>& ScalarField2D<T>::getDataTransfer() {
    return dataTransfer;
}

template<typename T>
ScalarFieldDataTransfer2D<T> const& ScalarField2D<T>::getDataTransfer() const {
    return dataTransfer;
}

template<typename T>
void ScalarField2D<T>::allocateMemory() {
    rawData = new T[(pluint)nx*(pluint)ny];
    field   = new T* [(pluint)nx];
    for (plint iX=0; iX<nx; ++iX) {
        field[iX] = rawData + (pluint)iX*(pluint)ny;
    }
}

template<typename T>
void ScalarField2D<T>::releaseMemory() {
    delete [] rawData; rawData = 0;
    delete [] field;
}

////////////////////// Class ScalarFieldDataTransfer2D /////////////////////////

template<typename T>
ScalarFieldDataTransfer2D<T>::ScalarFieldDataTransfer2D(ScalarField2D<T>& field_)
    : field(field_)
{ }

template<typename T>
plint ScalarFieldDataTransfer2D<T>::sizeOfCell() const {
    return 1;
}

template<typename T>
void ScalarFieldDataTransfer2D<T>::send(Box2D domain, T* buffer) const {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            *(buffer+iData) = field.get(iX,iY);
            ++iData;
        }
    }
}

template<typename T>
void ScalarFieldDataTransfer2D<T>::receive(Box2D domain, T const* buffer) {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            field.get(iX,iY) = *(buffer+iData);
            ++iData;
        }
    }
}

template<typename T>
void ScalarFieldDataTransfer2D<T>::attribute (
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from )
{
    PLB_PRECONDITION (typeid(from) == typeid(ScalarField2D<T> const&));
    PLB_PRECONDITION( contained(toDomain, field.getBoundingBox()) );
    ScalarField2D<T> const& fromField = (ScalarField2D<T> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            field.get(iX,iY) = fromField.get(iX+deltaX,iY+deltaY);
        }
    }
}


//////// Class TensorField2D //////////////////////////////////

template<typename T, int nDim>
TensorField2D<T,nDim>::TensorField2D(plint nx_, plint ny_)
    : nx(nx_), ny(ny_),
      dataTransfer(*this)
{
    allocateMemory();
}

template<typename T, int nDim>
TensorField2D<T,nDim>::~TensorField2D() {
    releaseMemory();
}

template<typename T, int nDim>
TensorField2D<T,nDim>::TensorField2D(TensorField2D<T,nDim> const& rhs) 
    : AtomicBlock2D<T>(rhs),
      nx(rhs.nx),
      ny(rhs.ny),
      dataTransfer(*this)
{
    allocateMemory();
    for (plint iData=0; iData<nx*ny; ++iData) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[iData][iDim] = rhs[iData][iDim];
        }
    }
}

template<typename T, int nDim>
TensorField2D<T,nDim>& TensorField2D<T,nDim>::operator=(TensorField2D<T,nDim> const& rhs) {
    TensorField2D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::swap(TensorField2D<T,nDim>& rhs) {
    AtomicBlock2D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::reset() {
    for (plint index=0; index<nx*ny; ++index) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template<typename T, int nDim>
identifiers::BlockId TensorField2D<T,nDim>::getBlockId() const {
    return identifiers::getTensorId<T,nDim>();
}

template<typename T, int nDim>
Box2D TensorField2D<T,nDim>::getBoundingBox() const {
    return Box2D(0, nx-1, 0, ny-1);
}

template<typename T, int nDim>
TensorFieldDataTransfer2D<T,nDim>& TensorField2D<T,nDim>::getDataTransfer() {
    return dataTransfer;
}

template<typename T, int nDim>
TensorFieldDataTransfer2D<T,nDim> const& TensorField2D<T,nDim>::getDataTransfer() const {
    return dataTransfer;
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::allocateMemory() {
    rawData = new Array<T,nDim>[(pluint)nx*(pluint)ny];
    field   = new Array<T,nDim>* [(pluint)nx];
    for (plint iX=0; iX<nx; ++iX) {
        field[iX] = rawData + (pluint)iX*(pluint)ny;
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::releaseMemory() {
    delete [] rawData; rawData = 0;
    delete [] field;
}

////////////////////// Class TensorFieldDataTransfer2D /////////////////////////

template<typename T, int nDim>
TensorFieldDataTransfer2D<T,nDim>::TensorFieldDataTransfer2D(TensorField2D<T,nDim>& field_)
    : field(field_)
{ }

template<typename T, int nDim>
plint TensorFieldDataTransfer2D<T,nDim>::sizeOfCell() const {
    return nDim;
}

template<typename T, int nDim>
void TensorFieldDataTransfer2D<T,nDim>::send(Box2D domain, T* buffer) const {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                *(buffer+iData) = field.get(iX,iY)[iDim];
                ++iData;
            }
        }
    }
}

template<typename T, int nDim>
void TensorFieldDataTransfer2D<T,nDim>::receive(Box2D domain, T const* buffer) {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                field.get(iX,iY)[iDim] = *(buffer+iData);
                ++iData;
            }
        }
    }
}

template<typename T, int nDim>
void TensorFieldDataTransfer2D<T,nDim>::attribute (
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from )
{
    PLB_PRECONDITION (typeid(from) == typeid(TensorField2D<T,nDim> const&));
    PLB_PRECONDITION( contained(toDomain, field.getBoundingBox()) );
    TensorField2D<T,nDim> const& fromField = (TensorField2D<T,nDim> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                field.get(iX,iY)[iDim] = fromField.get(iX+deltaX,iY+deltaY)[iDim];
            }
        }
    }
}

}  // namespace plb

#endif  
