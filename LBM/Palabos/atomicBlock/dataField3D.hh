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
 * Scalar, vector and tensor fields for 3D data analysis -- generic implementation.
 */
#ifndef DATA_FIELD_3D_HH
#define DATA_FIELD_3D_HH

#include "atomicBlock/dataField3D.h"
#include "core/block3D.hh"
#include "atomicBlock/atomicBlock3D.hh"
#include <algorithm>
#include <typeinfo>

namespace plb {

/////// Class ScalarField3D //////////////////////////////////

template<typename T>
ScalarField3D<T>::ScalarField3D(plint nx_, plint ny_, plint nz_)
    : nx(nx_), ny(ny_), nz(nz_),
      dataTransfer(*this)
{
    allocateMemory();
}

template<typename T>
ScalarField3D<T>::~ScalarField3D() {
    releaseMemory();
}

template<typename T>
ScalarField3D<T>::ScalarField3D(ScalarField3D<T> const& rhs)
    : AtomicBlock3D<T>(rhs),
      nx(rhs.nx), ny(rhs.ny), nz(rhs.nz),
      dataTransfer(*this)
{
    allocateMemory();
    for (pluint iData=0; iData<getSize(); ++iData) {
        (*this)[iData] = rhs[iData];
    }
}

template<typename T>
ScalarField3D<T>& ScalarField3D<T>::operator=(ScalarField3D<T> const& rhs) {
    ScalarField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
void ScalarField3D<T>::swap(ScalarField3D<T>& rhs) {
    AtomicBlock3D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template<typename T>
Box3D ScalarField3D<T>::getBoundingBox() const {
    return Box3D(0, nx-1, 0, ny-1, 0, nz-1);
}

template<typename T>
void ScalarField3D<T>::reset() {
    for (plint index=0; index<nx*ny; ++index) {
        (*this)[index] = T();
    }
}

template<typename T>
identifiers::BlockId ScalarField3D<T>::getBlockId() const {
    return identifiers::getScalarId<T>();
}

template<typename T>
ScalarFieldDataTransfer3D<T>& ScalarField3D<T>::getDataTransfer() {
    return dataTransfer;
}

template<typename T>
ScalarFieldDataTransfer3D<T> const& ScalarField3D<T>::getDataTransfer() const {
    return dataTransfer;
}

template<typename T>
void ScalarField3D<T>::allocateMemory() {
    rawData = new T [(pluint)nx*(pluint)ny*(pluint)nz];
    field   = new T** [(pluint)nx];
    for (plint iX=0; iX<nx; ++iX) {
        field[iX] = new T* [(pluint)ny];
        for (plint iY=0; iY<ny; ++iY) {
            field[iX][iY] = rawData + (pluint)nz*((pluint)iY+(pluint)ny*(pluint)iX);
        }
    }
}

template<typename T>
void ScalarField3D<T>::releaseMemory() {
    delete [] rawData; rawData = 0;
    for (plint iX=0; iX<nx; ++iX) {
      delete [] field[iX];
    }
    delete [] field;
}

////////////////////// Class ScalarFieldDataTransfer3D /////////////////////////

template<typename T>
ScalarFieldDataTransfer3D<T>::ScalarFieldDataTransfer3D(ScalarField3D<T>& field_)
    : field(field_)
{ }

template<typename T>
plint ScalarFieldDataTransfer3D<T>::sizeOfCell() const {
    return 1;
}

template<typename T>
void ScalarFieldDataTransfer3D<T>::send(Box3D domain, T* buffer) const {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *(buffer+iData) = field.get(iX,iY,iZ);
                ++iData;
            }
        }
    }
}

template<typename T>
void ScalarFieldDataTransfer3D<T>::receive(Box3D domain, T const* buffer) {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field.get(iX,iY,iZ) = *(buffer+iData);
                ++iData;
            }
        }
    }
}

template<typename T>
void ScalarFieldDataTransfer3D<T>::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D<T> const& from )
{
    PLB_PRECONDITION (typeid(from) == typeid(ScalarField3D<T> const&));
    PLB_PRECONDITION( contained(toDomain, field.getBoundingBox()) );
    ScalarField3D<T> const& fromField = (ScalarField3D<T> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                field.get(iX,iY,iZ) = fromField.get(iX+deltaX,iY+deltaY,iZ+deltaZ);
            }
        }
    }
}



//////// Class TensorField3D //////////////////////////////////

template<typename T, int nDim>
TensorField3D<T,nDim>::TensorField3D(plint nx_, plint ny_, plint nz_)
    : nx(nx_), ny(ny_), nz(nz_),
      dataTransfer(*this)
{
    allocateMemory();
}

template<typename T, int nDim>
TensorField3D<T,nDim>::~TensorField3D() {
    releaseMemory();
}

template<typename T, int nDim>
TensorField3D<T,nDim>::TensorField3D(TensorField3D<T,nDim> const& rhs)
    : AtomicBlock3D<T>(rhs),
      nx(rhs.nx), ny(rhs.ny), nz(rhs.nz),
      dataTransfer(*this)
{
    allocateMemory();
    for (plint iData=0; iData<nx*ny*nz; ++iData) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[iData][iDim] = rhs[iData][iDim];
        }
    }
}

template<typename T, int nDim>
TensorField3D<T,nDim>& TensorField3D<T,nDim>::operator=(TensorField3D<T,nDim> const& rhs) {
    TensorField3D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::swap(TensorField3D<T,nDim>& rhs) {
    AtomicBlock3D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template<typename T, int nDim>
Box3D TensorField3D<T,nDim>::getBoundingBox() const {
    return Box3D(0, nx-1, 0, ny-1, 0, nz-1);
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::reset() {
    for (plint index=0; index<nx*ny; ++index) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template<typename T, int nDim>
identifiers::BlockId TensorField3D<T,nDim>::getBlockId() const {
    return identifiers::getTensorId<T,nDim>();
}

template<typename T, int nDim>
TensorFieldDataTransfer3D<T,nDim>& TensorField3D<T,nDim>::getDataTransfer() {
    return dataTransfer;
}

template<typename T, int nDim>
TensorFieldDataTransfer3D<T,nDim> const& TensorField3D<T,nDim>::getDataTransfer() const {
    return dataTransfer;
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::allocateMemory() {
    rawData = new Array<T,nDim>   [(pluint)nx*(pluint)ny*(pluint)nz];
    field   = new Array<T,nDim>** [(pluint)nx];
    for (plint iX=0; iX<nx; ++iX) {
        field[iX] = new Array<T,nDim>* [(pluint)ny];
        for (plint iY=0; iY<ny; ++iY) {
            field[iX][iY] = rawData + (pluint)nz*((pluint)iY+(pluint)ny*(pluint)iX);
        }
    }
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::releaseMemory() {
    delete [] rawData; rawData = 0;
    for (plint iX=0; iX<nx; ++iX) {
        delete [] field[iX];
    }
    delete [] field;
}


////////////////////// Class TensorFieldDataTransfer3D /////////////////////////

template<typename T, int nDim>
TensorFieldDataTransfer3D<T,nDim>::TensorFieldDataTransfer3D(TensorField3D<T,nDim>& field_)
    : field(field_)
{ }

template<typename T, int nDim>
plint TensorFieldDataTransfer3D<T,nDim>::sizeOfCell() const {
    return nDim;
}

template<typename T, int nDim>
void TensorFieldDataTransfer3D<T,nDim>::send(Box3D domain, T* buffer) const {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    *(buffer+iData) = field.get(iX,iY,iZ)[iDim];
                    ++iData;
                }
            }
        }
    }
}

template<typename T, int nDim>
void TensorFieldDataTransfer3D<T,nDim>::receive(Box3D domain, T const* buffer) {
    PLB_PRECONDITION( contained(domain, field.getBoundingBox()) );
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    field.get(iX,iY,iZ)[iDim] = *(buffer+iData);
                    ++iData;
                }
            }
        }
    }
}

template<typename T, int nDim>
void TensorFieldDataTransfer3D<T,nDim>::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D<T> const& from )
{
    PLB_PRECONDITION (typeid(from) == typeid(TensorField3D<T,nDim> const&));
    PLB_PRECONDITION( contained(toDomain, field.getBoundingBox()) );
    TensorField3D<T,nDim> const& fromField = (TensorField3D<T,nDim> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    field.get(iX,iY,iZ)[iDim] = fromField.get(iX+deltaX,iY+deltaY,iZ+deltaZ)[iDim];
                }
            }
        }
    }
}

}  // namespace plb

#endif  
