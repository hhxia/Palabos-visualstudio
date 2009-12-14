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
 * Serializer and UnSerializer for atomic blocks -- generic implementation.
 */
#ifndef ATOMIC_BLOCK_SERIALIZER_3D_HH
#define ATOMIC_BLOCK_SERIALIZER_3D_HH

#include "atomicBlock/atomicBlockSerializer3D.h"
#include "core/plbDebug.h"

namespace plb {

////////// class AtomicBlockSerializer3D ////////////////////////////

template<typename T>
AtomicBlockSerializer3D<T>::AtomicBlockSerializer3D (
        AtomicBlock3D<T> const& block_, IndexOrdering::OrderingT ordering_ )
    : block(block_), ordering(ordering_),
      domain(block.getBoundingBox()),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0)
{ }

template<typename T>
AtomicBlockSerializer3D<T>::AtomicBlockSerializer3D (
        AtomicBlock3D<T> const& block_,
        Box3D domain_,
        IndexOrdering::OrderingT ordering_ )
    : block(block_), ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0)
{ }

template<typename T>
AtomicBlockSerializer3D<T>* AtomicBlockSerializer3D<T>::clone() const {
    return new AtomicBlockSerializer3D(*this);
}

template<typename T>
pluint AtomicBlockSerializer3D<T>::getSize() const {
    return domain.nCells() * block.getDataTransfer().sizeOfCell();
}

template<typename T>
const T* AtomicBlockSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) const {
    PLB_PRECONDITION( !isEmpty() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        bufferSize = domain.getNz() * block.getDataTransfer().sizeOfCell();
        buffer.reallocate(bufferSize);
        block.getDataTransfer().send(Box3D(iX,iX,iY,iY, domain.z0, domain.z1), buffer.get());
        ++iY;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iX;
        }
    }
    else {
        bufferSize = domain.getNx() * block.getDataTransfer().sizeOfCell();
        buffer.reallocate(bufferSize);
        block.getDataTransfer().send(Box3D(domain.x0, domain.x1, iY,iY,iZ,iZ), buffer.get());
        ++iY;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iZ;
        }
    }
    return buffer.get();
}

template<typename T>
bool AtomicBlockSerializer3D<T>::isEmpty() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iZ > domain.z1;
    }
}


////////// class AtomicBlockUnSerializer3D ////////////////////////////

template<typename T>
AtomicBlockUnSerializer3D<T>::AtomicBlockUnSerializer3D (
        AtomicBlock3D<T>& block_, IndexOrdering::OrderingT ordering_ )
    : block(block_), ordering(ordering_),
      domain(block.getBoundingBox()),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0)
{ }

template<typename T>
AtomicBlockUnSerializer3D<T>::AtomicBlockUnSerializer3D (
        AtomicBlock3D<T>& block_,
        Box3D domain_,
        IndexOrdering::OrderingT ordering_ )
    : block(block_), ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0)
{ }

template<typename T>
AtomicBlockUnSerializer3D<T>* AtomicBlockUnSerializer3D<T>::clone() const {
    return new AtomicBlockUnSerializer3D(*this);
}

template<typename T>
pluint AtomicBlockUnSerializer3D<T>::getSize() const {
    return domain.nCells() * block.getDataTransfer().sizeOfCell();
}

template<typename T>
T* AtomicBlockUnSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) {
    PLB_PRECONDITION( !isFull() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        bufferSize = domain.getNz() * block.getDataTransfer().sizeOfCell();
    }
    else {
        bufferSize = domain.getNx() * block.getDataTransfer().sizeOfCell();
    }
    buffer.reallocate(bufferSize);
    return buffer.get();
}

template<typename T>
void AtomicBlockUnSerializer3D<T>::commitData() {
    PLB_PRECONDITION( !isFull() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        block.getDataTransfer().receive(Box3D(iX,iX,iY,iY, domain.z0, domain.z1), buffer.get());
        ++iY;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iX;
        }
    }
    else {
        block.getDataTransfer().receive(Box3D(domain.x0, domain.x1, iY,iY,iZ,iZ), buffer.get());
        ++iY;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iZ;
        }
    }
}

template<typename T>
bool AtomicBlockUnSerializer3D<T>::isFull() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iZ > domain.z1;
    }
}

}  //  namespace plb

#endif  // ATOMIC_BLOCK_SERIALIZER_3D_HH
