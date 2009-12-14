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
 * Serializer and UnSerializer for atomic blocks -- header file.
 */
#ifndef ATOMIC_BLOCK_SERIALIZER_3D_H
#define ATOMIC_BLOCK_SERIALIZER_3D_H

#include "core/globalDefs.h"
#include "core/util.h"
#include "atomicBlock/atomicBlock3D.h"
#include "core/serializer.h"

namespace plb {

template<typename T>
class AtomicBlockSerializer3D : public DataSerializer<T> {
public:
    AtomicBlockSerializer3D(AtomicBlock3D<T> const& block_,
                            IndexOrdering::OrderingT ordering_);
    AtomicBlockSerializer3D(AtomicBlock3D<T> const& block_,
                            Box3D domain_,
                            IndexOrdering::OrderingT ordering_);
    virtual AtomicBlockSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    AtomicBlock3D<T> const& block;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable util::Buffer<T> buffer;
    mutable plint iX, iY, iZ;
};

template<typename T>
class AtomicBlockUnSerializer3D : public DataUnSerializer<T> {
public:
    AtomicBlockUnSerializer3D(AtomicBlock3D<T>& block_,
                              IndexOrdering::OrderingT ordering_);
    AtomicBlockUnSerializer3D(AtomicBlock3D<T>& block_,
                              Box3D domain_,
                              IndexOrdering::OrderingT ordering_);
    virtual AtomicBlockUnSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual T* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    AtomicBlock3D<T>& block;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable util::Buffer<T> buffer;
    mutable plint iX, iY, iZ;
};

}  //  namespace plb

#endif
