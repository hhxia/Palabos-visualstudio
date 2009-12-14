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
 * Serializer and UnSerializer for multi blocks -- header file.
 */
#ifndef MULTI_BLOCK_SERIALIZER_3D_H
#define MULTI_BLOCK_SERIALIZER_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "core/serializer.h"

namespace plb {

template<typename T>
class MultiBlockSerializer3D : public DataSerializer<T> {
public:
    MultiBlockSerializer3D(MultiBlock3D<T> const& multiBlock_,
                           IndexOrdering::OrderingT ordering_);
    MultiBlockSerializer3D(MultiBlock3D<T> const& multiBlock_,
                           Box3D domain_,
                           IndexOrdering::OrderingT ordering_);
    virtual MultiBlockSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
    void computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongZ(plint nextBlockId, plint nextChunkSize) const;
    void communicateBuffer(plint bufferSize, plint fromProc, bool isAllocated) const;
    void fillBufferWithZeros(plint nextChunkSize) const;
private:
    MultiBlock3D<T> const& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint iX, iY, iZ;
    mutable std::vector<T> buffer;
};

template<typename T>
class MultiBlockUnSerializer3D : public DataUnSerializer<T> {
public:
    MultiBlockUnSerializer3D(MultiBlock3D<T>& multiBlock_,
                             IndexOrdering::OrderingT ordering_);
    MultiBlockUnSerializer3D(MultiBlock3D<T>& multiBlock_,
                             Box3D domain_,
                             IndexOrdering::OrderingT ordering_);
    virtual MultiBlockUnSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual T* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
    void fillBufferAlongX(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongY(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongZ(plint nextBlockId, plint nextChunkSize);
    void communicateBuffer(plint bufferSize, plint toBlock, bool isAllocated) const;
private:
    MultiBlock3D<T>& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box3D domain;
    mutable plint iX, iY, iZ;
    mutable std::vector<T> buffer;
};

template<typename T>
class RawMultiBlockSerializer3D : public DataSerializer<T> {
public:
    RawMultiBlockSerializer3D(MultiBlock3D<T> const& multiBlock_);
    virtual RawMultiBlockSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    void communicateBuffer(plint bufferSize, plint fromBlock, bool isAllocated) const;
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
private:
    MultiBlock3D<T> const& multiBlock;
    mutable plint currentBlockId;
    mutable plint iX;
    mutable std::vector<T> buffer;
};

template<typename T>
class RawMultiBlockUnSerializer3D : public DataUnSerializer<T> {
public:
    RawMultiBlockUnSerializer3D(MultiBlock3D<T>& multiBlock_);
    virtual RawMultiBlockUnSerializer3D<T>* clone() const;
    virtual pluint getSize() const;
    virtual T* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    void communicateBuffer(plint bufferSize, plint toBlock, bool isAllocated) const;
    MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
private:
    MultiBlock3D<T>& multiBlock;
    mutable plint currentBlockId;
    mutable plint iX;
    mutable std::vector<T> buffer;
};

}  //  namespace plb

#endif  // MULTI_BLOCK_SERIALIZER_3D_H
