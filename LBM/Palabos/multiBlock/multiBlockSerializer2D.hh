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
#ifndef MULTI_BLOCK_SERIALIZER_2D_HH
#define MULTI_BLOCK_SERIALIZER_2D_HH

#include "parallelism/mpiManager.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"

namespace plb {

////////// class MultiBlockSerializer2D ////////////////////////////

template<typename T>
MultiBlockSerializer2D<T>::MultiBlockSerializer2D (
        MultiBlock2D<T> const& multiBlock_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(multiBlock.getBoundingBox()),
      iX(domain.x0), iY(domain.y0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockSerializer2D<T>::MultiBlockSerializer2D (
        MultiBlock2D<T> const& multiBlock_,
        Box2D domain_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_), ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockSerializer2D<T>* MultiBlockSerializer2D<T>::clone() const {
    return new MultiBlockSerializer2D<T>(*this);
}

template<typename T>
pluint MultiBlockSerializer2D<T>::getSize() const {
    if (ordering==IndexOrdering::memorySaving) {
        return
            getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
    }
    else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

template<typename T>
const T* MultiBlockSerializer2D<T>::getNextDataBuffer(pluint& bufferSize) const {
    PLB_PRECONDITION( !isEmpty() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkY(iX,iY, nextBlockId, nextChunkSize);
        if (iY+nextChunkSize>domain.y1+1) {
            nextChunkSize = domain.y1-iY+1;
        }
        if (hasData) {
            computeBufferAlongY(nextBlockId, nextChunkSize);
            bufferSize = nextChunkSize*multiBlock.sizeOfCell();
        }
        else {
            if (ordering==IndexOrdering::forward) {
                fillBufferWithZeros(nextChunkSize);
                bufferSize = nextChunkSize*multiBlock.sizeOfCell();
            }
            else {
                bufferSize = 0;
            }
        }
        iY += nextChunkSize;
        if (iY > domain.y1) {
            PLB_ASSERT(iY == domain.y1+1);
            iY = domain.y0;
            ++iX;
        }
    }
    else {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkX(iX,iY, nextBlockId, nextChunkSize);
        if (iX+nextChunkSize>domain.x1+1) {
            nextChunkSize = domain.x1-iX+1;
        }
        if (hasData) {
            computeBufferAlongX(nextBlockId, nextChunkSize);
            bufferSize = nextChunkSize*multiBlock.sizeOfCell();
        }
        else {
            fillBufferWithZeros(nextChunkSize);
            bufferSize = nextChunkSize*multiBlock.sizeOfCell();
        }
        iX += nextChunkSize;
        if (iX > domain.x1) {
            PLB_ASSERT(iX == domain.x1+1);
            iX = domain.x0;
            ++iY;
        }
    }
    return &buffer[0];
}

template<typename T>
bool MultiBlockSerializer2D<T>::isEmpty() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iY > domain.y1;
    }
}

template<typename T>
MultiBlockDistribution2D const& MultiBlockSerializer2D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool MultiBlockSerializer2D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}

template<typename T>
void MultiBlockSerializer2D<T>::computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters2D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        AtomicBlock2D<T> const& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(Box2D(localX,localX+nextChunkSize-1, localY, localY), &buffer[0]);
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
}

template<typename T>
void MultiBlockSerializer2D<T>::computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters2D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        AtomicBlock2D<T> const& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send(Box2D(localX,localX, localY, localY+nextChunkSize-1), &buffer[0]);
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
}

template<typename T>
void MultiBlockSerializer2D<T>::communicateBuffer (
        plint bufferSize, plint fromProc, bool isAllocated ) const
{
#ifdef PLB_MPI_PARALLEL
    // Exchanging a dummy message ahead of the main data has the effect
    //   of synchronizing sender and receiver and therefore avoiding network
    //   jams: given that only processor 0 is receiving and treating the data,
    //   it does a poor job handling all the requests from the other processors
    //   at the same time.
    plint dummyMessage;
    if (isAllocated && !global::mpi().isMainProcessor()) {
        global::mpi().receive(&dummyMessage, 1, 0);
        global::mpi().send(&buffer[0], bufferSize, 0);
    }
    if (!isAllocated && global::mpi().isMainProcessor()) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        global::mpi().send(&dummyMessage, 1, fromProc);
        global::mpi().receive(&buffer[0], bufferSize, fromProc);
    }
#endif
}

template<typename T>
void MultiBlockSerializer2D<T>::fillBufferWithZeros(plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
#ifdef PLB_MPI_PARALLEL
    if (global::mpi().isMainProcessor()) {
#endif
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        for (plint iBuffer=0; iBuffer<bufferSize+1; ++iBuffer) {
            buffer[iBuffer] = T();
        }
#ifdef PLB_MPI_PARALLEL
    }
#endif
}


////////// class MultiBlockUnSerializer2D ////////////////////////////

template<typename T>
MultiBlockUnSerializer2D<T>::MultiBlockUnSerializer2D (
        MultiBlock2D<T>& multiBlock_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(multiBlock.getBoundingBox()),
      iX(domain.x0), iY(domain.y0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockUnSerializer2D<T>::MultiBlockUnSerializer2D (
        MultiBlock2D<T>& multiBlock_,
        Box2D domain_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockUnSerializer2D<T>* MultiBlockUnSerializer2D<T>::clone() const {
    return new MultiBlockUnSerializer2D<T>(*this);
}


template<typename T>
pluint MultiBlockUnSerializer2D<T>::getSize() const {
    if (ordering==IndexOrdering::memorySaving) {
        return getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
    }
    else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

template<typename T>
T* MultiBlockUnSerializer2D<T>::getNextDataBuffer(pluint& bufferSize) {
    PLB_PRECONDITION( !isFull() );
    plint nextBlockId, nextChunkSize;
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        bool hasData = getMultiBlockDistribution().getNextChunkY(iX,iY, nextBlockId, nextChunkSize);
        if (iY+nextChunkSize>domain.y1+1) {
            nextChunkSize = domain.y1-iY+1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        if (ordering==IndexOrdering::memorySaving && !hasData) {
            bufferSize = 0;
        }
    }
    else {
        getMultiBlockDistribution().getNextChunkX(iX,iY, nextBlockId, nextChunkSize);
        if (iX+nextChunkSize>domain.x1+1) {
            nextChunkSize = domain.x1-iX+1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
    }
#ifdef PLB_MPI_PARALLEL
    if (global::mpi().isMainProcessor()) {
#endif
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
#ifdef PLB_MPI_PARALLEL
    }
#endif
    return &buffer[0];
}

template<typename T>
void MultiBlockUnSerializer2D<T>::commitData() {
    PLB_PRECONDITION( !isFull() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkY(iX,iY, nextBlockId, nextChunkSize);
        if (iY+nextChunkSize>domain.y1+1) {
            nextChunkSize = domain.y1-iY+1;
        }
        if (hasData) {
            fillBufferAlongY(nextBlockId, nextChunkSize);
        }
        iY += nextChunkSize;
        if (iY > domain.y1) {
            iY = domain.y0;
            ++iX;
        }
    }
    else {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkX(iX,iY, nextBlockId, nextChunkSize);
        if (iX+nextChunkSize>domain.x1+1) {
            nextChunkSize = domain.x1-iX+1;
        }
        if (hasData) {
            fillBufferAlongX(nextBlockId, nextChunkSize);
        }
        iX += nextChunkSize;
        if (iX > domain.x1) {
            iX = domain.x0;
            ++iY;
        }
    }
    // At the end of unserialization, duplicate overlaps.
    if (isFull()) {
        multiBlock.getBlockCommunicator().duplicateOverlaps(multiBlock);
    }
}

template<typename T>
bool MultiBlockUnSerializer2D<T>::isFull() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iY > domain.y1;
    }
}

template<typename T>
MultiBlockDistribution2D const& MultiBlockUnSerializer2D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool MultiBlockUnSerializer2D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}

template<typename T>
void MultiBlockUnSerializer2D<T>::fillBufferAlongX(plint nextBlockId, plint nextChunkSize) {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters2D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
    if (blockIsLocal) {
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        AtomicBlock2D<T>& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(Box2D(localX,localX+nextChunkSize-1, localY, localY), &buffer[0]);
    }
}

template<typename T>
void MultiBlockUnSerializer2D<T>::fillBufferAlongY(plint nextBlockId, plint nextChunkSize) {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters2D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
    if (blockIsLocal) {
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        AtomicBlock2D<T>& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive(Box2D(localX,localX, localY, localY+nextChunkSize-1), &buffer[0]);
    }
}

template<typename T>
void MultiBlockUnSerializer2D<T>::communicateBuffer (
        plint bufferSize, plint toProc, bool isAllocated ) const
{
#ifdef PLB_MPI_PARALLEL
    if (isAllocated && !global::mpi().isMainProcessor()) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        global::mpi().receive(&buffer[0], bufferSize, 0);
    }
    if (!isAllocated && global::mpi().isMainProcessor()) {
        global::mpi().send(&buffer[0], bufferSize, toProc);
    }
#endif
}

}  //  namespace plb

#endif  // MULTI_BLOCK_SERIALIZER_2D_HH
