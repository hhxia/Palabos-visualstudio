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
#ifndef MULTI_BLOCK_SERIALIZER_3D_HH
#define MULTI_BLOCK_SERIALIZER_3D_HH

#include "parallelism/mpiManager.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "core/plbDebug.h"

namespace plb {

////////// class MultiBlockSerializer3D ////////////////////////////

template<typename T>
MultiBlockSerializer3D<T>::MultiBlockSerializer3D (
        MultiBlock3D<T> const& multiBlock_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(multiBlock.getBoundingBox()),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockSerializer3D<T>::MultiBlockSerializer3D (
        MultiBlock3D<T> const& multiBlock_,
        Box3D domain_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_), ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockSerializer3D<T>* MultiBlockSerializer3D<T>::clone() const {
    return new MultiBlockSerializer3D<T>(*this);
}

template<typename T>
pluint MultiBlockSerializer3D<T>::getSize() const {
    if (ordering==IndexOrdering::memorySaving) {
        return
            getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
    }
    else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

template<typename T>
const T* MultiBlockSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) const {
    PLB_PRECONDITION( !isEmpty() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkZ(iX,iY,iZ, nextBlockId, nextChunkSize);
        if (iZ+nextChunkSize>domain.z1+1) {
            nextChunkSize = domain.z1-iZ+1;
        }
        if (hasData) {
            computeBufferAlongZ(nextBlockId, nextChunkSize);
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
        iZ += nextChunkSize;
        if (iZ > domain.z1) {
            PLB_ASSERT(iZ==domain.z1+1);
            iZ = domain.z0;
            ++iY;
            if (iY > domain.y1) {
                PLB_ASSERT(iY==domain.y1+1);
                iY = domain.y0;
                ++iX;
            }
        }
    }
    else {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkX(iX,iY,iZ, nextBlockId, nextChunkSize);
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
            PLB_ASSERT(iX==domain.x1+1);
            iX = domain.x0;
            ++iY;
            if (iY > domain.y1) {
                PLB_ASSERT(iY==domain.y1+1);
                iY = domain.y0;
                ++iZ;
            }
        }
    }
    return &buffer[0];
}

template<typename T>
bool MultiBlockSerializer3D<T>::isEmpty() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iZ > domain.z1;
    }
}

template<typename T>
MultiBlockDistribution3D const& MultiBlockSerializer3D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool MultiBlockSerializer3D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}

template<typename T>
void MultiBlockSerializer3D<T>::computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T> const& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send (
                Box3D(localX,localX+nextChunkSize-1, localY, localY, localZ, localZ),
                &buffer[0] );
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
}

template<typename T>
void MultiBlockSerializer3D<T>::computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T> const& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send (
                Box3D(localX,localX, localY, localY+nextChunkSize-1, localZ, localZ),
                &buffer[0] );
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
}

template<typename T>
void MultiBlockSerializer3D<T>::computeBufferAlongZ(plint nextBlockId, plint nextChunkSize) const {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T> const& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().send (
                Box3D(localX,localX, localY,localY, localZ, localZ+nextChunkSize-1),
                &buffer[0] );
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
}

template<typename T>
void MultiBlockSerializer3D<T>::communicateBuffer (
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
void MultiBlockSerializer3D<T>::fillBufferWithZeros(plint nextChunkSize) const {
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


////////// class MultiBlockUnSerializer3D ////////////////////////////

template<typename T>
MultiBlockUnSerializer3D<T>::MultiBlockUnSerializer3D (
        MultiBlock3D<T>& multiBlock_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(multiBlock.getBoundingBox()),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockUnSerializer3D<T>::MultiBlockUnSerializer3D (
        MultiBlock3D<T>& multiBlock_,
        Box3D domain_,
        IndexOrdering::OrderingT ordering_ )
    : multiBlock(multiBlock_),
      ordering(ordering_),
      domain(domain_),
      iX(domain.x0), iY(domain.y0), iZ(domain.z0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
MultiBlockUnSerializer3D<T>* MultiBlockUnSerializer3D<T>::clone() const {
    return new MultiBlockUnSerializer3D<T>(*this);
}


template<typename T>
pluint MultiBlockUnSerializer3D<T>::getSize() const {
    if (ordering==IndexOrdering::memorySaving) {
        return getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
    }
    else {
        return domain.nCells() * multiBlock.sizeOfCell();
    }
}

template<typename T>
T* MultiBlockUnSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) {
    PLB_PRECONDITION( !isFull() );
    plint nextBlockId, nextChunkSize;
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        bool hasData = getMultiBlockDistribution().getNextChunkZ(iX,iY,iZ, nextBlockId, nextChunkSize);
        if (iZ+nextChunkSize>domain.z1+1) {
            nextChunkSize = domain.z1-iZ+1;
        }
        bufferSize = nextChunkSize * multiBlock.sizeOfCell();
        if (ordering==IndexOrdering::memorySaving && !hasData) {
            bufferSize = 0;
        }
    }
    else {
        getMultiBlockDistribution().getNextChunkX(iX,iY,iZ, nextBlockId, nextChunkSize);
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
void MultiBlockUnSerializer3D<T>::commitData() {
    PLB_PRECONDITION( !isFull() );
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkZ(iX,iY,iZ, nextBlockId, nextChunkSize);
        if (iZ+nextChunkSize>domain.z1+1) {
            nextChunkSize = domain.z1-iZ+1;
        }
        if (hasData) {
            fillBufferAlongZ(nextBlockId, nextChunkSize);
        }
        iZ += nextChunkSize;
        if (iZ > domain.z1) {
            iZ = domain.z0;
            ++iY;
            if (iY > domain.y1) {
                iY = domain.y0;
                ++iX;
            }
        }
    }
    else {
        plint nextBlockId, nextChunkSize;
        bool hasData = getMultiBlockDistribution().getNextChunkX(iX,iY,iZ, nextBlockId, nextChunkSize);
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
            if (iY > domain.y1) {
                iY = domain.y0;
                ++iZ;
            }
        }
    }
    // At the end of unserialization, duplicate overlaps.
    if (isFull()) {
        multiBlock.getBlockCommunicator().duplicateOverlaps(multiBlock);
    }
}

template<typename T>
bool MultiBlockUnSerializer3D<T>::isFull() const {
    if (ordering==IndexOrdering::forward || ordering==IndexOrdering::memorySaving) {
        return iX > domain.x1;
    }
    else {
        return iZ > domain.z1;
    }
}

template<typename T>
MultiBlockDistribution3D const& MultiBlockUnSerializer3D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool MultiBlockUnSerializer3D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}

template<typename T>
void MultiBlockUnSerializer3D<T>::fillBufferAlongX(plint nextBlockId, plint nextChunkSize) {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
    if (blockIsLocal) {
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T>& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive (
                Box3D(localX,localX+nextChunkSize-1, localY, localY, localZ, localZ),
                &buffer[0] );
    }
}

template<typename T>
void MultiBlockUnSerializer3D<T>::fillBufferAlongY(plint nextBlockId, plint nextChunkSize) {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
    if (blockIsLocal) {
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T>& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive (
                Box3D(localX,localX, localY, localY+nextChunkSize-1, localZ,localZ), &buffer[0]);
    }
}

template<typename T>
void MultiBlockUnSerializer3D<T>::fillBufferAlongZ(plint nextBlockId, plint nextChunkSize) {
    plint bufferSize = nextChunkSize*multiBlock.sizeOfCell();
    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(nextBlockId);
    bool blockIsLocal = isLocal(nextBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);
    if (blockIsLocal) {
        plint localX = param.toLocalX(iX);
        plint localY = param.toLocalY(iY);
        plint localZ = param.toLocalZ(iZ);
        AtomicBlock3D<T>& nextBlock = multiBlock.getComponent(nextBlockId);
        nextBlock.getDataTransfer().receive (
                Box3D(localX,localX, localY,localY, localZ, localZ+nextChunkSize-1),
                &buffer[0] );
    }
}

template<typename T>
void MultiBlockUnSerializer3D<T>::communicateBuffer (
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


////////// class RawMultiBlockSerializer3D ////////////////////////////

template<typename T>
RawMultiBlockSerializer3D<T>::RawMultiBlockSerializer3D (
        MultiBlock3D<T> const& multiBlock_ )
    : multiBlock(multiBlock_),
      currentBlockId(0),
      iX(0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
RawMultiBlockSerializer3D<T>* RawMultiBlockSerializer3D<T>::clone() const {
    return new RawMultiBlockSerializer3D<T>(*this);
}

template<typename T>
pluint RawMultiBlockSerializer3D<T>::getSize() const {
    return
        getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
}

template<typename T>
const T* RawMultiBlockSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) const {
    PLB_PRECONDITION( !isEmpty() );

    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(currentBlockId);
    plint ev = param.getEnvelopeWidth();
    plint nx = param.getBulk().getNx();
    plint ny = param.getBulk().getNy();
    plint nz = param.getBulk().getNz();
    plint chunkSize = ny*nz;
    bufferSize = chunkSize*multiBlock.sizeOfCell();

    bool blockIsLocal = isLocal(currentBlockId);
    if (blockIsLocal) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
        AtomicBlock3D<T> const& currentBlock = multiBlock.getComponent(currentBlockId);
        currentBlock.getDataTransfer().send (
                Box3D(iX+ev, iX+ev, ev, ny-1+ev, ev, nz-1+ev), &buffer[0] );
    }
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);

    ++iX;
    if (iX==nx) {
        iX=0;
        ++currentBlockId;
    }
    return &buffer[0];
}

template<typename T>
void RawMultiBlockSerializer3D<T>::communicateBuffer (
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
bool RawMultiBlockSerializer3D<T>::isEmpty() const {
    return currentBlockId >= getMultiBlockDistribution().getNumBlocks();
}

template<typename T>
MultiBlockDistribution3D const& RawMultiBlockSerializer3D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool RawMultiBlockSerializer3D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}


////////// class RawMultiBlockUnSerializer3D ////////////////////////////

template<typename T>
RawMultiBlockUnSerializer3D<T>::RawMultiBlockUnSerializer3D(MultiBlock3D<T>& multiBlock_)
    : multiBlock(multiBlock_),
      currentBlockId(0),
      iX(0),
      buffer(1) // this avoids buffer of size 0 which one cannot point to
{ }

template<typename T>
RawMultiBlockUnSerializer3D<T>* RawMultiBlockUnSerializer3D<T>::clone() const {
    return new RawMultiBlockUnSerializer3D<T>(*this);
}


template<typename T>
pluint RawMultiBlockUnSerializer3D<T>::getSize() const {
    return
        getMultiBlockDistribution().getNumAllocatedBulkCells() * multiBlock.sizeOfCell();
}

template<typename T>
T* RawMultiBlockUnSerializer3D<T>::getNextDataBuffer(pluint& bufferSize) {
    PLB_PRECONDITION( !isFull() );

    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(currentBlockId);
    plint ny = param.getBulk().getNy();
    plint nz = param.getBulk().getNz();
    plint chunkSize = ny*nz;
    bufferSize = chunkSize*multiBlock.sizeOfCell();

    if (isLocal(currentBlockId)) {
        // The +1 avoids having a buffer of size 0 which one cannot point to
        buffer.resize(bufferSize+1);
    }
    return &buffer[0];
}

template<typename T>
void RawMultiBlockUnSerializer3D<T>::commitData() {
    PLB_PRECONDITION( !isFull() );

    BlockParameters3D const& param = getMultiBlockDistribution().getBlockParameters(currentBlockId);
    plint ev = param.getEnvelopeWidth();
    plint nx = param.getBulk().getNx();
    plint ny = param.getBulk().getNy();
    plint nz = param.getBulk().getNz();
    plint chunkSize = ny*nz;
    plint bufferSize = chunkSize*multiBlock.sizeOfCell();

    bool blockIsLocal = isLocal(currentBlockId);
    communicateBuffer(bufferSize, param.getProcId(), blockIsLocal);

    if (blockIsLocal) {
        AtomicBlock3D<T>& currentBlock = multiBlock.getComponent(currentBlockId);
        currentBlock.getDataTransfer().receive (
                Box3D(iX+ev, iX+ev, ev, ny-1+ev, ev, nz-1+ev), &buffer[0] );
    }

    ++iX;
    if (iX==nx) {
        iX=0;
        ++currentBlockId;
    }

    // At the end of unserialization, duplicate overlaps.
    if (isFull()) {
        multiBlock.getBlockCommunicator().duplicateOverlaps(multiBlock);
    }
}

template<typename T>
bool RawMultiBlockUnSerializer3D<T>::isFull() const {
    return currentBlockId >= getMultiBlockDistribution().getNumBlocks();
}

template<typename T>
MultiBlockDistribution3D const& RawMultiBlockUnSerializer3D<T>::getMultiBlockDistribution() const {
    return multiBlock.getMultiBlockManagement().getMultiBlockDistribution();
}

template<typename T>
bool RawMultiBlockUnSerializer3D<T>::isLocal(plint whichBlock) const {
    return multiBlock.getMultiBlockManagement().getThreadAttribution().isLocal (
            getMultiBlockDistribution().getBlockParameters(whichBlock).getProcId() );
}

template<typename T>
void RawMultiBlockUnSerializer3D<T>::communicateBuffer (
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

#endif  // MULTI_BLOCK_SERIALIZER_3D_HH
