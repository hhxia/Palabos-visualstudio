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
#ifndef MULTI_BLOCK_SERIALIZER_2D_H
#define MULTI_BLOCK_SERIALIZER_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"
#include "core/serializer.h"

namespace plb {

template<typename T>
class MultiBlockSerializer2D : public DataSerializer<T> {
public:
    MultiBlockSerializer2D(MultiBlock2D<T> const& multiBlock_,
                           IndexOrdering::OrderingT ordering_);
    MultiBlockSerializer2D(MultiBlock2D<T> const& multiBlock_,
                           Box2D domain_,
                           IndexOrdering::OrderingT ordering_);
    virtual MultiBlockSerializer2D<T>* clone() const;
    virtual pluint getSize() const;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    MultiBlockDistribution2D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
    void computeBufferAlongX(plint nextBlockId, plint nextChunkSize) const;
    void computeBufferAlongY(plint nextBlockId, plint nextChunkSize) const;
    void communicateBuffer(plint bufferSize, plint fromProc, bool isAllocated) const;
    void fillBufferWithZeros(plint nextChunkSize) const;
private:
    MultiBlock2D<T> const& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint iX, iY;
    mutable std::vector<T> buffer;
};

template<typename T>
class MultiBlockUnSerializer2D : public DataUnSerializer<T> {
public:
    MultiBlockUnSerializer2D(MultiBlock2D<T>& multiBlock_,
                             IndexOrdering::OrderingT ordering_);
    MultiBlockUnSerializer2D(MultiBlock2D<T>& multiBlock_,
                             Box2D domain_,
                             IndexOrdering::OrderingT ordering_);
    virtual MultiBlockUnSerializer2D<T>* clone() const;
    virtual pluint getSize() const;
    virtual T* getNextDataBuffer(pluint& bufferSize);
    virtual void commitData();
    virtual bool isFull() const;
private:
    MultiBlockDistribution2D const& getMultiBlockDistribution() const;
    bool isLocal(plint whichBlock) const;
    void fillBufferAlongX(plint nextBlockId, plint nextChunkSize);
    void fillBufferAlongY(plint nextBlockId, plint nextChunkSize);
    void communicateBuffer(plint bufferSize, plint toBlock, bool isAllocated) const;
private:
    MultiBlock2D<T>& multiBlock;
    IndexOrdering::OrderingT ordering;
    Box2D domain;
    mutable plint iX, iY;
    mutable std::vector<T> buffer;
};

}  //  namespace plb

#endif  // MULTI_BLOCK_SERIALIZER_2D_H
