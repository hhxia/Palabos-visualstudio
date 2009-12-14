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


#ifndef SERIALIZER_HH
#define SERIALIZER_HH

#include "parallelism/mpiManager.h"
#include "core/serializer.h"
#include "core/plbDebug.h"
#include <algorithm>

using namespace std;

namespace plb {

////////// class ScalingSerializer ////////////////////////////

template<typename T>
ScalingSerializer<T>::ScalingSerializer(DataSerializer<T> const* baseSerializer_, T scalingFactor_)
    : baseSerializer(baseSerializer_),
      scalingFactor(scalingFactor_)
{ }

template<typename T>
ScalingSerializer<T>::ScalingSerializer(ScalingSerializer<T> const& rhs)
    : baseSerializer(rhs.baseSerializer->clone()),
      scalingFactor(rhs.scalingFactor)
{ }

template<typename T>
ScalingSerializer<T>& ScalingSerializer<T>::operator=(ScalingSerializer<T> const& rhs)
{
    delete baseSerializer; baseSerializer = rhs.baseSerializer->clone();
    scalingFactor = rhs.scalingFactor;
    return *this;
}

template<typename T>
ScalingSerializer<T>* ScalingSerializer<T>::clone() const
{
    return new ScalingSerializer<T>(*this);
}

template<typename T>
ScalingSerializer<T>::~ScalingSerializer()
{
    delete baseSerializer;
}


template<typename T>
pluint ScalingSerializer<T>::getSize() const {
    return baseSerializer->getSize();
}

template<typename T>
const T* ScalingSerializer<T>::getNextDataBuffer(pluint& bufferSize) const {
    const T* unscaledBuffer = baseSerializer->getNextDataBuffer(bufferSize);
    scaledBuffer.reallocate(bufferSize);
    for (pluint iBuffer=0; iBuffer<bufferSize; ++iBuffer) {
        scaledBuffer[iBuffer] = unscaledBuffer[iBuffer] * scalingFactor;
    }
    return scaledBuffer.get();
}

template<typename T>
bool ScalingSerializer<T>::isEmpty() const {
    return baseSerializer->isEmpty();
}


////////// class TypeConversionSerializer ////////////////////////////

template<typename T, typename TConv>
TypeConversionSerializer<T,TConv>::TypeConversionSerializer (
        DataSerializer<T> const* baseSerializer_)
    : baseSerializer(baseSerializer_)
{ }

template<typename T, typename TConv>
TypeConversionSerializer<T,TConv>::TypeConversionSerializer (
        TypeConversionSerializer<T,TConv> const& rhs)
    : baseSerializer(rhs.baseSerializer->clone())
{ }

template<typename T, typename TConv>
TypeConversionSerializer<T,TConv>&
    TypeConversionSerializer<T,TConv>::operator=(TypeConversionSerializer<T,TConv> const& rhs)
{
    delete baseSerializer; baseSerializer = rhs.baseSerializer->clone();
    return *this;
}

template<typename T, typename TConv>
TypeConversionSerializer<T,TConv>*
    TypeConversionSerializer<T,TConv>::clone() const
{
    return new TypeConversionSerializer<T,TConv>(*this);
}

template<typename T, typename TConv>
TypeConversionSerializer<T,TConv>::~TypeConversionSerializer<T,TConv>() {
    delete baseSerializer;
}

template<typename T, typename TConv>
pluint TypeConversionSerializer<T,TConv>::getSize() const {
    return baseSerializer->getSize();
}

template<typename T, typename TConv>
const TConv* TypeConversionSerializer<T,TConv>::getNextDataBuffer(pluint& bufferSize) const {
    const T* originalBuffer = baseSerializer->getNextDataBuffer(bufferSize);
    convBuffer.reallocate(bufferSize);
    for (pluint iBuffer=0; iBuffer<bufferSize; ++iBuffer) {
        convBuffer[iBuffer] = static_cast<TConv>( originalBuffer[iBuffer] );
    }
    return convBuffer.get();
}

template<typename T, typename TConv>
bool TypeConversionSerializer<T,TConv>::isEmpty() const {
    return baseSerializer->isEmpty();
}


////////// Free functions ////////////////////////////

template<typename T>
void serializerToUnSerializer(DataSerializer<T> const* serializer, DataUnSerializer<T>* unSerializer) {
    PLB_PRECONDITION( serializer->getSize() == unSerializer->getSize() );
    pluint writePos = 0, readPos = 0;
    pluint serializerBufferSize =0, unSerializerBufferSize =0;
    const T* serializerBuffer =0;
    T* unSerializerBuffer =0;
    while (!unSerializer->isFull()) {
        if (readPos==serializerBufferSize) {
            serializerBuffer = serializer->getNextDataBuffer(serializerBufferSize);
            readPos = 0;
        }
        if (writePos==unSerializerBufferSize) {
            unSerializerBuffer = unSerializer->getNextDataBuffer(unSerializerBufferSize);
            writePos = 0;
        }

        pluint remainToRead = (plint)serializerBufferSize - (plint)readPos;
        pluint remainToWrite = (plint)unSerializerBufferSize - (plint)writePos;
        pluint nextChunk = min(remainToRead, remainToWrite);
        for (pluint iChunk=0; iChunk<nextChunk; ++iChunk, ++readPos, ++writePos) {
            if (global::mpi().isMainProcessor()) {
                unSerializerBuffer[writePos] = serializerBuffer[readPos];
            }
        }
        if (writePos==unSerializerBufferSize) {
            unSerializer->commitData();
        }
    }
    delete serializer;
    delete unSerializer;
}

template<typename T>
void serializerToSink(DataSerializer<T> const* serializer, SerializedWriter<T>* sink) {
    if (global::mpi().isMainProcessor()) {
        pluint dataSize = serializer->getSize();
        sink->writeHeader(dataSize);
    }
    while (!serializer->isEmpty()) {
        pluint bufferSize;
        const T* dataBuffer = serializer->getNextDataBuffer(bufferSize);
        if (global::mpi().isMainProcessor()) {
            sink->writeData(dataBuffer, bufferSize);
        }
    }
    delete sink;
    delete serializer;
}

template<typename T>
void sourceToUnSerializer(SerializedReader<T> const* source, DataUnSerializer<T>* unSerializer) {
    if (global::mpi().isMainProcessor()) {
        pluint dataSize = unSerializer->getSize();
        source->readHeader(dataSize);
    }
    while (!unSerializer->isFull()) {
        pluint bufferSize = 0;
        T* dataBuffer = unSerializer->getNextDataBuffer(bufferSize);
        if (global::mpi().isMainProcessor()) {
            source->readData(dataBuffer, bufferSize);
        }
        unSerializer->commitData();
    }
    delete source;
    delete unSerializer;
}

}  // namespace plb

#endif


