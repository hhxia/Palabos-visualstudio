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


#ifndef SERIALIZER_H
#define SERIALIZER_H

#include "core/globalDefs.h"
#include "core/util.h"

namespace plb {

template<typename T>
struct DataSerializer {
    virtual ~DataSerializer() { }
    virtual DataSerializer<T>* clone() const =0;
    virtual pluint getSize() const =0;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const =0;
    virtual bool isEmpty() const =0;
};

template<typename T>
struct DataUnSerializer {
    virtual ~DataUnSerializer() { }
    virtual DataUnSerializer<T>* clone() const =0;
    virtual pluint getSize() const =0;
    virtual T* getNextDataBuffer(pluint& bufferSize) =0;
    virtual void commitData() =0;
    virtual bool isFull() const =0;
};

template<typename T>
struct SerializedWriter {
    virtual ~SerializedWriter() { }
    virtual SerializedWriter<T>* clone() const =0;
    virtual void writeHeader(pluint dataSize) =0;
    virtual void writeData(T const* dataBuffer, pluint bufferSize) =0;
};


template<typename T>
struct SerializedReader {
    virtual ~SerializedReader() { }
    virtual SerializedReader<T>* clone() const =0;
    virtual void readHeader(pluint dataSize) const =0;
    virtual void readData(T* dataBuffer, pluint bufferSize) const =0;
};

template<typename T>
class ScalingSerializer : public DataSerializer<T> {
public:
    ScalingSerializer(DataSerializer<T> const* baseSerializer_, T scalingFactor_);
    ScalingSerializer(ScalingSerializer<T> const& rhs);
    ScalingSerializer<T>& operator=(ScalingSerializer<T> const& rhs);
    ScalingSerializer<T>* clone() const;
    ~ScalingSerializer<T>();
    virtual pluint getSize() const;
    virtual const T* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    DataSerializer<T> const* baseSerializer;
    T scalingFactor;
    mutable util::Buffer<T> scaledBuffer;
};

template<typename T, typename TConv>
class TypeConversionSerializer : public DataSerializer<TConv> {
public:
    TypeConversionSerializer(DataSerializer<T> const* baseSerializer_);
    TypeConversionSerializer(TypeConversionSerializer<T,TConv> const& rhs);
    TypeConversionSerializer<T,TConv>& operator=(TypeConversionSerializer<T,TConv> const& rhs);
    TypeConversionSerializer<T, TConv>* clone() const;
    ~TypeConversionSerializer();
    virtual pluint getSize() const;
    virtual const TConv* getNextDataBuffer(pluint& bufferSize) const;
    virtual bool isEmpty() const;
private:
    DataSerializer<T> const* baseSerializer;
    mutable util::Buffer<TConv> convBuffer;
};

template<typename T>
void serializerToUnSerializer(DataSerializer<T> const* serializer, DataUnSerializer<T>* unSerializer);

template<typename T>
void serializerToSink(DataSerializer<T> const* serializer, SerializedWriter<T>* sink);

template<typename T>
void sourceToUnSerializer(SerializedReader<T> const* source, DataUnSerializer<T>* unSerializer);


} // namespace plb

#endif
