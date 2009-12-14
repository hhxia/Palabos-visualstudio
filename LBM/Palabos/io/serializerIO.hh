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

#ifndef SERIALIZER_IO_HH
#define SERIALIZER_IO_HH

#include "io/serializerIO.h"
#include "io/base64.h"
#include "io/endianness.h"
#include "core/plbDebug.h"
#include "core/globalDefs.h"
#include <vector>
#include <limits>
#include <iomanip>
#include <istream>
#include <ostream>
#include <fstream>

namespace plb {

/* *************** Class Base64Writer ******************************** */

template<typename T>
class Base64Writer : public SerializedWriter<T> {
public:
    Base64Writer(std::ostream& ostr_, bool enforceUint_, bool switchEndianness_);
    Base64Writer(Base64Writer<T> const& rhs);
    virtual Base64Writer<T>* clone() const;
    ~Base64Writer();
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(T const* dataBuffer, pluint bufferSize);
private:
    std::ostream& ostr;
    bool enforceUint;
    bool switchEndianness;
    Base64Encoder<T>* dataEncoder;
};

template<typename T>
Base64Writer<T>::Base64Writer(std::ostream& ostr_, bool enforceUint_, bool switchEndianness_)
    : ostr(ostr_),
      enforceUint(enforceUint_),
      switchEndianness(switchEndianness_),
      dataEncoder(0)
{ }

template<typename T>
Base64Writer<T>::Base64Writer(Base64Writer<T> const& rhs)
    : ostr(rhs.ostr),
      enforceUint(rhs.enforceUint),
      dataEncoder(0)
{
    if (rhs.dataEncoder) {
        dataEncoder = new Base64Encoder<T>(*rhs.dataEncoder);
    }
}

template<typename T>
Base64Writer<T>* Base64Writer<T>::clone() const {
    return new Base64Writer<T>(*this);
}

template<typename T>
Base64Writer<T>::~Base64Writer() {
    delete dataEncoder;
}

template<typename T>
void Base64Writer<T>::writeHeader(pluint dataSize) {
    pluint binarySize = (pluint) (dataSize * sizeof(T));
    if (enforceUint) {
        Base64Encoder<unsigned int> sizeEncoder(ostr, 1);
        PLB_PRECONDITION(binarySize <= std::numeric_limits<unsigned int>::max());
        unsigned int uintBinarySize = (unsigned int)binarySize;
        if (switchEndianness) {
            endianByteSwap(uintBinarySize);
        }
        sizeEncoder.encode(&uintBinarySize, 1);
    }
    else {
        Base64Encoder<pluint> sizeEncoder(ostr, 1);
        if (switchEndianness) {
            endianByteSwap(binarySize);
        }
        sizeEncoder.encode(&binarySize, 1);
    }
    dataEncoder = new Base64Encoder<T>(ostr, dataSize);
}

template<typename T>
void Base64Writer<T>::writeData(T const* dataBuffer, pluint bufferSize)
{
    static std::vector<T> tmpData;
    if (switchEndianness) {
        tmpData.resize(bufferSize);
        for (pluint iData=0; iData<bufferSize; ++iData) {
            endianByteSwap(dataBuffer[iData], tmpData[iData]);
        }
        dataEncoder->encode(&tmpData[0], bufferSize);
    }
    else {
        dataEncoder->encode(dataBuffer, bufferSize);
    }
}


/* *************** Class Base64Reader ******************************** */

template<typename T>
class Base64Reader : public SerializedReader<T> {
public:
    Base64Reader(std::istream& istr_, bool enforceUint_, bool switchEndianness_);
    Base64Reader(Base64Reader<T> const& rhs);
    virtual Base64Reader<T>* clone() const;
    ~Base64Reader();
    virtual void readHeader(pluint dataSize) const;
    virtual void readData(T* dataBuffer, pluint bufferSize) const;
private:
    std::istream& istr;
    bool enforceUint;
    bool switchEndianness;
    mutable Base64Decoder<T>* dataDecoder;
};

template<typename T>
Base64Reader<T>::Base64Reader(std::istream& istr_, bool enforceUint_, bool switchEndianness_)
    : istr(istr_),
      enforceUint(enforceUint_),
      switchEndianness(switchEndianness_),
      dataDecoder(0)
{ }

template<typename T>
Base64Reader<T>::Base64Reader(Base64Reader<T> const& rhs)
    : istr(rhs.istr),
      enforceUint(rhs.enforceUint),
      dataDecoder(0)
{
    if (rhs.dataDecoder) {
        dataDecoder = new Base64Decoder<T>(*rhs.dataDecoder);
    }
}

template<typename T>
Base64Reader<T>* Base64Reader<T>::clone() const {
    return new Base64Reader<T>(*this);
}

template<typename T>
Base64Reader<T>::~Base64Reader() {
    delete dataDecoder;
}

template<typename T>
void Base64Reader<T>::readHeader(pluint dataSize) const {
    pluint binarySize = 0;
    if (enforceUint) {
        unsigned int uintBinarySize;
        Base64Decoder<unsigned int> sizeDecoder(istr, 1);
        sizeDecoder.decode(&uintBinarySize, 1);
        if (switchEndianness) {
            endianByteSwap(uintBinarySize);
        }
        binarySize = uintBinarySize;
    }
    else {
        Base64Decoder<pluint> sizeDecoder(istr, 1);
        sizeDecoder.decode(&binarySize, 1);
        if (switchEndianness) {
            endianByteSwap(binarySize);
        }
    }
#ifdef PLB_DEBUG
    pluint sizeFromInputFile = binarySize / sizeof(T);
#endif
    PLB_PRECONDITION(sizeFromInputFile == dataSize);

    dataDecoder = new Base64Decoder<T>(istr, dataSize);
}

template<typename T>
void Base64Reader<T>::readData(T* dataBuffer, pluint bufferSize) const
{
    static std::vector<T> tmpData;
    if (switchEndianness) {
        tmpData.resize(bufferSize);
        dataDecoder->decode(&tmpData[0], bufferSize);
        for (pluint iData=0; iData<bufferSize; ++iData) {
            endianByteSwap(tmpData[iData], dataBuffer[iData]);
        }
    }
    else {
        dataDecoder->decode(dataBuffer, bufferSize);
    }
}


/* *************** Class AsciiWriter ******************************** */

template<typename T>
class AsciiWriter : public SerializedWriter<T> {
public:
    AsciiWriter(std::ostream& ostr_, plint numDigits_);
    virtual AsciiWriter<T>* clone() const;
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(T const* dataBuffer, pluint bufferSize);
private:
    std::ostream& ostr;
    plint numDigits;
};


template<typename T>
AsciiWriter<T>::AsciiWriter(std::ostream& ostr_, plint numDigits_)
    : ostr(ostr_),
      numDigits(numDigits_)
{ }

template<typename T>
AsciiWriter<T>* AsciiWriter<T>::clone() const {
    return new AsciiWriter<T>(*this);
}

template<typename T>
void AsciiWriter<T>::writeHeader(pluint dataSize)
{ }

template<typename T>
void AsciiWriter<T>::writeData(T const* dataBuffer, pluint bufferSize) {
    // Write the output data.
    if (numDigits==0) {
        for (pluint iData=0; iData<bufferSize; ++iData) {
            ostr << dataBuffer[iData] << " ";
        }
    }
    else {
        for (pluint iData=0; iData<bufferSize; ++iData) {
            ostr << std::setprecision(numDigits) << dataBuffer[iData] << " ";
        }
    }
}

/* *************** Class AsciiReader ******************************** */

template<typename T>
class AsciiReader : public SerializedReader<T> {
public:
    AsciiReader(std::istream& istr_);
    virtual AsciiReader<T>* clone() const;
    virtual void readHeader(pluint dataSize) const;
    virtual void readData(T* dataBuffer, pluint bufferSize) const;
private:
    std::istream& istr;
};


template<typename T>
AsciiReader<T>::AsciiReader(std::istream& istr_)
    : istr(istr_)
{ }

template<typename T>
AsciiReader<T>* AsciiReader<T>::clone() const {
    return new AsciiReader<T>(*this);
}

template<typename T>
void AsciiReader<T>::readHeader(pluint dataSize) const
{ }

template<typename T>
void AsciiReader<T>::readData(T* dataBuffer, pluint bufferSize) const
{
    // Read the input data.
    for (pluint iData=0; iData<bufferSize; ++iData) {
        istr >> dataBuffer[iData];
    }
}



/* *************** Free functions ************************************ */

template<typename T>
    void serializerToBase64Stream(DataSerializer<T> const* serializer, std::ostream& ostr, bool enforceUint)
{
    serializerToSink (
            serializer,
            new Base64Writer<T>(ostr, enforceUint,
                                global::IOpolicy().getEndianSwitchOnBase64out()) );
}

template<typename T>
void base64StreamToUnSerializer(std::istream& istr, DataUnSerializer<T>* unSerializer, bool enforceUint) {
    sourceToUnSerializer (
            new Base64Reader<T>(istr, enforceUint,
                                global::IOpolicy().getEndianSwitchOnBase64in()),
            unSerializer);
}

template<typename T>
void serializerToAsciiStream(DataSerializer<T> const* serializer, std::ostream& ostr, plint numDigits)
{
    serializerToSink(serializer, new AsciiWriter<T>(ostr, numDigits));
}

template<typename T>
void asciiStreamToUnSerializer(std::istream& istr, DataUnSerializer<T>* unSerializer)
{
    sourceToUnSerializer(new AsciiReader<T>(istr), unSerializer);
}

} // namespace plb

#endif
