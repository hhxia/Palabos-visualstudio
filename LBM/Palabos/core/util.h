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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef UTIL_H
#define UTIL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

namespace util {

/// Compute square of a scalar value
template<typename T>
T sqr(T arg) {
    return arg*arg;
}

/// Test equality between two boolean values.
/** As we were testing with GCC 4.3.3, we observed that the
 *  equality operator (==) has an unexpected behavior. It
 *  compares the bools as if they were ints, and sometimes
 *  returns false, even if both bools are true. We don't know
 *  if this is a weird definition of the C++ language or a compiler
 *  bug, but in both cases, the present function offers a workaround.
 */
inline bool boolIsEqual(bool val1, bool val2) {
    return ( (val1 && val2) || !(val1 || val2) );
}

/// Round to next plint value
template<typename T>
plint roundToInt(T value) {
    return value>0 ?
               static_cast<plint>(value+(T)0.5) :
               static_cast<plint>(value-(T)0.5);
}

/// Round a signed integer to the next larger value which is divisible by step.
inline plint roundUp(plint value, plint step) {
    plint modulo = value % step;
    plint result = value-modulo;
    if (modulo>0) {
        result += step;
    }
    return result;
}

/// Round a signed integer to the next smaller value which is divisible by step.
inline plint roundDown(plint value, plint step) {
    plint modulo = value % step;
    plint result = value-modulo;
    if (modulo<0) {
        result -= step;
    }
    return result;
}


/// A simple class for handling buffer memory
/** This class can be seen as a replacement of the std::vector
 *  template. It is less powerful, but at least, unlike std::vector,
 *  it treats the type bool appropriately (this is the only reason
 *  why class Buffer was written). Note that class Buffer is not
 *  STL compatible (it is not an STL container), but it is thread-
 *  safe.
 */
template<typename T>
class Buffer {
public:
    /// Default constructor allcoates no memory.
    Buffer()
        : bufferSize(0),
          data(0)
    { }
    /// Constructor with buffer size; buffer elements are not default-initialized.
    Buffer(pluint bufferSize_)
        : bufferSize(bufferSize_),
          data(new T[bufferSize])
    { }
    Buffer(Buffer<T> const& rhs)
        : bufferSize(rhs.bufferSize),
          data(new T[bufferSize])
    {
        for (pluint iData=0; iData<bufferSize; ++iData) {
            data[iData] = rhs.data[iData];
        }
    }
    ~Buffer() {
        delete [] data;
    }
    /// Assignment-operator is thread-safe.
    Buffer<T>& operator=(Buffer<T> const& rhs) {
        Buffer<T> newBuffer(rhs);
        swap(newBuffer);
        return *this;
    }
    /// Swap with other buffer; this operation cannot throw.
    void swap(Buffer<T>& rhs) {
        std::swap(bufferSize, rhs.bufferSize);
        std::swap(data, rhs.data);
    }
    /// Resize the buffer, and keep values from before.
    void resize(pluint newBufferSize) {
        if (newBufferSize > bufferSize) {
            Buffer<T> newBuffer(newBufferSize);
            for (int iData=0; iData<bufferSize; ++iData) {
                newBuffer.data[iData] = data[iData];
            }
            swap(newBuffer);
        }
    }
    /// Resize the buffer, but don't keep values from before.
    void reallocate(pluint newBufferSize) {
        if (newBufferSize > bufferSize) {
            Buffer<T> newBuffer(newBufferSize);
            swap(newBuffer);
        }
    }
    /// Pointer to raw data.
    T* get() {
        return data;
    }
    /// Const-pointer to raw data.
    T const* get() const {
        return data;
    }
    /// Element access.
    T& operator[](pluint index) {
        PLB_PRECONDITION( index<bufferSize );
        return data[index];
    }
    /// Const-element access.
    T const& operator[](pluint index) const {
        PLB_PRECONDITION( index<bufferSize );
        return data[index];
    }
private:
    pluint bufferSize;
    T* data;
};

}  // namespace util

}  // namespace plb

#endif
