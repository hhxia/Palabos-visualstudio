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


#ifndef ENDIANNESS_H
#define ENDIANNESS_H

#include "core/globalDefs.h"
#include <algorithm>

namespace plb {

template<int nBytes> inline void endianByteSwapImpl(char const* src, char* dest);
template<int nBytes> inline void endianByteSwapImpl(char* value);

template<typename T> inline void endianByteSwap(T const& src, T& dest) {
    endianByteSwapImpl<sizeof(T)> (
            reinterpret_cast< char const* > (&src), reinterpret_cast< char* > (&dest) );
}

template<typename T> inline void endianByteSwap(T& value) {
    endianByteSwapImpl<sizeof(T)> ( reinterpret_cast< char* > (&value) );
}

// Specialization for 1-byte types.
template<>
inline void endianByteSwapImpl<1>(char const* src, char* dest)
{
    dest[0] = src[0];
}

// Specialization for 1-byte types.
template<>
inline void endianByteSwapImpl<1>(char* value)
{ }

// Specialization for 2-byte types.
template<>
inline void endianByteSwapImpl<2>(char const* src, char* dest)
{
    dest[0] = src[1];
    dest[1] = src[0];
}

// Specialization for 2-byte types.
template<>
inline void endianByteSwapImpl<2>(char* value)
{
    std::swap(value[0], value[1]);
}

// Specialization for 4-byte types.
template<>
inline void endianByteSwapImpl<4>(char const* src, char* dest)
{
    dest[0] = src[3];
    dest[1] = src[2];
    dest[2] = src[1];
    dest[3] = src[0];
}

// Specialization for 4-byte types.
template<>
inline void endianByteSwapImpl<4>(char* value)
{
    std::swap(value[0], value[3]);
    std::swap(value[1], value[2]);
}

// Specialization for 8-byte types.
template<>
inline void endianByteSwapImpl<8>(char const* src, char* dest)
    {
    dest[0] = src[7];
    dest[1] = src[6];
    dest[2] = src[5];
    dest[3] = src[4];
    dest[4] = src[3];
    dest[5] = src[2];
    dest[6] = src[1];
    dest[7] = src[0];
}

// Specialization for 8-byte types.
template<>
inline void endianByteSwapImpl<8>(char* value)
{
    std::swap(value[0], value[7]);
    std::swap(value[1], value[6]);
    std::swap(value[2], value[5]);
    std::swap(value[3], value[4]);
}


// Specialization for 12-byte types.
template<>
inline void endianByteSwapImpl<12>(char const* src, char* dest)
    {
    dest[0] = src[11];
    dest[1] = src[10];
    dest[2] = src[9];
    dest[3] = src[8];
    dest[4] = src[7];
    dest[5] = src[6];
    dest[6] = src[5];
    dest[7] = src[4];
    dest[8] = src[3];
    dest[9] = src[2];
    dest[10] = src[1];
    dest[11] = src[0];
}

// Specialization for 12-byte types.
template<>
inline void endianByteSwapImpl<12>(char* value)
{
    std::swap(value[0], value[11]);
    std::swap(value[1], value[10]);
    std::swap(value[2], value[9]);
    std::swap(value[3], value[8]);
    std::swap(value[4], value[7]);
    std::swap(value[5], value[6]);
}

} // namespace plb

#endif  // ENDIANNESS_H
