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

/* Acknowledgment: The strategy adopted here to encode
 * and decode Base64, and in particular the expression of the
 * arrays Base64Encoder::enc64 and Base64Decoder::dec64, 
 * are inspired by the open source library b64 by Bob Trower,
 * which is distributed under an MIT license at the address
 * http://base64.sourceforge.net/b64.c
 */


#ifndef BASE64_H
#define BASE64_H

#include "core/globalDefs.h"
#include <iosfwd>

namespace plb {

template<typename T>
class Base64Encoder {
public:
    Base64Encoder(std::ostream& ostr_, pluint fullLength_);
    void encode(const T* data, pluint length);
private:
    void fillOverflow(const unsigned char* charData, pluint charLength, pluint& pos);
    void flushOverflow();
    void writeSize();
    void encodeBlock( const unsigned char* data);
    void encodeUnfinishedBlock( const unsigned char* data, plint length);
private:
    static const char enc64[65];
private:
    std::ostream& ostr;
    pluint charFullLength;
    pluint numWritten;
    plint numOverflow;
    unsigned char overflow[3];
};

template<typename T>
class Base64Decoder {
public:
    Base64Decoder(std::istream& istr_, pluint fullLength_);
    void decode(T* data, pluint length);
private:
    void flushOverflow(unsigned char* charData, pluint charLength, pluint& pos);
    unsigned char getNext();
    void decodeBlock(unsigned char* data);
private:
    static const char dec64[82];
private:
    std::istream& istr;
    pluint charFullLength;
    pluint numRead;
    plint posOverflow;
    unsigned char overflow[3];
};

} // namespace plb

#endif
