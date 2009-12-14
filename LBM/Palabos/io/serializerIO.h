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

#ifndef SERIALIZER_IO_H
#define SERIALIZER_IO_H

#include "core/globalDefs.h"
#include "core/serializer.h"
#include <iosfwd>

namespace plb {

/// Take a Serializer, convert into Base64 format (ASCII based binary representation), and stream into output stream.
/** Ahead of the data, an integer value is encoded which stands for the total size of
 *  serialized data. For compatibility with the VTK file format (as of this writing),
 *  you can enforce that the type of this variable is converted to "unsigned int".
 *  Note that this may lead to errors on 64-bit platforms, if the total amount of
 *  data exceeds 2 GB.
 */
template<typename T>
void serializerToBase64Stream(DataSerializer<T> const* serializer, std::ostream& ostr, bool enforceUint=false);

/// Take an input stream with Base64 encoded binary content, and stream into an unSerializer
/** If the integer value which indicates the amount of data to be unSerialized is of type
 *  "unsigned int", this fact can be enforced with the flag enforceUplint to ensure
 *  compatibility between 32-bit and 64-bit platforms.
 */
template<typename T>
void base64StreamToUnSerializer(std::istream& istr, DataUnSerializer<T>* unSerializer, bool enforceUint=false);

/// Take a Serializer, convert and stream into output in ASCII format.
/** Number of digits in the ASCII representation of numbers is given by the variable numDigits.
 */
template<typename T>
void serializerToAsciiStream(DataSerializer<T> const* serializer, std::ostream& ostr, plint numDigits=8);


/// Take an UnSerializer and fill it with data from an ASCII-format input stream.
/** Number of digits in the ASCII representation of numbers is given by the variable numDigits.
 */
template<typename T>
void asciiStreamToUnSerializer(std::istream& istr, DataUnSerializer<T>* unSerializer);

} // namespace plb

#endif
