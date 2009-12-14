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
 * Input/Output based on serializers -- template instantiation.
 */
#include "io/serializerIO.h"
#include "io/serializerIO.hh"

namespace plb {

template class Base64Writer<double>;
template class Base64Reader<double>;
template class AsciiWriter<double>;
template class AsciiReader<double>;

template
void serializerToBase64Stream<double>(DataSerializer<double> const* serializer, std::ostream& ostr, bool enforceUint);

template
void base64StreamToUnSerializer<double>(std::istream& istr, DataUnSerializer<double>* unSerializer, bool enforceUint);

template
void serializerToAsciiStream<double>(DataSerializer<double> const* serializer, std::ostream& ostr, plint numDigits);

template
void asciiStreamToUnSerializer<double>(std::istream& istr, DataUnSerializer<double>* unSerializer);

}
