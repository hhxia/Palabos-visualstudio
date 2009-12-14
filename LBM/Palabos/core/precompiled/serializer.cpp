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
 * serializer and unserializer -- template instantiation.
 */

#include "parallelism/mpiManager.h"
#include "core/serializer.h"
#include "core/serializer.hh"

namespace plb {

template class ScalingSerializer<double>;
template class ScalingSerializer<int>;
template class TypeConversionSerializer<double,double>;
template class TypeConversionSerializer<double,float>;

template void serializerToUnSerializer<double>(DataSerializer<double> const* serializer,
                                               DataUnSerializer<double>* unSerializer);

template void serializerToSink<double>(DataSerializer<double> const* serializer,
                                       SerializedWriter<double>* sink);

template void sourceToUnSerializer<double>(SerializedReader<double> const* source,
                                           DataUnSerializer<double>* unSerializer);

}  // namespace plb
