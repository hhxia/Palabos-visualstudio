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
 * Scalar, vector and tensor fields for 2D data analysis -- template instantiation.
 */

#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataField2D.hh"

namespace plb {

template class MultiScalarField2D<double>;
template class MultiTensorField2D<double, 2>; // for vectors
template class MultiTensorField2D<double, 3>; // for symmetric tensors
template class MultiTensorField2D<double, 4>; // for non-symmetric tensors
template class MultiTensorField2D<double, 9>; // for full populations

}  // namespace plb
