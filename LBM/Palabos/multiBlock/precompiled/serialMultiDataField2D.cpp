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
 * Helper classes for serial 2D multiblock lattice -- template instantiation.
 */

#include "multiBlock/serialMultiDataField2D.h"
#include "multiBlock/serialMultiDataField2D.hh"

namespace plb {

    template class SerialScalarAccess2D<double>;
    template class SerialTensorAccess2D<double,2>;
    template class SerialTensorAccess2D<double,3>;
    template class SerialTensorAccess2D<double,4>;
    template class SerialTensorAccess2D<double,9>;

}
