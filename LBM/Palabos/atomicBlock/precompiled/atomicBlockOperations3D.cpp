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
 * A 3D multiblock -- template instantiation.
 */

#include "atomicBlock/atomicBlockOperations3D.h"
#include "atomicBlock/atomicBlockOperations3D.hh"

namespace plb {

template void executeDataProcessor<double> (
        DataProcessorGenerator3D<double> const& generator,
        std::vector<AtomicBlock3D<double>*> objects );

template void executeDataProcessor<double> (
        ReductiveDataProcessorGenerator3D<double>& generator,
        std::vector<AtomicBlock3D<double>*> objects );

template void addInternalProcessor<double> (
        DataProcessorGenerator3D<double> const& generator,
        std::vector<AtomicBlock3D<double>*> objects, plint level );

}
