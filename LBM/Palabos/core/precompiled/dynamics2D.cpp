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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "core/dynamics.h"
#include "core/dynamics.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

    template class Dynamics<double, descriptors::D2Q9Descriptor>;
    template class BasicBulkDynamics<double, descriptors::D2Q9Descriptor>;
    template class CompositeDynamics<double, descriptors::D2Q9Descriptor>;
    template class PreparePopulationsDynamics<double, descriptors::D2Q9Descriptor>;
    template class BulkCompositeDynamics<double, descriptors::D2Q9Descriptor>;
    template class BounceBack<double, descriptors::D2Q9Descriptor>;
    template class NoDynamics<double, descriptors::D2Q9Descriptor>;

}
