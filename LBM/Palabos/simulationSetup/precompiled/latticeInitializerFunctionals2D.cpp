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

#include "simulationSetup/latticeInitializerFunctionals2D.h"
#include "simulationSetup/latticeInitializerFunctionals2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

    template class OneCellFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class OneCellIndexedFunctional2D<double, descriptors::D2Q9Descriptor>;

    template class GenericLatticeFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class GenericIndexedLatticeFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class InstantiateDynamicsFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class InstantiateComplexDomainDynamicsFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class DynamicsFromMaskFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class InstantiateDotDynamicsFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class SetConstBoundaryVelocityFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class SetConstBoundaryDensityFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class IniConstEquilibriumFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class StripeOffDensityOffsetFunctional2D<double, descriptors::D2Q9Descriptor>;
    template class InstantiateCompositeDynamicsFunctional2D<double, descriptors::D2Q9Descriptor>;

}  // namespace plb
