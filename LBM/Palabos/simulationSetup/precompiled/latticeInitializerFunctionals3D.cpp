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

#include "simulationSetup/latticeInitializerFunctionals3D.h"
#include "simulationSetup/latticeInitializerFunctionals3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

    template class OneCellFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class OneCellIndexedFunctional3D<double, descriptors::D3Q19Descriptor>;

    template class GenericLatticeFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class GenericIndexedLatticeFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class InstantiateDynamicsFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class InstantiateComplexDomainDynamicsFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class DynamicsFromMaskFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class InstantiateDotDynamicsFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class SetConstBoundaryVelocityFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class SetConstBoundaryDensityFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class IniConstEquilibriumFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class StripeOffDensityOffsetFunctional3D<double, descriptors::D3Q19Descriptor>;
    template class InstantiateCompositeDynamicsFunctional3D<double, descriptors::D3Q19Descriptor>;

}  // namespace plb
