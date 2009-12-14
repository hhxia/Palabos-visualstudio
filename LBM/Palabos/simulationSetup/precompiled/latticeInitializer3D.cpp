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
#include "simulationSetup/latticeInitializer3D.h"
#include "simulationSetup/latticeInitializer3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

    template void apply<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain,
            OneCellFunctional3D<double, descriptors::D3Q19Descriptor>* f );
    template void applyIndexed<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain,
            OneCellIndexedFunctional3D<double, descriptors::D3Q19Descriptor>* f );
    template void defineDynamics<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain,
            Dynamics<double,descriptors::D3Q19Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice,
            Box3D boundingBox, DomainFunctional3D* functional,
            Dynamics<double,descriptors::D3Q19Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice, DotList3D const& dotList,
            Dynamics<double,descriptors::D3Q19Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice, plint iX, plint iY, plint iZ,
            Dynamics<double,descriptors::D3Q19Descriptor>* dynamics);
    template void setBoundaryVelocity<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain, Array<double,3> velocity );
    template void setBoundaryDensity<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain, double rho );
    template void initializeAtEquilibrium<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain, double rho, Array<double,3> velocity );
    template void stripeOffDensityOffset<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain, double deltaRho );
    template void setCompositeDynamics<double, descriptors::D3Q19Descriptor> (
            BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, Box3D domain,
            CompositeDynamics<double, descriptors::D3Q19Descriptor>* compositeDynamics );

}  // namespace plb
