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
#include "simulationSetup/latticeInitializer2D.h"
#include "simulationSetup/latticeInitializer2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

    template void apply<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain,
            OneCellFunctional2D<double, descriptors::D2Q9Descriptor>* f );
    template void applyIndexed<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain,
            OneCellIndexedFunctional2D<double, descriptors::D2Q9Descriptor>* f );
    template void defineDynamics<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain,
            Dynamics<double,descriptors::D2Q9Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice,
            Box2D boundingBox, DomainFunctional2D* functional,
            Dynamics<double,descriptors::D2Q9Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice, DotList2D const& dotList,
            Dynamics<double,descriptors::D2Q9Descriptor>* dynamics);
    template void defineDynamics<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice, plint iX, plint iY,
            Dynamics<double,descriptors::D2Q9Descriptor>* dynamics);
    template void setBoundaryVelocity<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain, Array<double,2> velocity );
    template void setBoundaryDensity<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain, double rho );
    template void initializeAtEquilibrium<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain, double rho, Array<double,2> velocity );
    template void stripeOffDensityOffset<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain, double deltaRho );
    template void setCompositeDynamics<double, descriptors::D2Q9Descriptor> (
            BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain,
            CompositeDynamics<double, descriptors::D2Q9Descriptor>* compositeDynamics );

}  // namespace plb
