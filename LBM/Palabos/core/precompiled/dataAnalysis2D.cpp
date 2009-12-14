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
#include "core/dataAnalysis2D.h"
#include "core/dataAnalysis2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

template double computeAverage(ScalarFieldBase2D<double>& scalarField, Box2D domain);
template double computeAverage(ScalarFieldBase2D<double>& scalarField);

template double computeMin(ScalarFieldBase2D<double>& scalarField, Box2D domain);
template double computeMin(ScalarFieldBase2D<double>& scalarField);

template double computeMax(ScalarFieldBase2D<double>& scalarField, Box2D domain);
template double computeMax(ScalarFieldBase2D<double>& scalarField);

template double computeBoundedAverage(ScalarFieldBase2D<double>& scalarField, Box2D domain);
template double computeBoundedAverage(ScalarFieldBase2D<double>& scalarField);

template double computeAverageDensity<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain );
template double computeAverageDensity<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice );

template double computeAverageRhoBar<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain );
template double computeAverageRhoBar<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice );

template double computeAverageEnergy<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, Box2D domain );
template double computeAverageEnergy<double, descriptors::D2Q9Descriptor> (
                    BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice );

}  // namespace plb
