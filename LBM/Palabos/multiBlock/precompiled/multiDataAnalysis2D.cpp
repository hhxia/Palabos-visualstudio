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

#include "multiBlock/multiDataAnalysis2D.h"
#include "multiBlock/multiDataAnalysis2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {
    
/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template void extractSubDomain<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& density, Box2D domain );

template std::auto_ptr<MultiBlockLattice2D<double,descriptors::D2Q9Descriptor> >
    extractSubDomain<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template void computeDensity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& density, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeDensity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeDensity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeRhoBar<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& rhoBar, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeRhoBar<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeRhoBar<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeKineticEnergy<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& energy, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeKineticEnergy<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeKineticEnergy<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeVelocityNorm<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& velocityNorm, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeVelocityNorm<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> > computeVelocityNorm<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeVelocityComponent<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& velocityComponent, Box2D domain, plint iComponent );

template std::auto_ptr<MultiScalarField2D<double> > computeVelocityComponent<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain, plint iComponent );

template std::auto_ptr<MultiScalarField2D<double> > computeVelocityComponent<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, plint iComponent );


template void computeVelocity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiTensorField2D<double,2>& velocity, Box2D domain );

template std::auto_ptr<MultiTensorField2D<double,2> > computeVelocity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain );

template std::auto_ptr<MultiTensorField2D<double,2> > computeVelocity<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computePopulation<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        MultiScalarField2D<double>& population, Box2D domain, plint iPop );

template std::auto_ptr<MultiScalarField2D<double> > computePopulation<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, Box2D domain, plint iPop );

template std::auto_ptr<MultiScalarField2D<double> > computePopulation<double, descriptors::D2Q9Descriptor> (
        MultiBlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, plint iPop );


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

template void extractSubDomain<double> (
        MultiScalarField2D<double>& field,
        MultiScalarField2D<double>& density, Box2D domain );

template std::auto_ptr<MultiScalarField2D<double> >
    extractSubDomain<double> ( MultiScalarField2D<double>& field, Box2D domain );


/* *************** MultiScalarField - Scalar operations *************** */

template void add(double scalar, MultiScalarField2D<double>& field, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(double scalar, MultiScalarField2D<double>& field, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(double scalar, MultiScalarField2D<double>& field);

template void add(MultiScalarField2D<double>& field, double scalar, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(MultiScalarField2D<double>& field, double scalar);

template void subtract(double scalar, MultiScalarField2D<double>& field, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(double scalar, MultiScalarField2D<double>& field, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(double scalar, MultiScalarField2D<double>& field);

template void subtract(MultiScalarField2D<double>& field, double scalar, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(MultiScalarField2D<double>& field, double scalar);

template void multiply(double scalar, MultiScalarField2D<double>& field, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(double scalar, MultiScalarField2D<double>& field, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(double scalar, MultiScalarField2D<double>& field);

template void multiply(MultiScalarField2D<double>& field, double scalar, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(MultiScalarField2D<double>& field, double scalar);

template void divide(double scalar, MultiScalarField2D<double>& field, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(double scalar, MultiScalarField2D<double>& field, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(double scalar, MultiScalarField2D<double>& field);

template void divide(MultiScalarField2D<double>& field, double scalar, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(MultiScalarField2D<double>& field, double scalar);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template void addInPlace(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template void addInPlace(MultiScalarField2D<double>& field, double scalar);
template void subtractInPlace(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template void subtractInPlace(MultiScalarField2D<double>& field, double scalar);
template void multiplyInPlace(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template void multiplyInPlace(MultiScalarField2D<double>& field, double scalar);
template void divideInPlace(MultiScalarField2D<double>& field, double scalar, Box2D domain);
template void divideInPlace(MultiScalarField2D<double>& field, double scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template void add(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > add(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void subtract(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > subtract(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void multiply(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > multiply(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void divide(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, MultiScalarField2D<double>& result, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template std::auto_ptr<MultiScalarField2D<double> > divide(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template void addInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template void addInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void subtractInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template void subtractInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void multiplyInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template void multiplyInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);

template void divideInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B, Box2D domain);
template void divideInPlace(MultiScalarField2D<double>& A, MultiScalarField2D<double>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template void extractSubDomain<double,2> (
        MultiTensorField2D<double,2>& field,
        MultiTensorField2D<double,2>& density, Box2D domain );

template std::auto_ptr<MultiTensorField2D<double,2> >
    extractSubDomain<double,2> ( MultiTensorField2D<double,2>& field, Box2D domain );


template void extractComponent<double,2> (
        MultiTensorField2D<double,2>& tensorField, MultiScalarField2D<double>& component, Box2D domain, int iComponent );

template std::auto_ptr<MultiScalarField2D<double> > extractComponent<double,2> (
        MultiTensorField2D<double,2>& tensorField, Box2D domain, int iComponent);

template std::auto_ptr<MultiScalarField2D<double> > extractComponent<double,2> (
        MultiTensorField2D<double,2>& tensorField, int iComponent);


template void computeNorm<double,2> (
        MultiTensorField2D<double,2>& tensorField, MultiScalarField2D<double>& component, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeNorm<double,2> (
        MultiTensorField2D<double,2>& tensorField, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeNorm<double,2> (
        MultiTensorField2D<double,2>& tensorField);


template void computeNormSqr<double,2> (
        MultiTensorField2D<double,2>& tensorField, MultiScalarField2D<double>& component, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeNormSqr<double,2> (
        MultiTensorField2D<double,2>& tensorField, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeNormSqr<double,2> (
        MultiTensorField2D<double,2>& tensorField);


template void computeVorticity<double>(MultiTensorField2D<double,2>& velocity, MultiScalarField2D<double>& vorticity, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeVorticity<double>(MultiTensorField2D<double,2>& velocity, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeVorticity<double>(MultiTensorField2D<double,2>& velocity);


template void computeBulkVorticity<double>(MultiTensorField2D<double,2>& velocity, MultiScalarField2D<double>& vorticity, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeBulkVorticity<double>(MultiTensorField2D<double,2>& velocity, Box2D domain);

template std::auto_ptr<MultiScalarField2D<double> > computeBulkVorticity<double>(MultiTensorField2D<double,2>& velocity);


/* *************** MultiTensorField - MultiTensorField operations *************** */

template void add(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, MultiTensorField2D<double,2>& result, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > add(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > add(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void subtract(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, MultiTensorField2D<double,2>& result, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > subtract(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > subtract(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void multiply(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, MultiTensorField2D<double,2>& result, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > multiply(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > multiply(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void divide(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, MultiTensorField2D<double,2>& result, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > divide(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template std::auto_ptr<MultiTensorField2D<double,2> > divide(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template void addInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template void addInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void subtractInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template void subtractInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void multiplyInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template void multiplyInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);

template void divideInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B, Box2D domain);
template void divideInPlace(MultiTensorField2D<double,2>& A, MultiTensorField2D<double,2>& B);


}  // namespace plb
