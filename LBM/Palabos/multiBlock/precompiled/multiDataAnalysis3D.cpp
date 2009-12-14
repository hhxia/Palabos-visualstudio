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

#include "multiBlock/multiDataAnalysis3D.h"
#include "multiBlock/multiDataAnalysis3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {
    
/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template void extractSubDomain<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& extractedLattice, Box3D domain );

template std::auto_ptr<MultiBlockLattice3D<double,descriptors::D3Q19Descriptor> >
    extractSubDomain<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );


template void computeDensity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& density, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeDensity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeDensity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeRhoBar<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& rhoBar, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeRhoBar<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeRhoBar<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeKineticEnergy<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& energy, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeKineticEnergy<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeKineticEnergy<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeVelocityNorm<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& velocityNorm, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeVelocityNorm<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> > computeVelocityNorm<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeVelocityComponent<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& velocityComponent, Box3D domain, plint iComponent );

template std::auto_ptr<MultiScalarField3D<double> > computeVelocityComponent<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain, plint iComponent );

template std::auto_ptr<MultiScalarField3D<double> > computeVelocityComponent<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, plint iComponent );


template void computeVelocity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiTensorField3D<double,3>& velocity, Box3D domain );

template std::auto_ptr<MultiTensorField3D<double,3> > computeVelocity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain );

template std::auto_ptr<MultiTensorField3D<double,3> > computeVelocity<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computePopulation<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& population, Box3D domain, plint iPop );

template std::auto_ptr<MultiScalarField3D<double> > computePopulation<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, Box3D domain, plint iPop );

template std::auto_ptr<MultiScalarField3D<double> > computePopulation<double, descriptors::D3Q19Descriptor> (
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, plint iPop );



/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

template void extractSubDomain<double> (
        MultiScalarField3D<double>& field,
        MultiScalarField3D<double>& density, Box3D domain );

template std::auto_ptr<MultiScalarField3D<double> >
    extractSubDomain<double> ( MultiScalarField3D<double>& field, Box3D domain );


/* *************** MultiScalarField - Scalar operations *************** */

template void add(double scalar, MultiScalarField3D<double>& field, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(double scalar, MultiScalarField3D<double>& field, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(double scalar, MultiScalarField3D<double>& field);

template void add(MultiScalarField3D<double>& field, double scalar, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(MultiScalarField3D<double>& field, double scalar);

template void subtract(double scalar, MultiScalarField3D<double>& field, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(double scalar, MultiScalarField3D<double>& field, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(double scalar, MultiScalarField3D<double>& field);

template void subtract(MultiScalarField3D<double>& field, double scalar, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(MultiScalarField3D<double>& field, double scalar);

template void multiply(double scalar, MultiScalarField3D<double>& field, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(double scalar, MultiScalarField3D<double>& field, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(double scalar, MultiScalarField3D<double>& field);

template void multiply(MultiScalarField3D<double>& field, double scalar, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(MultiScalarField3D<double>& field, double scalar);

template void divide(double scalar, MultiScalarField3D<double>& field, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(double scalar, MultiScalarField3D<double>& field, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(double scalar, MultiScalarField3D<double>& field);

template void divide(MultiScalarField3D<double>& field, double scalar, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(MultiScalarField3D<double>& field, double scalar);

/* *************** MultiScalarField - Scalar inplace operations *************** */

template void addInPlace(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template void addInPlace(MultiScalarField3D<double>& field, double scalar);
template void subtractInPlace(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template void subtractInPlace(MultiScalarField3D<double>& field, double scalar);
template void multiplyInPlace(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template void multiplyInPlace(MultiScalarField3D<double>& field, double scalar);
template void divideInPlace(MultiScalarField3D<double>& field, double scalar, Box3D domain);
template void divideInPlace(MultiScalarField3D<double>& field, double scalar);

/* *************** MultiScalarField - MultiScalarField operations *************** */

template void add(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > add(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void subtract(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > subtract(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void multiply(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > multiply(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void divide(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, MultiScalarField3D<double>& result, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template std::auto_ptr<MultiScalarField3D<double> > divide(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template void addInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template void addInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void subtractInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template void subtractInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void multiplyInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template void multiplyInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

template void divideInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B, Box3D domain);
template void divideInPlace(MultiScalarField3D<double>& A, MultiScalarField3D<double>& B);

                     
/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */


template void extractSubDomain<double,3> (
        MultiTensorField3D<double,3>& field,
        MultiTensorField3D<double,3>& extractedField, Box3D domain );

template std::auto_ptr<MultiTensorField3D<double,3> >
    extractSubDomain<double,3> ( MultiTensorField3D<double,3>& field, Box3D domain );

                      
template void extractComponent<double,3> (
        MultiTensorField3D<double,3>& tensorField, MultiScalarField3D<double>& component, Box3D domain, int iComponent );

template std::auto_ptr<MultiScalarField3D<double> > extractComponent<double,3> (
        MultiTensorField3D<double,3>& tensorField, Box3D domain, int iComponent);

template std::auto_ptr<MultiScalarField3D<double> > extractComponent<double,3> (
        MultiTensorField3D<double,3>& tensorField, int iComponent);


template void computeNorm<double,3> (
        MultiTensorField3D<double,3>& tensorField, MultiScalarField3D<double>& component, Box3D domain);

template std::auto_ptr<MultiScalarField3D<double> > computeNorm<double,3> (
        MultiTensorField3D<double,3>& tensorField, Box3D domain);

template std::auto_ptr<MultiScalarField3D<double> > computeNorm<double,3> (
        MultiTensorField3D<double,3>& tensorField);


template void computeNormSqr<double,3> (
        MultiTensorField3D<double,3>& tensorField, MultiScalarField3D<double>& component, Box3D domain);

template std::auto_ptr<MultiScalarField3D<double> > computeNormSqr<double,3> (
        MultiTensorField3D<double,3>& tensorField, Box3D domain);

template std::auto_ptr<MultiScalarField3D<double> > computeNormSqr<double,3> (
        MultiTensorField3D<double,3>& tensorField);


template void computeVorticity<double>(MultiTensorField3D<double,3>& velocity, MultiTensorField3D<double,3>& vorticity, Box3D domain);

template std::auto_ptr<MultiTensorField3D<double,3> > computeVorticity<double>(MultiTensorField3D<double,3>& velocity, Box3D domain);

template std::auto_ptr<MultiTensorField3D<double,3> > computeVorticity<double>(MultiTensorField3D<double,3>& velocity);


template void computeBulkVorticity<double>(MultiTensorField3D<double,3>& velocity, MultiTensorField3D<double,3>& vorticity, Box3D domain);

template std::auto_ptr<MultiTensorField3D<double,3> > computeBulkVorticity<double>(MultiTensorField3D<double,3>& velocity, Box3D domain);

template std::auto_ptr<MultiTensorField3D<double,3> > computeBulkVorticity<double>(MultiTensorField3D<double,3>& velocity);


/* *************** MultiTensorField - MultiTensorField operations *************** */

template void add(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, MultiTensorField3D<double,3>& result, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > add(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > add(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void subtract(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, MultiTensorField3D<double,3>& result, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > subtract(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > subtract(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void multiply(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, MultiTensorField3D<double,3>& result, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > multiply(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > multiply(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void divide(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, MultiTensorField3D<double,3>& result, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > divide(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template std::auto_ptr<MultiTensorField3D<double,3> > divide(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template void addInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template void addInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void subtractInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template void subtractInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void multiplyInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template void multiplyInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

template void divideInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B, Box3D domain);
template void divideInPlace(MultiTensorField3D<double,3>& A, MultiTensorField3D<double,3>& B);

}  // namespace plb
