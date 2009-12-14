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

#include "atomicBlock/atomicDataAnalysis2D.h"
#include "atomicBlock/atomicDataAnalysis2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {
    
/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template void computeDensity<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, ScalarField2D<double>& density);

template std::auto_ptr<ScalarField2D<double> > computeDensity<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeRhoBar<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, ScalarField2D<double>& density);

template std::auto_ptr<ScalarField2D<double> > computeRhoBar<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeKineticEnergy<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, ScalarField2D<double>& energy);

template std::auto_ptr<ScalarField2D<double> > computeKineticEnergy<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeVelocityNorm<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, ScalarField2D<double>& velocityNorm);

template std::auto_ptr<ScalarField2D<double> > computeVelocityNorm<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computeVelocityComponent<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& velocityComponent, plint iComponent);

template std::auto_ptr<ScalarField2D<double> > computeVelocityComponent<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, plint iComponent );


template void computeVelocity<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, TensorField2D<double,2>& velocity);

template std::auto_ptr<TensorField2D<double,2> > computeVelocity<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice );


template void computePopulation<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, ScalarField2D<double>& population, plint iPop);

template std::auto_ptr<ScalarField2D<double> > computePopulation<double, descriptors::D2Q9Descriptor> (
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice, plint iPop );


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** ScalarField - Scalar operations *************** */

template void add(double scalar, ScalarField2D<double>& field, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > add(double scalar, ScalarField2D<double>& field);
template void add(ScalarField2D<double>& field, double scalar, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > add(ScalarField2D<double>& field, double scalar);

template void subtract(double scalar, ScalarField2D<double>& field, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > subtract(double scalar, ScalarField2D<double>& field);

template void subtract(ScalarField2D<double>& field, double scalar, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > subtract(ScalarField2D<double>& field, double scalar);

template void multiply(double scalar, ScalarField2D<double>& field, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > multiply(double scalar, ScalarField2D<double>& field);

template void multiply(ScalarField2D<double>& field, double scalar, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > multiply(ScalarField2D<double>& field, double scalar);

template void divide(double scalar, ScalarField2D<double>& field, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > divide(double scalar, ScalarField2D<double>& field);

template void divide(ScalarField2D<double>& field, double scalar, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > divide(ScalarField2D<double>& field, double scalar);


/* *************** ScalarField - Scalar inplace operations *************** */

template void addInPlace(ScalarField2D<double>& field, double scalar);
template void subtractInPlace(ScalarField2D<double>& field, double scalar);
template void multiplyInPlace(ScalarField2D<double>& field, double scalar);
template void divideInPlace(ScalarField2D<double>& field, double scalar);


/* *************** ScalarField - ScalarField operations *************** */

template void add(ScalarField2D<double>& A, ScalarField2D<double>& B, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > add(ScalarField2D<double>& A, ScalarField2D<double>& B);

template void subtract(ScalarField2D<double>& A, ScalarField2D<double>& B, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > subtract(ScalarField2D<double>& A, ScalarField2D<double>& B);

template void multiply(ScalarField2D<double>& A, ScalarField2D<double>& B, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > multiply(ScalarField2D<double>& A, ScalarField2D<double>& B);

template void divide(ScalarField2D<double>& A, ScalarField2D<double>& B, ScalarField2D<double>& result);
template std::auto_ptr<ScalarField2D<double> > divide(ScalarField2D<double>& A, ScalarField2D<double>& B);


/* *************** ScalarField - ScalarField inplace operations *************** */

template void addInPlace(ScalarField2D<double>& A, ScalarField2D<double>& B);
template void subtractInPlace(ScalarField2D<double>& A, ScalarField2D<double>& B);
template void multiplyInPlace(ScalarField2D<double>& A, ScalarField2D<double>& B);
template void divideInPlace(ScalarField2D<double>& A, ScalarField2D<double>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template void extractComponent<double,2> (
        TensorField2D<double,2>& tensorField, ScalarField2D<double>& component, int iComponent );

template std::auto_ptr<ScalarField2D<double> > extractComponent<double,2> (
        TensorField2D<double,2>& tensorField, int iComponent);


template void computeNorm<double,2> (
        TensorField2D<double,2>& tensorField, ScalarField2D<double>& component);

template std::auto_ptr<ScalarField2D<double> > computeNorm<double,2> (
        TensorField2D<double,2>& tensorField);


template void computeNormSqr<double,2> (
        TensorField2D<double,2>& tensorField, ScalarField2D<double>& component);

template std::auto_ptr<ScalarField2D<double> > computeNormSqr<double,2> (
        TensorField2D<double,2>& tensorField);

template void computeVorticity<double>(TensorField2D<double,2>& velocity, ScalarField2D<double>& vorticity);
template std::auto_ptr<ScalarField2D<double> > computeVorticity<double>(TensorField2D<double,2>& velocity);

template void computeBulkVorticity<double>(TensorField2D<double,2>& velocity, ScalarField2D<double>& vorticity);
template std::auto_ptr<ScalarField2D<double> > computeBulkVorticity<double>(TensorField2D<double,2>& velocity);


/* *************** TensorField - TensorField operations *************** */

template void add(TensorField2D<double,2>& A, TensorField2D<double,2>& B, TensorField2D<double,2>& result);
template std::auto_ptr<TensorField2D<double,2> > add(TensorField2D<double,2>& A, TensorField2D<double,2>& B);

template void subtract(TensorField2D<double,2>& A, TensorField2D<double,2>& B, TensorField2D<double,2>& result);
template std::auto_ptr<TensorField2D<double,2> > subtract(TensorField2D<double,2>& A, TensorField2D<double,2>& B);

template void multiply(TensorField2D<double,2>& A, TensorField2D<double,2>& B, TensorField2D<double,2>& result);
template std::auto_ptr<TensorField2D<double,2> > multiply(TensorField2D<double,2>& A, TensorField2D<double,2>& B);

template void divide(TensorField2D<double,2>& A, TensorField2D<double,2>& B, TensorField2D<double,2>& result);
template std::auto_ptr<TensorField2D<double,2> > divide(TensorField2D<double,2>& A, TensorField2D<double,2>& B);


/* *************** TensorField - TensorField inplace operations *************** */

template void addInPlace(TensorField2D<double,2>& A, TensorField2D<double,2>& B);
template void subtractInPlace(TensorField2D<double,2>& A, TensorField2D<double,2>& B);
template void multiplyInPlace(TensorField2D<double,2>& A, TensorField2D<double,2>& B);
template void divideInPlace(TensorField2D<double,2>& A, TensorField2D<double,2>& B);


}  // namespace plb
