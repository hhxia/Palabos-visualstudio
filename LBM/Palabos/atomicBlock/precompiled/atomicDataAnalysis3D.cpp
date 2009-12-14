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

#include "atomicBlock/atomicDataAnalysis3D.h"
#include "atomicBlock/atomicDataAnalysis3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {
    
/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template void computeDensity<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, ScalarField3D<double>& density);

template std::auto_ptr<ScalarField3D<double> > computeDensity<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeRhoBar<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, ScalarField3D<double>& density);

template std::auto_ptr<ScalarField3D<double> > computeRhoBar<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeKineticEnergy<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, ScalarField3D<double>& energy);

template std::auto_ptr<ScalarField3D<double> > computeKineticEnergy<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeVelocityNorm<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, ScalarField3D<double>& velocityNorm);

template std::auto_ptr<ScalarField3D<double> > computeVelocityNorm<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computeVelocityComponent<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& velocityComponent, plint iComponent);

template std::auto_ptr<ScalarField3D<double> > computeVelocityComponent<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, plint iComponent );


template void computeVelocity<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, TensorField3D<double,3>& velocity);

template std::auto_ptr<TensorField3D<double,3> > computeVelocity<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice );


template void computePopulation<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, ScalarField3D<double>& population, plint iPop);

template std::auto_ptr<ScalarField3D<double> > computePopulation<double, descriptors::D3Q19Descriptor> (
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice, plint iPop );


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** ScalarField - Scalar operations *************** */

template void add(double scalar, ScalarField3D<double>& field, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > add(double scalar, ScalarField3D<double>& field);
template void add(ScalarField3D<double>& field, double scalar, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > add(ScalarField3D<double>& field, double scalar);

template void subtract(double scalar, ScalarField3D<double>& field, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > subtract(double scalar, ScalarField3D<double>& field);

template void subtract(ScalarField3D<double>& field, double scalar, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > subtract(ScalarField3D<double>& field, double scalar);

template void multiply(double scalar, ScalarField3D<double>& field, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > multiply(double scalar, ScalarField3D<double>& field);

template void multiply(ScalarField3D<double>& field, double scalar, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > multiply(ScalarField3D<double>& field, double scalar);

template void divide(double scalar, ScalarField3D<double>& field, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > divide(double scalar, ScalarField3D<double>& field);

template void divide(ScalarField3D<double>& field, double scalar, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > divide(ScalarField3D<double>& field, double scalar);


/* *************** ScalarField - Scalar inplace operations *************** */

template void addInPlace(ScalarField3D<double>& field, double scalar);
template void subtractInPlace(ScalarField3D<double>& field, double scalar);
template void multiplyInPlace(ScalarField3D<double>& field, double scalar);
template void divideInPlace(ScalarField3D<double>& field, double scalar);


/* *************** ScalarField - ScalarField operations *************** */

template void add(ScalarField3D<double>& A, ScalarField3D<double>& B, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > add(ScalarField3D<double>& A, ScalarField3D<double>& B);

template void subtract(ScalarField3D<double>& A, ScalarField3D<double>& B, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > subtract(ScalarField3D<double>& A, ScalarField3D<double>& B);

template void multiply(ScalarField3D<double>& A, ScalarField3D<double>& B, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > multiply(ScalarField3D<double>& A, ScalarField3D<double>& B);

template void divide(ScalarField3D<double>& A, ScalarField3D<double>& B, ScalarField3D<double>& result);
template std::auto_ptr<ScalarField3D<double> > divide(ScalarField3D<double>& A, ScalarField3D<double>& B);


/* *************** ScalarField - ScalarField inplace operations *************** */

template void addInPlace(ScalarField3D<double>& A, ScalarField3D<double>& B);
template void subtractInPlace(ScalarField3D<double>& A, ScalarField3D<double>& B);
template void multiplyInPlace(ScalarField3D<double>& A, ScalarField3D<double>& B);
template void divideInPlace(ScalarField3D<double>& A, ScalarField3D<double>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */
                      
template void extractComponent<double,3> (
        TensorField3D<double,3>& tensorField, ScalarField3D<double>& component, int iComponent );

template std::auto_ptr<ScalarField3D<double> > extractComponent<double,3> (
        TensorField3D<double,3>& tensorField, int iComponent);


template void computeNorm<double,3> (
        TensorField3D<double,3>& tensorField, ScalarField3D<double>& component);

template std::auto_ptr<ScalarField3D<double> > computeNorm<double,3> (
        TensorField3D<double,3>& tensorField);


template void computeNormSqr<double,3> (
        TensorField3D<double,3>& tensorField, ScalarField3D<double>& component);

template std::auto_ptr<ScalarField3D<double> > computeNormSqr<double,3> (
        TensorField3D<double,3>& tensorField);

template void computeVorticity<double>(TensorField3D<double,3>& velocity, TensorField3D<double,3>& vorticity);
template std::auto_ptr<TensorField3D<double,3> > computeVorticity<double>(TensorField3D<double,3>& velocity);

template void computeBulkVorticity<double>(TensorField3D<double,3>& velocity, TensorField3D<double,3>& vorticity);
template std::auto_ptr<TensorField3D<double,3> > computeBulkVorticity<double>(TensorField3D<double,3>& velocity);


/* *************** TensorField - TensorField operations *************** */

template void add(TensorField3D<double,3>& A, TensorField3D<double,3>& B, TensorField3D<double,3>& result);
template std::auto_ptr<TensorField3D<double,3> > add(TensorField3D<double,3>& A, TensorField3D<double,3>& B);

template void subtract(TensorField3D<double,3>& A, TensorField3D<double,3>& B, TensorField3D<double,3>& result);
template std::auto_ptr<TensorField3D<double,3> > subtract(TensorField3D<double,3>& A, TensorField3D<double,3>& B);

template void multiply(TensorField3D<double,3>& A, TensorField3D<double,3>& B, TensorField3D<double,3>& result);
template std::auto_ptr<TensorField3D<double,3> > multiply(TensorField3D<double,3>& A, TensorField3D<double,3>& B);

template void divide(TensorField3D<double,3>& A, TensorField3D<double,3>& B, TensorField3D<double,3>& result);
template std::auto_ptr<TensorField3D<double,3> > divide(TensorField3D<double,3>& A, TensorField3D<double,3>& B);


/* *************** TensorField - TensorField inplace operations *************** */

template void addInPlace(TensorField3D<double,3>& A, TensorField3D<double,3>& B);
template void subtractInPlace(TensorField3D<double,3>& A, TensorField3D<double,3>& B);
template void multiplyInPlace(TensorField3D<double,3>& A, TensorField3D<double,3>& B);
template void divideInPlace(TensorField3D<double,3>& A, TensorField3D<double,3>& B);


}  // namespace plb
