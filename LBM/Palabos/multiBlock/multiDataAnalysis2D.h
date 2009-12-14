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
 * Helper functions for domain initialization -- header file.
 */
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "core/dataAnalysisFunctionals2D.h"
#include <memory>


#ifndef MULTI_DATA_ANALYSIS_2D_H
#define MULTI_DATA_ANALYSIS_2D_H

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Extract Sub-Lattice ******************************* */

template<typename T, template<typename U> class Descriptor>
void extractSubDomain( MultiBlockLattice2D<T,Descriptor>& lattice,
                       MultiBlockLattice2D<T,Descriptor>& extractedLattice,
                       Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiBlockLattice2D<T,Descriptor> > extractSubDomain(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);


/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& density, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** RhoBar ******************************************** */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& rhoBar, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& energy, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& velocityNorm, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& velocityComponent,
                              Box2D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice,
                                                               Box2D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice, plint iComponent);


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiTensorField2D<T,Descriptor<T>::d>& velocity, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,Descriptor<T>::d> >
   computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,Descriptor<T>::d> >
    computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                             MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
   computeDeviatoricStress(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice,
                                 MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>& S, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
   computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice);


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& population,
                       Box2D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice,
                                                        Box2D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice, plint iPop);


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Generic operations *************** */

template<typename T, class Function>
void apply(Function f, MultiScalarField2D<T>& field, Box2D domain) {
    applyProcessingFunctional (
            new ApplyScalarFunctional2D<T,Function>(f), domain, field);
}

template<typename T, class Function>
void apply(Function f, MultiScalarField2D<T>& field) {
    apply(f, field, field.getBoundingBox());
}


template<typename T, class Function>
void evaluate(Function f, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain) {
    applyProcessingFunctional (
            new EvaluateScalarFunctional2D<T,Function>(f), domain, field, result);
}

template<typename T, class Function>
std::auto_ptr<MultiScalarField2D<T> > evaluate(Function f, MultiScalarField2D<T>& field, Box2D domain) {
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T, class Function>
std::auto_ptr<MultiScalarField2D<T> > evaluate(Function f, MultiScalarField2D<T>& field) {
    return evaluate(f, field, field.getBoundingBox());
}


/* *************** Extract Sub-ScalarField *************************** */

template<typename T>
void extractSubDomain( MultiScalarField2D<T>& field,
                       MultiScalarField2D<T>& extractedField,
                       Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > extractSubDomain(MultiScalarField2D<T>& field, Box2D domain);


/* *************** MultiScalarField - Scalar operations *************** */

template<typename T>
void add(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T>& field, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T>& field);


template<typename T>
void add(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& field, T scalar);


template<typename T>
void subtract(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T>& field, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T>& field);


template<typename T>
void subtract(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& field, T scalar);


template<typename T>
void multiply(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T>& field, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T>& field);


template<typename T>
void multiply(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& field, T scalar);


template<typename T>
void divide(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T>& field, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T>& field);


template<typename T>
void divide(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& field, T scalar);


/* *************** MultiScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
void addInPlace(MultiScalarField2D<T>& field, T scalar);

template<typename T>
void subtractInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
void subtractInPlace(MultiScalarField2D<T>& field, T scalar);

template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& field, T scalar);

template<typename T>
void divideInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain);

template<typename T>
void divideInPlace(MultiScalarField2D<T>& field, T scalar);


/* *************** MultiScalarField - MultiScalarField operations *************** */

template<typename T>
void add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
void addInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void subtractInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
void subtractInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


template<typename T>
void divideInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain);

template<typename T>
void divideInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiTensorField2D<T,nDim>& field,
                       MultiTensorField2D<T,nDim>& extractedField,
                       Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > extractSubDomain(MultiTensorField2D<T,nDim>& field, Box2D domain);


/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& component, Box2D domain, int iComponent);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > extractComponent(MultiTensorField2D<T,nDim>& tensorField, Box2D domain, int iComponent);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > extractComponent(MultiTensorField2D<T,nDim>& tensorField, int iComponent);


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& norm, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T,nDim>& tensorField, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T,nDim>& tensorField);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& normSqr, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T,nDim>& tensorField, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T,nDim>& tensorField);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& norm, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& normSqr, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& trace, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiTensorField2D<T,2>& velocity, MultiScalarField2D<T>& vorticity, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T,2>& velocity, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T,2>& velocity);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiTensorField2D<T,2>& velocity, MultiScalarField2D<T>& vorticity, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T,2>& velocity, Box2D domain);

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T,2>& velocity);


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiTensorField2D<T,2>& velocity, MultiTensorField2D<T,3>& S, Box2D domain);

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeStrainRate(MultiTensorField2D<T,2>& velocity, Box2D domain);

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeStrainRate(MultiTensorField2D<T,2>& velocity);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiTensorField2D<T,2>& velocity, MultiTensorField2D<T,3>& S, Box2D domain);

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeBulkStrainRate(MultiTensorField2D<T,2>& velocity, Box2D domain);

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeBulkStrainRate(MultiTensorField2D<T,2>& velocity);


/* *************** MultiTensorField - MultiTensorField operations *************** */

template<typename T, int nDim>
void add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
void addInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);


template<typename T, int nDim>
void divideInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain);

template<typename T, int nDim>
void divideInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B);

}  // namespace plb

#endif  // MULTI_DATA_ANALYSIS_2D_H
