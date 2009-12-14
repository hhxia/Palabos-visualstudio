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
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "core/dataAnalysisFunctionals3D.h"
#include <memory>


#ifndef MULTI_DATA_ANALYSIS_3D_H
#define MULTI_DATA_ANALYSIS_3D_H

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Extract Sub-Lattice ******************************* */

template<typename T, template<typename U> class Descriptor>
void extractSubDomain( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiBlockLattice3D<T,Descriptor>& extractedLattice,
                       Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiBlockLattice3D<T,Descriptor> > extractSubDomain(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);


/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& density, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** RhoBar ******************************************** */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& energy, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityNorm, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityComponent,
                              Box3D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                               Box3D domain, plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, plint iComponent);


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiTensorField3D<T,Descriptor<T>::d>& velocity, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >
   computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >
    computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                             MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
   computeDeviatoricStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice,
                                 MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
   computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice);


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& population,
                       Box3D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                        Box3D domain, plint iPop);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, plint iPop);



/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */


/* *************** Generic operations *************** */

template<typename T, class Function>
void apply(Function f, MultiScalarField3D<T>& field, Box3D domain) {
    applyProcessingFunctional (
            new ApplyScalarFunctional3D<T,Function>(f), domain, field);
}

template<typename T, class Function>
void apply(Function f, MultiScalarField3D<T>& field) {
    apply(f, field, field.getBoundingBox());
}


template<typename T, class Function>
void evaluate(Function f, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain) {
    applyProcessingFunctional (
            new EvaluateScalarFunctional3D<T,Function>(f), domain, field, result);
}

template<typename T, class Function>
std::auto_ptr<MultiScalarField3D<T> > evaluate(Function f, MultiScalarField3D<T>& field, Box3D domain) {
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    evaluate(f, field, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T, class Function>
std::auto_ptr<MultiScalarField3D<T> > evaluate(Function f, MultiScalarField3D<T>& field) {
    return evaluate(f, field, field.getBoundingBox());
}


/* *************** Extract Sub-ScalarField *************************** */

template<typename T>
void extractSubDomain( MultiScalarField3D<T>& field,
                       MultiScalarField3D<T>& extractedField,
                       Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > extractSubDomain(MultiScalarField3D<T>& field, Box3D domain);


/* *************** MultiScalarField - Scalar operations *************** */

template<typename T>
void add(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field);


template<typename T>
void add(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar);


template<typename T>
void subtract(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field);


template<typename T>
void subtract(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar);


template<typename T>
void multiply(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field);


template<typename T>
void multiply(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar);


template<typename T>
void divide(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field);


template<typename T>
void divide(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar);


/* *************** MultiScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar);

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar);

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar);

template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain);

template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar);


/* *************** MultiScalarField - MultiScalarField operations *************** */

template<typename T>
void add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


/* *************** MultiScalarField - MultiScalarField inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain);

template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiTensorField3D<T,nDim>& field,
                       MultiTensorField3D<T,nDim>& extractedField,
                       Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > extractSubDomain(MultiTensorField3D<T,nDim>& field, Box3D domain);


/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& component, Box3D domain, int iComponent);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, Box3D domain, int iComponent);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, int iComponent);


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& norm, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& norm, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& trace, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, Box3D domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity);


/* *************** Strain Rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity);


/* *************** MultiTensorField - MultiTensorField operations *************** */

template<typename T, int nDim>
void add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


/* *************** MultiTensorField - MultiTensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);


template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain);

template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B);

}  // namespace plb

#endif  // MULTI_DATA_ANALYSIS_3D_H
