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
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/dataAnalysisFunctionals2D.h"
#include <memory>

#ifndef ATOMIC_DATA_ANALYSIS_2D_H
#define ATOMIC_DATA_ANALYSIS_2D_H

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& density);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeDensity(BlockLattice2D<T,Descriptor>& lattice);


/* *************** RhoBar ******************************************** */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& rhoBar);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeRhoBar(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& energy);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeKineticEnergy(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& velocityNorm);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeVelocityNorm(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Velocity Component ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(BlockLattice2D<T,Descriptor>& lattice,
                              ScalarField2D<T>& velocityComponent,
                              plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeVelocityComponent(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(BlockLattice2D<T,Descriptor>& lattice,
                     TensorField2D<T,Descriptor<T>::d>& velocity);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,Descriptor<T>::d> >
    computeVelocity(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(BlockLattice2D<T,Descriptor>& lattice,
                             TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(BlockLattice2D<T,Descriptor>& lattice,
                                 TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& S);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(BlockLattice2D<T,Descriptor>& lattice);


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(BlockLattice2D<T,Descriptor>& lattice,
                       ScalarField2D<T>& population, plint iPop );

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> >
    computePopulation(BlockLattice2D<T,Descriptor>& lattice, plint iPop);


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Generic operations *************** */

template<typename T, class Function>
void apply(Function f, ScalarField2D<T>& field) {
    applyProcessingFunctional (
            new ApplyScalarFunctional2D<T,Function>(f), field.getBoundingBox(), field);
}

template<typename T, class Function>
void evaluate(Function f, ScalarField2D<T>& field, ScalarField2D<T>& result) {
    applyProcessingFunctional (
            new EvaluateScalarFunctional2D<T,Function>(f), field.getBoundingBox(), field, result);
}

template<typename T, class Function>
std::auto_ptr<ScalarField2D<T> > evaluate(Function f, ScalarField2D<T>& field) {
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    evaluate(f, field, result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


/* *************** ScalarField - Scalar operations *************** */

template<typename T>
void add(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(T scalar, ScalarField2D<T>& field);


template<typename T>
void add(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(ScalarField2D<T>& field, T scalar);


template<typename T>
void subtract(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(T scalar, ScalarField2D<T>& field);


template<typename T>
void subtract(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(ScalarField2D<T>& field, T scalar);


template<typename T>
void multiply(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(T scalar, ScalarField2D<T>& field);


template<typename T>
void multiply(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(ScalarField2D<T>& field, T scalar);


template<typename T>
void divide(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(T scalar, ScalarField2D<T>& field);


template<typename T>
void divide(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(ScalarField2D<T>& field, T scalar);



/* *************** ScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(ScalarField2D<T>& field, T scalar);

template<typename T>
void subtractInPlace(ScalarField2D<T>& field, T scalar);

template<typename T>
void multiplyInPlace(ScalarField2D<T>& field, T scalar);

template<typename T>
void divideInPlace(ScalarField2D<T>& field, T scalar);



/* *************** ScalarField - ScalarField operations *************** */

template<typename T>
void add(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(ScalarField2D<T>& A, ScalarField2D<T>& B);


template<typename T>
void subtract(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(ScalarField2D<T>& A, ScalarField2D<T>& B);


template<typename T>
void multiply(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(ScalarField2D<T>& A, ScalarField2D<T>& B);


template<typename T>
void divide(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result);

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(ScalarField2D<T>& A, ScalarField2D<T>& B);



/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B);

template<typename T>
void subtractInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B);

template<typename T>
void multiplyInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B);

template<typename T>
void divideInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& component, int iComponent);

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > extractComponent(TensorField2D<T,nDim>& tensorField, int iComponent);


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& norm);

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > computeNorm(TensorField2D<T,nDim>& norm);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& normSqr);

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > computeNormSqr(TensorField2D<T,nDim>& normSqr);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(TensorField2D<T,3>& tensorField, ScalarField2D<T>& norm);

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorNorm(TensorField2D<T,3>& norm);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(TensorField2D<T,3>& tensorField, ScalarField2D<T>& normSqr);

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorNormSqr(TensorField2D<T,3>& normSqr);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(TensorField2D<T,3>& tensorField, ScalarField2D<T>& trace);

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorTrace(TensorField2D<T,3>& trace);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(TensorField2D<T,2>& velocity, ScalarField2D<T>& vorticity);

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeVorticity(TensorField2D<T,2>& velocity);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(TensorField2D<T,2>& velocity, ScalarField2D<T>& vorticity);

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeBulkVorticity(TensorField2D<T,2>& velocity);


/* *************** Strain Rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(TensorField2D<T,2>& velocity, TensorField2D<T,3>& S);

template<typename T>
std::auto_ptr<TensorField2D<T,3> > computeStrainRate(TensorField2D<T,2>& velocity);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(TensorField2D<T,2>& velocity, TensorField2D<T,3>& S);

template<typename T>
std::auto_ptr<TensorField2D<T,3> > computeBulkStrainRate(TensorField2D<T,2>& velocity);


/* *************** TensorField - TensorField operations *************** */

template<typename T, int nDim>
void add(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > add(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);


template<typename T, int nDim>
void subtract(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > subtract(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);


template<typename T, int nDim>
void multiply(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > multiply(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);


template<typename T, int nDim>
void divide(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > divide(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);



/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);

template<typename T, int nDim>
void subtractInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);

template<typename T, int nDim>
void multiplyInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);

template<typename T, int nDim>
void divideInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);

}  // namespace plb

#endif  // ATOMIC_DATA_ANALYSIS_2D_H
