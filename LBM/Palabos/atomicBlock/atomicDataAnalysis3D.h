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
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/dataAnalysisFunctionals3D.h"
#include <memory>

#ifndef ATOMIC_DATA_ANALYSIS_3D_H
#define ATOMIC_DATA_ANALYSIS_3D_H

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& density);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeDensity(BlockLattice3D<T,Descriptor>& lattice);


/* *************** RhoBar ******************************************** */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& rhoBar);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeRhoBar(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& energy);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& velocityNorm);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Velocity Component ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(BlockLattice3D<T,Descriptor>& lattice,
                              ScalarField3D<T>& velocityComponent,
                              plint iComponent);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityComponent(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(BlockLattice3D<T,Descriptor>& lattice,
                     TensorField3D<T,Descriptor<T>::d>& velocity);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,Descriptor<T>::d> >
    computeVelocity(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(BlockLattice3D<T,Descriptor>& lattice,
                             TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(BlockLattice3D<T,Descriptor>& lattice);


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice,
                                 TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S);

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice);



/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(BlockLattice3D<T,Descriptor>& lattice,
                       ScalarField3D<T>& population, plint iPop );

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> >
    computePopulation(BlockLattice3D<T,Descriptor>& lattice, plint iPop);


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Generic operations *************** */

template<typename T, class Function>
void apply(Function f, ScalarField3D<T>& field) {
    applyProcessingFunctional (
            new ApplyScalarFunctional3D<T,Function>(f), field.getBoundingBox(), field);
}

template<typename T, class Function>
void evaluate(Function f, ScalarField3D<T>& field, ScalarField3D<T>& result) {
    applyProcessingFunctional (
            new EvaluateScalarFunctional3D<T,Function>(f), field.getBoundingBox(), field, result);
}

template<typename T, class Function>
std::auto_ptr<ScalarField3D<T> > evaluate(Function f, ScalarField3D<T>& field) {
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy());
    evaluate(f, field, result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


/* *************** ScalarField - Scalar operations *************** */

template<typename T>
void add(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(T scalar, ScalarField3D<T>& field);


template<typename T>
void add(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& field, T scalar);


template<typename T>
void subtract(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(T scalar, ScalarField3D<T>& field);


template<typename T>
void subtract(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& field, T scalar);


template<typename T>
void multiply(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(T scalar, ScalarField3D<T>& field);


template<typename T>
void multiply(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& field, T scalar);


template<typename T>
void divide(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(T scalar, ScalarField3D<T>& field);


template<typename T>
void divide(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& field, T scalar);



/* *************** ScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& field, T scalar);

template<typename T>
void subtractInPlace(ScalarField3D<T>& field, T scalar);

template<typename T>
void multiplyInPlace(ScalarField3D<T>& field, T scalar);

template<typename T>
void divideInPlace(ScalarField3D<T>& field, T scalar);



/* *************** ScalarField - ScalarField operations *************** */

template<typename T>
void add(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& A, ScalarField3D<T>& B);


template<typename T>
void subtract(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& A, ScalarField3D<T>& B);


template<typename T>
void multiply(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& A, ScalarField3D<T>& B);


template<typename T>
void divide(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result);

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& A, ScalarField3D<T>& B);



/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B);

template<typename T>
void subtractInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B);

template<typename T>
void multiplyInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B);

template<typename T>
void divideInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B);


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component, int iComponent);

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > extractComponent(TensorField3D<T,nDim>& tensorField, int iComponent);


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& norm);

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNorm(TensorField3D<T,nDim>& norm);


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& norm);

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNormSqr(TensorField3D<T,nDim>& norm);


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(TensorField3D<T,6>& tensorField, ScalarField3D<T>& norm);

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNorm(TensorField3D<T,6>& norm);


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(TensorField3D<T,6>& tensorField, ScalarField3D<T>& normSqr);

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNormSqr(TensorField3D<T,6>& normSqr);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField, ScalarField3D<T>& trace);

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T,6>& trace);


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField, ScalarField3D<T>& trace);

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T,6>& trace);


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity);

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeVorticity(TensorField3D<T,3>& velocity);


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity);

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeBulkVorticity(TensorField3D<T,3>& velocity);


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S);

template<typename T>
std::auto_ptr<TensorField3D<T,6> > computeStrainRate(TensorField3D<T,3>& velocity);


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S);

template<typename T>
std::auto_ptr<TensorField3D<T,6> > computeBulkStrainRate(TensorField3D<T,3>& velocity);


/* *************** TensorField - TensorField operations *************** */

template<typename T, int nDim>
void add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);


template<typename T, int nDim>
void subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);


template<typename T, int nDim>
void multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);


template<typename T, int nDim>
void divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result);

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);



/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);

template<typename T, int nDim>
void subtractInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);

template<typename T, int nDim>
void multiplyInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);

template<typename T, int nDim>
void divideInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);

}  // namespace plb

#endif  // ATOMIC_DATA_ANALYSIS_3D_H
