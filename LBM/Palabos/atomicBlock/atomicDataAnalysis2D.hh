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
#include "atomicBlock/atomicDataAnalysis2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataCouplingWrapper2D.h"
#include "core/dataAnalysisFunctionals2D.h"


#ifndef ATOMIC_DATA_ANALYSIS_2D_HH
#define ATOMIC_DATA_ANALYSIS_2D_HH

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& density)
{
    applyProcessingFunctional (
            new BoxDensityFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeDensity(BlockLattice2D<T,Descriptor>& lattice)
{
    ScalarField2D<T>* density = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeDensity(lattice, *density);
    return std::auto_ptr<ScalarField2D<T> >(density);
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& rhoBar)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeRhoBar(BlockLattice2D<T,Descriptor>& lattice)
{
    ScalarField2D<T>* rhoBar = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeRhoBar(lattice, *rhoBar);
    return std::auto_ptr<ScalarField2D<T> >(rhoBar);
}


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& energy)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeKineticEnergy(BlockLattice2D<T,Descriptor>& lattice)
{
    ScalarField2D<T>* energy = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeKineticEnergy(lattice, *energy);
    return std::auto_ptr<ScalarField2D<T> >(energy);
}


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& velocityNorm)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeVelocityNorm(BlockLattice2D<T,Descriptor>& lattice)
{
    ScalarField2D<T>* velocityNorm = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeVelocityNorm(lattice, *velocityNorm);
    return std::auto_ptr<ScalarField2D<T> >(velocityNorm);
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(BlockLattice2D<T,Descriptor>& lattice,
                              ScalarField2D<T>& velocityComponent,
                              plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional2D<T,Descriptor>(iComponent),
            lattice.getBoundingBox(), lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computeVelocityComponent (
        BlockLattice2D<T,Descriptor>& lattice, plint iComponent )
{
    ScalarField2D<T>* velocityComponent = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computeVelocityComponent(lattice, *velocityComponent, iComponent);
    return std::auto_ptr<ScalarField2D<T> >(velocityComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(BlockLattice2D<T,Descriptor>& lattice,
                     TensorField2D<T,Descriptor<T>::d>& velocity)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,Descriptor<T>::d> > computeVelocity(BlockLattice2D<T,Descriptor>& lattice)
{
    TensorField2D<T,Descriptor<T>::d>* velocity
        = new TensorField2D<T,Descriptor<T>::d>(lattice.getNx(), lattice.getNy());
    computeVelocity(lattice, *velocity);
    return std::auto_ptr<TensorField2D<T,Descriptor<T>::d> >(velocity);
}


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(BlockLattice2D<T,Descriptor>& lattice,
                             TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq)
{
    applyProcessingFunctional (
            new BoxDeviatoricStressFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeDeviatoricStress(BlockLattice2D<T,Descriptor>& lattice)
{
    TensorField2D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new TensorField2D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy());
    computeDeviatoricStress(lattice, *PiNeq);
    return std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(BlockLattice2D<T,Descriptor>& lattice,
                             TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& S)
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional2D<T,Descriptor>, lattice.getBoundingBox(), lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRateFromStress(BlockLattice2D<T,Descriptor>& lattice)
{
    TensorField2D<T,SymmetricTensor<T,Descriptor>::n>* S
        = new TensorField2D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy());
    computeStrainRateFromStress(lattice, *S);
    return std::auto_ptr<TensorField2D<T,SymmetricTensor<T,Descriptor>::n> >(S);
}


/* *************** Population *************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& population, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional2D<T,Descriptor>(iPop), lattice.getBoundingBox(), lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField2D<T> > computePopulation(BlockLattice2D<T,Descriptor>& lattice, plint iPop)
{
    ScalarField2D<T>* population = new ScalarField2D<T>(lattice.getNx(), lattice.getNy());
    computePopulation(lattice, *population, iPop);
    return std::auto_ptr<ScalarField2D<T> >(population);
}


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** ScalarField - Scalar operations *************** */

template<typename T>
void add(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(ScalarField2D<T>& field, T scalar)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    add(field, scalar, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void add(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(T scalar, ScalarField2D<T>& field)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    add(scalar, field, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void subtract(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(ScalarField2D<T>& field, T scalar)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    subtract(field, scalar, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void subtract(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(T scalar, ScalarField2D<T>& field)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    subtract(scalar, field, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void multiply(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(ScalarField2D<T>& field, T scalar)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    multiply(field, scalar, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void multiply(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(T scalar, ScalarField2D<T>& field)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    multiply(scalar, field, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void divide(ScalarField2D<T>& field, T scalar, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(ScalarField2D<T>& field, T scalar)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    divide(field, scalar, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void divide(T scalar, ScalarField2D<T>& field, ScalarField2D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional2D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(T scalar, ScalarField2D<T>& field)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(field.getNx(), field.getNy());
    divide(scalar, field, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}

/* *************** ScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(ScalarField2D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void subtractInPlace(ScalarField2D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void multiplyInPlace(ScalarField2D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void divideInPlace(ScalarField2D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional2D<T>(scalar), field.getBoundingBox(), field);
}


/* *************** ScalarField - ScalarField operations *************** */

template<typename T>
void add(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result)
{
    std::vector<ScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional2D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > add(ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(A.getNx(), A.getNy());
    add(A, B, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void subtract(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result)
{
    std::vector<ScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional2D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > subtract(ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(A.getNx(), A.getNy());
    subtract(A, B, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void multiply(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result)
{
    std::vector<ScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional2D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > multiply(ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(A.getNx(), A.getNy());
    multiply(A, B, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


template<typename T>
void divide(ScalarField2D<T>& A, ScalarField2D<T>& B, ScalarField2D<T>& result)
{
    std::vector<ScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional2D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > divide(ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    ScalarField2D<T>* result = new ScalarField2D<T>(A.getNx(), A.getNy());
    divide(A, B, *result);
    return std::auto_ptr<ScalarField2D<T> >(result);
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional2D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void subtractInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional2D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void multiplyInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional2D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void divideInPlace(ScalarField2D<T>& A, ScalarField2D<T>& B) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional2D<T>, A.getBoundingBox(), A, B );
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& component, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional2D<T,nDim>(iComponent), tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > extractComponent(TensorField2D<T,nDim>& tensorField, int iComponent)
{
    ScalarField2D<T>* component = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    extractComponent(tensorField, *component, iComponent);
    return std::auto_ptr<ScalarField2D<T> >(component);
}

/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& norm)
{
    applyProcessingFunctional (
            new ComputeNormFunctional2D<T,nDim>, tensorField.getBoundingBox(), norm, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > computeNorm(TensorField2D<T,nDim>& tensorField)
{
    ScalarField2D<T>* norm = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeNorm(tensorField, *norm);
    return std::auto_ptr<ScalarField2D<T> >(norm);
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(TensorField2D<T,nDim>& tensorField, ScalarField2D<T>& normSqr)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional2D<T,nDim>, tensorField.getBoundingBox(), normSqr, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField2D<T> > computeNormSqr(TensorField2D<T,nDim>& tensorField)
{
    ScalarField2D<T>* normSqr = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeNormSqr(tensorField, *normSqr);
    return std::auto_ptr<ScalarField2D<T> >(normSqr);
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(TensorField2D<T,3>& tensorField, ScalarField2D<T>& norm)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional2D<T>, tensorField.getBoundingBox(), norm, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorNorm(TensorField2D<T,3>& tensorField)
{
    ScalarField2D<T>* norm = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorNorm(tensorField, *norm);
    return std::auto_ptr<ScalarField2D<T> >(norm);
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(TensorField2D<T,3>& tensorField, ScalarField2D<T>& normSqr)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional2D<T>, tensorField.getBoundingBox(), normSqr, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorNormSqr(TensorField2D<T,3>& tensorField)
{
    ScalarField2D<T>* normSqr = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorNormSqr(tensorField, *normSqr);
    return std::auto_ptr<ScalarField2D<T> >(normSqr);
}


/* *************** Trace of each symmetric tensor of a field ************* */

template<typename T>
void computeSymmetricTensorTrace(TensorField2D<T,3>& tensorField, ScalarField2D<T>& trace)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional2D<T>, tensorField.getBoundingBox(), trace, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeSymmetricTensorTrace(TensorField2D<T,3>& tensorField)
{
    ScalarField2D<T>* trace = new ScalarField2D<T>(tensorField.getNx(), tensorField.getNy());
    computeSymmetricTensorTrace(tensorField, *trace);
    return std::auto_ptr<ScalarField2D<T> >(trace);
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(TensorField2D<T,2>& velocity, ScalarField2D<T>& vorticity)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional2D<T,2>, velocity.getBoundingBox(), vorticity, velocity, envelopeWidth );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeVorticity(TensorField2D<T,2>& velocity)
{
    ScalarField2D<T>* vorticity = new ScalarField2D<T>(velocity.getNx(), velocity.getNy());
    computeVorticity(velocity, *vorticity);
    return std::auto_ptr<ScalarField2D<T> >(vorticity);
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(TensorField2D<T,2>& velocity, ScalarField2D<T>& vorticity)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional2D<T,2>, velocity.getBoundingBox(), vorticity, velocity );
}

template<typename T>
std::auto_ptr<ScalarField2D<T> > computeBulkVorticity(TensorField2D<T,2>& velocity)
{
    ScalarField2D<T>* vorticity = new ScalarField2D<T>(velocity.getNx(), velocity.getNy());
    computeBulkVorticity(velocity, *vorticity);
    return std::auto_ptr<ScalarField2D<T> >(vorticity);
}


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(TensorField2D<T,2>& velocity, TensorField2D<T,3>& S)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional2D<T,2>, velocity.getBoundingBox(), velocity, S, envelopeWidth );
}

template<typename T>
std::auto_ptr<TensorField2D<T,3> > computeStrainRate(TensorField2D<T,2>& velocity)
{
    TensorField2D<T,3>* S = new TensorField2D<T,3>(velocity.getNx(), velocity.getNy());
    computeStrainRate(velocity, *S);
    return std::auto_ptr<TensorField2D<T,3> >(S);
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(TensorField2D<T,2>& velocity, TensorField2D<T,3>& S)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional2D<T,2>, velocity.getBoundingBox(), velocity, S );
}

template<typename T>
std::auto_ptr<TensorField2D<T,3> > computeBulkVorticity(TensorField2D<T,2>& velocity)
{
    TensorField2D<T,3>* S = new TensorField2D<T,3>(velocity.getNx(), velocity.getNy());
    computeBulkStrainRate(velocity, *S);
    return std::auto_ptr<TensorField2D<T,3> >(S);
}


/* *************** TensorField - TensorField operations *************** */

template<typename T, int nDim>
void add(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result)
{
    std::vector<TensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional2D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > add(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    TensorField2D<T,nDim>* result = new TensorField2D<T,nDim>(A.getNx(), A.getNy());
    add(A, B, *result);
    return std::auto_ptr<TensorField2D<T,nDim> >(result);
}


template<typename T, int nDim>
void subtract(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result)
{
    std::vector<TensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional2D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > subtract(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    TensorField2D<T,nDim>* result = new TensorField2D<T,nDim>(A.getNx(), A.getNy());
    subtract(A, B, *result);
    return std::auto_ptr<TensorField2D<T,nDim> >(result);
}


template<typename T, int nDim>
void multiply(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result)
{
    std::vector<TensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional2D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > multiply(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    TensorField2D<T,nDim>* result = new TensorField2D<T,nDim>(A.getNx(), A.getNy());
    multiply(A, B, *result);
    return std::auto_ptr<TensorField2D<T,nDim> >(result);
}


template<typename T, int nDim>
void divide(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B, TensorField2D<T,nDim>& result)
{
    std::vector<TensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional2D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField2D<T,nDim> > divide(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    TensorField2D<T,nDim>* result = new TensorField2D<T,nDim>(A.getNx(), A.getNy());
    divide(A, B, *result);
    return std::auto_ptr<TensorField2D<T,nDim> >(result);
}


/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional2D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void subtractInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional2D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void multiplyInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional2D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void divideInPlace(TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>, A.getBoundingBox(), A, B );
}

}  // namespace plb

#endif  // ATOMIC_DATA_ANALYSIS_2D_HH
