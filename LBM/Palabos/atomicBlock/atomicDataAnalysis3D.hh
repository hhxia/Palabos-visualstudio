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
#include "atomicBlock/atomicDataAnalysis3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataCouplingWrapper3D.h"
#include "core/dataAnalysisFunctionals3D.h"


#ifndef ATOMIC_DATA_ANALYSIS_3D_HH
#define ATOMIC_DATA_ANALYSIS_3D_HH

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& density)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeDensity(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* density
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeDensity(lattice, *density);
    return std::auto_ptr<ScalarField3D<T> >(density);
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& rhoBar)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeRhoBar(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* rhoBar
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeRhoBar(lattice, *rhoBar);
    return std::auto_ptr<ScalarField3D<T> >(rhoBar);
}


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& energy)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeKineticEnergy(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* energy
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeKineticEnergy(lattice, *energy);
    return std::auto_ptr<ScalarField3D<T> >(energy);
}


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& velocityNorm)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityNorm(BlockLattice3D<T,Descriptor>& lattice)
{
    ScalarField3D<T>* velocityNorm
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityNorm(lattice, *velocityNorm);
    return std::auto_ptr<ScalarField3D<T> >(velocityNorm);
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(BlockLattice3D<T,Descriptor>& lattice,
                              ScalarField3D<T>& velocityComponent,
                              plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent),
            lattice.getBoundingBox(), lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computeVelocityComponent (
        BlockLattice3D<T,Descriptor>& lattice, plint iComponent )
{
    ScalarField3D<T>* velocityComponent
        = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocityComponent(lattice, *velocityComponent, iComponent);
    return std::auto_ptr<ScalarField3D<T> >(velocityComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(BlockLattice3D<T,Descriptor>& lattice,
                     TensorField3D<T,Descriptor<T>::d>& velocity)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,Descriptor<T>::d> > computeVelocity(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,Descriptor<T>::d>* velocity
        = new TensorField3D<T,Descriptor<T>::d>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeVelocity(lattice, *velocity);
    return std::auto_ptr<TensorField3D<T,Descriptor<T>::d> >(velocity);
}



/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress(BlockLattice3D<T,Descriptor>& lattice,
                             TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq)
{
    applyProcessingFunctional (
            new BoxDeviatoricStressFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeDeviatoricStress(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new TensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeDeviatoricStress(lattice, *PiNeq);
    return std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice,
                             TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S)
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional3D<T,Descriptor>, lattice.getBoundingBox(), lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRateFromStress(BlockLattice3D<T,Descriptor>& lattice)
{
    TensorField3D<T,SymmetricTensor<T,Descriptor>::n>* S
        = new TensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computeStrainRateFromStress(lattice, *S);
    return std::auto_ptr<TensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(S);
}



/* *************** Population *************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& population, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional3D<T,Descriptor>(iPop), lattice.getBoundingBox(), lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<ScalarField3D<T> > computePopulation(BlockLattice3D<T,Descriptor>& lattice, plint iPop)
{
    ScalarField3D<T>* population = new ScalarField3D<T>(lattice.getNx(), lattice.getNy(), lattice.getNz());
    computePopulation(lattice, *population, iPop);
    return std::auto_ptr<ScalarField3D<T> >(population);
}


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** ScalarField - Scalar operations *************** */

template<typename T>
void add(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void add(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    add(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    subtract(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    multiply(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(ScalarField3D<T>& field, T scalar, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& field, T scalar)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(field, scalar, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(T scalar, ScalarField3D<T>& field, ScalarField3D<T>& result)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional3D<T>(scalar), field.getBoundingBox(), field, result );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(T scalar, ScalarField3D<T>& field)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(field.getNx(), field.getNy(), field.getNz());
    divide(scalar, field, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}

/* *************** ScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void subtractInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void multiplyInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}

template<typename T>
void divideInPlace(ScalarField3D<T>& field, T scalar) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional3D<T>(scalar), field.getBoundingBox(), field);
}


/* *************** ScalarField - ScalarField operations *************** */

template<typename T>
void add(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > add(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void subtract(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > subtract(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void multiply(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > multiply(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


template<typename T>
void divide(ScalarField3D<T>& A, ScalarField3D<T>& B, ScalarField3D<T>& result)
{
    std::vector<ScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional3D<T>, A.getBoundingBox(), fields );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > divide(ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    ScalarField3D<T>* result = new ScalarField3D<T>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::auto_ptr<ScalarField3D<T> >(result);
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void subtractInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void multiplyInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}

template<typename T>
void divideInPlace(ScalarField3D<T>& A, ScalarField3D<T>& B) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional3D<T>, A.getBoundingBox(), A, B );
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional3D<T,nDim>(iComponent), tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > extractComponent(TensorField3D<T,nDim>& tensorField, int iComponent)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    extractComponent(tensorField, *component, iComponent);
    return std::auto_ptr<ScalarField3D<T> >(component);
}


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component)
{
    applyProcessingFunctional (
            new ComputeNormFunctional3D<T,nDim>, tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNorm(TensorField3D<T,nDim>& tensorField)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNorm(tensorField, *component);
    return std::auto_ptr<ScalarField3D<T> >(component);
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(TensorField3D<T,nDim>& tensorField, ScalarField3D<T>& component)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional3D<T,nDim>, tensorField.getBoundingBox(), component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<ScalarField3D<T> > computeNormSqr(TensorField3D<T,nDim>& tensorField)
{
    ScalarField3D<T>* component = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeNormSqr(tensorField, *component);
    return std::auto_ptr<ScalarField3D<T> >(component);
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(TensorField3D<T,6>& tensorField, ScalarField3D<T>& norm)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional3D<T>, tensorField.getBoundingBox(), norm, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNorm(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* norm = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNorm(tensorField, *norm);
    return std::auto_ptr<ScalarField3D<T> >(norm);
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(TensorField3D<T,6>& tensorField, ScalarField3D<T>& normSqr)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional3D<T>, tensorField.getBoundingBox(), normSqr, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorNormSqr(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* normSqr = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorNormSqr(tensorField, *normSqr);
    return std::auto_ptr<ScalarField3D<T> >(normSqr);
}


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField, ScalarField3D<T>& trace)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional3D<T>, tensorField.getBoundingBox(), trace, tensorField );
}

template<typename T>
std::auto_ptr<ScalarField3D<T> > computeSymmetricTensorTrace(TensorField3D<T,6>& tensorField)
{
    ScalarField3D<T>* trace = new ScalarField3D<T>(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    computeSymmetricTensorTrace(tensorField, *trace);
    return std::auto_ptr<ScalarField3D<T> >(trace);
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional3D<T,3>, velocity.getBoundingBox(), velocity, vorticity, envelopeWidth );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeVorticity(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,3>* vorticity = new TensorField3D<T,3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeVorticity(velocity, *vorticity);
    return std::auto_ptr<TensorField3D<T,3> >(vorticity);
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(TensorField3D<T,3>& velocity, TensorField3D<T,3>& vorticity)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional3D<T,3>, velocity.getBoundingBox(), velocity, vorticity );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeBulkVorticity(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,3>* vorticity = new TensorField3D<T,3>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkVorticity(velocity, *vorticity);
    return std::auto_ptr<TensorField3D<T,3> >(vorticity);
}


/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional3D<T,3>, velocity.getBoundingBox(), velocity, S, envelopeWidth );
}

template<typename T>
std::auto_ptr<TensorField3D<T,3> > computeStrainRate(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,6>* S = new TensorField3D<T,6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeStrainRate(velocity, *S);
    return std::auto_ptr<TensorField3D<T,6> >(S);
}


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(TensorField3D<T,3>& velocity, TensorField3D<T,6>& S)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional3D<T,3>, velocity.getBoundingBox(), velocity, S );
}

template<typename T>
std::auto_ptr<TensorField3D<T,6> > computeBulkStrainRate(TensorField3D<T,3>& velocity)
{
    TensorField3D<T,6>* S = new TensorField3D<T,6>(velocity.getNx(), velocity.getNy(), velocity.getNz());
    computeBulkStrainRate(velocity, *S);
    return std::auto_ptr<TensorField3D<T,6> >(S);
}


/* *************** TensorField - TensorField operations *************** */

template<typename T, int nDim>
void add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > add(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    add(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > subtract(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    subtract(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > multiply(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    multiply(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


template<typename T, int nDim>
void divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B, TensorField3D<T,nDim>& result)
{
    std::vector<TensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional3D<T,nDim>, A.getBoundingBox(), fields );
}

template<typename T, int nDim>
std::auto_ptr<TensorField3D<T,nDim> > divide(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    TensorField3D<T,nDim>* result = new TensorField3D<T,nDim>(A.getNx(), A.getNy(), A.getNz());
    divide(A, B, *result);
    return std::auto_ptr<TensorField3D<T,nDim> >(result);
}


/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void subtractInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void multiplyInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

template<typename T, int nDim>
void divideInPlace(TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>, A.getBoundingBox(), A, B );
}

}  // namespace plb

#endif  // ATOMIC_DATA_ANALYSIS_3D_HH
