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
#include "multiBlock/multiDataAnalysis3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiDataCouplingWrapper3D.h"
#include "core/dataAnalysisFunctionals3D.h"

#ifndef MULTI_DATA_ANALYSIS_3D_HH
#define MULTI_DATA_ANALYSIS_3D_HH

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Extract Sub-Lattice ******************************* */

template<typename T, template<typename U> class Descriptor>
void extractSubDomain( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiBlockLattice3D<T,Descriptor>& extractedLattice,
                       Box3D domain)
{
    applyProcessingFunctional (
            new ExtractLatticeSubDomainFunctional3D<T,Descriptor>, domain, lattice, extractedLattice );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiBlockLattice3D<T,Descriptor> > extractSubDomain(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiBlockLattice3D<T,Descriptor>* extractedLattice = new MultiBlockLattice3D<T,Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return std::auto_ptr<MultiBlockLattice3D<T,Descriptor> >(extractedLattice);
}


/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& density, Box3D domain)
{
    applyProcessingFunctional (
            new BoxDensityFunctional3D<T,Descriptor>, domain, lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiScalarField3D<T>* density = new MultiScalarField3D<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(density);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeDensity(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeDensity(lattice, lattice.getBoundingBox());
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar, Box3D domain)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional3D<T,Descriptor>, domain, lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiScalarField3D<T>* rhoBar = new MultiScalarField3D<T>(lattice, domain);
    computeRhoBar(lattice, *rhoBar, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeRhoBar(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeRhoBar(lattice, lattice.getBoundingBox());
}


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& energy, Box3D domain)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional3D<T,Descriptor>, domain, lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiScalarField3D<T>* energy = new MultiScalarField3D<T>(lattice, domain);
    computeKineticEnergy(lattice, *energy, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(energy);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeKineticEnergy(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityNorm, Box3D domain)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional3D<T,Descriptor>, domain, lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiScalarField3D<T>* velocityNorm = new MultiScalarField3D<T>(lattice, domain);
    computeVelocityNorm(lattice, *velocityNorm, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(velocityNorm);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(MultiBlockLattice3D<T,Descriptor>& lattice) {
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& velocityComponent,
                              Box3D domain, plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent), domain, lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                               Box3D domain, plint iComponent)
{
    MultiScalarField3D<T>* velocityComponent = new MultiScalarField3D<T>(lattice, domain);
    computeVelocityComponent(lattice, *velocityComponent, domain, iComponent);
    return std::auto_ptr<MultiScalarField3D<T> >(velocityComponent);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(MultiBlockLattice3D<T,Descriptor>& lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiTensorField3D<T,Descriptor<T>::d>& velocity, Box3D domain)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> > computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiTensorField3D<T,Descriptor<T>::d>* velocity
        = new MultiTensorField3D<T,Descriptor<T>::d>(lattice, domain);
    computeVelocity(lattice, *velocity, domain);
    return std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >(velocity);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,Descriptor<T>::d> >
    computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq,
                              Box3D domain )
{
    applyProcessingFunctional (
            new BoxDeviatoricStressFunctional3D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeDeviatoricStress(lattice, *PiNeq, domain);
    return std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeDeviatoricStress(lattice, lattice.getBoundingBox());
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S,
                                  Box3D domain )
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional3D<T,Descriptor>, domain, lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain)
{
    MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>* S
        = new MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeStrainRateFromStress(lattice, *S, domain);
    return std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >(S);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice3D<T,Descriptor>& lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& population,
                       Box3D domain, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional3D<T,Descriptor>(iPop), domain, lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice,
                                                        Box3D domain, plint iPop)
{
    MultiScalarField3D<T>* population = new MultiScalarField3D<T>(lattice, domain);
    computePopulation(lattice, *population, domain, iPop);
    return std::auto_ptr<MultiScalarField3D<T> >(population);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computePopulation(MultiBlockLattice3D<T,Descriptor>& lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-ScalarField *************************** */

template<typename T>
void extractSubDomain( MultiScalarField3D<T>& field,
                       MultiScalarField3D<T>& extractedField,
                       Box3D domain)
{
    applyProcessingFunctional (
            new ExtractScalarSubDomainFunctional3D<T>, domain, field, extractedField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > extractSubDomain(MultiScalarField3D<T>& field, Box3D domain)
{
    MultiScalarField3D<T>* extractedField = new MultiScalarField3D<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(extractedField);
}

/* *************** MultiScalarField - Scalar operations *************** */

template<typename T>
void add(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    add(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}


template<typename T>
void add(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    add(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(T scalar, MultiScalarField3D<T>& field)
{
    return add(scalar, field, field.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtract(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(T scalar, MultiScalarField3D<T>& field)
{
    return subtract(scalar, field, field.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiply(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(T scalar, MultiScalarField3D<T>& field)
{
    return multiply(scalar, field, field.getBoundingBox());
}


template<typename T>
void divide(MultiScalarField3D<T>& field, T scalar, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    divide(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}


template<typename T>
void divide(T scalar, MultiScalarField3D<T>& field, MultiScalarField3D<T>& result, Box3D domain)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional3D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(field, domain);
    divide(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(T scalar, MultiScalarField3D<T>& field)
{
    return divide(scalar, field, field.getBoundingBox());
}


/* *************** MultiScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void addInPlace(MultiScalarField3D<T>& field, T scalar) {
    addInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& field, T scalar) {
    subtractInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& field, T scalar) {
    multiplyInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar, Box3D domain) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional3D<T>(scalar), domain, field);
}

template<typename T>
void divideInPlace(MultiScalarField3D<T>& field, T scalar) {
    divideInPlace(field, scalar, field.getBoundingBox());
}


/* *************** MultiScalarField - MultiScalarField operations *************** */

template<typename T>
void add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(A, domain);
    add(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > add(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(A, domain);
    subtract(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > subtract(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(A, domain);
    multiply(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > multiply(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return multiply(A, B, A.getBoundingBox());
}

template<typename T>
void divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<T>& result, Box3D domain)
{
    std::vector<MultiScalarField3D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional3D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain)
{
    MultiScalarField3D<T>* result = new MultiScalarField3D<T>(A, domain);
    divide(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > divide(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B)
{
    return divide(A, B, A.getBoundingBox());
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void addInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void subtractInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void multiplyInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, Box3D domain) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional3D<T>, domain, A, B );
}

template<typename T>
void divideInPlace(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiTensorField3D<T,nDim>& field,
                       MultiTensorField3D<T,nDim>& extractedField,
                       Box3D domain)
{
    applyProcessingFunctional (
            new ExtractTensorSubDomainFunctional3D<T,nDim>, domain, field, extractedField );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > extractSubDomain(MultiTensorField3D<T,nDim>& field, Box3D domain)
{
    MultiTensorField3D<T,nDim>* extractedField = new MultiTensorField3D<T,nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::auto_ptr<MultiTensorField3D<T,nDim> >(extractedField);
}


/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& component, Box3D domain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional3D<T,nDim>(iComponent), domain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, Box3D domain, int iComponent)
{
    MultiScalarField3D<T>* component = new MultiScalarField3D<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return std::auto_ptr<MultiScalarField3D<T> >(component);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > extractComponent(MultiTensorField3D<T,nDim>& tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& norm, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNormFunctional3D<T,nDim>, domain, norm, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField, Box3D domain)
{
    MultiScalarField3D<T>* norm = new MultiScalarField3D<T>(tensorField, domain);
    computeNorm(tensorField, *norm, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(norm);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNorm(MultiTensorField3D<T,nDim>& tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional3D<T,nDim>, domain, normSqr, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField, Box3D domain)
{
    MultiScalarField3D<T>* normSqr = new MultiScalarField3D<T>(tensorField, domain);
    computeNormSqr(tensorField, *normSqr, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(normSqr);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField3D<T> > computeNormSqr(MultiTensorField3D<T,nDim>& tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& norm, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional3D<T>, domain, norm, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    MultiScalarField3D<T>* norm = new MultiScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(norm);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNorm(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& normSqr, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional3D<T>, domain, normSqr, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    MultiScalarField3D<T>* normSqr = new MultiScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(normSqr);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorNormSqr(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, MultiScalarField3D<T>& trace, Box3D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional3D<T>, domain, trace, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField, Box3D domain)
{
    MultiScalarField3D<T>* trace = new MultiScalarField3D<T>(tensorField, domain);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return std::auto_ptr<MultiScalarField3D<T> >(trace);
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeSymmetricTensorTrace(MultiTensorField3D<T,6>& tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional3D<T,3>, domain, velocity, vorticity, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    MultiTensorField3D<T,3>* vorticity = new MultiTensorField3D<T,3>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return std::auto_ptr<MultiTensorField3D<T,3> >(vorticity);
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(MultiTensorField3D<T,3>& velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}


/* *************** Vorticity, without boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,3>& vorticity, Box3D domain)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional3D<T,3>, domain, velocity, vorticity );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    MultiTensorField3D<T,3>* vorticity = new MultiTensorField3D<T,3>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return std::auto_ptr<MultiTensorField3D<T,3> >(vorticity);
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,3> > computeBulkVorticity(MultiTensorField3D<T,3>& velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}



/* *************** Strain rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional3D<T,3>, domain, velocity, S, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    MultiTensorField3D<T,6>* S = new MultiTensorField3D<T,6>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return std::auto_ptr<MultiTensorField3D<T,6> >(S);
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeStrainRate(MultiTensorField3D<T,3>& velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}


/* *************** Str. rate, without boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, MultiTensorField3D<T,6>& S, Box3D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional3D<T,6>, domain, velocity, S );
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity, Box3D domain)
{
    MultiTensorField3D<T,6>* S = new MultiTensorField3D<T,6>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return std::auto_ptr<MultiTensorField3D<T,6> >(S);
}

template<typename T>
std::auto_ptr<MultiTensorField3D<T,6> > computeBulkStrainRate(MultiTensorField3D<T,3>& velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}


/* *************** MultiTensorField - MultiTensorField operations *************** */

template<typename T, int nDim>
void add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    MultiTensorField3D<T,nDim>* result = new MultiTensorField3D<T,nDim>(A, domain);
    add(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField3D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > add(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    MultiTensorField3D<T,nDim>* result = new MultiTensorField3D<T,nDim>(A, domain);
    subtract(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField3D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > subtract(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    MultiTensorField3D<T,nDim>* result = new MultiTensorField3D<T,nDim>(A, domain);
    multiply(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField3D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > multiply(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return multiply(A, B, A.getBoundingBox());
}

template<typename T, int nDim>
void divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, MultiTensorField3D<T,nDim>& result, Box3D domain)
{
    std::vector<MultiTensorField3D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional3D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain)
{
    MultiTensorField3D<T,nDim>* result = new MultiTensorField3D<T,nDim>(A, domain);
    divide(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField3D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField3D<T,nDim> > divide(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B)
{
    return divide(A, B, A.getBoundingBox());
}


/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void addInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void subtractInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B, Box3D domain) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void divideInPlace(MultiTensorField3D<T,nDim>& A, MultiTensorField3D<T,nDim>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}

}  // namespace plb

#endif  // MULTI_DATA_ANALYSIS_3D_HH
