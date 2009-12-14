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
#include "multiBlock/multiDataAnalysis2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataCouplingWrapper2D.h"
#include "core/dataAnalysisFunctionals2D.h"

#ifndef MULTI_DATA_ANALYSIS_2D_HH
#define MULTI_DATA_ANALYSIS_2D_HH

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Extract Sub-Lattice ******************************* */

template<typename T, template<typename U> class Descriptor>
void extractSubDomain( MultiBlockLattice2D<T,Descriptor>& lattice,
                       MultiBlockLattice2D<T,Descriptor>& extractedLattice,
                       Box2D domain)
{
    applyProcessingFunctional (
            new ExtractLatticeSubDomainFunctional2D<T,Descriptor>, domain, lattice, extractedLattice );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiBlockLattice2D<T,Descriptor> > extractSubDomain(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiBlockLattice2D<T,Descriptor>* extractedLattice = new MultiBlockLattice2D<T,Descriptor>(lattice, domain);
    extractSubDomain(lattice, *extractedLattice, domain);
    return std::auto_ptr<MultiBlockLattice2D<T,Descriptor> >(extractedLattice);
}


/* *************** Density ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& density, Box2D domain)
{
    applyProcessingFunctional (
            new BoxDensityFunctional2D<T,Descriptor>, domain, lattice, density );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiScalarField2D<T>* density = new MultiScalarField2D<T>(lattice, domain);
    computeDensity(lattice, *density, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(density);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeDensity(MultiBlockLattice2D<T,Descriptor>& lattice) {
    return computeDensity(lattice, lattice.getBoundingBox());
}


/* *************** RhoBar ******************************************* */

template<typename T, template<typename U> class Descriptor>
void computeRhoBar(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& rhoBar, Box2D domain)
{
    applyProcessingFunctional (
            new BoxRhoBarFunctional2D<T,Descriptor>, domain, lattice, rhoBar );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeRhoBar(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiScalarField2D<T>* rhoBar = new MultiScalarField2D<T>(lattice, domain);
    computeRhoBar(lattice, *rhoBar, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeRhoBar(MultiBlockLattice2D<T,Descriptor>& lattice) {
    return computeRhoBar(lattice, lattice.getBoundingBox());
}


/* *************** Kinetic Energy ************************************ */

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& energy, Box2D domain)
{
    applyProcessingFunctional (
            new BoxKineticEnergyFunctional2D<T,Descriptor>, domain, lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiScalarField2D<T>* energy = new MultiScalarField2D<T>(lattice, domain);
    computeKineticEnergy(lattice, *energy, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(energy);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeKineticEnergy(MultiBlockLattice2D<T,Descriptor>& lattice) {
    return computeKineticEnergy(lattice, lattice.getBoundingBox());
}


/* *************** Velocity Norm ************************************* */

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& velocityNorm, Box2D domain)
{
    applyProcessingFunctional (
            new BoxVelocityNormFunctional2D<T,Descriptor>, domain, lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiScalarField2D<T>* velocityNorm = new MultiScalarField2D<T>(lattice, domain);
    computeVelocityNorm(lattice, *velocityNorm, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(velocityNorm);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm(MultiBlockLattice2D<T,Descriptor>& lattice) {
    return computeVelocityNorm(lattice, lattice.getBoundingBox());
}


/* *************** Velocity Component ******************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& velocityComponent,
                              Box2D domain, plint iComponent)
{
    applyProcessingFunctional (
            new BoxVelocityComponentFunctional2D<T,Descriptor>(iComponent), domain, lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice,
                                                               Box2D domain, plint iComponent)
{
    MultiScalarField2D<T>* velocityComponent = new MultiScalarField2D<T>(lattice, domain);
    computeVelocityComponent(lattice, *velocityComponent, domain, iComponent);
    return std::auto_ptr<MultiScalarField2D<T> >(velocityComponent);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(MultiBlockLattice2D<T,Descriptor>& lattice, plint iComponent)
{
    return computeVelocityComponent(lattice, lattice.getBoundingBox(), iComponent);
}


/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiTensorField2D<T,Descriptor<T>::d>& velocity, Box2D domain)
{
    applyProcessingFunctional (
            new BoxVelocityFunctional2D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,Descriptor<T>::d> > computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiTensorField2D<T,Descriptor<T>::d>* velocity
        = new MultiTensorField2D<T,Descriptor<T>::d>(lattice, domain);
    computeVelocity(lattice, *velocity, domain);
    return std::auto_ptr<MultiTensorField2D<T,Descriptor<T>::d> >(velocity);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,Descriptor<T>::d> >
    computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice)
{
    return computeVelocity(lattice, lattice.getBoundingBox());
}


/* *************** Deviatoric Stress ********************************* */

template<typename T, template<typename U> class Descriptor>
void computeDeviatoricStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq,
                              Box2D domain )
{
    applyProcessingFunctional (
            new BoxDeviatoricStressFunctional2D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>* PiNeq
        = new MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeDeviatoricStress(lattice, *PiNeq, domain);
    return std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >(PiNeq);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeDeviatoricStress(MultiBlockLattice2D<T,Descriptor>& lattice)
{
    return computeDeviatoricStress(lattice, lattice.getBoundingBox());
}


/* *************** Strain Rate from Stress *************************** */

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                                  MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>& S,
                                  Box2D domain )
{
    applyProcessingFunctional (
            new BoxStrainRateFromStressFunctional2D<T,Descriptor>, domain, lattice, S );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain)
{
    MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>* S
        = new MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n>(lattice, domain);
    computeStrainRateFromStress(lattice, *S, domain);
    return std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >(S);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> >
    computeStrainRateFromStress(MultiBlockLattice2D<T,Descriptor>& lattice)
{
    return computeStrainRateFromStress(lattice, lattice.getBoundingBox());
}


/* *************** Population **************************************** */

template<typename T, template<typename U> class Descriptor>
void computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& population,
                       Box2D domain, plint iPop)
{
    applyProcessingFunctional (
            new BoxPopulationFunctional2D<T,Descriptor>(iPop), domain, lattice, population );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice,
                                                        Box2D domain, plint iPop)
{
    MultiScalarField2D<T>* population = new MultiScalarField2D<T>(lattice, domain);
    computePopulation(lattice, *population, domain, iPop);
    return std::auto_ptr<MultiScalarField2D<T> >(population);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField2D<T> > computePopulation(MultiBlockLattice2D<T,Descriptor>& lattice, plint iPop)
{
    return computePopulation(lattice, lattice.getBoundingBox(), iPop);
}


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-ScalarField *************************** */

template<typename T>
void extractSubDomain( MultiScalarField2D<T>& field,
                       MultiScalarField2D<T>& extractedField,
                       Box2D domain)
{
    applyProcessingFunctional (
            new ExtractScalarSubDomainFunctional2D<T>, domain, field, extractedField );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > extractSubDomain(MultiScalarField2D<T>& field, Box2D domain)
{
    MultiScalarField2D<T>* extractedField = new MultiScalarField2D<T>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(extractedField);
}

/* *************** MultiScalarField - Scalar operations *************** */

template<typename T>
void add(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& field, T scalar, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    add(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& field, T scalar)
{
    return add(field, scalar, field.getBoundingBox());
}


template<typename T>
void add(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_plus_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T>& field, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    add(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(T scalar, MultiScalarField2D<T>& field)
{
    return add(scalar, field, field.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_minus_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& field, T scalar, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    subtract(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& field, T scalar)
{
    return subtract(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtract(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new Alpha_minus_A_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T>& field, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    subtract(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(T scalar, MultiScalarField2D<T>& field)
{
    return subtract(scalar, field, field.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& field, T scalar, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    multiply(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& field, T scalar)
{
    return multiply(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiply(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_times_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T>& field, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    multiply(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(T scalar, MultiScalarField2D<T>& field)
{
    return multiply(scalar, field, field.getBoundingBox());
}


template<typename T>
void divide(MultiScalarField2D<T>& field, T scalar, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new A_dividedBy_alpha_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& field, T scalar, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    divide(field, scalar, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& field, T scalar)
{
    return divide(field, scalar, field.getBoundingBox());
}


template<typename T>
void divide(T scalar, MultiScalarField2D<T>& field, MultiScalarField2D<T>& result, Box2D domain)
{
    applyProcessingFunctional (
            new Alpha_dividedBy_A_functional2D<T>(scalar), domain, field, result );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T>& field, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(field, domain);
    divide(scalar, field, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(T scalar, MultiScalarField2D<T>& field)
{
    return divide(scalar, field, field.getBoundingBox());
}


/* *************** MultiScalarField - Scalar inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain) {
    applyProcessingFunctional (
            new A_plus_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template<typename T>
void addInPlace(MultiScalarField2D<T>& field, T scalar) {
    addInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain) {
    applyProcessingFunctional (
            new A_minus_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template<typename T>
void subtractInPlace(MultiScalarField2D<T>& field, T scalar) {
    subtractInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain) {
    applyProcessingFunctional (
            new A_times_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& field, T scalar) {
    multiplyInPlace(field, scalar, field.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField2D<T>& field, T scalar, Box2D domain) {
    applyProcessingFunctional (
            new A_dividedBy_alpha_inplace_functional2D<T>(scalar), domain, field);
}

template<typename T>
void divideInPlace(MultiScalarField2D<T>& field, T scalar) {
    divideInPlace(field, scalar, field.getBoundingBox());
}


/* *************** MultiScalarField - MultiScalarField operations *************** */

template<typename T>
void add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain)
{
    std::vector<MultiScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_plus_B_functional2D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(A, domain);
    add(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > add(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T>
void subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain)
{
    std::vector<MultiScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_minus_B_functional2D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(A, domain);
    subtract(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > subtract(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T>
void multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain)
{
    std::vector<MultiScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_times_B_functional2D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(A, domain);
    multiply(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > multiply(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B)
{
    return multiply(A, B, A.getBoundingBox());
}

template<typename T>
void divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, MultiScalarField2D<T>& result, Box2D domain)
{
    std::vector<MultiScalarField2D<T>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new A_dividedBy_B_functional2D<T>, domain, fields );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain)
{
    MultiScalarField2D<T>* result = new MultiScalarField2D<T>(A, domain);
    divide(A, B, *result, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(result);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > divide(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B)
{
    return divide(A, B, A.getBoundingBox());
}


/* *************** ScalarField - ScalarField inplace operations *************** */

template<typename T>
void addInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_plus_B_inplace_functional2D<T>, domain, A, B );
}

template<typename T>
void addInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void subtractInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_minus_B_inplace_functional2D<T>, domain, A, B );
}

template<typename T>
void subtractInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_times_B_inplace_functional2D<T>, domain, A, B );
}

template<typename T>
void multiplyInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T>
void divideInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B, Box2D domain) {
    applyProcessingFunctional (
            new A_dividedBy_B_inplace_functional2D<T>, domain, A, B );
}

template<typename T>
void divideInPlace(MultiScalarField2D<T>& A, MultiScalarField2D<T>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiTensorField2D<T,nDim>& field,
                       MultiTensorField2D<T,nDim>& extractedField,
                       Box2D domain)
{
    applyProcessingFunctional (
            new ExtractTensorSubDomainFunctional2D<T,nDim>, domain, field, extractedField );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > extractSubDomain(MultiTensorField2D<T,nDim>& field, Box2D domain)
{
    MultiTensorField2D<T,nDim>* extractedField = new MultiTensorField2D<T,nDim>(field, domain);
    extractSubDomain(field, *extractedField, domain);
    return std::auto_ptr<MultiTensorField2D<T,nDim> >(extractedField);
}


/* *************** Component (scalar-field) out of a tensor-field ****** */

template<typename T, int nDim>
void extractComponent(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& component, Box2D domain, int iComponent)
{
    applyProcessingFunctional (
            new ExtractTensorComponentFunctional2D<T,nDim>(iComponent), domain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > extractComponent(MultiTensorField2D<T,nDim>& tensorField, Box2D domain, int iComponent)
{
    MultiScalarField2D<T>* component = new MultiScalarField2D<T>(tensorField, domain);
    extractComponent(tensorField, *component, domain, iComponent);
    return std::auto_ptr<MultiScalarField2D<T> >(component);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > extractComponent(MultiTensorField2D<T,nDim>& tensorField, int iComponent)
{
    return extractComponent(tensorField, tensorField.getBoundingBox(), iComponent);
}


/* *************** Vector-norm of each cell in the field *************** */

template<typename T, int nDim>
void computeNorm(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& component, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeNormFunctional2D<T,nDim>, domain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T,nDim>& tensorField, Box2D domain)
{
    MultiScalarField2D<T>* component = new MultiScalarField2D<T>(tensorField, domain);
    computeNorm(tensorField, *component, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(component);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNorm(MultiTensorField2D<T,nDim>& tensorField)
{
    return computeNorm(tensorField, tensorField.getBoundingBox());
}


/* *************** Squared vector-norm of each cell in the field ******** */

template<typename T, int nDim>
void computeNormSqr(MultiTensorField2D<T,nDim>& tensorField, MultiScalarField2D<T>& component, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeNormSqrFunctional2D<T,nDim>, domain, component, tensorField );
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T,nDim>& tensorField, Box2D domain)
{
    MultiScalarField2D<T>* component = new MultiScalarField2D<T>(tensorField, domain);
    computeNormSqr(tensorField, *component, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(component);
}

template<typename T, int nDim>
std::auto_ptr<MultiScalarField2D<T> > computeNormSqr(MultiTensorField2D<T,nDim>& tensorField)
{
    return computeNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Tensor-norm of each symmetric tensor of a field ***** */

template<typename T>
void computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& norm, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormFunctional2D<T>, domain, norm, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField, Box2D domain)
{
    MultiScalarField2D<T>* norm = new MultiScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorNorm(tensorField, *norm, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(norm);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNorm(MultiTensorField2D<T,3>& tensorField)
{
    return computeSymmetricTensorNorm(tensorField, tensorField.getBoundingBox());
}


/* *************** Squared Tensor-norm of each symmetric tensor of a field*/

template<typename T>
void computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& normSqr, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorNormSqrFunctional2D<T>, domain, normSqr, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField, Box2D domain)
{
    MultiScalarField2D<T>* normSqr = new MultiScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorNormSqr(tensorField, *normSqr, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(normSqr);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorNormSqr(MultiTensorField2D<T,3>& tensorField)
{
    return computeSymmetricTensorNormSqr(tensorField, tensorField.getBoundingBox());
}


/* *************** Trace of each symmetric tensor of a field ************ */

template<typename T>
void computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField, MultiScalarField2D<T>& trace, Box2D domain)
{
    applyProcessingFunctional (
            new ComputeSymmetricTensorTraceFunctional2D<T>, domain, trace, tensorField );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField, Box2D domain)
{
    MultiScalarField2D<T>* trace = new MultiScalarField2D<T>(tensorField, domain);
    computeSymmetricTensorTrace(tensorField, *trace, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(trace);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeSymmetricTensorTrace(MultiTensorField2D<T,3>& tensorField)
{
    return computeSymmetricTensorTrace(tensorField, tensorField.getBoundingBox());
}


/* *************** Vorticity from Velocity field *********************** */

template<typename T>
void computeVorticity(MultiTensorField2D<T,2>& velocity, MultiScalarField2D<T>& vorticity, Box2D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxVorticityFunctional2D<T,2>, domain, vorticity, velocity, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T,2>& velocity, Box2D domain)
{
    MultiScalarField2D<T>* vorticity = new MultiScalarField2D<T>(velocity, domain);
    computeVorticity(velocity, *vorticity, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(vorticity);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeVorticity(MultiTensorField2D<T,2>& velocity)
{
    return computeVorticity(velocity, velocity.getBoundingBox());
}


/* *************** Vorticity, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkVorticity(MultiTensorField2D<T,2>& velocity, MultiScalarField2D<T>& vorticity, Box2D domain)
{
    applyProcessingFunctional (
            new BoxBulkVorticityFunctional2D<T,2>, domain, vorticity, velocity );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T,2>& velocity, Box2D domain)
{
    MultiScalarField2D<T>* vorticity = new MultiScalarField2D<T>(velocity, domain);
    computeBulkVorticity(velocity, *vorticity, domain);
    return std::auto_ptr<MultiScalarField2D<T> >(vorticity);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > computeBulkVorticity(MultiTensorField2D<T,2>& velocity)
{
    return computeBulkVorticity(velocity, velocity.getBoundingBox());
}



/* *************** Strain Rate from Velocity field ********************* */

template<typename T>
void computeStrainRate(MultiTensorField2D<T,2>& velocity, MultiTensorField2D<T,3>& S, Box2D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxStrainRateFunctional2D<T,2>, domain, S, velocity, envelopeWidth );
}

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeStrainRate(MultiTensorField2D<T,2>& velocity, Box2D domain)
{
    MultiTensorField2D<T,3>* S = new MultiTensorField2D<T,3>(velocity, domain);
    computeStrainRate(velocity, *S, domain);
    return std::auto_ptr<MultiTensorField2D<T,3> >(S);
}

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeStrainRate(MultiTensorField2D<T,2>& velocity)
{
    return computeStrainRate(velocity, velocity.getBoundingBox());
}


/* *************** Str. rate, witout boundary treatment, from Velocity field  */

template<typename T>
void computeBulkStrainRate(MultiTensorField2D<T,2>& velocity, MultiTensorField2D<T,3>& S, Box2D domain)
{
    plint envelopeWidth=1;
    applyProcessingFunctional (
            new BoxBulkStrainRateFunctional2D<T,2>, domain, velocity, S );
}

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeBulkStrainRate(MultiTensorField2D<T,2>& velocity, Box2D domain)
{
    MultiTensorField2D<T,3>* S = new MultiTensorField2D<T,3>(velocity, domain);
    computeBulkStrainRate(velocity, *S, domain);
    return std::auto_ptr<MultiTensorField2D<T,3> >(S);
}

template<typename T>
std::auto_ptr<MultiTensorField2D<T,3> > computeBulkStrainRate(MultiTensorField2D<T,2>& velocity)
{
    return computeBulkStrainRate(velocity, velocity.getBoundingBox());
}


/* *************** MultiTensorField - MultiTensorField operations *************** */

template<typename T, int nDim>
void add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain)
{
    std::vector<MultiTensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_plus_B_functional2D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain)
{
    MultiTensorField2D<T,nDim>* result = new MultiTensorField2D<T,nDim>(A, domain);
    add(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField2D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > add(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B)
{
    return add(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain)
{
    std::vector<MultiTensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_minus_B_functional2D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain)
{
    MultiTensorField2D<T,nDim>* result = new MultiTensorField2D<T,nDim>(A, domain);
    subtract(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField2D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > subtract(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B)
{
    return subtract(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain)
{
    std::vector<MultiTensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_times_B_functional2D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain)
{
    MultiTensorField2D<T,nDim>* result = new MultiTensorField2D<T,nDim>(A, domain);
    multiply(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField2D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > multiply(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B)
{
    return multiply(A, B, A.getBoundingBox());
}

template<typename T, int nDim>
void divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, MultiTensorField2D<T,nDim>& result, Box2D domain)
{
    std::vector<MultiTensorField2D<T,nDim>* > fields;
    fields.push_back(&A);
    fields.push_back(&B);
    fields.push_back(&result);
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_functional2D<T,nDim>, domain, fields );
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain)
{
    MultiTensorField2D<T,nDim>* result = new MultiTensorField2D<T,nDim>(A, domain);
    divide(A, B, *result, domain);
    return std::auto_ptr<MultiTensorField2D<T,nDim> >(result);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > divide(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B)
{
    return divide(A, B, A.getBoundingBox());
}


/* *************** TensorField - TensorField inplace operations *************** */

template<typename T, int nDim>
void addInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain) {
    applyProcessingFunctional (
            new Tensor_A_plus_B_inplace_functional2D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void addInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B) {
    addInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain) {
    applyProcessingFunctional (
            new Tensor_A_minus_B_inplace_functional2D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void subtractInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B) {
    subtractInPlace(A, B, A.getBoundingBox());
}


template<typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain) {
    applyProcessingFunctional (
            new Tensor_A_times_B_inplace_functional2D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void multiplyInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B) {
    multiplyInPlace(A,B, A.getBoundingBox());
}


template<typename T, int nDim>
void divideInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B, Box2D domain) {
    applyProcessingFunctional (
            new Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>, domain, A, B );
}

template<typename T, int nDim>
void divideInPlace(MultiTensorField2D<T,nDim>& A, MultiTensorField2D<T,nDim>& B) {
    divideInPlace(A, B, A.getBoundingBox());
}

}  // namespace plb

#endif  // MULTI_DATA_ANALYSIS_2D_HH
