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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef DYNAMICS_HH
#define DYNAMICS_HH

#include "core/dynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "core/latticeStatistics.h"
#include <algorithm>
#include <limits>

namespace plb {

/* *************** Class Dynamics ******************************************* */

/* By default, this method forwards the call to clone().
 */
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* Dynamics<T,Descriptor>::cloneComposeable() const {
    return clone();
}

template<typename T, template<typename U> class Descriptor>
T Dynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    if (whichParameter == dynamicParams::omega_shear) {
        return getOmega();
    }
    return 0.;
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    if (whichParameter == dynamicParams::omega_shear) {
        setOmega(value);
    }
}


template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const {
    f = cell.getRawPopulations();
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::getExternalField (
        Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext) const
{
    PLB_PRECONDITION(pos+size <= Descriptor<T>::ExternalField::numScalars);
    T const* externalData = cell.getExternal(pos);
    for (plint iExt=0; iExt<size; ++iExt) {
        ext[iExt] = externalData[iExt];
    }
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f)
{
    cell.getRawPopulations() = f;
}

template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::setExternalField (
        Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext)
{
    PLB_PRECONDITION(pos+size <= Descriptor<T>::ExternalField::numScalars);
    T* externalData = cell.getExternal(pos);
    for (plint iExt=0; iExt<size; ++iExt) {
        externalData[iExt] = ext[iExt];
    }
}

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineDensity(Cell<T,Descriptor>& cell, T density) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineTemperature(Cell<T,Descriptor>& cell, T temperature) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q) { }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineDeviatoricStress(Cell<T,Descriptor>& cell, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq)
{ }

/** This method does nothing by default, unless overloaded in
 *  a deriving class. This is OK, because it's the behavior
 *  users will expect: if they use "defineSomething" on a cell,
 *  they want it to define the boundary condition on boundary
 *  nodes, but nothing to happen in bulk nodes. In this way,
 *  they can be lazy and loop over the whole domain instead
 *  of tracking boundary nodes explicitly.
 **/
template<typename T, template<typename U> class Descriptor>
void Dynamics<T,Descriptor>::defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value)
{ }


/* *************** Class BasicBulkDynamics *************************** */

template<typename T, template<typename U> class Descriptor>
BasicBulkDynamics<T,Descriptor>::BasicBulkDynamics(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
    {
    return momentTemplates<T,Descriptor>::compute_rho(cell);
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const
{
    return Descriptor<T>::cs2 * computeDensity(cell) * computeTemperature(cell);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    momentTemplates<T,Descriptor>::compute_uLb(cell, u);
}

/** Defaults to 1.
 */
template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const
{
    return (T)1;
}

/** Defaults to zero.
 */
template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

/** Defaults to zero.
 */
template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeMoment(
        Cell<T,Descriptor> const& cell, plint momentId, T* moment ) const
{
    PLB_PRECONDITION( false );
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::getOmega() const {
    return omega;
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::setOmega(T omega_) {
    omega = omega_;
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::get_rhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j ) const
{
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
}

template<typename T, template<typename U> class Descriptor>
void BasicBulkDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
T BasicBulkDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return momentTemplates<T,Descriptor>::get_eBar(cell);
}


/* *************** Class CompositeDynamics *************************** */

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::CompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : baseDynamics(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::CompositeDynamics(CompositeDynamics<T,Descriptor> const& rhs)
    : baseDynamics(rhs.baseDynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>::~CompositeDynamics() {
    delete baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>& CompositeDynamics<T,Descriptor>::operator=(CompositeDynamics<T,Descriptor> const& rhs)
{
    delete baseDynamics;
    baseDynamics = rhs.baseDynamics->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
CompositeDynamics<T,Descriptor>* CompositeDynamics<T,Descriptor>::cloneWithNewBase (
        Dynamics<T,Descriptor>* baseDynamics_) const
{
    // First, create a clone, based on the dynamic type of CompositeDynamics
    CompositeDynamics<T,Descriptor>* newDynamics = clone();
    // Then, replace its original baseDynamics by the one we've received as parameter
    newDynamics->replaceBaseDynamics(baseDynamics_);
    return newDynamics;
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::replaceBaseDynamics(Dynamics<T,Descriptor>* newBaseDynamics) {
    delete baseDynamics;
    baseDynamics = newBaseDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& CompositeDynamics<T,Descriptor>::getBaseDynamics() {
    return *baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& CompositeDynamics<T,Descriptor>::getBaseDynamics() const {
    return *baseDynamics;
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell, BlockStatistics<T>& statistics )
{
    prepareCollision(cell);
    baseDynamics -> collide(cell, statistics);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeEquilibrium (
        plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr, T thetaBar) const
{
    return baseDynamics -> computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::regularize (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{
    baseDynamics -> regularize(cell, rhoBar, j, jSqr, PiNeq, thetaBar);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeDensity(cell);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computePressure(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeVelocity(Cell<T,Descriptor> const& cell,
                                                          Array<T,Descriptor<T>::d>& u ) const
{
    this->getBaseDynamics().computeVelocity(cell, u);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeTemperature(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    return this->getBaseDynamics().computeDeviatoricStress(cell, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeHeatFlux(Cell<T,Descriptor> const& cell,
                                                          Array<T,Descriptor<T>::d>& q ) const
{
    this->getBaseDynamics().computeHeatFlux(cell, q);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment ) const
{
    this->getBaseDynamics().computeMoment(cell, momentId, moment);
}

template<typename T, template<typename U> class Descriptor>
plint CompositeDynamics<T,Descriptor>::numDecomposedVariables(plint order) const {
    return baseDynamics->numDecomposedVariables(order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    baseDynamics->decompose(cell, rawData, order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    baseDynamics->recompose(cell, rawData, order);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{
    baseDynamics->rescale(rawData, xDxInv, xDt, order);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return this->getBaseDynamics().computeRhoBar(cell);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j ) const
{
    this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    this->getBaseDynamics().computeRhoBarJPiNeq(cell, rhoBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const {
    return this->getBaseDynamics().computeEbar(cell);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::getOmega() const {
    return baseDynamics -> getOmega();
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setOmega(T omega_) {
    baseDynamics -> setOmega(omega_);
}

template<typename T, template<typename U> class Descriptor>
T CompositeDynamics<T,Descriptor>::getParameter(plint whichParameter) const {
    return baseDynamics -> getParameter(whichParameter);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setParameter(plint whichParameter, T value) {
    baseDynamics -> setParameter(whichParameter, value);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const {
    baseDynamics -> getPopulations(cell, f);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::getExternalField(Cell<T,Descriptor> const& cell,
                                                    plint pos, plint size, T* ext) const
{
    baseDynamics -> getExternalField(cell, pos, size, ext);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f)
{
    baseDynamics -> setPopulations(cell, f);
}

template<typename T, template<typename U> class Descriptor>
void CompositeDynamics<T,Descriptor>::setExternalField(Cell<T,Descriptor>& cell, plint pos,
                                                    plint size, const T* ext)
{
    baseDynamics -> setExternalField(cell, pos, size, ext);
}

/* *************** Class PreparePopulationsDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
PreparePopulationsDynamics<T,Descriptor>::PreparePopulationsDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
void PreparePopulationsDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell) {
    completePopulations(cell);
}


/* *************** Class BulkCompositeDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
BulkCompositeDynamics<T,Descriptor>::BulkCompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : PreparePopulationsDynamics<T,Descriptor>(baseDynamics_)
{ }

/* *************** Class BounceBack ********************************** */

template<typename T, template<typename U> class Descriptor>
BounceBack<T,Descriptor>::BounceBack(T rho_)
    :rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
BounceBack<T,Descriptor>* BounceBack<T,Descriptor>::clone() const {
    return new BounceBack<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{
    for (plint iPop=1; iPop <= Descriptor<T>::q/2; ++iPop) {
        std::swap(cell[iPop], cell[iPop+Descriptor<T>::q/2]);
    }
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                               T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return rho;
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::getOmega() const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T BounceBack<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint BounceBack<T,Descriptor>::numDecomposedVariables(plint order) const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    cell.serialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    PLB_PRECONDITION( (plint)rawData.size() == numDecomposedVariables(order) );
    cell.unSerialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void BounceBack<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }

/* *************** Class NoDynamics ********************************** */

template<typename T, template<typename U> class Descriptor>
NoDynamics<T,Descriptor>* NoDynamics<T,Descriptor>::clone() const {
    return new NoDynamics<T,Descriptor>(*this);
}
 
template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics<T>& statistics )
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                            T jSqr, T thetaBar) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::regularize(
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
        T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar ) const
{ }


template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    return (T)1;
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u) const
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeHeatFlux (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& q) const
{
    q.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment) const
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::getOmega() const {
    return T();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::setOmega(T omega_)
{ }

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return (T)1 - Descriptor<T>::SkordosFactor();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar = (T)1 - Descriptor<T>::SkordosFactor();
    j.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    rhoBar = (T)1 - Descriptor<T>::SkordosFactor();
    j.resetToZero();
    PiNeq.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
T NoDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const
{
    return T();
}

template<typename T, template<typename U> class Descriptor>
plint NoDynamics<T,Descriptor>::numDecomposedVariables(plint order) const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::decompose (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order ) const
{
    rawData.resize(numDecomposedVariables(order));
    cell.serialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::recompose (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order ) const
{
    PLB_PRECONDITION( (plint)rawData.size() == numDecomposedVariables(order) );
    cell.unSerialize(&rawData[0]);
}

template<typename T, template<typename U> class Descriptor>
void NoDynamics<T,Descriptor>::rescale (
        std::vector<T>& rawData, T xDxInv, T xDt, plint order ) const
{ }

}  // namespace plb

#endif  // DYNAMICS_HH
