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
#ifndef BOUNDARY_DYNAMICS_HH
#define BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/boundaryDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

/* *************** Class BoundaryCompositeDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
BoundaryCompositeDynamics<T,Descriptor>::BoundaryCompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : PreparePopulationsDynamics<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
BoundaryCompositeDynamics<T,Descriptor>* BoundaryCompositeDynamics<T,Descriptor>::clone() const
{
    return new BoundaryCompositeDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T BoundaryCompositeDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const {
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computeDensity(tmpCell);
}

template<typename T, template<typename U> class Descriptor>
T BoundaryCompositeDynamics<T,Descriptor>::computePressure(Cell<T,Descriptor> const& cell) const {
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computePressure(tmpCell);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeVelocity(Cell<T,Descriptor> const& cell,
                                                              Array<T,Descriptor<T>::d>& u ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    this->getBaseDynamics().computeVelocity(tmpCell, u);
}

template<typename T, template<typename U> class Descriptor>
T BoundaryCompositeDynamics<T,Descriptor>::computeTemperature(Cell<T,Descriptor> const& cell) const {
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computeTemperature(tmpCell);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computeDeviatoricStress(tmpCell, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeHeatFlux(Cell<T,Descriptor> const& cell,
                                                              Array<T,Descriptor<T>::d>& q ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    this->getBaseDynamics().computeHeatFlux(tmpCell, q);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeMoment (
        Cell<T,Descriptor> const& cell, plint momentId, T* moment ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    this->getBaseDynamics().computeMoment(tmpCell, momentId, moment);
}

template<typename T, template<typename U> class Descriptor>
T BoundaryCompositeDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computeRhoBar(tmpCell);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    this->getBaseDynamics().computeRhoBarJ(tmpCell, rhoBar, j);
}

template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::computeRhoBarJPiNeq (
        Cell<T,Descriptor> const& cell, T& rhoBar, Array<T,Descriptor<T>::d>& j,
        Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
{
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    this->getBaseDynamics().computeRhoBarJPiNeq(tmpCell, rhoBar, j, PiNeq);
}

template<typename T, template<typename U> class Descriptor>
T BoundaryCompositeDynamics<T,Descriptor>::computeEbar(Cell<T,Descriptor> const& cell) const {
    Cell<T,Descriptor> tmpCell(cell);
    this -> completePopulations(tmpCell);
    return this->getBaseDynamics().computeEbar(tmpCell);
}

/** Do nothing inside this functions. This defaults to a behavior where
 *  the dynamics of BoundaryCompositeDynamics is identical to baseDynamics.
 *  More interesting behavior is achieved in derived classes which overload
 *  method completePopulations().
 */
template<typename T, template<typename U> class Descriptor>
void BoundaryCompositeDynamics<T,Descriptor>::completePopulations(Cell<T,Descriptor>& cell) const
{ }


/* *************** Class StoreDensityDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
StoreDensityDynamics<T,Descriptor>::StoreDensityDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : BoundaryCompositeDynamics<T,Descriptor>(baseDynamics_)
{
    rhoBar = Descriptor<T>::rhoBar((T)1);
}

template<typename T, template<typename U> class Descriptor>
StoreDensityDynamics<T,Descriptor>* StoreDensityDynamics<T,Descriptor>::clone() const
{
    return new StoreDensityDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T StoreDensityDynamics<T,Descriptor>::computeDensity(Cell<T,Descriptor> const& cell) const
{
    return Descriptor<T>::fullRho(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityDynamics<T,Descriptor>::defineDensity (
        Cell <T,Descriptor>& cell, T rho_ )
{
    rhoBar = Descriptor<T>::rhoBar(rho_);
}

template<typename T, template<typename U> class Descriptor>
T StoreDensityDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const
{
    return rhoBar;
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell, T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar_ = rhoBar;
    T rho = Descriptor<T>::fullRho(rhoBar);
    this->computeVelocity(cell, j);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] *= rho;
    }
}


/* *************** Class StoreVelocityDynamics *********************** */

template<typename T, template<typename U> class Descriptor>
StoreVelocityDynamics<T,Descriptor>::StoreVelocityDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : BoundaryCompositeDynamics<T,Descriptor>(baseDynamics_)
{
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
StoreVelocityDynamics<T,Descriptor>* StoreVelocityDynamics<T,Descriptor>::clone() const
{
    return new StoreVelocityDynamics<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void StoreVelocityDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& u_ ) const
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u_[iD] = u[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
void StoreVelocityDynamics<T,Descriptor>::defineVelocity (
        Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u_ )
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u[iD] = u_[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
T StoreVelocityDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    T rho = this->computeDensity(cell);
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void StoreVelocityDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell,
        T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    T rho = this->computeDensity(cell);
    rhoBar_ = Descriptor<T>::rhoBar(rho);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }
}


/* *************** Class StoreDensityAndVelocityDynamics ************* */

template<typename T, template<typename U> class Descriptor>
StoreDensityAndVelocityDynamics<T,Descriptor>::StoreDensityAndVelocityDynamics (
        Dynamics<T,Descriptor>* baseDynamics_ )
    : BoundaryCompositeDynamics<T,Descriptor>(baseDynamics_)
{
    rhoBar = Descriptor<T>::rhoBar((T)1);
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
StoreDensityAndVelocityDynamics<T,Descriptor>* StoreDensityAndVelocityDynamics<T,Descriptor>::clone() const
{
    return new StoreDensityAndVelocityDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T StoreDensityAndVelocityDynamics<T,Descriptor>::computeDensity (
        Cell<T,Descriptor> const& cell) const 
{
    return Descriptor<T>::fullRho(rhoBar);
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T,Descriptor>::computeVelocity (
        Cell<T,Descriptor> const& cell,
        Array<T,Descriptor<T>::d>& u_ ) const
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u_[iD] = u[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T,Descriptor>::defineDensity (
        Cell<T,Descriptor>& cell, T rho_)
{
    rhoBar = Descriptor<T>::rhoBar(rho_);
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T,Descriptor>::defineVelocity (
        Cell<T,Descriptor>& cell,
        Array<T,Descriptor<T>::d> const& u_ ) 
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u[iD] = u_[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
T StoreDensityAndVelocityDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    return rhoBar;
}

template<typename T, template<typename U> class Descriptor>
void StoreDensityAndVelocityDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell,
        T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar_ = rhoBar;
    T rho = Descriptor<T>::fullRho(rhoBar);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }
}


/* *************** Class StoreTemperatureAndVelocityDynamics ************* */

template<typename T, template<typename U> class Descriptor>
StoreTemperatureAndVelocityDynamics<T,Descriptor>::StoreTemperatureAndVelocityDynamics (
        Dynamics<T,Descriptor>* baseDynamics_ )
    : BoundaryCompositeDynamics<T,Descriptor>(baseDynamics_)
{
    thetaBar = T();
    u.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
StoreTemperatureAndVelocityDynamics<T,Descriptor>* StoreTemperatureAndVelocityDynamics<T,Descriptor>::clone() const {
    return new StoreTemperatureAndVelocityDynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T,Descriptor>::computeVelocity(
        Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::d>& u_ ) const
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u_[iD] = u[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
T StoreTemperatureAndVelocityDynamics<T,Descriptor>::computeTemperature (
        Cell<T,Descriptor> const& cell) const
{
    return thetaBar + (T)1;
}

template<typename T, template<typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T,Descriptor>::defineVelocity (
        Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u_ )
{
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u[iD] = u_[iD];
    }
}

template<typename T, template<typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T,Descriptor>::defineTemperature (
        Cell<T,Descriptor>& cell, T theta_ )
{
    thetaBar = theta_ - (T)1;
}

template<typename T, template<typename U> class Descriptor>
T StoreTemperatureAndVelocityDynamics<T,Descriptor>::computeRhoBar(Cell<T,Descriptor> const& cell) const {
    T rho = this->computeDensity(cell);
    return Descriptor<T>::rhoBar(rho);
}

template<typename T, template<typename U> class Descriptor>
void StoreTemperatureAndVelocityDynamics<T,Descriptor>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell,
        T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    T rho = this->computeDensity(cell);
    rhoBar_ = Descriptor<T>::rhoBar(rho);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = rho * u[iD];
    }
}


/* *************** Class VelocityDirichletBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    VelocityDirichletBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : StoreVelocityDynamics<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>*
    VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
T VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    computeRhoBar(Cell<T,Descriptor> const& cell) const 
{
    std::vector<int> const& onWallIndices
        = indexTemplates::subIndex<Descriptor<T>, direction, 0>();

    std::vector<int> const& normalIndices
        = indexTemplates::subIndex<Descriptor<T>, direction, orientation>();

    T rhoOnWall = T();
    for (pluint fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
        rhoOnWall += cell[onWallIndices[fIndex]];
    }

    T rhoNormal = T();
    for (pluint fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
        rhoNormal += cell[normalIndices[fIndex]];
    }

    T uNormal = (T)orientation * this->u[direction];
    T rhoBar =((T)2*rhoNormal+rhoOnWall-Descriptor<T>::SkordosFactor()*uNormal)
                  / ((T)1+uNormal);

    return rhoBar;
}


template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
T VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    computeDensity(Cell<T,Descriptor> const& cell) const 
{
    return Descriptor<T>::fullRho(computeRhoBar(cell));
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell,
        T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar_ = computeRhoBar(cell);
    T rho = Descriptor<T>::fullRho(rhoBar_);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = this->u[iD] * rho;
    }
}


/* *************** Class DensityDirichletBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    DensityDirichletBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : StoreDensityDynamics<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>*
    DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    computeJ( Cell<T,Descriptor> const& cell,
              Array<T,Descriptor<T>::d>& j_ ) const
{
    // All velocity components parallel to the wall are zero by definition.
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j_[iD] = T();
    }
    T rhoBar = this -> rhoBar;

    std::vector<int> const& onWallIndices
        = indexTemplates::subIndex<Descriptor<T>, direction, 0>();

    std::vector<int> const& normalIndices
        = indexTemplates::subIndex<Descriptor<T>, direction, orientation>();

    T rhoOnWall = T();
    for (pluint fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
        rhoOnWall += cell[onWallIndices[fIndex]];
    }

    T rhoNormal = T();
    for (pluint fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
        rhoNormal += cell[normalIndices[fIndex]];
    }

    j_[direction] = (T)orientation*((T)2*rhoNormal+rhoOnWall-rhoBar);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::
    computeVelocity( Cell<T,Descriptor> const& cell,
                     Array<T,Descriptor<T>::d>& u_ ) const
{
    this->computeJ(cell, u_);
    T invRho = Descriptor<T>::invRho(this->rhoBar);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        u_[iD] *= invRho;
    }
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>::computeRhoBarJ (
        Cell<T,Descriptor> const& cell,
        T& rhoBar_, Array<T,Descriptor<T>::d>& j) const
{
    rhoBar_ = this->rhoBar;
    this->computeJ(cell, j);
}

}  // namespace plb

#endif
