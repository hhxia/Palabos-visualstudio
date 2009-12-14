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
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef LATTICE_INITIALIZER_FUNCTIONALS_3D_HH
#define LATTICE_INITIALIZER_FUNCTIONALS_3D_HH

#include "simulationSetup/latticeInitializer3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
OneCellFunctional3D<T,Descriptor>::~OneCellFunctional3D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}


template<typename T, template<class U> class Descriptor>
void OneCellFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale)
{ }


template<typename T, template<class U> class Descriptor>
OneCellIndexedFunctional3D<T,Descriptor>::~OneCellIndexedFunctional3D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale)
{ }


template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::GenericLatticeFunctional3D (
        OneCellFunctional3D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::GenericLatticeFunctional3D (
        GenericLatticeFunctional3D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::~GenericLatticeFunctional3D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>& GenericLatticeFunctional3D<T,Descriptor>::operator= (
        GenericLatticeFunctional3D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                f->execute(lattice.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>*
    GenericLatticeFunctional3D<T,Descriptor>::clone() const
{
    return new GenericLatticeFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    f->getModificationPattern(isWritten);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericLatticeFunctional3D<T,Descriptor>::appliesTo() const {
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale) {
    f->rescale(dxScale, dtScale);
}



template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::GenericIndexedLatticeFunctional3D (
        OneCellIndexedFunctional3D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::GenericIndexedLatticeFunctional3D (
        GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::~GenericIndexedLatticeFunctional3D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>&
    GenericIndexedLatticeFunctional3D<T,Descriptor>::operator= (
        GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                f->execute ( iX+relativeOffset.x,
                             iY+relativeOffset.y,
                             iZ+relativeOffset.z,
                             lattice.get(iX,iY,iZ) );
            }
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>*
    GenericIndexedLatticeFunctional3D<T,Descriptor>::clone() const
{
    return new GenericIndexedLatticeFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    f->getModificationPattern(isWritten);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericIndexedLatticeFunctional3D<T,Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    f->rescale(dxScale, dtScale);
}


/* *************** Class InstantiateDynamicsFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::InstantiateDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::InstantiateDynamicsFunctional3D (
        InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>&
    InstantiateDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::~InstantiateDynamicsFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>*
    InstantiateDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class InstantiateComplexDomainDynamicsFunctional3D ** */

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, DomainFunctional3D* domain_ )
    : dynamics(dynamics_),
      domain(domain_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional3D (
        InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      domain(rhs.domain->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>&
    InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    delete domain; domain = rhs.domain->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::~InstantiateComplexDomainDynamicsFunctional3D()
{
    delete dynamics;
    delete domain;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::process (
        Box3D boundingBox, BlockLattice3D<T,Descriptor>& lattice )
{

    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
        for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
            for (plint iZ=boundingBox.z0; iZ<=boundingBox.z1; ++iZ) {
                if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y,iZ+relativeOffset.z)) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>*
    InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class InstantiateDotDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::InstantiateDotDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::InstantiateDotDynamicsFunctional3D (
        InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>&
    InstantiateDotDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::~InstantiateDotDynamicsFunctional3D()
{
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T,Descriptor>::process (
        DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
        Dot3D const& dot = dotList.getDot(iDot);
        lattice.attributeDynamics(dot.x,dot.y,dot.z, dynamics->clone());
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDotDynamicsFunctional3D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>*
    InstantiateDotDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateDotDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromMaskFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::DynamicsFromMaskFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, bool whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::DynamicsFromMaskFunctional3D (
        DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>&
    DynamicsFromMaskFunctional3D<T,Descriptor>::operator= (
        DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::~DynamicsFromMaskFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      ScalarField3D<T>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {
                bool flag = (bool) util::roundToInt (
                               mask.get(iX+offset.x, iY+offset.y, iZ+offset.z) );
                if ( util::boolIsEqual(flag, whichFlag) ) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromMaskFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>*
    DynamicsFromMaskFunctional3D<T,Descriptor>::clone() const 
{
    return new DynamicsFromMaskFunctional3D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromIntMaskFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::DynamicsFromIntMaskFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, int whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::DynamicsFromIntMaskFunctional3D (
        DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>&
    DynamicsFromIntMaskFunctional3D<T,Descriptor>::operator= (
        DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::~DynamicsFromIntMaskFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      ScalarField3D<T>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {
                bool flag = (bool) util::roundToInt (
                               mask.get(iX+offset.x, iY+offset.y, iZ+offset.z) );
                if ( flag == whichFlag ) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromIntMaskFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>*
    DynamicsFromIntMaskFunctional3D<T,Descriptor>::clone() const 
{
    return new DynamicsFromIntMaskFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetConstBoundaryVelocityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T,Descriptor>::SetConstBoundaryVelocityFunctional3D (
        Array<T,Descriptor<T>::d> velocity )
    : u(velocity)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).defineVelocity(u);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T,Descriptor>*
    SetConstBoundaryVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityFunctional3D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    T latticeVelocityScale = dtScale / dxScale;
    u[0] *= latticeVelocityScale;
    u[1] *= latticeVelocityScale;
    u[2] *= latticeVelocityScale;
}


/* ************* Class SetConstBoundaryDensityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T,Descriptor>::SetConstBoundaryDensityFunctional3D(T rho_)
    : rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).defineDensity(rho);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T,Descriptor>*
    SetConstBoundaryDensityFunctional3D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryDensityFunctional3D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale )
{
    // Nothing to do: density is scale invariant.
}


/* ************* Class IniConstEquilibriumFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T,Descriptor>::IniConstEquilibriumFunctional3D (
        T density_, Array<T,Descriptor<T>::d> velocity )
    : rhoBar(Descriptor<T>::rhoBar(density_)),
      j     (density_*velocity[0], density_*velocity[1], density_*velocity[2]),
      jSqr  (VectorTemplate<T,Descriptor>::normSqr(j))
{ }

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    lattice.get(iX,iY,iZ)[iPop] =
                        lattice.get(iX,iY,iZ).computeEquilibrium(iPop, rhoBar, j, jSqr);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T,Descriptor>*
    IniConstEquilibriumFunctional3D<T,Descriptor>::clone() const
{
    return new IniConstEquilibriumFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumFunctional3D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    T latticeVelocityScale = dtScale / dxScale;
    j[0] *= latticeVelocityScale;
    j[1] *= latticeVelocityScale;
    j[2] *= latticeVelocityScale;
    jSqr *= latticeVelocityScale*latticeVelocityScale;
}


/* ************* Class StripeOffDensityOffsetFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T,Descriptor>::StripeOffDensityOffsetFunctional3D(T deltaRho_)
    : deltaRho(deltaRho_)
{ }

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell         = lattice.get(iX,iY,iZ);
                Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
                plint orderOfDecomposition = 0;
                std::vector<T> rawData;
                dynamics.decompose(cell, rawData, orderOfDecomposition);
                T& rhoBar = rawData[0];
                rhoBar -= deltaRho;
                dynamics.recompose(cell, rawData, orderOfDecomposition);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T,Descriptor>*
    StripeOffDensityOffsetFunctional3D<T,Descriptor>::clone() const
{
    return new StripeOffDensityOffsetFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT StripeOffDensityOffsetFunctional3D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

// Density is scale invariant.
template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T,Descriptor>::rescale(T dxScale, T dtScale) { }


/* ************* Class InstantiateCompositeDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::InstantiateCompositeDynamicsFunctional3D (
        CompositeDynamics<T,Descriptor>* compositeDynamics_ )
    : compositeDynamics(compositeDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::InstantiateCompositeDynamicsFunctional3D (
        InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs )
    : compositeDynamics(rhs.compositeDynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>&
    InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete compositeDynamics; compositeDynamics = rhs.compositeDynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::~InstantiateCompositeDynamicsFunctional3D()
{
    delete compositeDynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Dynamics<T,Descriptor>* baseDynamics = lattice.get(iX,iY,iZ).getDynamics().cloneComposeable();
                Dynamics<T,Descriptor>* dynamics = compositeDynamics->cloneWithNewBase(baseDynamics);
                lattice.attributeDynamics(iX,iY,iZ, dynamics);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::appliesTo() const
{
    // Composite dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>*
    InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateCompositeDynamicsFunctional3D<T,Descriptor>(*this);
}

}  // namespace plb

#endif  // LATTICE_INITIALIZER_FUNCTIONALS_3D_HH
