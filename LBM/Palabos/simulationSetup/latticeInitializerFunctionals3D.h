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
 * Functionals for domain initialization -- header file.
 */
#ifndef LATTICE_INITIALIZER_FUNCTIONALS_3D_H
#define LATTICE_INITIALIZER_FUNCTIONALS_3D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
struct OneCellFunctional3D {
    virtual ~OneCellFunctional3D();
    virtual OneCellFunctional3D<T,Descriptor>* clone() const =0;
    virtual void execute(Cell<T,Descriptor>& cell) const =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
};

template<typename T, template<class U> class Descriptor>
struct OneCellIndexedFunctional3D {
    virtual ~OneCellIndexedFunctional3D();
    virtual OneCellIndexedFunctional3D<T,Descriptor>* clone() const =0;
    virtual void execute(plint iX, plint iY, plint iZ, Cell<T,Descriptor>& cell) const =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
};

struct DomainFunctional3D {
    virtual ~DomainFunctional3D() { }
    virtual bool operator() (plint iX, plint iY, plint iZ) const =0;
    virtual DomainFunctional3D* clone() const =0;
};

template<typename T, template<class U> class Descriptor>
class GenericLatticeFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    GenericLatticeFunctional3D(OneCellFunctional3D<T,Descriptor>* f_);
    GenericLatticeFunctional3D(GenericLatticeFunctional3D<T,Descriptor> const& rhs);
    ~GenericLatticeFunctional3D();
    GenericLatticeFunctional3D<T,Descriptor>& operator=(GenericLatticeFunctional3D<T,Descriptor> const& rhs);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual GenericLatticeFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
private:
    OneCellFunctional3D<T,Descriptor>* f;
};

template<typename T, template<class U> class Descriptor>
class GenericIndexedLatticeFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    GenericIndexedLatticeFunctional3D(OneCellIndexedFunctional3D<T,Descriptor>* f_);
    GenericIndexedLatticeFunctional3D(GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs);
    ~GenericIndexedLatticeFunctional3D();
    GenericIndexedLatticeFunctional3D<T,Descriptor>& operator=(GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual GenericIndexedLatticeFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
private:
    OneCellIndexedFunctional3D<T,Descriptor>* f;
};


/* *************** Class InstantiateDynamicsFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateDynamicsFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    InstantiateDynamicsFunctional3D(Dynamics<T,Descriptor>* dynamics_);
    InstantiateDynamicsFunctional3D(InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs);
    InstantiateDynamicsFunctional3D<T,Descriptor>& operator= (
            InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs );
    ~InstantiateDynamicsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateDynamicsFunctional3D<T,Descriptor>* clone() const ;
private:
    Dynamics<T,Descriptor>* dynamics;
};


/* ************* Class InstantiateComplexDomainDynamicsFunctional3D ** */

template<typename T, template<typename U> class Descriptor>
class InstantiateComplexDomainDynamicsFunctional3D
    : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    InstantiateComplexDomainDynamicsFunctional3D( Dynamics<T,Descriptor>* dynamics_,
                                                  DomainFunctional3D* domain_ );
    InstantiateComplexDomainDynamicsFunctional3D (
            InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs );
    InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>& operator= (
            InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs );
    ~InstantiateComplexDomainDynamicsFunctional3D();
    virtual void process(Box3D boundingBox, BlockLattice3D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    DomainFunctional3D* domain;
};


/* ************* Class InstantiateDotDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateDotDynamicsFunctional3D : public DotProcessingFunctional3D_L<T,Descriptor> {
public:
    InstantiateDotDynamicsFunctional3D(Dynamics<T,Descriptor>* dynamics_);
    InstantiateDotDynamicsFunctional3D(InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs);
    InstantiateDotDynamicsFunctional3D<T,Descriptor>& operator= (
            InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs );
    ~InstantiateDotDynamicsFunctional3D();
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateDotDynamicsFunctional3D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
};


/* ************* Class DynamicsFromMaskFunctional3D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
/** Note that the boolean mask is of type MultiScalarField3D<T> instead of MultiScalarField3D<bool>,
 *  because coupling lattices of different data types is not possible.
 */
template<typename T, template<typename U> class Descriptor>
class DynamicsFromMaskFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor> {
public:
    DynamicsFromMaskFunctional3D(Dynamics<T,Descriptor>* dynamics_, bool whichFlag_);
    DynamicsFromMaskFunctional3D(DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs);
    DynamicsFromMaskFunctional3D<T,Descriptor>& operator= (
            DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs );
    ~DynamicsFromMaskFunctional3D();
    virtual void process (
            Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                          ScalarField3D<T>& mask );
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DynamicsFromMaskFunctional3D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    bool whichFlag;
};


/* ************* Class DynamicsFromIntMaskFunctional3D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
/** Note that the boolean mask is of type MultiScalarField3D<T> instead of MultiScalarField3D<bool>,
 *  because coupling lattices of different data types is not possible.
 */
template<typename T, template<typename U> class Descriptor>
class DynamicsFromIntMaskFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor> {
public:
    DynamicsFromIntMaskFunctional3D(Dynamics<T,Descriptor>* dynamics_, int whichFlag_);
    DynamicsFromIntMaskFunctional3D(DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs);
    DynamicsFromIntMaskFunctional3D<T,Descriptor>& operator= (
            DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs );
    ~DynamicsFromIntMaskFunctional3D();
    virtual void process (
            Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                          ScalarField3D<T>& mask );
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DynamicsFromIntMaskFunctional3D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    int whichFlag;
};


/* ************* Class SetConstBoundaryVelocityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class SetConstBoundaryVelocityFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    SetConstBoundaryVelocityFunctional3D(Array<T,Descriptor<T>::d> velocity);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual SetConstBoundaryVelocityFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    Array<T,Descriptor<T>::d> u;
};


/* ************* Class SetConstBoundaryDensityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class SetConstBoundaryDensityFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    SetConstBoundaryDensityFunctional3D(T rho_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual SetConstBoundaryDensityFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T rho;
};


/* ************* Class IniConstEquilibriumFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class IniConstEquilibriumFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    IniConstEquilibriumFunctional3D(T density_, Array<T,Descriptor<T>::d> velocity);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual IniConstEquilibriumFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    T jSqr;
};

/* ************* Class StripeOffDensityOffsetFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class StripeOffDensityOffsetFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    StripeOffDensityOffsetFunctional3D(T deltaRho_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual StripeOffDensityOffsetFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T deltaRho;
};


/* ************* Class InstantiateCompositeDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateCompositeDynamicsFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    InstantiateCompositeDynamicsFunctional3D(CompositeDynamics<T,Descriptor>* compositeDynamics_);
    InstantiateCompositeDynamicsFunctional3D(InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs);
    InstantiateCompositeDynamicsFunctional3D<T,Descriptor>& operator= (
            InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs );
    ~InstantiateCompositeDynamicsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateCompositeDynamicsFunctional3D<T,Descriptor>* clone() const;
private:
    CompositeDynamics<T,Descriptor>* compositeDynamics;
};

}  // namespace plb

#endif  // LATTICE_INITIALIZER_FUNCTIONALS_3D_H
