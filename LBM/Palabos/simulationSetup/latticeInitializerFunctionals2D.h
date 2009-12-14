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
#ifndef LATTICE_INITIALIZER_FUNCTIONALS_2D_H
#define LATTICE_INITIALIZER_FUNCTIONALS_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
struct OneCellFunctional2D {
    virtual ~OneCellFunctional2D();
    virtual OneCellFunctional2D<T,Descriptor>* clone() const =0;
    virtual void execute(Cell<T,Descriptor>& cell) const =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
};

template<typename T, template<class U> class Descriptor>
struct OneCellIndexedFunctional2D {
    virtual ~OneCellIndexedFunctional2D();
    virtual OneCellIndexedFunctional2D<T,Descriptor>* clone() const =0;
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
};

struct DomainFunctional2D {
    virtual ~DomainFunctional2D() { }
    virtual bool operator() (plint iX, plint iY) const =0;
    virtual DomainFunctional2D* clone() const =0;
};

template<typename T, template<class U> class Descriptor>
class GenericLatticeFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    GenericLatticeFunctional2D(OneCellFunctional2D<T,Descriptor>* f_);
    GenericLatticeFunctional2D(GenericLatticeFunctional2D<T,Descriptor> const& rhs);
    ~GenericLatticeFunctional2D();
    GenericLatticeFunctional2D<T,Descriptor>& operator=(GenericLatticeFunctional2D<T,Descriptor> const& rhs);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual GenericLatticeFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
private:
    OneCellFunctional2D<T,Descriptor>* f;
};

template<typename T, template<class U> class Descriptor>
class GenericIndexedLatticeFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    GenericIndexedLatticeFunctional2D(OneCellIndexedFunctional2D<T,Descriptor>* f_);
    GenericIndexedLatticeFunctional2D(GenericIndexedLatticeFunctional2D<T,Descriptor> const& rhs);
    ~GenericIndexedLatticeFunctional2D();
    GenericIndexedLatticeFunctional2D<T,Descriptor>& operator=(GenericIndexedLatticeFunctional2D<T,Descriptor> const& rhs);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual GenericIndexedLatticeFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
private:
    OneCellIndexedFunctional2D<T,Descriptor>* f;
};


/* *************** Class InstantiateDynamicsFunctional2D ************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateDynamicsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    InstantiateDynamicsFunctional2D(Dynamics<T,Descriptor>* dynamics_);
    InstantiateDynamicsFunctional2D(InstantiateDynamicsFunctional2D<T,Descriptor> const& rhs);
    InstantiateDynamicsFunctional2D<T,Descriptor>& operator= (
            InstantiateDynamicsFunctional2D<T,Descriptor> const& rhs );
    ~InstantiateDynamicsFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateDynamicsFunctional2D<T,Descriptor>* clone() const ;
private:
    Dynamics<T,Descriptor>* dynamics;
};

/* ************* Class InstantiateComplexDomainDynamicsFunctional2D ** */

template<typename T, template<typename U> class Descriptor>
class InstantiateComplexDomainDynamicsFunctional2D
    : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    InstantiateComplexDomainDynamicsFunctional2D( Dynamics<T,Descriptor>* dynamics_,
                                                  DomainFunctional2D* domain_ );
    InstantiateComplexDomainDynamicsFunctional2D (
            InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor> const& rhs );
    InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>& operator= (
            InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor> const& rhs );
    ~InstantiateComplexDomainDynamicsFunctional2D();
    virtual void process(Box2D boundingBox, BlockLattice2D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    DomainFunctional2D* domain;
};


/* ************* Class InstantiateDotDynamicsFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateDotDynamicsFunctional2D : public DotProcessingFunctional2D_L<T,Descriptor> {
public:
    InstantiateDotDynamicsFunctional2D(Dynamics<T,Descriptor>* dynamics_);
    InstantiateDotDynamicsFunctional2D(InstantiateDotDynamicsFunctional2D<T,Descriptor> const& rhs);
    InstantiateDotDynamicsFunctional2D<T,Descriptor>& operator= (
            InstantiateDotDynamicsFunctional2D<T,Descriptor> const& rhs );
    ~InstantiateDotDynamicsFunctional2D();
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateDotDynamicsFunctional2D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
};


/* ************* Class DynamicsFromMaskFunctional2D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
/** Note that the boolean mask is of type MultiScalarField2D<T> instead of MultiScalarField2D<bool>,
 *  because coupling lattices of different data types is not possible.
 */
template<typename T, template<typename U> class Descriptor>
class DynamicsFromMaskFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor> {
public:
    DynamicsFromMaskFunctional2D(Dynamics<T,Descriptor>* dynamics_, bool whichFlag_);
    DynamicsFromMaskFunctional2D(DynamicsFromMaskFunctional2D<T,Descriptor> const& rhs);
    DynamicsFromMaskFunctional2D<T,Descriptor>& operator= (
            DynamicsFromMaskFunctional2D<T,Descriptor> const& rhs );
    ~DynamicsFromMaskFunctional2D();
    virtual void process (
            Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                          ScalarField2D<T>& mask );
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DynamicsFromMaskFunctional2D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    bool whichFlag;
};


/* ************* Class DynamicsFromIntMaskFunctional2D ************************ */

/// Assign dynamics to nodes specified by a boolean mask.
/** Note that the boolean mask is of type MultiScalarField2D<T> instead of MultiScalarField2D<bool>,
 *  because coupling lattices of different data types is not possible.
 */
template<typename T, template<typename U> class Descriptor>
class DynamicsFromIntMaskFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor> {
public:
    DynamicsFromIntMaskFunctional2D(Dynamics<T,Descriptor>* dynamics_, int whichFlag_);
    DynamicsFromIntMaskFunctional2D(DynamicsFromIntMaskFunctional2D<T,Descriptor> const& rhs);
    DynamicsFromIntMaskFunctional2D<T,Descriptor>& operator= (
            DynamicsFromIntMaskFunctional2D<T,Descriptor> const& rhs );
    ~DynamicsFromIntMaskFunctional2D();
    virtual void process (
            Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                          ScalarField2D<T>& mask );
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DynamicsFromIntMaskFunctional2D<T,Descriptor>* clone() const;
private:
    Dynamics<T,Descriptor>* dynamics;
    int whichFlag;
};


/* ************* Class SetConstBoundaryVelocityFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class SetConstBoundaryVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    SetConstBoundaryVelocityFunctional2D(Array<T,Descriptor<T>::d> velocity);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual SetConstBoundaryVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    Array<T,Descriptor<T>::d> u;
};


/* ************* Class SetConstBoundaryDensityFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class SetConstBoundaryDensityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    SetConstBoundaryDensityFunctional2D(T rho_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual SetConstBoundaryDensityFunctional2D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T rho;
};


/* ************* Class IniConstEquilibriumFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class IniConstEquilibriumFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    IniConstEquilibriumFunctional2D(T density_, Array<T,Descriptor<T>::d> velocity);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual IniConstEquilibriumFunctional2D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    T jSqr;
};

/* ************* Class StripeOffDensityOffsetFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class StripeOffDensityOffsetFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    StripeOffDensityOffsetFunctional2D(T deltaRho_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual StripeOffDensityOffsetFunctional2D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual void rescale(T dxScale, T dtScale);
private:
    T deltaRho;
};


/* ************* Class InstantiateCompositeDynamicsFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
class InstantiateCompositeDynamicsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    InstantiateCompositeDynamicsFunctional2D(CompositeDynamics<T,Descriptor>* compositeDynamics_);
    InstantiateCompositeDynamicsFunctional2D(InstantiateCompositeDynamicsFunctional2D<T,Descriptor> const& rhs);
    InstantiateCompositeDynamicsFunctional2D<T,Descriptor>& operator= (
            InstantiateCompositeDynamicsFunctional2D<T,Descriptor> const& rhs );
    ~InstantiateCompositeDynamicsFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual InstantiateCompositeDynamicsFunctional2D<T,Descriptor>* clone() const;
private:
    CompositeDynamics<T,Descriptor>* compositeDynamics;
};

}  // namespace plb

#endif  // LATTICE_INITIALIZER_FUNCTIONALS_2D_H
