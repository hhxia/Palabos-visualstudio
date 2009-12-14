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
 * Utilities to help users handle data processors -- header file.
 */
#ifndef DATA_PROCESSOR_WRAPPER_3D_H
#define DATA_PROCESSOR_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "core/blockSurface3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include <vector>

namespace plb {

template<typename T> class Block3D;
template<typename T> class AtomicBlock3D;
template<typename T, template<typename U> class Descriptor> class BlockLatticeBase3D;
template<typename T, template<typename U> class Descriptor> class BlockLattice3D;
template<typename T> class ScalarFieldBase3D;
template<typename T> class ScalarField3D;
template<typename T, int nDim> class TensorFieldBase3D;
template<typename T, int nDim> class TensorField3D;


/* *************** All flavors of Box processing functionals ********* */

/// Easy instantiation of boxed data processor (general case)
template<typename T>
struct BoxProcessingFunctional3D {
    virtual ~BoxProcessingFunctional3D() { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoxProcessingFunctional3D<T>* clone() const =0;
};

/// A Boxed data processor, automatically generated from a BoxProcessingFunctional3D
template<typename T>
class BoxProcessor3D : public DataProcessor3D<T> {
public:
    BoxProcessor3D(BoxProcessingFunctional3D<T>* functional_,
                   Box3D domain_, std::vector<AtomicBlock3D<T>*> atomicBlocks_);
    BoxProcessor3D(BoxProcessor3D<T> const& rhs);
    BoxProcessor3D<T>& operator=(BoxProcessor3D<T> const& rhs);
    ~BoxProcessor3D();
    Box3D getDomain() const;
    virtual void process();
    virtual BoxProcessor3D<T>* clone() const;
private:
    BoxProcessingFunctional3D<T>* functional;
    Box3D domain;
    std::vector<AtomicBlock3D<T>*> atomicBlocks;
};

/// An automatically created generator for the BoxProcessor3D
template<typename T>
class BoxProcessorGenerator3D : public BoxedDataProcessorGenerator3D<T> {
public:
    BoxProcessorGenerator3D(BoxProcessingFunctional3D<T>* functional_, Box3D domain);
    ~BoxProcessorGenerator3D();
    BoxProcessorGenerator3D(BoxProcessorGenerator3D<T> const& rhs);
    BoxProcessorGenerator3D<T>& operator=(BoxProcessorGenerator3D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> atomicBlocks) const;
    virtual BoxProcessorGenerator3D<T>* clone() const;
private:
    BoxProcessingFunctional3D<T>* functional;
};


/// Easy instantiation of boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoxProcessingFunctional3D_L : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template<typename T>
struct BoxProcessingFunctional3D_S : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoxProcessingFunctional3D_T : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoxProcessingFunctional3D_LL : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor1>& lattice1,
                                       BlockLattice3D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoxProcessingFunctional3D_SS : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field1,
                                       ScalarField3D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoxProcessingFunctional3D_TT : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, TensorField3D<T,nDim1>& field1,
                                       TensorField3D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoxProcessingFunctional3D_ST : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field1,
                                       TensorField3D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoxProcessingFunctional3D_LS : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoxProcessingFunctional3D_LT : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template<typename T, template<typename U> class Descriptor>
struct LatticeBoxProcessingFunctional3D : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template<typename T>
struct ScalarFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template<typename T, int nDim>
struct TensorFieldBoxProcessingFunctional3D : public BoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional3D_L<T,Descriptor>* functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional3D_L<T,Descriptor>* functional,
                                   Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice, plint level=0);

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D_S<T>* functional,
                               Box3D domain, ScalarFieldBase3D<T>& field);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D_S<T>* functional,
                                   Box3D domain, ScalarFieldBase3D<T>& field, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_T<T,nDim>* functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field);
template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_T<T,nDim>* functional,
                                   Box3D domain, TensorFieldBase3D<T,nDim>& field, plint level=0);



/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of dotted data processor (general case)
template<typename T>
struct DotProcessingFunctional3D {
    virtual ~DotProcessingFunctional3D() { }
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DotProcessingFunctional3D<T>* clone() const =0;
};

/// A Dotted data processor, automatically generated from a DotProcessingFunctional3D
template<typename T>
class DotProcessor3D : public DataProcessor3D<T> {
public:
    DotProcessor3D(DotProcessingFunctional3D<T>* functional_,
                   DotList3D const& dotList_, std::vector<AtomicBlock3D<T>*> atomicBlocks_);
    DotProcessor3D(DotProcessor3D<T> const& rhs);
    DotProcessor3D<T>& operator=(DotProcessor3D<T> const& rhs);
    ~DotProcessor3D();
    virtual void process();
    virtual DotProcessor3D<T>* clone() const;
    DotList3D const& getDotList() const;
private:
    DotProcessingFunctional3D<T>* functional;
    DotList3D dotList;
    std::vector<AtomicBlock3D<T>*> atomicBlocks;
};

/// An automatically created generator for the DotProcessor3D
template<typename T>
class DotProcessorGenerator3D : public DottedDataProcessorGenerator3D<T> {
public:
    DotProcessorGenerator3D(DotProcessingFunctional3D<T>* functional_, DotList3D const& dotList);
    ~DotProcessorGenerator3D();
    DotProcessorGenerator3D(DotProcessorGenerator3D<T> const& rhs);
    DotProcessorGenerator3D<T>& operator=(DotProcessorGenerator3D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> atomicBlocks) const;
    virtual DotProcessorGenerator3D<T>* clone() const;
private:
    DotProcessingFunctional3D<T>* functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct DotProcessingFunctional3D_L : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template<typename T>
struct DotProcessingFunctional3D_S : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template<typename T, int nDim>
struct DotProcessingFunctional3D_T : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct DotProcessingFunctional3D_LL : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor1>& lattice1,
                                                   BlockLattice3D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
struct DotProcessingFunctional3D_SS : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field1,
                                                   ScalarField3D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct DotProcessingFunctional3D_TT : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, TensorField3D<T,nDim1>& field1,
                                                   TensorField3D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct DotProcessingFunctional3D_ST : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field1,
                                                   TensorField3D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct DotProcessingFunctional3D_LS : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice,
                         ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct DotProcessingFunctional3D_LT : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice,
                                                   TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template<typename T, template<typename U> class Descriptor>
struct LatticeDotProcessingFunctional3D : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, std::vector<BlockLattice3D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template<typename T>
struct ScalarFieldDotProcessingFunctional3D : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, std::vector<ScalarField3D<T>*> scalarFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template<typename T, int nDim>
struct TensorFieldDotProcessingFunctional3D : public DotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, std::vector<TensorField3D<T,nDim>*> tensorFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional3D_L<T,Descriptor>* functional,
                               DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional3D_L<T,Descriptor>* functional,
                                   DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice, plint level=0);

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D_S<T>* functional,
                               DotList3D const& dotList, ScalarFieldBase3D<T>& field);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D_S<T>* functional,
                                   DotList3D const& dotList, ScalarFieldBase3D<T>& field, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_T<T,nDim>* functional,
                               DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field);
template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_T<T,nDim>* functional,
                                   DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field, plint level=0);


/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed data processor special boundary treatment (general case)
template<typename T>
class BoundedBoxProcessingFunctional3D {
public:
    virtual ~BoundedBoxProcessingFunctional3D() { }
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks ) =0;
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks ) =0;
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks ) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoundedBoxProcessingFunctional3D<T>* clone() const =0;
    BoxProcessingFunctional3D<T>* getBulkProcessor() const; 
    BoxProcessingFunctional3D<T>* getPlaneProcessor(int direction, int orientation) const; 
    BoxProcessingFunctional3D<T>* getEdgeProcessor(int plane, int normal1, int normal2) const; 
    BoxProcessingFunctional3D<T>* getCornerProcessor(int normalX, int normalY, int normalZ) const; 
public:
    class BulkWrapperFunctional : public BoxProcessingFunctional3D<T> {
    public:
        BulkWrapperFunctional(BoundedBoxProcessingFunctional3D<T>* boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const& rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional& operator=(BulkWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual BulkWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional;
    };
    class PlaneWrapperFunctional : public BoxProcessingFunctional3D<T> {
    public:
        PlaneWrapperFunctional(BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
                               int direction_, int orientation_);
        PlaneWrapperFunctional(PlaneWrapperFunctional const& rhs);
        ~PlaneWrapperFunctional();
        PlaneWrapperFunctional& operator=(PlaneWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual PlaneWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional;
        int direction, orientation;
    };
    class EdgeWrapperFunctional : public BoxProcessingFunctional3D<T> {
    public:
        EdgeWrapperFunctional(BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
                              int plane_, int normal1_, int normal2_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const& rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional& operator=(EdgeWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual EdgeWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional;
        int plane, normal1, normal2;
    };
    class CornerWrapperFunctional : public BoxProcessingFunctional3D<T> {
    public:
        CornerWrapperFunctional(BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
                                int normalX_, int normalY_, int normalZ_);
        CornerWrapperFunctional(CornerWrapperFunctional const& rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional& operator=(CornerWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual CornerWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional;
        int normalX, normalY, normalZ;
    };
};

/// Generic implementation of "apply" and "integrate" for Bounded BoxProcessingFunctionals
template<typename T>
class BoundedOneBlockProcessingFunctionalOperation3D {
public:
    BoundedOneBlockProcessingFunctionalOperation3D(Box3D const& domain, plint boundaryWidth_);
    void apply(BoundedBoxProcessingFunctional3D<T>* functional, Block3D<T>& block);
    void integrate(BoundedBoxProcessingFunctional3D<T>* functional, Block3D<T>& block, plint level);
private:
    BlockSurface3D surf;
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoundedBoxProcessingFunctional3D_L : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              BlockLattice3D<T,Descriptor>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                BlockLattice3D<T,Descriptor>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for a single scalar field
template<typename T>
struct BoundedBoxProcessingFunctional3D_S : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, ScalarField3D<T>& lattice) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& lattice ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              ScalarField3D<T>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                ScalarField3D<T>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoundedBoxProcessingFunctional3D_T : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, TensorField3D<T,nDim>& lattice) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& lattice ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoundedBoxProcessingFunctional3D_LL : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, BlockLattice3D<T,Descriptor1>& lattice1,
                                            BlockLattice3D<T,Descriptor2>& lattice2 ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2 ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              BlockLattice3D<T,Descriptor1>& lattice1,
                              BlockLattice3D<T,Descriptor2>& lattice2 ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                BlockLattice3D<T,Descriptor1>& lattice1,
                                BlockLattice3D<T,Descriptor2>& lattice2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoundedBoxProcessingFunctional3D_SS : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, ScalarField3D<T>& field1,
                                            ScalarField3D<T>& field2 ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2 ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              ScalarField3D<T>& field1,
                              ScalarField3D<T>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                ScalarField3D<T>& field1,
                                ScalarField3D<T>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoundedBoxProcessingFunctional3D_TT : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, TensorField3D<T,nDim1>& field1,
                                            TensorField3D<T,nDim2>& field2 ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2 ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim1>& field1,
                              TensorField3D<T,nDim2>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim1>& field1,
                                TensorField3D<T,nDim2>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoundedBoxProcessingFunctional3D_ST : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, ScalarField3D<T>& field1,
                                            TensorField3D<T,nDim>& field2 ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2 ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              ScalarField3D<T>& field1,
                              TensorField3D<T,nDim>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                ScalarField3D<T>& field1,
                                TensorField3D<T,nDim>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoundedBoxProcessingFunctional3D_LS : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                            ScalarField3D<T>& field ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              BlockLattice3D<T,Descriptor>& lattice,
                              ScalarField3D<T>& field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                BlockLattice3D<T,Descriptor>& lattice,
                                ScalarField3D<T>& field ) =0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoundedBoxProcessingFunctional3D_LT : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk( Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                            TensorField3D<T,nDim>& field ) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              BlockLattice3D<T,Descriptor>& lattice,
                              TensorField3D<T,nDim>& field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                BlockLattice3D<T,Descriptor>& lattice,
                                TensorField3D<T,nDim>& field ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct BoundedLatticeBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              std::vector<BlockLattice3D<T,Descriptor>*> lattices ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                std::vector<BlockLattice3D<T,Descriptor>*> lattices ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple scalar fields.
template<typename T>
struct BoundedScalarFieldBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, std::vector<ScalarField3D<T>*> field) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               std::vector<ScalarField3D<T>*> field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              std::vector<ScalarField3D<T>*> field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                std::vector<ScalarField3D<T>*> field ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct BoundedTensorFieldBoxProcessingFunctional3D : public BoundedBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, std::vector<TensorField3D<T,nDim>*> field) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              std::vector<TensorField3D<T,nDim>*> field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                std::vector<TensorField3D<T,nDim>*> field ) =0;
};



template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_L<T,Descriptor>* functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_L<T,Descriptor>* functional,
                                   Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_S<T>* functional,
                               Box3D domain, ScalarFieldBase3D<T>& field, plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_S<T>* functional,
                                   Box3D domain, ScalarFieldBase3D<T>& field,
                                   plint boundaryWidth, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_T<T,nDim>* functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field, plint boundaryWidth);
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_T<T,nDim>* functional,
                                   Box3D domain, TensorFieldBase3D<T,nDim>& field,
                                   plint boundaryWidth, plint level=0);


}  // namespace plb

#endif  // DATA_PROCESSOR_WRAPPER_3D_H
