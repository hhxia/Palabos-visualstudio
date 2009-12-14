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
#ifndef DATA_PROCESSOR_WRAPPER_2D_H
#define DATA_PROCESSOR_WRAPPER_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/blockSurface2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include <vector>

namespace plb {

template<typename T> class Block2D;
template<typename T> class AtomicBlock2D;
template<typename T, template<typename U> class Descriptor> class BlockLatticeBase2D;
template<typename T, template<typename U> class Descriptor> class BlockLattice2D;
template<typename T> class ScalarFieldBase2D;
template<typename T> class ScalarField2D;
template<typename T, int nDim> class TensorFieldBase2D;
template<typename T, int nDim> class TensorField2D;


/* *************** All flavors of Box processing functionals ********* */

/// Easy instantiation of boxed data processor (general case)
template<typename T>
struct BoxProcessingFunctional2D {
    virtual ~BoxProcessingFunctional2D() { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoxProcessingFunctional2D<T>* clone() const =0;
};

/// A Boxed data processor, automatically generated from a BoxProcessingFunctional2D
template<typename T>
class BoxProcessor2D : public DataProcessor2D<T> {
public:
    BoxProcessor2D(BoxProcessingFunctional2D<T>* functional_,
                   Box2D domain_, std::vector<AtomicBlock2D<T>*> atomicBlocks_);
    BoxProcessor2D(BoxProcessor2D<T> const& rhs);
    BoxProcessor2D<T>& operator=(BoxProcessor2D<T> const& rhs);
    ~BoxProcessor2D();
    Box2D getDomain() const;
    virtual void process();
    virtual BoxProcessor2D<T>* clone() const;
private:
    BoxProcessingFunctional2D<T>* functional;
    Box2D domain;
    std::vector<AtomicBlock2D<T>*> atomicBlocks;
};

/// An automatically created generator for the BoxProcessor2D
template<typename T>
class BoxProcessorGenerator2D : public BoxedDataProcessorGenerator2D<T> {
public:
    BoxProcessorGenerator2D(BoxProcessingFunctional2D<T>* functional_, Box2D domain);
    ~BoxProcessorGenerator2D();
    BoxProcessorGenerator2D(BoxProcessorGenerator2D<T> const& rhs);
    BoxProcessorGenerator2D<T>& operator=(BoxProcessorGenerator2D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> atomicBlocks) const;
    virtual BoxProcessorGenerator2D<T>* clone() const;
private:
    BoxProcessingFunctional2D<T>* functional;
};


/// Easy instantiation of boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoxProcessingFunctional2D_L : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template<typename T>
struct BoxProcessingFunctional2D_S : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoxProcessingFunctional2D_T : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoxProcessingFunctional2D_LL : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor1>& lattice1,
                                       BlockLattice2D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoxProcessingFunctional2D_SS : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field1,
                                       ScalarField2D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoxProcessingFunctional2D_TT : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, TensorField2D<T,nDim1>& field1,
                                       TensorField2D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoxProcessingFunctional2D_ST : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field1,
                                       TensorField2D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoxProcessingFunctional2D_LS : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoxProcessingFunctional2D_LT : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices
template<typename T, template<typename U> class Descriptor>
struct LatticeBoxProcessingFunctional2D : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalarFields
template<typename T>
struct ScalarFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> scalarFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensorFields
template<typename T, int nDim>
struct TensorFieldBoxProcessingFunctional2D : public BoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> tensorFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional2D_L<T,Descriptor>* functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional2D_L<T,Descriptor>* functional,
                                   Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, plint level=0);

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D_S<T>* functional,
                               Box2D domain, ScalarFieldBase2D<T>& field);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D_S<T>* functional,
                                   Box2D domain, ScalarFieldBase2D<T>& field, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_T<T,nDim>* functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field);
template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_T<T,nDim>* functional,
                                   Box2D domain, TensorFieldBase2D<T,nDim>& field, plint level=0);



/* *************** All flavors of Dot processing functionals ********* */

/// Easy instantiation of dotted data processor (general case)
template<typename T>
struct DotProcessingFunctional2D {
    virtual ~DotProcessingFunctional2D() { }
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DotProcessingFunctional2D<T>* clone() const =0;
};

/// A Dotted data processor, automatically generated from a DotProcessingFunctional2D
template<typename T>
class DotProcessor2D : public DataProcessor2D<T> {
public:
    DotProcessor2D(DotProcessingFunctional2D<T>* functional_,
                   DotList2D const& dotList_, std::vector<AtomicBlock2D<T>*> atomicBlocks_);
    DotProcessor2D(DotProcessor2D<T> const& rhs);
    DotProcessor2D<T>& operator=(DotProcessor2D<T> const& rhs);
    ~DotProcessor2D();
    virtual void process();
    virtual DotProcessor2D<T>* clone() const;
    DotList2D const& getDotList() const;
private:
    DotProcessingFunctional2D<T>* functional;
    DotList2D dotList;
    std::vector<AtomicBlock2D<T>*> atomicBlocks;
};

/// An automatically created generator for the DotProcessor2D
template<typename T>
class DotProcessorGenerator2D : public DottedDataProcessorGenerator2D<T> {
public:
    DotProcessorGenerator2D(DotProcessingFunctional2D<T>* functional_, DotList2D const& dotList);
    ~DotProcessorGenerator2D();
    DotProcessorGenerator2D(DotProcessorGenerator2D<T> const& rhs);
    DotProcessorGenerator2D<T>& operator=(DotProcessorGenerator2D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> atomicBlocks) const;
    virtual DotProcessorGenerator2D<T>* clone() const;
private:
    DotProcessingFunctional2D<T>* functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct DotProcessingFunctional2D_L : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template<typename T>
struct DotProcessingFunctional2D_S : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template<typename T, int nDim>
struct DotProcessingFunctional2D_T : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct DotProcessingFunctional2D_LL : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor1>& lattice1,
                                                   BlockLattice2D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
struct DotProcessingFunctional2D_SS : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field1,
                                                   ScalarField2D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct DotProcessingFunctional2D_TT : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, TensorField2D<T,nDim1>& field1,
                                                   TensorField2D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct DotProcessingFunctional2D_ST : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field1,
                                                   TensorField2D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct DotProcessingFunctional2D_LS : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice,
                         ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct DotProcessingFunctional2D_LT : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice,
                                                   TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices
template<typename T, template<typename U> class Descriptor>
struct LatticeDotProcessingFunctional2D : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, std::vector<BlockLattice2D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalarFields
template<typename T>
struct ScalarFieldDotProcessingFunctional2D : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, std::vector<ScalarField2D<T>*> scalarFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensorFields
template<typename T, int nDim>
struct TensorFieldDotProcessingFunctional2D : public DotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, std::vector<TensorField2D<T,nDim>*> tensorFields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional2D_L<T,Descriptor>* functional,
                               DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional2D_L<T,Descriptor>* functional,
                                   DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice, plint level=0);

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D_S<T>* functional,
                               DotList2D const& dotList, ScalarFieldBase2D<T>& field);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D_S<T>* functional,
                                   DotList2D const& dotList, ScalarFieldBase2D<T>& field, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_T<T,nDim>* functional,
                               DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field);
template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_T<T,nDim>* functional,
                                   DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field, plint level=0);


/* *************** All flavors of Bounded Box processing functionals ********* */

/// Easy instantiation of boxed data processor special boundary treatment (general case)
template<typename T>
class BoundedBoxProcessingFunctional2D {
public:
    virtual ~BoundedBoxProcessingFunctional2D() { }
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks ) =0;
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks ) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoundedBoxProcessingFunctional2D<T>* clone() const =0;
    BoxProcessingFunctional2D<T>* getBulkProcessor() const; 
    BoxProcessingFunctional2D<T>* getEdgeProcessor(int direction, int orientation) const; 
    BoxProcessingFunctional2D<T>* getCornerProcessor(int normalX, int normalY) const; 
public:
    class BulkWrapperFunctional : public BoxProcessingFunctional2D<T> {
    public:
        BulkWrapperFunctional(BoundedBoxProcessingFunctional2D<T>* boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const& rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional& operator=(BulkWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual BulkWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional;
    };
    class EdgeWrapperFunctional : public BoxProcessingFunctional2D<T> {
    public:
        EdgeWrapperFunctional(BoundedBoxProcessingFunctional2D<T>* boundedFunctional_, int direction_, int orientation_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const& rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional& operator=(EdgeWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual EdgeWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional;
        int direction, orientation;
    };
    class CornerWrapperFunctional : public BoxProcessingFunctional2D<T> {
    public:
        CornerWrapperFunctional(BoundedBoxProcessingFunctional2D<T>* boundedFunctional_, int normalX_, int normalY_);
        CornerWrapperFunctional(CornerWrapperFunctional const& rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional& operator=(CornerWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual CornerWrapperFunctional* clone() const;
    private:
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional;
        int normalX, normalY;
    };
};

/// Generic implementation of "apply" and "integrate" for Bounded BoxProcessingFunctionals
template<typename T>
class BoundedOneBlockProcessingFunctionalOperation2D {
public:
    BoundedOneBlockProcessingFunctionalOperation2D(Box2D const& domain, plint boundaryWidth_);
    void apply(BoundedBoxProcessingFunctional2D<T>* functional, Block2D<T>& block);
    void integrate(BoundedBoxProcessingFunctional2D<T>* functional, Block2D<T>& block, plint level);
private:
    BlockSurface2D surf;
};

/// Easy instantiation of bounded boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoundedBoxProcessingFunctional2D_L : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              BlockLattice2D<T,Descriptor>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                BlockLattice2D<T,Descriptor>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for a single scalar field
template<typename T>
struct BoundedBoxProcessingFunctional2D_S : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, ScalarField2D<T>& lattice) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoundedBoxProcessingFunctional2D_T : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, TensorField2D<T,nDim>& lattice) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              TensorField2D<T,nDim>& lattice ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                TensorField2D<T,nDim>& lattice ) =0;
};

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoundedBoxProcessingFunctional2D_LL : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, BlockLattice2D<T,Descriptor1>& lattice1,
                                            BlockLattice2D<T,Descriptor2>& lattice2 ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              BlockLattice2D<T,Descriptor1>& lattice1,
                              BlockLattice2D<T,Descriptor2>& lattice2 ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                BlockLattice2D<T,Descriptor1>& lattice1,
                                BlockLattice2D<T,Descriptor2>& lattice2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoundedBoxProcessingFunctional2D_SS : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, ScalarField2D<T>& field1,
                                            ScalarField2D<T>& field2 ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& field1,
                              ScalarField2D<T>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& field1,
                                ScalarField2D<T>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoundedBoxProcessingFunctional2D_TT : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, TensorField2D<T,nDim1>& field1,
                                            TensorField2D<T,nDim2>& field2 ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              TensorField2D<T,nDim1>& field1,
                              TensorField2D<T,nDim2>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                TensorField2D<T,nDim1>& field1,
                                TensorField2D<T,nDim2>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoundedBoxProcessingFunctional2D_ST : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, ScalarField2D<T>& field1,
                                            TensorField2D<T,nDim>& field2 ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& field1,
                              TensorField2D<T,nDim>& field2 ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& field1,
                                TensorField2D<T,nDim>& field2 ) =0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoundedBoxProcessingFunctional2D_LS : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                            ScalarField2D<T>& field ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              BlockLattice2D<T,Descriptor>& lattice,
                              ScalarField2D<T>& field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                BlockLattice2D<T,Descriptor>& lattice,
                                ScalarField2D<T>& field ) =0;
};

/// Easy instantiation of bounded boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoundedBoxProcessingFunctional2D_LT : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk( Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                            TensorField2D<T,nDim>& field ) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              BlockLattice2D<T,Descriptor>& lattice,
                              TensorField2D<T,nDim>& field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                BlockLattice2D<T,Descriptor>& lattice,
                                TensorField2D<T,nDim>& field ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct BoundedLatticeBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              std::vector<BlockLattice2D<T,Descriptor>*> lattices ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                std::vector<BlockLattice2D<T,Descriptor>*> lattices ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple scalar fields.
template<typename T>
struct BoundedScalarFieldBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T>*> field) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              std::vector<ScalarField2D<T>*> field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                std::vector<ScalarField2D<T>*> field ) =0;
};

/// Easy instantiation of bounded boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct BoundedTensorFieldBoxProcessingFunctional2D : public BoundedBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, std::vector<TensorField2D<T,nDim>*> field) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              std::vector<TensorField2D<T,nDim>*> field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                std::vector<TensorField2D<T,nDim>*> field ) =0;
};



template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_L<T,Descriptor>* functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_L<T,Descriptor>* functional,
                                   Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_S<T>* functional,
                               Box2D domain, ScalarFieldBase2D<T>& field, plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_S<T>* functional,
                                   Box2D domain, ScalarFieldBase2D<T>& field,
                                   plint boundaryWidth, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_T<T,nDim>* functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field, plint boundaryWidth);
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_T<T,nDim>* functional,
                                   Box2D domain, TensorFieldBase2D<T,nDim>& field,
                                   plint boundaryWidth, plint level=0);


}  // namespace plb

#endif  // DATA_PROCESSOR_WRAPPER_2D_H
