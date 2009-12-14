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
#ifndef REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_H
#define REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_H

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


/* *************** All flavors of reductive box processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
template<typename T>
class ReductiveBoxProcessingFunctional2D {
public:
    virtual ~ReductiveBoxProcessingFunctional2D() { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual ReductiveBoxProcessingFunctional2D<T>* clone() const =0;
    virtual BlockStatistics<T> const& getStatistics() const =0;
    virtual BlockStatistics<T>& getStatistics() =0;
};

/// ReductiveBoxProcessingFunctional2D which instantiates its own statistics object
template<typename T>
class PlainReductiveBoxProcessingFunctional2D : public ReductiveBoxProcessingFunctional2D<T> {
public:
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
};

/// A reductive boxed data processor, automatically generated from a ReductiveBoxProcessingFunctional2D
template<typename T>
class ReductiveBoxProcessor2D : public DataProcessor2D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessor2D(ReductiveBoxProcessingFunctional2D<T>* functional_,
                            Box2D domain_, std::vector<AtomicBlock2D<T>*> atomicBlocks_);
    Box2D getDomain() const;
    virtual void process();
    virtual ReductiveBoxProcessor2D<T>* clone() const;
private:
    ReductiveBoxProcessingFunctional2D<T>* functional;
    Box2D domain;
    std::vector<AtomicBlock2D<T>*> atomicBlocks;
};

/// An automatically created generator for the ReductiveBoxProcessor2D
template<typename T>
class ReductiveBoxProcessorGenerator2D : public BoxedReductiveDataProcessorGenerator2D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessorGenerator2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessorGenerator2D(ReductiveBoxProcessingFunctional2D<T>* functional_, Box2D domain);
    ~ReductiveBoxProcessorGenerator2D();
    ReductiveBoxProcessorGenerator2D(ReductiveBoxProcessorGenerator2D<T> const& rhs);
    ReductiveBoxProcessorGenerator2D<T>& operator=(ReductiveBoxProcessorGenerator2D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual ReductiveBoxProcessorGenerator2D<T>* clone() const;
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
    ReductiveBoxProcessingFunctional2D<T> const& getFunctional() const;
private:
    ReductiveBoxProcessingFunctional2D<T>* functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional2D_L : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template<typename T>
struct ReductiveBoxProcessingFunctional2D_S : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template<typename T, int nDim>
struct ReductiveBoxProcessingFunctional2D_T : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct ReductiveBoxProcessingFunctional2D_LL : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor1>& lattice1,
                                       BlockLattice2D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct ReductiveBoxProcessingFunctional2D_SS : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field1,
                                       ScalarField2D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct ReductiveBoxProcessingFunctional2D_TT : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, TensorField2D<T,nDim1>& field1,
                                       TensorField2D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct ReductiveBoxProcessingFunctional2D_ST : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, ScalarField2D<T>& field1,
                                       TensorField2D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional2D_LS : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       ScalarField2D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct ReductiveBoxProcessingFunctional2D_LT : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       TensorField2D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct ReductiveLatticeBoxProcessingFunctional2D : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalar fields.
template<typename T>
struct ReductiveScalarFieldBoxProcessingFunctional2D : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct ReductiveTensorFieldBoxProcessingFunctional2D : public PlainReductiveBoxProcessingFunctional2D<T> {
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};



template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_L<T,Descriptor>& functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_S<T>& functional,
                               Box2D domain, ScalarFieldBase2D<T>& field);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_T<T,nDim>& functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field);


/* *************** All flavors of reductive dot processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
template<typename T>
class ReductiveDotProcessingFunctional2D {
public:
    virtual ~ReductiveDotProcessingFunctional2D() { }
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual ReductiveDotProcessingFunctional2D<T>* clone() const =0;
    virtual BlockStatistics<T> const& getStatistics() const =0;
    virtual BlockStatistics<T>& getStatistics() =0;
};

/// ReductiveDotProcessingFunctional2D which instantiates its own statistics object
template<typename T>
class PlainReductiveDotProcessingFunctional2D : public ReductiveDotProcessingFunctional2D<T> {
public:
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
};


/// A ReductiveDotted data processor, automatically generated from a ReductiveDotProcessingFunctional2D
template<typename T>
class ReductiveDotProcessor2D : public DataProcessor2D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor2D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveDotProcessor2D(ReductiveDotProcessingFunctional2D<T>* functional_,
                            DotList2D const& dotList_, std::vector<AtomicBlock2D<T>*> atomicBlocks_);
    DotList2D const& getDotList() const;
    virtual void process();
    virtual ReductiveDotProcessor2D<T>* clone() const;
private:
    ReductiveDotProcessingFunctional2D<T>* functional;
    DotList2D dotList;
    std::vector<AtomicBlock2D<T>*> atomicBlocks;
};

/// An automatically created generator for the ReductiveDotProcessor2D
template<typename T>
class ReductiveDotProcessorGenerator2D : public DottedReductiveDataProcessorGenerator2D<T> {
public:
    ReductiveDotProcessorGenerator2D(ReductiveDotProcessingFunctional2D<T>* functional_, DotList2D const& dotList);
    ~ReductiveDotProcessorGenerator2D();
    ReductiveDotProcessorGenerator2D(ReductiveDotProcessorGenerator2D<T> const& rhs);
    ReductiveDotProcessorGenerator2D<T>& operator=(ReductiveDotProcessorGenerator2D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual ReductiveDotProcessorGenerator2D<T>* clone() const;
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
    ReductiveDotProcessingFunctional2D<T> const& getFunctional() const;
private:
    ReductiveDotProcessingFunctional2D<T>* functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct ReductiveDotProcessingFunctional2D_L : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice) =0;
    virtual ReductiveDotProcessingFunctional2D_L<T,Descriptor>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template<typename T>
struct ReductiveDotProcessingFunctional2D_S : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field) =0;
    virtual ReductiveDotProcessingFunctional2D_S<T>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template<typename T, int nDim>
struct ReductiveDotProcessingFunctional2D_T : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, TensorField2D<T,nDim>& field) =0;
    virtual ReductiveDotProcessingFunctional2D_T<T,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct ReductiveDotProcessingFunctional2D_LL : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor1>& lattice1,
                                                   BlockLattice2D<T,Descriptor2>& lattice2) =0;
    virtual ReductiveDotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
struct ReductiveDotProcessingFunctional2D_SS : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field1,
                                                   ScalarField2D<T>& field2) =0;
    virtual ReductiveDotProcessingFunctional2D_SS<T>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct ReductiveDotProcessingFunctional2D_TT : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, TensorField2D<T,nDim1>& field1,
                                                   TensorField2D<T,nDim2>& field2) =0;
    virtual ReductiveDotProcessingFunctional2D_TT<T,nDim1,nDim2>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct ReductiveDotProcessingFunctional2D_ST : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, ScalarField2D<T>& field1,
                                                   TensorField2D<T,nDim>& field2) =0;
    virtual ReductiveDotProcessingFunctional2D_ST<T,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct ReductiveDotProcessingFunctional2D_LS : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice,
                         ScalarField2D<T>& field) =0;
    virtual ReductiveDotProcessingFunctional2D_LS<T,Descriptor>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct ReductiveDotProcessingFunctional2D_LT : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice,
                                                   TensorField2D<T,nDim>& field) =0;
    virtual ReductiveDotProcessingFunctional2D_LT<T,Descriptor,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct ReductiveLatticeDotProcessingFunctional2D : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotlist, std::vector<BlockLattice2D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotlist, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalar fields.
template<typename T>
struct ReductiveScalarFieldDotProcessingFunctional2D : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotlist, std::vector<ScalarField2D<T>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotlist, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensor fields.
template<typename T, int nDim>
struct ReductiveTensorFieldDotProcessingFunctional2D : public PlainReductiveDotProcessingFunctional2D<T> {
    virtual void process(DotList2D const& dotlist, std::vector<TensorField2D<T,nDim>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList2D const& dotlist, std::vector<AtomicBlock2D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_L<T,Descriptor>& functional,
                               DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice);
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_S<T>& functional,
                               DotList2D const& dotList, ScalarFieldBase2D<T>& field);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_T<T,nDim>& functional,
                               DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field);

/* *************** All flavors of Bounded Reductive Box processing functionals ********* */

/// Easy instantiation of bounded reductive boxed data processor (general case)
template<typename T>
class BoundedReductiveBoxProcessingFunctional2D {
public:
    virtual ~BoundedReductiveBoxProcessingFunctional2D() { }
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoundedReductiveBoxProcessingFunctional2D<T>* clone() const =0;
    ReductiveBoxProcessingFunctional2D<T>* getBulkProcessor() const; 
    ReductiveBoxProcessingFunctional2D<T>* getEdgeProcessor(int direction, int orientation) const; 
    ReductiveBoxProcessingFunctional2D<T>* getCornerProcessor(int normalX, int normalY) const; 
    BlockStatistics<T> const& getStatistics() const;
    BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
public:
    class BulkWrapperFunctional : public ReductiveBoxProcessingFunctional2D<T> {
    public:
        BulkWrapperFunctional(BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const& rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional& operator=(BulkWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual BulkWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional;
    };
    class EdgeWrapperFunctional : public ReductiveBoxProcessingFunctional2D<T> {
    public:
        EdgeWrapperFunctional(BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_,
                              int direction_, int orientation_);
        EdgeWrapperFunctional(EdgeWrapperFunctional const& rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional& operator=(EdgeWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual EdgeWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional;
        int direction, orientation;
    };
    class CornerWrapperFunctional : public ReductiveBoxProcessingFunctional2D<T> {
    public:
        CornerWrapperFunctional(BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_,
                                int normalX_, int normalY_);
        CornerWrapperFunctional(CornerWrapperFunctional const& rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional& operator=(CornerWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual CornerWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional;
        int normalX, normalY;
    };
};

/// Generic implementatio1 of "apply" and "integrate" for Bounded Reductive BoxProcessingFunctionals
template<typename T>
class BoundedReductiveOneBlockProcessingFunctionalOperation2D {
public:
    BoundedReductiveOneBlockProcessingFunctionalOperation2D(Box2D const& domain, plint boundaryWidth_);
    void apply(BoundedReductiveBoxProcessingFunctional2D<T>& functional, Block2D<T>& block);
private:
    BlockSurface2D surf;
};

/// Easy instantiation of bounded reductive boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional2D_L : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for a single scalar field
template<typename T>
struct BoundedReductiveBoxProcessingFunctional2D_S : public BoundedReductiveBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, ScalarField2D<T>& field) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& field ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_T : public BoundedReductiveBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, TensorField2D<T,nDim>& field) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              TensorField2D<T,nDim>& field ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                TensorField2D<T,nDim>& field ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoundedReductiveBoxProcessingFunctional2D_LL : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoundedReductiveBoxProcessingFunctional2D_SS : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoundedReductiveBoxProcessingFunctional2D_TT : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_ST : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional2D_LS : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoundedReductiveBoxProcessingFunctional2D_LT : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveLatticeBoxProcessingFunctional2D : public BoundedReductiveBoxProcessingFunctional2D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for multiple scalar fields.
template<typename T>
struct BoundedReductiveScalarFieldBoxProcessingFunctional2D : public BoundedReductiveBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T>*> fields) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              std::vector<ScalarField2D<T>*> fields ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                std::vector<ScalarField2D<T>*> fields ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct BoundedReductiveTensorFieldBoxProcessingFunctional2D : public BoundedReductiveBoxProcessingFunctional2D<T> {
    virtual void processBulkGeneric(Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);
    virtual void processEdgeGeneric( int direction, int orientation, Box2D domain,
                                     std::vector<AtomicBlock2D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, Box2D domain,
                                       std::vector<AtomicBlock2D<T>*> atomicBlocks );

    virtual void processBulk(Box2D domain, std::vector<TensorField2D<T,nDim>*> fields) =0;
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              std::vector<TensorField2D<T,nDim>*> fields ) =0;
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                std::vector<TensorField2D<T,nDim>*> fields ) =0;
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_L<T,Descriptor>& functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_S<T>& functional,
                               Box2D domain, ScalarFieldBase2D<T>& field, plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_T<T,nDim>& functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field, plint boundaryWidth);



}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_H
