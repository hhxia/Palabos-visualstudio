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
#ifndef REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_H
#define REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_H

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


/* *************** All flavors of reductive box processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
template<typename T>
class ReductiveBoxProcessingFunctional3D {
public:
    virtual ~ReductiveBoxProcessingFunctional3D() { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual ReductiveBoxProcessingFunctional3D<T>* clone() const =0;
    virtual BlockStatistics<T> const& getStatistics() const =0;
    virtual BlockStatistics<T>& getStatistics() =0;
};

/// ReductiveBoxProcessingFunctional3D which instantiates its own statistics object
template<typename T>
class PlainReductiveBoxProcessingFunctional3D : public ReductiveBoxProcessingFunctional3D<T> {
public:
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
};

/// A reductive boxed data processor, automatically generated from a ReductiveBoxProcessingFunctional3D
template<typename T>
class ReductiveBoxProcessor3D : public DataProcessor3D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessor3D(ReductiveBoxProcessingFunctional3D<T>* functional_,
                            Box3D domain_, std::vector<AtomicBlock3D<T>*> atomicBlocks_);
    Box3D getDomain() const;
    virtual void process();
    virtual ReductiveBoxProcessor3D<T>* clone() const;
private:
    ReductiveBoxProcessingFunctional3D<T>* functional;
    Box3D domain;
    std::vector<AtomicBlock3D<T>*> atomicBlocks;
};

/// An automatically created generator for the ReductiveBoxProcessor3D
template<typename T>
class ReductiveBoxProcessorGenerator3D : public BoxedReductiveDataProcessorGenerator3D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessorGenerator3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveBoxProcessorGenerator3D(ReductiveBoxProcessingFunctional3D<T>* functional_, Box3D domain);
    ~ReductiveBoxProcessorGenerator3D();
    ReductiveBoxProcessorGenerator3D(ReductiveBoxProcessorGenerator3D<T> const& rhs);
    ReductiveBoxProcessorGenerator3D<T>& operator=(ReductiveBoxProcessorGenerator3D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual ReductiveBoxProcessorGenerator3D<T>* clone() const;
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
    ReductiveBoxProcessingFunctional3D<T> const& getFunctional() const;
private:
    ReductiveBoxProcessingFunctional3D<T>* functional;
};

/// Easy instantiation of boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional3D_L : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single scalar field
template<typename T>
struct ReductiveBoxProcessingFunctional3D_S : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for a single tensor field
template<typename T, int nDim>
struct ReductiveBoxProcessingFunctional3D_T : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct ReductiveBoxProcessingFunctional3D_LL : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor1>& lattice1,
                                       BlockLattice3D<T,Descriptor2>& lattice2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct ReductiveBoxProcessingFunctional3D_SS : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field1,
                                       ScalarField3D<T>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct ReductiveBoxProcessingFunctional3D_TT : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, TensorField3D<T,nDim1>& field1,
                                       TensorField3D<T,nDim2>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct ReductiveBoxProcessingFunctional3D_ST : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, ScalarField3D<T>& field1,
                                       TensorField3D<T,nDim>& field2) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct ReductiveBoxProcessingFunctional3D_LS : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct ReductiveBoxProcessingFunctional3D_LT : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,nDim>& field) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct ReductiveLatticeBoxProcessingFunctional3D : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple scalar fields.
template<typename T>
struct ReductiveScalarFieldBoxProcessingFunctional3D : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct ReductiveTensorFieldBoxProcessingFunctional3D : public PlainReductiveBoxProcessingFunctional3D<T> {
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};



template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_L<T,Descriptor>& functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice);
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_S<T>& functional,
                               Box3D domain, ScalarFieldBase3D<T>& field);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_T<T,nDim>& functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field);


/* *************** All flavors of reductive dot processing functionals ********* */

/// Easy instantiation of reductive boxed data processor (general case)
template<typename T>
class ReductiveDotProcessingFunctional3D {
public:
    virtual ~ReductiveDotProcessingFunctional3D() { }
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual ReductiveDotProcessingFunctional3D<T>* clone() const =0;
    virtual BlockStatistics<T> const& getStatistics() const =0;
    virtual BlockStatistics<T>& getStatistics() =0;
};

/// ReductiveDotProcessingFunctional3D which instantiates its own statistics object
template<typename T>
class PlainReductiveDotProcessingFunctional3D : public ReductiveDotProcessingFunctional3D<T> {
public:
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
};


/// A ReductiveDotted data processor, automatically generated from a ReductiveDotProcessingFunctional3D
template<typename T>
class ReductiveDotProcessor3D : public DataProcessor3D<T> {
public:
    /** \param functional_ The functional is not owned by the ReductiveBoxProcessor3D,
     *                     i.e. it is not deleted in the destructor.
     */
    ReductiveDotProcessor3D(ReductiveDotProcessingFunctional3D<T>* functional_,
                            DotList3D const& dotList_, std::vector<AtomicBlock3D<T>*> atomicBlocks_);
    DotList3D const& getDotList() const;
    virtual void process();
    virtual ReductiveDotProcessor3D<T>* clone() const;
private:
    ReductiveDotProcessingFunctional3D<T>* functional;
    DotList3D dotList;
    std::vector<AtomicBlock3D<T>*> atomicBlocks;
};

/// An automatically created generator for the ReductiveDotProcessor3D
template<typename T>
class ReductiveDotProcessorGenerator3D : public DottedReductiveDataProcessorGenerator3D<T> {
public:
    ReductiveDotProcessorGenerator3D(ReductiveDotProcessingFunctional3D<T>* functional_, DotList3D const& dotList);
    ~ReductiveDotProcessorGenerator3D();
    ReductiveDotProcessorGenerator3D(ReductiveDotProcessorGenerator3D<T> const& rhs);
    ReductiveDotProcessorGenerator3D<T>& operator=(ReductiveDotProcessorGenerator3D<T> const& rhs);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual ReductiveDotProcessorGenerator3D<T>* clone() const;
    virtual BlockStatistics<T> const& getStatistics() const;
    virtual BlockStatistics<T>& getStatistics();
    ReductiveDotProcessingFunctional3D<T> const& getFunctional() const;
private:
    ReductiveDotProcessingFunctional3D<T>* functional;
};

/// Easy instantiation of dotted data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct ReductiveDotProcessingFunctional3D_L : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice) =0;
    virtual ReductiveDotProcessingFunctional3D_L<T,Descriptor>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single scalar field
template<typename T>
struct ReductiveDotProcessingFunctional3D_S : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field) =0;
    virtual ReductiveDotProcessingFunctional3D_S<T>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for a single tensor field
template<typename T, int nDim>
struct ReductiveDotProcessingFunctional3D_T : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, TensorField3D<T,nDim>& field) =0;
    virtual ReductiveDotProcessingFunctional3D_T<T,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct ReductiveDotProcessingFunctional3D_LL : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor1>& lattice1,
                                                   BlockLattice3D<T,Descriptor2>& lattice2) =0;
    virtual ReductiveDotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
struct ReductiveDotProcessingFunctional3D_SS : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field1,
                                                   ScalarField3D<T>& field2) =0;
    virtual ReductiveDotProcessingFunctional3D_SS<T>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct ReductiveDotProcessingFunctional3D_TT : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, TensorField3D<T,nDim1>& field1,
                                                   TensorField3D<T,nDim2>& field2) =0;
    virtual ReductiveDotProcessingFunctional3D_TT<T,nDim1,nDim2>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct ReductiveDotProcessingFunctional3D_ST : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, ScalarField3D<T>& field1,
                                                   TensorField3D<T,nDim>& field2) =0;
    virtual ReductiveDotProcessingFunctional3D_ST<T,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct ReductiveDotProcessingFunctional3D_LS : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice,
                         ScalarField3D<T>& field) =0;
    virtual ReductiveDotProcessingFunctional3D_LS<T,Descriptor>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct ReductiveDotProcessingFunctional3D_LT : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice,
                                                   TensorField3D<T,nDim>& field) =0;
    virtual ReductiveDotProcessingFunctional3D_LT<T,Descriptor,nDim>* clone() const =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct ReductiveLatticeDotProcessingFunctional3D : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotlist, std::vector<BlockLattice3D<T,Descriptor>*> lattices) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotlist, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple scalar fields.
template<typename T>
struct ReductiveScalarFieldDotProcessingFunctional3D : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotlist, std::vector<ScalarField3D<T>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotlist, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};

/// Easy instantiation of dotted data processor for multiple tensor fields.
template<typename T, int nDim>
struct ReductiveTensorFieldDotProcessingFunctional3D : public PlainReductiveDotProcessingFunctional3D<T> {
    virtual void process(DotList3D const& dotlist, std::vector<TensorField3D<T,nDim>*> fields) =0;
    /// Invoke parent-method "processGenericBlocks" through a type-cast
    virtual void processGenericBlocks(DotList3D const& dotlist, std::vector<AtomicBlock3D<T>*> atomicBlocks);
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_L<T,Descriptor>& functional,
                               DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice);
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_S<T>& functional,
                               DotList3D const& dotList, ScalarFieldBase3D<T>& field);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_T<T,nDim>& functional,
                               DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field);

/* *************** All flavors of Bounded Reductive Box processing functionals ********* */

/// Easy instantiation of bounded reductive boxed data processor (general case)
template<typename T>
class BoundedReductiveBoxProcessingFunctional3D {
public:
    virtual ~BoundedReductiveBoxProcessingFunctional3D() { }
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks) =0;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void rescale(T dxScale, T dtScale);
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BoundedReductiveBoxProcessingFunctional3D<T>* clone() const =0;
    ReductiveBoxProcessingFunctional3D<T>* getBulkProcessor() const; 
    ReductiveBoxProcessingFunctional3D<T>* getPlaneProcessor(int direction, int orientation) const; 
    ReductiveBoxProcessingFunctional3D<T>* getEdgeProcessor(int plane, int normal1, int normal2) const; 
    ReductiveBoxProcessingFunctional3D<T>* getCornerProcessor(int normalX, int normalY, int normalZ) const; 
    BlockStatistics<T> const& getStatistics() const;
    BlockStatistics<T>& getStatistics();
private:
    BlockStatistics<T> statistics;
public:
    class BulkWrapperFunctional : public ReductiveBoxProcessingFunctional3D<T> {
    public:
        BulkWrapperFunctional(BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_);
        BulkWrapperFunctional(BulkWrapperFunctional const& rhs);
        ~BulkWrapperFunctional();
        BulkWrapperFunctional& operator=(BulkWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual BulkWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional;
    };
    class PlaneWrapperFunctional : public ReductiveBoxProcessingFunctional3D<T> {
    public:
        PlaneWrapperFunctional(BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
                              int direction_, int orientation_);
        PlaneWrapperFunctional(PlaneWrapperFunctional const& rhs);
        ~PlaneWrapperFunctional();
        PlaneWrapperFunctional& operator=(PlaneWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual PlaneWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional;
        int direction, orientation;
    };
    class EdgeWrapperFunctional : public ReductiveBoxProcessingFunctional3D<T> {
    public:
        EdgeWrapperFunctional(BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
                              int plane, int normal1, int normal2);
        EdgeWrapperFunctional(EdgeWrapperFunctional const& rhs);
        ~EdgeWrapperFunctional();
        EdgeWrapperFunctional& operator=(EdgeWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual EdgeWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional;
        int plane, normal1, normal2;
    };
    class CornerWrapperFunctional : public ReductiveBoxProcessingFunctional3D<T> {
    public:
        CornerWrapperFunctional(BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
                                int normalX_, int normalY_, int normalZ_);
        CornerWrapperFunctional(CornerWrapperFunctional const& rhs);
        ~CornerWrapperFunctional();
        CornerWrapperFunctional& operator=(CornerWrapperFunctional const& rhs);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
        virtual BlockDomain::DomainT appliesTo() const;
        virtual void rescale(T dxScale, T dtScale);
        virtual void getModificationPattern(std::vector<bool>& isWritten) const;
        virtual CornerWrapperFunctional* clone() const;
        virtual BlockStatistics<T> const& getStatistics() const;
        virtual BlockStatistics<T>& getStatistics();
    private:
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional;
        int normalX, normalY, normalZ;
    };
};

/// Generic implementatio1 of "apply" and "integrate" for Bounded Reductive BoxProcessingFunctionals
template<typename T>
class BoundedReductiveOneBlockProcessingFunctionalOperation3D {
public:
    BoundedReductiveOneBlockProcessingFunctionalOperation3D(Box3D const& domain, plint boundaryWidth_);
    void apply(BoundedReductiveBoxProcessingFunctional3D<T>& functional, Block3D<T>& block);
private:
    BlockSurface3D surf;
};

/// Easy instantiation of bounded reductive boxed data processor for a single lattice
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional3D_L : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for a single scalar field
template<typename T>
struct BoundedReductiveBoxProcessingFunctional3D_S : public BoundedReductiveBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, ScalarField3D<T>& field) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              ScalarField3D<T>& field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                ScalarField3D<T>& field ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for a single tensor field
template<typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_T : public BoundedReductiveBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, TensorField3D<T,nDim>& field) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& field ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& field ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& field ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
struct BoundedReductiveBoxProcessingFunctional3D_LL : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for ScalarField-ScalarField coupling
template<typename T>
struct BoundedReductiveBoxProcessingFunctional3D_SS : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
struct BoundedReductiveBoxProcessingFunctional3D_TT : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_ST : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveBoxProcessingFunctional3D_LS : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
struct BoundedReductiveBoxProcessingFunctional3D_LT : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for multiple lattices.
template<typename T, template<typename U> class Descriptor>
struct BoundedReductiveLatticeBoxProcessingFunctional3D : public BoundedReductiveBoxProcessingFunctional3D<T> {
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

/// Easy instantiation of bounded reductive boxed data processor for multiple scalar fields.
template<typename T>
struct BoundedReductiveScalarFieldBoxProcessingFunctional3D : public BoundedReductiveBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, std::vector<ScalarField3D<T>*> fields) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               std::vector<ScalarField3D<T>*> fields ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              std::vector<ScalarField3D<T>*> fields ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                std::vector<ScalarField3D<T>*> fields ) =0;
};

/// Easy instantiation of bounded reductive boxed data processor for multiple tensor fields.
template<typename T, int nDim>
struct BoundedReductiveTensorFieldBoxProcessingFunctional3D : public BoundedReductiveBoxProcessingFunctional3D<T> {
    virtual void processBulkGeneric(Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);
    virtual void processPlaneGeneric( int direction, int orientation, Box3D domain,
                                      std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processEdgeGeneric( int plane, int normal1, int normal2, Box3D domain,
                                     std::vector<AtomicBlock3D<T>*> atomicBlocks );
    virtual void processCornerGeneric( int normalX, int normalY, int normalZ, Box3D domain,
                                       std::vector<AtomicBlock3D<T>*> atomicBlocks );

    virtual void processBulk(Box3D domain, std::vector<TensorField3D<T,nDim>*> fields) =0;
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> fields ) =0;
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              std::vector<TensorField3D<T,nDim>*> fields ) =0;
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                std::vector<TensorField3D<T,nDim>*> fields ) =0;
};


template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>& functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_S<T>& functional,
                               Box3D domain, ScalarFieldBase3D<T>& field, plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>& functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field, plint boundaryWidth);



}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_H
