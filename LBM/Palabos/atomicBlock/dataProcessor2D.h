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
 * Interface for dataProcessing steps -- header file.
 */
#ifndef DATA_PROCESSOR_2D_H
#define DATA_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/blockStatistics.h"
#include <vector>
#include <algorithm>

namespace plb {

// Forward declarations
template<typename T> class AtomicBlock2D;

/// DataProcessors are used to run extended operations on a lattice or data field
template<typename T>
struct DataProcessor2D {
    virtual ~DataProcessor2D() { }
    /// Execute processing operation
    virtual void process() =0;
    /// Clone Data Processor, on its dynamic type
    virtual DataProcessor2D<T>* clone() const =0;
    /// Extent of application area (0 for purely local operations)
    virtual plint extent() const;
    /// Extent of application area along a direction (0 or 1)
    virtual plint extent(int direction) const;
};

/// This is a factory class generating LatticeProcessors
/** The LatticeProcessorGenerator can be tailored (shifted/reduced) to
 *  a sublattice, after which the LatticeProcessor is generated. The
 *  LatticeProcessor itself is static, i.e. the coordinates of the
 *  sublattice to which it refers cannot be changed after construction.
 *  Instead, a new LatticeProcessor must be generated with the generator.
 */
template<typename T>
struct DataProcessorGenerator2D {
    virtual ~DataProcessorGenerator2D();
    /// Shift the domain of application of this data processor.
    virtual void shift(plint deltaX, plint deltaY) =0;
    /// Multiply coordinates of the domain of application of this data processor.
    virtual void multiply(plint scale) =0;
    /// Divide coordinates of the domain of application of this data processor.
    virtual void divide(plint scale) =0;
    /// Extract a subdomain (in-place operation).
    /** \return True if original domain of application and subDomain intersect.
     */
    virtual bool extract(Box2D subDomain) =0;
    /// Generate the data processor.
    virtual DataProcessor2D<T>* generate (
            std::vector<AtomicBlock2D<T>*> atomicBlocks ) const =0;
    /// Clone DataProcessorGenerator, based on its dynamics type.
    virtual DataProcessorGenerator2D<T>* clone() const =0;
    /// Indicates whether data processor should be applied on envelope or not. Defaults to false.
    virtual BlockDomain::DomainT appliesTo() const;
    /// Rescale the physical units of the data processor. Defaults to no rescaling.
    virtual void rescale(T dxScale, T dtScale);
    /// Tell which blocks are modified (written) when the processor is applied on them.
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
};

template<typename T>
class BoxedDataProcessorGenerator2D : public DataProcessorGenerator2D<T> {
public:
    BoxedDataProcessorGenerator2D(Box2D domain_);
    virtual void shift(plint deltaX, plint deltaY);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box2D subDomain);
    Box2D getDomain() const;
private:
    Box2D domain;
};

template<typename T>
class DottedDataProcessorGenerator2D : public DataProcessorGenerator2D<T> {
public:
    DottedDataProcessorGenerator2D(DotList2D const& dots_);
    virtual void shift(plint deltaX, plint deltaY);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box2D subDomain);
    DotList2D const& getDotList() const;
private:
    DotList2D dots;
};

template<typename T>
class ReductiveDataProcessorGenerator2D {
public:
    virtual ~ReductiveDataProcessorGenerator2D();
    /// Shift the domain of application of this data processor.
    virtual void shift(plint deltaX, plint deltaY) =0;
    /// Multiply coordinates of the domain of application of this data processor.
    virtual void multiply(plint scale) =0;
    /// Divide coordinates of the domain of application of this data processor.
    virtual void divide(plint scale) =0;
    /// Extract a subdomain (in-place operation).
    /** \return True if original domain of application and subDomain intersect.
     */
    virtual bool extract(Box2D subDomain) =0;
    /// Generate the data processor.
    virtual DataProcessor2D<T>* generate (
            std::vector<AtomicBlock2D<T>*> atomicBlocks ) =0;
    /// Clone ReductiveDataProcessorGenerator, based on its dynamics type.
    virtual ReductiveDataProcessorGenerator2D<T>* clone() const =0;
    /// Const handle to statistics object.
    virtual BlockStatistics<T> const& getStatistics() const =0;
    /// Non-const handle to statistics object.
    virtual BlockStatistics<T>& getStatistics() =0;
    /// Indicates whether data processor should be applied on envelope or not. Defaults to false.
    virtual BlockDomain::DomainT appliesTo() const;
    /// Rescale the physical units of the data processor.
    virtual void rescale(T dxScale, T dtScale);
    /// Tell which blocks are modified (written) when the processor is applied on them.
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
};

template<typename T>
class BoxedReductiveDataProcessorGenerator2D : public ReductiveDataProcessorGenerator2D<T> {
public:
    BoxedReductiveDataProcessorGenerator2D(Box2D domain_);
    virtual void shift(plint deltaX, plint deltaY);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box2D subDomain);
    Box2D getDomain() const;
private:
    Box2D domain;
};

template<typename T>
class DottedReductiveDataProcessorGenerator2D : public ReductiveDataProcessorGenerator2D<T> {
public:
    DottedReductiveDataProcessorGenerator2D(DotList2D const& dots_);
    virtual void shift(plint deltaX, plint deltaY);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box2D subDomain);
    DotList2D const& getDotList() const;
private:
    DotList2D dots;
};

}  // namespace plb

#endif  // DATA_PROCESSOR_2D_H
