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

#ifndef DATA_PROCESSOR_2D_HH
#define DATA_PROCESSOR_2D_HH

#include "atomicBlock/dataProcessor2D.h"
#include "core/util.h"

namespace plb {

////////////////////// Class DataProcessor2D /////////////////

/** This method returns the maximum extent of the processor, over any
 *  direction. By extent one means the the size of the neighborhood
 *  for non-local accesses. For example, if the processor implements
 *  a second-order accurate laplace operator, the value returned
 *  by extent() is 1, because the second derivative is evaluated
 *  with help of the -1 2 -1 stencil, requiring one left and one
 *  right neighbor.
 *  The default implementation of this method returns an extent of 1,
 *  based on the assumption that most LB stuff is somehow based on
 *  nearest-neighbor interaction.
 *  This is a bit dangerous though, as one easily forgets to override
 *  the method in case of larger-than-nearest-neighbor Processors.
 *  Still, LatticeProcessors are much easier to write when there's only
 *  one method to override, so we'll let it be this way.
 */
template<typename T>
plint DataProcessor2D<T>::extent() const {
    return 1;
}

/** By default, this method assumes a symmetric neighborhood relation
 *  and refers to the non-directed version of extent().
 */
template<typename T>
plint DataProcessor2D<T>::extent(int direction) const {
    return extent();
}

////////////////////// Class DataProcessorGenerator2D /////////////////

template<typename T>
DataProcessorGenerator2D<T>::~DataProcessorGenerator2D()
{ }

template<typename T>
BlockDomain::DomainT DataProcessorGenerator2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/** Default implementation does nothing, to reflect the fact that many
 *  processor generators are invariant under change of scale. It is
 *  however important to remember to redefine this method for data
 *  processors which dependent on physical units, because otherwise
 *  the processor yields erroneous results in a MultiGrid environment.
 *
 *  \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 */
template<typename T>
void DataProcessorGenerator2D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void DataProcessorGenerator2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}


////////////////////// Class BoxedDataProcessorGenerator2D /////////////////

template<typename T>
BoxedDataProcessorGenerator2D<T>::BoxedDataProcessorGenerator2D(Box2D domain_)
    : domain(domain_)
{
}

template<typename T>
void BoxedDataProcessorGenerator2D<T>::shift(plint deltaX, plint deltaY) {
    domain = domain.shift(deltaX, deltaY);
}

template<typename T>
void BoxedDataProcessorGenerator2D<T>::multiply(plint scale) {
    domain = domain.multiply(scale);
}

template<typename T>
void BoxedDataProcessorGenerator2D<T>::divide(plint scale) {
    domain = domain.divide(scale);
}

template<typename T>
bool BoxedDataProcessorGenerator2D<T>::extract(Box2D subDomain) {
    Box2D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
Box2D BoxedDataProcessorGenerator2D<T>::getDomain() const {
    return domain;
}


////////////////////// Class DottedDataProcessorGenerator2D /////////////////

template<typename T>
DottedDataProcessorGenerator2D<T>::DottedDataProcessorGenerator2D (
        DotList2D const& dots_)
    : dots(dots_)
{ }

template<typename T>
void DottedDataProcessorGenerator2D<T>::shift(plint deltaX, plint deltaY) {
    dots = dots.shift(deltaX,deltaY);
}

template<typename T>
void DottedDataProcessorGenerator2D<T>::multiply(plint scale) {
    dots = dots.multiply(scale);
}

template<typename T>
void DottedDataProcessorGenerator2D<T>::divide(plint scale) {
    dots = dots.divide(scale);
}

template<typename T>
bool DottedDataProcessorGenerator2D<T>::extract(Box2D subDomain) {
    DotList2D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
DotList2D const& DottedDataProcessorGenerator2D<T>::getDotList() const {
    return dots;
}


////////////////////// Class ReductiveDataProcessorGenerator2D /////////////////

template<typename T>
ReductiveDataProcessorGenerator2D<T>::~ReductiveDataProcessorGenerator2D()
{ }

template<typename T>
BlockDomain::DomainT ReductiveDataProcessorGenerator2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/** Default implementation does nothing, to reflect the fact that many
 *  processor generators are invariant under change of scale. It is
 *  however important to remember to redefine this method for data
 *  processors which dependent on physical units, because otherwise
 *  the processor yields erroneous results in a MultiGrid environment.
 *
 *  \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 */
template<typename T>
void ReductiveDataProcessorGenerator2D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveDataProcessorGenerator2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}


////////////////////// Class BoxedReductiveDataProcessorGenerator2D /////////////////

template<typename T>
BoxedReductiveDataProcessorGenerator2D<T>::BoxedReductiveDataProcessorGenerator2D(Box2D domain_)
    : domain(domain_)
{ }

template<typename T>
void BoxedReductiveDataProcessorGenerator2D<T>::shift(plint deltaX, plint deltaY) {
    domain = domain.shift(deltaX, deltaY);
}

template<typename T>
void BoxedReductiveDataProcessorGenerator2D<T>::multiply(plint scale) {
    domain = domain.multiply(scale);
}

template<typename T>
void BoxedReductiveDataProcessorGenerator2D<T>::divide(plint scale) {
    domain = domain.divide(scale);
}

template<typename T>
bool BoxedReductiveDataProcessorGenerator2D<T>::extract(Box2D subDomain) {
    Box2D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
Box2D BoxedReductiveDataProcessorGenerator2D<T>::getDomain() const {
    return domain;
}


////////////////////// Class DottedReductiveDataProcessorGenerator2D /////////////////

template<typename T>
DottedReductiveDataProcessorGenerator2D<T>::DottedReductiveDataProcessorGenerator2D (
        DotList2D const& dots_)
    : dots(dots_)
{ }

template<typename T>
void DottedReductiveDataProcessorGenerator2D<T>::shift(plint deltaX, plint deltaY) {
    dots = dots.shift(deltaX,deltaY);
}

template<typename T>
void DottedReductiveDataProcessorGenerator2D<T>::multiply(plint scale) {
    dots = dots.multiply(scale);
}

template<typename T>
void DottedReductiveDataProcessorGenerator2D<T>::divide(plint scale) {
    dots = dots.divide(scale);
}

template<typename T>
bool DottedReductiveDataProcessorGenerator2D<T>::extract(Box2D subDomain) {
    DotList2D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
DotList2D const& DottedReductiveDataProcessorGenerator2D<T>::getDotList() const {
    return dots;
}


}  // namespace plb

#endif  // DATA_PROCESSOR_2D_HH
