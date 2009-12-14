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

#ifndef DATA_PROCESSOR_3D_HH
#define DATA_PROCESSOR_3D_HH

#include "atomicBlock/dataProcessor3D.h"
#include "core/util.h"

namespace plb {

////////////////////// Class DataProcessor3D /////////////////

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
plint DataProcessor3D<T>::extent() const {
    return 1;
}

/** By default, this method assumes a symmetric neighborhood relation
 *  and refers to the non-directed version of extent().
 */
template<typename T>
plint DataProcessor3D<T>::extent(int direction) const {
    return extent();
}

////////////////////// Class DataProcessorGenerator3D /////////////////

template<typename T>
DataProcessorGenerator3D<T>::~DataProcessorGenerator3D()
{ }

template<typename T>
BlockDomain::DomainT DataProcessorGenerator3D<T>::appliesTo() const
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
void DataProcessorGenerator3D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void DataProcessorGenerator3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}


////////////////////// Class BoxedDataProcessorGenerator3D /////////////////

template<typename T>
BoxedDataProcessorGenerator3D<T>::BoxedDataProcessorGenerator3D(Box3D domain_)
    : domain(domain_)
{
}

template<typename T>
void BoxedDataProcessorGenerator3D<T>::shift(plint deltaX, plint deltaY, plint deltaZ) {
    domain = domain.shift(deltaX, deltaY, deltaZ);
}

template<typename T>
void BoxedDataProcessorGenerator3D<T>::multiply(plint scale) {
    domain = domain.multiply(scale);
}

template<typename T>
void BoxedDataProcessorGenerator3D<T>::divide(plint scale) {
    domain = domain.divide(scale);
}

template<typename T>
bool BoxedDataProcessorGenerator3D<T>::extract(Box3D subDomain) {
    Box3D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
Box3D BoxedDataProcessorGenerator3D<T>::getDomain() const {
    return domain;
}


////////////////////// Class DottedDataProcessorGenerator3D /////////////////

template<typename T>
DottedDataProcessorGenerator3D<T>::DottedDataProcessorGenerator3D (
        DotList3D const& dots_)
    : dots(dots_)
{ }

template<typename T>
void DottedDataProcessorGenerator3D<T>::shift(plint deltaX, plint deltaY, plint deltaZ) {
    dots = dots.shift(deltaX,deltaY,deltaZ);
}

template<typename T>
void DottedDataProcessorGenerator3D<T>::multiply(plint scale) {
    dots = dots.multiply(scale);
}

template<typename T>
void DottedDataProcessorGenerator3D<T>::divide(plint scale) {
    dots = dots.divide(scale);
}

template<typename T>
bool DottedDataProcessorGenerator3D<T>::extract(Box3D subDomain) {
    DotList3D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
DotList3D const& DottedDataProcessorGenerator3D<T>::getDotList() const {
    return dots;
}


////////////////////// Class ReductiveDataProcessorGenerator3D /////////////////

template<typename T>
ReductiveDataProcessorGenerator3D<T>::~ReductiveDataProcessorGenerator3D()
{ }

template<typename T>
BlockDomain::DomainT ReductiveDataProcessorGenerator3D<T>::appliesTo() const
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
void ReductiveDataProcessorGenerator3D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveDataProcessorGenerator3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}


////////////////////// Class BoxedReductiveDataProcessorGenerator3D /////////////////

template<typename T>
BoxedReductiveDataProcessorGenerator3D<T>::BoxedReductiveDataProcessorGenerator3D(Box3D domain_)
    : domain(domain_)
{ }

template<typename T>
void BoxedReductiveDataProcessorGenerator3D<T>::shift(plint deltaX, plint deltaY, plint deltaZ) {
    domain = domain.shift(deltaX, deltaY, deltaZ);
}

template<typename T>
void BoxedReductiveDataProcessorGenerator3D<T>::multiply(plint scale) {
    domain = domain.multiply(scale);
}

template<typename T>
void BoxedReductiveDataProcessorGenerator3D<T>::divide(plint scale) {
    domain = domain.divide(scale);
}

template<typename T>
bool BoxedReductiveDataProcessorGenerator3D<T>::extract(Box3D subDomain) {
    Box3D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
Box3D BoxedReductiveDataProcessorGenerator3D<T>::getDomain() const {
    return domain;
}


////////////////////// Class DottedReductiveDataProcessorGenerator3D /////////////////

template<typename T>
DottedReductiveDataProcessorGenerator3D<T>::DottedReductiveDataProcessorGenerator3D (
        DotList3D const& dots_)
    : dots(dots_)
{ }

template<typename T>
void DottedReductiveDataProcessorGenerator3D<T>::shift(plint deltaX, plint deltaY, plint deltaZ) {
    dots = dots.shift(deltaX,deltaY,deltaZ);
}

template<typename T>
void DottedReductiveDataProcessorGenerator3D<T>::multiply(plint scale) {
    dots = dots.multiply(scale);
}

template<typename T>
void DottedReductiveDataProcessorGenerator3D<T>::divide(plint scale) {
    dots = dots.divide(scale);
}

template<typename T>
bool DottedReductiveDataProcessorGenerator3D<T>::extract(Box3D subDomain) {
    DotList3D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
DotList3D const& DottedReductiveDataProcessorGenerator3D<T>::getDotList() const {
    return dots;
}


}  // namespace plb

#endif  // DATA_PROCESSOR_3D_HH
