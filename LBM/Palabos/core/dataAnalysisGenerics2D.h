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
 * Generic algorithms for data analysis.
 */
#ifndef DATA_ANALYSIS_GENERICS_2D_H
#define DATA_ANALYSIS_GENERICS_2D_H

namespace plb {

/* ******** CountLatticeElementsFunctional2D **************************** */

template<typename T, template<typename U> class Descriptor, class BoolMask> 
CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>::CountLatticeElementsFunctional2D (
        BoolMask boolMask_)
  : countId(this->getStatistics().subscribeIntSum()),
    boolMask(boolMask_)
{ }

template<typename T, template<typename U> class Descriptor, class BoolMask> 
void CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            if (boolMask(cell)) {
                statistics.gatherIntSum(countId, 1);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>*
    CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>::clone() const
{
    return new CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>(*this);
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}


/* ******** CountScalarElementsFunctional2D **************************** */

template<typename T, class BoolMask> 
CountScalarElementsFunctional2D<T,BoolMask>::CountScalarElementsFunctional2D (
        BoolMask boolMask_)
  : countId(this->getStatistics().subscribeIntSum()),
    boolMask(boolMask_)
{ }

template<typename T, class BoolMask> 
void CountScalarElementsFunctional2D<T,BoolMask>::process (
        Box2D domain, ScalarField2D<T>& field )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T value = field.get(iX,iY);
            if (boolMask(value)) {
                statistics.gatherIntSum(countId, 1);
            }
        }
    }
}

template<typename T, class BoolMask> 
CountScalarElementsFunctional2D<T,BoolMask>*
    CountScalarElementsFunctional2D<T,BoolMask>::clone() const
{
    return new CountScalarElementsFunctional2D<T,BoolMask>(*this);
}

template<typename T, class BoolMask> 
plint CountScalarElementsFunctional2D<T,BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}


/* ******** CountTensorElementsFunctional2D **************************** */

template<typename T, int nDim, class BoolMask> 
CountTensorElementsFunctional2D<T,nDim,BoolMask>::CountTensorElementsFunctional2D (
        BoolMask boolMask_)
  : countId(this->getStatistics().subscribeIntSum()),
    boolMask(boolMask_)
{ }

template<typename T, int nDim, class BoolMask> 
void CountTensorElementsFunctional2D<T,nDim,BoolMask>::process (
        Box2D domain, TensorField2D<T,nDim>& field )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,nDim> const& value = field.get(iX,iY);
            if (boolMask(value)) {
                statistics.gatherIntSum(countId, 1);
            }
        }
    }
}

template<typename T, int nDim, class BoolMask> 
CountTensorElementsFunctional2D<T,nDim,BoolMask>*
    CountTensorElementsFunctional2D<T,nDim,BoolMask>::clone() const
{
    return new CountTensorElementsFunctional2D<T,nDim,BoolMask>(*this);
}

template<typename T, int nDim, class BoolMask> 
plint CountTensorElementsFunctional2D<T,nDim,BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}


/* ******** ApplyScalarFunctional2D ************************************* */

template<typename T, class Function>
ApplyScalarFunctional2D<T,Function>::ApplyScalarFunctional2D(Function f_)
    : f(f_)
{ }

template<typename T,class Function>
void ApplyScalarFunctional2D<T,Function>::process (
        Box2D domain, ScalarField2D<T>& field )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            field.get(iX,iY) = f(field.get(iX,iY));
        }
    }
}

template<typename T,class Function>
ApplyScalarFunctional2D<T,Function>* ApplyScalarFunctional2D<T,Function>::clone() const {
    return new ApplyScalarFunctional2D<T,Function>(*this);
}

template<typename T,class Function>
void ApplyScalarFunctional2D<T,Function>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T,class Function>
BlockDomain::DomainT ApplyScalarFunctional2D<T,Function>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** EvaluateScalarFunctional2D ************************************* */

template<typename T, class Function>
EvaluateScalarFunctional2D<T,Function>::EvaluateScalarFunctional2D(Function f_)
    : f(f_)
{ }

template<typename T, class Function>
void EvaluateScalarFunctional2D<T,Function>::process (
        Box2D domain, ScalarField2D<T>& field, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(field, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = f(field.get(iX,iY));
        }
    }
}

template<typename T, class Function>
EvaluateScalarFunctional2D<T,Function>* EvaluateScalarFunctional2D<T,Function>::clone() const {
    return new EvaluateScalarFunctional2D<T,Function>(*this);
}

template<typename T, class Function>
void EvaluateScalarFunctional2D<T,Function>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, class Function>
BlockDomain::DomainT EvaluateScalarFunctional2D<T,Function>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_GENERICS_2D_H
