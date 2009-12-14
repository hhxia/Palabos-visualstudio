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
 * Scalar, vector and tensor fields for 2D data fields -- generic implementation.
 */

#ifndef MULTI_DATA_FIELD_2D_HH
#define MULTI_DATA_FIELD_2D_HH

#include <algorithm>
#include <limits>
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"

namespace plb {

/////// Class MultiScalarField2D //////////////////////////////////

template<typename T>
MultiScalarField2D<T>::MultiScalarField2D (
        MultiBlockManagement2D const& multiBlockManagement_,
        BlockCommunicator2D<T>* blockCommunicator_,
        CombinedStatistics<T>* combinedStatistics_,
        MultiScalarAccess2D<T>* multiScalarAccess_ )
    : MultiBlock2D<T>(multiBlockManagement_, blockCommunicator_, combinedStatistics_ ),
      multiScalarAccess(multiScalarAccess_)
{
    allocateFields();
}

template<typename T>
MultiScalarField2D<T>::MultiScalarField2D(plint nx, plint ny)
    : MultiBlock2D<T>(nx,ny),
      multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField2D<T>::~MultiScalarField2D() {
    deAllocateFields();
    delete multiScalarAccess;
}

template<typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiScalarField2D<T> const& rhs)
    : ScalarFieldBase2D<T>(rhs),
      MultiBlock2D<T>(rhs),
      multiScalarAccess(rhs.multiScalarAccess->clone())
{
    allocateFields();
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        *(fields[iBlock]) = *(rhs.fields[iBlock]);
    }
}

template<typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiBlock2D<T> const& rhs)
    : MultiBlock2D<T>(rhs),
      multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField2D<T>::MultiScalarField2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop)
    : MultiBlock2D<T>(rhs, subDomain, crop),
      multiScalarAccess(defaultMultiBlockPolicy2D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField2D<T>& MultiScalarField2D<T>::operator=(MultiScalarField2D<T> const& rhs) {
    MultiScalarField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
void MultiScalarField2D<T>::swap(MultiScalarField2D<T>& rhs) {
    MultiBlock2D<T>::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiScalarAccess, rhs.multiScalarAccess);
}

template<typename T>
void MultiScalarField2D<T>::reset() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        fields[iBlock] -> reset();
    }
}

template<typename T>
void MultiScalarField2D<T>::allocateFields() 
{
    for (plint iBlock=0; iBlock < getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            Box2D const& envelope = this->getMultiBlockManagement().getEnvelope(iBlock);
            ScalarField2D<T>* newField =
                new ScalarField2D<T> (
                        envelope.getNx(), envelope.getNy() );
            newField -> setLocation(Dot2D(envelope.x0, envelope.y0));
            fields.push_back(newField);
        }
        else {
            fields.push_back( 0 );
        }
    }
}

template<typename T>
void MultiScalarField2D<T>::deAllocateFields() 
{
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        delete fields[iBlock];
    }
}

template<typename T>
T& MultiScalarField2D<T>::get(plint iX, plint iY) {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    return multiScalarAccess->getDistributedScalar(iX,iY, this->getMultiBlockManagement(), fields);
}

template<typename T>
T const& MultiScalarField2D<T>::get(plint iX, plint iY) const {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    return multiScalarAccess->getDistributedScalar(iX,iY, this->getMultiBlockManagement(), fields);
}

template<typename T>
AtomicBlock2D<T>& MultiScalarField2D<T>::getComponent(plint iBlock) {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T>
AtomicBlock2D<T> const& MultiScalarField2D<T>::getComponent(plint iBlock) const {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T>
plint MultiScalarField2D<T>::sizeOfCell() const {
    return 1;
}

template<typename T>
identifiers::BlockId MultiScalarField2D<T>::getBlockId() const {
    return identifiers::getScalarId<T>();
}

template<typename T>
BlockParameters2D const& MultiScalarField2D<T>::getParameters(plint iParam) const {
    return getMultiBlockDistribution().getBlockParameters(iParam);
}

template<typename T>
MultiBlockDistribution2D const& MultiScalarField2D<T>::getMultiBlockDistribution() const {
    return this->getMultiBlockManagement().getMultiBlockDistribution();
}




//////// Class MultiTensorField2D //////////////////////////////////

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::MultiTensorField2D (
        MultiBlockManagement2D const& multiBlockManagement_,
        BlockCommunicator2D<T>* blockCommunicator_,
        CombinedStatistics<T>* combinedStatistics_,
        MultiTensorAccess2D<T,nDim>* multiTensorAccess_ )
    : MultiBlock2D<T>(multiBlockManagement_, blockCommunicator_, combinedStatistics_ ),
      multiTensorAccess(multiTensorAccess_)
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::MultiTensorField2D(plint nx, plint ny)
    : MultiBlock2D<T>(nx,ny),
      multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::~MultiTensorField2D() {
    deAllocateFields();
    delete multiTensorAccess;
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::MultiTensorField2D(MultiTensorField2D<T,nDim> const& rhs)
    : TensorFieldBase2D<T,nDim>(rhs),
      MultiBlock2D<T>(rhs),
      multiTensorAccess(rhs.multiTensorAccess->clone())
{
    allocateFields();
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        *(fields[iBlock]) = *(rhs.fields[iBlock]);
    }
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::MultiTensorField2D(MultiBlock2D<T> const& rhs)
    : MultiBlock2D<T>(rhs),
      multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>::MultiTensorField2D(MultiBlock2D<T> const& rhs, Box2D subDomain, bool crop)
    : MultiBlock2D<T>(rhs, subDomain, crop),
      multiTensorAccess(defaultMultiBlockPolicy2D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>& MultiTensorField2D<T,nDim>::operator=(MultiTensorField2D<T,nDim> const& rhs) {
    MultiTensorField2D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, int nDim>
void MultiTensorField2D<T,nDim>::swap(MultiTensorField2D<T,nDim>& rhs) {
    MultiBlock2D<T>::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiTensorAccess, rhs.multiTensorAccess);
}

template<typename T, int nDim>
void MultiTensorField2D<T,nDim>::reset() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        fields[iBlock] -> reset();
    }
}

template<typename T, int nDim>
void MultiTensorField2D<T,nDim>::allocateFields() 
{
    for (plint iBlock=0; iBlock < getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            Box2D const& envelope = this->getMultiBlockManagement().getEnvelope(iBlock);
            TensorField2D<T,nDim>* newField =
                new TensorField2D<T,nDim> (
                        envelope.getNx(), envelope.getNy() );
            newField -> setLocation(Dot2D(envelope.x0, envelope.y0));
            fields.push_back(newField);
        }
        else {
            fields.push_back( 0 );
        }
    }
}


template<typename T, int nDim>
void MultiTensorField2D<T,nDim>::deAllocateFields() 
{
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        delete fields[iBlock];
    }
}

template<typename T, int nDim>
inline Array<T,nDim>& 
MultiTensorField2D<T,nDim>::get(plint iX, plint iY) {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    return multiTensorAccess->getDistributedTensor(iX,iY, this->getMultiBlockManagement(), fields);
}

template<typename T, int nDim>
inline Array<T,nDim> const& 
MultiTensorField2D<T,nDim>::get(plint iX, plint iY) const {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    return multiTensorAccess->getDistributedTensor(iX,iY, this->getMultiBlockManagement(), fields);
}

template<typename T, int nDim>
AtomicBlock2D<T>& MultiTensorField2D<T,nDim>::getComponent(plint iBlock) {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T, int nDim>
AtomicBlock2D<T> const& MultiTensorField2D<T,nDim>::getComponent(plint iBlock) const {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T, int nDim>
plint MultiTensorField2D<T,nDim>::sizeOfCell() const {
    return nDim;
}

template<typename T, int nDim>
identifiers::BlockId MultiTensorField2D<T,nDim>::getBlockId() const {
    return identifiers::getTensorId<T,nDim>();
}

template<typename T, int nDim>
BlockParameters2D const& MultiTensorField2D<T,nDim>::getParameters(plint iParam) const {
    return getMultiBlockDistribution().getBlockParameters(iParam);
}

template<typename T, int nDim>
MultiBlockDistribution2D const& MultiTensorField2D<T,nDim>::getMultiBlockDistribution() const {
    return this->getMultiBlockManagement().getMultiBlockDistribution();
}

}  // namespace plb

#endif  
