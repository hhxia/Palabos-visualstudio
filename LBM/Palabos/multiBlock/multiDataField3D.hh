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
 * Scalar, vector and tensor fields for 3D data fields -- generic implementation.
 */

#ifndef MULTI_DATA_FIELD_3D_HH
#define MULTI_DATA_FIELD_3D_HH

#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include <algorithm>
#include <limits>

namespace plb {

/////// Class MultiScalarField3D //////////////////////////////////

template<typename T>
MultiScalarField3D<T>::MultiScalarField3D (
        MultiBlockManagement3D const& multiBlockManagement_,
        BlockCommunicator3D<T>* blockCommunicator_,
        CombinedStatistics<T>* combinedStatistics_,
        MultiScalarAccess3D<T>* multiScalarAccess_ )
    : MultiBlock3D<T>(multiBlockManagement_, blockCommunicator_, combinedStatistics_ ),
      multiScalarAccess(multiScalarAccess_)
{
    allocateFields();
}

template<typename T>
MultiScalarField3D<T>::MultiScalarField3D(plint nx, plint ny, plint nz)
    : MultiBlock3D<T>(nx,ny,nz),
      multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField3D<T>::~MultiScalarField3D() {
    deAllocateFields();
    delete multiScalarAccess;
}

template<typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiScalarField3D<T> const& rhs)
    : ScalarFieldBase3D<T>(rhs),
      MultiBlock3D<T>(rhs),
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
MultiScalarField3D<T>::MultiScalarField3D(MultiBlock3D<T> const& rhs)
    : MultiBlock3D<T>(rhs),
      multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField3D<T>::MultiScalarField3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop)
    : MultiBlock3D<T>(rhs, subDomain, crop),
      multiScalarAccess(defaultMultiBlockPolicy3D().getMultiScalarAccess<T>())
{
    allocateFields();
}

template<typename T>
MultiScalarField3D<T>& MultiScalarField3D<T>::operator=(MultiScalarField3D<T> const& rhs) {
    MultiScalarField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
void MultiScalarField3D<T>::swap(MultiScalarField3D<T>& rhs) {
    MultiBlock3D<T>::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiScalarAccess, rhs.multiScalarAccess);
}

template<typename T>
void MultiScalarField3D<T>::reset() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        fields[iBlock] -> reset();
    }
}

template<typename T>
void MultiScalarField3D<T>::allocateFields() 
{
    for (plint iBlock=0; iBlock < getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            Box3D const& envelope = this->getMultiBlockManagement().getEnvelope(iBlock);
            ScalarField3D<T>* newField =
                new ScalarField3D<T> (
                        envelope.getNx(), envelope.getNy(), envelope.getNz() );
            newField -> setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
            fields.push_back(newField);
        }
        else {
            fields.push_back( 0 );
        }
    }
}

template<typename T>
void MultiScalarField3D<T>::deAllocateFields() 
{
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        delete fields[iBlock];
    }
}

template<typename T>
inline T& MultiScalarField3D<T>::get(plint iX, plint iY, plint iZ) {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
    return multiScalarAccess->getDistributedScalar(iX,iY,iZ, this->getMultiBlockManagement(), fields);
}

template<typename T>
inline T const& MultiScalarField3D<T>::get(plint iX, plint iY, plint iZ) const {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
    return multiScalarAccess->getDistributedScalar(iX,iY,iZ, this->getMultiBlockManagement(), fields);
}

template<typename T>
AtomicBlock3D<T>& MultiScalarField3D<T>::getComponent(plint iBlock) {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T>
AtomicBlock3D<T> const& MultiScalarField3D<T>::getComponent(plint iBlock) const {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T>
plint MultiScalarField3D<T>::sizeOfCell() const {
    return 1;
}

template<typename T>
identifiers::BlockId MultiScalarField3D<T>::getBlockId() const {
    return identifiers::getScalarId<T>();
}

template<typename T>
BlockParameters3D const& MultiScalarField3D<T>::getParameters(plint iParam) const {
    return getMultiBlockDistribution().getBlockParameters(iParam);
}

template<typename T>
MultiBlockDistribution3D const& MultiScalarField3D<T>::getMultiBlockDistribution() const {
    return this->getMultiBlockManagement().getMultiBlockDistribution();
}



//////// Class MultiTensorField3D //////////////////////////////////

template<typename T, int nDim>
MultiTensorField3D<T,nDim>::MultiTensorField3D (
        MultiBlockManagement3D const& multiBlockManagement_,
        BlockCommunicator3D<T>* blockCommunicator_,
        CombinedStatistics<T>* combinedStatistics_,
        MultiTensorAccess3D<T,nDim>* multiTensorAccess_ )
    : MultiBlock3D<T>(multiBlockManagement_, blockCommunicator_, combinedStatistics_ ),
      multiTensorAccess(multiTensorAccess_)
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField3D<T,nDim>::MultiTensorField3D(plint nx, plint ny, plint nz)
    : MultiBlock3D<T>(nx,ny,nz),
      multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField3D<T,nDim>::~MultiTensorField3D() {
    deAllocateFields();
    delete multiTensorAccess;
}

template<typename T, int nDim>
MultiTensorField3D<T,nDim>::MultiTensorField3D(MultiTensorField3D<T,nDim> const& rhs)
    : TensorFieldBase3D<T,nDim>(rhs),
      MultiBlock3D<T>(rhs),
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
MultiTensorField3D<T,nDim>::MultiTensorField3D(MultiBlock3D<T> const& rhs)
    : MultiBlock3D<T>(rhs),
      multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField3D<T,nDim>::MultiTensorField3D(MultiBlock3D<T> const& rhs, Box3D subDomain, bool crop)
    : MultiBlock3D<T>(rhs, subDomain, crop),
      multiTensorAccess(defaultMultiBlockPolicy3D().getMultiTensorAccess<T,nDim>())
{
    allocateFields();
}

template<typename T, int nDim>
MultiTensorField3D<T,nDim>& MultiTensorField3D<T,nDim>::operator=(MultiTensorField3D<T,nDim> const& rhs) {
    MultiTensorField3D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, int nDim>
void MultiTensorField3D<T,nDim>::swap(MultiTensorField3D<T,nDim>& rhs) {
    MultiBlock3D<T>::swap(rhs);
    fields.swap(rhs.fields);
    std::swap(multiTensorAccess, rhs.multiTensorAccess);
}

template<typename T, int nDim>
void MultiTensorField3D<T,nDim>::reset() {
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        fields[iBlock] -> reset();
    }
}

template<typename T, int nDim>
void MultiTensorField3D<T,nDim>::allocateFields() 
{
    for (plint iBlock=0; iBlock < getMultiBlockDistribution().getNumBlocks(); ++iBlock) {
        if ( this->getMultiBlockManagement().getThreadAttribution().isLocal (
                  getParameters(iBlock).getProcId()) )
        {
            Box3D const& envelope = this->getMultiBlockManagement().getEnvelope(iBlock);
            TensorField3D<T,nDim>* newField =
                new TensorField3D<T,nDim> (
                        envelope.getNx(), envelope.getNy(), envelope.getNz() );
            newField -> setLocation(Dot3D(envelope.x0, envelope.y0, envelope.z0));
            fields.push_back(newField);
        }
        else {
            fields.push_back( 0 );
        }
    }
}

template<typename T, int nDim>
void MultiTensorField3D<T,nDim>::deAllocateFields() 
{
    std::vector<plint> const& relevantBlocks = this->getRelevantBlocks();
    for (plint rBlock=0; rBlock < this->getNumRelevantBlocks(); ++rBlock) {
        plint iBlock = relevantBlocks[rBlock];
        delete fields[iBlock];
    }
}

template<typename T, int nDim>
inline Array<T,nDim>& 
MultiTensorField3D<T,nDim>::get(plint iX, plint iY, plint iZ) {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
    return multiTensorAccess->getDistributedTensor(iX,iY,iZ, this->getMultiBlockManagement(), fields);
}

template<typename T, int nDim>
inline Array<T,nDim> const& 
MultiTensorField3D<T,nDim>::get(plint iX, plint iY, plint iZ) const {
    PLB_PRECONDITION(iX>=0 && iX<this->getNx());
    PLB_PRECONDITION(iY>=0 && iY<this->getNy());
    PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
    return multiTensorAccess->getDistributedTensor(iX,iY,iZ, this->getMultiBlockManagement(), fields);
}

template<typename T, int nDim>
AtomicBlock3D<T>& MultiTensorField3D<T,nDim>::getComponent(plint iBlock) {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T, int nDim>
AtomicBlock3D<T> const& MultiTensorField3D<T,nDim>::getComponent(plint iBlock) const {
    PLB_PRECONDITION( iBlock<(plint)fields.size() );
    return *fields[iBlock];
}

template<typename T, int nDim>
plint MultiTensorField3D<T,nDim>::sizeOfCell() const {
    return nDim;
}

template<typename T, int nDim>
identifiers::BlockId MultiTensorField3D<T,nDim>::getBlockId() const {
    return identifiers::getTensorId<T,nDim>();
}

template<typename T, int nDim>
BlockParameters3D const& MultiTensorField3D<T,nDim>::getParameters(plint iParam) const {
    return getMultiBlockDistribution().getBlockParameters(iParam);
}

template<typename T, int nDim>
MultiBlockDistribution3D const& MultiTensorField3D<T,nDim>::getMultiBlockDistribution() const {
    return this->getMultiBlockManagement().getMultiBlockDistribution();
}

}  // namespace plb

#endif  
