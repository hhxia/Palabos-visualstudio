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
 * Serial version of the 3D block communicator -- generic implementation.
 */
#ifndef SERIAL_BLOCK_COMMUNICATOR_3D_HH
#define SERIAL_BLOCK_COMMUNICATOR_3D_HH

#include "multiBlock/serialBlockCommunicator3D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "multiBlock/multiBlock3D.h"
#include "core/plbDebug.h"

namespace plb {

////////////////////// Class SerialBlockCommunicator3D /////////////////////

template<typename T>
SerialBlockCommunicator3D<T>::SerialBlockCommunicator3D()
{ }

template<typename T>
SerialBlockCommunicator3D<T>* SerialBlockCommunicator3D<T>::clone() const {
    return new SerialBlockCommunicator3D<T>;
}

template<typename T>
void SerialBlockCommunicator3D<T>::broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const {
    // Nothing to do in the serial case
}

template<typename T>
void SerialBlockCommunicator3D<T>::broadCastVector(T* vect, plint length, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const {
    // Nothing to do in the serial case
}

template<typename T>
void SerialBlockCommunicator3D<T>::copyOverlap(Overlap3D const& overlap, MultiBlock3D<T>& multiBlock) const
{
    MultiBlockManagement3D const& multiBlockManagement = multiBlock.getMultiBlockManagement();
    MultiBlockDistribution3D const& multiBlockDistribution = multiBlockManagement.getMultiBlockDistribution();
    plint originalId = overlap.getOriginalId();
    plint overlapId  = overlap.getOverlapId();
    BlockParameters3D const& originalParameters = multiBlockDistribution.getBlockParameters(originalId);
    BlockParameters3D const& overlapParameters  = multiBlockDistribution.getBlockParameters(overlapId);

    Box3D originalCoords(originalParameters.toLocal(overlap.getOriginalCoordinates()));
    Box3D overlapCoords(overlapParameters.toLocal(overlap.getOverlapCoordinates()));

    PLB_PRECONDITION(originalCoords.x1-originalCoords.x0 == overlapCoords.x1-overlapCoords.x0);
    PLB_PRECONDITION(originalCoords.y1-originalCoords.y0 == overlapCoords.y1-overlapCoords.y0);
    PLB_PRECONDITION(originalCoords.z1-originalCoords.z0 == overlapCoords.z1-overlapCoords.z0);

    AtomicBlock3D<T>* originalBlock = &multiBlock.getComponent(originalId);
    AtomicBlock3D<T>* overlapBlock = &multiBlock.getComponent(overlapId);
    plint deltaX = originalCoords.x0 - overlapCoords.x0;
    plint deltaY = originalCoords.y0 - overlapCoords.y0;
    plint deltaZ = originalCoords.z0 - overlapCoords.z0;

    overlapBlock -> getDataTransfer().attribute(overlapCoords, deltaX, deltaY, deltaZ, *originalBlock);
}

template<typename T>
void SerialBlockCommunicator3D<T>::duplicateOverlaps(MultiBlock3D<T>& multiBlock) const
{
    MultiBlockManagement3D const& multiBlockManagement = multiBlock.getMultiBlockManagement();

    // Non-periodic communication
    MultiBlockDistribution3D const& multiBlockDistribution = multiBlockManagement.getMultiBlockDistribution();
    for (plint iOverlap=0; iOverlap<multiBlockDistribution.getNumNormalOverlaps(); ++iOverlap) {
        copyOverlap(multiBlockDistribution.getNormalOverlap(iOverlap), multiBlock);
    }

    // Periodic communication
    PeriodicitySwitch3D<T> const& periodicity = multiBlock.periodicity();
    for (plint iOverlap=0; iOverlap<multiBlockDistribution.getNumPeriodicOverlaps(); ++iOverlap) {
        PeriodicOverlap3D const& pOverlap = multiBlockDistribution.getPeriodicOverlap(iOverlap);
        if (periodicity.get(pOverlap.normalX, pOverlap.normalY, pOverlap.normalZ)) {
            copyOverlap(pOverlap.overlap, multiBlock);
        }
    }
}

template<typename T>
void SerialBlockCommunicator3D<T>::signalPeriodicity() const
{ }


}  // namespace plb

#endif
