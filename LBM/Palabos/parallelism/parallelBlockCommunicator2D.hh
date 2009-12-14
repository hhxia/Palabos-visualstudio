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
 * Helper classes for parallel 2D multiblock lattice -- generic implementation.
 */
#ifndef PARALLEL_BLOCK_COMMUNICATOR_2D_HH
#define PARALLEL_BLOCK_COMMUNICATOR_2D_HH

#include "parallelism/parallelBlockCommunicator2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include <cmath>

namespace plb {

#ifdef PLB_MPI_PARALLEL

////////////////////// Class ParallelBlockCommunicator2D /////////////////////

template<typename T>
ParallelBlockCommunicator2D<T>::ParallelBlockCommunicator2D()
    : needsUpdate(true)
{ }

template<typename T>
ParallelBlockCommunicator2D<T>::ParallelBlockCommunicator2D (
        ParallelBlockCommunicator2D<T> const& rhs )
  : needsUpdate(true)
{ }

template<typename T>
ParallelBlockCommunicator2D<T>* ParallelBlockCommunicator2D<T>::clone() const {
    return new ParallelBlockCommunicator2D<T>;
}

template<typename T>
void ParallelBlockCommunicator2D<T>::broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const {
    plint fromProc = multiBlockManagement.getMultiBlockDistribution().getBlockParameters(fromBlock).getProcId();
    global::mpi().bCast(&scalar, 1, fromProc);
}

template<typename T>
void ParallelBlockCommunicator2D<T>::broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const {
    plint fromProc = multiBlockManagement.getMultiBlockDistribution().getBlockParameters(fromBlock).getProcId();
    global::mpi().bCast(data, length, fromProc);
}

template<typename T>
void ParallelBlockCommunicator2D<T>::setupCommunicationStructure (
        MultiBlockManagement2D const& multiBlockManagement, plint sizeOfCell,
        PeriodicitySwitch2D<T> const& periodicity ) const
{
    if (!needsUpdate) {
        return;
    }
    needsUpdate = false;

    sendPackage.clear();
    recvPackage.clear();
    sendRecvPackage.clear();

    SendRecvPool sendPool, recvPool;

    RelevantIndexes2D const& relevantIndexes = multiBlockManagement.getRelevantIndexes();
    MultiBlockDistribution2D const& dataDistribution = multiBlockManagement.getMultiBlockDistribution();
    for (pluint iRelevant=0; iRelevant<relevantIndexes.getNormalOverlaps().size(); ++iRelevant) {
        plint iOverlap = relevantIndexes.getNormalOverlaps()[iRelevant];
        subscribeOverlap (
            dataDistribution.getNormalOverlap(iOverlap), multiBlockManagement,
            sendPool, recvPool, sizeOfCell );
    }
    for (pluint iRelevant=0; iRelevant<relevantIndexes.getPeriodicOverlaps().size(); ++iRelevant) {
        plint iOverlap = relevantIndexes.getPeriodicOverlaps()[iRelevant];
        PeriodicOverlap2D const& pOverlap = dataDistribution.getPeriodicOverlap(iOverlap);
        if (periodicity.get(pOverlap.normalX,pOverlap.normalY)) {
            subscribeOverlap (
                pOverlap.overlap, multiBlockManagement,
                sendPool, recvPool, sizeOfCell );
        }
    }
    sendComm = SendPoolCommunicator<T>(sendPool);
    recvComm = RecvPoolCommunicator<T>(recvPool);
}

template<typename T>
void ParallelBlockCommunicator2D<T>::subscribeOverlap (
        Overlap2D const& overlap, MultiBlockManagement2D const& multiBlockManagement,
        SendRecvPool& sendPool, SendRecvPool& recvPool, plint sizeOfCell ) const
{
    CommunicationInfo2D info;

    info.fromBlockId = overlap.getOriginalId();
    info.toBlockId   = overlap.getOverlapId();

    MultiBlockDistribution2D const& dataDistribution = multiBlockManagement.getMultiBlockDistribution();
    BlockParameters2D const& originalParameters = dataDistribution.getBlockParameters(info.fromBlockId);
    BlockParameters2D const& overlapParameters  = dataDistribution.getBlockParameters(info.toBlockId);

    info.fromProcessId = originalParameters.getProcId();
    info.toProcessId   = overlapParameters.getProcId();

    info.fromDomain = originalParameters.toLocal(overlap.getOriginalCoordinates());
    info.toDomain   = overlapParameters.toLocal(overlap.getOverlapCoordinates());

    plint lx = info.fromDomain.x1-info.fromDomain.x0+1;
    plint ly = info.fromDomain.y1-info.fromDomain.y0+1;
    PLB_PRECONDITION(lx == info.toDomain.x1-info.toDomain.x0+1);
    PLB_PRECONDITION(ly == info.toDomain.y1-info.toDomain.y0+1);

    plint numberOfCells = lx*ly;

    ThreadAttribution const& attribution = multiBlockManagement.getThreadAttribution();
    if (attribution.isLocal(info.fromProcessId) && attribution.isLocal(info.toProcessId)) {
        sendRecvPackage.push_back(info);
    }
    else if (attribution.isLocal(info.fromProcessId)) {
        sendPackage.push_back(info);
        sendPool.subscribeMessage(info.toProcessId, numberOfCells*sizeOfCell);
    }
    else if (attribution.isLocal(info.toProcessId)) {
        recvPackage.push_back(info);
        recvPool.subscribeMessage(info.fromProcessId, numberOfCells*sizeOfCell);
    }
}

template<typename T>
void ParallelBlockCommunicator2D<T>::duplicateOverlaps(MultiBlock2D<T>& multiBlock) const
{
    MultiBlockManagement2D const& multiBlockManagement = multiBlock.getMultiBlockManagement();
    PeriodicitySwitch2D<T> const& periodicity          = multiBlock.periodicity();

    setupCommunicationStructure(multiBlockManagement, multiBlock.sizeOfCell(), periodicity);

    // 1. Non-blocking receives.
    recvComm.startBeingReceptive();

    // 2. Non-blocking sends.
    for (unsigned iSend=0; iSend<sendPackage.size(); ++iSend) {
        CommunicationInfo2D const& info = sendPackage[iSend];
        AtomicBlock2D<T>& fromBlock = multiBlock.getComponent(info.fromBlockId);
        fromBlock.getDataTransfer().send (
                info.fromDomain, sendComm.getSendBuffer(info.toProcessId) );
        sendComm.acceptMessage(info.toProcessId);
    }

    // 3. Local copies which require no communication.
    for (unsigned iSendRecv=0; iSendRecv<sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo2D const& info = sendRecvPackage[iSendRecv];
        AtomicBlock2D<T>& fromBlock = multiBlock.getComponent(info.fromBlockId);
        AtomicBlock2D<T>& toBlock = multiBlock.getComponent(info.toBlockId);
        plint deltaX = info.fromDomain.x0 - info.toDomain.x0;
        plint deltaY = info.fromDomain.y0 - info.toDomain.y0;
        toBlock.getDataTransfer().attribute (
                info.toDomain, deltaX, deltaY, fromBlock );
    }

    // 4. Finalize the receives.
    for (unsigned iRecv=0; iRecv<recvPackage.size(); ++iRecv) {
        CommunicationInfo2D const& info = recvPackage[iRecv];
        AtomicBlock2D<T>& toBlock = multiBlock.getComponent(info.toBlockId);
        toBlock.getDataTransfer().receive (
                info.toDomain, recvComm.receiveMessage(info.fromProcessId) );
    }

    // 5. Finalize the sends.
    sendComm.finalize();
}

template<typename T>
void ParallelBlockCommunicator2D<T>::signalPeriodicity() const {
    needsUpdate = true;
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif
