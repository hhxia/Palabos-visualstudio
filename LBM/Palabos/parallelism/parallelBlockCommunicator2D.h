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
 * Helper classes for parallel 2D multiblock lattice -- header file.
 */

#ifndef PARALLEL_BLOCK_COMMUNICATOR_2D_H
#define PARALLEL_BLOCK_COMMUNICATOR_2D_H

#include "core/globalDefs.h"
#include "core/periodicity2D.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "parallelism/sendRecvPool.h"
#include "parallelism/communicationPackage2D.h"
#include <vector>

namespace plb {

#ifdef PLB_MPI_PARALLEL

template<typename T> class DataTransmittor2D;

template<typename T>
class ParallelBlockCommunicator2D : public BlockCommunicator2D<T> {
public:
    ParallelBlockCommunicator2D();
    ParallelBlockCommunicator2D(ParallelBlockCommunicator2D<T> const& rhs);
    virtual ParallelBlockCommunicator2D<T>* clone() const;
    virtual void broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const;
    virtual void broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const;
    virtual void duplicateOverlaps(MultiBlock2D<T>& multiBlock) const;
    virtual void signalPeriodicity() const;
private:
    ParallelBlockCommunicator2D<T>& operator= (
            ParallelBlockCommunicator2D<T> const& rhs ) { return *this; }
    void setupCommunicationStructure (
            MultiBlockManagement2D const& multiBlockManagement, plint sizeOfCell,
            PeriodicitySwitch2D<T> const& periodicity ) const;
    void subscribeOverlap (
        Overlap2D const& overlap, MultiBlockManagement2D const& multiBlockManagement,
        SendRecvPool& sendPool, SendRecvPool& recvPool, plint sizeOfCell ) const;
private:
    mutable bool needsUpdate;
    mutable CommunicationPackage2D sendPackage;
    mutable CommunicationPackage2D recvPackage;
    mutable CommunicationPackage2D sendRecvPackage;
    mutable SendPoolCommunicator<T> sendComm;
    mutable RecvPoolCommunicator<T> recvComm;
};

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif
