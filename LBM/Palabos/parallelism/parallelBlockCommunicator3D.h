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
 * Helper classes for parallel 3D multiblock lattice -- header file.
 */

#ifndef PARALLEL_BLOCK_COMMUNICATOR_3D_H
#define PARALLEL_BLOCK_COMMUNICATOR_3D_H

#include "core/globalDefs.h"
#include "core/periodicity3D.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "parallelism/sendRecvPool.h"
#include "parallelism/communicationPackage3D.h"
#include <vector>

namespace plb {

#ifdef PLB_MPI_PARALLEL

template<typename T> class DataTransmittor3D;

template<typename T>
class ParallelBlockCommunicator3D : public BlockCommunicator3D<T> {
public:
    ParallelBlockCommunicator3D();
    ParallelBlockCommunicator3D(ParallelBlockCommunicator3D<T> const& rhs);
    virtual ParallelBlockCommunicator3D<T>* clone() const;
    virtual void broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const;
    virtual void broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const;
    virtual void duplicateOverlaps(MultiBlock3D<T>& multiBlock) const;
    virtual void signalPeriodicity() const;
private:
    ParallelBlockCommunicator3D<T>& operator= (
            ParallelBlockCommunicator3D<T> const& rhs ) { return *this; }
    void setupCommunicationStructure (
            MultiBlockManagement3D const& multiBlockManagement, plint sizeOfCell,
            PeriodicitySwitch3D<T> const& periodicity ) const;
    void subscribeOverlap (
        Overlap3D const& overlap, MultiBlockManagement3D const& multiBlockManagement,
        SendRecvPool& sendPool, SendRecvPool& recvPool, plint sizeOfCell ) const;
private:
    mutable bool needsUpdate;
    mutable CommunicationPackage3D sendPackage;
    mutable CommunicationPackage3D recvPackage;
    mutable CommunicationPackage3D sendRecvPackage;
    mutable SendPoolCommunicator<T> sendComm;
    mutable RecvPoolCommunicator<T> recvComm;
};

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif
