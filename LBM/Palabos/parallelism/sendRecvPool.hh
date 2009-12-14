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
 * Wrapper functions that simplify the use of MPI
 */

#ifndef SEND_RECV_POOL_HH
#define SEND_RECV_POOL_HH

#include "parallelism/mpiManager.h"
#include "parallelism/sendRecvPool.h"
#include "core/plbDebug.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

template <typename T>
SendPoolCommunicator<T>::SendPoolCommunicator(SendRecvPool const& pool)
    : subscriptions(pool.begin(), pool.end())
{ }

template <typename T>
T* SendPoolCommunicator<T>::getSendBuffer(int toProc) {
    CommunicatorEntry<T>& entry = subscriptions[toProc];
    PLB_PRECONDITION( entry.currentMessage < (int)entry.lengths.size() );
    return &entry.data[entry.currentPos];
}

template <typename T>
void SendPoolCommunicator<T>::acceptMessage(int toProc) {
    CommunicatorEntry<T>& entry = subscriptions[toProc];
    PLB_PRECONDITION( entry.currentMessage < (int)entry.lengths.size() );
    int length = entry.lengths[entry.currentMessage];
    entry.currentMessage++;
    entry.currentPos += length;

    if (entry.currentMessage==(int)entry.lengths.size()) {
        startCommunication(toProc);
        entry.reset();
    }
}

template <typename T>
void SendPoolCommunicator<T>::finalize() {
    typename std::map<int, CommunicatorEntry<T> >::iterator iter = subscriptions.begin();
    for (; iter != subscriptions.end(); ++iter) {
        CommunicatorEntry<T>& entry = iter->second;
        global::mpi().wait(&entry.request, &entry.status);
    }
}

template <typename T>
void SendPoolCommunicator<T>::startCommunication(int toProc)
{
    CommunicatorEntry<T>& entry = subscriptions[toProc];
    global::mpi().iSend(&entry.data[0], entry.data.size(), toProc, &entry.request);
}

template <typename T>
RecvPoolCommunicator<T>::RecvPoolCommunicator(SendRecvPool const& pool)
    : subscriptions(pool.begin(), pool.end())
{ }

template <typename T>
void RecvPoolCommunicator<T>::startBeingReceptive() {
    typename std::map<int, CommunicatorEntry<T> >::iterator iter = subscriptions.begin();
    for (; iter != subscriptions.end(); ++iter) {
        int fromProc = iter->first;
        CommunicatorEntry<T>& entry = iter->second;
        global::mpi().iRecv(&entry.data[0], entry.data.size(), fromProc, &entry.request);
    }
}

template <typename T>
T const* RecvPoolCommunicator<T>::receiveMessage(int fromProc) {
    CommunicatorEntry<T>& entry = subscriptions[fromProc];
    PLB_PRECONDITION( entry.currentMessage < (int)entry.lengths.size() );
    if (entry.currentMessage==0) {
        finalize(fromProc);
    }
    T const* message = &entry.data[entry.currentPos];
    int length = entry.lengths[entry.currentMessage];
    entry.currentMessage++;
    entry.currentPos += length;
    if (entry.currentMessage==(int)entry.lengths.size()) {
        entry.reset();
    }
    return message;
}

template <typename T>
void RecvPoolCommunicator<T>::finalize(int fromProc) {
    CommunicatorEntry<T>& entry = subscriptions[fromProc];
    global::mpi().wait(&entry.request, &entry.status);
}

#endif // PLB_MPI_PARALLEL

}  // namespace plb


#endif  // SEND_RECV_POOL_HH
