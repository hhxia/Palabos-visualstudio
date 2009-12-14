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

#ifndef SEND_RECV_POOL_H
#define SEND_RECV_POOL_H

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include <map>
#include <vector>

namespace plb {

#ifdef PLB_MPI_PARALLEL

struct PoolEntry {
    PoolEntry()
        : cumDataLength(0), lengths()
    { }
    int cumDataLength;
    std::vector<int> lengths;
};

class SendRecvPool {
public:
    typedef std::map<int, PoolEntry> SubsT;
public:
    void subscribeMessage(int proc, int numData) {
        PoolEntry& entry = subscriptions[proc];
        entry.lengths.push_back(numData);
        entry.cumDataLength += numData;
    }
    void clear() {
        subscriptions.clear();
    }
    SubsT::const_iterator begin() const { return subscriptions.begin(); }
    SubsT::const_iterator end() const { return subscriptions.end(); }
private:
    SubsT subscriptions;
};

template <typename T>
struct CommunicatorEntry {
    CommunicatorEntry() { } 
    CommunicatorEntry(PoolEntry const& poolEntry)
        : lengths(poolEntry.lengths),
          data(poolEntry.cumDataLength),
          currentMessage(0),
          currentPos(0)
    { }
    void reset() {
        currentMessage=0;
        currentPos=0;
    }
    std::vector<int> lengths;
    std::vector<T> data;
    int currentMessage;
    int currentPos;
    MPI_Request request;
    MPI_Status  status;
};

template <typename T>
class SendPoolCommunicator {
public:
    SendPoolCommunicator() { }
    SendPoolCommunicator(SendRecvPool const& pool);
    T* getSendBuffer(int toProc);
    void acceptMessage(int toProc);
    void finalize();
private:
    void startCommunication(int toProc);
private:
    std::map<int, CommunicatorEntry<T> > subscriptions;
};

template <typename T>
class RecvPoolCommunicator {
public:
    RecvPoolCommunicator() { }
    RecvPoolCommunicator(SendRecvPool const& pool);
    void startBeingReceptive();
    T const* receiveMessage(int fromProc);
private:
    void finalize(int fromProc);
private:
    std::map<int, CommunicatorEntry<T> > subscriptions;
};

#endif  // PLB_MPI_PARALLEL

}  // namespace plb

#endif  // SEND_RECV_POOL_H
