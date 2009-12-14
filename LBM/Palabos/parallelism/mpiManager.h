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

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#include "core/globalDefs.h"

#ifdef PLB_MPI_PARALLEL
#include "mpi.h"
#include <vector>
#endif


namespace plb {

namespace global {

#ifdef PLB_MPI_PARALLEL

/// Wrapper functions that simplify the use of MPI
class MpiManager {
public:
    /// Initializes the mpi manager
    void init(int *argc, char ***argv, bool verbous=false);
    /// Returns the number of processes
    int getSize() const;
    /// Returns the process ID
    int getRank() const;
    /// Returns process ID of main processor
    int bossId() const;
    /// Tells whether current processor is main processor
    bool isMainProcessor() const;
    /// Returns universal MPI-time in seconds
    double getTime() const;

    /// Synchronizes the processes
    void barrier(MPI_Comm comm = MPI_COMM_WORLD);

    /// Sends data at *buf, blocking
    template <typename T>
    void send(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Sends data at *buf, non blocking
    template <typename T>
    void iSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Sends data at *buf, non blocking and request free
    template <typename T>
    void iSendRequestFree(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Receives data at *buf, blocking
    template <typename T>
    void receive(T *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Receives data at *buf, non blocking
    template <typename T>
    void iRecv(T *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Send and receive data between two partners
    template <typename T>
    void sendRecv(T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = 0,
                  MPI_Comm comm = MPI_COMM_WORLD);

    /// Sends data to master processor
    template <typename T>
    void sendToMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

    /// Scatter data from one processor over multiple processors
    template <typename T>
    void scatterV(T *sendBuf, T *recvBuf, int* sendCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Gather data from multiple processors to one processor
    template <typename T>
    void gatherV(T* sendBuf, T* recvBuf, int *recvCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);


    /// Broadcast data from one processor to multiple processors
    template <typename T>
    void bCast(T* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Broadcast data when root is unknown to other processors
    template <typename T>
    void bCastThroughMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

    /// Reduction operation toward one processor
    template <typename T>
    void reduce(T sendVal, T& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);

    /// Element-per-element reduction of a vector of data
    template <typename T>
    void reduceVect(std::vector<T>& sendVal, std::vector<T>& recvVal,
                    MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Reduction operation, followed by a broadcast
    template <typename T>
    void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

    /// Complete a non-blocking MPI operation
    void wait(MPI_Request* request, MPI_Status* status);

private:
    /// Implementation code for Scatter
    template <typename T>
    void scatterv_impl(T *sendBuf, int* sendCounts, int* displs,
                       T* recvBuf, int recvCount, int root, MPI_Comm comm);

    /// Implementation code for Gather
    template <typename T>
    void gatherv_impl(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs,
                      int root, MPI_Comm comm);
private:
    MpiManager();
    ~MpiManager();
private:
    int numTasks, taskId;
    bool ok;

friend MpiManager& mpi();
};

#else  // #ifdef PLB_MPI_PARALLEL

class MpiManager {
public:
    /// Initializes the mpi manager
    void init(int *argc, char ***argv, bool verbous=false) { }
    /// Returns the number of processes
    int getSize() const { return 1; }
    /// Returns the process ID
    int getRank() const { return 0; }
    /// Returns process ID of main processor
    int bossId() const { return 0; }
    /// Tells whether current processor is main processor
    bool isMainProcessor() const { return true; }
    /// Broadcast data from one processor to multiple processors
    template <typename T>
    void bCast(T* sendBuf, int sendCount, int root = 0) { }

friend MpiManager& mpi();
};

#endif  // PLB_MPI_PARALLEL

inline MpiManager& mpi() {
    static MpiManager instance;
    return instance;
}

}  // namespace global

}  // namespace plb


#endif  // MPI_MANAGER_H
