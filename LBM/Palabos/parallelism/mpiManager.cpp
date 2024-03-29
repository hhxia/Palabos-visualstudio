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
 * Wrapper functions that simplify the use of MPI, template instatiations
 */

#ifdef PLB_MPI_PARALLEL

#include "parallelism/mpiManager.h"
#include <iostream>

namespace plb {

namespace global {

MpiManager::MpiManager() : ok(false)
{ }

MpiManager::~MpiManager() {
    if (ok) {
        MPI_Finalize();
        ok = false;
    }
}

void MpiManager::init(int *argc, char ***argv, bool verbous) {
    if (verbous) {
        std::cerr << "Constructing an MPI thread" << std::endl;
    }
    int ok1 = MPI_Init(argc, argv);
    int ok2 = MPI_Comm_rank(MPI_COMM_WORLD,&taskId);
    int ok3 = MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
    ok = (ok1==0 && ok2==0 && ok3==0);
}

int MpiManager::getSize() const {
    return numTasks;
}

int MpiManager::getRank() const {
    return taskId;
}

int MpiManager::bossId() const {
    return 0;
}

bool MpiManager::isMainProcessor() const {
    return bossId() == getRank();
}

double MpiManager::getTime() const {
    if (!ok) return 0.;
    return MPI_Wtime();
}

void MpiManager::barrier(MPI_Comm comm) {
    if (!ok) return;
    MPI_Barrier(comm);
}

template <>
void MpiManager::send<char>(char *buf, int count, int dest, int tag, MPI_Comm comm) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm);
}

template <>
void MpiManager::send<int>(int *buf, int count, int dest, int tag, MPI_Comm comm) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm);
}

template <>
void MpiManager::send<long>(long *buf, int count, int dest, int tag, MPI_Comm comm) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_LONG, dest, tag, comm);
}

template <>
void MpiManager::send<float>(float *buf, int count, int dest, int tag, MPI_Comm comm) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm);
}

template <>
void MpiManager::send<double>(double *buf, int count, int dest, int tag, MPI_Comm comm) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm);
}

template <>
void MpiManager::iSend<char>
    (char *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, request);
    }
}

template <>
void MpiManager::iSend<int>
    (int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, request);
    }
}

template <>
void MpiManager::iSend<float>
    (float *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, request);
    }
}

template <>
void MpiManager::iSend<double>
    (double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, request);
    }
}

template <>
void MpiManager::iSendRequestFree<char>
    (char *buf, int count, int dest, int tag, MPI_Comm comm)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<int>
    (int *buf, int count, int dest, int tag, MPI_Comm comm)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<float>
    (float *buf, int count, int dest, int tag, MPI_Comm comm)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<double>
    (double *buf, int count, int dest, int tag, MPI_Comm comm)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::receive<char>(char *buf, int count, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::receive<int>(int *buf, int count, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<long>(long *buf, int count, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_LONG, source, tag, comm, &status);
}

template <>
void MpiManager::receive<float>(float *buf, int count, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<double>(double *buf, int count, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::sendToMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<long>(long* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}


template <>
void MpiManager::iRecv<char>(char *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, request);
    }
}

template <>
void MpiManager::iRecv<int>(int *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, request);
    }
}

template <>
void MpiManager::iRecv<long>(long *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_LONG, source, tag, comm, request);
    }
}

template <>
void MpiManager::iRecv<float>(float *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, request);
    }
}

template <>
void MpiManager::iRecv<double>(double *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, request);
    }
}

template <>
void MpiManager::sendRecv<char>
    (char *sendBuf, char *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<int>
    (int *sendBuf, int *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_INT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<float>
    (float *sendBuf, float *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_FLOAT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<double>
    (double *sendBuf, double *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_DOUBLE, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::scatterv_impl<char>(char* sendBuf, int* sendCounts, int* displs,
                                     char* recvBuf, int recvCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_CHAR,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::scatterv_impl<int>(int *sendBuf, int* sendCounts, int* displs,
                                int* recvBuf, int recvCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_INT,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_INT, root, comm);
}

template <>
void MpiManager::scatterv_impl<long>(long *sendBuf, int* sendCounts, int* displs,
                                     long* recvBuf, int recvCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_LONG,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_LONG, root, comm);
}

template <>
void MpiManager::scatterv_impl<float>(float *sendBuf, int* sendCounts, int* displs,
                                      float* recvBuf, int recvCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_FLOAT,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::scatterv_impl<double>(double *sendBuf, int* sendCounts, int* displs,
                                       double* recvBuf, int recvCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Scatterv(static_cast<void*>(sendBuf),
                 sendCounts, displs, MPI_DOUBLE,
                 static_cast<void*>(recvBuf),
                 recvCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::gatherv_impl<char>(char* sendBuf, int sendCount,
                                    char* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_CHAR,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_CHAR,
                root, comm);
}

template <>
void MpiManager::gatherv_impl<int>(int* sendBuf, int sendCount,
                                   int* recvBuf, int* recvCounts, int* displs,
                                   int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_INT,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_INT,
                root, comm);
}

template <>
void MpiManager::gatherv_impl<long>(long* sendBuf, int sendCount,
                                    long* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_LONG,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_LONG,
                root, comm);
}

template <>
void MpiManager::gatherv_impl<float>(float* sendBuf, int sendCount,
                                 float* recvBuf, int* recvCounts, int* displs,
                                 int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_FLOAT,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_FLOAT,
                root, comm);
}

template <>
void MpiManager::gatherv_impl<double>(double* sendBuf, int sendCount,
                                  double* recvBuf, int* recvCounts, int* displs,
                                  int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_DOUBLE,
                static_cast<void*>(recvBuf), recvCounts, displs, MPI_DOUBLE,
                root, comm);
}

template <>
void MpiManager::bCast<char>(char* sendBuf, int sendCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::bCast<int>(int* sendBuf, int sendCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_INT, root, comm);
}

template <>
void MpiManager::bCast<long>(long* sendBuf, int sendCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_LONG, root, comm);
}

template <>
void MpiManager::bCast<float>(float* sendBuf, int sendCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::bCast<double>(double* sendBuf, int sendCount, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Bcast(static_cast<void*>(sendBuf),
              sendCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::bCastThroughMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<long>(long* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
    bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::reduce<char>(char sendVal, char& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduce<int>(int sendVal, int& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduce<long>(long sendVal, long& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_LONG, op, root, comm);
}

template <>
void MpiManager::reduce<float>(float sendVal, float& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduce<double>(double sendVal, double& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&sendVal),
               static_cast<void*>(&recvVal), 1, MPI_DOUBLE, op, root, comm);
}

template <>
void MpiManager::reduceVect<char>(std::vector<char>& sendVal, std::vector<char>& recvVal,
                                 MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&(sendVal[0])),
               static_cast<void*>(&(recvVal[0])),
               sendVal.size(), MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduceVect<int>(std::vector<int>& sendVal, std::vector<int>& recvVal,
                                 MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&(sendVal[0])),
               static_cast<void*>(&(recvVal[0])),
               sendVal.size(), MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduceVect<long>(std::vector<long>& sendVal, std::vector<long>& recvVal,
                                  MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&(sendVal[0])),
               static_cast<void*>(&(recvVal[0])),
               sendVal.size(), MPI_LONG, op, root, comm);
}

template <>
void MpiManager::reduceVect<float>(std::vector<float>& sendVal, std::vector<float>& recvVal,
                               MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&(sendVal[0])),
               static_cast<void*>(&(recvVal[0])),
               sendVal.size(), MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduceVect<double>(std::vector<double>& sendVal, std::vector<double>& recvVal,
                                MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    MPI_Reduce(static_cast<void*>(&(sendVal[0])),
               static_cast<void*>(&(recvVal[0])),
               sendVal.size(), MPI_DOUBLE, op, root, comm);
}

template <>
void MpiManager::reduceAndBcast<char>(char& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    char recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_CHAR, op, root, comm);
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_CHAR, root, comm);

}

template <>
void MpiManager::reduceAndBcast<int>(int& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    int recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_INT, op, root, comm);
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_INT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<long>(long& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    int recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG, op, root, comm);
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_LONG, root, comm);

}

template <>
void MpiManager::reduceAndBcast<float>(float& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    float recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_FLOAT, op, root, comm);
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_FLOAT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<double>(double& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
    if (!ok) return;
    double recvVal;
    MPI_Reduce(&reductVal, &recvVal, 1, MPI_DOUBLE, op, root, comm);
    reductVal = recvVal;
    MPI_Bcast(&reductVal, 1, MPI_DOUBLE, root, comm);

}

void MpiManager::wait(MPI_Request* request, MPI_Status* status)
{
    if (!ok) return;
    MPI_Wait(request, status);
}

}  // namespace global


}

#endif
