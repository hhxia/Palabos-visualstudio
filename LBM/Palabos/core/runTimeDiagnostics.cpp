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

#include "core/runTimeDiagnostics.h"
#include "parallelism/mpiManager.h"

namespace plb {

DiagnosticFileSingleton::DiagnosticFileSingleton(std::string fName)
    : ofile ( ( std::string(global::directories().getLogOutDir() + fName) ).c_str() )
{ }

void DiagnosticFileSingleton::write(std::string message) {
    ofile << message << std::endl;
}

DiagnosticFileSingleton& warningFile() {
    static DiagnosticFileSingleton warningFileSingleton("plbWarning.dat");
    return warningFileSingleton;
}

void plbWarning(bool issueWarning, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int warningCounter = issueWarning ? 1 : 0;
    int warningSum = 0;
    global::mpi().reduce(warningCounter, warningSum, MPI_SUM);
    issueWarning = warningSum > 0;

    if (global::mpi().isMainProcessor()) {
#endif
        if (issueWarning) {
            warningFile().write(message);
        }
#ifdef PLB_MPI_PARALLEL
    }
#endif
}

void plbMainProcWarning(bool issueWarning, std::string message) {
    if (issueWarning) {
        warningFile().write(message);
    }
}


DiagnosticFileSingleton& errorFile() {
    static DiagnosticFileSingleton errorFileSingleton("plbError.dat");
    return errorFileSingleton;
}

void plbMemoryError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw XflMemoryException(message);
    }
}

void plbIOError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw XflIOException(message);
    }
}

void plbNetworkError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw XflNetworkException(message);
    }
}

void plbLogicError(bool issueError, std::string message) {
#ifdef PLB_MPI_PARALLEL
    int errorCounter = issueError ? 1 : 0;
    global::mpi().reduceAndBcast(errorCounter,  MPI_SUM);
    issueError = errorCounter > 0;
#endif

    if (issueError) {
        errorFile().write(message);
        throw XflLogicException(message);
    }
}


XflMemoryException::XflMemoryException(std::string message_) throw()
    : message(std::string("Palabos Memory error: ") + message_)
{ }

const char* XflMemoryException::what() const throw() {
    return message.c_str();
}


XflIOException::XflIOException(std::string message_) throw()
    : message(std::string("Palabos I/O error: ") + message_)
{ }

const char* XflIOException::what() const throw() {
    return message.c_str();
}


XflNetworkException::XflNetworkException(std::string message_) throw()
    : message(std::string("Palabos Network error: ") + message_)
{ }

const char* XflNetworkException::what() const throw() {

    return message.c_str();
}

XflLogicException::XflLogicException(std::string message_) throw()
    : message(std::string("Palabos logic error: ") + message_)
{ }

const char* XflLogicException::what() const throw() {
    return message.c_str();
}

}  // namespace plb
