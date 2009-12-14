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

#ifndef RUN_TIME_DIAGNOSTICS_H
#define RUN_TIME_DIAGNOSTICS_H

#include "core/globalDefs.h"
#include <string>
#include <fstream>
#include <exception>

namespace plb {

class DiagnosticFileSingleton {
public:
    DiagnosticFileSingleton(std::string fName);
    void write(std::string message);
private:
    DiagnosticFileSingleton();
    std::ofstream ofile;
};

DiagnosticFileSingleton& warningFile();
void plbWarning(bool issueWarning, std::string message);
void plbMainProcWarning(bool issueWarning, std::string message);

DiagnosticFileSingleton& errorFile();
void plbMemoryError(bool issueError, std::string message);
void plbIOError(bool issueError, std::string message);
void plbNetworkError(bool issueError, std::string message);
void plbLogicError(bool issueError, std::string message);

class XflException : public std::exception
{ };

class XflMemoryException : public XflException {
public:
    XflMemoryException(std::string message_) throw();
    virtual ~XflMemoryException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class XflIOException : public XflException {
public:
    XflIOException(std::string message_) throw();
    virtual ~XflIOException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class XflNetworkException : public XflException {
public:
    XflNetworkException(std::string message_) throw();
    virtual ~XflNetworkException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

class XflLogicException : public XflException {
public:
    XflLogicException(std::string message_) throw();
    virtual ~XflLogicException() throw() { }
    virtual const char* what() const throw();
private:
    std::string message;
};

}  // namespace plb

#endif  // RUN_TIME_DIAGNOSTICS_H
