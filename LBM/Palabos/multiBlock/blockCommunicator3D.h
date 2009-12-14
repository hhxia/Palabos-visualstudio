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
 * Block Communicator -- Abstract base class.
 */
#ifndef BLOCK_COMMUNICATOR_3D_H
#define BLOCK_COMMUNICATOR_3D_H

#include "core/globalDefs.h"

namespace plb {

template<typename T> class MultiBlock3D;
class MultiBlockManagement3D;

template<typename T>
struct BlockCommunicator3D {
    virtual ~BlockCommunicator3D() { }
    virtual BlockCommunicator3D<T>* clone() const =0;
    virtual void broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const =0;
    virtual void broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const =0;
    virtual void duplicateOverlaps(MultiBlock3D<T>& multiBlock) const =0;
    virtual void signalPeriodicity() const =0;
};

}  // namespace plb

#endif
