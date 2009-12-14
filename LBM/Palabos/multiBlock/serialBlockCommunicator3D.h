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
 * Serial version of the 3D block communicator -- header file.
 */
#ifndef SERIAL_BLOCK_COMMUNICATOR_3D_H
#define SERIAL_BLOCK_COMMUNICATOR_3D_H

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

class MultiBlockManagement3D;

template<typename T>
class SerialBlockCommunicator3D : public BlockCommunicator3D<T> {
public:
    SerialBlockCommunicator3D();
    virtual SerialBlockCommunicator3D<T>* clone() const;
    virtual void broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const;
    virtual void broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement3D const& multiBlockManagement) const;
    virtual void duplicateOverlaps(MultiBlock3D<T>& multiBlock) const;
    virtual void signalPeriodicity() const;
private:
    void copyOverlap(Overlap3D const& overlap, MultiBlock3D<T>& multiBlock) const;
};

}  // namespace plb

#endif
