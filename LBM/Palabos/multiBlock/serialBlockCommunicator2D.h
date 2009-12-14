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
 * Serial version of the 2D block communicator -- header file.
 */
#ifndef SERIAL_BLOCK_COMMUNICATOR_2D_H
#define SERIAL_BLOCK_COMMUNICATOR_2D_H

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

class MultiBlockManagement2D;

template<typename T>
class SerialBlockCommunicator2D : public BlockCommunicator2D<T> {
public:
    SerialBlockCommunicator2D();
    virtual SerialBlockCommunicator2D<T>* clone() const;
    virtual void broadCastScalar(T& scalar, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const;
    virtual void broadCastVector(T* data, plint length, plint fromBlock, MultiBlockManagement2D const& multiBlockManagement) const;
    virtual void duplicateOverlaps(MultiBlock2D<T>& multiBlock) const;
    virtual void signalPeriodicity() const;
private:
    void copyOverlap(Overlap2D const& overlap, MultiBlock2D<T>& multiBlock) const;
};

}  // namespace plb

#endif
