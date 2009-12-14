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
 * 3D Blocks -- generic implementation.
 */
#ifndef BLOCK_3D_HH
#define BLOCK_3D_HH

#include "core/block3D.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

template<typename T>
Block3D<T>::Block3D()
    : internalStatistics(),
      periodicitySwitch(*this)
{ }

template<typename T>
Block3D<T>::Block3D(Block3D<T> const& rhs)
  : internalStatistics(rhs.internalStatistics),
    periodicitySwitch(*this, rhs.periodicitySwitch)
{ }

template<typename T>
void Block3D<T>::swap(Block3D<T>& rhs) {
    std::swap(internalStatistics, rhs.internalStatistics);
    std::swap(periodicitySwitch, rhs.periodicitySwitch);
}

template<typename T>
BlockStatistics<T>& Block3D<T>::getInternalStatistics() {
    return internalStatistics;
}

template<typename T>
PeriodicitySwitch3D<T> const& Block3D<T>::periodicity() const {
    return periodicitySwitch;
}

template<typename T>
PeriodicitySwitch3D<T>& Block3D<T>::periodicity() {
    return periodicitySwitch;
}

template<typename T>
plint Block3D<T>::getNx() const {
    return getBoundingBox().getNx();
}

template<typename T>
plint Block3D<T>::getNy() const {
    return getBoundingBox().getNy();
}

template<typename T>
plint Block3D<T>::getNz() const {
    return getBoundingBox().getNz();
}

template<typename T>
BlockStatistics<T> const& Block3D<T>::getInternalStatistics() const {
    return internalStatistics;
}

template<typename T>
void Block3D<T>::signalPeriodicity()
{ }

template <typename T>
void copySerializedBlock(Block3D<T> const& from, Block3D<T>& to, IndexOrdering::OrderingT ordering) {
    PLB_PRECONDITION( from.getBoundingBox().nCells() == to.getBoundingBox().nCells() );
    serializerToUnSerializer( from.getBlockSerializer(from.getBoundingBox(), ordering),
                              to.getBlockUnSerializer(to.getBoundingBox(), ordering) );
}

}  // namespace plb

#endif
