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
 * Classes and methods for handling periodicity in 3D -- implementation.
 */

#ifndef PERIODICITY_3D_HH
#define PERIODICITY_3D_HH

#include "core/periodicity3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"

namespace plb {

template<typename T>
PeriodicitySwitch3D<T>::PeriodicitySwitch3D(Block3D<T>& block_)
    : block(block_)
{
    // Default all boundaries to non-periodic.
    toggleAll(false);
}

template<typename T>
PeriodicitySwitch3D<T>::PeriodicitySwitch3D(Block3D<T>& block_, PeriodicitySwitch3D<T> const& rhs)
    : periodicity(rhs.periodicity),
      block(block_)
{ }

template<typename T>
PeriodicitySwitch3D<T>& PeriodicitySwitch3D<T>::operator=(PeriodicitySwitch3D<T> const& rhs)
{
    // Don't modify the Block3D reference, as it is immutable.
    periodicity = rhs.periodicity;
    return *this;
}

template<typename T>
void PeriodicitySwitch3D<T>::toggle(plint direction, bool periodic) {
    PLB_PRECONDITION( direction==0 || direction==1 || direction==2 );
    periodicity[direction] = periodic;
    block.signalPeriodicity();
}

template<typename T>
void PeriodicitySwitch3D<T>::toggleAll(bool periodic) {
    toggle(0, periodic);
    toggle(1, periodic);
    toggle(2, periodic);
    block.signalPeriodicity();
}

template<typename T>
bool PeriodicitySwitch3D<T>::get(plint direction) const {
    PLB_PRECONDITION( direction==0 || direction==1 || direction==2 );

    return periodicity[direction];
}

template<typename T>
bool PeriodicitySwitch3D<T>::get(plint normalX, plint normalY, plint normalZ) const {
    bool testX = normalX != 0;
    bool testY = normalY != 0;
    bool testZ = normalZ != 0;
    return ( (!testX || (testX && periodicity[0]) ) &&
             (!testY || (testY && periodicity[1]) ) &&
             (!testZ || (testZ && periodicity[2]) ) );
}

template<typename T>
Box3D PeriodicitySwitch3D<T>::getPeriodicEnvelope(Box3D const& bulk, plint envelopeWidth) const {
    Box3D envelope(bulk);
    if (get(0)) { // If periodic in x, extend bulk in x-direction
        envelope.x0 -= envelopeWidth;
        envelope.x1 += envelopeWidth;
    }
    if (get(1)) { // If periodic in y, extend bulk in y-direction
        envelope.y0 -= envelopeWidth;
        envelope.y1 += envelopeWidth;
    }
    if (get(2)) { // If periodic in y, extend bulk in y-direction
        envelope.z0 -= envelopeWidth;
        envelope.z1 += envelopeWidth;
    }
    return envelope;
}

} // namespace plb

#endif  //PERIODICITY_3D_HH
