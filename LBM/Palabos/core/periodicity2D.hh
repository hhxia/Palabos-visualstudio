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
 * Classes and methods for handling periodicity in 2D -- implementation.
 */

#ifndef PERIODICITY_2D_HH
#define PERIODICITY_2D_HH

#include "core/periodicity2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"

namespace plb {

template<typename T>
PeriodicitySwitch2D<T>::PeriodicitySwitch2D(Block2D<T>& block_)
    : block(block_)
{
    // Default all boundaries to non-periodic.
    toggleAll(false);
}

template<typename T>
PeriodicitySwitch2D<T>::PeriodicitySwitch2D(Block2D<T>& block_, PeriodicitySwitch2D<T> const& rhs)
    : periodicity(rhs.periodicity),
      block(block_)
{ }

template<typename T>
PeriodicitySwitch2D<T>& PeriodicitySwitch2D<T>::operator=(PeriodicitySwitch2D<T> const& rhs)
{
    // Don't modify the Block2D reference, as it is immutable.
    periodicity = rhs.periodicity;
    return *this;
}

template<typename T>
void PeriodicitySwitch2D<T>::toggle(plint direction, bool periodic) {
    PLB_PRECONDITION( direction==0 || direction==1 );
    periodicity[direction] = periodic;
    block.signalPeriodicity();
}

template<typename T>
void PeriodicitySwitch2D<T>::toggleAll(bool periodic) {
    toggle(0, periodic);
    toggle(1, periodic);
    block.signalPeriodicity();
}

template<typename T>
bool PeriodicitySwitch2D<T>::get(plint direction) const {
    PLB_PRECONDITION( direction==0 || direction==1 );

    return periodicity[direction];
}

template<typename T>
bool PeriodicitySwitch2D<T>::get(plint normalX, plint normalY) const {
    bool testX = normalX != 0;
    bool testY = normalY != 0;
    return ( (!testX || (testX && periodicity[0]) ) &&
             (!testY || (testY && periodicity[1]) ) );
}

template<typename T>
Box2D PeriodicitySwitch2D<T>::getPeriodicEnvelope(Box2D const& bulk, plint envelopeWidth) const {
    Box2D envelope(bulk);
    if (get(0)) { // If periodic in x, extend bulk in x-direction
        envelope.x0 -= envelopeWidth;
        envelope.x1 += envelopeWidth;
    }
    if (get(1)) { // If periodic in y, extend bulk in y-direction
        envelope.y0 -= envelopeWidth;
        envelope.y1 += envelopeWidth;
    }
    return envelope;
}

} // namespace plb

#endif  //PERIODICITY_2D_HH
