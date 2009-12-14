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
 * Classes and methods for handling periodicity in 2D -- header file.
 */
#ifndef PERIODICITY_2D_H
#define PERIODICITY_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/array.h"

namespace plb {

template<typename T> class Block2D;

/// Storage of the information whether a particular direction in a block is periodic.
template<typename T>
class PeriodicitySwitch2D {
public:
    /// The constructor defaults all directions to false (i.e., non-periodic).
    PeriodicitySwitch2D(Block2D<T>& block_);
    PeriodicitySwitch2D(Block2D<T>& block_, PeriodicitySwitch2D<T> const& rhs);
    PeriodicitySwitch2D<T>& operator=(PeriodicitySwitch2D<T> const& rhs);

    /// Set periodicity status of a direction (direction=0 means x-direction etc.)
    void toggle(plint direction, bool periodic);
    /// Set periodicity status of all directions synchronously
    void toggleAll(bool periodic);
    /// Get periodicity status of a direction;
    bool get(plint direction) const;
    /// Get periodicity along a general direction;
    bool get(plint normalX, plint normalY) const;
    /// Extend the bulk in each periodic direction, and return the result.
    Box2D getPeriodicEnvelope(Box2D const& bulk, plint envelopeWidth) const;
private:
    Array<bool,2> periodicity;
    Block2D<T>& block;
};

} // namespace plb

#endif  // PERIODICITY_2D_H
