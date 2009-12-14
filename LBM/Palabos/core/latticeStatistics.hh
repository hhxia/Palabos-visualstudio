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


#ifndef LATTICE_STATISTICS_HH
#define LATTICE_STATISTICS_HH

#include "core/latticeStatistics.h"
#include "core/block2D.h"
#include "core/block3D.h"

namespace plb {

template<typename T>
inline void gatherStatistics(BlockStatistics<T>& statistics, T rhoBar, T uSqr) {
    statistics.gatherAverage(LatticeStatistics::avRhoBar, rhoBar);
    statistics.gatherAverage(LatticeStatistics::avUSqr, uSqr);
    statistics.gatherMax(LatticeStatistics::maxUSqr, uSqr);
    statistics.incrementStats();
}

} // namespace plb

#endif
