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

#include "parallelism/mpiManager.h"
#include "core/plbTimer.h"
#include <map>

namespace plb {

namespace global {

XflTimer::XflTimer()
    : cumulativeTime(0.),
      isOn(false)
{ }

void XflTimer::start() {
#ifdef PLB_MPI_PARALLEL
    startTime = mpi().getTime();
#else
    startClock = clock();
#endif
    isOn = true;
}

void XflTimer::restart() {
    reset();
    start();
}

double XflTimer::stop() {
    cumulativeTime = getTime();
    isOn = false;
    return cumulativeTime;
}

void XflTimer::reset() {
    cumulativeTime = 0.;
}

double XflTimer::getTime() const {
    if (isOn) {
#ifdef PLB_MPI_PARALLEL
        return cumulativeTime + mpi().getTime()-startTime;
#else
        return cumulativeTime + (double)(clock()-startClock)
                              / (double)CLOCKS_PER_SEC;
#endif
    }
    else {
        return cumulativeTime;
    }
}

XflTimer& timer(std::string nameOfTimer) {
    static std::map<std::string, XflTimer> timerCollection;
    return timerCollection[nameOfTimer];
}

}  // namespace global

}  // namespace plb
