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
 * A timer class for benchmarking program parts -- header file.
 */
#ifndef PLB_TIMER_H
#define PLB_TIMER_H

#include <ctime>

namespace plb {

namespace global {

/// A cumulative timer for benchmarking program parts and summing up the times.
/** In serial, the C function clock() is used. In parallel, the MPI
 *  function is used.
 */
class XflTimer {
public:
    XflTimer();
    /// Proceed with time measurement.
    void start();
    /// Reset timer to zero and start time measurement.
    void restart();
    /// Interrupt time measurement ( you can still proceed with start() ).
    /* \return Current cumulative time.
     */
    double stop();
    /// Reset timer to zero.
    void reset();
    /// Get current cumulative time.
    double getTime() const;
private:
    double cumulativeTime;
    bool   isOn;
#ifdef PLB_MPI_PARALLEL
    double startTime;
#else
    clock_t startClock;
#endif
friend XflTimer& timer(std::string nameOfTimer);
};

XflTimer& timer(std::string nameOfTimer);

}  // namespace global

}  // namespace plb

#endif
