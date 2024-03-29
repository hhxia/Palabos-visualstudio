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

#ifndef BENCHMARK_UTIL_HH
#define BENCHMARK_UTIL_HH

#include "core/globalDefs.h"
#include "io/parallelIO.h"
#include "algorithm/benchmarkUtil.h"
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <float.h>
#include <cassert>


namespace plb {

namespace util {


/////////// Class ValueTracer ////////////////////////

template<typename T>
ValueTracer<T>::ValueTracer(T u, T L, T _epsilon)
    : deltaT((plint)(L/u/2.)),
      epsilon(_epsilon),
      t(0),
      converged(false)
{ }

template<typename T>
plint ValueTracer<T>::getDeltaT() const {
    return deltaT;
}

template<typename T>
void ValueTracer<T>::takeValue(T val, bool doPrint) {
    values.push_back(val);
    if ((plint)values.size() > abs(deltaT)) {
        values.erase(values.begin());
        if (doPrint && t%deltaT==0) {
            T average = computeAverage();
            T stdDev = computeStdDev(average);
            pcout << "average=" << average << "; stdDev/average=" << stdDev/average << std::endl;
        }
    }
    ++t;
}

template<typename T>
void ValueTracer<T>::resetScale(T u, T L) {
    t = t%deltaT;
    deltaT = (plint) (L/u/2.);
    if ( (plint)values.size() > abs(deltaT) ) {
        values.erase(values.begin(), values.begin() + (values.size()-deltaT) );
    }
}

template<typename T>
void ValueTracer<T>::resetValues() {
    t = 0;
    if ((plint)values.size() > 0) {
        values.erase(values.begin(), values.begin() + values.size() );
    }
}

template<typename T>
bool ValueTracer<T>::hasConverged() const {
    if ((plint)values.size() < abs(deltaT)) {
        return false;
    }
    else {
        T average = computeAverage();
        T stdDev = computeStdDev(average);
		if (!_isnan(stdDev/average))
            return stdDev/average < epsilon;
        else {
            pcout << "simulation diverged.\n";
            return true;
        }
    }
}

template<typename T>
bool ValueTracer<T>::hasConvergedMinMax() const {
    if ((plint)values.size() < abs(deltaT)) {
        return false;
    }
    else {
        T minEl = *min_element(values.begin(), values.end());
        T maxEl = *max_element(values.begin(), values.end());
        T average = computeAverage();
        return (maxEl - minEl)/average < epsilon;
    }
}

template<typename T>
T ValueTracer<T>::computeAverage() const {
    return accumulate(values.begin(), values.end(), 0.) / values.size();
}

template<typename T>
T ValueTracer<T>::computeStdDev(T average) const {
    plint n = values.size();
    T sqrDev = 0.;
    for (plint i=0; i<n; ++i) {
        sqrDev += (values[i]-average)*(values[i]-average);
    }
    return sqrt(sqrDev/(n-1));
}

template<typename T>
void ValueTracer<T>::setEpsilon(T epsilon_) {
    epsilon = epsilon_;
}



/////////// Class BisectStepper ////////////////////////

template<typename T>
BisectStepper<T>::BisectStepper(T _iniVal, T _step)
    : iniVal(_iniVal), step(_step), state(first)
{
    if (step==0.) {
        step = iniVal/5.;
    }
    assert(step>0.);
}

template<typename T>
T BisectStepper<T>::getVal(bool stable, bool doPrint) {
    std::stringstream message;
    switch(state) {
    case first:
        if(stable) {
            currentVal = iniVal+step;
            state = up;
            message << "[" << iniVal << ",infty]";
        }
        else {
            currentVal = iniVal-step;
            state = down;
            message << "[-infty," << iniVal << "]";
        }
        break;
    case up:
        if (stable) {
            message << "[" << currentVal << ",infty]";
            currentVal += step;
        }
        else {
            lowerVal = currentVal-step;
            upperVal = currentVal;
            currentVal = 0.5*(lowerVal+upperVal);
            state = bisect;
            message << "[" << lowerVal << "," << upperVal << "]";
        }
        break;
    case down:
        if (!stable) {
            message << "[-infty," << currentVal << "]";
            currentVal -= step;
        }
        else {
            lowerVal = currentVal;
            upperVal = currentVal+step;
            currentVal = 0.5*(lowerVal+upperVal);
            state = bisect;
            message << "[" << lowerVal << "," << upperVal << "]";
        }
        break;
    case bisect:
        if (stable) {
            lowerVal = currentVal;
        }
        else {
            upperVal = currentVal;
        }
        currentVal = 0.5*(lowerVal+upperVal);
        message << "[" << lowerVal << "," << upperVal << "]";
        break;
    }
    if (doPrint) {
        pcout << "Value in range " << message.str() << std::endl;
    }
    return currentVal;
}

template<typename T>
bool BisectStepper<T>::hasConverged(T epsilon) const {
  return (state==bisect) && ((upperVal-lowerVal)/lowerVal < epsilon);
}

} // namespace util

} // namespace plb

#endif
