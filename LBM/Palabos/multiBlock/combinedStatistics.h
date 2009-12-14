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
 * The CombinedStatistics class -- header file.
 */
#ifndef COMBINED_STATISTICS_H
#define COMBINED_STATISTICS_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include <vector>

namespace plb {

template<typename T>
class CombinedStatistics {
public:
    virtual ~CombinedStatistics();
    virtual CombinedStatistics<T>* clone() const =0;
    void combine (
            std::vector<BlockStatistics<T> const*>& individualStatistics,
            BlockStatistics<T>& result ) const;
protected:
    virtual void reduceStatistics (
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights,
            std::vector<T>& sumObservables,
            std::vector<T>& maxObservables,
            std::vector<plint>& intSumObservables ) const =0;
private:
    void computeLocalAverage (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights ) const;
    void computeLocalSum (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<T>& sumObservables ) const;
    void computeLocalMax (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<T>& maxObservables ) const;
    void computeLocalIntSum (
            std::vector<BlockStatistics<T> const*> const& individualStatistics,
            std::vector<plint>& intSumObservables ) const;
};

template<typename T>
class SerialCombinedStatistics : public CombinedStatistics<T> {
public:
    virtual SerialCombinedStatistics<T>* clone() const;
protected:
    virtual void reduceStatistics (
            std::vector<T>& averageObservables,
            std::vector<T>& sumWeights,
            std::vector<T>& sumObservables,
            std::vector<T>& maxObservables,
            std::vector<plint>& intSumObservables ) const;
};

}  // namespace plb

#endif
