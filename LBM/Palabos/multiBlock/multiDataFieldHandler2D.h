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
 * Handler for 2D multi data field structure -- header file.
 */

#ifndef MULTI_DATA_FIELD_HANDLER_2D_H
#define MULTI_DATA_FIELD_HANDLER_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement2D.h"
#include <vector>

namespace plb {

template<typename T>
struct MultiDataFieldHandler2D {
    virtual ~MultiDataFieldHandler2D() { }
    virtual plint getNx() const =0;
    virtual plint getNy() const =0;
    virtual MultiBlockDistribution2D const& getMultiBlockDistribution() const =0;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly) const =0;
    virtual T reduceSum(T localSum) const =0;
    virtual T reduceAverage(T localAverage, T localWeight) const =0;
    virtual T reduceMin(T localMin) const =0;
    virtual T reduceMax(T localMax) const =0;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const =0;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const =0;
};

template<typename T>
class SerialMultiDataFieldHandler2D : public MultiDataFieldHandler2D<T> {
public:
    SerialMultiDataFieldHandler2D(MultiBlockDistribution2D const& multiBlockDistribution_);
    virtual plint getNx() const;
    virtual plint getNy() const;
    virtual MultiBlockDistribution2D const& getMultiBlockDistribution() const;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const;
private:
    MultiBlockDistribution2D multiBlockDistribution;
};


#ifdef PLB_MPI_PARALLEL
template<typename T>
class ParallelMultiDataFieldHandler2D : public MultiDataFieldHandler2D<T> {
public:
    ParallelMultiDataFieldHandler2D(MultiBlockDistribution2D const& multiBlockDistribution_);
    virtual plint getNx() const;
    virtual plint getNy() const;
    virtual MultiBlockDistribution2D const& getMultiBlockDistribution() const;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const;
private:
    MultiBlockDistribution2D multiBlockDistribution;
};
#endif

}  // namespace plb

#endif
