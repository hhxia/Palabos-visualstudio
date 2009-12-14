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
 * Handler for 3D multi data field structure -- header file.
 */

#ifndef MULTI_DATA_FIELD_HANDLER_3D_H
#define MULTI_DATA_FIELD_HANDLER_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockManagement3D.h"
#include <vector>

namespace plb {

template<typename T>
struct MultiDataFieldHandler3D {
    virtual ~MultiDataFieldHandler3D() { }
    virtual plint getNx() const =0;
    virtual plint getNy() const =0;
    virtual plint getNz() const =0;
    virtual MultiBlockDistribution3D const& getMultiBlockDistribution() const =0;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly, plint& lz) const =0;
    virtual T reduceSum(T localSum) const =0;
    virtual T reduceAverage(T localAverage, T localWeight) const =0;
    virtual T reduceMin(T localMin) const =0;
    virtual T reduceMax(T localMax) const =0;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const =0;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const =0;
};

template<typename T>
class SerialMultiDataFieldHandler3D : public MultiDataFieldHandler3D<T> {
public:
    SerialMultiDataFieldHandler3D(MultiBlockDistribution3D const& multiBlockDistribution_);
    virtual plint getNx() const;
    virtual plint getNy() const;
    virtual plint getNz() const;
    virtual MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly, plint& lz) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const;
private:
    MultiBlockDistribution3D multiBlockDistribution;
};


#ifdef PLB_MPI_PARALLEL
template<typename T>
class ParallelMultiDataFieldHandler3D : public MultiDataFieldHandler3D<T> {
public:
    ParallelMultiDataFieldHandler3D(MultiBlockDistribution3D const& multiBlockDistribution_);
    virtual plint getNx() const;
    virtual plint getNy() const;
    virtual plint getNz() const;
    virtual MultiBlockDistribution3D const& getMultiBlockDistribution() const;
    virtual bool getLocalEnvelope(plint iBlock, plint& lx, plint& ly, plint& lz) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, plint fromBlock) const;
    virtual void broadCastVector(T* vect, plint size, plint fromBlock) const;
private:
    MultiBlockDistribution3D multiBlockDistribution;
};
#endif

}  // namespace plb

#endif
