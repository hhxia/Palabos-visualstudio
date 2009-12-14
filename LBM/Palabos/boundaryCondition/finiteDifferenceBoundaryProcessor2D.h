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

#ifndef BOUNDARY_PROCESSOR_2D_H
#define BOUNDARY_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/blockLattice2D.h"

namespace plb {

/**
* This class computes the Skordos BC
* on a flat wall in 2D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class StraightFdBoundaryProcessor2D : public DataProcessor2D<T> {
public:
    StraightFdBoundaryProcessor2D(Box2D domain_, BlockLattice2D<T,Descriptor>& lattice_);
    virtual void process();
    virtual StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>* clone() const;
private:
    template<int deriveDirection>
    void interpolateGradients (
            BlockLattice2D<T,Descriptor> const& lattice,
            Array<T,Descriptor<T>::d>& velDeriv, plint iX, plint iY ) const;
private:
    Box2D domain;
    BlockLattice2D<T,Descriptor>& lattice;
};

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class StraightFdBoundaryProcessorGenerator2D : public BoxedDataProcessorGenerator2D<T>
{
public:
    StraightFdBoundaryProcessorGenerator2D(Box2D domain);
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> objects) const;
    virtual DataProcessorGenerator2D<T>* clone() const;
};

/**
* This class computes the Skordos BC in 2D on a convex
* corner but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
class OuterVelocityCornerProcessor2D : public DataProcessor2D<T> {
public:
    OuterVelocityCornerProcessor2D(plint x_, plint y_, BlockLattice2D<T,Descriptor>& lattice);
    virtual plint extent() const { return 2; }
    virtual plint extent(int whichDirection) const { return 2; }
    virtual void process();
    virtual OuterVelocityCornerProcessor2D<T,Descriptor,xNormal,yNormal>* clone() const;
private:
    plint x, y;
    BlockLattice2D<T,Descriptor>& lattice;
};

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
class OuterVelocityCornerProcessorGenerator2D : public BoxedDataProcessorGenerator2D<T>
{
public:
    OuterVelocityCornerProcessorGenerator2D(plint x_, plint y_);
    virtual DataProcessor2D<T>* generate(std::vector<AtomicBlock2D<T>*> objects) const;
    virtual DataProcessorGenerator2D<T>*  clone() const;
};

}

#endif  // BOUNDARY_PROCESSOR_2D_H
