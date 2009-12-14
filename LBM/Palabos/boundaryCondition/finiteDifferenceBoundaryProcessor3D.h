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

#ifndef BOUNDARY_PROCESSOR_3D_H
#define BOUNDARY_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

/**
* This class computes the Skordos BC
* on a plane wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class PlaneFdBoundaryProcessor3D : public DataProcessor3D<T>
{
public:
    PlaneFdBoundaryProcessor3D(Box3D domain_, BlockLattice3D<T,Descriptor>& lattice_);
    virtual void process();
    virtual PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>* clone() const;
private:
    template<int deriveDirection>
    void interpolateGradients (
            BlockLattice3D<T,Descriptor> const& block,
            Array<T,Descriptor<T>::d>& velDeriv, plint iX, plint iY, plint iZ ) const;
private:
    Box3D domain;
    BlockLattice3D<T,Descriptor>& lattice;
};

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class PlaneFdBoundaryProcessorGenerator3D : public BoxedDataProcessorGenerator3D<T>
{
public:
    PlaneFdBoundaryProcessorGenerator3D(Box3D domain);
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> objects) const;
    virtual DataProcessorGenerator3D<T>*  clone() const;
};

/**
* This class computes the Skordos BC
* on a convex edge wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Descriptor,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessor3D : public DataProcessor3D<T> {
public:
    enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
    OuterVelocityEdgeProcessor3D(Box3D domain_, BlockLattice3D<T,Descriptor>& lattice_);
    virtual plint extent() const { return 2; }
    virtual plint extent(int whichDirection) const { return 2; }
    virtual void process();
    virtual OuterVelocityEdgeProcessor3D<T,Descriptor,plane,normal1,normal2>* clone() const;
private:
    T getNeighborRho(plint x, plint y, plint z, plint step1, plint step2,
                     BlockLattice3D<T,Descriptor> const& lattice);
    template<int deriveDirection, int orientation>
    void interpolateGradients (
            BlockLattice3D<T,Descriptor> const& lattice,
            Array<T,Descriptor<T>::d>& velDeriv, plint iX, plint iY, plint iZ ) const;
private:
    Box3D domain;
    BlockLattice3D<T,Descriptor>& lattice;
};

template<typename T, template<typename U> class Descriptor,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessorGenerator3D
    : public BoxedDataProcessorGenerator3D<T>
{
public:
    OuterVelocityEdgeProcessorGenerator3D(Box3D domain);
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> objects) const;
    virtual DataProcessorGenerator3D<T>* clone() const;
};

template<typename T, template<typename U> class Descriptor,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessor3D : public DataProcessor3D<T> {
public:
    OuterVelocityCornerProcessor3D(plint x_, plint y_, plint z_, BlockLattice3D<T,Descriptor>& lattice_);
    virtual plint extent() const { return 2; }
    virtual plint extent(int whichDirection) const { return 2; }
    virtual void process();
    virtual OuterVelocityCornerProcessor3D<T,Descriptor,xNormal,yNormal,zNormal>* clone() const;
private:
    plint x,y,z;
    BlockLattice3D<T,Descriptor>& lattice;
};

template<typename T, template<typename U> class Descriptor,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessorGenerator3D
    : public BoxedDataProcessorGenerator3D<T>
{
public:
    OuterVelocityCornerProcessorGenerator3D(plint x_, plint y_, plint z_);
    virtual DataProcessor3D<T>* generate(std::vector<AtomicBlock3D<T>*> objects) const;
    virtual DataProcessorGenerator3D<T>* clone() const;
};

}

#endif  // BOUNDARY_PROCESSOR_3D_H
