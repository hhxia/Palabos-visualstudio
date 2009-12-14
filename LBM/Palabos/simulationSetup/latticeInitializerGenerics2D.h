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
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef LATTICE_INITIALIZER_GENERICS_2D_H
#define LATTICE_INITIALIZER_GENERICS_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, class VelocityFunction>
class SetCustomBoundaryVelocityFunctional2D : public OneCellIndexedFunctional2D<T,Descriptor>
{
public:
    SetCustomBoundaryVelocityFunctional2D(VelocityFunction f_)
        : f(f_),
          velocityScale( (T)1 )
    { }
    virtual SetCustomBoundaryVelocityFunctional2D<T,Descriptor,VelocityFunction>* clone() const {
        return new SetCustomBoundaryVelocityFunctional2D<T,Descriptor,VelocityFunction>(*this);
    }
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const {
        Array<T,Descriptor<T>::d> u;
        f(iX, iY, u);
        u[0] *= velocityScale;
        u[1] *= velocityScale;
        cell.defineVelocity(u);
    }
    virtual void rescale(T dxScale, T dtScale) {
        velocityScale *= dtScale / dxScale;
    }
private:
    VelocityFunction f;
    T velocityScale;
};

template<typename T, template<typename U> class Descriptor, class DensityFunction>
class SetCustomBoundaryDensityFunctional2D : public OneCellIndexedFunctional2D<T,Descriptor>
{
public:
    SetCustomBoundaryDensityFunctional2D(DensityFunction f_) : f(f_)
    { }
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const {
        // No rescaling needed: rho is scale invariant.
        T rho = f(iX, iY);
        cell.defineDensity(rho);
    }
    virtual SetCustomBoundaryDensityFunctional2D<T,Descriptor,DensityFunction>* clone() const {
        return new SetCustomBoundaryDensityFunctional2D<T,Descriptor,DensityFunction>(*this);
    }
private:
    DensityFunction f;
};

template<typename T, template<typename U> class Descriptor, class RhoUFunction>
class IniCustomEquilibriumFunctional2D : public OneCellIndexedFunctional2D<T,Descriptor>
{
public:
    IniCustomEquilibriumFunctional2D(RhoUFunction f_)
        : f(f_),
          velocityScale( (T)1 )
    { }
    virtual void execute(plint iX, plint iY, Cell<T,Descriptor>& cell) const {
        Array<T,Descriptor<T>::d> j;
        T rho;
        f(iX, iY, rho, j);
        for (int iD=0; iD<Descriptor<T>::d; ++iD) {
            j[iD] *= rho;
            j[iD] *= velocityScale;
        }
        T rhoBar = Descriptor<T>::rhoBar(rho);
        T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
        for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
            cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
        }
    }
    virtual IniCustomEquilibriumFunctional2D<T,Descriptor,RhoUFunction>* clone() const {
        return new IniCustomEquilibriumFunctional2D<T,Descriptor,RhoUFunction>(*this);
    }
    virtual void rescale(T dxScale, T dtScale) {
        velocityScale *= dtScale / dxScale;
    }
private:
    RhoUFunction f;
    T velocityScale;
};

template<typename T, template<class U> class Descriptor, class VelocityFunction>
void setBoundaryVelocity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, VelocityFunction f) {
    applyIndexed (
            lattice, domain,
            new SetCustomBoundaryVelocityFunctional2D<T,Descriptor,VelocityFunction> (f) );
}

template<typename T, template<class U> class Descriptor, class DensityFunction>
void setBoundaryDensity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, DensityFunction f) {
    applyIndexed (
            lattice, domain,
            new SetCustomBoundaryDensityFunctional2D<T,Descriptor,DensityFunction> (f) );
}

template<typename T, template<class U> class Descriptor, class RhoUFunction>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, RhoUFunction f) {
    applyIndexed (
            lattice, domain,
            new IniCustomEquilibriumFunctional2D<T,Descriptor,RhoUFunction> (f) );
}

}  // namespace plb

#endif  // LATTICE_INITIALIZER_GENERICS_2D_H
