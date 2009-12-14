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

#ifndef BOUNDARY_PROCESSOR_2D_HH
#define BOUNDARY_PROCESSOR_2D_HH

#include "boundaryCondition/finiteDifferenceBoundaryProcessor2D.h"
#include "finiteDifference/finiteDifference2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <typeinfo>

namespace plb {

///////////  StraightFdBoundaryProcessor2D ///////////////////////////////////

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>::
    StraightFdBoundaryProcessor2D(Box2D domain_, BlockLattice2D<T,Descriptor>& lattice_)
    : domain(domain_), lattice(lattice_)
{
    PLB_PRECONDITION(domain.x0==domain.x1 || domain.y0==domain.y1);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>::process() {
    typedef SymmetricTensorImpl<T,Descriptor<T>::d> S;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        Array<T,Descriptor<T>::d> dx_u, dy_u;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor>& cell = lattice.get(iX,iY);
            Dynamics<T,Descriptor>& dynamics = cell.getDynamics();

            T rho = cell.computeDensity();
            Array<T,Descriptor<T>::d> u;
            cell.computeVelocity(u);
            interpolateGradients<0>(lattice, dx_u, iX, iY);
            interpolateGradients<1>(lattice, dy_u, iX, iY);
            T dx_ux = dx_u[0];
            T dy_ux = dy_u[0];
            T dx_uy = dx_u[1];
            T dy_uy = dy_u[1];
            T omega = cell.getDynamics().getOmega();
            T sToPi = - rho / Descriptor<T>::invCs2 / omega;
            Array<T,SymmetricTensor<T,Descriptor>::n> pi;
            pi[S::xx] = (T)2 * dx_ux * sToPi;
            pi[S::yy] = (T)2 * dy_uy * sToPi;
            pi[S::xy] = (dx_uy + dy_ux) * sToPi;

            Array<T,Descriptor<T>::d> j;
            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                j[iD] = rho*u[iD];
            }
            // Computation of the particle distribution functions
            // according to the regularized formula
            T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
              cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
                               offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>*
    StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>::clone() const
{
    return new StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
template<int deriveDirection>
void StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation>::
    interpolateGradients(BlockLattice2D<T,Descriptor> const& lattice,
                         Array<T,Descriptor<T>::d>& velDeriv, plint iX, plint iY) const
{
    fd::DirectedGradients2D<T,Descriptor,direction,orientation,direction==deriveDirection>::
        o1_velocityDerivative(velDeriv, lattice, iX, iY);
}

////////  StraightFdBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
StraightFdBoundaryProcessorGenerator2D<T,Descriptor, direction,orientation>::
    StraightFdBoundaryProcessorGenerator2D(Box2D domain)
    : BoxedDataProcessorGenerator2D<T>(domain)
{ }

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
DataProcessor2D<T>*
    StraightFdBoundaryProcessorGenerator2D<T,Descriptor,direction,orientation>::generate(std::vector<AtomicBlock2D<T>*> objects) const
{
    PLB_PRECONDITION( typeid(*objects[0]) == typeid(BlockLattice2D<T,Descriptor>) );
    return new StraightFdBoundaryProcessor2D<T,Descriptor,direction,orientation> (
            this->getDomain(),
            *dynamic_cast<BlockLattice2D<T,Descriptor>*> (objects[0]) );
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
DataProcessorGenerator2D<T>*
    StraightFdBoundaryProcessorGenerator2D<T,Descriptor,direction,orientation>::clone() const
{
    return new
        StraightFdBoundaryProcessorGenerator2D<T,Descriptor,direction,orientation> (this->getDomain());
}
 
/////////// OuterVelocityCornerProcessor2D /////////////////////////////////////

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>::
    OuterVelocityCornerProcessor2D(plint x_, plint y_, BlockLattice2D<T,Descriptor>& lattice_)
    : x(x_), y(y_),
      lattice(lattice_)
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
void OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>::process() {
    typedef SymmetricTensorImpl<T,Descriptor<T>::d> S;

    T rho10 = lattice.get(x-1*xNormal, y-0*yNormal).computeDensity();
    T rho01 = lattice.get(x-0*xNormal, y-1*yNormal).computeDensity();

    T rho = (T)1/(T)2*(rho01+rho10);
 
    Array<T,Descriptor<T>::d> dx_u, dy_u;
    fd::DirectedGradients2D<T, Descriptor, 0, xNormal, true>::o1_velocityDerivative(dx_u, lattice, x,y);
    fd::DirectedGradients2D<T, Descriptor, 1, yNormal, true>::o1_velocityDerivative(dy_u, lattice, x,y);
    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];

    Cell<T,Descriptor>& cell = lattice.get(x,y);
    Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
    T omega = dynamics.getOmega();

    T sToPi = - rho / Descriptor<T>::invCs2 / omega;
    Array<T,SymmetricTensor<T,Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;

    // Computation of the particle distribution functions
    // according to the regularized formula
    Array<T,Descriptor<T>::d> u, j;
    lattice.get(x,y).computeVelocity(u);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = rho*u[iD];
    }
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] =
            dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
                offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
    }
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>*
    OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>::clone() const
{
    return new OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>(*this);
}


////////  OuterVelocityCornerProcessorGenerator2D ////////////////////////////

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
OuterVelocityCornerProcessorGenerator2D<T, Descriptor, xNormal, yNormal>::
    OuterVelocityCornerProcessorGenerator2D(plint x_, plint y_)
        : BoxedDataProcessorGenerator2D<T>(Box2D(x_, x_, y_, y_))
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
DataProcessor2D<T>*
OuterVelocityCornerProcessorGenerator2D<T, Descriptor, xNormal, yNormal>::generate (
        std::vector<AtomicBlock2D<T>*> objects ) const
{
    PLB_PRECONDITION( typeid(*objects[0]) == typeid(BlockLattice2D<T,Descriptor>) );
    return new OuterVelocityCornerProcessor2D<T, Descriptor, xNormal, yNormal>
                   ( this->getDomain().x0, this->getDomain().y0,
                     *dynamic_cast<BlockLattice2D<T,Descriptor>*> (objects[0]) );
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
OuterVelocityCornerProcessorGenerator2D<T, Descriptor, xNormal, yNormal>::clone() const
{
    return new OuterVelocityCornerProcessorGenerator2D<T, Descriptor, xNormal, yNormal>
               ( this->getDomain().x0, this->getDomain().y0);
}

}  // namespace plb

#endif  // BOUNDARY_PROCESSOR_2D_HH
