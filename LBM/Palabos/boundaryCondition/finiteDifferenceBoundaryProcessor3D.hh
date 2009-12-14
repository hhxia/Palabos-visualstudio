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

#ifndef BOUNDARY_PROCESSOR_3D_HH
#define BOUNDARY_PROCESSOR_3D_HH

#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.h"
#include "finiteDifference/finiteDifference3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <typeinfo>

namespace plb {

////////  PlaneFdBoundaryProcessor3D ///////////////////////////////////

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>::
    PlaneFdBoundaryProcessor3D(Box3D domain_, BlockLattice3D<T,Descriptor>& lattice_)
    : domain(domain_),
      lattice(lattice_)
{
    PLB_ASSERT(domain.x0==domain.x1 || domain.y0==domain.y1 ||
               domain.z0==domain.z1);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>::process()
{
    typedef SymmetricTensorImpl<T,Descriptor<T>::d> S;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        Array<T,Descriptor<T>::d> dx_u, dy_u, dz_u;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
                T rho = cell.computeDensity();
                Array<T,Descriptor<T>:: d> u;
                cell.computeVelocity(u);

                interpolateGradients<0> ( lattice, dx_u, iX, iY, iZ );
                interpolateGradients<1> ( lattice, dy_u, iX, iY, iZ );
                interpolateGradients<2> ( lattice, dz_u, iX, iY, iZ );
                T dx_ux = dx_u[0];
                T dy_ux = dy_u[0];
                T dz_ux = dz_u[0];
                T dx_uy = dx_u[1];
                T dy_uy = dy_u[1];
                T dz_uy = dz_u[1];
                T dx_uz = dx_u[2];
                T dy_uz = dy_u[2];
                T dz_uz = dz_u[2];
                T omega = cell.getDynamics().getOmega();
                T sToPi = - rho / Descriptor<T>::invCs2 / omega;
                Array<T,SymmetricTensor<T,Descriptor>::n> pi;
                pi[S::xx] = (T)2 * dx_ux * sToPi;
                pi[S::yy] = (T)2 * dy_uy * sToPi;
                pi[S::zz] = (T)2 * dz_uz * sToPi;
                pi[S::xy] = (dx_uy + dy_ux) * sToPi;
                pi[S::xz] = (dx_uz + dz_ux) * sToPi;
                pi[S::yz] = (dy_uz + dz_uy) * sToPi;

                Array<T,Descriptor<T>::d> j;
                for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                    j[iD] = rho*u[iD];
                }
                T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

                // Computation of the particle distribution functions
                // according to the regularized formula
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
                    cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
                                     offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>*
     PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>::clone() const
{
    return new PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
template<int deriveDirection>
void PlaneFdBoundaryProcessor3D<T,Descriptor,direction,orientation>::
    interpolateGradients(BlockLattice3D<T,Descriptor> const& lattice, Array<T,Descriptor<T>::d>& velDeriv,
                         plint iX, plint iY, plint iZ) const
{
    fd::DirectedGradients3D<T, Descriptor, direction, orientation, deriveDirection, direction==deriveDirection>::
        o1_velocityDerivative(velDeriv, lattice, iX, iY, iZ);
}


////////  PlaneFdBoundaryProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
PlaneFdBoundaryProcessorGenerator3D<T,Descriptor,direction,orientation>::
    PlaneFdBoundaryProcessorGenerator3D(Box3D domain)
    : BoxedDataProcessorGenerator3D<T>(domain)
{ }

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
DataProcessor3D<T>* PlaneFdBoundaryProcessorGenerator3D<T,Descriptor,direction,orientation>::
        generate(std::vector<AtomicBlock3D<T>*> objects) const
{
    PLB_PRECONDITION( typeid(*objects[0]) == typeid(BlockLattice3D<T,Descriptor>) );
    return new PlaneFdBoundaryProcessor3D<T,Descriptor, direction,orientation>
        ( this->getDomain(),
          *dynamic_cast<BlockLattice3D<T,Descriptor>*> (objects[0]) );
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
DataProcessorGenerator3D<T>*
    PlaneFdBoundaryProcessorGenerator3D<T,Descriptor,direction,orientation>::clone() const
{
    return new PlaneFdBoundaryProcessorGenerator3D<T,Descriptor,direction,orientation> (this->getDomain());
}
 

////////  OuterVelocityEdgeProcessor3D ///////////////////////////////////

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>::
    OuterVelocityEdgeProcessor3D(Box3D domain_, BlockLattice3D<T,Descriptor>& lattice_)
    : domain(domain_),
      lattice(lattice_)
{
    PLB_ASSERT (
            (plane==2 && domain.x0==domain.x1 && domain.y0==domain.y1) ||
            (plane==1 && domain.x0==domain.x1 && domain.z0==domain.z1) ||
            (plane==0 && domain.y0==domain.y1 && domain.z0==domain.z1)     );

}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
void OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>::process() {
    typedef SymmetricTensorImpl<T,Descriptor<T>::d> S;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                Dynamics<T,Descriptor>& dynamics = cell.getDynamics();

                T rho10 = getNeighborRho(iX,iY,iZ,1,0, lattice);
                T rho01 = getNeighborRho(iX,iY,iZ,0,1, lattice);
                T rho20 = getNeighborRho(iX,iY,iZ,2,0, lattice);
                T rho02 = getNeighborRho(iX,iY,iZ,0,2, lattice);
                T rho = (T)2/(T)3*(rho01+rho10)-(T)1/(T)6*(rho02+rho20);

                std::vector<Array<T,3> > dA_uB_(3);
                interpolateGradients<plane,0>            ( lattice, dA_uB_[0], iX, iY, iZ );
                interpolateGradients<direction1,normal1> ( lattice, dA_uB_[1], iX, iY, iZ );
                interpolateGradients<direction2,normal2> ( lattice, dA_uB_[2], iX, iY, iZ );
                std::vector<Array<T,3> > dA_uB(3);
                for (int iBeta=0; iBeta<3; ++iBeta) {
                    dA_uB[plane][iBeta]      = dA_uB_[0][iBeta];
                    dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
                    dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
                }
                T omega = dynamics.getOmega();
                T sToPi = - rho / Descriptor<T>::invCs2 / omega;
                Array<T,SymmetricTensor<T,Descriptor>::n> pi;
                pi[S::xx] = (T)2 * dA_uB[0][0] * sToPi;
                pi[S::yy] = (T)2 * dA_uB[1][1] * sToPi;
                pi[S::zz] = (T)2 * dA_uB[2][2] * sToPi;
                pi[S::xy] = (dA_uB[0][1]+dA_uB[1][0]) * sToPi;
                pi[S::xz] = (dA_uB[0][2]+dA_uB[2][0]) * sToPi;
                pi[S::yz] = (dA_uB[1][2]+dA_uB[2][1]) * sToPi;

                // Computation of the particle distribution functions
                // according to the regularized formula
                Array<T,Descriptor<T>::d> u, j;
                cell.computeVelocity(u);

                for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                    j[iD] = rho*u[iD];
                }
                T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
                                     offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>*
    OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>::clone() const
{
    return new OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>(*this);
}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
T OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>::
    getNeighborRho(plint x, plint y, plint z, plint step1, plint step2, BlockLattice3D<T,Descriptor> const& lattice)
{
    Array<int,3> coords (x, y, z);
    coords[direction1] += -normal1*step1;
    coords[direction2] += -normal2*step2;
    return lattice.get(coords[0], coords[1], coords[2]).computeDensity();
}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
template<int deriveDirection, int orientation>
void OuterVelocityEdgeProcessor3D<T,Descriptor, plane,normal1,normal2>::
    interpolateGradients(BlockLattice3D<T,Descriptor> const& lattice,
                         Array<T,Descriptor<T>::d>& velDeriv,
                         plint iX, plint iY, plint iZ) const
{
    fd::DirectedGradients3D<T,Descriptor,deriveDirection,orientation,deriveDirection,deriveDirection!=plane>::
        o1_velocityDerivative(velDeriv, lattice, iX, iY, iZ);
}

////////  OuterVelocityEdgeProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>::
    OuterVelocityEdgeProcessorGenerator3D(Box3D domain)
    : BoxedDataProcessorGenerator3D<T>(domain)
{ }

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
DataProcessor3D<T>*
    OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>::
         generate(std::vector<AtomicBlock3D<T>*> objects) const
{
    PLB_PRECONDITION( typeid(*objects[0]) == typeid(BlockLattice3D<T,Descriptor>) );
    return new OuterVelocityEdgeProcessor3D < T,Descriptor, plane,normal1,normal2 >
        ( this->getDomain(),
          *dynamic_cast<BlockLattice3D<T,Descriptor>*> (objects[0]) );
}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>::clone() const
{
    return new OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2 > (this->getDomain());
}

/////////// OuterVelocityCornerProcessor3D /////////////////////////////////////

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessor3D<T, Descriptor, xNormal, yNormal, zNormal>::
    OuterVelocityCornerProcessor3D ( plint x_, plint y_, plint z_, BlockLattice3D<T,Descriptor>& lattice_ )
    : x(x_), y(y_), z(z_),
      lattice(lattice_)
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerProcessor3D<T, Descriptor, xNormal, yNormal, zNormal>::process() {
    typedef SymmetricTensorImpl<T,Descriptor<T>::d> S;
    Cell<T,Descriptor>& cell = lattice.get(x,y,z);
    Dynamics<T,Descriptor>& dynamics = cell.getDynamics();

    T rho100 = lattice.get(x - 1*xNormal, y - 0*yNormal, z - 0*zNormal).computeDensity();
    T rho010 = lattice.get(x - 0*xNormal, y - 1*yNormal, z - 0*zNormal).computeDensity();
    T rho001 = lattice.get(x - 0*xNormal, y - 0*yNormal, z - 1*zNormal).computeDensity();
    T rho = (T)1/(T)3 * (rho001 + rho010 + rho100);

    Array<T,Descriptor<T>::d> dx_u, dy_u, dz_u;
    fd::DirectedGradients3D<T, Descriptor, 0, xNormal, 0, true>::o1_velocityDerivative(dx_u, lattice, x,y,z);
    fd::DirectedGradients3D<T, Descriptor, 1, yNormal, 0, true>::o1_velocityDerivative(dy_u, lattice, x,y,z);
    fd::DirectedGradients3D<T, Descriptor, 2, zNormal, 0, true>::o1_velocityDerivative(dz_u, lattice, x,y,z);

    T dx_ux = dx_u[0];
    T dy_ux = dy_u[0];
    T dz_ux = dz_u[0];
    T dx_uy = dx_u[1];
    T dy_uy = dy_u[1];
    T dz_uy = dz_u[1];
    T dx_uz = dx_u[2];
    T dy_uz = dy_u[2];
    T dz_uz = dz_u[2];
    T omega = dynamics.getOmega();
    T sToPi = - rho / Descriptor<T>::invCs2 / omega;
    Array<T,SymmetricTensor<T,Descriptor>::n> pi;
    pi[S::xx] = (T)2 * dx_ux * sToPi;
    pi[S::yy] = (T)2 * dy_uy * sToPi;
    pi[S::zz] = (T)2 * dz_uz * sToPi;
    pi[S::xy] = (dx_uy + dy_ux) * sToPi;
    pi[S::xz] = (dx_uz + dz_ux) * sToPi;
    pi[S::yz] = (dy_uz + dz_uy) * sToPi;

    // Computation of the particle distribution functions
    // according to the regularized formula
    Array<T,Descriptor<T>::d> u, j;
    cell.computeVelocity(u);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] = rho*u[iD];
    }
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = dynamics.computeEquilibrium(iPop,Descriptor<T>::rhoBar(rho),j,jSqr) +
                         offEquilibriumTemplates<T,Descriptor>::fromPiToFneq(iPop, pi);
    }
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessor3D<T, Descriptor, xNormal, yNormal, zNormal>*
     OuterVelocityCornerProcessor3D<T, Descriptor, xNormal, yNormal, zNormal>::clone() const
{
    return new OuterVelocityCornerProcessor3D<T, Descriptor, xNormal, yNormal, zNormal>(*this);
}


////////  OuterVelocityCornerProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal>::
    OuterVelocityCornerProcessorGenerator3D(plint x_, plint y_, plint z_)
    : BoxedDataProcessorGenerator3D<T>(Box3D(x_,x_, y_,y_, z_,z_))
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
DataProcessor3D<T>*
    OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal>::
        generate(std::vector<AtomicBlock3D<T>*> objects) const
{
    PLB_PRECONDITION( typeid(*objects[0]) == typeid(BlockLattice3D<T,Descriptor>) );
    return new OuterVelocityCornerProcessor3D<T,Descriptor, xNormal,yNormal,zNormal>
                   ( this->getDomain().x0, this->getDomain().y0, this->getDomain().z0,
                     *dynamic_cast<BlockLattice3D<T,Descriptor>*> (objects[0]) );
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal>::clone() const
{
    return new OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal, yNormal, zNormal>
                   (this->getDomain().x0, this->getDomain().y0, this->getDomain().z0);
}

}  // namespace plb

#endif  // BOUNDARY_PROCESSOR_3D_HH
