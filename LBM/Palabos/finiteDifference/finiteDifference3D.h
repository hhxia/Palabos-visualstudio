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

#ifndef FINITE_DIFFERENCE_3D_H
#define FINITE_DIFFERENCE_3D_H

#include "core/globalDefs.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

namespace fd {

    template<typename T, template<typename U> class Descriptor,
             int direction, int orientation, int deriveDirection,
             bool orthogonal>
    struct DirectedGradients3D {
        /// Nearest-neighbor evaluation of velocity derivative, first-order accurate
        ///   only along the boundary normal.
        static void o1_velocityDerivative( Array<T,Descriptor<T>::d>& velDeriv,
                                           BlockLattice3D<T,Descriptor> const& blockLattice,
                                           plint iX, plint iY, plint iZ );
        /// Nearest-neighbor evaluation of density derivative, first-order accurate
        ///   only along the boundary normal.
        static void o1_densityDerivative( T& rhoDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ );
        /// Next-to-nearest-neibhbor, second-order accurate evaluation of velocity
        ///   derivative.
        static void o2_velocityDerivative( Array<T,Descriptor<T>::d>& velDeriv,
                                           BlockLattice3D<T,Descriptor> const& blockLattice,
                                           plint iX, plint iY, plint iZ );
        /// Next-to-nearest-neibhbor, second-order accurate evaluation of density
        ///   derivative.
        static void o2_densityDerivative( T& rhoDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ );
    };

    // Implementation for orthogonal==true; i.e. the derivative is along the boundary normal.
    template<typename T, template<typename U> class Descriptor,
             int direction, int orientation, int deriveDirection>
    struct DirectedGradients3D<T, Descriptor, direction, orientation,
                               deriveDirection, true>
    {
        static void o1_velocityDerivative(Array<T,Descriptor<T>::d>& velDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ)
        {
            Array<T,Descriptor<T>::d> u0, u1;
            
            blockLattice.get(iX,iY,iZ).computeVelocity(u0);
            blockLattice.get (
                iX+(direction==0 ? (-orientation):0),
                iY+(direction==1 ? (-orientation):0),
                iZ+(direction==2 ? (-orientation):0)  ).computeVelocity(u1);

            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                velDeriv[iD] = -orientation * fd::o1_fwd_diff(u0[iD], u1[iD]);
            }
        }

        static void o1_densityDerivative(T& rhoDeriv,
                                         BlockLattice3D<T,Descriptor> const& blockLattice,
                                         plint iX, plint iY, plint iZ)
        {
            // note that the derivative runs along direction.
            T rho0 = blockLattice.get(iX,iY,iZ).computeRho();
            T rho1 = blockLattice.get (
                        iX+(direction==0 ? (-orientation):0),
                        iY+(direction==1 ? (-orientation):0),
                        iZ+(direction==2 ? (-orientation):0)  ).computeRho();
            rhoDeriv = -orientation * fd::o1_fwd_diff(rho0, rho1);
        }

        static void o2_velocityDerivative(Array<T,Descriptor<T>::d>& velDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ)
        {
            Array<T,Descriptor<T>::d> u0, u1, u2;
            
            blockLattice.get(iX,iY,iZ).computeVelocity(u0);
            blockLattice.get (
                iX+(direction==0 ? (-orientation):0),
                iY+(direction==1 ? (-orientation):0),
                iZ+(direction==2 ? (-orientation):0)  ).computeVelocity(u1);
            blockLattice.get (
                iX+(direction==0 ? (-2*orientation):0),
                iY+(direction==1 ? (-2*orientation):0),
                iZ+(direction==2 ? (-2*orientation):0) ).computeVelocity(u2);

            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                velDeriv[iD] = -orientation * fd::fwd_diff(u0[iD], u1[iD], u2[iD]);
            }
        }

        static void o2_densityDerivative(T& rhoDeriv,
                                         BlockLattice3D<T,Descriptor> const& blockLattice,
                                         plint iX, plint iY, plint iZ)
        {
            // note that the derivative runs along direction.
            T rho0 = blockLattice.get(iX,iY,iZ).computeRho();
            T rho1 = blockLattice.get (
                        iX+(direction==0 ? (-orientation):0),
                        iY+(direction==1 ? (-orientation):0),
                        iZ+(direction==2 ? (-orientation):0)  ).computeRho();
            T rho2 = blockLattice.get (
                        iX+(direction==0 ? (-2*orientation):0),
                        iY+(direction==1 ? (-2*orientation):0),
                        iZ+(direction==2 ? (-2*orientation):0) ).computeRho();

            rhoDeriv = -orientation * fd::fwd_diff(rho0, rho1, rho2);
        }
    };

    // Implementation for orthogonal==false; i.e. the derivative is aligned with the boundary.
    template<typename T, template<typename U> class Descriptor,
             int direction, int orientation, int deriveDirection>
    struct DirectedGradients3D<T, Descriptor, direction, orientation,
                               deriveDirection, false>
    {
        static void  o1_velocityDerivative(Array<T,Descriptor<T>::d>& velDeriv,
                                           BlockLattice3D<T,Descriptor> const& blockLattice,
                                           plint iX, plint iY, plint iZ)
        {
            // Along the boundary, second-order accuracy is achieved with a nearest-
            //   neighbor scheme.
            o2_velocityDerivative(velDeriv, blockLattice, iX,iY,iZ);
        }
        static void  o1_densityDerivative(T& rhoDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ)
        {
            // Along the boundary, second-order accuracy is achieved with a nearest-
            //   neighbor scheme.
            o2_densityDerivative(rhoDeriv, blockLattice, iX,iY,iZ);
        }
        static void  o2_velocityDerivative(Array<T,Descriptor<T>::d>& velDeriv,
                                           BlockLattice3D<T,Descriptor> const& blockLattice,
                                           plint iX, plint iY, plint iZ)
        {
            Array<T,Descriptor<T>::d> u_p1, u_m1;
            
            blockLattice.get (
                iX+(deriveDirection==0 ? 1:0),
                iY+(deriveDirection==1 ? 1:0),
                iZ+(deriveDirection==2 ? 1:0) ).computeVelocity(u_p1);

            blockLattice.get (
                iX+(deriveDirection==0 ? (-1):0),
                iY+(deriveDirection==1 ? (-1):0),
                iZ+(deriveDirection==2 ? (-1):0) ).computeVelocity(u_m1);

            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                velDeriv[iD] = fd::ctl_diff(u_p1[iD],u_m1[iD]);
            }
        }

        static void  o2_densityDerivative(T& rhoDeriv,
                                          BlockLattice3D<T,Descriptor> const& blockLattice,
                                          plint iX, plint iY, plint iZ)
        {
            typedef Descriptor<T> L;

            T rho_p1 = blockLattice.get (
                        iX+(deriveDirection==0 ? 1:0),
                        iY+(deriveDirection==1 ? 1:0),
                        iZ+(deriveDirection==2 ? 1:0) ).computeRho();

            T rho_m1 = blockLattice.get (
                        iX+(deriveDirection==0 ? (-1):0),
                        iY+(deriveDirection==1 ? (-1):0),
                        iZ+(deriveDirection==2 ? (-1):0) ).computeRho();

            rhoDeriv = fd::ctl_diff(rho_p1, rho_m1);
        }
    };

}  // namespace fd

}  // namespace plb


#endif
