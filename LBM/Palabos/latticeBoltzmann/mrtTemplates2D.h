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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * Template specializations for some computationally intensive LB 
 * functions of the header file mrtTemplates.h, for the D2Q9 grid.
 */

#ifndef MRT_TEMPLATES_2D_H
#define MRT_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D2Q9 lattice
template<typename T>
struct mrtTemplates<T, descriptors::MRTD2Q9Descriptor> 
{
    typedef descriptors::MRTD2Q9Descriptor<T> Descriptor; 

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                    T rho, Array<T,2> const& u, T uSqr )
    {
//         momentsEq[0] = T();
        momentsEq[1] = -(T)2*rho + (T)3*uSqr;
        momentsEq[2] = rho - (T)3*uSqr;
//         momentsEq[3] = T();
        momentsEq[4] = -u[0];
//         momentsEq[5] = T();
        momentsEq[6] = -u[1];
        momentsEq[7] = u[0]*u[0] - u[1]*u[1];
        momentsEq[8] = u[0]*u[1];
    }
    
    /// Computation of all moments (specialized for d2q9)
    static void computeMoments(Array<T,Descriptor::q>& moments, Cell<T,descriptors::MRTD2Q9Descriptor>& cell)
    {
//         moments[0] = cell[0] + cell[1] + cell[2] + cell[3] + 
//                 cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + (T)1;

        moments[1] = -(T)4*cell[0] +(T)2*cell[1] - cell[2] + (T)2*cell[3] - cell[4] + 
                (T)2*cell[5] - cell[6] + (T)2*cell[7] - cell[8] - (T)2;

        moments[2] = (T)4*cell[0] + cell[1] - (T)2*cell[2] + cell[3] - (T)2*cell[4] + 
                cell[5] - (T)2*cell[6] + cell[7] - (T)2*cell[8] + (T)1;

//         moments[3] = - cell[1] - cell[2] - cell[3] +
//                 cell[5] + cell[6] + cell[7];

        moments[4] = - cell[1] + (T)2*cell[2] - cell[3] +
                cell[5] - (T)2*cell[6] + cell[7];

//         moments[5] = cell[1] - cell[3] - cell[4] - 
//                 cell[5] + cell[7] + cell[8];

        moments[6] = cell[1] - cell[3] + (T)2*cell[4] - 
                cell[5] + cell[7] - (T)2*cell[8];

        moments[7] = cell[2] - cell[4] + cell[6] - cell[8];

        moments[8] = - cell[1] + cell[3] - cell[5] + cell[7];
    }
    
    /// MRT collision step
    static T mrtCollision( Cell<T,descriptors::MRTD2Q9Descriptor>& cell,
                           T rho, Array<T,2> const& u,
                           T invM_S[9][9] )
    {
        typedef descriptors::MRTD2Q9Descriptor<T> L;
        T uSqr = VectorTemplateImpl<T,2>::normSqr(u);
        Array<T,9> moments, momentsEq;

        computeMoments(moments,cell);
        computeEquilibrium(momentsEq,rho,u,uSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom4 = moments[4] - momentsEq[4];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        
        cell[0] -= invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        cell[1] -= invM_S[1][1]*mom1 +
                   invM_S[1][2]*mom2 +
                   invM_S[1][4]*mom4 + 
                   invM_S[1][6]*mom6 +
                   invM_S[1][8]*mom8;
        
        cell[2] -= invM_S[2][1]*mom1 +
                   invM_S[2][2]*mom2 +
                   invM_S[2][4]*mom4 +
                   invM_S[2][7]*mom7;
        
        cell[3] -= invM_S[3][1]*mom1 +
                   invM_S[3][2]*mom2 +
                   invM_S[3][4]*mom4 +
                   invM_S[3][6]*mom6 +
                   invM_S[3][8]*mom8;
        
        cell[4] -= invM_S[4][1]*mom1 +
                   invM_S[4][2]*mom2 +
                   invM_S[4][6]*mom6 +
                   invM_S[4][7]*mom7;
        
        cell[5] -= invM_S[5][1]*mom1 +
                   invM_S[5][2]*mom2 +
                   invM_S[5][4]*mom4 +
                   invM_S[5][6]*mom6 +
                   invM_S[5][8]*mom8;
        
        cell[6] -= invM_S[6][1]*mom1 +
                   invM_S[6][2]*mom2 +
                   invM_S[6][4]*mom4 +
                   invM_S[6][7]*mom7;
        
        cell[7] -= invM_S[7][1]*mom1 +
                   invM_S[7][2]*mom2 +
                   invM_S[7][4]*mom4 +
                   invM_S[7][6]*mom6 +
                   invM_S[7][8]*mom8;
        
        cell[8] -= invM_S[8][1]*mom1 +
                   invM_S[8][2]*mom2 +
                   invM_S[8][6]*mom6 +
                   invM_S[8][7]*mom7;

        return uSqr;
    }

};


}  // namespace plb

#endif
