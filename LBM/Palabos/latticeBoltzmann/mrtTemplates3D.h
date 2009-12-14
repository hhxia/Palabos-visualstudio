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

#ifndef MRT_TEMPLATES_3D_H
#define MRT_TEMPLATES_3D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D2Q9 lattice
template<typename T>
struct mrtTemplates<T, descriptors::MRTD3Q19Descriptor> 
{
    typedef descriptors::MRTD3Q19Descriptor<T> Descriptor; 

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                    T rho, Array<T,3> const& u, T uSqr )
    {
        momentsEq[1] = rho*(-(T)11 + (T)19*uSqr);
        momentsEq[2] = rho*((T)3 - (T)11/(T)2 * uSqr);
        
        momentsEq[4] = -(T)2/(T)3*rho*u[0];
        momentsEq[6] = -(T)2/(T)3*rho*u[1];
        momentsEq[8] = -(T)2/(T)3*rho*u[2];
        
        momentsEq[9] = rho*((T)2*u[0]*u[0] - u[1]*u[1]-u[2]*u[2]);
        momentsEq[10] = -momentsEq[9]/(T)2;
        momentsEq[11] = rho*(u[1]*u[1]-u[2]*u[2]);
        momentsEq[12] = -momentsEq[11]/(T)2;
        
        momentsEq[13] = rho*u[0]*u[1];
        momentsEq[14] = rho*u[1]*u[2];
        momentsEq[15] = rho*u[0]*u[2];
        
    }
    
    /// Computation of all moments (specialized for d3q19)
    static void computeMoments(Array<T,Descriptor::q>& moments, Cell<T,descriptors::MRTD3Q19Descriptor>& cell)
    {
        moments[1] =
                -(T)30*cell[0]-(T)11*cell[1]-(T)11*cell[2]-(T)11*cell[3]
                +(T)8*cell[4]+(T)8*cell[5]+(T)8*cell[6]+(T)8*cell[7]
                +(T)8*cell[8]+(T)8*cell[9]-(T)11*cell[10]-(T)11*cell[11]
                -(T)11*cell[12]+(T)8*cell[13]+(T)8*cell[14]+(T)8*cell[15]
                +(T)8*cell[16]+(T)8*cell[17]+(T)8*cell[18] - (T)11;
        
        moments[2] =
                (T)12*cell[0]-(T)4*cell[1]-(T)4*cell[2]-(T)4*cell[3]
                +cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]
                -(T)4*cell[10]-(T)4*cell[11]-(T)4*cell[12]+cell[13]
                +cell[14]+cell[15]+cell[16]+cell[17]+cell[18] + (T)3;
        
        moments[4] =
                (T)4*cell[1]-cell[4]-cell[5]-cell[6]-cell[7]
                -(T)4*cell[10]+cell[13]+cell[14]+cell[15]+cell[16];
        
        moments[6] =
                (T)4*cell[2]-cell[4]+cell[5]-cell[8]-cell[9]
                -(T)4*cell[11]+cell[13]-cell[14]+cell[17]+cell[18];
        
        moments[8] =
                (T)4*cell[3]-cell[6]+cell[7]-cell[8]+cell[9]
                -(T)4*cell[12]+cell[15]-cell[16]+cell[17]-cell[18];
        
        moments[9] =
                ((T)2*cell[1]-cell[2]-cell[3]+cell[4]+cell[5]
                +cell[6]+cell[7]-(T)2*cell[8]-(T)2*cell[9]+(T)2*cell[10]
                -cell[11]-cell[12]+cell[13]+cell[14]+cell[15]+cell[16]
                -(T)2*cell[17]-(T)2*cell[18]) /*/ (T)3*/; // M.f[9]=3*pxx;
        
        moments[10] =
                -(T)4*cell[1]+(T)2*cell[2]+(T)2*cell[3]+cell[4]
                +cell[5]+cell[6]+cell[7]-(T)2*cell[8]-(T)2*cell[9]
                -(T)4*cell[10]+(T)2*cell[11]+(T)2*cell[12]+cell[13]
                +cell[14]+cell[15]+cell[16]-(T)2*cell[17]-(T)2*cell[18];
        
        moments[11] =
                (cell[2]-cell[3]+cell[4]+cell[5]-cell[6]-cell[7]
                +cell[11]-cell[12]+cell[13]+cell[14]-cell[15]-cell[16])/*/(T)3*/;// M.f[9]=3*pxx;
        
        moments[12] =
                -(T)2*cell[2]+(T)2*cell[3]+cell[4]+cell[5]
                -cell[6]-cell[7]-(T)2*cell[11]+(T)2*cell[12]
                +cell[13]+cell[14]-cell[15]-cell[16];
        
        moments[13] =
                cell[4]-cell[5]+cell[13]-cell[14];
        
        moments[14] =
                cell[8]-cell[9]+cell[17]-cell[18];
        
        moments[15] =
                cell[6]-cell[7]+cell[15]-cell[16];
        
        moments[16] =
                -cell[4]-cell[5]+cell[6]+cell[7]+cell[13]+cell[14]-cell[15]-cell[16];
        
        moments[17] =
                cell[4]-cell[5]-cell[8]-cell[9]-cell[13]+cell[14]+cell[17]+cell[18];
        
        moments[18] =
                -cell[6]+cell[7]+cell[8]-cell[9]+cell[15]-cell[16]-cell[17]+cell[18];
        
    }
    
    /// MRT collision step
    static T mrtCollision( Cell<T,descriptors::MRTD3Q19Descriptor>& cell,
                           T rho, Array<T,3> const& u, T invM_S[19][19] )
    {
        typedef descriptors::MRTD3Q19Descriptor<T> L;
        T uSqr = VectorTemplateImpl<T,3>::normSqr(u);
        Array<T,19> moments, momentsEq;

        computeMoments(moments,cell);
        computeEquilibrium(momentsEq,rho,u,uSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom4 = moments[4] - momentsEq[4];
        T mom6 = moments[6] - momentsEq[6];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];

        
        cell[0] -= invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        cell[1] -= invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        cell[2] -= invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        cell[3] -= invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        cell[4] -= invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][4]*mom4
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13
                +invM_S[4][16]*mom16
                +invM_S[4][17]*mom17;
        
        cell[5] -= invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][4]*mom4
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13
                +invM_S[5][16]*mom16
                +invM_S[5][17]*mom17;
        
        cell[6] -= invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][4]*mom4
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15
                +invM_S[6][16]*mom16
                +invM_S[6][18]*mom18;
        
        cell[7] -= invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][4]*mom4
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15
                +invM_S[7][16]*mom16
                +invM_S[7][18]*mom18;
        
        cell[8] -= invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][6]*mom6
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14
                +invM_S[8][17]*mom17
                +invM_S[8][18]*mom18;
        
        cell[9] -= invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][6]*mom6
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14
                +invM_S[9][17]*mom17
                +invM_S[9][18]*mom18;
        
        cell[10] -= invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        cell[11] -= invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        cell[12] -= invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        cell[13] -= invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][4]*mom4
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13
                +invM_S[13][16]*mom16
                +invM_S[13][17]*mom17;
        
        cell[14] -= invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][4]*mom4
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13
                +invM_S[14][16]*mom16
                +invM_S[14][17]*mom17;
        
        cell[15] -= invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][4]*mom4
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15
                +invM_S[15][16]*mom16
                +invM_S[15][18]*mom18;
        
        cell[16] -= invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][4]*mom4
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15
                +invM_S[16][16]*mom16
                +invM_S[16][18]*mom18;
        
        cell[17] -= invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][6]*mom6
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14
                +invM_S[17][17]*mom17
                +invM_S[17][18]*mom18;
        
        cell[18] -= invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][6]*mom6
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14
                +invM_S[18][17]*mom17
                +invM_S[18][18]*mom18;
        
        return uSqr;
    }

};


}  // namespace plb

#endif
