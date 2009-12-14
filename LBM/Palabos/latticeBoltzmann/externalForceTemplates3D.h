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

/* Main author: Orestis Malaspinas
 */

/** \file
 * 3D specialization of externalForceTemplates functions
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_3D_H
#define EXTERNAL_FORCE_TEMPLATES_3D_H

namespace plb {

template<typename T>
struct externalForceTemplates<T, descriptors::ForcedD3Q19Descriptor> 
{
static void addGuoForce(
                Cell<T,descriptors::ForcedD3Q19Descriptor>& cell,
                Array<T,descriptors::ForcedD3Q19Descriptor<T>::d> const& u, T omega, T amplitude)
{
    static const int forceBeginsAt
        = descriptors::ForcedD3Q19Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = cell.getExternal(forceBeginsAt);
    T mu = amplitude*((T)1-omega/(T)2);
    
    static const T oneOver6 = (T)1/(T)6;
    static const T oneOver12 = (T)1/(T)12;
    
    cell[0] += -mu*(force[0]*u[0]+force[1]*u[1]+force[2]*u[2]);
    
    cell[1] += oneOver6*mu*(force[0]*(-(T)1+2*u[0])-force[1]*u[1]-force[2]*u[2]);
    
    cell[2] += -oneOver6*mu*(force[0]*u[0]+force[1]*((T)1-2*u[1])+force[2]*u[2]);
    
    cell[3] += -oneOver6*mu*(force[0]*u[0]
                           + force[1]*u[1]
                           + force[2]*((T)1-2*u[2]));
    
    cell[4] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]+3*u[1])
                           + force[1]*(-(T)1+2*u[1]+3*u[0])
                           - force[2]*u[2]);
    
    cell[5] += oneOver12*mu*( force[0]*(-(T)1+2*u[0]-3*u[1])
                            + force[1]*((T)1+2*u[1]-3*u[0])
                            - force[2]*u[2]);
    
    cell[6] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]+3*u[2])
                           - force[1]*u[1]
                           + force[2]*(-(T)1+2*u[2]+3*u[0]));
    
    cell[7] += oneOver12*mu*(force[0]*(-(T)1+2*u[0]-3*u[2])
                           - force[1]*u[1]
                           + force[2]*((T)1+2*u[2]-3*u[0]));
    
    cell[8] += -oneOver12*mu*(force[0]*u[0]
                            + force[1]*((T)1-2*u[1]-3*u[2])
                            + force[2]*((T)1-2*u[2]-3*u[1]));
    
    cell[9] += -oneOver12*mu*(force[0]*u[0]
                            + force[1]*((T)1-2*u[1]+3*u[2])
                            + force[2]*(-(T)1-2*u[2]+3*u[1]));
    
    cell[10] += oneOver6*mu*(force[0]*((T)1+2*u[0])
                            -force[1]*u[1]
                            -force[2]*u[2]);
    
    cell[11] += -oneOver6*mu*(force[0]*u[0]
                            +force[1]*(-(T)1-2*u[1])
                            +force[2]*u[2]);
    
    cell[12] += -oneOver6*mu*(force[0]*u[0]
                             +force[1]*u[1]
                             +force[2]*(-(T)1-2*u[2]));
    
    cell[13] += oneOver12*mu*(force[0]*((T)1+2*u[0]+3*u[1])
                             +force[1]*((T)1+2*u[1]+3*u[0])
                             -force[2]*u[2]);
    
    cell[14] += oneOver12*mu*(force[0]*((T)1+2*u[0]-3*u[1])
                             +force[1]*(-(T)1+2*u[1]-3*u[0])
                             -force[2]*u[2]);
    
    cell[15] += oneOver12*mu*(force[0]*((T)1+2*u[0]+3*u[2])
                             -force[1]*u[1]
                             +force[2]*((T)1+2*u[2]+3*u[0]));
    
    cell[16] += oneOver12*mu*(force[0]*((T)1+2*u[0]-3*u[2])
                             -force[1]*u[1]
                             +force[2]*(-(T)1+2*u[2]-3*u[0]));
    
    cell[17] += -oneOver12*mu*(force[0]*u[0]
                              +force[1]*(-(T)1-2*u[1]-3*u[2])
                              +force[2]*(-(T)1-2*u[2]-3*u[1]));
    
    cell[18] += -oneOver12*mu*(force[0]*u[0]
                              +force[1]*(-(T)1-2*u[1]+3*u[2])
                              +force[2]*((T)1-2*u[2]+3*u[1]));
}

};

}  // namespace plb

#endif
