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
 * 2D specialization of externalForceTemplates functions
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_2D_H
#define EXTERNAL_FORCE_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {
    
template<typename T>
struct externalForceTemplates<T, descriptors::ForcedD2Q9Descriptor> 
{
static void addGuoForce(
                Cell<T,descriptors::ForcedD2Q9Descriptor>& cell,
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> const& u, T omega, T amplitude)
{
    static const int forceBeginsAt
        = descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = cell.getExternal(forceBeginsAt);
    T mu = amplitude*((T)1-omega/(T)2);
    
    static const T oneOver3 = (T)1/(T)3;
    static const T oneOver12 = (T)1/(T)12;
    static const T fourOver3 = (T)4/(T)3;
    
    const T twoUx   = (T)2*u[0];
    const T threeUx = (T)3*u[0];
    
    const T twoUy   = (T)2*u[1];
    const T threeUy = (T)3*u[1];

    cell[0] += -fourOver3*mu*(force[0]*u[0]+force[1]*u[1]);
    
    cell[1] += oneOver12*mu*(force[0]*(-(T)1+twoUx-threeUy)+force[1]*((T)1+twoUy-threeUx));
    
    cell[2] += oneOver3*mu*(force[0]*(-(T)1+twoUx)-force[1]*u[1]);
    
    cell[3] += oneOver12*mu*(force[0]*(-(T)1+twoUx+threeUy)+force[1]*(-(T)1+twoUy+threeUx));
    
    cell[4] += -oneOver3*mu*(force[0]*u[0]+force[1]*((T)1-twoUy));
    
    cell[5] += oneOver12*mu*(force[0]*((T)1+twoUx-threeUy)+force[1]*(-(T)1+twoUy-threeUx));
    
    cell[6] += oneOver3*mu*(force[0]*((T)1+twoUx)-force[1]*u[1]);
    
    cell[7] += oneOver12*mu*(force[0]*((T)1+twoUx+threeUy)+force[1]*((T)1+twoUy+threeUx));
    
    cell[8] += -oneOver3*mu*(force[0]*u[0]+force[1]*(-(T)1-twoUy));
}

};

}  // namespace plb

#endif
