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

#ifndef FD_STENCILS_1D_H
#define FD_STENCILS_1D_H

#include "core/globalDefs.h"

namespace plb {

namespace fd {

    /// Second-order central gradient (u_p1 = u(x+1))
    template<typename T>
    T ctl_diff(T u_p1, T u_m1) {
        return (u_p1 - u_m1) / (T)2;
    }

    /// Second-order forward gradient (u_1 = u(x+1))
    template<typename T>
    T fwd_diff(T u_0, T u_1, T u_2) {
        return (-(T)3*u_0 + (T)4*u_1 - (T)1*u_2) / (T)2;
    }

    /// Second-order backward gradient (u_1 = u(x-1))
    template<typename T>
    T bwd_diff(T u_0, T u_1, T u_2) {
        return -fwd_diff(u_0, u_1, u_2);
    }

    /// First-order forward gradient (u_1 = u(x+1))
    template<typename T>
    T o1_fwd_diff(T u_0, T u_1) {
        return (-u_0 + u_1);
    }

    /// Value at u_0 for which asymmetric gradient is zero (u_1 = u(x+1))
    template<typename T>
    T boundaryZeroGradient(T u_1, T u_2) {
        return (T)4/(T)3*u_1 - (T)1/(T)3*u_2;
    }

    /// Linear interpolation (yields u0 at pos=0 and u1 at pos=1)
    template<typename T>
    T linearInterpolate(T u_0, T u_1, T pos) {
        return ((T)1-pos)*u_0 + pos*u_1;
    }

}  // namespace fd

}  // namespace plb


#endif
