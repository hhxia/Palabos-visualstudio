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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_HELPERS_H
#define ENTROPIC_LB_HELPERS_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct entropicLbTemplates 
{
    /// Computation of equilibrium distribution
    static T equilibrium( plint iPop, T rho, Array<T,Descriptor<T>::d> const& u)
    {
        typedef Descriptor<T> L;
        const T invCs = sqrt(L::invCs2);
        const T sqt3 = sqrt(3.0);
        T prod = (T)1;
        for (int iD=0; iD < Descriptor<T>::d; ++iD)
        {
            T uc = u[iD] * invCs; // u[iD] / c_s

            prod *= ((T)2 - sqrt(1.0+uc*uc)) * 
                    pow((2.0 / sqt3 * uc + 
                    sqrt(1.0+uc*uc))/(1.0-uc/sqt3),
                        L::c[iPop][iD]/sqt3*invCs);
        }
        return rho*L::t[iPop]*prod-L::SkordosFactor()*L::t[iPop];
    }
    
    /// Computation of equilibrium distribution
    static T equilibriumApprox( plint iPop, T rho, Array<T,Descriptor<T>::d> const& u)
    {
        typedef Descriptor<T> L;

        T uSqr = VectorTemplate<T,Descriptor>::normSqr(u);
        T cu = T();
        for (int iD=0; iD < Descriptor<T>::d; ++iD)
        {
            cu += L::c[iPop][iD]*u[iD];
        }
        
        return rho * L::t[iPop] * (1.0 +
                cu*L::invCs2 - 0.5 * uSqr*L::invCs2 + 0.5*pow(L::invCs2,2)*cu*cu
                - 0.5*pow(L::invCs2,2)*cu*uSqr + pow(cu,3)*pow(L::invCs2,3)/6.0
                + 0.125*uSqr*uSqr*pow(L::invCs2,2) - 0.25*cu*cu*uSqr*pow(L::invCs2,3)
                + pow(cu,4)*pow(L::invCs2,4)/24.0)-L::SkordosFactor()*L::t[iPop];
    }
};

}

#include "latticeBoltzmann/entropicLbTemplates2D.h"
#include "latticeBoltzmann/entropicLbTemplates3D.h"

#endif
