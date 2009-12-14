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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef MRT_TEMPLATES_H
#define MRT_TEMPLATES_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/mrtLattices.h"

namespace plb {

/// All helper functions are inside this structure
template<typename T, template<typename U> class Descriptor>
struct mrtTemplates {

    /// Computation of equilibrium distribution (in moments space)
    static T equilibrium( plint iPop, T rho, Array<T,Descriptor<T>::d> const& u,
                          const T uSqr )
    {
        T equ = T();
        for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop)
        {
            equ += Descriptor<T>::M[iPop][jPop] * 
                    (dynamicsTemplates<T,Descriptor>::equilibrium(jPop,rho,u,uSqr) + 
                     Descriptor<T>::t[jPop]);
        }
        
        return equ;
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor<T>::q>& momentsEq, 
                                    T rho, Array<T,Descriptor<T>::d> const& u,
                                    const T uSqr )
    {
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        {
            momentsEq[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop)
            {
                momentsEq[iPop] += Descriptor<T>::M[iPop][jPop] * 
                        (dynamicsTemplates<T,Descriptor>::equilibrium(jPop,rho,u,uSqr) + 
                        Descriptor<T>::t[jPop]);
            }
        }
    }
    
    static void computeMoments(T moments[Descriptor<T>::q], Cell<T,Descriptor> &cell)
    {
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop)
        {
            moments[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop)
            {
                moments[iPop] += Descriptor<T>::M[iPop][jPop] * 
                        (cell[jPop] + Descriptor<T>::t[jPop]);
            }
        }
    }
    
    /// MRT collision step
    static T mrtCollision( Cell<T,Descriptor>& cell,
                           T rho, Array<T,Descriptor<T>::d> const& u,
                           T invM_S[Descriptor<T>::q][Descriptor<T>::q])
    {
        T uSqr = VectorTemplate<T,Descriptor>::normSqr(u);
        Array<T,Descriptor<T>::q>& momentsEq;
        Array<T,Descriptor<T>::q>& moments;
        
        computeMoments(moments,cell);
        computeEquilibrium(momentsEq,rho,u,uSqr);
    
        for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) 
        {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor<T>::q; ++jPop)
            {
                collisionTerm += invM_S[iPop][jPop] * 
                        (moments[jPop] - momentsEq[jPop]);
            }
            cell[iPop] -= collisionTerm;
        }
        
        return uSqr;
    }

};  // struct mrtHelpers

}  // namespace plb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "latticeBoltzmann/mrtTemplates2D.h"
#include "latticeBoltzmann/mrtTemplates3D.h"

#endif
