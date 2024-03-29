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
 * Helper functions for the implementation of non-Newtonian Carreau
 * dynamics. 
 */
#ifndef CARREAU_DYNAMICS_TEMPLATES_H
#define CARREAU_DYNAMICS_TEMPLATES_H

#include "core/util.h"

namespace plb {

template<typename T,int N>
struct ImplicitOmega
{
	T operator()(T alpha, T nu0_nuInfOverCs2, T nuInfOverCs2, T nMinusOneOverTwo, T omega0)
    {
        ImplicitOmega<T,N-1> omega;
		T omSqr = omega(alpha,nu0_nuInfOverCs2,nuInfOverCs2,nMinusOneOverTwo,omega0);
        omSqr *= omSqr;
        
		return (T)2/((T)1+(T)2*nu0_nuInfOverCs2*
                   pow((T)1+alpha*omSqr,nMinusOneOverTwo)+nuInfOverCs2);
    }
};

template<typename T>
struct ImplicitOmega<T,0>
{
	T operator()(T alpha, T nu0_nuInfOverCs2, T nuInfOverCs2, T nMinusOneOverTwo, T omega0)
    {
		return (T)2/((T)1+(T)2*nu0_nuInfOverCs2*pow((T)1+
				alpha*omega0*omega0,nMinusOneOverTwo)+nuInfOverCs2);
    }
};

/// This structure forwards the calls to the appropriate helper class
template<typename T, int N>
struct carreauDynamicsTemplates {
	static T fromPiAndRhoToOmega(T alpha, T nu0_nuInfOverCs2, T nuInfOverCs2, T nMinusOneOverTwo, T omega0) 
{
    ImplicitOmega<T,N> omega;
    
	return omega(alpha,nu0_nuInfOverCs2,nuInfOverCs2,nMinusOneOverTwo,omega0);
}

};  // struct carreauDynamicsTemplates

}  // namespace plb

#endif  // CARREAU_DYNAMICS_TEMPLATES_H
