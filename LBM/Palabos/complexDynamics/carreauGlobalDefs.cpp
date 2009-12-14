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
 * Global Definitions -- implementation file.
 */

#include "complexDynamics/carreauGlobalDefs.h"

namespace plb {

namespace global {

void CarreauParametersClass::setNu0(double nu0_)
{
    nu0 = nu0_;
}

void CarreauParametersClass::setNuInf(double nuInf_)
{
	nuInf = nuInf_;
}

void CarreauParametersClass::setLambda(double lambda_)
{
    lambda = lambda_;
}


void CarreauParametersClass::setExponent(double n_)
{
    n = n_;
}


double CarreauParametersClass::getNu0() const
{
    return nu0;
}

double CarreauParametersClass::getNuInf() const
{
	return nuInf;
}


double CarreauParametersClass::getLambda() const
{
    return lambda;
}


double CarreauParametersClass::getExponent() const
{
    return n;
}

}  // namespace global

}  // namespace plb
