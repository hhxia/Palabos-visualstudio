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

#ifndef ZOU_HE_BOUNDARY_3D_H
#define ZOU_HE_BOUNDARY_3D_H

#include "core/globalDefs.h"
#include "boundaryCondition/boundaryCondition3D.h"

namespace plb {

////////// Factory function for Zou/He BC ///////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createZouHeBoundaryCondition3D();

}


#endif
