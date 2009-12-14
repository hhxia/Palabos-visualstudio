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

/** \file
 * Groups all the generic 3D template files in the boundaryConditions directory.
 */
#include "boundaryCondition/boundaryDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics3D.hh"
#include "boundaryCondition/equilibriumBoundaryDynamics.hh"
#include "boundaryCondition/inamuroAnalyticalDynamics.hh"
#include "boundaryCondition/zouHeBoundary3D.hh"
#include "boundaryCondition/zouHeDynamics.hh"
#include "boundaryCondition/boundaryCondition3D.hh"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.hh"
#include "boundaryCondition/neumannCondition3D.hh"
#include "boundaryCondition/bounceBackModels.hh"
#include "boundaryCondition/bounceBackModels3D.hh"
