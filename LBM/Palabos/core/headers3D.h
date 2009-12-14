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
 * Groups all the include files for basic 3D dynamics.
 */
#include "core/globalDefs.h"
#include "core/plbInit.h"
#include "core/geometry3D.h"
#include "core/identifiers.h"
#include "core/units.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "core/blockStatistics.h"
#include "core/dataFieldBase3D.h"
#include "core/serializer.h"
#include "core/block3D.h"
#include "core/blockLatticeBase3D.h"
#include "core/latticeStatistics.h"
#include "core/dataAnalysis3D.h"
#include "core/dataAnalysisFunctionals3D.h"
#include "core/plbTimer.h"
#include "core/periodicity3D.h"
