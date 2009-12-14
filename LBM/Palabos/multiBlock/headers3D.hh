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
 * Groups all the template instantiations for the 3D multiBlock.
 */

#include "multiBlock/multiBlock3D.hh"
#include "multiBlock/multiBlockLattice3D.hh"
#include "multiBlock/multiDataFieldHandler3D.hh"
#include "multiBlock/multiDataField3D.hh"
#include "multiBlock/serialMultiBlockLattice3D.hh"
#include "multiBlock/serialMultiDataField3D.hh"
#include "multiBlock/serialBlockCommunicator3D.hh"
#include "multiBlock/multiDataCouplingWrapper3D.hh"
#include "multiBlock/reductiveMultiDataCouplingWrapper3D.hh"
#include "multiBlock/multiBlockOperations3D.hh"
#include "multiBlock/multiDataAnalysis3D.hh"
#include "multiBlock/multiBlockSerializer3D.hh"
#include "multiBlock/multiLatticeInitializer3D.hh"
#include "multiBlock/combinedStatistics.hh"
