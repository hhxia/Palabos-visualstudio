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
 * #include "core/globalDefs.h"
 * Groups all the include files for the 3D multiBlock.
 */

#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataFieldHandler3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/serialMultiBlockLattice3D.h"
#include "multiBlock/serialMultiDataField3D.h"
#include "multiBlock/serialBlockCommunicator3D.h"
#include "multiBlock/staticRepartitions3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiDataCouplingWrapper3D.h"
#include "multiBlock/reductiveMultiDataCouplingWrapper3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiDataAnalysis3D.h"
#include "multiBlock/multiLatticeInitializer3D.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlockInfo3D.h"

