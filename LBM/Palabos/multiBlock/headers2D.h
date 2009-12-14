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
 * Groups all the include files for the 2D multiBlock.
 */

#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataFieldHandler2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/serialMultiBlockLattice2D.h"
#include "multiBlock/serialMultiDataField2D.h"
#include "multiBlock/serialBlockCommunicator2D.h"
#include "multiBlock/staticRepartitions2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/multiDataCouplingWrapper2D.h"
#include "multiBlock/reductiveMultiDataCouplingWrapper2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiDataAnalysis2D.h"
#include "multiBlock/multiLatticeInitializer2D.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include "multiBlock/combinedStatistics.h"
