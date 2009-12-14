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
 * Utilities for 2D multi data distributions -- header file.
 */

#ifndef STATIC_REPARTITIONS_2D_H
#define STATIC_REPARTITIONS_2D_H

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

/// A 2D field of scalar values used to indicate the type of the cells.
/// Any positive value indicates an active (bulk, boundary) cell, 
/// while zero indicates a non-active (no-dynamics) cell
typedef ScalarField2D<unsigned char> CellTypeField2D;

/// Create a nx-by-ny data distribution
MultiBlockDistribution2D createRegularMultiBlockDistribution2D (
        plint nx, plint ny, plint numBlocksX, plint numBlocksY,
        plint envelopeWidth );

/// Create a data distribution with regular blocks, as evenly distributed as possible
MultiBlockDistribution2D createRegularMultiBlockDistribution2D(plint nx, plint ny, plint envelopeWidth,
                                                               int numProc = global::mpi().getSize());

/// Create a data distribution by slicing the domain (a block of nX*nY cells
/// as defined by cellTypeField) into numBlocks blocks along the x-direction. 
/// The x-extent of the blocks is chosen such as to obtain an approximately 
/// equal number of active cells in each block.
MultiBlockDistribution2D createXSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth );

/// cf above.
MultiBlockDistribution2D createYSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth );

/// Create x-sliced data distribution, balancing the number of active cells between blocks,
/// implicitly setting numBlocks = #processors
MultiBlockDistribution2D createXSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField, plint envelopeWidth=1);

/// cf above
MultiBlockDistribution2D createYSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField, plint envelopeWidth=1);

}  // namespace plb


#endif
