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
 * Utilities for 3D multi data distributions -- header file.
 */

#ifndef STATIC_REPARTITIONS_3D_H
#define STATIC_REPARTITIONS_3D_H

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

/// A 3D field of scalar values used to indicate the type of the cells.
/// Any positive value indicates an active (bulk, boundary) cell, 
/// while zero indicates a non-active (no-dynamics) cell
typedef ScalarField3D<unsigned char> CellTypeField3D;

/// Create a nx-by-ny-by-nz data distribution
MultiBlockDistribution3D createRegularMultiBlockDistribution3D (
        plint nx, plint ny, plint nz, plint numBlocksX, plint numBlocksY, plint numBlocksZ,
        plint envelopeWidth );

/// Create a data distribution with regular blocks, as evenly distributed as possible
MultiBlockDistribution3D createRegularMultiBlockDistribution3D(plint nx, plint ny, plint nz, plint envelopeWidth,
                                                               int numProc = global::mpi().getSize());

/// Create a data distribution by slicing the domain (a block of nX*nY*nZ cells
/// as defined by cellTypeField) into numBlocks blocks along the x-direction. 
/// The x-extent of the blocks is chosen such as to obtain an approximately 
/// equal number of active cells in each block.
MultiBlockDistribution3D createXSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth );

/// cf above.
MultiBlockDistribution3D createYSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth );
        
/// cf above.
MultiBlockDistribution3D createZSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth );

/// Create x-sliced data distribution, balancing the number of active cells between blocks,
/// implicitly setting numBlocks = #processors
MultiBlockDistribution3D createXSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth=1);

/// cf above
MultiBlockDistribution3D createYSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth=1);

/// cf above
MultiBlockDistribution3D createZSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth=1);

}  // namespace plb


#endif
