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
 * Utilities for 2D multi data distributions -- implementation.
 */

#include "multiBlock/staticRepartitions2D.h"
#include "atomicBlock/dataField2D.hh"
#include "core/blockStatistics.hh"
#include "algorithm/basicAlgorithms.h"

namespace plb {

////////////////////// function createRegularMultiBlockDistribution2D /////////////////////

MultiBlockDistribution2D createRegularMultiBlockDistribution2D (
        plint nx, plint ny, plint numBlocksX, plint numBlocksY,
        plint envelopeWidth )
{
    MultiBlockDistribution2D dataGeometry(nx, ny);
    plint posX = 0;
    plint iBlock = 0;
    for (plint iBlockX=0; iBlockX<numBlocksX; ++iBlockX) {
        plint lx = nx / numBlocksX;
        if (iBlockX < nx%numBlocksX) ++lx;
        plint posY = 0;
        for (plint iBlockY=0; iBlockY<numBlocksY; ++iBlockY) {
             plint ly = ny / numBlocksY;
             if (iBlockY < ny%numBlocksY) ++ly;
             dataGeometry.addBlock(Box2D(posX, posX+lx-1, posY, posY+ly-1), envelopeWidth, iBlock);
             iBlock++;
             posY += ly;
        }
        posX += lx;
    }
    return dataGeometry;
}

MultiBlockDistribution2D createRegularMultiBlockDistribution2D (
        plint nx, plint ny, plint envelopeWidth, int numProc) {
    std::vector<int> repartition = algorithm::evenRepartition(numProc, 2);
    return createRegularMultiBlockDistribution2D(nx, ny, repartition[0], repartition[1], envelopeWidth);
}

MultiBlockDistribution2D createXSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth )
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    
    std::vector<int> numActivePerSlice;
    plint numActiveTotal = 0;
    for(plint iX=0; iX<nX; iX++) {
        plint numActiveCurrentSlice = 0;
        for (plint iY=0; iY<nY; iY++) {
            if (cellTypeField.get(iX,iY) > 0) numActiveCurrentSlice++;
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;
    
    MultiBlockDistribution2D dataGeometry(nX, nY);
    
    plint iX=0;
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        plint posX = iX;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock<numActivePerBlock && iX<nX) {
            numActiveCurrentBlock += numActivePerSlice[iX];
            iX++;
        }
        dataGeometry.addBlock(Box2D(posX, iX-1, 0, nY-1), envelopeWidth, iBlock);
    }
    return dataGeometry;
}

MultiBlockDistribution2D createYSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth )
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    
    std::vector<int> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iY=0; iY<nY; iY++) {
        plint numActiveCurrentSlice = 0;
        for(plint iX=0; iX<nX; iX++) {
            if (cellTypeField.get(iX,iY) > 0) numActiveCurrentSlice++;
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;
    
    MultiBlockDistribution2D dataGeometry(nX, nY);
    
    plint iY=0;
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        plint posY = iY;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock<numActivePerBlock && iY<nY) {
            numActiveCurrentBlock += numActivePerSlice[iY];
            iY++;
        }
        dataGeometry.addBlock(Box2D(0, nX-1, posY, iY-1), envelopeWidth, iBlock);
    }
    return dataGeometry;
}

MultiBlockDistribution2D createXSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField, plint envelopeWidth)
{
    return createXSlicedMultiBlockDistribution2D(cellTypeField, global::mpi().getSize(), envelopeWidth);
}

MultiBlockDistribution2D createYSlicedMultiBlockDistribution2D (
        CellTypeField2D const& cellTypeField, plint envelopeWidth)
{
    return createYSlicedMultiBlockDistribution2D(cellTypeField, global::mpi().getSize(), envelopeWidth);
}

}  // namespace plb
