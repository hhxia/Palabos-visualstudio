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
 * Utilities for 3D multi data distributions -- implementation.
 */

#include "multiBlock/staticRepartitions3D.h"
#include "atomicBlock/dataField3D.hh"
#include "core/blockStatistics.hh"
#include "algorithm/basicAlgorithms.h"

namespace plb {

////////////////////// function createRegularMultiBlockDistribution3D /////////////////////

MultiBlockDistribution3D createRegularMultiBlockDistribution3D (
        plint nx, plint ny, plint nz, plint numBlocksX, plint numBlocksY, plint numBlocksZ,
        plint envelopeWidth )
{
    MultiBlockDistribution3D dataGeometry(nx, ny, nz);
    plint posX = 0;
    plint iBlock = 0;
    for (plint iBlockX=0; iBlockX<numBlocksX; ++iBlockX) {
        plint lx = nx / numBlocksX;
        if (iBlockX < nx%numBlocksX) ++lx;
        plint posY = 0;
        for (plint iBlockY=0; iBlockY<numBlocksY; ++iBlockY) {
            plint ly = ny / numBlocksY;
            if (iBlockY < ny%numBlocksY) ++ly;
            plint posZ = 0;
            for (plint iBlockZ=0; iBlockZ<numBlocksZ; ++iBlockZ) {
                plint lz = nz / numBlocksZ;
                if (iBlockZ < nz%numBlocksZ) ++lz;
                dataGeometry.addBlock(Box3D(posX, posX+lx-1, posY, posY+ly-1, posZ, posZ+lz-1),
                                      envelopeWidth, iBlock++);
                posZ += lz;
            }
            posY += ly;
        }
        posX += lx;
    }
    return dataGeometry;
}

MultiBlockDistribution3D createRegularMultiBlockDistribution3D (
        plint nx, plint ny, plint nz, plint envelopeWidth, int numProc) {
    std::vector<int> repartition = algorithm::evenRepartition(numProc, 3);
    return createRegularMultiBlockDistribution3D ( nx, ny, nz,
                                         repartition[0], repartition[1], repartition[2],
                                         envelopeWidth );
}

MultiBlockDistribution3D createXSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth )
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();
    
    std::vector<int> numActivePerSlice;
    plint numActiveTotal = 0;
    for(plint iX=0; iX<nX; iX++) {
        plint numActiveCurrentSlice = 0;
        for (plint iY=0; iY<nY; iY++) {
            for (plint iZ=0; iZ<nZ; iZ++) {
                if (cellTypeField.get(iX,iY,iZ) > 0) numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;
    
    MultiBlockDistribution3D dataGeometry(nX, nY, nZ);
    
    plint iX=0;
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        plint posX = iX;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock<numActivePerBlock && iX<nX) {
            numActiveCurrentBlock += numActivePerSlice[iX];
            iX++;
        }
        dataGeometry.addBlock(Box3D(posX, iX-1, 0, nY-1, 0, nZ-1), envelopeWidth, iBlock);
    }
    return dataGeometry;
}

MultiBlockDistribution3D createYSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth )
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();
    
    std::vector<int> numActivePerSlice;
    plint numActiveTotal = 0;
    for (plint iY=0; iY<nY; iY++) {
        plint numActiveCurrentSlice = 0;
        for(plint iX=0; iX<nX; iX++) {
            for(plint iZ=0; iZ<nZ; iZ++) {
                if (cellTypeField.get(iX,iY,iZ) > 0) numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;
    
    MultiBlockDistribution3D dataGeometry(nX, nY, nZ);
    
    plint iY=0;
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        plint posY = iY;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock<numActivePerBlock && iY<nY) {
            numActiveCurrentBlock += numActivePerSlice[iY];
            iY++;
        }
        dataGeometry.addBlock(Box3D(0, nX-1, posY, iY-1, 0, nZ-1), envelopeWidth, iBlock);
    }
    return dataGeometry;
}

MultiBlockDistribution3D createZSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks,
        plint envelopeWidth )
{
    plint nX = cellTypeField.getNx();
    plint nY = cellTypeField.getNy();
    plint nZ = cellTypeField.getNz();
    
    std::vector<int> numActivePerSlice;
    plint numActiveTotal = 0;
    for(plint iZ=0; iZ<nZ; iZ++) {
        plint numActiveCurrentSlice = 0;
        for(plint iX=0; iX<nX; iX++) {
            for (plint iY=0; iY<nY; iY++) {
                if (cellTypeField.get(iX,iY,iZ) > 0) numActiveCurrentSlice++;
            }
        }
        numActivePerSlice.push_back(numActiveCurrentSlice);
        numActiveTotal += numActiveCurrentSlice;
    }
    plint numActivePerBlock = numActiveTotal / numBlocks;    
    
    MultiBlockDistribution3D dataGeometry(nX, nY, nZ);
    
    plint iZ=0;
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        plint posZ = iZ;
        plint numActiveCurrentBlock = 0;
        while (numActiveCurrentBlock<numActivePerBlock && iZ<nZ) {
            numActiveCurrentBlock += numActivePerSlice[iZ];
            iZ++;
        }
        dataGeometry.addBlock(Box3D(0, nX-1, 0, nY-1, posZ, iZ-1), envelopeWidth, iBlock);
    }
    return dataGeometry;
}

MultiBlockDistribution3D createXSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth)
{
    return createXSlicedMultiBlockDistribution3D(cellTypeField, global::mpi().getSize(), envelopeWidth);
}

MultiBlockDistribution3D createYSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth)
{
    return createYSlicedMultiBlockDistribution3D(cellTypeField, global::mpi().getSize(), envelopeWidth);
}

MultiBlockDistribution3D createZSlicedMultiBlockDistribution3D (
        CellTypeField3D const& cellTypeField, plint envelopeWidth)
{
    return createZSlicedMultiBlockDistribution3D(cellTypeField, global::mpi().getSize(), envelopeWidth);
}

}  // namespace plb
