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
 * Helper functions for domain initialization -- header file.
 */
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include <string>


#ifndef MULTI_BLOCK_INFO_3D_H
#define MULTI_BLOCK_INFO_3D_H

namespace plb {

template<typename T>
std::string getMultiBlockInfo(MultiBlock3D<T> const& multiBlock);

template<typename T>
std::string getMultiBlockInfo(MultiBlock3D<T> const& multiBlock) {
    plint nx = multiBlock.getNx();
    plint ny = multiBlock.getNy();
    plint nz = multiBlock.getNz();

    MultiBlockManagement3D const& management = multiBlock.getMultiBlockManagement();
    MultiBlockDistribution3D const& distribution = management.getMultiBlockDistribution();
    if (distribution.getNumBlocks()==0) {
        return std::string("Empty multi-block\n");
    }
    plint maxNumCells = distribution.getBlockParameters(0).getBulk().nCells();
    plint largestBlock = 0;
    plint minNumCells = distribution.getBlockParameters(0).getBulk().nCells();
    plint smallestBlock = 0;
    plint numAllocatedCells = 0;
    for (pluint iBlock=0; iBlock<distribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const&  parameters = distribution.getBlockParameters(iBlock);
        plint numCells = parameters.getBulk().nCells();
        numAllocatedCells += numCells;
        if (numCells>maxNumCells) {
            maxNumCells = numCells;
            largestBlock = iBlock;
        }
        if (numCells<minNumCells) {
            minNumCells = numCells;
            smallestBlock = iBlock;
        }
    }

    BlockParameters3D const& smallest = distribution.getBlockParameters(smallestBlock);
    BlockParameters3D const& largest = distribution.getBlockParameters(largestBlock);
    std::stringstream blockInfo;
    blockInfo << "Size of the multi-block:     " << nx << "-by-" << ny << "-by-" << nz << "\n";
    blockInfo << "Number of atomic-blocks:     " << distribution.getNumBlocks() << "\n";
    blockInfo << "Smallest atomic-block:       " << smallest.getBulk().getNx() << "-by-"
                                                 << smallest.getBulk().getNy() << "-by-"
                                                 << smallest.getBulk().getNz() << "\n";
    blockInfo << "Largest atomic-block:        "  << largest.getBulk().getNx()  << "-by-"
                                                  << largest.getBulk().getNy() << "-by-"
                                                  << largest.getBulk().getNz() << "\n";
    blockInfo << "Number of allocated cells:   " << (T)numAllocatedCells/1.e6 << " million\n";
    blockInfo << "Fraction of allocated cells: " << (T)numAllocatedCells/(T)(nx*ny*nz)*100 << " percent\n";

    return blockInfo.str();
}

}  // namespace plb

#endif  // MULTI_BLOCK_INFO_3D_H
