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

/* Main author: Daniel Lagrava
 */

#ifdef PLB_USE_CVMLCPP

#include "parallelism/mpiManager.h"
#include "libraryInterfaces/CVMLCPP_treeReader.h"
#include "libraryInterfaces/CVMLCPP_treeReader.hh"
#include "libraryInterfaces/CVMLCPP_treeReader3D.h"
#include "libraryInterfaces/CVMLCPP_treeReader3D.hh"
#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/DTree"
#include "cvmlcpp/volume/Geometry"
#include "cvmlcpp/volume/VolumeIO"
#include "cvmlcpp/volume/Voxelizer"

using namespace std;

namespace plb {

/// decoding a table of height indexes ranging between 1 and 7 
Box3D decodeIndexes3D(std::vector<int>& indexes, int height)
{
    PLB_PRECONDITION( indexes.size() >= height+1 );

    int sizeDomaine = 1;
    int maxSize = (1  << height);
    int x0, y0, z0, x1, y1, z1;
    x0 = y0 = x1 = y1 = z0 = z1 = 0;
    int blockSize = 1; // we start from the smaller block size
    for (int i = height; i > 0; --i){
        x0 +=  blockSize*(indexes[i] & 1); // first bit
        y0 += blockSize*((indexes[i] & 2) >> 1); // second bit
        z0 += blockSize*((indexes[i] & 4) >> 2); // third bit 
        blockSize *= 2;
    }
    x1 = x0 + maxSize/blockSize;
    y1 = y0 + maxSize/blockSize;
    z1 = z0 + maxSize/blockSize;
        
    return  Box3D(x0,x1,y0,y1,z0,z1);
}

void computeDomainParameters3D(std::string fName, plint N,
                               double& voxelSize, Array<double,3>& offset)
{
    cvmlcpp::Geometry<float> geometry;
    cvmlcpp::readSTL(geometry, fName);

    double geometrySize = 0.;
    for(int d = 0; d < 3; ++d) {
        geometrySize = max (
                geometrySize,
                double(geometry.max(d)) - double(geometry.min(d))
        );
        offset[d] = geometry.min(d);
    }
    PLB_ASSERT(geometrySize > 0.0);
    voxelSize = geometrySize / (double) N;
}
    
TreeReader3D<int>* createUniformTreeReader3D(std::string fName, plint N, int limitDepth)
{
    cvmlcpp::Geometry<float> geometry;
    cvmlcpp::readSTL(geometry, fName);
    
    // Voxelization.
    double geometrySize = 0.;
    for(int d = 0; d < 3; ++d) {
        geometrySize = max (
                geometrySize,
                double(geometry.max(d)) - double(geometry.min(d))
        );
    }
    PLB_ASSERT(geometrySize > 0.0);
    cvmlcpp::DTree<int,3u> octree(0);
    double voxelSize = geometrySize / (double) N;
    // Create the octree.
    cvmlcpp::voxelize(geometry, octree, voxelSize);
   
    bool equalSizeBlocks = global::mpi().getSize()==1 ? false : true;
    // Create the reader.
    return new TreeReader3D<int>(octree,limitDepth,equalSizeBlocks);
}

TreeReader3D<int>* createUniformTreeReader3D(cvmlcpp::DTree<int,3u> &octree, plint N, int limitDepth)
{
    bool equalSizeBlocks = global::mpi().getSize()==1 ? false : true;
    // Create the reader.
    return new TreeReader3D<int>(octree,limitDepth,equalSizeBlocks);
}




}  // namespace cfl

#endif  // PLB_USE_CVMLCPP
