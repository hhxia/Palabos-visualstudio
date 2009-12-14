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

#ifndef TREE_READER_3D_H
#define TREE_READER_3D_H

#include "core/geometry3D.h"
#include "atomicBlock/dataField3D.h"
#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/DTree"
#include <vector>

namespace plb {

template <typename T>
void treeUnboundedDecomposition3D(cvmlcpp::DTree<T,3> tree, std::vector< Box3D>& domains);

template <typename T>
class TreeReader3D {
    public:
        typedef cvmlcpp::DTreeProxy<T,3> DNode3D;
    private:
        cvmlcpp::DTree<T,3> tree;
        plint limitHeight;
        bool equalSizeBlocks;
        std::vector<Box3D> domains;
        std::vector<DNode3D> cutOffNodes;
        int sideLength;
        
        void analyzeTree();
        void addBox3DAndNode(DNode3D &node, plint height);
        
    public:
        // Constructor
        TreeReader3D(cvmlcpp::DTree<T,3> tree_, 
                     plint limitHeight_, 
                     bool equalSizeBlocks_=false) : tree(tree_),
                                                    limitHeight(limitHeight_),
                                                    equalSizeBlocks(equalSizeBlocks_)
        { analyzeTree();}
        
        // Destructor
        ~TreeReader3D(){}
        
        std::vector<Box3D> getDomains() const;
        bool isEntirelyLiquid(plint domainNumber) const;
        template <class U> void createBooleanMask(plint domainNumber, ScalarField3D<U> &mask, plint envelopeWidth) const;

        // Getters
        plint getNx(){return sideLength;}
        plint getNy(){return sideLength;}
        plint getNz(){return sideLength;}
};

void computeDomainParameters3D(std::string fName, plint N,
                               double& voxelSize, Array<double,3>& offset);

TreeReader3D<int>* createUniformTreeReader3D(std::string fName, plint N, int limitDepth);

TreeReader3D<int>* createUniformTreeReader3D(cvmlcpp::DTree<int,3u> &octree, plint N, int limitDepth);


Box3D decodeIndexes3D(std::vector<int>& indexes, int height);

}  // namespace plb

#endif  // TREE_READER_3D_H

#endif  // PLB_USE_CVMLCPP

#include "libraryInterfaces/CVMLCPP_treeReader3D.hh"
