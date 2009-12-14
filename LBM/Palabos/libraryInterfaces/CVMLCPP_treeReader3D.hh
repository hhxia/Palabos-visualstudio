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

#ifndef TREE_READER_3D_HH
#define TREE_READER_3D_HH

#include <vector>
#include <queue>

#include "core/plbDebug.h"
#include "io/parallelIO.h"
#include "simulationSetup/dataFieldInitializer3D.h"
#include "libraryInterfaces/CVMLCPP_treeReader.h"
#include "libraryInterfaces/CVMLCPP_treeReader3D.h"

#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/Geometry"
#include "cvmlcpp/volume/VolumeIO"
#include "cvmlcpp/volume/Voxelizer"
#include "cvmlcpp/volume/DTree"


namespace plb {

/// 3D version of the tree to matrix convertion
template <typename T, typename U>
void convertToMatrix3D(cvmlcpp::DTree<T,3> tree, ScalarField3D<U> &result, plint envelopeWidth){
    typedef cvmlcpp::DTreeProxy<T,3> DNode3D;
    int height = treeDepth<T,3>(tree);
    
    // we suppose a field of size 0..2^height-1 x 0..2^height-1 x 0..2^height-1
    PLB_PRECONDITION(result.getNx() == (1 << height)+2*envelopeWidth);
    PLB_PRECONDITION(result.getNy() == (1 << height)+2*envelopeWidth);
    PLB_PRECONDITION(result.getNz() == (1 << height)+2*envelopeWidth);
    // initialize to false (solid nodes)
    for (plint iX = 0; iX < result.getNx(); ++iX){
        for (plint iY = 0; iY < result.getNy(); ++iY){
            for (plint iZ = 0; iZ < result.getNz(); ++iZ){
                result.get(iX, iY,iZ) = (U)false;
            }
        }
    }
    
    // container for the results
    std::vector< Box3D> domains;
    // decompose the whole tree in Box3D which contain 1
    treeUnboundedDecomposition3D<T>(tree,domains);
    
    // iteration over the results
    std::vector< Box3D >::iterator it = domains.begin();
    while (it != domains.end()){
        Box3D box = *it;
        for (plint iX = box.x0; iX < box.x1; ++iX){
            for (plint iY = box.y0; iY < box.y1; ++iY){
                for (plint iZ = box.z0; iZ < box.z1; ++iZ){
                    result.get(iX+envelopeWidth, iY+envelopeWidth,iZ+envelopeWidth) = (U)true;
                }
            }
        }
        ++it;
    }
}


/// Convert the DTree structure into a set of Box3D without any depth bound
template <typename T>
void treeUnboundedDecomposition3D(cvmlcpp::DTree<T,3> tree, std::vector< Box3D>& decomposedDomains){
    typedef cvmlcpp::DTreeProxy<T,3> DNode3D;
    int currentDepth = 0;
    
    std::queue< DNode3D > q;
    q.push(tree.root());
    // retrieve the height of the tree for further computations
    int height = treeDepth<T,3>(tree);
    
    // table that contains the current index list 
    std::vector<int> indexes(height+1);
    for (int i = 0; i < height+1; ++i) indexes[i] = 0;
    
    while (!q.empty()){
        DNode3D node = q.front();
        q.pop();
        
        // multiply all domains by two if we go deeper
        if (currentDepth < node.depth()){
            currentDepth = node.depth();
            // multiply every vector element by 3
            std::vector<Box3D>::iterator it = decomposedDomains.begin();
            while (it != decomposedDomains.end()){
                (*it) = (*it).multiply(2);
                ++it;
            }
        }
    
        // if the node is a branch
        if (!node.isLeaf()){			
            // we examine all of its children if we have not reached a certain depth
            for (int i = 0; i < (1 << 3); ++i){
                    q.push(node[i]);
            }
        }
        else 	{
            if (node() == 1){
                findPath(tree, node, indexes);
                Box3D domain = decodeIndexes3D(indexes, node.depth());
                decomposedDomains.push_back(domain);
            }
        }
    }	
}

template <typename T>
std::vector<Box3D> TreeReader3D<T>::getDomains() const {
    return domains;
}

template <typename T>
bool TreeReader3D<T>::isEntirelyLiquid(plint domainNumber) const {
    return cutOffNodes[domainNumber].isLeaf();
}

/// convert a node in ScalarField3D
template <typename T>
template <typename U>
void TreeReader3D<T>::createBooleanMask(plint domainNumber, ScalarField3D<U> &mask, plint envelopeWidth) const {
    DNode3D node = cutOffNodes[domainNumber];
    if (isEntirelyLiquid(domainNumber)){
        // We return a mask full of false (fluid-only)
        setToConstant(mask, mask.getBoundingBox(), (U) true);
    }
    else {
        convertToMatrix3D(node.clone(),mask, envelopeWidth);
    }
}

template <typename T>
void TreeReader3D<T>::addBox3DAndNode(DNode3D &node, plint height){
    std::vector<int> indexes(height+1);
    for (int i = 0; i < height+1; ++i) indexes[i] = 0;
    findPath(tree,node,indexes);
    Box3D domain = decodeIndexes3D(indexes, node.depth());
    domain = domain.multiply((1<< height)/(1 << node.depth()));
    // Subdivide all the domain in case it is a leaf node before the limit depth
    plint dividedDomainNumber = 1 << 3*(limitHeight-node.depth());
    plint divisionFactor = 1 << (limitHeight-node.depth()); // Cubic root of the dividedDomainNumber
    if (equalSizeBlocks && node.isLeaf() && dividedDomainNumber > 1){
        plint lx = domain.getNx()-1;
        plint ly = domain.getNy()-1;
        plint lz = domain.getNz()-1;
        for (plint iX = 0; iX < divisionFactor; ++iX){
            for (plint iY = 0; iY < divisionFactor; ++iY){
                for (plint iZ = 0; iZ < divisionFactor; ++iZ){
                    Box3D newDomain(domain.x0 + iX*(lx/divisionFactor), domain.x0 + (iX+1)*lx/divisionFactor,
                                    domain.y0 + iY*(ly/divisionFactor), domain.y0 + (iY+1)*ly/divisionFactor,
                                    domain.z0 + iZ*(lz/divisionFactor), domain.z0 + (iZ+1)*lz/divisionFactor );
                    cutOffNodes.push_back(node);
                    domains.push_back(newDomain);
                }
            }
        } 
    }
    else {
        cutOffNodes.push_back(node);
        domains.push_back(domain);
    }
}

template <typename T>
void TreeReader3D<T>::analyzeTree(){
    plint currentDepth = 0;

    std::queue< DNode3D > q;
    q.push(tree.root());
    // retrieve the height of the tree for further computations
    int height = treeDepth<T,3>(tree);
    // compute the length size of the cube that contains the domain
    sideLength = 1 << height;
    // table that contains the current index list 
    
    
    while (!q.empty()){
        DNode3D node = q.front();
        q.pop();
        
        if (currentDepth < node.depth()){
            currentDepth = node.depth();
        }
    
        // if the node is a branch
        if (!node.isLeaf()){
            // we examine all of its children if we have not reached a certain depth
            if (node.depth() < limitHeight){
                for (int i = 0; i < (1 << 3); ++i){
                    q.push(node[i]);
                }
            }
            else {
                if (node.isLeaf()){
                    if (node() == 1){
                        addBox3DAndNode(node,height);
                    }
                }
                else {
                    // reconstruction of the indexes up to this node
                    addBox3DAndNode(node,height);
                }
            }
            
        }
        // if we have a leaf node that CONTAINS only liquid before the intended depth
        else {
            if (node() == 1){
                addBox3DAndNode(node,height);
            }
        }
    }

    for (pluint iDomain=0; iDomain<domains.size(); ++iDomain) {
        domains[iDomain].x1--;
        domains[iDomain].y1--;
        domains[iDomain].z1--;
    }
}

}  // namespace plb

#endif  // TREE_READER_3D_HH

#endif  // PLB_USE_CVMLCPP
