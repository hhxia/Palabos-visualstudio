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

#ifndef TREE_READER_2D_HH
#define TREE_READER_2D_HH

#include <vector>
#include <queue>

#include "core/plbDebug.h"
#include "simulationSetup/dataFieldInitializer2D.h"
#include "libraryInterfaces/CVMLCPP_treeReader.h"
#include "libraryInterfaces/CVMLCPP_treeReader2D.h"

#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/Geometry"
#include "cvmlcpp/volume/VolumeIO"
#include "cvmlcpp/volume/Voxelizer"
#include "cvmlcpp/volume/DTree"


namespace plb {

/****************************************************
*********2D partial template specialization**********
*****************************************************/
template < typename T, typename U>
void convertToMatrix2D(cvmlcpp::DTree<T,2> tree, ScalarField2D<U> &result, plint envelopeWidth){
    typedef cvmlcpp::DTreeProxy<T,2> DNode2D;
    int height = treeDepth<T,2>(tree);

    // we suppose a field of size 0..2^height-1 x 0..2^height-1
    PLB_PRECONDITION(result.getNx() == (1 << height)+2*envelopeWidth);
    PLB_PRECONDITION(result.getNy() == (1 << height)+2*envelopeWidth);
    // initialize the mask
    for (plint iX = 0; iX < result.getNx(); ++iX){
        for (plint iY = 0; iY < result.getNy(); ++iY){
            result.get(iX, iY) = (U)false;
        }
    }

    // container for the results
    std::vector< Box2D > domains;
    
    // call to the treeUnboundedDecomposition2D function with big depth
    treeUnboundedDecomposition2D<T>(tree,domains);
    
    // iteration over the results
    std::vector< Box2D >::iterator it = domains.begin();
    while (it != domains.end()){
        Box2D newBox = *it;
        ++it;
        for (int iX = newBox.x0; iX < newBox.x1; ++iX){
            for (int iY = newBox.y0; iY < newBox.y1; ++iY){
                result.get(iX, iY) = true;
            }
        }
    }
}

/// Convert a table of indexes (0..2^D, D the tree dimension) into a Box2D 
Box2D decodeIndexes2D(int* indexes, int height){
    int sizeDomaine = 1;
    int maxSize = (1  << height);
    int x0, y0, x1, y1;
    x0 = y0 = x1 = y1 = 0;
    int blockSize = 1; // we start from the smaller block size
    for (int i = height; i > 0; --i){
        x0 +=  blockSize*(indexes[i] & 1);
        y0 += blockSize*((indexes[i] & 2) >> 1);
        blockSize *= 2;
    }
    x1 = x0 + maxSize/blockSize;
    y1 = y0 + maxSize/blockSize;
    return Box2D(x0,x1,y0,y1);
}

/// Adding a Box2D and the node that corresponds to it
template <typename T>
void TreeReader2D<T>::addBox2DAndNode(DNode2D &node, plint height){
    int* indexes = new int[height+1];
    for (int i = 0; i < height+1; ++i) indexes[i] = 0;
    findPath(tree,node,indexes);
    Box2D domain = decodeIndexes2D(indexes, node.depth());
    domain = domain.multiply((1<< height)/(1 << node.depth()));
    cutOffNodes.push_back(node);
    domains.push_back(domain);
    delete indexes;
}

template <typename T>
void treeUnboundedDecomposition2D(cvmlcpp::DTree<T,2> tree, std::vector< Box2D>& domains){
			
    typedef cvmlcpp::DTreeProxy<T,2> DNode2D;
    int currentDepth = 0;

    std::queue< DNode2D > q;
    q.push(tree.root());
    // retrieve the height of the tree for further computations
    int height = treeDepth<T,2>(tree);

    // table that contains the current index list 
    int* indexes = new int[height+1];
    for (int i = 0; i < height+1; ++i) indexes[i] = 0;
    
    while (!q.empty()){
        DNode2D node = q.front();
        q.pop();
        
        if (currentDepth < node.depth()){
            currentDepth = node.depth();
        }
        
        // reconstruction of the indexes up to this node
        findPath(tree,node,indexes);
                
        // if the node is a branch
        if (!node.isLeaf()){			
            // we examine all of its children if we have not reached a certain depth
            for (int i = 0; i < (1 << 2); ++i){
                q.push(node[i]);
            }
        }	
        else {
            if (node() == 1){
                Box2D domain = decodeIndexes2D(indexes, node.depth());
                domain = domain.multiply((1<< height)/(1 << node.depth()));
                domains.push_back(domain);
            }
        }
    }	
}

/// convert a node in ScalarField3D
template <typename T>
template <typename U>
void TreeReader2D<T>::createBooleanMask (
        plint domainNumber, ScalarField2D<U> &mask, plint envelopeWidth) const
{
  DNode2D node = cutOffNodes[domainNumber];
  if (isEntirelyLiquid(domainNumber)){
    // we return a mask full of true
    initializeAtConstant(mask, mask.getBoundingBox(), (U) true);
  }
  else {
    convertToMatrix2D(node.clone(),mask, envelopeWidth);
  }
}



/// Convert the DTree structure into a set of Box2D bounded by the depth
template <typename T>
void TreeReader2D<T>::analyzeTree(){
    plint currentDepth = 0;

    std::queue< DNode2D > q;
    q.push(tree.root());
    // retrieve the height of the tree for further computations
    int height = treeDepth<T,2>(tree);
    // compute the length size of the cube that contains the domain
    sideLength = 1 << height;
    // table that contains the current index list 
    
    
    while (!q.empty()){
        DNode2D node = q.front();
        q.pop();
        
        if (currentDepth < node.depth()){
            currentDepth = node.depth();
        }
    
        // if the node is a branch
        if (!node.isLeaf()){
            // we examine all of its children if we have not reached a certain depth
            if (node.depth() < limitHeight){
                for (int i = 0; i < (1 << 2); ++i){
                    q.push(node[i]);
                }
            }
            else {
                if (node.isLeaf()){
                    if (node() == 1){
                        addBox2DAndNode(node,height);
                    }
                }
                else {
                    // reconstruction of the indexes up to this node
                    addBox2DAndNode(node,height);
                }
            }
        }
        // if we have a leaf node that CONTAINS only liquid before the intended depth
        else {
            if (node() == 1){
                addBox2DAndNode(node,height);
            }
        }
    }

    for (pluint iDomain=0; iDomain<domains.size(); ++iDomain) {
        domains[iDomain].x1--;
        domains[iDomain].y1--;
    }
}

}  // namespace plb

#endif  // TREE_READER_2D_HH

#endif  // PLB_USE_CVMLCPP
