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

#ifndef TREE_READER_2D_H
#define TREE_READER_2D_H

#include "core/geometry2D.h"
#include "atomicBlock/dataField2D.h"
#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/DTree"

namespace plb {

template <typename T>
void treeUnboundedDecomposition2D(cvmlcpp::DTree<T,2> tree,std::vector< Box2D>& domains);

/// class that converts a DTree into a vector of Box2D and scalar fields
template <typename T>
class TreeReader2D {
    public:
        typedef cvmlcpp::DTreeProxy<T,2> DNode2D;
    private:
        cvmlcpp::DTree<T,2> tree;
        plint limitHeight;
        std::vector<Box2D> domains;
        std::vector<DNode2D> cutOffNodes;
        int sideLength;

        void analyzeTree();
        void addBox2DAndNode(DNode2D &node, plint height);
    public:
        // constructors
        TreeReader2D(cvmlcpp::DTree<T,2> tree_, int limitHeight_) :
            tree(tree_), limitHeight(limitHeight_)
        { analyzeTree();}
        // destructor
        ~TreeReader2D(){}

        // operations
        std::vector<Box2D> getDomains() const;
        bool isEntirelyLiquid(plint domainNumber) const;
        template <class U> void createBooleanMask(plint domainNumber, ScalarField2D<U> &mask, plint envelopeWidth) const;

        // getters
        plint getNx(){return sideLength;}
        plint getNy(){return sideLength;}
};

}  // namespace plb

#endif  // TREE_READER_2D_H

#endif  // PLB_USE_CVMLCPP

#include "libraryInterfaces/CVMLCPP_treeReader2D.hh"
