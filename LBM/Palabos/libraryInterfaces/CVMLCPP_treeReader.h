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

#ifndef TREE_READER_H
#define TREE_READER_H

#include "cvmlcpp/base/Matrix"
#include "cvmlcpp/volume/DTree"
#include <vector>

namespace plb {

/// Use BFS to find the height of a DTree
template <typename T, std::size_t D> int treeDepth(cvmlcpp::DTree<T,D> tree);

/// using DFS to find a given node and getting the index collection that created it
template <typename T, std::size_t D>
void findPath(cvmlcpp::DTree<T,D> tree, cvmlcpp::DTreeProxy<T, D> &wantedNode,
              std::vector<int>& indexes);

}

#endif  // TREE_READER_H

#endif  // PLB_USE_CVMLCPP
