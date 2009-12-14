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
 * Templates for finding indexes for a specific subset of the neighborhood
 *  -- header file
 */
#ifndef INDEX_TEMPLATES_H
#define INDEX_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/array.h"
#include <algorithm>
#include <vector>

namespace plb {

namespace indexTemplates {

/// Compute the opposite of a given direction
template <typename Descriptor> inline plint opposite(plint iPop) {
    if (iPop==0) return 0;
    if (iPop<=Descriptor::q/2) return iPop + Descriptor::q/2;
    return iPop - Descriptor::q/2;
}

template <typename Descriptor>
plint findVelocity(Array<int,Descriptor::d> const& v) {
    for (plint iPop=0; iPop<Descriptor::q; ++iPop) {
        bool fit = true;
        for (int iD=0; iD<Descriptor::d; ++iD) {
            if (Descriptor::c[iPop][iD] != v[iD]) {
                fit = false;
                break;
            }
        }
        if (fit) return iPop;
    }
    return Descriptor::q;
}

template <typename Descriptor, plint index, plint value>
class SubIndex {
private:
    SubIndex() {
        for (int iVel=0; iVel<Descriptor::q; ++iVel) {
            if (Descriptor::c[iVel][index]==value) {
                indices.push_back(iVel);
            }
        }
    }

    std::vector<int> indices;

    template <typename Descriptor_, plint index_, plint value_>
    friend std::vector<int> const& subIndex();
};

template <typename Descriptor, plint index, plint value>
std::vector<int> const& subIndex() {
    static SubIndex<Descriptor, index, value> subIndexSingleton;
    return subIndexSingleton.indices;
}

/**
* finds distributions incoming into the wall
* but we want the ones outgoing from the wall,
* therefore we have to take the opposite ones.
*/
template <typename Descriptor, int direction, int orientation>
class SubIndexOutgoing {
private:
    SubIndexOutgoing() // finds the indexes outgoing from the walls
    {
        indices = indexTemplates::subIndex<Descriptor,direction,orientation>();

        for (pluint iPop = 0; iPop < indices.size(); ++iPop) {
            indices[iPop] = indexTemplates::opposite<Descriptor>(indices[iPop]);
        }

    }

    std::vector<int> indices;

    template <typename Descriptor_, int direction_, int orientation_>
    friend std::vector<int> const& subIndexOutgoing();
};

template <typename Descriptor, int direction, int orientation>
std::vector<int> const& subIndexOutgoing() {
    static SubIndexOutgoing<Descriptor, direction, orientation> subIndexOutgoingSingleton;
    return subIndexOutgoingSingleton.indices;
}

///finds all the remaining indexes of a lattice given some other indexes
template <typename Descriptor>
std::vector<int> remainingIndexes(const std::vector<int> &indices)
{
    std::vector<int> remaining;
    for (plint iPop = 0; iPop < Descriptor::q; ++iPop)
    {
        bool found = false;
        for (pluint jPop = 0; jPop < indices.size(); ++jPop)
        {
            if (indices[jPop] == iPop)
            {
                found = true;
            }
        }
        if (!found)
        {
            remaining.push_back(iPop);
        }
    }
    return remaining;
}

/// finds the indexes outgoing from a 2D corner
template <typename Descriptor, int xNormal, int yNormal>
class SubIndexOutgoingCorner2D {
private:
    SubIndexOutgoingCorner2D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect ( xNormal, yNormal );
        std::vector<int> knownIndexes;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = xNormal; vect[1] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0; vect[1] = yNormal;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0; vect[1] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        indices = indexTemplates::remainingIndexes<L>(knownIndexes);
    }

    std::vector<int> indices;

    template <typename Descriptor_, int direction_, int orientation_>
    friend std::vector<int> const& subIndexOutgoingCorner2D();
};

template <typename Descriptor, int xNormal, int yNormal>
std::vector<int> const& subIndexOutgoingCorner2D() {
    static SubIndexOutgoingCorner2D<Descriptor, xNormal, yNormal> subIndexOutgoingCorner2DSingleton;
    return subIndexOutgoingCorner2DSingleton.indices;
}

}  // namespace indexTemplates

}  // namespace plb

#endif
