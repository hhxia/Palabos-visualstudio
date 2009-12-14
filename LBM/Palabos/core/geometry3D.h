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


#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include <algorithm>
#include <iterator>
#include <vector>
using namespace std;

namespace plb {

/// Coordinates of a 3D Box
struct Box3D {
    Box3D() : x0(), x1(), y0(), y1(), z0(), z1() { }
    Box3D(plint x0_, plint x1_, plint y0_, plint y1_, plint z0_, plint z1_)
        : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
    { }
    /// Return same box, shifted by (deltaX,deltaY,deltaZ)
    Box3D shift(plint deltaX, plint deltaY, plint deltaZ) const {
        return Box3D(x0+deltaX, x1+deltaX, y0+deltaY, y1+deltaY, z0+deltaZ, z1+deltaZ);
    }
    /// Return same box, rescaled by a factor scaling
    Box3D multiply(plint scaling) const {
        return Box3D(scaling*x0, scaling*x1, scaling*y0,
                     scaling*y1, scaling*z0, scaling*z1);
    }
    /// Return same box, rescaled by a factor 1/scaling
    Box3D divide(plint scaling) const {
        return Box3D(x0/scaling, x1/scaling, y0/scaling,
                     y1/scaling, z0/scaling, z1/scaling);
    }
    /// Rescale by 1/scaling and make sure to fit into higher-level box by trimming excess.
    /** This function should be used with grid refinement if the 
     *  coordinates of a fine grid are rescaled to fit on a coarse grid.
     *  The lower bounds of the resulting coarse box are increased by one
     *  if their fine-grid original version was not divisible by scaling;
     *  and the upper bounds are decreased in the same way.
     *  This makes sure that the coarse grid result does not exceed the
     *  bounds of the fine grid.
     */
    Box3D divideAndFitSmaller(plint scaling) const {
        return Box3D (util::roundUp(x0,scaling)/scaling, util::roundDown(x1,scaling)/scaling,
                      util::roundUp(y0,scaling)/scaling, util::roundDown(y1,scaling)/scaling,
                      util::roundUp(z0,scaling)/scaling, util::roundDown(z1,scaling)/scaling);
    }

    /// Rescale by 1/scaling and make sure to fit into higher-level box by snapping to larger size.
    /** This function should be used with grid refinement if the 
     *  coordinates of a fine grid are rescaled to fit on a coarse grid.
     *  The lower bounds of the resulting coarse box are decreased by one
     *  if their fine-grid original version was not divisible by scaling;
     *  and the upper bounds are increased in the same way.
     *  This makes sure that the coarse grid result contains the
     *  bounds of the fine grid.
     */
    Box3D divideAndFitLarger(plint scaling) const {
        return Box3D (util::roundDown(x0,scaling)/scaling, util::roundUp(x1,scaling)/scaling,
                      util::roundDown(y0,scaling)/scaling, util::roundUp(y1,scaling)/scaling,
                      util::roundDown(z0,scaling)/scaling, util::roundUp(z1,scaling)/scaling);
    }

    /// Add a border of nCells cells to the box
    Box3D enlarge(plint nCells) const {
        return Box3D(x0-nCells, x1+nCells, y0-nCells, y1+nCells, z0-nCells, z1+nCells);
    }
    /// Number of cells in x-direction
    plint getNx()  const { return (x1-x0+1); }
    /// Number of cells in y-direction
    plint getNy()  const { return (y1-y0+1); }
    /// Number of cells in z-direction
    plint getNz()  const { return (z1-z0+1); }
    /// Total number of cells in the box
    plint nCells() const { return getNx()*getNy()*getNz(); }
    /// Return the maximum of getNx(), getNy(), and getNz()
    plint getMaxWidth() const { return max(max(getNx(), getNy()), getNz()); }

    plint x0, x1, y0, y1, z0, z1;
};

/// Coordinates of a 3D point
struct Dot3D {
    Dot3D() : x(), y(), z() { }
    Dot3D(plint x_, plint y_, plint z_) : x(x_), y(y_), z(z_) { }
    plint x, y, z;
};

/// List of 3D points, used to describe a subdomain
struct DotList3D {
    DotList3D() { }
    DotList3D(std::vector<Dot3D> const& dots_) : dots(dots_) { }
    /// Add one more poplint to the list
    void addDot(Dot3D dot) {
        dots.push_back(dot);
    }
    /// Add more points to the list
    void addDots(std::vector<Dot3D> const& dots_) {
        dots.insert(dots.end(), dots_.begin(), dots_.end());
    }
    /// Get const reference to one of the dots
    Dot3D const& getDot(plint whichDot) const {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Get non-const reference to one of the dots
    Dot3D& getDot(plint whichDot) {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Return same poplint list, shifted by (deltaX,deltaY,deltaZ)
    DotList3D shift(plint deltaX, plint deltaY, plint deltaZ) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x+deltaX, it->y+deltaY, it->z+deltaZ) );
        }
        return DotList3D(newDots);
    }
    /// Return same box, rescaled by a factor scaling
    DotList3D multiply(plint scaling) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x*scaling, it->y*scaling, it->z*scaling) );
        }
        return DotList3D(newDots);
    }
    /// Return same box, rescaled by a factor 1/scaling
    DotList3D divide(plint scaling) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x/scaling, it->y/scaling, it->z/scaling) );
        }
        return DotList3D(newDots);

    }
    /// Get total number of points
    plint getN() const {
        return dots.size();
    }

    std::vector<Dot3D> dots;
};

/// Decide if lattice poplint is contained in 3D box, boundaries inclusive
inline bool contained(plint x, plint y, plint z, Box3D const& box) {
    return x>=box.x0 && x<=box.x1 &&
           y>=box.y0 && y<=box.y1 &&
           z>=box.z0 && z<=box.z1;
}

/// Decide if lattice poplint is contained in 3D box, boundaries inclusive
inline bool contained(Dot3D dot, Box3D const& box) {
    return contained(dot.x, dot.y, dot.z, box);
}

/// Decide if 3D box1 is contained in 3D box2, boundaries inclusive
inline bool contained(Box3D const& box1, Box3D const& box2) {
    return box1.x0>=box2.x0 && box1.x1<=box2.x1 &&
           box1.y0>=box2.y0 && box1.y1<=box2.y1 &&
           box1.z0>=box2.z0 && box1.z1<=box2.z1;
}

/// Compute intersection between two 3D boxes
/** \return false if boxes don't intersect
 */
using namespace std;
inline bool intersect(Box3D const& box1, Box3D const& box2, Box3D& inters) {
    inters.x0 = max(box1.x0,box2.x0);
    inters.y0 = max(box1.y0,box2.y0);
    inters.z0 = max(box1.z0,box2.z0);

    inters.x1 = min(box1.x1,box2.x1);
    inters.y1 = min(box1.y1,box2.y1);
    inters.z1 = min(box1.z1,box2.z1);

    return inters.x1>=inters.x0 && inters.y1>=inters.y0 && inters.z1>=inters.z0;
}

/// Compute intersection between a 3D box and a 3D DotList
/** \return false if the two don't intersect
 */
inline bool intersect(Box3D const& box, DotList3D const& dotlist, DotList3D& inters) {
    inters = DotList3D();
    std::vector<Dot3D>::const_iterator it = dotlist.dots.begin();
    for (; it != dotlist.dots.end(); ++it) {
        if ( contained(it->x,it->y,it->z, box) ) {
            inters.addDot(*it);
        }
    }
    return !inters.dots.empty();
}

/// Except the domain of box "toExcept" form the domain of box "originalBox"
/** The result consists of three boxes, which are added to the vector "result"
 */
inline void except(Box3D const& originalBox, Box3D const& toExcept, std::vector<Box3D>& result)
{
    Box3D intersection;
    bool doesIntersect = intersect(originalBox, toExcept, intersection);
    if (!doesIntersect) {
        result.push_back(originalBox);
    }
    else {
        if (intersection.x0 != originalBox.x0) {
            result.push_back(Box3D(originalBox.x0, intersection.x0-1, originalBox.y0, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.x1 != originalBox.x1) {
            result.push_back(Box3D(intersection.x1+1, originalBox.x1, originalBox.y0, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.y0 != originalBox.y0) {
            result.push_back(Box3D(intersection.x0, intersection.x1, originalBox.y0, intersection.y0-1, originalBox.z0, originalBox.z1));
        }
        if (intersection.y1 != originalBox.y1) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y1+1, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.z0 != originalBox.z0) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y0, intersection.y1, originalBox.z0, intersection.z0-1));
        }
        if (intersection.z1 != originalBox.z1) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y0, intersection.y1, intersection.z1+1, originalBox.z1));
        }
    }
}


} // namespace plb

#endif
