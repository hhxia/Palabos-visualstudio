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
 * Automatic creation of a 3D non-sparse multi-block -- generatic implementation
 */

#ifdef PLB_USE_CVMLCPP

#ifndef FULL_MULTI_BLOCK_3D_HH
#define FULL_MULTI_BLOCK_3D_HH

#include "libraryInterfaces/CVMLCPP_fullMultiBlock3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"
#include <cvmlcpp/volume/Geometry>
#include <cvmlcpp/volume/Voxelizer>

using namespace std;

namespace plb {

template<typename T, template<typename U> class Descriptor>
class DynamicsFromBoolMaskFunctional3D
    : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    DynamicsFromBoolMaskFunctional3D(Dynamics<T,Descriptor>* dynamics_)
        : dynamics(dynamics_)
    { }
    ~DynamicsFromBoolMaskFunctional3D() {
        delete dynamics;
    }
    DynamicsFromBoolMaskFunctional3D(DynamicsFromBoolMaskFunctional3D<T,Descriptor> const& rhs)
        : dynamics(rhs.dynamics->clone())
    { }
    DynamicsFromBoolMaskFunctional3D<T,Descriptor>& operator= (
            DynamicsFromBoolMaskFunctional3D<T,Descriptor> const& rhs )
    {
        delete dynamics; dynamics = rhs.dynamics->clone();
        return *this;
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& geometry);
    virtual DynamicsFromBoolMaskFunctional3D<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    Dynamics<T,Descriptor>* dynamics;
};

template<typename T, template<typename U> class Descriptor>
void DynamicsFromBoolMaskFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& geometry )
{
    // Take the defensive point of view that lattice and geometry are possibly
    //   shifted in space, and correct their relative displacement.
    Dot3D offset = computeRelativeDisplacement(lattice, geometry);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Round the real-valued geometry flag to the nearest integer, then cast to bool.
                bool wallNode = (bool) util::roundToInt(geometry.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                // If required, instantiate bounce-back dynamics on current node.
                if (wallNode) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromBoolMaskFunctional3D<T,Descriptor>* DynamicsFromBoolMaskFunctional3D<T,Descriptor>::clone() const
{
    return new DynamicsFromBoolMaskFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromBoolMaskFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
ScalarField3D<T>* boolMaskFromSTL(std::string stlFile, plint N)
{
    cvmlcpp::Geometry<float> geometry;
    cvmlcpp::readSTL(geometry, stlFile);
    
    // Find dimensions
    double geometrySize = 0.;
    for(int d = 0; d < 3; ++d) {
        geometrySize = max (
                geometrySize,
                double(geometry.max(d)) - double(geometry.min(d))
        );
    }
    PLB_ASSERT( geometrySize > 0.0 );
    cvmlcpp::Matrix<int,3u> matrix;
    double voxelSize = geometrySize / (double)N;
    cvmlcpp::voxelize(geometry, matrix, voxelSize);
    
    plint nx = *matrix.extents();
    plint ny = *(matrix.extents()+1);
    plint nz = *(matrix.extents()+2);

    ScalarField3D<T>* boolMask = new ScalarField3D<T>(nx,ny,nz);
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                boolMask->get(iX,iY,iZ) = (T)matrix[iX][iY][iZ];
            }
        }
    }

    return boolMask;
}

template<typename T, template <typename U> class Descriptor>
void dynamicsFromBoolMask (
        ScalarField3D<T>& boolMask,
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Dynamics<T,Descriptor>* dynamics, bool flag )
{
    MultiScalarField3D<T> multiBoolMask(lattice);
    copySerializedBlock(boolMask, multiBoolMask);

    defineDynamics(lattice, multiBoolMask, dynamics, flag);
}

}  // namespace plb

#endif  // FULL_MULTI_BLOCK_3D_HH

#endif  // PLB_USE_CVMLCPP
