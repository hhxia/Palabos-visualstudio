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
 * Neumann boundary conditions -- generic implementation.
 */
#ifndef NEUMANN_CONDITION_3D_HH
#define NEUMANN_CONDITION_3D_HH

#include "boundaryCondition/neumannCondition3D.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    std::vector<int> const& unknownIndices = indexTemplates::subIndex<Descriptor<T>, direction, -orientation>();
    enum {
        normalX = direction==0 ? orientation : 0,
        normalY = direction==1 ? orientation : 0,
        normalZ = direction==2 ? orientation : 0
    };
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (pluint fIndex=0; fIndex<unknownIndices.size(); ++fIndex) {
                    plint iPop = unknownIndices[fIndex];
                    lattice.get(iX,iY,iZ)[iPop] = lattice.get(iX-normalX, iY-normalY, iZ-normalZ)[iPop];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>*
    CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>::clone() const
{
    return new CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
BlockDomain::DomainT CopyUnknownPopulationsFunctional3D<T,Descriptor,direction,orientation>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    lattice.get(iX,iY,iZ)[iPop] = lattice.get(iX-normalX, iY-normalY, iZ-normalZ)[iPop];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>*
    CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>::clone() const
{
    return new CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
BlockDomain::DomainT CopyAllPopulationsFunctional3D<T,Descriptor,normalX,normalY,normalZ>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> u;
                lattice.get(iX-normalX, iY-normalY, iZ-normalZ).computeVelocity(u);
                lattice.get(iX, iY, iZ).defineVelocity(u);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>*
    CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::clone() const
{
    return new CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
BlockDomain::DomainT CopyVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> u;
                lattice.get(iX-normalX, iY-normalY, iZ-normalZ).computeVelocity(u);
                if (normalX!=0) {
                    u[0] = T();
                }
                if (normalY!=0) {
                    u[1] = T();
                }
                if (normalZ!=0) {
                    u[2] = T();
                }
                lattice.get(iX, iY, iZ).defineVelocity(u);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>*
    CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::clone() const
{
    return new CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
BlockDomain::DomainT CopyTangentialVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> u;
                lattice.get(iX-normalX, iY-normalY, iZ-normalZ).computeVelocity(u);
                if (normalX==0) {
                    u[0] = T();
                }
                if (normalY==0) {
                    u[1] = T();
                }
                if (normalZ==0) {
                    u[2] = T();
                }
                lattice.get(iX, iY, iZ).defineVelocity(u);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>*
    CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::clone() const
{
    return new CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
BlockDomain::DomainT CopyNormalVelocityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX, iY, iZ).defineDensity (
                        lattice.get(iX-normalX, iY-normalY, iZ-normalZ).computeDensity() );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>*
    CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::clone() const
{
    return new CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>(*this);
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
void CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor, int normalX, int normalY, int normalZ> 
BlockDomain::DomainT CopyDensityFunctional3D<T,Descriptor,normalX,normalY,normalZ>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // NEUMANN_CONDITION_3D_HH
