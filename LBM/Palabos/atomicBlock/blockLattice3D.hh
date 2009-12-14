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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_3D_HH
#define BLOCK_LATTICE_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/util.h"
#include <algorithm>
#include <typeinfo>
using namespace std;
namespace plb {

// Class BlockLattice3D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::BlockLattice3D (
        plint nx_, plint ny_, plint nz_,
        Dynamics<T,Descriptor>* backgroundDynamics_ )
    : nx(nx_), ny(ny_), nz(nz_),
      backgroundDynamics(backgroundDynamics_),
      dataTransfer(*this)
{
    // Allocate memory, and initialize dynamics.
    allocateMemory();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                grid[iX][iY][iZ].attributeDynamics(backgroundDynamics);
            }
        }
    }
    // Attribute default value to the standard statistics (average uSqr,
    //   max uSqr, average rho). These have previously been subscribed
    //   in the constructor of BlockLatticeBase3D.
    std::vector<T> average, sum, max;
    std::vector<plint> intSum;
    average.push_back(Descriptor<T>::rhoBar((T)1));
                             // default average rho to 1, to avoid division by
                             // zero in constRhoBGK and related models
    average.push_back(T());  // default average uSqr to 0
    max.push_back(T());      // default max uSqr to 0
    plint numCells = 1;       // pretend fictitious cell to evaluate statistics
    this->getInternalStatistics().evaluate (average, sum, max, intSum, numCells);
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::~BlockLattice3D()
{
    releaseMemory();
}

/** The whole data of the lattice is duplicated. This includes
 * both particle distribution function and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>::BlockLattice3D(BlockLattice3D<T,Descriptor> const& rhs)
    : BlockLatticeBase3D<T,Descriptor>(rhs),
      AtomicBlock3D<T>(rhs),
      nx(rhs.nx),
      ny(rhs.ny),
      nz(rhs.nz),
      backgroundDynamics(rhs.backgroundDynamics->clone()),
      dataTransfer(*this)
{
    allocateMemory();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                Cell<T,Descriptor>& cell = grid[iX][iY][iZ];
                // Assign cell from rhs
                cell = rhs.grid[iX][iY][iZ];
                // Get an independent clone of the dynamics,
                //   or assign backgroundDynamics
                if (&cell.getDynamics()==rhs.backgroundDynamics) {
                    cell.attributeDynamics(backgroundDynamics);
                }
                else {
                    cell.attributeDynamics(cell.getDynamics().clone());
                }
            }
        }
    }
}

/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated. This includes both particle distribution function
 * and external fields. 
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice3D<T,Descriptor>& BlockLattice3D<T,Descriptor>::operator= (
        BlockLattice3D<T,Descriptor> const& rhs )
{
    BlockLattice3D<T,Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the 
 * lattice are copied, and not the lattice itself.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::swap(BlockLattice3D& rhs) {
    BlockLatticeBase3D<T,Descriptor>::swap(rhs);
    AtomicBlock3D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
}

/// For an AtomicBlock, the lower left corner is always at the origin
template<typename T, template<typename U> class Descriptor>
Box3D BlockLattice3D<T,Descriptor>::getBoundingBox() const {
    return Box3D(0, nx-1, 0, ny-1, 0, nz-1);
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::initialize() {
    this->executeInternalProcessors();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::specifyStatisticsStatus(Box3D domain, bool status) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collide(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                grid[iX][iY][iZ].collide(this->getInternalStatistics());
                grid[iX][iY][iZ].revert();
            }
        }
    }
}

/** \sa collide(int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collide() {
    collide(this->getBoundingBox());
}

/** The distribution function never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should 
 * be streamed outside are simply left untouched.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::stream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    bulkStream(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                     domain.y0+vicinity,domain.y1-vicinity,
                     domain.z0+vicinity,domain.z1-vicinity) );

    boundaryStream(domain, Box3D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x1-vicinity-1,domain.x1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z0,domain.z0+vicinity-1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z1-vicinity+1,domain.z1) );
}

/** At the end of this method, the methods finalizeIteration()
 * and executeInternalProcessors() are automatically invoked.
 * \sa stream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::stream()
{
    stream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
 * because memory is traversed only once instead of twice.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D3Q19)
    collide(Box3D(domain.x0,domain.x0+vicinity-1,
                  domain.y0,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x1-vicinity+1,domain.x1,
                  domain.y0,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0,domain.y0+vicinity-1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y1-vicinity+1,domain.y1,
                  domain.z0,domain.z1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0+vicinity,domain.y1-vicinity,
                  domain.z0,domain.z0+vicinity-1) );
    collide(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                  domain.y0+vicinity,domain.y1-vicinity,
                  domain.z1-vicinity+1,domain.z1) );

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(Box3D(domain.x0+vicinity,domain.x1-vicinity,
                               domain.y0+vicinity,domain.y1-vicinity,
                               domain.z0+vicinity,domain.z1-vicinity) );

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(domain, Box3D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1,
                                 domain.z0,domain.z1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z0,domain.z0+vicinity-1) );
    boundaryStream(domain, Box3D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0+vicinity,domain.y1-vicinity,
                                 domain.z1-vicinity+1,domain.z1) );
}

/** At the end of this method, finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa collideAndStream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::collideAndStream() {
    collideAndStream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::incrementTime() {
    this->getTimeCounter().incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::allocateMemory() {
    rawData = new Cell<T,Descriptor> [nx*ny*nz];
    grid    = new Cell<T,Descriptor>** [nx];
    for (plint iX=0; iX<nx; ++iX) {
        grid[iX] = new Cell<T,Descriptor>* [ny];
        for (plint iY=0; iY<ny; ++iY) {
            grid[iX][iY] = rawData + nz*(iY+ny*iX);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::releaseMemory() {
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            for (plint iZ=0; iZ<nz; ++iZ) {
                Dynamics<T,Descriptor>* dynamics = &grid[iX][iY][iZ].getDynamics();
                if (dynamics != backgroundDynamics) {
                    delete dynamics;
                }
            }
        }
    }
    delete backgroundDynamics;
    delete [] rawData;
    for (plint iX=0; iX<nx; ++iX) {
        delete [] grid[iX];
    }
    delete [] grid;
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::attributeDynamics (
        plint iX, plint iY, plint iZ, Dynamics<T,Descriptor>* dynamics )
{
    Dynamics<T,Descriptor>* previousDynamics = &grid[iX][iY][iZ].getDynamics();
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX][iY][iZ].attributeDynamics(dynamics);
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& BlockLattice3D<T,Descriptor>::getBackgroundDynamics() {
    return *backgroundDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& BlockLattice3D<T,Descriptor>::getBackgroundDynamics() const {
    return *backgroundDynamics;
}

/** This method is slower than bulkStream(int,int,int,int), because it must
 * be verified which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::boundaryStream(Box3D bound, Box3D domain) {
    // Make sure bound is contained within current lattice
    PLB_PRECONDITION( contained(bound, this->getBoundingBox()) );
    // Make sure domain is contained within bound
    PLB_PRECONDITION( contained(domain, bound) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if ( nextX>=bound.x0 && nextX<=bound.x1 &&
                         nextY>=bound.y0 && nextY<=bound.y1 &&
                         nextZ>=bound.z0 && nextZ<=bound.z1 )
                    {
                        std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                                  grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

/** This method is faster than boundaryStream(int,int,int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::bulkStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    std::swap(grid[iX][iY][iZ][iPop+Descriptor<T>::q/2],
                              grid[nextX][nextY][nextZ][iPop]);
                }
            }
        }
    }
}

/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::bulkCollideAndStream(Box3D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    // For cache efficiency, memory is traversed block-wise. The three outer loops enumerate
    //   the blocks, whereas the three inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    // Outer loops.
    for (plint outerX=domain.x0; outerX<=domain.x1; outerX+=blockSize) {
        for (plint outerY=domain.y0; outerY<=domain.y1+blockSize-1; outerY+=blockSize) {
            for (plint outerZ=domain.z0; outerZ<=domain.z1+2*(blockSize-1); outerZ+=blockSize) {
                // Inner loops.
                plint dx = 0;
                for (plint innerX=outerX;
                     innerX <= min(outerX+blockSize-1, domain.x1);
                     ++innerX, ++dx)
                {
                    // Y-index is shifted in negative direction at each x-increment. to ensure
                    //   that only post-collision cells are accessed during the swap-operation
                    //   of the streaming.
                    plint minY = outerY-dx;
                    plint maxY = minY+blockSize-1;
                    plint dy = 0;
                    for (plint innerY=max(minY,domain.y0);
                         innerY <= min(maxY, domain.y1);
                         ++innerY, ++dy)
                    {
                        // Z-index is shifted in negative direction at each x-increment. and at each
                        //    y-increment, to ensure that only post-collision cells are accessed during
                        //    the swap-operation of the streaming.
                        plint minZ = outerZ-dx-dy;
                        plint maxZ = minZ+blockSize-1;
                        for (plint innerZ=max(minZ,domain.z0);
                             innerZ <= min(maxZ, domain.z1);
                             ++innerZ)
                        {
                            // Collide the cell.
                            grid[innerX][innerY][innerZ].collide (
                                    this->getInternalStatistics() );
                            // Swap the populations on the cell, and then with post-collision
                            //   neighboring cell, to perform the streaming step.
                            latticeTemplates<T,Descriptor>::swapAndStream3D (
                                    grid, innerX, innerY, innerZ );
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::implementPeriodicity() {
    static const plint vicinity = Descriptor<T>::vicinity;
    plint maxX = nx-1;
    plint maxY = ny-1;
    plint maxZ = nz-1;
    bool periodicX = this->periodicity().get(0);
    bool periodicY = this->periodicity().get(1);
    bool periodicZ = this->periodicity().get(2);

    // Periodicity of planes orthogonal to x-axis.
    if (periodicX) {
        periodicDomain(Box3D(-vicinity,-1, 0,maxY, 0,maxZ));
    }
    // Periodicity of planes orthogonal to y-axis.
    if (periodicY) {
        periodicDomain(Box3D(0,maxX, -vicinity,-1, 0,maxZ));
    }
    // Periodicity of planes orthogonal to z-axis.
    if (periodicZ) {
        periodicDomain(Box3D(0,maxX, 0,maxY, -vicinity,-1));
    }

    // Periodicity of edges in y-z plane.
    if (periodicY && periodicZ) {
        // Periodicity between (-y,-z) and (+y,+z) edge.
        periodicDomain(Box3D(0,maxX, -vicinity,-1, -vicinity,-1));
        // Periodicity between (-y,+z) and (+y,-z) edge.
        periodicDomain(Box3D(0,maxX, -vicinity,-1, maxZ+1,maxZ+vicinity));
    }
    // Periodicity of edges in x-z plane.
    if (periodicX && periodicZ) {
        // Periodicity between (-x,-z) and (+x,+z) edge.
        periodicDomain(Box3D(-vicinity,-1, 0,maxY, -vicinity,-1));
        // Periodicity between (-x,+z) and (+x,-z) edge.
        periodicDomain(Box3D(maxX+1,maxX+vicinity, 0,maxY, -vicinity,-1));
    }
    // Periodicity of edges in x-y plane.
    if (periodicX && periodicY) {
        // Periodicity between (-x,-y) and (+x,+y) edge.
        periodicDomain(Box3D(-vicinity,-1, -vicinity,-1, 0,maxZ));
        // Periodicity between (-x,+y) and (+x,-y) edge.
        periodicDomain(Box3D(-vicinity,-1, maxY+1,maxY+vicinity, 0,maxZ));
    }

    // Periodicity of corners.
    if (periodicX && periodicY && periodicZ) {
        // Periodicity of (+1,+1,+1) and (-1,-1,-1)
        periodicDomain(Box3D(maxX+1,maxX+vicinity, maxY+1,maxY+vicinity, maxZ+1,maxZ+vicinity));
        // Periodicity of (+1,+1,-1) and (-1,-1,+1)
        periodicDomain(Box3D(maxX+1,maxX+vicinity, maxY+1,maxY+vicinity, -vicinity,-1));
        // Periodicity of (+1,-1,+1) and (-1,+1,-1)
        periodicDomain(Box3D(maxX+1,maxX+vicinity, -vicinity,-1, maxZ+1,maxZ+vicinity));
        // Periodicity of (+1,-1,-1) and (-1,+1,+1)
        periodicDomain(Box3D(maxX+1,maxX+vicinity, -vicinity,-1, -vicinity,-1));
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice3D<T,Descriptor>::periodicDomain(Box3D domain) {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                    plint prevX = iX - Descriptor<T>::c[iPop][0];
                    plint prevY = iY - Descriptor<T>::c[iPop][1];
                    plint prevZ = iZ - Descriptor<T>::c[iPop][2];

                    if ( (prevX>=0 && prevX<nx) &&
                         (prevY>=0 && prevY<ny) &&
                         (prevZ>=0 && prevZ<nz) )
                    {
                        plint nextX = (iX+nx)%nx;
                        plint nextY = (iY+ny)%ny;
                        plint nextZ = (iZ+nz)%nz;
                        std::swap (
                            grid[prevX][prevY][prevZ][indexTemplates::opposite<Descriptor<T> >(iPop)],
                            grid[nextX][nextY][nextZ][iPop] );
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
identifiers::BlockId BlockLattice3D<T,Descriptor>::getBlockId() const {
    return identifiers::getLatticeId<T,Descriptor>();
}

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer3D<T,Descriptor>& BlockLattice3D<T,Descriptor>::getDataTransfer() {
    return dataTransfer;
}

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer3D<T,Descriptor> const& BlockLattice3D<T,Descriptor>::getDataTransfer() const {
    return dataTransfer;
}


////////////////////// Class BlockLatticeDataTransfer3D /////////////////////////

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer3D<T,Descriptor>::BlockLatticeDataTransfer3D(BlockLattice3D<T,Descriptor>& lattice_)
    : lattice(lattice_)
{ }

template<typename T, template<typename U> class Descriptor>
plint BlockLatticeDataTransfer3D<T,Descriptor>::sizeOfCell() const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::send(Box3D domain, T* buffer) const {
    PLB_PRECONDITION(contained(domain, lattice.getBoundingBox()));
    plint cellSize = sizeOfCell();
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).serialize(buffer+iData);
                iData += cellSize;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::receive(Box3D domain, T const* buffer) {
    PLB_PRECONDITION(contained(domain, lattice.getBoundingBox()));
    plint cellSize = sizeOfCell();
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).unSerialize(buffer+iData);
                iData += cellSize;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer3D<T,Descriptor>::attribute (
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D<T> const& from)
{
    PLB_PRECONDITION (typeid(from) == typeid(BlockLattice3D<T,Descriptor> const&));
    PLB_PRECONDITION(contained(toDomain, lattice.getBoundingBox()));
    BlockLattice3D<T,Descriptor> const& fromLattice = (BlockLattice3D<T,Descriptor> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).attributeValues (
                        fromLattice.get(iX+deltaX,iY+deltaY,iZ+deltaZ) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CachePolicy3D& BlockLattice3D<T,Descriptor>::cachePolicy() {
    static CachePolicy3D cachePolicySingleton(30);
    return cachePolicySingleton;
}

}  // namespace plb

#endif
