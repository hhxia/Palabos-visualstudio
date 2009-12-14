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
 * The dynamics of a 2D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_2D_HH
#define BLOCK_LATTICE_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "core/dynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/latticeTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "core/util.h"
#include <algorithm>
#include <typeinfo>
using namespace std;

namespace plb {

////////////////////// Class BlockLattice2D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice2D<T,Descriptor>::BlockLattice2D (
        plint nx_, plint ny_,
        Dynamics<T,Descriptor>* backgroundDynamics_ )
    : nx(nx_), ny(ny_),
      backgroundDynamics(backgroundDynamics_),
      dataTransfer(*this)
{
    // Allocate memory and attribute dynamics.
    allocateMemory();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            grid[iX][iY].attributeDynamics(backgroundDynamics);
        }
    }
    // Attribute default value to the standard statistics (average uSqr,
    //   max uSqr, average rho). These have previously been subscribed
    //   in the constructor of BlockLatticeBase2D.
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
BlockLattice2D<T,Descriptor>::~BlockLattice2D()
{
    releaseMemory();
}

/** The whole data of the lattice is duplicated. This includes
 * both particle distribution function and external fields.
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice2D<T,Descriptor>::BlockLattice2D(BlockLattice2D<T,Descriptor> const& rhs)
    : BlockLatticeBase2D<T,Descriptor>(rhs),
      AtomicBlock2D<T>(rhs),
      nx(rhs.nx),
      ny(rhs.ny),
      backgroundDynamics(rhs.backgroundDynamics->clone()),
      dataTransfer(*this)
{
    allocateMemory();
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            Cell<T,Descriptor>& cell = grid[iX][iY];
            // Assign cell from rhs
            cell = rhs.grid[iX][iY];
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

/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated. This includes both particle distribution function
 * and external fields. 
 * \warning The dynamics objects and internalProcessors are not copied
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Descriptor>
BlockLattice2D<T,Descriptor>& BlockLattice2D<T,Descriptor>::operator= (
        BlockLattice2D<T,Descriptor> const& rhs )
{
    BlockLattice2D<T,Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the 
 * lattice are copied, and not the lattice itself.
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::swap(BlockLattice2D& rhs) {
    BlockLatticeBase2D<T,Descriptor>::swap(rhs);
    AtomicBlock2D<T>::swap(rhs);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(backgroundDynamics, rhs.backgroundDynamics);
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
}

/// For an AtomicBlock, the lower left corner is always at the origin
template<typename T, template<typename U> class Descriptor>
Box2D BlockLattice2D<T,Descriptor>::getBoundingBox() const {
    return Box2D(0, nx-1, 0, ny-1);
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::initialize() {
    this->executeInternalProcessors();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::specifyStatisticsStatus (
        Box2D domain, bool status )
{
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            grid[iX][iY].specifyStatisticsStatus(status);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::collide(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            grid[iX][iY].collide(this->getInternalStatistics());
            grid[iX][iY].revert();
        }
    }
}

/** \sa collide(int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::collide() {
    collide(this->getBoundingBox());
}

/** The distribution functions never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should
 * be streamed outside are simply left untouched.
 * The finalization of an iteration step is not automatically executed,
 * as it is in the method stream(). If you want it to be executed, you
 * must explicitly call the methods finalizeIteration() and
 * executeInternalProcessors().
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::stream(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    bulkStream( Box2D(domain.x0+vicinity,domain.x1-vicinity,
                      domain.y0+vicinity,domain.y1-vicinity) );

    boundaryStream(domain, Box2D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1));
}

/** At the end of this method, the methods finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa stream(int,int,int,int)
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::stream()
{
    stream(this->getBoundingBox());

    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int) and stream(int,int,int,int), because memory
 * is traversed only once instead of twice.
 * The finalization of an iteration step is not automatically invoked by this
 * method, as it is in the method stream(). If you want it to be executed, you
 * must explicitly call the methods finalizeIteration() and
 * executeInternalProcessors().
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::collideAndStream(Box2D domain)
{
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    static const plint vicinity = Descriptor<T>::vicinity;

    // First, do the collision on cells within a boundary envelope of width
    // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
    collide(Box2D(domain.x0,domain.x0+vicinity-1, domain.y0,domain.y1));
    collide(Box2D(domain.x1-vicinity+1,domain.x1, domain.y0,domain.y1));
    collide(Box2D(domain.x0+vicinity,domain.x1-vicinity, domain.y0,domain.y0+vicinity-1));
    collide(Box2D(domain.x0+vicinity,domain.x1-vicinity, domain.y1-vicinity+1,domain.y1));

    // Then, do the efficient collideAndStream algorithm in the bulk,
    // excluding the envelope (this is efficient because there is no
    // if-then-else statement within the loop, given that the boundary
    // region is excluded)
    bulkCollideAndStream(Box2D(domain.x0+vicinity,domain.x1-vicinity,
                               domain.y0+vicinity,domain.y1-vicinity));

    // Finally, do streaming in the boundary envelope to conclude the
    // collision-stream cycle
    boundaryStream(domain, Box2D(domain.x0,domain.x0+vicinity-1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x1-vicinity+1,domain.x1,
                                 domain.y0,domain.y1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y0,domain.y0+vicinity-1));
    boundaryStream(domain, Box2D(domain.x0+vicinity,domain.x1-vicinity,
                                 domain.y1-vicinity+1,domain.y1));
}

/** At the end of this method, the methods finalizeIteration() and
 * executeInternalProcessors() are automatically invoked.
 * \sa collideAndStream(int,int,int,int) */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::collideAndStream() {
    collideAndStream(this->getBoundingBox());
    
    implementPeriodicity();

    this->executeInternalProcessors();
    this->evaluateStatistics();
    this->incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::incrementTime() {
    this->getTimeCounter().incrementTime();
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::allocateMemory() {
    rawData = new Cell<T,Descriptor> [nx*ny];
    grid    = new Cell<T,Descriptor>* [nx];
    for (plint iX=0; iX<nx; ++iX) {
        grid[iX] = rawData + iX*ny;
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::releaseMemory() {
    for (plint iX=0; iX<nx; ++iX) {
        for (plint iY=0; iY<ny; ++iY) {
            Dynamics<T,Descriptor>* dynamics = &grid[iX][iY].getDynamics();
            if (dynamics != backgroundDynamics) {
                delete dynamics;
            }
        }
    }
    delete backgroundDynamics;
    delete [] rawData;
    delete [] grid;
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::attributeDynamics(plint iX, plint iY, Dynamics<T,Descriptor>* dynamics) {
    Dynamics<T,Descriptor>* previousDynamics = &grid[iX][iY].getDynamics();
    if (previousDynamics != backgroundDynamics) {
        delete previousDynamics;
    }
    grid[iX][iY].attributeDynamics(dynamics);
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& BlockLattice2D<T,Descriptor>::getBackgroundDynamics() {
    return *backgroundDynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& BlockLattice2D<T,Descriptor>::getBackgroundDynamics() const {
    return *backgroundDynamics;
}

/** This method is slower than bulkStream(int,int,int,int), because one needs
 * to verify which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::boundaryStream(Box2D bound, Box2D domain) {
    // Make sure bound is contained within current lattice
    PLB_PRECONDITION( contained(bound, this->getBoundingBox()) );
    // Make sure domain is contained within bound
    PLB_PRECONDITION( contained(domain, bound) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                plint nextX = iX + Descriptor<T>::c[iPop][0];
                plint nextY = iY + Descriptor<T>::c[iPop][1];
                if (nextX>=bound.x0 && nextX<=bound.x1 && nextY>=bound.y0 && nextY<=bound.y1) {
                    std::swap(grid[iX][iY][iPop+Descriptor<T>::q/2],
                              grid[nextX][nextY][iPop]);
                }
            }
        }
    }
}

/** This method is faster than boundaryStream(int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::bulkStream(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
                plint nextX = iX + Descriptor<T>::c[iPop][0];
                plint nextY = iY + Descriptor<T>::c[iPop][1];
                std::swap(grid[iX][iY][iPop+Descriptor<T>::q/2],
                          grid[nextX][nextY][iPop]);
            }
        }
    }
}

/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int)
 * \sa collideAndStream()
 */
/*
 *template<typename T, template<typename U> class Descriptor>
 *void BlockLattice2D<T,Descriptor>::bulkCollideAndStream(Box2D domain) {
 *    // Make sure domain is contained within current lattice
 *    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );
 *
 *    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
 *        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
 *            grid[iX][iY].collide(this->getInternalStatistics());
 *            latticeTemplates<T,Descriptor>::swapAndStream2D(grid, iX, iY);
 *        }
 *    }
 *}
 */

/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::bulkCollideAndStream(Box2D domain) {
    // Make sure domain is contained within current lattice
    PLB_PRECONDITION( contained(domain, this->getBoundingBox()) );

    // For cache efficiency, memory is traversed block-wise. The two outer loops enumerate
    //   the blocks, whereas the two inner loops enumerate the cells inside each block.
    const plint blockSize = cachePolicy().getBlockSize();
    // Outer loops.
    for (plint outerX=domain.x0; outerX<=domain.x1; outerX+=blockSize) {
        for (plint outerY=domain.y0; outerY<=domain.y1+blockSize-1; outerY+=blockSize) {
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
                for (plint innerY=max(minY,domain.y0);
                     innerY <= min(maxY, domain.y1);
                     ++innerY)
                {
                    // Collide the cell.
                    grid[innerX][innerY].collide (
                            this->getInternalStatistics() );
                    // Swap the populations on the cell, and then with post-collision
                    //   neighboring cell, to perform the streaming step.
                    latticeTemplates<T,Descriptor>::swapAndStream2D (
                            grid, innerX, innerY );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::periodicDomain(Box2D domain) {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                plint prevX = iX - Descriptor<T>::c[iPop][0];
                plint prevY = iY - Descriptor<T>::c[iPop][1];
                if ( (prevX>=0 && prevX<nx) &&
                     (prevY>=0 && prevY<ny) )
                {
                    plint nextX = (iX+nx)%nx;
                    plint nextY = (iY+ny)%ny;
                    std::swap (
                        grid[prevX][prevY][indexTemplates::opposite<Descriptor<T> >(iPop)],
                        grid[nextX][nextY][iPop] );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLattice2D<T,Descriptor>::implementPeriodicity() {
    static const plint vicinity = Descriptor<T>::vicinity;
    plint maxX = nx-1;
    plint maxY = ny-1;
    bool periodicX = this->periodicity().get(0);
    bool periodicY = this->periodicity().get(1);
    if (periodicX) {  // Periodicity of edges in x-direction.
        periodicDomain(Box2D(-vicinity,-1,0,maxY));
    }
    if (periodicY) {  // Periodicity of edges in y-direction.
        periodicDomain(Box2D(0,maxX,-vicinity,-1));
    }
    if (periodicX && periodicY) {
        // Periodicity between (-1,-1) and (+1,+1) corner.
        periodicDomain(Box2D(-vicinity,-1,-vicinity,-1));
        // Periodicity between (-1,+1) and (+1,-1) corner.
        periodicDomain(Box2D(-vicinity,-1,maxY+1,maxY+vicinity));
    }
}

template<typename T, template<typename U> class Descriptor>
identifiers::BlockId BlockLattice2D<T,Descriptor>::getBlockId() const {
    return identifiers::getLatticeId<T,Descriptor>();
}

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer2D<T,Descriptor>& BlockLattice2D<T,Descriptor>::getDataTransfer() {
    return dataTransfer;
}

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer2D<T,Descriptor> const& BlockLattice2D<T,Descriptor>::getDataTransfer() const {
    return dataTransfer;
}


////////////////////// Class BlockLatticeDataTransfer2D /////////////////////////

template<typename T, template<typename U> class Descriptor>
BlockLatticeDataTransfer2D<T,Descriptor>::BlockLatticeDataTransfer2D(BlockLattice2D<T,Descriptor>& lattice_)
    : lattice(lattice_)
{ }

template<typename T, template<typename U> class Descriptor>
plint BlockLatticeDataTransfer2D<T,Descriptor>::sizeOfCell() const {
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer2D<T,Descriptor>::send(Box2D domain, T* buffer) const {
    PLB_PRECONDITION(contained(domain, lattice.getBoundingBox()));
    plint cellSize = sizeOfCell();
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).serialize(buffer+iData);
            iData += cellSize;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer2D<T,Descriptor>::receive(Box2D domain, T const* buffer) {
    PLB_PRECONDITION(contained(domain, lattice.getBoundingBox()));
    plint cellSize = sizeOfCell();
    plint iData=0;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).unSerialize(buffer+iData);
            iData += cellSize;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void BlockLatticeDataTransfer2D<T,Descriptor>::attribute (
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D<T> const& from )
{
    PLB_PRECONDITION (typeid(from) == typeid(BlockLattice2D<T,Descriptor> const&));
    PLB_PRECONDITION(contained(toDomain, lattice.getBoundingBox()));
    BlockLattice2D<T,Descriptor> const& fromLattice = (BlockLattice2D<T,Descriptor> const&) from;
    for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
        for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
            lattice.get(iX,iY).attributeValues(fromLattice.get(iX+deltaX,iY+deltaY));
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CachePolicy2D& BlockLattice2D<T,Descriptor>::cachePolicy() {
    static CachePolicy2D cachePolicySingleton(200);
    return cachePolicySingleton;
}

}  // namespace plb

#endif
