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
 * Helper classes for parallel 2D multiblock lattice -- header file.
 */
#ifndef PARALLEL_MULTI_BLOCK_LATTICE_2D_H
#define PARALLEL_MULTI_BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "parallelism/parallelBlockCommunicator2D.h"
#include "multiBlock/multiBlockLattice2D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ParallelCellAccess2D : public MultiCellAccess2D<T,Descriptor> {
public:
    ParallelCellAccess2D();
    ~ParallelCellAccess2D();
    virtual Cell<T,Descriptor>& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<BlockLattice2D<T,Descriptor>*>& lattices );
    virtual Cell<T,Descriptor> const& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::vector<BlockLattice2D<T,Descriptor>*> const& lattices ) const;
    virtual void broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                               MultiBlockManagement2D const& multiBlockManagement) const;
    ParallelCellAccess2D<T,Descriptor>* clone() const;
private:
    mutable Cell<T,Descriptor> distributedCell;
    mutable std::vector<Cell<T,Descriptor>*> baseCells;
    mutable std::vector<Cell<T,Descriptor> const*> constBaseCells;
    mutable Dynamics<T,Descriptor>* parallelDynamics;
};

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif
