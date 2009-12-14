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
 * Default policy for instantiating serial/parallel multi-blocks
 */
#ifndef DEFAULT_MULTI_BLOCK_POLICY_2D_H
#define DEFAULT_MULTI_BLOCK_POLICY_2D_H

#include "core/globalDefs.h"
#include "multiBlock/serialBlockCommunicator2D.h"
#include "parallelism/parallelBlockCommunicator2D.h"
#include "multiBlock/combinedStatistics.h"
#include "parallelism/parallelStatistics.h"
#include "multiBlock/serialMultiBlockLattice2D.h"
#include "parallelism/parallelMultiBlockLattice2D.h"
#include "multiBlock/serialMultiDataField2D.h"
#include "parallelism/parallelMultiDataField2D.h"
#include "multiBlock/threadAttribution.h"
#include "multiBlock/staticRepartitions2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

class DefaultMultiBlockPolicy2D {
public:
    template<typename T>
    BlockCommunicator2D<T>* getBlockCommunicator() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelBlockCommunicator2D<T>();
#else
        return new SerialBlockCommunicator2D<T>();
#endif
    }

    template<typename T>
    CombinedStatistics<T>* getCombinedStatistics() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelCombinedStatistics<T>();
#else
        return new SerialCombinedStatistics<T>();
#endif
    }

    template<typename T, template<typename U> class Descriptor>
    MultiCellAccess2D<T,Descriptor>* getMultiCellAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelCellAccess2D<T,Descriptor>();
#else
        return new SerialCellAccess2D<T,Descriptor>();
#endif
    }

    template<typename T>
    MultiScalarAccess2D<T>* getMultiScalarAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelScalarAccess2D<T>();
#else
        return new SerialScalarAccess2D<T>();
#endif
    }

    template<typename T, int nDim>
    MultiTensorAccess2D<T,nDim>* getMultiTensorAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelTensorAccess2D<T,nDim>();
#else
        return new SerialTensorAccess2D<T,nDim>();
#endif
    }

    ThreadAttribution* getThreadAttribution() {
#ifdef PLB_MPI_PARALLEL
        return new OneToOneThreadAttribution();
#else
        return new SerialThreadAttribution();
#endif
    }

    MultiBlockManagement2D getMultiBlockManagement(plint nx, plint ny) {
        plint envelopeWidth=1;
        return MultiBlockManagement2D( createRegularMultiBlockDistribution2D(nx,ny, envelopeWidth, numProcesses),
                                       getThreadAttribution() );
    }

    void setNumProcesses(int numProcesses_) {
        numProcesses = numProcesses_;
    }

    int getNumProcesses() const {
        return numProcesses;
    }
private:
    DefaultMultiBlockPolicy2D()
        : numProcesses(global::mpi().getSize())
    { }
    friend DefaultMultiBlockPolicy2D& defaultMultiBlockPolicy2D();
private:
    int numProcesses;
};


inline DefaultMultiBlockPolicy2D& defaultMultiBlockPolicy2D() {
    static DefaultMultiBlockPolicy2D singleton;
    return singleton;
}

}  // namespace plb

#endif
