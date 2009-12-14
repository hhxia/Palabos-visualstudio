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
#ifndef DEFAULT_MULTI_BLOCK_POLICY_3D_H
#define DEFAULT_MULTI_BLOCK_POLICY_3D_H

#include "core/globalDefs.h"
#include "multiBlock/serialBlockCommunicator3D.h"
#include "parallelism/parallelBlockCommunicator3D.h"
#include "multiBlock/combinedStatistics.h"
#include "parallelism/parallelStatistics.h"
#include "multiBlock/serialMultiBlockLattice3D.h"
#include "parallelism/parallelMultiBlockLattice3D.h"
#include "multiBlock/serialMultiDataField3D.h"
#include "parallelism/parallelMultiDataField3D.h"
#include "multiBlock/threadAttribution.h"
#include "multiBlock/staticRepartitions3D.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

class DefaultMultiBlockPolicy3D {
public:
    template<typename T>
    BlockCommunicator3D<T>* getBlockCommunicator() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelBlockCommunicator3D<T>();
#else
        return new SerialBlockCommunicator3D<T>();
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
    MultiCellAccess3D<T,Descriptor>* getMultiCellAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelCellAccess3D<T,Descriptor>();
#else
        return new SerialCellAccess3D<T,Descriptor>();
#endif
    }

    template<typename T>
    MultiScalarAccess3D<T>* getMultiScalarAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelScalarAccess3D<T>();
#else
        return new SerialScalarAccess3D<T>();
#endif
    }

    template<typename T, int nDim>
    MultiTensorAccess3D<T,nDim>* getMultiTensorAccess() {
#ifdef PLB_MPI_PARALLEL
        return new ParallelTensorAccess3D<T,nDim>();
#else
        return new SerialTensorAccess3D<T,nDim>();
#endif
    }

    ThreadAttribution* getThreadAttribution() {
#ifdef PLB_MPI_PARALLEL
        return new OneToOneThreadAttribution();
#else
        return new SerialThreadAttribution();
#endif
    }

    MultiBlockManagement3D getMultiBlockManagement(plint nx, plint ny, plint nz) {
        plint envelopeWidth = 1;
        return MultiBlockManagement3D( createRegularMultiBlockDistribution3D(nx,ny,nz, envelopeWidth, numProcesses),
                                       getThreadAttribution() );
    }

    void setNumProcesses(int numProcesses_) {
        numProcesses = numProcesses_;
    }

    int getNumProcesses() const {
        return numProcesses;
    }
private:
    DefaultMultiBlockPolicy3D()
        : numProcesses(global::mpi().getSize())
    { }
    friend DefaultMultiBlockPolicy3D& defaultMultiBlockPolicy3D();
private:
    int numProcesses;
};

inline DefaultMultiBlockPolicy3D& defaultMultiBlockPolicy3D() {
    static DefaultMultiBlockPolicy3D singleton;
    return singleton;
}

}  // namespace plb

#endif
