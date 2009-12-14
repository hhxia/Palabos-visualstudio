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
 * Operations on the 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_OPERATIONS_2D_H
#define MULTI_BLOCK_OPERATIONS_2D_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include <vector>

namespace plb {

template<typename T> class MultiBlock2D;
template<typename T> class DataProcessorGenerator2D;
template<typename T> class ReductiveDataProcessorGenerator2D;

template<typename T, class OriginalGenerator, class MutableGenerator>
class MultiProcessing2D {
public:
    MultiProcessing2D( OriginalGenerator& generator_,
                       std::vector<MultiBlock2D<T>*> multiBlocks_ );
    ~MultiProcessing2D();
    void extractProcessorsOnFirstBlock(BlockDomain::DomainT appliesTo);
    void intersectWithRemainingBlocks(BlockDomain::DomainT appliesTo);
    void subdivideGenerator();
    void adjustCoordinates();
    std::vector<MutableGenerator*> const& getRetainedGenerators() const;
    std::vector<std::vector<plint> > const& getAtomicBlockNumbers() const;
    std::vector<MultiBlock2D<T>*> multiBlocksWhichRequireUpdate() const;
    void updateEnvelopesWhereRequired();
private:
    void extractGeneratorOnBlocks(std::vector<Box2D> const& finalDomains,
                                  std::vector<std::vector<plint> > const& finalIds,
                                  plint shiftX=0, plint shiftY=0);
private:
    OriginalGenerator&              generator;
    std::vector<MultiBlock2D<T>*>   multiBlocks;
    MultiBlock2D<T>*                firstMultiBlock;
    std::vector<MutableGenerator*>  retainedGenerators;
    std::vector<std::vector<plint> > atomicBlockNumbers;
};

template<typename T>
void executeDataProcessor( DataProcessorGenerator2D<T> const& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks );

template<typename T>
void executeDataProcessor( ReductiveDataProcessorGenerator2D<T>& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks );

template<typename T>
void addInternalProcessor( DataProcessorGenerator2D<T> const& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks, plint level=0 );

} // namespace plb

#endif
