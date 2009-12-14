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
 * Operations on the 2D multiblock -- generic implementation.
 */
#ifndef MULTI_BLOCK_OPERATIONS_2D_HH
#define MULTI_BLOCK_OPERATIONS_2D_HH

#include "multiBlock/multiBlockOperations2D.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/domainManipulation2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/atomicBlockOperations2D.h"
#include "multiGrid/multiScale.h"
#include "core/plbDebug.h"


namespace plb {

/* *************** Class MultiProcessing2D *************************** */

template<typename T, class OriginalGenerator, class MutableGenerator>
MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::MultiProcessing2D (
        OriginalGenerator& generator_,
        std::vector<MultiBlock2D<T>*> multiBlocks_ )
    : generator(generator_),
      multiBlocks(multiBlocks_)
{
    PLB_PRECONDITION( multiBlocks.size()>=1 );
    firstMultiBlock = multiBlocks[0];

    // Subdivide the original generator into smaller generators which act
    //   on the intersection of all implied blocks. At this stage, all coordinates
    //   are global, thus relative to the multi-block, not to the individual
    //   atomic-blocks.
    subdivideGenerator();

    // Then, convert coordinates from a global representation to a local one,
    //   relative to the atomoic-blocks of the first multi-block.
    adjustCoordinates();
}

template<typename T, class OriginalGenerator, class MutableGenerator>
void MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::subdivideGenerator()
{
    // To start with, determine which multi-blocks are read and which are written
    std::vector<bool> isWritten(multiBlocks.size());
    generator.getModificationPattern(isWritten);
    PLB_ASSERT( isWritten.size() == multiBlocks.size() );

    // The reference block (the one for which the envelope is included if
    //   the domain generator.appliesTo() include the envelope) is either the
    //   multi-block which is written, or the first multi-block if all are read-only.
    pluint referenceBlock = 0;
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        if (isWritten[iBlock]) {
            referenceBlock = iBlock;
            break;
        }
    }

    // In debug mode, make sure that a most one multi-block is written when envelope is included.
#ifdef PLB_DEBUG
    if ( BlockDomain::usesEnvelope(generator.appliesTo()) ) {
        plint numWritten = 0;
        for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
            if (isWritten[iBlock]) {
                ++numWritten;
            }
        }
        PLB_ASSERT( numWritten <= 1 );
    }
#endif
    
    // The first step is to access the domains of the the atomic blocks, as well
    //   as their IDs in each of the coupled multi blocks. The domain corresponds
    //   to the bulk and/or to the envelope, depending on the value of generator.appliesTo().
    std::vector<std::vector<DomainAndId2D> > domainsWithId(multiBlocks.size());
    for (pluint iMulti=0; iMulti<multiBlocks.size(); ++iMulti) {
        std::vector<plint> const& relevantBlocks
            = multiBlocks[iMulti] -> getMultiBlockManagement().getRelevantIndexes().getBlocks();
        for (pluint rBlock=0; rBlock < relevantBlocks.size(); ++rBlock) {
            plint iBlock = relevantBlocks[rBlock];
            BlockParameters2D const& params = multiBlocks[iMulti] -> getParameters(iBlock);
            switch (generator.appliesTo()) {
                case BlockDomain::bulk:
                    domainsWithId[iMulti].push_back(DomainAndId2D(params.getBulk(),iBlock));
                    break;
                case BlockDomain::bulkAndEnvelope:
                    domainsWithId[iMulti].push_back(DomainAndId2D(params.getEnvelope(),iBlock));
                    break;
                case BlockDomain::envelope:
                    // For the reference block, we restrict ourselves to the envelope, because
                    //   that's the desired domain of application.
                    if (iMulti==referenceBlock) {
                        std::vector<Box2D> envelopeOnly;
                        except(params.getEnvelope(), params.getBulk(), envelopeOnly);
                        for (pluint iEnvelope=0; iEnvelope<envelopeOnly.size(); ++iEnvelope) {
                            domainsWithId[iMulti].push_back(DomainAndId2D(envelopeOnly[iEnvelope], iBlock));
                        }
                    }
                    // For the other blocks, we need to take bulk and envelope, because all these domains
                    //   potentially intersect with the envelope of the reference block.
                    else {
                        domainsWithId[iMulti].push_back(DomainAndId2D(params.getEnvelope(),iBlock));
                    }
                    break;
            }
        }
    }

    // If the multi-blocks are not at the same level of grid refinement, the level
    //   of the first block is taken as reference, and the coordinates of the other
    //   blocks are rescaled accordingly.
    plint firstLevel = multiBlocks[0]->getMultiBlockManagement().getRefinementLevel();
    for (pluint iMulti=1; iMulti<multiBlocks.size(); ++iMulti) {
        plint relativeLevel = firstLevel -
                             multiBlocks[iMulti]->getMultiBlockManagement().getRefinementLevel();
        if (relativeLevel != 0) {
            for (pluint iBlock=0; iBlock<domainsWithId[iMulti].size(); ++iBlock) {
                domainsWithId[iMulti][iBlock].domain =
                    global::getDefaultMultiScaleManager<T>().scaleBox (
                            domainsWithId[iMulti][iBlock].domain, relativeLevel );
            }
        }
    }

    // If the envelopes are included as well, it is assumed that at most one of
    //   the multi blocks has write-access. All others (those that have read-only
    //   access) need to be non-overlaping, to avoid multiple writes on the cells
    //   of the write-access-multi-block. Thus, overlaps are now eliminitated in
    //   the read-access-multi-blocks.
    if ( BlockDomain::usesEnvelope(generator.appliesTo()) ) {
        for (pluint iMulti=0; iMulti<multiBlocks.size(); ++iMulti) {
            if (!isWritten[iMulti]) {
                std::vector<DomainAndId2D> nonOverlapBlocks(getNonOverlapingBlocks(domainsWithId[iMulti]));
                domainsWithId[iMulti].swap(nonOverlapBlocks);
            }
        }
    }

    // This is the heart of the whole procedure: intersecting atomic blocks
    //   between all coupled multi blocks are identified.
    std::vector<Box2D> finalDomains;
    std::vector<std::vector<plint> > finalIds;
    intersectDomainsAndIds(domainsWithId, finalDomains, finalIds);

    // And, to end with, re-create processor generators adapted to the
    //   computed domains of intersection.
    if ( BlockDomain::usesEnvelope(generator.appliesTo()) ) {
        // In case the envelope is included, periodicity must be explicitly treated.
        //   Indeed, the user indicates the domain of applicability with respect to
        //   bulk nodes only. The generator is therefore shifted in all space directions
        //   to englobe periodic boundary nodes as well.
        plint shiftX = firstMultiBlock->getNx();
        plint shiftY = firstMultiBlock->getNy();
        PeriodicitySwitch2D<T> const& periodicity = firstMultiBlock->periodicity();
        for (plint orientX=-1; orientX<=+1; ++orientX) {
            for (plint orientY=-1; orientY<=+1; ++orientY) {
                if (periodicity.get(orientX,orientY)) {
                    extractGeneratorOnBlocks(finalDomains, finalIds, orientX*shiftX, orientY*shiftY);
                }
            }
        }
    }
    else {
        extractGeneratorOnBlocks(finalDomains, finalIds);
    }
}

template<typename T, class OriginalGenerator, class MutableGenerator>
void MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::extractGeneratorOnBlocks (
        std::vector<Box2D> const& finalDomains,
        std::vector<std::vector<plint> > const& finalIds,
        plint shiftX, plint shiftY)
{
    MutableGenerator* originalGenerator=generator.clone();
    originalGenerator->shift(shiftX,shiftY);
    for (pluint iDomain=0; iDomain<finalDomains.size(); ++iDomain) {
        MutableGenerator* extractedGenerator = originalGenerator->clone();
        if (extractedGenerator->extract(finalDomains[iDomain]) ) {
            retainedGenerators.push_back(extractedGenerator);
            atomicBlockNumbers.push_back(finalIds[iDomain]);
        }
        else {
            delete extractedGenerator;
        }
    }
    delete originalGenerator;
}


template<typename T, class OriginalGenerator, class MutableGenerator>
MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::~MultiProcessing2D() {
    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        delete retainedGenerators[iGenerator];
    }
}

template<typename T, class OriginalGenerator, class MutableGenerator>
void MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::adjustCoordinates() {
    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        // The generator is adjusted to local coordinates with respect to the first block.
        //   If the other blocks have a relative displacement wrt. the first block, this
        //   must be explicitly coded in the data processor.
        plint firstNumber = atomicBlockNumbers[iGenerator][0];
        BlockParameters2D const& params = firstMultiBlock->getParameters(firstNumber);
        Box2D const& envelope = params.getEnvelope();
        retainedGenerators[iGenerator]->shift(-envelope.x0, -envelope.y0);
    }
}


template<typename T, class OriginalGenerator, class MutableGenerator>
std::vector<MutableGenerator*> const&
    MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::getRetainedGenerators() const
{
    return retainedGenerators;
}

template<typename T, class OriginalGenerator, class MutableGenerator>
std::vector<std::vector<plint> > const&
    MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::getAtomicBlockNumbers() const
{
    return atomicBlockNumbers;
}

template<typename T, class OriginalGenerator, class MutableGenerator>
std::vector<MultiBlock2D<T>*> 
    MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::multiBlocksWhichRequireUpdate() const
{
    std::vector<MultiBlock2D<T>*> multiBlocksModifiedByProcessor;
    // If the generator includes envelopes, the envelopes need no update in any case.
    if ( ! BlockDomain::usesEnvelope(generator.appliesTo()) ) {
        // Otherwise, all blocks which have been modified by the processor must
        //   be updated.
        std::vector<bool> isWritten(multiBlocks.size());
        generator.getModificationPattern(isWritten);
        for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
            if (isWritten[iBlock]) {
                multiBlocksModifiedByProcessor.push_back(multiBlocks[iBlock]);
            }
        }
    }
    return multiBlocksModifiedByProcessor;
}

template<typename T, class OriginalGenerator, class MutableGenerator>
void MultiProcessing2D<T,OriginalGenerator,MutableGenerator>::updateEnvelopesWhereRequired()
{
    std::vector<MultiBlock2D<T>*> multiBlocks = multiBlocksWhichRequireUpdate();
    for (pluint iBlock=0; iBlock<multiBlocks.size(); ++iBlock) {
        multiBlocks[iBlock]->getBlockCommunicator().
            duplicateOverlaps(*multiBlocks[iBlock]);
    }
}


template<typename T>
void executeDataProcessor( DataProcessorGenerator2D<T> const& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks )
{
    MultiProcessing2D<T, DataProcessorGenerator2D<T> const, DataProcessorGenerator2D<T> >
        multiProcessing(generator, multiBlocks);
    std::vector<DataProcessorGenerator2D<T>*> const& retainedGenerators = multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const& atomicBlockNumbers = multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D<T>*> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock=0; iBlock<extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] = &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock version" of executeDataProcessor.
        plb::executeDataProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks);
    }
    // In the "executeProcessor" version, envelopes are updated right here, because the processor
    //   has already been executed. This behavior is unlike the behavior of the "addInternalProcessor" version,
    //   where envelopes are updated from within the multi-block, after the execution of internal processors.
    multiProcessing.updateEnvelopesWhereRequired();
}

template<typename T>
void executeDataProcessor( ReductiveDataProcessorGenerator2D<T>& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks )
{
    MultiProcessing2D<T, ReductiveDataProcessorGenerator2D<T>, ReductiveDataProcessorGenerator2D<T> >
        multiProcessing(generator, multiBlocks);
    std::vector<ReductiveDataProcessorGenerator2D<T>*> const& retainedGenerators = multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const& atomicBlockNumbers = multiProcessing.getAtomicBlockNumbers();

    std::vector<BlockStatistics<T> const*> individualStatistics(retainedGenerators.size());
    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D<T>*> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock=0; iBlock<extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] = &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock Reductive version" of executeDataProcessor.
        plb::executeDataProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks);
        individualStatistics[iGenerator] = &(retainedGenerators[iGenerator]->getStatistics());
    }
    multiBlocks[0]->getCombinedStatistics().combine(individualStatistics, generator.getStatistics());
    // In the "executeProcessor" version, envelopes are updated right here, because the processor
    //   has already been executed. This behavior is unlike the behavior of the "addInternalProcessor" version,
    //   where envelopes are updated from within the multi-block, after the execution of internal processors.
    multiProcessing.updateEnvelopesWhereRequired();
}

template<typename T>
void addInternalProcessor( DataProcessorGenerator2D<T> const& generator,
                           std::vector<MultiBlock2D<T>*> multiBlocks, plint level )
{
    MultiProcessing2D<T, DataProcessorGenerator2D<T> const, DataProcessorGenerator2D<T> >
        multiProcessing(generator, multiBlocks);
    std::vector<DataProcessorGenerator2D<T>*> const& retainedGenerators = multiProcessing.getRetainedGenerators();
    std::vector<std::vector<plint> > const& atomicBlockNumbers = multiProcessing.getAtomicBlockNumbers();

    for (pluint iGenerator=0; iGenerator<retainedGenerators.size(); ++iGenerator) {
        std::vector<AtomicBlock2D<T>*> extractedAtomicBlocks(multiBlocks.size());
        for (pluint iBlock=0; iBlock<extractedAtomicBlocks.size(); ++iBlock) {
            extractedAtomicBlocks[iBlock] = &multiBlocks[iBlock]->getComponent(atomicBlockNumbers[iGenerator][iBlock]);
        }
        // Delegate to the "AtomicBlock version" of addInternal.
        plb::addInternalProcessor(*retainedGenerators[iGenerator], extractedAtomicBlocks, level);
    }
    // Subscribe the processor in the multi-block. This guarantees that the multi-block is aware
    //   of the maximal current processor level, and it instantiates the communication pattern
    //   for an update of envelopes after processor execution.
    multiBlocks[0]->subscribeProcessor (
            level,
            multiProcessing.multiBlocksWhichRequireUpdate(),
            BlockDomain::usesEnvelope(generator.appliesTo()) );
}

}  // namespace plb

#endif  // MULTI_BLOCK_OPERATIONS_2D_HH
