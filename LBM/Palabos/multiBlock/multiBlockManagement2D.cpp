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
 * Geometry specifications for 2D multiblocks -- implementation.
 */

#include "multiBlock/multiBlockManagement2D.h"
#include "core/plbDebug.h"
#include <algorithm>


namespace plb {

////////////////////// Class BlockParameters2D /////////////////////////////

BlockParameters2D::BlockParameters2D(Box2D bulk_, plint envelopeWidth_, plint procId_,
                                     bool leftX, bool rightX, bool leftY, bool rightY)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(bulk_),
      envelope(bulk.x0-envelopeWidth, bulk.x1+envelopeWidth,
               bulk.y0-envelopeWidth, bulk.y1+envelopeWidth),
      nonPeriodicEnvelope(envelope)
{
    if (leftX) {
        nonPeriodicEnvelope.x0 += envelopeWidth;
    }
    if (rightX) {
        nonPeriodicEnvelope.x1 -= envelopeWidth;
    }

    if (leftY) {
        nonPeriodicEnvelope.y0 += envelopeWidth;
    }
    if (rightY) {
        nonPeriodicEnvelope.y1 -= envelopeWidth;
    }
}

BlockParameters2D::BlockParameters2D (
        Box2D bulk_, Box2D envelope_, Box2D nonPeriodicEnvelope_,
        plint envelopeWidth_, plint procId_)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(bulk_),
      envelope(envelope_),
      nonPeriodicEnvelope(nonPeriodicEnvelope_)
{ }


////////////////////// Struct PeriodicOverlap2D /////////////////////

PeriodicOverlap2D::PeriodicOverlap2D(Overlap2D const& overlap_, plint normalX_, plint normalY_)
    : overlap(overlap_),
      normalX(normalX_),
      normalY(normalY_)
{ }

////////////////////// Class MultiBlockDistribution2D /////////////////////

MultiBlockDistribution2D::MultiBlockDistribution2D(plint nx_, plint ny_)
    : boundingBox(0, nx_-1, 0, ny_-1),
      periodicEnvelopeSize(0)
{ }

MultiBlockDistribution2D::MultiBlockDistribution2D(Box2D boundingBox_)
    : boundingBox(boundingBox_),
      periodicEnvelopeSize(0)
{ }

MultiBlockDistribution2D& MultiBlockDistribution2D::operator=(MultiBlockDistribution2D const& rhs) {
    boundingBox = rhs.boundingBox;
    periodicEnvelopeSize = rhs.periodicEnvelopeSize;
    blocks = rhs.blocks;
    normalOverlaps = rhs.normalOverlaps;
    periodicOverlaps = rhs.periodicOverlaps;
    return *this;
}


BlockParameters2D const& MultiBlockDistribution2D::getBlockParameters(plint whichBlock) const {
    PLB_PRECONDITION( whichBlock < getNumBlocks() );
    return blocks[whichBlock];
}

Overlap2D const& MultiBlockDistribution2D::getNormalOverlap(plint whichOverlap) const {
    PLB_PRECONDITION( whichOverlap < getNumNormalOverlaps() );
    return normalOverlaps[whichOverlap];
}

PeriodicOverlap2D const& MultiBlockDistribution2D::getPeriodicOverlap(plint whichOverlap) const {
    PLB_PRECONDITION( whichOverlap < getNumPeriodicOverlaps() );
    return periodicOverlaps[whichOverlap];
}

void MultiBlockDistribution2D::addBlock(Box2D bulk, plint envelopeWidth, plint procId) {
    PLB_PRECONDITION( contained(bulk, getBoundingBox()) );

    BlockParameters2D newBlock (
            bulk, envelopeWidth, procId,
            bulk.x0==boundingBox.x0, bulk.x1==boundingBox.x1,
            bulk.y0==boundingBox.y0, bulk.y1==boundingBox.y1 );
    computeNormalOverlaps(newBlock);
    blocks.push_back(newBlock);
    computePeriodicOverlaps();
}

plint MultiBlockDistribution2D::locate(plint x, plint y, plint guess) const {
    PLB_PRECONDITION( contained(x,y, getBoundingBox()) );
    PLB_PRECONDITION( guess < getNumBlocks() );

    for (plint iBlock=0; iBlock<(plint)blocks.size(); ++iBlock, guess = (guess+1)%blocks.size()) {
        Box2D const& coord = blocks[guess].getBulk();
        if (contained(x, y, coord)) {
            return guess;
        }
    }
    return -1;
}

void MultiBlockDistribution2D::computeNormalOverlaps(BlockParameters2D const& newBlock) {
    neighbors.resize(getNumBlocks()+1);
    Box2D intersection;
    plint iNew = getNumBlocks();
    for (plint iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
        if (intersect(blocks[iBlock].getBulk(), newBlock.getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap2D(iBlock, iNew, intersection));
            neighbors[iBlock].push_back(iNew);
        }
        if (intersect(newBlock.getBulk(), blocks[iBlock].getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap2D(iNew, iBlock, intersection));
            neighbors[iNew].push_back(iBlock);
        }
    }
}

void MultiBlockDistribution2D::computePeriodicOverlaps() {
    // It is assumed that a new block has been added; overlaps between
    //   this block and all existing ones will be identified
    plint iNew = getNumBlocks()-1;
    BlockParameters2D const& newBlock = blocks[iNew];
    Box2D intersection;
    for (plint dx=-1; dx<=+1; dx+=1) {
        for (plint dy=-1; dy<=+1; dy+=1) {
            if (dx!=0 || dy!=0) {
                // The new block is shifted by the length of the full multi block in each space
                //   direction. Consequently, overlaps between the original multi block and the
                //   shifted new block are identified as periodic overlaps.
                plint shiftX = dx*getBoundingBox().getNx();
                plint shiftY = dy*getBoundingBox().getNy();
                Box2D newBulk(newBlock.getBulk().shift(shiftX,shiftY));
                Box2D newEnvelope(newBlock.getEnvelope().shift(shiftX,shiftY));
                // Check overlap which each existing block, including with the newly added one.
                for (plint iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
                    // Does the envelope of the shifted new block overlap with the bulk of a previous
                    //   block? If yes, add an overlap, in which the previous block has the "original
                    //   position", and the new block has the "overlap position".
                    if (intersect(blocks[iBlock].getBulk(), newEnvelope, intersection)) {
                        periodicOverlaps.push_back (
                                PeriodicOverlap2D (
                                    Overlap2D(iBlock, iNew, intersection, shiftX, shiftY),
                                    dx, dy ) );
                        neighbors[iBlock].push_back(iNew);
                        updatePeriodicEnvelopeSize(intersection.getMaxWidth());
                    }
                    // Does the bulk of the shifted new block overlap with the envelope of a previous
                    //   block? If yes, add an overlap, in which the new block has the "original position",
                    //   and the previous block has the "overlap position".
                    //   If we are in the situation in which the newly added block is periodic with itself,
                    //   this step must be skipped, because otherwise the overlap is counted twice.
                    if (!(iBlock==iNew) &&
                        intersect(newBulk, blocks[iBlock].getEnvelope(), intersection))
                    {
                        intersection = intersection.shift(-shiftX,-shiftY);
                        periodicOverlaps.push_back (
                                PeriodicOverlap2D (
                                    Overlap2D(iNew, iBlock, intersection, -shiftX, -shiftY),
                                    -dx, -dy ) );
                        neighbors[iNew].push_back(iBlock);
                        updatePeriodicEnvelopeSize(intersection.getMaxWidth());
                    }
                }
            }
        }
    }
}

void MultiBlockDistribution2D::updatePeriodicEnvelopeSize(plint newSize) {
    if (newSize > periodicEnvelopeSize) {
        periodicEnvelopeSize = newSize;
    }
}

pluint MultiBlockDistribution2D::getNumAllocatedBulkCells() const {
    pluint numCells = 0;
    for (pluint iBlock=0; iBlock<blocks.size(); ++iBlock) {
        numCells += (pluint)blocks[iBlock].getBulkLx() * (pluint)blocks[iBlock].getBulkLy();
    }
    return numCells;
}

bool MultiBlockDistribution2D::getNextChunkX(plint iX, plint iY, plint& nextLattice, plint& nextChunkSize) const {
    nextLattice = locate(iX,iY);
    if (nextLattice == -1) {
        plint exploreX = iX+1;
        while(exploreX<getBoundingBox().getNx() && locate(exploreX,iY)==-1) {
            ++exploreX;
        }
        nextChunkSize = exploreX-iX;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().x1-iX+1;
        return true;
    }
}

bool MultiBlockDistribution2D::getNextChunkY(plint iX, plint iY, plint& nextLattice, plint& nextChunkSize) const {
    nextLattice = locate(iX,iY);
    if (nextLattice == -1) {
        plint exploreY = iY+1;
        while(exploreY<getBoundingBox().getNy() && locate(iX,exploreY)==-1) {
            ++exploreY;
        }
        nextChunkSize = exploreY-iY;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().y1-iY+1;
        return true;
    }
}

RelevantIndexes2D::RelevantIndexes2D( MultiBlockDistribution2D const& multiBlockDistribution,
                                      ThreadAttribution const& attribution )
{
    if (attribution.allBlocksAreLocal()) {
        listAllIndexes (
                multiBlockDistribution.getNumBlocks(),
                multiBlockDistribution.getNumNormalOverlaps(),
                multiBlockDistribution.getNumPeriodicOverlaps() );
        boundingBox = multiBlockDistribution.getBoundingBox();
    }
    else {
        computeRelevantIndexesInParallel(multiBlockDistribution, attribution);
    }
}

void RelevantIndexes2D::listAllIndexes (
        plint numBlocks, plint numNormalOverlaps, plint numPeriodicOverlaps )
{
    myBlocks.resize(numBlocks);
    nearbyBlocks.resize(numBlocks);
    for (plint iBlock=0; iBlock<numBlocks; ++iBlock) {
        myBlocks[iBlock]     = iBlock;
        nearbyBlocks[iBlock] = iBlock;
    }
    normalOverlaps.resize(numNormalOverlaps);
    for (plint iOverlap=0; iOverlap<numNormalOverlaps; ++iOverlap) {
        normalOverlaps[iOverlap] = iOverlap;
    }
    periodicOverlaps.resize(numPeriodicOverlaps);
    periodicOverlapWithRemoteData.resize(numPeriodicOverlaps);
    for (plint iOverlap=0; iOverlap<numPeriodicOverlaps; ++iOverlap) {
        periodicOverlaps[iOverlap] = iOverlap;
        periodicOverlapWithRemoteData[iOverlap] = iOverlap;
    }
}

void RelevantIndexes2D::computeRelevantIndexesInParallel(MultiBlockDistribution2D const& multiBlockDistribution,
                                                         ThreadAttribution const& attribution)
{
    for (plint iBlock=0; iBlock<multiBlockDistribution.getNumBlocks(); ++iBlock) {
        if (attribution.isLocal( multiBlockDistribution.getBlockParameters(iBlock).getProcId() ) ) {
            Box2D const& newBlock = multiBlockDistribution.getBlockParameters(iBlock).getEnvelope();
            if (myBlocks.empty()) {
                boundingBox = newBlock;
            }
            else {
                if (newBlock.x0 < boundingBox.x0) boundingBox.x0 = newBlock.x0;
                if (newBlock.x1 > boundingBox.x1) boundingBox.x1 = newBlock.x1;
                if (newBlock.y0 < boundingBox.y0) boundingBox.y0 = newBlock.y0;
                if (newBlock.y1 > boundingBox.y1) boundingBox.y1 = newBlock.y1;
            }
            myBlocks.push_back(iBlock);
            nearbyBlocks.push_back(iBlock);
        }
    }
    for (plint iOverlap=0; iOverlap<multiBlockDistribution.getNumNormalOverlaps(); ++iOverlap) {
        Overlap2D const& overlap = multiBlockDistribution.getNormalOverlap(iOverlap);
        plint originalProc = multiBlockDistribution.getBlockParameters(overlap.getOriginalId()).getProcId();
        plint overlapProc = multiBlockDistribution.getBlockParameters(overlap.getOverlapId()).getProcId();
        if (attribution.isLocal(originalProc)) {
            nearbyBlocks.push_back( overlap.getOverlapId() );
        }
        if (attribution.isLocal(originalProc) || attribution.isLocal(overlapProc)) {
            normalOverlaps.push_back(iOverlap);
        }
    }
    for (plint iOverlap=0; iOverlap<multiBlockDistribution.getNumPeriodicOverlaps(); ++iOverlap) {
        Overlap2D const& overlap = multiBlockDistribution.getPeriodicOverlap(iOverlap).overlap;
        plint originalProc = multiBlockDistribution.getBlockParameters(overlap.getOriginalId()).getProcId();
        plint overlapProc = multiBlockDistribution.getBlockParameters(overlap.getOverlapId()).getProcId();
        if (attribution.isLocal(originalProc)) {
            nearbyBlocks.push_back( overlap.getOverlapId() );
        }
        if (attribution.isLocal(overlapProc)) {
            periodicOverlapWithRemoteData.push_back(iOverlap);
        }
        if (attribution.isLocal(originalProc) || attribution.isLocal(overlapProc)) {
            periodicOverlaps.push_back(iOverlap);
        }
    }
    // Erase duplicates in nearbyBlocks
    std::sort(nearbyBlocks.begin(), nearbyBlocks.end());
    std::vector<plint>::iterator newEnd = unique(nearbyBlocks.begin(), nearbyBlocks.end());
    nearbyBlocks.erase(newEnd, nearbyBlocks.end());
}

////////////////////// Class MultiBlockManagement2D /////////////////////////////

MultiBlockManagement2D::MultiBlockManagement2D( MultiBlockDistribution2D const& multiBlockDistribution_,
                                                ThreadAttribution* threadAttribution_, plint refinementLevel_ )
    : multiBlockDistribution(multiBlockDistribution_),
      threadAttribution(threadAttribution_),
      relevantIndexes(multiBlockDistribution, getThreadAttribution()),
      refinementLevel(refinementLevel_)
{ }

MultiBlockManagement2D::MultiBlockManagement2D(MultiBlockManagement2D const& rhs)
    : multiBlockDistribution(rhs.multiBlockDistribution),
      threadAttribution(rhs.threadAttribution->clone()),
      relevantIndexes(rhs.relevantIndexes),
      refinementLevel(rhs.refinementLevel)
{ }

MultiBlockManagement2D& MultiBlockManagement2D::operator=(MultiBlockManagement2D const& rhs) {
    MultiBlockManagement2D newBlockManagement(rhs);
    newBlockManagement.swap(*this);
    return *this;
}

void MultiBlockManagement2D::swap(MultiBlockManagement2D& rhs) {
    std::swap(threadAttribution, rhs.threadAttribution);
    std::swap(multiBlockDistribution, rhs.multiBlockDistribution);
    std::swap(relevantIndexes, rhs.relevantIndexes);
}

MultiBlockManagement2D::~MultiBlockManagement2D() {
    delete threadAttribution;
}

Box2D const& MultiBlockManagement2D::getEnvelope(plint iBlock) const {
    return multiBlockDistribution.getBlockParameters(iBlock).getEnvelope();
}

ThreadAttribution const& MultiBlockManagement2D::getThreadAttribution() const {
    return *threadAttribution;
}


bool MultiBlockManagement2D::findInLocalBulk (
            plint iX, plint iY, plint& foundId, plint& localX, plint& localY, plint guess ) const
{
    foundId = multiBlockDistribution.locate(iX,iY, guess);
    BlockParameters2D const& parameters = multiBlockDistribution.getBlockParameters(foundId);
    localX = parameters.toLocalX(iX);
    localY = parameters.toLocalY(iY);
    return foundId >= 0;
}

bool MultiBlockManagement2D::findAllLocalRepresentations (
            plint iX, plint iY, std::vector<plint>& foundId,
            std::vector<plint>& foundX, std::vector<plint>& foundY ) const
{
    bool hasBulkCell = false;
    // First, search in all blocks which are local to the current processor, including in the envelopes.
    //   These blocks are confined within the boundingBox, so checking for inclusion in the boundingBox
    //   eliminates most queries and thus enhances efficiency.
    if (contained(iX,iY, relevantIndexes.getBoundingBox())) {
        for (pluint iBlock=0; iBlock < relevantIndexes.getBlocks().size(); ++iBlock) {
            BlockParameters2D const& parameters
                = multiBlockDistribution.getBlockParameters(relevantIndexes.getBlocks()[iBlock]);
            Box2D const& envelope = parameters.getEnvelope();
            if (contained(iX, iY, envelope)) {
                Box2D const& bulk = parameters.getBulk();
                if (contained(iX, iY, bulk)) {
                    hasBulkCell = true;
                    foundId.insert(foundId.begin(), relevantIndexes.getBlocks()[iBlock]);
                    foundX.insert(foundX.begin(), iX-envelope.x0);
                    foundY.insert(foundY.begin(), iY-envelope.y0);
                }
                else {
                    foundId.push_back(relevantIndexes.getBlocks()[iBlock]);
                    foundX.push_back(iX-envelope.x0);
                    foundY.push_back(iY-envelope.y0);
                }
            }
        }
    }
    // Here's a subtlety: with periodic boundary conditions, one may need to take into account
    //   a cell which is not inside the boundingBox, because it's at the opposite boundary.
    //   Therefore, this loop checks all blocks which overlap with the current one by periodicity.
    for (pluint iRelevant=0; iRelevant<relevantIndexes.getPeriodicOverlapWithRemoteData().size(); ++iRelevant) {
        plint iOverlap = relevantIndexes.getPeriodicOverlapWithRemoteData()[iRelevant];
        Overlap2D const& overlap = multiBlockDistribution.getPeriodicOverlap(iOverlap).overlap;
        if (contained(iX,iY, overlap.getOriginalCoordinates())) {
            plint overlapId = overlap.getOverlapId();
            foundId.push_back(overlapId);
            BlockParameters2D const& parameters = multiBlockDistribution.getBlockParameters(overlapId);
            foundX.push_back(parameters.toLocalX(iX-overlap.getShiftX()));
            foundY.push_back(parameters.toLocalY(iY-overlap.getShiftY()));
        }
    }
    return hasBulkCell;
}

plint MultiBlockManagement2D::getRefinementLevel() const {
    return refinementLevel;
}

MultiBlockDistribution2D extractMultiBlockDistribution (
        MultiBlockDistribution2D const& originalDistribution,
        Box2D subDomain, bool crop )
{
    MultiBlockDistribution2D extractedDistribution(crop ? subDomain : originalDistribution.getBoundingBox());
    for (plint iBlock=0; iBlock<originalDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters2D const& originalParameters = originalDistribution.getBlockParameters(iBlock);
        Box2D intersection;
        if (intersect(subDomain, originalParameters.getBulk(), intersection) ) {
            extractedDistribution.addBlock(intersection,
                                           originalParameters.getEnvelopeWidth(),
                                           originalParameters.getProcId() );
        }
    }
    return extractedDistribution;
}

MultiBlockManagement2D extractMultiBlockManagement (
        MultiBlockManagement2D const& originalManagement,
        Box2D subDomain, bool crop )
{
    return MultiBlockManagement2D (
            extractMultiBlockDistribution (
                originalManagement.getMultiBlockDistribution(), subDomain, crop ),
            originalManagement.getThreadAttribution().clone(),
            originalManagement.getRefinementLevel() );

}

}  // namespace plb
