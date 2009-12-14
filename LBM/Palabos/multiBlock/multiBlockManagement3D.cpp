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
 * Geometry specifications for 3D multiblocks -- implementation.
 */

#include "multiBlock/multiBlockManagement3D.h"
#include "core/plbDebug.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

////////////////////// Class BlockParameters3D /////////////////////////////

BlockParameters3D::BlockParameters3D(Box3D bulk_, plint envelopeWidth_, plint procId_,
                                     bool leftX, bool rightX, bool leftY, bool rightY, bool leftZ, bool rightZ)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(bulk_),
      envelope(bulk.x0-envelopeWidth, bulk.x1+envelopeWidth, bulk.y0-envelopeWidth,
               bulk.y1+envelopeWidth, bulk.z0-envelopeWidth, bulk.z1+envelopeWidth),
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
    if (leftZ) {
        nonPeriodicEnvelope.z0 += envelopeWidth;
    }
    if (rightZ) {
        nonPeriodicEnvelope.z1 -= envelopeWidth;
    }
}

BlockParameters3D::BlockParameters3D (
        Box3D bulk_, Box3D envelope_, Box3D nonPeriodicEnvelope_,
        plint envelopeWidth_, plint procId_)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(bulk_),
      envelope(envelope_),
      nonPeriodicEnvelope(nonPeriodicEnvelope_)
{ }


////////////////////// Struct PeriodicOverlap3D /////////////////////

PeriodicOverlap3D::PeriodicOverlap3D(Overlap3D const& overlap_, plint normalX_, plint normalY_, plint normalZ_)
    : overlap(overlap_),
      normalX(normalX_),
      normalY(normalY_),
      normalZ(normalZ_)
{ }



////////////////////// Class MultiBlockDistribution3D /////////////////////

MultiBlockDistribution3D::MultiBlockDistribution3D(plint nx_, plint ny_, plint nz_)
    : boundingBox(0, nx_-1, 0, ny_-1, 0, nz_-1),
      periodicEnvelopeSize(0)
{ }

MultiBlockDistribution3D::MultiBlockDistribution3D(Box3D boundingBox_)
    : boundingBox(boundingBox_),
      periodicEnvelopeSize(0)
{ }

MultiBlockDistribution3D& MultiBlockDistribution3D::operator=(MultiBlockDistribution3D const& rhs) {
    boundingBox = rhs.boundingBox;
    periodicEnvelopeSize = rhs.periodicEnvelopeSize;
    blocks = rhs.blocks;
    normalOverlaps = rhs.normalOverlaps;
    periodicOverlaps = rhs.periodicOverlaps;
    return (*this);
}


BlockParameters3D const& MultiBlockDistribution3D::getBlockParameters(plint whichBlock) const {
    PLB_PRECONDITION( whichBlock < getNumBlocks() );
    return blocks[whichBlock];
}
Overlap3D const& MultiBlockDistribution3D::getNormalOverlap(plint whichOverlap) const {
    PLB_PRECONDITION( whichOverlap < getNumNormalOverlaps() );
    return normalOverlaps[whichOverlap];
}
PeriodicOverlap3D const& MultiBlockDistribution3D::getPeriodicOverlap(plint whichOverlap) const {
    PLB_PRECONDITION( whichOverlap < getNumPeriodicOverlaps() );
    return periodicOverlaps[whichOverlap];
}

void MultiBlockDistribution3D::addBlock(Box3D bulk, plint envelopeWidth, plint procId)
{
    PLB_PRECONDITION( contained(bulk, getBoundingBox()) );

    BlockParameters3D newBlock (
            bulk, envelopeWidth, procId,
            bulk.x0==boundingBox.x0, bulk.x1==boundingBox.x1,
            bulk.y0==boundingBox.y0, bulk.y1==boundingBox.y1,
            bulk.z0==boundingBox.z0, bulk.z1==boundingBox.z1 );

    computeNormalOverlaps(newBlock);
    blocks.push_back(newBlock);
    computePeriodicOverlaps();
}

plint MultiBlockDistribution3D::locate(plint x, plint y, plint z, plint guess) const {
    PLB_PRECONDITION( contained(x,y,z, getBoundingBox()) );
    PLB_PRECONDITION( guess < getNumBlocks() );

    for (plint iBlock=0; iBlock<(plint)getNumBlocks(); ++iBlock, guess = (guess+1)%blocks.size()) {
        Box3D const& coord = blocks[guess].getBulk();
        if (contained(x, y, z, coord)) {
            return guess;
        }
    }
    return -1;
}

void MultiBlockDistribution3D::computeNormalOverlaps(BlockParameters3D const& newBlock) {
    neighbors.resize(getNumBlocks()+1);
    Box3D intersection;
    plint iNew = getNumBlocks();
    for (plint iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
        if (intersect(blocks[iBlock].getBulk(), newBlock.getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(iBlock, iNew, intersection));
            neighbors[iBlock].push_back(iNew);
        }
        if (intersect(newBlock.getBulk(), blocks[iBlock].getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(iNew, iBlock, intersection));
            neighbors[iNew].push_back(iBlock);
        }
    }
}

void MultiBlockDistribution3D::computePeriodicOverlaps() {
    // It is assumed that a new block has been added; overlaps between
    //   this block and all existing ones will be identified
    plint iNew = getNumBlocks()-1;
    BlockParameters3D const& newBlock = blocks[iNew];
    Box3D intersection;
    for (plint dx=-1; dx<=+1; dx+=1) {
        for (plint dy=-1; dy<=+1; dy+=1) {
            for (plint dz=-1; dz<=+1; dz+=1) {
                if (dx!=0 || dy!=0 || dz!=0) {
                    // The new block is shifted by the length of the full multi block in each space
                    //   direction. Consequently, overlaps between the original multi block and the
                    //   shifted new block are identified as periodic overlaps.
                    plint shiftX = dx*getBoundingBox().getNx();
                    plint shiftY = dy*getBoundingBox().getNy();
                    plint shiftZ = dz*getBoundingBox().getNz();
                    Box3D newBulk(newBlock.getBulk().shift(shiftX,shiftY,shiftZ));
                    Box3D newEnvelope(newBlock.getEnvelope().shift(shiftX,shiftY,shiftZ));
                    // Check overlap which each existing block, including with the newly added one.
                    for (plint iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
                        // Does the envelope of the shifted new block overlap with the bulk of a previous
                        //   block? If yes, add an overlap, in which the previous block has the "original
                        //   position", and the new block has the "overlap position".
                        if (intersect(blocks[iBlock].getBulk(), newEnvelope, intersection)) {
                            periodicOverlaps.push_back (
                                    PeriodicOverlap3D (
                                        Overlap3D(iBlock, iNew, intersection, shiftX, shiftY, shiftZ),
                                        dx, dy, dz ) );
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
                            intersection = intersection.shift(-shiftX,-shiftY, -shiftZ);
                            periodicOverlaps.push_back (
                                    PeriodicOverlap3D (
                                        Overlap3D(iNew, iBlock, intersection, -shiftX, -shiftY, -shiftZ),
                                        -dx, -dy, -dz ) );
                            neighbors[iNew].push_back(iBlock);
                            updatePeriodicEnvelopeSize(intersection.getMaxWidth());
                        }
                    }
                }
            }
        }
    }
}

void MultiBlockDistribution3D::updatePeriodicEnvelopeSize(plint newSize) {
    if (newSize > periodicEnvelopeSize) {
        periodicEnvelopeSize = newSize;
    }
}

pluint MultiBlockDistribution3D::getNumAllocatedBulkCells() const {
    plint numCells = 0;
    for (pluint iBlock=0; iBlock<blocks.size(); ++iBlock) {
        numCells += (pluint)blocks[iBlock].getBulkLx() *
                    (pluint)blocks[iBlock].getBulkLy() *
                    (pluint)blocks[iBlock].getBulkLz();
    }
    return numCells;
}

bool MultiBlockDistribution3D::getNextChunkX(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        plint exploreX = iX+1;
        while(exploreX<getBoundingBox().getNx() && locate(exploreX,iY,iZ)==-1) {
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

bool MultiBlockDistribution3D::getNextChunkY(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        plint exploreY = iY+1;
        while(exploreY<getBoundingBox().getNy() && locate(iX,exploreY,iZ)==-1) {
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

bool MultiBlockDistribution3D::getNextChunkZ(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        plint exploreZ = iZ+1;
        while(exploreZ<getBoundingBox().getNz() && locate(iX,iY,exploreZ)==-1) {
            ++exploreZ;
        }
        nextChunkSize = exploreZ-iZ;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().z1-iZ+1;
        return true;
    }
}

RelevantIndexes3D::RelevantIndexes3D( MultiBlockDistribution3D const& multiBlockDistribution,
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

void RelevantIndexes3D::listAllIndexes (
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

void RelevantIndexes3D::computeRelevantIndexesInParallel(MultiBlockDistribution3D const& multiBlockDistribution,
                                                         ThreadAttribution const& attribution)
{
    for (plint iBlock=0; iBlock<multiBlockDistribution.getNumBlocks(); ++iBlock) {
        if (attribution.isLocal( multiBlockDistribution.getBlockParameters(iBlock).getProcId() ) ) {
            Box3D const& newBlock = multiBlockDistribution.getBlockParameters(iBlock).getEnvelope();
            if (myBlocks.empty()) {
                boundingBox = newBlock;
            }
            else {
                if (newBlock.x0 < boundingBox.x0) boundingBox.x0 = newBlock.x0;
                if (newBlock.x1 > boundingBox.x1) boundingBox.x1 = newBlock.x1;
                if (newBlock.y0 < boundingBox.y0) boundingBox.y0 = newBlock.y0;
                if (newBlock.y1 > boundingBox.y1) boundingBox.y1 = newBlock.y1;
                if (newBlock.z0 < boundingBox.z0) boundingBox.z0 = newBlock.z0;
                if (newBlock.z1 > boundingBox.z1) boundingBox.z1 = newBlock.z1;
            }
            myBlocks.push_back(iBlock);
            nearbyBlocks.push_back(iBlock);
        }
    }
    for (plint iOverlap=0; iOverlap<multiBlockDistribution.getNumNormalOverlaps(); ++iOverlap) {
        Overlap3D const& overlap = multiBlockDistribution.getNormalOverlap(iOverlap);
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
        Overlap3D const& overlap = multiBlockDistribution.getPeriodicOverlap(iOverlap).overlap;
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


MultiBlockManagement3D::MultiBlockManagement3D( MultiBlockDistribution3D const& multiBlockDistribution_,
                                                ThreadAttribution* threadAttribution_, plint refinementLevel_ )
    : multiBlockDistribution(multiBlockDistribution_),
      threadAttribution(threadAttribution_),
      relevantIndexes(multiBlockDistribution, getThreadAttribution()),
      refinementLevel(refinementLevel_)
{ }

MultiBlockManagement3D::MultiBlockManagement3D(MultiBlockManagement3D const& rhs)
    : multiBlockDistribution(rhs.multiBlockDistribution),
      threadAttribution(rhs.threadAttribution->clone()),
      relevantIndexes(rhs.relevantIndexes),
      refinementLevel(rhs.refinementLevel)
{ }

MultiBlockManagement3D& MultiBlockManagement3D::operator=(MultiBlockManagement3D const& rhs) {
    MultiBlockManagement3D newBlockManagement(rhs);
    newBlockManagement.swap(*this);
    return *this;
}

void MultiBlockManagement3D::swap(MultiBlockManagement3D& rhs) {
    std::swap(threadAttribution, rhs.threadAttribution);
    std::swap(multiBlockDistribution, rhs.multiBlockDistribution);
    std::swap(relevantIndexes, rhs.relevantIndexes);
}

MultiBlockManagement3D::~MultiBlockManagement3D() {
    delete threadAttribution;
}

Box3D const& MultiBlockManagement3D::getEnvelope(plint iBlock) const {
    return multiBlockDistribution.getBlockParameters(iBlock).getEnvelope();
}

ThreadAttribution const& MultiBlockManagement3D::getThreadAttribution() const {
    return *threadAttribution;
}


bool MultiBlockManagement3D::findInLocalBulk (
            plint iX, plint iY, plint iZ, plint& foundId, plint& localX, plint& localY, plint& localZ, plint guess ) const
{
    foundId = multiBlockDistribution.locate(iX,iY,iZ, guess);
    BlockParameters3D const& parameters = multiBlockDistribution.getBlockParameters(foundId);
    localX = parameters.toLocalX(iX);
    localY = parameters.toLocalY(iY);
    localZ = parameters.toLocalZ(iZ);
    return foundId >= 0;
}

bool MultiBlockManagement3D::findAllLocalRepresentations (
            plint iX, plint iY, plint iZ, std::vector<plint>& foundId,
            std::vector<plint>& foundX, std::vector<plint>& foundY, std::vector<plint>& foundZ ) const
{
    bool hasBulkCell = false;
    // First, search in all blocks which are local to the current processor, including in the envelopes.
    //   These blocks are confined within the boundingBox, so checking for inclusion in the boundingBox
    //   eliminates most queries and thus enhances efficiency.
    if (contained(iX,iY,iZ, relevantIndexes.getBoundingBox())) {
        for (pluint iBlock=0; iBlock < relevantIndexes.getBlocks().size(); ++iBlock) {
            BlockParameters3D const& parameters
                = multiBlockDistribution.getBlockParameters(relevantIndexes.getBlocks()[iBlock]);
            Box3D const& envelope = parameters.getEnvelope();
            if (contained(iX, iY, iZ, envelope)) {
                Box3D const& bulk = parameters.getBulk();
                if (contained(iX, iY, iZ, bulk)) {
                    hasBulkCell = true;
                    foundId.insert(foundId.begin(), relevantIndexes.getBlocks()[iBlock]);
                    foundX.insert(foundX.begin(), iX-envelope.x0);
                    foundY.insert(foundY.begin(), iY-envelope.y0);
                    foundZ.insert(foundZ.begin(), iZ-envelope.z0);
                }
                else {
                    foundId.push_back(relevantIndexes.getBlocks()[iBlock]);
                    foundX.push_back(iX-envelope.x0);
                    foundY.push_back(iY-envelope.y0);
                    foundZ.push_back(iZ-envelope.z0);
                }
            }
        }
    }
    // Here's a subtlety: with periodic boundary conditions, one may need to take into account
    //   a cell which is not inside the boundingBox, because it's at the opposite boundary.
    //   Therefore, this loop checks all blocks which overlap with the current one by periodicity.
    for (pluint iRelevant=0; iRelevant<relevantIndexes.getPeriodicOverlapWithRemoteData().size(); ++iRelevant) {
        plint iOverlap = relevantIndexes.getPeriodicOverlapWithRemoteData()[iRelevant];
        Overlap3D const& overlap = multiBlockDistribution.getPeriodicOverlap(iOverlap).overlap;
        if (contained(iX,iY,iZ, overlap.getOriginalCoordinates())) {
            plint overlapId = overlap.getOverlapId();
            foundId.push_back(overlapId);
            BlockParameters3D const& parameters = multiBlockDistribution.getBlockParameters(overlapId);
            foundX.push_back(parameters.toLocalX(iX-overlap.getShiftX()));
            foundY.push_back(parameters.toLocalY(iY-overlap.getShiftY()));
            foundZ.push_back(parameters.toLocalZ(iZ-overlap.getShiftZ()));
        }
    }
    return hasBulkCell;
}

plint MultiBlockManagement3D::getRefinementLevel() const {
    return refinementLevel;
}


MultiBlockDistribution3D extractMultiBlockDistribution (
        MultiBlockDistribution3D const& originalDistribution,
        Box3D subDomain, bool crop )
{
    MultiBlockDistribution3D extractedDistribution(crop ? subDomain : originalDistribution.getBoundingBox());
    for (plint iBlock=0; iBlock<originalDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& originalParameters = originalDistribution.getBlockParameters(iBlock);
        Box3D intersection;
        if (intersect(subDomain, originalParameters.getBulk(), intersection) ) {
            extractedDistribution.addBlock(intersection,
                                           originalParameters.getEnvelopeWidth(),
                                           originalParameters.getProcId() );
        }
    }
    return extractedDistribution;
}

MultiBlockManagement3D extractMultiBlockManagement (
        MultiBlockManagement3D const& originalManagement,
        Box3D subDomain, bool crop )
{
    return MultiBlockManagement3D (
            extractMultiBlockDistribution (
                originalManagement.getMultiBlockDistribution(), subDomain, crop ),
            originalManagement.getThreadAttribution().clone(),
            originalManagement.getRefinementLevel() );
}

}  // namespace plb
