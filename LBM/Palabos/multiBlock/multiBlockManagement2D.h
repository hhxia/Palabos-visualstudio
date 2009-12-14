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
 * Geometry specifications for 2D multiblocks -- header file.
 */

#ifndef MULTI_BLOCK_ARRANGEMENT_2D_H
#define MULTI_BLOCK_ARRANGEMENT_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "multiBlock/threadAttribution.h"
#include <vector>

namespace plb {

class BlockParameters2D {
public:
    BlockParameters2D(Box2D bulk_, plint envelopeWidth_, plint procId_,
                      bool leftX, bool rightX, bool leftY, bool rightY);
    BlockParameters2D(Box2D bulk_, Box2D envelope_, Box2D nonPeriodicEnvelope_,
                      plint envelopeWidth_, plint procId_);
    plint getEnvelopeWidth()    const { return envelopeWidth; }
    plint getProcId()           const { return procId; }
    Box2D const& getBulk()     const { return bulk; }
    Box2D const& getEnvelope() const { return envelope; }
    Box2D const& getNonPeriodicEnvelope() const { return nonPeriodicEnvelope; }
    plint getBulkLx()           const { return bulk.getNx(); }
    plint getBulkLy()           const { return bulk.getNy(); }
    plint getEnvelopeLx()       const { return envelope.getNx(); }
    plint getEnvelopeLy()       const { return envelope.getNy(); }
    plint toLocalX(plint iX)     const { return iX-envelope.x0; }
    plint toLocalY(plint iY)     const { return iY-envelope.y0; }
    Box2D toLocal(Box2D const& coord) const {
        return Box2D(coord.shift(-envelope.x0, -envelope.y0));
    }
private:
    plint envelopeWidth, procId;
    Box2D bulk, envelope, nonPeriodicEnvelope;
};

class Overlap2D {
public:
    Overlap2D(plint originalId_, plint overlapId_, Box2D const& intersection_)
        : originalId(originalId_), overlapId(overlapId_),
          originalRegion(intersection_),
          overlapRegion(intersection_)
    { }
    Overlap2D(plint originalId_, plint overlapId_,
              Box2D const& originalRegion_,
              plint shiftX, plint shiftY)
        : originalId(originalId_), overlapId(overlapId_),
          originalRegion(originalRegion_),
          overlapRegion(originalRegion.shift(-shiftX, -shiftY))
    { }
    plint getOriginalId() const { return originalId; }
    plint getOverlapId()  const { return overlapId; }
    Box2D const& getOriginalCoordinates() const { return originalRegion; }
    Box2D const& getOverlapCoordinates() const  { return overlapRegion; }
    plint getShiftX() const { return originalRegion.x0 - overlapRegion.x0; }
    plint getShiftY() const { return originalRegion.y0 - overlapRegion.y0; }
private: 
    plint originalId, overlapId;
    Box2D originalRegion, overlapRegion;
};

/// This structure holds both overlap information and orientation of the boundary.
/** In case of periodic overlaps, it is important to know the orientation of the
 *  boundary, additionally to the coordinates of the overlap region. This is
 *  required when the communication step within a multi block is executed. Given
 *  that the user can selectively swith on/off periodicity, the multi block
 *  must be able to decide which periodic overlaps to communicate and which not.
 */
struct PeriodicOverlap2D {
    PeriodicOverlap2D(Overlap2D const& overlap_, plint normalX_, plint normalY_);
    Overlap2D overlap;
    plint      normalX;
    plint      normalY;
};

class MultiBlockDistribution2D {
public:
    MultiBlockDistribution2D(plint nx_, plint ny_);
    MultiBlockDistribution2D(Box2D boundingBox_);
    MultiBlockDistribution2D& operator=(MultiBlockDistribution2D const& rhs);
    Box2D const& getBoundingBox()  const { return boundingBox; }
    plint getNumBlocks()            const { return blocks.size(); }
    plint getNumNormalOverlaps()    const { return normalOverlaps.size(); }
    plint getNumPeriodicOverlaps()  const { return periodicOverlaps.size(); }
    plint getPeriodicEnvelopeSize() const { return periodicEnvelopeSize; }
    void addBlock(Box2D domain, plint envelopeWidth, plint procId=0);
    BlockParameters2D const& getBlockParameters(plint whichBlock) const;
    Overlap2D   const& getNormalOverlap(plint whichOverlap) const;
    PeriodicOverlap2D const& getPeriodicOverlap(plint whichOverlap) const;
    plint locate(plint iX, plint iY, plint guess=0) const;
    pluint getNumAllocatedBulkCells() const;
    bool getNextChunkX(plint iX, plint iY, plint& nextLattice, plint& nextChunkSize) const;
    bool getNextChunkY(plint iX, plint iY, plint& nextLattice, plint& nextChunkSize) const;
private:
    void computeNormalOverlaps(BlockParameters2D const& newBlock);
    void computePeriodicOverlaps();
    void updatePeriodicEnvelopeSize(plint newSize);
private:
    Box2D boundingBox;
    plint periodicEnvelopeSize;
    std::vector<BlockParameters2D> blocks;
    std::vector<Overlap2D> normalOverlaps;
    std::vector<PeriodicOverlap2D> periodicOverlaps;
    std::vector<std::vector<plint> > neighbors;
};

/// Indexes of Blocks and Overlaps which are relevant in the parallel case
class RelevantIndexes2D {
public:
    RelevantIndexes2D(MultiBlockDistribution2D const& multiBlockDistribution, ThreadAttribution const& attribution);
    /// Index of all blocks local to current processor
    std::vector<plint> const& getBlocks()                const { return myBlocks; }
    /// Index of all blocks with which current processor has communication
    std::vector<plint> const& getNearbyBlocks()          const { return nearbyBlocks; }
    /// Index of all overlaps for which original or overlap data are on current processor
    std::vector<plint> const& getNormalOverlaps()        const { return normalOverlaps; }
    /// Index of all periodic overlaps for which original or overlap data are on current processor
    std::vector<plint> const& getPeriodicOverlaps()      const { return periodicOverlaps; }
    /// Index of all periodic overlaps for which overlap data are on current processor
    std::vector<plint> const& getPeriodicOverlapWithRemoteData() const { return periodicOverlapWithRemoteData; }
    /// Bounding box for the envelope of all blocks which are on current processor
    Box2D const& getBoundingBox()                       const { return boundingBox; }
private:
    void listAllIndexes(plint numBlocks, plint numNormalOverlaps, plint numPeriodicOverlaps);
    void computeRelevantIndexesInParallel(MultiBlockDistribution2D const& multiBlockDistribution,
                                          ThreadAttribution const& attribution);
private:
    std::vector<plint> myBlocks;
    std::vector<plint> nearbyBlocks;
    std::vector<plint> normalOverlaps;
    std::vector<plint> periodicOverlaps;
    std::vector<plint> periodicOverlapWithRemoteData;
    Box2D boundingBox;
};

class MultiBlockManagement2D {
public:
    MultiBlockManagement2D( MultiBlockDistribution2D const& multiBlockDistribution_,
                            ThreadAttribution* threadAttribution_,
                            plint refinementLevel_ =0);
    MultiBlockManagement2D(MultiBlockManagement2D const& rhs);
    MultiBlockManagement2D& operator=(MultiBlockManagement2D const& rhs);
    void swap(MultiBlockManagement2D& rhs);
    ~MultiBlockManagement2D();
    Box2D const& getBoundingBox() const { return multiBlockDistribution.getBoundingBox(); }
    MultiBlockDistribution2D const& getMultiBlockDistribution() const { return multiBlockDistribution; }
    RelevantIndexes2D const& getRelevantIndexes() const { return relevantIndexes; }
    Box2D const& getEnvelope(plint iBlock) const;
    ThreadAttribution const& getThreadAttribution() const;
    bool findInLocalBulk (
            plint iX, plint iY, plint& foundId, plint& localX, plint& localY, plint guess=0 ) const;
    bool findAllLocalRepresentations (
            plint iX, plint iY, std::vector<plint>& foundId,
            std::vector<plint>& foundX, std::vector<plint>& foundY ) const;
    plint getRefinementLevel() const;
private:
    MultiBlockDistribution2D multiBlockDistribution;
    ThreadAttribution*       threadAttribution;
    RelevantIndexes2D        relevantIndexes;
    plint                     refinementLevel;
};

/// Create a new distribution, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new distribution is equal to
 *  the specified sub-domain. If crop is false, the bounding-box is the same as
 *  the bounding-box of the original distribution.
 */
MultiBlockDistribution2D extractMultiBlockDistribution (
        MultiBlockDistribution2D const& originalDistribution,
        Box2D subDomain, bool crop );

/// Create a new block-management, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new block-management is
 *  equal to the specified sub-domain. If crop is false, the bounding-box is the same
 *  as the bounding-box of the original block-management.
 */
MultiBlockManagement2D extractMultiBlockManagement (
        MultiBlockManagement2D const& originalBlockManagement,
        Box2D subDomain, bool crop );

}  // namespace plb

#endif
