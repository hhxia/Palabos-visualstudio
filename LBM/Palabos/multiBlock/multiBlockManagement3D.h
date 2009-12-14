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
 * Geometry specifications for 3D multiblocks -- header file.
 */

#ifndef MULTI_BLOCK_ARRANGEMENT_3D_H
#define MULTI_BLOCK_ARRANGEMENT_3D_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "multiBlock/threadAttribution.h"
#include <vector>

namespace plb {

class BlockParameters3D {
public:
    BlockParameters3D(Box3D bulk_, plint envelopeWidth_, plint procId_,
                      bool leftX, bool rightX, bool leftY,
                      bool rightY, bool leftZ, bool rightZ);
    BlockParameters3D(Box3D bulk_, Box3D envelope_, Box3D nonPeriodicEnvelope_,
                      plint envelopeWidth_, plint procId_);
    plint getEnvelopeWidth()     const { return envelopeWidth; }
    plint getProcId()            const { return procId; }
    Box3D const& getBulk()      const { return bulk; }
    Box3D const& getEnvelope( ) const { return envelope; }
    Box3D const& getNonPeriodicEnvelope() const { return nonPeriodicEnvelope; }
    plint getBulkLx()            const { return bulk.getNx(); }
    plint getBulkLy()            const { return bulk.getNy(); }
    plint getBulkLz()            const { return bulk.getNz(); }
    plint getEnvelopeLx()        const { return envelope.getNx(); }
    plint getEnvelopeLy()        const { return envelope.getNy(); }
    plint getEnvelopeLz()        const { return envelope.getNz(); }
    plint toLocalX(plint iX)      const { return iX-envelope.x0; }
    plint toLocalY(plint iY)      const { return iY-envelope.y0; }
    plint toLocalZ(plint iZ)      const { return iZ-envelope.z0; }
    Box3D toLocal(Box3D const& coord) const {
        return Box3D(coord.shift(-envelope.x0, -envelope.y0, -envelope.z0));
    }
private:
    plint envelopeWidth, procId;
    Box3D bulk, envelope, nonPeriodicEnvelope;
};

class Overlap3D {
public:
    Overlap3D(plint originalId_, plint overlapId_, Box3D const& intersection_)
        : originalId(originalId_), overlapId(overlapId_),
          originalRegion(intersection_),
          overlapRegion(intersection_)
    { }
    Overlap3D(plint originalId_, plint overlapId_,
              Box3D const& originalRegion_,
              plint shiftX, plint shiftY, plint shiftZ)
        : originalId(originalId_), overlapId(overlapId_),
          originalRegion(originalRegion_),
          overlapRegion(originalRegion.shift(-shiftX, -shiftY, -shiftZ))
    { }
    plint getOriginalId() const { return originalId; }
    plint getOverlapId()  const { return overlapId; }
    Box3D const& getOriginalCoordinates() const { return originalRegion; }
    Box3D const& getOverlapCoordinates() const  { return overlapRegion; }
    plint getShiftX() const { return originalRegion.x0 - overlapRegion.x0; }
    plint getShiftY() const { return originalRegion.y0 - overlapRegion.y0; }
    plint getShiftZ() const { return originalRegion.z0 - overlapRegion.z0; }
private: 
    plint originalId, overlapId;
    Box3D originalRegion, overlapRegion;
};

/// This structure holds both overlap information and orientation of the boundary.
/** In case of periodic overlaps, it is important to know the orientation of the
 *  boundary, additionally to the coordinates of the overlap region. This is
 *  required when the communication step within a multi block is executed. Given
 *  that the user can selectively swith on/off periodicity, the multi block
 *  must be able to decide which periodic overlaps to communicate and which not.
 */
struct PeriodicOverlap3D {
    PeriodicOverlap3D(Overlap3D const& overlap_, plint normalX_, plint normalY_, plint normalZ_);
    Overlap3D overlap;
    plint      normalX;
    plint      normalY;
    plint      normalZ;
};

class MultiBlockDistribution3D {
public:
    MultiBlockDistribution3D(plint nx_, plint ny_, plint nz_);
    MultiBlockDistribution3D(Box3D boundingBox_);
    MultiBlockDistribution3D& operator=(MultiBlockDistribution3D const& rhs);
    Box3D const& getBoundingBox() const { return boundingBox; }
    plint getNumBlocks()            const { return blocks.size(); }
    plint getNumNormalOverlaps()    const { return normalOverlaps.size(); }
    plint getNumPeriodicOverlaps()  const { return periodicOverlaps.size(); }
    plint getPeriodicEnvelopeSize() const { return periodicEnvelopeSize; }
    void addBlock(Box3D bulk, plint envelopeWidth, plint procId=0);
    BlockParameters3D const& getBlockParameters(plint whichBlock) const;
    Overlap3D   const& getNormalOverlap(plint whichOverlap) const;
    PeriodicOverlap3D const& getPeriodicOverlap(plint whichOverlap) const;
    plint locate(plint iX, plint iY, plint iZ, plint guess=0) const;
    pluint getNumAllocatedBulkCells() const;
    bool getNextChunkX(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const;
    bool getNextChunkY(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const;
    bool getNextChunkZ(plint iX, plint iY, plint iZ, plint& nextLattice, plint& nextChunkSize) const;
private:
    void computeNormalOverlaps(BlockParameters3D const& newBlock);
    void computePeriodicOverlaps();
    void updatePeriodicEnvelopeSize(plint newSize);
private:
    Box3D boundingBox;
    plint periodicEnvelopeSize;
    std::vector<BlockParameters3D> blocks;
    std::vector<Overlap3D> normalOverlaps;
    std::vector<PeriodicOverlap3D> periodicOverlaps;
    std::vector<std::vector<plint> > neighbors;
};

/// Indexes of Blocks and Overlaps which are relevant in the parallel case
class RelevantIndexes3D {
public:
    RelevantIndexes3D(MultiBlockDistribution3D const& multiBlockDistribution, ThreadAttribution const& attribution);
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
    Box3D const& getBoundingBox()                       const { return boundingBox; }
private:
    void listAllIndexes(plint numBlocks, plint numNormalOverlaps, plint numPeriodicOverlaps);
    void computeRelevantIndexesInParallel(MultiBlockDistribution3D const& multiBlockDistribution,
                                          ThreadAttribution const& attribution);
private:
    std::vector<plint> myBlocks;
    std::vector<plint> nearbyBlocks;
    std::vector<plint> normalOverlaps;
    std::vector<plint> periodicOverlaps;
    std::vector<plint> periodicOverlapWithRemoteData;
    Box3D boundingBox;
};

class MultiBlockManagement3D {
public:
    MultiBlockManagement3D( MultiBlockDistribution3D const& multiBlockDistribution_,
                            ThreadAttribution* threadAttribution_,
                            plint refinementLevel_ =0);
    MultiBlockManagement3D(MultiBlockManagement3D const& rhs);
    MultiBlockManagement3D& operator=(MultiBlockManagement3D const& rhs);
    void swap(MultiBlockManagement3D& rhs);
    ~MultiBlockManagement3D();
    Box3D const& getBoundingBox() const { return multiBlockDistribution.getBoundingBox(); }
    MultiBlockDistribution3D const& getMultiBlockDistribution() const { return multiBlockDistribution; }
    RelevantIndexes3D const& getRelevantIndexes() const { return relevantIndexes; }
    Box3D const& getEnvelope(plint iBlock) const;
    ThreadAttribution const& getThreadAttribution() const;
    bool findInLocalBulk (
            plint iX, plint iY, plint iZ, plint& foundId, plint& localX, plint& localY, plint& localZ, plint guess=0 ) const;
    bool findAllLocalRepresentations (
            plint iX, plint iY, plint iZ, std::vector<plint>& foundId,
            std::vector<plint>& foundX, std::vector<plint>& foundY, std::vector<plint>& foundZ ) const;
    plint getRefinementLevel() const;
private:
    MultiBlockDistribution3D multiBlockDistribution;
    ThreadAttribution*  threadAttribution;
    RelevantIndexes3D   relevantIndexes;
    plint                refinementLevel;
};

/// Create a new distribution, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new distribution is equal to
 *  the specified sub-domain. If crop is false, the bounding-box is the same as
 *  the bounding-box of the original distribution.
 */
MultiBlockDistribution3D extractMultiBlockDistribution (
        MultiBlockDistribution3D const& originalDistribution,
        Box3D subDomain, bool crop );

/// Create a new block-management, corresponding to a sub-domain of the old one.
/** If the parameter crop is true, the bounding-box of the new block-management is
 *  equal to the specified sub-domain. If crop is false, the bounding-box is the same
 *  as the bounding-box of the original block-management.
 */
MultiBlockManagement3D extractMultiBlockManagement (
        MultiBlockManagement3D const& originalManagement,
        Box3D subDomain, bool crop );

}  // namespace plb

#endif
