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
 * Global Definitions -- header file.
 */
#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

#ifdef PLB_MPI_PARALLEL
#include "mpi.h"
#endif

#include <string>

namespace plb {

/// Integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using flplint instead of plint ensures
 *  64-bit compatibility of the code.
 */
typedef ptrdiff_t plint;

/// Unsigned integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using fluplint instead of unsigned plint 
 *  ensures 64-bit compatibility of the code.
 **/
typedef size_t pluint;

/// Ordering of indices when a BlockXD is converted into a serial data stream.
/** Signification of constants:
 *    - forward:  Right-most index (y in 2D and z in 3D) is contiguous in memory.
 *                For non-allocated parts of the Block, the value 0 is produced for
 *                output, and values are ignored during input.
 *    - backward: Left-most index (x) is contiguous in memory.
 *                For non-allocated parts of the Block, the value 0 is produced for
 *                output, and values are ignored during input.
 *    - memorySaving: Ordering is forward (this respects the natural ordering in
 *                    Palabos). Non-allocated parts of the Block are neither written
 *                    or read: memory savings in the program are reflected by memory
 *                    savings on the disk.
 **/
namespace IndexOrdering {
    enum OrderingT {forward, backward, memorySaving};
}

/// Sub-domain of an atomic-block, on which for example a data processor is executed.
/** Signification of constants:
 *      - bulk: Refers to bulk-nodes, without envelope.
 *      - envelope: Refers to nodes on the envelope, without the bulk.
 *      - bulkAndEnvelope: Refers to the full domain.
 **/
namespace BlockDomain {

    enum DomainT {bulk, envelope, bulkAndEnvelope};
    inline bool usesEnvelope(DomainT domain) {
        return domain==envelope || domain==bulkAndEnvelope;
    }

}

namespace global {

class IOpolicyClass {
public:
    void setIndexOrderingForStreams(IndexOrdering::OrderingT streamOrdering_);
    IndexOrdering::OrderingT getIndexOrderingForStreams() const;
    void setEndianSwitchOnBase64out(bool doSwitch);
    bool getEndianSwitchOnBase64out();
    void setEndianSwitchOnBase64in(bool doSwitch);
    bool getEndianSwitchOnBase64in();
private:
    IOpolicyClass();
private:
    IndexOrdering::OrderingT streamOrdering;
    bool endianSwitchOnBase64out;
    bool endianSwitchOnBase64in;
    friend IOpolicyClass& IOpolicy();
};
    
class Directories {
public:
    /// Collectively set all output directories to the value specified by outputDir.
    void setOutputDir(std::string outputDir);
    /// Specify output directory for logfile.
    void setLogOutDir(std::string logOutDir_);
    /// Specify output directory for the ImageWriter.
    void setImageOutDir(std::string imageOutDir_);
    /// Specify output directory for the VtkDataWriter.
    void setVtkOutDir(std::string inputDir_);
    /// Specify directory for input files.
    void setInputDir(std::string inputDir);

    /// Get output directory for logfile.
    std::string getLogOutDir() const;
    /// Get output directory for the ImageWriter.
    std::string getImageOutDir() const;
    /// Get output directory for the VtkDataWriter.
    std::string getVtkOutDir() const;
    /// Get directory for input files.
    std::string getInputDir() const;
private:
    Directories();
private:
    std::string logOutDir;
    std::string imageOutDir;
    std::string vtkOutDir;
    std::string inputDir;
friend Directories& directories();
};

inline Directories& directories() {
    static Directories singleton;
    return singleton;
}

inline IOpolicyClass& IOpolicy() {
    static IOpolicyClass singleton;
    return singleton;
}

}  // namespace global

}  // namespace plb

#endif
