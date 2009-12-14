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

#ifndef SERIALIZER_IO_3D_H
#define SERIALIZER_IO_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "multiBlock/multiBlock3D.h"
#include "io/serializerIO.h"
#include "io/parallelIO.h"

namespace plb {

/// Save the content of a Block3D into a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 */
template<typename T>
void saveBinaryBlock(Block3D<T> const& block, std::string fName, bool enforceUint=false);

/// Load the content of a Block3D from a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 */
template<typename T>
void loadBinaryBlock(Block3D<T>& block, std::string fName, bool enforceUint=false);

/// Save the raw content of a MultiBlock3D efficiently into a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is read, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This function is more efficient than the function saveBinaryBlock, but it is
 *  dependent on the structure of the multi-block. If data is written with saveRawMultiBlock()
 *  in a parallel program, the program which loads the data with loadRawMultiBlock() must
 *  have exactly the same number of cores.
 */
template<typename T>
void saveRawMultiBlock(MultiBlock3D<T> const& block, std::string fName, bool enforceUint=false);

/// Load the raw content of a MultiBlock3D efficiently from a Base64 encoded binary file.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is read, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This function is more efficient than the function loadBinaryBlock, but it is
 *  dependent on the structure of the multi-block. If data is written with saveRawMultiBlock()
 *  in a parallel program, the program which loads the data with loadRawMultiBlock() must
 *  have exactly the same number of cores.
 */
template<typename T>
void loadRawMultiBlock(MultiBlock3D<T>& block, std::string fName, bool enforceUint=false);

/// Flush the content of a Block3D into a generic C++ stream with space-separated ASCII words.
/** The content includes external scalars in the case of a BlockLattice3D. Only raw
 *  data is written, and structural information such as length and width of the
 *  block is lost.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This file format is not exact and should be used for data post-processing only,
 *  and not for checkpointing.
 */
template<typename T>
std::ostream& operator<<(std::ostream& ostr, Block3D<T> const& block);

/// Flush the content of a generic C++ stream with ASCII content into a Block3D.
/** The content includes external scalars in the case of a BlockLattice3D.
 *
 *  Index-ordering for this operation can be chosen through a call to
 *  global::IOpolicy().setIndexOrderingForStreams(IndexOrdering::OrderingT).
 *
 *  This file format is not exact and should be used for data post-processing only,
 *  and not for checkpointing.
 */
template<typename T>
std::istream& operator>>(std::istream& istr, Block3D<T>& block);

} // namespace plb

#endif  // SERIALIZER_IO_3D_H
