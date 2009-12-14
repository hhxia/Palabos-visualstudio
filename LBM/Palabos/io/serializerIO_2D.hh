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

#ifndef SERIALIZER_IO_2D_HH
#define SERIALIZER_IO_2D_HH

#include "io/serializerIO_2D.h"
#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "parallelism/mpiManager.h"
#include "multiBlock/multiBlockSerializer2D.h"
#include <istream>
#include <ostream>
#include <fstream>

namespace plb {

template<typename T>
void saveBinaryBlock(Block2D<T> const& block, std::string fName, bool enforceUint) {
    std::ofstream* ostr = 0;
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fName.c_str());
        plbMainProcWarning( !(*ostr), std::string("Could not open binary file ")+
                                      fName+std::string(" for saving") );
    }
    serializerToBase64Stream<T>( block.getBlockSerializer(block.getBoundingBox(),
                                 global::IOpolicy().getIndexOrderingForStreams()), *ostr, enforceUint );
    delete ostr;
}

template<typename T>
void loadBinaryBlock(Block2D<T>& block, std::string fName, bool enforceUint) {
    std::ifstream* istr = 0;
    if (global::mpi().isMainProcessor()) {
        istr = new std::ifstream(fName.c_str());
        plbMainProcWarning( !(*istr), std::string("Could not open binary file ")+
                                      fName+std::string(" for reading") );
    }
    base64StreamToUnSerializer<T> (
            *istr, block.getBlockUnSerializer (
                block.getBoundingBox(),
                global::IOpolicy().getIndexOrderingForStreams()),
            enforceUint );
    delete istr;
}

template<typename T>
std::ostream& operator<<(std::ostream& ostr, Block2D<T> const& block) {
    plint numDigits = 0; // Number of digits is handled by iostream manipulators,
                        // as opposed to being imposed here.
    serializerToAsciiStream<T>( block.getBlockSerializer (
                                    block.getBoundingBox(),
                                    global::IOpolicy().getIndexOrderingForStreams() ),
                                ostr, numDigits );
    return ostr;
}

template<typename T>
std::istream& operator>>(std::istream& istr, Block2D<T>& block) {
    asciiStreamToUnSerializer<T>( istr,
                                  block.getBlockUnSerializer (
                                      block.getBoundingBox(),
                                      global::IOpolicy().getIndexOrderingForStreams()
                                  ) );
    return istr;
}

} // namespace plb

#endif
