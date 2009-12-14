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
 * Definition of external fields for a lattice descriptor -- header file
 */
#ifndef EXTERNAL_FIELDS_H
#define EXTERNAL_FIELDS_H

#include "core/globalDefs.h"

namespace plb {

namespace descriptors {

struct NoExternalField {
    static const int numScalars = 0;
    static const int numSpecies = 0;
    static const int forceBeginsAt = 0;
    static const int sizeOfForce   = 0;
};

struct NoExternalFieldBase {
    typedef NoExternalField ExternalField;
};

struct Force2dDescriptor {
    static const int numScalars = 2;
    static const int numSpecies = 1;
    static const int forceBeginsAt = 0;
    static const int sizeOfForce   = 2;
};

struct Force2dDescriptorBase {
    typedef Force2dDescriptor ExternalField;
};

struct Force3dDescriptor {
    static const int numScalars = 3;
    static const int numSpecies = 1;
    static const int forceBeginsAt = 0;
    static const int sizeOfForce   = 3;
};

struct Force3dDescriptorBase {
    typedef Force3dDescriptor ExternalField;
};

}  // namespace descriptors

}  // namespace plb

#endif
