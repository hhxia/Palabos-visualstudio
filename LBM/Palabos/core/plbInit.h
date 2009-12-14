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
 * LB initialisation routine -- header file.
 */
#ifndef PLB_INIT_H
#define PLB_INIT_H

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"

namespace plb {

inline void plbInit(int *argc, char ***argv, bool verbous=false) {
    global::mpi().init(argc, argv, verbous);
}

}  // namespace plb

#endif
