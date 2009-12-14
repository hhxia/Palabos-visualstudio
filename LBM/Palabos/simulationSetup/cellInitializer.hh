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
 * Helper functions for domain initialization -- header file.
 */
#ifndef CELL_INITIALIZER_HH
#define CELL_INITIALIZER_HH

#include "simulationSetup/cellInitializer.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(Cell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity) {
    Array<T,Descriptor<T>::d> j;
    VectorTemplate<T,Descriptor>::multiplyByScalar(velocity, density, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(density);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

}  // namespace plb

#endif  // CELL_INITIALIZER_HH
