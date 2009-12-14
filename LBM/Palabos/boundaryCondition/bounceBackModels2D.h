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
 * BounceBack dynamics models in 2D -- header file.
 */
#ifndef BOUNCE_BACK_MODELS_2D_H
#define BOUNCE_BACK_MODELS_2D_H

#include "boundaryCondition/bounceBackModels.h"
#include "simulationSetup/latticeInitializer2D.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase2D<T,Descriptor>& lattice, Box2D boundingBox,
        DomainFunctional2D* domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase2D<T,Descriptor>& lattice, DotList2D const& dotList );

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_2D_H
