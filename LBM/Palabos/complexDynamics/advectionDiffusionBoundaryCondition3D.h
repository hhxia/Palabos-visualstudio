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

/* Main author: Orestis Malaspinas
 */

/** \file A helper for initialising 3D boundaries -- header file.  */

#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/blockLatticeBase3D.h"
#include "complexDynamics/advectionDiffusionDynamics.h"

#include <vector>
#include <list>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class OnLatticeAdvectionDiffusionBoundaryCondition3D {
public:
    virtual ~OnLatticeAdvectionDiffusionBoundaryCondition3D() { }

    // 3D boundary condition for temperature:
    virtual void addTemperatureBoundary0N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureBoundary0P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureBoundary1N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureBoundary1P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureBoundary2N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureBoundary2P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;

    virtual void addTemperatureEdge0NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge0NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge0PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge0PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge1NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge1NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge1PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge1PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge2NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge2NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge2PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureEdge2PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;

    virtual void addTemperatureCornerNNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerNNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerNPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerNPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerPNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerPNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerPPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;
    virtual void addTemperatureCornerPPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice ) =0;

    void setTemperatureConditionOnBlockBoundaries(BlockLatticeBase3D<T,Descriptor>& lattice);
};


//////  Factory function for Regularized Thermal BC

template<typename T, template<typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>*
    createLocalAdvectionDiffusionBoundaryCondition3D(
            BlockLatticeBase3D<T,Descriptor>& block);

} //namespace plb


#endif
