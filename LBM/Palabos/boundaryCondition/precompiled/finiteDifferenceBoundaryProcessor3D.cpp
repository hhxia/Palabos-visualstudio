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

#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.hh"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/dataProcessor3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

using namespace descriptors;

template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          0,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 0,-1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          1,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 1,-1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,1>;
template class PlaneFdBoundaryProcessor3D<double, descriptors::D3Q19Descriptor,          2,-1>;
template class PlaneFdBoundaryProcessorGenerator3D<double, descriptors::D3Q19Descriptor, 2,-1>;

template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 0, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 0, 1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 1, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 1, 1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2,-1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2,-1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2,-1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2,-1, 1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2, 1,-1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2, 1,-1>;
template class OuterVelocityEdgeProcessor3D          <double, descriptors::D3Q19Descriptor, 2, 1, 1>;
template class OuterVelocityEdgeProcessorGenerator3D <double, descriptors::D3Q19Descriptor, 2, 1, 1>;

template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1,-1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1,-1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1,-1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1,-1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1, 1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1, 1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor, -1, 1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor, -1, 1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1,-1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1,-1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1,-1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1,-1, 1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1, 1,-1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1, 1,-1>;
template class OuterVelocityCornerProcessor3D          <double, descriptors::D3Q19Descriptor,  1, 1, 1>;
template class OuterVelocityCornerProcessorGenerator3D <double, descriptors::D3Q19Descriptor,  1, 1, 1>;

}
