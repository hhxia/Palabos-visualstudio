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

#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** ReductiveBoxed Data Processor functionals ****************** */

template class ReductiveBoxProcessor3D<double>;
template class ReductiveBoxProcessorGenerator3D<double>;


template class ReductiveBoxProcessingFunctional3D<double>;

template class ReductiveBoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveBoxProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice );

template class ReductiveBoxProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        ReductiveBoxProcessingFunctional3D_S<double>& functional,
        Box3D domain, ScalarFieldBase3D<double>& field );

template class ReductiveBoxProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        ReductiveBoxProcessingFunctional3D_T<double,3>& functional,
        Box3D domain, TensorFieldBase3D<double,3>& field );

template class ReductiveBoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class ReductiveBoxProcessingFunctional3D_SS<double>;
template class ReductiveBoxProcessingFunctional3D_TT<double,3,3>;
template class ReductiveBoxProcessingFunctional3D_ST<double,3>;
template class ReductiveBoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class ReductiveBoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class ReductiveLatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class ReductiveScalarFieldBoxProcessingFunctional3D<double>;
template class ReductiveTensorFieldBoxProcessingFunctional3D<double,3>;


/* *************** ReductiveDotted Data Processor functionals ***************** */

template class ReductiveDotProcessor3D<double>;
template class ReductiveDotProcessorGenerator3D<double>;

template class ReductiveDotProcessingFunctional3D<double>;

template class ReductiveDotProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveDotProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>& functional,
        DotList3D const& domain, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice );

template class ReductiveDotProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        ReductiveDotProcessingFunctional3D_S<double>& functional,
        DotList3D const& domain, ScalarFieldBase3D<double>& field );

template class ReductiveDotProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        ReductiveDotProcessingFunctional3D_T<double,3>& functional,
        DotList3D const& domain, TensorFieldBase3D<double,3>& field );

template class ReductiveDotProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class ReductiveDotProcessingFunctional3D_SS<double>;
template class ReductiveDotProcessingFunctional3D_TT<double,3,3>;
template class ReductiveDotProcessingFunctional3D_ST<double,3>;
template class ReductiveDotProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class ReductiveDotProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class ReductiveLatticeDotProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class ReductiveScalarFieldDotProcessingFunctional3D<double>;
template class ReductiveTensorFieldDotProcessingFunctional3D<double,3>;

/* *************** Bounded ReductiveBoxed Data Processor functionals ****************** */

template class BoundedReductiveBoxProcessingFunctional3D<double>;
template class BoundedReductiveOneBlockProcessingFunctionalOperation3D<double>;

template class BoundedReductiveBoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedReductiveBoxProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        BoundedReductiveBoxProcessingFunctional3D_S<double>& functional,
        Box3D domain, ScalarFieldBase3D<double>& field,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        BoundedReductiveBoxProcessingFunctional3D_T<double,3>& functional,
        Box3D domain, TensorFieldBase3D<double,3>& field,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class BoundedReductiveBoxProcessingFunctional3D_SS<double>;
template class BoundedReductiveBoxProcessingFunctional3D_TT<double,3,3>;
template class BoundedReductiveBoxProcessingFunctional3D_ST<double,3>;
template class BoundedReductiveBoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class BoundedReductiveBoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class BoundedReductiveLatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoundedReductiveScalarFieldBoxProcessingFunctional3D<double>;
template class BoundedReductiveTensorFieldBoxProcessingFunctional3D<double,3>;

}  // namespace plb
