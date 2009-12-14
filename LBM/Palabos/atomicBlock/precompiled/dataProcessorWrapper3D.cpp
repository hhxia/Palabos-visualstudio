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

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/dataProcessorWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessor3D<double>;
template class BoxProcessorGenerator3D<double>;


template class BoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain, BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, plint level );

template class BoxProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        BoxProcessingFunctional3D_S<double>* functional,
        Box3D domain, ScalarFieldBase3D<double>& field );
template
void integrateProcessingFunctional<double> (
        BoxProcessingFunctional3D_S<double>* functional,
        Box3D domain, ScalarFieldBase3D<double>& field, plint level );

template class BoxProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        BoxProcessingFunctional3D_T<double,3>* functional,
        Box3D domain, TensorFieldBase3D<double,3>& field );
template
void integrateProcessingFunctional<double,3> (
        BoxProcessingFunctional3D_T<double,3>* functional,
        Box3D domain, TensorFieldBase3D<double,3>& field, plint level );

template class BoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class BoxProcessingFunctional3D_SS<double>;
template class BoxProcessingFunctional3D_TT<double,3,3>;
template class BoxProcessingFunctional3D_ST<double,3>;
template class BoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class BoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class LatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class ScalarFieldBoxProcessingFunctional3D<double>;
template class TensorFieldBoxProcessingFunctional3D<double,3>;


/* *************** Dotted Data Processor functionals ***************** */

template class DotProcessor3D<double>;
template class DotProcessorGenerator3D<double>;

template class DotProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList, BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, plint level );

template class DotProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        DotProcessingFunctional3D_S<double>* functional,
        DotList3D const& dotList, ScalarFieldBase3D<double>& field );
template
void integrateProcessingFunctional<double> (
        DotProcessingFunctional3D_S<double>* functional,
        DotList3D const& dotList, ScalarFieldBase3D<double>& field, plint level );

template class DotProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        DotProcessingFunctional3D_T<double,3>* functional,
        DotList3D const& dotList, TensorFieldBase3D<double,3>& field );
template
void integrateProcessingFunctional<double,3> (
        DotProcessingFunctional3D_T<double,3>* functional,
        DotList3D const& dotList, TensorFieldBase3D<double,3>& field, plint level );

template class DotProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class DotProcessingFunctional3D_SS<double>;
template class DotProcessingFunctional3D_TT<double,3,3>;
template class DotProcessingFunctional3D_ST<double,3>;
template class DotProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class DotProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class LatticeDotProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class ScalarFieldDotProcessingFunctional3D<double>;
template class TensorFieldDotProcessingFunctional3D<double,3>;


/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional3D<double>;

template class BoundedBoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_L<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain, BlockLatticeBase3D<double,descriptors::D3Q19Descriptor>& lattice, plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_L<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain, BlockLatticeBase3D<double, descriptors::D3Q19Descriptor>& lattice, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional3D_S<double>;
template
void applyProcessingFunctional<double> (
        BoundedBoxProcessingFunctional3D_S<double>* functional,
        Box3D domain, ScalarFieldBase3D<double>& field, plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedBoxProcessingFunctional3D_S<double>* functional,
        Box3D domain, ScalarFieldBase3D<double>& field, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional3D_T<double,3>;
template
void applyProcessingFunctional<double,3> (
        BoundedBoxProcessingFunctional3D_T<double,3>* functional,
        Box3D domain, TensorFieldBase3D<double,3>& field, plint boundaryWidth );
template
void integrateProcessingFunctional<double,3> (
        BoundedBoxProcessingFunctional3D_T<double,3>* functional,
        Box3D domain, TensorFieldBase3D<double,3>& field, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>;
template class BoundedBoxProcessingFunctional3D_SS<double>;
template class BoundedBoxProcessingFunctional3D_TT<double,3,3>;
template class BoundedBoxProcessingFunctional3D_ST<double,3>;
template class BoundedBoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>;
template class BoundedBoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>;
template class BoundedLatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoundedScalarFieldBoxProcessingFunctional3D<double>;
template class BoundedTensorFieldBoxProcessingFunctional3D<double,3>;



}  // namespace plb
