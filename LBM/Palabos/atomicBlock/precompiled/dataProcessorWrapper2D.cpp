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

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/dataProcessorWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessor2D<double>;
template class BoxProcessorGenerator2D<double>;


template class BoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain, BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, plint level );

template class BoxProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        BoxProcessingFunctional2D_S<double>* functional,
        Box2D domain, ScalarFieldBase2D<double>& field );
template
void integrateProcessingFunctional<double> (
        BoxProcessingFunctional2D_S<double>* functional,
        Box2D domain, ScalarFieldBase2D<double>& field, plint level );

template class BoxProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        BoxProcessingFunctional2D_T<double,2>* functional,
        Box2D domain, TensorFieldBase2D<double,2>& field );
template
void integrateProcessingFunctional<double,2> (
        BoxProcessingFunctional2D_T<double,2>* functional,
        Box2D domain, TensorFieldBase2D<double,2>& field, plint level );

template class BoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class BoxProcessingFunctional2D_SS<double>;
template class BoxProcessingFunctional2D_TT<double,2,2>;
template class BoxProcessingFunctional2D_ST<double,2>;
template class BoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class BoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class LatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class ScalarFieldBoxProcessingFunctional2D<double>;
template class TensorFieldBoxProcessingFunctional2D<double,2>;


/* *************** Dotted Data Processor functionals ***************** */

template class DotProcessor2D<double>;
template class DotProcessorGenerator2D<double>;

template class DotProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList, BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, plint level );

template class DotProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        DotProcessingFunctional2D_S<double>* functional,
        DotList2D const& dotList, ScalarFieldBase2D<double>& field );
template
void integrateProcessingFunctional<double> (
        DotProcessingFunctional2D_S<double>* functional,
        DotList2D const& dotList, ScalarFieldBase2D<double>& field, plint level );

template class DotProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        DotProcessingFunctional2D_T<double,2>* functional,
        DotList2D const& dotList, TensorFieldBase2D<double,2>& field );
template
void integrateProcessingFunctional<double,2> (
        DotProcessingFunctional2D_T<double,2>* functional,
        DotList2D const& dotList, TensorFieldBase2D<double,2>& field, plint level );

template class DotProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class DotProcessingFunctional2D_SS<double>;
template class DotProcessingFunctional2D_TT<double,2,2>;
template class DotProcessingFunctional2D_ST<double,2>;
template class DotProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class DotProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class LatticeDotProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class ScalarFieldDotProcessingFunctional2D<double>;
template class TensorFieldDotProcessingFunctional2D<double,2>;


/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional2D<double>;

template class BoundedBoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice, plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain, BlockLatticeBase2D<double, descriptors::D2Q9Descriptor>& lattice, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        BoundedBoxProcessingFunctional2D_S<double>* functional,
        Box2D domain, ScalarFieldBase2D<double>& field, plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedBoxProcessingFunctional2D_S<double>* functional,
        Box2D domain, ScalarFieldBase2D<double>& field, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        BoundedBoxProcessingFunctional2D_T<double,2>* functional,
        Box2D domain, TensorFieldBase2D<double,2>& field, plint boundaryWidth );
template
void integrateProcessingFunctional<double,2> (
        BoundedBoxProcessingFunctional2D_T<double,2>* functional,
        Box2D domain, TensorFieldBase2D<double,2>& field, plint boundaryWidth, plint level );

template class BoundedBoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class BoundedBoxProcessingFunctional2D_SS<double>;
template class BoundedBoxProcessingFunctional2D_TT<double,2,2>;
template class BoundedBoxProcessingFunctional2D_ST<double,2>;
template class BoundedBoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class BoundedBoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class BoundedLatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class BoundedScalarFieldBoxProcessingFunctional2D<double>;
template class BoundedTensorFieldBoxProcessingFunctional2D<double,2>;



}  // namespace plb
