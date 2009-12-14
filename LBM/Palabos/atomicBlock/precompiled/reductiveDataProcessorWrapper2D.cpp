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

#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

/* *************** ReductiveBoxed Data Processor functionals ****************** */

template class ReductiveBoxProcessor2D<double>;
template class ReductiveBoxProcessorGenerator2D<double>;


template class ReductiveBoxProcessingFunctional2D<double>;

template class ReductiveBoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveBoxProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice );

template class ReductiveBoxProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        ReductiveBoxProcessingFunctional2D_S<double>& functional,
        Box2D domain, ScalarFieldBase2D<double>& field );

template class ReductiveBoxProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        ReductiveBoxProcessingFunctional2D_T<double,2>& functional,
        Box2D domain, TensorFieldBase2D<double,2>& field );

template class ReductiveBoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class ReductiveBoxProcessingFunctional2D_SS<double>;
template class ReductiveBoxProcessingFunctional2D_TT<double,2,2>;
template class ReductiveBoxProcessingFunctional2D_ST<double,2>;
template class ReductiveBoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class ReductiveBoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class ReductiveLatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class ReductiveScalarFieldBoxProcessingFunctional2D<double>;
template class ReductiveTensorFieldBoxProcessingFunctional2D<double,2>;


/* *************** ReductiveDotted Data Processor functionals ***************** */

template class ReductiveDotProcessor2D<double>;
template class ReductiveDotProcessorGenerator2D<double>;

template class ReductiveDotProcessingFunctional2D<double>;

template class ReductiveDotProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveDotProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>& functional,
        DotList2D const& domain, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice );

template class ReductiveDotProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        ReductiveDotProcessingFunctional2D_S<double>& functional,
        DotList2D const& domain, ScalarFieldBase2D<double>& field );

template class ReductiveDotProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        ReductiveDotProcessingFunctional2D_T<double,2>& functional,
        DotList2D const& domain, TensorFieldBase2D<double,2>& field );

template class ReductiveDotProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class ReductiveDotProcessingFunctional2D_SS<double>;
template class ReductiveDotProcessingFunctional2D_TT<double,2,2>;
template class ReductiveDotProcessingFunctional2D_ST<double,2>;
template class ReductiveDotProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class ReductiveDotProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class ReductiveLatticeDotProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class ReductiveScalarFieldDotProcessingFunctional2D<double>;
template class ReductiveTensorFieldDotProcessingFunctional2D<double,2>;

/* *************** Bounded ReductiveBoxed Data Processor functionals ****************** */

template class BoundedReductiveBoxProcessingFunctional2D<double>;
template class BoundedReductiveOneBlockProcessingFunctionalOperation2D<double>;

template class BoundedReductiveBoxProcessingFunctional2D_L<double, descriptors::D2Q9Descriptor>;
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedReductiveBoxProcessingFunctional2D_L<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain, BlockLatticeBase2D<double,descriptors::D2Q9Descriptor>& lattice,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional2D_S<double>;
template
void applyProcessingFunctional<double> (
        BoundedReductiveBoxProcessingFunctional2D_S<double>& functional,
        Box2D domain, ScalarFieldBase2D<double>& field,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional2D_T<double,2>;
template
void applyProcessingFunctional<double,2> (
        BoundedReductiveBoxProcessingFunctional2D_T<double,2>& functional,
        Box2D domain, TensorFieldBase2D<double,2>& field,
        plint boundaryWidth );

template class BoundedReductiveBoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>;
template class BoundedReductiveBoxProcessingFunctional2D_SS<double>;
template class BoundedReductiveBoxProcessingFunctional2D_TT<double,2,2>;
template class BoundedReductiveBoxProcessingFunctional2D_ST<double,2>;
template class BoundedReductiveBoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>;
template class BoundedReductiveBoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>;
template class BoundedReductiveLatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>;
template class BoundedReductiveScalarFieldBoxProcessingFunctional2D<double>;
template class BoundedReductiveTensorFieldBoxProcessingFunctional2D<double,2>;

}  // namespace plb
