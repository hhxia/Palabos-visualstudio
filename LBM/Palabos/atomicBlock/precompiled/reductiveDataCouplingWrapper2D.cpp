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

#include "atomicBlock/reductiveDataCouplingWrapper2D.h"
#include "atomicBlock/reductiveDataCouplingWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<double>(ReductiveBoxProcessingFunctional2D<double>& functional,
                                       Box2D domain, std::vector<AtomicBlock2D<double>*> objects);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        ReductiveBoxProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2 );
template
void applyProcessingFunctional<double> (
        ReductiveBoxProcessingFunctional2D_SS<double>& functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2 );
template
void applyProcessingFunctional<double,2,2> (
        ReductiveBoxProcessingFunctional2D_TT<double,2,2>& functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2 );
template
void applyProcessingFunctional<double,2> (
        ReductiveBoxProcessingFunctional2D_ST<double,2>& functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2 );
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveBoxProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field );
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        ReductiveBoxProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveLatticeBoxProcessingFunctional2D<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices );

template
void applyProcessingFunctional<double> (
        ReductiveScalarFieldBoxProcessingFunctional2D<double>& functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields );

template
void applyProcessingFunctional<double,2> (
        ReductiveTensorFieldBoxProcessingFunctional2D<double,2>& functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields );



/* *************** Dotted Data Processor functionals ***************** */

template
void applyProcessingFunctional<double>(ReductiveDotProcessingFunctional2D<double>& functional,
                                       DotList2D const& dotList, std::vector<AtomicBlock2D<double>*> objects);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        ReductiveDotProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>& functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2 );
template
void applyProcessingFunctional<double> (
        ReductiveDotProcessingFunctional2D_SS<double>& functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2 );
template
void applyProcessingFunctional<double,2,2> (
        ReductiveDotProcessingFunctional2D_TT<double,2,2>& functional,
        DotList2D const& dotList,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2 );
template
void applyProcessingFunctional<double,2> (
        ReductiveDotProcessingFunctional2D_ST<double,2>& functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2 );
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveDotProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>& functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field );
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        ReductiveDotProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>& functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        ReductiveLatticeDotProcessingFunctional2D<double,descriptors::D2Q9Descriptor>& functional,
        DotList2D const& dotList, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices );

template
void applyProcessingFunctional<double> (
        ReductiveScalarFieldDotProcessingFunctional2D<double>& functional,
        DotList2D const& dotList, std::vector<ScalarField2D<double>*> fields );

template
void applyProcessingFunctional<double,2> (
        ReductiveTensorFieldDotProcessingFunctional2D<double,2>& functional,
        DotList2D const& dotList, std::vector<TensorField2D<double,2>*> fields );



/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D<double>;

template
void applyProcessingFunctional<double>(BoundedReductiveBoxProcessingFunctional2D<double>& functional,
                                       Box2D domain, std::vector<AtomicBlock2D<double>*> objects,
                                       plint boundaryWidth);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        BoundedReductiveBoxProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double> (
        BoundedReductiveBoxProcessingFunctional2D_SS<double>& functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double,2,2> (
        BoundedReductiveBoxProcessingFunctional2D_TT<double,2,2>& functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double,2> (
        BoundedReductiveBoxProcessingFunctional2D_ST<double,2>& functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedReductiveBoxProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field,
        plint boundaryWidth);
template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        BoundedReductiveBoxProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>& functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field,
        plint boundaryWidth);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedReductiveLatticeBoxProcessingFunctional2D<double,descriptors::D2Q9Descriptor>& functional,
        Box2D domain, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices,
        plint boundaryWidth);

template
void applyProcessingFunctional<double> (
        BoundedReductiveScalarFieldBoxProcessingFunctional2D<double>& functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields, plint boundaryWidth );

template
void applyProcessingFunctional<double,2> (
        BoundedReductiveTensorFieldBoxProcessingFunctional2D<double,2>& functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields, plint boundaryWidth );


}  // namespace plb
