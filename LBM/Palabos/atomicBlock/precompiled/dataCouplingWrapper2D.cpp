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

#include "atomicBlock/dataCouplingWrapper2D.h"
#include "atomicBlock/dataCouplingWrapper2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<double>(BoxProcessingFunctional2D<double>* functional,
                                       Box2D domain, std::vector<AtomicBlock2D<double>*> atomicBlocks);
template
void integrateProcessingFunctional<double>(BoxProcessingFunctional2D<double>* functional,
                                           Box2D domain, std::vector<AtomicBlock2D<double>*> atomicBlocks, plint level);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor, descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2 );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice2, plint level );

template
void applyProcessingFunctional<double> (
        BoxProcessingFunctional2D_SS<double>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2 );
template
void integrateProcessingFunctional<double> (
        BoxProcessingFunctional2D_SS<double>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2, plint level );

template
void applyProcessingFunctional<double,2,2> (
        BoxProcessingFunctional2D_TT<double,2,2>* functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2 );
template
void integrateProcessingFunctional<double,2,2> (
        BoxProcessingFunctional2D_TT<double,2,2>* functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2, plint level );

template
void applyProcessingFunctional<double,2> (
        BoxProcessingFunctional2D_ST<double,2>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2 );
template
void integrateProcessingFunctional<double,2> (
        BoxProcessingFunctional2D_ST<double,2>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        BoxProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        BoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        LatticeBoxProcessingFunctional2D<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        LatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain, std::vector<BlockLattice2D<double, descriptors::D2Q9Descriptor>*> lattices, plint level );

template
void applyProcessingFunctional<double> (
        ScalarFieldBoxProcessingFunctional2D<double>* functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields );
template
void integrateProcessingFunctional<double> (
        ScalarFieldBoxProcessingFunctional2D<double>* functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields, plint level );

template
void applyProcessingFunctional<double,2> (
        TensorFieldBoxProcessingFunctional2D<double,2>* functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields );
template
void integrateProcessingFunctional<double,2> (
        TensorFieldBoxProcessingFunctional2D<double,2>* functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields, plint level );



/* *************** Dotted Data Processor functionals ***************** */

template
void applyProcessingFunctional<double>(DotProcessingFunctional2D<double>* functional,
                                       DotList2D const& dotList, std::vector<AtomicBlock2D<double>*> atomicBlocks);
template
void integrateProcessingFunctional<double>(DotProcessingFunctional2D<double>* functional,
                                           DotList2D const& dotList, std::vector<AtomicBlock2D<double>*> atomicBlocks, plint level);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2 );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice2, plint level );

template
void applyProcessingFunctional<double> (
        DotProcessingFunctional2D_SS<double>* functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2 );
template
void integrateProcessingFunctional<double> (
        DotProcessingFunctional2D_SS<double>* functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2, plint level );

template
void applyProcessingFunctional<double,2,2> (
        DotProcessingFunctional2D_TT<double,2,2>* functional,
        DotList2D const& dotList,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2 );
template
void integrateProcessingFunctional<double,2,2> (
        DotProcessingFunctional2D_TT<double,2,2>* functional,
        DotList2D const& dotList,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2, plint level );

template
void applyProcessingFunctional<double,2> (
        DotProcessingFunctional2D_ST<double,2>* functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2 );
template
void integrateProcessingFunctional<double,2> (
        DotProcessingFunctional2D_ST<double,2>* functional,
        DotList2D const& dotList,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        DotProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        DotProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        DotProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>* functional,
        DotList2D const& dotList,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        LatticeDotProcessingFunctional2D<double,descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        LatticeDotProcessingFunctional2D<double, descriptors::D2Q9Descriptor>* functional,
        DotList2D const& dotList, std::vector<BlockLattice2D<double, descriptors::D2Q9Descriptor>*> lattices, plint level );

template
void applyProcessingFunctional<double> (
        ScalarFieldDotProcessingFunctional2D<double>* functional,
        DotList2D const& dotList, std::vector<ScalarField2D<double>*> fields );
template
void integrateProcessingFunctional<double> (
        ScalarFieldDotProcessingFunctional2D<double>* functional,
        DotList2D const& dotList, std::vector<ScalarField2D<double>*> fields, plint level );

template
void applyProcessingFunctional<double,2> (
        TensorFieldDotProcessingFunctional2D<double,2>* functional,
        DotList2D const& dotList, std::vector<TensorField2D<double,2>*> fields );
template
void integrateProcessingFunctional<double,2> (
        TensorFieldDotProcessingFunctional2D<double,2>* functional,
        DotList2D const& dotList, std::vector<TensorField2D<double,2>*> fields, plint level );


/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedCoupledBlocksProcessingFunctionalOperation2D<double>;

template
void applyProcessingFunctional<double>(BoundedBoxProcessingFunctional2D<double>* functional,
                                       Box2D domain, std::vector<AtomicBlock2D<double>*> atomicBlocks,
                                       plint boundaryWidth);
template
void integrateProcessingFunctional<double>(BoundedBoxProcessingFunctional2D<double>* functional,
                                           Box2D domain, std::vector<AtomicBlock2D<double>*> atomicBlocks,
                                           plint boundaryWidth, plint level);

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_LL<double,descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_LL<double, descriptors::D2Q9Descriptor,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice1,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double> (
        BoundedBoxProcessingFunctional2D_SS<double>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedBoxProcessingFunctional2D_SS<double>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        ScalarField2D<double>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,2,2> (
        BoundedBoxProcessingFunctional2D_TT<double,2,2>* functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double,2,2> (
        BoundedBoxProcessingFunctional2D_TT<double,2,2>* functional,
        Box2D domain,
        TensorField2D<double,2>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,2> (
        BoundedBoxProcessingFunctional2D_ST<double,2>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double,2> (
        BoundedBoxProcessingFunctional2D_ST<double,2>* functional,
        Box2D domain,
        ScalarField2D<double>& field1,
        TensorField2D<double,2>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_LS<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedBoxProcessingFunctional2D_LS<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        ScalarField2D<double>& field,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        BoundedBoxProcessingFunctional2D_LT<double,descriptors::D2Q9Descriptor,2>* functional,
        Box2D domain,
        BlockLattice2D<double,descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor,2> (
        BoundedBoxProcessingFunctional2D_LT<double, descriptors::D2Q9Descriptor,2>* functional,
        Box2D domain,
        BlockLattice2D<double, descriptors::D2Q9Descriptor>& lattice,
        TensorField2D<double,2>& field,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedLatticeBoxProcessingFunctional2D<double,descriptors::D2Q9Descriptor>* functional,
        Box2D domain, std::vector<BlockLattice2D<double,descriptors::D2Q9Descriptor>*> lattices,
        plint boundaryWidth);
template
void integrateProcessingFunctional<double, descriptors::D2Q9Descriptor> (
        BoundedLatticeBoxProcessingFunctional2D<double, descriptors::D2Q9Descriptor>* functional,
        Box2D domain, std::vector<BlockLattice2D<double, descriptors::D2Q9Descriptor>*> lattices,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double> (
        BoundedScalarFieldBoxProcessingFunctional2D<double>* functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields, plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedScalarFieldBoxProcessingFunctional2D<double>* functional,
        Box2D domain, std::vector<ScalarField2D<double>*> fields, plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,2> (
        BoundedTensorFieldBoxProcessingFunctional2D<double,2>* functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields, plint boundaryWidth );
template
void integrateProcessingFunctional<double,2> (
        BoundedTensorFieldBoxProcessingFunctional2D<double,2>* functional,
        Box2D domain, std::vector<TensorField2D<double,2>*> fields, plint boundaryWidth, plint level );


}  // namespace plb
