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

#include "atomicBlock/dataCouplingWrapper3D.h"
#include "atomicBlock/dataCouplingWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<double>(BoxProcessingFunctional3D<double>* functional,
                                       Box3D domain, std::vector<AtomicBlock3D<double>*> atomicBlocks);
template
void integrateProcessingFunctional<double>(BoxProcessingFunctional3D<double>* functional,
                                           Box3D domain, std::vector<AtomicBlock3D<double>*> atomicBlocks, plint level);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2 );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice2, plint level );

template
void applyProcessingFunctional<double> (
        BoxProcessingFunctional3D_SS<double>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2 );
template
void integrateProcessingFunctional<double> (
        BoxProcessingFunctional3D_SS<double>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2, plint level );

template
void applyProcessingFunctional<double,3,3> (
        BoxProcessingFunctional3D_TT<double,3,3>* functional,
        Box3D domain,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2 );
template
void integrateProcessingFunctional<double,3,3> (
        BoxProcessingFunctional3D_TT<double,3,3>* functional,
        Box3D domain,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2, plint level );

template
void applyProcessingFunctional<double,3> (
        BoxProcessingFunctional3D_ST<double,3>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2 );
template
void integrateProcessingFunctional<double,3> (
        BoxProcessingFunctional3D_ST<double,3>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        BoxProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        BoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        LatticeBoxProcessingFunctional3D<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain, std::vector<BlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        LatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain, std::vector<BlockLattice3D<double, descriptors::D3Q19Descriptor>*> lattices, plint level );

template
void applyProcessingFunctional<double> (
        ScalarFieldBoxProcessingFunctional3D<double>* functional,
        Box3D domain, std::vector<ScalarField3D<double>*> fields );
template
void integrateProcessingFunctional<double> (
        ScalarFieldBoxProcessingFunctional3D<double>* functional,
        Box3D domain, std::vector<ScalarField3D<double>*> fields, plint level );

template
void applyProcessingFunctional<double,3> (
        TensorFieldBoxProcessingFunctional3D<double,3>* functional,
        Box3D domain, std::vector<TensorField3D<double,3>*> fields );
template
void integrateProcessingFunctional<double,3> (
        TensorFieldBoxProcessingFunctional3D<double,3>* functional,
        Box3D domain, std::vector<TensorField3D<double,3>*> fields, plint level );



/* *************** Dotted Data Processor functionals ***************** */

template
void applyProcessingFunctional<double>(DotProcessingFunctional3D<double>* functional,
                                       DotList3D const& dotList, std::vector<AtomicBlock3D<double>*> atomicBlocks);
template
void integrateProcessingFunctional<double>(DotProcessingFunctional3D<double>* functional,
                                           DotList3D const& dotList, std::vector<AtomicBlock3D<double>*> atomicBlocks, plint level);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2 );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice2, plint level );

template
void applyProcessingFunctional<double> (
        DotProcessingFunctional3D_SS<double>* functional,
        DotList3D const& dotList,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2 );
template
void integrateProcessingFunctional<double> (
        DotProcessingFunctional3D_SS<double>* functional,
        DotList3D const& dotList,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2, plint level );

template
void applyProcessingFunctional<double,3,3> (
        DotProcessingFunctional3D_TT<double,3,3>* functional,
        DotList3D const& dotList,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2 );
template
void integrateProcessingFunctional<double,3,3> (
        DotProcessingFunctional3D_TT<double,3,3>* functional,
        DotList3D const& dotList,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2, plint level );

template
void applyProcessingFunctional<double,3> (
        DotProcessingFunctional3D_ST<double,3>* functional,
        DotList3D const& dotList,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2 );
template
void integrateProcessingFunctional<double,3> (
        DotProcessingFunctional3D_ST<double,3>* functional,
        DotList3D const& dotList,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        DotProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        DotProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        DotProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>* functional,
        DotList3D const& dotList,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        LatticeDotProcessingFunctional3D<double,descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList, std::vector<BlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        LatticeDotProcessingFunctional3D<double, descriptors::D3Q19Descriptor>* functional,
        DotList3D const& dotList, std::vector<BlockLattice3D<double, descriptors::D3Q19Descriptor>*> lattices, plint level );

template
void applyProcessingFunctional<double> (
        ScalarFieldDotProcessingFunctional3D<double>* functional,
        DotList3D const& dotList, std::vector<ScalarField3D<double>*> fields );
template
void integrateProcessingFunctional<double> (
        ScalarFieldDotProcessingFunctional3D<double>* functional,
        DotList3D const& dotList, std::vector<ScalarField3D<double>*> fields, plint level );

template
void applyProcessingFunctional<double,3> (
        TensorFieldDotProcessingFunctional3D<double,3>* functional,
        DotList3D const& dotList, std::vector<TensorField3D<double,3>*> fields );
template
void integrateProcessingFunctional<double,3> (
        TensorFieldDotProcessingFunctional3D<double,3>* functional,
        DotList3D const& dotList, std::vector<TensorField3D<double,3>*> fields, plint level );


/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedCoupledBlocksProcessingFunctionalOperation3D<double>;

template
void applyProcessingFunctional<double>(BoundedBoxProcessingFunctional3D<double>* functional,
                                       Box3D domain, std::vector<AtomicBlock3D<double>*> atomicBlocks,
                                       plint boundaryWidth);
template
void integrateProcessingFunctional<double>(BoundedBoxProcessingFunctional3D<double>* functional,
                                           Box3D domain, std::vector<AtomicBlock3D<double>*> atomicBlocks,
                                           plint boundaryWidth, plint level);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_LL<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice1,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double> (
        BoundedBoxProcessingFunctional3D_SS<double>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedBoxProcessingFunctional3D_SS<double>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        ScalarField3D<double>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,3,3> (
        BoundedBoxProcessingFunctional3D_TT<double,3,3>* functional,
        Box3D domain,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double,3,3> (
        BoundedBoxProcessingFunctional3D_TT<double,3,3>* functional,
        Box3D domain,
        TensorField3D<double,3>& field1,
        TensorField3D<double,3>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,3> (
        BoundedBoxProcessingFunctional3D_ST<double,3>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double,3> (
        BoundedBoxProcessingFunctional3D_ST<double,3>* functional,
        Box3D domain,
        ScalarField3D<double>& field1,
        TensorField3D<double,3>& field2,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedBoxProcessingFunctional3D_LS<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        ScalarField3D<double>& field,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        BoundedBoxProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>* functional,
        Box3D domain,
        BlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field,
        plint boundaryWidth );
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        BoundedBoxProcessingFunctional3D_LT<double, descriptors::D3Q19Descriptor,3>* functional,
        Box3D domain,
        BlockLattice3D<double, descriptors::D3Q19Descriptor>& lattice,
        TensorField3D<double,3>& field,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedLatticeBoxProcessingFunctional3D<double,descriptors::D3Q19Descriptor>* functional,
        Box3D domain, std::vector<BlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices,
        plint boundaryWidth);
template
void integrateProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedLatticeBoxProcessingFunctional3D<double, descriptors::D3Q19Descriptor>* functional,
        Box3D domain, std::vector<BlockLattice3D<double, descriptors::D3Q19Descriptor>*> lattices,
        plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double> (
        BoundedScalarFieldBoxProcessingFunctional3D<double>* functional,
        Box3D domain, std::vector<ScalarField3D<double>*> fields, plint boundaryWidth );
template
void integrateProcessingFunctional<double> (
        BoundedScalarFieldBoxProcessingFunctional3D<double>* functional,
        Box3D domain, std::vector<ScalarField3D<double>*> fields, plint boundaryWidth, plint level );

template
void applyProcessingFunctional<double,3> (
        BoundedTensorFieldBoxProcessingFunctional3D<double,3>* functional,
        Box3D domain, std::vector<TensorField3D<double,3>*> fields, plint boundaryWidth );
template
void integrateProcessingFunctional<double,3> (
        BoundedTensorFieldBoxProcessingFunctional3D<double,3>* functional,
        Box3D domain, std::vector<TensorField3D<double,3>*> fields, plint boundaryWidth, plint level );


}  // namespace plb
