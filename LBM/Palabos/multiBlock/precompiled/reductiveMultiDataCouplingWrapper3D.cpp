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

#include "multiBlock/reductiveMultiDataCouplingWrapper3D.h"
#include "multiBlock/reductiveMultiDataCouplingWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** Boxed Data Processor functionals ****************** */

template
void applyProcessingFunctional<double>(ReductiveBoxProcessingFunctional3D<double>& functional,
                                       Box3D domain,
                                       std::vector<MultiBlock3D<double>*> multiBlocks);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        ReductiveBoxProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2 );

template
void applyProcessingFunctional<double> (
        ReductiveBoxProcessingFunctional3D_SS<double>& functional,
        Box3D domain,
        MultiScalarField3D<double>& field1,
        MultiScalarField3D<double>& field2 );

template
void applyProcessingFunctional<double,3,3> (
        ReductiveBoxProcessingFunctional3D_TT<double,3,3>& functional,
        Box3D domain,
        MultiTensorField3D<double,3>& field1,
        MultiTensorField3D<double,3>& field2 );

template
void applyProcessingFunctional<double,3> (
        ReductiveBoxProcessingFunctional3D_ST<double,3>& functional,
        Box3D domain,
        MultiScalarField3D<double>& field1,
        MultiTensorField3D<double,3>& field2 );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveBoxProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& field );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        ReductiveBoxProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiTensorField3D<double,3>& field );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveLatticeBoxProcessingFunctional3D<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain, std::vector<MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices );

template
void applyProcessingFunctional<double> (
        ReductiveScalarFieldBoxProcessingFunctional3D<double>& functional,
        Box3D domain, std::vector<MultiScalarField3D<double>*> fields );

template
void applyProcessingFunctional<double,3> (
        ReductiveTensorFieldBoxProcessingFunctional3D<double,3>& functional,
        Box3D domain, std::vector<MultiTensorField3D<double,3>*> fields );





/* *************** Dotted Data Processor functionals ***************** */

template
void applyProcessingFunctional<double>(ReductiveDotProcessingFunctional3D<double>& functional,
                                       DotList3D const& dotList,
                                       std::vector<MultiBlock3D<double>*> multiBlocks);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        ReductiveDotProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>& functional,
        DotList3D const& dotList,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2 );

template
void applyProcessingFunctional<double> (
        ReductiveDotProcessingFunctional3D_SS<double>& functional,
        DotList3D const& dotList,
        MultiScalarField3D<double>& field1,
        MultiScalarField3D<double>& field2 );

template
void applyProcessingFunctional<double,3,3> (
        ReductiveDotProcessingFunctional3D_TT<double,3,3>& functional,
        DotList3D const& dotList,
        MultiTensorField3D<double,3>& field1,
        MultiTensorField3D<double,3>& field2 );

template
void applyProcessingFunctional<double,3> (
        ReductiveDotProcessingFunctional3D_ST<double,3>& functional,
        DotList3D const& dotList,
        MultiScalarField3D<double>& field1,
        MultiTensorField3D<double,3>& field2 );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveDotProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>& functional,
        DotList3D const& dotList,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& field );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        ReductiveDotProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>& functional,
        DotList3D const& dotList,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiTensorField3D<double,3>& field );

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        ReductiveLatticeDotProcessingFunctional3D<double,descriptors::D3Q19Descriptor>& functional,
        DotList3D const& dotList, std::vector<MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices );

template
void applyProcessingFunctional<double> (
        ReductiveScalarFieldDotProcessingFunctional3D<double>& functional,
        DotList3D const& dotList, std::vector<MultiScalarField3D<double>*> fields );

template
void applyProcessingFunctional<double,3> (
        ReductiveTensorFieldDotProcessingFunctional3D<double,3>& functional,
        DotList3D const& dotList, std::vector<MultiTensorField3D<double,3>*> fields );



/* *************** Bounded Boxed Multi Data Processor functionals ****************** */

template class BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D<double>;

template
void applyProcessingFunctional<double>(BoundedReductiveBoxProcessingFunctional3D<double>& functional,
                                       Box3D domain, std::vector<MultiBlock3D<double>*> objects,
                                       plint boundaryWidth);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor> (
        BoundedReductiveBoxProcessingFunctional3D_LL<double,descriptors::D3Q19Descriptor,descriptors::D3Q19Descriptor>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice1,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double> (
        BoundedReductiveBoxProcessingFunctional3D_SS<double>& functional,
        Box3D domain,
        MultiScalarField3D<double>& field1,
        MultiScalarField3D<double>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double,3,3> (
        BoundedReductiveBoxProcessingFunctional3D_TT<double,3,3>& functional,
        Box3D domain,
        MultiTensorField3D<double,3>& field1,
        MultiTensorField3D<double,3>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double,3> (
        BoundedReductiveBoxProcessingFunctional3D_ST<double,3>& functional,
        Box3D domain,
        MultiScalarField3D<double>& field1,
        MultiTensorField3D<double,3>& field2,
        plint boundaryWidth);
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedReductiveBoxProcessingFunctional3D_LS<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiScalarField3D<double>& field,
        plint boundaryWidth);
template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor,3> (
        BoundedReductiveBoxProcessingFunctional3D_LT<double,descriptors::D3Q19Descriptor,3>& functional,
        Box3D domain,
        MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>& lattice,
        MultiTensorField3D<double,3>& field,
        plint boundaryWidth);

template
void applyProcessingFunctional<double, descriptors::D3Q19Descriptor> (
        BoundedReductiveLatticeBoxProcessingFunctional3D<double,descriptors::D3Q19Descriptor>& functional,
        Box3D domain, std::vector<MultiBlockLattice3D<double,descriptors::D3Q19Descriptor>*> lattices,
        plint boundaryWidth);

template
void applyProcessingFunctional<double> (
        BoundedReductiveScalarFieldBoxProcessingFunctional3D<double>& functional,
        Box3D domain, std::vector<MultiScalarField3D<double>*> fields, plint boundaryWidth );

template
void applyProcessingFunctional<double,3> (
        BoundedReductiveTensorFieldBoxProcessingFunctional3D<double,3>& functional,
        Box3D domain, std::vector<MultiTensorField3D<double,3>*> fields, plint boundaryWidth );


}  // namespace plb
