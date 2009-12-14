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

/** \file
 * Utilities to help users handle data processors -- header file.
 */
#ifndef MULTI_DATA_COUPLING_WRAPPER_2D_H
#define MULTI_DATA_COUPLING_WRAPPER_2D_H

#include "core/globalDefs.h"
#include "core/blockSurface2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"

namespace plb {

// Forward declarations
template<typename T> class MultiBlock2D;
template<typename T, template<typename U> class Descriptor> class MultiBlockLattice2D;
template<typename T> class MultiScalarField2D;
template<typename T, int nDim> class MultiTensorField2D;

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks, plint level=0);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2);
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2, plint level=0);

/// Easy instantiation of boxed data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D_SS<T>* functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D_SS<T>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2, plint level=0);

/// Easy instantiation of boxed data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                               Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2);
template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                                   Box2D domain,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2, plint level=0);

/// Easy instantiation of boxed data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_ST<T,nDim>* functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2);
template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_ST<T,nDim>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2, plint level=0);

/// Easy instantiation of boxed data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional2D_LS<T,Descriptor>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LS<T,Descriptor>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field, plint level=0);

/// Easy instantiation of boxed data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                               Box2D domain, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                                   Box2D domain, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint level=0);

template<typename T>
void applyProcessingFunctional(ScalarFieldBoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiScalarField2D<T>*> fields);
template<typename T>
void integrateProcessingFunctional(ScalarFieldBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiScalarField2D<T>*> fields, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                               Box2D domain, std::vector<MultiTensorField2D<T,nDim>*> fields);
template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                                   Box2D domain, std::vector<MultiTensorField2D<T,nDim>*> fields, plint level=0);



/// Easy instantiation of dotted data processor (general case)
template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D<T>* functional,
                               DotList2D const& dotList, std::vector<MultiBlock2D<T>*> multiBlocks);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D<T>* functional,
                                   DotList2D const& dotList, std::vector<MultiBlock2D<T>*> multiBlocks, plint level=0);

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(DotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2);
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(DotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2, plint level=0);

/// Easy instantiation of dotted data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D_SS<T>* functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D_SS<T>* functional,
                                   DotList2D const& dotList,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2, plint level=0);

/// Easy instantiation of dotted data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(DotProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                               DotList2D const& dotList,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2);
template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(DotProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                                   DotList2D const& dotList,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2, plint level=0);

/// Easy instantiation of dotted data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_ST<T,nDim>* functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2);
template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_ST<T,nDim>* functional,
                                   DotList2D const& dotList,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2, plint level=0);

/// Easy instantiation of dotted data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional2D_LS<T,Descriptor>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional2D_LS<T,Descriptor>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field, plint level=0);

/// Easy instantiation of dotted data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeDotProcessingFunctional2D<T,Descriptor>* functional,
                               DotList2D const& dotList, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeDotProcessingFunctional2D<T,Descriptor>* functional,
                                   DotList2D const& dotList, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint level=0);

template<typename T>
void applyProcessingFunctional(ScalarFieldDotProcessingFunctional2D<T>* functional,
                               DotList2D const& dotList, std::vector<MultiScalarField2D<T>*> fields);
template<typename T>
void integrateProcessingFunctional(ScalarFieldDotProcessingFunctional2D<T>* functional,
                                   DotList2D const& dotList, std::vector<MultiScalarField2D<T>*> fields, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldDotProcessingFunctional2D<T,nDim>* functional,
                               DotList2D const& dotList, std::vector<MultiTensorField2D<T,nDim>*> fields);
template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldDotProcessingFunctional2D<T,nDim>* functional,
                                   DotList2D const& dotList, std::vector<MultiTensorField2D<T,nDim>*> fields, plint level=0);


/// Generic implementation of "apply" and "integrate" for Bounded BoxProcessingFunctionals into MultiBlocks
template<typename T>
class BoundedCoupledMultiBlocksProcessingFunctionalOperation2D {
public:
    BoundedCoupledMultiBlocksProcessingFunctionalOperation2D(Box2D const& domain, plint boundaryWidth_);
    void apply(BoundedBoxProcessingFunctional2D<T>* functional, std::vector<MultiBlock2D<T>*> multiBlocks);
    void integrate(BoundedBoxProcessingFunctional2D<T>* functional, std::vector<MultiBlock2D<T>*> multiBlocks, plint level);
private:
    BlockSurface2D surf;
};

/// Easy instantiation of bounded boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks,
                               plint boundaryWidth );
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2,
                               plint boundaryWidth = Descriptor1<T>::boundaryWidth );
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2,
                                   plint boundaryWidth = Descriptor1<T>::boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_SS<T>* functional, Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2,
                               plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_SS<T>* functional, Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional, Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2,
                               plint boundaryWidth);
template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional, Box2D domain,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_ST<T,nDim>* functional, Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2,
                               plint boundaryWidth );
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_ST<T,nDim>* functional, Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LS<T,Descriptor>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LS<T,Descriptor>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                               Box2D domain, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                                   Box2D domain, std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T>
void applyProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiScalarField2D<T>*> fields,
                               plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiScalarField2D<T>*> fields,
                                   plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                               Box2D domain, std::vector<MultiTensorField2D<T,nDim>*> fields,
                               plint boundaryWidth);
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                                   Box2D domain, std::vector<MultiTensorField2D<T,nDim>*> fields,
                                   plint boundaryWidth);

}  // namespace plb

#endif  // MULTI_DATA_COUPLING_WRAPPER_2D_H
