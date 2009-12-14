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
#ifndef MULTI_DATA_COUPLING_WRAPPER_3D_H
#define MULTI_DATA_COUPLING_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "core/blockSurface2D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"

namespace plb {

// Forward declarations
template<typename T> class MultiBlock3D;
template<typename T, template<typename U> class Descriptor> class MultiBlockLattice3D;
template<typename T> class MultiScalarField3D;
template<typename T, int nDim> class MultiTensorField3D;

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks, plint level=0);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2);
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2, plint level=0);

/// Easy instantiation of boxed data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2);
template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2, plint level=0);

/// Easy instantiation of boxed data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                               Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2);
template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   Box3D domain,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2, plint level=0);

/// Easy instantiation of boxed data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2);
template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2, plint level=0);

/// Easy instantiation of boxed data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field, plint level=0);

/// Easy instantiation of boxed data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint level=0);

template<typename T>
void applyProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiScalarField3D<T>*> fields);
template<typename T>
void integrateProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiScalarField3D<T>*> fields, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                               Box3D domain, std::vector<MultiTensorField3D<T,nDim>*> fields);
template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain, std::vector<MultiTensorField3D<T,nDim>*> fields, plint level=0);



/// Easy instantiation of dotted data processor (general case)
template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList, std::vector<MultiBlock3D<T>*> multiBlocks);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList, std::vector<MultiBlock3D<T>*> multiBlocks, plint level=0);

/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2);
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2, plint level=0);

/// Easy instantiation of dotted data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2);
template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional,
                                   DotList3D const& dotList,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2, plint level=0);

/// Easy instantiation of dotted data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                               DotList3D const& dotList,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2);
template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   DotList3D const& dotList,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2, plint level=0);

/// Easy instantiation of dotted data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2);
template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                                   DotList3D const& dotList,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2, plint level=0);

/// Easy instantiation of dotted data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field, plint level=0);

/// Easy instantiation of dotted data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                               DotList3D const& dotList, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                                   DotList3D const& dotList, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint level=0);

template<typename T>
void applyProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList, std::vector<MultiScalarField3D<T>*> fields);
template<typename T>
void integrateProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList, std::vector<MultiScalarField3D<T>*> fields, plint level=0);

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                               DotList3D const& dotList, std::vector<MultiTensorField3D<T,nDim>*> fields);
template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                                   DotList3D const& dotList, std::vector<MultiTensorField3D<T,nDim>*> fields, plint level=0);


/// Generic implementation of "apply" and "integrate" for Bounded BoxProcessingFunctionals into MultiBlocks
template<typename T>
class BoundedCoupledMultiBlocksProcessingFunctionalOperation3D {
public:
    BoundedCoupledMultiBlocksProcessingFunctionalOperation3D(Box3D const& domain, plint boundaryWidth_);
    void apply(BoundedBoxProcessingFunctional3D<T>* functional, std::vector<MultiBlock3D<T>*> multiBlocks);
    void integrate(BoundedBoxProcessingFunctional3D<T>* functional, std::vector<MultiBlock3D<T>*> multiBlocks, plint level);
private:
    BlockSurface3D surf;
};

/// Easy instantiation of bounded boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks,
                               plint boundaryWidth );
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth = Descriptor1<T>::boundaryWidth );
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2,
                                   plint boundaryWidth = Descriptor1<T>::boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiScalarField-MultiScalarField coupling
template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional, Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2,
                               plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional, Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiTensorField-MultiTensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2,
                               plint boundaryWidth);

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, Box3D domain,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for MultiScalarField-MultiTensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional, Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2,
                               plint boundaryWidth );
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional, Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2,
                                   plint boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for Lattice-MultiScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

/// Easy instantiation of bounded boxed data processor for Lattice-MultiTensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);
template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain, std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices,
                                   plint boundaryWidth = Descriptor<T>::boundaryWidth, plint level=0);

template<typename T>
void applyProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiScalarField3D<T>*> fields,
                               plint boundaryWidth);
template<typename T>
void integrateProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiScalarField3D<T>*> fields,
                                   plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                               Box3D domain, std::vector<MultiTensorField3D<T,nDim>*> fields,
                               plint boundaryWidth);
template<typename T, int nDim>
void integrateProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain, std::vector<MultiTensorField3D<T,nDim>*> fields,
                                   plint boundaryWidth);

}  // namespace plb

#endif  // MULTI_DATA_COUPLING_WRAPPER_3D_H
