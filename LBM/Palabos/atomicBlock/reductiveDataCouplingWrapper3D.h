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
#ifndef REDUCTIVE_DATA_COUPLING_WRAPPER_3D_H
#define REDUCTIVE_DATA_COUPLING_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "core/blockSurface3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"

namespace plb {

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2);

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_SS<T>& functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2);

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional, Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2);

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2);

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field);

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices);
template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<ScalarField3D<T>*> fields);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain, std::vector<TensorField3D<T,nDim>*> fields);



/// Easy instantiation of dotted data processor (general case)
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks);
/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional, DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2);

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_SS<T>& functional, DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2);

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_TT<T,nDim1,nDim2>& functional, DotList3D const& dotList,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2);

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_ST<T,nDim>& functional, DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2);

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LS<T,Descriptor>& functional, DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field);

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional3D<T,Descriptor>& functional,
                               DotList3D const& dotList, std::vector<BlockLattice3D<T,Descriptor>*> lattices);
template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList, std::vector<ScalarField3D<T>*> fields);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional3D<T,nDim>& functional,
                               DotList3D const& dotList, std::vector<TensorField3D<T,nDim>*> fields);



/// Generic implementation of "apply" and "integrate" for Bounded Reductive BoxProcessingFunctionals
template<typename T>
class BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D {
public:
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D (
            Box3D const& domain, plint boundaryWidth_ );
    void apply(BoundedReductiveBoxProcessingFunctional3D<T>& functional,
               std::vector<AtomicBlock3D<T>*> atomicBlocks );
private:
    BlockSurface3D surf;
};

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_SS<T>& functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional, Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field,
                               plint boundaryWidth);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<ScalarField3D<T>*> field, plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain, std::vector<TensorField3D<T,nDim>*> field, plint boundaryWidth);



}  // namespace plb

#endif  // REDUCTIVE_DATA_COUPLING_WRAPPER_3D_H
