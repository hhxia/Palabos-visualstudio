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
#ifndef REDUCTIVE_DATA_COUPLING_WRAPPER_2D_H
#define REDUCTIVE_DATA_COUPLING_WRAPPER_2D_H

#include "core/globalDefs.h"
#include "core/blockSurface2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"

namespace plb {

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2);

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_SS<T>& functional, Box2D domain,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2);

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional, Box2D domain,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2);

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional, Box2D domain,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2);

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field);

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices);
template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<ScalarField2D<T>*> fields);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain, std::vector<TensorField2D<T,nDim>*> fields);



/// Easy instantiation of dotted data processor (general case)
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks);
/// Easy instantiation of dotted data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional, DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2);

/// Easy instantiation of dotted data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_SS<T>& functional, DotList2D const& dotList,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2);

/// Easy instantiation of dotted data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_TT<T,nDim1,nDim2>& functional, DotList2D const& dotList,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2);

/// Easy instantiation of dotted data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_ST<T,nDim>& functional, DotList2D const& dotList,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2);

/// Easy instantiation of dotted data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LS<T,Descriptor>& functional, DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field);

/// Easy instantiation of dotted data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional2D<T,Descriptor>& functional,
                               DotList2D const& dotList, std::vector<BlockLattice2D<T,Descriptor>*> lattices);
template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList, std::vector<ScalarField2D<T>*> fields);
template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional2D<T,nDim>& functional,
                               DotList2D const& dotList, std::vector<TensorField2D<T,nDim>*> fields);



/// Generic implementation of "apply" and "integrate" for Bounded Reductive BoxProcessingFunctionals
template<typename T>
class BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D {
public:
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D (
            Box2D const& domain, plint boundaryWidth_ );
    void apply(BoundedReductiveBoxProcessingFunctional2D<T>& functional,
               std::vector<AtomicBlock2D<T>*> atomicBlocks );
private:
    BlockSurface2D surf;
};

/// Easy instantiation of boxed data processor (general case)
template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for lattice-lattice coupling
template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for ScalarField-ScalarField coupling
template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_SS<T>& functional, Box2D domain,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for TensorField-TensorField coupling
template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional, Box2D domain,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for ScalarField-TensorField coupling
template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional, Box2D domain,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for Lattice-ScalarField coupling
template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field,
                               plint boundaryWidth);

/// Easy instantiation of boxed data processor for Lattice-TensorField coupling
template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional, Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field,
                               plint boundaryWidth);

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattice,
                               plint boundaryWidth = Descriptor<T>::boundaryWidth);

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<ScalarField2D<T>*> field, plint boundaryWidth);

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain, std::vector<TensorField2D<T,nDim>*> field, plint boundaryWidth);



}  // namespace plb

#endif  // REDUCTIVE_DATA_COUPLING_WRAPPER_2D_H
