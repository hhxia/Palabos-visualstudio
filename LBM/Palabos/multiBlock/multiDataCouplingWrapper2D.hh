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


#ifndef MULTI_DATA_COUPLING_WRAPPER_2D_HH
#define MULTI_DATA_COUPLING_WRAPPER_2D_HH

#include "multiBlock/multiDataCouplingWrapper2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing2D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks)
{
    executeDataProcessor (
            BoxProcessorGenerator2D<T>(functional, domain),
            multiBlocks );
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks, plint level)
{
    addInternalProcessor (
            BoxProcessorGenerator2D<T>(functional, domain),
            multiBlocks, level );
}

/* *************** BoxProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D_SS<T>* functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D_SS<T>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                               Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                                   Box2D domain,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_ST<T,nDim>* functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_ST<T,nDim>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional2D_LS<T,Descriptor>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LS<T,Descriptor>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                                   Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** LatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                               Box2D domain,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                                   Box2D domain,
                                   std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** ScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain,
                                   std::vector<MultiScalarField2D<T>*> fields)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


template<typename T>
void integrateProcessingFunctional(ScalarFieldBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain,
                                   std::vector<MultiScalarField2D<T>*> fields, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** TensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                                   Box2D domain,
                                   std::vector<MultiTensorField2D<T,nDim>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                                   Box2D domain,
                                   std::vector<MultiTensorField2D<T,nDim>*> fields, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}



/* *************** DotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D<T>* functional,
                               DotList2D const& dotList, std::vector<MultiBlock2D<T>*> multiBlocks)
{
    executeDataProcessor (
            DotProcessorGenerator2D<T>(functional, dotList),
            multiBlocks );
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D<T>* functional,
                                   DotList2D const& dotList, std::vector<MultiBlock2D<T>*> multiBlocks, plint level)
{
    addInternalProcessor (
            DotProcessorGenerator2D<T>(functional, dotList),
            multiBlocks, level );
}

/* *************** DotProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(DotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(DotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D_SS<T>* functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D_SS<T>* functional,
                                   DotList2D const& dotList,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(DotProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                               DotList2D const& dotList,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(DotProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                                   DotList2D const& dotList,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_ST<T,nDim>* functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_ST<T,nDim>* functional,
                                   DotList2D const& dotList,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional2D_LS<T,Descriptor>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional2D_LS<T,Descriptor>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_LT<T,Descriptor,nDim>* functional,
                                   DotList2D const& dotList,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}

/* *************** LatticeDotProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeDotProcessingFunctional2D<T,Descriptor>* functional,
                               DotList2D const& dotList,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeDotProcessingFunctional2D<T,Descriptor>* functional,
                                   DotList2D const& dotList,
                                   std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** ScalarFieldDotProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldDotProcessingFunctional2D<T>* functional,
                               DotList2D const& dotList,
                               std::vector<MultiScalarField2D<T>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


template<typename T>
void integrateProcessingFunctional(ScalarFieldDotProcessingFunctional2D<T>* functional,
                                   DotList2D const& dotList,
                                   std::vector<MultiScalarField2D<T>*> fields, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** TensorFieldDotProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldDotProcessingFunctional2D<T,nDim>* functional,
                               DotList2D const& dotList,
                               std::vector<MultiTensorField2D<T,nDim>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldDotProcessingFunctional2D<T,nDim>* functional,
                                   DotList2D const& dotList,
                                   std::vector<MultiTensorField2D<T,nDim>*> fields, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}



/* *************** Class BoundedCoupledMultiBlocksProcessingFunctionalOperation2D ******** */

template<typename T>
BoundedCoupledMultiBlocksProcessingFunctionalOperation2D<T>::BoundedCoupledMultiBlocksProcessingFunctionalOperation2D (
        Box2D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedCoupledMultiBlocksProcessingFunctionalOperation2D<T>::apply (
        BoundedBoxProcessingFunctional2D<T>* functional, std::vector<MultiBlock2D<T>*> multiBlocks )
{
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getBulkProcessor(), surf.bulk()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,-1), surf.edge0N()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,+1), surf.edge0P()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,-1), surf.edge1N()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,+1), surf.edge1P()), multiBlocks);

    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,-1), surf.cornerNN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,-1), surf.cornerPN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,+1), surf.cornerNP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,+1), surf.cornerPP()), multiBlocks);
    delete functional;
}

template<typename T>
void BoundedCoupledMultiBlocksProcessingFunctionalOperation2D<T>::integrate (
        BoundedBoxProcessingFunctional2D<T>* functional, std::vector<MultiBlock2D<T>*> multiBlocks, plint level )
{
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getBulkProcessor(), surf.bulk()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,-1), surf.edge0N()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,+1), surf.edge0P()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,-1), surf.edge1N()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,+1), surf.edge1P()), multiBlocks, level);

    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,-1), surf.cornerNN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,-1), surf.cornerPN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,+1), surf.cornerNP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,+1), surf.cornerPP()), multiBlocks, level);
    delete functional;
}

/* *************** BoundedBoxProcessing2D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D<T>* functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks,
                               plint boundaryWidth )
{
    BoundedCoupledMultiBlocksProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks,
                                   plint boundaryWidth, plint level)
{
    BoundedCoupledMultiBlocksProcessingFunctionalOperation2D<T>(domain, boundaryWidth).integrate(functional, multiBlocks, level);
}

/* *************** BoundedBoxProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice2D<T,Descriptor2>& lattice2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_SS<T>* functional, Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2,
                               plint boundaryWidth)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_SS<T>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiScalarField2D<T>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional, Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2,
                               plint boundaryWidth)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>* functional,
                                   Box2D domain,
                                   MultiTensorField2D<T,nDim1>& field1,
                                   MultiTensorField2D<T,nDim2>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_ST<T,nDim>* functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_ST<T,nDim>* functional,
                                   Box2D domain,
                                   MultiScalarField2D<T>& field1,
                                   MultiTensorField2D<T,nDim>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LS<T,Descriptor>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LS<T,Descriptor>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiScalarField2D<T>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional, Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>* functional, Box2D domain,
                                   MultiBlockLattice2D<T,Descriptor>& lattice,
                                   MultiTensorField2D<T,nDim>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedLatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                               Box2D domain,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>* functional,
                                   Box2D domain,
                                   std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional2D<T>* functional,
                               Box2D domain,
                               std::vector<MultiScalarField2D<T>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional2D<T>* functional,
                                   Box2D domain,
                                   std::vector<MultiScalarField2D<T>*> fields, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedTensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                               Box2D domain,
                               std::vector<MultiTensorField2D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>* functional,
                                   Box2D domain,
                                   std::vector<MultiTensorField2D<T,nDim>*> fields, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_DATA_COUPLING_WRAPPER_2D_HH
