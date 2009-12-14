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


#ifndef MULTI_DATA_COUPLING_WRAPPER_3D_HH
#define MULTI_DATA_COUPLING_WRAPPER_3D_HH

#include "multiBlock/multiDataCouplingWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks)
{
    executeDataProcessor (
            BoxProcessorGenerator3D<T>(functional, domain),
            multiBlocks );
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks, plint level)
{
    addInternalProcessor (
            BoxProcessorGenerator3D<T>(functional, domain),
            multiBlocks, level );
}

/* *************** BoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                               Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   Box3D domain,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** BoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                                   Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** LatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain,
                                   std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** ScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain,
                                   std::vector<MultiScalarField3D<T>*> fields)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


template<typename T>
void integrateProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain,
                                   std::vector<MultiScalarField3D<T>*> fields, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}


/* *************** TensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain,
                                   std::vector<MultiTensorField3D<T,nDim>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain,
                                   std::vector<MultiTensorField3D<T,nDim>*> fields, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, level);
}



/* *************** DotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList, std::vector<MultiBlock3D<T>*> multiBlocks)
{
    executeDataProcessor (
            DotProcessorGenerator3D<T>(functional, dotList),
            multiBlocks );
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList, std::vector<MultiBlock3D<T>*> multiBlocks, plint level)
{
    addInternalProcessor (
            DotProcessorGenerator3D<T>(functional, dotList),
            multiBlocks, level );
}

/* *************** DotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional,
                                   DotList3D const& dotList,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                               DotList3D const& dotList,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   DotList3D const& dotList,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                                   DotList3D const& dotList,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** DotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional,
                                   DotList3D const& dotList,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}

/* *************** LatticeDotProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                               DotList3D const& dotList,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                                   DotList3D const& dotList,
                                   std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** ScalarFieldDotProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList,
                               std::vector<MultiScalarField3D<T>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


template<typename T>
void integrateProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList,
                                   std::vector<MultiScalarField3D<T>*> fields, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}


/* *************** TensorFieldDotProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                               DotList3D const& dotList,
                               std::vector<MultiTensorField3D<T,nDim>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                                   DotList3D const& dotList,
                                   std::vector<MultiTensorField3D<T,nDim>*> fields, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, multiBlocks, level);
}



/* *************** Class BoundedCoupledMultiBlocksProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedCoupledMultiBlocksProcessingFunctionalOperation3D<T>::BoundedCoupledMultiBlocksProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedCoupledMultiBlocksProcessingFunctionalOperation3D<T>::apply (
        BoundedBoxProcessingFunctional3D<T>* functional, std::vector<MultiBlock3D<T>*> multiBlocks )
{
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()), multiBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()), multiBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()), multiBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1, -1), surf.cornerNNN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1,  1), surf.cornerNNP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1, -1), surf.cornerNPN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1,  1), surf.cornerNPP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1, -1), surf.cornerPNN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1,  1), surf.cornerPNP()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1, -1), surf.cornerPPN()), multiBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1,  1), surf.cornerPPP()), multiBlocks);

    delete functional;

}

template<typename T>
void BoundedCoupledMultiBlocksProcessingFunctionalOperation3D<T>::integrate (
        BoundedBoxProcessingFunctional3D<T>* functional, std::vector<MultiBlock3D<T>*> multiBlocks, plint level )
{
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()), multiBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()), multiBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()), multiBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1, -1), surf.cornerNNN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1,  1), surf.cornerNNP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1, -1), surf.cornerNPN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1,  1), surf.cornerNPP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1, -1), surf.cornerPNN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1,  1), surf.cornerPNP()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1, -1), surf.cornerPPN()), multiBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1,  1), surf.cornerPPP()), multiBlocks, level);

    delete functional;
}

/* *************** BoundedBoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks,
                               plint boundaryWidth )
{
    BoundedCoupledMultiBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, multiBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks,
                                   plint boundaryWidth, plint level)
{
    BoundedCoupledMultiBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).integrate(functional, multiBlocks, level);
}

/* *************** BoundedBoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor1>& lattice1,
                                   MultiBlockLattice3D<T,Descriptor2>& lattice2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional, Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2,
                               plint boundaryWidth)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiScalarField3D<T>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2,
                               plint boundaryWidth)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   Box3D domain,
                                   MultiTensorField3D<T,nDim1>& field1,
                                   MultiTensorField3D<T,nDim2>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional,
                                   Box3D domain,
                                   MultiScalarField3D<T>& field1,
                                   MultiTensorField3D<T,nDim>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiScalarField3D<T>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                                   MultiBlockLattice3D<T,Descriptor>& lattice,
                                   MultiTensorField3D<T,nDim>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain,
                                   std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                               Box3D domain,
                               std::vector<MultiScalarField3D<T>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain,
                                   std::vector<MultiScalarField3D<T>*> fields, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}


/* *************** BoundedTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                               Box3D domain,
                               std::vector<MultiTensorField3D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain,
                                   std::vector<MultiTensorField3D<T,nDim>*> fields, plint boundaryWidth, plint level )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, multiBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // MULTI_DATA_COUPLING_WRAPPER_3D_HH
