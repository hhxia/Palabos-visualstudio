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


#ifndef DATA_COUPLING_WRAPPER_3D_HH
#define DATA_COUPLING_WRAPPER_3D_HH

#include "atomicBlock/dataCouplingWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/atomicBlockOperations3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    executeDataProcessor( BoxProcessorGenerator3D<T>(functional, domain),
                          atomicBlocks );
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks, plint level)
{
    addInternalProcessor( BoxProcessorGenerator3D<T>(functional, domain),
                          atomicBlocks, level );
}

/* *************** BoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor1>& lattice1,
                                   BlockLattice3D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** BoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D_SS<T>* functional,
                                   Box3D domain,
                                   ScalarField3D<T>& field1,
                                   ScalarField3D<T>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** BoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   Box3D domain,
                                   TensorField3D<T,nDim1>& field1,
                                   TensorField3D<T,nDim2>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** BoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                               Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_ST<T,nDim>* functional,
                                   Box3D domain,
                                   ScalarField3D<T>& field1,
                                   TensorField3D<T,nDim>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** BoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   ScalarField3D<T>& field, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** BoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   TensorField3D<T,nDim>& field, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}

/* *************** LatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain,
                                   std::vector<BlockLattice3D<T,Descriptor>*> lattices, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** ScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                               Box3D domain,
                               std::vector<ScalarField3D<T>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T>
void integrateProcessingFunctional(ScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain,
                                   std::vector<ScalarField3D<T>*> fields, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** TensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                               Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain,
                                   std::vector<TensorField3D<T,nDim>*> fields, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, level);
}


/* *************** DotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    executeDataProcessor( DotProcessorGenerator3D<T>(functional, dotList),
                          atomicBlocks );
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks, plint level)
{
    addInternalProcessor( DotProcessorGenerator3D<T>(functional, dotList),
                          atomicBlocks, level );
}

/* *************** DotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, DotList3D const& dotList,
                                   BlockLattice3D<T,Descriptor1>& lattice1,
                                   BlockLattice3D<T,Descriptor2>& lattice2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** DotProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional, DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D_SS<T>* functional,
                                   DotList3D const& dotList,
                                   ScalarField3D<T>& field1,
                                   ScalarField3D<T>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** DotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, DotList3D const& dotList,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(DotProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   DotList3D const& dotList,
                                   TensorField3D<T,nDim1>& field1,
                                   TensorField3D<T,nDim2>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** DotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                               DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_ST<T,nDim>* functional,
                                   DotList3D const& dotList,
                                   ScalarField3D<T>& field1,
                                   TensorField3D<T,nDim>& field2, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** DotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional, DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional3D_LS<T,Descriptor>* functional, DotList3D const& dotList,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   ScalarField3D<T>& field, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** DotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, DotList3D const& dotList,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   TensorField3D<T,nDim>& field, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}

/* *************** LatticeDotProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                               DotList3D const& dotList,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(LatticeDotProcessingFunctional3D<T,Descriptor>* functional,
                                   DotList3D const& dotList,
                                   std::vector<BlockLattice3D<T,Descriptor>*> lattices, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** ScalarFieldDotProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                               DotList3D const& dotList,
                               std::vector<ScalarField3D<T>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


template<typename T>
void integrateProcessingFunctional(ScalarFieldDotProcessingFunctional3D<T>* functional,
                                   DotList3D const& dotList,
                                   std::vector<ScalarField3D<T>*> fields, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** TensorFieldDotProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                               DotList3D const& dotList,
                               std::vector<TensorField3D<T,nDim>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

template<typename T, int nDim>
void integrateProcessingFunctional(TensorFieldDotProcessingFunctional3D<T,nDim>* functional,
                                   DotList3D const& dotList,
                                   std::vector<TensorField3D<T,nDim>*> fields, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, dotList, atomicBlocks, level);
}


/* *************** Class BoundedCoupledBlocksProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedCoupledBlocksProcessingFunctionalOperation3D<T>::BoundedCoupledBlocksProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedCoupledBlocksProcessingFunctionalOperation3D<T>::apply (
        BoundedBoxProcessingFunctional3D<T>* functional, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()), atomicBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()), atomicBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()), atomicBlocks);

    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1, -1), surf.cornerNNN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1,  1), surf.cornerNNP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1, -1), surf.cornerNPN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1,  1), surf.cornerNPP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1, -1), surf.cornerPNN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1,  1), surf.cornerPNP()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1, -1), surf.cornerPPN()), atomicBlocks);
    executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1,  1), surf.cornerPPP()), atomicBlocks);

    delete functional;
}

template<typename T>
void BoundedCoupledBlocksProcessingFunctionalOperation3D<T>::integrate (
        BoundedBoxProcessingFunctional3D<T>* functional, std::vector<AtomicBlock3D<T>*> atomicBlocks, plint level )
{
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()), atomicBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()), atomicBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()), atomicBlocks, level);

    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1, -1, -1), surf.cornerNNN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1, -1,  1), surf.cornerNNP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1,  1, -1), surf.cornerNPN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1,  1,  1), surf.cornerNPP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1, -1, -1), surf.cornerPNN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1, -1,  1), surf.cornerPNP()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1,  1, -1), surf.cornerPPN()), atomicBlocks, level);
    addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1,  1,  1), surf.cornerPPP()), atomicBlocks, level);

    delete functional;

}

/* *************** BoundedBoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks,
                               plint boundaryWidth )
{
    BoundedCoupledBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, atomicBlocks);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks,
                                   plint boundaryWidth, plint level)
{
    BoundedCoupledBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).integrate(functional, atomicBlocks, level);
}

/* *************** BoundedBoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor1>& lattice1,
                                   BlockLattice3D<T,Descriptor2>& lattice2,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional, Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_SS<T>* functional,
                                   Box3D domain,
                                   ScalarField3D<T>& field1,
                                   ScalarField3D<T>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional, Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, int nDim1, int nDim2>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>* functional,
                                   Box3D domain,
                                   TensorField3D<T,nDim1>& field1,
                                   TensorField3D<T,nDim2>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional,
                               Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_ST<T,nDim>* functional,
                                   Box3D domain,
                                   ScalarField3D<T>& field1,
                                   TensorField3D<T,nDim>& field2,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LS<T,Descriptor>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   ScalarField3D<T>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>* functional, Box3D domain,
                                   BlockLattice3D<T,Descriptor>& lattice,
                                   TensorField3D<T,nDim>& field,
                                   plint boundaryWidth, plint level)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                               Box3D domain,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>* functional,
                                   Box3D domain,
                                   std::vector<BlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                               Box3D domain,
                               std::vector<ScalarField3D<T>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T>
void integrateProcessingFunctional(BoundedScalarFieldBoxProcessingFunctional3D<T>* functional,
                                   Box3D domain,
                                   std::vector<ScalarField3D<T>*> fields, plint boundaryWidth, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}


/* *************** BoundedTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                               Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>* functional,
                                   Box3D domain,
                                   std::vector<TensorField3D<T,nDim>*> fields, plint boundaryWidth, plint level )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    integrateProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth, level);
}

}  // namespace plb

#endif  // DATA_COUPLING_WRAPPER_3D_HH
