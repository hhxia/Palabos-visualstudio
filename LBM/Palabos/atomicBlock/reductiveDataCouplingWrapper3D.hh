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


#ifndef REDUCTIVE_DATA_COUPLING_WRAPPER_3D_HH
#define REDUCTIVE_DATA_COUPLING_WRAPPER_3D_HH

#include "atomicBlock/reductiveDataCouplingWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/atomicBlockOperations3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"
#include "multiBlock/combinedStatistics.h"

namespace plb {

/* *************** ReductiveBoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D<T> generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks );
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing3D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_SS<T>& functional,
                               Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional,
                               Box3D domain,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

/* *************** ReductiveLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}


/* *************** ReductiveScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain,
                               std::vector<ScalarField3D<T>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}


/* *************** ReductiveTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveDotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D<T> generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_SS<T>& functional,
                               DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               DotList3D const& dotList,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_ST<T,nDim>& functional,
                               DotList3D const& dotList,
                               ScalarField3D<T>& field1,
                               TensorField3D<T,nDim>& field2 )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LS<T,Descriptor>& functional,
                               DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               DotList3D const& dotList,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

/* *************** ReductiveLatticeDotProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional3D<T,Descriptor>& functional,
                               DotList3D const& dotList,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


/* *************** ReductiveScalarFieldDotProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList,
                               std::vector<ScalarField3D<T>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


/* *************** ReductiveTensorFieldDotProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional3D<T,nDim>& functional,
                               DotList3D const& dotList,
                               std::vector<TensorField3D<T,nDim>*> fields )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** Class BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D<T>::
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D<T>::apply (
        BoundedReductiveBoxProcessingFunctional3D<T>& functional, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ReductiveBoxProcessorGenerator3D<T>*> generators;

    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getBulkProcessor(), surf.bulk()) );

    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(0,-1), surf.surface0N()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(0,+1), surf.surface0P()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(1,-1), surf.surface1N()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(1,+1), surf.surface1P()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(2,-1), surf.surface2N()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getPlaneProcessor(2,+1), surf.surface2P()) );

    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(0, -1, -1), surf.edge0NN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(0, -1,  1), surf.edge0NP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(0,  1, -1), surf.edge0PN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(0,  1,  1), surf.edge0PP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(1, -1, -1), surf.edge1NN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(1, -1,  1), surf.edge1NP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(1,  1, -1), surf.edge1PN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(1,  1,  1), surf.edge1PP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(2, -1, -1), surf.edge2NN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(2, -1,  1), surf.edge2NP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(2,  1, -1), surf.edge2PN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getEdgeProcessor(2,  1,  1), surf.edge2PP()) );

    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor(-1, -1, -1), surf.cornerNNN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor(-1, -1,  1), surf.cornerNNP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor(-1,  1, -1), surf.cornerNPN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor(-1,  1,  1), surf.cornerNPP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor( 1, -1, -1), surf.cornerPNN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor( 1, -1,  1), surf.cornerPNP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor( 1,  1, -1), surf.cornerPPN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator3D<T> (
                functional.getCornerProcessor( 1,  1,  1), surf.cornerPPP()) );

    std::vector<BlockStatistics<T> const*> individualStatistics(generators.size());

    for (pluint iGenerator=0; iGenerator<generators.size(); ++iGenerator) {
            executeDataProcessor(*generators[iGenerator], atomicBlocks);
        individualStatistics[iGenerator] = &(generators[iGenerator]->getStatistics());
    }

    SerialCombinedStatistics<T>().combine(individualStatistics, functional.getStatistics());

    for (pluint iGenerator=0; iGenerator<generators.size(); ++iGenerator) {
        delete generators[iGenerator];
    }
}


/* *************** BoundedReductiveBoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks,
                               plint boundaryWidth )
{
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).
        apply(functional, atomicBlocks);
}


/* *************** BoundedReductiveBoxProcessing3D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor1>& lattice1,
                               BlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_SS<T>& functional,
                               Box3D domain,
                               ScalarField3D<T>& field1,
                               ScalarField3D<T>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               Box3D domain,
                               TensorField3D<T,nDim1>& field1,
                               TensorField3D<T,nDim2>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional,
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



/* *************** BoundedReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               ScalarField3D<T>& field,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               Box3D domain,
                               BlockLattice3D<T,Descriptor>& lattice,
                               TensorField3D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock3D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain,
                               std::vector<BlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain,
                               std::vector<ScalarField3D<T>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain,
                               std::vector<TensorField3D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock3D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


}  // namespace plb

#endif  // REDUCTIVE_DATA_COUPLING_WRAPPER_3D_HH
