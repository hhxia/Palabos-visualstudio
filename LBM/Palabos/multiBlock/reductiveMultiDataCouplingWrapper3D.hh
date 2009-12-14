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


#ifndef REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_3D_HH
#define REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_3D_HH

#include "multiBlock/reductiveMultiDataCouplingWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** ReductiveBoxProcessing3D, general case *************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D<T>& functional,
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D<T> generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_SS<T>& functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain,
                               std::vector<MultiScalarField3D<T>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain,
                               std::vector<MultiTensorField3D<T,nDim>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}




/* *************** ReductiveDotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList, std::vector<MultiBlock3D<T>*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D<T> generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_SS<T>& functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               DotList3D const& dotList,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_ST<T,nDim>& functional,
                               DotList3D const& dotList,
                               MultiScalarField3D<T>& field1,
                               MultiTensorField3D<T,nDim>& field2 )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LS<T,Descriptor>& functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               DotList3D const& dotList,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveLatticeDotProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional3D<T,Descriptor>& functional,
                               DotList3D const& dotList,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveScalarFieldDotProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional3D<T>& functional,
                               DotList3D const& dotList,
                               std::vector<MultiScalarField3D<T>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveTensorFieldDotProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional3D<T,nDim>& functional,
                               DotList3D const& dotList,
                               std::vector<MultiTensorField3D<T,nDim>*> fields )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}



/* *************** Class BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D<T>::
    BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D<T>::apply (
        BoundedReductiveBoxProcessingFunctional3D<T>& functional, std::vector<MultiBlock3D<T>*> multiBlocks )
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
        executeDataProcessor(*generators[iGenerator], multiBlocks);
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
                               Box3D domain, std::vector<MultiBlock3D<T>*> multiBlocks,
                               plint boundaryWidth )
{
    BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation3D<T>(domain, boundaryWidth).
        apply(functional, multiBlocks);
}


/* *************** BoundedReductiveBoxProcessing3D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor1>& lattice1,
                               MultiBlockLattice3D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_SS<T>& functional,
                               Box3D domain,
                               MultiScalarField3D<T>& field1,
                               MultiScalarField3D<T>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>& functional,
                               Box3D domain,
                               MultiTensorField3D<T,nDim1>& field1,
                               MultiTensorField3D<T,nDim2>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>& functional,
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



/* *************** BoundedReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiScalarField3D<T>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>& functional,
                               Box3D domain,
                               MultiBlockLattice3D<T,Descriptor>& lattice,
                               MultiTensorField3D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock3D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock3D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock3D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveLatticeBoxProcessing3D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>& functional,
                               Box3D domain,
                               std::vector<MultiBlockLattice3D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock3D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing3D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>& functional,
                               Box3D domain,
                               std::vector<MultiScalarField3D<T>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing3D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>& functional,
                               Box3D domain,
                               std::vector<MultiTensorField3D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock3D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock3D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


}  // namespace plb

#endif  // REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_3D_HH
