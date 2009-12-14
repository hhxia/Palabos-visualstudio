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


#ifndef REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_2D_HH
#define REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_2D_HH

#include "multiBlock/reductiveMultiDataCouplingWrapper2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** ReductiveBoxProcessing2D, general case *************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator2D<T> generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveBoxProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_SS<T>& functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveLatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain,
                               std::vector<MultiScalarField2D<T>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}


/* *************** ReductiveTensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain,
                               std::vector<MultiTensorField2D<T,nDim>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks);
}




/* *************** ReductiveDotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList, std::vector<MultiBlock2D<T>*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator2D<T> generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** ReductiveDotProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_SS<T>& functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               DotList2D const& dotList,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_ST<T,nDim>& functional,
                               DotList2D const& dotList,
                               MultiScalarField2D<T>& field1,
                               MultiTensorField2D<T,nDim>& field2 )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LS<T,Descriptor>& functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               DotList2D const& dotList,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveLatticeDotProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional2D<T,Descriptor>& functional,
                               DotList2D const& dotList,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveScalarFieldDotProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList,
                               std::vector<MultiScalarField2D<T>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}


/* *************** ReductiveTensorFieldDotProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional2D<T,nDim>& functional,
                               DotList2D const& dotList,
                               std::vector<MultiTensorField2D<T,nDim>*> fields )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, multiBlocks);
}



/* *************** Class BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation2D ******** */

template<typename T>
BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation2D<T>::
    BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation2D (
        Box2D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation2D<T>::apply (
        BoundedReductiveBoxProcessingFunctional2D<T>& functional, std::vector<MultiBlock2D<T>*> multiBlocks )
{
    std::vector<ReductiveBoxProcessorGenerator2D<T>*> generators;
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getBulkProcessor(), surf.bulk()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getEdgeProcessor(0,-1), surf.edge0N()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getEdgeProcessor(0,+1), surf.edge0P()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getEdgeProcessor(1,-1), surf.edge1N()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getEdgeProcessor(1,+1), surf.edge1P()) );

    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getCornerProcessor(-1,-1), surf.cornerNN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getCornerProcessor(+1,-1), surf.cornerPN()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getCornerProcessor(-1,+1), surf.cornerNP()) );
    generators.push_back(new ReductiveBoxProcessorGenerator2D<T> (
                functional.getCornerProcessor(+1,+1), surf.cornerPP()) );

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


/* *************** BoundedReductiveBoxProcessing2D, general case *************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<MultiBlock2D<T>*> multiBlocks,
                               plint boundaryWidth )
{
    BoundedReductiveCoupledMultiBlocksProcessingFunctionalOperation2D<T>(domain, boundaryWidth).
        apply(functional, multiBlocks);
}


/* *************** BoundedReductiveBoxProcessing2D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor1>& lattice1,
                               MultiBlockLattice2D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_SS<T>& functional,
                               Box2D domain,
                               MultiScalarField2D<T>& field1,
                               MultiScalarField2D<T>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               Box2D domain,
                               MultiTensorField2D<T,nDim1>& field1,
                               MultiTensorField2D<T,nDim2>& field2,
                               plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&field1);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional,
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



/* *************** BoundedReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiScalarField2D<T>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               Box2D domain,
                               MultiBlockLattice2D<T,Descriptor>& lattice,
                               MultiTensorField2D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<MultiBlock2D<T>*> multiBlocks(2);
    multiBlocks[0] = dynamic_cast<MultiBlock2D<T>*>(&lattice);
    multiBlocks[1] = dynamic_cast<MultiBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveLatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain,
                               std::vector<MultiBlockLattice2D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        multiBlocks[iLattice] = dynamic_cast<MultiBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain,
                               std::vector<MultiScalarField2D<T>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain,
                               std::vector<MultiTensorField2D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<MultiBlock2D<T>*> multiBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        multiBlocks[iField] = dynamic_cast<MultiBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, multiBlocks, boundaryWidth);
}


}  // namespace plb

#endif  // REDUCTIVE_MULTI_DATA_COUPLING_WRAPPER_2D_HH
