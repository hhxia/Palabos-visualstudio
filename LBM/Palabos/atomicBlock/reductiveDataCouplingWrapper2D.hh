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


#ifndef REDUCTIVE_DATA_COUPLING_WRAPPER_2D_HH
#define REDUCTIVE_DATA_COUPLING_WRAPPER_2D_HH

#include "atomicBlock/reductiveDataCouplingWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/atomicBlockOperations2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"
#include "multiBlock/combinedStatistics.h"

namespace plb {

/* *************** ReductiveBoxProcessing2D, general case *************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D<T>& functional,
                               Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator2D<T> generator(functional.clone(), domain);
    executeDataProcessor(generator, atomicBlocks );
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing2D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_SS<T>& functional,
                               Box2D domain,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               Box2D domain,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional,
                               Box2D domain,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks);
}

/* *************** ReductiveLatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain,
                               std::vector<BlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}


/* *************** ReductiveScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain,
                               std::vector<ScalarField2D<T>*> fields )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}


/* *************** ReductiveTensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain,
                               std::vector<TensorField2D<T,nDim>*> fields )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks);
}



/* *************** ReductiveDotProcessing, general case ***************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator2D<T> generator(functional.clone(), dotList);
    executeDataProcessor(generator, atomicBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_SS<T>& functional,
                               DotList2D const& dotList,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               DotList2D const& dotList,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_ST<T,nDim>& functional,
                               DotList2D const& dotList,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2 )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LS<T,Descriptor>& functional,
                               DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               DotList2D const& dotList,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}

/* *************** ReductiveLatticeDotProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveLatticeDotProcessingFunctional2D<T,Descriptor>& functional,
                               DotList2D const& dotList,
                               std::vector<BlockLattice2D<T,Descriptor>*> lattices )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


/* *************** ReductiveScalarFieldDotProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(ReductiveScalarFieldDotProcessingFunctional2D<T>& functional,
                               DotList2D const& dotList,
                               std::vector<ScalarField2D<T>*> fields )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}


/* *************** ReductiveTensorFieldDotProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveTensorFieldDotProcessingFunctional2D<T,nDim>& functional,
                               DotList2D const& dotList,
                               std::vector<TensorField2D<T,nDim>*> fields )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, dotList, atomicBlocks);
}



/* *************** Class BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D ******** */

template<typename T>
BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D<T>::
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D (
        Box2D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D<T>::apply (
        BoundedReductiveBoxProcessingFunctional2D<T>& functional, std::vector<AtomicBlock2D<T>*> atomicBlocks )
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
        executeDataProcessor(*generators[iGenerator], atomicBlocks);
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
                               Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks,
                               plint boundaryWidth )
{
    BoundedReductiveCoupledBlocksProcessingFunctionalOperation2D<T>(domain, boundaryWidth).
        apply(functional, atomicBlocks);
}


/* *************** BoundedReductiveBoxProcessing2D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor1>& lattice1,
                               BlockLattice2D<T,Descriptor2>& lattice2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&lattice2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_SS<T>& functional,
                               Box2D domain,
                               ScalarField2D<T>& field1,
                               ScalarField2D<T>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>& functional,
                               Box2D domain,
                               TensorField2D<T,nDim1>& field1,
                               TensorField2D<T,nDim2>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>& functional,
                               Box2D domain,
                               ScalarField2D<T>& field1,
                               TensorField2D<T,nDim>& field2,
                               plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&field1);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field2);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               ScalarField2D<T>& field,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}



/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>& functional,
                               Box2D domain,
                               BlockLattice2D<T,Descriptor>& lattice,
                               TensorField2D<T,nDim>& field,
                               plint boundaryWidth)
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(2);
    atomicBlocks[0] = dynamic_cast<AtomicBlock2D<T>*>(&lattice);
    atomicBlocks[1] = dynamic_cast<AtomicBlock2D<T>*>(&field);
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveLatticeBoxProcessing2D ****************************************** */

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>& functional,
                               Box2D domain,
                               std::vector<BlockLattice2D<T,Descriptor>*> lattices, plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(lattices.size());
    for (pluint iLattice=0; iLattice<lattices.size(); ++iLattice) {
        atomicBlocks[iLattice] = dynamic_cast<AtomicBlock2D<T>*>(lattices[iLattice]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing2D ****************************************** */

template<typename T>
void applyProcessingFunctional(BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>& functional,
                               Box2D domain,
                               std::vector<ScalarField2D<T>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing2D ****************************************** */

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>& functional,
                               Box2D domain,
                               std::vector<TensorField2D<T,nDim>*> fields, plint boundaryWidth )
{
    std::vector<AtomicBlock2D<T>*> atomicBlocks(fields.size());
    for (pluint iField=0; iField<fields.size(); ++iField) {
        atomicBlocks[iField] = dynamic_cast<AtomicBlock2D<T>*>(fields[iField]);
    }
    applyProcessingFunctional(functional, domain, atomicBlocks, boundaryWidth);
}


}  // namespace plb

#endif  // REDUCTIVE_DATA_COUPLING_WRAPPER_2D_HH
