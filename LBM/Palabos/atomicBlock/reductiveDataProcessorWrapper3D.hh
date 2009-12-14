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


#ifndef REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_HH
#define REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_HH

#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"
#include "multiBlock/combinedStatistics.h"

namespace plb {

/* *************** Class ReductiveBoxProcessingFunctional3D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT ReductiveBoxProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void ReductiveBoxProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveBoxProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}



/* *************** Class PlainReductiveBoxProcessingFunctional3D ************************* */

template<typename T>
BlockStatistics<T> const& PlainReductiveBoxProcessingFunctional3D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& PlainReductiveBoxProcessingFunctional3D<T>::getStatistics() {
    return statistics;
}


/* *************** Class ReductiveBoxProcessor3D ************************************ */

template<typename T>
ReductiveBoxProcessor3D<T>::ReductiveBoxProcessor3D (
        ReductiveBoxProcessingFunctional3D<T>* functional_,
        Box3D domain_, std::vector<AtomicBlock3D<T>*> atomicBlocks_)
    : functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
Box3D ReductiveBoxProcessor3D<T>::getDomain() const {
    return domain;
}

template<typename T>
void ReductiveBoxProcessor3D<T>::process() {
    functional -> processGenericBlocks(domain, atomicBlocks);
    functional -> getStatistics().evaluate();
}

template<typename T>
ReductiveBoxProcessor3D<T>* ReductiveBoxProcessor3D<T>::clone() const {
    return new ReductiveBoxProcessor3D<T>(*this);
}


/* *************** Class ReductiveBoxProcessorGenerator3D *************************** */

template<typename T>
ReductiveBoxProcessorGenerator3D<T>::ReductiveBoxProcessorGenerator3D (
        ReductiveBoxProcessingFunctional3D<T>* functional_,
        Box3D domain )
    : BoxedReductiveDataProcessorGenerator3D<T>(domain),
      functional(functional_)
{ }

template<typename T>
ReductiveBoxProcessorGenerator3D<T>::~ReductiveBoxProcessorGenerator3D() {
    delete functional;
}

template<typename T>
ReductiveBoxProcessorGenerator3D<T>::ReductiveBoxProcessorGenerator3D(ReductiveBoxProcessorGenerator3D<T> const& rhs)
    : BoxedReductiveDataProcessorGenerator3D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
ReductiveBoxProcessorGenerator3D<T>& ReductiveBoxProcessorGenerator3D<T>::operator= (
        ReductiveBoxProcessorGenerator3D<T> const& rhs )
{
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT ReductiveBoxProcessorGenerator3D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void ReductiveBoxProcessorGenerator3D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void ReductiveBoxProcessorGenerator3D<T>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor3D<T>* ReductiveBoxProcessorGenerator3D<T>::generate (
        std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    // Don't clone functional. Given that the functional contains the BlockStatistics object,
    //   everybody must poplint to the same instance.
    return new ReductiveBoxProcessor3D<T>(functional, this->getDomain(), atomicBlocks);
}

template<typename T>
ReductiveBoxProcessorGenerator3D<T>* ReductiveBoxProcessorGenerator3D<T>::clone() const {
    return new ReductiveBoxProcessorGenerator3D<T>(*this);
}

template<typename T>
BlockStatistics<T> const& ReductiveBoxProcessorGenerator3D<T>::getStatistics() const {
    return functional->getStatistics();
}

template<typename T>
BlockStatistics<T>& ReductiveBoxProcessorGenerator3D<T>::getStatistics() {
    return functional->getStatistics();
}

template<typename T>
ReductiveBoxProcessingFunctional3D<T> const& ReductiveBoxProcessorGenerator3D<T>::getFunctional() const {
    return *functional;
}


/* *************** ReductiveBoxProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveBoxProcessingFunctional3D_L<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_L<T,Descriptor>& functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D<T> generator(functional.clone(), domain);
    lattice.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing3D_S ******************************************* */

template<typename T>
void ReductiveBoxProcessingFunctional3D_S<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    process(domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_S<T>& functional,
                               Box3D domain, ScalarFieldBase3D<T>& field)
{

    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D<T> generator(functional.clone(), domain);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing3D_T ******************************************* */

template<typename T, int nDim>
void ReductiveBoxProcessingFunctional3D_T<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D_T<T,nDim>& functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D<T> generator(functional.clone(), domain);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void ReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void ReductiveBoxProcessingFunctional3D_SS<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void ReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void ReductiveBoxProcessingFunctional3D_ST<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void ReductiveBoxProcessingFunctional3D_LS<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void ReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveLatticeBoxProcessing3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ReductiveScalarFieldBoxProcessing3D ******************************************* */

template<typename T>
void ReductiveScalarFieldBoxProcessingFunctional3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks)
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveTensorFieldBoxProcessing3D ******************************************* */

template<typename T, int nDim>
void ReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** Class ReductiveDotProcessingFunctional3D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT ReductiveDotProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void ReductiveDotProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveDotProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}



/* *************** Class PlainReductiveDotProcessingFunctional3D ************************* */

template<typename T>
BlockStatistics<T> const& PlainReductiveDotProcessingFunctional3D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& PlainReductiveDotProcessingFunctional3D<T>::getStatistics() {
    return statistics;
}


/* *************** Class ReductiveDotProcessor3D ************************************ */

template<typename T>
ReductiveDotProcessor3D<T>::ReductiveDotProcessor3D(ReductiveDotProcessingFunctional3D<T>* functional_,
               DotList3D const& dotList_, std::vector<AtomicBlock3D<T>*> atomicBlocks_)
    : functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
DotList3D const& ReductiveDotProcessor3D<T>::getDotList() const {
    return dotList;
}

template<typename T>
void ReductiveDotProcessor3D<T>::process() {
    functional -> processGenericBlocks(dotList, atomicBlocks);
}

template<typename T>
ReductiveDotProcessor3D<T>* ReductiveDotProcessor3D<T>::clone() const {
    return new ReductiveDotProcessor3D<T>(*this);
}


/* *************** Class ReductiveDotProcessorGenerator3D *************************** */

template<typename T>
ReductiveDotProcessorGenerator3D<T>::ReductiveDotProcessorGenerator3D (
        ReductiveDotProcessingFunctional3D<T>* functional_,
        DotList3D const& dotList )
    : DottedReductiveDataProcessorGenerator3D<T>(dotList),
      functional(functional_)
{ }

template<typename T>
ReductiveDotProcessorGenerator3D<T>::~ReductiveDotProcessorGenerator3D() {
    delete functional;
}

template<typename T>
ReductiveDotProcessorGenerator3D<T>::ReductiveDotProcessorGenerator3D (
        ReductiveDotProcessorGenerator3D<T> const& rhs )
    : DottedReductiveDataProcessorGenerator3D<T>(*this),
      functional(rhs.functional->clone())
{ }

template<typename T>
ReductiveDotProcessorGenerator3D<T>& ReductiveDotProcessorGenerator3D<T>::operator= (
        ReductiveDotProcessorGenerator3D<T> const& rhs )
{
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT ReductiveDotProcessorGenerator3D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void ReductiveDotProcessorGenerator3D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void ReductiveDotProcessorGenerator3D<T>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    functional->getModificationPattern(isWritten);
}


template<typename T>
DataProcessor3D<T>* ReductiveDotProcessorGenerator3D<T>::generate (
        std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    return new ReductiveDotProcessor3D<T>(functional, this->getDotList(), atomicBlocks);
}

template<typename T>
ReductiveDotProcessorGenerator3D<T>* ReductiveDotProcessorGenerator3D<T>::clone() const {
    return new ReductiveDotProcessorGenerator3D<T>(*this);
}

template<typename T>
BlockStatistics<T> const& ReductiveDotProcessorGenerator3D<T>::getStatistics() const {
    return functional->getStatistics();
}

template<typename T>
BlockStatistics<T>& ReductiveDotProcessorGenerator3D<T>::getStatistics() {
    return functional->getStatistics();
}

template<typename T>
ReductiveDotProcessingFunctional3D<T> const& ReductiveDotProcessorGenerator3D<T>::getFunctional() const {
    return *functional;
}


/* *************** ReductiveDotProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveDotProcessingFunctional3D_L<T,Descriptor>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_L<T,Descriptor>& functional,
                               DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D<T> generator(functional.clone(), dotList);
    lattice.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing3D_S ******************************************* */

template<typename T>
void ReductiveDotProcessingFunctional3D_S<T>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_S<T>& functional,
                               DotList3D const& dotList, ScalarFieldBase3D<T>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D<T> generator(functional.clone(), dotList);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing3D_T ******************************************* */

template<typename T, int nDim>
void ReductiveDotProcessingFunctional3D_T<T,nDim>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional3D_T<T,nDim>& functional,
                               DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D<T> generator(functional.clone(), dotList);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void ReductiveDotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing3D_SS ****************************************** */

template<typename T>
void ReductiveDotProcessingFunctional3D_SS<T>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void ReductiveDotProcessingFunctional3D_TT<T,nDim1,nDim2>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void ReductiveDotProcessingFunctional3D_ST<T,nDim>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void ReductiveDotProcessingFunctional3D_LS<T,Descriptor>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void ReductiveDotProcessingFunctional3D_LT<T,Descriptor,nDim>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveLatticeDotProcessing3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveLatticeDotProcessingFunctional3D<T,Descriptor>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ReductiveScalarFieldDotProcessing3D ******************************************* */

template<typename T>
void ReductiveScalarFieldDotProcessingFunctional3D<T>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveTensorFieldDotProcessing3D ******************************************* */

template<typename T, int nDim>
void ReductiveTensorFieldDotProcessingFunctional3D<T,nDim>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}


/* *************** Class BoundedReductiveBoxProcessingFunctional3D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}

template<typename T>
ReductiveBoxProcessingFunctional3D<T>* BoundedReductiveBoxProcessingFunctional3D<T>::getBulkProcessor() const {
    return new BulkWrapperFunctional(this->clone());
}

template<typename T>
ReductiveBoxProcessingFunctional3D<T>* BoundedReductiveBoxProcessingFunctional3D<T>::getPlaneProcessor (
                                  int direction, int orientation ) const
{
    return new PlaneWrapperFunctional(this->clone(), direction, orientation);
}

template<typename T>
ReductiveBoxProcessingFunctional3D<T>* BoundedReductiveBoxProcessingFunctional3D<T>::getEdgeProcessor (
                                  int plane, int normal1, int normal2 ) const
{
    return new EdgeWrapperFunctional(this->clone(), plane, normal1, normal2);
}

template<typename T>
ReductiveBoxProcessingFunctional3D<T>* BoundedReductiveBoxProcessingFunctional3D<T>::getCornerProcessor (
                                  int normalX, int normalY, int normalZ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY, normalZ);
}

template<typename T>
BlockStatistics<T> const& BoundedReductiveBoxProcessingFunctional3D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& BoundedReductiveBoxProcessingFunctional3D<T>::getStatistics() {
    return statistics;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_ )
    : boundedFunctional(boundedFunctional_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BulkWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone())
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::~BulkWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::operator= (
            BulkWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional3D<T>::BulkWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::PlaneWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
        int direction_, int orientation_)
    : boundedFunctional(boundedFunctional_),
      direction(direction_),
      orientation(orientation_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::PlaneWrapperFunctional (
        PlaneWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      direction(rhs.direction),
      orientation(rhs.orientation)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::~PlaneWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::operator= (
            PlaneWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processPlaneGeneric(direction, orientation, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::clone() const
{
    return new PlaneWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
        int plane_, int normal1_, int normal2_)
    : boundedFunctional(boundedFunctional_),
      plane(plane_),
      normal1(normal1_),
      normal2(normal2_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        EdgeWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      plane(rhs.plane),
      normal1(rhs.normal1),
      normal2(rhs.normal2)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::~EdgeWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::operator= (
            EdgeWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    plane = rhs.plane;
    normal1 = rhs.normal1;
    normal2 = rhs.normal2;
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processEdgeGeneric(plane, normal1, normal2, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional3D<T>* boundedFunctional_,
        int normalX_,
        int normalY_,
        int normalZ_)
    : boundedFunctional(boundedFunctional_),
      normalX(normalX_),
      normalY(normalY_),
      normalZ(normalZ_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        CornerWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      normalX(rhs.normalX),
      normalY(rhs.normalY),
      normalZ(rhs.normalZ)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::~CornerWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::operator= (
            CornerWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    normalZ = rhs.normalZ;
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processCornerGeneric(normalX, normalY, normalZ, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional3D<T>::CornerWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveOneBlockProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedReductiveOneBlockProcessingFunctionalOperation3D<T>::BoundedReductiveOneBlockProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedReductiveOneBlockProcessingFunctionalOperation3D<T>::apply (
        BoundedReductiveBoxProcessingFunctional3D<T>& functional, Block3D<T>& block )
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
        block.executeDataProcessor(*generators[iGenerator]);
        individualStatistics[iGenerator] = &(generators[iGenerator]->getStatistics());
    }

    SerialCombinedStatistics<T>().combine(individualStatistics, functional.getStatistics());

    for (pluint iGenerator=0; iGenerator<generators.size(); ++iGenerator) {
        delete generators[iGenerator];
    }
}


/* *************** BoundedReductiveBoxProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane, normal1, normal2, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_L<T,Descriptor>& functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}


/* *************** BoundedReductiveBoxProcessing3D_S ******************************************* */

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane, normal1, normal2, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_S<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_S<T>& functional,
                               Box3D domain, ScalarFieldBase3D<T>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}


/* *************** BoundedReductiveBoxProcessing3D_T ******************************************* */

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane, normal1, normal2, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D_T<T,nDim>& functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}

/* *************** BoundedReductiveBoxProcessing3D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                   dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing3D_SS ****************************************** */

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_SS<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_SS<T>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_SS<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional3D_SS<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_ST<T,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional3D_LS<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane, normal1, normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveLatticeBoxProcessing3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processPlane(direction, orientation, domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processEdge(plane, normal1, normal2, domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional3D<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, normalZ, domain, lattices);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing3D ******************************************* */

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processPlane(direction, orientation, domain, fields);
}

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processEdge(plane, normal1, normal2, domain, fields);
}

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional3D<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing3D ******************************************* */

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processPlane(direction, orientation, domain, fields);
}

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processEdge(plane, normal1, normal2, domain, fields);
}

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional3D<T,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSOR_WRAPPER_3D_HH
