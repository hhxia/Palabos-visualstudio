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


#ifndef REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH
#define REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH

#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"
#include "multiBlock/combinedStatistics.h"

namespace plb {

/* *************** Class ReductiveBoxProcessingFunctional2D ************************* */

/** Operation is not applied to envelope by default. **/
template<typename T>
BlockDomain::DomainT ReductiveBoxProcessingFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void ReductiveBoxProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveBoxProcessingFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}

/* *************** Class PlainReductiveBoxProcessingFunctional2D ************************* */

template<typename T>
BlockStatistics<T> const& PlainReductiveBoxProcessingFunctional2D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& PlainReductiveBoxProcessingFunctional2D<T>::getStatistics() {
    return statistics;
}


/* *************** Class ReductiveBoxProcessor2D ************************************ */

template<typename T>
ReductiveBoxProcessor2D<T>::ReductiveBoxProcessor2D (
        ReductiveBoxProcessingFunctional2D<T>* functional_,
        Box2D domain_, std::vector<AtomicBlock2D<T>*> atomicBlocks_)
    : functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
Box2D ReductiveBoxProcessor2D<T>::getDomain() const {
    return domain;
}

template<typename T>
void ReductiveBoxProcessor2D<T>::process() {
    functional -> processGenericBlocks(domain, atomicBlocks);
    functional -> getStatistics().evaluate();
}

template<typename T>
ReductiveBoxProcessor2D<T>* ReductiveBoxProcessor2D<T>::clone() const {
    return new ReductiveBoxProcessor2D<T>(*this);
}


/* *************** Class ReductiveBoxProcessorGenerator2D *************************** */

template<typename T>
ReductiveBoxProcessorGenerator2D<T>::ReductiveBoxProcessorGenerator2D (
        ReductiveBoxProcessingFunctional2D<T>* functional_,
        Box2D domain )
    : BoxedReductiveDataProcessorGenerator2D<T>(domain),
      functional(functional_)
{ }

template<typename T>
ReductiveBoxProcessorGenerator2D<T>::~ReductiveBoxProcessorGenerator2D() {
    delete functional;
}

template<typename T>
ReductiveBoxProcessorGenerator2D<T>::ReductiveBoxProcessorGenerator2D(ReductiveBoxProcessorGenerator2D<T> const& rhs)
    : BoxedReductiveDataProcessorGenerator2D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
ReductiveBoxProcessorGenerator2D<T>& ReductiveBoxProcessorGenerator2D<T>::operator= (
        ReductiveBoxProcessorGenerator2D<T> const& rhs )
{
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT ReductiveBoxProcessorGenerator2D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void ReductiveBoxProcessorGenerator2D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void ReductiveBoxProcessorGenerator2D<T>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor2D<T>* ReductiveBoxProcessorGenerator2D<T>::generate (
        std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    // Don't clone functional. Given that the functional contains the BlockStatistics object,
    //   everybody must point to the same instance.
    return new ReductiveBoxProcessor2D<T>(functional, this->getDomain(), atomicBlocks);
}

template<typename T>
ReductiveBoxProcessorGenerator2D<T>* ReductiveBoxProcessorGenerator2D<T>::clone() const {
    return new ReductiveBoxProcessorGenerator2D<T>(*this);
}

template<typename T>
BlockStatistics<T> const& ReductiveBoxProcessorGenerator2D<T>::getStatistics() const {
    return functional->getStatistics();
}

template<typename T>
BlockStatistics<T>& ReductiveBoxProcessorGenerator2D<T>::getStatistics() {
    return functional->getStatistics();
}

template<typename T>
ReductiveBoxProcessingFunctional2D<T> const& ReductiveBoxProcessorGenerator2D<T>::getFunctional() const {
    return *functional;
}


/* *************** ReductiveBoxProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveBoxProcessingFunctional2D_L<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_L<T,Descriptor>& functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator2D<T> generator(functional.clone(), domain);
    lattice.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing2D_S ******************************************* */

template<typename T>
void ReductiveBoxProcessingFunctional2D_S<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks)
{
    process(domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_S<T>& functional,
                               Box2D domain, ScalarFieldBase2D<T>& field)
{

    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator2D<T> generator(functional.clone(), domain);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing2D_T ******************************************* */

template<typename T, int nDim>
void ReductiveBoxProcessingFunctional2D_T<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveBoxProcessingFunctional2D_T<T,nDim>& functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator2D<T> generator(functional.clone(), domain);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveBoxProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void ReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void ReductiveBoxProcessingFunctional2D_SS<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void ReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void ReductiveBoxProcessingFunctional2D_ST<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void ReductiveBoxProcessingFunctional2D_LS<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void ReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveLatticeBoxProcessing2D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ReductiveScalarFieldBoxProcessing2D ******************************************* */

template<typename T>
void ReductiveScalarFieldBoxProcessingFunctional2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks)
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** ReductiveTensorFieldBoxProcessing2D ******************************************* */

template<typename T, int nDim>
void ReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** Class ReductiveDotProcessingFunctional2D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT ReductiveDotProcessingFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void ReductiveDotProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void ReductiveDotProcessingFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}


/* *************** Class PlainReductiveDotProcessingFunctional2D ************************* */

template<typename T>
BlockStatistics<T> const& PlainReductiveDotProcessingFunctional2D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& PlainReductiveDotProcessingFunctional2D<T>::getStatistics() {
    return statistics;
}


/* *************** Class ReductiveDotProcessor2D ************************************ */

template<typename T>
ReductiveDotProcessor2D<T>::ReductiveDotProcessor2D(ReductiveDotProcessingFunctional2D<T>* functional_,
               DotList2D const& dotList_, std::vector<AtomicBlock2D<T>*> atomicBlocks_)
    : functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
DotList2D const& ReductiveDotProcessor2D<T>::getDotList() const {
    return dotList;
}

template<typename T>
void ReductiveDotProcessor2D<T>::process() {
    functional -> processGenericBlocks(dotList, atomicBlocks);
}

template<typename T>
ReductiveDotProcessor2D<T>* ReductiveDotProcessor2D<T>::clone() const {
    return new ReductiveDotProcessor2D<T>(*this);
}


/* *************** Class ReductiveDotProcessorGenerator2D *************************** */

template<typename T>
ReductiveDotProcessorGenerator2D<T>::ReductiveDotProcessorGenerator2D (
        ReductiveDotProcessingFunctional2D<T>* functional_,
        DotList2D const& dotList )
    : DottedReductiveDataProcessorGenerator2D<T>(dotList),
      functional(functional_)
{ }

template<typename T>
ReductiveDotProcessorGenerator2D<T>::~ReductiveDotProcessorGenerator2D() {
    delete functional;
}

template<typename T>
ReductiveDotProcessorGenerator2D<T>::ReductiveDotProcessorGenerator2D (
        ReductiveDotProcessorGenerator2D<T> const& rhs )
    : DottedReductiveDataProcessorGenerator2D<T>(*this),
      functional(rhs.functional->clone())
{ }

template<typename T>
ReductiveDotProcessorGenerator2D<T>& ReductiveDotProcessorGenerator2D<T>::operator= (
        ReductiveDotProcessorGenerator2D<T> const& rhs )
{
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT ReductiveDotProcessorGenerator2D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void ReductiveDotProcessorGenerator2D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void ReductiveDotProcessorGenerator2D<T>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor2D<T>* ReductiveDotProcessorGenerator2D<T>::generate (
        std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    return new ReductiveDotProcessor2D<T>(functional, this->getDotList(), atomicBlocks);
}

template<typename T>
ReductiveDotProcessorGenerator2D<T>* ReductiveDotProcessorGenerator2D<T>::clone() const {
    return new ReductiveDotProcessorGenerator2D<T>(*this);
}

template<typename T>
BlockStatistics<T> const& ReductiveDotProcessorGenerator2D<T>::getStatistics() const {
    return functional->getStatistics();
}

template<typename T>
BlockStatistics<T>& ReductiveDotProcessorGenerator2D<T>::getStatistics() {
    return functional->getStatistics();
}

template<typename T>
ReductiveDotProcessingFunctional2D<T> const& ReductiveDotProcessorGenerator2D<T>::getFunctional() const {
    return *functional;
}


/* *************** ReductiveDotProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveDotProcessingFunctional2D_L<T,Descriptor>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_L<T,Descriptor>& functional,
                               DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator2D<T> generator(functional.clone(), dotList);
    lattice.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing2D_S ******************************************* */

template<typename T>
void ReductiveDotProcessingFunctional2D_S<T>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_S<T>& functional,
                               DotList2D const& dotList, ScalarFieldBase2D<T>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator2D<T> generator(functional.clone(), dotList);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing2D_T ******************************************* */

template<typename T, int nDim>
void ReductiveDotProcessingFunctional2D_T<T,nDim>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(dotList, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(ReductiveDotProcessingFunctional2D_T<T,nDim>& functional,
                               DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator2D<T> generator(functional.clone(), dotList);
    field.executeDataProcessor(generator);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** ReductiveDotProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void ReductiveDotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing2D_SS ****************************************** */

template<typename T>
void ReductiveDotProcessingFunctional2D_SS<T>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void ReductiveDotProcessingFunctional2D_TT<T,nDim1,nDim2>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void ReductiveDotProcessingFunctional2D_ST<T,nDim>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void ReductiveDotProcessingFunctional2D_LS<T,Descriptor>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** ReductiveDotProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void ReductiveDotProcessingFunctional2D_LT<T,Descriptor,nDim>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** ReductiveLatticeDotProcessing2D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void ReductiveLatticeDotProcessingFunctional2D<T,Descriptor>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ReductiveScalarFieldDotProcessing2D ******************************************* */

template<typename T>
void ReductiveScalarFieldDotProcessingFunctional2D<T>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** ReductiveTensorFieldDotProcessing2D ******************************************* */

template<typename T, int nDim>
void ReductiveTensorFieldDotProcessingFunctional2D<T,nDim>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}


/* *************** Class BoundedReductiveBoxProcessingFunctional2D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** By default, it is assumed that none of the blocks are written, they
 *  are only read.
 */
template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = false;
    }
}


template<typename T>
ReductiveBoxProcessingFunctional2D<T>* BoundedReductiveBoxProcessingFunctional2D<T>::getBulkProcessor() const {
    return new BulkWrapperFunctional(this->clone());
}

template<typename T>
ReductiveBoxProcessingFunctional2D<T>* BoundedReductiveBoxProcessingFunctional2D<T>::getEdgeProcessor (
                                  int direction, int orientation ) const
{
    return new EdgeWrapperFunctional(this->clone(), direction, orientation);
}

template<typename T>
ReductiveBoxProcessingFunctional2D<T>* BoundedReductiveBoxProcessingFunctional2D<T>::getCornerProcessor (
                                  int normalX, int normalY ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY);
}

template<typename T>
BlockStatistics<T> const& BoundedReductiveBoxProcessingFunctional2D<T>::getStatistics() const {
    return statistics;
}

template<typename T>
BlockStatistics<T>& BoundedReductiveBoxProcessingFunctional2D<T>::getStatistics() {
    return statistics;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_ )
    : boundedFunctional(boundedFunctional_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BulkWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone())
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::~BulkWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::operator= (
            BulkWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional2D<T>::BulkWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_,
        int direction_, int orientation_)
    : boundedFunctional(boundedFunctional_),
      direction(direction_),
      orientation(orientation_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        EdgeWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      direction(rhs.direction),
      orientation(rhs.orientation)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::~EdgeWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::operator= (
            EdgeWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processEdgeGeneric(direction, orientation, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional ** */

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        BoundedReductiveBoxProcessingFunctional2D<T>* boundedFunctional_,
        int normalX_,
        int normalY_ )
    : boundedFunctional(boundedFunctional_),
      normalX(normalX_),
      normalY(normalY_)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        CornerWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      normalX(rhs.normalX),
      normalY(rhs.normalY)
{ }

template<typename T>
BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::~CornerWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional&
    BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::operator= (
            CornerWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    return *this;
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processCornerGeneric(normalX, normalY, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::
         getModificationPattern(std::vector<bool>& isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional*
    BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

template<typename T>
BlockStatistics<T> const&
    BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::getStatistics() const
{
    return boundedFunctional -> getStatistics();
}

template<typename T>
BlockStatistics<T>&
    BoundedReductiveBoxProcessingFunctional2D<T>::CornerWrapperFunctional::getStatistics()
{
    return boundedFunctional -> getStatistics();
}


/* *************** Class BoundedReductiveOneBlockProcessingFunctionalOperation2D ******** */

template<typename T>
BoundedReductiveOneBlockProcessingFunctionalOperation2D<T>::BoundedReductiveOneBlockProcessingFunctionalOperation2D (
        Box2D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
    { }

template<typename T>
void BoundedReductiveOneBlockProcessingFunctionalOperation2D<T>::apply (
        BoundedReductiveBoxProcessingFunctional2D<T>& functional, Block2D<T>& block )
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
        block.executeDataProcessor(*generators[iGenerator]);
        individualStatistics[iGenerator] = &(generators[iGenerator]->getStatistics());
    }

    SerialCombinedStatistics<T>().combine(individualStatistics, functional.getStatistics());

    for (pluint iGenerator=0; iGenerator<generators.size(); ++iGenerator) {
        delete generators[iGenerator];
    }
}


/* *************** BoundedReductiveBoxProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T,Descriptor>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_L<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_L<T,Descriptor>& functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}


/* *************** BoundedReductiveBoxProcessing2D_S ******************************************* */

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_S<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_S<T>& functional,
                               Box2D domain, ScalarFieldBase2D<T>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}


/* *************** BoundedReductiveBoxProcessing2D_T ******************************************* */

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_T<T,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional2D_T<T,nDim>& functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& lattice, plint boundaryWidth)
{
    BoundedReductiveOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}

/* *************** BoundedReductiveBoxProcessing2D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedReductiveBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                   dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing2D_SS ****************************************** */

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_SS<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_SS<T>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedReductiveBoxProcessingFunctional2D_SS<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedReductiveBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_ST<T,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveBoxProcessingFunctional2D_LS<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedReductiveBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedReductiveLatticeBoxProcessing2D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processEdge(direction, orientation, domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedReductiveLatticeBoxProcessingFunctional2D<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, domain, lattices);
}


/* *************** BoundedReductiveScalarFieldBoxProcessing2D ******************************************* */

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    processEdge(direction, orientation, domain, fields);
}

template<typename T>
void BoundedReductiveScalarFieldBoxProcessingFunctional2D<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}


/* *************** BoundedReductiveTensorFieldBoxProcessing2D ******************************************* */

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    processEdge(direction, orientation, domain, fields);
}

template<typename T, int nDim>
void BoundedReductiveTensorFieldBoxProcessingFunctional2D<T,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}

}  // namespace plb

#endif  // REDUCTIVE_DATA_PROCESSOR_WRAPPER_2D_HH
