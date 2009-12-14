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


#ifndef DATA_PROCESSOR_WRAPPER_2D_HH
#define DATA_PROCESSOR_WRAPPER_2D_HH

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class BoxProcessingFunctional2D ************************* */

/** Operations are not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT BoxProcessingFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoxProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void BoxProcessingFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}

/* *************** Class BoxProcessor2D ************************************ */

template<typename T>
BoxProcessor2D<T>::BoxProcessor2D(BoxProcessingFunctional2D<T>* functional_,
               Box2D domain_, std::vector<AtomicBlock2D<T>*> atomicBlocks_)
    : functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
BoxProcessor2D<T>::BoxProcessor2D(BoxProcessor2D<T> const& rhs)
    : functional(rhs.functional->clone()),
      domain(rhs.domain), atomicBlocks(rhs.atomicBlocks)
{ }

template<typename T>
BoxProcessor2D<T>& BoxProcessor2D<T>::operator=(BoxProcessor2D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    domain = rhs.domain;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

template<typename T>
BoxProcessor2D<T>::~BoxProcessor2D() {
    delete functional;
}

template<typename T>
Box2D BoxProcessor2D<T>::getDomain() const {
    return domain;
}

template<typename T>
void BoxProcessor2D<T>::process() {
    functional -> processGenericBlocks(domain, atomicBlocks);
}

template<typename T>
BoxProcessor2D<T>* BoxProcessor2D<T>::clone() const {
    return new BoxProcessor2D<T>(*this);
}


/* *************** Class BoxProcessorGenerator2D *************************** */

template<typename T>
BoxProcessorGenerator2D<T>::BoxProcessorGenerator2D(BoxProcessingFunctional2D<T>* functional_, Box2D domain)
    : BoxedDataProcessorGenerator2D<T>(domain),
      functional(functional_)
{ }

template<typename T>
BoxProcessorGenerator2D<T>::~BoxProcessorGenerator2D() {
    delete functional;
}

template<typename T>
BoxProcessorGenerator2D<T>::BoxProcessorGenerator2D(BoxProcessorGenerator2D<T> const& rhs)
    : BoxedDataProcessorGenerator2D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
BoxProcessorGenerator2D<T>& BoxProcessorGenerator2D<T>::operator=(BoxProcessorGenerator2D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT BoxProcessorGenerator2D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void BoxProcessorGenerator2D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void BoxProcessorGenerator2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor2D<T>* BoxProcessorGenerator2D<T>::generate(std::vector<AtomicBlock2D<T>*> atomicBlocks) const {
    return new BoxProcessor2D<T>(functional->clone(), this->getDomain(), atomicBlocks);
}

template<typename T>
BoxProcessorGenerator2D<T>* BoxProcessorGenerator2D<T>::clone() const {
    return new BoxProcessorGenerator2D<T>(*this);
}


/* *************** BoxProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoxProcessingFunctional2D_L<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional2D_L<T,Descriptor>* functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    lattice.executeDataProcessor(BoxProcessorGenerator2D<T>(functional, domain));
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional2D_L<T,Descriptor>* functional,
                                   Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, plint level)
{
    lattice.addInternalProcessor(BoxProcessorGenerator2D<T>(functional, domain), level);
}


/* *************** BoxProcessing2D_S ******************************************* */

template<typename T>
void BoxProcessingFunctional2D_S<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional2D_S<T>* functional,
                               Box2D domain, ScalarFieldBase2D<T>& field)
{
    field.executeDataProcessor(BoxProcessorGenerator2D<T>(functional, domain));
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional2D_S<T>* functional,
                                   Box2D domain, ScalarFieldBase2D<T>& field, plint level)
{
    field.addInternalProcessor(BoxProcessorGenerator2D<T>(functional, domain), level);
}


/* *************** BoxProcessing2D_T ******************************************* */

template<typename T, int nDim>
void BoxProcessingFunctional2D_T<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional2D_T<T,nDim>* functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& field)
{
    field.executeDataProcessor(BoxProcessorGenerator2D<T>(functional, domain));
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional2D_T<T,nDim>* functional,
                                   Box2D domain, TensorFieldBase2D<T,nDim>& field, plint level)
{
    field.addInternalProcessor(BoxProcessorGenerator2D<T>(functional, domain), level);
}


/* *************** BoxProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing2D_SS ****************************************** */

template<typename T>
void BoxProcessingFunctional2D_SS<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void BoxProcessingFunctional2D_ST<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoxProcessingFunctional2D_LS<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** LatticeBoxProcessing2D ************************************** */

template<typename T, template<typename U> class Descriptor>
void LatticeBoxProcessingFunctional2D<T,Descriptor>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ScalarFieldBoxProcessing2D *********************************** */

template<typename T>
void ScalarFieldBoxProcessingFunctional2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** TensorFieldBoxProcessing2D *********************************** */

template<typename T, int nDim>
void TensorFieldBoxProcessingFunctional2D<T,nDim>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}


/* *************** Class DotProcessingFunctional2D ************************* */

/** Operation is not applied to envelope by default. **/
template<typename T>
BlockDomain::DomainT DotProcessingFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void DotProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void DotProcessingFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}



/* *************** Class DotProcessor2D ************************************ */

template<typename T>
DotProcessor2D<T>::DotProcessor2D(DotProcessingFunctional2D<T>* functional_,
               DotList2D const& dotList_, std::vector<AtomicBlock2D<T>*> atomicBlocks_)
    : functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
DotProcessor2D<T>::DotProcessor2D(DotProcessor2D<T> const& rhs)
    : functional(rhs.functional->clone()),
      dotList(rhs.dotList), atomicBlocks(rhs.atomicBlocks)
{ }

template<typename T>
DotProcessor2D<T>& DotProcessor2D<T>::operator=(DotProcessor2D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    dotList = rhs.dotList;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

template<typename T>
DotProcessor2D<T>::~DotProcessor2D() {
    delete functional;
}

template<typename T>
void DotProcessor2D<T>::process() {
    functional -> processGenericBlocks(dotList, atomicBlocks);
}

template<typename T>
DotProcessor2D<T>* DotProcessor2D<T>::clone() const {
    return new DotProcessor2D<T>(*this);
}

template<typename T>
DotList2D const& DotProcessor2D<T>::getDotList() const {
    return dotList;
}


/* *************** Class DotProcessorGenerator2D *************************** */

template<typename T>
DotProcessorGenerator2D<T>::DotProcessorGenerator2D(DotProcessingFunctional2D<T>* functional_, DotList2D const& dotList)
    : DottedDataProcessorGenerator2D<T>(dotList),
      functional(functional_)
{ }

template<typename T>
DotProcessorGenerator2D<T>::~DotProcessorGenerator2D() {
    delete functional;
}

template<typename T>
DotProcessorGenerator2D<T>::DotProcessorGenerator2D(DotProcessorGenerator2D<T> const& rhs)
    : DottedDataProcessorGenerator2D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
DotProcessorGenerator2D<T>& DotProcessorGenerator2D<T>::operator=(DotProcessorGenerator2D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT DotProcessorGenerator2D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void DotProcessorGenerator2D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void DotProcessorGenerator2D<T>::getModificationPattern(std::vector<bool>& isWritten) const
{
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor2D<T>* DotProcessorGenerator2D<T>::generate(std::vector<AtomicBlock2D<T>*> atomicBlocks) const {
    return new DotProcessor2D<T>(functional->clone(), this->getDotList(), atomicBlocks);
}

template<typename T>
DotProcessorGenerator2D<T>* DotProcessorGenerator2D<T>::clone() const {
    return new DotProcessorGenerator2D<T>(*this);
}


/* *************** DotProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void DotProcessingFunctional2D_L<T,Descriptor>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional2D_L<T,Descriptor>* functional,
                               DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    lattice.executeDataProcessor(DotProcessorGenerator2D<T>(functional, dotList));
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional2D_L<T,Descriptor>* functional,
                                   DotList2D const& dotList, BlockLatticeBase2D<T,Descriptor>& lattice, plint level)
{
    lattice.addInternalProcessor(DotProcessorGenerator2D<T>(functional, dotList), level);
}


/* *************** DotProcessing2D_S ******************************************* */

template<typename T>
void DotProcessingFunctional2D_S<T>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional2D_S<T>* functional,
                               DotList2D const& dotList, ScalarFieldBase2D<T>& field)
{
    field.executeDataProcessor(DotProcessorGenerator2D<T>(functional, dotList));
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional2D_S<T>* functional,
                                   DotList2D const& dotList, ScalarFieldBase2D<T>& field, plint level)
{
    field.addInternalProcessor(DotProcessorGenerator2D<T>(functional, dotList), level);
}


/* *************** DotProcessing2D_T ******************************************* */

template<typename T, int nDim>
void DotProcessingFunctional2D_T<T,nDim>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional2D_T<T,nDim>* functional,
                               DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field)
{
    field.executeDataProcessor(DotProcessorGenerator2D<T>(functional, dotList));
}

template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional2D_T<T,nDim>* functional,
                                   DotList2D const& dotList, TensorFieldBase2D<T,nDim>& field, plint level)
{
    field.addInternalProcessor(DotProcessorGenerator2D<T>(functional, dotList), level);
}


/* *************** DotProcessing2D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void DotProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks ) {
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing2D_SS ****************************************** */

template<typename T>
void DotProcessingFunctional2D_SS<T>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void DotProcessingFunctional2D_TT<T,nDim1,nDim2>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void DotProcessingFunctional2D_ST<T,nDim>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void DotProcessingFunctional2D_LS<T,Descriptor>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void DotProcessingFunctional2D_LT<T,Descriptor,nDim>::processGenericBlocks(DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** LatticeDotProcessing2D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void LatticeDotProcessingFunctional2D<T,Descriptor>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ScalarFieldDotProcessing2D ******************************************* */

template<typename T>
void ScalarFieldDotProcessingFunctional2D<T>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** TensorFieldDotProcessing2D ******************************************* */

template<typename T, int nDim>
void TensorFieldDotProcessingFunctional2D<T,nDim>::processGenericBlocks (
        DotList2D const& dotList, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}


/* *************** Class BoundedBoxProcessingFunctional2D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoundedBoxProcessingFunctional2D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void BoundedBoxProcessingFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}

template<typename T>
BoxProcessingFunctional2D<T>* BoundedBoxProcessingFunctional2D<T>::getBulkProcessor() const {
    return new BulkWrapperFunctional(this->clone());
}

template<typename T>
BoxProcessingFunctional2D<T>* BoundedBoxProcessingFunctional2D<T>::getEdgeProcessor (
                                  int direction, int orientation ) const
{
    return new EdgeWrapperFunctional(this->clone(), direction, orientation);
}

template<typename T>
BoxProcessingFunctional2D<T>* BoundedBoxProcessingFunctional2D<T>::getCornerProcessor (
                                  int normalX, int normalY ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY);
}

/* *************** Class BoundedBoxProcessingFunctional2D::BulkWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional_ )
    : boundedFunctional(boundedFunctional_)
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BulkWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone())
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::~BulkWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional&
    BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::operator= (
            BulkWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional*
    BoundedBoxProcessingFunctional2D<T>::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional_,
        int direction_, int orientation_)
    : boundedFunctional(boundedFunctional_),
      direction(direction_),
      orientation(orientation_)
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        EdgeWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      direction(rhs.direction),
      orientation(rhs.orientation)
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::~EdgeWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional&
    BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::operator= (
            EdgeWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processEdgeGeneric(direction, orientation, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional*
    BoundedBoxProcessingFunctional2D<T>::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional2D::CornerWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        BoundedBoxProcessingFunctional2D<T>* boundedFunctional_,
        int normalX_,
        int normalY_ )
    : boundedFunctional(boundedFunctional_),
      normalX(normalX_),
      normalY(normalY_)
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        CornerWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      normalX(rhs.normalX),
      normalY(rhs.normalY)
{ }

template<typename T>
BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::~CornerWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional&
    BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::operator= (
            CornerWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    return *this;
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    boundedFunctional->processCornerGeneric(normalX, normalY, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional*
    BoundedBoxProcessingFunctional2D<T>::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}


/* *************** Class BoundedOneBlockProcessingFunctionalOperation2D ******** */

template<typename T>
BoundedOneBlockProcessingFunctionalOperation2D<T>::BoundedOneBlockProcessingFunctionalOperation2D (
        Box2D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedOneBlockProcessingFunctionalOperation2D<T>::apply (
        BoundedBoxProcessingFunctional2D<T>* functional, Block2D<T>& block )
{
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getBulkProcessor(), surf.bulk()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,-1), surf.edge0N()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,+1), surf.edge0P()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,-1), surf.edge1N()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,+1), surf.edge1P()));

    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,-1), surf.cornerNN()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,-1), surf.cornerPN()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,+1), surf.cornerNP()));
    block.executeDataProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,+1), surf.cornerPP()));
    delete functional;
}

template<typename T>
void BoundedOneBlockProcessingFunctionalOperation2D<T>::integrate (
        BoundedBoxProcessingFunctional2D<T>* functional, Block2D<T>& block, plint level )
{
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getBulkProcessor(), surf.bulk()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,-1), surf.edge0N()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(0,+1), surf.edge0P()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,-1), surf.edge1N()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getEdgeProcessor(1,+1), surf.edge1P()), level);

    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,-1), surf.cornerNN()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,-1), surf.cornerPN()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(-1,+1), surf.cornerNP()), level);
    block.addInternalProcessor(BoxProcessorGenerator2D<T>(functional->getCornerProcessor(+1,+1), surf.cornerPP()), level);
    delete functional;
}


/* *************** BoundedBoxProcessing2D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_L<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_L<T,Descriptor>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_L<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_L<T,Descriptor>* functional,
                               Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_L<T,Descriptor>* functional,
                                   Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}


/* *************** BoundedBoxProcessing2D_S ******************************************* */

template<typename T>
void BoundedBoxProcessingFunctional2D_S<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedBoxProcessingFunctional2D_S<T>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedBoxProcessingFunctional2D_S<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_S<T>* functional,
                               Box2D domain, ScalarFieldBase2D<T>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_S<T>* functional,
                                   Box2D domain, ScalarFieldBase2D<T>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}


/* *************** BoundedBoxProcessing2D_T ******************************************* */

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_T<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_T<T,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge(direction, orientation, domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_T<T,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, domain, dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional2D_T<T,nDim>* functional,
                               Box2D domain, TensorFieldBase2D<T,nDim>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D_T<T,nDim>* functional,
                                   Box2D domain, TensorFieldBase2D<T,nDim>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation2D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}

/* *************** BoundedBoxProcessing2D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional2D_LL<T,Descriptor1,Descriptor2>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor1>&>(*atomicBlocks[0]),
                   dynamic_cast<BlockLattice2D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing2D_SS ****************************************** */

template<typename T>
void BoundedBoxProcessingFunctional2D_SS<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedBoxProcessingFunctional2D_SS<T>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedBoxProcessingFunctional2D_SS<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing2D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional2D_TT<T,nDim1,nDim2>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<TensorField2D<T,nDim1>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing2D_ST ****************************************** */

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_ST<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_ST<T,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional2D_ST<T,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing2D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_LS<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_LS<T,Descriptor>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional2D_LS<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField2D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing2D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processEdgeGeneric (
        int direction, int orientation,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processEdge( direction, orientation, domain,
                 dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional2D_LT<T,Descriptor,nDim>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, domain,
                   dynamic_cast<BlockLattice2D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField2D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedLatticeBoxProcessing2D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>::processEdgeGeneric (
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
void BoundedLatticeBoxProcessingFunctional2D<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<BlockLattice2D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice2D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, domain, lattices);
}

/* *************** BoundedScalarFieldBoxProcessing2D ******************************************* */

template<typename T>
void BoundedScalarFieldBoxProcessingFunctional2D<T>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T>
void BoundedScalarFieldBoxProcessingFunctional2D<T>::processEdgeGeneric (
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
void BoundedScalarFieldBoxProcessingFunctional2D<T>::processCornerGeneric (
        int normalX, int normalY,
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<ScalarField2D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, domain, fields);
}


/* *************** BoundedTensorFieldBoxProcessing2D ******************************************* */

template<typename T, int nDim>
void BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>::processBulkGeneric (
        Box2D domain, std::vector<AtomicBlock2D<T>*> atomicBlocks )
{
    std::vector<TensorField2D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField2D<T,nDim>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T, int nDim>
void BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>::processEdgeGeneric (
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
void BoundedTensorFieldBoxProcessingFunctional2D<T,nDim>::processCornerGeneric (
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

#endif  // DATA_PROCESSOR_WRAPPER_2D_HH
