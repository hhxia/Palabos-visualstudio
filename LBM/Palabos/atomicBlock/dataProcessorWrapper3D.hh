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


#ifndef DATA_PROCESSOR_WRAPPER_3D_HH
#define DATA_PROCESSOR_WRAPPER_3D_HH

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class BoxProcessingFunctional3D ************************* */

/** Operation is not applied to envelope by default. **/
template<typename T>
BlockDomain::DomainT BoxProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoxProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void BoxProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}


/* *************** Class BoxProcessor3D ************************************ */

template<typename T>
BoxProcessor3D<T>::BoxProcessor3D(BoxProcessingFunctional3D<T>* functional_,
               Box3D domain_, std::vector<AtomicBlock3D<T>*> atomicBlocks_)
    : functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
BoxProcessor3D<T>::BoxProcessor3D(BoxProcessor3D<T> const& rhs)
    : functional(rhs.functional->clone()),
      domain(rhs.domain), atomicBlocks(rhs.atomicBlocks)
{ }

template<typename T>
BoxProcessor3D<T>& BoxProcessor3D<T>::operator=(BoxProcessor3D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    domain = rhs.domain;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

template<typename T>
BoxProcessor3D<T>::~BoxProcessor3D() {
    delete functional;
}

template<typename T>
Box3D BoxProcessor3D<T>::getDomain() const {
    return domain;
}

template<typename T>
void BoxProcessor3D<T>::process() {
    functional -> processGenericBlocks(domain, atomicBlocks);
}

template<typename T>
BoxProcessor3D<T>* BoxProcessor3D<T>::clone() const {
    return new BoxProcessor3D<T>(*this);
}


/* *************** Class BoxProcessorGenerator3D *************************** */

template<typename T>
BoxProcessorGenerator3D<T>::BoxProcessorGenerator3D(BoxProcessingFunctional3D<T>* functional_, Box3D domain)
    : BoxedDataProcessorGenerator3D<T>(domain),
      functional(functional_)
{ }

template<typename T>
BoxProcessorGenerator3D<T>::~BoxProcessorGenerator3D() {
    delete functional;
}

template<typename T>
BoxProcessorGenerator3D<T>::BoxProcessorGenerator3D(BoxProcessorGenerator3D<T> const& rhs)
    : BoxedDataProcessorGenerator3D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
BoxProcessorGenerator3D<T>& BoxProcessorGenerator3D<T>::operator=(BoxProcessorGenerator3D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT BoxProcessorGenerator3D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void BoxProcessorGenerator3D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void BoxProcessorGenerator3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor3D<T>* BoxProcessorGenerator3D<T>::generate(std::vector<AtomicBlock3D<T>*> atomicBlocks) const {
    return new BoxProcessor3D<T>(functional->clone(), this->getDomain(), atomicBlocks);
}

template<typename T>
BoxProcessorGenerator3D<T>* BoxProcessorGenerator3D<T>::clone() const {
    return new BoxProcessorGenerator3D<T>(*this);
}


/* *************** BoxProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoxProcessingFunctional3D_L<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoxProcessingFunctional3D_L<T,Descriptor>* functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    lattice.executeDataProcessor(BoxProcessorGenerator3D<T>(functional, domain));
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoxProcessingFunctional3D_L<T,Descriptor>* functional,
                                   Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice, plint level)
{
    lattice.addInternalProcessor(BoxProcessorGenerator3D<T>(functional, domain), level);
}


/* *************** BoxProcessing3D_S ******************************************* */

template<typename T>
void BoxProcessingFunctional3D_S<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoxProcessingFunctional3D_S<T>* functional,
                               Box3D domain, ScalarFieldBase3D<T>& field)
{
    field.executeDataProcessor(BoxProcessorGenerator3D<T>(functional, domain));
}

template<typename T>
void integrateProcessingFunctional(BoxProcessingFunctional3D_S<T>* functional,
                                   Box3D domain, ScalarFieldBase3D<T>& field, plint level)
{
    field.addInternalProcessor(BoxProcessorGenerator3D<T>(functional, domain), level);
}


/* *************** BoxProcessing3D_T ******************************************* */

template<typename T, int nDim>
void BoxProcessingFunctional3D_T<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process(domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoxProcessingFunctional3D_T<T,nDim>* functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& field)
{
    field.executeDataProcessor(BoxProcessorGenerator3D<T>(functional, domain));
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoxProcessingFunctional3D_T<T,nDim>* functional,
                                   Box3D domain, TensorFieldBase3D<T,nDim>& field, plint level)
{
    field.addInternalProcessor(BoxProcessorGenerator3D<T>(functional, domain), level);
}


/* *************** BoxProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing3D_SS ****************************************** */

template<typename T>
void BoxProcessingFunctional3D_SS<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void BoxProcessingFunctional3D_ST<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** BoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoxProcessingFunctional3D_LS<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    process ( domain,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** LatticeBoxProcessing3D ************************************** */

template<typename T, template<typename U> class Descriptor>
void LatticeBoxProcessingFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(domain, lattices);
}

/* *************** ScalarFieldBoxProcessing3D *********************************** */

template<typename T>
void ScalarFieldBoxProcessingFunctional3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}

/* *************** TensorFieldBoxProcessing3D *********************************** */

template<typename T, int nDim>
void TensorFieldBoxProcessingFunctional3D<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(domain, fields);
}


/* *************** Class DotProcessingFunctional3D ************************* */

/** Operation is not executed on envelope by default. **/
template<typename T>
BlockDomain::DomainT DotProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void DotProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void DotProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}



/* *************** Class DotProcessor3D ************************************ */

template<typename T>
DotProcessor3D<T>::DotProcessor3D(DotProcessingFunctional3D<T>* functional_,
               DotList3D const& dotList_, std::vector<AtomicBlock3D<T>*> atomicBlocks_)
    : functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

template<typename T>
DotProcessor3D<T>::DotProcessor3D(DotProcessor3D<T> const& rhs)
    : functional(rhs.functional->clone()),
      dotList(rhs.dotList), atomicBlocks(rhs.atomicBlocks)
{ }

template<typename T>
DotProcessor3D<T>& DotProcessor3D<T>::operator=(DotProcessor3D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    dotList = rhs.dotList;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

template<typename T>
DotProcessor3D<T>::~DotProcessor3D() {
    delete functional;
}

template<typename T>
void DotProcessor3D<T>::process() {
    functional -> processGenericBlocks(dotList, atomicBlocks);
}

template<typename T>
DotProcessor3D<T>* DotProcessor3D<T>::clone() const {
    return new DotProcessor3D<T>(*this);
}

template<typename T>
DotList3D const& DotProcessor3D<T>::getDotList() const {
    return dotList;
}


/* *************** Class DotProcessorGenerator3D *************************** */

template<typename T>
DotProcessorGenerator3D<T>::DotProcessorGenerator3D(DotProcessingFunctional3D<T>* functional_, DotList3D const& dotList)
    : DottedDataProcessorGenerator3D<T>(dotList),
      functional(functional_)
{ }

template<typename T>
DotProcessorGenerator3D<T>::~DotProcessorGenerator3D() {
    delete functional;
}

template<typename T>
DotProcessorGenerator3D<T>::DotProcessorGenerator3D(DotProcessorGenerator3D<T> const& rhs)
    : DottedDataProcessorGenerator3D<T>(rhs),
      functional(rhs.functional->clone())
{ }

template<typename T>
DotProcessorGenerator3D<T>& DotProcessorGenerator3D<T>::operator=(DotProcessorGenerator3D<T> const& rhs) {
    delete functional; functional = rhs.functional->clone();
    return *this;
}

template<typename T>
BlockDomain::DomainT DotProcessorGenerator3D<T>::appliesTo() const {
    return functional->appliesTo();
}

template<typename T>
void DotProcessorGenerator3D<T>::rescale(T dxScale, T dtScale) {
    functional->rescale(dxScale, dtScale);
}

template<typename T>
void DotProcessorGenerator3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    functional->getModificationPattern(isWritten);
}

template<typename T>
DataProcessor3D<T>* DotProcessorGenerator3D<T>::generate(std::vector<AtomicBlock3D<T>*> atomicBlocks) const {
    return new DotProcessor3D<T>(functional->clone(), this->getDotList(), atomicBlocks);
}

template<typename T>
DotProcessorGenerator3D<T>* DotProcessorGenerator3D<T>::clone() const {
    return new DotProcessorGenerator3D<T>(*this);
}


/* *************** DotProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void DotProcessingFunctional3D_L<T,Descriptor>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(DotProcessingFunctional3D_L<T,Descriptor>* functional,
                               DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    lattice.executeDataProcessor(DotProcessorGenerator3D<T>(functional, dotList));
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(DotProcessingFunctional3D_L<T,Descriptor>* functional,
                                   DotList3D const& dotList, BlockLatticeBase3D<T,Descriptor>& lattice, plint level)
{
    lattice.addInternalProcessor(DotProcessorGenerator3D<T>(functional, dotList), level);
}


/* *************** DotProcessing3D_S ******************************************* */

template<typename T>
void DotProcessingFunctional3D_S<T>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(DotProcessingFunctional3D_S<T>* functional,
                               DotList3D const& dotList, ScalarFieldBase3D<T>& field)
{
    field.executeDataProcessor(DotProcessorGenerator3D<T>(functional, dotList));
}

template<typename T>
void integrateProcessingFunctional(DotProcessingFunctional3D_S<T>* functional,
                                   DotList3D const& dotList, ScalarFieldBase3D<T>& field, plint level)
{
    field.addInternalProcessor(DotProcessorGenerator3D<T>(functional, dotList), level);
}


/* *************** DotProcessing3D_T ******************************************* */

template<typename T, int nDim>
void DotProcessingFunctional3D_T<T,nDim>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process(dotList, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(DotProcessingFunctional3D_T<T,nDim>* functional,
                               DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field)
{
    field.executeDataProcessor(DotProcessorGenerator3D<T>(functional, dotList));
}

template<typename T, int nDim>
void integrateProcessingFunctional(DotProcessingFunctional3D_T<T,nDim>* functional,
                                   DotList3D const& dotList, TensorFieldBase3D<T,nDim>& field, plint level)
{
    field.addInternalProcessor(DotProcessorGenerator3D<T>(functional, dotList), level);
}


/* *************** DotProcessing3D_LL******************************************* */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void DotProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks ) {
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
              dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing3D_SS ****************************************** */

template<typename T>
void DotProcessingFunctional3D_SS<T>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void DotProcessingFunctional3D_TT<T,nDim1,nDim2>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void DotProcessingFunctional3D_ST<T,nDim>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void DotProcessingFunctional3D_LS<T,Descriptor>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}


/* *************** DotProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void DotProcessingFunctional3D_LT<T,Descriptor,nDim>::processGenericBlocks(DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks) {
    process ( dotList,
              dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
              dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** LatticeDotProcessing3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void LatticeDotProcessingFunctional3D<T,Descriptor>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    process(dotList, lattices);
}

/* *************** ScalarFieldDotProcessing3D ******************************************* */

template<typename T>
void ScalarFieldDotProcessingFunctional3D<T>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}

/* *************** TensorFieldDotProcessing3D ******************************************* */

template<typename T, int nDim>
void TensorFieldDotProcessingFunctional3D<T,nDim>::processGenericBlocks (
        DotList3D const& dotList, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    process(dotList, fields);
}


/* *************** Class BoundedBoxProcessingFunctional3D ************************* */

/** Operation is not applied to envelope by default. **/
template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
template<typename T>
void BoundedBoxProcessingFunctional3D<T>::rescale(T dxScale, T dtScale)
{ }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
template<typename T>
void BoundedBoxProcessingFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}


template<typename T>
BoxProcessingFunctional3D<T>* BoundedBoxProcessingFunctional3D<T>::getBulkProcessor() const {
    return new BulkWrapperFunctional(this->clone());
}

template<typename T>
BoxProcessingFunctional3D<T>* BoundedBoxProcessingFunctional3D<T>::getPlaneProcessor (
                                  int direction, int orientation ) const
{
    return new PlaneWrapperFunctional(this->clone(), direction, orientation);
}

template<typename T>
BoxProcessingFunctional3D<T>* BoundedBoxProcessingFunctional3D<T>::getEdgeProcessor (
                                  int plane, int normal1, int normal2 ) const
{
    return new EdgeWrapperFunctional(this->clone(), plane, normal1, normal2);
}

template<typename T>
BoxProcessingFunctional3D<T>* BoundedBoxProcessingFunctional3D<T>::getCornerProcessor (
                                  int normalX, int normalY, int normalZ ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY, normalZ);
}

/* *************** Class BoundedBoxProcessingFunctional3D::BulkWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional_ )
    : boundedFunctional(boundedFunctional_)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::BulkWrapperFunctional (
        BulkWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone())
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::~BulkWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional&
    BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::operator= (
            BulkWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional*
    BoundedBoxProcessingFunctional3D<T>::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::PlaneWrapperFunctional (
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
        int direction_, int orientation_)
    : boundedFunctional(boundedFunctional_),
      direction(direction_),
      orientation(orientation_)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::PlaneWrapperFunctional (
        PlaneWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      direction(rhs.direction),
      orientation(rhs.orientation)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::~PlaneWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional&
    BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::operator= (
            PlaneWrapperFunctional const& rhs )
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processPlaneGeneric(direction, orientation, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}


template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional*
    BoundedBoxProcessingFunctional3D<T>::PlaneWrapperFunctional::clone() const
{
    return new PlaneWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
        int plane_, int normal1_, int normal2_ )
    : boundedFunctional(boundedFunctional_),
      plane(plane_),
      normal1(normal1_),
      normal2(normal2_)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::EdgeWrapperFunctional (
        EdgeWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      plane(rhs.plane),
      normal1(rhs.normal1),
      normal2(rhs.normal2)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::~EdgeWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional&
    BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::operator= (
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
void BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processEdgeGeneric(plane, normal1, normal2, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}


template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional*
    BoundedBoxProcessingFunctional3D<T>::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::CornerWrapperFunctional ** */

template<typename T>
BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        BoundedBoxProcessingFunctional3D<T>* boundedFunctional_,
        int normalX_,
        int normalY_,
        int normalZ_ )
    : boundedFunctional(boundedFunctional_),
      normalX(normalX_),
      normalY(normalY_),
      normalZ(normalZ_)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::CornerWrapperFunctional (
        CornerWrapperFunctional const& rhs )
    : boundedFunctional(rhs.boundedFunctional->clone()),
      normalX(rhs.normalX),
      normalY(rhs.normalY),
      normalZ(rhs.normalZ)
{ }

template<typename T>
BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::~CornerWrapperFunctional() {
    delete boundedFunctional;
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional&
    BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::operator= (
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
void BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    boundedFunctional->processCornerGeneric(normalX, normalY, normalZ, domain, atomicBlocks);
}

template<typename T>
BlockDomain::DomainT BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::rescale (
        T dxScale, T dtScale )
{
    boundedFunctional->rescale(dxScale, dtScale);
}

template<typename T>
void BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

template<typename T>
typename BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional*
    BoundedBoxProcessingFunctional3D<T>::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}


/* *************** Class BoundedOneBlockProcessingFunctionalOperation3D ******** */

template<typename T>
BoundedOneBlockProcessingFunctionalOperation3D<T>::BoundedOneBlockProcessingFunctionalOperation3D (
        Box3D const& domain, plint boundaryWidth_ )
    : surf(domain, boundaryWidth_)
{ }

template<typename T>
void BoundedOneBlockProcessingFunctionalOperation3D<T>::apply (
        BoundedBoxProcessingFunctional3D<T>* functional, Block3D<T>& block )
{
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()));

    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()));

    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()));

    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1, -1), surf.cornerNNN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1, -1,  1), surf.cornerNNP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1, -1), surf.cornerNPN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor(-1,  1,  1), surf.cornerNPP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1, -1), surf.cornerPNN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1, -1,  1), surf.cornerPNP()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1, -1), surf.cornerPPN()));
    block.executeDataProcessor(BoxProcessorGenerator3D<T>(functional->getCornerProcessor( 1,  1,  1), surf.cornerPPP()));

    delete functional;
}

template<typename T>
void BoundedOneBlockProcessingFunctionalOperation3D<T>::integrate (
        BoundedBoxProcessingFunctional3D<T>* functional, Block3D<T>& block, plint level )
{
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getBulkProcessor(), surf.bulk()), level);

    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,-1), surf.surface0N()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(0,+1), surf.surface0P()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,-1), surf.surface1N()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(1,+1), surf.surface1P()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,-1), surf.surface2N()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getPlaneProcessor(2,+1), surf.surface2P()), level);

    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1, -1), surf.edge0NN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0, -1,  1), surf.edge0NP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1, -1), surf.edge0PN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(0,  1,  1), surf.edge0PP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1, -1), surf.edge1NN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1, -1,  1), surf.edge1NP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1, -1), surf.edge1PN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(1,  1,  1), surf.edge1PP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1, -1), surf.edge2NN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2, -1,  1), surf.edge2NP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1, -1), surf.edge2PN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(2,  1,  1), surf.edge2PP()), level);

    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1, -1, -1), surf.cornerNNN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1, -1,  1), surf.cornerNNP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1,  1, -1), surf.cornerNPN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor(-1,  1,  1), surf.cornerNPP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1, -1, -1), surf.cornerPNN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1, -1,  1), surf.cornerPNP()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1,  1, -1), surf.cornerPPN()), level);
    block.addInternalProcessor(BoxProcessorGenerator3D<T>(functional->getEdgeProcessor( 1,  1,  1), surf.cornerPPP()), level);

    delete functional;
}


/* *************** BoundedBoxProcessing3D_L ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_L<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_L<T,Descriptor>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_L<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane,normal1,normal2, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_L<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]));
}

template<typename T, template<typename U> class Descriptor>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_L<T,Descriptor>* functional,
                               Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T, template<typename U> class Descriptor>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_L<T,Descriptor>* functional,
                                   Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}


/* *************** BoundedBoxProcessing3D_S ******************************************* */

template<typename T>
void BoundedBoxProcessingFunctional3D_S<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedBoxProcessingFunctional3D_S<T>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedBoxProcessingFunctional3D_S<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane,normal1,normal2, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void BoundedBoxProcessingFunctional3D_S<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]));
}

template<typename T>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_S<T>* functional,
                               Box3D domain, ScalarFieldBase3D<T>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_S<T>* functional,
                                   Box3D domain, ScalarFieldBase3D<T>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}


/* *************** BoundedBoxProcessing3D_T ******************************************* */

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_T<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk(domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_T<T,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane(direction, orientation, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_T<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge(plane,normal1,normal2, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_T<T,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner(normalX, normalY, normalZ, domain, dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[0]));
}

template<typename T, int nDim>
void applyProcessingFunctional(BoundedBoxProcessingFunctional3D_T<T,nDim>* functional,
                               Box3D domain, TensorFieldBase3D<T,nDim>& lattice, plint boundaryWidth)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).apply(functional, lattice);
}

template<typename T, int nDim>
void integrateProcessingFunctional(BoundedBoxProcessingFunctional3D_T<T,nDim>* functional,
                                   Box3D domain, TensorFieldBase3D<T,nDim>& lattice,
                                   plint boundaryWidth, plint level)
{
    BoundedOneBlockProcessingFunctionalOperation3D<T>(domain, boundaryWidth).integrate(functional, lattice, level);
}

/* *************** BoundedBoxProcessing3D_LL ****************************************** */

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                 dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U1> class Descriptor1,
                     template<typename U2> class Descriptor2>
void BoundedBoxProcessingFunctional3D_LL<T,Descriptor1,Descriptor2>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    PLB_PRECONDITION( atomicBlocks.size() == 2 );
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor1>&>(*atomicBlocks[0]),
                   dynamic_cast<BlockLattice3D<T,Descriptor2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing3D_SS ****************************************** */

template<typename T>
void BoundedBoxProcessingFunctional3D_SS<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedBoxProcessingFunctional3D_SS<T>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedBoxProcessingFunctional3D_SS<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T>
void BoundedBoxProcessingFunctional3D_SS<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing3D_TT ****************************************** */

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim1, int nDim2>
void BoundedBoxProcessingFunctional3D_TT<T,nDim1,nDim2>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<TensorField3D<T,nDim1>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim2>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing3D_ST ****************************************** */

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_ST<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_ST<T,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_ST<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, int nDim>
void BoundedBoxProcessingFunctional3D_ST<T,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing3D_LS ****************************************** */

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_LS<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_LS<T,Descriptor>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_LS<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor>
void BoundedBoxProcessingFunctional3D_LS<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<ScalarField3D<T>&>(*atomicBlocks[1]) );
}

/* *************** BoundedBoxProcessing3D_LT ****************************************** */

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processBulk( domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processPlaneGeneric (
        int direction, int orientation,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processPlane( direction, orientation, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processEdge( plane,normal1,normal2, domain,
                 dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                 dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void BoundedBoxProcessingFunctional3D_LT<T,Descriptor,nDim>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    processCorner( normalX, normalY, normalZ, domain,
                   dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]),
                   dynamic_cast<TensorField3D<T,nDim>&>(*atomicBlocks[1]) );
}

/* *************** BoundedLatticeBoxProcessing3D ******************************************* */

template<typename T, template<typename U> class Descriptor>
void BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processBulk(domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>::processPlaneGeneric (
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
void BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processEdge(plane,normal1,normal2, domain, lattices);
}

template<typename T, template<typename U> class Descriptor>
void BoundedLatticeBoxProcessingFunctional3D<T,Descriptor>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<BlockLattice3D<T,Descriptor>*> lattices(atomicBlocks.size());
    for (pluint iLattice=0; iLattice<atomicBlocks.size(); ++iLattice) {
        lattices[iLattice] = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[iLattice]);
    }
    processCorner(normalX, normalY, normalZ, domain, lattices);
}

/* *************** BoundedScalarFieldBoxProcessing3D ******************************************* */

template<typename T>
void BoundedScalarFieldBoxProcessingFunctional3D<T>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T>
void BoundedScalarFieldBoxProcessingFunctional3D<T>::processPlaneGeneric (
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
void BoundedScalarFieldBoxProcessingFunctional3D<T>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processEdge(plane,normal1,normal2, domain, fields);
}

template<typename T>
void BoundedScalarFieldBoxProcessingFunctional3D<T>::processCornerGeneric (
        int normalX, int normalY, int normalZ,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<ScalarField3D<T>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[iField]);
    }
    processCorner(normalX, normalY, normalZ, domain, fields);
}


/* *************** BoundedTensorFieldBoxProcessing3D ******************************************* */

template<typename T, int nDim>
void BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>::processBulkGeneric (
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processBulk(domain, fields);
}

template<typename T, int nDim>
void BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>::processPlaneGeneric (
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
void BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>::processEdgeGeneric (
        int plane, int normal1, int normal2,
        Box3D domain, std::vector<AtomicBlock3D<T>*> atomicBlocks )
{
    std::vector<TensorField3D<T,nDim>*> fields(atomicBlocks.size());
    for (pluint iField=0; iField<atomicBlocks.size(); ++iField) {
        fields[iField] = dynamic_cast<TensorField3D<T,nDim>*>(atomicBlocks[iField]);
    }
    processEdge(plane,normal1,normal2, domain, fields);
}

template<typename T, int nDim>
void BoundedTensorFieldBoxProcessingFunctional3D<T,nDim>::processCornerGeneric (
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

#endif  // DATA_PROCESSOR_WRAPPER_3D_HH
