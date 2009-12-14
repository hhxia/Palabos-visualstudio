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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef DATA_ANALYSIS_3D_HH
#define DATA_ANALYSIS_3D_HH

#include "core/dataAnalysis3D.h"
#include "core/dataAnalysisFunctionals3D.h"
#include "core/blockStatistics.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/dataCouplingWrapper3D.h"
#include "atomicBlock/reductiveDataCouplingWrapper3D.h"

namespace plb {

/* *************** Data-analysis functions for ScalarField *********** */

template<typename T>
T computeAverage(ScalarFieldBase3D<T>& scalarField, Box3D domain) {
    BoxScalarSumFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getSumScalar() / (T) domain.nCells();
}

template<typename T>
T computeAverage(ScalarFieldBase3D<T>& scalarField) {
    return computeAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeMin(ScalarFieldBase3D<T>& scalarField, Box3D domain) {
    BoxScalarMinFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMinScalar();
}

template<typename T>
T computeMin(ScalarFieldBase3D<T>& scalarField) {
    return computeMin(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeMax(ScalarFieldBase3D<T>& scalarField, Box3D domain) {
    BoxScalarMaxFunctional3D<T> functional;
    applyProcessingFunctional(functional, domain, scalarField);
    return functional.getMaxScalar();
}

template<typename T>
T computeMax(ScalarFieldBase3D<T>& scalarField) {
    return computeMax(scalarField, scalarField.getBoundingBox());
}


template<typename T>
T computeBoundedAverage(ScalarFieldBase3D<T>& scalarField, Box3D domain) {
    BoundedBoxScalarSumFunctional3D<T> functional;
    plint envelopeWidth=1;
    applyProcessingFunctional(functional, domain, scalarField, envelopeWidth);
    return functional.getSumScalar() /
             (T) ( (domain.getNx()-1)*(domain.getNy()-1) );
}

template<typename T>
T computeBoundedAverage(ScalarFieldBase3D<T>& scalarField) {
    return computeBoundedAverage(scalarField, scalarField.getBoundingBox());
}


template<typename T, class BoolMask> 
plint count(ScalarFieldBase3D<T>& field, Box3D domain, BoolMask boolMask)
{
    CountScalarElementsFunctional3D<T,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, class BoolMask> 
plint count(ScalarFieldBase3D<T>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}


/* *************** Data-analysis functions for TensorField ********** */

template<typename T, plint nDim, class BoolMask> 
plint count(TensorFieldBase3D<T,nDim>& field, Box3D domain, BoolMask boolMask)
{
    CountTensorElementsFunctional3D<T,nDim,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

template<typename T, plint nDim, class BoolMask> 
plint count(TensorFieldBase3D<T,nDim>& field, BoolMask boolMask)
{
    return count(field, field.getBoundingBox(), boolMask);
}

/* *************** Data-analysis functions for BlockLattice ********** */

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return Descriptor<T>::fullRho( functional.getSumRhoBar() / (T) domain.nCells() );
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase3D<T,Descriptor>& lattice) {
    return computeAverageDensity(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain) {
    BoxSumRhoBarFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumRhoBar() / (T) domain.nCells();
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase3D<T,Descriptor>& lattice) {
    return computeAverageRhoBar(lattice, lattice.getBoundingBox());
}


template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain) 
{
    BoxSumEnergyFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSumEnergy() / (T) domain.nCells();;
}

template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase3D<T,Descriptor>& lattice) {
    return computeAverageEnergy(lattice, lattice.getBoundingBox());
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, BoolMask boolMask)
{
    CountLatticeElementsFunctional3D<T,Descriptor,BoolMask> functional(boolMask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase3D<T,Descriptor>& lattice, BoolMask boolMask)
{
    return count(lattice, lattice.getBoundingBox(), boolMask);
}


}  // namespace plb

#endif  // DATA_ANALYSIS_3D_HH
