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
#ifndef DATA_ANALYSIS_FUNCTIONALS_2D_HH
#define DATA_ANALYSIS_FUNCTIONALS_2D_HH

#include "core/dataAnalysisFunctionals2D.h"
#include "core/blockStatistics.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "atomicBlock/dataCouplingWrapper2D.h"
#include "atomicBlock/reductiveDataCouplingWrapper2D.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "core/plbDebug.h"
#include "finiteDifference/fdStencils1D.h"
#include <cmath>

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template<typename T, template<typename U> class Descriptor> 
BoxSumRhoBarFunctional2D<T,Descriptor>::BoxSumRhoBarFunctional2D()
    : sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumRhoBarFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            statistics.gatherSum (
                    sumRhoBarId, cell.getDynamics().computeRhoBar(cell)
            );
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumRhoBarFunctional2D<T,Descriptor>*
    BoxSumRhoBarFunctional2D<T,Descriptor>::clone() const
{
    return new BoxSumRhoBarFunctional2D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumRhoBarFunctional2D<T,Descriptor>::getSumRhoBar() const {
    return this->getStatistics().getSum(sumRhoBarId);
}


template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional2D<T,Descriptor>::BoxSumEnergyFunctional2D()
    : sumEnergyId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumEnergyFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            T uNormSqr = VectorTemplate<T,Descriptor>::normSqr(velocity);
            statistics.gatherSum(sumEnergyId, uNormSqr);
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional2D<T,Descriptor>*
    BoxSumEnergyFunctional2D<T,Descriptor>::clone() const
{
    return new BoxSumEnergyFunctional2D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumEnergyFunctional2D<T,Descriptor>::getSumEnergy() const {
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}


/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
void ExtractLatticeSubDomainFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice1,
                      BlockLattice2D<T,Descriptor>& lattice2 )
{
    Dot2D offset = computeRelativeDisplacement(lattice1, lattice2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice2.attributeDynamics(iX+offset.x,iY+offset.y,
                                       lattice1.get(iX,iY).getDynamics().clone());
            lattice2.get(iX+offset.x,iY+offset.y).
                attributeValues(lattice1.get(iX,iY));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
ExtractLatticeSubDomainFunctional2D<T,Descriptor>* ExtractLatticeSubDomainFunctional2D<T,Descriptor>::clone() const
{
    return new ExtractLatticeSubDomainFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void ExtractLatticeSubDomainFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT ExtractLatticeSubDomainFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            scalarField.get(iX+offset.x,iY+offset.y)
                = lattice.get(iX,iY).computeDensity();
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxDensityFunctional2D<T,Descriptor>* BoxDensityFunctional2D<T,Descriptor>::clone() const
{
    return new BoxDensityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxDensityFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            scalarField.get(iX+offset.x,iY+offset.y)
                = cell.getDynamics().computeRhoBar(cell);
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxRhoBarFunctional2D<T,Descriptor>* BoxRhoBarFunctional2D<T,Descriptor>::clone() const
{
    return new BoxRhoBarFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxRhoBarFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            scalarField.get(iX+offset.x,iY+offset.y)
                = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxKineticEnergyFunctional2D<T,Descriptor>* BoxKineticEnergyFunctional2D<T,Descriptor>::clone() const
{
    return new BoxKineticEnergyFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxKineticEnergyFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            scalarField.get(iX+offset.x,iY+offset.y)
                = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityNormFunctional2D<T,Descriptor>* BoxVelocityNormFunctional2D<T,Descriptor>::clone() const
{
    return new BoxVelocityNormFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityNormFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional2D<T,Descriptor>::BoxVelocityComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            scalarField.get(iX+offset.x,iY+offset.y) = velocity[iComponent];
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional2D<T,Descriptor>* BoxVelocityComponentFunctional2D<T,Descriptor>::clone() const
{
    return new BoxVelocityComponentFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityComponentFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, TensorField2D<T,Descriptor<T>::d>& tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computeVelocity(tensorField.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityFunctional2D<T,Descriptor>* BoxVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new BoxVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxDeviatoricStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq )
{
    Dot2D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computeDeviatoricStress(PiNeq.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxDeviatoricStressFunctional2D<T,Descriptor>* BoxDeviatoricStressFunctional2D<T,Descriptor>::clone() const
{
    return new BoxDeviatoricStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxDeviatoricStressFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxDeviatoricStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        TensorField2D<T,SymmetricTensor<T,Descriptor>::n>& S )
{
    Dot2D offset = computeRelativeDisplacement(lattice, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            Array<T,SymmetricTensor<T,Descriptor>::n>& element = S.get(iX+offset.x,iY+offset.y);
            cell.computeDeviatoricStress(element);
            T omega     = cell.getDynamics().getOmega();
            T rhoBar    = cell.getDynamics().computeRhoBar(cell);
            T prefactor = - omega * Descriptor<T>::invCs2 *
                            Descriptor<T>::invRho(rhoBar) / (T)2;
            for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                element[iTensor] *= prefactor;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxStrainRateFromStressFunctional2D<T,Descriptor>* BoxStrainRateFromStressFunctional2D<T,Descriptor>::clone() const
{
    return new BoxStrainRateFromStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxStrainRateFromStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional2D<T,Descriptor>::BoxPopulationFunctional2D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            scalarField.get(iX+offset.x,iY+offset.y) = lattice.get(iX,iY)[iComponent];
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional2D<T,Descriptor>* BoxPopulationFunctional2D<T,Descriptor>::clone() const
{
    return new BoxPopulationFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxPopulationFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for ScalarField ******** */

template<typename T>
BoxScalarSumFunctional2D<T>::BoxScalarSumFunctional2D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoxScalarSumFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            statistics.gatherSum(sumScalarId, scalarField.get(iX,iY));
        }
    }
}

template<typename T>
BoxScalarSumFunctional2D<T>* BoxScalarSumFunctional2D<T>::clone() const
{
    return new BoxScalarSumFunctional2D<T>(*this);
}

template<typename T>
T BoxScalarSumFunctional2D<T>::getSumScalar() const {
    return this->getStatistics().getSum(sumScalarId);
}


template<typename T>
BoxScalarMinFunctional2D<T>::BoxScalarMinFunctional2D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMinFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // BlockStatistics computes only maximum, no minimum. Therefore,
            //   the relation min(x) = -max(-x) is used.
            statistics.gatherMax(maxScalarId, -scalarField.get(iX,iY));
        }
    }
}

template<typename T>
BoxScalarMinFunctional2D<T>* BoxScalarMinFunctional2D<T>::clone() const
{
    return new BoxScalarMinFunctional2D<T>(*this);
}

template<typename T>
T BoxScalarMinFunctional2D<T>::getMinScalar() const {
    // The minus sign accounts for the relation min(x) = -max(-x).
    return  - this->getStatistics().getMax(maxScalarId);
}


template<typename T>
BoxScalarMaxFunctional2D<T>::BoxScalarMaxFunctional2D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMaxFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            statistics.gatherMax(maxScalarId, scalarField.get(iX,iY));
        }
    }
}

template<typename T>
BoxScalarMaxFunctional2D<T>* BoxScalarMaxFunctional2D<T>::clone() const
{
    return new BoxScalarMaxFunctional2D<T>(*this);
}

template<typename T>
T BoxScalarMaxFunctional2D<T>::getMaxScalar() const {
    return this->getStatistics().getMax(maxScalarId);
}


template<typename T>
BoundedBoxScalarSumFunctional2D<T>::BoundedBoxScalarSumFunctional2D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoundedBoxScalarSumFunctional2D<T>::processBulk (
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            statistics.gatherSum(sumScalarId, scalarField.get(iX,iY));
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional2D<T>::processEdge (
        int direction, int orientation,
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Edge nodes have a weight of 0.5, because only 50% of the
            //   cell centered at the node is inside the computational domain.
            statistics.gatherSum(sumScalarId, scalarField.get(iX,iY) / (T)2);
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional2D<T>::processCorner (
        int normalX, int normalY,
        Box2D domain, ScalarField2D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Corner nodes have a weight of 0.25, because only 25% of the
            //   cell centered at the node is inside the computational domain.
            statistics.gatherSum(sumScalarId, scalarField.get(iX,iY) / (T)4);
        }
    }
}

template<typename T>
BoundedBoxScalarSumFunctional2D<T>* BoundedBoxScalarSumFunctional2D<T>::clone() const
{
    return new BoundedBoxScalarSumFunctional2D<T>(*this);
}

template<typename T>
T BoundedBoxScalarSumFunctional2D<T>::getSumScalar() const {
    return this->getStatistics().getSum(sumScalarId);
}


/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
void ExtractScalarSubDomainFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& field1,
                      ScalarField2D<T>& field2 )
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            field2.get(iX+offset.x,iY+offset.y) = field1.get(iX,iY);
        }
    }
}

template<typename T>
ExtractScalarSubDomainFunctional2D<T>* ExtractScalarSubDomainFunctional2D<T>::clone() const
{
    return new ExtractScalarSubDomainFunctional2D<T>(*this);
}

template<typename T>
void ExtractScalarSubDomainFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT ExtractScalarSubDomainFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_alpha_functional2D ************************************* */

template<typename T>
A_plus_alpha_functional2D<T>::A_plus_alpha_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = A.get(iX,iY) + alpha;
        }
    }
}

template<typename T>
A_plus_alpha_functional2D<T>* A_plus_alpha_functional2D<T>::clone() const {
    return new A_plus_alpha_functional2D<T>(*this);
}

template<typename T>
void A_plus_alpha_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_functional2D ************************************** */

template<typename T>
A_minus_alpha_functional2D<T>::A_minus_alpha_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = A.get(iX,iY) - alpha;
        }
    }
}

template<typename T>
A_minus_alpha_functional2D<T>* A_minus_alpha_functional2D<T>::clone() const {
    return new A_minus_alpha_functional2D<T>(*this);
}

template<typename T>
void A_minus_alpha_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_minus_A_functional2D ************************************* */

template<typename T>
Alpha_minus_A_functional2D<T>::Alpha_minus_A_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_minus_A_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = alpha - A.get(iX,iY);
        }
    }
}

template<typename T>
Alpha_minus_A_functional2D<T>* Alpha_minus_A_functional2D<T>::clone() const {
    return new Alpha_minus_A_functional2D<T>(*this);
}

template<typename T>
void Alpha_minus_A_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT Alpha_minus_A_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_times_alpha_functional2D ************************************* */

template<typename T>
A_times_alpha_functional2D<T>::A_times_alpha_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = A.get(iX,iY) * alpha;
        }
    }
}

template<typename T>
A_times_alpha_functional2D<T>* A_times_alpha_functional2D<T>::clone() const {
    return new A_times_alpha_functional2D<T>(*this);
}

template<typename T>
void A_times_alpha_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_functional2D ************************************* */

template<typename T>
A_dividedBy_alpha_functional2D<T>::A_dividedBy_alpha_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = A.get(iX,iY) / alpha;
        }
    }
}

template<typename T>
A_dividedBy_alpha_functional2D<T>* A_dividedBy_alpha_functional2D<T>::clone() const {
    return new A_dividedBy_alpha_functional2D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_dividedBy_A_functional2D ************************************* */

template<typename T>
Alpha_dividedBy_A_functional2D<T>::Alpha_dividedBy_A_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_dividedBy_A_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& result )
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offset.x,iY+offset.y) = alpha / A.get(iX,iY);
        }
    }
}

template<typename T>
Alpha_dividedBy_A_functional2D<T>* Alpha_dividedBy_A_functional2D<T>::clone() const {
    return new Alpha_dividedBy_A_functional2D<T>(*this);
}

template<typename T>
void Alpha_dividedBy_A_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT Alpha_dividedBy_A_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_plus_alpha_inplace_functional2D ************************************* */

template<typename T>
A_plus_alpha_inplace_functional2D<T>::A_plus_alpha_inplace_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) += alpha;
        }
    }
}

template<typename T>
A_plus_alpha_inplace_functional2D<T>* A_plus_alpha_inplace_functional2D<T>::clone() const {
    return new A_plus_alpha_inplace_functional2D<T>(*this);
}

template<typename T>
void A_plus_alpha_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_inplace_functional2D ************************************** */

template<typename T>
A_minus_alpha_inplace_functional2D<T>::A_minus_alpha_inplace_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) -= alpha;
        }
    }
}

template<typename T>
A_minus_alpha_inplace_functional2D<T>* A_minus_alpha_inplace_functional2D<T>::clone() const {
    return new A_minus_alpha_inplace_functional2D<T>(*this);
}

template<typename T>
void A_minus_alpha_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_alpha_inplace_functional2D ************************************* */

template<typename T>
A_times_alpha_inplace_functional2D<T>::A_times_alpha_inplace_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) *= alpha;
        }
    }
}

template<typename T>
A_times_alpha_inplace_functional2D<T>* A_times_alpha_inplace_functional2D<T>::clone() const {
    return new A_times_alpha_inplace_functional2D<T>(*this);
}

template<typename T>
void A_times_alpha_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_inplace_functional2D ************************************* */

template<typename T>
A_dividedBy_alpha_inplace_functional2D<T>::A_dividedBy_alpha_inplace_functional2D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) /= alpha;
        }
    }
}

template<typename T>
A_dividedBy_alpha_inplace_functional2D<T>* A_dividedBy_alpha_inplace_functional2D<T>::clone() const {
    return new A_dividedBy_alpha_inplace_functional2D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_functional2D ****************************************** */

template<typename T>
void A_plus_B_functional2D<T>::process (
        Box2D domain, std::vector<ScalarField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField2D<T>& A = *fields[0];
    ScalarField2D<T>& B = *fields[1];
    ScalarField2D<T>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) + B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T>
A_plus_B_functional2D<T>* A_plus_B_functional2D<T>::clone() const {
    return new A_plus_B_functional2D<T>(*this);
}

template<typename T>
void A_plus_B_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_B_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_functional2D ****************************************** */

template<typename T>
void A_minus_B_functional2D<T>::process (
        Box2D domain, std::vector<ScalarField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField2D<T>& A = *fields[0];
    ScalarField2D<T>& B = *fields[1];
    ScalarField2D<T>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) - B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T>
A_minus_B_functional2D<T>* A_minus_B_functional2D<T>::clone() const {
    return new A_minus_B_functional2D<T>(*this);
}

template<typename T>
void A_minus_B_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_B_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_functional2D ****************************************** */

template<typename T>
void A_times_B_functional2D<T>::process (
        Box2D domain, std::vector<ScalarField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField2D<T>& A = *fields[0];
    ScalarField2D<T>& B = *fields[1];
    ScalarField2D<T>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) * B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T>
A_times_B_functional2D<T>* A_times_B_functional2D<T>::clone() const {
    return new A_times_B_functional2D<T>(*this);
}

template<typename T>
void A_times_B_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_times_B_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_functional2D ****************************************** */

template<typename T>
void A_dividedBy_B_functional2D<T>::process (
        Box2D domain, std::vector<ScalarField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField2D<T>& A = *fields[0];
    ScalarField2D<T>& B = *fields[1];
    ScalarField2D<T>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) / B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T>
A_dividedBy_B_functional2D<T>* A_dividedBy_B_functional2D<T>::clone() const {
    return new A_dividedBy_B_functional2D<T>(*this);
}

template<typename T>
void A_dividedBy_B_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_inplace_functional2D ****************************************** */

template<typename T>
void A_plus_B_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) += B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T>
A_plus_B_inplace_functional2D<T>* A_plus_B_inplace_functional2D<T>::clone() const {
    return new A_plus_B_inplace_functional2D<T>(*this);
}

template<typename T>
void A_plus_B_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_plus_B_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_inplace_functional2D ****************************************** */

template<typename T>
void A_minus_B_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) -= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T>
A_minus_B_inplace_functional2D<T>* A_minus_B_inplace_functional2D<T>::clone() const {
    return new A_minus_B_inplace_functional2D<T>(*this);
}

template<typename T>
void A_minus_B_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_minus_B_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_inplace_functional2D ****************************************** */

template<typename T>
void A_times_B_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) *= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T>
A_times_B_inplace_functional2D<T>* A_times_B_inplace_functional2D<T>::clone() const {
    return new A_times_B_inplace_functional2D<T>(*this);
}

template<typename T>
void A_times_B_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_times_B_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_inplace_functional2D ****************************************** */

template<typename T>
void A_dividedBy_B_inplace_functional2D<T>::process (
        Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) -= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T>
A_dividedBy_B_inplace_functional2D<T>* A_dividedBy_B_inplace_functional2D<T>::clone() const {
    return new A_dividedBy_B_inplace_functional2D<T>(*this);
}

template<typename T>
void A_dividedBy_B_inplace_functional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_functional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/// Finite Difference operations on data-fields
namespace fdDataField {

template<typename T, int nDim>
inline T bulkXderiv (
        TensorField2D<T,nDim> const& velocity, plint iX, plint iY, int iD )
{
    T dxu = fd::ctl_diff( velocity.get(iX+1,iY)[iD],
                          velocity.get(iX-1,iY)[iD] );
    return dxu;
}

template<typename T, int nDim>
inline T bulkYderiv (
        TensorField2D<T,nDim> const& velocity, plint iX, plint iY, int iD )
{
    T dyu = fd::ctl_diff( velocity.get(iX,iY+1)[iD],
                          velocity.get(iX,iY-1)[iD] );
    return dyu;
}

template<typename T, int nDim>
inline T edgeXderiv (
        TensorField2D<T,nDim> const& velocity, int direction, int orientation,
        plint iX, plint iY, int iD )
{
    if (direction==0) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX              ,iY)[iD],
                             velocity.get(iX-1*orientation,iY)[iD] );
    }
    else {
        return bulkXderiv(velocity, iX,iY, iD);
    }
}

template<typename T, int nDim>
inline T edgeYderiv (
        TensorField2D<T,nDim> const& velocity, int direction, int orientation,
        plint iX, plint iY, int iD )
{
    if (direction==1) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY              )[iD],
                              velocity.get(iX,iY-1*orientation)[iD] );
    }
    else {
        return bulkYderiv(velocity, iX,iY, iD);
    }
}

template<typename T, int nDim>
inline T cornerXderiv (
        TensorField2D<T,nDim> const& velocity,
        int normalX, int normalY,
        plint iX, plint iY, int iD )
{
    int orientation = normalX;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX              ,iY)[iD],
                         velocity.get(iX-1*orientation,iY)[iD] );
}

template<typename T, int nDim>
inline T cornerYderiv (
        TensorField2D<T,nDim> const& velocity,
        int normalX, int normalY,
        plint iX, plint iY, int iD )
{
    int orientation = normalY;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX,iY              )[iD],
                         velocity.get(iX,iY-1*orientation)[iD] );
}

template<typename T, int nDim>
inline T bulkVorticity(TensorField2D<T,nDim> const& velocity, plint iX, plint iY)
{
    T dxuy = bulkXderiv(velocity, iX,iY, 1);
    T dyux = bulkYderiv(velocity, iX,iY, 0);
    return dxuy - dyux;
}

template<typename T, int nDim>
inline T edgeVorticity( TensorField2D<T,nDim> const& velocity, int direction, int orientation,
                        plint iX, plint iY )
{
    T dxuy = edgeXderiv(velocity, direction, orientation, iX,iY, 1);
    T dyux = edgeYderiv(velocity, direction, orientation, iX,iY, 0);
    return dxuy - dyux;
}

template<typename T, int nDim>
inline T cornerVorticity( TensorField2D<T,nDim> const& velocity, int normalX, int normalY,
                          plint iX, plint iY )
{
    T dxuy = cornerXderiv(velocity, normalX, normalY, iX,iY, 1);
    T dyux = cornerYderiv(velocity, normalX, normalY, iX,iY, 0);
    return dxuy - dyux;
}

}  // fdDataField



template<typename T, int nDim>
void ExtractTensorSubDomainFunctional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& field1,
                      TensorField2D<T,nDim>& field2 )
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                field2.get(iX+offset.x,iY+offset.y)[iDim] = field1.get(iX,iY)[iDim];
            }
        }
    }
}

template<typename T, int nDim>
ExtractTensorSubDomainFunctional2D<T,nDim>* ExtractTensorSubDomainFunctional2D<T,nDim>::clone() const
{
    return new ExtractTensorSubDomainFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorSubDomainFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorSubDomainFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
ExtractTensorComponentFunctional2D<T,nDim>::ExtractTensorComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, int nDim>
void ExtractTensorComponentFunctional2D<T,nDim>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            scalarField.get(iX,iY) = tensorField.get(iX+offset.x,iY+offset.y)[iComponent];
        }
    }
}

template<typename T, int nDim>
ExtractTensorComponentFunctional2D<T,nDim>* ExtractTensorComponentFunctional2D<T,nDim>::clone() const
{
    return new ExtractTensorComponentFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorComponentFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorComponentFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormFunctional2D<T,nDim>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            scalarField.get(iX,iY) = std::sqrt( VectorTemplateImpl<T,nDim>::normSqr (
                                                    tensorField.get(iX+offset.x,iY+offset.y) ) );
        }
    }
}

template<typename T, int nDim>
ComputeNormFunctional2D<T,nDim>* ComputeNormFunctional2D<T,nDim>::clone() const
{
    return new ComputeNormFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormSqrFunctional2D<T,nDim>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,nDim>& tensorField )
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            scalarField.get(iX,iY) = VectorTemplateImpl<T,nDim>::normSqr (
                                         tensorField.get(iX+offset.x,iY+offset.y) );
        }
    }
}

template<typename T, int nDim>
ComputeNormSqrFunctional2D<T,nDim>* ComputeNormSqrFunctional2D<T,nDim>::clone() const
{
    return new ComputeNormSqrFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormSqrFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormSqrFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,3>& tensorField )
{
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Array<T,3>& el = tensorField.get(iX+offset.x,iY+offset.y);
                scalarField.get(iX,iY) = std::sqrt ( 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                        // .. and off-diagonal component twice, due to symmetry.
                        (T)2 * util::sqr(el[tensor::xy]) );
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormFunctional2D<T>* ComputeSymmetricTensorNormFunctional2D<T>::clone() const
{
    return new ComputeSymmetricTensorNormFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormSqrFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,3>& tensorField )
{
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,3>& el = tensorField.get(iX+offset.x,iY+offset.y);
            scalarField.get(iX,iY) = 
                    // Count diagonal components once ...
                            util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                    // .. and off-diagonal components twice, due to symmetry.
                    (T)2 * util::sqr(el[tensor::xy]);
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormSqrFunctional2D<T>* ComputeSymmetricTensorNormSqrFunctional2D<T>::clone() const
{
    return new ComputeSymmetricTensorNormSqrFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormSqrFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormSqrFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void ComputeSymmetricTensorTraceFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& scalarField,
                      TensorField2D<T,3>& tensorField )
{
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,3>& el = tensorField.get(iX+offset.x,iY+offset.y);
            scalarField.get(iX,iY) = el[tensor::xx] + el[tensor::yy];
        }
    }
}

template<typename T>
ComputeSymmetricTensorTraceFunctional2D<T>* ComputeSymmetricTensorTraceFunctional2D<T>::clone() const
{
    return new ComputeSymmetricTensorTraceFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorTraceFunctional2D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorTraceFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void BoxBulkVorticityFunctional2D<T,nDim>::process (
        Box2D domain, ScalarField2D<T>& vorticity,
                      TensorField2D<T,nDim>& velocity )
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            vorticity.get(iX,iY) =
                fdDataField::bulkVorticity(velocity, iX2,iY2);
        }
    }
}

template<typename T, int nDim>
BoxBulkVorticityFunctional2D<T,nDim>* BoxBulkVorticityFunctional2D<T,nDim>::clone() const
{
    return new BoxBulkVorticityFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkVorticityFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, int nDim>
void BoxVorticityFunctional2D<T,nDim>::processBulk (
        Box2D domain, ScalarField2D<T>& vorticity,
                      TensorField2D<T,nDim>& velocity )
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            vorticity.get(iX,iY) =
                fdDataField::bulkVorticity(velocity, iX2,iY2);
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional2D<T,nDim>::processEdge (
        int direction, int orientation, Box2D domain,
        ScalarField2D<T>& vorticity, TensorField2D<T,nDim>& velocity )
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            vorticity.get(iX,iY) =
                fdDataField::edgeVorticity(velocity,direction,orientation, iX2,iY2);
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional2D<T,nDim>::processCorner (
        int normalX, int normalY,  Box2D domain,
        ScalarField2D<T>& vorticity, TensorField2D<T,nDim>& velocity )
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            vorticity.get(iX,iY) =
                fdDataField::cornerVorticity(velocity,normalX,normalY, iX2,iY2);
        }
    }
}

template<typename T, int nDim>
BoxVorticityFunctional2D<T,nDim>* BoxVorticityFunctional2D<T,nDim>::clone() const
{
    return new BoxVorticityFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxVorticityFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}


template<typename T, int nDim>
BlockDomain::DomainT BoxVorticityFunctional2D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxBulkStrainRateFunctional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& velocity,
                      TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2);
            el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, 1) +
                               fdDataField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T, int nDim>
BoxBulkStrainRateFunctional2D<T,nDim>* BoxBulkStrainRateFunctional2D<T,nDim>::clone() const
{
    return new BoxBulkStrainRateFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkStrainRateFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkStrainRateFunctional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void BoxStrainRateFunctional2D<T,nDim>::processBulk (
        Box2D domain, TensorField2D<T,nDim>& velocity,
                      TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2);
            el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, 1) +
                               fdDataField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional2D<T,nDim>::processEdge (
        int direction, int orientation, Box2D domain,
        TensorField2D<T,nDim>& velocity, TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2);
            el[tensor::xx] = fdDataField::edgeXderiv(velocity, direction,orientation, iX, iY, 0);
            el[tensor::yy] = fdDataField::edgeYderiv(velocity, direction,orientation, iX, iY, 1);
            el[tensor::xy] = ( fdDataField::edgeXderiv(velocity, direction,orientation, iX, iY, 1) +
                               fdDataField::edgeYderiv(velocity, direction,orientation, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional2D<T,nDim>::processCorner (
        int normalX, int normalY,  Box2D domain,
        TensorField2D<T,nDim>& velocity, TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2);
            el[tensor::xx] = fdDataField::cornerXderiv(velocity, normalX,normalY, iX, iY, 0);
            el[tensor::yy] = fdDataField::cornerYderiv(velocity, normalX,normalY, iX, iY, 1);
            el[tensor::xy] = ( fdDataField::cornerXderiv(velocity, normalX,normalY, iX, iY, 1) +
                               fdDataField::cornerYderiv(velocity, normalX,normalY, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T, int nDim>
BoxStrainRateFunctional2D<T,nDim>* BoxStrainRateFunctional2D<T,nDim>::clone() const
{
    return new BoxStrainRateFunctional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxStrainRateFunctional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}


template<typename T, int nDim>
BlockDomain::DomainT BoxStrainRateFunctional2D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


/* ******** Tensor_A_plus_B_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_plus_B_functional2D<T,nDim>::process (
        Box2D domain, std::vector<TensorField2D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField2D<T,nDim>& A = *fields[0];
    TensorField2D<T,nDim>& B = *fields[1];
    TensorField2D<T,nDim>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) + B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_functional2D<T,nDim>* Tensor_A_plus_B_functional2D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_functional2D<T,nDim>::process (
        Box2D domain, std::vector<TensorField2D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField2D<T,nDim>& A = *fields[0];
    TensorField2D<T,nDim>& B = *fields[1];
    TensorField2D<T,nDim>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) - B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_functional2D<T,nDim>* Tensor_A_minus_B_functional2D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_functional2D<T,nDim>::process (
        Box2D domain, std::vector<TensorField2D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField2D<T,nDim>& A = *fields[0];
    TensorField2D<T,nDim>& B = *fields[1];
    TensorField2D<T,nDim>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) * B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_functional2D<T,nDim>* Tensor_A_times_B_functional2D<T,nDim>::clone() const {
    return new Tensor_A_times_B_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_dividedBy_B_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional2D<T,nDim>::process (
        Box2D domain, std::vector<TensorField2D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField2D<T,nDim>& A = *fields[0];
    TensorField2D<T,nDim>& B = *fields[1];
    TensorField2D<T,nDim>& result = *fields[2];
    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            result.get(iX+offsetResult.x,iY+offsetResult.y)
                = A.get(iX,iY) / B.get(iX+offsetB.x,iY+offsetB.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_functional2D<T,nDim>* Tensor_A_dividedBy_B_functional2D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_plus_B_inplace_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) += B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_inplace_functional2D<T,nDim>* Tensor_A_plus_B_inplace_functional2D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_inplace_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_inplace_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_inplace_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) -= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_inplace_functional2D<T,nDim>* Tensor_A_minus_B_inplace_functional2D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_inplace_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_inplace_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_inplace_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) *= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_inplace_functional2D<T,nDim>* Tensor_A_times_B_inplace_functional2D<T,nDim>::clone() const {
    return new Tensor_A_times_B_inplace_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_inplace_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_dividedBy_B_inplace_functional2D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>::process (
        Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B)
{
    Dot2D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            A.get(iX,iY) -= B.get(iX+offset.x,iY+offset.y);
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>* Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONALS_2D_HH
