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
#ifndef DATA_ANALYSIS_FUNCTIONALS_3D_HH
#define DATA_ANALYSIS_FUNCTIONALS_3D_HH

#include "core/dataAnalysisFunctionals3D.h"
#include "core/blockStatistics.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/dataCouplingWrapper3D.h"
#include "atomicBlock/reductiveDataCouplingWrapper3D.h"
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
BoxSumRhoBarFunctional3D<T,Descriptor>::BoxSumRhoBarFunctional3D()
    : sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumRhoBarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                statistics.gatherSum (
                        sumRhoBarId, cell.getDynamics().computeRhoBar(cell)
                );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumRhoBarFunctional3D<T,Descriptor>*
    BoxSumRhoBarFunctional3D<T,Descriptor>::clone() const
{
    return new BoxSumRhoBarFunctional3D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumRhoBarFunctional3D<T,Descriptor>::getSumRhoBar() const {
    return this->getStatistics().getSum(sumRhoBarId);
}


template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional3D<T,Descriptor>::BoxSumEnergyFunctional3D()
    : sumEnergyId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                T uNormSqr = VectorTemplate<T,Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional3D<T,Descriptor>*
    BoxSumEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new BoxSumEnergyFunctional3D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumEnergyFunctional3D<T,Descriptor>::getSumEnergy() const {
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}


/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
void ExtractLatticeSubDomainFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice1,
                      BlockLattice3D<T,Descriptor>& lattice2 )
{
    Dot3D offset = computeRelativeDisplacement(lattice1, lattice2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice2.attributeDynamics(iX+offset.x,iY+offset.y,iZ+offset.z,
                                           lattice1.get(iX,iY,iZ).getDynamics().clone());
                lattice2.get(iX+offset.x,iY+offset.y,iZ+offset.z).
                    attributeValues(lattice1.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
ExtractLatticeSubDomainFunctional3D<T,Descriptor>* ExtractLatticeSubDomainFunctional3D<T,Descriptor>::clone() const
{
    return new ExtractLatticeSubDomainFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void ExtractLatticeSubDomainFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT ExtractLatticeSubDomainFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = lattice.get(iX,iY,iZ).computeDensity();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxDensityFunctional3D<T,Descriptor>* BoxDensityFunctional3D<T,Descriptor>::clone() const
{
    return new BoxDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxDensityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = cell.getDynamics().computeRhoBar(cell);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxRhoBarFunctional3D<T,Descriptor>* BoxRhoBarFunctional3D<T,Descriptor>::clone() const
{
    return new BoxRhoBarFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxRhoBarFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxKineticEnergyFunctional3D<T,Descriptor>* BoxKineticEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new BoxKineticEnergyFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxKineticEnergyFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityNormFunctional3D<T,Descriptor>* BoxVelocityNormFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityNormFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityNormFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional3D<T,Descriptor>::BoxVelocityComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = velocity[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional3D<T,Descriptor>* BoxVelocityComponentFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityComponentFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityComponentFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::d>& tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeVelocity (
                        tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityFunctional3D<T,Descriptor>* BoxVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}


template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxDeviatoricStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq )
{
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeDeviatoricStress (
                        PiNeq.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxDeviatoricStressFunctional3D<T,Descriptor>* BoxDeviatoricStressFunctional3D<T,Descriptor>::clone() const
{
    return new BoxDeviatoricStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxDeviatoricStressFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxDeviatoricStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S )
{
    Dot3D offset = computeRelativeDisplacement(lattice, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                Array<T,SymmetricTensor<T,Descriptor>::n>& element
                    = S.get(iX+offset.x,iY+offset.y,iZ+offset.z);
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
}

template<typename T, template<typename U> class Descriptor> 
BoxStrainRateFromStressFunctional3D<T,Descriptor>* BoxStrainRateFromStressFunctional3D<T,Descriptor>::clone() const
{
    return new BoxStrainRateFromStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxStrainRateFromStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional3D<T,Descriptor>::BoxPopulationFunctional3D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = lattice.get(iX,iY,iZ)[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional3D<T,Descriptor>* BoxPopulationFunctional3D<T,Descriptor>::clone() const
{
    return new BoxPopulationFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxPopulationFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for scalar-field ******* */

template<typename T>
BoxScalarSumFunctional3D<T>::BoxScalarSumFunctional3D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoxScalarSumFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarSumFunctional3D<T>* BoxScalarSumFunctional3D<T>::clone() const
{
    return new BoxScalarSumFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarSumFunctional3D<T>::getSumScalar() const {
    return this->getStatistics().getSum(sumScalarId);
}


template<typename T>
BoxScalarMinFunctional3D<T>::BoxScalarMinFunctional3D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMinFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                statistics.gatherMax(maxScalarId, -scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarMinFunctional3D<T>* BoxScalarMinFunctional3D<T>::clone() const
{
    return new BoxScalarMinFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarMinFunctional3D<T>::getMinScalar() const {
    // The minus sign accounts for the relation min(x) = -max(-x).
    return  - this->getStatistics().getMax(maxScalarId);
}


template<typename T>
BoxScalarMaxFunctional3D<T>::BoxScalarMaxFunctional3D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMaxFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherMax(maxScalarId, scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarMaxFunctional3D<T>* BoxScalarMaxFunctional3D<T>::clone() const
{
    return new BoxScalarMaxFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarMaxFunctional3D<T>::getMaxScalar() const {
    return this->getStatistics().getMax(maxScalarId);
}


template<typename T>
BoundedBoxScalarSumFunctional3D<T>::BoundedBoxScalarSumFunctional3D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processBulk (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processPlane (
        int direction, int orientation,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Plane boundary nodes have a weight of 0.5, because only 50% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, scalarField.get(iX,iY,iZ) / (T)2);
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Edge nodes have a weight of 0.25, because only 25% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, scalarField.get(iX,iY,iZ) / (T)4);
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics<T>& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Corner nodes have a weight of 0.125, because only 1/8 of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, scalarField.get(iX,iY,iZ) / (T)8);
            }
        }
    }
}

template<typename T>
BoundedBoxScalarSumFunctional3D<T>* BoundedBoxScalarSumFunctional3D<T>::clone() const
{
    return new BoundedBoxScalarSumFunctional3D<T>(*this);
}

template<typename T>
T BoundedBoxScalarSumFunctional3D<T>::getSumScalar() const {
    return this->getStatistics().getSum(sumScalarId);
}


/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
void ExtractScalarSubDomainFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& field1,
                      ScalarField3D<T>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field2.get(iX+offset.x,iY+offset.y,iZ+offset.z) = field1.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
ExtractScalarSubDomainFunctional3D<T>* ExtractScalarSubDomainFunctional3D<T>::clone() const
{
    return new ExtractScalarSubDomainFunctional3D<T>(*this);
}

template<typename T>
void ExtractScalarSubDomainFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT ExtractScalarSubDomainFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_plus_alpha_functional3D ************************************* */

template<typename T>
A_plus_alpha_functional3D<T>::A_plus_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) + alpha;
            }
        }
    }
}

template<typename T>
A_plus_alpha_functional3D<T>* A_plus_alpha_functional3D<T>::clone() const {
    return new A_plus_alpha_functional3D<T>(*this);
}

template<typename T>
void A_plus_alpha_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_functional3D ************************************** */

template<typename T>
A_minus_alpha_functional3D<T>::A_minus_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) - alpha;
            }
        }
    }
}

template<typename T>
A_minus_alpha_functional3D<T>* A_minus_alpha_functional3D<T>::clone() const {
    return new A_minus_alpha_functional3D<T>(*this);
}

template<typename T>
void A_minus_alpha_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_minus_A_functional3D ************************************* */

template<typename T>
Alpha_minus_A_functional3D<T>::Alpha_minus_A_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_minus_A_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = alpha - A.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
Alpha_minus_A_functional3D<T>* Alpha_minus_A_functional3D<T>::clone() const {
    return new Alpha_minus_A_functional3D<T>(*this);
}

template<typename T>
void Alpha_minus_A_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT Alpha_minus_A_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_times_alpha_functional3D ************************************* */

template<typename T>
A_times_alpha_functional3D<T>::A_times_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) * alpha;
            }
        }
    }
}

template<typename T>
A_times_alpha_functional3D<T>* A_times_alpha_functional3D<T>::clone() const {
    return new A_times_alpha_functional3D<T>(*this);
}

template<typename T>
void A_times_alpha_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_functional3D ************************************* */

template<typename T>
A_dividedBy_alpha_functional3D<T>::A_dividedBy_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) / alpha;
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_functional3D<T>* A_dividedBy_alpha_functional3D<T>::clone() const {
    return new A_dividedBy_alpha_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_dividedBy_A_functional3D ************************************* */

template<typename T>
Alpha_dividedBy_A_functional3D<T>::Alpha_dividedBy_A_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_dividedBy_A_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = alpha / A.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
Alpha_dividedBy_A_functional3D<T>* Alpha_dividedBy_A_functional3D<T>::clone() const {
    return new Alpha_dividedBy_A_functional3D<T>(*this);
}

template<typename T>
void Alpha_dividedBy_A_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T>
BlockDomain::DomainT Alpha_dividedBy_A_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_plus_alpha_inplace_functional3D ************************************* */

template<typename T>
A_plus_alpha_inplace_functional3D<T>::A_plus_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += alpha;
            }
        }
    }
}

template<typename T>
A_plus_alpha_inplace_functional3D<T>* A_plus_alpha_inplace_functional3D<T>::clone() const {
    return new A_plus_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_plus_alpha_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_inplace_functional3D ************************************** */

template<typename T>
A_minus_alpha_inplace_functional3D<T>::A_minus_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= alpha;
            }
        }
    }
}

template<typename T>
A_minus_alpha_inplace_functional3D<T>* A_minus_alpha_inplace_functional3D<T>::clone() const {
    return new A_minus_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_minus_alpha_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_alpha_inplace_functional3D ************************************* */

template<typename T>
A_times_alpha_inplace_functional3D<T>::A_times_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= alpha;
            }
        }
    }
}

template<typename T>
A_times_alpha_inplace_functional3D<T>* A_times_alpha_inplace_functional3D<T>::clone() const {
    return new A_times_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_times_alpha_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_inplace_functional3D ************************************* */

template<typename T>
A_dividedBy_alpha_inplace_functional3D<T>::A_dividedBy_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) /= alpha;
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_inplace_functional3D<T>* A_dividedBy_alpha_inplace_functional3D<T>::clone() const {
    return new A_dividedBy_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_functional3D ****************************************** */

template<typename T>
void A_plus_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_plus_B_functional3D<T>* A_plus_B_functional3D<T>::clone() const {
    return new A_plus_B_functional3D<T>(*this);
}

template<typename T>
void A_plus_B_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_plus_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_functional3D ****************************************** */

template<typename T>
void A_minus_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_minus_B_functional3D<T>* A_minus_B_functional3D<T>::clone() const {
    return new A_minus_B_functional3D<T>(*this);
}

template<typename T>
void A_minus_B_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_minus_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_functional3D ****************************************** */

template<typename T>
void A_times_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_times_B_functional3D<T>* A_times_B_functional3D<T>::clone() const {
    return new A_times_B_functional3D<T>(*this);
}

template<typename T>
void A_times_B_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_times_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_functional3D ****************************************** */

template<typename T>
void A_dividedBy_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_dividedBy_B_functional3D<T>* A_dividedBy_B_functional3D<T>::clone() const {
    return new A_dividedBy_B_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_inplace_functional3D ****************************************** */

template<typename T>
void A_plus_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_plus_B_inplace_functional3D<T>* A_plus_B_inplace_functional3D<T>::clone() const {
    return new A_plus_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_plus_B_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_plus_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_inplace_functional3D ****************************************** */

template<typename T>
void A_minus_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_minus_B_inplace_functional3D<T>* A_minus_B_inplace_functional3D<T>::clone() const {
    return new A_minus_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_minus_B_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_minus_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_inplace_functional3D ****************************************** */

template<typename T>
void A_times_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_times_B_inplace_functional3D<T>* A_times_B_inplace_functional3D<T>::clone() const {
    return new A_times_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_times_B_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_times_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_inplace_functional3D ****************************************** */

template<typename T>
void A_dividedBy_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_dividedBy_B_inplace_functional3D<T>* A_dividedBy_B_inplace_functional3D<T>::clone() const {
    return new A_dividedBy_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_inplace_functional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

/// Finite Difference operations on data-fields
namespace fdDataField {

template<typename T, int nDim>
inline T bulkXderiv (
        TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ, int iD )
{
    T dxu = fd::ctl_diff( velocity.get(iX+1,iY,iZ)[iD],
                          velocity.get(iX-1,iY,iZ)[iD] );
    return dxu;
}

template<typename T, int nDim>
inline T bulkYderiv (
        TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ, int iD )
{
    T dyu = fd::ctl_diff( velocity.get(iX,iY+1,iZ)[iD],
                          velocity.get(iX,iY-1,iZ)[iD] );
    return dyu;
}

template<typename T, int nDim>
inline T bulkZderiv (
        TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ, int iD )
{
    T dzu = fd::ctl_diff( velocity.get(iX,iY,iZ+1)[iD],
                          velocity.get(iX,iY,iZ-1)[iD] );
    return dzu;
}

template<typename T, int nDim>
inline T planeXderiv (
        TensorField3D<T,nDim> const& velocity, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==0) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX              ,iY,iZ)[iD],
                             velocity.get(iX-1*orientation,iY,iZ)[iD] );
    }
    else {
        return bulkXderiv(velocity, iX,iY,iZ, iD);
    }
}

template<typename T, int nDim>
inline T planeYderiv (
        TensorField3D<T,nDim> const& velocity, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==1) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY              ,iZ)[iD],
                             velocity.get(iX,iY-1*orientation,iZ)[iD] );
    }
    else {
        return bulkYderiv(velocity, iX,iY,iZ, iD);
    }
}

template<typename T, int nDim>
inline T planeZderiv (
        TensorField3D<T,nDim> const& velocity, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==2) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY,iZ              )[iD],
                             velocity.get(iX,iY,iZ-1*orientation)[iD] );
    }
    else {
        return bulkZderiv(velocity, iX,iY,iZ, iD);
    }
}

template<typename T, int nDim>
inline T edgeXderiv (
        TensorField3D<T,nDim> const& velocity,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==0) {
        return bulkXderiv(velocity, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==1 ? direction2 : direction1;
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX              ,iY,iZ)[iD],
                             velocity.get(iX-1*orientation,iY,iZ)[iD] );
    }
}

template<typename T, int nDim>
inline T edgeYderiv (
        TensorField3D<T,nDim> const& velocity,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==1) {
        return bulkYderiv(velocity, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==0 ? direction1 : direction2;
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY              ,iZ)[iD],
                             velocity.get(iX,iY-1*orientation,iZ)[iD] );
    }
}

template<typename T, int nDim>
inline T edgeZderiv (
        TensorField3D<T,nDim> const& velocity,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==2) {
        return bulkZderiv(velocity, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==0 ? direction2 : direction1;
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY,iZ              )[iD],
                             velocity.get(iX,iY,iZ-1*orientation)[iD] );
    }
}

template<typename T, int nDim>
inline T cornerXderiv (
        TensorField3D<T,nDim> const& velocity,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalX;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX              ,iY,iZ)[iD],
                         velocity.get(iX-1*orientation,iY,iZ)[iD] );
}

template<typename T, int nDim>
inline T cornerYderiv (
        TensorField3D<T,nDim> const& velocity,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalY;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX,iY              ,iZ)[iD],
                         velocity.get(iX,iY-1*orientation,iZ)[iD] );
}

template<typename T, int nDim>
inline T cornerZderiv (
        TensorField3D<T,nDim> const& velocity,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalZ;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX,iY,iZ              )[iD],
                         velocity.get(iX,iY,iZ-1*orientation)[iD] );
}


template<typename T, int nDim>
inline T bulkVorticityX(TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ )
{
    T dyuz = fdDataField::bulkYderiv(velocity, iX,iY,iZ, 2);
    T dzuy = fdDataField::bulkZderiv(velocity, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T, int nDim>
inline T bulkVorticityY(TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ )
{
    T dzux = fdDataField::bulkZderiv(velocity, iX,iY,iZ, 0);
    T dxuz = fdDataField::bulkXderiv(velocity, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T, int nDim>
inline T bulkVorticityZ(TensorField3D<T,nDim> const& velocity, plint iX, plint iY, plint iZ )
{
    T dxuy = fdDataField::bulkXderiv(velocity, iX,iY,iZ, 1);
    T dyux = fdDataField::bulkYderiv(velocity, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T, int nDim>
inline T planeVorticityX( TensorField3D<T,nDim> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dyuz = fdDataField::planeYderiv(velocity,direction,orientation, iX,iY,iZ, 2);
    T dzuy = fdDataField::planeZderiv(velocity,direction,orientation, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T, int nDim>
inline T planeVorticityY( TensorField3D<T,nDim> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dzux = fdDataField::planeZderiv(velocity,direction,orientation, iX,iY,iZ, 0);
    T dxuz = fdDataField::planeXderiv(velocity,direction,orientation, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T, int nDim>
inline T planeVorticityZ( TensorField3D<T,nDim> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dxuy = fdDataField::planeXderiv(velocity,direction,orientation, iX,iY,iZ, 1);
    T dyux = fdDataField::planeYderiv(velocity,direction,orientation, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T, int nDim>
inline T edgeVorticityX( TensorField3D<T,nDim> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dyuz = fdDataField::edgeYderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 2);
    T dzuy = fdDataField::edgeZderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T, int nDim>
inline T edgeVorticityY( TensorField3D<T,nDim> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dzux = fdDataField::edgeZderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 0);
    T dxuz = fdDataField::edgeXderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T, int nDim>
inline T edgeVorticityZ( TensorField3D<T,nDim> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dxuy = fdDataField::edgeXderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 1);
    T dyux = fdDataField::edgeYderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T, int nDim>
inline T cornerVorticityX( TensorField3D<T,nDim> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dyuz = fdDataField::cornerYderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 2);
    T dzuy = fdDataField::cornerZderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T, int nDim>
inline T cornerVorticityY( TensorField3D<T,nDim> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dzux = fdDataField::cornerZderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 0);
    T dxuz = fdDataField::cornerXderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T, int nDim>
inline T cornerVorticityZ( TensorField3D<T,nDim> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dxuy = fdDataField::cornerXderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 1);
    T dyux = fdDataField::cornerYderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 0);

    return dxuy - dyux;
}

}  // fdDataField


template<typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& field1,
                      TensorField3D<T,nDim>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    field2.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                        = field1.get(iX,iY,iZ)[iDim];
                }
            }
        }
    }
}

template<typename T, int nDim>
ExtractTensorSubDomainFunctional3D<T,nDim>* ExtractTensorSubDomainFunctional3D<T,nDim>::clone() const
{
    return new ExtractTensorSubDomainFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorSubDomainFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
ExtractTensorComponentFunctional3D<T,nDim>::ExtractTensorComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, int nDim>
void ExtractTensorComponentFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ)
                    = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iComponent];
            }
        }
    }
}

template<typename T, int nDim>
ExtractTensorComponentFunctional3D<T,nDim>* ExtractTensorComponentFunctional3D<T,nDim>::clone() const
{
    return new ExtractTensorComponentFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorComponentFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorComponentFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ) = std::sqrt( VectorTemplateImpl<T,nDim>::normSqr (
                                                           tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) ) );
            }
        }
    }
}

template<typename T, int nDim>
ComputeNormFunctional3D<T,nDim>* ComputeNormFunctional3D<T,nDim>::clone() const
{
    return new ComputeNormFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormSqrFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ) = VectorTemplateImpl<T,nDim>::normSqr (
                                            tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, int nDim>
ComputeNormSqrFunctional3D<T,nDim>* ComputeNormSqrFunctional3D<T,nDim>::clone() const
{
    return new ComputeNormSqrFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormSqrFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormSqrFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = std::sqrt ( 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal components twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) +util::sqr(el[tensor::yz])) );
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormFunctional3D<T>* ComputeSymmetricTensorNormFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorNormFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal components twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) +util::sqr(el[tensor::yz]));
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormSqrFunctional3D<T>* ComputeSymmetricTensorNormSqrFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = el[tensor::xx] + el[tensor::yy] + el[tensor::zz];
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorTraceFunctional3D<T>* ComputeSymmetricTensorTraceFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorTraceFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorTraceFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, int nDim>
void BoxBulkVorticityFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& velocity,
                      TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkVorticityX(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::bulkVorticityY(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::bulkVorticityZ(velocity, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
BoxBulkVorticityFunctional3D<T,nDim>* BoxBulkVorticityFunctional3D<T,nDim>::clone() const
{
    return new BoxBulkVorticityFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkVorticityFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processBulk (
        Box3D domain, TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkVorticityX(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::bulkVorticityY(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::bulkVorticityZ(velocity, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processPlane (
        int direction, int orientation, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::planeVorticityX(velocity,direction,orientation, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::planeVorticityY(velocity,direction,orientation, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::planeVorticityZ(velocity,direction,orientation, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::edgeVorticityX(velocity,plane,normal1,normal2, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::edgeVorticityY(velocity,plane,normal1,normal2, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::edgeVorticityZ(velocity,plane,normal1,normal2, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::cornerVorticityX(velocity,normalX,normalY,normalZ, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::cornerVorticityY(velocity,normalX,normalY,normalZ, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::cornerVorticityZ(velocity,normalX,normalY,normalZ, iX,iY,iZ);
            }
        }
    }
}


template<typename T, int nDim>
BoxVorticityFunctional3D<T,nDim>* BoxVorticityFunctional3D<T,nDim>::clone() const
{
    return new BoxVorticityFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}


template<typename T, int nDim>
BlockDomain::DomainT BoxVorticityFunctional3D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& velocity,
                      TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
BoxBulkStrainRateFunctional3D<T,nDim>* BoxBulkStrainRateFunctional3D<T,nDim>::clone() const
{
    return new BoxBulkStrainRateFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkStrainRateFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processBulk (
        Box3D domain, TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processPlane (
        int direction, int orientation, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 1) +
                               fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) +
                               fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{

    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) +
                                   fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2);
            }
        }
    }
}


template<typename T, int nDim>
BoxStrainRateFunctional3D<T,nDim>* BoxStrainRateFunctional3D<T,nDim>::clone() const
{
    return new BoxStrainRateFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = true;
}


template<typename T, int nDim>
BlockDomain::DomainT BoxStrainRateFunctional3D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


/* ******** Tensor__B_functional3D ****************************************** */

template<typename T, int nDim>
void Tensor_A_plus_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_functional3D<T,nDim>* Tensor_A_plus_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_functional3D<T,nDim>* Tensor_A_minus_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_functional3D<T,nDim>* Tensor_A_times_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_times_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_dividedBy_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_functional3D<T,nDim>* Tensor_A_dividedBy_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false;
    isWritten[1] = false;
    isWritten[2] = true;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_plus_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_inplace_functional3D<T,nDim>* Tensor_A_plus_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_inplace_functional3D<T,nDim>* Tensor_A_minus_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_inplace_functional3D<T,nDim>* Tensor_A_times_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_times_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_dividedBy_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>* Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[2] = true;
    isWritten[1] = false;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONALS_3D_HH
