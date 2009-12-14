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

/* The original version of this file was written by Orestis Malaspinas
 * and Andrea Parmigiani.
 */

#ifndef SHAN_CHEN_PROCESSOR_2D_HH
#define SHAN_CHEN_PROCESSOR_2D_HH

#include "multiPhysics/shanChenProcessor2D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference2D.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {


// The two classes and the function which follow implement a template mechanism
//   to compute the force according to the following rule:
//   - If there is a force in the external fields (ExternalField::sizeOfForce=3), return it
//   - If there is no force in the external fields (ExternalField::sizeOfForce=0), return 0

/// Default implementation of ExternalForceAccess2D: return 0
template<typename T, template<typename U> class Descriptor, plint numForceComponents>
struct ExternalForceAccess2D {
    static T getComponent(Cell<T,Descriptor> const& cell, plint iD) {
        return T();
    }
};

/// Specialization of ExternalForceAccess2D: return force from external scalar if available
template<typename T, template<typename U> class Descriptor>
struct ExternalForceAccess2D<T,Descriptor,2> {
    static T getComponent(Cell<T,Descriptor> const& cell, plint iD) {
        return *(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt+iD));
    }
};


/// Automatic instantiation of ExternalForceAccess2D, depending on the Descriptor
template<typename T, template<typename U> class Descriptor>
T getExternalForceComponent2D(Cell<T,Descriptor> const& cell, plint iD) {
    return ExternalForceAccess2D<T, Descriptor, Descriptor<T>::ExternalField::sizeOfForce>::getComponent(cell,iD);
}


/* *************** ShanChenMultiComponentProcessor2D ***************** */

template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor2D <T,Descriptor>::ShanChenMultiComponentProcessor2D(T G_)
    : G(G_)
{ }

template<typename T, template<typename U> class Descriptor>
void ShanChenMultiComponentProcessor2D<T,Descriptor>::process (
        Box2D domain,
        std::vector<BlockLattice2D<T,Descriptor>*> lattices )
{
    // Number of species (or components) which are coupled in this Shan/Chen multi-component fluid.
    plint numSpecies = (plint) lattices.size();
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;
    // Handle to external scalars
    enum {
        densityOffset  = D::ExternalField::densityBeginsAt,
        momentumOffset = D::ExternalField::momentumBeginsAt,
    };
    
    // Compute per-lattice density  and momentum on every site and on each
    //   lattice, and store result in external scalars;  envelope cells are included,
    //   because they are needed to compute the interaction potential in the following.
    //   Note that the per-lattice value of the momentum is stored temporarily only, as
    //   it is corrected later on, based on the common fluid velocity.
    for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
        for (plint iX=domain.x0-1; iX<=domain.x1+1; ++iX) {
            for (plint iY=domain.y0-1; iY<=domain.y1+1; ++iY) {
                // Get "intelligent" value of density through cell object, to account
                //   for the fact that the density value can be user-defined, for example
                //   on boundaries.
                Cell<T,Descriptor>& cell = lattices[iSpecies]->get(iX,iY);
                // And store the result into the corresponding external scalar.
                *cell.getExternal(densityOffset) = cell.computeDensity();
                // Compute momentum through direct access to particle populations, and store
                //   result in corresponding external scalars. Note that Cell::computeVelocity
                //   cannot be used, because it returns the velocity of the external scalars,
                //   not the velocity computed from the particle populations.
                Array<T,Descriptor<T>::d> j;
                momentTemplates<T,Descriptor>::get_j(cell,j);
                for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                    *(cell.getExternal(momentumOffset)+iD) = j[iD];
                }
            }
        }
    }

    // Temporary the relaxation parameters omega at a place where they are less
    //   verbous to access.
    std::vector<T> omega(numSpecies);

    // Compute the interaction force between the species, and store it by
    //   means of a velocity correction in the external velocity field.
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Computation of the common density over all populations, weighted by
            //   the relaxation parameters omega.
            T weightedDensity = T();
            for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                Cell<T,Descriptor> const& cell = lattices[iSpecies]->get(iX,iY);
                // Take this opportunity to store relaxation parameters in vector omega.
                omega[iSpecies] = cell.getDynamics().getOmega();
                weightedDensity += omega[iSpecies] * (*cell.getExternal(densityOffset));
            }
            // Computation of the common velocity, shared among all populations.
            Array<T,Descriptor<T>::d> uTot;
            for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                uTot[iD] = T();
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    T *momentum = lattices[iSpecies]->get(iX,iY).getExternal(momentumOffset);
                    uTot[iD] += momentum[iD] * omega[iSpecies];
                }
                uTot[iD] /= weightedDensity;
            }

            // Computation of the interaction potential.
            // 1. Allocate memory and initialize to zero.
            std::vector<Array<T,D::d> > rhoContribution(numSpecies);
            for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                rhoContribution[iSpecies].resetToZero();
            }
            // 2. Compute the term \sum_i ( t_i rho(x+c_i,t) c_i )
            for (plint iPop = 0; iPop < D::q; ++iPop) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor> const& cell = lattices[iSpecies]->get(nextX,nextY);
                    T rho = *cell.getExternal(densityOffset);
                    for (int iD = 0; iD < D::d; ++iD) {
                        rhoContribution[iSpecies][iD] += D::t[iPop] * rho * D::c[iPop][iD];
                    }
                }
            }

            // Computation and storage of the final velocity, consisting
            //   of uTot plus the momentum difference due to interaction
            //   potential and external force
            for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                Cell<T,Descriptor>& cell = lattices[iSpecies]->get(iX,iY);
                T *momentum = cell.getExternal(momentumOffset);
                for (int iD = 0; iD < D::d; ++iD) {
                    momentum[iD] = uTot[iD];
                    // Initialize force contribution with force from external fields if there
                    //   is any, or with zero otherwise.
                    T forceContribution = getExternalForceComponent2D(cell, iD);
                    // Then, add a contribution from the potential of all other species.
                    for (plint iPartnerSpecies=0; iPartnerSpecies<numSpecies; ++iPartnerSpecies) {
                        if (iPartnerSpecies != iSpecies) {
                            forceContribution -= G * rhoContribution[iPartnerSpecies][iD];
                        }
                    }
                    momentum[iD] += (T)1/omega[iSpecies]*forceContribution;
                    // Multiply by rho to covnert from velocity to momentum.
                    momentum[iD] *= *cell.getExternal(densityOffset);
                }
            }
        }
    }
}


template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor2D<T,Descriptor>*
    ShanChenMultiComponentProcessor2D<T,Descriptor>::clone() const
{
    return new ShanChenMultiComponentProcessor2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ShanChenMultiComponentProcessor2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const
{
    // All blocks are modified by the Shan/Chen processor.
    for (pluint iBlock=0; iBlock<isWritten.size(); ++iBlock) {
        isWritten[iBlock] = true;
    }
}


/* *************** ShanChenSingleComponentProcessor2D ***************** */

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor2D<T,Descriptor>::ShanChenSingleComponentProcessor2D (
        T G_, interparticlePotential::PsiFunction<T>* Psi_ )
    : G(G_),
      Psi(Psi_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor2D<T,Descriptor>::~ShanChenSingleComponentProcessor2D() {
    // Pointer to Psi function is owned; delete it in the destructor.
    delete Psi;
}

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor2D<T,Descriptor>::ShanChenSingleComponentProcessor2D (
        ShanChenSingleComponentProcessor2D<T,Descriptor> const& rhs )
    : G(rhs.G),
      Psi(rhs.Psi->clone())
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor2D<T,Descriptor>&
    ShanChenSingleComponentProcessor2D<T,Descriptor>::operator= (
        ShanChenSingleComponentProcessor2D<T,Descriptor> const& rhs )
{
    G = rhs.G;
    delete Psi; Psi = rhs.Psi->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor2D<T,Descriptor>*
    ShanChenSingleComponentProcessor2D<T,Descriptor>::clone() const
{
    return new ShanChenSingleComponentProcessor2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ShanChenSingleComponentProcessor2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;
    // Handle to external scalars
    enum {
        densityOffset  = D::ExternalField::densityBeginsAt,
        momentumOffset = D::ExternalField::momentumBeginsAt,
    };

    plint nx = domain.getNx() + 2;  // Include a one-cell boundary
    plint ny = domain.getNy() + 2;  // Include a one-cell boundary
    plint offsetX = domain.x0-1;
    plint offsetY = domain.y0-1;
    ScalarField2D<T> psiField(nx,ny);
    
    // Compute density and momentum on every site and store result in external scalars;
    //   furthermore, evaluate the interaction potential Psi and store it into a ScalarField.
    //   Envelope cells are included, because they are needed to compute the interaction potential
    //   in the following. Note that the value of the momentum is stored temporarily only, as
    //   it is corrected later on to include corrections due to the interaction potential.
    for (plint iX=domain.x0-1; iX<=domain.x1+1; ++iX) {
        for (plint iY=domain.y0-1; iY<=domain.y1+1; ++iY) {
            // Get "intelligent" value of density through cell object, to account
            //   for the fact that the density value can be user-defined, for example
            //   on boundaries.
            Cell<T,Descriptor>& cell = lattice.get(iX,iY);
            T rho = cell.computeDensity();
            // Evaluate potential function psi.
            psiField.get(iX-offsetX, iY-offsetY) = Psi->compute(rho);
            // Store density into the corresponding external scalar.
            *cell.getExternal(densityOffset) = rho;
            // Compute momentum through direct access to particle populations, and store
            //   result in corresponding external scalars. Note that Cell::computeVelocity
            //   cannot be used, because it returns the velocity of the external scalars,
            //   not the velocity computed from the particle populations.
            Array<T,Descriptor<T>::d> j;
            momentTemplates<T,Descriptor>::get_j(cell,j);
            for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                *(cell.getExternal(momentumOffset)+iD) = j[iD];
            }
        }
    }

    // Compute the interparticle forces, and store they by means of a 
    //   velocity correction in the external velocity field.
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,D::d> rhoContribution;
            rhoContribution.resetToZero();
            // Compute the term \sum_i ( t_i psi(x+c_i,t) c_i )
            for (plint iPop = 0; iPop < D::q; ++iPop) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                T psi = psiField.get(nextX-offsetX, nextY-offsetY);
                for (int iD = 0; iD < D::d; ++iD) {
                    rhoContribution[iD] += D::t[iPop] * psi * D::c[iPop][iD];
                }
            }

            // Computation and storage of the final momentum, including tho momentum
            //   difference due to interaction potential and the external force.
            Cell<T,Descriptor>& cell = lattice.get(iX,iY);
            T *momentum = cell.getExternal(momentumOffset);
            for (int iD = 0; iD < D::d; ++iD) {
                // Initialize force contribution with force from external fields if there
                //   is any, or with zero otherwise.
                T forceContribution = getExternalForceComponent2D(cell, iD);
                // Add interaction term.
                T psi = psiField.get(iX-offsetX, iY-offsetY);
                forceContribution -= G * psi * rhoContribution[iD];
                // Include into total momentum.
                momentum[iD] += (T)1/cell.getDynamics().getOmega()*forceContribution;
            }
        }
    }
}


}  // namespace plb

#endif  // SHAN_CHEN_PROCESSOR_2D_HH
