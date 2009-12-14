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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_HH
#define REGULARIZED_BOUNDARY_DYNAMICS_HH

#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "core/cell.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
struct boundaryTemplates {
    static void compute_PiNeq (
         Dynamics<T,Descriptor> const& dynamics,
         Cell<T,Descriptor> const& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr,
         Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq )
    {
        typedef Descriptor<T> L;

        std::vector<int> const& onWallIndices = indexTemplates::subIndex<L, direction, 0>();
        std::vector<int> const& normalIndices = indexTemplates::subIndex<L, direction, orientation>();

        // Compute off-equilibrium for known particle populations.
        Array<T,Descriptor<T>::q> fNeq;
        for (pluint fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
            plint iPop = onWallIndices[fIndex];
            fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
        }
        for (pluint fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
            plint iPop = normalIndices[fIndex];
            if (iPop == 0) {
                fNeq[iPop] = T();  // fNeq[0] will not be used anyway
            }
            else {
                fNeq[iPop] = cell[iPop] - dynamics.computeEquilibrium(iPop, rhoBar, j, jSqr);
            }
        }

        // Compute PiNeq from fNeq, by using "bounce-back of off-equilibrium part" rule.
        int iPi = 0;
        for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
            for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
                PiNeq[iPi] = T();
                for (pluint fIndex=0; fIndex<onWallIndices.size(); ++fIndex)
                {
                    const plint iPop = onWallIndices[fIndex];
                    PiNeq[iPi] += L::c[iPop][iAlpha]*L::c[iPop][iBeta]*fNeq[iPop];
                }
                for (pluint fIndex=0; fIndex<normalIndices.size(); ++fIndex)
                {
                    const plint iPop = normalIndices[fIndex];
                    PiNeq[iPi] += (T)2 * L::c[iPop][iAlpha]*L::c[iPop][iBeta]* fNeq[iPop];
                }
                ++iPi;
            }
        }
    }

};  // struct boundaryTemplates



/* *************** Class RegularizedVelocityBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::
    RegularizedVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>*
    RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>::
    completePopulations(Cell<T,Descriptor>& cell) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    this -> computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    boundaryTemplates<T,Descriptor,direction,orientation>::compute_PiNeq (
            this->getBaseDynamics(), cell, rhoBar, j, jSqr, PiNeq );

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}


/* *************** Class RegularizedDensityBoundaryDynamics ************* */

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>::
    RegularizedDensityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_)
        : DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>*
    RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>::clone() const
{
    return new RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
void RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>::
    completePopulations(Cell<T,Descriptor>& cell) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    this -> computeRhoBarJ(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
    boundaryTemplates<T,Descriptor,direction,orientation>::compute_PiNeq (
            this->getBaseDynamics(), cell, rhoBar, j, jSqr, PiNeq );

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

}  // namespace plb

#endif
