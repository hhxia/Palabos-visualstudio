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

/* Main author: Orestis Malaspinas
 */

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_HH
#define ADVECTION_DIFFUSION_BOUNDARIES_HH

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "core/util.h"
#include "complexDynamics/utilAdvectionDiffusion.h"
#include "latticeBoltzmann/advectionDiffusionLattices.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void DensityClosure ( Cell<T,Descriptor>& cell, Dynamics<T,Descriptor> const& dynamics )
{
    typedef Descriptor<T> D;
    
    T rho    = dynamics.computeDensity(cell);
    T rhoBar = D::rhoBar(rho);
    
    plint missingNormal = 0;
    std::vector<int> missingDiagonal = indexTemplates::subIndexOutgoing<D,direction,orientation>();
    std::vector<int> knownIndexes   = indexTemplates::remainingIndexes<D>(missingDiagonal);
   // here I know all missing and non missing f_i
    for (pluint iPop = 0; iPop < missingDiagonal.size(); ++iPop)
    {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < D::d; ++iDim)
            numOfNonNullComp += abs(D::c[missingDiagonal[iPop]][iDim]);

        if (numOfNonNullComp == 1)
        {
            missingNormal = missingDiagonal[iPop];
            missingDiagonal.erase(missingDiagonal.begin()+iPop);
            break;
        }
    }
    
    T sum = T();
    for (pluint iPop = 0; iPop < knownIndexes.size(); ++iPop)
    {
        sum += cell[knownIndexes[iPop]];
    }
    cell[missingNormal] = rhoBar - sum;
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void RegularizedClosure ( Cell<T,Descriptor>& cell, Dynamics<T,Descriptor> const& dynamics )
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T,Descriptor> adTempl;
    
    T rhoBar = dynamics.computeRhoBar(cell);
    T* u = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    Array<T,D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD)
    {
        jEq[iD] = Descriptor<T>::fullRho(rhoBar) * u[iD];
    }
    plint missingNormal = 0;
    std::vector<int> missingDiagonal = indexTemplates::subIndexOutgoing<D,direction,orientation>();
   // here I know all missing and non missing f_i
    for (pluint iPop = 0; iPop < missingDiagonal.size(); ++iPop)
    {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < D::d; ++iDim)
            numOfNonNullComp += abs(D::c[missingDiagonal[iPop]][iDim]);

        if (numOfNonNullComp == 1)
        {
            missingNormal = missingDiagonal[iPop];
            missingDiagonal.erase(missingDiagonal.begin()+iPop);
            break;
        }
    }
    
    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices 
    
    cell[missingNormal] = 
        adTempl::bgk_ma1_equilibrium(missingNormal, rhoBar, jEq) 
        -(cell[indexTemplates::opposite<D>(missingNormal)]
        - adTempl::bgk_ma1_equilibrium(
                indexTemplates::opposite<D>(missingNormal), rhoBar, jEq) ) ;
}

// ============= flat wall standard boundary ==================//

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
AdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>::AdvectionDiffusionBoundaryDynamics (
        Dynamics<T,Descriptor>* baseDynamics )
    : StoreDensityDynamics<T,Descriptor>(baseDynamics)
{ }

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
AdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>*
    AdvectionDiffusionBoundaryDynamics<T,Descriptor, direction, orientation>::clone() const
{
    return new AdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void AdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>::completePopulations(Cell<T,Descriptor>& cell) const
{
    DensityClosure<T,Descriptor,direction,orientation>(cell, *this);
}

// ============= flat wall regularized boundary ==================//

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>::
        RegularizedAdvectionDiffusionBoundaryDynamics (Dynamics<T,Descriptor>* baseDynamics )
    : StoreDensityDynamics<T,Descriptor>(baseDynamics)
{ }

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>*
    RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor, direction, orientation>::clone() const
{
    return new RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>(*this);
}

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
void RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>::
        completePopulations(Cell<T,Descriptor>& cell) const
{
    RegularizedClosure<T,Descriptor,direction,orientation>(cell, *this);
}

// =============== 2D corners ===================//


template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T,Descriptor,xNormal,yNormal>::AdvectionDiffusionCornerDynamics2D (
        Dynamics<T,Descriptor>* baseDynamics )
    : StoreDensityDynamics<T,Descriptor>(baseDynamics)
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T,Descriptor,xNormal,yNormal>*
    AdvectionDiffusionCornerDynamics2D<T,Descriptor, xNormal,yNormal>::clone() const
{
    return new AdvectionDiffusionCornerDynamics2D<T,Descriptor,xNormal,yNormal>(*this);
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,Descriptor,xNormal,yNormal>::completePopulations(Cell<T,Descriptor>& cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T,Descriptor> adTempl;
    
    T rhoBar = this->computeRhoBar(cell);
    T* u = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    Array<T,D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD)
    {
        jEq[iD] = Descriptor<T>::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = utilAdvDiff::subIndexOutgoing2DonCorners<D,xNormal,yNormal>();
    // here I know all missing and non missing f_i
    
    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices 
    
    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
    {
        cell[unknownIndexes[iPop]] = 
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq) 
                -(cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
                - adTempl::bgk_ma1_equilibrium(
                        indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq) ) ;
    }
}

// =============== 3D corners ===================//

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>::AdvectionDiffusionCornerDynamics3D (
        Dynamics<T,Descriptor>* baseDynamics )
    : StoreDensityDynamics<T,Descriptor>(baseDynamics)
{ }

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>*
    AdvectionDiffusionCornerDynamics3D<T,Descriptor, xNormal,yNormal,zNormal>::clone() const
{
    return new AdvectionDiffusionCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(*this);
}

template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>::completePopulations(Cell<T,Descriptor>& cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T,Descriptor> adTempl;
    
    T rhoBar = this->computeRhoBar(cell);
    T* u = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    Array<T,D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD)
    {
        jEq[iD] = Descriptor<T>::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = utilAdvDiff::subIndexOutgoing3DonCorners<D,xNormal,yNormal,zNormal>();
    // here I know all missing and non missing f_i
    
    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices 
    
    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
    {
        cell[unknownIndexes[iPop]] = 
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq) 
                -(cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
                - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq) ) ;
    }
}

// =============== 3D edges ===================//

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>::AdvectionDiffusionEdgeDynamics3D (
        Dynamics<T,Descriptor>* baseDynamics )
    : StoreDensityDynamics<T,Descriptor>(baseDynamics)
{ }

template<typename T, template<typename U> class Descriptor,  int plane, int normal1, int normal2>
AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>*
    AdvectionDiffusionEdgeDynamics3D<T,Descriptor, plane,normal1, normal2>::clone() const
{
    return new AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(*this);
}

template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
void AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1, normal2>::completePopulations(Cell<T,Descriptor>& cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T,Descriptor> adTempl;
    
    T rhoBar = this->computeRhoBar(cell);
    T* u = cell.getExternal(Descriptor<T>::ExternalField::velocityBeginsAt);
    Array<T,D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD)
    {
        jEq[iD] = D::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = utilAdvDiff::subIndexOutgoing3DonEdges<D,plane,normal1, normal2>();
    // here I know all missing and non missing f_i
    
    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices 
    
    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop)
    {
        cell[unknownIndexes[iPop]] = 
                adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq) 
                -(cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
                - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq) ) ;
    }
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_BOUNDARIES_HH
