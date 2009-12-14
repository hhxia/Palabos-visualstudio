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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H

#include "core/globalDefs.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "boundaryCondition/boundaryDynamics.h"

namespace plb {

/// Advection-diffusion dynamics on flat boundaries
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class AdvectionDiffusionBoundaryDynamics : public StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    AdvectionDiffusionBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object, based on its dynamic type
    virtual AdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
};

/// Advection-diffusion dynamics on flat boundaries
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class RegularizedAdvectionDiffusionBoundaryDynamics : public StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    RegularizedAdvectionDiffusionBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object, based on its dynamic type
    virtual RegularizedAdvectionDiffusionBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
};

/// Advection-diffusion dynamics on 2D corners
template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal>
class AdvectionDiffusionCornerDynamics2D : public StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    AdvectionDiffusionCornerDynamics2D(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>* clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
};


/// Advection-diffusion dynamics on 3D edges
template<typename T, template<typename U> class Descriptor, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgeDynamics3D : public StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    AdvectionDiffusionEdgeDynamics3D(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object, based on its dynamic type
    virtual AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>* clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;

};


/// Advection-diffusion dynamics on 3D corners
template<typename T, template<typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
class AdvectionDiffusionCornerDynamics3D : public StoreDensityDynamics<T,Descriptor>
{
public:
    /// Constructor
    AdvectionDiffusionCornerDynamics3D(Dynamics<T,Descriptor>* baseDynamics);

    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>* clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
};

}  // namespace plb

#endif
