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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef VARIABLE_OMEGA_DYNAMICS_H
#define VARIABLE_OMEGA_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/// A dynamics for space-dependent relaxation parameter, generic with respect to base dynamics.
template<typename T, template<typename U> class Descriptor>
class VariableOmegaDynamics : public CompositeDynamics<T,Descriptor> {
public:
    VariableOmegaDynamics(Dynamics<T,Descriptor>* baseDynamics_);
    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual T getOmegaFromCell(Cell<T,Descriptor> const& cell) const =0;
};

/// A dynamics for relaxation parameter dependent on off-equilibrium stress, generic with respect to base dynamics.
template<typename T, template<typename U> class Descriptor>
class OmegaFromPiDynamics : public VariableOmegaDynamics<T,Descriptor> {
public:
    OmegaFromPiDynamics(Dynamics<T,Descriptor>* baseDynamics_);
    virtual T getOmegaFromCell(Cell<T,Descriptor> const& cell) const;
    virtual T getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const =0;
};

} // namespace plb

#endif  // VARIABLE_OMEGA_DYNAMICS_H
