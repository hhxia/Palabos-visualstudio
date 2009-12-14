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
 * Definition of a LB cell -- generic implementation.
 */
#ifndef CELL_HH
#define CELL_HH

#include <algorithm>
#include "core/cell.h"
#include "core/util.h"

namespace plb {

////////////////////////// Class Cell /////////////////////////////

/** The possibility to default construct Cell objects facilitates
 * their use in various types of containers. However, they can not
 * be used directly after construction; the method attributeDynamics()
 * must be used first.
 */
template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>::Cell()
    : takesStat(true), dynamics(0)
{
    iniPop();
    iniExternal();
}

/** This constructor initializes the dynamics, but not the values
 * of the distribution functions. Remember that the dynamics is not
 * owned by the Cell object, the user must ensure its proper
 * destruction and a sufficient life time.
 */
template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>::Cell(Dynamics<T,Descriptor>* dynamics_)
    : takesStat(true), dynamics(dynamics_)
{
    iniPop();
    iniExternal();
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::attributeDynamics(Dynamics<T,Descriptor>* dynamics_) {
    dynamics = dynamics_;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& Cell<T,Descriptor>::getDynamics() const {
    PLB_PRECONDITION(dynamics);
    return *dynamics;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>& Cell<T,Descriptor>::getDynamics() {
    PLB_PRECONDITION(dynamics);
    return *dynamics;
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::revert() {
    for (plint iPop=1; iPop<=Descriptor<T>::q/2; ++iPop) {
        std::swap(f[iPop],f[iPop+Descriptor<T>::q/2]);
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::iniPop() {
    f.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::iniExternal() {
    for (plint iData=0; iData<Descriptor<T>::ExternalField::numScalars; ++iData) {
        *external.get(iData) = T();
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::serialize(T* data) const {
    const plint q = Descriptor<T>::q;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;
    for (plint iPop=0; iPop<q; ++iPop) {
        data[iPop] = f[iPop];
    }
    for (plint iExternal=0; iExternal < numExt; ++iExternal) {
        data[iExternal+q] = *external.get(iExternal);
    }
}

template<typename T, template<typename U> class Descriptor>
void Cell<T,Descriptor>::unSerialize(T const* data) {
    const plint q = Descriptor<T>::q;
    const plint numExt = Descriptor<T>::ExternalField::numScalars;
    for (plint iPop=0; iPop<q; ++iPop) {
        f[iPop] = data[iPop];
    }
    for (plint iExternal=0; iExternal < numExt; ++iExternal) {
        *external.get(iExternal) = data[iExternal+q];
    }
}

}  // namespace plb

#endif  // CELL_HH
