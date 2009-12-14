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

#ifndef SHAN_CHEN_PROCESSOR_3D_H
#define SHAN_CHEN_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiPhysics/interparticlePotential.h"

namespace plb {

/// Shan-Chen coupling for multi-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenMultiComponentProcessor3D : public LatticeBoxProcessingFunctional3D<T,Descriptor> {
public:
    ShanChenMultiComponentProcessor3D(T G_);
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices );
    virtual ShanChenMultiComponentProcessor3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    T G;
};

/// Shan-Chen coupling for single-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenSingleComponentProcessor3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    ShanChenSingleComponentProcessor3D(T G_, interparticlePotential::PsiFunction<T>* Psi_);
    ~ShanChenSingleComponentProcessor3D();
    ShanChenSingleComponentProcessor3D(ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs);
    ShanChenSingleComponentProcessor3D& operator=(ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice );
    virtual ShanChenSingleComponentProcessor3D<T,Descriptor>* clone() const;
private:
    T G;
    interparticlePotential::PsiFunction<T>* Psi;
};

}

#endif  // SHAN_CHEN_LATTICES_3D_H
