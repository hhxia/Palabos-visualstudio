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

#ifndef SHAN_CHEN_PROCESSOR_2D_H
#define SHAN_CHEN_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "multiPhysics/interparticlePotential.h"

namespace plb {

/// Shan-Chen coupling for multi-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenMultiComponentProcessor2D : public LatticeBoxProcessingFunctional2D<T,Descriptor> {
public:
    ShanChenMultiComponentProcessor2D(T G_);
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T,Descriptor>*> lattices );
    virtual ShanChenMultiComponentProcessor2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
private:
    T G;
};

/// Shan-Chen coupling for single-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenSingleComponentProcessor2D : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    ShanChenSingleComponentProcessor2D(T G_, interparticlePotential::PsiFunction<T>* Psi_);
    ~ShanChenSingleComponentProcessor2D();
    ShanChenSingleComponentProcessor2D(ShanChenSingleComponentProcessor2D<T,Descriptor> const& rhs);
    ShanChenSingleComponentProcessor2D& operator=(ShanChenSingleComponentProcessor2D<T,Descriptor> const& rhs);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice );
    virtual ShanChenSingleComponentProcessor2D<T,Descriptor>* clone() const;
private:
    T G;
    interparticlePotential::PsiFunction<T>* Psi;
};

}

#endif  // SHAN_CHEN_LATTICES_2D_H
