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
 * BounceBack dynamics models in 3D -- generic implementation.
 */
#ifndef BOUNCE_BACK_MODELS_3D_HH
#define BOUNCE_BACK_MODELS_3D_HH

#include "boundaryCondition/bounceBackModels3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include <typeinfo>

namespace plb {

/* *************** Class InitializeMomentumExchangeFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
class InitializeMomentumExchangeFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
                    if (typeid(dynamics) == typeid(MomentumExchangeBounceBack<T,Descriptor>&)) {
                        std::vector<plint> fluidDirections;
                        for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                            plint nextX = iX + Descriptor<T>::c[0];
                            plint nextY = iY + Descriptor<T>::c[1];
                            plint nextZ = iZ + Descriptor<T>::c[2];
                            Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY).getDynamics();
                            if ( typeid(partner) ==
                                 typeid(MomentumExchangeBounceBack<T,Descriptor>&) )
                            {
                                fluidDirections.push_back(iPop);
                            }
                        }

                        if (!fluidDirections.empty()) {
                            MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                                dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                            bounceBackDynamics.setFluidDirections(fluidDirections);
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual InitializeMomentumExchangeFunctional3D<T,Descriptor>* clone() const 
    {
        return new InitializeMomentumExchangeFunctional3D<T,Descriptor>(*this);
    }
};

/* ************* Class MomentumExchangeComplexDomainFunctional3D ** */

template<typename T, template<typename U> class Descriptor>
class MomentumExchangeComplexDomainFunctional3D
    : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    MomentumExchangeComplexDomainFunctional3D( DomainFunctional3D* domain_ )
        : domain(domain_)
    { }
    MomentumExchangeComplexDomainFunctional3D (
            MomentumExchangeComplexDomainFunctional3D<T,Descriptor> const& rhs )
        : domain(rhs.domain->clone())
    { }
    MomentumExchangeComplexDomainFunctional3D<T,Descriptor>& operator= (
            MomentumExchangeComplexDomainFunctional3D<T,Descriptor> const& rhs )
    {
        delete domain; domain = rhs.domain->clone();
    }

    ~MomentumExchangeComplexDomainFunctional3D() {
        delete domain;
    }
    virtual void process(Box3D boundingBox, BlockLattice3D<T,Descriptor>& lattice) {
        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
            for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
                for (plint iZ=boundingBox.z0; iZ<=boundingBox.z1; ++iZ) {
                    if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y,iZ+relativeOffset.z)) {
                        Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
                        if (typeid(dynamics) == typeid(MomentumExchangeBounceBack<T,Descriptor>&)) {
                            std::vector<plint> fluidDirections;
                            for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                                plint nextX = iX + Descriptor<T>::c[iPop][0];
                                plint nextY = iY + Descriptor<T>::c[iPop][1];
                                plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                                Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY,nextZ).getDynamics();
                                if ( typeid(partner) ==
                                     typeid(MomentumExchangeBounceBack<T,Descriptor>&) )
                                {
                                    fluidDirections.push_back(iPop);
                                }
                            }

                            if (!fluidDirections.empty()) {
                                MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                                    dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                                bounceBackDynamics.setFluidDirections(fluidDirections);
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual MomentumExchangeComplexDomainFunctional3D<T,Descriptor>* clone() const 
    {
        return new MomentumExchangeComplexDomainFunctional3D<T,Descriptor>(*this);
    }
private:
    DomainFunctional3D* domain;
};

template<typename T, template<typename U> class Descriptor>
class InitializeDotMomentumExchangeFunctional3D : public DotProcessingFunctional3D_L<T,Descriptor> {
public:
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice) {
        for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
            Dot3D const& dot = dotList.getDot(iDot);
            plint iX = dot.x;
            plint iY = dot.y;
            plint iZ = dot.z;
            Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
            if (typeid(dynamics) == typeid(MomentumExchangeBounceBack<T,Descriptor>&)) {
                std::vector<plint> fluidDirections;
                for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY,nextZ).getDynamics();
                    if ( typeid(partner) ==
                         typeid(MomentumExchangeBounceBack<T,Descriptor>&) )
                    {
                        fluidDirections.push_back(iPop);
                    }
                }

                if (!fluidDirections.empty()) {
                    MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                        dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                    bounceBackDynamics.setFluidDirections(fluidDirections);
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual InitializeDotMomentumExchangeFunctional3D<T,Descriptor>* clone() const 
    {
        return new InitializeDotMomentumExchangeFunctional3D<T,Descriptor>(*this);
    }
};


template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain )
{
    applyProcessingFunctional (
        new InitializeMomentumExchangeFunctional3D<T,Descriptor>(), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase3D<T,Descriptor>& lattice, Box3D boundingBox,
        DomainFunctional3D* domain )
{
    applyProcessingFunctional (
            new MomentumExchangeComplexDomainFunctional3D<T,Descriptor>(domain), boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLatticeBase3D<T,Descriptor>& lattice, DotList3D const& dotList )
{
    applyProcessingFunctional (
        new InitializeDotMomentumExchangeFunctional3D<T,Descriptor>(), dotList, lattice );
}

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_3D_H
