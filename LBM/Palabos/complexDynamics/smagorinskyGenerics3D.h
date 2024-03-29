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

#ifndef SMAGORINSKY_GENERICS_3D_H
#define SMAGORINSKY_GENERICS_3D_H

namespace plb {

template<typename T, template<typename U> class Descriptor, class SmagoFunction>
class StaticSmagorinskyFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    StaticSmagorinskyFunctional3D(SmagoFunction smagoFunction_, T cSmago0_)
        : smagoFunction(smagoFunction_),
          cSmago0(cSmago0_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    T cSmago = smagoFunction (
                            iX+relativeOffset.x,
                            iY+relativeOffset.y,
                            iZ+relativeOffset.z,
                            cSmago0 );
                    Dynamics<T,Descriptor>* baseDynamics = lattice.get(iX,iY,iZ).getDynamics().cloneComposeable();
                    T omega0 = baseDynamics->getOmega();
                    lattice.attributeDynamics ( iX,iY,iZ, 
                            new SmagorinskyDynamics<T,Descriptor> (
                                baseDynamics, omega0, cSmago ) );
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Composite dynamics needs to be instantiated everywhere, including envelope.
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual StaticSmagorinskyFunctional3D<T,Descriptor,SmagoFunction>* clone() const 
    {
        return new StaticSmagorinskyFunctional3D<T,Descriptor,SmagoFunction>(*this);
    }
private:
    SmagoFunction smagoFunction;
    T cSmago0;
};

template<typename T>
T constOmegaFromOmega0(plint iX, plint iY, plint iZ, T omega0) {
    return omega0;
}

template<typename T, template<typename U> class Descriptor>
void instantiateStaticSmagorinsky(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, T cSmago) {
    instantiateStaticSmagorinsky(lattice, domain, constOmegaFromOmega0<T>, cSmago);
}

template<typename T, template<typename U> class Descriptor>
void instantiateStaticSmagorinsky(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, T cSmago) {
    instantiateStaticSmagorinsky(lattice, domain, constOmegaFromOmega0<T>, cSmago);
}

template<typename T, template<typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(BlockLattice3D<T,Descriptor>& lattice, Box3D domain, SmagoFunction smagoFunction, T cSmago0)
{
    applyProcessingFunctional (
            new StaticSmagorinskyFunctional3D<T,Descriptor,SmagoFunction>(smagoFunction, cSmago0),
            domain, lattice );
}

template<typename T, template<typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain, SmagoFunction smagoFunction, T cSmago0)
{
    applyProcessingFunctional (
            new StaticSmagorinskyFunctional3D<T,Descriptor,SmagoFunction>(smagoFunction, cSmago0),
            domain, lattice );
}

} // namespace plb

#endif  // SMAGORINSKY_GENERICS_3D_H
