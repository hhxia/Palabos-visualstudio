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
 * Helper functions for data field initialization -- header file.
 */
#ifndef DATA_FIELD_INITIALIZER_2D_HH
#define DATA_FIELD_INITIALIZER_2D_HH

#include "simulationSetup/dataFieldInitializer2D.h"
#include "core/array.h"

namespace plb {

template<typename T>
class IniConstScalarFunctional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    IniConstScalarFunctional2D(T value_)
        : value(value_)
    { }
    virtual void process(Box2D domain, ScalarField2D<T>& field) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                field.get(iX,iY) = value;
            }
        }
    }
    virtual IniConstScalarFunctional2D<T>* clone() const {
        return new IniConstScalarFunctional2D<T>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Include boundary right away, to avoid need for envelope update.
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void rescale(T dxScale, T dtScale) { }
private:
    T value;
};


template<typename T>
void setToConstant(ScalarFieldBase2D<T>& field, Box2D domain, T value) {
    applyProcessingFunctional(new IniConstScalarFunctional2D<T>(value), domain, field);
}


template<typename T, int nDim>
class IniConstTensorFunctional2D : public BoxProcessingFunctional2D_T<T,nDim>
{
public:
    IniConstTensorFunctional2D(Array<T,nDim> const& value_)
        : value(value_)
    { }
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                field.get(iX,iY) = value;
            }
        }
    }
    virtual IniConstTensorFunctional2D<T,nDim>* clone() const {
        return new IniConstTensorFunctional2D<T,nDim>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Include boundary right away, to avoid need for envelope update.
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void rescale(T dxScale, T dtScale) { }
private:
    Array<T,nDim> value;
};

template<typename T, int nDim>
void setToConstant( TensorFieldBase2D<T,nDim>& field, Box2D domain,
                    Array<T,nDim> const& value )
{
    applyProcessingFunctional (
            new IniConstTensorFunctional2D<T,nDim>(value), domain, field );
}


template<typename T>
class SetToCoordinateFunctional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    SetToCoordinateFunctional2D(plint index_)
        : index(index_)
    {
        PLB_ASSERT( index >= 0 && index <=1 );
    }
    virtual void process(Box2D domain, ScalarField2D<T>& field) {
        Dot2D relativeOffset = field.getLocation();
        Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
        Array<plint,2> pos;
        for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
            for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
                field.get(pos[0],pos[1]) = (T) (pos[index]+ofs[index]);
            }
        }
    }
    virtual SetToCoordinateFunctional2D<T>* clone() const {
        return new SetToCoordinateFunctional2D<T>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Boundary cannot be included, because periodic boundaries would get the wrong value.
        return BlockDomain::bulk;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void rescale(T dxScale, T dtScale) { }
private:
    plint index;
};


template<typename T>
void setToCoordinate(ScalarFieldBase2D<T>& field, Box2D domain, plint index) {
    applyProcessingFunctional(new SetToCoordinateFunctional2D<T>(index), domain, field);
}

template<typename T>
class SetToCoordinatesFunctional2D : public BoxProcessingFunctional2D_T<T,2>
{
public:
    SetToCoordinatesFunctional2D()
    { }
    virtual void process(Box2D domain, TensorField2D<T,2>& field) {
        Dot2D relativeOffset = field.getLocation();
        Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
        Array<plint,2> pos;
        for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
            for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
                Array<T,2>& cell = field.get(pos[0], pos[1]);
                cell[0] = (T) (pos[0]+ofs[0]);
                cell[1] = (T) (pos[1]+ofs[1]);
            }
        }
    }
    virtual SetToCoordinatesFunctional2D<T>* clone() const {
        return new SetToCoordinatesFunctional2D<T>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Boundary cannot be included, because periodic boundaries would get the wrong value.
        return BlockDomain::bulk;
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void rescale(T dxScale, T dtScale) { }
};


template<typename T>
void setToCoordinates(TensorFieldBase2D<T,2>& field, Box2D domain) {
    applyProcessingFunctional(new SetToCoordinatesFunctional2D<T>, domain, field);
}

}  // namespace plb

#endif  // DATA_FIELD_INITIALIZER_2D_HH

