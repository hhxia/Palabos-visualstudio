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
 * Helper functions for data field initialization -- generic implementation.
 */
#ifndef DATA_FIELD_INITIALIZER_GENERICS_2D_H
#define DATA_FIELD_INITIALIZER_GENERICS_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"

namespace plb {

template<typename T, class Function>
class SetToScalarFunctionFunctional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    SetToScalarFunctionFunctional2D(Function f_)
        : f(f_)
    { }
    virtual void process(Box2D domain, ScalarField2D<T>& field) {
        Dot2D relativeOffset = field.getLocation();
        Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
        Array<plint,2> pos;
        for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
            for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
                field.get(pos[0],pos[1]) = 
                    f(pos[0]+ofs[0], pos[1]+ofs[1]);
            }
        }
    }
    virtual SetToScalarFunctionFunctional2D<T,Function>* clone() const {
        return new SetToScalarFunctionFunctional2D<T,Function>(*this);
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
    Function f;
};

template<typename T, class Function>
void setToFunction(ScalarFieldBase2D<T>& field, Box2D domain, Function f) {
    applyProcessingFunctional (
            new SetToScalarFunctionFunctional2D<T,Function>(f), domain, field );
}

template<typename T, int nDim, class Function>
class SetToTensorFunctionFunctional2D : public BoxProcessingFunctional2D_T<T,nDim>
{
public:
    SetToTensorFunctionFunctional2D(Function f_)
        : f(f_)
    { }
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field) {
        Dot2D relativeOffset = field.getLocation();
        Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
        Array<plint,2> pos;
        Array<T,nDim> value;
        for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
            for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
                f(pos[0]+ofs[0], pos[1]+ofs[1], value);
                field.get(pos[0],pos[1]) = value;
            }
        }
    }
    virtual SetToTensorFunctionFunctional2D<T,nDim,Function>* clone() const {
        return new SetToTensorFunctionFunctional2D<T,nDim,Function>(*this);
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
    Function f;
};

template<typename T, int nDim, class Function>
void setToFunction(TensorFieldBase2D<T,nDim>& field, Box2D domain, Function f) {
    applyProcessingFunctional (
            new SetToTensorFunctionFunctional2D<T,nDim,Function>(f), domain, field );
}

}  // namespace plb

#endif  // DATA_FIELD_INITIALIZER_GENERICS_2D_H
