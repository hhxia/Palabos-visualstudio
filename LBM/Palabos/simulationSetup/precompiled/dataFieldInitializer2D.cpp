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

#include "simulationSetup/dataFieldInitializer2D.h"
#include "simulationSetup/dataFieldInitializer2D.hh"

namespace plb {

    template class IniConstScalarFunctional2D<double>;
    template void setToConstant(ScalarFieldBase2D<double>& field, Box2D domain, double value);
    template class IniConstTensorFunctional2D<double,2>;
    template void setToConstant(TensorFieldBase2D<double,2>& field, Box2D domain, Array<double,2> const& value);
    template class SetToCoordinateFunctional2D<double>;
    template void setToCoordinate(ScalarFieldBase2D<double>& field, Box2D domain, plint index);
    template class SetToCoordinatesFunctional2D<double>;
    template void setToCoordinates(TensorFieldBase2D<double,2>& field, Box2D domain);

}  // namespace plb
