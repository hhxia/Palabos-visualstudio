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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "core/dataAnalysisFunctionals3D.h"
#include "core/dataAnalysisFunctionals3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

template class ExtractLatticeSubDomainFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxSumRhoBarFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxSumEnergyFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxDensityFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxRhoBarFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxKineticEnergyFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxVelocityComponentFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxVelocityNormFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxVelocityFunctional3D<double, descriptors::D3Q19Descriptor>;
template class BoxPopulationFunctional3D<double, descriptors::D3Q19Descriptor>;


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

template class BoxScalarSumFunctional3D<double>;
template class BoundedBoxScalarSumFunctional3D<double>;

template class BoxScalarMinFunctional3D<double>;
template class BoxScalarMaxFunctional3D<double>;

template class ExtractScalarSubDomainFunctional3D<double>;

template class A_plus_alpha_functional3D<double>;
template class A_minus_alpha_functional3D<double>;
template class Alpha_minus_A_functional3D<double>;
template class A_times_alpha_functional3D<double>;
template class A_dividedBy_alpha_functional3D<double>;
template class Alpha_dividedBy_A_functional3D<double>;

template class A_plus_alpha_inplace_functional3D<double>;
template class A_minus_alpha_inplace_functional3D<double>;
template class A_times_alpha_inplace_functional3D<double>;
template class A_dividedBy_alpha_inplace_functional3D<double>;

template class A_plus_B_functional3D<double>;
template class A_minus_B_functional3D<double>;
template class A_times_B_functional3D<double>;
template class A_dividedBy_B_functional3D<double>;

template class A_plus_B_inplace_functional3D<double>;
template class A_minus_B_inplace_functional3D<double>;
template class A_times_B_inplace_functional3D<double>;
template class A_dividedBy_B_inplace_functional3D<double>;


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template class ExtractTensorSubDomainFunctional3D<double,3>;
template class ExtractTensorComponentFunctional3D<double,3>;
template class ComputeNormFunctional3D<double,3>;
template class ComputeNormSqrFunctional3D<double,3>;
template class BoxBulkVorticityFunctional3D<double,3>;
template class BoxVorticityFunctional3D<double,3>;

template class Tensor_A_plus_B_functional3D<double,3>;
template class Tensor_A_minus_B_functional3D<double,3>;
template class Tensor_A_times_B_functional3D<double,3>;
template class Tensor_A_dividedBy_B_functional3D<double,3>;

template class Tensor_A_plus_B_inplace_functional3D<double,3>;
template class Tensor_A_minus_B_inplace_functional3D<double,3>;
template class Tensor_A_times_B_inplace_functional3D<double,3>;
template class Tensor_A_dividedBy_B_inplace_functional3D<double,3>;

}
