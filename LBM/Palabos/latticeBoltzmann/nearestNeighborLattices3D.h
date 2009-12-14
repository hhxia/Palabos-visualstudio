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
 * Descriptor for nearest-neighbor 3D lattices. In principle, thanks
 * to static genericity of the code, definition of the lattice and
 * implementation of the dynamics are independent. Of course, there
 * remain a few dependencies. For example, a dynamics which uses higher-
 * order moments such as heat flux needs a lattice with an extended
 * neighborhood. Another example is the D3Q13 lattice, which works only
 * with its very own dynamics.
 *
 *  -- header file
 */
#ifndef NEAREST_NEIGHBOR_LATTICES_3D_H
#define NEAREST_NEIGHBOR_LATTICES_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/externalFields.h"
#include "latticeBoltzmann/roundOffPolicy.h"
#include <vector>

namespace plb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {

     /// D3Q13 lattice
    template <typename T> struct D3Q13Constants
    {
        enum { d = 3, q = 13 };       ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
        static const T lambda_e;      ///< relaxation parameter for the bulk stress
        static const T lambda_h;      ///< additional relaxation parameter
    };

    template <typename T> struct D3Q13DescriptorBase
        : public D3Q13Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D3Q13DescriptorBase<T> BaseDescriptor;
    };

    /// D3Q15 lattice
    template <typename T> struct D3Q15Constants
    {
        enum { d = 3, q = 15 };       ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
    };

    template <typename T> struct D3Q15DescriptorBase
        : public D3Q15Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D3Q15DescriptorBase<T> BaseDescriptor;
    };

    /// D3Q19 lattice
    template <typename T> struct D3Q19Constants
    {
        enum { d = 3, q = 19 };       ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
    };

    template <typename T> struct D3Q19DescriptorBase
        : public D3Q19Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D3Q19DescriptorBase<T> BaseDescriptor;
    };

    /// D3Q27 lattice
    template <typename T> struct D3Q27Constants
    {
        enum { d = 3, q = 27 };       ///< number of dimensions/distr. functions
        static const T invD;          ///< 1 / (number of dimensions)
        static const int vicinity;    ///< size of neighborhood
        static const int c[q][d];     ///< lattice directions
        static const int cNormSqr[q]; ///< norm-square of the vector c
        static const T t[q];          ///< lattice weights
        static const T cs2;           ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;        ///< 1 / cs2
    };

    template <typename T> struct D3Q27DescriptorBase
        : public D3Q27Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D3Q27DescriptorBase<T> BaseDescriptor;
    };

    template <typename T> struct D3Q13Descriptor
        : public D3Q13DescriptorBase<T>, public NoExternalFieldBase
    { };

    template <typename T> struct ForcedD3Q13Descriptor
        : public D3Q13DescriptorBase<T>, public Force3dDescriptorBase
    { };

    template <typename T> struct D3Q15Descriptor
        : public D3Q15DescriptorBase<T>, public NoExternalFieldBase
    { };

    template <typename T> struct ForcedD3Q15Descriptor
        : public D3Q15DescriptorBase<T>, public Force3dDescriptorBase
    { };

    template <typename T> struct D3Q19Descriptor
        : public D3Q19DescriptorBase<T>, public NoExternalFieldBase
    { };

    template <typename T> struct ForcedD3Q19Descriptor
        : public D3Q19DescriptorBase<T>, public Force3dDescriptorBase
    { };

    template <typename T> struct D3Q27Descriptor
        : public D3Q27DescriptorBase<T>, public NoExternalFieldBase
    { };

    template <typename T> struct ForcedD3Q27Descriptor
        : public D3Q27DescriptorBase<T>, public Force3dDescriptorBase
    { };

}  // namespace descriptors

}  // namespace plb

#endif
