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
 * Templates for common geometric operations (scalar product, vector-
 * matric operations etc.).
 *  -- header file
 */
#ifndef GEOMETRIC_OPERATION_TEMPLATES_H
#define GEOMETRIC_OPERATION_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/array.h"
#include <algorithm>
#include <vector>

namespace plb {

template <typename T, int d> struct VectorTemplateImpl;

template <typename T, template<typename U> class Descriptor>
struct VectorTemplate {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Compute scalar product between two vectors
    static T scalarProduct(Array<T,d> const& u1, Array<T,d> const& u2) {
        return VectorTemplateImpl<T,d>::scalarProduct(u1,u2);
    }
    /// Compute scalar product between two a c-array and a plb-array
    static T scalarProduct(const T u1[d], Array<T,d> const& u2) {
        return VectorTemplateImpl<T,d>::scalarProduct(u1,u2);
    }
    /// Compute norm-square of a vector
    static T normSqr(Array<T,d> const& u) {
        return VectorTemplateImpl<T,d>::normSqr(u);
    }
    /// Multiply vector elements component-wise by a scalar
    static void multiplyByScalar(Array<T,d>& u, T scalar) {
        VectorTemplateImpl<T,d>::multiplyByScalar(u,scalar);
    }
    /// Multiply vector elements component-wise by a scalar and store in second vector
    static void multiplyByScalar(Array<T,d> const& u, T scalar, Array<T,d>& result) {
        VectorTemplateImpl<T,d>::multiplyByScalar(u,scalar,result);
    }
};

template <typename T, int d>
struct VectorTemplateImpl {
    static T scalarProduct(Array<T,d> const& u1, Array<T,d> const& u2) {
        T result = T();
        for (int iD=0; iD<d; ++iD) {
            return result += u1[iD]*u2[iD];
        }
        return result;
    }
    static T scalarProduct(const T u1[d], Array<T,d> const& u2) {
        T result = T();
        for (int iD=0; iD<d; ++iD) {
            return result += u1[iD]*u2[iD];
        }
        return result;
    }
    static T normSqr(Array<T,d> const& u) {
        return scalarProduct(u,u);
    }
    static void multiplyByScalar(Array<T,d>& u, T scalar) {
        for (int iD=0; iD<d; ++iD) {
            u[iD] *= scalar;
        }
    }
    static void multiplyByScalar(Array<T,d> const& u, T scalar, Array<T,d>& result) {
        for (int iD=0; iD<d; ++iD) {
            result[iD] = u[iD]*scalar;
        }
    }
};


template <typename T>
struct VectorTemplateImpl<T,2> {
    static T scalarProduct(Array<T,2> const& u1, Array<T,2> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1];
    }
    static T scalarProduct(const T u1[2], Array<T,2> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1];
    }
    static T normSqr(Array<T,2> const& u) {
        return u[0]*u[0] + u[1]*u[1];
    }
    static void multiplyByScalar(Array<T,2>& u, T scalar) {
        u[0] *= scalar;
        u[1] *= scalar;
    }
    static void multiplyByScalar(Array<T,2> const& u, T scalar, Array<T,2>& result) {
        result[0] = u[0]*scalar;
        result[1] = u[1]*scalar;
    }
};


template <typename T>
struct VectorTemplateImpl<T,3> {
    static T scalarProduct(Array<T,3> const& u1, Array<T,3> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
    }
    static T scalarProduct(const T u1[3], Array<T,3> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
    }
    static T normSqr(Array<T,3> const& u) {
        return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    }
    static void multiplyByScalar(Array<T,3>& u, T scalar) {
        u[0] *= scalar;
        u[1] *= scalar;
        u[2] *= scalar;
    }
    static void multiplyByScalar(Array<T,3> const& u, T scalar, Array<T,3>& result) {
        result[0] = u[0]*scalar;
        result[1] = u[1]*scalar;
        result[2] = u[2]*scalar;
    }
};


template <typename T, int d> struct SymmetricTensorImpl { };

template<typename T> struct SymmetricTensorImpl<T,2> {
    static const int n = 3;
    enum Indices { xx=0, xy=1, yy=2 };
    static void matVectMult(Array<T,n> const& mat, Array<T,2> const& vect, Array<T,2>& result) {
        result[0] = mat[xx]*vect[0] + mat[xy]*vect[1];
        result[1] = mat[xy]*vect[0] + mat[yy]*vect[1];
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return mat[xx]*mat[xx]+mat[yy]*mat[yy] + (T)2*mat[xy]*mat[xy];
    }
};

template<typename T> struct SymmetricTensorImpl<T,3> {
    static const int n = 6;
    enum Indices { xx=0, xy=1, xz=2, yy=3, yz=4, zz=5 };
    static void matVectMult(Array<T,n> const& mat, Array<T,3> const& vect, Array<T,3>& result) {
        result[0] = mat[xx]*vect[0] + mat[xy]*vect[1] + mat[xz]*vect[2];
        result[1] = mat[xy]*vect[0] + mat[yy]*vect[1] + mat[yz]*vect[2];
        result[2] = mat[xz]*vect[0] + mat[yz]*vect[1] + mat[zz]*vect[2];
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return mat[xx]*mat[xx]+mat[yy]*mat[yy]+mat[zz]*mat[zz]
               + (T)2*mat[xy]*mat[xy] + (T)2*mat[xz]*mat[xz] + (T)2*mat[yz]*mat[yz];
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template<typename U> class Descriptor>
struct SymmetricTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricTensorImpl<T,d>::n;
    static void matVectMult(Array<T,n> const& mat, Array<T,d> const& vect, Array<T,d>& result) {
        SymmetricTensorImpl<T,d>::matVectMult(mat, vect, result);
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return SymmetricTensorImpl<T,d>::tensorNormSqr(mat);
    }
};

}  // namespace plb

#endif
