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
 * Class Array, a fixed-width vector type for general use in Palabos.
 */

#ifndef ARRAY_H
#define ARRAY_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include <algorithm>

using namespace std;

namespace plb {

/** A simple array class, which is slightly more convenient
 *  to use than the pure C-array. It initializes the data to
 *  zero in the default constructor, can be assigned from an
 *  Array of different type, contains bound verifications for
 *  debugging, and is specialized for the 2D and 3D case.
 *  
 *  ATTENTION: Values are not reset to zero in the constructor,
 *  for efficiency reason. If you need them to default to
 *  zero, you need to invoke method resetToZero() explicitly.
 */
template<typename T, pluint size>
class Array {
public:
    Array() { }
    template<typename U, pluint uSize>
    Array(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min(size,uSize), data);
    }
    template<typename U, pluint uSize>
    Array<T,size>& operator=(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min(size,uSize), data);
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<size);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<size);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        std::copy(cArray, cArray+size, data);
    }
    void to_cArray(T* cArray) const {
        std::copy(data, data+size, cArray);
    }
    void resetToZero() {
        std::fill_n(data, size, T());
    }
    Array<T,size>& operator += (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] += b[i];
        }
        return *this;
    }
    Array<T,size>& operator += (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] += alpha;
        }
        return *this;
    }
    Array<T,size>& operator -= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] -= b[i];
        }
        return *this;
    }
    Array<T,size>& operator -= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] -= alpha;
        }
        return *this;
    }
    Array<T,size>& operator *= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] *= b[i];
        }
        return *this;
    }
    Array<T,size>& operator *= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] *= alpha;
        }
        return *this;
    }
    Array<T,size>& operator /= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] /= b[i];
        }
        return *this;
    }
    Array<T,size>& operator /= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] /= alpha;
        }
        return *this;
    }
private:
    T data[size];
};

template<typename T, pluint size>
Array<T,size> operator+(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator+(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] + alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator+(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha + a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] - alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha - a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] * alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha * a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator/(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator/(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] / alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator/(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha / a[i];
    }
    return result;
}


template<typename T>
class Array<T,2> {
public:
    Array() { }
    Array(T x, T y) {
        data[0] = x;
        data[1] = y;
    }
    template<typename U, pluint uSize>
    Array(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min((pluint)2,uSize), data);
    }
    template<typename U, pluint uSize>
    Array<T,2>& operator=(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min((pluint)2,uSize), data);
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<2);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<2);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
    }
    Array<T,2>& operator+=(Array<T,2> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        return *this;
    }
    Array<T,2>& operator+=(T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        return *this;
    }
    Array<T,2>& operator-=(Array<T,2> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        return *this;
    }
    Array<T,2>& operator-=(T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        return *this;
    }
    Array<T,2>& operator*=(Array<T,2> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        return *this;
    }
    Array<T,2>& operator*=(T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        return *this;
    }
    Array<T,2>& operator/=(Array<T,2> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        return *this;
    }
    Array<T,2>& operator/=(T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        return *this;
    }
private:
    T data[2];
};

template<typename T>
Array<T,2> operator+(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]+b[0], a[1]+b[1]);
}

template<typename T>
Array<T,2> operator+(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]+alpha, a[1]+alpha);
}

template<typename T>
Array<T,2> operator+(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha+a[0], alpha+a[1]);
}

template<typename T>
Array<T,2> operator-(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]-b[0], a[1]-b[1]);
}

template<typename T>
Array<T,2> operator-(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]-alpha, a[1]-alpha);
}

template<typename T>
Array<T,2> operator-(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha-a[0], alpha-a[1]);
}

template<typename T>
Array<T,2> operator*(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]*b[0], a[1]*b[1]);
}

template<typename T>
Array<T,2> operator*(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]*alpha, a[1]*alpha);
}

template<typename T>
Array<T,2> operator*(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha*a[0], alpha*a[1]);
}

template<typename T>
Array<T,2> operator/(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]/b[0], a[1]/b[1]);
}

template<typename T>
Array<T,2> operator/(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]/alpha, a[1]/alpha);
}

template<typename T>
Array<T,2> operator/(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha/a[0], alpha/a[1]);
}


template<typename T>
class Array<T,3> {
public:
    Array() { }
    Array(T x, T y, T z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    template<typename U, pluint uSize>
    Array(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min((pluint)3,uSize), data);
    }
    template<typename U, pluint uSize>
    Array<T,3>& operator=(Array<U,uSize> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+min((pluint)3,uSize), data);
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<3);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<3);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
        data[2] = cArray[2];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
        cArray[2] = data[2];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
        data[2] = T();
    }
    Array<T,3>& operator += (Array<T,3> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
        return *this;
    }
    Array<T,3>& operator += (T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        data[2] += alpha;
        return *this;
    }
    Array<T,3>& operator -= (Array<T,3> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        data[2] -= b[2];
        return *this;
    }
    Array<T,3>& operator -= (T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        data[2] -= alpha;
        return *this;
    }
    Array<T,3>& operator *= (Array<T,3> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        data[2] *= b[2];
        return *this;
    }
    Array<T,3>& operator *= (T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        data[2] *= alpha;
        return *this;
    }
    Array<T,3>& operator /= (Array<T,3> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        data[2] /= b[2];
        return *this;
    }
    Array<T,3>& operator /= (T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        data[2] /= alpha;
        return *this;
    }
private:
    T data[3];
};

template<typename T>
Array<T,3> operator+(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

template<typename T>
Array<T,3> operator+(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]+alpha, a[1]+alpha, a[2]+alpha);
}

template<typename T>
Array<T,3> operator+(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha+a[0], alpha+a[1], alpha+a[2]);
}

template<typename T>
Array<T,3> operator-(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

template<typename T>
Array<T,3> operator-(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]-alpha, a[1]-alpha, a[2]-alpha);
}

template<typename T>
Array<T,3> operator-(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha-a[0], alpha-a[1], alpha-a[2]);
}

template<typename T>
Array<T,3> operator*(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

template<typename T>
Array<T,3> operator*(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]*alpha, a[1]*alpha, a[2]*alpha);
}

template<typename T>
Array<T,3> operator*(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha*a[0], alpha*a[1], alpha*a[2]);
}

template<typename T>
Array<T,3> operator/(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}

template<typename T>
Array<T,3> operator/(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]/alpha, a[1]/alpha, a[2]/alpha);
}

template<typename T>
Array<T,3> operator/(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha/a[0], alpha/a[1], alpha/a[2]);
}

} // end namespace plb

#endif
