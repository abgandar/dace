/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2016 Politecnico di Milano (2014 Dinamica Srl)                    *
* Licensed under the Apache License, Version 2.0 (the "License");             *
* you may not use this file except in compliance with the License.            *
* You may obtain a copy of the License at                                     *
*                                                                             *
*    http://www.apache.org/licenses/LICENSE-2.0                               *
*                                                                             *
* Unless required by applicable law or agreed to in writing, software         *
* distributed under the License is distributed on an "AS IS" BASIS,           *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    *
* See the License for the specific language governing permissions and         *
* limitations under the License.                                              *
*                                                                             *
*******************************************************************************/

/*
 * AlgebraicVector_t.h
 *
 *  Created on: Sep. 10, 2014
 *      Author: Dinamica Srl
 */

/*  Templated function definitions for AlgebraicVector class.

    This header file contains the definition of all functions in the (templated) AlgebraicVector class.
*/

#ifndef DINAMICA_ALGEBRAICVECTOR_T_H_
#define DINAMICA_ALGEBRAICVECTOR_T_H_

// C++ stdlib classes used only internally in this implementation
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// DACE classes
#include "dace/PromotionTrait.h"
#include "dace/MathExtension.h"
#include "dace/compiledDA.h"
#include "dace/AlgebraicVector.h"
#ifdef WITH_ALGEBRAICMATRIX
#include "dace/AlgebraicMatrix.h"
#include "dace/AlgebraicMatrix_t.h"
#endif /* WITH_ALGEBRAICMATRIX */

namespace DACE {

/********************************************************************************
*     Element access
*********************************************************************************/
/** Extracts elements from AlgebraicVector.
    @note This routine performs range checking and throws an error if the indices are out of range.
    @param[in] first The index of the first element to be extracted.
    @param[in] last  The index of the last element to be extracted.
    @return A new AlgebraicVector with elements from position first to last.
    @throw std::runtime_error
*/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::extract(const size_t first, const size_t last) const {
    if(first>=this->size() || last>=this->size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::extract: Indices out of bounds.");

    return AlgebraicVector<T>(*this, first, last);
}

/** Append an AlgebraicVector to the end of the current one and return the new vector.
    @param[in] obj The AlgebraicVector to be appended.
    @return A new AlgebraicVector containing the elements of both vectors, cast upwards if necessary.
*/
template<typename T> template<typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> AlgebraicVector<T>::concat(const std::vector<U> &obj) const {
    const size_t size1 = this->size();
    const size_t size2 = obj.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> res(size1+size2);

    for(size_t i=0; i<size1; i++)
        res[i] = (*this)[i];
    for(size_t i=0; i<size2; i++)
        res[i+size1] = obj[i];

    return res;
}

/** Append elements of vector @e obj to the end of ourself.
    Converts the type @e U to match our @e T if necessary.
    @param[in] obj A vector of elements to append.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator<<(const std::vector<U> &obj) {
    const size_t size = obj.size();
    for(size_t i=0; i<size; i++) {
        (*this).push_back((T)obj[i]);
    }
    return *this;
}

/***********************************************************************************
*     Coefficient access routines
************************************************************************************/
/** Return the constant parts of each element.
    @return AlgebraicVector<double> containing the constant part of each element.
 */
template<typename T> AlgebraicVector<double> AlgebraicVector<T>::cons() const {
    using DACE::cons;

    const size_t size = this->size();
    AlgebraicVector<double> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cons((*this)[i]);
    }
    return temp;
}

#ifdef WITH_ALGEBRAICMATRIX
/** Return the linear part of a polynomial map in AlgebraicVector<T>.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @return A @p AlgebraicMatrix<double> of dimension size by nvar, where size is the
    size of the @p AlgebraicVector<T> considered and nvar is the number of variables defined
    during the DACE initialization. Each row contains the linear part of the corresponding
    DA included in the original @p AlgebraicVector<T>.
 */
template<typename T> AlgebraicMatrix<double> AlgebraicVector<T>::linear() const {
    const size_t size = this->size();
    const int nvar = DA::getMaxVariables();

    AlgebraicMatrix<double> out(size, nvar);
    for(size_t i=0; i<size; i++) {
          out.setrow(i, (*this)[i].linear());
    }
    return out;
}
#else
/** Return the linear part of a polynomial map in AlgebraicVector<T>.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @return A std::vector<std::vector<double>>, where each std::vector<double> contains.
    the linear part of the corresponding element in the original AlgebraicVector<T>.
 */
template<typename T> std::vector<std::vector<double>> AlgebraicVector<T>::linear() const {
    const size_t size = this->size();

    std::vector< std::vector<double> > out(size);
    for(size_t i=0; i<size; i++) {
          out[i] = (*this)[i].linear();
    }
    return out;
}
#endif /* WITH_ALGEBRAICMATRIX */

/********************************************************************************
*     Assignments, Copying & Filtering
*********************************************************************************/
/** Add ourselves to the given AlgebraicVector componentwise.
    @param[in] obj An AlgebraicVector.
    @return A reference to ourselves.
    @throw std::runtime_error
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator+=(const AlgebraicVector<U> &obj) {
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator+=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] += obj[i];
    }
    return *this;
}

/** Add ourselves to the given scalar componentwise.
    @param[in] obj A scalar value.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator+=(const U &obj) {
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] += obj;
    }
    return *this;
}

/** Subtract the given AlgebraicVector from ourselves componentwise.
    @param[in] obj An AlgebraicVector.
    @return A reference to ourselves.
    @throw std::runtime_error
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator-=(const AlgebraicVector<U> &obj) {
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator-=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] -= obj[i];
    }
    return *this;
}

/** Subtract the given scalar from ourselves componentwise.
    @param[in] obj A scalar value.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator-=(const U &obj) {
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] -= obj;
    }
    return *this;
}

/** Multiply ourselves with the given AlgebraicVector componentwise.
    @param[in] obj An AlgebraicVector.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator*=(const AlgebraicVector<U> &obj) {
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator*=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] *= obj[i];
    }
    return *this;
}

/** Multiply ourselves with the given scalar.
    @param[in] obj A scalar value.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator*=(const U &obj) {
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] *= obj;
    }
    return *this;
}

/** Divide ourselves by the given AlgebraicVector componentwise.
    @param[in] obj An AlgebraicVector.
    @return A reference to ourselves.
    @throw std::runtime_error
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator/=(const AlgebraicVector<U> &obj) {
    const size_t size = this->size();
    if(size != obj.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator/=: Vectors must have the same length.");

    for(size_t i=0; i<size; i++) {
        (*this)[i] /= obj[i];
    }
    return *this;
}

/** Divide ourselves by the given scalar.
    @param[in] obj A scalar value.
    @return A reference to ourselves.
 */
template<typename T> template<typename U> AlgebraicVector<T>& AlgebraicVector<T>::operator/=(const U &obj) {
    const size_t size = this->size();
    for(size_t i=0; i<size; i++) {
        (*this)[i] /= obj;
    }
    return *this;
}

/** Returns an AlgebraicVector<T> with all monomials of order less than @e min and greater
    than @e max removed (trimmed).
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] min The minimum order to be preserved.
    @param[in] max The maximum order to be preserved.
    @return A new AlgebraicVector<T>.
*/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::trim(const unsigned int min, const unsigned int max) const {
    AlgebraicVector<T> tmp(this->size());

    for(size_t i=0; i<this->size(); i++) {
        tmp[i] = (*this)[i].trim(min, max);
    }

    return tmp;
}

/********************************************************************************
*     Basic arithmetic operations
*********************************************************************************/
/** Returns the additive inverse of the vector.
    @return A new AlgebraicVector with the opposite sign.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::operator-() const {
    return -1.0*(*this);
}

/** Compute the derivative of a AlgebraicVector<T> with respect to variable @e p.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] p The independent variable number with respect to which the derivative is calculated.
    @return A new AlgebraicVector<T>.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::deriv(const unsigned int p) const {
    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].deriv(p);
    }

    return temp;
}

/** Compute the integral of a AlgebraicVector<T> with respect to variable @e p.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] p The independent variable number with respect to which the integral is calculated.
    @return A new AlgebraicVector<T>.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::integ(const unsigned int p) const {
    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].integ(p);
    }

    return temp;
}

/********************************************************************************
*     Intrinsic functions
*********************************************************************************/
/** Componentwise application of the absolute value function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::absolute() const {
    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = absolute((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the truncation function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::trunc() const {
    using std::trunc;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = trunc((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the round function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::round() const {
    using std::round;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = round((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the mod function.
    @param[in] p divisor.
    @return A new AlgebraicVector.
 */
template<typename T> template<typename U> AlgebraicVector<T> AlgebraicVector<T>::mod(const U &p) const {
    using DACE::mod;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = mod((*this)[i], p);
    }
    return temp;
}

/** Componentwise application of the integer power function.
    @param[in] p The power to raise each element to.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::pow(const int p) const {
    using std::pow;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = pow((*this)[i], p);
    }
    return temp;
}

/** Componentwise application of the double power function.
    @param[in] p The power to raise each element to.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::pow(const double p) const {
    using std::pow;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = pow((*this)[i], p);
    }
    return temp;
}

/** Componentwise application of the p-th root function.
    @param[in] p The root to be computed.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::root(const int p) const {
    using DACE::root;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = root((*this)[i],p);
    }
    return temp;
}

/** Componentwise application of the multiplicative inverse function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::minv() const {
    using DACE::minv;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = minv((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the square function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sqr() const {
    using DACE::sqr;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqr((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the square root function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sqrt() const {
    using std::sqrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqrt((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the inverse square root function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::isrt() const {
    using DACE::isrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = isrt((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the cube root function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cbrt() const {
    using std::cbrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cbrt((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the inverse cube root function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::icbrt() const {
    using DACE::icbrt;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sqrt((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hypotenuse function hypot(x,y).
    @param[in] Y The AlgebraicVector<T> second argument.
    @return A new AlgebraicVector.
*/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::hypot(const AlgebraicVector<T> &Y) const {
    using std::hypot;

    const size_t size = this->size();
    if(Y.size() != size)
        throw std::runtime_error("DACE::AlgebraicVector<T>::hypot(): Vectors must have the same length.");

    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = hypot((*this)[i], Y[i]);
    }
    return temp;
}

/** Componentwise application of the exponential function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::exp() const {
    using std::exp;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = exp((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the natural logarithm function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log() const {
    using std::log;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the logarithm function relative to given base.
    @param[in] b The base for the logarithm.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::logb(const double b) const {
    using DACE::logb;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = logb((*this)[i], b);}
    return temp;
}

/** Componentwise application of the decadic logarithm function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log10() const {
    using std::log10;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log10((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the binary logarithm function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::log2() const {
    using std::log2;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = log2((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the sine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sin() const {
    using std::sin;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sin((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the cosine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cos() const {
    using std::cos;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cos((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the tangent function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::tan() const {
    using std::tan;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = tan((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the arcsine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::asin() const {
    using std::asin;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = asin((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the arccosine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::acos() const {
    using std::acos;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = acos((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the arctangent function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atan() const {
    using std::atan;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atan((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the four-quadrant arctangent of @p Y/X where
    @p Y is the current object and @p X is the argument.
    @param[in] X The AlgebraicVector<T> representing X.
    @return A new AlgebraicVector with elements in @f$ [-pi, pi] @f$ .
    @throw std::runtime_error
*/
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atan2(const AlgebraicVector<T> &X) const {
    using std::atan2;

    const size_t size = this->size();
    if(X.size() != size)
        throw std::runtime_error("DACE::AlgebraicVector<T>::atan2(): Vectors must have the same length.");

    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atan2((*this)[i], X[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic sine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::sinh() const {
    using std::sinh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = sinh((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic cosine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::cosh() const {
    using std::cosh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = cosh((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic tangent function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::tanh() const {
    using std::tanh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = tanh((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic arcsine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::asinh() const {
    using std::asinh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = asinh((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic arccosine function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::acosh() const {
    using std::acosh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = acosh((*this)[i]);
    }
    return temp;
}

/** Componentwise application of the hyperbolic arctangent function.
    @return A new AlgebraicVector.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::atanh() const {
    using std::atanh;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = atanh((*this)[i]);
    }
    return temp;
}

/***********************************************************************************
*    Vector routines
************************************************************************************/
/** Compute the dot product with another AlgebraicVector.
    @param[in] obj The other AlgebraicVector.
    @return A scalar value representing dot (inner) product.
    @throw std::runtime_error
 */
template<typename T> template<typename U> typename PromotionTrait<T,U>::returnType AlgebraicVector<T>::dot(const AlgebraicVector<U> &obj) const {
    const size_t size = this->size();
    if(size != obj.size())
          throw std::runtime_error("DACE::AlgebraicVector<T>::dot(): Vectors must have the same length.");

    typename PromotionTrait<T,U>::returnType temp = 0.0;
    for(size_t i=0; i<size; i++) {
        temp += (*this)[i] * obj[i];
    }
    return temp;
}

/** Compute the cross product with another 3D AlgebraicVector.
    @param[in] obj The other AlgebraicVector.
    @return A new AlgebraicVector.
    @throw std::runtime_error
 */
template<typename T> template<typename U> AlgebraicVector<typename PromotionTrait<T,U>::returnType> AlgebraicVector<T>::cross(const AlgebraicVector<U> &obj) const {
    if((this->size() != 3) || (obj.size() != 3))
        throw std::runtime_error("DACE::AlgebraicVector<T>::cross(): Inputs must be 3 element AlgebraicVectors.");

    AlgebraicVector<typename PromotionTrait<T,U>::returnType> temp(3);

    temp[0] = ((*this)[1] * obj[2]) - ((*this)[2] * obj[1]);
    temp[1] = ((*this)[2] * obj[0]) - ((*this)[0] * obj[2]);
    temp[2] = ((*this)[0] * obj[1]) - ((*this)[1] * obj[0]);
    return temp;
}

/** Compute the length (Euclidean vector norm).
    @return Euclidean length of the vector.
 */
template<typename T> T AlgebraicVector<T>::length() const {
    using std::sqrt; using DACE::sqr;       // Implementational note: these using statements are very subtle and absolutely needed.
                                            // They force the compiler to perform argument dependent lookup (ADL) which then finds
                                            // the correct sqrt() and sqr() functions even if they are not in DACE:: or std::!
    const size_t size = this->size();
    T norm = 0.0;
    for(size_t i=0; i<size; i++) {
        norm = norm + sqr((*this)[i]);
    }

    return sqrt(norm);
}

/** Normalize the vector.
    @return An AlgebraicVector<T> of unit length.
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::normalize() const {
    using DACE::minv;

    AlgebraicVector<T> temp(*this);
    temp *= minv(this->length());

    return temp;
}

/** Invert the polynomial map represented by this AlgebraicVector<DA>.
    Map inversion requires that the linear part of this polynomial map is invertible.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @return The inverted polynomial map.
    @throw std::runtime_error
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::invert() const {
    const unsigned int ord = DA::getTO();
    const size_t nvar = this->size();

    if(nvar > DA::getMaxVariables())
        throw std::runtime_error("DACE::AlgebraicVector<T>::inverse: dimension of vector exceeds maximum number of DA variables.");

    // Create DA identity
    AlgebraicVector<T> DDA = AlgebraicVector<T>::identity(nvar);

    // Split map into constant part AC, non-constant part M, and non-linear part AN
    AlgebraicVector<double> AC = this->cons();
    AlgebraicVector<T> M = this->trim(1);
    AlgebraicVector<T> AN = M.trim(2);

#ifdef WITH_ALGEBRAICMATRIX
    // Extract the linear coefficients matrix
    AlgebraicMatrix<double> AL = M.linear();

    // Compute the inverse of linear coefficients matrix
    AlgebraicMatrix<double> AI = AL.inv();

    // Compute DA representation of the inverse of the linear part of the map and its composition with non-linear part AN
    compiledDA AIoAN(AI*AN);
    AlgebraicVector<DA> Linv = AI*DDA;
#else
    // Compute the inverse of linear coefficients matrix
    std::vector<std::vector<double>> AI = M.linear();
    matrix_inverse(AI);

    // Compute DA representation of the inverse of the linear part of the map and its composition with non-linear part AN
    AlgebraicVector<T> Linv(nvar);
    // Linv = AI*AN
    for(size_t i=0; i<nvar; i++) {
        Linv[i] = 0.0;
        for(size_t j=0; j<nvar; j++)
            Linv[i] += AI[i][j]*AN[j];
    }
    compiledDA AIoAN(Linv);
    // Linv = AI*DDA
    for(size_t i=0; i<nvar; i++) {
        Linv[i] = 0.0;
        for(size_t j=0; j<nvar; j++)
            Linv[i] += AI[i][j]*DDA[j];
    }
#endif /* WITH_ALGEBRAICMATRIX */

    // Iterate to obtain the inverse map
    AlgebraicVector<T> MI = Linv;
    for(unsigned int i=1; i<ord; i++) {
        DA::setTO(i+1);
        MI = Linv - AIoAN.eval(MI);
    }

    return MI.eval(DDA-AC);
}

/***********************************************************************************
*     Polynomial evaluation routines
************************************************************************************/
/** Evaluate each element of a vector of DA with any vector type @e U containing the arguments
    and return a vector of results of the same type @e U.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] args A vector (e.g. AlgebraicVector<> or std::vector<>) of arguments.
    @return A new vector of same type as argument args containing the results of the evaluation

    @see compiledDA
    @see AlgebraicVector::compile()
 */
template<typename T> template<typename U> U AlgebraicVector<T>::operator()(const U &args) const {
    return compiledDA(*this)(args);
}

/** Evaluate each element of a vector of DA with a braced initializer list of type @e U
    and return an AlgebraicVector<U> with the results.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] l A braced initializer list containing the arguments of type U.
    @return A new AlgebraicVector<U> containing the results of the evaluation.

    @see compiledDA
    @see AlgebraicVector::compile()
 */
template<typename T> template<typename U> AlgebraicVector<U> AlgebraicVector<T>::operator()(const std::initializer_list<U> l) const {
    return compiledDA(*this)(l);
}

/** Evaluate each element of a vector of DA with an array of arithmetic type @e T arguments.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] args A C array of arithmetic type @e T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array @e args.
    @return A new AlgebraicVector<U> containing the results of the evaluation.

    @see AlgebraicVector::compile()
    @see compiledDA
 */
template<typename T> template<typename U> AlgebraicVector<U> AlgebraicVector<T>::operator()(const U args[], const unsigned int length) const {
    return compiledDA(*this)(args, length);
}

/** Evaluate each element of a vector of DA with any vector type @e U containing the arguments
    and return a vector of results of the same type @e U.
    @deprecated Replaced by AlgebraicVector::operator()().
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] args A vector (e.g. AlgebraicVector<> or std::vector<>) of arguments.
    @return A new vector of same type as argument @e args containing the results of the evaluation

    @see AlgebraicVector::operator()()
    @see AlgebraicVector::compile()
    @see compiledDA
 */
template<typename T> template<typename U> U AlgebraicVector<T>::eval(const U &args) const {
    return compiledDA(*this)(args);
}

/** Evaluate each element of a vector of DA with a braced initializer list of type @e U
    and return an AlgebraicVector<U> with the results.
    @deprecated Replaced by AlgebraicVector::operator()().
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] l A braced initializer list containing the arguments of type @e U.
    @return A new AlgebraicVector<U> containing the results of the evaluation.

    @see AlgebraicVector::operator()()
    @see AlgebraicVector::compile()
    @see compiledDA
 */
template<typename T> template<typename U> AlgebraicVector<U> AlgebraicVector<T>::eval(const std::initializer_list<U> l) const {
    return compiledDA(*this)(l);
}

/** Evaluate each element of a vector of DA with an array of arithmetic type T arguments.
    @deprecated Replaced by AlgebraicVector::operator()().
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] args A C array of arithmetic type T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array args.
    @return A new AlgebraicVector<U> containing the results of the evaluation.

    @see AlgebraicVector::operator()()
    @see AlgebraicVector::compile()
    @see compiledDA
 */
template<typename T> template<typename U> AlgebraicVector<U> AlgebraicVector<T>::eval(const U args[], const unsigned int length) const {
    return compiledDA(*this)(args, length);
}

/** Evaluate each element of a vector of DA with a single arithmetic type U argument.
    @deprecated Replaced by AlgebraicVector::operator()() with braced initializer list (e.g. `da({arg})`).
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] arg A single variable of arithmetic type @e U representing the first independent DA variable.
    @return The result of the evaluation.

    @see AlgebraicVector::compile()
    @see compiledDA
 */
template<typename T> template<typename U> AlgebraicVector<U> AlgebraicVector<T>::evalScalar(const U &arg) const {
    return compiledDA(*this)({arg});
}

/** Compile vector of DA polynomials and create a compiledDA object.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @return The compiled DA object.
    @see DA::compile()
 */
template<typename T> compiledDA AlgebraicVector<T>::compile() const {
    return compiledDA(*this);
}

/** Partial evaluation of vector of DA polynomials. For each element of the vector,
    variable @e var is replaced by the value @e val. The resulting vector of DA is returned.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] var The independent variable number to be replaced.
    @param[in] val The value by which to replace the variable.
    @return A new AlgebraicVector<DA>.

    @see DA::plug()
 */
template<typename T> AlgebraicVector<T> AlgebraicVector<T>::plug(const unsigned int var, const double val) const {
    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = (*this)[i].plug(var,val);
    }

    return temp;
}

/***********************************************************************************
*     Norm routines
************************************************************************************/
/** Componentwise application of the norm function.
    @param[in] type The type of norm to be computed.
    @return A new AlgebraicVector.
    @see DA::norm()
 */
template<typename T> AlgebraicVector<double> AlgebraicVector<T>::norm(const unsigned int type) const {
    using DACE::norm;

    const size_t size = this->size();
    AlgebraicVector<T> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = norm((*this)[i], type);
    }
    return temp;
}

/***********************************************************************************
*     Input/Output routines
************************************************************************************/
/** Convert to string.
    @return String representing the AlgebraicVector.
 */
template<typename T> std::string AlgebraicVector<T>::toString() const {
    std::ostringstream strs;
    strs << *this << std::endl;

    return strs.str();
}

/********************************************************************************
*     Static creation routines
*********************************************************************************/
/** Return the DA identity of dimension @e n.
    @param[in] n The dimension.
    @return AlgebraicVector<DA> containing the DA identity in @e n dimensions.
 */
template<typename T> AlgebraicVector<DA> AlgebraicVector<T>::id(const size_t n) {
    AlgebraicVector<DA> temp(n);
    for(size_t i=0; i < n; i++) {
        temp[i] = DA::id(i+1);
    }

    return temp;
}

/** Return the DA identity of dimension @e n.
    Legacy alias for AlgebraicVector::id().
    @deprecated Replaced by AlgebraicVector::id().

    @param[in] n The dimension.
    @return AlgebraicVector<DA> containing the DA identity in @e n dimensions.
 */
template<typename T> AlgebraicVector<DA> AlgebraicVector<T>::identity(const size_t n) {
    return id(n);
}

/********************************************************************************
*     Private routines
*********************************************************************************/
#ifndef WITH_ALGEBRAICMATRIX
/** @cond */
/* Internal routine to compute a matrix inverse of a double precision matrix.
   Algorithm based on the Gauss elimination with full pivot (from the Numerical
   Cookbook) adapted for C++. This is the same algorithm but a different
   implementation than in the AlgebraicMatrix class.
   This is NOT intended for public use. Limited error checking is performed
   in accordance with the exclusive use of this routine in map inversion.
 */
template<typename T> void AlgebraicVector<T>::matrix_inverse(std::vector<std::vector<double>> &A) {
    using std::abs;

    const size_t n = A.size();
    std::vector<size_t> indexc(n), indexr(n), ipiv(n, 0);

    for (size_t i=0; i<n; i++) {
        size_t icol = 0, irow = 0;
        double big = 0.0;
        for (size_t j=0; j<n; j++)
            if (ipiv[j] == 0)
                for (size_t k=0; k<n; k++)
                    if (ipiv[k] == 0)
                        if (abs(A[j][k]) >= big) {
                            big = abs(A[j][k]);
                            irow = j;
                            icol = k;
                        }
        ipiv[icol] = 1;
        if (irow != icol)
            for (size_t l=0; l<n; l++) std::swap(A[irow][l], A[icol][l]);
        indexr[i] = irow;
        indexc[i] = icol;
        if (A[icol][icol] == 0.0) throw std::runtime_error("DACE::AlgebraicVector<T>::inverse: linear matrix inverse does not exist.");
        const double pivinv = 1.0/A[icol][icol];
        A[icol][icol] = 1.0;
        for (size_t l=0; l<n; l++) A[icol][l] *= pivinv;
        for (size_t ll=0; ll<n; ll++)
            if (ll != icol) {
                const double temp = A[ll][icol];
                A[ll][icol] = 0.0;
                for (size_t l=0; l<n; l++) A[ll][l] -= A[icol][l]*temp;
            }
    }

    for (size_t i=n; i>0; i--)
        if (indexr[i-1] != indexc[i-1])
            for (size_t k=0; k<n; k++)
                std::swap(A[k][indexr[i-1]], A[k][indexc[i-1]]);
}
/** @endcond */
#endif /* WITH_ALGEBRAICMATRIX */


/***********************************************************************************
*
*     Non-member functions
*
************************************************************************************/

/***********************************************************************************
*     Coefficient Access Functions
************************************************************************************/
/** Return the constant parts.
    @return An AlgebraicVector<double>.
    @see AlgebraicVector<T>::cons
 */
template<typename T> AlgebraicVector<double> cons(const AlgebraicVector<T> &obj) {
    return obj.cons();
}

#ifdef WITH_ALGEBRAICMATRIX
/** Return the linear part of a polynomial map in AlgebraicVector<T>.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] obj The AlgebraicVector<T> to extract linear part from.
    @return A AlgebraicMatrix<double> of dimension @e size by @e nvar, where @e size is the
    size of the AlgebraicVector<T> considered and @e nvar is the number of variables defined
    during the DACE initialization. Each row contains the linear part of the corresponding
    DA included in the original AlgebraicVector<T>.

    @see AlgebraicVector<T>::linear
 */
template<typename T> AlgebraicMatrix<double> linear(const AlgebraicVector<T> &obj) {
#else
/** Return the linear part of a polynomial map in AlgebraicVector<T>.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.

    @param[in] obj The AlgebraicVector<T> to extract linear part from.
    @return A std::vector<std::vector<double>>, where each std::vector<double> contains.
    the linear part of the corresponding element in the original AlgebraicVector<T>.

    @see AlgebraicVector<T>::linear
 */
template<typename T> std::vector<std::vector<double>> linear(const AlgebraicVector<T> &obj) {
#endif /* WITH_ALGEBRAICMATRIX */
    return obj.linear();
}

/***********************************************************************************
*     Basic Arithmetic Operators
************************************************************************************/
/** Componentwise addition between two AlgebraicVectors.
    @param[in] obj1 The first AlgebraicVector.
    @param[in] obj2 The second AlgebraicVector.
    @return A new AlgebraicVector.
    @throw std::runtime_error
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator+(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator+: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] + obj2[i];
    }
    return temp;
}

/** Componentwise addition between a AlgebraicVector and a scalar value.
    @param[in] obj1 An AlgebraicVector.
    @param[in] obj2 A scalar value.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator+(const AlgebraicVector<T> &obj1, const U &obj2) {
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] + obj2;
    }
    return temp;
}

/** Componentwise addition between a scalar value and a AlgebraicVector.
    @param[in] obj1 A scalar value.
    @param[in] obj2 A AlgebraicVector.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator+(const T &obj1, const AlgebraicVector<U> &obj2) {
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 + obj2[i];
    }
    return temp;
}

/** Componentwise subtraction between two AlgebraicVectors.
    @param[in] obj1 The first AlgebraicVector.
    @param[in] obj2 The second AlgebraicVector.
    @return A new AlgebraicVector.
    @throw std::runtime_error
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator-(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator-: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] - obj2[i];
    }
    return temp;
}

/** Componentwise subtraction between a AlgebraicVector and a scalar value.
    @param[in] obj1 A AlgebraicVector.
    @param[in] obj2 A scalar value.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator-(const AlgebraicVector<T> &obj1, const U &obj2) {
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] - obj2;
    }
    return temp;
}

/** Componentwise subtraction between a scalar value and a AlgebraicVector.
    @param[in] obj1 A scalar value.
    @param[in] obj2 A AlgebraicVector.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator-(const T &obj1, const AlgebraicVector<U> &obj2) {
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 - obj2[i];
    }
    return temp;
}

/** Componentwise Componentwise multiplication between two AlgebraicVectors.
    @param[in] obj1 The first AlgebraicVector.
    @param[in] obj2 The second AlgebraicVector.
    @return A new AlgebraicVector.
    @throw std::runtime_error
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator*(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator*: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] * obj2[i];
    }
    return temp;
}

/** Componentwise multiplication between a AlgebraicVector and a scalar value.
    @param[in] obj1 A AlgebraicVector.
    @param[in] obj2 A scalar value.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator*(const AlgebraicVector<T> &obj1, const U &obj2) {
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] * obj2;
    }
    return temp;
}

/** Componentwise multiplication between a scalar value and a AlgebraicVector.
    @param[in] obj1 A scalar value.
    @param[in] obj2 A AlgebraicVector.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator*(const T &obj1, const AlgebraicVector<U> &obj2) {
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 * obj2[i];
    }
    return temp;
}

/** Componentwise division between two AlgebraicVectors.
    @param[in] obj1 The first AlgebraicVector.
    @param[in] obj2 The second AlgebraicVector.
    @return A new AlgebraicVector.
    @throw std::runtime_error
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator/(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    if(obj1.size() != obj2.size())
        throw std::runtime_error("DACE::AlgebraicVector<T>::operator/: Vectors must have the same length.");

    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] / obj2[i];
    }
    return temp;
}

/** Componentwise division between a AlgebraicVector and a scalar value.
    @param[in] obj1 A AlgebraicVector.
    @param[in] obj2 A scalar value.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator/(const AlgebraicVector<T> &obj1, const U &obj2) {
    const size_t size = obj1.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1[i] / obj2;
    }
    return temp;
}

/** Componentwise division between a scalar value and a AlgebraicVector.
    @param[in] obj1 A scalar value..
    @param[in] obj2 A AlgebraicVector..
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T, U>::returnType> operator/(const T &obj1, const AlgebraicVector<U> &obj2) {
    const size_t size = obj2.size();
    AlgebraicVector<typename PromotionTrait<T, U>::returnType> temp(size);
    for(size_t i=0; i<size; i++) {
        temp[i] = obj1 / obj2[i];
    }
    return temp;
}

/***********************************************************************************
*     Basic Arithmetic Functions
************************************************************************************/
/** Compute the derivative of a AlgebraicVector<T> with respect to variable p.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The independent variable number with respect to which the derivative is calculated.
    @return A new AlgebraicVector<T>.

    @see AlgebraicVector<T>::deriv
 */
template<typename T> AlgebraicVector<T> deriv(const AlgebraicVector<T> &obj, const unsigned int p) {
    return obj.deriv(p);
}

/** Compute the integral of a AlgebraicVector<T> with respect to variable p.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The independent variable number with respect to which the integral is calculated.
    @return A new AlgebraicVector<T> containing the result of the integration.
    @see AlgebraicVector<T>::integ
 */
template<typename T> AlgebraicVector<T> integ(const AlgebraicVector<T> &obj, const unsigned int p) {
    return obj.integ(p);
}

/***********************************************************************************
*     Filtering Functions
************************************************************************************/
/** Returns an AlgebraicVector<T> with all monomials of order less than @e min and greater
    than @e max removed (trimmed).
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj The AlgebraicVector<T> to be trimmed.
    @param[in] min The minimum order to be preserved.
    @param[in] max The maximum order to be preserved.
    @return A new AlgebraicVector<T>.
    @see AlgebraicVector<T>::trim
*/
template<typename T> AlgebraicVector<T> trim(const AlgebraicVector<T> &obj, unsigned int min, unsigned int max) {
    return obj.trim(min, max);
}

/***********************************************************************************
*     Intrinsic Functions
************************************************************************************/
/** Componentwise application of the absolute value function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::absolute
 */
template<typename T> AlgebraicVector<T> absolute(const AlgebraicVector<T> &obj) {
    return obj.absolute();
}

/** Componentwise application of the truncation function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::trunc
 */
template<typename T> AlgebraicVector<T> trunc(const AlgebraicVector<T> &obj) {
    return obj.trunc();
}

/** Componentwise application of the round function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::round
 */
template<typename T> AlgebraicVector<T> round(const AlgebraicVector<T> &obj) {
    return obj.round();
}

/** Componentwise application of the modulo function.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The divisor.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::mod
 */
template<typename T, typename U> AlgebraicVector<T> mod(const AlgebraicVector<T> &obj, const U &p) {
    return obj.mod(p);
}

/** Componentwise application of the integer power function.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The power to raise each element to.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::pow
 */
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const int p) {
    return obj.pow(p);
}

/** Componentwise application of the double power function.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The power to raise each element to.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::pow
 */
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const double p) {
    return obj.pow(p);
}

/** Componentwise application of the root function.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] p The root to be computed.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::root
 */
template<typename T> AlgebraicVector<T> root(const AlgebraicVector<T> &obj, const int p) {
    return obj.root(p);
}

/** Componentwise application of the multiplicative inverse function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::minv
 */
template<typename T> AlgebraicVector<T> minv(const AlgebraicVector<T> &obj) {
    return obj.minv();
}

/** Componentwise application of the square function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::sqr
 */
template<typename T> AlgebraicVector<T> sqr(const AlgebraicVector<T> &obj) {
    return obj.sqr();
}

/** Componentwise application of the square root function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::sqrt
 */
template<typename T> AlgebraicVector<T> sqrt(const AlgebraicVector<T> &obj) {
    return obj.sqrt();
}

/** Componentwise application of the inverse square root function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::isrt
 */
template<typename T> AlgebraicVector<T> isrt(const AlgebraicVector<T> &obj) {
    return obj.isrt();
}

/** Componentwise application of the cube root function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::cbrt
 */
template<typename T> AlgebraicVector<T> cbrt(const AlgebraicVector<T> &obj) {
    return obj.cbrt();
}

/** Componentwise application of the inverse cube root function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::icbrt
 */
template<typename T> AlgebraicVector<T> icbrt(const AlgebraicVector<T> &obj) {
    return obj.icbrt();
}

/** Componentwise application of the hypotenuse function.
    @param[in] X The AlgebraicVector<T> containing X.
    @param[in] Y The AlgebraicVector<T> containing Y.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::hypot
 */
template<typename T> AlgebraicVector<T> hypot(const AlgebraicVector<T> &X, const AlgebraicVector<T> &Y) {
    return X.hypot(Y);
}

/** Componentwise application of the exponential function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::exp
 */
template<typename T> AlgebraicVector<T> exp(const AlgebraicVector<T> &obj) {
    return obj.exp();
}

/** Componentwise application of the natural logarithm function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::log
 */
template<typename T> AlgebraicVector<T> log(const AlgebraicVector<T> &obj) {
    return obj.log();
}

/** Componentwise application of the logarithm function relative to given base.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] b The base for the logarithm.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::logb
 */
template<typename T> AlgebraicVector<T> logb(const AlgebraicVector<T> &obj, const double b) {
    return obj.logb(b);
}

/** Componentwise application of the decadic logarithm function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::log10
 */
template<typename T> AlgebraicVector<T> log10(const AlgebraicVector<T> &obj) {
    return obj.log10();
}

/** Componentwise application of the binary logarithm function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::log2
 */
template<typename T> AlgebraicVector<T> log2(const AlgebraicVector<T> &obj) {
    return obj.log2();
}

/** Componentwise application of the sine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::sin
 */
template<typename T> AlgebraicVector<T> sin(const AlgebraicVector<T> &obj) {
    return obj.sin();
}

/** Componentwise application of the cosine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::cos
 */
template<typename T> AlgebraicVector<T> cos(const AlgebraicVector<T> &obj) {
    return obj.cos();
}

/** Componentwise application of the tangent function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::tan
 */
template<typename T> AlgebraicVector<T> tan(const AlgebraicVector<T> &obj) {
    return obj.tan();
}

/** Componentwise application of the arcsine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::asin
 */
template<typename T> AlgebraicVector<T> asin(const AlgebraicVector<T> &obj) {
    return obj.asin();
}

/** Componentwise application of the arccosine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::acos
 */
template<typename T> AlgebraicVector<T> acos(const AlgebraicVector<T> &obj) {
    return obj.acos();
}

/** Componentwise application of the arctangent function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::atan
 */
template<typename T> AlgebraicVector<T> atan(const AlgebraicVector<T> &obj) {
    return obj.atan();
}

/** Componentwise application of the four-quadrant arctangent of @p Y/X.
    @param[in] Y The AlgebraicVector<T> containing Y.
    @param[in] X The AlgebraicVector<T> containing X.
    @return A new AlgebraicVector with elements in @f$ [-pi, pi] @f$ .
    @see AlgebraicVector<T>::atan2
 */
template<typename T> AlgebraicVector<T> atan2(const AlgebraicVector<T> &Y, const AlgebraicVector<T> &X) {
    return Y.atan2(X);
}

/** Componentwise application of the hyperbolic sine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::sinh
 */
template<typename T> AlgebraicVector<T> sinh(const AlgebraicVector<T> &obj) {
    return obj.sinh();
}

/** Componentwise application of the hyperbolic cosine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::cosh
 */
template<typename T> AlgebraicVector<T> cosh(const AlgebraicVector<T> &obj) {
    return obj.cosh();
}

/** Componentwise application of the hyperbolic tangent function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::tanh
 */
template<typename T> AlgebraicVector<T> tanh(const AlgebraicVector<T> &obj) {
    return obj.tanh();
}

/** Componentwise application of the hyperbolic arcsine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::asinh
 */
template<typename T> AlgebraicVector<T> asinh(const AlgebraicVector<T> &obj) {
    return obj.asinh();
}

/** Componentwise application of the hyperbolic arccosine function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::acosh
 */
template<typename T> AlgebraicVector<T> acosh(const AlgebraicVector<T> &obj) {
    return obj.acosh();
}

/** Componentwise application of the hyperbolic arctangent function.
    @param[in] obj An AlgebraicVector<T>.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::atanh
 */
template<typename T> AlgebraicVector<T> atanh(const AlgebraicVector<T> &obj) {
    return obj.atanh();
}

/***********************************************************************************
*     Norm & Estimation Functions
************************************************************************************/
/** Componentwise application of the norm function.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] type The type of norm to be computed.
    @return A new AlgebraicVector.
    @see AlgebraicVector<T>::norm
    @see DA::norm
 */
template<typename T> AlgebraicVector<double> norm(const AlgebraicVector<T> &obj, const unsigned int type) {
    return obj.norm(type);
}

/***********************************************************************************
*     Evaluation Functions
************************************************************************************/
/** Evaluate each element of a vector of DA with any vector type @e U containing the arguments
    and return a vector of results of the same type @e U.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] args A vector (e.g. AlgebraicVector<> or std::vector<>) of arguments.
    @return A new vector of same type as argument args containing the.
    results of the evaluation
    @see AlgebraicVector<T>::eval()
 */
template<typename T, typename U> U eval(const AlgebraicVector<T> &obj, const U &args) {
    return obj.eval(args);
}

/** Evaluate each element of a vector of DA with a braced initializer list of type @e U
    and return an AlgebraicVector<U> with the results.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] l A braced initializer list containing the arguments.
    @return A new AlgebraicVector<U> containing the results of the evaluation.
    @see AlgebraicVector<T>::eval()
 */
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const std::initializer_list<U> l) {
    return obj.eval(l);
}

/** Evaluate each element of a vector of DA with an array of arithmetic type @e T arguments.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @param[in] obj An AlgebraicVector<T>.
    @param[in] args A C array of arithmetic type T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array args.
    @return A new AlgebraicVector<U> containing the results of the evaluation.
    @see AlgebraicVector<T>::eval()
 */
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const U args[], const unsigned int length) {
    return obj.eval(args, length);
}

/** Evaluate each element of a vector of DA with a single arithmetic type @e U argument.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @note For efficient repeated evaluation of the same AlgebraicVector use the corresponding method
    in class DACE::compiledDA.
    @param[in] obj The AlgebraicVector<T> to evaluate.
    @param[in] arg The argument of type @e T.
    @return The result of the evaluation.
    @see AlgebraicVector<T>::evalScalar()
 */
template<typename T, typename U> AlgebraicVector<U> evalScalar(const AlgebraicVector<T> &obj, const U &arg) {
    return obj.evalScalar(arg);
}

/** Compile vector of polynomials and create a compiledDA object.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj The AlgebraicVector to compile.
    @return The compiled DA object.
    @see AlgebraicVector<T>::compile()
 */
template<typename T> compiledDA compile(const AlgebraicVector<T> &obj) {
    return obj.compile();
}

/** Partial evaluation of vector of polynomials. In each element of the vector,
    variable @e var is replaced by the value @e val. The resulting vector of DAs
    is returned.
    @warning This function only works on @p AlgebraicVector<DA>. When called on other
    data types (e.g. double) a compiler error is issued.
    @param[in] obj The vector to partially evaluate.
    @param[in] var The independent variable number to be replaced.
    @param[in] val The value by which to replace the variable.
    @return A new DA object.
    @see AlgebraicVector<T>::plug()
 */
template<typename T> AlgebraicVector<T> plug(const AlgebraicVector<T> &obj, const unsigned int var, const double val) {
    return obj.plug(var, val);
}

/***********************************************************************************
*     Input/Output Functions
************************************************************************************/
/** Output a vector to a C++ output stream.
    @param[in] out A C++ output stream.
    @param[in] obj An AlgebraicVector to be written to the stream.
    @return Reference to output stream @e out.
 */
template<typename T> std::ostream& operator<<(std::ostream &out, const AlgebraicVector<T> &obj) {
    const size_t size = obj.size();

    out << "[[[ " << size << " vector" << std::endl;
    for(size_t i=0; i<size;i++) {
        out << obj[i] << std::endl;
    }
    out << "]]]" << std::endl;

    return out;
}

/** Read a vector from a C++ input stream.
    @param[in] in A C++ input stream.
    @param[in] obj An AlgebraicVector to be read from the stream.
    @return Reference to input stream @e in.
 */
template<typename T> std::istream& operator>>(std::istream &in, AlgebraicVector<T> &obj) {
    std::string init_line;
    size_t vec_size;

    // try to read the first line
    getline(in, init_line);
    if(in.good()) {
        // retrieve the size of the vector to be read
        std::size_t found = init_line.find_first_of(' ');
        std::string size_str(init_line,4,found-4);
        if(!(std::istringstream(size_str) >> vec_size)) vec_size = 0;

        // resize the object to meet the size of the vector to be read
        obj.resize(vec_size);

        // fill the AlgebraicVector
        for (size_t i = 0; in.good() && (i < vec_size); i++) {
            in >> obj[i];

            // check the next character
            if (in.peek() == '\n')       // the previous operator>> does not consume the \n character when an AlgebraicVector<T> (with T != DA) is considered
                in.ignore();            // ignore the next character
        }
        // skip the line at the end of a AlgebraicVector (containing ]]])
        getline(in, init_line);
    } else {
        obj.resize(0);
    }

    return in;
}

/***********************************************************************************
*     Vector Functions
************************************************************************************/
/** Compute the dot product between two AlgebraicVectors.
   @param[in] obj1 An AlgebraicVector.
   @param[in] obj2 An AlgebraicVector.
   @return A scalar value.
 */
template<typename T, typename U> typename PromotionTrait<T,U>::returnType dot(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    return obj1.dot(obj2);
}

/** Compute the cross product between two 3D AlgebraicVectors.
    @param[in] obj1 An AlgebraicVector.
    @param[in] obj2 An AlgebraicVector.
    @return A new AlgebraicVector.
 */
template<typename T, typename U> AlgebraicVector<typename PromotionTrait<T,U>::returnType> cross(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2) {
    return obj1.cross(obj2);
}

/** Compute vector length (Euclidean vector norm).
    @param[in] obj An AlgebraicVector<T>.
    @return Euclidean length of the vector.
 */
template<typename T> T length(const AlgebraicVector<T> &obj) {
    return obj.length();
}

/** Normalize an AlgebraicVector<T>.
    @param[in] obj An AlgebraicVector<T> to normalize.
    @return An AlgebraicVector<T> of unit length.
    @see AlgebraicVector<T>::normalize
 */
template<typename T> AlgebraicVector<T> normalize(const AlgebraicVector<T> &obj) {
    return obj.normalize();
}

}

#endif /* DINAMICA_ALGEBRAICVECTOR_T_H_ */
