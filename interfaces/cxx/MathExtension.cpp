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
 * MathExtension.cpp
 *
 *  Created on: Sep. 22, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes
#include <cmath>

// DACE classes
#include "dace/config.h"
#include "dace/MathExtension.h"

namespace DACE {

/** Absolute value.
    @param[in] x The function argument.
 */
double absolute(const double x) {
    return std::abs(x);
}

/** Constant part. For double type this is just x.
    @param[in] x The function argument.
 */
double cons(const double x) {
    return x;
}

/** Logarithm relative to base b.
    @param[in] x The function argument.
    @param[in] b The base of the logarithm (must be positive).
 */
double logb(const double x, const double b) {
    return std::log(x)/std::log(b);
}

/** Inverse square root 1/sqrt(x).
    @param[in] x The function argument.
 */
double isrt(const double x) {
    return 1.0/std::sqrt(x);
}

/** Inverse cube root 1/cbrt(x).
    @param[in] x The function argument.
 */
double icbrt(const double x) {
    return 1.0/std::cbrt(x);
}

/** Square of x.
    @param[in] x The function argument.
 */
double sqr(const double x) {
    return x*x;
}

/** Multiplicative inverse 1/x.
    @param[in] x The function argument.
 */
double minv(const double x) {
    return 1.0/x;
}

/** Modulo function (remainder of x/p).
    @param[in] x The dividend.
    @param[in] p The divisor.
 */
double mod(const double x, const double p) {
    return std::fmod(x, p);
}

/** p-th root of x.
    @param[in] x The function argument.
    @param[in] p The root to take.
 */
double root(const double x, const int p) {
    return std::pow(x, 1.0/p);
}

/** norm of x.
    @param[in] x The function argument.
    @param[in] type the type of norm (ignored for double).
 */
double norm(const double x, const int type) {
    return std::abs(x);
}

}
