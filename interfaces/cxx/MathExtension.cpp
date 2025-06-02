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
#include "dace/dacecore.h"

namespace DACE {

/** Absolute value of @e x.
    @param[in] x The function argument.
 */
double absolute(const double x) {
    return std::abs(x);
}

/** Constant part. For double type this is just @e x.
    @param[in] x The function argument.
 */
double cons(const double x) {
    return x;
}

/** Logarithm relative to base @e b.
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

/** Square of @e x.
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

/** p-th root of @e x.
    @param[in] x The function argument.
    @param[in] p The root to take.
 */
double root(const double x, const int p) {
    return std::pow(x, 1.0/p);
}

/** norm of @e x.
    @param[in] x The function argument.
    @param[in] type The type of norm (ignored for double).
 */
double norm(const double x, const int type) {
    return std::abs(x);
}

/** Bessel J function of @e x.
    @param[in] n The order of the Bessel function.
    @param[in] x The function argument.
 */
double BesselJFunction(const int n, const double x) {
    double res;
    BesselWrapper(x, n, n, -1, &res);
    return res;
}

/** Bessel Y function of @e x.
    @param[in] n The order of the Bessel function.
    @param[in] x The function argument.
 */
double BesselYFunction(const int n, const double x) {
    double res;
    BesselWrapper(x, n, n, 1, &res);
    return res;
}

/** Bessel I function of @e x.
    @param[in] n The order of the Bessel function.
    @param[in] x The function argument.
    @param[in] scaled If true the result is scaled by `exp(x)`.
 */
double BesselIFunction(const int n, const double x, const bool scaled) {
    double res;
    ModifiedBesselWrapper(x, n, n, scaled ? -2 : -1, &res);
    return res;
}

/** Bessel K function of @e x.
    @param[in] n The order of the Bessel function.
    @param[in] x The function argument.
    @param[in] scaled If true the result is scaled by `exp(x)`.
 */
double BesselKFunction(const int n, const double x, const bool scaled) {
    double res;
    ModifiedBesselWrapper(x, n, n, scaled ? 2 : 1, &res);
    return res;
}

/** Gamma function of @e x.
    @param[in] x The function argument.
 */
double GammaFunction(const double x) {
    return std::tgamma(x);
}

/** Log gamma function of @e x.
    @param[in] x The function argument.
 */
double LogGammaFunction(const double x) {
    return std::lgamma(x);
}

/** Legendre polynomial of degree @e n of @e x.
    @param[in] n The degree of the Legendre polynomial.
    @param[in] x The function argument.
 */
double LegendrePolynomial(const unsigned int n, const double x) {
//    return std::legendre(n, x);   // not yet widely available
    if(n == 0) return 1.0;
    if(n == 1) return x;

    double P[3] = {1.0, x, 0.0};
    for(unsigned int i = 2; i <= n; i++)
        P[i%3] = ((2*i-1)*x*P[(i-1)%3] - (i-1)*P[(i-2)%3])/i;
    return P[n%3];
}

/** Hermite polynomial of degree @e n of @e x.
    @param[in] n The degree of the Hermite polynomial.
    @param[in] x The function argument.
 */
double HermitePolynomial(const unsigned int n, const double x) {
//    return std::hermite(n, x);   // not yet widely available
    if(n == 0) return 1.0;
    if(n == 1) return 2.0*x;

    double P[3] = {1.0, 2.0*x, 0.0};
    for(unsigned int i = 2; i <= n; i++)
        P[i%3] = 2.0*x*P[(i-1)%3] - 2.0*(i-1)*P[(i-2)%3];
    return P[n%3];
}

/** Laguerre polynomial of degree @e n of @e x.
    @param[in] n The degree of the Laguerre polynomial.
    @param[in] x The function argument.
 */
double LaguerrePolynomial(const unsigned int n, const double x) {
//    return std::laguerre(n, x);   // not yet widely available
    if(n == 0) return 1.0;
    if(n == 1) return 1.0-x;

    double P[3] = {1.0, 1.0-x, 0.0};
    for(unsigned int i = 2; i <= n; i++)
        P[i%3] = ((2*i-1-x)*P[(i-1)%3] - (i-1)*P[(i-2)%3])/i;
    return P[n%3];
}

}
