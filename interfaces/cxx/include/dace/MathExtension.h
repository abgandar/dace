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
 * MathExtension.h
 *
 *  Created on: Sep. 22, 2014
 *      Author: Dinamica Srl
 */

/*  DACE extensions to standard math routines provided by the standard library

    This file provides double overloads of math functions provided by DACE for DA and not
    available in the standard library. These overloads provide sensible double versions,
    allowing using the same code for both DA and double data types (e.g. via templates).
 */

#ifndef DINAMICA_MATHEXTENSION_H_
#define DINAMICA_MATHEXTENSION_H_

namespace DACE {

/** @name Double Intrinsic Functions
 * @{
 */
DACE_API double absolute(const double x);
DACE_API double cons(const double x);
DACE_API double logb(const double x, const double b = 10.0);
DACE_API double isrt(const double x);
DACE_API double icbrt(const double x);
DACE_API double sqr(const double x);
DACE_API double minv(const double x);
DACE_API double mod(const double x, const double p);
DACE_API double root(const double x, const int p = 2);
DACE_API double norm(const double x, const int type = 0);
DACE_API double BesselJFunction(const int n, const double x);
DACE_API double BesselYFunction(const int n, const double x);
DACE_API double BesselIFunction(const int n, const double x, const bool scaled = false);
DACE_API double BesselKFunction(const int n, const double x, const bool scaled = false);
DACE_API double GammaFunction(const double x);
DACE_API double LogGammaFunction(const double x);
DACE_API double LegendrePolynomial(const unsigned int n, const double x);
DACE_API double HermitePolynomial(const unsigned int n, const double x);
DACE_API double LaguerrePolynomial(const unsigned int n, const double x);
/** @} */

}

#endif /* DINAMICA_MATHEXTENSION_H_ */

