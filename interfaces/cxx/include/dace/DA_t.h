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
 * DA_t.h
 *
 *  Created on: Apr 07, 2014
 *      Author: Dinamica Srl
 */

/*  Templated function definitions for main DA class.

    This header file contains the definition of templated functions in the DA class.
*/

#ifndef DINAMICA_DA_T_H_
#define DINAMICA_DA_T_H_

// DACE classes
#include "dace/compiledDA.h"
#include "dace/DA.h"

namespace DACE {

/********************************************************************************
*     DA polynomial evaluation routines
*********************************************************************************/
/** Generic evaluation of the DA with a vector of arithmetic type T arguments.
    @param[in] args std::vector<T> of arithmetic type T with which the DA vector is evaluated
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @note This function can be called with a braced initializer list. However, C++
    is not able to derive the type of elements of an initializer list automatically.
    That means operator() must be called explicitly as e.g. <double>({1.0, 2.0, 3.0}) when
    used with initializer lists.
    @see compiledDA::operator()
 */
template<class T> T DA::operator()(const std::vector<T> &args) const {
    return compiledDA(*this)(args)[0];
}

/** Generic evaluation of the DA with an array of arithmetic type T arguments.
    @param[in] args array of arithmetic type T with which the DA vector is evaluated
    @param[in] length number of elements in the array args
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @see compiledDA::operator()
 */
template<class T> T DA::operator()(const T args[], const unsigned int length) const {
    return compiledDA(*this)(args, length)[0];
}

/** Generic evaluation of the DA with a vector of arithmetic type T arguments.
    @param[in] args std::vector<T> of arithmetic type T with which the DA vector is evaluated
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @note This function can be called with a braced initializer list. However, C++
    is not able to derive the type of elements of an initializer list automatically.
    That means eval() must be called explicitly as e.g. eval<double>({1.0, 2.0, 3.0}) when
    used with initializer lists.
    @deprecated Replaced by DA::operator().
    @see DA::operator()
    @see compiledDA::operator()
 */
template<class T> T DA::eval(const std::vector<T> &args) const {
    return (*this)(args);
}

/** Generic evaluation of the DA with an array of arithmetic type T arguments.
    @param[in] args array of arithmetic type T with which the DA vector is evaluated
    @param[in] length number of elements in the array args
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @deprecated Replaced by DA::operator().
    @see DA::operator()
    @see compiledDA::operator()
 */
template<class T> T DA::eval(const T args[], const unsigned int length) const {
    return (*this)(args, length);
}

/** Generic evaluation of the DA with a single arithmetic type T argument.
    @param[in] arg single variable of arithmetic type T of the first independent DA variable
    @return The result of the evaluation
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @see compiledDA::evalScalar
 */
template<class T> T DA::evalScalar(const T &arg) const {
    return compiledDA(*this).evalScalar(arg)[0];
}

/** Generic evaluation of the DA with a vector of arithmetic type T arguments.
    @param[in] da a DA object
    @param[in] args std::vector<T> of arithmetic type T with which the DA vector is evaluated
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @note This function can be called with a braced initializer list. However, C++
    is not able to derive the type of elements of an initializer list automatically.
    That means eval() must be called explicitly as e.g. eval<double>(x, {1.0, 2.0, 3.0}) when
    used with initializer lists.
    @see compiledDA::operator()
 */
template<class T> T eval(const DA &da, const std::vector<T> &args) {
    return da(args);
}

/** Generic evaluation of the DA with an array of arithmetic type T arguments.
    @param[in] da a DA object
    @param[in] args array of arithmetic type T with which the DA vector is evaluated
    @param[in] length number of elements in the array args
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @see compiledDA::operator()
 */
template<class T> T eval(const DA &da, const T args[], const unsigned int length) {
    return da(args, length);
}

/** Generic evaluation of the DA with a single arithmetic type T arguments.
    @param[in] da a DA object
    @param[in] arg single variable of arithmetic type T of the first independent DA variable
    @return The result of the evaluation
    @throw DACE::DACEException
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    compiledDA.
    @see compiledDA::evalScalar
 */
template<class T> T evalScalar(const DA &da, const T &arg) {
    return da.evalScalar(arg);
}

}

#endif /* DINAMICA_DA_T_H_ */
