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
/** Generic evaluation of the DA with a vector-like type @e T of arguments.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] args A vector-like object of arithemtic type (e.g. std::vector<double>)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see compiledDA::operator()()
 */
template<class T> typename T::value_type DA::operator()(const T &args) const {
    return compiledDA(*this)(args)[0];
}

/** Generic evaluation of the DA with arguments in a braced initializer list of type @e T.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] l A braced initializer list of arithemtic type (e.g. `{1.0, 2.0, 3.0}`)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see compiledDA::operator()()
 */
template<class T> T DA::operator()(const std::initializer_list<T> l) const {
    return compiledDA(*this)(l)[0];
}

/** Generic evaluation of the DA with an array of arithmetic type @e T arguments.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] args An array of arithmetic type @e T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array args.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see compiledDA::operator()()
 */
template<class T> T DA::operator()(const T args[], const unsigned int length) const {
    return compiledDA(*this)(args, length)[0];
}

/** Generic evaluation of the DA with a vector-like type @e T of arguments.
    @deprecated Replaced by DA::operator()().
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] args A vector-like object of arithemtic type (e.g. std::vector<double>)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
    @see compiledDA::operator()()
 */
template<class T> typename T::value_type DA::eval(const T &args) const {
    return (*this)(args);
}

/** Generic evaluation of the DA with arguments in a braced initializer list of type @e T.
    @deprecated Replaced by DA::operator()().
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] l A braced initializer list of arithemtic type (e.g. `{1.0, 2.0, 3.0}`)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
    @see compiledDA::operator()()
 */
template<class T> T DA::eval(const std::initializer_list<T> l) const {
    return (*this)(l);
}

/** Generic evaluation of the DA with an array of arithmetic type @e T arguments.
    @deprecated Replaced by DA::operator()().
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] args A C array of arithmetic type @e T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array @e args.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
    @see compiledDA::operator()()
 */
template<class T> T DA::eval(const T args[], const unsigned int length) const {
    return (*this)(args, length);
}

/** Generic evaluation of the DA with a single arithmetic type @e T argument.
    @deprecated Replaced by DA::operator()() with braced initializer list (e.g. `da({arg})`).
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] arg A single variable of arithmetic type @e T of the first independent DA variable.
    @return The result of the evaluation.
    @see compiledDA::operator()()
 */
template<class T> T DA::evalScalar(const T &arg) const {
    return compiledDA(*this).evalScalar(arg)[0];
}

/** Generic evaluation of the DA with a vector-like type @e T of arguments.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] da A DA object.
    @param[in] args A vector-like object of arithemtic type (e.g. std::vector<double>)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see compiledDA::operator()()
 */
template<class T> typename T::value_type eval(const DA &da, const T &args) {
    return da(args);
}

/** Generic evaluation of the DA with arguments in a braced initializer list of type @e T.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] da A DA object.
    @param[in] l A braced initializer list of arithemtic type (e.g. `{1.0, 2.0, 3.0}`)
    with which the DA vector is evaluated.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
    @see compiledDA::operator()()
 */
template<class T> T eval(const DA &da, const std::initializer_list<T> l) {
    return da(l);
}

/** Generic evaluation of the DA with an array of arithmetic type @e T arguments.
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] da A DA object.
    @param[in] args An array of arithmetic type @e T with which the DA vector is evaluated.
    @param[in] length The number of elements in the array @e args.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
    @see compiledDA::operator()()
 */
template<class T> T eval(const DA &da, const T args[], const unsigned int length) {
    return da(args, length);
}

/** Generic evaluation of the DA with a single arithmetic type @e T arguments.
    @deprecated Replaced by eval() with braced initializer list (e.g. `eval(da, {arg})`).
    @note To be used only for single polynomial evaluation. For multiple
    evaluations of the same polynomial use the corresponding method in class
    DACE::compiledDA.
    @param[in] da A DA object.
    @param[in] arg A single variable of arithmetic type @e T of the first independent DA variable.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::operator()()
 */
template<class T> T evalScalar(const DA &da, const T &arg) {
    return da.evalScalar(arg);
}

}

#endif /* DINAMICA_DA_T_H_ */
