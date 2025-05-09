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
 * compiledDA_t.h
 *
 *  Created on: Apr 07, 2014
 *      Author: Dinamica Srl
 */

/*  Templated function definitions for compiledDA class.

    This header file contains the definition of templated functions in the compiledDA class.
*/

#ifndef DINAMICA_COMPILEDDA_T_H_
#define DINAMICA_COMPILEDDA_T_H_

// DACE classes
#include "dace/compiledDA.h"

/** @addtogroup DACECXX C++ Interface
 *  @{
 */

namespace DACE {

/********************************************************************************
*     compiledDA evaluation routines
*********************************************************************************/
/** Evaluate the compiled polynomial with a vector of any arithmetic type
    (such as DA or double) and return vector of results.
    @param[in] args the values of the independent DA variables to evaluate
    with. Must be a std::vector<> (or derived class) of an arithmetic
    type. If less than the number of independent DA variables defined
    during the DACE initialization are given, the missing entries are
    assumed to be zero.
    @return Vector with the result of the evaluation. The vector is of
    the same type as the argument args.
 */
template<class V> V compiledDA::eval(const V &args) const {
    V res(dim);
    eval(args, res);

    return res;
}

/** Evaluate the compiled polynomial with a braced initializer list of any arithmetic type
    (such as DA or double) and return vector of results.
    @param[in] l the values of the independent DA variables to evaluate
    with. Must be a braced initializer list of an arithmetic
    type. If less than the number of independent DA variables defined
    during the DACE initialization are given, the missing entries are
    assumed to be zero.
    @return std::vector with the result of the evaluation
    @note C++ is not able to derive the type of elements of an initializer list automatically.
    That means eval() must be called explicitly as e.g. eval<double>({1.0, 2.0, 3.0}) when
    used with initializer lists.
 */
template<class T> std::vector<T> compiledDA::eval(const std::initializer_list<T> l) const {
    std::vector<T> res(dim);
    eval(std::vector<T>(l), res);

    return res;
}

/** Evaluate the compiled polynomial with an array of any arithmetic type
    (such as DA or double) and return vector of results.
    @param[in] args array of the values of the independent DA variables to
    evaluate with
    @param[in] length Size of the array args[]. If less than the number of
    variables defined during the DACE initialization are given, the
    missing entries are assumed to be zero.
    @return Vector with the result of the evaluation. The vector is of
    type std::vector<V>.
 */
template<class T> std::vector<T> compiledDA::eval(const T args[], const unsigned int length) const {
    std::vector<T> arg(args,args+length);
    std::vector<T> res(dim);
    eval(arg, res);

    return res;
}

/** Evaluate the compiled polynomial with a single argument of any
    arithmetic type (such as DA or double) and return vector of results.
    @param[in] arg The value of the first independent DA variable to evaluate
    with. All remaining independent DA variables are assumed to be zero.
    @return Vector with the result of the evaluation. The vector is of
    type std::vector<V>.
 */
template<class T> std::vector<T> compiledDA::evalScalar(const T &arg) const {
    std::vector<T> args(1);
    std::vector<T> res(dim);
    args[0] = arg;
    eval(args, res);

    return res;
}

/** Evaluate the compiled polynomial with a vector of arithmetic type T
    (such as DA or double) and return the result in the vector res.
    @param[in] args the values of the independent DA variables to evaluate
    with. Must be a std::vector<> (or derived class) of an arithmetic
    type. If less than the number of independent DA variables defined
    during the DACE initialization are given, the missing entries are
    assumed to be zero.
    @param[out] res the vector containing the result of the evaluation.
    Must be at least of length dim.
    @tparam T Any arithmetic type. Must support at least:\n
     - default constructor    T::T()\n
     - assignment             T::operator=(const T& t)\n
     - assignment & addition  T::operator+=(const T& t)\n
     - double multiplication  T::operator*(const double d)\n
     - double addition        T::operator+(const double d)
 */
template<class T> void compiledDA::eval(const std::vector<T> &args, std::vector<T> &res) const {
    const unsigned int narg = args.size();
    unsigned int jlskip = ord+1;
    double *p = ac+2;
    T *xm = new T[ord+1];

    // prepare temporary powers
    xm[0] = args[0]*0.0 + 1.0;
    // constant part
    for(unsigned int i=0; i<dim; i++, p++)
        res[i] = args[0]*0.0 + (*p);
    // higher order terms
    for(unsigned int i=1; i<terms; i++){
        unsigned int jl = (unsigned int)(*p); p++;
        unsigned int jv = (unsigned int)(*p)-1; p++;
        if(jl > jlskip)
        {
            p += dim;
            continue;
        }
        if(jv >= narg)
        {
            jlskip = jl;
            p += dim;
            continue;
        }
        jlskip = ord+1;
        xm[jl] = xm[jl-1]*args[jv];
        for(unsigned int j=0; j<dim; j++, p++)
            if((*p)!=0.0) res[j] += xm[jl]*(*p);
    }

    delete[] xm;
}

}

#endif /* DINAMICA_COMPILEDDA_T_H_ */

/** @}*/
