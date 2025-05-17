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
 * compiledDA.h
 *
 *  Created on: Mar 01, 2014
 *      Author: Dinamica Srl
 */

/*  DA representation for efficient evaluation.

    This header file contains the compiledDA class for a representation of one
    or more DA objects that is prepared ("compiled") for efficient evaluation.
*/

#ifndef DINAMICA_COMPILEDDA_H_
#define DINAMICA_COMPILEDDA_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <initializer_list>

namespace DACE {

// forward declaration
class DA;

/** Representation of a precomputed DA polynomial for efficient evaluation.
    @ingroup DACECXX

    Use this class to convert a DACE::DA or vector of DACE::DA into an optimized
    representation for evaluation if you need to evaluate the polynomial repeatedly
    with different arguments. That way the compilation step is only carried out once.

    Example:
    @code
        #include <iostream>
        #include <chrono>
        #include <dace/dace.h>

        std::vector<double> v({0.1, 0.2});
        double arr[2] = {0.1, 0.2};

        DA::init(10, 2);

        // direct DA evaluations (all identical result)
        DA x = DA::random(0.7);
        std::cout << x({0.1, 0.2}) << std::endl;
        std::cout << x(v) << std::endl;
        std::cout << x(arr, 2) << std::endl;

        // compiledDA evaluation (all identical result)
        compiledDA cda(x);
        std::cout << cda({0.1, 0.2})[0] << std::endl;
        std::cout << cda(v)[0] << std::endl;
        std::cout << cda(arr, 2)[0] << std::endl;

        // compiledDA vector evaluation
        AlgebraicVector<DA> Y(3);
        Y[0] = DA::random(0.7);
        Y[1] = DA::random(0.7);
        Y[2] = DA::random(0.7);
        compiledDA cY(Y);
        std::vector<double> res = cY({0.1, 0.2});          // evaluate all 3 components at once
        std::cout << res[0] << "  " << res[1] << "  " << res[2] << std::endl;
        std::cout << Y[0](v) << "  " << Y[1](v) << "  " << Y[2](v) << std::endl;    // identical manual evaluation

        // performance comparison
        constexpr unsigned int N = 10000;
        std::vector<double> r(N);
        auto t0 = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i < N; i++)
        {
            v[0] = 0.1 + (double)i/N;
            r[i] = x(v);
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time DA (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << std::endl;

        t0 = std::chrono::high_resolution_clock::now();
        for(unsigned int i = 0; i < N; i++)
        {
            v[0] = 0.1 + (double)i/N;
            r[i] = cda(v)[0];
        }
        t1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time compiledDA (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << std::endl;
    @endcode
 */
class DACE_API compiledDA
{
private:
    double *ac;             //!< Compiled polynomial evaluation data
    unsigned int dim;       //!< Number of polynomials (dimension)
    unsigned int ord;       //!< Maximum order of the polynomial
    unsigned int vars;      //!< Number of variables in the polynomial
    unsigned int terms;     //!< Number of terms in the polynomial

public:
    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    compiledDA(const compiledDA &cda);
    compiledDA(const DA &da);
    compiledDA(const std::vector<DA> &da);
    ~compiledDA() throw();

    /********************************************************************************
    *     Assignments
    *********************************************************************************/
    compiledDA& operator=(const compiledDA &cda);

    /********************************************************************************
    *     Evaluation operators
    *********************************************************************************/
    template<class V> V operator()(const V &args) const;
    template<class T> std::vector<T> operator()(const std::initializer_list<T> l) const;
    template<class T> std::vector<T> operator()(const T args[], const unsigned int length) const;
    template<class T> void operator()(const std::vector<T> &args, std::vector<T> &res) const;

    /********************************************************************************
    *     Evaluation
    *********************************************************************************/
    template<class T> std::vector<T> evalScalar(const T &arg) const;
    template<class V> V eval(const V &args) const;
    template<class T> std::vector<T> eval(const std::initializer_list<T> l) const;
    template<class T> std::vector<T> eval(const T args[], const unsigned int length) const;
    template<class T> void eval(const std::vector<T> &args, std::vector<T> &res) const;

    /********************************************************************************
    *     Member access routines
    *********************************************************************************/
    const double* getAc() const;
    unsigned int getDim() const;
    unsigned int getOrd() const;
    unsigned int getVars() const;
    unsigned int getTerms() const;
};

// specializations for particularly efficient evaluation with double and DA arguments implemented in the library
/** @{ */
template<> DACE_API void compiledDA::operator()(const std::vector<DA> &args, std::vector<DA> &res) const;
template<> DACE_API void compiledDA::operator()(const std::vector<double> &args, std::vector<double> &res) const;
/** @} */
}

#endif /* DINAMICA_COMPILEDDA_H_ */
