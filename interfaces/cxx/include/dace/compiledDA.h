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

/*! \file

    \brief DA representation for efficient repeated evaluation.

    This header file contains the compiledDA class for a representation of one
    or more DA objects that is prepared ("compiled") for efficient repeated evaluation.
*/

#ifndef DINAMICA_COMPILEDDA_H_
#define DINAMICA_COMPILEDDA_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <initializer_list>

namespace DACE {

// forward declaration
class DA;

/*! Representation of a precomputed DA polynomial for efficient evaluation.
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
    *     Evaluation
    *********************************************************************************/
    template<class V> V eval(const V &args) const;
    template<class T> std::vector<T> eval(const std::initializer_list<T> l) const;
    template<class T> std::vector<T> eval(const T args[], const unsigned int length) const;
    template<class T> std::vector<T> evalScalar(const T &arg) const;
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
template<> DACE_API void compiledDA::eval(const std::vector<DA> &args, std::vector<DA> &res) const;
template<> DACE_API void compiledDA::eval(const std::vector<double> &args, std::vector<double> &res) const;

}

#endif /* DINAMICA_COMPILEDDA_H_ */
