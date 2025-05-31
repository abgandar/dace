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
 * compiledDA.cpp
 *
 *  Created on: Mar 01, 2014
 *      Author: Dinamica Srl
 */

// DACE classes
#include "dace/config.h"
#include "dace/compiledDA.h"
#include "dace/DA.h"
#include "dace/DACEException.h"

namespace DACE {

/********************************************************************************
*     Constructors & Destructors
*********************************************************************************/
/** Create a copy of a compiledDA object.
    @param[in] cda The compiled DA object to be copied.
 */
compiledDA::compiledDA(const compiledDA &cda) {
    dim = cda.dim;
    terms = cda.terms;
    ord = cda.ord;
    ac = new double[terms*(dim+2)];
    for(int i = terms*(dim+2)-1; i >= 0; i--) ac[i] = cda.ac[i];
}

/** Create a vector of compiledDA objects from a vector of DA objects.
    @param[in] da A vector of DA objects.
    @throw DACE::DACEException
 */
compiledDA::compiledDA(const std::vector<DA> &da) {
    dim = (unsigned int)da.size();
    if(dim < 1) DACEException(16, 04);

    ac = new double[DA::getMaxMonomials()*(dim+2)];
    const DACEDA **mb = new const DACEDA*[dim];
    for(unsigned int i = 0; i < dim; i++) mb[i] = &(da[i].m_index);
    daceEvalTree(mb, dim, ac, terms, ord);
    delete[] mb;
    if(daceGetError()) DACEException();
}

/** Create a compiledDA object from a DA object.
    @param[in] da A DA object.
    @throw DACE::DACEException
 */
compiledDA::compiledDA(const DA &da) {
    ac = new double[DA::getMaxMonomials()*3];
    dim = 1;
    const DACEDA *mb = &(da.m_index);
    daceEvalTree(&mb, dim, ac, terms, ord);
    if(daceGetError()) DACEException();
}

/** Destructor.
 */
compiledDA::~compiledDA() throw() {
    delete[] ac;
}

/********************************************************************************
*     Assignments
*********************************************************************************/
/** Copy the content of a given compiledDA object into the current
    compiledDA.
    @param[in] cda The compiledDA vector to be copied.
    @return The compiledDA object with the same content of the given compiledDA.
 */
compiledDA& compiledDA::operator=(const compiledDA &cda) {
    if(this != &cda)
    {
        dim = cda.dim;
        terms = cda.terms;
        ord = cda.ord;
        if(ac) delete[] ac;
        ac = new double[terms*(dim+2)];
        for(int i = terms*(dim+2)-1; i >= 0; i--) ac[i] = cda.ac[i];
    }

    return *this;
}

/********************************************************************************
*     Evaluation operator template specializations
*********************************************************************************/
// double evaluation
template<> void compiledDA::operator()(const std::vector<double> &args, std::vector<double> &res) const {
    res.reserve(dim);
    daceEvalTreeDouble(res.data(), dim, args.data(), args.size(), ac, terms, ord);
}

// DA evaluation
template<> void compiledDA::operator()(const std::vector<DA> &args, std::vector<DA> &res) const {
    DACEDA **r = new DACEDA*[dim];
    const DACEDA **a = new const DACEDA*[args.size()];
    res.reserve(dim);
    for(unsigned int i = 0; i < dim; i++) r[i] = &res[i].m_index;
    for(unsigned int i = 0; i < args.size(); i++) a[i] = &args[i].m_index;
    daceEvalTreeDA(r, dim, a, args.size(), ac, terms, ord);
    delete[] r;
    delete[] a;
}

/********************************************************************************
*     Member access routines
*********************************************************************************/
/** Return the coefficient array.
    @return The coefficient array.
 */
const double* compiledDA::getAc() const {
    return this->ac;
}

/** Return the number of DAs (dimension).
    @return The dimension.
 */
unsigned int compiledDA::getDim() const {
    return this->dim;
}

/** Return the maximum order.
    @return The maximum order.
 */
unsigned int compiledDA::getOrd() const {
    return this->ord;
}

/** Return the number of terms.
    @return the number of terms.
 */
unsigned int compiledDA::getTerms() const {
    return this->terms;
}

}
