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

/** @addtogroup DACECXX C++ Interface
    @{
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
    @param[in] cda compiled DA object to be copied
 */
compiledDA::compiledDA(const compiledDA &cda) {
    dim = cda.dim;
    terms = cda.terms;
    vars = cda.vars;
    ord = cda.ord;
    ac = new double[terms*(dim+2)];
    for(int i=terms*(dim+2)-1; i>=0; i--) ac[i] = cda.ac[i];
}

/** Create a vector of compiledDA objects from a vector of DA objects.
    @param[in] da vector of DA objects
    @throw DACE::DACEException
 */
compiledDA::compiledDA(const std::vector<DA> &da) {
    dim = (unsigned int)da.size();
    if(dim<1) DACEException(16,04);

    ac = new double[DA::getMaxMonomials()*(dim+2)];
    unsigned int nterms, nvars, nord;
    const DACEDA **mb = new const DACEDA*[dim];
    for(unsigned int i=0; i<dim; i++) mb[i] = &(da[i].m_index);
    daceEvalTree(mb,(int)dim,ac,nterms,nvars,nord);
    terms = (unsigned int) nterms;
    vars = (unsigned int) nvars;
    ord = (unsigned int) nord;
    delete[] mb;
    if(daceGetError()) DACEException();
}

/** Create a compiledDA object from a DA object.
    @param[in] da vector
    @throw DACE::DACEException
 */
compiledDA::compiledDA(const DA &da) {
    ac = new double[DA::getMaxMonomials()*3];
    dim = 1;
    unsigned int nterms, nvars, nord;
    const DACEDA *mb = &(da.m_index);
    daceEvalTree(&mb,(int)dim,ac,nterms,nvars,nord);
    terms = (unsigned int) nterms;
    vars = (unsigned int) nvars;
    ord = (unsigned int) nord;
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
    @param[in] cda compiledDA vector to be copied
    @return The compiledDA object with the same content of the given compiledDA
 */
compiledDA& compiledDA::operator=(const compiledDA &cda) {
    if(this != &cda){
        dim = cda.dim;
        terms = cda.terms;
        vars = cda.vars;
        ord = cda.ord;
        ac = new double[terms*(dim+2)];
        for(int i=terms*(dim+2)-1; i>=0; i--) ac[i] = cda.ac[i];}

    return *this;
}

/********************************************************************************
*     Evaluation overloads and template specialization
*********************************************************************************/
// double evaluation
template<> void compiledDA::eval(const std::vector<double> &args, std::vector<double> &res) const {
    const size_t narg = args.size();
    double *p = ac+2;
    double *xm = new double[ord+1];

    // prepare temporary powers
    xm[0] = 1.0;
    // constant part
    for(unsigned int i=0; i<dim; i++, p++)
        res[i] = (*p);
    // higher order terms
    for(unsigned int i=1; i<terms; i++){
        unsigned int jl = (unsigned int)(*p); p++;
        unsigned int jv = (unsigned int)(*p)-1; p++;
        if(jv < narg)
            xm[jl] = xm[jl-1]*args[jv];
        else
            xm[jl] = 0;
        for(unsigned int j=0; j<dim; j++, p++)
            res[j] += xm[jl]*(*p);
    }

    delete[] xm;
}

// DA evaluation
template<> void compiledDA::eval(const std::vector<DA> &args, std::vector<DA> &res) const {
    const size_t narg = args.size();
    unsigned int jlskip = ord+1;
    double *p = ac+2;
    DACEDA *xm = new DACEDA[ord+1];
    DACEDA tmp;

    // allocate temporary DA variables in dace
    for(unsigned int i=0; i<ord+1; i++) daceAllocateDA(xm[i],0);
    daceAllocateDA(tmp,0);
    // prepare temporary powers
    daceCreateConstant(xm[0],1.0);

    // constant part
    for(unsigned int i=0; i<dim; i++, p++)
        daceCreateConstant(res[i].m_index,*p);
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
        daceMultiply(xm[jl-1],args[jv].m_index,xm[jl]);
        for(unsigned int j=0; j<dim; j++, p++)
            if((*p)!=0.0)
            {
                daceMultiplyDouble(xm[jl],*p,tmp);
                daceAdd(res[j].m_index,tmp,res[j].m_index);
            }
    }
    // deallocate memory
    daceFreeDA(tmp);
    for(int i=ord; i>=0; i--) daceFreeDA(xm[i]);
    delete[] xm;

    if(daceGetError()) DACEException();
}

/********************************************************************************
*     Member access routines
*********************************************************************************/
/** Return the coefficient array.
    @return coefficient array
 */
const double* compiledDA::getAc() const {
    return this->ac;
}

/** Return the number of DAs (dimension).
    @return dimension
 */
unsigned int compiledDA::getDim() const {
    return this->dim;
}

/** Return the maximum order.
    @return maximum order
 */
unsigned int compiledDA::getOrd() const {
    return this->ord;
}

/** Return the maximum number of variables.
    @return maximum number of variables
 */
unsigned int compiledDA::getVars() const {
    return this->vars;
}

/** Return the number of terms.
    @return number of terms
 */
unsigned int compiledDA::getTerms() const {
    return this->terms;
}

}

/** @} */
