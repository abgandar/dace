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
 * DA.cpp
 *
 *  Created on: Feb 24, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <algorithm>
#include <exception>
#include <sstream>
#include <cmath>

// DACE classes
#include "dace/config.h"
#include "dace/compiledDA.h"
#include "dace/DACEException.h"
#include "dace/Monomial.h"
#include "dace/Interval.h"
#include "dace/DA.h"
#include "dace/AlgebraicVector.h"
#include "dace/AlgebraicVector_t.h"

namespace DACE {

/********************************************************************************
*     DACE static class variables
*********************************************************************************/
bool DA::initialized = false;           // not initialized yet
std::stack<unsigned int> DA::TOstack;   // truncation order stack, initially empty

/********************************************************************************
*     DACE Setup
*********************************************************************************/
/** Initialize the DACE and set the maximum order and the
    maximum number of variables.
    @warning MUST BE CALLED BEFORE ANY OTHER DA ROUTINE CAN BE USED.
    @note This routine performs a mandatory version check to compare the version
    of the C++ interface used to compile the program to the version of the
    DACE library that is linked dynamically at runtime.
    @param[in] ord The order of the Taylor polynomials.
    @param[in] nvar The number of independent variables.
    @see DA::checkVersion
 */
void DA::init(const unsigned int ord, const unsigned int nvar) {
    try {
        checkVersion();
    } catch(DACEException const& ex) {
        std::cerr << ex << std::endl;
        std::terminate();
    }
    daceInitialize(ord,nvar);
    if(daceGetError()) DACEException();
    initialized = true;
}

/** Returns the inizialisation status of the DACE.
    @return True if the DACE has previously been initialized by a call
    to DA::init, or false otherwise.
    @see DA::init
 */
bool DA::isInitialized() {
    return initialized;
}

/** Return the major and minor version number of the DACE Core along with the patch level of the library.
    @param[out] maj The major DACE version number.
    @param[out] min The minor DACE version number.
    @param[out] patch The patch level of DACE version.
    @see DA::checkVersion
 */
void DA::version(int &maj, int &min, int &patch) {
    daceGetVersion(maj, min, patch);
    if(daceGetError()) DACEException();
}

/** Check the DACE core library version.
    Check the DACE core library version linked to this C++ interface against the interface version
    and throw an exception if the versions don't match.
    @note This routine is called automatically by DA::init to ensure compatibility with the current
    runtime environment.
    @throw DACE::DACEException
 */
void DA::checkVersion() {
    int maj, min, patch;
    version(maj, min, patch);
    if((maj!=DACE_CPP_MAJOR)||(min!=DACE_CPP_MINOR)||(maj!=DACE_MAJOR_VERSION)||(min!=DACE_MINOR_VERSION)) DACEException(20,99);
}

/** Return the maximum order currently set for DA operations.
    @return The maximum order, or zero if undefined.
    @throw DACE::DACEException
 */
unsigned int DA::getMaxOrder() {
    const unsigned int ord = daceGetMaxOrder();
    if(daceGetError()) DACEException();

    return ord;
}

/** Set the cutoff value @e eps to a new value and return the previous value.
    @param[in] eps The new cutoff value.
    @return The previous cutoff value, or zero if undefined.
    @throw DACE::DACEException
 */
double DA::setEps(const double eps) {
    const double old_eps = daceSetEpsilon(eps);
    if(daceGetError()) DACEException();

    return old_eps;
}

/** Return the cutoff value currently set for DA operations.
    @return The cutoff value, or zero if undefined.
    @throw DACE::DACEException
 */
double DA::getEps() {
    const double eps = daceGetEpsilon();
    if(daceGetError()) DACEException();

    return eps;
}

/** Return the machine epsilon (pessimistic estimate).
    @return The machine epsilon, or zero if undefined.
    @throw DACE::DACEException
 */
double DA::getEpsMac() {
    const double epsmac = daceGetMachineEpsilon();
    if(daceGetError()) DACEException();

    return epsmac;
}

/** Return the maximum number of variables set for the computations.
    @return The maximum number of variables, or zero if undefined.
    @throw DACE::DACEException
 */
unsigned int DA::getMaxVariables() {
    const unsigned int nvar = daceGetMaxVariables();
    if(daceGetError()) DACEException();

    return nvar;
}

/** Return the maximum number of monomials possible for the current
    order and number of variables.
    @return The maximum number of monomials, or zero if undefined.
    @throw DACE::DACEException
 */
unsigned int DA::getMaxMonomials() {
    const unsigned int nmmax = daceGetMaxMonomials();
    if(daceGetError()) DACEException();

    return nmmax;
}

/** Set the truncation order to a new value and return the previous value.
    All terms larger than the truncation order are discarded in subsequent operations.
    @param[in] ot The new truncation order.
    @return The previous truncation order, or zero if undefined.
    @throw DACE::DACEException
    @see DA::getTO
    @see DA::pushTO
    @see DA::popTO
 */
unsigned int DA::setTO(const unsigned int ot) {
    const unsigned int old_no = daceSetTruncationOrder(ot);
    if(daceGetError()) DACEException();

    return old_no;
}

/** Return the truncation order currently set for the computations.
    All terms larger than the truncation order are discarded in subsequent operations.
    @return The current truncation order, or zero if undefined.
    @throw DACE::DACEException
    @see DA::setTO
    @see DA::pushTO
    @see DA::popTO
 */
unsigned int DA::getTO() {
    const unsigned int no = daceGetTruncationOrder();
    if(daceGetError()) DACEException();

    return no;
}

/** Set a new truncation order, saving the previous one on the truncation order stack.
    All terms larger than the truncation order are discarded in subsequent operations.
    @param[in] ot The new truncation order.
    @throw DACE::DACEException
    @see DA::getTO
    @see DA::setTO
    @see DA::popTO
 */
void DA::pushTO(const unsigned int ot) {
    TOstack.push(daceSetTruncationOrder(ot));
    if(daceGetError()) DACEException();
}

/** Restore the previous truncation order from the truncation order stack.
    All terms larger than the truncation order are discarded in subsequent operations.
    @throw DACE::DACEException
    @see DA::getTO
    @see DA::setTO
    @see DA::pushTO
 */
void DA::popTO() {
    if(!TOstack.empty()) {
        daceSetTruncationOrder(TOstack.top());
        TOstack.pop();
        if(daceGetError()) DACEException();
    }
}

/********************************************************************************
*     Constructors & Destructors
*********************************************************************************/
/** Create an empty DA object representing the constant zero function.
    @throw DACE::DACEException
 */
DA::DA() {
    daceAllocateDA(m_index, 0);
    if(daceGetError()) DACEException();
}

/** Create a copy of a DA object.
    @param[in] da The DA object to be copied.
    @throw DACE::DACEException
 */
DA::DA(const DA &da) {
    daceAllocateDA(m_index, 0);
    daceCopy(da.m_index, m_index);
    if(daceGetError()) DACEException();
}

/** Move all resources from a DA object to us without copying.
    @param[in] da The DA object to be moved.
    @throw DACE::DACEException
 */
DA::DA(DA &&da) {
    m_index = da.m_index;
    daceInvalidateDA(da.m_index);   // prevent other DA from freeing our memory
    if(daceGetError()) DACEException();
}

/** Create a DA object with the constant part equal to @e c.
    @warning This routine *MUST* be called with a floating point type as the first argument, e.g. `DA(1.0)`.
    Expressions involving integer data types such as `DA(1)` will be interpreted as the first independent
    DA variable instead of the constant DA object of the given value.
    @param[in] c The double value to set as constant part of the DA object.
    @throw DACE::DACEException
 */
DA::DA(const double c) {
    daceAllocateDA(m_index, 0);
    daceCreateConstant(m_index, c);
    if(daceGetError()) DACEException();
}

/** Create a DA object as @e c times the independent variable number @e i.
    @deprecated The single argument version (with default 1.0 for c) is replaced by DA::id().
    It is recommended to always use DA::id() for consistency.
    @note When used in its one argument form (with the default argument 1.0 for c), this routine *MUST*
    be called with an integer type as the first argument, e.g. `DA(1)`.
    Expressions involving floating point data types such as `DA(1.0)` will be interpreted as the constant
    DA of the given value instead of the first independent variable.
    @param[in] i The independent variable number (i=0 means the constant part).
    @param[in] c The coefficient for the given independent variable.
    By default, this value is set to 1.0.
    @throw DACE::DACEException
    @see DA::id
 */
DA::DA(const unsigned int i, const double c) {
    daceAllocateDA(m_index, 0);
    daceCreateVariable(m_index, i, c);
    if(daceGetError()) DACEException();
}

/** Create a DA object as @e c times the independent variable number @e i.
    @deprecated The single argument version (with default 1.0 for c) is replaced by DA::id().
    It is recommended to always use DA::id() for consistency.
    @note When used in its one argument form (with the default argument 1.0 for c), this routine *MUST*
    be called with an integer type as the first argument, e.g. DA(1).
    Expressions involving floating point data types such as DA(1.0) will be interpreted as the constant
    DA of the given value instead of the first independent variable.
    @param[in] i The independent variable number (i=0 means the constant part).
    @param[in] c The coefficient for the given independent variable.
    By default, this value is set to 1.0.
    @throw DACE::DACEException
    @see DA::id
 */
DA::DA(const int i, const double c) {
    daceAllocateDA(m_index, 0);
    daceCreateVariable(m_index, (unsigned int)i, c);
    if(daceGetError()) DACEException();
}

/** Destroy a DA object and free the associated object in the DACE core.
 */
DA::~DA() throw() {
    daceFreeDA(m_index);

    // Never throw an exception from a destructor. Instead, we clear the error and ignore it. There is not much the user could do in this case anyway.
    if(daceGetError()) daceClearError();
}

/********************************************************************************
*     Coefficient information, access, and extraction routines
*********************************************************************************/

/** Return the number of non-zero coefficients of a DA object.
    @return The number of non-zero coefficients stored in the DA object.
    @throw DACE::DACEException
 */
unsigned int DA::size() const {
    const unsigned int res = daceGetLength(m_index);
    if(daceGetError()) DACEException();

    return res;
}

/** Return the order of a DA object.
    @return Lowest order of the non-zero monomials or @p UINT_MAX if it is the zero DA.
    @throw DACE::DACEException
 */
unsigned int DA::order() const {
    const unsigned int res = daceGetOrder(m_index);
    if(daceGetError()) DACEException();

    return res;
}

/** Return the degree of a DA object.
    @return Highest order of the non-zero monomials or @p INT_MIN if it is the zero DA.
    @throw DACE::DACEException
 */
int DA::degree() const {
    const int res = daceGetDegree(m_index);
    if(daceGetError()) DACEException();

    return res;
}

/** Check if a DA object has any NAN coefficients.
    @return True is any coefficients of the DA object are NAN.
    @throw DACE::DACEException
 */
int DA::isnan() const {
    const int temp = daceIsNan(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Check if a DA object has any INF coefficients.
    @return True is any coefficients of the DA object are INF.
    @throw DACE::DACEException
 */
int DA::isinf() const {
    const int temp = daceIsInf(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Return the constant part of a DA object.
    @return A double corresponding to the constant part of the DA object.
    @throw DACE::DACEException
 */
double DA::cons() const {
    const double temp = daceGetConstant(m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Return the linear part of a DA object.
    @return An AlgebraicVector<dobule> containing the linear coefficients of
    each independent DA variable in the DA object.
    @throw DACE::DACEException
 */
AlgebraicVector<double> DA::linear() const {
     AlgebraicVector<double> v(daceGetMaxVariables());

    daceGetLinear(m_index, v.data());
    if(daceGetError()) DACEException();

    return v;
}

/** Compute the gradient of the DA object.
    @return An AlgebraicVector<DA> containing the derivatives
    of the DA object with respect to all independent DA variables.
    @throw DACE::DACEException
 */
AlgebraicVector<DA> DA::gradient() const {
    const unsigned int nvar = daceGetMaxVariables();
    AlgebraicVector<DA> temp(nvar);
    for(unsigned int i = 0; i < nvar; i++) {
        temp[i] = deriv(i+1);
    }

    return temp;
}

/** Return a specific coefficient of a DA object.
    @param[in] jj A vector of the exponents of the coefficient to retrieve.
    @return The coefficient of DA object corresponding to the given vector of exponents.
    @throw DACE::DACEException
 */
double DA::getCoefficient(const std::vector<unsigned int> &jj) const {
    double coeff;
    const unsigned int nvar = daceGetMaxVariables();
    if(jj.size() >= nvar) {
        coeff = daceGetCoefficient(m_index, jj.data());
    } else {
        std::vector<unsigned int> temp(jj);
        temp.resize(nvar, 0);
        coeff = daceGetCoefficient(m_index, temp.data());
    }
    if(daceGetError()) DACEException();

    return coeff;
}

/** Set a specific coefficient into a DA object.
    @param[in] jj A vector of the exponents of the coefficient to be set.
    @param[in] coeff The value to be set as coefficient.
    @throw DACE::DACEException
 */
void DA::setCoefficient(const std::vector<unsigned int> &jj, const double coeff) {
    // check arguments
    const unsigned int nvar = daceGetMaxVariables();

    if(jj.size() >= nvar) {
        daceSetCoefficient(m_index, jj.data(), coeff);
    } else {
        std::vector<unsigned int> temp(jj);
        temp.resize(nvar, 0);
        daceSetCoefficient(m_index, temp.data(), coeff);
    }

    if(daceGetError()) DACEException();
}

/** Return the Monomial corresponding to the non-zero coefficient at the given
    position in the DA object (monomials use one based indexing!).
    @param[in] npos The position within the DA object. The ordering of the Monomials
    within a DA object is implementation dependent and does not correspond
    to the order in which Monomials are listed in the ASCII output routines.
    @return A Monomial object containing both the coefficient and
    the vector of exponents corresponding to the given position.
    If the requested monomial is not present in the DA object, a Monomial
    with coefficient set to 0.0 is returned.
    @throw DACE::DACEException
    @see Monomial
    @see DA::getMonomial
    @see DA::getMonomials
 */
Monomial DA::getMonomial(const unsigned int npos) const {
    Monomial m;
    getMonomial(npos, m);
    return m;
}

/** Return the Monomial corresponding to the non-zero coefficient at the given
    position in the DA object (monomials use one based indexing!).
    @param[in] npos The position within the DA object. The ordering of the Monomials.
    within a DA object is implementation dependent and does not correspond
    to the order in which Monomials are listed in the ASCII output routines.
    @param[out] m The monomial object in which to store the corresponding monomial.
    @throw DACE::DACEException
    @see Monomial
    @see DA::getMonomial
    @see DA::getMonomials
 */
void DA::getMonomial(const unsigned int npos, Monomial &m) const {
    daceGetCoefficientAt(m_index, (int)npos, m.m_jj.data(), m.m_coeff);
    if(daceGetError()) DACEException();
}

/** Return a vector of all Monomials in the DA object (differently from
    getMonomial() where only a single Monomial, corresponding to a specified
    position in the DA object, is returned).
    @return A vector of Monomial objects containing both the coefficient and
    the exponents corresponding to each monomial in the DA object. The monomials
    are returned in the same order as in the DACE ASCII output (that is, they
    are sorted by order).
    @throw DACE::DACEException
    @see Monomial
    @see DA::getMonomial
 */
std::vector<Monomial> DA::getMonomials() const {
    const unsigned int nord = daceGetMaxOrder();
    const unsigned int s = size();
    std::vector<Monomial> res(s), out(s);

    for(unsigned int i = 0; i < s; i++)
    daceGetCoefficientAt(m_index, i+1, res[i].m_jj.data(), res[i].m_coeff);
    if(daceGetError()) DACEException();

    // compute the order of each monomial
    std::vector<unsigned int> sum(s);
    for(unsigned int i = 0; i < s; i++) {
        sum[i] = res[i].order();
    }

    // sort monomials by order
    unsigned int k = 0;
    for(unsigned int ord = 0; ord <= nord; ord++) {
        for(unsigned int i = 0; i < s; i++) {
            if( sum[i] == ord ) {
                out[k] = res[i];
                k++;
            }
        }
    }

    return out;
}

/********************************************************************************
*     Assignments
*********************************************************************************/
/** Move all resources from a DA object to us without copying.
    @param[in] da The DA object to be moved.
    @throw DACE::DACEException
 */
DA& DA::operator=(DA &&da) {
    // do a switch-a-roo: we get theirs, they get ours, and then they can even free it for us!
    std::swap(m_index, da.m_index);

    return *this;
}

/** Copy the content of a different DA object into the DA object.
    @param[in] da The DA object to be copied.
    @return The DA object with the same content of the given DA object.
    @throw DACE::DACEException
 */
DA& DA::operator=(const DA &da) {
    if(this != &da) {
        daceCopy(da.m_index, m_index);
        if(daceGetError()) DACEException();
    }

    return *this;
}

/** Copy a constant polynomial of value @e c into the DA.
    @param[in] c The constant value to be copied.
    @return The DA object representing the constant function with value @e c.
    @throw DACE::DACEException
 */
DA& DA::operator=(const double c) {
    daceCreateConstant(m_index, c);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the sum between the DA object and the given one.
    The result is directly computed into the current DA object.
    @param[in] da The DA object to be added.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator+=(const DA &da) {
    daceAdd(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the sum between the current DA object and a given constant.
    The result is directly computed into the current DA object.
    @param[in] c The constant value to be added.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator+=(const double c) {
    daceAddDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the difference between the current DA object and the given one.
    The result is directly computed into the current DA object.
    @param[in] da The DA object to be subtracted.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator-=(const DA &da) {
    daceSubtract(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the difference between the current DA object and a given constant.
    The result is directly computed into the current DA object.
    @param[in] c constant value to be subtracted.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator-=(const double c) {
    daceSubtractDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the product between the current DA object and the given one.
    The result is directly computed into the current DA object.
    @param[in] da The DA object to be multiplied.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator*=(const DA &da) {
    daceMultiply(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the product between the current DA object and a given constant.
    The result is directly computed into the current DA object.
    @param[in] c The constant value to be multiplied.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator*=(const double c) {
    daceMultiplyDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the division between the current DA object and the given one.
    The result is directly computed into the current DA object.
    @param[in] da The DA object by which the current DA is divided.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator/=(const DA &da) {
    daceDivide(m_index, da.m_index, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/** Compute the division between the current DA object and a given constant.
    The result is directly computed into the current DA object.
    @param[in] c The constant value by which the current DA is divided.
    @return The current DA object with modified contents.
    @throw DACE::DACEException
 */
DA& DA::operator/=(const double c) {
    daceDivideDouble(m_index, c, m_index);
    if(daceGetError()) DACEException();

    return *this;
}

/********************************************************************************
*     Algebraic operations
*********************************************************************************/
/** Compute the additive inverse of the given DA object.
    @param[in] da The DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator-(const DA &da) {
    DA temp;
    daceMultiplyDouble(da.m_index, -1.0, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the addition between two DA objects.
    @param[in] da1 The first DA object.
    @param[in] da2 The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator+(const DA &da1, const DA &da2) {
    DA temp;
    daceAdd(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the addition between a DA object and a given constant.
    @param[in] da The DA object.
    @param[in] c The constant.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator+(const DA &da, const double c) {
    DA temp;
    daceAddDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the addition between a given constant and a DA object.
    @param[in] c The constant.
    @param[in] da the DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator+(const double c, const DA &da) {
    DA temp;
    daceAddDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the subtraction between two DA objects.
    @param[in] da1 The first DA object.
    @param[in] da2 The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator-(const DA &da1, const DA &da2) {
    DA temp;
    daceSubtract(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the subtraction between a DA object and a given constant.
    @param[in] da The DA object.
    @param[in] c The constant.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator-(const DA &da, const double c) {
    DA temp;
    daceSubtractDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the subtraction between a given constant and a DA object.
    @param[in] c The constant.
    @param[in] da The DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator-(const double c, const DA &da) {
    DA temp;
    daceDoubleSubtract(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the multiplication between two DA objects.
    @param[in] da1 The first DA object.
    @param[in] da2 The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator*(const DA &da1, const DA &da2) {
    DA temp;
    daceMultiply(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the multiplication between a DA object and a given constant.
    @param[in] da The DA object.
    @param[in] c The constant.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator*(const DA &da, const double c) {
    DA temp;
    daceMultiplyDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the multiplication between a given constant and a DA object.
    @param[in] c The constant.
    @param[in] da The DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator*(const double c, const DA &da) {
    DA temp;
    daceMultiplyDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the division between two DA objects.
    @param[in] da1 The first DA object.
    @param[in] da2 The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator/(const DA &da1, const DA &da2) {
    DA temp;
    daceDivide(da1.m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the division between a DA object and a given constant.
    @param[in] da The DA object.
    @param[in] c The constant.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator/(const DA &da, const double c) {
    DA temp;
    daceDivideDouble(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the division between a given constant and a DA object.
    @param[in] c The constant.
    @param[in] da The DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA operator/(const double c, const DA &da) {
    DA temp;
    daceDoubleDivide(da.m_index, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*     Math routines
*********************************************************************************/
/** Multiply the DA vector with another DA vector monomial by monomial.
    This is the equivalent of coefficient-wise multiplication (like in DA addition).
    @param[in] da The DA vector to multiply with coefficient-wise.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::evalMonomial
 */
DA DA::multiplyMonomials(const DA &da) const {
    DA temp;
    daceMultiplyMonomials(m_index, da.m_index, temp.m_index);
    if (daceGetError()) DACEException();

    return temp;
}

/** Divide by independent variable @e var raised to power @e p.
    @param[in] var The independente variable number to divide by.
    @param[in] p The power of the independent variable.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::divide(const unsigned int var, const unsigned int p) const {
    DA temp;
    daceDivideByVariable(m_index, var, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the derivative of a DA object with respect to variable @e i.
    @param[in] i The indepednent variable number with respect to which the derivative is calculated.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::deriv(const unsigned int i) const {
    DA temp;
    daceDifferentiate(i, m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the derivative of a DA object with respect to variables @e ind.
    If @e ind has fewer entries than there are independent variables, the missing
    entries are assumed to be zero. If @e ind has more entries than there are independent
    variables, extra values are ignored.

    @param[in] ind A vector containing the number of derivatives to take for each
    independent variable.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::deriv(const std::vector<unsigned int> ind) const {
    DA temp(*this);
    const unsigned int size = std::min((unsigned int)ind.size(),daceGetMaxVariables());

    for(unsigned int i=0; i<size; i++)
        for(unsigned int j=0; j<ind[i]; j++)
            daceDifferentiate((i+1), temp.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the integral of a DA object with respect to variable @e i.
    @param[in] i The independent variable number with respect to which the integral is calculated.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::integ(const unsigned int i) const {
    DA temp;
    daceIntegrate(i, m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the integral of a DA object with respect to variables @e i.
    If @e ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If @e ind has more
    than there are independent variables, extra values are ignored.

    @param[in] ind A vector containing the number of integrals to take for each
    independent variable.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::integ(const std::vector<unsigned int> ind) const {
    DA temp(*this);
    const unsigned int size = std::min((unsigned int)ind.size(),daceGetMaxVariables());

    for(unsigned int i=0; i<size; i++)
        for(unsigned int j=0; j<ind[i]; j++)
            daceIntegrate((i+1), temp.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Returns a DA object with all monomials of order less than @e min and greater
    than @e max removed.
    @param[in] min The minimum order to keep in the DA object.
    @param[in] max The maximum order to keep in the DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::trim(const unsigned int min, const unsigned int max) const {
    DA temp;
    daceTrim(m_index, min, max, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Absolute value of a DA object.
    Returns either the DA or the negative of the DA based on the constant part.
    @return A new DA object.
    @throw DACE::DACEException

    @see DA::abs
    @see DA::norm
 */
DA DA::absolute() const {
    DA temp;
    daceAbsolute(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Truncate the constant part of a DA object to an integer.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::trunc() const {
    DA temp;
    daceTruncate(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Round the constant part of a DA object to an integer.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::round() const {
    DA temp;
    daceRound(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Take floating-point remainder of constant part divided by double p.
    @param[in] p The double divisor.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::mod(const double p) const {
    DA temp;
    daceModuloDouble(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Take remainder of DA divided by another DA.
    This computes `x - trunc(cons(x)/cons(y))*y`.
    @param[in] da The DA divisor.
    @return A new DA object.
    @throw DACE::DACEException

    @see daceModulo()
 */
DA DA::mod(const DA &da) const {
    DA temp;
    daceModulo(m_index, da.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Elevate a DA object to a given integer power.
    @param[in] p The power to which the DA object is raised.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::pow(const int p) const {
    DA temp;
    dacePower(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Elevate a DA object to a given real power. The constant part must be positive.
    @param[in] p The power to which the DA object is raised.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::pow(const double p) const {
    DA temp;
    dacePowerDouble(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the p-th root of a DA object.
    @param[in] p The root to be computed.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::root(const int p) const {
    DA temp;
    daceRoot(m_index, p, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the multiplicative inverse of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::minv() const {
    DA temp;
    daceMultiplicativeInverse(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the square of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::sqr() const {
    DA temp;
    daceSquare(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the square root of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::sqrt() const {
    DA temp;
    daceSquareRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the inverse square root of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::isrt() const {
    DA temp;
    daceInverseSquareRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the cubic root of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::cbrt() const {
    DA temp;
    daceCubicRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the inverse cubic root of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::icrt() const {
    DA temp;
    daceInverseCubicRoot(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hypotenuse (`sqrt(a*a+b*b)`) of a DA object and the given DA argument.
    @param[in] da The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::hypot(const DA &da) const {
    DA temp;
    daceHypotenuse(m_index, da.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the exponential of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::exp() const {
    DA temp;
    daceExponential(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the natural logarithm of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::log() const {
    DA temp;
    daceLogarithm(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the logarithm of a DA object with respect to a given base.
    @param[in] b The base with respect to which the logarithm is computed (default: base 10).
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::logb(const double b) const {
    DA temp;
    daceLogarithmBase(m_index, b, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the 10 based logarithm of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::log10() const {
    DA temp;
    daceLogarithm10(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the 2 based logarithm of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::log2() const {
    DA temp;
    daceLogarithm2(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the sine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::sin() const {
    DA temp;
    daceSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the cosine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::cos() const {
    DA temp;
    daceCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the tangent of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::tan() const {
    DA temp;
    daceTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the arcsine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::asin() const {
    DA temp;
    daceArcSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the arccosine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::acos() const {
    DA temp;
    daceArcCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the arctangent of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::atan() const {
    DA temp;
    daceArcTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the four-quadrant arctangent of Y/X. Y is the current DA object,
    whereas X is the given DA.
    @param[in] da The second DA object containing X.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::atan2(const DA &da) const {
    DA temp;
    daceArcTangent2(m_index, da.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic sine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::sinh() const {
    DA temp;
    daceHyperbolicSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic cosine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::cosh() const {
    DA temp;
    daceHyperbolicCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic tangent of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::tanh() const {
    DA temp;
    daceHyperbolicTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic arcsine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::asinh() const {
    DA temp;
    daceHyperbolicArcSine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic arccosine of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::acosh() const {
    DA temp;
    daceHyperbolicArcCosine(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the hyperbolic arctangent of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::atanh() const {
    DA temp;
    daceHyperbolicArcTangent(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the error function of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::erf() const {
    DA temp;
    daceErrorFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the complementary error function of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::erfc() const {
    DA temp;
    daceComplementaryErrorFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the @e n-th Bessel function of first type @f$J_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @note This function fails if the result is too large to be represented in double precision.
    @param[in] n The order of the Bessel function.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::BesselJFunction(const int n) const {
    DA temp;
    daceBesselJFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the @e n-th Bessel function of second type @f$Y_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @note This function fails if the result is too large to be represented in double precision.
    @param[in] n The order of the Bessel function.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::BesselYFunction(const int n) const {
    DA temp;
    daceBesselYFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the @e n-th modified Bessel function of first type @f$I_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @note This function fails if the result is too large to be represented in double precision.
    @param[in] n The order of the Bessel function.
    @param[in] scaled If true, the modified Bessel function is scaled
    by a factor `exp(-x)`, i.e. `exp(-x)I_n(x)` is returned.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::BesselIFunction(const int n, const bool scaled) const {
    DA temp;
    daceBesselIFunction(m_index, n, scaled, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the @e n-th modified Bessel function of second type @f$K_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @note This function fails if the result is too large to be represented in double precision.
    @param[in] n The order of the Bessel function.
    @param[in] scaled If true, the modified Bessel function is scaled.
    by a factor `exp(x)`, i.e. `exp(x)K_n(x)` is returned.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::BesselKFunction(const int n, const bool scaled) const {
    DA temp;
    daceBesselKFunction(m_index, n, scaled, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the Gamma function of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::GammaFunction() const {
    DA temp;
    daceGammaFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the Logarithmic Gamma function (i.e. the natural logarithm of Gamma) of a DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::LogGammaFunction() const {
    DA temp;
    daceLogGammaFunction(m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Compute the @e n-th order Psi function, i.e. the (n+1)st derivative of the Logarithmic Gamma
    function, of a DA object.
    @param[in] n The order of the Psi function.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::PsiFunction(const unsigned int n) const {
    DA temp;
    dacePsiFunction(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the Legendre polynomial of degree @e n.
    @param[in] n The degree of the Legendre polynomial.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::LegendrePolynomial(const unsigned int n) const {
    DA temp;
    daceLegendrePolynomial(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the associated Legendre polynomial of degree @e n and order @e m.
    @param[in] n The degree of the associated Legendre polynomial.
    @param[in] m The order of the associated Legendre polynomial.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::AssociatedLegendrePolynomial(const unsigned int n, const unsigned int m) const {
    DA temp;
    daceAssociatedLegendrePolynomial(m_index, n, m, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the Hermite polynomial of degree @e n.
    @param[in] n The degree of the Hermite polynomial.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::HermitePolynomial(const unsigned int n) const {
    DA temp;
    daceHermitePolynomial(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the Laguerre polynomial of degree @e n.
    @param[in] n The degree of the Laguerre polynomial.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::LaguerrePolynomial(const unsigned int n) const {
    DA temp;
    daceLaguerrePolynomial(m_index, n, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the associated Laguerre polynomial of degree @e n and order @e m.
    @param[in] n The degree of the associated Laguerre polynomial.
    @param[in] m The order of the associated Laguerre polynomial.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::AssociatedLaguerrePolynomial(const unsigned int n, const unsigned int m) const {
    DA temp;
    daceAssociatedLaguerrePolynomial(m_index, n, m, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluate the spherical harmonic of degree @e n and order @e m.
    This is the spherical harmonic @f$Y_n^m(\theta, \varphi)@f$ for @f$\varphi = 0@f$
    sometimes also called the spherical associated Legendre functions.
    To obtain the full complex spherical harmonic, multiply by @f$e^{im\varphi}@f$.
    @param[in] n The degree of the spherical harmonic.
    @param[in] m The order of the spherical harmonic.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::SphericalHarmonic(const unsigned int n, const unsigned int m) const {
    DA temp;
    daceSphericalHarmonic(m_index, n, m, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Beta function (Euler integral of first kind).
    @param[in] da2 A DA object.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::BetaFunction(const DA &da2) const {
    DA temp;
    daceBetaFunction(m_index, da2.m_index, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*    Norm and estimation routines
*********************************************************************************/
/** Compute different types of norms for a DA object.
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @return The resulting norm.
    @throw DACE::DACEException
 */
double DA::norm(const unsigned int type) const {
    const double c = daceNorm(m_index, type);
    if(daceGetError()) DACEException();

    return c;
}

/** Extract different types of order sorted norms from a DA object.
    @param[in] var The grouping\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @return The resulting norms for each group.
    @throw DACE::DACEException
 */
std::vector<double> DA::orderNorm(const unsigned int var, const unsigned int type) const {
    std::vector<double> v(daceGetMaxOrder()+1);
    daceOrderedNorm(m_index, var, type, v.data());
    if(daceGetError()) DACEException();

    return v;
}

/** Estimate different types of order sorted norms for terms of a DA object
    up to a specified order.
    @note If estimation is not possible, zero is returned for all requested orders.
    @param[in] var The grouping:\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable @e var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @param[in] nc The maximum order to be estimated (default: max order).
    @return The resulting norms and estimated norms for each group.
    @throw DACE::DACEException
 */
std::vector<double> DA::estimNorm(const unsigned int var, const unsigned int type, const unsigned int nc) const {
    std::vector<double> v(nc+1);
    daceEstimate(m_index, var, type, v.data(), NULL, nc);
    if(daceGetError()) DACEException();

    return v;
}

/** Estimate different types of order sorted norms for terms of a DA object
    up to a specified order with error estimates.
    @note If estimation is not possible, zero is returned for all requested orders.
    @param[out] err Returns the amount by which the estimate underestimates the actual ordered norm
    of the terms in the polynomial up to the minimum of @e nc or the maximum computation order.
    @param[in] var The grouping\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable @e var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @param[in] nc The maximum order to be estimated (default: max order).
    @return The resulting norms and estimated norms for each group.
    @throw DACE::DACEException
 */
std::vector<double> DA::estimNorm(std::vector<double> &err, const unsigned int var, const unsigned int type, const unsigned int nc) const {
    std::vector<double> v(nc+1);
    err.resize(std::min(nc, daceGetMaxOrder())+1);
    daceEstimate(m_index, var, type, v.data(), err.data(), nc);
    if(daceGetError()) DACEException();

    return v;
}

/** Compute lower and upper bounds of a DA object over the domain @f$[-1,1]^n@f$.
    @return An Interval object containing both the lower and the upper bound of the DA object.
    @throw DACE::DACEException
    @see DACE::Interval
 */
Interval DA::bound() const {
    Interval i;
    daceGetBounds(m_index, i.m_lb, i.m_ub);
    if(daceGetError()) DACEException();

    return i;
}

/** Estimate the convergence radius of the DA object.
    Based on the coefficients of the DA object, evaluating the DA with any values for the
    independent variables of norm less than the returned convergence radius are estimated
    to have a truncation error of less than @e eps.
    @warning This is an estimate based on the assumption of exponential convergence of the
    coefficients. It is not a rigorous bound.
    @param[in] eps The requested tolerance.
    @param[in] type The type of norm (default: sum norm).
    @return The estimated convergence radius.
    @throw DACE::DACEException
 */
double DA::convRadius(const double eps, const unsigned int type) const {
    const unsigned int ord = daceGetTruncationOrder();

    std::vector<double> res = estimNorm(0, type, ord+1);
    return std::pow(eps/res[ord+1],1.0/(ord+1));
}

/********************************************************************************
*     DACE polynomial evaluation routines
*********************************************************************************/
/** Compile current DA object and create a compiledDA object.
    @return The compiled DA object.
    @throw DACE::DACEException
 */
compiledDA DA::compile() const {
    return compiledDA(*this);
}

/** Partial evaluation of a DA object. In the DA object, variable @e var is
    replaced by the value @e val. The resulting DA object is returned.
    @param[in] var The independent variable number to be replaced.
    @param[in] val The value by which to replace the variable.
    @return A new DA object.
    @throw DACE::DACEException
 */
DA DA::plug(const unsigned int var, const double val) const {
    DA temp;
    daceEvalVariable(m_index,var,val,temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Evaluates the DA vector using the coefficients in argument values as the values for each monomial.
    This is equivalent to a monomial-wise dot product of two DA vectors.
    @param[in] values A DA vector containing the values of each monomial.
    @return The result of the evaluation.
    @throw DACE::DACEException
    @see DA::multiplyMonomial
 */
double DA::evalMonomials(const DA &values) const {
    const double res = daceEvalMonomials(m_index, values.m_index);
    if (daceGetError()) DACEException();

    return res;
}

/** Partial evaluation of a DA object. In the DA object, variable @e from is
    replaced by the value @e val times variable @e to.
    @param[in] from The independent variable number to be replaced.
    @param[in] to The independent variable number to be inserted instead.
    @param[in] val The value by which to scale the inserted variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::replaceVariable
 */
DA DA::replaceVariable(const unsigned int from, const unsigned int to, const double val) const {
    DA temp;
    daceReplaceVariable(m_index, from, to, val, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Scaling of an independent variable. In the DA object, variable @e var is
    replaced by the value @e val times var.
    @param[in] var The independent variable number to be scaled.
    @param[in] val The value by which to scale the variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::scaleVariable
 */
DA DA::scaleVariable(const unsigned int var, const double val) const {
    DA temp;
    daceScaleVariable(m_index, var, val, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/** Affine translation of an independent variable. In the DA object, variable @e var is
    replaced by <tt>a*var+c</tt>.
    @param[in] var The independent variable number to be translated.
    @param[in] a The value by which to scale the variable.
    @param[in] c The value by which to shift the variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::translateVariable
 */
DA DA::translateVariable(const unsigned int var, const double a, const double c) const {
    DA temp;
    daceTranslateVariable(m_index, var, a, c, temp.m_index);
    if(daceGetError()) DACEException();

    return temp;
}

/********************************************************************************
*     DACE input/output routines
*********************************************************************************/
/** Convert DA object to string.
    @return A string.
    @throw DACE::DACEException
 */
std::string DA::toString() const {
    // initialize 2D char array
    unsigned int nstr = daceGetMaxMonomials() + 2;
    char *ss = new char[nstr*DACE_STRLEN];

    daceWrite( m_index, ss, nstr );

    // copy from char array to string
    std::string s;
    for(unsigned int i=0; i<(unsigned int)nstr; i++) {
        ss[(i+1)*DACE_STRLEN-1] = '\0'; // should already be done by the Fortran code but just to be sure terminate string
        s.append(&ss[i*DACE_STRLEN]);
        s.append(1,'\n');
    }

    // delete 2D char array
    delete[] ss;
    if(daceGetError()) DACEException();

    return s;
}

/** Output operator.
    @param[in] out A C++ output stream.
    @param[in] da A DA object to be output to the stream.
    @return The C++ output stream.
    @throw DACE::DACEException
    @see DA::toString
 */
std::ostream& operator<<(std::ostream &out, const DA &da) {
    out << toString(da);

    return out;
}

/** Input operator. Reads both string and binary DA representations from a stream.
    @note When using binary IO operations, make sure the stream is opened in ios_base::binary mode!
    Some C++ libraries are known to mangle the input otherwise which will break the ability to read binary DA objects.
    Setting the binary flag for all IO (also text based) does not affect the reading and is recommended.
    @param[in] in A C++ input stream.
    @param[in] da The DA object to read to from the stream.
    @return The C++ input stream.
    @throw DACE::DACEException
    @see DA::fromString
 */
std::istream& operator>>(std::istream &in, DA &da) {
    storedDA sda(in);

    if(sda.isValid()) {
        da = sda;                                       // automatically cast sDA to DA
        return in;
    } else {
        const std::string endstr = "------------------------------------------------";  // from daceio.c
        std::string line = sda;                         // automatically cast sDA to string
        std::vector<std::string> strs;
        strs.reserve(5000);                             // estimate that 5000 lines will be an OK guess for many use cases

        if(line.length() > 0)
        {
            // parse the content of line into string array
            const std::string::iterator end = line.end();
            std::string::iterator p = line.begin();
            while(p != end)
            {
                const std::string::iterator p0 = p;
                while(p != end && *p != '\n') p++;
                strs.emplace_back(p0, p);
                if(p != end) p++;
            }

            // complete the last line, in case it was not a full line that was read
            if(*(end-1) != '\n')
            {
                getline(in, line);
                strs[strs.size()-1] += line;
            }
        }

        if(!strs.empty()) {
            if(!strs.back().empty()) {
                // check that last line is not the terminator line and in case remove it
                if(strs.back().compare(4, 31, endstr, 0, 31) == 0) {
                    strs.pop_back();
                } else {
                    // read the istream until end string is found and put each line in the string vector (end condition taken from daceio.c)
                    for(getline(in, line); in.good() && (line.compare(4, 31, endstr, 0, 31) != 0); getline(in, line))
                    strs.push_back(line);
                }
            }
            // convert string vector to DA
            da = DA::fromString(strs);
        }
    }

    return in;
}

/** Write a binary representation of the DA as a blob to output stream.
    @param[in] os A C++ output stream.
    @throw DACE::DACEException
 */
void DA::write(std::ostream &os) const {
    os << storedDA(*this);
}


/********************************************************************************
*     DACE static factory routines
*********************************************************************************/
/** Create a DA object and fill with random entries.
    For cm < 0, the DA object is filled with random numbers.
    For cm > 0, the DA object is filled with weighted decaying numbers.
    @param[in] cm The filling factor.
    @throw DACE::DACEException
 */
DA DA::random(const double cm) {
    DA temp;
    daceCreateRandom(temp.m_index, cm);
    if(daceGetError()) DACEException();

    return temp;
}

/** Create a DA object representing the identity function in independent DA
    variable number @e var.
    @param[in] var The independent DA variable number.
    @param[in] c The coefficient for the independent variable.
    @throw DACE::DACEException
 */
DA DA::id(const unsigned int var, const double c) {
    return DA(var, c);
}

/** Create a DA object representing the identity function in independent DA
    variable number @e var.
    Legacy alias for DA::id().
    @deprecated Replaced by DA::id().
    @param[in] var The independent DA variable number.
    @param[in] c The coefficient for the independent variable.
    @throw DACE::DACEException
    @see DA::id
 */
DA DA::identity(const unsigned int var, const double c) {
    return id(var, c);
}

/** Convert a string to DA object.
    @param[in] str A string.
    @return A DA object.
    @throw DACE::DACEException
 */
DA DA::fromString(const std::string &str) {
    std::istringstream ssin(str);
    DA temp;
    ssin >> temp;

    return temp;
}

/** Convert a vector of strings to DA object.
    @param[in] str A vector of strings, each representing one line of the input.
    @return A DA object.
    @throw DACE::DACEException
 */
DA DA::fromString(const std::vector<std::string> &str) {
    // create 2D char array
    const unsigned int nstr = (unsigned int)str.size();
    char *ss = new char[nstr*DACE_STRLEN];

    // fill the array with blanks
    for(unsigned int i=0; i<nstr*DACE_STRLEN; i++)
        ss[i] = ' ';

    // fill char array with rows of the string vector
    for(unsigned int i=0; i<nstr; i++) {
        str[i].copy(&ss[i*DACE_STRLEN], DACE_STRLEN);
    }

    DA da;
    daceRead(da.m_index, ss, nstr);

    // delete 2D char array
    delete[] ss;
    if(daceGetError()) DACEException();

    return da;
}

/** Read a binary representation of a DA from input stream.
    @note When using binary IO operations, make sure the stream is opened in ios_base::binary mode!
    Some C++ libraries are known to mangle the input otherwise which will break the ability to read binary DA objects.
    Setting the binary flag for all IO (also text based) does not affect the reading and is recommended.
    @param[in] is A C++ input stream.
    @throw DACE::DACEException
 */
DA DA::read(std::istream &is) {
    storedDA sda(is);
    return sda;         // automatically cast to DA
}

/********************************************************************************
*     DACE various routines
*********************************************************************************/
/** Print debugging information about the internal memory status of the DACE core.
 */
void DA::memdump() {
    daceMemoryDump();
}

/********************************************************************************
*     DACE non-member functions
*********************************************************************************/
/** Return the number of non-zero coefficients of a DA object.
    @param[in] da The DA object.
    @return The number of non-zero coefficients in the DA object.
    @throw DACE::DACEException
    @see DA::size
 */
unsigned int size(const DA &da) {
    return da.size();
}

/** Return the order of a DA object.
    @param[in] da The DA object.
    @return Lowest order of the non-zero monomials or @p UINT_MAX if it is the zero DA.
    @throw DACE::DACEException
    @see DA::order
 */
unsigned int order(const DA &da) {
    return da.order();
}

/** Return the degree of a DA object.
    @param[in] da The DA object.
    @return Highest order of the non-zero monomials or @p INT_MIN if it is the zero DA.
    @throw DACE::DACEException
    @see DA::degree
 */
int degree(const DA &da) {
    return da.degree();
}

/** Check if a DA object has any NAN coefficients.
    @param[in] da A DA object.
    @return True if any coefficients of the DA object are NAN.
    @throw DACE::DACEException
 */
int isnan(const DA &da) {
    return da.isnan();
}

/** Check if a DA object has any INF coefficients.
    @param[in] da A DA object.
    @return True if any coefficients of the DA object are INF.
    @throw DACE::DACEException
 */
int isinf(const DA &da) {
    return da.isinf();
}

/** Return the constant part of a DA object.
    @param[in] da A DA object.
    @return A double corresponding to the constant part of the DA object.
    @throw DACE::DACEException
 */
double cons(const DA &da) {
    return da.cons();
}

/** Return the linear part of a DA object.
    @param[in] da A DA object.
    @return An AlgebraicVector<dobule> containing the linear coefficients of
    each independent DA variable in the given DA object.
    @throw DACE::DACEException
 */
AlgebraicVector<double> linear(const DA &da) {
    return da.linear();
}

/** Compute the gradient of a DA object.
    @param[in] da the given DA object.
    @return A AlgebraicVector<DA> containing the derivatives of the DA object with respect to all
    independent DA variables.
    @throw DACE::DACEException
 */
AlgebraicVector<DA> gradient(const DA &da) {
    return da.gradient();
}

/** Divide by independent variable @e var raised to power @e p.
    @param[in] da A DA object.
    @param[in] var The independent variable number to divide by.
    @param[in] p The power of variable @e var to divide by.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::divide
 */
DA divide(const DA &da, const unsigned int var, const unsigned int p) {
    return da.divide(var, p);
}

/** Compute the derivative of a DA object with respect to variable i.
    @param[in] da A DA object.
    @param[in] i The independent variable number with respect to which the derivative is calculated.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::deriv
 */
DA deriv(const DA &da, const unsigned int i) {
    return da.deriv(i);
}

/** Compute the derivative of a DA object with respect to variables @e ind.
    @param[in] da A DA object.
    @param[in] ind A vector containing the number of derivatives to take for each
    independent variable. If @e ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If @e ind has more
    entries than there are independent variables, extra values are ignored.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::deriv
 */
DA deriv(const DA &da, const std::vector<unsigned int> ind) {
    return da.deriv(ind);
}

/** Compute the integral of a DA object with respect to variable @e i.
    @param[in] da A DA object.
    @param[in] i Independent variable number with respect to which the integral is calculated.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::integ
 */
DA integ(const DA &da, const unsigned int i) {
    return da.integ(i);
}

/** Compute the integral of a DA object with respect to variable @e i.
    @param[in] da A DA object.
    @param[in] ind A vector containing the number of derivatives to take for each
    independent variable. If @e ind has fewer entries than there are independent
    variables, the missing entries are assumed to be zero. If @e ind has more
    entries than there are independent variables, extra values are ignored.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::integ
 */
DA integ(const DA &da, const std::vector<unsigned int> ind) {
    return da.integ(ind);
}

/** Returns a DA object with all monomials of order less than @e min and greater
    than @e max removed (trimmed).
    @param[in] da A DA object.
    @param[in] min The minimum order to keep in the DA object.
    @param[in] max The maximum order to keep in the DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::trim
 */
DA trim(const DA &da, const unsigned int min, const unsigned int max) {
    return da.trim(min,max);
}

/** Absolute value of a DA object.
    Returns either the DA or the negative of the DA based on the constant part.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::absolute
    @see DA::abs
    @see DA::norm
 */
DA absolute(const DA &da) {
    return da.absolute();
}

/** Truncate the constant part of a DA object to an integer.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::trunc
 */
DA trunc(const DA &da) {
    return da.trunc();
}

/** Round the constant part of a DA object to an integer.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::round
 */
DA round(const DA &da) {
    return da.round();
}

/** Compute the floating-point remainder of c/p (c modulo p),
    where @e c is the constant part of the given DA object.
    @param[in] da A DA object.
    @param[in] p The costant with respect to which the modulo function is computed.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::mod
 */
DA mod(const DA &da, double p) {
    return da.mod(p);
}

/** Compute remainder of one DA divided by another DA.
    This computes `x - trunc(cons(x)/cons(y))*y`.
    @param[in] da1 A DA object.
    @param[in] da2 The DA object with respect to which the modulo function is computed.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::mod
    @see daceModulo()
 */
DA mod(const DA &da1, const DA &da2) {
    return da1.mod(da2);
}

/** Raise a DA object to a given integer power.
    @param[in] da A DA object.
    @param[in] p The power to which the DA object is raised.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::pow
 */
DA pow(const DA &da, int p) {
    return da.pow(p);
}

/** Raise a DA object to a given integer power. The constant part must be positive.
    @param[in] da A DA object.
    @param[in] p The power to which the DA object is raised.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::pow
 */
DA pow(const DA &da, double p) {
    return da.pow(p);
}

/** Compute the @e p-th root of a DA object.
    @param[in] da A DA object.
    @param[in] p The root to be computed.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::root
 */
DA root(const DA &da, int p) {
    return da.root(p);
}

/** Compute the multiplicative inverse of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::minv
 */
DA minv(const DA &da) {
    return da.minv();
}

/** Compute the square of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::sqr
 */
DA sqr(const DA &da) {
    return da.sqr();
}

/** Compute the square root of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::sqrt
 */
DA sqrt(const DA &da) {
    return da.sqrt();
}

/** Compute the inverse square root of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::isrt
 */
DA isrt(const DA &da) {
    return da.isrt();
}

/** Compute the cubic root of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::cbrt
 */
DA cbrt(const DA &da) {
    return da.cbrt();
}

/** Compute the inverse cubic root of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::icrt
 */
DA icrt(const DA &da) {
    return da.icrt();
}

/** Compute the hypotenuse (`sqrt(X*X + Y*Y)`) of two DA objects.
    @param[in] X The first DA object.
    @param[in] Y The second DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::hypot
 */
DA hypot(const DA &X, const DA &Y) {
    return X.hypot(Y);
}

/** Compute the exponential of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::exp
 */
DA exp(const DA &da) {
    return da.exp();
}

/** Compute the natural logarithm of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::log
 */
DA log(const DA &da) {
    return da.log();
}

/** Compute the logarithm of a DA object with respect to a given base.
    @param[in] da A DA object.
    @param[in] b The base with respect to which the logarithm is computed (default: base 10).
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::logb
 */
DA logb(const DA &da, const double b) {
    return da.logb(b);
}

/** Compute the 10 based logarithm of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::log10
 */
DA log10(const DA &da) {
    return da.log10();
}

/** Compute the 2 based logarithm of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::log2
 */
DA log2(const DA &da) {
    return da.log2();
}

/** Compute the sine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::sin
 */
DA sin(const DA &da) {
    return da.sin();
}

/** Compute the cosine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::cos
 */
DA cos(const DA &da) {
    return da.cos();
}

/** Compute the tangent of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::tan
 */
DA tan(const DA &da) {
    return da.tan();
}

/** Compute the arcsine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::asin
 */
DA asin(const DA &da) {
    return da.asin();
}

/** Compute the arccosine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::acos
 */
DA acos(const DA &da) {
    return da.acos();
}

/** Compute the arctangent of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::atan
 */
DA atan(const DA &da) {
    return da.atan();
}

/** Compute the four-quadrant arctangent of Y/X.
    @param[in] Y A DA object.
    @param[in] X A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::atan2
 */
DA atan2(const DA &Y, const DA &X) {
    return Y.atan2(X);
}

/** Compute the hyperbolic sine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::sinh
 */
DA sinh(const DA &da) {
    return da.sinh();
}

/** Compute the hyperbolic cosine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::cosh
 */
DA cosh(const DA &da) {
    return da.cosh();
}

/** Compute the hyperbolic tangent of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::tanh
 */
DA tanh(const DA &da) {
    return da.tanh();
}

/** Compute the hyperbolic arcsine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::asinh
 */
DA asinh(const DA &da) {
    return da.asinh();
}

/** Compute the hyperbolic arccosine of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::acosh
 */
DA acosh(const DA &da) {
    return da.acosh();
}

/** Compute the hyperbolic arctangent of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::atanh
 */
DA atanh(const DA &da) {
    return da.atanh();
}

/** Compute the error function of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::erf
 */
DA erf(const DA &da) {
    return da.erf();
}

/** Compute the complementary error function of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::erfc
 */
DA erfc(const DA &da) {
    return da.erfc();
}

/** Compute the @e n-th Bessel function of first type @f$J_n@f$ of a DA object.
    Alias of BesselJFunction for C compatible naming.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselJFunction
 */
DA jn(const int n, const DA &da) {
    return da.BesselJFunction(n);
}

/** Compute the @e n-th Bessel function of second type @f$Y_n@f$ of a DA object.
    Alias of BesselYFunction for C compatible naming.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselYFunction
 */
DA yn(const int n, const DA &da) {
    return da.BesselYFunction(n);
}

/** Compute the @e n-th Bessel function of first type @f$J_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselJFunction
 */
DA BesselJFunction(const int n, const DA &da) {
    return da.BesselJFunction(n);
}

/** Compute the @e n-th Bessel function of second type @f$Y_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselYFunction
 */
DA BesselYFunction(const int n, const DA &da) {
    return da.BesselYFunction(n);
}

/** Compute the @e n-th modified Bessel function of first type @f$I_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @param[in] scaled If true, the modified Bessel function is scaled.
    by a factor `exp(-x)`, i.e. `exp(-x)I_n(x)` is returned.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselIFunction
 */
DA BesselIFunction(const int n, const DA &da, const bool scaled) {
    return da.BesselIFunction(n, scaled);
}

/** Compute the @e n-th modified Bessel function of second type @f$K_n@f$ of a DA object.
    @note The DA must have non-negative constant part while the order is allowed to be negative.
    @param[in] n The order of the Bessel function.
    @param[in] da A DA object.
    @param[in] scaled If true, the modified Bessel function is scaled.
    by a factor `exp(-x)`, i.e. `exp(-x)K_n(x)` is returned.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BesselKFunction
 */
DA BesselKFunction(const int n, const DA &da, const bool scaled) {
    return da.BesselKFunction(n, scaled);
}

/** Compute the Gamma function of a DA object.
    Alias of GammaFunction() for C99 compatible naming.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::GammaFunction
 */
DA tgamma(const DA &da) {
    return da.GammaFunction();
}

/** Compute the Log Gamma function of a DA object.
    Alias of LogGammaFunction() for C99 compatible naming.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::LogGammaFunction
 */
DA lgamma(const DA &da) {
    return da.LogGammaFunction();
}

/** Compute the Gamma function of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::GammaFunction
 */
DA GammaFunction(const DA &da) {
    return da.GammaFunction();
}

/** Compute the Log Gamma function of a DA object.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::LogGammaFunction
 */
DA LogGammaFunction(const DA &da) {
    return da.LogGammaFunction();
}

/** Compute the @e n-th Psi function of a DA object.
    @param[in] n The order of the Psi function to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::PsiFunction
 */
DA PsiFunction(const unsigned int n, const DA &da) {
    return da.PsiFunction(n);
}

/** Evaluate the Legendre polynomial of degree @e n.
    @param[in] n The degree of the Legendre polynomial to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::LegendrePolynomial
 */
DA LegendrePolynomial(const unsigned int n, const DA &da) {
    return da.LegendrePolynomial(n);
}

/** Evaluate the associated Legendre polynomial of degree @e n and order @e m.
    @param[in] n The degree of the associated Legendre polynomial to compute.
    @param[in] m The order of the associated Legendre polynomial to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::AssociatedLegendrePolynomial
 */
DA AssociatedLegendrePolynomial(const unsigned int n, const unsigned int m, const DA &da) {
    return da.AssociatedLegendrePolynomial(n, m);
}

/** Evaluate the Hermite polynomial of degree @e n.
    @param[in] n The degree of the Hermite polynomial to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::HermitePolynomial
 */
DA HermitePolynomial(const unsigned int n, const DA &da) {
    return da.HermitePolynomial(n);
}

/** Evaluate the Laguerre polynomial of degree @e n.
    @param[in] n The degree of the Laguerre polynomial to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::LaguerrePolynomial
 */
DA LaguerrePolynomial(const unsigned int n, const DA &da) {
    return da.LaguerrePolynomial(n);
}

/** Evaluate the Laguerre polynomial of degree @e n and order @e m.
    @param[in] n The degree of the associated Laguerre polynomial to compute.
    @param[in] m The degree of the associated Laguerre polynomial to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::AssociatedLaguerrePolynomial
 */
DA AssociatedLaguerrePolynomial(const unsigned int n, const unsigned int m, const DA &da) {
    return da.AssociatedLaguerrePolynomial(n, m);
}

/** Evaluate the spherical harmonic of degree @e n and order @e m.
    @param[in] n The degree of the spherical harmonic to compute.
    @param[in] m The degree of the spherical harmonic to compute.
    @param[in] da A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::SphericalHarmonic
 */
DA SphericalHarmonic(const unsigned int n, const unsigned int m, const DA &da) {
    return da.SphericalHarmonic(n, m);
}

/** Beta function (Euler integral of first kind).
    @param[in] da1 A DA object.
    @param[in] da2 A DA object.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::BetaFunction
 */
DA BetaFunction(const DA &da1, const DA &da2) {
    return da1.BetaFunction(da2);
}

/** Compute different types of norms for a DA object.
    @param[in] da A DA object.
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @return The resulting norm.
    @throw DACE::DACEException
    @see DA::norm
 */
double norm(const DA &da, unsigned int type) {
    return da.norm(type);
}

/** Compute different types of order sorted norms for terms of a DA object.
    @param[in] da A DA object.
    @param[in] var The grouping\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @return The resulting norms for each group.
    @throw DACE::DACEException
    @see DA::onorm
 */
std::vector<double> orderNorm(const DA &da, const unsigned int var, const unsigned int type) {
    return da.orderNorm(var, type);
}

/** Estimate different types of order sorted norms for terms of a DA object
    up to a specified order.
    @note If estimation is not possible, zero is returned for all requested orders.
    @param[in] da A DA object.
    @param[in] var The grouping:\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @param[in] nc The maximum order (defualt: max order).
    @return The resulting norms and estimated norms for each group.
    @throw DACE::DACEException
    @see DA::estim
 */
std::vector<double> estimNorm(const DA &da, const unsigned int var, const unsigned int type, const unsigned int nc) {
    return da.estimNorm(var, type, nc);
}

/** Estimate different types of order sorted norms for terms of a DA object
    up to a specified order with error estimates.
    @note If estimation is not possible, zero is returned for all requested orders.
    @param[in] da A DA object.
    @param[out] err Returns the amount by which the estimate underestimates
    the actual ordered norm of the terms in the polynomial up to the minimum
    of nc or the maximum computation order
    @param[in] var The grouping:\n
     0: Terms are sorted by their order (Default)\n
    >0: Terms are sorted by the exponent of variable var
    @param[in] type The type of norm to be computed. Possible norms are:\n
     0: Max norm (Default)\n
     1: Sum norm\n
    >1: Vector norm of given type
    @param[in] nc The maximum order (default: max order).
    @return A double corresponding to the result of the operation.
    @throw DACE::DACEException
    @see DA::estim
 */
std::vector<double> estimNorm(const DA &da, std::vector<double> &err, const unsigned int var, const unsigned int type, const unsigned int nc) {
    return da.estimNorm(err, var, type, nc);
}

/** Compute lower and upper bounds of a DA object over the domain @f$[-1,1]^n@f$.
    @param[in] da A DA object.
    @return An Interval object containing both the lower and the upper bound
    of the DA object.
    @throw DACE::DACEException
    @see DA::bound
 */
Interval bound(const DA &da) {
    return da.bound();
}

/** Estimate the convergence radius of the DA object.
    Based on the coefficients of the DA object, evaluating the DA with any values for the
    independent variables of norm less than the returned convergence radius are estimated
    to have a truncation error of less than @e eps.
    @warning This is an estimate based on the assumption of exponential convergence of the
    coefficients. It is not a rigorous bound.
    @param[in] da the given DA object.
    @param[in] eps The requested tolerance.
    @param[in] type The type of norm (default: sum norm).
    @return The estimated convergence radius.
    @throw DACE::DACEException
    @see DA::conv_radius
 */
double convRadius(const DA &da, const double eps, const unsigned int type) {
    return da.convRadius(eps, type);
}

/** Partial evaluation of a DA object.
    In the DA object, variable @e var is replaced by the value @e val.
    @param[in] da A DA object.
    @param[in] var The independent variable number to be replaced.
    @param[in] val the value by which to replace the variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::plug
 */
DA plug(const DA &da, const unsigned int var, const double val) {
    return da.plug(var, val);
}

/** Partial evaluation of a DA object.
    In the DA object, variable @e from is replaced by the value @e val times variable @e to.
    @param[in] da A DA object.
    @param[in] from The independent variable number to be replaced.
    @param[in] to The independent variable number to be inserted instead.
    @param[in] val The value by which to scale the inserted variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::replaceVariable
 */
DA replaceVariable(const DA &da, const unsigned int from, const unsigned int to, const double val) {
    return da.replaceVariable(from, to, val);
}

/** Scaling of an independent variable.
    In the DA object, variable @e var is replaced by the value @e val times @e var.
    @param[in] da A DA object.
    @param[in] var The independent variable number to be scaled.
    @param[in] val the value by which to scale the variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::scaleVariable
 */
DA scaleVariable(const DA &da, const unsigned int var, const double val) {
    return da.scaleVariable(var, val);
}

/** Affine translation of an independent variable.
    In the DA object, variable @e var is replaced by `a*var + c`.
    @param[in] da A DA object.
    @param[in] var The independent variable number to be translated.
    @param[in] a The value by which to scale the variable.
    @param[in] c The value by which to shift the variable.
    @return A new DA object.
    @throw DACE::DACEException
    @see DA::translateVariable
 */
DA translateVariable(const DA &da, const unsigned int var, const double a, const double c) {
    return da.translateVariable(var, a, c);
}

/** Compile a given DA object and create a compiledDA object.
    @return A compiled DA object.
    @throw DACE::DACEException
    @see DA::compile
 */
compiledDA compile(const DA &da) {
    return da.compile();
}

/** Convert DA object to a string.
    @param[in] da A DA object.
    @return A string representing the DA.
    @throw DACE::DACEException
    @see DA::toString
 */
std::string toString(const DA &da) {
    return da.toString();
}

/** Write binary representation of DA object to given output stream.
    @param[in] da A DA object.
    @param[in] os A C++ output stream.
    @throw DACE::DACEException
    @see DA::write
 */
void write(const DA &da, std::ostream &os) {
    return da.write(os);
}

namespace comp_cons {
    bool operator<(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)<cons(rhs);
    }

    bool operator<(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs<cons(rhs);
    }

    bool operator<(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)<rhs;
    }

    bool operator>(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)>cons(rhs);
    }

    bool operator>(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs>cons(rhs);
    }

    bool operator>(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)>rhs;
    }

    bool operator<=(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)<=cons(rhs);
    }

    bool operator<=(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs<=cons(rhs);
    }

    bool operator<=(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)<=rhs;
    }

    bool operator>=(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)>=cons(rhs);
    }

    bool operator>=(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs>=cons(rhs);
    }

    bool operator>=(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)>=rhs;
    }

    bool operator==(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)==cons(rhs);
    }

    bool operator==(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs==cons(rhs);
    }

    bool operator==(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)==rhs;
    }

    bool operator!=(const DA &lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A DA object.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)!=cons(rhs);
    }

    bool operator!=(const double lhs, const DA &rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A double.
    @param[in] rhs A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return lhs!=cons(rhs);
    }

    bool operator!=(const DA &lhs, const double rhs) {
/** Compare constant part of DAs.
    @warning This operator only compares constant parts! It does not check any higher order terms
    so DAs comparing equal may still represent different polynomials.
    @param[in] lhs A DA object.
    @param[in] rhs A double.
    @throw DACE::DACEException
    @see DA::cons
 */
        return cons(lhs)!=rhs;
    }
}

namespace abs_cons {
    double abs(const DA &da) {
/** Absolute value of constant part.
    @param[in] da A DA object.
    @throw DACE::DACEException
    @see DA::cons
 */
        return std::abs(da.cons());
    }
}

namespace abs_max {
    double abs(const DA &da) {
/** Largest coefficient in absolute value.
    @param[in] da A DA object.
    @throw DACE::DACEException
    @see DA::norm
 */
        return da.norm(0);
    }
}

namespace abs_sum {
    double abs(const DA &da) {
/** Sum of absolute values of all coefficients.
    @param[in] da A DA object.
    @throw DACE::DACEException
    @see DA::norm
 */
        return da.norm(1);
    }
}

// static class variable with size of binary header as reported by DACE core
const unsigned int storedDA::headerSize = daceBlobSize(NULL);

/** Create new storedDA from an existing DA.
    @param[in] da A DA object.
 */
storedDA::storedDA(const DA &da) {
    unsigned int len;

    daceExportBlob(da.m_index, NULL, len);
    resize(len);
    daceExportBlob(da.m_index, data(), len);
    if(daceGetError()) DACEException();
}

/** Create new storedDA by copying from a buffer.
    @param[in] data A vector of bytes.
 */
storedDA::storedDA(const std::vector<char> &data) : std::vector<char>(data) {
}

/** Create a new storedDA by reading from a stream.
    @param[in] is A C++ input stream.
 */
storedDA::storedDA(std::istream &is) : std::vector<char>(storedDA::headerSize) {
    // read blob header
    is.read(data(), headerSize);
    if(is.gcount() != headerSize)
    {
        resize((size_t)is.gcount());
        return;
    }

    // check validity of blob header and read the rest
    const unsigned int len = daceBlobSize(data());
    if(len == 0)
    {
        return;
    }
    else if(len > headerSize)
    {
        resize(len);
        is.read(data()+headerSize, len-headerSize);
        if(is.gcount() != (len-headerSize))
        {
            resize((size_t)(headerSize+is.gcount()));
            return;
        }
    }
}

/** Return if this storedDA data appears to be valid.
    @return True if storedDA appears valid.
 */
bool storedDA::isValid() const {
    const size_t s1 = size();
    // check first that we have the minimum amount of data for the DACE blob header
    if(s1 < headerSize)
        return false;

    const unsigned int s2 = daceBlobSize(data());
    // is the blob valid
    if(s2 == 0)
        return false;

    // do we have the amount of data claimed in the header
    return s1 >= s2;
}

/** Convert stored data to string.
    @return A string of the data.
 */
storedDA::operator std::string() const {
    return std::string(data(), size());
}

/** Convert storedDA back to DA.
    @return The DA.
 */
storedDA::operator DA() const {
    DA da;

    if(isValid())
    {
        daceImportBlob(data(), da.m_index);
        if(daceGetError()) DACEException();
    }
    else
    {
        DACEException(15, 06);
    }

    return da;
}

/** storedDA stream output operator. Writes binary data to output stream.
    @param[in] out A C++ output stream.
    @param[in] sda The storedDA.
    @return The C++ output stream.
 */
std::ostream& operator<<(std::ostream &out, const storedDA &sda) {
    out.write(sda.data(), sda.size());

    return out;
}

}
