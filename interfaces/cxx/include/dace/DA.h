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
 * DA.h
 *
 *  Created on: Feb 24, 2014
 *      Author: Dinamica Srl
 */

/*  Main DA class and related utilities.

    This header file contains the DA class representing a single DA polynomial.
    It implements an object oriented interface to access all math function, as
    well as non-member functions for classical functional math notation.
    It also contains the storedDA class for serializing a DA object into an
    opaque binary representation that can be stored, transmitted, and converted
    back into a DA object without loss of precision.
*/

#ifndef DINAMICA_DA_H_
#define DINAMICA_DA_H_

// DACE C++ interface version (must match the version returned by DACEVER)
#define DACE_CPP_MAJOR (2)
#define DACE_CPP_MINOR (1)

// C++ stdlib classes used in this public interface
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <initializer_list>

#include "dace/dacecore.h"

/** Namespace containing the %DACE C++ interface.
    @ingroup DACECXX

    All components of the %DACE @ref DACECXX "C++ interface" live in this namespace.
    This includes the DA classes as well as the non-member functions acting on DA and
    AlgebraicVector objects.
 */
namespace DACE {

// forward declarations
class compiledDA;
class storedDA;
class DACEException;
class Monomial;
class Interval;
class DA;
template<typename T> class AlgebraicVector;

/** Basic DA class representing a single polynomial.
    @ingroup DACECXX

    Implements an object oriented interface to access all math function.

    Additionally, the DACE namespace provides non-member functions to access
    the math functions in classical functional math notation.

    @see DACE
    @see DACE::storedDA
 */
class DACE_API DA
{
    friend class compiledDA;
    friend class storedDA;

private:
    static bool initialized;                                                //!< Indicates if DA::init() was called
    static std::stack<unsigned int> TOstack;                                //!< Truncation order stack
    DACEDA m_index;                                                         //!< Index pointing to DA vector in DACE core library

public:
   typedef double value_type;                                               //!< Underlying type of DA coefficients (for compatibility with C++ std lib, boost, etc)

    /********************************************************************************
    *     DACE initialization and state related routines
    *********************************************************************************/
    /** @name Initialization & State
     * @{
     */
    static void init(const unsigned int ord, const unsigned int nvar);
    static bool isInitialized();
    static void version(int &maj, int &min, int &patch);
    static void checkVersion();
    static double setEps(const double eps);
    static double getEps();
    static double getEpsMac();
    static unsigned int getMaxOrder();
    static unsigned int getMaxVariables();
    static unsigned int getMaxMonomials();
    static unsigned int getTO();
    static unsigned int setTO(const unsigned int ot = DA::getMaxOrder());
    static void pushTO(const unsigned int ot = DA::getMaxOrder());
    static void popTO();
    /** @} */

    /********************************************************************************
    *     Constructors & Destructors (Variable Creation)
    *********************************************************************************/
    /** @name Constructors & Destructors (Variable Creation)
     * @{
     */
    DA();
    DA(const DA &da);
    DA(DA &&da);

    explicit DA(const int i, const double c = 1.0);
    explicit DA(const unsigned int i, const double c = 1.0);
    DA(const double c);

    ~DA() throw();
    /** @} */

    /********************************************************************************
    *     Coefficient information, access, and extraction
    *********************************************************************************/
    /** @name Coefficient Access
     * @{
     */
    unsigned int size() const;
    unsigned int order() const;
    int degree() const;
    int isnan() const;
    int isinf() const;
    double cons() const;
    AlgebraicVector<double> linear() const;
    AlgebraicVector<DA> gradient() const;
    double getCoefficient(const std::vector<unsigned int> &jj) const;
    void setCoefficient(const std::vector<unsigned int> &jj, const double coeff);
    Monomial getMonomial(const unsigned int npos) const;
    void getMonomial(const unsigned int npos, Monomial &m) const;
    std::vector<Monomial> getMonomials() const;
    /** @} */

    /********************************************************************************
    *     Assignments, Copying & Filtering
    *********************************************************************************/
    /** @name Assignment, Copying & Filtering
     * @{
     */
    DA& operator=(DA &&da);
    DA& operator=(const DA &da);
    DA& operator=(const double c);

    DA& operator+=(const DA &da);
    DA& operator+=(const double c);
    DA& operator-=(const DA &da);
    DA& operator-=(const double c);
    DA& operator*=(const DA &da);
    DA& operator*=(const double c);
    DA& operator/=(const DA &da);
    DA& operator/=(const double c);

    DA trim(const unsigned int min, const unsigned int max = DA::getMaxOrder()) const;
    /** @} */

    /********************************************************************************
    *     Basic arithmetic operations
    *********************************************************************************/
    /** @cond */
    friend DACE_API DA operator-(const DA &da);
    friend DACE_API DA operator+(const DA &da1, const DA &da2);
    friend DACE_API DA operator+(const DA &da, const double c);
    friend DACE_API DA operator+(const double c, const DA &da);
    friend DACE_API DA operator-(const DA &da1, const DA &da2);
    friend DACE_API DA operator-(const DA &da, const double c);
    friend DACE_API DA operator-(const double c, const DA &da);
    friend DACE_API DA operator*(const DA &da1, const DA &da2);
    friend DACE_API DA operator*(const DA &da, const double c);
    friend DACE_API DA operator*(const double c, const DA &da);
    friend DACE_API DA operator/(const DA &da1, const DA &da2);
    friend DACE_API DA operator/(const DA &da, const double c);
    friend DACE_API DA operator/(const double c, const DA &da);
    /** @endcond */

    /** @name Basic Arithmetic
     * See also: @ref DABasicArithmeticOperators "DA Basic Arithmetic Operators"
     * @{
     */
    DA multiplyMonomials(const DA &da) const;
    DA divide(const unsigned int var, const unsigned int p = 1) const;
    DA deriv(const unsigned int i) const;
    DA deriv(const std::vector<unsigned int> ind) const;
    DA integ(const unsigned int i) const;
    DA integ(const std::vector<unsigned int> ind) const;
    /** @} */

    /********************************************************************************
    *     Intrinsic functions
    *********************************************************************************/
    /** @name Intrinsics
     * @{
     */
    DA absolute() const;
    DA trunc() const;
    DA round() const;
    DA mod(const double p) const;
    DA pow(const int p) const;
    DA pow(const double p) const;
    DA root(const int p = 2) const;
    DA minv() const;
    DA sqr() const;
    DA sqrt() const;
    DA isrt() const;
    DA cbrt() const;
    DA icrt() const;
    DA hypot(const DA &da) const;
    DA exp() const;
    DA log() const;
    DA logb(const double b = 10.0) const;
    DA log10() const;
    DA log2() const;
    DA sin() const;
    DA cos() const;
    DA tan() const;
    DA asin() const;
    DA acos() const;
    DA atan() const;
    DA atan2(const DA &da) const;
    DA sinh() const;
    DA cosh() const;
    DA tanh() const;
    DA asinh() const;
    DA acosh() const;
    DA atanh() const;
    DA erf() const;
    DA erfc() const;
    DA BesselJFunction(const int n) const;
    DA BesselYFunction(const int n) const;
    DA BesselIFunction(const int n, const bool scaled = false) const;
    DA BesselKFunction(const int n, const bool scaled = false) const;
    DA GammaFunction() const;
    DA LogGammaFunction() const;
    DA PsiFunction(const unsigned int n) const;
    /** @} */

    /********************************************************************************
    *    Norm and estimation routines
    *********************************************************************************/
    /** @name Norm & Estimation
     * @{
     */
    double norm(const unsigned int type = 0) const;
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
    std::vector<double> estimNorm(std::vector<double> &err, const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
    Interval bound() const;
    double convRadius(const double eps, const unsigned int type = 1) const;
    /** @} */

    /********************************************************************************
    *     Polynomial evaluation routines & operators
    *********************************************************************************/
    /** @name Evaluation
     * @{
     */
    template<class T> typename T::value_type operator()(const T &args) const;
    template<class T> T operator()(const std::initializer_list<T> l) const;
    template<class T> T operator()(const T args[], const unsigned int length) const;
    template<class T> typename T::value_type eval(const T &args) const;
    template<class T> T eval(const std::initializer_list<T> l) const;
    template<class T> T eval(const T args[], const unsigned int length) const;
    template<class T> T evalScalar(const T &arg) const;
    compiledDA compile() const;
    DA plug(const unsigned int var, const double val = 0.0) const;
    double evalMonomials(const DA &values) const;
    DA replaceVariable(const unsigned int from = 0, const unsigned int to = 0, const double val = 1.0) const;
    DA scaleVariable(const unsigned int var = 0, const double val = 1.0) const;
    DA translateVariable(const unsigned int var = 0, const double a = 1.0, const double c = 0.0) const;
    /** @} */

    /********************************************************************************
    *     Input/output routines
    *********************************************************************************/
    /** @name Input/Output
     * @{
     */
    std::string toString() const;
    void write(std::ostream &os) const;
    /** @} */

    /********************************************************************************
    *     Static creation routines
    *********************************************************************************/
    /** @name Creation routines
     * @{
     */
    static DA random(const double cm);
    static DA id(const unsigned int var, const double c = 1.0);
    static DA identity(const unsigned int var, const double c = 1.0);
    static DA fromString(const std::string &str);
    static DA fromString(const std::vector<std::string> &str);
    static DA read(std::istream &is);
    /** @} */

    /********************************************************************************
    *     Miscellaneous routines
    *********************************************************************************/
    /** @name Miscellaneous
     * @{
     */
    static void memdump();
    /** @} */
};

/********************************************************************************
*     DA non-member functions
*********************************************************************************/
/** @name DA Coefficient Access Functions
 * @{
 */
DACE_API unsigned int size(const DA &da);
DACE_API unsigned int order(const DA &da);
DACE_API int degree(const DA &da);
DACE_API int isnan(const DA &da);
DACE_API int isinf(const DA &da);
DACE_API double cons(const DA &da);
DACE_API AlgebraicVector<double> linear(const DA &da);
DACE_API AlgebraicVector<DA> gradient(const DA &da);
/** @} */

/** @anchor DABasicArithmeticOperators
 * @name DA Basic Arithmetic Operators
 * @{
 */
DACE_API DA operator-(const DA &da);
DACE_API DA operator+(const DA &da1, const DA &da2);
DACE_API DA operator+(const DA &da, const double c);
DACE_API DA operator+(const double c, const DA &da);
DACE_API DA operator-(const DA &da1, const DA &da2);
DACE_API DA operator-(const DA &da, const double c);
DACE_API DA operator-(const double c, const DA &da);
DACE_API DA operator*(const DA &da1, const DA &da2);
DACE_API DA operator*(const DA &da, const double c);
DACE_API DA operator*(const double c, const DA &da);
DACE_API DA operator/(const DA &da1, const DA &da2);
DACE_API DA operator/(const DA &da, const double c);
DACE_API DA operator/(const double c, const DA &da);
/** @} */

/** @name DA Basic Arithmetic
 * @{
 */
DACE_API DA divide(const DA &da, const unsigned int var, const unsigned int p = 1);
DACE_API DA deriv(const DA &da, const unsigned int i);
DACE_API DA deriv(const DA &da, const std::vector<unsigned int> ind);
DACE_API DA integ(const DA &da, const unsigned int i);
DACE_API DA integ(const DA &da, const std::vector<unsigned int> ind);
/** @} */

/** @name DA Filtering Functions
 * @{
 */
DACE_API DA trim(const DA &da, const unsigned int min, const unsigned int max = DA::getMaxOrder());
/** @} */

/** @name DA Intrinsic Functions
 *  @{
 */
DACE_API DA absolute(const DA &da);
DACE_API DA trunc(const DA &da);
DACE_API DA round(const DA &da);
DACE_API DA mod(const DA &da, const double p);
DACE_API DA pow(const DA &da, const int p);
DACE_API DA pow(const DA &da, const double p);
DACE_API DA root(const DA &da, const int p = 2);
DACE_API DA minv(const DA &da);
DACE_API DA sqr(const DA &da);
DACE_API DA sqrt(const DA &da);
DACE_API DA isrt(const DA &da);
DACE_API DA cbrt(const DA &da);
DACE_API DA icrt(const DA &da);
DACE_API DA hypot(const DA &da1, const DA &da2);
DACE_API DA exp(const DA &da);
DACE_API DA log(const DA &da);
DACE_API DA logb(const DA &da, const double b = 10.0);
DACE_API DA log10(const DA &da);
DACE_API DA log2(const DA &da);
DACE_API DA sin(const DA &da);
DACE_API DA cos(const DA &da);
DACE_API DA tan(const DA &da);
DACE_API DA asin(const DA &da);
DACE_API DA acos(const DA &da);
DACE_API DA atan(const DA &da);
DACE_API DA atan2(const DA &da1, const DA &da2);
DACE_API DA sinh(const DA &da);
DACE_API DA cosh(const DA &da);
DACE_API DA tanh(const DA &da);
DACE_API DA asinh(const DA &da);
DACE_API DA acosh(const DA &da);
DACE_API DA atanh(const DA &da);
DACE_API DA erf(const DA &da);
DACE_API DA erfc(const DA &da);
DACE_API DA jn(const int n, const DA &da);
DACE_API DA yn(const int n, const DA &da);
DACE_API DA BesselJFunction(const int n, const DA &da);
DACE_API DA BesselYFunction(const int n, const DA &da);
DACE_API DA BesselIFunction(const int n, const DA &da, const bool scaled = false);
DACE_API DA BesselKFunction(const int n, const DA &da, const bool scaled = false);
DACE_API DA tgamma(const DA &da);
DACE_API DA lgamma(const DA &da);
DACE_API DA GammaFunction(const DA &da);
DACE_API DA LogGammaFunction(const DA &da);
DACE_API DA PsiFunction(const unsigned int n, const DA &da);
/** @} */

/** @name DA Norm & Estimation Functions
 * @{
 */
DACE_API double norm(const DA &da, unsigned int type = 0);
DACE_API std::vector<double> orderNorm(const DA &da, unsigned int var = 0, unsigned int type = 0);
DACE_API std::vector<double> estimNorm(const DA &da, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DA::getMaxOrder());
DACE_API std::vector<double> estimNorm(const DA &da, std::vector<double> &err, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DA::getMaxOrder());
DACE_API Interval bound(const DA &da);
DACE_API double convRadius(const DA &da, const double eps, const unsigned int type = 1);
/** @} */

/** @name DA Evaluation Functions and Operators
 * @{
 */
template<class T> typename T::value_type eval(const DA &da, const T &args);
template<class T> T eval(const DA &da, const std::initializer_list<T> l);
template<class T> T eval(const DA &da, const T args[], const unsigned int length);
template<class T> T evalScalar(const DA &da, const T &arg);
DACE_API compiledDA compile(const DA &da);
DACE_API DA plug(const DA &da, const unsigned int var, const double val = 0.0);
DACE_API DA replaceVariable(const DA &da, const unsigned int from = 0, const unsigned int to = 0, const double val = 1.0);
DACE_API DA scaleVariable(const DA &da, const unsigned int var = 0, const double val = 1.0);
DACE_API DA translateVariable(const DA &da, const unsigned int var = 0, const double a = 1.0, const double c = 0.0);
/** @} */

/** @name DA Input/Output Functions
 * @{
 */
DACE_API std::string toString(const DA &da);
DACE_API void write(const DA &da, std::ostream &os);
DACE_API std::ostream& operator<<(std::ostream &out, const DA &da);
DACE_API std::istream& operator>>(std::istream &in, DA &da);
/** @} */

/** Implementation of DA comparisons using only the constant part.
    To use these comparisons import them into the DACE namespace using
    @code
        namespace DACE { using namespace DACE::comp_cons; };
    @endcode
    @warning These comparisons only use the constant part. This may not be appropriate
    for all purposes. Read the Manual on @ref DACOMP for further information.
 */
namespace comp_cons {
    DACE_API bool operator<(const DA &lhs, const DA &rhs);
    DACE_API bool operator<(const double lhs, const DA &rhs);
    DACE_API bool operator<(const DA &lhs, const double rhs);
    DACE_API bool operator>(const DA &lhs, const DA &rhs);
    DACE_API bool operator>(const double lhs, const DA &rhs);
    DACE_API bool operator>(const DA &lhs, const double rhs);
    DACE_API bool operator<=(const DA &lhs, const DA &rhs);
    DACE_API bool operator<=(const double lhs, const DA &rhs);
    DACE_API bool operator<=(const DA &lhs, const double rhs);
    DACE_API bool operator>=(const DA &lhs, const DA &rhs);
    DACE_API bool operator>=(const double lhs, const DA &rhs);
    DACE_API bool operator>=(const DA &lhs, const double rhs);
    DACE_API bool operator==(const DA &lhs, const DA &rhs);
    DACE_API bool operator==(const double lhs, const DA &rhs);
    DACE_API bool operator==(const DA &lhs, const double rhs);
    DACE_API bool operator!=(const DA &lhs, const DA &rhs);
    DACE_API bool operator!=(const double lhs, const DA &rhs);
    DACE_API bool operator!=(const DA &lhs, const double rhs);
}

/** Implementation of abs(DA) using only the absolute constant part.
    To use this function import it into the DACE namespace using
    @code
        namespace DACE { using namespace DACE::abs_cons; };
    @endcode
 */
namespace abs_cons {
    DACE_API double abs(const DA &da);
}

/** Implementation of abs(DA) using the largest absolute coefficient.
    To use this function import it into the DACE namespace using
    @code
        namespace DACE { using namespace DACE::abs_max; };
    @endcode
 */
namespace abs_max {
    DACE_API double abs(const DA &da);
}

/** Implementation of abs(DA) using the sum of all absolute coefficients.
    To use this function import it into the DACE namespace using
    @code
        namespace DACE { using namespace DACE::abs_sum; };
    @endcode
 */
namespace abs_sum {
    DACE_API double abs(const DA &da);
}

/** Represents a DA vector in a binary, setup independent format.
    @ingroup DACECXX

    Use this class for serializing DA objects into an opaque binary representation
    that can be stored, transmitted, and converted back into DA object without
    loss of precision and across platforms.

    Example:
    @code
        DA x0 = sin(DA::id(1));
        DA y0 = cos(DA::id(1));
        std::fstream f("DA.bin", f.binary | f.trunc | f.in | f.out);
        f << storedDA(x0) << storedDA(y0);      // write as binary DAs
        f.seekp(0);
        DA x;
        f >> x;           // read directly into DA, handles both storedDA and text format
        storedDA ys(f);   // read storedDA from stream (uncommon)
        DA y = ys;        // convert storedDA back to DA
        std::cout << x0-x << y0-y;
    @endcode

    @see DACE::DA
 */
class DACE_API storedDA : std::vector<char>
{
private:
    static const unsigned int headerSize;           //!< size of the binary DA header (from DACE core)

public:
    storedDA(const DA &da);
    storedDA(const std::vector<char> &data);
    storedDA(std::istream &is);

    bool isValid() const;

    operator DA() const;
    operator std::string() const;

    friend DACE_API std::ostream& operator<<(std::ostream &out, const storedDA &sda);
};

/** @name Input/Output Functions
 * @{
 */
DACE_API std::ostream& operator<<(std::ostream &out, const storedDA &sda);
/** @} */

}

#endif /* DINAMICA_DA_H_ */
