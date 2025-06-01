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
 * AlgebraicVector.h
 *
 *  Created on: Sep. 10, 2014
 *      Author: Dinamica Srl
 */

/*  Main AlgebraicVector class.

    This header file contains the AlgebraicVector class to simplify math using vectors
    of any algebraic type, including DA and double. It derives from std::vector, so
    can be used in all library interfaces where a std::vector is accepted.
*/

#ifndef DINAMICA_ALGEBRAICVECTOR_H_
#define DINAMICA_ALGEBRAICVECTOR_H_

// C++ stdlib classes required for interface definition
#include <vector>
#include <initializer_list>
#include <type_traits>

// DACE classes required for interface definition (DA.h needed for DA::getMaxOrder(), DA::getMaxVariables() default arguments)
#include "dace/DA.h"

namespace DACE {

/** Vector of any algebraic type.
    @ingroup DACECXX

    Provides vector-vector, vector-scalar, and componentwise operations.

    This class is templated allowing it to store elements of any algebraic type T.
    However, some of the member functions (e.g. AlgebraicVector::eval) require the
    datatype T to be DACE::DA (or something that implements the same interface).
    Most basic operations that make sense for scalars and DAs will work for both.

    Example:
    @code
        #include <dace/dace.h>

        DA::init(5, 3);                         // order 5, 3 variables

        AlgebraicVector<double> dbl(3);         // vector of length 3
        AlgebraicVector<DA> da(3);

        dbl = {1.0, 2.0, 3.0};                  // braced initializer list notation
        da = {DA::id(1), DA::id(2), DA::id(3)};
        //da = AlgebraicVector<DA>::identity(); // same as above

        da = sin(da) + 3.1*dbl - 1.0;           // can mix scalars, vector types, result upcast to DA

        da = da.deriv(1) + vectorDA{0.0, DA::id(2), DA::id(3)};     // vectorDA is short for AlgebraicVector<DA>
        //dbl.deriv(1);                         // compiler error: no derivative of double

        std::cout << da[0] - cos(DA::id(1));        // prints a zero DA vector
    @endcode
 */
template<typename T> class AlgebraicVector : public std::vector<T>
{
public:
    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    /** @name Constructors & Destructors
     * @{
     */
    /** Default constructor to create empty AlgebraicVector.
    */
    AlgebraicVector() : std::vector<T>() {};

    /** Constructor with size to allocate a vector of the given size with elements initialized using their default constructor.
        @param[in] size The initial length of the AlgebraicVector.
    */
    explicit AlgebraicVector(const size_t size) : std::vector<T>(size) {};

    /** Constructor with size and value to allocate a vector of the given size with elements initialized as copies of d.
        @param[in] size The initial length of the AlgebraicVector.
        @param[in] d    The initial value for the elements.
    */
    AlgebraicVector(const size_t size, const T &d) : std::vector<T>(size, d) {};

    /** Copy constructor to create a copy of any existing vector.
        @param[in] v A vector to be copied into AlgebraicVector.
    */
    AlgebraicVector(const std::vector<T> &v) : std::vector<T>(v) {};

    /** Extraction constructor to copy only a given range of elements from vector v.
        @note The constructor does not perform any range checking for the extraction.
        @param[in] v A vector to be copied into AlgebraicVector.
        @param[in] first The index of the first element to be copied.
        @param[in] last The index of the last element to be copied.
        @see AlgebraicVector<T>::extract
    */
    AlgebraicVector(const std::vector<T> &v, const size_t first, const size_t last) : std::vector<T>(v.begin()+first, v.begin()+last+1) {};

    /** Constructor to create a vector from an initializer list.
        @param[in] l A braced initializer list to be copied into the AlgebraicVector.
    */
    AlgebraicVector(std::initializer_list<T> l) : std::vector<T>(l) {};
    /** @} */

    /********************************************************************************
    *     Element access
    *********************************************************************************/
    /** @name Element Access
     * @{
     */
    AlgebraicVector<T> extract(const size_t first, const size_t last) const;
    template<typename U> AlgebraicVector<typename std::common_type_t<T, U>> concat(const std::vector<U> &obj) const;
    template<typename U> AlgebraicVector<T>& operator<<(const std::vector<U> &obj);
    /** @} */

    /********************************************************************************
    *     Coefficient access
    *********************************************************************************/
    /** @name Coefficient Access
     * @{
     */
    AlgebraicVector<double> cons() const;
    std::vector<std::vector<double>> linear() const;
    /** @} */

    /********************************************************************************
    *     Assignments, Copying & Filtering
    *********************************************************************************/
    /** @name Assignment, Copying & Filtering
     * @{
     */
    template<typename U> AlgebraicVector<T>& operator+=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator+=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const U &obj);

    AlgebraicVector<T> trim(const unsigned int min, const unsigned int max = DA::getMaxOrder()) const;
    /** @} */

    /********************************************************************************
    *     Basic arithmetic operations
    *********************************************************************************/
    /** @name Basic Arithmetic
     * See also: @ref AlgebraicVectorBasicArithmeticOperators "Algebraic Vector Basic Arithmetic Operators"
     * @{
     */
    AlgebraicVector<T> operator-() const;
    AlgebraicVector<T> deriv(const unsigned int p) const;
    AlgebraicVector<T> integ(const unsigned int p) const;
    /** @} */

    /********************************************************************************
    *     Intrinsic functions
    *********************************************************************************/
    /** @name Intrinsics
     * @{
     */
    AlgebraicVector<T> absolute() const;
    AlgebraicVector<T> trunc() const;
    AlgebraicVector<T> round() const;
    template<typename U> AlgebraicVector<T> mod(const U &p) const;
    AlgebraicVector<T> pow(const int p) const;
    AlgebraicVector<T> pow(const double p) const;
    AlgebraicVector<T> root(const int p = 2) const;
    AlgebraicVector<T> minv() const;
    AlgebraicVector<T> sqr() const;
    AlgebraicVector<T> sqrt() const;
    AlgebraicVector<T> isrt() const;
    AlgebraicVector<T> cbrt() const;
    AlgebraicVector<T> icbrt() const;
    AlgebraicVector<T> hypot(const AlgebraicVector<T> &Y) const;
    AlgebraicVector<T> exp() const;
    AlgebraicVector<T> log() const;
    AlgebraicVector<T> logb(const double b = 10.0) const;
    AlgebraicVector<T> log10() const;
    AlgebraicVector<T> log2() const;
    AlgebraicVector<T> sin() const;
    AlgebraicVector<T> cos() const;
    AlgebraicVector<T> tan() const;
    AlgebraicVector<T> asin() const;
    AlgebraicVector<T> acos() const;
    AlgebraicVector<T> atan() const;
    AlgebraicVector<T> atan2(const AlgebraicVector<T> &X) const;
    AlgebraicVector<T> sinh() const;
    AlgebraicVector<T> cosh() const;
    AlgebraicVector<T> tanh() const;
    AlgebraicVector<T> asinh() const;
    AlgebraicVector<T> acosh() const;
    AlgebraicVector<T> atanh() const;
    AlgebraicVector<T> erf() const;
    AlgebraicVector<T> erfc() const;
    AlgebraicVector<T> BesselJFunction(const int n) const;
    AlgebraicVector<T> BesselYFunction(const int n) const;
    AlgebraicVector<T> BesselIFunction(const int n, const bool scaled = false) const;
    AlgebraicVector<T> BesselKFunction(const int n, const bool scaled = false) const;
    AlgebraicVector<T> GammaFunction() const;
    AlgebraicVector<T> LogGammaFunction() const;
    AlgebraicVector<T> PsiFunction(const unsigned int n) const;
    /** @} */

    /***********************************************************************************
    *    Vector routines
    ************************************************************************************/
    /** @name Vector routines
     * @{
     */
    template<typename U> typename std::common_type_t<T, U> dot(const AlgebraicVector<U> &obj) const;
    template<typename U> AlgebraicVector<typename std::common_type_t<T, U>> cross(const AlgebraicVector<U> &obj) const;
    T length() const;
    AlgebraicVector<T> normalize() const;
    AlgebraicVector<T> invert() const;
    // XXX: various Jacobians, gradients, curls, etc?
    /** @} */

    /********************************************************************************
    *     Norm and estimation routines
    *********************************************************************************/
    /** @name Norm & Estimation
     * @{
     */
    AlgebraicVector<double> norm(const unsigned int type = 0) const;
    /* XXX: define and add the norm estimation routines from DA including convergence radius estimation
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
    Interval bound() const;
    double convRadius(const double eps, const unsigned int type = 1) const;
    */
    /** @} */

    /********************************************************************************
    *     Polynomial evaluation operators and routines
    *********************************************************************************/
    /** @name Evaluation
     * @{
     */
    template<typename U> U operator()(const U &args) const;
    template<typename U> AlgebraicVector<U> operator()(const std::initializer_list<U> l) const;
    template<typename U> AlgebraicVector<U> operator()(const U args[], const unsigned int length) const;

    template<typename U> U eval(const U &args) const;
    template<typename U> AlgebraicVector<U> eval(const std::initializer_list<U> l) const;
    template<typename U> AlgebraicVector<U> eval(const U args[], const unsigned int length) const;

    template<typename U> AlgebraicVector<U> evalScalar(const U &arg) const;

    compiledDA compile() const;
    AlgebraicVector<T> plug(const unsigned int var, const double val = 0.0) const;
    /** @} */

    /***********************************************************************************
    *     Input/Output routines
    ************************************************************************************/
    /** @name Input/Output
     *  @{
     */
    std::string toString() const;
    /** @} */

    /********************************************************************************
    *     Static creation routines
    *********************************************************************************/
    /** @name Creation routines
     *  @{
     */
    static AlgebraicVector<DA> id(const size_t n = DA::getMaxVariables());
    static AlgebraicVector<DA> identity(const size_t n = DA::getMaxVariables());
    /** @} */

private:
    static void matrix_inverse(std::vector<std::vector<double>> &A);        // Private helper routine for double precision matrix inversion
};

/********************************************************************************
*     AlgebraicVector non-member functions
*********************************************************************************/
/** @name Vector Coefficient Access Functions
 *  @{
 */
template<typename T> AlgebraicVector<double> cons(const AlgebraicVector<T> &obj);
template<typename T> std::vector<std::vector<double>> linear(const AlgebraicVector<T> &obj);
/** @} */

/** @anchor AlgebraicVectorBasicArithmeticOperators
 * @name Vector Basic Arithmetic Operators
 * @{
 */
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator+( const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator+( const AlgebraicVector<T> &obj1, const U &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator+( const T &obj1, const AlgebraicVector<U> &obj2);

template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator-( const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator-( const AlgebraicVector<T> &obj1, const U &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator-( const T &obj1, const AlgebraicVector<U> &obj2);

template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator*( const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator*( const AlgebraicVector<T> &obj1, const U &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator*( const T &obj1, const AlgebraicVector<U> &obj2);

template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator/( const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator/( const AlgebraicVector<T> &obj1, const U &obj2);
template<typename T,typename U> AlgebraicVector<typename std::common_type_t<T, U>> operator/( const T &obj1, const AlgebraicVector<U> &obj2);
/** @} */

/** @name Vector Basic Arithmetic
 * @{
 */
template<typename T> AlgebraicVector<T> deriv(const AlgebraicVector<T> &obj, const unsigned int p);
template<typename T> AlgebraicVector<T> integ(const AlgebraicVector<T> &obj, const unsigned int p);
/** @} */

/** @name Vector Filtering Functions
 * @{
 */
template<typename T> AlgebraicVector<T> trim(const AlgebraicVector<T> &obj, unsigned int min, unsigned int max = DA::getMaxOrder());
/** @} */

 /** @name Vector Intrinsic Functions
 *  @{
 */
template<typename T> AlgebraicVector<T> absolute(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> trunc(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> round(const AlgebraicVector<T> &obj);
template<typename T, typename U> AlgebraicVector<T> mod(const AlgebraicVector<T> &obj, const U &p);
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const double p);
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const int p);
template<typename T> AlgebraicVector<T> root(const AlgebraicVector<T> &obj, const int p = 2);
template<typename T> AlgebraicVector<T> minv(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqr(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> isrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cbrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> icbrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> hypot(const AlgebraicVector<T> &X, const AlgebraicVector<T> &Y);
template<typename T> AlgebraicVector<T> exp(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> log(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> logb(const AlgebraicVector<T> &obj, const double b = 10.0);
template<typename T> AlgebraicVector<T> log10(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> log2(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sin(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cos(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> tan(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> asin(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> acos(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atan(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atan2(const AlgebraicVector<T> &Y, const AlgebraicVector<T> &X);
template<typename T> AlgebraicVector<T> sinh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cosh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> tanh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> asinh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> acosh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> atanh(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> erf(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> erfc(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> jn(const int n, const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> yn(const int n, const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cyl_bessel_i(const int n, const AlgebraicVector<T> &da);
template<typename T> AlgebraicVector<T> cyl_bessel_j(const int n, const AlgebraicVector<T> &da);
template<typename T> AlgebraicVector<T> cyl_bessel_k(const int n, const AlgebraicVector<T> &da);
template<typename T> AlgebraicVector<T> cyl_neumann(const int n, const AlgebraicVector<T> &da);
template<typename T> AlgebraicVector<T> BesselJFunction(const int n, const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> BesselYFunction(const int n, const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> BesselIFunction(const int n, const AlgebraicVector<T> &obj, const bool scaled = false);
template<typename T> AlgebraicVector<T> BesselKFunction(const int n, const AlgebraicVector<T> &obj, const bool scaled = false);
template<typename T> AlgebraicVector<T> GammaFunction(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> LogGammaFunction(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> tgamma(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> lgamma(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> PsiFunction(const AlgebraicVector<T> &obj, const unsigned int n);
/** @} */

/** @name Vector Norm & Estimation Functions
 * @{
 */
template<typename T> AlgebraicVector<double> norm(const AlgebraicVector<T> &obj, const unsigned int type = 0);
/** @} */

/** @name Vector Evaluation Functions
 *  @{
 */
template<typename T, typename U> U eval(const AlgebraicVector<T> &obj, const U &args);
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const std::initializer_list<U> l);
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const U args[], const unsigned int length);
template<typename T, typename U> AlgebraicVector<U> evalScalar(const AlgebraicVector<T> &obj, const U &arg);
template<typename T> compiledDA compile(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> plug(const AlgebraicVector<T> &obj, const unsigned int var, const double val = 0.0);
/** @} */

/** @name Vector Input/Output Functions
 * @{
 */
template<typename T> std::ostream& operator<<(std::ostream &out, const AlgebraicVector<T> &obj);
template<typename T> std::istream& operator>>(std::istream &in, AlgebraicVector<T> &obj);
/** @} */

/** @name Vector Functions
 *  @{
 */
template<typename T, typename U> typename std::common_type_t<T, U> dot(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T, typename U> AlgebraicVector<typename std::common_type_t<T, U>> cross(const AlgebraicVector<T> &obj1, const AlgebraicVector<U> &obj2);
template<typename T> T length(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> normalize(const AlgebraicVector<T> &obj);
/** @} */

// shortcuts for common vector types
typedef AlgebraicVector<DA> vectorDA;           //!< Short for AlgebraicVector<DA>.
typedef AlgebraicVector<double> vectordb;       //!< Short for AlgebraicVector<double>.

}

#endif /* DINAMICA_ALGEBRAICVECTOR_H_ */
