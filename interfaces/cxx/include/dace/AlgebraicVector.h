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

#ifndef DINAMICA_ALGEBRAICVECTOR_H_
#define DINAMICA_ALGEBRAICVECTOR_H_

// C++ stdlib classes required for interface definition
#include <vector>
#include <initializer_list>

// DACE classes required for interface definition (DA.h needed for DA::getMaxOrder(), DA::getMaxVariables() default arguments)
#include "dace/PromotionTrait.h"
#include "dace/DA.h"

namespace DACE {

// forward declarations
#ifdef WITH_ALGEBRAICMATRIX
template<typename T> class AlgebraicMatrix;
#endif

/*! Generic vector class to handle vectors of algebraic types and their algebraic operations. */
template<typename T> class AlgebraicVector : public std::vector<T>
{
public:
    /***********************************************************************************
    *     Constructors
    ************************************************************************************/
    /*! Default constructor to create empty AlgebraicVector
    */
    AlgebraicVector() : std::vector<T>() {};
    /*! Constructor with size to allocate a vector of the given size with elements initialized using their default constructor.
        \param[in] size initial length of the AlgebraicVector.
    */
    explicit AlgebraicVector(const size_t size) : std::vector<T>(size) {};
    /*! Constructor with size and value to allocate a vector of the given size with elements initialized as copies of d.
        \param[in] size initial length of the AlgebraicVector.
        \param[in] d    initial value for the elements
    */
    AlgebraicVector(const size_t size, const T &d) : std::vector<T>(size, d) {};
    /*! Copy constructor to create a copy of any existing vector.
        \param[in] v vector to be copied into AlgebraicVector
    */
    AlgebraicVector(const std::vector<T> &v) : std::vector<T>(v) {};
    /*! Extraction constructor to copy only a given range of elements from vector v.
        \param[in] v vector to be copied into AlgebraicVector
        \param[in] first index of the first element to be copied
        \param[in] last index of the last element to be copied
        \note The constructor does not perform any range checking for the extraction.
        \sa AlgebraicVector<T>::extract
    */
    AlgebraicVector(const std::vector<T> &v, const size_t first, const size_t last) : std::vector<T>(v.begin()+first, v.begin()+last+1) {};
    /*! Constructor to create a vector from an initializer list.
        \param[in] l braced initializer list to be copied into the AlgebraicVector
    */
    AlgebraicVector(std::initializer_list<T> l) : std::vector<T>(l) {};

    /***********************************************************************************
    *     Element and coefficient access / extraction routines
    ************************************************************************************/
    AlgebraicVector<T> extract(const size_t first, const size_t last) const;
    template<typename U> AlgebraicVector<typename PromotionTrait<T,U>::returnType> concat(const std::vector<U> &obj) const;
    AlgebraicVector<double> cons() const;
#ifdef WITH_ALGEBRAICMATRIX
    AlgebraicMatrix<double> linear() const;
#else
    std::vector< std::vector<double> > linear() const;
#endif /* WITH_ALGEBRAICMATRIX */

    /***********************************************************************************
    *     Operator overloads
    ************************************************************************************/
    AlgebraicVector<T> operator-() const;
    template<typename U> AlgebraicVector<T>& operator+=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator+=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator-=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator*=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const AlgebraicVector<U> &obj);
    template<typename U> AlgebraicVector<T>& operator/=(const U &obj);
    template<typename U> AlgebraicVector<T>& operator<<(const std::vector<U> &obj);

    /***********************************************************************************
    *     Math routines
    ************************************************************************************/
    AlgebraicVector<T> absolute() const;
    AlgebraicVector<T> trunc() const;
    AlgebraicVector<T> round() const;
    AlgebraicVector<T> mod(const double p) const;
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

    /***********************************************************************************
    *    Vector routines
    ************************************************************************************/
    template<typename V> typename PromotionTrait<T,V>::returnType dot(const AlgebraicVector<V> &obj) const;
    template<typename V> AlgebraicVector<typename PromotionTrait<T,V>::returnType> cross(const AlgebraicVector<V> &obj) const;
    T length() const;
    AlgebraicVector<T> normalize() const;
    // XXX: various Jacobians, gradients, curls, etc?

    /***********************************************************************************
    *     Special routines (DA related)
    ************************************************************************************/
    AlgebraicVector<T> deriv(const unsigned int p) const;
    AlgebraicVector<T> integ(const unsigned int p) const;
    template<typename V> V eval(const V &args) const;
    template<typename U> AlgebraicVector<U> eval(const std::initializer_list<U> l) const;
    template<typename U> AlgebraicVector<U> evalScalar(const U &arg) const;
    compiledDA compile() const;
    AlgebraicVector<T> plug(const unsigned int var, const double val = 0.0) const;
    AlgebraicVector<T> trim(const unsigned int min, const unsigned int max = DA::getMaxOrder()) const;
    AlgebraicVector<T> invert() const;

    /********************************************************************************
    *     DA norm routines
    *********************************************************************************/
    AlgebraicVector<double> norm(const unsigned int type = 0) const;
    /* XXX: define and add the norm estimation routines from DA including convergence radius estimation
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DA::getMaxOrder()) const;
    Interval bound() const;
    double convRadius(const double eps, const unsigned int type = 1) const;
    */

    /********************************************************************************
    *     Static factory routines
    *********************************************************************************/
    static AlgebraicVector<DA> identity(const size_t n = DA::getMaxVariables());

    /***********************************************************************************
    *     Input/Output routines
    ************************************************************************************/
    std::string toString() const;

private:
#ifndef WITH_ALGEBRAICMATRIX
    static void matrix_inverse(std::vector< std::vector<double> > &A);        // Private helper routine for double precision matrix inversion
#endif /* WITH_ALGEBRAICMATRIX */
};

// operators
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator+( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator-( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator*( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const AlgebraicVector<U> &obj1, const V &obj2);
template<typename U,typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> operator/( const U &obj1, const AlgebraicVector<V> &obj2);

template<typename U> std::ostream& operator<<(std::ostream &out, const AlgebraicVector<U> &obj);
template<typename U> std::istream& operator>>(std::istream &in, AlgebraicVector<U> &obj);

// Declaration of external functional style wrappers to access AlgebraicVector functions
template<typename T> AlgebraicVector<double> cons(const AlgebraicVector<T> &obj);
#ifdef WITH_ALGEBRAICMATRIX
template<typename T> AlgebraicMatrix<double> linear(const AlgebraicVector<T> &obj);
#else
template<typename T> std::vector< std::vector<double> > linear(const AlgebraicVector<T> &obj);
#endif /* WITH_ALGEBRAICMATRIX */
template<typename T> AlgebraicVector<T> deriv(const AlgebraicVector<T> &obj, const unsigned int p);
template<typename T> AlgebraicVector<T> integ(const AlgebraicVector<T> &obj, const unsigned int p);
template<typename T> AlgebraicVector<T> absolute(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> trunc(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> round(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> mod(const AlgebraicVector<T> &obj, const double p);
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const double p);
template<typename T> AlgebraicVector<T> pow(const AlgebraicVector<T> &obj, const int p);
template<typename T> AlgebraicVector<T> root(const AlgebraicVector<T> &obj, const int p = 2);
template<typename T> AlgebraicVector<T> minv(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqr(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> sqrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> isrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> cbrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> icbrt(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> hypot(const AlgebraicVector<T> &Y);
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
template<typename U, typename V> typename PromotionTrait<U,V>::returnType dot(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename U, typename V> AlgebraicVector<typename PromotionTrait<U,V>::returnType> cross(const AlgebraicVector<U> &obj1, const AlgebraicVector<V> &obj2);
template<typename T> T length(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> normalize(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> trim(const AlgebraicVector<T> &obj, unsigned int min, unsigned int max = DA::getMaxOrder());
template<typename T, typename V> V eval(const AlgebraicVector<T> &obj, const V &args);
template<typename T, typename U> AlgebraicVector<U> eval(const AlgebraicVector<T> &obj, const std::initializer_list<U> l);
template<typename T, typename U> AlgebraicVector<U> evalScalar(const AlgebraicVector<T> &obj, const U &arg);
template<typename T> compiledDA compile(const AlgebraicVector<T> &obj);
template<typename T> AlgebraicVector<T> plug(const AlgebraicVector<T> &obj, const unsigned int var, const double val = 0.0);
template<typename T> AlgebraicVector<double> norm(const AlgebraicVector<T> &obj, const unsigned int type = 0);

// specializations for various DA specific routines implemented and instantiated directly in the library instead of in a template
#ifdef WITH_ALGEBRAICMATRIX
template<> DACE_API AlgebraicMatrix<double> AlgebraicVector<DA>::linear() const;
template<> DACE_API AlgebraicMatrix<double> linear(const AlgebraicVector<DA> &obj);
#else
template<> DACE_API std::vector<std::vector<double>> AlgebraicVector<DA>::linear() const;
template<> DACE_API void AlgebraicVector<DA>::matrix_inverse(std::vector< std::vector<double> > &A);
template<> DACE_API std::vector< std::vector<double> > linear(const AlgebraicVector<DA> &obj);
#endif /* WITH_ALGEBRAICMATRIX */
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::trim(const unsigned int min, const unsigned int max) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::deriv(const unsigned int p) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::integ(const unsigned int p) const;
template<> DACE_API compiledDA AlgebraicVector<DA>::compile() const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::plug(const unsigned int var, const double val) const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::invert() const;
template<> DACE_API AlgebraicVector<DA> AlgebraicVector<DA>::identity(const size_t n);
template<> DACE_API AlgebraicVector<DA> trim(const AlgebraicVector<DA> &obj, unsigned int min, unsigned int max);
template<> DACE_API AlgebraicVector<DA> deriv(const AlgebraicVector<DA> &obj, const unsigned int p);
template<> DACE_API AlgebraicVector<DA> integ(const AlgebraicVector<DA> &obj, const unsigned int p);
template<> DACE_API compiledDA compile(const AlgebraicVector<DA> &obj);
template<> DACE_API AlgebraicVector<DA> plug(const AlgebraicVector<DA> &obj, const unsigned int var, const double val);

// shortcuts for common vector types
typedef AlgebraicVector<DA> vectorDA;       //!< Shorthand notation for AlgebraicVector<DA>.
typedef AlgebraicVector<double> vectordb;   //!< Shorthand notation for AlgebraicVector<double>.

}

#endif /* DINAMICA_ALGEBRAICVECTOR_H_ */
