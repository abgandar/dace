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
 * AlgebraicMatrix.h
 *
 *  Created on: July 17, 2014
 *      Author: Dinamica Srl
 */

/*  Experimental AlgebraicMatrix class.

    This header file contains the AlgebraicMatrix class. This is experimental and
    not supported by default. It may have bugs or not work as expected.
*/

#ifndef DINAMICA_ALGEBRAICMATRIX_H_
#define DINAMICA_ALGEBRAICMATRIX_H_

// C++ stdlib classes required for interface definition
#include <vector>
#include <iostream>
#include <type_traits>

namespace DACE {

// forward declarations
class DA;
template<typename T> class AlgebraicVector;

/** Matrix of any algebraic type. */
template <class T> class AlgebraicMatrix
{
public:
    /***********************************************************************************
    *     Constructors & Destructors
    ************************************************************************************/
    AlgebraicMatrix() : _nrows(0), _ncols(0) {};    //!< Default Constructor

    /** Constructor for square matrices.
        @param[in] size size of the matrix, i.e. the number of rows and columns.
     */
    explicit AlgebraicMatrix(const int size) : _nrows(size), _ncols(size), _data(size*size,0.0) { };

    /** Constructor for rectangular matrices.
        @param[in] nrows number of rows of the matrix
        @param[in] ncols number of columns of the matrix
     */
    AlgebraicMatrix(const int nrows, const int ncols) : _nrows(nrows), _ncols(ncols), _data(nrows*ncols,0.0) { };

    /** Constructor for rectangular matrices that allows to set all elements equal to a variable.
        @param[in] nrows number of rows of the matrix
        @param[in] ncols number of columns of the matrix
        @param[in] d     matrix elements value
     */
    AlgebraicMatrix(const int nrows, const int ncols, const T &d) : _nrows(nrows), _ncols(ncols), _data(nrows*ncols, d) { };

    /***********************************************************************************
    *     Output number of rows, columns, and size
    ************************************************************************************/
    /** Returns the number of columns of the matrix
        @return number of columns of the matrix.
     */
    unsigned int ncols() const { return this->_ncols; };

    /** Returns the number of rows of the matrix
        @return number of rows of the matrix.
     */
    unsigned int nrows() const { return this->_nrows; };

    /** Returns the number of elements of the matrix
        @return number of elements of the matrix.
     */
    unsigned int size() const { return this->_data.size(); };

    void resize(int size);
    void resize(int rows, int cols);

    /***********************************************************************************
    *     Element access routine
    ************************************************************************************/
    T& at(const unsigned int irow, const unsigned int icol);
    const T& at(const unsigned int irow, const unsigned int icol) const;

    std::vector<T> getrow(const unsigned int irow) const;
    std::vector<T> getcol(const unsigned int icol) const;

    void setrow(const unsigned int irow, const std::vector<T> &obj);
    void setcol(const unsigned int icol, const std::vector<T> &obj);

    AlgebraicMatrix<T> submat(const unsigned int first_row, const unsigned int first_col, const unsigned int last_row, const unsigned int last_col) const;
    AlgebraicMatrix<T> submat(const unsigned int last_row, const unsigned int last_col) const;

    /***********************************************************************************
    *     Matrix operations
    ************************************************************************************/
    AlgebraicMatrix<T> transpose() const;
    T det() const;
    AlgebraicMatrix<T> inv() const;

    /***********************************************************************************
    *     Coefficient access routines
    ************************************************************************************/
    AlgebraicMatrix<double> cons() const;

private:
    unsigned int _nrows;    //!< Number of rows of the matrix
    unsigned int _ncols;    //!< Number of columns of the matrix
    std::vector<T> _data;   //!< Container for storing the elements linearly

    static unsigned int pivot(unsigned int& k, const unsigned int ii, const AlgebraicMatrix<T>& A, std::vector<unsigned int>& P, std::vector<unsigned int>& R, std::vector<unsigned int>& C1, std::vector<unsigned int >& C2, T& det);
    static void eliminate(const unsigned int k, AlgebraicMatrix<T>& A, std::vector<unsigned int>& R);
};

/***********************************************************************************
 *     Operators
 ************************************************************************************/
template<typename U> std::ostream& operator<<(std::ostream &out, const AlgebraicMatrix<U> &obj);
template<> DACE_API  std::ostream& operator<<(std::ostream &out, const AlgebraicMatrix<DA> &obj);
template<typename U> std::istream& operator>>(std::istream &in, AlgebraicMatrix<U> &obj);
template<> DACE_API  std::istream& operator>>(std::istream &in, AlgebraicMatrix<DA> &obj);

template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator+( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator+( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator+( const U &obj1, const AlgebraicMatrix<V> &obj2 );

template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator-( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator-( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator-( const U &obj1, const AlgebraicMatrix<V> &obj2 );

template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator*( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename std::common_type_t<U, V>> operator*( const U &obj1, const AlgebraicMatrix<V> &obj2 );
template<typename U,typename V> AlgebraicVector<typename std::common_type_t<U, V>> operator*( const AlgebraicVector<U> &obj1, const AlgebraicMatrix<V> &obj2 );
template<typename U,typename V> AlgebraicVector<typename std::common_type_t<U, V>> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicVector<V> &obj2 );

/***********************************************************************************
 *     Functional style wrappers
 ************************************************************************************/
template<class T> AlgebraicMatrix<T> transpose(const AlgebraicMatrix<T> &obj);
template<class T> T det(const AlgebraicMatrix<T> &obj);
template<class T> AlgebraicMatrix<T> inv(const AlgebraicMatrix<T> &obj);
template<class T> AlgebraicMatrix<double> cons(const AlgebraicMatrix<T> &obj);

/***********************************************************************************
 *     Type definitions
 ************************************************************************************/
typedef AlgebraicMatrix<DA> matrixDA;           //!< Short for AlgebraicMatrix<DA>
typedef AlgebraicMatrix<double> matrixdb;       //!< Short for AlgebraicMatrix<double>

}

#endif /* DINAMICA_ALGEBRAICMATRIX_H_ */
