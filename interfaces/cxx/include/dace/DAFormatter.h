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
 * DAFormatter.h
 *
 *  Created on: Oct 18, 2014
 *      Author: Dinamica Srl
 */

/*  Formatters allow printing of DAs in a variety of formats other than the standard output

    The standard output format of DA via the C++ iostream interface provides a human-readable text
    representation of a DA vector. In some circumstances it is preferable to print a DA vector in
    a different format, e.g. to output the polynomial in a particular programming language or in
    LaTeX format.

    The abstract DAFormatter class provides an abstract interface for this functionality.
    The DASimpleFormatter class is an implementation of a simple formatter that can output DAs
    in a variety of pre-defined or user-supplied formats including Python, C, Matlab, Fortran,
    and LaTeX code.
*/

#ifndef DINAMICA_DAFORMATTER_H_
#define DINAMICA_DAFORMATTER_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <string>

/** @addtogroup DACECXX C++ Interface
    @{
 */

namespace DACE {

// forward declaration
class DA;

/** Abstract class providing an interface to output DA vectors in advanced formats.

    The standard output format of DA via the C++ iostream interface provides a human-readable text
    representation of a DA vector. In some circumstances it is preferable to print a DA vector in
    a different format, e.g. to output the polynomial in a particular programming language or in
    LaTeX format.

    This abstract class provides an interface for this functionality.
    The DACE::DADefaultFormatter and DACE::DASimpleFormatter classes are an implementation of
    simple formatters that can output DAs in various pre-defined or user-supplied formats
    including Python, C, Matlab, Fortran, and LaTeX code.

    @see DACE::DADefaultFormatter
    @see DACE::DASimpleFormatter
 */
class DACE_API DAFormatter
{
public:
    /** Format a single DA and return a string representation.
        @param da DA object
        @return Formatted string representation.
        @pure This function must be implemented by any derived DA formatter at a minimum.
     */
    virtual std::string operator()(const DA &da) = 0;

    virtual std::string operator()(const std::vector<DA> &da);

    /** Format a single DA and return a string representation.
        @param da DA object
        @return Formatted string representation.
        @deprecated Replaced by DAFormatter::operator().
        @see DAFormatter::operator()
     */
    std::string format(const DA &da) {
        return (*this)(da);
    };

    /** Format a vector of DAs and return a string representation.
        Usually this just formats each DA one after the other.
        @param da Vector of DA objects
        @return Formatted string representation.
        @deprecated Replaced by DAFormatter::operator().
        @see DAFormatter::operator()
     */
    std::string format(const std::vector<DA> &da) {
        return (*this)(da);
    };
};

/** DADefaultFormatter formats DAs using the built-in DA output function.
 */
class DACE_API DADefaultFormatter : public DAFormatter
{
public:
    std::string operator()(const DA &da);
};

/** Structure containing the elements of a simple format as used by the DASimpleFormatter.

    Each monomial consisting of coefficient C and exponents e_1, e_2, ... is formatted
    by outputting these values:

        pos/neg << abs(C)

    followed by this for each exponent e_i

        mul << pre_pow << var << pre_var << i+first_var << post_var << pow << e_i+first_pow << post_pow

    or if shortening is enabled and e_i == 1

        mul << var << pre_var << i+first_var << post_var

    @see DASimpleFormatter
*/
struct DASimpleFormat {
    std::string pos, neg, mul, pre_pow, var, pre_var, post_var, pow, post_pow, linebreak;
    int first_var,              ///< number of the first independent DA variable (e.g. start counting at 0 or 1)
        first_pow;              ///< offset added to the power of an independent DA variable
    unsigned int monperline;    ///< number of monomials formatted per line before a line break is output
    bool shorten;               ///< if true, exponents equal to one are output using the short format
};

/** Formats a DA vector using simple rules to output code suitable for various programming languages.

    Example:
    @code
        DASimpleFormatter sf(DASimpleFormatter::LATEX);
        DA x = sin(DA(1));
        std::string res = sf(x);
        // res now is "+1 \cdot x_{1} -0.1666666666666667 \cdot x_{1}^{3} +0.008333333333333333 \cdot x_{1}^{5}"
    @endcode

    @see DASimpleFormat
 */
class DACE_API DASimpleFormatter : public DAFormatter
{
public:
    static const DASimpleFormat C;
    static const DASimpleFormat C_POW;
    static const DASimpleFormat FORTRAN;
    static const DASimpleFormat FORTRAN_POW;
    static const DASimpleFormat MATLAB;
    static const DASimpleFormat MATLAB_POW;
    static const DASimpleFormat PYTHON;
    static const DASimpleFormat PYTHON_NP;
    static const DASimpleFormat PYTHON_POW;
    static const DASimpleFormat LATEX;

    DASimpleFormat sf;                                              //!< Format the instance uses

    DASimpleFormatter() : sf(C) {};                                 //!< Default constructor formats as C code
    DASimpleFormatter(const DASimpleFormat &isf) : sf(isf) {};      //!< Create formatter using the provided format

    std::string operator()(const DA &da);
};

}

#endif /* DINAMICA_DAFORMATTER_H_ */

/** @} */
