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
 * DAFormatter.cpp
 *
 *  Created on: Oct 20, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>

// DACE classes
#include "dace/config.h"
#include "dace/DAFormatter.h"
#include "dace/Monomial.h"
#include "dace/DA.h"

namespace DACE {

/// C code formatter. Expects an array p[var][n] containing the n-th power of independent variable var.
const DASimpleFormat DASimpleFormatter::C =           { "+",  "-",  "*",        "",     "p", "[",  "",  "][",  "]", " \\\n\t",     0, -1, 20, false };
/// C code formatter. Uses the pow function repeatedly and only expects an array x[var] containing the value of independent variable var.
const DASimpleFormat DASimpleFormatter::C_POW =       { "+",  "-",  "*",        "pow(", "x", "[",  "]", ",",   ")", " \\\n\t",     0,  0, 20, true  };
/// Fortran code formatter. Expects an array p(var,n) containing the n-th power of independent variable var.
const DASimpleFormat DASimpleFormatter::FORTRAN =     { "+",  "-",  "*",        "",     "p", "(",  "",  ",",   ")", " &\n     &",  1,  0, 20, false };
/// Fortran code formatter. Uses ** repeatedly and only expects an array x(var) containing the value of independent variable var.
const DASimpleFormat DASimpleFormatter::FORTRAN_POW = { "+",  "-",  "*",        "",     "x", "(",  ")", "**(", ")", " &\n     &",  1,  0, 20, true  };
/// Matlab code formatter. Expects an array p(var,n) containing the n-th power of independent variable var.
const DASimpleFormat DASimpleFormatter::MATLAB =      { "+",  "-",  ".*",       "",     "p", "(",  "",  ",",   ")", " ...\n\t",    1,  0, 20, false };
/// Matlab code formatter. Uses ^ repeatedly and only expects an array x(var) containing the value of independent variable var.
const DASimpleFormat DASimpleFormatter::MATLAB_POW =  { "+",  "-",  ".*",       "",     "x", "(",  ")", ".^(", ")", " ...\n\t",    1,  0, 20, true  };
/// Python code formatter. Expects a list p[var][n] containing the n-th power of independent variable var.
const DASimpleFormat DASimpleFormatter::PYTHON =      { "+",  "-",  "*",        "",     "p", "[",  "",  "][",  "]", " \\\n\t",     0, -1, 20, false };
/// Python code formatter. Expects a numpy array p[var,n] containing the n-th power of independent variable var.
const DASimpleFormat DASimpleFormatter::PYTHON_NP =   { "+",  "-",  "*",        "",     "p", "[",  "",  ",",   "]", " \\\n\t",     0, -1, 20, false };
/// Python code formatter. Uses ** repeatedly and only expects a list x[var] containing the value of independent variable var.
const DASimpleFormat DASimpleFormatter::PYTHON_POW =  { "+",  "-",  "*",        "",     "x", "[",  "]", "**(", ")", " \\\n\t",     0,  0, 20, true  };
/// LaTeX formatter. Outputs the polynomial as a nicely formatted equation.
const DASimpleFormat DASimpleFormatter::LATEX =       { " +", " -", " \\cdot ", "",     "x", "_{", "}", "^{",  "}", " \n\t",       1,  0, 20, true  };

/** Format a single DA and return a string representation.
    @param da DA object
    @return formatted string representation
 */
std::string DASimpleFormatter::format(const DA &da){
    const std::vector<Monomial> monomials = da.getMonomials();
    const size_t size = monomials.size();
    std::ostringstream res;

    res.precision(16);
    for(size_t i=0; i<size; i++) {
        if(monomials[i].m_coeff < 0)
            res << sf.neg << -monomials[i].m_coeff;
        else
            res << sf.pos << monomials[i].m_coeff;

        for(size_t j=0; j<monomials[i].m_jj.size(); j++) {
            if(monomials[i].m_jj[j] <= 0)
                continue;
            else if(sf.shorten && monomials[i].m_jj[j] == 1)
                res << sf.mul << sf.var << sf.pre_var << j+sf.first_var << sf.post_var;
            else
                res << sf.mul << sf.pre_pow << sf.var << sf.pre_var << j+sf.first_var << sf.post_var << sf.pow << monomials[i].m_jj[j]+sf.first_pow << sf.post_pow;
        }
        if((i+1)%sf.monperline == 0 && i+1<size)
            res << sf.linebreak;
    }

    return res.str();
}

/** Format a vector of DAs and return a string representation.
    This just formats each DA in the vector one after the other.
    @param da vector of DA objects
    @return formatted string representation
 */
std::string DASimpleFormatter::format(const std::vector<DA> &da){
    std::ostringstream res;

    for(unsigned int i=0; i<da.size(); i++)
        res << format(da[i]) << std::endl;

    return res.str();
}

}
