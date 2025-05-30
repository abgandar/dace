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
 * Monomial.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>
#include <iomanip>

// DACE classes
#include "dace/config.h"
#include "dace/Monomial.h"
#include "dace/DA.h"

namespace DACE {

/** Create a Monomial object large enough to hold all current monomials.
*/
Monomial::Monomial() : m_jj(DA::getMaxVariables()), m_coeff(0.0) {
}

/** Compute the order of the monomial.
    @return Order of the monomial.
 */
unsigned int Monomial::order() const {
    unsigned int ord = 0;

    for(unsigned int i = 0; i < m_jj.size(); i++)
        ord += m_jj[i];

    return ord;
}

/** Convert monomial to string.
    @return A string representing the monomial in human readable form.
 */
std::string Monomial::toString() const {
    std::ostringstream oss;

    oss << "     I  COEFFICIENT              ORDER EXPONENTS" << std::endl;

    // value and order
    oss << "     1  ";
    oss << std::setiosflags(std::ios::uppercase) << std::setprecision(16) << std::scientific << std::setw(24) << m_coeff;
    oss << std::setw(4) << order() << std::setw(1) << ' ';

    // exponents
    for(unsigned int i = 0; i < m_jj.size(); i++){
        oss << std::setw(1) << ' ' << std::setw(2) << m_jj[i];}

    oss << std::endl;
    oss << "------------------------------------------------" << std::endl;

    return oss.str();
}

/** Monomial stream output operator.
    @param[in] out C++ output stream.
    @param[in] m The Monomial to be printed to the stream.
    @return The C++ output stream.
 */
std::ostream& operator<<(std::ostream &out, const Monomial &m) {
    return out << m.toString();
}

}
