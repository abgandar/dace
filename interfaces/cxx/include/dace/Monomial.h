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
 * Monomial.h
 *
 *  Created on: Mar 10, 2014
 *      Author: Dinamica Srl
 */

/*  Simple storage class for a single monomial

    The Monomial class represents a single monomial, containing a vector of exponents
    and the corresponding coefficient. This class is only used for user convenience
    when extracting monomials from DAs and is not related to the internal representation.
*/

#ifndef DINAMICA_MONOMIAL_H_
#define DINAMICA_MONOMIAL_H_

// C++ stdlib classes used in this public interface
#include <vector>
#include <string>
#include <ostream>

/** @addtogroup DACECXX C++ Interface
    @{
 */

namespace DACE {

/** Represents a single monomial.

    This class represents a single monomial, containing a vector of exponents
    and the corresponding coefficient. This class is only used for user convenience
    when extracting monomials from DAs and is not related to the internal representation.
 */
class DACE_API Monomial
{
public:
    std::vector<unsigned int> m_jj;     //!< Vector of exponents
    double m_coeff;                     //!< Coefficient

    Monomial();
    unsigned int order() const;
    std::string toString() const;
};

DACE_API std::ostream& operator<< (std::ostream &out, const Monomial &m);

}

#endif /* DINAMICA_MONOMIAL_H_ */

/** @} */
