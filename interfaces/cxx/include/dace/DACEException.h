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
 * DACEException.h
 *
 *  Created on: Mar 11, 2014
 *      Author: Dinamica Srl
 */

/*  DACE Exceptions are thrown when a DA operations fails within the C++ interface.

    The DACEException class contains methods for error handling within the DACE C++ interface.
    Whenever an error occurs during a DA operation either nothing happens, a warning is printed,
    or this exception is thrown. What happens depends on the settings in this class.
*/

#ifndef DINAMICA_DACEEXCEPTION_H_
#define DINAMICA_DACEEXCEPTION_H_

// C++ stdlib classes used in this public interface
#include <exception>
#include <string>
#include <ostream>

namespace DACE {

/** Contains methods for error handling within the DACE C++ interface.

    Whenever an error occurs during a DA operation either nothing happens, a warning is printed,
    or this exception is thrown. What happens is determined via the static settings in this class.
 */
class DACE_API DACEException : public std::exception
{
private:
    int m_x;                    //!< Severity
    int m_yy;                   //!< Error code
    std::string msg;            //!< Error message
    static int severity;        //!< Default severity (errors equal or larger raise a DACEException)
    static bool warning;        //!< Default warning status (all errors print a warning if true)
    void execute() const;
    void updateMessage();

public:
    DACEException();
    DACEException(const int exc_sv, const int exc_id);
    ~DACEException() throw();

    const char* what() const throw();

    static void setSeverity(const int n);
    static void setWarning(const bool w);

    friend DACE_API std::ostream& operator<< (std::ostream &out, const DACEException &ex);
};

}

#endif /* DINAMICA_DACEEXCEPTION_H_ */
