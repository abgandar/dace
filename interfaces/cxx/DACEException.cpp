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
 * DACEException.cpp
 *
 *  Created on: Mar 11, 2014
 *      Author: Dinamica Srl
 */

// C++ stdlib classes used only internally in the implementation
#include <sstream>
#include <iostream>

// DACE classes
#include "dace/config.h"
#include "dace/DACEException.h"
#include "dace/dacecore.h"

namespace DACE {

    int DACEException::severity = 6;
    bool DACEException::warning = true;

    /********************************************************************************
    *     Constructors & Destructors
    *********************************************************************************/
    /** Create DACEException from current error code in DACE core.
        This also clears the current DACE core error code. It then executes the
        appropriate action for the exception based on current settings.
     */
    DACEException::DACEException() {
        m_x = daceGetErrorX();
        m_yy = daceGetErrorYY();
        updateMessage();
        daceClearError();
        execute();
    }

    /** Create and execute a DACEException object with given severity and ID codes.
        @param exc_sv The severity code of the error plus 10 to indicate it originated in the C++ interface.
        @param exc_id The ID code of the error.
     */
    DACEException::DACEException(const int exc_sv, const int exc_id) {
        m_x = exc_sv;
        m_yy = exc_id;
        updateMessage();
        execute();
    }

    /** Destructor.
     */
    DACEException::~DACEException() throw() {
        // nothing to do, just overloading the virtual destructor of the parent class
    }

    /********************************************************************************
    *     Private member functions
    *********************************************************************************/
    /** Update the error message of this exception based on its ID.
     */
    void DACEException::updateMessage() {
        struct errstrings {
            int ID;
            const char* msg;
        };

        // error codes known to the DACE C++ interface
        static const errstrings DACEerr[] = {
            { 000, "DACE: Unknown DACE C++ interface error. Contact DACE Developers for filing a bug report."},
            {1101, "DA::getCoeff: Not enough exponents, missing exponents treated as zero"},
            {1102, "DA::setCoeff: Not enough exponents, missing exponents treated as zero"},
            {1103, "DA::getCoeff: More exponents than variables, ignoring extra exponents"},
            {1104, "DA::setCoeff: More exponents than variables, ignoring extra exponents"},
            {1604, "compiledDA::compiledDA(): Dimension lower than 1"},
            {1605, "compiledDA::eval: Argument size lower than the number of variables in the polynomial"},
            {1506, "storedDA::operator DA(): Invalid data, can't convert to DA"},
            {2099, "DA::checkVersion: DACE C++ interface header file and DACE core library version mismatch"}
        };
        static const int length = sizeof(DACEerr)/sizeof(errstrings);

        const int id = m_x*100 + m_yy;
        std::stringstream s;

        // severity above 10 is C++ interface error (subtract 10 to get severity)
        if(m_x > 10)
        {
            int i;
            for(i = length-1; (i > 0) && (DACEerr[i].ID != id); i--);
            s << DACEerr[i].msg << " (ID: " << DACEerr[i].ID << ")" ;
        }
        else
        {
            s << daceGetErrorMessage() << " (ID: " << id << ")";
        }
        msg = s.str();
    }

    /** Execute this exception, i.e. throw or print warning based on current settings.
        @throw DACE::DACEException
     */
    void DACEException::execute() const {
        const int sev = m_x%11; // modulo 11 to handle both DACE core and C++ interface severity codes

        if(sev >= severity) {
            throw *this;
        } else if(warning) {
            std::cerr << "Warning: " << msg << std::endl;
        }
    }

    /********************************************************************************
    *     Public member functions
    *********************************************************************************/
    /** Return a human readable error string representing this exception.
        @return A C string containing the error message.
     */
    const char* DACEException::what() const throw() {
        return msg.c_str();
    }

    /********************************************************************************
    *     Static member functions
    *********************************************************************************/
    /** Set severity level.
        Errors with severity code greater or equal to this
        value will throw an exception.
        Severity levels are:\n
        0 = Warning:   Informative, no action required\n
        1 = Warning:   Serius, possible wrong implementation\n
        6 = Error:     Recoverable, assumptions have been made\n
        9 = Error:     Unrecoverable, new call to DA::init is required to
                       reinitialize DACE, interface objects are no longer valid\n
        10 = Critical: Crash in the DACE, just printing as much as possible
                       and dying.
        @param n The new severity level value.
     */
    void DACEException::setSeverity(const int n) {
        severity = n;
    }

    /** Set the current mode for printing of warnings.
        @param w The new warning status (prints warnings if true).
     */
    void DACEException::setWarning(const bool w) {
        warning = w;
    }

    /********************************************************************************
    *     Related functions
    *********************************************************************************/
    /** DACEException stream output operator.
        @param[in] out A C++ output stream.
        @param[in] ex The DACEException to be printed to the stream.
        @return The C++ output stream.
     */
    std::ostream& operator<<(std::ostream &out, const DACEException &ex) {
        return out << ex.what() << std::endl;
    }

}
