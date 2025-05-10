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
 *  daceerror.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

/** @addtogroup DACE Core
    @{
 */

// indicate that we want to use strXXX_s functions, if available
#define __STDC_WANT_LIB_EXT1__ 1
#include <string.h>
#include <stdio.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"

/********************************************************************************
 *     DACE error messages and their IDs
 ********************************************************************************/
/* Error ID and associated human readable error message */
typedef struct {
    int ID;             // Internal ID of the error
    const char* msg;    // Human readable error message
} errstrings;

/* List of all known errors and their ID (100*severity + error code) */
static const errstrings DACEerr[] = {
    {   0, "Unknown DACE error. Contact DACE developers for filing a bug report."},
    {1001, "Dynamic memory allocation failure"},
    {1002, "Out of memory"},
    {1003, "DACE has not been initialized"},
    {1004, "DA object not allocated"},
    {1005, "Incorrect number of monomials"},
    {1006, "Incorrect DA coding arrays"},
    {1007, "Requested length too long"},
    {1008, "Error in monomial evaluation tree construction"},
    {   9, ""},
    {  10, ""},
    { 911, "Order and/or variable too large"},
    {  12, ""},
    {  13, ""},
    {  14, ""},
    {  15, ""},
    {  16, ""},
    {  17, ""},
    {  18, ""},
    {  19, ""},
    {  20, ""},
    { 621, "Not enough storage"},
    { 622, "Order too large"},
    { 623, "Number of variables too high"},
    { 624, "Invalid independent variable"},
    { 625, "Invalid DA codes"},
    { 626, "Invalid encoded exponent"},
    {  27, ""},
    {  28, ""},
    {  29, ""},
    {  30, ""},
    { 631, "Invalid data"},
    { 632, "Unknown format"},
    { 633, "DA vector too long"},
    { 634, "Not enough lines to read"},
    {  35, ""},
    {  36, ""},
    {  37, ""},
    {  38, ""},
    {  39, ""},
    {  40, ""},
    { 641, "Dividing by zero"},
    { 642, "Inverse does not exists"},
    { 643, "Non-integer power of non-positive DA"},
    { 644, "Zero-th root does not exist"},
    { 645, "Even root of negative DA"},
    { 646, "Odd root of zero DA"},
    { 647, "Negative constant part in logarithm"},
    { 648, "Base of logarithm must be positive"},
    { 649, "Cosine is zero in tangent"},
    { 650, "Out of domain"},
    { 651, "No estimate is possible"},
    {  52, ""},
    {  53, ""},
    {  54, ""},
    {  55, ""},
    {  56, ""},
    {  57, ""},
    {  58, ""},
    {  59, ""},
    {  60, ""},
    { 161, "Free or invalid variable"},
    { 162, "Truncation order too high"},
    { 163, "Inacurate estimate"},
    { 164, "Numbering out of order"},
    { 165, "Too many variables"},
    { 166, "Duplicate monomial"},
    { 167, "Order increased to 1"},
    { 168, "Variable increased to 1"},
    {  69, ""},
    {  70, ""},
    {  71, ""},
    {  72, ""},
    {  73, ""},
    {  74, ""},
    {  75, ""},
    {  76, ""},
    {  77, ""},
    {  78, ""},
    {  79, ""},
    {  80, ""},
    {  81, ""},
    {  82, ""},
    {  83, ""},
    {  84, ""},
    {  85, ""},
    {  86, ""},
    {  87, ""},
    {  88, ""},
    {  89, ""},
    {  90, ""},
    {  91, ""},
    {  92, ""},
    {  93, ""},
    {  94, ""},
    {  95, ""},
    {  96, ""},
    {  97, ""},
    {  98, ""},
    {  99, ""}
};

/********************************************************************************
 *     DACE error state routine
 ********************************************************************************/

/** Return the current error ID.
    The error ID is 100*severity + error code, where severity is between 0 - 10
    and the error code between 0 - 99.
    @return The current error ID
*/
unsigned int daceGetError()
{
    return DACEDbg.ierr;
}

/** Return current error severity (between 0 - 10).
    @return The current error severity
*/
unsigned int daceGetErrorX()
{
    return DACEDbg.ixerr;
}

/** Return the current error code (between 0 - 99).
    @return The current error code
*/
unsigned int daceGetErrorYY()
{
    return DACEDbg.iyyerr;
}

/** Return the function name where current error occured.
    @return The function name.
*/
const char* daceGetErrorFunctionName()
{
    return DACEDbg.name;
}

/** Return the current error message.
    @return The human readable error message for current error.
*/
const char* daceGetErrorMessage()
{
    return DACEDbg.msg;
}

/** Clear the current DACE error code.
*/
void daceClearError()
{
    DACEDbg.ierr = 0;
    DACEDbg.ixerr = 0;
    DACEDbg.iyyerr = 0;
    *DACEDbg.name = '\0';
    *DACEDbg.msg = '\0';
}

/**   Set DACE error state for errors within the DACE.

      The error codes are defined as XYY with X indicating the severity and
      YY corresponding to the actual error code

      Severity Levels X

      1  = Info:  Informative, no action required

      3  = Warning:  Serious, possibly incorrect use of DACE routines

      6  = Error:    Recoverable, result may not be correct or assumptions have
                     been made

      9  = Error:    Unrecoverable, new call to daceInitialize is required to
                     reinitialize DACE, DACE objects are no longer valid

      10 = Critical: Crash in the DACE, just printing as much as possible
                     and die.

      @param[in] c name of function where the error happened
      @param[in] ix is the error severity code
      @param[in] iyy is the error code
 */

void daceSetError(const char *c, const unsigned int ix, const unsigned int iyy)
{
    // check if it is a critical error
    if(ix == DACE_PANIC)
    {
        fprintf(stderr, "DACE critical error %u in %s:\n%s\nbye bye!\n", DACEerr[iyy%100].ID, c, DACEerr[iyy%100].msg);
        exit(1);
    }
    else
    {
        if(ix > DACEDbg.ixerr)
        {
            DACEDbg.ierr = (ix%11)*100 + iyy%100;
            DACEDbg.ixerr = ix%11;
            DACEDbg.iyyerr = iyy%100;
#ifdef HAVE_SAFE_STRINGS
            strncpy_s(DACEDbg.name, ERROR_FUN_SIZE, c, ERROR_FUN_SIZE-1);
            strncpy_s(DACEDbg.msg, ERROR_MSG_SIZE, c, ERROR_MSG_SIZE-1);
            strncat_s(DACEDbg.msg, ERROR_MSG_SIZE, ": ", ERROR_MSG_SIZE-strnlen_s(DACEDbg.msg, ERROR_MSG_SIZE)-1);
            strncat_s(DACEDbg.msg, ERROR_MSG_SIZE, DACEerr[DACEDbg.iyyerr].msg, ERROR_MSG_SIZE-strnlen_s(DACEDbg.msg, ERROR_MSG_SIZE)-1);
#else
            strncpy(DACEDbg.name, c, ERROR_FUN_SIZE-1); DACEDbg.name[ERROR_FUN_SIZE-1] = '\0';
            strncpy(DACEDbg.msg, c, ERROR_MSG_SIZE-1); DACEDbg.msg[ERROR_MSG_SIZE-1] = '\0';
            strncat(DACEDbg.msg, ": ", ERROR_MSG_SIZE-strnlen(DACEDbg.msg, ERROR_MSG_SIZE)-1);
            strncat(DACEDbg.msg, DACEerr[DACEDbg.iyyerr].msg, ERROR_MSG_SIZE-strnlen(DACEDbg.msg, ERROR_MSG_SIZE)-1);
#endif
        }
    }
}

/** @} */
