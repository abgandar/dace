/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2019 University of Southampton                                    *
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
 *  dacecontrib.h
 *
 *  Created on: November 8, 2019
 *      Author: Alexander Wittig
 */

/// @cond CONTRIB

/*  Internal contributed functions to DACE core.

    This file contains all externally contributed routines used in the DACE core.
    It is never included publicly by DACE users or high level interfaces.
*/

#ifndef DINAMICA_DACECONTRIB_H_
#define DINAMICA_DACECONTRIB_H_

/** @addtogroup DACEContrib
    @{
 */
double zeta_(const double x, const double q, unsigned int *err);
double dgamma_(const double *x);
double psi_(const double *x);
int ribesl_(double *x, double *alpha, long int *nb, long int *ize, double *b, long int *ncalc);
int rjbesl_(double *x, double *alpha, long int *nb, double *b, long int *ncalc);
int rkbesl_(double *x, double *alpha, long int *nb, long int *ize, double *b, long int *ncalc);
int rybesl_(double *x, double *alpha, long int *nb, double *b, long int *ncalc);
/** @} */

#endif /* DINAMICA_DACECONTRIB_H_ */

/// @endcond
