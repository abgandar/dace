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
 * dace_s.h
 *
 *  Created on: Jul 20, 2018
 *      Author: University of Southampton
 */

/*! @file

     @brief Public DACE C++ header file for static linking.

     This header file just includes all other DACE C++ headers in the right order.
     It is identical to include/dace/dace.h except on Windows where it explicitly
     sets the linkage of all DACE functions to local. Use this when linking with
     the static version of the DACE library (instead of the DLL).
*/

#ifndef DINAMICA_DACE_S_H_
#define DINAMICA_DACE_S_H_

// explicitly define DACE_API to be empty for static linking
#define DACE_API
#include "dace/dace.h"

#endif /* DINAMICA_DACE_S_H_ */
