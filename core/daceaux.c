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
 *  daceaux.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"


/** Return the minimum between two unsigned integers.
    @return Minimum between @e a and @e b.
*/
unsigned int umin(const unsigned int a, const unsigned int b) { return (a > b)? b : a; }

/** Return the maximum between two unsigned integers.
    @return Maximum between @e a and @e b.
*/
unsigned int umax(const unsigned int a, const unsigned int b) { return (a < b)? b : a; }

/** Raise a double to positive integer power.
    @param[in] a The base value.
    @param[in] b The power.
    @return @e a raised to the power of @e b.
 */
double pown(double a, unsigned int b)
{
    double res = 1.0;
    while(b)
    {
        if(b & 1u)
            res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

/** Raise integer to positive integer power.
    @param[in] a The base value.
    @param[in] b The power.
    @return @e a raised to the power of @e b.
 */
int npown(int a, unsigned int b)
{
    int res = 1;
    while(b)
    {
        if(b & 1u)
            res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    /** Wrapper for C calloc function (allocate memory and zero it) with DACE error handling.
        @param[in] count The number of elements to allocate.
        @param[in] size The size (in bytes) of each element.
        @return A pointer to the newly allocated memory.
     */
    void* dacecalloc(size_t count, size_t size)
    {
        void* res = calloc(count, size);
        if(!res)
        {
            daceSetError(__func__, DACE_PANIC, 1);
            exit(1);
        }
        return res;
    }

    /** Wrapper for C malloc function (allocate memory) with DACE error handling.
        @param[in] size The size (in bytes) of memory.
        @return A pointer to the newly allocated memory.
     */
    void* dacemalloc(size_t size)
    {
        void* res = malloc(size);
        if(!res)
        {
            daceSetError(__func__, DACE_PANIC, 1);
            exit(1);
        }
        return res;
    }

    /** Wrapper for C malloc and memset functions (allocate memory and zero it) with DACE error handling.
        @param[in] size The size (in bytes) of memory.
        @return A pointer to the newly allocated memory.
     */
    void* dacemalloc0(size_t size)
    {
        void* res = dacemalloc(size);
        memset(res, 0, size);

        return res;
    }

    /** Wrapper for C realloc function (reallocate memory) with DACE error handling.
        @param[in] ptr A pointer to the currently allocated memory.
        @param[in] size The size (in bytes) of the new memory.
        @return A pointer to the newly reallocated memory.
     */
    void* dacerealloc(void* ptr, size_t size)
    {
        void* res = realloc(ptr, size);
        if(!res)
        {
            daceSetError(__func__, DACE_PANIC, 1);
            exit(1);
        }
        return res;
    }

    /** Wrapper for C free function (free memory).
        @param[in] ptr A pointer to the currently allocated memory.
     */
    void dacefree(void* ptr)
    {
        if(ptr) free(ptr);
    }
#endif      // DACE_MEMORY_MODEL != DACE_MEMORY_STATIC

/** Return a single integer containing all @e nv exponents in @e p with maximum order @e no.
    @note No checking is performed to see if this operation will overflow! This check is done in daceInit().
    @param[in] p A C array of @e nv exponents.
    @param[in] no The maximum computation order used for encoding. Must be strictly positive!
    @param[in] nv The number of exponents in @e p. May be zero, in which case zero is returned.
    @return The corresponding exponents combined in one integer.
*/
unsigned int daceEncodeExponents(const unsigned int p[], const unsigned int no, const unsigned int nv)
{
    if( nv == 0 ) return 0;
    const unsigned int base = no+1;
    unsigned int res = p[nv-1];
    // going backwards to be compatible with old code where res = \sum_{i=0}^{nv-1} p_i * base^i
    for( const unsigned int *ptr = p+nv-2; ptr >= p; ptr-- )
        res = res*base + *ptr;
    return res;
}

/** Decode a single integer containing @e nv exponents of maximum order @e no into @e p.
    @note No checking is performed to see if this operation will overflow! This check is done in daceInit().
    @param[in] ic The encoded integer to decode.
    @param[in] no The maximum computation order used for encoding. Must be strictly positive!
    @param[in] nv The number of exponents in @e p. May be zero, in which case zero is returned.
    @param[out] p A C array of @e nv exponents.
    @return The sum of the decoded exponents (i.e. the order).
*/
unsigned int daceDecodeExponents(unsigned int ic, const unsigned int no, const unsigned int nv, unsigned int p[])
{
    const unsigned int base = no+1;
    unsigned int order = 0;

    for(unsigned int i = 0; i < nv; i++)
    {
        p[i] = ic%base;
        ic /= base;
        order += p[i];
    }

    if(ic != 0)
    {
        daceSetError(__func__, DACE_ERROR, 26);
        for(unsigned int i = 0; i < nv; i++) p[i] = 0;
        return 0;
    }

    return order;
}

/** Enumerate all monomials with @e nv variables and up to maximum order @e no in @e p.
    When called the first time, @e p should be all zeros (i.e. the constant or zeroth
    order monomial). Then call it repeatedly with the previous state of @e p and it will
    put the next monomial in @e p as well as return the order of that monomial.

    When there are no more monomials left, the function will cycle back to the constant monomial
    of order zero.
    @note The sequence of monomials is deterministic but NOT sorted by order. It depends on both
    @e nv and @e no.
    @note This function handles the cases nv=0 and no=0 correctly (i.e. returning 0).
    @param[out] p A C array of @e nv exponents.
    @param[in] no The maximum order of monomials being enumerated.
    @param[in] nv The number of exponents in @e p.
    @return The order of the next monomial returned in @e p.
*/
unsigned int daceNextMonomial(unsigned int p[], const unsigned int no, const unsigned int nv)
{
    unsigned int o = 0;
    for( unsigned int i = 0; i < nv; i++ ) o += p[i];
    for( unsigned int *ptr = p; ptr < p+nv; ptr++ )
    {
        if( o < no )
        {
            *ptr = *ptr + 1;
            return o+1;
        }
        else
        {
            o -= *ptr;
            *ptr = 0;
        }
    }
    // no more monomials
    return 0;
}

/** Enumerate all monomials with nv variables and up to maximum order no in p.
    When called the first time, @e p should be all zeros (i.e. the constant or zeroth
    order monomial). Then call it repeatedly with the previous state of @e p and
    it will put the next monomial in @e p as well as return the order of that monomial.

    When there are no more monomials left, the funcion will cycle back to
    the constant monomial of order zero.
    @note The sequence of monomials is sorted by order, i.e. first all first order
    monomials are returned, then all second order and so on.
    @note This function handles the cases nv=0 and no=0 correctly (i.e. returning 0).
    @param[out] p A C array of @e nv exponents.
    @param[in] no The maximum order of monomials being enumerated.
    @param[in] nv The number of exponents in @e p.
    @return The order of the next monomial returned in @e p.
*/
unsigned int daceNextOrderedMonomial(unsigned int p[], const unsigned int no, const unsigned int nv)
{
    if( nv == 0 || no == 0 ) return 0;
    unsigned int o = 0;
    for( unsigned int i = 0; i < nv; i++ ) o += p[i];
    const unsigned int oo = daceNextMonomial(p+1, o, nv-1);
    if( oo == 0 )
        o = (o+1)%(no+1);        // jump to next order
    p[0] = o-oo;    // complete the monomial up to order o
    return o;
}

/** Compute number of monomials of order @e no in @e nv variables.
    @param[in] no The maximum order of monomials.
    @param[in] nv The number of variables in monomials.
    @return The total number of monomials.
*/
unsigned int daceCountMonomials(unsigned int no, unsigned int nv)
{
    double dnumda = 1.0;
    const unsigned int mm = umax(nv, no);

    for(unsigned int i = 1; i <= umin(nv, no); i++)
    {
        dnumda = (dnumda*(mm+i))/i;
    }
    return (unsigned int)dnumda;
}

/** Encode the given exponents in @e jj into a DA coding integer as stored in `monomial.ii`.
    @param[in] jj A C array of @e nvmax exponents to encode.
    @return The DA coding integer as stored in `monomial.ii`.
*/
unsigned int daceEncode(const unsigned int jj[])
{
    const unsigned int base = DACECom.nomax+1;
    unsigned int io = 0, ic1 = 0, ic2 = 0;

    for(int i = DACECom.nvmax-1; i >= (int)DACECom.nv1; i--)
    {
        ic2 = ic2*base + jj[i];
        io += jj[i];
    }

    for(int i = DACECom.nv1-1; i >= 0; i--)
    {
        ic1 = ic1*base + jj[i];
        io += jj[i];
    }

    if(io > DACECom.nomax)
    {
        daceSetError(__func__, DACE_ERROR, 22);
        return 0;
    }

    return DACECom.ia1[ic1] + DACECom.ia2[ic2];
}

/** Decode the given DA coding integer as stored in `monomial.ii` into @e jj.
    @param[in] jc A DA coding integer as stored in `monomial.ii` to decode.
    @param[out] jj A C array of @e nvmax exponents.
*/
void daceDecode(const unsigned int jc, unsigned int jj[])
{
    const unsigned int io = daceDecodeExponents(DACECom.ie1[jc], DACECom.nomax, DACECom.nv1, jj) +
                            daceDecodeExponents(DACECom.ie2[jc], DACECom.nomax, DACECom.nv2, jj+DACECom.nv1);

    if(io > DACECom.nomax)
    {
        daceSetError(__func__, DACE_ERROR, 25);
        for(unsigned int i = 0; i < DACECom.nvmax; i++) jj[i] = 0;
        return;
    }
}

/** Pack monomials in @e cc into DA object @e inc.
    @param[in] cc A C array of @e nmmax monomials.
    @param[in] inc A pointer to DA object to pack the monomials into.
*/
void dacePack(double * restrict cc, DACEDA * restrict inc)
{
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    monomial* ic = ipoc;
#ifdef DACE_FILTERING
    if(DACECom.lfi == 0)
#endif
    {
        // provide loop without error checking for big enough target DA
        if(ilmc >= DACECom.nmmax)
        {
            for(unsigned int i = 0; i < DACECom.nmmax; i++)
            {
                if(!(fabs(cc[i]) <= DACECom_t.eps))  //  && (DACECom.ieo[i] <= DACECom_t.nocut) is avoided for performance
                {
                    ic->ii = i;
                    ic->cc = cc[i];
                    ic++;
                }
                cc[i] = 0.0;
            }
        }
        else
        {
            for(unsigned int i = 0; i < DACECom.nmmax; i++)
            {
                if(!(fabs(cc[i]) <= DACECom_t.eps) && DACECom.ieo[i] <= DACECom_t.nocut)     // here we remove also cut orders to save as much space as possible
                {
                    if(ic >= ipoc+ilmc)
                    {
                        daceSetError(__func__, DACE_ERROR, 21);
                        for(unsigned int j = i; j < DACECom.nmmax; j++) cc[j] = 0.0;
                        break;
                    }
                    ic->ii = i;
                    ic->cc = cc[i];
                    ic++;
                }
                cc[i] = 0.0;
            }
        }
    }
#ifdef DACE_FILTERING
    else
    {
        for(unsigned int i = 0; i < DACECom.nmmax; i++)
        {
            if((DACECom.ifi[i] != 0) && !(fabs(cc[i]) <= DACECom_t.eps) && (DACECom.ieo[i] <= DACECom_t.nocut))
            {
                if(ic >= ipoc+ilmc)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    for(unsigned int j = 0; j < DACECom.nmmax; j++) cc[j] = 0.0;
                    break;
                }
                ic->ii = i;
                ic->cc = cc[i];
                ic++;
            }
            cc[i] = 0.0;
        }
    }
#endif

    daceSetLength(inc, ic-ipoc);
}

/** Return a pseudo-random number between 0.0 and 1.0.
    @return A pseudo-random number in @f$[0.0, 1.0]@f$.
*/
double daceRandom()
{
    return (double)rand()/(double)RAND_MAX;
}
