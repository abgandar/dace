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
 *  dacemath.c
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

// MS C library needs this to trigger it to define math constants
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "dace/config.h"
#include "dace/dacebase.h"
#include "dace/daceaux.h"
#include "dacecontrib.h"

// define various math constants in case they have not been defined by math.h
// these are non-standard C, but most C libraries have them
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_PI_2
#define M_PI_2 (1.57079632679489661923)
#endif

/********************************************************************************
 *     Basic DACE arithmetic operations
 *********************************************************************************/

/** Perform addition of two DA objects.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceAdd(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    if(!daceIsSameObject(ina, inc) && !daceIsSameObject(inb, inc))
    {
        daceWeightedSum(ina, 1.0, inb, 1.0, inc);
    }
    else
    {
        DACEDA idaadd;
        daceAllocateDA(&idaadd, 0);
        daceWeightedSum(ina, 1.0, inb, 1.0, &idaadd);
        daceCopy(&idaadd, inc);
        daceFreeDA(&idaadd);
    }
}

/** Perform subtraction of two DA objects.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceSubtract(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    if(!daceIsSameObject(ina, inc) && !daceIsSameObject(inb, inc))
    {
        daceWeightedSum(ina, 1.0, inb, -1.0, inc);
    }
    else
    {
        DACEDA idasub;
        daceAllocateDA(&idasub, 0);
        daceWeightedSum(ina, 1.0, inb, -1.0, &idasub);
        daceCopy(&idasub, inc);
        daceFreeDA(&idasub);
    }
}

/** Perform multiplication of two DA objects.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceMultiply(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
// These should use thread local storage (TLS) for multithread safe implementations
// see https://en.wikipedia.org/wiki/Thread-local_storage
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    static DACE_THREAD_LOCAL double cc[DACE_STATIC_NMMAX] = {0};
    static DACE_THREAD_LOCAL extended_monomial emb[DACE_STATIC_NMMAX];
    static DACE_THREAD_LOCAL extended_monomial *ipbeg[DACE_STATIC_NOMAX+1];
    static DACE_THREAD_LOCAL extended_monomial *ipend[DACE_STATIC_NOMAX+1];
	static DACE_THREAD_LOCAL unsigned int nomax = 0;
	static DACE_THREAD_LOCAL unsigned int nvmax = 0;

    // make sure static memory is correctly allocated
    if(nomax != DACECom.nomax || nvmax != DACECom.nvmax)
    {
		nomax = DACECom.nomax;
		nvmax = DACECom.nvmax;
		ipbeg[0] = &emb[0];
		for(unsigned int i = 1; i <= DACECom.nomax; i++)
			ipbeg[i] = emb + daceCountMonomials(i - 1, DACECom.nvmax);
    }
#else
    static DACE_THREAD_LOCAL double *restrict cc = NULL;
    static DACE_THREAD_LOCAL extended_monomial *emb = NULL;
    static DACE_THREAD_LOCAL extended_monomial **ipbeg = NULL;
    static DACE_THREAD_LOCAL extended_monomial **ipend = NULL;
    static DACE_THREAD_LOCAL unsigned int nomax = 0;
	static DACE_THREAD_LOCAL unsigned int nvmax = 0;

    // make sure static memory is correctly allocated
	if(nomax != DACECom.nomax || nvmax != DACECom.nvmax)
	{
		nomax = DACECom.nomax;
		nvmax = DACECom.nvmax;
		dacefree(cc);
        dacefree(emb);
        dacefree(ipbeg);
        dacefree(ipend);
        cc = (double*) dacecalloc(DACECom.nmmax, sizeof(double));
        emb = (extended_monomial*) dacecalloc(DACECom.nmmax, sizeof(extended_monomial));
        ipbeg = (extended_monomial**) dacecalloc(DACECom.nomax+1, sizeof(extended_monomial*));
        ipend = (extended_monomial**) dacecalloc(DACECom.nomax+1, sizeof(extended_monomial*));
		ipbeg[0] = &emb[0];
		for(unsigned int i = 1; i <= DACECom.nomax; i++)
			ipbeg[i] = emb + daceCountMonomials(i - 1, DACECom.nvmax);
	}
#endif

    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    // sort so that ina is the short DA vector
    if(illa > illb)
    {
        unsigned int t1;
        t1 = illb; illb = illa; illa = t1;
        t1 = ilmb; ilmb = ilma; ilma = t1;
        monomial* t2;
        t2 = ipoa; ipoa = ipob; ipob = t2;
    }

    for(unsigned int i = 0; i <= DACECom_t.nocut; i++) ipend[i] = ipbeg[i];

    // sort vector b by order
    for(const monomial *ib = ipob; ib < ipob+illb; ib++)
    {
        const unsigned int noib = DACECom.ieo[ib->ii];
        if(noib > DACECom_t.nocut) continue;
        ipend[noib]->i1 = DACECom.ie1[ib->ii];
        ipend[noib]->i2 = DACECom.ie2[ib->ii];
        ipend[noib]->cc = ib->cc;
        ipend[noib]++;
    }

    // perform actual multiplication
    for(const monomial *ia = ipoa; ia < ipoa+illa; ia++)
    {
        const unsigned int i1ia = DACECom.ie1[ia->ii];
        const unsigned int i2ia = DACECom.ie2[ia->ii];
        const double ccia = ia->cc;
        // Note: all of these inner loops can safely be run in parallel
        //#pragma omp parallel for
        for(int noib = DACECom_t.nocut-DACECom.ieo[ia->ii]; noib >= 0; noib--)
        {
            for(const extended_monomial *ib = ipbeg[noib]; ib < ipend[noib]; ib++)
            {
                const unsigned int ic = DACECom.ia1[i1ia+ib->i1] + DACECom.ia2[i2ia+ib->i2];
                cc[ic] += ccia*ib->cc;
            }
        }
    }

    dacePack(cc, inc);
}

/** Multiply two DA vectors component-wise, i.e. each monomial of ina with the corresponding monomial of inb
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
    @see daceEvalMonomials
 */
void daceMultiplyMonomials(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
	monomial *ipoa; unsigned int ilma, illa;
	monomial *ipob; unsigned int ilmb, illb;
	monomial *ipoc; unsigned int ilmc, illc;

	daceVariableInformation(ina, &ipoa, &ilma, &illa);
	daceVariableInformation(inb, &ipob, &ilmb, &illb);
	daceVariableInformation(inc, &ipoc, &ilmc, &illc);

	monomial *ib = ipob, *ic = ipoc;
	monomial *const ibmax = ipob + ilmb, *const icmax = ipoc + ilmc;

	for(monomial *i = ipoa; i < ipoa + illa; i++)
	{
		while(ib->ii < i->ii && ib < ibmax)
			ib++;
		if(ib == ibmax) break;
		if(ib->ii == i->ii)
		{
			if(ic >= icmax)
			{
				daceSetError(__func__, DACE_ERROR, 21);
				break;
			}
			ic->cc = i->cc*ib->cc;
			ic->ii = i->ii;
			ic++;
		}
	}
}

/** Perform division of two DA objects.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceDivide(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA idadiv;
    const double cons = daceGetConstant(ina)/daceGetConstant(inb);

    daceAllocateDA(&idadiv, 0);
    daceMultiplicativeInverse(inb, &idadiv);
    daceMultiply(ina, &idadiv, inc);
    // set constant part explicitly to ensure correct FP rounding
    if(daceGetConstant(inc) != 0.0)
        daceSetCoefficient0(inc, 0, cons);
    daceFreeDA(&idadiv);
}

/** Square a DA object.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the DA object to square.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceSquare(const DACEDA *ina, DACEDA *inb)
{
    daceMultiply(ina, ina, inb);
}

/** Add constant to a DA object.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to add.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceAddDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    if(!daceIsSameObject(ina, inb))
        daceCopy(ina, inb);

    daceSetCoefficient0(inb, 0, daceGetConstant(inb)+ckon);
}

/** Subtract DA object from constant.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to subtract from.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceDoubleSubtract(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    daceMultiplyDouble(ina, -1.0, inb);
    daceSetCoefficient0(inb, 0, daceGetConstant(inb)+ckon);
}

/** Subtract constant from a DA object.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to subtract.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceSubtractDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    daceAddDouble(ina, -ckon, inb);
}

/** Multiply constant and DA object.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to multiply by.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceMultiplyDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    monomial *ib = ipob;

    if(illa <= ilmb)
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc*ckon;
            if(!(fabs(c) <= DACECom_t.eps))
            {
                ib->cc = c;
                ib->ii = ia->ii;
                ib++;
            }
        }
    }
    else
    {
        monomial *const ibmax = ipob+ilmb;
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc*ckon;
            if(!(fabs(c) <= DACECom_t.eps))
            {
                if(ib >= ibmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    break;
                }
                ib->cc = c;
                ib->ii = ia->ii;
                ib++;
            }
        }
    }

    daceSetLength(inb, ib-ipob);
}

/** Divide DA object by a constant.
    @note This routine is aliasing safe, i.e. @e inb can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to divide by.
    @param[out] inb A pointer to the DA object to store the result in.
 */
void daceDivideDouble(const DACEDA *ina, const double ckon, DACEDA *inb)
{
    if(ckon == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 41);
        daceCreateConstant(inb, 0.0);
        return;
    }

    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);

    monomial *ib = ipob;

    if(illa <= ilmb)
    {
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc/ckon;
            if(!(fabs(c) <= DACECom_t.eps))
            {
                ib->cc = c;
                ib->ii = ia->ii;
                ib++;
            }
        }
    }
    else
    {
        monomial *const ibmax = ipob+ilmb;
        for(monomial *ia = ipoa; ia < ipoa+illa; ia++)
        {
            if(DACECom.ieo[ia->ii] > DACECom_t.nocut)
                continue;

            const double c = ia->cc/ckon;
            if(!(fabs(c) <= DACECom_t.eps))
            {
                if(ib >= ibmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    break;
                }
                ib->cc = c;
                ib->ii = ia->ii;
                ib++;
            }
        }
    }

    daceSetLength(inb, ib-ipob);
}

/** Divide constant by DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] ckon The constant value to divide.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceDoubleDivide(const DACEDA *ina, const double ckon, DACEDA *inc)
{
    const double cons = ckon/daceGetConstant(ina);
    daceMultiplicativeInverse(ina, inc);
    daceMultiplyDouble(inc, ckon, inc);
    // set constant part explicitly to ensure correct FP rounding
    if(daceGetConstant(inc) != 0.0)
        daceSetCoefficient0(inc, 0, cons);
}

/** Divide a DA vector by a single variable to some power, if possible.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] var The number of the independent variable by which to divide.
    @param[in] p The power of the independent variable.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceDivideByVariable(const DACEDA *ina, const unsigned int var, const unsigned int p, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(var < 1 || var > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    // treat a few special cases
    if(p == 0)
    {
        // dividing by 1
        daceCopy(ina, inc);
        return;
    }
    else if(illa == 0)
    {
        // dividing 0 by anything
        daceCreateConstant(inc, 0.0);
        return;
    }
    else if(p > DACECom.nomax)
    {
        // dividing non-zero DA by too high a power
        daceSetError(__func__, DACE_ERROR, 42);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = var-1;
    if(var > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(var > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            if(ipow < p)
            {
                daceSetError(__func__, DACE_ERROR, 42);
                daceCreateConstant(inc, 0.0);
                return;
            }
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2-p*idiv];
            ic->cc = i->cc;
            ic++;
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            if(ipow < p)
            {
                daceSetError(__func__, DACE_ERROR, 42);
                daceCreateConstant(inc, 0.0);
                return;
            }
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1-p*idiv] + DACECom.ia2[ic2];
            ic->cc = i->cc;
            ic++;
        }
    }

    daceSetLength(inc, ic-ipoc);
}

/** Derivative of DA object with respect to a given independent variable.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] idif The number of the independent variable with respect to which the
    derivative is taken.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceDifferentiate(const unsigned int idif, const DACEDA *ina, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(idif < 1 || idif > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = idif-1;
    if(idif > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(idif > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            if(ipow == 0 || DACECom.ieo[i->ii] > DACECom_t.nocut+1)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2-idiv];
            ic->cc = i->cc*ipow;
            ic++;
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            if(ipow == 0 || DACECom.ieo[i->ii] > DACECom_t.nocut+1)
                continue;
            if(ic >= icmax)
            {
                daceSetError(__func__, DACE_ERROR, 21);
                break;
            }
            ic->ii = DACECom.ia1[ic1-idiv] + DACECom.ia2[ic2];
            ic->cc = i->cc*ipow;
            ic++;
        }
    }

    daceSetLength(inc, ic-ipoc);
}

/** Integral of DA object with respect to a given independent variable.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] iint The number of the independent variable with respect to which the
    integral is taken.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceIntegrate(const unsigned int iint, const DACEDA *ina, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    if(iint < 1 || iint > DACECom.nvmax)
    {
        daceSetError(__func__, DACE_ERROR, 24);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const unsigned int ibase = DACECom.nomax+1;
    unsigned int j = iint-1;
    if(iint > DACECom.nv1)
        j = j-DACECom.nv1;
    const unsigned int idiv = npown(ibase, j);

    monomial *ic = ipoc;
    monomial *const icmax = ipoc+ilmc;

    if(iint > DACECom.nv1)
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            if(DACECom.ieo[i->ii] >= DACECom_t.nocut)
                continue;
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic2/idiv)%ibase;
            const double ccc = i->cc/(ipow+1);
            if(!(fabs(ccc) <= DACECom_t.eps))
            {
                if(ic >= icmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    break;
                }
                ic->ii = DACECom.ia1[ic1] + DACECom.ia2[ic2+idiv];
                ic->cc = ccc;
                ic = ic+1;
            }
        }
    }
    else
    {
        for(monomial *i = ipoa; i < ipoa+illa; i++)
        {
            if(DACECom.ieo[i->ii] >= DACECom_t.nocut)
                continue;
            const unsigned int ic1 = DACECom.ie1[i->ii];
            const unsigned int ic2 = DACECom.ie2[i->ii];
            const unsigned int ipow = (ic1/idiv)%ibase;
            const double ccc = i->cc/(ipow+1);
            if(!(fabs(ccc) <= DACECom_t.eps))
            {
                if(ic >= icmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    break;
                }
                ic->ii = DACECom.ia1[ic1+idiv] + DACECom.ia2[ic2];
                ic->cc = ccc;
                ic = ic+1;
            }
        }
    }

    daceSetLength(inc, ic-ipoc);
}

/********************************************************************************
 *     DACE intrinsic function routines
 *********************************************************************************/

/** Absolute value of a DA object.
    Returns either the DA or the negative of the DA based on the constant part.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
    @see daceNorm
 */
void daceAbsolute(const DACEDA *ina, DACEDA *inc)
{
    daceCopy(ina, inc);
    const double a0 = daceGetConstant(inc);
    if(a0 < 0.0)
        daceMultiplyDouble(inc, -1.0, inc);
}

/** Truncate the constant part of a DA object to an integer.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceTruncate(const DACEDA *ina, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, rint(daceGetConstant(inc)));
}

/** Round the constant part of a DA object to an integer.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceRound(const DACEDA *ina, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, round(daceGetConstant(inc)));
}

/** Modulo the constant part of a DA object by a double.
    The constant part of the result is `fmod(cons(a), p)`.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] p Value with respect to which to compute the modulo.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceModuloDouble(const DACEDA *ina, const double p, DACEDA *inc)
{
    daceCopy(ina, inc);
    daceSetCoefficient0(inc, 0, fmod(daceGetConstant(inc), p));
}

/** Modulo a DA object by another DA object.
    This function calculates `a - trunc(cons(a)/cons(b))*b`. The constant part of the result
    is therefore the same as `fmod(cons(a), cons(b))`.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the second DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceModulo(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA temp;
    const double f = trunc(daceGetConstant(ina)/daceGetConstant(inb));
    daceAllocateDA(&temp, 0);
    daceMultiplyDouble(inb, f, &temp);
    daceSubtract(ina, &temp, inc);
    daceFreeDA(&temp);
}

/** Raise a DA object to the @e p-th power.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] p The power to which to raise the DA object.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void dacePowerDouble(const DACEDA *ina, const double p, DACEDA *inc)
{
    // check simple cases
    if(p == 0.0)
    {
        daceCreateConstant(inc, 1.0);
        return;
    }
    else if(p == (int)p)
    {
        dacePower(ina, (int)p, inc);
        return;
    }

    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 43);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = pow(a0, p);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = xf[i-1]/i*(p-(i-1));

    daceDivideDouble(ina, a0, inc);     // more accurate than including a0 in series (uses non-linear part in EvaluateSeries)
    daceEvaluateSeries(inc, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
   dacefree(xf);
#endif
}

/** Raise a DA object to the @e p-th integer power.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] np The power to which to raise the DA object.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void dacePower(const DACEDA *ina, const int np, DACEDA *inc)
{
    DACEDA itemp;

    // handle some common simple cases directly
    switch(np)
    {
        case 0:
            daceCreateConstant(inc, 1.0);
            return;

        case 1:
            daceCopy(ina, inc);
            return;

        case -1:
            daceMultiplicativeInverse(ina, inc);
            return;
    }

    // handle all other cases, again with common special cases hard coded
    switch(abs(np))
    {
        case 2:
            daceSquare(ina, inc);
            break;

        case 3:
            daceAllocateDA(&itemp, 0);
            daceSquare(ina, &itemp);
            daceMultiply(ina, &itemp, inc);
            daceFreeDA(&itemp);
            break;

        case 4:
            daceAllocateDA(&itemp, 0);
            daceSquare(ina, &itemp);
            daceSquare(&itemp, inc);
            daceFreeDA(&itemp);
            break;

        default:
            daceAllocateDA(&itemp, 0);
            daceCopy(ina, &itemp);
            daceCreateConstant(inc, 1.0);
            unsigned int inp = abs(np);
            while(inp)
            {
                if(inp & 1u)
                    daceMultiply(inc, &itemp, inc);
                inp >>= 1;
                if(inp)
                    daceSquare(&itemp, &itemp);
            }
            daceFreeDA(&itemp);
    }

    if(np < 0)
        daceMultiplicativeInverse(inc, inc);
}

/** Take the @e np-th root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] np The root to take of the DA object.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceRoot(const DACEDA *ina, const int np, DACEDA *inc)
{
    if(np == 0)
    {
        daceSetError(__func__, DACE_ERROR, 44);
        daceCreateConstant(inc, 0.0);
        return;
    }

    const double a0 = daceGetConstant(ina);
    const unsigned int iodd = abs(np) & 1u;

    if((iodd == 0) && (a0 <= 0.0))
    {
        daceSetError(__func__, DACE_ERROR, 45);
        daceCreateConstant(inc, 0.0);
        return;
    }
    else if((iodd == 1) && (a0 == 0.0))
    {
        daceSetError(__func__, DACE_ERROR, 46);
        daceCreateConstant(inc, 0.0);
        return;
    }

    double cr = 1.0/np;
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = copysign(pow(fabs(a0), cr), a0);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-1]/i*cr;
        cr--;
    }

    daceDivideDouble(ina, a0, inc);     // more accurate than including a0 in series (uses non-linear part in EvaluateSeries)
    daceEvaluateSeries(inc, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
   dacefree(xf);
#endif
}

/** Compute the multiplicative inverse of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceMultiplicativeInverse(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);

    if(a0 == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 41);
        daceCreateConstant(inc, 0.0);
        return;
    }

    if(DACECom_t.nocut < 5)
    {
        // lower orders: compute series directly
        daceMultiplicativeInverse0(ina, inc, a0);
    }
    else
    {
        // higher orders: use iteration
        const unsigned int nocut = DACECom_t.nocut;
        DACECom_t.nocut = 2;
        daceMultiplicativeInverse0(ina, inc, a0);
        DACEDA temp;
        daceAllocateDA(&temp, 0);
        for(unsigned int ord = 3; ord <= nocut; ord *= 2)
        {
            DACECom_t.nocut = umin(nocut, 2*ord-1);
            daceMultiply(ina, inc, &temp);
            daceDoubleSubtract(&temp, 2.0, &temp);
            daceMultiply(inc, &temp, inc);
        }
        daceFreeDA(&temp);
    }
}

/** Compute the multiplicative inverse of a DA object using series expansion.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
    @param[in] a0 The constant part of @e ina.
 */
void daceMultiplicativeInverse0(const DACEDA *ina, DACEDA *inc, const double a0)
{
    daceDivideDouble(ina, a0, inc);

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double *xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = 1.0/a0;
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = -xf[i-1];

    daceEvaluateSeries(inc, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the square root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceSquareRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, 2, inc);
}

/** Compute the inverse square root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceInverseSquareRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, -2, inc);
}

/** Compute the cubic root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceCubicRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, 3, inc);
}

/** Compute the inverse cubic root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceInverseCubicRoot(const DACEDA *ina, DACEDA *inc)
{
    daceRoot(ina, -3, inc);
}

/** Compute the hypotenuse of two DA objects.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.or inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the second DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHypotenuse(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA itemp1, itemp2;

    daceAllocateDA(&itemp1, 0);
    daceAllocateDA(&itemp2, 0);
    daceSquare(ina, &itemp1);
    daceSquare(inb, &itemp2);
    daceAdd(&itemp1, &itemp2, inc);
    daceRoot(inc, 2, inc);
    daceFreeDA(&itemp2);
    daceFreeDA(&itemp1);
}

/** Compute the exponential of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceExponential(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = exp(daceGetConstant(ina));
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        xf[i] = xf[i-1]/i;

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the natural logarithm root of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogarithm(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0)
    {
        daceSetError(__func__, DACE_ERROR, 47);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    daceDivideDouble(ina, a0, inc);
    xf[0] = log(a0);
    xf[1] = 1.0;
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-1]/i*(i-1);
    }

    daceEvaluateSeries(inc, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the logarithm with respect to base @e b of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] b The base of the logarithm to use.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogarithmBase(const DACEDA *ina, const double b, DACEDA *inc)
{
    if(b <= 0)
    {
        daceSetError(__func__, DACE_ERROR, 48);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceLogarithm(ina, inc);
    daceDivideDouble(inc, log(b), inc);
}

/** Compute the decadic logarithm of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogarithm10(const DACEDA *ina, DACEDA *inc)
{
    daceLogarithmBase(ina, 10.0, inc);
}

/** Compute the binary logarithm of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogarithm2(const DACEDA *ina, DACEDA *inc)
{
    daceLogarithmBase(ina, 2.0, inc);
}

/** Compute the sine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceSine(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    xf[0] = sin(a0);
    xf[1] = cos(a0);
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the cosine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceCosine(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    xf[0] = cos(a0);
    xf[1] = -sin(a0);
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = -xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the tangent of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(cos(daceGetConstant(ina)) == 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 49);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceAllocateDA(&itemp, 0);
    daceSine(ina, &itemp);
    daceCosine(ina, inc);
    daceDivide(&itemp, inc, inc);
    daceFreeDA(&itemp);
}

/** Compute the arcsine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceArcSine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceAllocateDA(&itemp, 0);
    daceSquare(ina, &itemp);
    daceDoubleSubtract(&itemp, 1.0, &itemp);
    daceSquareRoot(&itemp, &itemp);
    daceDivide(ina, &itemp, inc);
    daceArcTangent(inc, inc);
    daceFreeDA(&itemp);
}

/** Compute the arccosine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceArcCosine(const DACEDA *ina, DACEDA *inc)
{
    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }
    daceArcSine(ina, inc);
    daceDoubleSubtract(inc, M_PI_2, inc);
}

/** Compute the arctangent of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceArcTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA iarg;
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1] = {0};
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif
    const double a0 = daceGetConstant(ina);

    daceAllocateDA(&iarg, 0);
    daceMultiplyDouble(ina, a0, &iarg);
    daceAddDouble(&iarg, 1.0, &iarg);
    daceSubtractDouble(ina, a0, inc);
    daceDivide(inc, &iarg, &iarg);

    double s = 1.0;
    xf[0] = atan(a0);
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i+=2)
    {
        xf[i] = s/i;
        s = -s;
    }

    daceEvaluateSeries(&iarg, xf, inc);
    daceFreeDA(&iarg);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Arctangent of @e ina / @e inb with proper sign in @f$[-pi, pi]@f$.
    This function follows the C standard `atan2(y,x)` function syntax.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] inb A pointer to the second DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceArcTangent2(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    const double cx = daceGetConstant(inb);
    const double cy = daceGetConstant(ina);

    if(cx == 0.0 && cy == 0.0)
    {
        daceCreateConstant(inc, 0.0);
    }
    else
    {
        if(fabs(cy) > fabs(cx))
        {
            daceDivide(inb, ina, inc);
            daceArcTangent(inc, inc);
            if(cy < 0.0)
            {
                daceDoubleSubtract(inc, -M_PI_2, inc);
            }
            else
            {
                daceDoubleSubtract(inc, M_PI_2, inc);
            }
        }
        else
        {
            daceDivide(ina, inb, inc);
            daceArcTangent(inc, inc);
            if(cx < 0.0)
            {
                if(cy > 0.0)
                {
                    daceAddDouble(inc, M_PI, inc);
                }
                else
                {
                    daceAddDouble(inc, -M_PI, inc);
                }
            }
        }
    }
}

/** Compute the hyperbolic sine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicSine(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    xf[0] = sinh(a0);
    xf[1] = cosh(a0);

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the hyperbolic cosine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicCosine(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    xf[0] = cosh(a0);
    xf[1] = sinh(a0);

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = xf[i-2]/(i*(i-1));
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the hyperbolic tangent of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;
    const double a0 = daceGetConstant(ina);

    daceAllocateDA(&itemp, 0);
    if(a0 > 0.0)
    {
        daceMultiplyDouble(ina, -2.0, &itemp);
        daceExponential(&itemp, &itemp);
        daceAddDouble(&itemp, 1.0, inc);
        daceDoubleSubtract(&itemp, 1.0, &itemp);
        daceDivide(&itemp, inc, inc);
    }
    else
    {
        daceMultiplyDouble(ina, 2.0, &itemp);
        daceExponential(&itemp, &itemp);
        daceAddDouble(&itemp, 1.0, inc);
        daceAddDouble(&itemp, -1.0, &itemp);
        daceDivide(&itemp, inc, inc);
    }
    daceFreeDA(&itemp);
}

/** Compute the hyperbolic arcsince of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicArcSine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    daceAllocateDA(&itemp, 0);
    daceSquare(ina, inc);
    daceAddDouble(inc, 1.0, &itemp);
    daceSquareRoot(&itemp, inc);
    daceAdd(ina, inc, &itemp);
    daceLogarithm(&itemp, inc);
    daceFreeDA(&itemp);
}

/** Compute the hyperbolic arccosine of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicArcCosine(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(daceGetConstant(ina) <= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceAllocateDA(&itemp, 0);
    daceSquare(ina, inc);
    daceSubtractDouble(inc, 1.0, &itemp);
    daceSquareRoot(&itemp, inc);
    daceAdd(ina, inc, &itemp);
    daceLogarithm(&itemp, inc);
    daceFreeDA(&itemp);
}

/** Compute the hyperbolic arctangent of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHyperbolicArcTangent(const DACEDA *ina, DACEDA *inc)
{
    DACEDA itemp;

    if(fabs(daceGetConstant(ina)) >= 1.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceAllocateDA(&itemp, 0);
    daceAddDouble(ina, 1.0, &itemp);
    daceDoubleSubtract(ina, 1.0, inc);
    daceDivide(&itemp, inc, inc);
    daceLogarithm(inc, &itemp);
    daceMultiplyDouble(&itemp, 0.5, inc);
    daceFreeDA(&itemp);
}

/** Compute the error function of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceErrorFunction(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    double factor = 2.0*exp(-a0*a0)/sqrt(M_PI);
    xf[0] = erf(a0);
    xf[1] = factor;
    double Hi2 = 1.0;       // Hermite polynomial H_{i-2} = H_0
    double Hi1 = 2.0*a0;    // Hermite polynomial H_{i-1} = H_1

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        factor /= -((double)i);
        xf[i] = factor*Hi1;
        const double temp = 2.0*a0*Hi1 - 2.0*(i-1)*Hi2;     // recursion relation: H_i = 2*x*H_{i-1} - 2*(i-1)*H_{i-2}
        Hi2 = Hi1;
        Hi1 = temp;
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the complementary error function of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceComplementaryErrorFunction(const DACEDA *ina, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    const double a0 = daceGetConstant(ina);
    double factor = -2.0*exp(-a0*a0)/sqrt(M_PI);
    xf[0] = erfc(a0);
    xf[1] = factor;
    double Hi2 = 1.0;       // Hermite polynomial H_{i-2} = H_0
    double Hi1 = 2.0*a0;    // Hermite polynomial H_{i-1} = H_1

    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
		factor /= -((double)i);
		xf[i] = factor*Hi1;
        const double temp = 2.0*a0*Hi1 - 2.0*(i-1)*Hi2;     // recursion relation: H_i = 2*x*H_{i-1} - 2*(i-1)*H_{i-2}
        Hi2 = Hi1;
        Hi1 = temp;
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}


/// @cond

// Wrappers for contributed netlib Bessel functions (not for public use)

/** Compute value of Bessel functions @f$J_n, Y_n@f$ for @f$n \in [n0, n1]@f$.
    @param[in] x The function argument (non-negative).
    @param[in] n0 The lowest order of the Bessel functions to calculate (n0 <= n1).
    @param[in] n1 The highest order of the Bessel functions to calculate (n0 <= n1).
    @param[in] type The type of function to evaluate:\n
              -1: Bessel J function\n
               1: Bessel Y function
    @param[out] bz A C array of size `n1-n0+1` containing the values of @f$B_{n0}, B_{n0+1}, ..., B_{n1}@f$ .
    @return Returns 0 if all values are calculated accurately, -1 if x is too large
           to calculate the result or another error occured, or +1 if some of the
           results are of reduced accuracy.
 */
int BesselWrapper(const double x, const int n0, const int n1, const int type, double *bz)
{
    long int nb = (abs(n0) > abs(n1) ? abs(n0) : abs(n1))+1, ncalc;
	double xx = x, alpha = 0.0;

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    #define DACE_STATIC_MAX_BESSEL_ORDER 100
    if( DACE_STATIC_MAX_BESSEL_ORDER < nb ) return -1;
    double b[DACE_STATIC_MAX_BESSEL_ORDER];
#else
    double* b = (double*) dacecalloc(nb, sizeof(double));
#endif

	if(type < 0)
        rjbesl_(&xx, &alpha, &nb, b, &ncalc);
    else
        rybesl_(&xx, &alpha, &nb, b, &ncalc);

	// discombobulate results
    if(ncalc >= 0)
    {
        ncalc = (ncalc == nb ? 0 : 1);
        double s = (n0%2 == 0 ? 1.0 : -1.0);
        for(int i = n0; i <= n1; i++)
        {
            if(i >= 0)
                *(bz++) = b[i];
            else
            {
                *(bz++) = s*b[-i];    // for integer orders considered here, (-1)^n J_n = J_{-n}, and (-1)^n Y_n = Y_{-n}
                s *= -1.0;
            }
        }
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(b);
#endif
    return ncalc < 0 ? -1 : ncalc;
}

/** Compute value of modified Bessel functions @f$I_n@f$, @f$K_n@f$ for @f$n \in [n0, n1]@f$.
    @param[in] x The function argument (non-negative).
    @param[in] n0 The lowest order of the Bessel functions to calculate (n0 <= n1).
    @param[in] n1 The highest order of the Bessel functions to calculate (n0 <= n1).
    @param[in] type The type of function to evaluate:\n
              -2: Bessel I function, scaled (i.e. `exp(-x)*I_n(x)`)\n
              -1: Bessel I function\n
               1: Bessel K function\n
               2: Bessel K function, scaled (i.e. `exp(x)*K_n(x)`)
    @param[out] bz Array of size `n1-n0+1` containing the values of @f$B_{n0}, B_{n0+1}, ..., B_{n1}@f$.
    @return Returns 0 if all values are calculated accurately, -1 if x is too large
           to calculate the result or another error occured, or +1 if some of the
           results are of reduced accuracy.
 */
int ModifiedBesselWrapper(const double x, const int n0, const int n1, const int type, double *bz)
{
    long int nb = (abs(n0) > abs(n1) ? abs(n0) : abs(n1))+1, ize = abs(type), ncalc;
	double xx = x, alpha = 0.0;

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    #define DACE_STATIC_MAX_BESSEL_ORDER 100
    if( DACE_STATIC_MAX_BESSEL_ORDER < nb ) return -1;
    double b[DACE_STATIC_MAX_BESSEL_ORDER];
#else
    double* b = (double*) dacecalloc(nb, sizeof(double));
#endif

	if(type < 0)
        ribesl_(&xx, &alpha, &nb, &ize, b, &ncalc);
    else
        rkbesl_(&xx, &alpha, &nb, &ize, b, &ncalc);

	// discombobulate results
    if(ncalc >= 0)
    {
        ncalc = (ncalc == nb ? 0 : 1);
        for(int i = n0; i <= n1; i++)
            *(bz++) = b[abs(i)];    // for integer orders considered here, I_n = I_{-n}, and for all orders K_n = K_{-n}
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(b);
#endif
    return ncalc < 0 ? -1 : ncalc;
}

/// @endcond

/** Compute the modified Bessel function @f$I_n@f$ of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part >= 0).
    @param[in] n The order of the Bessel function.
    @param[in] scaled If true, the scaled Bessel function is computed (i.e. `exp(-x)*I_n(x)`).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceBesselIFunction(const DACEDA *ina, const int n, const bool scaled, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    const int res = ModifiedBesselWrapper(a0, n-DACECom_t.nocut, n+DACECom_t.nocut, scaled ? -2 : -1, bz);
    if(res >= 0)
    {
        if(scaled)
            daceEvaluateScaledModifiedBesselFunction(ina, bz, 1.0, inc);
        else
            daceEvaluateBesselFunction(ina, bz, 1.0, 1.0, inc);
    }
    else
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(bz);
#endif
}

/** Compute the modified Bessel function @f$K_n@f$ of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part >= 0).
    @param[in] n The order of the Bessel function.
    @param[in] scaled If true, the scaled Bessel function is computed (i.e. `exp(x)*K_n(x)`).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceBesselKFunction(const DACEDA *ina, const int n, const bool scaled, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    const int res = ModifiedBesselWrapper(a0, n-DACECom_t.nocut, n+DACECom_t.nocut, scaled ? 2 : 1, bz);
    if(res >= 0)
    {
        if(scaled)
            daceEvaluateScaledModifiedBesselFunction(ina, bz, -1.0, inc);
        else
            daceEvaluateBesselFunction(ina, bz, 1.0, -1.0, inc);
    }
    else
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(bz);
#endif
}

/** Compute the Bessel function @f$J_n@f$ of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part >= 0).
    @param[in] n The order of the Bessel function.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceBesselJFunction(const DACEDA *ina, const int n, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    const int res = BesselWrapper(a0, n-DACECom_t.nocut, n+DACECom_t.nocut, -1, bz);
    if(res >= 0)
       daceEvaluateBesselFunction(ina, bz, -1.0, 1.0, inc);
    else
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(bz);
#endif
}

/** Compute the Bessel function @f$Y_n@f$ of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part >= 0).
    @param[in] n The order of the Bessel function.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceBesselYFunction(const DACEDA *ina, const int n, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double bz[2*DACE_STATIC_NOMAX+1];
#else
    double* bz = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    const int res = BesselWrapper(a0, n-DACECom_t.nocut, n+DACECom_t.nocut, 1, bz);
    if(res >= 0)
        daceEvaluateBesselFunction(ina, bz, -1.0, 1.0, inc);
    else
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
    }

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(bz);
#endif
}

/** Evaluate a Bessel function with coefficients @e bz with the non-constant part of @e ina.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] bz A C array of `2*nocut+1` elements containing Bessel functions of orders `n-nocut, ..., n+nocut`.
    @param[in] type Either -1.0 for normal Bessel functions, or +1.0 for modified Bessel functions.
    @param[in] ktype Either -1.0 for modified Bessel K function, or +1.0 for all other Bessel functions.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceEvaluateBesselFunction(const DACEDA *ina, const double bz[], const double type, const double ktype, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
    double binomial[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
    double* binomial = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = bz[DACECom_t.nocut];
    binomial[0] = 1.0;
    double factor = 1.0;
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
    {
        factor *= ktype*0.5/i;
        // calculate binomial coefficients i choose j based on previously calculated i-1 choose j.
        binomial[i] = 1.0;
        for(unsigned int j = i-1; j > 0; j--)
            binomial[j] += binomial[j-1];
        // Calculate n-th derivative of Bessel function C, see http://dlmf.nist.gov/10.6
        // bz contains values of C_{n-o} to C_{n+o} of constant part of ina
        double sign = 1.0, c = 0.0;
        xf[i] = 0.0;
        for(unsigned int j = 0; j <= i; j++)
        {
            // use Kahan summation, since signs oscillate and magnitudes can also vary greatly
            const double y = binomial[j]*sign*bz[DACECom_t.nocut-i+2*j] - c;
            const double t = xf[i] + y;
            c = (t - xf[i]) - y;
            xf[i] = t;
            // in infinite precision the above is equivalent to:
            // xf[i] += binomial[j]*sign*bz[DACECom_t.nocut-i+2*j];
            sign *= type;
        }
        xf[i] *= factor;
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(binomial);
    dacefree(xf);
#endif
}

/** Evaluate a scaled modified Bessel function with coefficients @e bz with the non-constant part of @e ina.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] bz A C array of `2*nocut+1` elements containing modified Bessel functions of orders `n-nocut, ..., n+nocut`.
    @param[in] ktype Either -1.0 for scaled Bessel K function, or +1.0 for scaled Bessel I function.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceEvaluateScaledModifiedBesselFunction(const DACEDA *ina, const double bz[], const double ktype, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
    double binomial[2*DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
    double* binomial = (double*) dacecalloc(2*DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = bz[DACECom_t.nocut];
    binomial[0] = 1.0;
    double factor = 1.0;
    for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
    {
        factor *= ktype*0.5/i;
        // calculate binomial coefficients 2*i-1 choose j based on previously calculated 2*i-2 choose j.
        binomial[2*i-1] = 1.0;
        for(unsigned int j = 2*i-2; j > 0; j--)
            binomial[j] += binomial[j-1];
        // calculate binomial coefficients 2*i choose j based on previously calculated 2*i-1 choose j.
        binomial[2*i] = 1.0;
        for(unsigned int j = 2*i-1; j > 0; j--)
            binomial[j] += binomial[j-1];
        // Calculate n-th derivative of Bessel function C
        // bz contains values of C_{n-o} to C_{n+o} of constant part of ina
        double sign = 1.0, c = 0.0;
        xf[i] = 0.0;
        for(unsigned int j = 0; j <= 2*i; j++)
        {
            // use Kahan summation, since signs oscillate and magnitudes can also vary greatly
            const double y = binomial[j]*sign*bz[DACECom_t.nocut-i+j] - c;
            const double t = xf[i] + y;
            c = (t - xf[i]) - y;
            xf[i] = t;
            // in infinite precision the above is equivalent to:
            // xf[i] += binomial[j]*sign*bz[DACECom_t.nocut-i+j];
            sign *= -1.0;
        }
        xf[i] *= factor;
    }

    daceEvaluateSeries(ina, xf, inc);
#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(binomial);
    dacefree(xf);
#endif
}

/** Compute the partial Logarithmic Gamma function of a DA object (without constant part).
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @note No argument checking is performed to ensure values are within allowable range.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] a0 The constant part of @e ina.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogGammaFunction0(const DACEDA *ina, const double a0, DACEDA *inc)
{
#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    xf[0] = 0.0;
    xf[1] = psi_(&a0);
    double s = 1.0;
    for(unsigned int i = 2; i < DACECom_t.nocut+1; i++)
    {
        xf[i] = (s/i)*zeta_(i, a0, NULL);
        s *= -1.0;
    }

    daceEvaluateSeries(ina, xf, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the Logarithmic Gamma function of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part != 0, -1, -2, ...).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLogGammaFunction(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0 && trunc(a0) == a0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceLogGammaFunction0(ina, a0, inc);
    daceSetCoefficient0(inc, 0, log(dgamma_(&a0)));
}

/** Compute the Gamma function of a DA object.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part != 0, -1, -2, ...).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceGammaFunction(const DACEDA *ina, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0 && trunc(a0) == a0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceLogGammaFunction0(ina, a0, inc);
    daceExponential(inc, inc);
    daceMultiplyDouble(inc, dgamma_(&a0), inc);
}

/** Compute the @e n-th Psi function (the @e n+1 derivative of the logarithmic gamma function) of a DA object.
    For @e n = 0 this is the digamma function.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on (constant part != 0, -1, -2, ...).
    @param[in] n The order of the Psi function ( @e n >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void dacePsiFunction(const DACEDA *ina, const unsigned int n, DACEDA *inc)
{
    const double a0 = daceGetConstant(ina);
    if(a0 <= 0.0 && trunc(a0) == a0)
    {
        daceSetError(__func__, DACE_ERROR, 50);
        daceCreateConstant(inc, 0.0);
        return;
    }

#if DACE_MEMORY_MODEL == DACE_MEMORY_STATIC
    double xf[DACE_STATIC_NOMAX+1];
#else
    double* xf = (double*) dacecalloc(DACECom_t.nocut+1, sizeof(double));
#endif

    if(n == 0)
    {
        xf[0] = psi_(&a0);
        double s = 1.0;
        for(unsigned int i = 1; i < DACECom_t.nocut+1; i++)
        {
            xf[i] = s*zeta_(i+1, a0, NULL);
            s *= -1.0;
        }
    }
    else
    {
        double fac = (n%2 ? 1.0 : -1.0);
        for(unsigned int i = 2; i <= n; i++) fac *= i;
        for(unsigned int i = 0; i < DACECom_t.nocut+1; i++)
        {
            xf[i] = fac*zeta_(n+i+1, a0, NULL);
            fac = -(fac/(i+1))*(n+i+1);
        }
    }

    daceEvaluateSeries(ina, xf, inc);

#if DACE_MEMORY_MODEL != DACE_MEMORY_STATIC
    dacefree(xf);
#endif
}

/** Compute the Legendre polynomial of degree @e n.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] n The degree of the Legendre polynomial (@e n >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceLegendrePolynomial(const DACEDA *ina, const unsigned int n, DACEDA *inc)
{
    if(n == 0)
    {
        daceCreateConstant(inc, 1.0);
    }
    else if(n == 1)
    {
        daceCopy(ina, inc);
    }
    else
    {
        DACEDA P[3], itemp;
        for(unsigned int i = 0; i < 3; i++)
            daceAllocateDA(&P[i], 0);
        daceAllocateDA(&itemp, 0);

        daceCreateConstant(&P[0], 1.0);
        daceCopy(ina, &P[1]);
        for(unsigned int i = 2; i <= n; i++)
        {
            daceMultiply(ina, &P[(i-1)%3], &itemp);
            daceWeightedSum(&itemp, (2*i-1)/(double)i, &P[(i-2)%3], -(double)(i-1)/i, &P[i%3]);
        }
        daceCopy(&P[n%3], inc);

        for(unsigned int i = 0; i < 3; i++)
            daceFreeDA(&P[i]);
        daceFreeDA(&itemp);
    }
}

/* Double factorial of @e n.
 */
inline double ffact(const int n)
{
    double res = 1.0;
    for(int i = n; i >= 2; i -= 2)
        res *= i;
    return res;
}

/** Compute the associated Legendre polynomial of degree @e n and order @e m.
    This function is only a polynomial if @m is even.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] n The degree of the associated Legendre polynomial (@e n >= 0).
    @param[in] m The order of the associated Legendre polynomial (@e n >= @e m >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceAssociatedLegendrePolynomial(const DACEDA *ina, const unsigned int n, const unsigned int m, DACEDA *inc)
{
    if(m > n)
    {
        daceCreateConstant(inc, 0.0);
    }
    else if(n == 0)
    {
        daceCreateConstant(inc, 1.0);
    }
    else
    {
        DACEDA P[3], itemp;
        for(unsigned int i = 0; i < 3; i++)
            daceAllocateDA(&P[i], 0);
        daceAllocateDA(&itemp, 0);

        // calculate P_m,m
        daceSquare(ina, &itemp);
        daceDoubleSubtract(&itemp, 1.0, &itemp);
        dacePowerDouble(&itemp, 0.5*m, &itemp);
        daceMultiplyDouble(&itemp, (m%2 ? -1.0 : 1.0)*ffact(2*(int)m-1), &P[m%3]);
        if(m == n)
        {
            daceCopy(&P[m%3], inc);
            for(unsigned int i = 0; i < 3; i++)
                daceFreeDA(&P[i]);
            daceFreeDA(&itemp);
            return;
        }

        // calculate P_m+1,m
        daceMultiply(ina, &P[m%3], &P[(m+1)%3]);
        daceMultiplyDouble(&P[(m+1)%3], 2*m+1, &P[(m+1)%3]);
        if(m+1 == n)
        {
            daceCopy(&P[(m+1)%3], inc);
            for(unsigned int i = 0; i < 3; i++)
                daceFreeDA(&P[i]);
            daceFreeDA(&itemp);
            return;
        }

        // iterate to P_n,m
        for(unsigned int i = m+2; i <= n; i++)
        {
            daceMultiply(ina, &P[(i-1)%3], &itemp);
            daceWeightedSum(&itemp, (2*i-1)/((double)i-m), &P[(i-2)%3], -(double)(i+m-1)/((double)i-m), &P[i%3]);
        }
        daceCopy(&P[n%3], inc);

        for(unsigned int i = 0; i < 3; i++)
            daceFreeDA(&P[i]);
        daceFreeDA(&itemp);
    }
}

/** Compute the (physicist's) Hermite polynomial of degree @e n.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] n The degree of the Hermite polynomial (@e n >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceHermitePolynomial(const DACEDA *ina, const unsigned int n, DACEDA *inc)
{
    if(n == 0)
    {
        daceCreateConstant(inc, 1.0);
    }
    else if(n == 1)
    {
        daceMultiplyDouble(ina, 2.0, inc);
    }
    else
    {
        DACEDA P[3], itemp;
        for(unsigned int i = 0; i < 3; i++)
            daceAllocateDA(&P[i], 0);
        daceAllocateDA(&itemp, 0);

        daceCreateConstant(&P[0], 1.0);
        daceMultiplyDouble(ina, 2.0, &P[1]);
        for(unsigned int i = 2; i <= n; i++)
        {
            daceMultiply(ina, &P[(i-1)%3], &itemp);
            daceWeightedSum(&itemp, 2.0, &P[(i-2)%3], -2.0*(i-1), &P[i%3]);
        }
        daceCopy(&P[n%3], inc);

        for(unsigned int i = 0; i < 3; i++)
            daceFreeDA(&P[i]);
        daceFreeDA(&itemp);
    }
}

/** Compute the Laguerre polynomial of degree @e n.
    This is the same as the associated Laguerre polynomial of order m = 0.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] n The degree of the Laguerre polynomial (@e n >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
    @see daceAssociatedLaguerrePolynomial
*/
void daceLaguerrePolynomial(const DACEDA *ina, const unsigned int n, DACEDA *inc)
{
    daceAssociatedLaguerrePolynomial(ina, n, 0, inc);
}

/** Compute the associated Laguerre polynomial of degree @e n and order @e m.
    For @e m = 0 this yields the regular Laguerre polynomials.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] n The degree of the associated Laguerre polynomial (@e n >= 0).
    @param[in] m The order of the associated Laguerre polynomial (@e m >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
    @see daceLaguerrePolynomial
 */
void daceAssociatedLaguerrePolynomial(const DACEDA *ina, const unsigned int n, const unsigned int m, DACEDA *inc)
{
    if(n == 0)
    {
        daceCreateConstant(inc, 1.0);
    }
    else if(n == 1)
    {
        daceDoubleSubtract(ina, 1.0+m, inc);
    }
    else
    {
        DACEDA P[3], itemp;
        for(unsigned int i = 0; i < 3; i++)
            daceAllocateDA(&P[i], 0);
        daceAllocateDA(&itemp, 0);

        daceCreateConstant(&P[0], 1.0);
        daceDoubleSubtract(ina, 1.0+m, &P[1]);
        for(unsigned int i = 2; i <= n; i++)
        {
            daceDoubleSubtract(ina, 2*i+m-1, &itemp);
            daceMultiply(&itemp, &P[(i-1)%3], &itemp);
            daceWeightedSum(&itemp, 1.0/i, &P[(i-2)%3], -(double)(i+m-1)/i, &P[i%3]);
        }
        daceCopy(&P[n%3], inc);

        for(unsigned int i = 0; i < 3; i++)
            daceFreeDA(&P[i]);
        daceFreeDA(&itemp);
    }
}

/** Compute the spherical harmonic of degree @e n and order @e m.
    This is the spherical harmonic @f$Y_n^m(\theta, \varphi)@f$ for @f$\varphi = 0@f$
    sometimes also called the spherical associated Legendre functions.
    To obtain the full complex spherical harmonic, multiply by @f$e^{im\varphi}@f$.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object containing the polar angle @f$\theta@f$.
    @param[in] n The degree of the spherical harmonic (@e n >= 0).
    @param[in] m The order of the spherical harmonic (@e n >= @e m >= 0).
    @param[out] inc A pointer to the DA object to store the result in.
    @see daceLegendrePolynomial
 */
void daceSphericalHarmonic(const DACEDA *ina, const unsigned int n, const unsigned int m, DACEDA *inc)
{
    if(m > n)
    {
        daceCreateConstant(inc, 0.0);
        return;
    }

    daceCosine(ina, inc);
    daceAssociatedLegendrePolynomial(inc, n, m, inc);
    double fact = (m%2 ? -1.0 : 1.0)*(2*n+1)/(4*M_PI);
    for(unsigned int i = n+m; i > n-m; i--)
        fact /= i;
    daceMultiplyDouble(inc, sqrt(fact), inc);
}

/** Compute the Beta function (Euler integral of the first kind).
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina or @e inb.
    @param[in] inb A pointer to the first DA object to operate on.
    @param[in] ina A pointer to the second DA object to operate on.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceBetaFunction(const DACEDA *ina, const DACEDA *inb, DACEDA *inc)
{
    DACEDA itemp1, itemp2, itemp3;

    daceAllocateDA(&itemp1, 0);
    daceAllocateDA(&itemp2, 0);
    daceAllocateDA(&itemp3, 0);

    // exp(LogGamma(a) + LogGamma(b) - LogGamma(a+b))
    daceLogGammaFunction(ina, &itemp1);
    daceLogGammaFunction(inb, &itemp2);
    daceAdd(&itemp1, &itemp2, &itemp3);
    daceAdd(ina, inb, &itemp1);
    daceLogGammaFunction(&itemp1, &itemp1);
    daceSubtract(&itemp3, &itemp1, &itemp2);
    daceExponential(&itemp2, inc);

    daceFreeDA(&itemp3);
    daceFreeDA(&itemp2);
    daceFreeDA(&itemp1);
}

/** Evaluate a polynomial with coefficients @e xf with the non-constant part of @e ina.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] xf A C array of `nocut+1` elements containing the coefficients of the polynomial.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceEvaluateSeries(const DACEDA *ina, const double xf[], DACEDA *inc)
{
    daceEvaluateSeries0(ina, xf, DACECom_t.nocut, inc);
}

/** Evaluate a polynomial with coefficients @e xf with the non-constant part of @e ina.
    @note This routine is aliasing safe, i.e. @e inc can be the same as @e ina.
    @param[in] ina A pointer to the DA object to operate on.
    @param[in] xf A C array of @e n+1 elements containing the coefficients of the polynomial.
    @param[in] n The order of the polynomial with coefficients in @e xf.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceEvaluateSeries0(const DACEDA *ina, const double xf[], const unsigned int n, DACEDA *inc)
{
    DACEDA inon;
    const unsigned int nocut = DACECom_t.nocut;
    const unsigned int nmax = umin(nocut, n);

    daceAllocateDA(&inon, 0);
    daceCopy(ina, &inon);
    daceSetCoefficient0(&inon, 0, 0.0);

    DACECom_t.nocut = nocut-nmax+1;
    daceMultiplyDouble(&inon, xf[nmax], inc);
    daceAddDouble(inc, xf[nmax-1], inc);

    // evaluate series
    for(int i = nmax-2; i >= 0; i--)
    {
        DACECom_t.nocut = nocut-i;
        daceMultiply(&inon, inc, inc);
        daceAddDouble(inc, xf[i], inc);
    }

    DACECom_t.nocut = nocut;
    daceFreeDA(&inon);
}

/** Compute the weighted sum of two DA objects.
    @warning This routine is NOT aliasing safe! So @e inc *MUST BE DIFFERENT* from @e ina and @e inb.
    @param[in] ina A pointer to the first DA object to operate on.
    @param[in] afac The weighting factor to multiply @e ina by.
    @param[in] inb A pointer to the second DA object to operate on.
    @param[in] bfac The weighting factor to multiply @e inb by.
    @param[out] inc A pointer to the DA object to store the result in.
 */
void daceWeightedSum(const DACEDA *ina, const double afac, const DACEDA *inb, const double bfac, DACEDA *inc)
{
    monomial *ipoa; unsigned int ilma, illa;
    monomial *ipob; unsigned int ilmb, illb;
    monomial *ipoc; unsigned int ilmc, illc;

    daceVariableInformation(ina, &ipoa, &ilma, &illa);
    daceVariableInformation(inb, &ipob, &ilmb, &illb);
    daceVariableInformation(inc, &ipoc, &ilmc, &illc);

    monomial *ia = ipoa, *ib = ipob, *ic = ipoc;
    monomial *const iamax = ipoa+illa, *const ibmax = ipob+illb, *const icmax = ipoc+ilmc;

    if(illa > 0 && illb > 0)
    {
        // both polynomials have coefficients, merge until one runs out
        unsigned int ja = ia->ii;
        unsigned int jb = ib->ii;
        while(true)
        {
            if(ja == jb)
            {
                // add the two terms
                if(DACECom.ieo[ja] <= DACECom_t.nocut)
                {
                    const double ccc = ia->cc*afac + ib->cc*bfac;
                    if(!(fabs(ccc) <= DACECom_t.eps))
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ia->ii;
                        ic++;
                    }
                }
                ia++; ib++;
                if(ia >= iamax || ib >= ibmax) break;
                ja = ia->ii;
                jb = ib->ii;
            }
            else if(ja < jb)
            {
                // store term a
                if(DACECom.ieo[ja] <= DACECom_t.nocut)
                {
                    const double ccc = ia->cc*afac;
                    if(!(fabs(ccc) <= DACECom_t.eps))
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ia->ii;
                        ic++;
                    }
                }
                ia++;
                if(ia >= iamax) break;
                ja = ia->ii;
            }
            else
            {
                // store term b
                if(DACECom.ieo[jb] <= DACECom_t.nocut)
                {
                    const double ccc = ib->cc*bfac;
                    if(!(fabs(ccc) <= DACECom_t.eps))
                    {
                        if(ic >= icmax)
                        {
                            daceSetError(__func__, DACE_ERROR, 21);
                            daceSetLength(inc, ilmc);
                            return;
                        }
                        ic->cc = ccc;
                        ic->ii = ib->ii;
                        ic++;
                    }
                }
                ib++;
                if(ib >= ibmax) break;
                jb = ib->ii;
            }
        }
    }

    // copy any remaining terms from either ina or inb
    monomial *ismin, *ismax;
    double fac;
    if(ia < iamax)
    {
        ismin = ia;
        ismax = iamax;
        fac = afac;
    }
    else
    {
        ismin = ib;
        ismax = ibmax;
        fac = bfac;
    }

    for(monomial *is = ismin; is < ismax; is++)
    {
        if(DACECom.ieo[is->ii] <= DACECom_t.nocut)
        {
            const double ccc = is->cc*fac;
            if(!(fabs(ccc) <= DACECom_t.eps))
            {
                if(ic >= icmax)
                {
                    daceSetError(__func__, DACE_ERROR, 21);
                    daceSetLength(inc, ilmc);
                    return;
                }
                ic->cc = ccc;
                ic->ii = is->ii;
                ic++;
            }
        }
    }

    daceSetLength(inc, ic-ipoc);
}
