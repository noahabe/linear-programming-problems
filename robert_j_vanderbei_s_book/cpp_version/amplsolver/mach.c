/****************************************************************
Copyright (C) AT&T 1992
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of AT&T or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

AT&T DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL AT&T OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

/* machine-dependent initializations */

#include "arith.h"
#include "math.h"

#ifndef Long
#define Long long
#endif

double Infinity, negInfinity;
#ifdef __cplusplus
extern "C" void Mach();
#endif
#ifdef KR_headers
#define Void /*void*/
#else
#define Void void
#endif

 void
Mach(Void)
{
#ifdef IEEE_8087
	Long *L = (Long *)&Infinity;
	L[1] = 0x7ff00000;
	L[0] = 0;
/* would use #elif defined, but MIPS doesn't understand #elif */
#else
#ifdef IEEE_MC68k
	Long *L = (Long *)&Infinity;
	L[0] = 0x7ff00000;
	L[1] = 0;
#else
#ifdef HUGE_VAL
	Infinity = HUGE_VAL;
#else
#ifdef HUGE
#ifdef CRAY
	Infinity = 0.5 * HUGE;
	/* brain-damaged machine dies testing if (HUGE > -HUGE) */
#else
	Infinity = HUGE;
#endif
#else
	Cannot define Infinity!!!
#endif
#endif
#endif
#endif
	negInfinity = -Infinity;
	}
