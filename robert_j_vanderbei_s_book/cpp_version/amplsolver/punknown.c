/****************************************************************
Copyright (C) AT&T 1992, 1993, 1994
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

/* Print and skip past unknown keywords, possibly followed by = value. */

#include "stdio.h"

 char *
#ifdef KR_headers
pr_unknown(f, s) FILE *f; char *s;
#else
pr_unknown(FILE *f, char *s)
#endif
{
	char *s1;

	for(s1 = s; *s1 > ' ' && *s1 != '='; s1++);
	fprintf(f, "Unknown keyword \"%.*s\"\n", s1-s, s);
	while(*s1 <= ' ' && *s1)
		s1++;
	if (*s1 == '=') {
		while(*++s1)
			if (*s1 > ' ') {
				while(*++s1 > ' ');
				break;
			}
		}
	return s1;
	}
