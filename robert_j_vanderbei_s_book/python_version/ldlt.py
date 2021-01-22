
from linalg import *

# Converted from C++ to Python by Peiyuan Xu (peiyuanxu71@gmail.com) Southlake TX November 2018
# /*********************************************************************/
# /***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
# /***    All Rights Reserved                                        ***/
# /*********************************************************************/

# include <stdlib.h>
# include <math.h>
# include <string.h>
# //#include <unistd.h>
# include <sys/types.h>
# include <time.h>
# include <stdlib.h>

# ifdef QuadPrec
# include "Quad.h"
# define double Quad
# else
# define high(x) (x)
# endif

# include "lp.h"
# include "myalloc.h"
# include "ldlt.h"
# include "linalg.h"
# LP *openlp(void);

# define _EPS 1.0e-8
# define _EPSSOL 1.0e-6  /* Zero tolerance for consistent eqns w/dep rows */
# define _EPSNUM 0.0     /* Initial zero tolerance for dependent rows */
# define _EPSCDN 1.0e-12 /* Zero tolerance for ill-conditioning test */
# define _EPSDIAG 1.0e-14 /* diagonal perturbation */
# define _STABLTY 1.0    /* 0=fast/unstable, 1=slow/stable */
# define _NOREORD 0
# define _MD  1
# define _MLF 2
# define _DENSE -1
# define _UNSET 0
# define _PRIMAL 1
# define _DUAL 2

_EPS = 1.0e-8
_EPSSOL = 1.0e-6  # /* Zero tolerance for consistent eqns w/dep rows */
_EPSNUM = 0.0  # /* Initial zero tolerance for dependent rows */
_EPSCDN = 1.0e-12  # /* Zero tolerance for ill-conditioning test */
_EPSDIAG = 1.0e-14  # /* diagonal perturbation */
_STABLTY = 1.0  # /* 0=fast/unstable, 1=slow/stable */
_NOREORD = 0
_MD = 1
_MLF = 2
_DENSE = -1
_UNSET = 0
_PRIMAL = 1
_DUAL = 2

# /* Prototype static functions */

# static void inv_sym(
# int m,
# int n,
# int *kA,
# int *iA,
# int *kAt,
# int *iAt,
# int *kQ,
# int *iQ,
# int *bndmark,
# int *rngmark,
# int *tier,
# int verbose
# );

# static void swap(
# int     *v,
# int     i,
# int     j
# );

# static void qksort(
# int     *v,
# int     left,
# int     right
# );

# static void hfall(
# int     heapnum,
# int     *key,
# int     *iheap,
# int     *heap,
# int     cur
# );

# static void hrise(
# int     *key,
# int     *iheap,
# int     *heap,
# int     cur
# );

# static int rawsolve(
# int m,
# int n,
# double *z
# );

# static void lltsym(
# int     m,
# int     *degree,
# int     **nbrs,
# int     *bndmark,
# int     verbose
# );

# static void lltnum(
# int     m,
# int 	n,
# double  *AAt,
# int     *iAAt
# );

# /* Declare static variables */

# static  int     *kAAt=NULL, *iAAt=NULL, *mark=NULL, *perm=NULL, *iperm=NULL;
# static  double  *AAt =NULL, *diag=NULL;

# kAAt=None
# iAAt=None
# mark=None
# perm=None
# iperm=None
# # static  double  *AAt =NULL, *diag=NULL;
# AAt =None
diag = None;

# static  double  epssol, epsnum, epscdn, stablty, epsdiag;
# static  int     method, dense;
# static  int     ndep;
# static  int     denwin;
# static  int     pdf;


# static  double *y_k=NULL, *x_k=NULL,
# *r=NULL, *s=NULL, *z=NULL, *Qx=NULL;

y_k = None
x_k = None
r = None
s = None
z = None
Qx = None


# static  LP *lp=NULL;

# /* Define functions */

# void ldltfac(
# int     m,      /* rows */
# int     n,      /* columns */
# int     *kA,    /* constraint matrix in three linear arrays */
# int     *iA,
# double  *A,
# double  *dn,    /* diagonal matrix for upper-left  corner */
# double  *dm,    /* diagonal matrix for lower-right corner */
# int     *kAt,   /* A^T in three linear arrays */
# int     *iAt,
# double  *At,
# int     verbose /* verbosity */
# )
def ldltfac(m, n, kA, iA, A, dn, dm, kAt, iAt, At, verbose):
    global lp
    # {
    # int i,j;

    # if (lp == None or lp.m <= 0): #{
    lp = LP();
    lp.Q = MALLOC_DOUBLE(0)  # MALLOC( lp->Q,    0, double );
    lp.iQ = MALLOC_INT(0)  # MALLOC( lp->iQ,   0, int );
    lp.kQ = MALLOC_INT(n + 1)  # CALLOC( lp->kQ, n+1, int );
    lp.bndmark = MALLOC_INT(n)  # MALLOC( lp->bndmark, n, int );
    lp.rngmark = MALLOC_INT(m)  # MALLOC( lp->rngmark, m, int );

    for j in range(0, n):  # for (j=0; j<n; j++) { lp->bndmark[j] = BDD_BELOW; }
        lp.bndmark[j] = BDD_BELOW;
    for i in range(0, m):  # for (i=0; i<m; i++) { lp->rngmark[i] = INFINITE; }
        lp.rngmark[i] = INFINITE;

    lp.m = m;
    lp.n = n;
    lp.kA = kA;
    lp.iA = iA;
    lp.A = A;
    lp.kAt = kAt;
    lp.iAt = iAt;
    lp.At = At;
    lp.verbose = verbose;
    # }
    inv_num(lp, dn, dm);


# }

# void inv_num(
# LP 	*lp,	/* pointer to the linear program's data */
# double  *dn,    /* diagonal matrix for upper-left  corner */
# double  *dm     /* diagonal matrix for lower-right corner */
# )

def inv_num(lp, dn, dm):
    global diag
    global epsdiag
    global kAAt
    global iAAt
    global AAt
    global iperm
    global mark
    global perm
    global ndep
    global epssol
    global epsnum
    global epscdn
    global epsdiag
    global stablty
    global method
    global dense
    global pdf

    # {
    # int     m       =lp->m;
    # int     n       =lp->n;
    # int     *kA     =lp->kA;
    # int     *iA     =lp->iA;
    # double  *A      =lp->A;
    # int     *kAt    =lp->kAt;
    # int     *iAt    =lp->iAt;
    # double  *At     =lp->At;
    # int     *kQ     =lp->kQ;
    # int     *iQ     =lp->iQ;
    # double  *Q      =lp->Q;
    # int     *bndmark=lp->bndmark;
    # int     *rngmark=lp->rngmark;
    # int	*tier	=lp->tier;
    # int     verbose =lp->verbose;
    # int     max     =lp->max;

    m = lp.m;
    n = lp.n;
    kA = lp.kA;
    iA = lp.iA;
    A = lp.A;
    kAt = lp.kAt;
    iAt = lp.iAt;
    At = lp.At;
    kQ = lp.kQ;
    iQ = lp.iQ;
    Q = lp.Q;
    bndmark = lp.bndmark;
    rngmark = lp.rngmark;
    tier = lp.tier;
    verbose = lp.verbose;
    max = lp.max;

    # int i, j, k, row, col, *iwork;

    # /*----------------------------------------------+
    # | Set up                                        |
    # |            /       2      t \                 |
    # |            | -( Q+D  )   A  |                 |
    # |     K  =   |       n        |                 |
    # |            |              2 |                 |
    # |            |     A       D  |                 |
    # |            \              m /                 |
    # |                                               |
    # | in AAt, iAAt, kAAt, and diag. Only the lower  |
    # | triangle of a permutation of K is stored.     |
    # |                                              */

    # /*----------------------------------------------+
    # | If data structures are not set up, read       |
    # | $HOME/.syseq to set parameters and then       |
    # | do a symbolic factorization.  The symmetric   |
    # | reordering is put into perm[] and its         |
    # | inverse is put into iperm[]:                  |
    # |      new_index = iperm[ old_index ]           |
    # |                                              */

    if (diag == None or len(diag) == 0):  # if (diag == NULL) {
        epssol = _EPSSOL;
        epsnum = _EPSNUM;
        epscdn = _EPSCDN;
        epsdiag = _EPSDIAG;
        stablty = _STABLTY;
        method = _MD;
        dense = _DENSE;
        pdf = _UNSET;

        inv_sym(m, n, kA, iA, kAt, iAt, kQ, iQ, \
                bndmark, rngmark, tier, verbose);
    # }

    # /*----------------------------------------------+
    # | Get memory for integer work array.            |
    # |                                              */

    iwork = MALLOC_INT(m + n)  # MALLOC( iwork, m+n, int );

    # /*----------------------------------------------+
    # | Store the diagonal of K in diag[].            |
    # |                                              */
    # print("epsdiag:{:15.15f}".format(epsdiag))  #TODO
    for j in range(0, n):  # for (j=0; j<n; j++) { diag[iperm[j]]   = -MAX(dn[j],epsdiag); }
        diag[iperm[j]] = (-1) * MAX(dn[j], epsdiag);
    for i in range(0, m):  # for (i=0; i<m; i++) { diag[iperm[n+i]] =  MAX(dm[i],epsdiag); }
        diag[iperm[n + i]] = MAX(dm[i], epsdiag);

        # /*----------------------------------------------+
        # | Store lower triangle of permutation of K      |
        # | in AAt[], iAAt[], kAAt[].                     |
        # |                                              */

    for j in range(0, n):  # for (j=0; j<n; j++) {
        col = iperm[j];  # /* col is a new_index */
        for k in range(kAAt[col], kAAt[col + 1]):  # for (k=kAAt[col]; k<kAAt[col+1]; k++) {
            iwork[iAAt[k]] = k;
            AAt[k] = 0.0;
        # }
        for k in range(kA[j], kA[j + 1]):  # for (k=kA[j]; k<kA[j+1]; k++) {
            row = iperm[n + iA[k]];  # /* row is a new_index */
            if (row > col):
                AAt[iwork[row]] = A[k];
        # }
        for k in range(kQ[j], kQ[j + 1]):  # for (k=kQ[j]; k<kQ[j+1]; k++) {
            row = iperm[iQ[k]];  # /* row is a new_index */
            if (row > col):
                AAt[iwork[row]] = -max * Q[k];
            elif (row == col):
                diag[row] = diag[row] - max * Q[k];
        # }
    # }
    for i in range(0, m):  # for (i=0; i<m; i++) {
        col = iperm[n + i];
        for k in range(kAAt[col], kAAt[col + 1]):  # for (k=kAAt[col]; k<kAAt[col+1]; k++) {
            iwork[iAAt[k]] = k;
            AAt[k] = 0.0;
        # }
        for k in range(kAt[i], kAt[i + 1]):  # for (k=kAt[i]; k<kAt[i+1]; k++) {
            row = iperm[iAt[k]];
            if (row > col):
                AAt[iwork[row]] = At[k];
        # }
    # }

    iwork = None  # FREE(iwork);

    # /*----------------------------------------------+
    # | Going into lltnum, any row/column for which   |
    # | mark[] is set to FALSE will be treated as     |
    # | a non-existing row/column.  On return,        |
    # | any dependent rows are also marked FALSE.     |
    # |                                              */

    for i in range(0, m + n):  # for (i=0; i<m+n; i++) mark[i] = TRUE;
        mark[i] = TRUE;

    lltnum(m, n, AAt, iAAt);

    for i in range(0, m + n):  # for (i=0; i<m+n; i++) {
        if ((perm[i] < n and diag[i] > 0.0 and lp.verbose > 2) or \
                (perm[i] >= n and diag[i] < 0.0 and lp.verbose > 2) \
                ):  # {
            print("nonconvex subproblem: diag[{:4d}] = {:10.2e} ".format(i, diag[i]));
            # }
    # }

    # {
    # double mindiag = HUGE_VAL;
    mindiag = HUGE_VAL;

    for i in range(0, m + n):  # for (i=0; i<m+n; i++) {
        if (ABS(diag[i]) < mindiag):
            mindiag = ABS(diag[i]);
    # }

    if (mindiag < 1.0e-14):  # {
        epsdiag = epsdiag * 10;
        if (verbose > 2):
            print("mindiag = {:10.2e}, epsdiag = {:10.2e}".format(high(mindiag), high(epsdiag)));
    # }
    # }

    if (verbose > 1 and ndep > 0):
        print("     dependent  rows:    {:d}".format(ndep));


# }

# void forwardbackward(
# double  *Dn,    /* diagonal matrix for upper-left  corner */
# double  *Dm,    /* diagonal matrix for lower-right corner */
# double	*dx,
# double	*dy
# )
def forwardbackward(Dn, Dm, dx, dy):
    # {
    global lp
    solve(lp, Dn, Dm, dx, dy);


# }

# /*----------------------------------------------+
# | The routine solve() uses rawsolve() together  |
# | with iterative refinement to produce the best |
# | possible solutions for the given              |
# | factorization.                               */

# int solve(
# LP 	*lp,
# double	*Dn,
# double	*Dm,
# double  *c,
# double  *b
# )
def solve(lp, Dn, Dm, c, b):
    # {
    global y_k
    global x_k
    global r
    global s
    global z
    global Qx
    global iperm

    m = lp.m;
    n = lp.n;
    kA = lp.kA;
    iA = lp.iA;
    A = lp.A;
    kAt = lp.kAt;
    iAt = lp.iAt;
    At = lp.At;
    kQ = lp.kQ;
    iQ = lp.iQ;
    Q = lp.Q;
    max = lp.max;

    # int i, j, pass=0, consistent = TRUE;
    passPY = 0;
    consistent = TRUE;

    # int m2;
    # double maxrs, oldmaxrs, maxbc;

    m2 = m + n;

    if (y_k == None or len(y_k) == 0):
        y_k = MALLOC_DOUBLE(m);
    else:
        REALLOC_DOUBLE(y_k, m);
    # if (y_k == NULL)   {  MALLOC( y_k, m,  double );}
    # else	           { REALLOC( y_k, m,  double );}
    if (x_k == None or len(x_k) == 0):
        x_k = MALLOC_DOUBLE(n);
    else:
        REALLOC_DOUBLE(x_k, n);
    # if (x_k == NULL)   {  MALLOC( x_k, n,  double );}
    # else	           { REALLOC( x_k, n,  double );}
    if (r == None or len(r) == 0):
        r = MALLOC_DOUBLE(m);
    else:
        REALLOC_DOUBLE(r, m);
    # if (r == NULL)     {  MALLOC( r,   m,  double );}
    # else	           { REALLOC( r,   m,  double );}
    if (s == None or len(s) == 0):
        s = MALLOC_DOUBLE(n);
    else:
        REALLOC_DOUBLE(s, n);
    # if (s == NULL)     {  MALLOC( s,   n,  double );}
    # else	           { REALLOC( s,   n,  double );}
    if (z == None or len(z) == 0):
        z = MALLOC_DOUBLE(m2);
    else:
        REALLOC_DOUBLE(z, m2);
        # if (z == NULL)     {  MALLOC( z,   m2, double );}
    # else	           { REALLOC( z,   m2, double );}
    if (Qx == None or len(Qx) == 0):
        Qx = MALLOC_DOUBLE(n);
    else:
        REALLOC_DOUBLE(Qx, n);
        # if (Qx == NULL)    {  MALLOC( Qx,  n,  double );}
    # else	           { REALLOC( Qx,  n,  double );}

    maxbc = MAX(maxv(b, m), maxv(c, n)) + 1;

    maxrs = HUGE_VAL;
    while True:  # do {
        if (passPY == 0):  # {
            for j in range(0, n):  # for (j=0; j<n; j++) z[iperm[j]]   = c[j];
                z[iperm[j]] = c[j];
            for i in range(0, m):  # for (i=0; i<m; i++) z[iperm[n+i]] = b[i];
                z[iperm[n + i]] = b[i];
            # }
        else:  # {
            for j in range(0, n):  # for (j=0; j<n; j++) z[iperm[j]]   = s[j];
                z[iperm[j]] = s[j];
            for i in range(0, m):  # for (i=0; i<m; i++) z[iperm[n+i]] = r[i];
                z[iperm[n + i]] = r[i];
            # }

        consistent = rawsolve(m, n, z);

        if (passPY == 0):  # {
            for j in range(0, n):  # for (j=0; j<n; j++) x_k[j] = z[iperm[j]];
                x_k[j] = z[iperm[j]];
            for i in range(0, m):  # for (i=0; i<m; i++) y_k[i] = z[iperm[n+i]];
                y_k[i] = z[iperm[n + i]];
        # }
        else:  # {
            for j in range(0, n):  # for (j=0; j<n; j++) x_k[j] = x_k[j] + z[iperm[j]];
                x_k[j] = x_k[j] + z[iperm[j]];
            for i in range(0, m):  # for (i=0; i<m; i++) y_k[i] = y_k[i] + z[iperm[n+i]];
                y_k[i] = y_k[i] + z[iperm[n + i]];
            # }

        smx(m, n, A, kA, iA, x_k, r);
        smx(n, m, At, kAt, iAt, y_k, s);
        smx(n, n, Q, kQ, iQ, x_k, Qx);

        for j in range(0, n):  # for (j=0; j<n; j++) {
            s[j] = c[j] - (s[j] - Dn[j] * x_k[j] - max * Qx[j]);
            # }
        for i in range(0, m):  # for (i=0; i<m; i++) {
            r[i] = b[i] - (r[i] + Dm[i] * y_k[i]);
            # }

        oldmaxrs = maxrs;
        maxrs = MAX(maxv(r, m), maxv(s, n));

        # /* --- for tuning purposes --- */
        if (lp.verbose > 2 and passPY > 0):  # {
            print("refinement({:3d}): {:8.2e} {:8.2e} {:8.2e} ".format( \
                passPY, high(maxv(s, n)), high(maxv(r, m)), high(maxrs / maxbc)));
            # }
            # /* */

        passPY = passPY + 1;
        if (maxrs > 1.0e-10 * maxbc and maxrs < oldmaxrs / 2):
            pass
        else:
            break;
    # } while( maxrs > 1.0e-10*maxbc && maxrs < oldmaxrs/2 );

    if (maxrs > oldmaxrs and passPY > 1):  # {
        for j in range(0, n):  # for (j=0; j<n; j++) x_k[j] = x_k[j] - z[iperm[j]];
            x_k[j] = x_k[j] - z[iperm[j]];
        for i in range(0, m):  # for (i=0; i<m; i++) y_k[i] = y_k[i] - z[iperm[n+i]];
            y_k[i] = y_k[i] - z[iperm[n + i]];
    # }

    # /*----------------------------------------------------------
    # | overwrite input with output                             */

    for j in range(0, n):  # for (j=0; j<n; j++) c[j] = x_k[j];
        c[j] = x_k[j];
    for i in range(0, m):  # for (i=0; i<m; i++) b[i] = y_k[i];
        b[i] = y_k[i];

    return (consistent);


# }

# /*----------------------------------------------+
# | The routine rawsolve() does the forward,      |
# | diagonal, and backward substitions to solve   |
# | systems of equations involving the known      |
# | factorization.                               */

# static int rawsolve(
# int m,
# int n,
# double *z
# )

def rawsolve(m, n, z):
    global ndep
    global epssol
    global mark
    global kAAt
    global iAAt
    global AAt
    global diag

    # {
    # register int i, consistent = TRUE;
    # int m2, k, row;
    # register double eps = 0.0;
    # double beta;
    consistent = TRUE;
    eps = 0.0;

    m2 = m + n;

    if (ndep == 1):
        eps = epssol * maxv(z, m);

        # /*------------------------------------------------------+
        # |                                                       |
        # |               -1                                      |
        # |       z  <-  L  z                                     |
        # |                                                      */

    for i in range(0, m2):  # for (i=0; i<m2; i++) {
        if (mark[i] == 1):  # {
            beta = z[i];
            for k in range(kAAt[i], kAAt[i + 1]):  # for (k=kAAt[i]; k<kAAt[i+1]; k++) {
                row = iAAt[k];
                z[row] = z[row] - AAt[k] * beta;
            # }
        # }
        elif (fabs(z[i]) > eps):  # {
            consistent = FALSE;
        # }
        else:  # {
            z[i] = 0.0;
        # }
    # }

    # /*------------------------------------------------------+
    # |                                                       |
    # |               -1                                      |
    # |       z  <-  D  z                                     |
    # |                                                      */

    for i in range(m2 - 1, -1, -1):  # for (i=m2-1; i>=0; i--) {
        if (mark[i] == 1):  # {
            z[i] = z[i] / diag[i];
        # }
        elif (fabs(z[i]) > eps):  # {
            consistent = FALSE;
        # }
        else:  # {
            z[i] = 0.0;
        # }
    # }

    # /*------------------------------------------------------+
    # |                                                       |
    # |                t -1                                   |
    # |       z  <-  (L )  z                                  |
    # |                                                      */

    for i in range(m2 - 1, -1, -1):  # for (i=m2-1; i>=0; i--) {
        if (mark[i] == 1):  # {
            beta = z[i];
            for k in range(kAAt[i], kAAt[i + 1]):  # for (k=kAAt[i]; k<kAAt[i+1]; k++) {
                beta = beta - AAt[k] * z[iAAt[k]];
            # }
            z[i] = beta;
        # }
        elif (fabs(z[i]) > eps):  # {
            consistent = FALSE;
        # }
        else:  # {
            z[i] = 0.0;
        # }
    # }

    return (consistent);


# }

def inv_clo():
    global perm
    global iperm
    global iAAt
    global kAAt
    global mark
    global AAt
    global diag
    global y_k
    global x_k
    global r
    global s
    global z
    global Qx

    perm = None
    iperm = None
    iAAt = None
    kAAt = None
    mark = None
    AAt = None
    diag = None
    y_k = None
    x_k = None
    r = None
    s = None
    z = None
    Qx = None


# {
# FREE( perm );  FREE( iperm ); FREE( iAAt ); FREE( kAAt );
# FREE( mark );  FREE( AAt );   FREE( diag );
# FREE( y_k );  FREE( x_k );
# FREE( r );     FREE( s );     FREE( z );    FREE( Qx );
# }


# /* Define static functions */

# static void lltnum(
# int     m,
# int     n,
# double  *AAt,
# int     *iAAt
# )
def lltnum(m, n, AAt, iAAt):
    global denwin
    global ndep
    global diag
    global perm
    global kAAt
    global epsnum
    global mark
    # {
    # register int    kk, k_end;
    # int    *first, *link, i, j, newj, k, k_bgn, row;
    # int    sgn_diagi;
    # int    m2 = m+n;
    m2 = m + n;
    # register double lij_dj, *temp;
    # double lij;
    # double diagi, maxoffdiag, maxdiag=0.0;
    maxdiag = 0.0;

    # /*------------------------------------------------------+
    # |                                                       |
    # | the input is a symmetric matrix  A  with lower        |
    # |       triangle stored sparsely in                     |
    # |       kAAt[], iAAt[], AAt[] and with the diagonal     |
    # |       stored in dia[].                                |
    # | the output is the lower triangular matrix  L          |
    # |       stored sparsely in  kAAt,iAAt,AAt  and          |
    # |       a diagonal matrix  D  stored in the diag.       |
    # |                  t                                    |
    # |          A  = LDL                                     |
    # |                                                       |
    # +------------------------------------------------------*/

    # /*------------------------------------------------------+
    # | initialize constants                                 */

    temp = MALLOC_DOUBLE(m2);  # CALLOC( temp,  m2, double );
    first = MALLOC_INT(m2);  # MALLOC( first, m2, int );
    link = MALLOC_INT(m2);  # MALLOC( link,  m2, int );
    for i in range(0, m2):  # for (i=0; i<m2; i++) link[i] = -1;
        link[i] = -1;

    for i in range(0, m2):  # for (i=0; i<m2; i++) {
        if (ABS(diag[i]) > maxdiag):
            maxdiag = ABS(diag[i]);
    # }
    ndep = 0;

    # /*------------------------------------------------------+
    # | begin main loop - this code is taken from George and  |
    # | Liu's book, pg. 155, modified to do LDLt instead      |
    # | of LLt factorization.                                */
    for i in range(0, m2):  # for (i=0; i<m2; i++) {
        diagi = diag[i];
        if perm[i] < n:  # sgn_diagi = perm[i] < n ? -1 : 1;
            sgn_diagi = -1;
        else:
            sgn_diagi = 1;
        # for (j=link[i]; j != -1; j=newj) {
        j = link[i];
        while j != -1:
            newj = link[j];
            k = first[j];
            lij = AAt[k];
            lij_dj = lij * diag[j];
            diagi = diagi - lij * lij_dj;
            k_bgn = k + 1;
            k_end = kAAt[j + 1];
            if (k_bgn < k_end):  # {
                first[j] = k_bgn;
                row = iAAt[k_bgn];
                link[j] = link[row];
                link[row] = j;
                if (j < denwin):  # {
                    for kk in range(k_bgn, k_end):  # for (kk=k_bgn; kk<k_end; kk++)
                        temp[iAAt[kk]] = temp[iAAt[kk]] + lij_dj * AAt[kk];
                # }
                else:  # {
                    # double *ptr;
                    # ptr = &temp[row];
                    # for (kk=k_bgn; kk<k_end; kk++) {
                    #    *ptr += lij_dj*AAt[kk];
                    #    ptr++;
                    # }
                    ptr = 0;
                    for kk in range(k_bgn, k_end):  # for (kk=k_bgn; kk < k_end; kk++) {
                        temp[row + ptr] = temp[row + ptr] + lij_dj * AAt[kk];
                        ptr = ptr + 1;
                    # }

                # }
            # }
            j = newj;
        # }
        k_bgn = kAAt[i];
        k_end = kAAt[i + 1];
        for kk in range(k_bgn, k_end):  # for (kk=k_bgn; kk<k_end; kk++) {
            row = iAAt[kk];
            AAt[kk] = AAt[kk] - temp[row];
        # }
        if (fabs(diagi) <= epsnum * maxdiag or mark[i] == FALSE):  # {
            # /*
            #    if (sgn_diagi*diagi <= epsnum*maxdiag || mark[i] == FALSE) {
            # */
            ndep = ndep + 1;
            maxoffdiag = 0.0;
            for kk in range(k_bgn, k_end):  # for (kk=k_bgn; kk<k_end; kk++) {
                maxoffdiag = MAX(maxoffdiag, ABS(AAt[kk]));
            # }
            if (maxoffdiag < 1.0e+6 * _EPS):  # {
                mark[i] = FALSE;
            # }
            else:  # {
                diagi = sgn_diagi * _EPS;
            # }
        # }
        diag[i] = diagi;
        if (k_bgn < k_end):  # {
            first[i] = k_bgn;
            row = iAAt[k_bgn];
            link[i] = link[row];
            link[row] = i;
            for kk in range(k_bgn, k_end):  # for (kk=k_bgn; kk<k_end; kk++) {
                row = iAAt[kk];
                if (mark[i]):  # {
                    AAt[kk] = AAt[kk] / diagi;
                # }
                else:  # {
                    AAt[kk] = 0.0;
                # }
                temp[row] = 0.0;
            # }
        # }
    # }

    link = None  # FREE( link );
    first = None  # FREE( first );
    temp = None  # FREE( temp  );


# }

# static void inv_sym(
# int m,
# int n,
# int *kA,
# int *iA,
# int *kAt,
# int *iAt,
# int *kQ,
# int *iQ,
# int *bndmark,
# int *rngmark,
# int *tier,
# int verbose
# )

def inv_sym(m, n, kA, iA, kAt, iAt, kQ, iQ, bndmark, rngmark, tier, verbose):
    # {
    global pdf
    global mark
    global dense
    # int i, j, k, n1, separable=TRUE;
    separable = TRUE;

    # int     *degree, **nbrs;

    # /*----------------------------------------------+
    # | Set up adjacency structure for                |
    # |                                               |
    # |            /       2      t \                 |
    # |            | -( Q+D  )   A  |                 |
    # |     K  =   |       n        |                 |
    # |            |              2 |                 |
    # |            |     A       D  |                 |
    # |            \              m /                 |
    # |                                               |
    # +----------------------------------------------*/

    degree = MALLOC_INT(n + m);  # MALLOC( degree, n+m, int   );
    nbrs = MALLOC2D(n + m);  # MALLOC( nbrs,   n+m, int * );

    # /*-----------------------------------------------------+
    # | First check to see if the problem is separable.     */

    for j in range(0, n):  # for (j=0; j<n; j++) {
        for k in range(kQ[j], kQ[j + 1]):  # for (k=kQ[j]; k<kQ[j+1]; k++) {
            if (iQ[k] != j):  # {
                separable = FALSE;
                break;
            # }
            # }
    # }

    # /*----------------------------------------------------+
    # | Select ordering priority (primal or dual)          */

    # {
    # double dense, fraction, pfillin, dfillin;

    fraction = 1.0e0;
    for j in range(0, n):  # for (j=0; j<n; j++) {
        denseLocal = float(kA[j + 1] - kA[j]) / float(m + 1);  # (double) (kA[j+1]-kA[j])/(m+1);
        fraction = fraction * (1.0e0 - denseLocal * denseLocal);
    # }
    pfillin = 0.5 * m * m * (1.0e0 - fraction);
    if (verbose > 2):
        print("primal fillin estimate: {:10.0f}".format( \
            high(pfillin)));

    fraction = 1.0e0;
    for i in range(0, m):  # for (i=0; i<m; i++) {
        denseLocal = float(kAt[i + 1] - kAt[i]) / float(n + 1);  # (double) (kAt[i+1]-kAt[i])/(n+1);
        fraction = fraction * (1.0e0 - denseLocal * denseLocal);
    # }
    dfillin = 0.5 * n * n * (1.0e0 - fraction);
    if (verbose > 2):
        print("dual   fillin estimate: {:10.0f}\n".format( \
            high(dfillin)));

    if (pdf == _UNSET):  # {
        if (3 * pfillin <= dfillin and separable == TRUE):  # {
            pdf = _PRIMAL;
            if (verbose > 2):
                print("Ordering priority favors PRIMAL");
        # }
        else:  # {
            pdf = _DUAL;
            if (verbose > 2): print("Ordering priority favors DUAL");
        # }
    # }
    # }

    # /*----------------------------------------------+
    # | Initialize nbrs so that nbrs[col][k] con-     |
    # | tains the row index of the k_th nonzero in    |
    # | column col.                                   |
    # | Initialize degree so that degree[col] con-    |
    # | tains the number of nonzeros in column col.   |
    # |                                              */

    for j in range(0, n):  # for (j=0; j<n; j++) {
        # int ne;
        ne = kA[j + 1] - kA[j] + kQ[j + 1] - kQ[j];
        nbrs[j] = MALLOC_INT(ne);  # MALLOC( nbrs[j],    ne, int );

        ne = 0;
        for k in range(kA[j], kA[j + 1]):  # for (k=kA[j]; k<kA[j+1]; k++) {
            nbrs[j][ne] = n + iA[k];
            ne = ne + 1;
        # }
        for k in range(kQ[j], kQ[j + 1]):  # for (k=kQ[j]; k<kQ[j+1]; k++) {
            if (iQ[k] != j):  # {
                nbrs[j][ne] = iQ[k];
                ne = ne + 1;
            # }
        # }

        degree[j] = ne;

    # }
    for i in range(0, m):  # for (i=0; i<m; i++) {
        # int ne;
        ne = kAt[i + 1] - kAt[i];
        nbrs[n + i] = MALLOC_INT(ne);  # MALLOC( nbrs[n+i],    ne, int );

        degree[n + i] = ne;

        ne = 0;
        for k in range(kAt[i], kAt[i + 1]):  # for (k=kAt[i]; k<kAt[i+1]; k++) {
            nbrs[n + i][ne] = iAt[k];
            ne = ne + 1;
        # }
    # }

    # /*----------------------------------------------+
    # | Initialize tier to contain the ordering       |
    # | priority scheme.                              |
    # |                                              */

    if (tier == None or len(tier) == 0):  # {
        tier = MALLOC_INT(n + m);  # MALLOC( tier, n+m, int );
        n1 = 0;
        if (pdf == _PRIMAL):  # {
            for j in range(0, n):  # for (j=0; j<n; j++) {
                if (bndmark[j] != FREEVAR):  # {
                    tier[j] = 0;  # /* 0 */
                # }
                else:  # {
                    tier[j] = 1;  # /* 2 */
                # }
            # }
            for i in range(0, m):  # for (i=0; i<m; i++) {
                if (rngmark[i] == UNCONST):  # {
                    tier[n + i] = 1;  # /* 4 */
                    n1 = n1 + 1;
                # }
                elif (rngmark[i] == INFINITE):  # {
                    tier[n + i] = 1;  # /* 1 */
                # }
                else:  # {
                    tier[n + i] = 1;  # /* 3 */
                    n1 = n1 + 1;
                # }
            # }
        # }
        else:  # {
            for j in range(0, n):  # for (j=0; j<n; j++) {
                if (bndmark[j] != FREEVAR):  # {
                    tier[j] = 1;  # /* 1 */
                # }
                else:  # {
                    tier[j] = 1;  # /* 3 */
                    n1 = n1 + 1;
                # }
            # }
            for i in range(0, m):  # for (i=0; i<m; i++) {
                if (rngmark[i] == UNCONST):  # {
                    tier[n + i] = 1;  # /* 4 */
                # }
                elif (rngmark[i] == INFINITE):  # {
                    tier[n + i] = 0;  # /* 0 */
                # }
                else:  # {
                    tier[n + i] = 1;  # /* 2 */
                # }
            # }
            # }
    # }

    # /*---------------------------------------------------------+
    # | compute maximum column degree of tier zero columns      */

    if (dense < 0):  # {
        # int *colhisto, tot, max, cnt;
        # float denfac=3;
        denfac = 3;

        colhisto = MALLOC_INT(n + m + 1);  # CALLOC( colhisto, n+m+1, int );

        for i in range(0, n + m):  # for (i=0; i<n+m; i++)
            if (tier[i] == 0):  # {
                colhisto[degree[i]] = colhisto[degree[i]] + 1;
            # }

        tot = 0;
        max = n1;
        for i in range(0, n + m + 1):  # for (i=0; i<=n+m; i++) {
            tot = tot + colhisto[i];
            if (tot >= max):
                break;
        # }
        i = i + 1;
        tot = 0;
        cnt = 0;
        for j in range(0, n + m):  # for j in range(0,n+m+1):     #for (j=0; j<=n+m; j++) {
            if (tier[j] == 0):  # {
                tot = tot + degree[j];
                cnt = cnt + 1;
            # }
        # }
        dense = int(denfac * i);
        # /*
        # dense = (int)(denfac*MAX(i,tot/cnt));
        # printf("i = %d, n = %d, m = %d, n1 = %d \n", i,n,m,n1);
        # printf("tot = %d, cnt = %d\n", tot, cnt);
        # */

        colhisto = None  # FREE( colhisto );
    # }

    if (verbose > 2):
        print("dense:                 {:5d}\n".format(dense));

        # /*----------------------------------------------+
        # | Get memory for mark[].                       */

    mark = MALLOC_INT(m + n);  # MALLOC( mark, m+n, int );

    lltsym(n + m, degree, nbrs, tier, verbose);

    degree = None;
    nbrs = None;
    tier = None;

    #    FREE( degree); FREE( nbrs ); FREE( tier );


# }

# static void lltsym(
# int     m,
# int     *degree,
# int     **nbrs,
# int     *tier,
# int     verbose
# )
def lltsym(m, degree, nbrs, tier, verbose):
    global stablty
    global method
    global diag
    global perm
    global iperm
    global kAAt
    global iAAt
    global AAt
    global denwin
    global dense
    # {
    # register int    kk, tag=0, deg, nbr_deg, nbr, nbr2, nbr3;
    tag = 0

    # int    i, ii, i2=0, k, kkk, row, node, cnt,
    # ni, nd, nz, aatnz, lnz, okey,
    # /* memfree=0, */
    # heapnum, cur,
    # maxcolkey=0, penalty;
    i2 = 0
    maxcolkey = 0

    # register int    *spc, *locfil, *iwork, *iwork2;

    # int    *node_nbrs, *nbr_nbrs, *dst,
    # *hkey, *heap, *iheap;

    # double narth;

    for i in range(0, m):  # for (i=0; i<m; i++)
        if (tier[i] == 0):
            maxcolkey = MAX(degree[i], maxcolkey);

    if (verbose > 2):
        print("max tier zero degree:  {:5d}".format(maxcolkey));

    penalty = stablty * m;
    if (method == _MLF):
        penalty = penalty * m;

    if (verbose > 2):
        print("ordering penalty:      {:5.0f}".format(high(penalty)));

        # /*---------------------------------------------------------+
        # | allocate space for perm and iperm.                      */

    perm = MALLOC_INT(m);  # MALLOC( perm,  m, int );
    iperm = MALLOC_INT(m);  # MALLOC( iperm, m, int );

    # /*---------------------------------------------------------+
    # | allocate space for work arrays.                         */

    dst = MALLOC_INT(m);  # MALLOC( dst,     m, int   );
    spc = MALLOC_INT(m);  # MALLOC( spc,     m, int   );
    locfil = MALLOC_INT(m);  # MALLOC( locfil,  m, int   );
    hkey = MALLOC_INT(m);  # MALLOC( hkey,    m, int   );
    heap = MALLOC_INT(m + 1);  # MALLOC( heap,    m, int   );
    iheap = MALLOC_INT(m);  # MALLOC( iheap,   m, int   );
    iwork = MALLOC_INT(m);  # MALLOC( iwork,   m, int   );
    iwork2 = MALLOC_INT(m);  # MALLOC( iwork2,  m, int   );

    # heap--;         /* so that indexing starts from 1 */  TODO

    # /*---------------------------------------------------------+
    # | calculate the number of nonzeros in A.                  */

    aatnz = 0;
    for i in range(0, m):  # for (i=0; i<m; i++) aatnz += degree[i];
        aatnz = aatnz + degree[i];
    aatnz = aatnz / 2;
    lnz = int(aatnz);

    # /*---------------------------------------------------------+
    # | allocate enough space to store symbolic structure of L   |
    # | (without any fillin).                                   */

    kAAt = MALLOC_INT(m + 1);  # MALLOC( kAAt, m+1, int );
    iAAt = MALLOC_INT(lnz);  # MALLOC( iAAt, lnz, int );

    # /*----------------------------------------------+
    # | To reduce the number of REALLOC's, a separate |
    # | array spc[] is set up which tells how much    |
    # | space has been allocated for each node so far |
    # | (hence, spc[i] will always be >= degree[i]). */

    for i in range(0, m):  # for (i=0; i<m; i++) { spc[i] = degree[i]; }
        spc[i] = degree[i];

        # /*---------------------------------------------------------+
        # | miscellaneous initializations.                          */

    for i in range(0, m):  # for (i=0; i<m; i++) {
        perm[i] = -1;
        iperm[i] = -1;
        iwork[i] = 0;
        iwork2[i] = -1;
        # }

        # /*---------------------------------------------------------+
        # | compute local fill for each node                        */

    if (method == _MLF):  # {
        for node in range(0, m):  # for (node=0; node<m; node++) {
            #        int lf = 0;
            lf = 0;
            deg = degree[node];
            node_nbrs = nbrs[node];
            for k in range(0, deg):  # for (k=0; k<deg; k++) {
                nbr = node_nbrs[k];
                nbr_nbrs = nbrs[nbr];
                nbr_deg = degree[nbr];
                tag = tag + 1;
                for kk in range(0, nbr_deg):  # for (kk=0; kk<nbr_deg; kk++) {
                    iwork[nbr_nbrs[kk]] = tag;
                # }
                for kk in range(k + 1, deg):  # for (kk=k+1; kk<deg; kk++) {
                    if (iwork[node_nbrs[kk]] != tag):
                        lf = lf + 1;
                # }
            # }
            locfil[node] = lf;
        # }

    # }

    # /*---------------------------------------------------------+
    # | Hkey determines the basis for our ordering heuristic.    |
    # | For example, if hkey[node] = degree[node], then the code |
    # | will generate a minimum degree ordering.  Implicit in    |
    # | hkey is the tie-breaking rule - currently the rule       |
    # | is somewhat random.  To make it first occuring minimum,  |
    # | change the formula to:                                   |
    # |       hkey[node] = degree[node]*m + node;                |
    # | warning: with this definition of hkey, there is the      |
    # | possibility of integer overflow on moderately large      |
    # | problems.                                                |
    # |                                                          |
    # | Nodes from the last m are assigned a penalty, as are     |
    # | nodes corresponding to dense columns.                    |
    # |                                                         */

    for node in range(0, m):  # for (node=0; node<m; node++) {
        if (method == _MD):  # {
            hkey[node] = degree[node];
        # }
        elif (method == _MLF):  # {
            hkey[node] = locfil[node];
        # }
        else:  # {
            hkey[node] = node;
        # }
    # }
    for node in range(0, m):  # for (node=0; node<m; node++) {
        if (degree[node] > dense and tier[node] == 0):  # {
            tier[node] = 1;
        # }
        hkey[node] = hkey[node] + int(tier[node] * penalty);
    # }

    # /*---------------------------------------------------------+
    # | set up heap structure for quickly accessing minimum.    */

    heapnum = m;
    for node in range(m - 1, -1, -1):  # for (node=m-1; node>=0; node--) {
        cur = node + 1;
        iheap[node] = cur;
        heap[cur] = node;
        hfall(heapnum, hkey, iheap, heap, cur);
    # }

    # /*---------------------------------------------------------+
    # | the min degree ordering loop                             |
    # |                                                         */

    tag = 0;
    nz = 0;
    i = 0;
    kAAt[0] = 0;

    denwin = m;
    while (i < m):  # {

        # /* compute min hkey and find node achieving the min */

        node = heap[1];
        deg = degree[node];
        node_nbrs = nbrs[node];

        if (deg >= m - 1 - i):
            denwin = i;

            # /* mark nodes for elimination: the min hkey node
            #   and any node 'indistinguishable' from it */

        perm[i] = node;
        iperm[node] = i;

        nd = 0;
        i2 = i + 1;
        for k in range(0, deg):  # for (k=0; k<deg; k++) iperm[node_nbrs[k]] = i;
            iperm[node_nbrs[k]] = i;
        for k in range(0, deg):  # for (k=0; k<deg; k++) {
            nbr = node_nbrs[k];
            if (degree[nbr] == deg and tier[nbr] == tier[node]):  # {
                nbr_nbrs = nbrs[nbr];
                for kk in range(0, deg):  # for (kk=0; kk<deg; kk++)
                    if (iperm[nbr_nbrs[kk]] < i):
                        break;
                if (kk == deg):  # {
                    perm[i2] = nbr;
                    iperm[nbr] = i2;
                    i2 = i2 + 1;
                # }
                else:  # {
                    dst[nd] = nbr;
                    nd = nd + 1;
                # }
            # }
            else:  # {
                dst[nd] = nbr;
                nd = nd + 1;
            # }
        # }

        # /* reallocate space for iAAt as necessary */

        ni = i2 - i;  # /* number of indistinguishable nodes */

        cnt = nz + (deg * (deg + 1) - (deg - ni) * (deg - ni + 1)) / 2;
        if (cnt > lnz):  # {
            lnz = MAX(cnt, 2 * lnz);
            REALLOC_INT(iAAt, lnz);  # REALLOC( iAAt, lnz, int );
        # }

        # /* copy adjacency lists in iAAt, kAAt */

        for ii in range(i, i2):  # for (ii=i; ii<i2; ii++) {
            node = perm[ii];
            node_nbrs = nbrs[node];

            kAAt[ii + 1] = kAAt[ii] + deg;

            for k in range(0, degree[node]):  # for (k=0; k<degree[node]; k++) {
                nbr = node_nbrs[k];
                row = iperm[nbr];
                if (row > ii):  # {
                    iAAt[nz] = nbr;
                    nz = nz + 1;
                # }
                elif (row == i):  # {
                    if (nbr != perm[i]):  # {
                        iAAt[nz] = nbr;
                        nz = nz + 1;
                    # }
                # }
            # }
            deg = deg - 1;
            # }

            # /* decrement degrees for each distinguishable
            #   neighbor of 'node' corresponding to the removal of
            #   'node' and the indistinguishable neighbors */

        node = perm[i];
        for k in range(0, nd):  # for (k=0; k<nd; k++) {
            nbr = dst[k];
            nbr_nbrs = nbrs[nbr];
            degree[nbr] = degree[nbr] - 1;
            nbr_deg = degree[nbr];
            # for (kk=0; nbr_nbrs[kk] != node; kk++) ;
            kk = 0;
            while nbr_nbrs[kk] != node:
                kk = kk + 1;

            while kk < nbr_deg:  # for ( ; kk<nbr_deg; kk++) {
                nbr_nbrs[kk] = nbr_nbrs[kk + 1];
                kk = kk + 1;
            # }
        # }
        if (i2 > i + 1):  # {
            for k in range(0, nd):  # for (k=0; k<nd; k++) {
                nbr = dst[k];
                nbr_nbrs = nbrs[nbr];
                nbr_deg = degree[nbr];
                cnt = 0;
                for kk in range(0, nbr_deg):  # for (kk=0; kk<nbr_deg; kk++) {
                    if (iperm[nbr_nbrs[kk]] > i):  # {
                        cnt = cnt + 1;
                    # }
                    else:  # {
                        nbr_nbrs[kk - cnt] = nbr_nbrs[kk];
                    # }
                # }
                degree[nbr] = degree[nbr] - cnt;
            # }
        # }

        for ii in range(i, i2):  # for (ii=i; ii<i2; ii++) {
            node = perm[ii];

            cur = iheap[node];
            okey = hkey[heap[cur]];
            heap[cur] = heap[heapnum];
            iheap[heap[cur]] = cur;
            heapnum = heapnum - 1;
            if (okey < hkey[heap[cur]]):
                hfall(heapnum, hkey, iheap, heap, cur);
            else:
                hrise(hkey, iheap, heap, cur);
        # }

        # /* add links: between each pair of distinguishable
        # nodes adjacent to min-deg node which don't
        # already have a link */

        if (method != _MLF or locfil[perm[i]] > 0):  # {
            # /*
            # if (nd > 1) cnt = memfree / (nd*(nd-1));
            # */
            for k in range(0, nd):  # for (k=0; k<nd; k++) {
                nbr = dst[k];
                nbr_deg = degree[nbr];
                nbr_nbrs = nbrs[nbr];
                tag = tag + 1;
                for kk in range(0, nbr_deg):  # for (kk=0; kk<nbr_deg; kk++) {
                    iwork[nbr_nbrs[kk]] = tag;
                # }
                for kk in range(k + 1, nd):  # for (kk=k+1; kk<nd; kk++) {
                    nbr2 = dst[kk];
                    if (iwork[nbr2] != tag):  # {
                        # int ne, cnt2;

                        if (method == _MLF):  # {
                            cnt2 = 0;
                            ne = degree[nbr2];
                            for kkk in range(0, ne):  # for (kkk=0; kkk<ne; kkk++) {
                                nbr3 = nbrs[nbr2][kkk];
                                if (iwork[nbr3] == tag):  # {
                                    locfil[nbr3] = locfil[nbr3] - 1;
                                    hkey[nbr3] = locfil[nbr3];
                                    if (tier[nbr3] != 0):
                                        hkey[nbr3] = hkey[nbr3] + tier[nbr3] * penalty;
                                    hrise(hkey, iheap, heap, iheap[nbr3]);
                                    cnt2 = cnt2 + 1;
                                # }
                            # }
                            locfil[nbr] = locfil[nbr] + (degree[nbr] - cnt2);
                            locfil[nbr2] = locfil[nbr2] + (degree[nbr2] - cnt2);
                        # }

                        ne = degree[nbr];
                        if (ne >= spc[nbr]):  # {
                            # /*
                            # spc[nbr] = MAX(spc[nbr]+cnt+1,2*spc[nbr]);
                            # memfree -= cnt;
                            # */
                            spc[nbr] = spc[nbr] * 2;
                            REALLOC_INT(nbrs[nbr], spc[nbr]);  # REALLOC( nbrs[nbr],  spc[nbr], int );
                        # }
                        nbrs[nbr][ne] = nbr2;
                        degree[nbr] = degree[nbr] + 1;
                        if (method == _MLF):
                            iwork[nbr2] = tag;

                        ne = degree[nbr2];
                        if (ne >= spc[nbr2]):  # {
                            # /*
                            # spc[nbr2] = MAX(spc[nbr2]+cnt+1,2*spc[nbr2]);
                            # memfree -= cnt;
                            # */
                            spc[nbr2] = spc[nbr2] * 2;
                            REALLOC_INT(nbrs[nbr2], spc[nbr2]);  # REALLOC( nbrs[nbr2], spc[nbr2], int  );
                        # }
                        nbrs[nbr2][ne] = nbr;
                        degree[nbr2] = degree[nbr2] + 1;
                    # }
                # }
            # }
        # }

        # /* adjust heap */

        for k in range(0, nd):  # for (k=0; k<nd; k++) {
            nbr = dst[k];
            if (method == _MD):  # {
                hkey[nbr] = degree[nbr];
            # }
            elif (method == _MLF):  # {
                locfil[nbr] = locfil[nbr] - ni * (degree[nbr] - nd + 1);
                hkey[nbr] = locfil[nbr];
            # }
            else:  # {
                hkey[nbr] = nbr;
            # }
            if (tier[nbr] != 0):
                hkey[nbr] = hkey[nbr] + int(tier[nbr] * penalty);
            hrise(hkey, iheap, heap, iheap[nbr]);
            hfall(heapnum, hkey, iheap, heap, iheap[nbr]);
            # }

        for ii in range(i, i2):  # for (ii=i; ii<i2; ii++) {
            node = perm[ii];
            nbrs[node] = None  # FREE(nbrs[node]);  /* memfree += spc[node]; */
        # }

        i = i2;
    # }
    if (verbose > 2):
        print("size of dense window = {:d} ".format(m - denwin));

    # heap++;
    dst = None  # FREE(dst);
    # FREE(spc); FREE(locfil); FREE(hkey); FREE(heap); FREE(iheap);
    spc = None
    locfil = None
    hkey = None
    heap = None
    iheap = None
    # FREE(iwork); FREE(iwork2);
    iwork = None
    iwork2 = None

    for k in range(0, kAAt[m]):  # for (k=0; k<kAAt[m]; k++) iAAt[k] = iperm[iAAt[k]];
        iAAt[k] = iperm[iAAt[k]];

    for i in range(0, m):  # for (i=0; i<m; i++) qksort(iAAt, kAAt[i], kAAt[i+1]-1);
        qksort(iAAt, kAAt[i], kAAt[i + 1] - 1);

        # /*---------------------------------------------------------+
        # | calculate and print statistics.                         */

    narth = 0.0e0;
    for i in range(0, m):  # for (i=0; i<m; i++) {
        k = kAAt[i + 1] - kAAt[i];
        narth = narth + k * k;  # (double) k*k;
    # }
    narth = narth + 3 * kAAt[m] + m;

    lnz = kAAt[m];

    REALLOC_INT(iAAt, lnz);  # REALLOC( iAAt, lnz,     int );

    if (verbose > 1):  # {
        print("nonzeros:    L {:10d},  arith_ops {:18.0f}".format( \
            lnz, high(narth)));
    # }

    # /*---------------------------------------------------------+
    # | allocate remaining storage.                             */

    AAt = MALLOC_DOUBLE(lnz);  # MALLOC( AAt,  lnz, double );
    diag = MALLOC_DOUBLE(m);  # MALLOC( diag,   m, double );


# }

# /*  *** qksort: sort v[left]...v[right] into increasing order ***
# reference: The C Programming Language,
# Kernighan and Ritchie
# 2nd. edition, page 87.  */

# static void qksort(
# int     *v,
# int     left,
# int     right
# )
def qksort(v, left, right):
    # {
    # int i, last;

    if (left >= right):
        return;  # /* do nothing if array contains
        # fewer than two elements */

    swap(v, left, int((left + right) / 2));  # /* move partition elem */
    last = left;  # /* to v[left] */
    for i in range(left + 1, right + 1):  # for (i = left+1; i <= right; i++)       /* partition */
        if (v[i] < v[left]):
            # swap(v, ++last, i);
            last = last + 1;
            swap(v, last, i);
    swap(v, left, last);  # /* restore partition elem */
    qksort(v, left, last - 1);
    qksort(v, last + 1, right);


# }

# /*  *** swap: interchange v[i] and v[j] */

# static void swap(
# int     *v,
# int     i,
# int     j
# )
def swap(v, i, j):
    # {
    # int temp;

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;


# }

# static void hfall(
# int     heapnum,
# int     *key,
# int     *iheap,
# int     *heap,
# int     cur
# )
def hfall(heapnum, key, iheap, heap, cur):
    # {
    #        int child;

    child = 2 * cur;
    while (child <= heapnum):  # {
        if (child < heapnum and \
                key[heap[child + 1]] < key[heap[child]]):
            child = child + 1;
        if (key[heap[cur]] > key[heap[child]]):  # {
            swap(heap, cur, child);
            swap(iheap, heap[cur], heap[child]);
            cur = child;
            child = 2 * cur;
        # }
        else:
            break;
    # }


# }

# static void hrise(
# int     *key,
# int     *iheap,
# int     *heap,
# int     cur
# )
def hrise(key, iheap, heap, cur):
    # {
    # int parent;

    parent = int(cur / 2);
    while (parent > 0):  # {
        if (key[heap[parent]] > key[heap[cur]]):  # {
            swap(heap, cur, parent);
            swap(iheap, heap[cur], heap[parent]);
            cur = parent;
            parent = int(cur / 2);
        # }
        else:
            break;
    # }
# }
