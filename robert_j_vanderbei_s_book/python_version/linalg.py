from utils import *

#include <stdlib.h>

#ifdef QuadPrec
#include "Quad.h"
#define double Quad
#else
#define high(x) (x)
#endif

#include "linalg.h"
#include "myalloc.h"
#include "macros.h"

#/*---------------------------------------------------------------+
#|  inner product between n-vectors x and y                      */

#double dotprod( double *x, double *y, int n)
def dotprod(x, y, n):
#{
        #int i;
        #double dotprod=0.0e0;
    dotprod=0.0e0;

    for i in range(0, n): #for (i=0; i<n; i++) dotprod += x[i]*y[i];
        dotprod = dotprod + x[i]*y[i]

    return (dotprod);
#}

def sdotprod(c, x_B, basics, m):
    prod = 0.0;

    for i in range(0, m):
        prod = prod + c[basics[i]] * x_B[i];

    return prod;



#/*---------------------------------------------------------------+
#|  y = basis submatrix of (a,ka,ia) times x                     */

#void bmx( int m, double *a, int *ka, int *ia, int *basis,
#          double *x, double *y)
def bmx(  m,  a,  ka,  ia,  basis, x,  y):
#{
#        int i,j,k;

    for i in range(0, m):    #for (i=0; i<m; i++) y[i] = 0.0e0;
        y[i] = 0.0e0;
    for i in range(0, m): #for (i=0; i<m; i++) {
        j = basis[i];
        for k in range(ka[j], ka[j+1]): #for (k=ka[j]; k<ka[j+1]; k++)
            y[ia[k]] = y[ia[k]] + a[k]*x[i];
    #}
#}

#/*---------------------------------------------------------------+
#|  y = basis submatrix of (a,ka,ia) transpose times x           */

#void btmx( int m, double *a, int *ka, int *ia, int *basis,
#          double *x, double *y)
def btmx(  m,  a,  ka,  ia,  basis, x,  y):
#{
#        int i,j,k;
    for i in range(0, m): #for (i=0; i<m; i++) y[i] = 0.0e0;
        y[i] = 0.0e0;
    for i in range(0, m): #for (i=0; i<m; i++) {
        j = basis[i];
        for k in range(ka[j], ka[j+1]): #for (k=ka[j]; k<ka[j+1]; k++)
            y[i] = y[i] + a[k]*x[ia[k]];
    #}
#}

#/*---------------------------------------------------------------+
#|  y = sparse matrix (a,ka,ia) times x                          */

#void smx( int m, int n, double *a, int *ka, int *ia, double *x, double *y)
def smx(  m,  n,  a,  ka,  ia,  x,  y):
#{
    #int i,j,k;

    for i in range(0, m): #for (i=0; i<m; i++) y[i] = 0.0e0;
        y[i] = 0.0e0;
    for j in range(0, n): #for (j=0; j<n; j++)
        for k in range(ka[j], ka[j+1]): #for (k=ka[j]; k<ka[j+1]; k++)
            y[ia[k]] = y[ia[k]] + a[k]*x[j];
#}

#/*---------------------------------------------------------------+
#|  (kat,iat,at) = transpose of (ka,ia,a)                        */

#void atnum( int m, int n, int *ka, int *ia, double *a,
#        int *kat,int *iat, double *at
#        )
def atnum(  m,  n,  ka,  ia,  a, kat, iat,  at):
#{
#        int i,j,k,row,addr;
#        int *iwork;

    iwork = MALLOC_INT(m)#        CALLOC( iwork, m, int );

    for k in range(0, ka[n]):  #for (k=0; k<ka[n]; k++) {
        row = ia[k];
        iwork[row]=iwork[row]+1;
    #}
    kat[0] = 0;
    for i in range(0, m): #for (i=0; i<m; i++) {
        kat[i+1] = kat[i] + iwork[i];
        iwork[i] = 0;
    #}
    for j in range(0, n): #for (j=0; j<n; j++) {
        for k in range(ka[j], ka[j+1]): #for (k=ka[j]; k<ka[j+1]; k++) {
            row = ia[k];
            addr = kat[row] +iwork[row];
            iwork[row]=iwork[row]+1;
            iat[addr] = j;
            at[addr]  = a[k];
        #}
    #}
    iwork = None #FREE( iwork );
#}

# +---------------------------------------------------------------+
# |  compute componentwise maximum of n-vector x                  |

def maxv(x, n):
    maxv = 0.0e0;

    for i in range(0, n):
        maxv = MAX(maxv, ABS(x[i]));
    return maxv;
