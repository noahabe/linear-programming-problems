from ldlt import *

# Converted from C++ to Python by Peiyuan Xu (peiyuanxu71@gmail.com) Southlake TX December 2018
# /*****************************************************************************
#
#                Implementation of the Primal-Dual Interior Point Method
#                              Robert J. Vanderbei
#                                28 November 1994
#
# ******************************************************************************/
# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# ifdef QuadPrec
# include "Quad.h"
# define double Quad
# else
# define high(x) (x)
# endif

# include "linalg.h"
# include "ldlt.h"
# include "myalloc.h"
# include "macros.h"
# /*
# define ABS(x)    ((x) > 0   ? (x) : -(x))
# define SGN(x)    ((x) > 0   ? (1.0) : (-1.0))
# define MAX(x,y)  ((x) > (y) ? (x) : (y))
# define MIN(x,y)  ((x) < (y) ? (x) : (y))
# */

# define EPS 1.0e-6
# define MAX_ITER 200

EPS = 1.0e-6
MAX_ITER = 200


# int solver(int m,int n,int nz,int *iA, int *kA,
#		double *A, double *b, double *c, double f,
#		double *x, double *y, double *w, double *z)
def solver(m, n, nz, iA, kA, A, b, c, f, x, y, w, z):
    global r
    # {
    # double  *dx, *dw, *dy, *dz;                          /* step directions */
    # double  *rho, *sigma, normr0, norms0;	 	 /* infeasibilites */
    # double  *D, *E;			                 /* diagonal matrices */
    # double  gamma, delta, mu, theta, r;                  /* parameters */

    # double  *At;			 /* arrays for A^T */
    # int     *iAt, *kAt;

    # int     i,j,iter,v=1,status=5;
    v = 1;
    status = 5;

    # float   primal_obj, dual_obj, normr, norms;

    # /*******************************************************************
    # * Allocate memory for arrays.
    # *******************************************************************/

    dx = MALLOC_DOUBLE(n);  # MALLOC(    dx, n,   double );
    dw = MALLOC_DOUBLE(m);  # MALLOC(    dw, m,   double );
    dy = MALLOC_DOUBLE(m);  # MALLOC(    dy, m,   double );
    dz = MALLOC_DOUBLE(n);  # MALLOC(    dz, n,   double );
    rho = MALLOC_DOUBLE(m);  # MALLOC(   rho, m,   double );
    sigma = MALLOC_DOUBLE(n);  # MALLOC( sigma, n,   double );
    D = MALLOC_DOUBLE(n);  # MALLOC(     D, n,   double );
    E = MALLOC_DOUBLE(m);  # MALLOC(     E, m,   double );

    At = MALLOC_DOUBLE(nz);  # MALLOC(   At,  nz,  double );
    iAt = MALLOC_INT(nz);  # MALLOC(  iAt,  nz,  int );
    kAt = MALLOC_INT(m + 1);  # MALLOC(  kAt, m+1,  int );

    # /****************************************************************
    # *  Verify input.                				    *
    # ****************************************************************/

    if (m < 20 and n < 20):  # {
        # int k;
        AA = MALLOC2D(20)  # double AA[20][20];
        for i in range(0, 20):
            AA[i] = MALLOC_DOUBLE(20)
        # for (j=0; j<n; j++)
        #	for (i=0; i<m; i++)
        #		AA[i][j] = 0;
        for j in range(0, n):
            for i in range(0, m):
                AA[i][j] = 0;

        # for (j=0; j<n; j++) {
        #	for (k=kA[j]; k<kA[j+1]; k++) {
        #		AA[iA[k]][j] = A[k];
        #	}
        # }
        for j in range(0, n):
            for k in range(kA[j], kA[j + 1]):
                AA[iA[k]][j] = A[k];

        print("A <= b: ");
        # for (i=0; i<m; i++) {
        #	for (j=0; j<n; j++) {
        #		printf(" %5.1f", AA[i][j]);
        #	}
        #	printf("<= %5.1f \n", b[i]);
        # }

        for i in range(0, m):
            for j in range(0, n):
                print(" {:5.1f}".format(AA[i][j]), end="")
            print("<= {:5.1f} ".format(b[i]))

        print("")  # printf("\n");

        print("c: ");
        # for (j=0; j<n; j++) printf(" %5.1f", c[j]);
        for j in range(0, n):
            print(" {:5.1f}".format(c[j]), end="");
        print("");
    # }

    # /****************************************************************
    # *  Initialization.              				    *
    # ****************************************************************/

    for j in range(0, n):  # for (j=0; j<n; j++) {
        x[j] = 1000.0;
        z[j] = 1000.0;
    # }

    for i in range(0, m):  # for (i=0; i<m; i++) {
        w[i] = 1000.0;
        y[i] = 1000.0;
    # }

    atnum(m, n, kA, iA, A, kAt, iAt, At);

    delta = 0.02;
    r = 0.9;

    normr0 = HUGE_VAL;
    norms0 = HUGE_VAL;

    # /****************************************************************
    # * 	Display Banner.
    # ****************************************************************/

    print("m = {:d},n = {:d},nz = {:d}".format(m, n, nz));
    print( \
        "------------------------------------------------------------------\n" \
        "         |           Primal          |            Dual           |\n" \
        "  Iter   |  Obj Value       Infeas   |  Obj Value       Infeas   |\n" \
        "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n" \
        );
    # fflush(stdout); TODO

    # /****************************************************************
    # * 	Iteration.
    # ****************************************************************/

    for iter in range(0, MAX_ITER):  # for (iter=0; iter<MAX_ITER; iter++) {

        # /*************************************************************
        # * STEP 1: Compute infeasibilities.
        # *************************************************************/

        smx(m, n, A, kA, iA, x, rho);

        for i in range(0, m):  # for (i=0; i<m; i++) {
            rho[i] = b[i] - rho[i] - w[i];
        # }

        normr = sqrt(dotprod(rho, rho, m));

        smx(n, m, At, kAt, iAt, y, sigma);
        for j in range(0, n):  # for (j=0; j<n; j++) {
            sigma[j] = c[j] - sigma[j] + z[j];
        # }
        norms = sqrt(dotprod(sigma, sigma, n));

        # /*************************************************************
        # * STEP 2: Compute duality gap.
        # *************************************************************/

        gamma = dotprod(z, x, n) + dotprod(y, w, m);

        # /*************************************************************
        # * Print statistics.
        # *************************************************************/

        primal_obj = dotprod(c, x, n) + f;
        dual_obj = dotprod(b, y, m) + f;
        print("{:8d}   {:14.7e}  {:8.1e}    {:14.7e}  {:8.1e} ".format(iter, primal_obj, normr, dual_obj, norms));
        # fflush(stdout);

        # /*************************************************************
        # * STEP 2.5: Check stopping rule.
        # *************************************************************/

        if (normr < EPS and norms < EPS and gamma < EPS):  # {
            status = 0;
            break;  # /* OPTIMAL */
        # }
        if (normr > 10 * normr0):  # {
            status = 2;
            break;  # /* PRIMAL INFEASIBLE (unreliable) */
        # }
        if (norms > 10 * norms0):  # {
            status = 4;
            break;  # /* DUAL INFEASIBLE (unreliable) */
        # }

        # /*************************************************************
        # * STEP 3: Compute central path parameter.
        # *************************************************************/

        mu = delta * gamma / (n + m);

        # /*************************************************************
        # * STEP 4: Compute step directions.
        # *************************************************************/

        for j in range(0, n):  # for (j=0; j<n; j++) { D[j] = z[j]/x[j]; }
            D[j] = z[j] / x[j];
        for i in range(0, m):  # for (i=0; i<m; i++) { E[i] = w[i]/y[i]; }
            E[i] = w[i] / y[i];

        ldltfac(n, m, kAt, iAt, At, E, D, kA, iA, A, v);

        for j in range(0, n):  # for (j=0; j<n; j++) { dx[j] = sigma[j] - z[j] + mu/x[j]; }
            dx[j] = sigma[j] - z[j] + mu / x[j];
        for i in range(0, m):  # for (i=0; i<m; i++) { dy[i] = rho[i]   + w[i] - mu/y[i]; }
            dy[i] = rho[i] + w[i] - mu / y[i];

        forwardbackward(E, D, dy, dx);

        for j in range(0, n):  # for (j=0; j<n; j++) { dz[j] = mu/x[j] - z[j] - D[j]*dx[j]; }
            dz[j] = mu / x[j] - z[j] - D[j] * dx[j];
        for i in range(0, m):  # for (i=0; i<m; i++) { dw[i] = mu/y[i] - w[i] - E[i]*dy[i]; }
            dw[i] = mu / y[i] - w[i] - E[i] * dy[i];

        # /*************************************************************
        # * STEP 5: Ratio test to find step length.
        # *************************************************************/

        theta = 0.0;
        for j in range(0, n):  # for (j=0; j<n; j++) {
            if (theta < -dx[j] / x[j]):  # {
                theta = -dx[j] / x[j];
            # }
            if (theta < -dz[j] / z[j]):  # {
                theta = -dz[j] / z[j];
            # }
        # }
        for i in range(0, m):  # for (i=0; i<m; i++) {
            if (theta < -dy[i] / y[i]):  # {
                theta = -dy[i] / y[i];
            # }
            if (theta < -dw[i] / w[i]):  # {
                theta = -dw[i] / w[i];
            # }
        # }
        theta = MIN(r / theta, 1.0);

        # /*************************************************************
        # * STEP 6: Step to new point
        # *************************************************************/

        for j in range(0, n):  # for (j=0; j<n; j++) {
            x[j] = x[j] + theta * dx[j];
            z[j] = z[j] + theta * dz[j];
        # }
        for i in range(0, m):  # for (i=0; i<m; i++) {
            y[i] = y[i] + theta * dy[i];
            w[i] = w[i] + theta * dw[i];
        # }

        normr0 = normr;
        norms0 = norms;
    # }

    # /****************************************************************
    # * 	Free work space                                             *
    # ****************************************************************/

    w = None  # FREE(     w );
    z = None  # FREE(     z );
    dx = None  # FREE(    dx );
    dw = None  # FREE(    dw );
    dy = None  # FREE(    dy );
    dz = None  # FREE(    dz );
    rho = None  # FREE(   rho );
    sigma = None  # FREE( sigma );
    D = None  # FREE(     D );
    E = None  # FREE(     E );

    At = None  # FREE(   At );
    iAt = None  # FREE(  iAt );
    kAt = None  # FREE(  kAt );

    return status;

# }   /* End of solver */
