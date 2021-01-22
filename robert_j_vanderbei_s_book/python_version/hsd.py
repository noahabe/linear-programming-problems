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

# define EPS 1.0e-12
# define MAX_ITER 200

EPS = 1.0e-12
MAX_ITER = 200


# int solver(int m,int n,int nz,int *iA, int *kA,
#		double *A, double *b, double *c, double f,
#		double *x, double *y, double *w, double *z)
def solver(m, n, nz, iA, kA, A, b, c, f, x, y, w, z):
    # {
    # double  *dx, *dw, *dy, *dz;                          /* step directions */
    # double  *fx, *fy, *gx, *gy;
    # double  phi, psi, dphi, dpsi;
    # double  *rho, *sigma, normr, norms;	 		 /* infeasibilites */
    # double  *D, *E;			                 /* diagonal matrices */
    # double  gamma, delta, mu, theta;                     /* parameters */

    # double  *At;			 /* arrays for A^T */
    # int     *iAt, *kAt;

    # int     i,j,iter,v=1,status=5;
    v = 1;
    status = 5;

    # double  primal_obj, dual_obj;

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
    fx = MALLOC_DOUBLE(n);  # MALLOC(    fx, n,   double );
    fy = MALLOC_DOUBLE(m);  # MALLOC(    fy, m,   double );
    gx = MALLOC_DOUBLE(n);  # MALLOC(    gx, n,   double );
    gy = MALLOC_DOUBLE(m);  # MALLOC(    gy, m,   double );

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
        x[j] = 1.0;
        z[j] = 1.0;
    # }

    for i in range(0, m):  # for (i=0; i<m; i++) {
        w[i] = 1.0;
        y[i] = 1.0;
    # }

    phi = 1.0;
    psi = 1.0;

    atnum(m, n, kA, iA, A, kAt, iAt, At);

    # /****************************************************************
    # * 	Display Banner.
    # ****************************************************************/

    print("m = {:d},n = {:d},nz = {:d}".format(m, n, nz));
    print( \
        "--------------------------------------------------------------------------\n" \
        "         |           Primal          |            Dual           |       |\n" \
        "  Iter   |  Obj Value       Infeas   |  Obj Value       Infeas   |  mu   |\n" \
        "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n" \
        );
    # fflush(stdout);

    # /****************************************************************
    # * 	Iteration.
    # ****************************************************************/

    for iter in range(0, MAX_ITER):  # for (iter=0; iter<MAX_ITER; iter++) {

        # /*************************************************************
        # * STEP 1: Compute mu and centering parameter delta.
        # *************************************************************/

        mu = (dotprod(z, x, n) + dotprod(w, y, m) + phi * psi) / (n + m + 1);

        if (iter % 2 == 0):  # {
            delta = 0.0;
        # }
        else:  # {
            delta = 1.0;
        # }

        # /*************************************************************
        # * STEP 1: Compute primal and dual objective function values.
        # *************************************************************/

        primal_obj = dotprod(c, x, n);
        dual_obj = dotprod(b, y, m);

        # /*************************************************************
        # * STEP 2: Check stopping rule.
        # *************************************************************/

        if (mu < EPS):  # {
            if (phi > psi):  # {
                status = 0;
                break;  # /* OPTIMAL */
            # }
            elif (dual_obj < 0.0):  # {
                status = 2;
                break;  # /* PRIMAL INFEASIBLE */
            # }
            elif (primal_obj > 0.0):  # {
                status = 4;
                break;  # /* DUAL INFEASIBLE */
            # }
            else:
                # {
                print("Trouble in river city ");
                status = 4;
                break;
        # }
        # }

        # /*************************************************************
        # * STEP 3: Compute infeasibilities.
        # *************************************************************/

        smx(m, n, A, kA, iA, x, rho);
        for i in range(0, m):  # for (i=0; i<m; i++) {
            rho[i] = rho[i] - b[i] * phi + w[i];
        # }
        normr = sqrt(dotprod(rho, rho, m)) / phi;
        for i in range(0, m):  # for (i=0; i<m; i++) {
            rho[i] = -(1 - delta) * rho[i] + w[i] - delta * mu / y[i];
        # }

        smx(n, m, At, kAt, iAt, y, sigma);
        for j in range(0, n):  # for (j=0; j<n; j++) {
            sigma[j] = -sigma[j] + c[j] * phi + z[j];
        # }
        norms = sqrt(dotprod(sigma, sigma, n)) / phi;
        for j in range(0, n):  # for (j=0; j<n; j++) {
            sigma[j] = -(1 - delta) * sigma[j] + z[j] - delta * mu / x[j];
        # }

        gamma = -(1 - delta) * (dual_obj - primal_obj + psi) + psi - delta * mu / phi;

        # /*************************************************************
        # * Print statistics.
        # *************************************************************/

        print("{:8d}   {:14.7e}  {:8.1e}    {:14.7e}  {:8.1e}  {:8.1e} ".format( \
            iter, high(primal_obj / phi + f), high(normr), \
            high(dual_obj / phi + f), high(norms), high(mu)));
        # fflush(stdout);

        # /*************************************************************
        # * STEP 4: Compute step directions.
        # *************************************************************/

        # for (j=0; j<n; j++) { D[j] = z[j]/x[j]; }
        for j in range(0, n):
            D[j] = z[j] / x[j];
        # for (i=0; i<m; i++) { E[i] = w[i]/y[i]; }
        for i in range(0, m):
            E[i] = w[i] / y[i];

        ldltfac(n, m, kAt, iAt, At, E, D, kA, iA, A, v);

        # for (j=0; j<n; j++) { fx[j] = -sigma[j]; }
        for j in range(0, n):
            fx[j] = -sigma[j];
        # for (i=0; i<m; i++) { fy[i] =  rho[i]; }
        for i in range(0, m):
            fy[i] = rho[i];

        forwardbackward(E, D, fy, fx);

        # for (j=0; j<n; j++) { gx[j] = -c[j]; }
        for j in range(0, n):
            gx[j] = -c[j];
        # for (i=0; i<m; i++) { gy[i] = -b[i]; }
        for i in range(0, m):
            gy[i] = -b[i];

        forwardbackward(E, D, gy, gx);

        dphi = (dotprod(c, fx, n) - dotprod(b, fy, m) + gamma) / \
               (dotprod(c, gx, n) - dotprod(b, gy, m) - psi / phi);

        # for (j=0; j<n; j++) { dx[j] = fx[j] - gx[j]*dphi; }
        for j in range(0, n):
            dx[j] = fx[j] - gx[j] * dphi;
        # for (i=0; i<m; i++) { dy[i] = fy[i] - gy[i]*dphi; }
        for i in range(0, m):
            dy[i] = fy[i] - gy[i] * dphi;

        # for (j=0; j<n; j++) { dz[j] = delta*mu/x[j] - z[j] - D[j]*dx[j]; }
        for j in range(0, n):
            dz[j] = delta * mu / x[j] - z[j] - D[j] * dx[j];

        # for (i=0; i<m; i++) { dw[i] = delta*mu/y[i] - w[i] - E[i]*dy[i]; }
        for i in range(0, m):
            dw[i] = delta * mu / y[i] - w[i] - E[i] * dy[i];
        dpsi = delta * mu / phi - psi - (psi / phi) * dphi;

        # /*************************************************************
        # * STEP 5: Compute step length.
        # *************************************************************/

        if (iter % 2 == 0):  # {
            pass
        # }
        else:  # {
            theta = 1.0;
        # }
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
        if (theta < -dphi / phi):  # {
            theta = -dphi / phi;
        # }
        if (theta < -dpsi / psi):  # {
            theta = -dpsi / psi;
        # }
        theta = MIN(0.95 / theta, 1.0);

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
        phi = phi + theta * dphi;
        psi = psi + theta * dpsi;
    # }

    for j in range(0, n):  # for (j=0; j<n; j++) {
        x[j] /= phi;
        z[j] /= phi;
    # }
    for i in range(0, m):  # for (i=0; i<m; i++) {
        y[i] /= phi;
        w[i] /= phi;
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
    fx = None  # FREE(    fx );
    fx = None  # FREE(    fy );
    gx = None  # FREE(    gx );
    gy = None  # FREE(    gy );

    At = None  # FREE(   At );
    iAt = None  # FREE(  iAt );
    kAt = None  # FREE(  kAt );

    return status;

# }   /* End of solver */
