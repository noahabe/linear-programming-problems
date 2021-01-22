import random
import math
import array
import time
import sys
from valind import valind
from tree import *
from lurefac import *

# Converted to Java using C++ to Java Converter - Premium Edition (version 18.10.10)
# Converted from Java to Python by Peiyuan Xu (peiyuanxu71@gmail.com) Southlake TX October 2018
# ******************************************************************************
#
#               Implementation of the
#		Primal Dual (i.e. Self Dual) Simplex Method
#		R. Vanderbei, 3 October 1994
#
# Solves problems in the form:
#
#	     T
#	max c x
#
#	A x  = b
#	  x >= 0
#
# A is an m by N matrix (it is convenient to reserve n for
# the difference N-m).  Artificial variables have been
# added (hence N > m).  One may assume that the last
# m columns of A are invertible and hence can be used as
# a starting basis.

HUGE_VAL = 10000000.0  
EPS = 1.0e-8  
EPS1 = 1.0e-8  
EPS2 = 1.0e-12  
EPS3 = 1.0e-10  
MAX_ITER = 1000000  

a = []  
tag = []  
link = []  
currtag = 1  

a_Nt = []  
tag_Nt = []  
link_Nt = []  

# lurefac functions

EPSSOL = 1.0e-5;  
EPSNUM = 1.0e-9;  
NOREORD = 0;  
MD = 1;  


rank = 0

colperm = []
icolperm = []
rowperm = []
irowperm = []


degL = []
degLt = []

L = []
Lt = []

degU = []
degUt = []

U = []
Ut = []

diag = []

newcol = []

inewcol = []
nnewcol = 0;  
nr = 0;  

perm = []
iperm = []
row_list = []

ngauss = []

gauss = []

rows = []

col_outs = []

imaxs = []

cumtime = 0.0;  
ocumtime = 0.0;  

TRUE = 1;  
FALSE = 0;  

y = []
yy = []

iwork = []
dwork = []
pwork = []

call = 0;  
eps = 0.0

HUGE_VAL = 10000000.0  
EPS = 1.0e-8  
EPS1 = 1.0e-8  
EPS2 = 1.0e-12  
EPS3 = 1.0e-10  
MAX_ITER = 1000000  

a = []  
tag = []  
link = []  
currtag = 1  

a_Nt = []  
tag_Nt = []  
link_Nt = []  

# lurefac functions

EPSSOL = 1.0e-5;  
EPSNUM = 1.0e-9;  
NOREORD = 0;  
MD = 1;  


rank = 0

colperm = []
icolperm = []
rowperm = []
irowperm = []


degL = []
degLt = []

L = []
Lt = []

degU = []
degUt = []

U = []
Ut = []

diag = []

newcol = []

inewcol = []
nnewcol = 0;  
nr = 0;  

perm = []
iperm = []
row_list = []

ngauss = []

gauss = []

rows = []

col_outs = []

imaxs = []

cumtime = 0.0;  
ocumtime = 0.0;  

TRUE = 1;  
FALSE = 0;  

y = []
yy = []

iwork = []
dwork = []
pwork = []

call = 0;  
eps = 0.0

#random.seed(100)

def solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z):

    global v

    x_B = []; 
    y_N = []; 

    xbar_B = []; 
    ybar_N = []; 

    dy_N = []; 
    idy_N = []; 
    ndy_N = 0; 

    dx_B = []; 
    idx_B = []; 

    at = []; 
    iat = [];
    kat = [];

    rscale = [];
    cscale = [];

    basics = []; 
    nonbasics = []; 
    basicflag = []; 

    iter = 0; 
    v = False;

    mu = HUGE_VAL

    vec = [];
    ivec = [];

    status = 0;

    # ********************************************************************
    # * For convenience, we put...
    # ********************************************************************

    N = n + m;

    # ********************************************************************
    # * Add slack variables.  We assume the calling routine allocated
    # * enough memory.
    # ********************************************************************

    i = 0;
    k = ka[n];
    for j in range(n, N):  
        c[j] = 0.0;
        a[k] = 1.0;
        ia[k] = i;
        i = i + 1;
        k = k + 1;
        ka[j + 1] = k;
    nz = k;

    # ********************************************************************
    # * Read in the Data and initialize the common memory sites.
    # ********************************************************************

    x_B = MALLOC_DOUBLE(m)  
    xbar_B = MALLOC_DOUBLE(m)  
    dx_B = MALLOC_DOUBLE(m)  
    y_N = MALLOC_DOUBLE(n)  
    ybar_N = MALLOC_DOUBLE(n)  
    dy_N = MALLOC_DOUBLE(n)  
    vec = MALLOC_DOUBLE(m)  
    nonbasics = MALLOC_INT(n)  
    basics = MALLOC_INT(m)  
    basicflag = MALLOC_INT(N)  
    idx_B = MALLOC_INT(m)  
    idy_N = MALLOC_INT(n)  
    ivec = MALLOC_INT(m)  
    at = MALLOC_DOUBLE(nz)  
    iat = MALLOC_INT(nz)  
    kat = MALLOC_INT(m + 1)  
    rscale = MALLOC_DOUBLE(m)  
    cscale = MALLOC_DOUBLE(n)  

    # ********************************************************************
    # * Add slack variables.  We assume the calling routine allocated
    # * enough memory.
    # ********************************************************************

    i = 0;
    k = ka[n];
    for j in range(n, N):  
        c[j] = 0.0;
        a[k] = 1.0;
        ia[k] = i;
        i = i + 1;
        k = k + 1;
        ka[j + 1] = k;
    nz = k;
    atnum(m, N, ka, ia, a, kat, iat, at);
    for j in range(0, n):  
        for k in range(ka[j], ka[j + 1]):  
            Asq = a[k] * a[k];
            rscale[ia[k]] += Asq;
            cscale[j] += Asq;

    for j in range(0, n):  
        cscale[j] = math.sqrt(cscale[j])

    for i in range(0, m):  
        rscale[i] = math.sqrt(rscale[i])

    for j in range(0, n):
        nonbasics[j] = j;
        basicflag[j] = -j - 1;
        y_N[j] = -c[j];
        ybar_N[j] = drand48() + cscale[j];

    for i in range(0, m):
        basics[i] = n + i;
        basicflag[n + i] = i;
        x_B[i] = b[i];
        xbar_B[i] = drand48() + rscale[i];

    lufac(m, ka, ia, a, basics, 0);

    print('')
    print('m = {:d}, n = {:d}, nz = {:d}'.format(m, N, nz))  

    '''
    print( \
        "---------------------------------------------------------------------------\n" + \
        "          |   Primal      |        |                           arithmetic  \n" + \
        "  Iter    |  Obj Value    |   mu   |   nonz(L)     nonz(U)     operations  \n" + \
        "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" \
        , flush=True);
    '''
    print( \
        "-----------------------------------\n" + \
        "          |   Primal      |        \n" + \
        "  Iter    |  Obj Value    |   mu   \n" + \
        "- - - - - - - - - - - - - - - - - -\n" \
        , flush=True);

    # *****************************************************************
    # * 	Main loop                                             *
    # *****************************************************************

    for iter in range(0, MAX_ITER):  
        # **************************************************************
        # * STEP 1: Find mu                                            *
        # **************************************************************

        mu = -HUGE_VAL;
        col_in = -1;
        for j in range(0, n):  
            if ybar_N[j] > EPS2:
                if mu < -y_N[j] / ybar_N[j]:
                    mu = -y_N[j] / ybar_N[j];
                    col_in = j;

        col_out = -1;
        for i in range(0, m):
            if xbar_B[i] > EPS2:
                if mu < -x_B[i] / xbar_B[i]:
                    mu = -x_B[i] / xbar_B[i];
                    col_out = i;
                    col_in = -1;

        if mu <= EPS3:  
            status = 0;
            break;

        if col_out >= 0:  

            # **************************************************************
            # *                          -1  T                             *
            # * STEP 2: Compute dy  = -(B  N) e                            *
            # *                   N            i			   *
            # *         where i = col_out                                  *
            # **************************************************************

            vec[0] = -1.0;
            ivec[0] = col_out;
            nvec = 1;

            nvec = btsolve(m, vec, ivec, nvec);

            ndy_N = Nt_times_y(N, at, iat, kat, basicflag, vec, ivec, nvec, dy_N, idy_N, ndy_N);

            # **************************************************************
            # * STEP 3: Ratio test to find entering column                 *
            # **************************************************************

            col_in = ratio_test(dy_N, idy_N, ndy_N, y_N, ybar_N, mu);

            if col_in == -1:  
                status = 2;
                break;

                # **************************************************************
                # *                        -1                                  *
                # * STEP 4: Compute dx  = B  N e                               *
                # *                   B         j                              *
                # *                                                            *
                # **************************************************************

            j = nonbasics[col_in];
            i = 0
            k = ka[j]
            while k < ka[j + 1]:  
                dx_B[i] = a[k];
                idx_B[i] = ia[k];
                i = i + 1
                k = k + 1
            ndx_B = i;
            ndx_B = bsolve(m, dx_B, idx_B, ndx_B);
        else:

            # **************************************************************
            # *                        -1                                  *
            # * STEP 2: Compute dx  = B  N e                               *
            # *                   B         j                              *
            # **************************************************************
        
            j = nonbasics[col_in];
            i = 0
            k = ka[j]
            while k < ka[j + 1]:  
                dx_B[i] = a[k];
                idx_B[i] = ia[k];
                i = i + 1
                k = k + 1
            ndx_B = i;
            ndx_B = bsolve(m, dx_B, idx_B, ndx_B);

            # **************************************************************
            # * STEP 3: Ratio test to find leaving column                  *
            # **************************************************************

            col_out = ratio_test(dx_B, idx_B, ndx_B, x_B, xbar_B, mu);

            if col_out == -1:  
                status = 1;
                break;

            # **************************************************************
            # *                          -1  T                             *
            # * STEP 4: Compute dy  = -(B  N) e                            *
            # *                   N            i			   *
            # *                                                            *
            # **************************************************************

            vec[0] = -1.0;
            ivec[0] = col_out;
            nvec = 1;
            nvec = btsolve(m, vec, ivec, nvec);

            ndy_N = Nt_times_y(N, at, iat, kat, basicflag, vec, ivec, nvec, dy_N, idy_N, ndy_N);

        # **************************************************************
        # *                                                            *
        # * STEP 5: Put       t = x /dx                                *
        # *                        i   i                               *
        # *                   _   _                                    *
        # *                   t = x /dx                                *
        # *                        i   i                               *
        # *                   s = y /dy                                *
        # *                        j   j                               *
        # *                   _   _                                    *
        # *                   s = y /dy                                *
        # *                        j   j                               *
        # **************************************************************

        # /* this is inefficient - it should be fixed |
        for k in range(0, ndx_B):  
            if idx_B[k] == col_out:
                break;

        t = x_B[col_out] / dx_B[k];
        tbar = xbar_B[col_out] / dx_B[k];

        # /* this is inefficient - it should be fixed |
        for k in range(0, ndy_N):  
            if idy_N[k] == col_in:
                break;

        s = y_N[col_in] / dy_N[k];
        sbar = ybar_N[col_in] / dy_N[k];

        # **************************************************************
        # *                                _    _    _                 *
        # * STEP 7: Set y  = y  - s dy     y  = y  - s dy              *
        # *              N    N       N     N    N       N             *
        # *                                _    _                      *
        # *             y  = s             y  = s                      *
        # *              i                  i                          *
        # *             _    _    _                                    *
        # *             x  = x  - t dx     x  = x  - t dx              *
        # *              B    B       B     B    B       B             *
        # *             _    _                                         *
        # *             x  = t             x  = t                      *
        # *              j                  j                          *
        # **************************************************************

        for k in range(0, ndy_N):  
            j = idy_N[k];
            y_N[j] -= s * dy_N[k];
            ybar_N[j] -= sbar * dy_N[k];

        y_N[col_in] = s;
        ybar_N[col_in] = sbar;

        for k in range(0, ndx_B):  
            i = idx_B[k];
            x_B[i] -= t * dx_B[k];
            xbar_B[i] -= tbar * dx_B[k];

        x_B[col_out] = t;
        xbar_B[col_out] = tbar;

        # *************************************************************
        # * STEP 8: Update basis                                      *
        # *************************************************************

        i = basics[col_out];
        j = nonbasics[col_in];
        basics[col_out] = j;
        nonbasics[col_in] = i;
        basicflag[i] = -col_in - 1;
        basicflag[j] = col_out;

        # *************************************************************
        # * STEP 9: Refactor basis and print statistics               *
        # *************************************************************

        from_scratch = refactor(m, ka, ia, a, basics, col_out, v);

        if from_scratch == TRUE:
            primal_obj = sdotprod(c, x_B, basics, m) + f;
            print('{:8d}   {:14.7e} {:9.2e} '.format(iter, high(primal_obj), high(mu)));

    primal_obj = sdotprod(c, x_B, basics, m) + f;
    print("");
    print('{:8d}   {:14.7e} {:9.2e} '.format(iter, high(primal_obj), high(mu)));

    # ******************************************************************
    # * 	Transcribe solution to x vector and dual solution to y *
    # ******************************************************************

    for j in range(0, N):
        x[j] = 0.0;
    for i in range(0, m):
        x[basics[i]] = x_B[i];
    for j in range(0, N):
        y[j] = 0.0;
    for i in range(0, n):
        y[nonbasics[i]] = y_N[i];

    # ****************************************************************
    # * 	Split out slack variables and shift dual variables.  *
    # ****************************************************************

    for j in range(0, n):
        z[j] = y[j];
    for i in range(0, m):
        y[i] = y[n + i];
        w[i] = x[n + i];

    # ****************************************************************
    # * 	Free work space                                      *
    # ****************************************************************

    vec = []  
    ivec = []  
    x_B = []  
    y_N = []  
    dx_B = []  
    idx_B = []  
    dy_N = []  
    idy_N = []  
    nonbasics = []  
    basics = []  

    return status;

# the main function

'''
ka = [0, 3, 6, 8, 0, 0, 0, 0, 0, 0]

ia = [0, 1, 3, 0, 2, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

a = [-1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

b = [-5.0, -10.0, 7.0, 5.0, 0, 0]

c = [-1.0, -4.0, -9.0, 0, 0, 0, 0, 0, 0]


m = 4
n = 3
nz = 8

print("ka: ")
for i in range(0,len(ka)):
    print("i={} ka[{}]={}".format(i, i, ka[i]))

print("A: ")
for j  in range(0, n):
    for k in range(ka[j], ka[j+1]):
        print("j:{:d} k:{:d} iA[k]:{:5d} A[k]:{:10.5f} ".format(j, k, ia[k], a[k]))
    print("")

print("")
print("b: ")
for i in range(0, len(b)):
    print("{:10.5f} ".format(b[i]))
print("")
print("c: ")
for i in range(0, len(c)):
    print("{:10.5f} ".format(c[i]))
print("")



f = 0
x = MALLOC_DOUBLE(nz+1)
y = MALLOC_DOUBLE(nz+1)
w = MALLOC_DOUBLE(nz+1)
z = MALLOC_DOUBLE(nz+1)

solver(m, n, nz, ia, ka, a, b, c, f, x, y, w, z)
'''

#pd_reset()

