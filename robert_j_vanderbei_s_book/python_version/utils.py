import random
import math
import array
import time
from valind import *

class LP:
    m = 0  # int m;		/* number of rows */
    n = 0  # int n;		/* number of columns */
    nz = 0  # int nz;		/* number of nonzeros */
    A = []  # double *A;	/* pointer to array of nonzero values in A */
    iA = []  # int *iA;	/* pointer to array of corresponding row indices */
    kA = []  # int *kA;	/* pointer to array of indices into A (and iA)
    # indicating where each new column of A begins */
    b = []  # double *b;	/* pointer to array containing right-hand side */
    c = []  # double *c;	/* pointer to array containing objective function */
    f = 0.  # double f;	/* fixed adjustment to objective function */
    r = []  # double *r;	/* pointer to array containing range vector */
    l = []  # double *l;	/* pointer to array containing lower bounds */
    u = []  # double *u;	/* pointer to array containing upper bounds */
    varsgn = []  # int *varsgn;	/* array indicating which variables were declared to
    # be non-positive	*/
    rowlab = []  # char **rowlab;	/* array of strings containing row labels */
    collab = []  # char **collab;	/* array of strings containing column labels */

    qnz = 0  # int qnz;	/* number of nonzeros in lower triangle	of Q */
    Q = []  # double *Q;	/* pointer to array of nonzero values of Q */
    iQ = []  # int *iQ;	/* pointer to array of corresponding row indices */
    kQ = []  # int *kQ;	/* pointer to array of indices into Q (and iQ)
    # indicating where each new column of Q begins */

    At = []  # double *At;	/* pointer to array of nonzero values in At */
    iAt = []  # int *iAt;	/* pointer to array of corresponding row indices */
    kAt = []  # int *kAt;	/* pointer to array of indices into At (and iAt) */

    bndmark = []  # int *bndmark;	/* pointer to array of bound marks */
    rngmark = []  # int *rngmark;	/* pointer to array of range marks */

    w = []  # double *w;	/* pointer to array containing primal surpluses	*/
    x = []  # double *x;	/* pointer to array containing primal solution */
    y = []  # double *y;	/* pointer to array containing dual solution */
    z = []  # double *z;	/* pointer to array containing dual slacks */
    p = []  # double *p;	/* pointer to array containing range slacks */
    q = []  # double *q;	/* pointer to array containing dual range slacks */
    s = []  # double *s;	/* pointer to array containing dual for	ub slacks */
    t = []  # double *t;	/* pointer to array containing upper bound slacks */
    v = []  # double *v;	/* pointer to array containing dual for	range (w) */
    g = []  # double *g;	/* pointer to array containing lower bound slacks */
    ub = []  # double *ub;	/* pointer to array containing shifted upper bounds */

    max = 0  # int max;	/* max = -1, min = 1 */
    inftol = 0.  # double inftol;	/* infeasibility tolerance */
    inftol2 = 0.  # double inftol2;	/* infeasibility for stopping rule */
    sf_req = 0  # int sf_req;	/* significant figures requested */
    itnlim = 0  # int itnlim;	/* iteration limit */
    timlim = 0.  # double timlim;	/* time limit */
    verbose = 0  # int verbose;	/* level of verbosity */
    epssol = 0.  # double epssol;	/* epsilon tolerance in f/b solve */
    epsnum = 0.  # double epsnum;	/* epsilon tolerance in num fact */
    epscdn = 0.  # double epscdn;	/* epsilon tolerance for conditioning */
    stablty = 0.  # double stablty;	/* mixing factor for stability */
    method = 0  # int method;	/* reordering method */
    dense = 0  # int dense;	/* threshold for dense columns/rows */
    pdf = 0  # int pdf;	/* order to favor primal (ADA^T) or dual (A^TDA) */
    name = " "  # char name[15];	/* string containing problem name */
    obj = " "  # char obj[11];	/* string containing objective function	name */
    rhs = " "  # char rhs[11];	/* string containing right-hand	side name */
    ranges = " "  # char ranges[11];/* string containing range set name */
    bounds = " "  # char bounds[11];/* string containing bound set name */

    tier = []  # int *tier;	/* tier for factorization priorities */

    param = []  # char **param;	/* array of strings containing user parameters */
    np = 0  # int np;		/* number of user parameters */

    # int  (*stopping_rule)(void *);/* pointer to stopping	rule fcn */
    # void (*init_vars)(void *);    /* pointer to initialization fcn */

    # void (*h_init)(void *);       /* pointer to initialization hook fcn */
    # void (*h_update)(void *);     /* pointer to update hook fcn */
    # void (*h_step)(void *);       /* pointer to step hook fcn */

    iter = 0  # int    iter;	    /* current iteration number	*/
    elaptime = 0.  # double elaptime;    /* elapsed time */
    pres = 0.  # double pres;	    /* primal residual (i.e. infeasibility) */
    dres = 0.  # double dres;	    /* dual   residual (i.e. infeasibility) */
    sigfig = 0  # int    sigfig;	    /* significant figures */
    primal_obj = 0.  # double primal_obj;  /* primal objective	value */
    dual_obj = 0.  # double dual_obj;    /* dual   objective	value */


# } LP;

# lp = LP()

# pdf = 0

# from ldlt.py
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

# epssol = _EPSSOL;
# epsnum = _EPSNUM;
# epscdn = _EPSCDN;
# epsdiag= _EPSDIAG;
# stablty= _STABLTY;
# method = _MD;
# dense  = _DENSE;
# pdf    = _UNSET;

# denwin = 0

UNSET = 0
PRIMAL = 1
DUAL = 2

FINITE = 0x1
INFINITE = 0x2
UNCONST = 0x4
FREEVAR = 0x1
BDD_BELOW = 0x2
BDD_ABOVE = 0x4
BOUNDED = 0x8

HUGE_VAL = 10000000000000.0
EPS = 1.0e-8
EPS1 = 1.0e-8
EPS2 = 1.0e-12
EPS3 = 1.0e-10
MAX_ITER = 1000000

# a = []
# tag = []
# link = []
# currtag = 1

a_Nt = []
tag_Nt = []
link_Nt = []
currtag = 1

# lurefac functions

EPSSOL = 1.0e-5;
EPSNUM = 1.0e-9;
NOREORD = 0;
MD = 1;

# rank = 0

# colperm = []
# icolperm = []
# rowperm = []
# irowperm = []


# degL = []
# degLt = []

# L = []
# Lt = []

# degU = []
# degUt = []

# U = []
# Ut = []

# #diag = []

# newcol = []

# inewcol = []
# nnewcol = 0;
# nr = 0;

# perm = []
# iperm = []
# row_list = []

# ngauss = []

# gauss = []

# rows = []

# col_outs = []

# imaxs = []

cumtime = 0.0;
ocumtime = 0.0;

TRUE = 1;
FALSE = 0;


# y = []
# yy = []

# iwork = []
# dwork = []
# pwork = []

# call = 0;
# eps = 0.0


def pd_reset():
    global HUGE_VAL
    global EPS
    global EPS1
    global EPS2
    global EPS3
    global MAX_ITER
    global a
    global tag
    global link
    global currtag
    global a_Nt
    global tag_Nt
    global link_Nt
    global EPSSOL
    global EPSNUM
    global NOREORD
    global MD
    global rank
    global colperm
    global icolperm
    global rowperm
    global irowperm
    global degL
    global degLt
    global L
    global Lt
    global degU
    global degUt
    global U
    global Ut
    global diag
    global newcol
    global inewcol
    global nnewcol
    global nr
    global perm
    global iperm
    global row_list
    global ngauss
    global gauss
    global rows
    global col_outs
    global imaxs
    global cumtime
    global ocumtime
    global TRUE
    global FALSE
    global y
    global yy
    global iwork
    global dwork
    global pwork
    global call
    global eps
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



def REALLOC_INT(x, size):
    if (len(x) == size):
        return x
    y = array.array('i', [0] * size)
    if len(x) > 0:
        for i in range(0, len(x)):
            if (i >= size):
                break;
            y[i] = x[i];
    x.clear()
    for i in range(0, len(y)):
        x.append(y[i])


def REALLOC_DOUBLE(x, size):
    if (len(x) == size):
        return x
    y = array.array('d', [0.0] * size)
    if len(x) > 0:
        for i in range(0, len(x)):
            if (i >= size):
                break;
            y[i] = x[i];
    x.clear()
    for i in range(0, len(y)):
        x.append(y[i])

def REALLOC_VALIND(x, size):  # int[] REALLOC(int[] x, int size)
    if (len(x) == size):
        return x
    y = []
    for i in range(0, size):
        if (i < len(x)):
            y.append(x[i])
        else:
            y.append(valind())
    x.clear()
    for i in range(0, len(y)):
        x.append(y[i])

def REALLOC2D(x, size):
    if (len(x) == size):
        return x
    y = []
    for i in range(0, size):
        if i < len(x):
            y.append(x[i])
        else:
            y.append([])
    x.clear()
    for i in range(0, len(y)):
        x.append(y[i])


def MALLOC_DOUBLE(size):
    y = []
    for i in range(0, size):
        y.append(0.0)
    return y


def MALLOC_INT(size):
    y = []
    for i in range(0, size):
        y.append(0)
    return y

def MALLOC_VALIND(size):
    v1 = []
    for i in range(0, size):
        v1.append(valind())

    return v1

def MALLOC2D(size):
    v1 = []
    for i in range(0, size):
        v1.append([])
    return v1




random.seed(100)


#def drand48():
#    return random.uniform(0, 1)
def drand48():
	return 0.

def high(v):
    return v


def sqrt(v):
    return math.sqrt(v);


def fabs(v):
    return math.fabs(v);


def Display_Solution(m, basics, x):
    print('SOLUTION:\n\n')

    for i in range(0, m):
        print('  X[:d] = :f\n', basics[i], high(x[i]));


def Nt_times_y(n, at, iat, kat, basicflag, y, iy, ny, yN, iyN, pnyN):
    global a_Nt
    global tag_Nt
    global link_Nt
    global currtag

    if len(a_Nt) == 0:
        a_Nt = MALLOC_DOUBLE(n)
    if len(tag_Nt) == 0:
        tag_Nt = MALLOC_INT(n)
    if len(link_Nt) == 0:
        link_Nt = MALLOC_INT(n + 2 + 1)

    jj = -1;
    for k in range(0, ny):

        i = iy[k];
        for kk in range(kat[i], kat[i + 1]):

            j = iat[kk];
            if (basicflag[j] < 0):

                if (tag_Nt[j] != currtag):
                    a_Nt[j] = 0.0;
                    tag_Nt[j] = currtag;
                    link_Nt[jj + 1] = j;
                    jj = j;

                a_Nt[j] += y[k] * at[kk];

    link_Nt[jj + 1] = n;
    currtag = currtag + 1;

    k = 0;
    jj = link_Nt[-1 + 1]
    while (jj < n):

        if (ABS(a_Nt[jj]) > EPS1):
            yN[k] = a_Nt[jj];
            iyN[k] = -basicflag[jj] - 1;
            k = k + 1;

        jj = link_Nt[jj + 1]

    pnyN = k;
    return pnyN;


def ratio_test(dy, idy, ndy, y, ybar, mu):
    global HUGE_VAL
    jj = -1
    min = HUGE_VAL
    for k in range(0, ndy):
        if (dy[k] > EPS1):
            j = idy[k];
            if ((y[j] + mu * ybar[j]) / dy[k] < min):
                min = (y[j] + mu * ybar[j]) / dy[k];
                jj = j;
                kk = k;
    return jj




def cycperm(start, end, perm, iperm):

    if (start < end):
        k = perm[start];
        for j in range(start, end):
            perm[j] = perm[j + 1];
        perm[end] = k;
        for j in range(start, end + 1):
            iperm[perm[j]] = j;
    elif (start > end):
        k = perm[start];
        for j in range(start, end, -1):
            perm[j] = perm[j - 1];
        perm[end] = k;
        for j in range(start, end - 1, -1):
            iperm[perm[j]] = j;



def clock():
    return time.clock()


def ABS(v):
    return abs(v)


def MAX(v1, v2):
    return max(v1, v2);

def MIN(v1, v2):
    return min(v1, v2);

def Bswap(v, i, j):
    temp = v[i];
    v[i] = v[j];
    v[j] = temp;


statmsg = []
statmsg.append("optimal solution")  # /* 0 */
statmsg.append("primal unbounded")  # /* 1 */
statmsg.append("primal infeasible")  # /* 2 */
statmsg.append("dual unbounded")  # /* 3 */
statmsg.append("dual infeasible")  # /* 4 */
statmsg.append("iteration limit")  # /* 5 */
statmsg.append("infinite lower bounds - not implemented")  # /* 6 */
statmsg.append("suboptimal solution")  # /* 7 */
# \};
