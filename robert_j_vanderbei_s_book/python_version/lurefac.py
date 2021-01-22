from utils import *
from ldlt import *
from tree import *

nr = 0;

colperm = []
icolperm = []
rowperm = []
irowperm = []
diag = []
degL = []
L = []
degUt = []
Ut = []
degLt = []
Lt = []
degU = []
U = []
newcol = []
inewcol = []
cumtime = 0.0;
y = []
yy = []
tag = []
currtag=1;

tag_bs = []
y_bs = []
yy_bs = []
currtag_bs=1;

iwork = []
dwork = []
pwork = []

perm = []
iperm = []
rows = []
col_outs = []
imaxs = []
row_list = []
ngauss = []
gauss = []

call=0;

#public void lufac( int m, int[] kA, int[] iA, double[] A, int[] basis, int v )

def lufac(m, kA, iA, A, basis, v):
    global HUGE_VAL
    global EPS
    global EPS1
    global EPS2
    global EPS3
    global MAX_ITER
    global a
    global link
    global currtag
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

    # for debugging only...
    global col2, kk, B, U, Ut, row, kkk

    method=MD

    hkey = []
    heap = []
    iheap = []
    iwork = []
    iwork2 = []

    degB = []
    degBt = []

    B = []
    Bt = []

    starttime = clock();

    col = 0;

        #+---------------------------------------------------------+
        #| Set/reset number of refactorizations to 0.              |


    if (nr > 0):
        nr = 0;

        #+---------------------------------------------------------+
        #| allocate space for perm, iperm, and diag.               |

    if len(colperm) == 0:
        colperm = MALLOC_INT(m)
    else:
        REALLOC_INT(colperm, m)
    if len(icolperm) == 0:
        icolperm = MALLOC_INT(m)
    else:
        REALLOC_INT(icolperm, m)

    if len(rowperm) == 0:
        rowperm = MALLOC_INT(m)
    else:
        REALLOC_INT(rowperm, m)


    if len(irowperm) == 0:
        irowperm = MALLOC_INT(m)
    else:
        REALLOC_INT(irowperm, m)


    if len(diag) == 0:
        diag = MALLOC_DOUBLE(m)
    else:
        REALLOC_DOUBLE(diag, m)

    #+---------------------------------------------------------+
    #| allocate space for work arrays.                         |

    degB = MALLOC_INT(m)
    degBt = MALLOC_INT(m)
    hkey = MALLOC_INT(m)
    heap = MALLOC_INT(m)
    iheap = MALLOC_INT(m)
    iwork = MALLOC_INT(m)
    iwork2 = MALLOC_INT(m)

    heap = MALLOC_INT(m+1);

    #+---------------------------------------------------------+
    #| calculate degrees in B and Bt                           |

    for i in range(0,m):
        degBt[i] = 0;
    for i in range(0,m):
        degB[i] = kA[ basis[i]+1 ] - kA[ basis[i] ];
        for k in range(kA[ basis[i] ], kA[ basis[i]+1 ]):
            degBt[ iA[k] ] = degBt[ iA[k] ] + 1;

    #+---------------------------------------------------------+
    #| initialize B and Bt                                     |

    B = MALLOC2D(m);
    Bt = MALLOC2D(m);
    for i in range(0,m):
        B[i] = MALLOC_VALIND(degB[i]);
        Bt[i] = MALLOC_VALIND(degBt[i]);

    for i in range(0,m):
        iwork[i] = 0;
    for j in range(0,m):
        kkk = 0;
        for k in range(kA[ basis[j] ], kA[ basis[j]+1 ]):
            row = iA[k];
            kk  = iwork[row];
            B[j][kkk].i = row;
            B[j][kkk].d = A[k];
            Bt[row][kk].i = j;
            Bt[row][kk].d = A[k];
            iwork[row] = iwork[row] + 1;
            kkk = kkk + 1;

    #+---------------------------------------------------------+
    #| miscellaneous initializations.                          |

    for i in range(0, m):
        icolperm[i] = -1;
        irowperm[i] = -1;
        iwork[i] = 0;
        iwork2[i] = -1;

    rank = m;
    tag = 0;
    Bnz = 0;
    Btnz = 0;

    #+---------------------------------------------------------+
    #| hkey encodes the tie-breaking rule - currently the rule |
    #| is somewhat random.  to make it first occuring minimum, |
    #| change the formula to:                                  |
    #|       hkey[node] = degree[node]*m + node;               |
    #| warning: with this definition of hkey, there is the     |
    #| possibility of integer overflow on moderately large     |
    #| problems.                                               |
    #|                                                         |

    for j in range(0,m):
        if (method == MD):
            hkey[j] = degBt[j];
        else:
            hkey[j] = j;

        if (method == MD and hkey[j]==0):
            hkey[j]=m+1;

    #+---------------------------------------------------------+
    #| set up heap structure for quickly accessing minimum.    |

    heapnum = m;
    for j in range(m-1, -1, -1):
        cur = j+1;
        iheap[j] = cur;
        heap[cur] = j;
        hfall( heapnum, hkey, iheap, heap, cur );

    #+---------------------------------------------------------+
    #| the min degree ordering loop                            |

    for i in range(0, m):
        #+------------------------------------------------+
        #|  select row with min column degree             |

        row    = heap[1];
        rowdeg = degBt[row];

        if (rowdeg == 0):
            print("singular matrix. rank deficiency = " +  str(m-i) + "\n");
            rank = i;
            break;

        #+------------------------------------------------+
        #|  select pivot element from this row by         |
        #|  choosing nonzero whose col is of minimal      |
        #|  degree                                        |

        if (method == MD):
            coldeg = m+1;
            for k in range(0, rowdeg):
                if ( degB[ Bt[row][k].i ] < coldeg and ABS(Bt[row][k].d ) > EPSNUM ):
                    col    = Bt[row][k].i;
                    coldeg = degB[col];
            if (coldeg == m+1):
                hkey[row]=m+2;
                hfall( heapnum, hkey, iheap, heap, iheap[row] );
                if (hkey[heap[1]] == m+2):
                    print("singular matrix. rank deficiency = {:d}\n".format(m-i));
                    rank = i;
                    break;
                else:
                    i = i-1
                    continue;
        else:
            col    = Bt[row][degBt[row]-1].i;
            coldeg = degB[col];
            for k in range(0, rowdeg):
                if ( Bt[row][k].i == row ):
                    col    = row;
                    coldeg = degB[col];
                    break;

        #+------------------------------------------------+
        #|  update permutation information                |

        colperm[i] = col;
        icolperm[col] = i;

        rowperm[i] = row;
        irowperm[row] = i;

        #+------------------------------------------------+
        #|  remove eliminated elements from B and Bt      |

        for k in range(0, coldeg):
            row2 = B[col][k].i;
            kk = 0
            while Bt[row2][kk].i != col:
                kk = kk + 1

            if (row2 != row):
                degBt[row2] = degBt[row2] - 1;
                Bswap( Bt[row2], degBt[row2], kk );

        for k in range(0, rowdeg):
            col2 = Bt[row][k].i;
            kk = 0
            while B[col2][kk].i != row:
                kk = kk + 1;

            degB[col2] = degB[col2] - 1;
            Bswap( B[col2], degB[col2],  kk );

        kk=0
        while(Bt[row][kk].i != col):
            kk = kk + 1;

        diag[i] = Bt[row][kk].d;

        degBt[row] = degBt[row] - 1;
        Bswap( Bt[row], degBt[row], kk );

        #+------------------------------------------------+
        #|  update heap                                   |

        okey = hkey[col];
        heap[1] = heap[heapnum];
        iheap[heap[1]] = 1;
        heapnum = heapnum - 1;
        if (okey < hkey[heap[1]]):
            hfall(heapnum, hkey, iheap, heap, 1);

        #+------------------------------------------------+
        #|  generate fillin and update elements           |

        for k in range(0, degB[col]):
            row2 = B[col][k].i;
            tag = tag + 1;
            for kk in range(0, degBt[row2]):
                col2 = Bt[row2][kk].i;
                iwork[ col2] = tag;
                iwork2[col2] = kk;
            for kk in range(0, degBt[row]):
                col2 = Bt[row][kk].i;
                if ( iwork[col2] == tag ):
                    kkk = iwork2[col2];
                    Bt[row2][kkk].d = Bt[row2][kkk].d - B[col][k].d * Bt[row][kk].d / diag[i];
                    if ( ABS(Bt[row2][kkk].d) <= 1.0e-12):
                        degBt[row2] = degBt[row2] - 1;
                        col3 = Bt[row2][degBt[row2]].i;
                        iwork [col3] = iwork [col2];
                        iwork2[col3] = iwork2[col2];
                        Bswap( Bt[row2], degBt[row2], kkk );
                else:
                    deg = degBt[row2];
                    REALLOC_VALIND( Bt[row2], deg+1);
                    Bt[row2][deg].i = col2;
                    Bt[row2][deg].d = - B[col][k].d * Bt[row][kk].d / diag[i];
                    degBt[row2] = degBt[row2]+1;

        for k in range(0, degBt[row]):
            col2 = Bt[row][k].i;
            tag = tag + 1;
            for kk in range(0, degB[col2]):
                row2 = B[col2][kk].i;
                iwork[ row2] = tag;
                iwork2[row2] = kk;
            for kk in range(0, degB[col]):
                row2 = B[col][kk].i;
                if ( iwork[row2] == tag ):
                    kkk = iwork2[row2];
                    B[col2][kkk].d = B[col2][kkk].d - B[col][kk].d * Bt[row][k].d / diag[i];
                    if ( ABS(B[col2][kkk].d) <= 1.0e-12):
                        degB[col2] = degB[col2] - 1;
                        row3 = B[col2][degB[col2]].i;
                        iwork [row3] = iwork [row2];
                        iwork2[row3] = iwork2[row2];
                        Bswap( B[col2], degB[col2], kkk );
                else:
                    deg = degB[col2];
                    REALLOC_VALIND(B[col2], deg+1);
                    B[col2][deg].i = row2;
                    B[col2][deg].d = (-1) * B[col][kk].d * Bt[row][k].d/diag[i];
                    degB[col2] = degB[col2]+1;

        #+------------------------------------------------+
        #|  adjust heap                                   |

        for k in range(0, degB[col]):
            row2 = B[col][k].i;
            if (method == MD):
                hkey[row2] = degBt[row2];
                if (hkey[row2]==0):
                    hkey[row2]=m+1;
            else:
                hkey[row2] = row2;
            hrise( hkey, iheap, heap, iheap[row2] );
            hfall( heapnum, hkey, iheap, heap, iheap[row2] );

    #+------------------------------------------------+
    #|  process dependent rows/cols                   |

    i = rank;
    for col in range(0, m):
        if (icolperm[col] == -1):
            colperm[i] = col;
            icolperm[col] = i;
            degB[col] = 0;
            i = i + 1;

    i = rank;
    for row in range(0,m):
        if (irowperm[row] == -1):
            rowperm[i] = row;
            irowperm[row] = i;
            degBt[row] = 0;
            i = i + 1;

    for i in range(rank, m):
        diag[i] = 0.0;

    #+------------------------------------------------+
    #|  divide each column of L by diagonal           |

    for col in range(0, m):
        for k in range(0, degB[col]):
            i = icolperm[col];
            B[col][k].d = B[col][k].d / diag[i];

    #+---------------------------------------------------------+
    #| calculate and print statistics.                         |

    narth = 0.0e0;
    for i in range(0, m):
        k = degB[i];
        narth += k*k;
        Bnz  += k;
        k = degBt[i];
        narth += k*k;
        Btnz += k;
    narth = narth + 3*Bnz + 3*Btnz + 2*m;

    if (v == TRUE):
        print("{:9d}   {:9d} {:15.0f}".format( Bnz, Btnz, narth));

    if len(degL) == 0:
        degL = MALLOC_INT(m);
    else:
        REALLOC_INT(degL, m);

    if len(L) == 0:
        L = MALLOC2D(m);
    else:
        REALLOC2D(L, m);

    if len(degUt) == 0:
        degUt = MALLOC_INT(m)
    else:
        REALLOC_INT(degUt, m);

    if len(Ut) == 0:
        Ut = MALLOC2D(m);
    else:
        REALLOC2D(Ut, m);
    for i in range(0, m):
        col = colperm[i];
        degL[i] = degB[col];
        REALLOC_VALIND( L[i], degL[i]);
        for k in range(0, degL[i]):
            L[i][k].d = B[col][k].d;
            L[i][k].i = irowperm[ B[col][k].i ];
    for i in range(0,m):
        row = rowperm[i];
        degUt[i] = degBt[row];
        REALLOC_VALIND( Ut[i], degUt[i]);
        for k in range(0, degUt[i]):
            Ut[i][k].d =           Bt[row][k].d;
            Ut[i][k].i = icolperm[ Bt[row][k].i ];
    if len(degLt) == 0:
        degLt = MALLOC_INT(m)
    else:
        REALLOC_INT(degLt, m);
    if len(Lt) == 0:
        Lt = MALLOC2D(m);
    else:
        REALLOC2D(Lt, m);

    if len(degU) == 0:
        degU = MALLOC_INT(m)
    else:
        REALLOC_INT(degU, m);


    if len(U) == 0:
        U = MALLOC2D(m);
    else:
        REALLOC2D(U, m);
    ratnum(m, degL , L , degLt, Lt);
    ratnum(m, degUt, Ut, degU , U );

    endtime = clock()
    cumtime = cumtime + (endtime - starttime);


#+-----------------------------------------------------------------+
#| LU factorization.                                               |
#| Input:                                                          |
#|    m          number of rows (= number of columns)              |
#|    kA, iA, A  three array sparse representation of m by n       |
#|               matrix A                                          |
#|    basis      array of m out of n of the column indices of A    |
#|               specify a submatrix B of A                        |
#| Output:                                                         |
#|    static global variables (only visible to other routines in   |
#|    this file:                                                   |
#|                                                                 |
#|    rank       rank of B                                         |
#|    B, degB    ragged array representation of L                  |
#|    Bt, degBt  ragged array representation of U transpose        |
#|               without its diagonal                              |
#|    diag       diagonal entries of U                             |
#|    colperm, icolperm, rowperm, irowperm                         |
#|               column and row permutations and their inverses    |

def refactor(m, kA, iA, A, basics, col_out, v):
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

    Utnz = 0;

    rffactor = 1.0;

    #+------------------------------------------------------+
    #| Check if it is time to refactor from scratch         |

    call = call + 1;
    if (col_out < 0 or call <= 1):
        ocumtime = 0.0;
        cumtime  = 0.0;
        lufac( m, kA, iA, A, basics, v );
        cumtime  = cumtime * rffactor;
        from_scratch = TRUE;
        return from_scratch;
    if ( call > 3 and cumtime/call >= ocumtime/(call-1) ):
        ocumtime = 0.0;
        cumtime  = 0.0;
        call = 1;
        lufac( m, kA, iA, A, basics, v );
        cumtime  = cumtime * rffactor;
        from_scratch = TRUE;
        return from_scratch;

    ocumtime  = cumtime;
    starttime = clock();

    #+------------------------------------------------------+
    #| Allocate storage for work arrays                     |

    if len(iwork) == 0:
        iwork = MALLOC_INT(m)
    if len(dwork) == 0:
        dwork = MALLOC_DOUBLE(m)
    if len(pwork) == 0:
        pwork = MALLOC2D(m);

    #+------------------------------------------------------+
    #| Put col_out into `new' indices                       |

    col_out = icolperm[col_out];

    #+------------------------------------------------------+
    #| Compute largest row index for new column             |

    if len(inewcol) == 0:
        print('ERROR: refactoring before bsolving \n');
        sys.exit(0);

    imax=0;
    for k in range(0, nnewcol):
	    imax = MAX(imax, inewcol[k]);

    if (imax < col_out):
	    print("singular matrix \n");
	    from_scratch = FALSE;
	    return from_scratch;

    #+------------------------------------------------------+
    #| Insert newcol into col_out column of U (and Ut)      |
    #|                                                      |
    #|             0 1 2 3 4 5                              |
    #|                                                      |
    #|          0  x * x x x x                              |
    #|          1    * x x x x                              |
    #|    U  =  2    * x x x x (here col_out=1 and imax=4)  |
    #|          3    *   x x x                              |
    #|          4    *     x x                              |
    #|          5            x                              |

    #+ first remove oldcol from Ut |
    for k in range(0, degU[col_out]):
        row = U[col_out][k].i;
        for kk in range(0, degUt[row]):
            if Ut[row][kk].i == col_out: break
                 #+ INEFFICIENT |
        if (kk < degUt[row]):
            degUt[row] = degUt[row] - 1;
            Bswap( Ut[row], degUt[row], kk );

    degU[col_out] = nnewcol;
    REALLOC_VALIND(U[col_out], nnewcol);
    kkkk = 0;
    diag[col_out] = 0.0;

    for k in range(0, nnewcol):
        row = inewcol[k];
        val =  newcol[k];
        if (row != col_out):
            U[col_out][kkkk].i = row;
            U[col_out][kkkk].d = val;
            kkkk = kkkk + 1;

            kkk = degUt[row];
            degUt[row] = degUt[row] + 1;
            REALLOC_VALIND(Ut[row], degUt[row]);
            Ut[row][kkk].i = col_out;
            Ut[row][kkk].d = val;
        else:
            diag[row] = val;
    degU[col_out] = kkkk;

    #+------------------------------------------------------+
    #| Allocate storage for permutation arrays and shift    |
    #| so that indexing begins at col_out                   |

    REALLOC2D(perm, nr+1);
    REALLOC2D(iperm, nr+1);
    perm[nr] = MALLOC_INT(imax+1);
    iperm[nr] = MALLOC_INT(imax+1);


    REALLOC_INT(rows, nr+1);
    REALLOC_INT(col_outs, nr+1);
    REALLOC_INT(imaxs, nr+1);

    #+------------------------------------------------------+
    #| Initialize permutation arrays so that col_out is     |
    #| cyclically permuted to imax.  After permutation:     |
    #|                                                      |
    #|             0 2 3 4 1 5                              |
    #|                                                      |
    #|          0  x x x x * x                              |
    #|    U  =  2    x x x * x (here col_out=1 and imax=4)  |
    #|          3      x x * x                              |
    #|          4        x * x                              |
    #|          1    x x x * x                              |
    #|          5            x                              |
    #|                                                      |

    for j in range(col_out, imax):
        perm[nr][j]   = j+1;
        iperm[nr][j+1] = j;
    perm[nr][imax]    = col_out;
    iperm[nr][col_out] = imax;

    #+------------------------------------------------------+
    #| Look for singleton columns/rows and permute columns  |
    #| to upper-left and rows to lower-right position in    |
    #| bump.  Don't forget that the diagonal is stored      |
    #| separately in diag[] and that this contributes one   |
    #| nonzero to each column/row investigated.             |

    bumpstart = col_out;
    bumpend   = imax;
    changes = 1
    while changes > 0:  #{ changed from do..while
        changes = 0;

        # /*-----------------------------------------------------+
        # | First look for columns.                              |
        # |                                                      |
        # |       0 1 2 3 4 5          0 3 1 2 4 5               |
        # |                                                      |
        # |    0  x x x x * x       0  x x x x * x               |
        # |    1    x x   * x       3    x     * x               |
        # |    2      x   * x  -->  1      x x * x               |
        # |    3        x * x       2        x * x               |
        # |    4    x x   * x       4      x x * x               |
        # |    5            x       5            x               |
        # |                                                      |

        for j in range(bumpstart, bumpend):
            col = perm[nr][j];
            cnt = 0;
            for k in range(0, degU[col]):
                Ui = U[col][k].i;
                if Ui >= col_out and Ui <= imax:
                    row = iperm[nr][ Ui ];
                    if (bumpstart <= row and row <= bumpend):
                        cnt = cnt + 1;

            if (cnt == 0):
                cycperm(j, bumpstart, perm[nr], iperm[nr]);
                bumpstart = bumpstart + 1;
                changes = changes + 1;
        #+------------------------------------------------------+
        #| Now look for rows.                                   |
        #|                                                      |
        #|       0 1 2 3 4 5          0 2 3 4 1 5               |
        #|                                                      |
        #|    0  x x x x * x       0  x x x * x x               |
        #|    1    x       x       2    x x *   x               |
        #|    2      x x * x  -->  3      x *   x               |
        #|    3        x * x       4    x x * x x               |
        #|    4    x x x * x       1          x x               |
        #|    5            x       5            x               |
        #|                                                      |

        for i in range(bumpend-1, bumpstart-1, -1):
            row = perm[nr][i];
            cnt = 0;
            for k in range(0, degUt[row]):
                Uti = Ut[row][k].i;
                if (Uti >= col_out and Uti <= imax):
                    col = iperm[nr][ Uti ];
                    if (bumpstart <= col and col <= bumpend):
                        cnt = cnt + 1;
            if (cnt == 0):
                cycperm(i, bumpend, perm[nr], iperm[nr]);
                bumpend = bumpend - 1;
                changes = changes + 1;

        #+------------------------------------------------------+
        #| Permute rows/columns of U and Ut.                    |

        #+------------------------------------------------------+
        #| Permute columns of U and diag.                       |

    for j in range(col_out, imax+1):
        dwork[j] = diag[j];
        iwork[j] = degU[j];
        pwork[j] = U[j];

    for j in range(col_out, imax+1):
	    diag[j]  = dwork[perm[nr][j]];
	    degU[j]  = iwork[perm[nr][j]];
	    U[j]     = pwork[perm[nr][j]];

    #+------------------------------------------------------+
    #| Permute rows of U.                                   |

    for j in range(col_out, m):
        for k in range(0, degU[j]):
            row = U[j][k].i;
            if (col_out <= row and row <= imax):
                U[j][k].i = iperm[nr][row];

    #+------------------------------------------------------+
    #| Permute rows of Ut.                                  |

    for i in range(col_out, imax+1):
        iwork[i] = degUt[i];
        pwork[i] = Ut[i];
    for i in range(col_out, imax+1):
	    degUt[i]  = iwork[perm[nr][i]];
	    Ut[i]     = pwork[perm[nr][i]];

    #+------------------------------------------------------+
    #| Permute columns of Ut.                               |

    for i in range(0, imax+1):
        for k in range(0, degUt[i]):
            col = Ut[i][k].i;
            if (col_out <= col and col <= imax):
                Ut[i][k].i = iperm[nr][col];

    #+------------------------------------------------------+
    #| Record bump row for later use.                       |

    row          = bumpend;
    rows[nr]     = row;
    col_outs[nr] = col_out;
    imaxs[nr]    = imax;

    if len(y) == 0:
        MALLOC_DOUBLE(y, m)
    if len(tag) == 0:
        MALLOC_INT(tag, m)

    #+------------------------------------------------------+
    #| Scatter bump row into a dense vector.                |

    for k in range(0, degUt[row]):
	    col = Ut[row][k].i;
	    y[col] = Ut[row][k].d;
	    tag[col] = currtag;
	    addtree(col);
    y[row] = diag[row];
    tag[row] = currtag;
    addtree(row);

    #+------------------------------------------------------+
    #| Remove bump row from U.                              |

    for k in range(0,degUt[row]):
        col = Ut[row][k].i;
        for kk in range(0, degU[col]):
            if (U[col][kk].i == row):
                break;
        if (kk < degU[col]):
            degU[col] = degU[col] - 1;
            Bswap(U[col], degU[col], kk);

    #+------------------------------------------------------+
    #| Do Gaussian elimination on scatter vector.           |

    REALLOC2D( row_list, nr+1);
    row_list[nr] = MALLOC_INT(m);
    REALLOC_INT(ngauss, nr+1);
    REALLOC2D( gauss, nr+1);
    gauss[nr] = MALLOC_DOUBLE(m);

    k=0;
    col=getfirst()
    while col<bumpend:
        row2 = col;
        row_list[nr][k] = row2;
        gauss[nr][k] = y[col] / diag[row2];
        for kk in range(0, degUt[row2]):
            col2 = Ut[row2][kk].i;
            if (tag[col2] != currtag):
                y[col2] = 0.0;
                tag[col2] = currtag;
                addtree(col2);
            y[col2] = y[col2] - gauss[nr][k] * Ut[row2][kk].d;
        k = k + 1;
        col = getnext()
    if (col != bumpend):
        print("ERROR: col != bumpend \n");
    ngauss[nr] = k;
    REALLOC_DOUBLE(gauss[nr], k);
    REALLOC_INT(row_list[nr], k);

    #+------------------------------------------------------+
    #| Add eliminated row to U.  kk counts nonzeros in      |
    #| eliminated row.                                      |

    diag[col] = y[col];

    kk = 0;
    col=getnext()
    while col != -1:
        if ( ABS(y[col])>EPS ):
            k = degU[col];
            degU[col] = degU[col] + 1;
            REALLOC_VALIND( U[col], degU[col]);
            U[col][k].i = row;
            U[col][k].d = y[col];
            kk = kk + 1;
        col = getnext()

    REALLOC_VALIND( Ut[row], kk);

    #+------------------------------------------------------+
    #| Remove bump row from Ut and replace with eliminated  |
    #| row.                                                 |

    k = 0;
    col=getfirst()
    while col != -1:
        if ( col>bumpend and ABS(y[col]) > EPS ):
            Ut[row][k].d = y[col];
            Ut[row][k].i = col;
            k = k + 1;
        col = getnext()
    degUt[row] = k;

    if (k != kk):
        print("ERROR: alloc'ed wrong size for Ut\n");

    currtag = currtag + 1;
    killtree();

    #+------------------------------------------------------+
    #| Apply permutation to colperm and icolperm            |

    for j in range(col_out, imax+1):
        iwork[j] = colperm[j];
    for j in range(col_out, imax+1):
        icolperm[ colperm[ perm[nr][j]]] = j;
    for j in range(col_out, imax+1):
        colperm[iperm[nr][j]] = iwork[j];

    #+------------------------------------------------------+
    #| Increment number of refactorizations.                |

    nr = nr + 1;

    for i in range(0, m):
        k = degUt[i];
        Utnz = Utnz + k;

    if (v == TRUE):
        print('            {:9d} '.format(Utnz), flush=True);

    endtime = clock();
    cumtime = cumtime + (endtime - starttime);

    from_scratch = FALSE;
    return from_scratch;




#+-----------------------------------------------------------------+
#| Forward/backward solve using LU factorization                   |
#| Input:                                                          |
#|    m          dimension of array y                              |
#|    y          array containing right-hand side                  |
#|                                                                 |
#|    static global variables (assumed setup by lufac()):          |
#|                                                                 |
#|    rank       rank of B                                         |
#|    L, degL    ragged array representation of L                  |
#|    Ut, degUt  ragged array representation of U transpose        |
#|               without its diagonal                              |
#|    diag       diagonal entries of U                             |
#|    colperm, icolperm, rowperm, irowperm                         |
#|               column and row permutations and their inverses    |
#| Output:                                                         |
#|                                          -1                     |
#|    y          array containing solution B  y                    |
#|                                                                 |
#|    integer flag indicating whether system is consistent         |

def bsolve(m, sy, iy, pny):
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

    ny = pny;
    consistent = TRUE;
    eps = EPSSOL;

    starttime = clock();

    if (len(y) == 0):
        y = MALLOC_DOUBLE(m)
    if (len(yy) == 0):
        yy = MALLOC_DOUBLE(m)
    if (len(tag) == 0):
        tag = MALLOC_INT(m)
    if (len(newcol) == 0):
        newcol = MALLOC_DOUBLE(m)
    if (len(inewcol) == 0):
        inewcol = MALLOC_INT(m)

    for k in range(0, ny):
        i = irowperm[iy[k]];
        y[i] = sy[k];
        tag[i] = currtag;
        addtree(i);

    if (rank < m):
        eps = EPSSOL * maxv(sy,ny);

    #+------------------------------------------------------+
    #|               -1                                     |
    #|       y  <-  L  y                                    |

    i = getfirst()
    while i < rank and i != -1:
        beta = y[i];
        for k in range(0, degL[i]):
            row = L[i][k].i;
            if (tag[row] != currtag):
                y[row] = 0.0;
                tag[row] = currtag;
                addtree(row);
            y[row] -= L[i][k].d * beta;
        i = getnext()
    #+------------------------------------------------------+
    #| Apply refactorization row operations.                |

    for jr in range(0,nr):
        #+--------------------------------------------------+
        #| Gather sparse vector.                            |
        k=0;
        for j in range(col_outs[jr], imaxs[jr]+1):
            if (tag[j] == currtag):
                sy[k] = y[j];
                iy[k] =   j;
                k = k+1;
                tag[j] = tag[j]-1;
                deltree(j);
        ny = k;

        #+--------------------------------------------------+
        #| Scatter and permute.                             |

        for k in range(0, ny):
            i = iperm[jr][iy[k]];
            y[i] = sy[k];
            tag[i] = currtag;
            addtree(i);

        #+--------------------------------------------------+
        #| Apply row operations.                            |

        row = rows[jr];
        for k in range(0,ngauss[jr]):
            row2 = row_list[jr][k];
            if (tag[row] != currtag):
                y[row] = 0.0;
                tag[row] = currtag;
                addtree(row);
            if (tag[row2] == currtag):
                y[row] -= gauss[jr][k] * y[row2];

    #+------------------------------------------------------+
    #|                                       -1             |
    #| Set aside sparse intermediate vector L  P a  for     |
    #|                                            j         |
    #| refactorization routine.                             |

    nnewcol = 0;
    i=getfirst()
    while (i != -1):
        if ( ABS(y[i]) > EPS ):
            newcol [nnewcol] = y[i];
            inewcol[nnewcol] = i;
            nnewcol = nnewcol + 1;
        i = getnext()

    #+------------------------------------------------------+
    #|               -1                                     |
    #|       y  <-  U  y                                    |

    i = getlast()
    while i >= rank and i != -1:
        if ( ABS( y[i] ) > eps ):
            consistent = FALSE;
        y[i] = 0.0;
        i=getprev()

    while i >= 0:
        beta = y[i]/diag[i];
        for k in range(0, degU[i]):
            row = U[i][k].i;
            if (tag[row] != currtag):
                y[row] = 0.0;
                tag[row] = currtag;
                addtree(row);
            y[row] = y[row] - U[i][k].d * beta;
        y[i] = beta;
        i=getprev()

    ny = 0;
    i = getfirst()
    while i != -1:
        if ( ABS(y[i]) > EPS ):
	        sy[ny] = y[i];
	        iy[ny] = colperm[i];
	        ny = ny + 1;
        i = getnext()

    currtag = currtag + 1;
    killtree();

    endtime = clock();
    cumtime = cumtime + (endtime - starttime);

    return ny


#int     btsolve(
#        int m,
#        double[] sy,
#        int[] iy,
#        int pny //int[] pny
#)

def btsolve(m, sy, iy, pny):
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

    ny = pny;
    consistent = TRUE;
    eps = EPSSOL;

    starttime = clock();

    if len(y) == 0:
        y = MALLOC_DOUBLE(m)
    if len(tag) == 0:
        tag = MALLOC_DOUBLE(m)

    for k in range(0, ny):
        i = icolperm[iy[k]];
        y[i] = sy[k];
        tag[i] = currtag;
        addtree(i);

    if (rank < m):
        eps = EPSSOL * maxv(sy,ny);

    #+------------------------------------------------------+
    #|               -T                                     |
    #|       y  <-  U  y                                    |

    i = getfirst()
    while (i < rank and i != -1):
        beta = y[i]/diag[i];
        for k in range(0,degUt[i]):
            row = Ut[i][k].i;
            if (tag[row] != currtag):
                y[row] = 0.0;
                tag[row] = currtag;
                addtree(row);
            y[row] -= Ut[i][k].d * beta;
        y[i] = beta;
        i=getnext();

    #+------------------------------------------------------+
    #| Apply refactorization row operations.                |

    for jr in range(nr-1, -1, -1):

        #+--------------------------------------------------+
        #| Apply row operations.                            |

        row = rows[jr];
        for k in range(ngauss[jr]-1, -1, -1):
            row2 = row_list[jr][k];
            if (tag[row2] != currtag):
                y[row2] = 0.0;
                tag[row2] = currtag;
                addtree(row2);
            if (tag[row] == currtag):
                y[row2] = y[row2] - gauss[jr][k] * y[row];

        #+--------------------------------------------------+
        #| Gather sparse vector.                            |

        k=0;
        for j in range(col_outs[jr], imaxs[jr]+1):
            if (tag[j] == currtag):
                sy[k] = y[j];
                iy[k] =   j;
                k = k + 1;
                tag[j] = tag[j] - 1;
                deltree(j);
        ny = k;

        #+--------------------------------------------------+
        #| Scatter and permute.                             |

        for k in range(0, ny):
            i = perm[jr][iy[k]];
            y[i] = sy[k];
            tag[i] = currtag;
            addtree(i);

    #+------------------------------------------------------+
    #|               -T                                     |
    #|       y  <-  L  y                                    |

    i=getlast()
    while (i >= rank and i != -1):
        if ( ABS( y[i] ) > eps ):
            consistent = FALSE;
        y[i] = 0.0;
        i=getprev()

    while(i>=0):
        beta = y[i];
        for k in range(0, degLt[i]):
            row = Lt[i][k].i;
            if (tag[row] != currtag):
                y[row] = 0.0;
                tag[row] = currtag;
                addtree(row);
            y[row] -= Lt[i][k].d * beta;
        i = getprev()

    ny = 0;
    i=getfirst();
    while (i != -1):
        if ABS(y[i]) > EPS:
            sy[ny] = y[i]
            iy[ny] = rowperm[i]
            ny = ny + 1

        i = getnext()
    currtag = currtag + 1;
    killtree();

    endtime = clock();
    cumtime = cumtime + (endtime - starttime);

    return ny

def ratnum(m, degB, B, degBt, Bt):
    iwork = []

    for i in range(0, m):
        degBt[i] = 0;

    for j in range(0, m):
        for k in range(0, degB[j]):
            i = B[j][k].i;
            degBt[i] = degBt[i] + 1;

    iwork = MALLOC_INT(m)
    for i in range(0, m):
        REALLOC_VALIND(Bt[i], degBt[i]);

    for j in range(0, m):
        for k in range(0, degB[j]):
            i = B[j][k].i;
            Bt[i][iwork[i]].i = j;
            Bt[i][iwork[i]].d = B[j][k].d;
            iwork[i] = iwork[i] + 1;

    iwork = []


# /*-----------------------------------------------------------------+
# | Forward/backward solve using LU factorization                    |
# | Input:                                                           |
# |    m          dimension of array y                               |
# |    y          array containing right-hand side                   |
# |                                                                  |
# |    static global variables (assumed setup by lufac()):           |
# |                                                                  |
# |    rank       rank of B                                          |
# |    L, degL    ragged array representation of L                   |
# |    Ut, degUt  ragged array representation of U transpose         |
# |               without its diagonal                               |
# |    diag       diagonal entries of U                              |
# |    colperm, icolperm, rowperm, irowperm                          |
# |               column and row permutations and their inverses     |
# | Output:                                                          |
# |                                          -1                      |
# |    y          array containing solution B  y                     |
# |                                                                  |
# |    integer flag indicating whether system is consistent         */

#int     dbsolve(int m, double *y)
def dbsolve(m, y):
    global newcol
    global inewcol
    global cumtime
#{
        # int i;
        # int k, row, consistent=TRUE;
        # double beta, *dwork;
        # double eps;
    consistent=TRUE;

	#double starttime, endtime;

    starttime = clock() #(double) clock();

    dwork = MALLOC_DOUBLE(m) #MALLOC( dwork, m, double );

    eps = 0.0; #TODO
    if (rank < m):
        eps = EPSSOL * maxv(y,m);
    for i in range(0,m): #for (i=0; i<m; i++)
        dwork[i] = y[i];
    for i in range(0, m): #for (i=0; i<m; i++)
        y[irowperm[i]] = dwork[i];

        #/*------------------------------------------------------+
        #|               -1                                      |
        #|       y  <-  L  y                                    */

    for i in range(0, rank): #for (i=0; i<rank; i++) {
        beta = y[i];
        for k in range(0, degL[i]): #for (k=0; k<degL[i]; k++) {
            row = L[i][k].i;
            y[row] = y[row] - L[i][k].d * beta;
        #}
    #}

    #    /*------------------------------------------------------+
	#|                                       -1              |
    #    | Set aside sparse intermediate vector L  P a  for      |
	#|                                            j          |
    #    | refactorization routine.                             */

    #if ( newcol == NULL) MALLOC(  newcol, m, double );
    if (newcol == None or len(newcol) == 0):
        newcol = MALLOC_DOUBLE(m);
	
	#if (inewcol == NULL) MALLOC( inewcol, m, int );
    if (inewcol == None or len(inewcol) == 0):
        inewcol = MALLOC_INT(m);

    nnewcol = 0;
    for i in range(0, m): #for (i=0; i<m; i++) {
        if ( ABS(y[i]) > EPS ): #{
            newcol [nnewcol] = y[i];
            inewcol[nnewcol] = i;
            nnewcol = nnewcol + 1;
	    #}
	#}

        #/*------------------------------------------------------+
        #|               -1                                      |
        #|       y  <-  U  y                                    */

    for i in range(m-1, rank-1, -1): #for (i=m-1; i>=rank; i--) {
        if ( ABS( y[i] ) > eps ):
            consistent = FALSE;
        y[i] = 0.0;
    #}
    for i in range(rank-1, -1, -1): #for (i=rank-1; i>=0; i--) {
        beta = y[i]/diag[i];
        for k in range(0, degU[i]): #for (k=0; k<degU[i]; k++) {
            row = U[i][k].i;
            y[row] = y[row] - U[i][k].d * beta;
		#}
        y[i] = beta;
    #}

    for i in range(0, m): #for (i=0; i<m; i++)
        dwork[i] = y[i];
    for i in range(0, m): #for (i=0; i<m; i++)
        y[colperm[i]] = dwork[i];

    dwork = None #FREE( dwork );

    endtime = clock(); #(double) clock();
    cumtime = cumtime + (endtime - starttime);

    return y; #consistent;
#}

