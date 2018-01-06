#line 1 "sgelsd.f"
/* sgelsd.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "sgelsd.f"
/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b81 = 0.;

/* > \brief <b> SGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b
> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGELSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, */
/*                          RANK, WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELSD computes the minimum-norm solution to a real linear least */
/* > squares problem: */
/* >     minimize 2-norm(| b - A*x |) */
/* > using the singular value decomposition (SVD) of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > */
/* > The problem is solved in three steps: */
/* > (1) Reduce the coefficient matrix A to bidiagonal form with */
/* >     Householder transformations, reducing the original problem */
/* >     into a "bidiagonal least squares problem" (BLS) */
/* > (2) Solve the BLS using a divide and conquer approach. */
/* > (3) Apply back all the Householder tranformations to solve */
/* >     the original least squares problem. */
/* > */
/* > The effective rank of A is determined by treating as zero those */
/* > singular values which are less than RCOND times the largest singular */
/* > value. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, B is overwritten by the N-by-NRHS solution */
/* >          matrix X.  If m >= n and RANK = n, the residual */
/* >          sum-of-squares for the solution in the i-th column is given */
/* >          by the sum of squares of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,max(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (min(M,N)) */
/* >          The singular values of A in decreasing order. */
/* >          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          RCOND is used to determine the effective rank of A. */
/* >          Singular values S(i) <= RCOND*S(1) are treated as zero. */
/* >          If RCOND < 0, machine precision is used instead. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >          The effective rank of A, i.e., the number of singular values */
/* >          which are greater than RCOND*S(1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK must be at least 1. */
/* >          The exact minimum amount of workspace needed depends on M, */
/* >          N and NRHS. As long as LWORK is at least */
/* >              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2, */
/* >          if M is greater than or equal to N or */
/* >              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2, */
/* >          if M is less than N, the code will execute correctly. */
/* >          SMLSIZ is returned by ILAENV and is equal to the maximum */
/* >          size of the subproblems at the bottom of the computation */
/* >          tree (usually about 25), and */
/* >             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the array WORK and the */
/* >          minimum size of the array IWORK, and returns these values as */
/* >          the first entries of the WORK and IWORK arrays, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          LIWORK >= max(1, 3*MINMN*NLVL + 11*MINMN), */
/* >          where MINMN = MIN( M,N ). */
/* >          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  the algorithm for computing the SVD failed to converge; */
/* >                if INFO = i, i off-diagonal elements of an intermediate */
/* >                bidiagonal form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realGEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int sgelsd_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *iwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer ie, il, mm;
    static doublereal eps, anrm, bnrm;
    static integer itau, nlvl, iascl, ibscl;
    static doublereal sfmin;
    static integer minmn, maxmn, itaup, itauq, mnthr, nwork;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *), sgebrd_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slalsd_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), slascl_(char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);
    static integer wlalsd;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    slacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), slaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int sormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer liwork, minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int sormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical lquery;
    static integer smlsiz;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 256 "sgelsd.f"
    /* Parameter adjustments */
#line 256 "sgelsd.f"
    a_dim1 = *lda;
#line 256 "sgelsd.f"
    a_offset = 1 + a_dim1;
#line 256 "sgelsd.f"
    a -= a_offset;
#line 256 "sgelsd.f"
    b_dim1 = *ldb;
#line 256 "sgelsd.f"
    b_offset = 1 + b_dim1;
#line 256 "sgelsd.f"
    b -= b_offset;
#line 256 "sgelsd.f"
    --s;
#line 256 "sgelsd.f"
    --work;
#line 256 "sgelsd.f"
    --iwork;
#line 256 "sgelsd.f"

#line 256 "sgelsd.f"
    /* Function Body */
#line 256 "sgelsd.f"
    *info = 0;
#line 257 "sgelsd.f"
    minmn = min(*m,*n);
#line 258 "sgelsd.f"
    maxmn = max(*m,*n);
#line 259 "sgelsd.f"
    lquery = *lwork == -1;
#line 260 "sgelsd.f"
    if (*m < 0) {
#line 261 "sgelsd.f"
	*info = -1;
#line 262 "sgelsd.f"
    } else if (*n < 0) {
#line 263 "sgelsd.f"
	*info = -2;
#line 264 "sgelsd.f"
    } else if (*nrhs < 0) {
#line 265 "sgelsd.f"
	*info = -3;
#line 266 "sgelsd.f"
    } else if (*lda < max(1,*m)) {
#line 267 "sgelsd.f"
	*info = -5;
#line 268 "sgelsd.f"
    } else if (*ldb < max(1,maxmn)) {
#line 269 "sgelsd.f"
	*info = -7;
#line 270 "sgelsd.f"
    }

/*     Compute workspace. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 279 "sgelsd.f"
    if (*info == 0) {
#line 280 "sgelsd.f"
	minwrk = 1;
#line 281 "sgelsd.f"
	maxwrk = 1;
#line 282 "sgelsd.f"
	liwork = 1;
#line 283 "sgelsd.f"
	if (minmn > 0) {
#line 284 "sgelsd.f"
	    smlsiz = ilaenv_(&c__9, "SGELSD", " ", &c__0, &c__0, &c__0, &c__0,
		     (ftnlen)6, (ftnlen)1);
#line 285 "sgelsd.f"
	    mnthr = ilaenv_(&c__6, "SGELSD", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
/* Computing MAX */
#line 286 "sgelsd.f"
	    i__1 = (integer) (log((doublereal) minmn / (doublereal) (smlsiz + 
		    1)) / log(2.)) + 1;
#line 286 "sgelsd.f"
	    nlvl = max(i__1,0);
#line 288 "sgelsd.f"
	    liwork = minmn * 3 * nlvl + minmn * 11;
#line 289 "sgelsd.f"
	    mm = *m;
#line 290 "sgelsd.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns. */

#line 295 "sgelsd.f"
		mm = *n;
/* Computing MAX */
#line 296 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "SGEQRF", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 296 "sgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 298 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "SORMQR", 
			"LT", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 298 "sgelsd.f"
		maxwrk = max(i__1,i__2);
#line 300 "sgelsd.f"
	    }
#line 301 "sgelsd.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined. */

/* Computing MAX */
#line 305 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, 
			"SGEBRD", " ", &mm, n, &c_n1, &c_n1, (ftnlen)6, (
			ftnlen)1);
#line 305 "sgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 307 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, "SORMBR"
			, "QLT", &mm, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 307 "sgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 309 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, 
			"SORMBR", "PLN", n, nrhs, n, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 309 "sgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing 2nd power */
#line 311 "sgelsd.f"
		i__1 = smlsiz + 1;
#line 311 "sgelsd.f"
		wlalsd = *n * 9 + (*n << 1) * smlsiz + (*n << 3) * nlvl + *n *
			 *nrhs + i__1 * i__1;
/* Computing MAX */
#line 313 "sgelsd.f"
		i__1 = maxwrk, i__2 = *n * 3 + wlalsd;
#line 313 "sgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 314 "sgelsd.f"
		i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,
			i__2), i__2 = *n * 3 + wlalsd;
#line 314 "sgelsd.f"
		minwrk = max(i__1,i__2);
#line 315 "sgelsd.f"
	    }
#line 316 "sgelsd.f"
	    if (*n > *m) {
/* Computing 2nd power */
#line 317 "sgelsd.f"
		i__1 = smlsiz + 1;
#line 317 "sgelsd.f"
		wlalsd = *m * 9 + (*m << 1) * smlsiz + (*m << 3) * nlvl + *m *
			 *nrhs + i__1 * i__1;
#line 319 "sgelsd.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                           than rows. */

#line 324 "sgelsd.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "SGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 326 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			    ilaenv_(&c__1, "SGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 326 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 328 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * 
			    ilaenv_(&c__1, "SORMBR", "QLT", m, nrhs, m, &c_n1,
			     (ftnlen)6, (ftnlen)3);
#line 328 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 330 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			    ilaenv_(&c__1, "SORMBR", "PLN", m, nrhs, m, &c_n1,
			     (ftnlen)6, (ftnlen)3);
#line 330 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 332 "sgelsd.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 333 "sgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 333 "sgelsd.f"
			maxwrk = max(i__1,i__2);
#line 334 "sgelsd.f"
		    } else {
/* Computing MAX */
#line 335 "sgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 335 "sgelsd.f"
			maxwrk = max(i__1,i__2);
#line 336 "sgelsd.f"
		    }
/* Computing MAX */
#line 337 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m + *nrhs * ilaenv_(&c__1, "SORMLQ"
			    , "LT", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)2);
#line 337 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 339 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + wlalsd;
#line 339 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/*     XXX: Ensure the Path 2a case below is triggered.  The workspace */
/*     calculation should use queries for all routines eventually. */
/* Computing MAX */
/* Computing MAX */
#line 342 "sgelsd.f"
		    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), 
			    i__3 = max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 342 "sgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 2) + *m * *m + max(i__3,i__4)
			    ;
#line 342 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 344 "sgelsd.f"
		} else {

/*                 Path 2 - remaining underdetermined cases. */

#line 348 "sgelsd.f"
		    maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, "SGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 350 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, 
			    "SORMBR", "QLT", m, nrhs, n, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 350 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 352 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * ilaenv_(&c__1, "SORM"\
			    "BR", "PLN", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)
			    3);
#line 352 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 354 "sgelsd.f"
		    i__1 = maxwrk, i__2 = *m * 3 + wlalsd;
#line 354 "sgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 355 "sgelsd.f"
		}
/* Computing MAX */
#line 356 "sgelsd.f"
		i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *m, i__1 = max(i__1,
			i__2), i__2 = *m * 3 + wlalsd;
#line 356 "sgelsd.f"
		minwrk = max(i__1,i__2);
#line 357 "sgelsd.f"
	    }
#line 358 "sgelsd.f"
	}
#line 359 "sgelsd.f"
	minwrk = min(minwrk,maxwrk);
#line 360 "sgelsd.f"
	work[1] = (doublereal) maxwrk;
#line 361 "sgelsd.f"
	iwork[1] = liwork;

#line 363 "sgelsd.f"
	if (*lwork < minwrk && ! lquery) {
#line 364 "sgelsd.f"
	    *info = -12;
#line 365 "sgelsd.f"
	}
#line 366 "sgelsd.f"
    }

#line 368 "sgelsd.f"
    if (*info != 0) {
#line 369 "sgelsd.f"
	i__1 = -(*info);
#line 369 "sgelsd.f"
	xerbla_("SGELSD", &i__1, (ftnlen)6);
#line 370 "sgelsd.f"
	return 0;
#line 371 "sgelsd.f"
    } else if (lquery) {
#line 372 "sgelsd.f"
	return 0;
#line 373 "sgelsd.f"
    }

/*     Quick return if possible. */

#line 377 "sgelsd.f"
    if (*m == 0 || *n == 0) {
#line 378 "sgelsd.f"
	*rank = 0;
#line 379 "sgelsd.f"
	return 0;
#line 380 "sgelsd.f"
    }

/*     Get machine parameters. */

#line 384 "sgelsd.f"
    eps = slamch_("P", (ftnlen)1);
#line 385 "sgelsd.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 386 "sgelsd.f"
    smlnum = sfmin / eps;
#line 387 "sgelsd.f"
    bignum = 1. / smlnum;
#line 388 "sgelsd.f"
    slabad_(&smlnum, &bignum);

/*     Scale A if max entry outside range [SMLNUM,BIGNUM]. */

#line 392 "sgelsd.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 393 "sgelsd.f"
    iascl = 0;
#line 394 "sgelsd.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 398 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 399 "sgelsd.f"
	iascl = 1;
#line 400 "sgelsd.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 404 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 405 "sgelsd.f"
	iascl = 2;
#line 406 "sgelsd.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 410 "sgelsd.f"
	i__1 = max(*m,*n);
#line 410 "sgelsd.f"
	slaset_("F", &i__1, nrhs, &c_b81, &c_b81, &b[b_offset], ldb, (ftnlen)
		1);
#line 411 "sgelsd.f"
	slaset_("F", &minmn, &c__1, &c_b81, &c_b81, &s[1], &c__1, (ftnlen)1);
#line 412 "sgelsd.f"
	*rank = 0;
#line 413 "sgelsd.f"
	goto L10;
#line 414 "sgelsd.f"
    }

/*     Scale B if max entry outside range [SMLNUM,BIGNUM]. */

#line 418 "sgelsd.f"
    bnrm = slange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 419 "sgelsd.f"
    ibscl = 0;
#line 420 "sgelsd.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 424 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 425 "sgelsd.f"
	ibscl = 1;
#line 426 "sgelsd.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 430 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 431 "sgelsd.f"
	ibscl = 2;
#line 432 "sgelsd.f"
    }

/*     If M < N make sure certain entries of B are zero. */

#line 436 "sgelsd.f"
    if (*m < *n) {
#line 436 "sgelsd.f"
	i__1 = *n - *m;
#line 436 "sgelsd.f"
	slaset_("F", &i__1, nrhs, &c_b81, &c_b81, &b[*m + 1 + b_dim1], ldb, (
		ftnlen)1);
#line 436 "sgelsd.f"
    }

/*     Overdetermined case. */

#line 441 "sgelsd.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined. */

#line 445 "sgelsd.f"
	mm = *m;
#line 446 "sgelsd.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns. */

#line 450 "sgelsd.f"
	    mm = *n;
#line 451 "sgelsd.f"
	    itau = 1;
#line 452 "sgelsd.f"
	    nwork = itau + *n;

/*           Compute A=Q*R. */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 457 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 457 "sgelsd.f"
	    sgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);

/*           Multiply B by transpose(Q). */
/*           (Workspace: need N+NRHS, prefer N+NRHS*NB) */

#line 463 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 463 "sgelsd.f"
	    sormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R. */

#line 468 "sgelsd.f"
	    if (*n > 1) {
#line 469 "sgelsd.f"
		i__1 = *n - 1;
#line 469 "sgelsd.f"
		i__2 = *n - 1;
#line 469 "sgelsd.f"
		slaset_("L", &i__1, &i__2, &c_b81, &c_b81, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 470 "sgelsd.f"
	    }
#line 471 "sgelsd.f"
	}

#line 473 "sgelsd.f"
	ie = 1;
#line 474 "sgelsd.f"
	itauq = ie + *n;
#line 475 "sgelsd.f"
	itaup = itauq + *n;
#line 476 "sgelsd.f"
	nwork = itaup + *n;

/*        Bidiagonalize R in A. */
/*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

#line 481 "sgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 481 "sgelsd.f"
	sgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R. */
/*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

#line 488 "sgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 488 "sgelsd.f"
	sormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 493 "sgelsd.f"
	slalsd_("U", &smlsiz, n, nrhs, &s[1], &work[ie], &b[b_offset], ldb, 
		rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)1);
#line 495 "sgelsd.f"
	if (*info != 0) {
#line 496 "sgelsd.f"
	    goto L10;
#line 497 "sgelsd.f"
	}

/*        Multiply B by right bidiagonalizing vectors of R. */

#line 501 "sgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 501 "sgelsd.f"
	sormbr_("P", "L", "N", n, nrhs, n, &a[a_offset], lda, &work[itaup], &
		b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

#line 504 "sgelsd.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 504 "sgelsd.f"
	i__1 = *m, i__2 = (*m << 1) - 4, i__1 = max(i__1,i__2), i__1 = max(
		i__1,*nrhs), i__2 = *n - *m * 3, i__1 = max(i__1,i__2);
#line 504 "sgelsd.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__1,wlalsd)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm. */

#line 510 "sgelsd.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 511 "sgelsd.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 511 "sgelsd.f"
	    i__1 = (*m << 2) + *m * *lda + max(i__3,i__4), i__2 = *m * *lda + 
		    *m + *m * *nrhs, i__1 = max(i__1,i__2), i__2 = (*m << 2) 
		    + *m * *lda + wlalsd;
#line 511 "sgelsd.f"
	    if (*lwork >= max(i__1,i__2)) {
#line 511 "sgelsd.f"
		ldwork = *lda;
#line 511 "sgelsd.f"
	    }
#line 513 "sgelsd.f"
	    itau = 1;
#line 514 "sgelsd.f"
	    nwork = *m + 1;

/*        Compute A=L*Q. */
/*        (Workspace: need 2*M, prefer M+M*NB) */

#line 519 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 519 "sgelsd.f"
	    sgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);
#line 521 "sgelsd.f"
	    il = nwork;

/*        Copy L to WORK(IL), zeroing out above its diagonal. */

#line 525 "sgelsd.f"
	    slacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 526 "sgelsd.f"
	    i__1 = *m - 1;
#line 526 "sgelsd.f"
	    i__2 = *m - 1;
#line 526 "sgelsd.f"
	    slaset_("U", &i__1, &i__2, &c_b81, &c_b81, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 528 "sgelsd.f"
	    ie = il + ldwork * *m;
#line 529 "sgelsd.f"
	    itauq = ie + *m;
#line 530 "sgelsd.f"
	    itaup = itauq + *m;
#line 531 "sgelsd.f"
	    nwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL). */
/*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

#line 536 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 536 "sgelsd.f"
	    sgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L. */
/*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 543 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 543 "sgelsd.f"
	    sormbr_("Q", "L", "T", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 549 "sgelsd.f"
	    slalsd_("U", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)
		    1);
#line 551 "sgelsd.f"
	    if (*info != 0) {
#line 552 "sgelsd.f"
		goto L10;
#line 553 "sgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of L. */

#line 557 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 557 "sgelsd.f"
	    sormbr_("P", "L", "N", m, nrhs, m, &work[il], &ldwork, &work[
		    itaup], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Zero out below first M rows of B. */

#line 563 "sgelsd.f"
	    i__1 = *n - *m;
#line 563 "sgelsd.f"
	    slaset_("F", &i__1, nrhs, &c_b81, &c_b81, &b[*m + 1 + b_dim1], 
		    ldb, (ftnlen)1);
#line 564 "sgelsd.f"
	    nwork = itau + *m;

/*        Multiply transpose(Q) by B. */
/*        (Workspace: need M+NRHS, prefer M+NRHS*NB) */

#line 569 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 569 "sgelsd.f"
	    sormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 572 "sgelsd.f"
	} else {

/*        Path 2 - remaining underdetermined cases. */

#line 576 "sgelsd.f"
	    ie = 1;
#line 577 "sgelsd.f"
	    itauq = ie + *m;
#line 578 "sgelsd.f"
	    itaup = itauq + *m;
#line 579 "sgelsd.f"
	    nwork = itaup + *m;

/*        Bidiagonalize A. */
/*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 584 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 584 "sgelsd.f"
	    sgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors. */
/*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

#line 591 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 591 "sgelsd.f"
	    sormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 596 "sgelsd.f"
	    slalsd_("L", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)
		    1);
#line 598 "sgelsd.f"
	    if (*info != 0) {
#line 599 "sgelsd.f"
		goto L10;
#line 600 "sgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of A. */

#line 604 "sgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 604 "sgelsd.f"
	    sormbr_("P", "L", "N", n, nrhs, m, &a[a_offset], lda, &work[itaup]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

#line 607 "sgelsd.f"
	}
#line 607 "sgelsd.f"
    }

/*     Undo scaling. */

#line 611 "sgelsd.f"
    if (iascl == 1) {
#line 612 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 613 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 615 "sgelsd.f"
    } else if (iascl == 2) {
#line 616 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 617 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 619 "sgelsd.f"
    }
#line 620 "sgelsd.f"
    if (ibscl == 1) {
#line 621 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 622 "sgelsd.f"
    } else if (ibscl == 2) {
#line 623 "sgelsd.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 624 "sgelsd.f"
    }

#line 626 "sgelsd.f"
L10:
#line 627 "sgelsd.f"
    work[1] = (doublereal) maxwrk;
#line 628 "sgelsd.f"
    iwork[1] = liwork;
#line 629 "sgelsd.f"
    return 0;

/*     End of SGELSD */

} /* sgelsd_ */

