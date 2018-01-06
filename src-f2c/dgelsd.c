#line 1 "dgelsd.f"
/* dgelsd.f -- translated by f2c (version 20100827).
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

#line 1 "dgelsd.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b82 = 0.;

/* > \brief <b> DGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b
> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGELSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELSD computes the minimum-norm solution to a real linear least */
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
/* > (3) Apply back all the Householder transformations to solve */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A in decreasing order. */
/* >          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN), */
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

/* > \date December 2016 */

/* > \ingroup doubleGEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int dgelsd_(integer *m, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebrd_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlalsd_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), dlascl_(char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), dgeqrf_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer wlalsd;
    extern /* Subroutine */ int dormlq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer liwork, minwrk, maxwrk;
    static doublereal smlnum;
    static logical lquery;
    static integer smlsiz;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 255 "dgelsd.f"
    /* Parameter adjustments */
#line 255 "dgelsd.f"
    a_dim1 = *lda;
#line 255 "dgelsd.f"
    a_offset = 1 + a_dim1;
#line 255 "dgelsd.f"
    a -= a_offset;
#line 255 "dgelsd.f"
    b_dim1 = *ldb;
#line 255 "dgelsd.f"
    b_offset = 1 + b_dim1;
#line 255 "dgelsd.f"
    b -= b_offset;
#line 255 "dgelsd.f"
    --s;
#line 255 "dgelsd.f"
    --work;
#line 255 "dgelsd.f"
    --iwork;
#line 255 "dgelsd.f"

#line 255 "dgelsd.f"
    /* Function Body */
#line 255 "dgelsd.f"
    *info = 0;
#line 256 "dgelsd.f"
    minmn = min(*m,*n);
#line 257 "dgelsd.f"
    maxmn = max(*m,*n);
#line 258 "dgelsd.f"
    mnthr = ilaenv_(&c__6, "DGELSD", " ", m, n, nrhs, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 259 "dgelsd.f"
    lquery = *lwork == -1;
#line 260 "dgelsd.f"
    if (*m < 0) {
#line 261 "dgelsd.f"
	*info = -1;
#line 262 "dgelsd.f"
    } else if (*n < 0) {
#line 263 "dgelsd.f"
	*info = -2;
#line 264 "dgelsd.f"
    } else if (*nrhs < 0) {
#line 265 "dgelsd.f"
	*info = -3;
#line 266 "dgelsd.f"
    } else if (*lda < max(1,*m)) {
#line 267 "dgelsd.f"
	*info = -5;
#line 268 "dgelsd.f"
    } else if (*ldb < max(1,maxmn)) {
#line 269 "dgelsd.f"
	*info = -7;
#line 270 "dgelsd.f"
    }

#line 272 "dgelsd.f"
    smlsiz = ilaenv_(&c__9, "DGELSD", " ", &c__0, &c__0, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)1);

/*     Compute workspace. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 281 "dgelsd.f"
    minwrk = 1;
#line 282 "dgelsd.f"
    liwork = 1;
#line 283 "dgelsd.f"
    minmn = max(1,minmn);
/* Computing MAX */
#line 284 "dgelsd.f"
    i__1 = (integer) (log((doublereal) minmn / (doublereal) (smlsiz + 1)) / 
	    log(2.)) + 1;
#line 284 "dgelsd.f"
    nlvl = max(i__1,0);

#line 287 "dgelsd.f"
    if (*info == 0) {
#line 288 "dgelsd.f"
	maxwrk = 0;
#line 289 "dgelsd.f"
	liwork = minmn * 3 * nlvl + minmn * 11;
#line 290 "dgelsd.f"
	mm = *m;
#line 291 "dgelsd.f"
	if (*m >= *n && *m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns. */

#line 295 "dgelsd.f"
	    mm = *n;
/* Computing MAX */
#line 296 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, 
		    n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 296 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 298 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "DORMQR", "LT", 
		    m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 298 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
#line 300 "dgelsd.f"
	}
#line 301 "dgelsd.f"
	if (*m >= *n) {

/*           Path 1 - overdetermined or exactly determined. */

/* Computing MAX */
#line 305 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, "DGEBRD"
		    , " ", &mm, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 305 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 307 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, "DORMBR", 
		    "QLT", &mm, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 307 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 309 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, "DORMBR",
		     "PLN", n, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 309 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
/* Computing 2nd power */
#line 311 "dgelsd.f"
	    i__1 = smlsiz + 1;
#line 311 "dgelsd.f"
	    wlalsd = *n * 9 + (*n << 1) * smlsiz + (*n << 3) * nlvl + *n * *
		    nrhs + i__1 * i__1;
/* Computing MAX */
#line 312 "dgelsd.f"
	    i__1 = maxwrk, i__2 = *n * 3 + wlalsd;
#line 312 "dgelsd.f"
	    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 313 "dgelsd.f"
	    i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,i__2), 
		    i__2 = *n * 3 + wlalsd;
#line 313 "dgelsd.f"
	    minwrk = max(i__1,i__2);
#line 314 "dgelsd.f"
	}
#line 315 "dgelsd.f"
	if (*n > *m) {
/* Computing 2nd power */
#line 316 "dgelsd.f"
	    i__1 = smlsiz + 1;
#line 316 "dgelsd.f"
	    wlalsd = *m * 9 + (*m << 1) * smlsiz + (*m << 3) * nlvl + *m * *
		    nrhs + i__1 * i__1;
#line 317 "dgelsd.f"
	    if (*n >= mnthr) {

/*              Path 2a - underdetermined, with many more columns */
/*              than rows. */

#line 322 "dgelsd.f"
		maxwrk = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, 
			&c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 323 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1, (
			ftnlen)6, (ftnlen)1);
#line 323 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 325 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * ilaenv_(&
			c__1, "DORMBR", "QLT", m, nrhs, m, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 325 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 327 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			ilaenv_(&c__1, "DORMBR", "PLN", m, nrhs, m, &c_n1, (
			ftnlen)6, (ftnlen)3);
#line 327 "dgelsd.f"
		maxwrk = max(i__1,i__2);
#line 329 "dgelsd.f"
		if (*nrhs > 1) {
/* Computing MAX */
#line 330 "dgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 330 "dgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 331 "dgelsd.f"
		} else {
/* Computing MAX */
#line 332 "dgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 332 "dgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 333 "dgelsd.f"
		}
/* Computing MAX */
#line 334 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m + *nrhs * ilaenv_(&c__1, "DORMLQ", 
			"LT", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)2);
#line 334 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 336 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + wlalsd;
#line 336 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/*     XXX: Ensure the Path 2a case below is triggered.  The workspace */
/*     calculation should use queries for all routines eventually. */
/* Computing MAX */
/* Computing MAX */
#line 339 "dgelsd.f"
		i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 =
			 max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 339 "dgelsd.f"
		i__1 = maxwrk, i__2 = (*m << 2) + *m * *m + max(i__3,i__4);
#line 339 "dgelsd.f"
		maxwrk = max(i__1,i__2);
#line 341 "dgelsd.f"
	    } else {

/*              Path 2 - remaining underdetermined cases. */

#line 345 "dgelsd.f"
		maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, "DGEBRD", " ", m,
			 n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 347 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, "DORMBR"
			, "QLT", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)3);
#line 347 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 349 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORMBR", 
			"PLN", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)3);
#line 349 "dgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 351 "dgelsd.f"
		i__1 = maxwrk, i__2 = *m * 3 + wlalsd;
#line 351 "dgelsd.f"
		maxwrk = max(i__1,i__2);
#line 352 "dgelsd.f"
	    }
/* Computing MAX */
#line 353 "dgelsd.f"
	    i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *m, i__1 = max(i__1,i__2), 
		    i__2 = *m * 3 + wlalsd;
#line 353 "dgelsd.f"
	    minwrk = max(i__1,i__2);
#line 354 "dgelsd.f"
	}
#line 355 "dgelsd.f"
	minwrk = min(minwrk,maxwrk);
#line 356 "dgelsd.f"
	work[1] = (doublereal) maxwrk;
#line 357 "dgelsd.f"
	iwork[1] = liwork;
#line 359 "dgelsd.f"
	if (*lwork < minwrk && ! lquery) {
#line 360 "dgelsd.f"
	    *info = -12;
#line 361 "dgelsd.f"
	}
#line 362 "dgelsd.f"
    }

#line 364 "dgelsd.f"
    if (*info != 0) {
#line 365 "dgelsd.f"
	i__1 = -(*info);
#line 365 "dgelsd.f"
	xerbla_("DGELSD", &i__1, (ftnlen)6);
#line 366 "dgelsd.f"
	return 0;
#line 367 "dgelsd.f"
    } else if (lquery) {
#line 368 "dgelsd.f"
	goto L10;
#line 369 "dgelsd.f"
    }

/*     Quick return if possible. */

#line 373 "dgelsd.f"
    if (*m == 0 || *n == 0) {
#line 374 "dgelsd.f"
	*rank = 0;
#line 375 "dgelsd.f"
	return 0;
#line 376 "dgelsd.f"
    }

/*     Get machine parameters. */

#line 380 "dgelsd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 381 "dgelsd.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 382 "dgelsd.f"
    smlnum = sfmin / eps;
#line 383 "dgelsd.f"
    bignum = 1. / smlnum;
#line 384 "dgelsd.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A if max entry outside range [SMLNUM,BIGNUM]. */

#line 388 "dgelsd.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 389 "dgelsd.f"
    iascl = 0;
#line 390 "dgelsd.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 394 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 395 "dgelsd.f"
	iascl = 1;
#line 396 "dgelsd.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 400 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 401 "dgelsd.f"
	iascl = 2;
#line 402 "dgelsd.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 406 "dgelsd.f"
	i__1 = max(*m,*n);
#line 406 "dgelsd.f"
	dlaset_("F", &i__1, nrhs, &c_b82, &c_b82, &b[b_offset], ldb, (ftnlen)
		1);
#line 407 "dgelsd.f"
	dlaset_("F", &minmn, &c__1, &c_b82, &c_b82, &s[1], &c__1, (ftnlen)1);
#line 408 "dgelsd.f"
	*rank = 0;
#line 409 "dgelsd.f"
	goto L10;
#line 410 "dgelsd.f"
    }

/*     Scale B if max entry outside range [SMLNUM,BIGNUM]. */

#line 414 "dgelsd.f"
    bnrm = dlange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 415 "dgelsd.f"
    ibscl = 0;
#line 416 "dgelsd.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 420 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 421 "dgelsd.f"
	ibscl = 1;
#line 422 "dgelsd.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 426 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 427 "dgelsd.f"
	ibscl = 2;
#line 428 "dgelsd.f"
    }

/*     If M < N make sure certain entries of B are zero. */

#line 432 "dgelsd.f"
    if (*m < *n) {
#line 432 "dgelsd.f"
	i__1 = *n - *m;
#line 432 "dgelsd.f"
	dlaset_("F", &i__1, nrhs, &c_b82, &c_b82, &b[*m + 1 + b_dim1], ldb, (
		ftnlen)1);
#line 432 "dgelsd.f"
    }

/*     Overdetermined case. */

#line 437 "dgelsd.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined. */

#line 441 "dgelsd.f"
	mm = *m;
#line 442 "dgelsd.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns. */

#line 446 "dgelsd.f"
	    mm = *n;
#line 447 "dgelsd.f"
	    itau = 1;
#line 448 "dgelsd.f"
	    nwork = itau + *n;

/*           Compute A=Q*R. */
/*           (Workspace: need 2*N, prefer N+N*NB) */

#line 453 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 453 "dgelsd.f"
	    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);

/*           Multiply B by transpose(Q). */
/*           (Workspace: need N+NRHS, prefer N+NRHS*NB) */

#line 459 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 459 "dgelsd.f"
	    dormqr_("L", "T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R. */

#line 464 "dgelsd.f"
	    if (*n > 1) {
#line 465 "dgelsd.f"
		i__1 = *n - 1;
#line 465 "dgelsd.f"
		i__2 = *n - 1;
#line 465 "dgelsd.f"
		dlaset_("L", &i__1, &i__2, &c_b82, &c_b82, &a[a_dim1 + 2], 
			lda, (ftnlen)1);
#line 466 "dgelsd.f"
	    }
#line 467 "dgelsd.f"
	}

#line 469 "dgelsd.f"
	ie = 1;
#line 470 "dgelsd.f"
	itauq = ie + *n;
#line 471 "dgelsd.f"
	itaup = itauq + *n;
#line 472 "dgelsd.f"
	nwork = itaup + *n;

/*        Bidiagonalize R in A. */
/*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

#line 477 "dgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 477 "dgelsd.f"
	dgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R. */
/*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

#line 484 "dgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 484 "dgelsd.f"
	dormbr_("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 489 "dgelsd.f"
	dlalsd_("U", &smlsiz, n, nrhs, &s[1], &work[ie], &b[b_offset], ldb, 
		rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)1);
#line 491 "dgelsd.f"
	if (*info != 0) {
#line 492 "dgelsd.f"
	    goto L10;
#line 493 "dgelsd.f"
	}

/*        Multiply B by right bidiagonalizing vectors of R. */

#line 497 "dgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 497 "dgelsd.f"
	dormbr_("P", "L", "N", n, nrhs, n, &a[a_offset], lda, &work[itaup], &
		b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

#line 500 "dgelsd.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 500 "dgelsd.f"
	i__1 = *m, i__2 = (*m << 1) - 4, i__1 = max(i__1,i__2), i__1 = max(
		i__1,*nrhs), i__2 = *n - *m * 3, i__1 = max(i__1,i__2);
#line 500 "dgelsd.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__1,wlalsd)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm. */

#line 506 "dgelsd.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 507 "dgelsd.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 507 "dgelsd.f"
	    i__1 = (*m << 2) + *m * *lda + max(i__3,i__4), i__2 = *m * *lda + 
		    *m + *m * *nrhs, i__1 = max(i__1,i__2), i__2 = (*m << 2) 
		    + *m * *lda + wlalsd;
#line 507 "dgelsd.f"
	    if (*lwork >= max(i__1,i__2)) {
#line 507 "dgelsd.f"
		ldwork = *lda;
#line 507 "dgelsd.f"
	    }
#line 509 "dgelsd.f"
	    itau = 1;
#line 510 "dgelsd.f"
	    nwork = *m + 1;

/*        Compute A=L*Q. */
/*        (Workspace: need 2*M, prefer M+M*NB) */

#line 515 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 515 "dgelsd.f"
	    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);
#line 517 "dgelsd.f"
	    il = nwork;

/*        Copy L to WORK(IL), zeroing out above its diagonal. */

#line 521 "dgelsd.f"
	    dlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 522 "dgelsd.f"
	    i__1 = *m - 1;
#line 522 "dgelsd.f"
	    i__2 = *m - 1;
#line 522 "dgelsd.f"
	    dlaset_("U", &i__1, &i__2, &c_b82, &c_b82, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 524 "dgelsd.f"
	    ie = il + ldwork * *m;
#line 525 "dgelsd.f"
	    itauq = ie + *m;
#line 526 "dgelsd.f"
	    itaup = itauq + *m;
#line 527 "dgelsd.f"
	    nwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL). */
/*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

#line 532 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 532 "dgelsd.f"
	    dgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L. */
/*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 539 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 539 "dgelsd.f"
	    dormbr_("Q", "L", "T", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 545 "dgelsd.f"
	    dlalsd_("U", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)
		    1);
#line 547 "dgelsd.f"
	    if (*info != 0) {
#line 548 "dgelsd.f"
		goto L10;
#line 549 "dgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of L. */

#line 553 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 553 "dgelsd.f"
	    dormbr_("P", "L", "N", m, nrhs, m, &work[il], &ldwork, &work[
		    itaup], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Zero out below first M rows of B. */

#line 559 "dgelsd.f"
	    i__1 = *n - *m;
#line 559 "dgelsd.f"
	    dlaset_("F", &i__1, nrhs, &c_b82, &c_b82, &b[*m + 1 + b_dim1], 
		    ldb, (ftnlen)1);
#line 560 "dgelsd.f"
	    nwork = itau + *m;

/*        Multiply transpose(Q) by B. */
/*        (Workspace: need M+NRHS, prefer M+NRHS*NB) */

#line 565 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 565 "dgelsd.f"
	    dormlq_("L", "T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 568 "dgelsd.f"
	} else {

/*        Path 2 - remaining underdetermined cases. */

#line 572 "dgelsd.f"
	    ie = 1;
#line 573 "dgelsd.f"
	    itauq = ie + *m;
#line 574 "dgelsd.f"
	    itaup = itauq + *m;
#line 575 "dgelsd.f"
	    nwork = itaup + *m;

/*        Bidiagonalize A. */
/*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

#line 580 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 580 "dgelsd.f"
	    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors. */
/*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

#line 587 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 587 "dgelsd.f"
	    dormbr_("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 592 "dgelsd.f"
	    dlalsd_("L", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &iwork[1], info, (ftnlen)
		    1);
#line 594 "dgelsd.f"
	    if (*info != 0) {
#line 595 "dgelsd.f"
		goto L10;
#line 596 "dgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of A. */

#line 600 "dgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 600 "dgelsd.f"
	    dormbr_("P", "L", "N", n, nrhs, m, &a[a_offset], lda, &work[itaup]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

#line 603 "dgelsd.f"
	}
#line 603 "dgelsd.f"
    }

/*     Undo scaling. */

#line 607 "dgelsd.f"
    if (iascl == 1) {
#line 608 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 609 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 611 "dgelsd.f"
    } else if (iascl == 2) {
#line 612 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 613 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 615 "dgelsd.f"
    }
#line 616 "dgelsd.f"
    if (ibscl == 1) {
#line 617 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 618 "dgelsd.f"
    } else if (ibscl == 2) {
#line 619 "dgelsd.f"
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 620 "dgelsd.f"
    }

#line 622 "dgelsd.f"
L10:
#line 623 "dgelsd.f"
    work[1] = (doublereal) maxwrk;
#line 624 "dgelsd.f"
    iwork[1] = liwork;
#line 625 "dgelsd.f"
    return 0;

/*     End of DGELSD */

} /* dgelsd_ */

