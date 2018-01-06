#line 1 "zgelsd.f"
/* zgelsd.f -- translated by f2c (version 20100827).
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

#line 1 "zgelsd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b80 = 0.;

/* > \brief <b> ZGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b
> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGELSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELSD computes the minimum-norm solution to a real linear least */
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
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. N >= 0. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, B is overwritten by the N-by-NRHS solution matrix X. */
/* >          If m >= n and RANK = n, the residual sum-of-squares for */
/* >          the solution in the i-th column is given by the sum of */
/* >          squares of the modulus of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M,N). */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK must be at least 1. */
/* >          The exact minimum amount of workspace needed depends on M, */
/* >          N and NRHS. As long as LWORK is at least */
/* >              2*N + N*NRHS */
/* >          if M is greater than or equal to N or */
/* >              2*M + M*NRHS */
/* >          if M is less than N, the code will execute correctly. */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the array WORK and the */
/* >          minimum sizes of the arrays RWORK and IWORK, and returns */
/* >          these values as the first entries of the WORK, RWORK and */
/* >          IWORK arrays, and no error message related to LWORK is issued */
/* >          by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK)) */
/* >          LRWORK >= */
/* >             10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + */
/* >             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ) */
/* >          if M is greater than or equal to N or */
/* >             10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS + */
/* >             MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ) */
/* >          if M is less than N, the code will execute correctly. */
/* >          SMLSIZ is returned by ILAENV and is equal to the maximum */
/* >          size of the subproblems at the bottom of the computation */
/* >          tree (usually about 25), and */
/* >             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) */
/* >          On exit, if INFO = 0, RWORK(1) returns the minimum LRWORK. */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
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

/* > \ingroup complex16GEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int zgelsd_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info)
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
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zgebrd_(integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlalsd_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     doublecomplex *, doublereal *, integer *, integer *, ftnlen), 
	    zlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *);
    static integer ldwork;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer liwork, minwrk, maxwrk;
    static doublereal smlnum;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen);
    static integer lrwork;
    static logical lquery;
    static integer nrwork, smlsiz;
    extern /* Subroutine */ int zunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

#line 275 "zgelsd.f"
    /* Parameter adjustments */
#line 275 "zgelsd.f"
    a_dim1 = *lda;
#line 275 "zgelsd.f"
    a_offset = 1 + a_dim1;
#line 275 "zgelsd.f"
    a -= a_offset;
#line 275 "zgelsd.f"
    b_dim1 = *ldb;
#line 275 "zgelsd.f"
    b_offset = 1 + b_dim1;
#line 275 "zgelsd.f"
    b -= b_offset;
#line 275 "zgelsd.f"
    --s;
#line 275 "zgelsd.f"
    --work;
#line 275 "zgelsd.f"
    --rwork;
#line 275 "zgelsd.f"
    --iwork;
#line 275 "zgelsd.f"

#line 275 "zgelsd.f"
    /* Function Body */
#line 275 "zgelsd.f"
    *info = 0;
#line 276 "zgelsd.f"
    minmn = min(*m,*n);
#line 277 "zgelsd.f"
    maxmn = max(*m,*n);
#line 278 "zgelsd.f"
    lquery = *lwork == -1;
#line 279 "zgelsd.f"
    if (*m < 0) {
#line 280 "zgelsd.f"
	*info = -1;
#line 281 "zgelsd.f"
    } else if (*n < 0) {
#line 282 "zgelsd.f"
	*info = -2;
#line 283 "zgelsd.f"
    } else if (*nrhs < 0) {
#line 284 "zgelsd.f"
	*info = -3;
#line 285 "zgelsd.f"
    } else if (*lda < max(1,*m)) {
#line 286 "zgelsd.f"
	*info = -5;
#line 287 "zgelsd.f"
    } else if (*ldb < max(1,maxmn)) {
#line 288 "zgelsd.f"
	*info = -7;
#line 289 "zgelsd.f"
    }

/*     Compute workspace. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 298 "zgelsd.f"
    if (*info == 0) {
#line 299 "zgelsd.f"
	minwrk = 1;
#line 300 "zgelsd.f"
	maxwrk = 1;
#line 301 "zgelsd.f"
	liwork = 1;
#line 302 "zgelsd.f"
	lrwork = 1;
#line 303 "zgelsd.f"
	if (minmn > 0) {
#line 304 "zgelsd.f"
	    smlsiz = ilaenv_(&c__9, "ZGELSD", " ", &c__0, &c__0, &c__0, &c__0,
		     (ftnlen)6, (ftnlen)1);
#line 305 "zgelsd.f"
	    mnthr = ilaenv_(&c__6, "ZGELSD", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
/* Computing MAX */
#line 306 "zgelsd.f"
	    i__1 = (integer) (log((doublereal) minmn / (doublereal) (smlsiz + 
		    1)) / log(2.)) + 1;
#line 306 "zgelsd.f"
	    nlvl = max(i__1,0);
#line 308 "zgelsd.f"
	    liwork = minmn * 3 * nlvl + minmn * 11;
#line 309 "zgelsd.f"
	    mm = *m;
#line 310 "zgelsd.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns. */

#line 315 "zgelsd.f"
		mm = *n;
/* Computing MAX */
#line 316 "zgelsd.f"
		i__1 = maxwrk, i__2 = *n * ilaenv_(&c__1, "ZGEQRF", " ", m, n,
			 &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 316 "zgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 318 "zgelsd.f"
		i__1 = maxwrk, i__2 = *nrhs * ilaenv_(&c__1, "ZUNMQR", "LC", 
			m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 318 "zgelsd.f"
		maxwrk = max(i__1,i__2);
#line 320 "zgelsd.f"
	    }
#line 321 "zgelsd.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined. */

/* Computing MAX */
/* Computing 2nd power */
#line 325 "zgelsd.f"
		i__3 = smlsiz + 1;
#line 325 "zgelsd.f"
		i__1 = i__3 * i__3, i__2 = *n * (*nrhs + 1) + (*nrhs << 1);
#line 325 "zgelsd.f"
		lrwork = *n * 10 + (*n << 1) * smlsiz + (*n << 3) * nlvl + 
			smlsiz * 3 * *nrhs + max(i__1,i__2);
/* Computing MAX */
#line 327 "zgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (mm + *n) * ilaenv_(&c__1, 
			"ZGEBRD", " ", &mm, n, &c_n1, &c_n1, (ftnlen)6, (
			ftnlen)1);
#line 327 "zgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 329 "zgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + *nrhs * ilaenv_(&c__1, 
			"ZUNMBR", "QLC", &mm, nrhs, n, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 329 "zgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 331 "zgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"ZUNMBR", "PLN", n, nrhs, n, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 331 "zgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 333 "zgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + *n * *nrhs;
#line 333 "zgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 334 "zgelsd.f"
		i__1 = (*n << 1) + mm, i__2 = (*n << 1) + *n * *nrhs;
#line 334 "zgelsd.f"
		minwrk = max(i__1,i__2);
#line 335 "zgelsd.f"
	    }
#line 336 "zgelsd.f"
	    if (*n > *m) {
/* Computing MAX */
/* Computing 2nd power */
#line 337 "zgelsd.f"
		i__3 = smlsiz + 1;
#line 337 "zgelsd.f"
		i__1 = i__3 * i__3, i__2 = *n * (*nrhs + 1) + (*nrhs << 1);
#line 337 "zgelsd.f"
		lrwork = *m * 10 + (*m << 1) * smlsiz + (*m << 3) * nlvl + 
			smlsiz * 3 * *nrhs + max(i__1,i__2);
#line 339 "zgelsd.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                           than rows. */

#line 344 "zgelsd.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "ZGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 346 "zgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			    ilaenv_(&c__1, "ZGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 346 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 348 "zgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * 
			    ilaenv_(&c__1, "ZUNMBR", "QLC", m, nrhs, m, &c_n1,
			     (ftnlen)6, (ftnlen)3);
#line 348 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 350 "zgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			    ilaenv_(&c__1, "ZUNMLQ", "LC", n, nrhs, m, &c_n1, 
			    (ftnlen)6, (ftnlen)2);
#line 350 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 352 "zgelsd.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 353 "zgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 353 "zgelsd.f"
			maxwrk = max(i__1,i__2);
#line 354 "zgelsd.f"
		    } else {
/* Computing MAX */
#line 355 "zgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 355 "zgelsd.f"
			maxwrk = max(i__1,i__2);
#line 356 "zgelsd.f"
		    }
/* Computing MAX */
#line 357 "zgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *m * *nrhs;
#line 357 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
/*     XXX: Ensure the Path 2a case below is triggered.  The workspace */
/*     calculation should use queries for all routines eventually. */
/* Computing MAX */
/* Computing MAX */
#line 360 "zgelsd.f"
		    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), 
			    i__3 = max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 360 "zgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 2) + *m * *m + max(i__3,i__4)
			    ;
#line 360 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 362 "zgelsd.f"
		} else {

/*                 Path 2 - underdetermined. */

#line 366 "zgelsd.f"
		    maxwrk = (*m << 1) + (*n + *m) * ilaenv_(&c__1, "ZGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 368 "zgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *nrhs * ilaenv_(&c__1, 
			    "ZUNMBR", "QLC", m, nrhs, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 368 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 370 "zgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "ZUNMBR", "PLN", n, nrhs, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 370 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 372 "zgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * *nrhs;
#line 372 "zgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 373 "zgelsd.f"
		}
/* Computing MAX */
#line 374 "zgelsd.f"
		i__1 = (*m << 1) + *n, i__2 = (*m << 1) + *m * *nrhs;
#line 374 "zgelsd.f"
		minwrk = max(i__1,i__2);
#line 375 "zgelsd.f"
	    }
#line 376 "zgelsd.f"
	}
#line 377 "zgelsd.f"
	minwrk = min(minwrk,maxwrk);
#line 378 "zgelsd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 379 "zgelsd.f"
	iwork[1] = liwork;
#line 380 "zgelsd.f"
	rwork[1] = (doublereal) lrwork;

#line 382 "zgelsd.f"
	if (*lwork < minwrk && ! lquery) {
#line 383 "zgelsd.f"
	    *info = -12;
#line 384 "zgelsd.f"
	}
#line 385 "zgelsd.f"
    }

#line 387 "zgelsd.f"
    if (*info != 0) {
#line 388 "zgelsd.f"
	i__1 = -(*info);
#line 388 "zgelsd.f"
	xerbla_("ZGELSD", &i__1, (ftnlen)6);
#line 389 "zgelsd.f"
	return 0;
#line 390 "zgelsd.f"
    } else if (lquery) {
#line 391 "zgelsd.f"
	return 0;
#line 392 "zgelsd.f"
    }

/*     Quick return if possible. */

#line 396 "zgelsd.f"
    if (*m == 0 || *n == 0) {
#line 397 "zgelsd.f"
	*rank = 0;
#line 398 "zgelsd.f"
	return 0;
#line 399 "zgelsd.f"
    }

/*     Get machine parameters. */

#line 403 "zgelsd.f"
    eps = dlamch_("P", (ftnlen)1);
#line 404 "zgelsd.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 405 "zgelsd.f"
    smlnum = sfmin / eps;
#line 406 "zgelsd.f"
    bignum = 1. / smlnum;
#line 407 "zgelsd.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A if max entry outside range [SMLNUM,BIGNUM]. */

#line 411 "zgelsd.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 412 "zgelsd.f"
    iascl = 0;
#line 413 "zgelsd.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 417 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 418 "zgelsd.f"
	iascl = 1;
#line 419 "zgelsd.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 423 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 424 "zgelsd.f"
	iascl = 2;
#line 425 "zgelsd.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 429 "zgelsd.f"
	i__1 = max(*m,*n);
#line 429 "zgelsd.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 430 "zgelsd.f"
	dlaset_("F", &minmn, &c__1, &c_b80, &c_b80, &s[1], &c__1, (ftnlen)1);
#line 431 "zgelsd.f"
	*rank = 0;
#line 432 "zgelsd.f"
	goto L10;
#line 433 "zgelsd.f"
    }

/*     Scale B if max entry outside range [SMLNUM,BIGNUM]. */

#line 437 "zgelsd.f"
    bnrm = zlange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 438 "zgelsd.f"
    ibscl = 0;
#line 439 "zgelsd.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 443 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 444 "zgelsd.f"
	ibscl = 1;
#line 445 "zgelsd.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 449 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 450 "zgelsd.f"
	ibscl = 2;
#line 451 "zgelsd.f"
    }

/*     If M < N make sure B(M+1:N,:) = 0 */

#line 455 "zgelsd.f"
    if (*m < *n) {
#line 455 "zgelsd.f"
	i__1 = *n - *m;
#line 455 "zgelsd.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, (
		ftnlen)1);
#line 455 "zgelsd.f"
    }

/*     Overdetermined case. */

#line 460 "zgelsd.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined. */

#line 464 "zgelsd.f"
	mm = *m;
#line 465 "zgelsd.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 469 "zgelsd.f"
	    mm = *n;
#line 470 "zgelsd.f"
	    itau = 1;
#line 471 "zgelsd.f"
	    nwork = itau + *n;

/*           Compute A=Q*R. */
/*           (RWorkspace: need N) */
/*           (CWorkspace: need N, prefer N*NB) */

#line 477 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 477 "zgelsd.f"
	    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);

/*           Multiply B by transpose(Q). */
/*           (RWorkspace: need N) */
/*           (CWorkspace: need NRHS, prefer NRHS*NB) */

#line 484 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 484 "zgelsd.f"
	    zunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R. */

#line 489 "zgelsd.f"
	    if (*n > 1) {
#line 490 "zgelsd.f"
		i__1 = *n - 1;
#line 490 "zgelsd.f"
		i__2 = *n - 1;
#line 490 "zgelsd.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 492 "zgelsd.f"
	    }
#line 493 "zgelsd.f"
	}

#line 495 "zgelsd.f"
	itauq = 1;
#line 496 "zgelsd.f"
	itaup = itauq + *n;
#line 497 "zgelsd.f"
	nwork = itaup + *n;
#line 498 "zgelsd.f"
	ie = 1;
#line 499 "zgelsd.f"
	nrwork = ie + *n;

/*        Bidiagonalize R in A. */
/*        (RWorkspace: need N) */
/*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB) */

#line 505 "zgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 505 "zgelsd.f"
	zgebrd_(&mm, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], &
		work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R. */
/*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB) */

#line 512 "zgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 512 "zgelsd.f"
	zunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 517 "zgelsd.f"
	zlalsd_("U", &smlsiz, n, nrhs, &s[1], &rwork[ie], &b[b_offset], ldb, 
		rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1], info, (
		ftnlen)1);
#line 520 "zgelsd.f"
	if (*info != 0) {
#line 521 "zgelsd.f"
	    goto L10;
#line 522 "zgelsd.f"
	}

/*        Multiply B by right bidiagonalizing vectors of R. */

#line 526 "zgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 526 "zgelsd.f"
	zunmbr_("P", "L", "N", n, nrhs, n, &a[a_offset], lda, &work[itaup], &
		b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

#line 529 "zgelsd.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 529 "zgelsd.f"
	i__1 = *m, i__2 = (*m << 1) - 4, i__1 = max(i__1,i__2), i__1 = max(
		i__1,*nrhs), i__2 = *n - *m * 3;
#line 529 "zgelsd.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__1,i__2)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm. */

#line 535 "zgelsd.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 536 "zgelsd.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 536 "zgelsd.f"
	    i__1 = (*m << 2) + *m * *lda + max(i__3,i__4), i__2 = *m * *lda + 
		    *m + *m * *nrhs;
#line 536 "zgelsd.f"
	    if (*lwork >= max(i__1,i__2)) {
#line 536 "zgelsd.f"
		ldwork = *lda;
#line 536 "zgelsd.f"
	    }
#line 538 "zgelsd.f"
	    itau = 1;
#line 539 "zgelsd.f"
	    nwork = *m + 1;

/*        Compute A=L*Q. */
/*        (CWorkspace: need 2*M, prefer M+M*NB) */

#line 544 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 544 "zgelsd.f"
	    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);
#line 546 "zgelsd.f"
	    il = nwork;

/*        Copy L to WORK(IL), zeroing out above its diagonal. */

#line 550 "zgelsd.f"
	    zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 551 "zgelsd.f"
	    i__1 = *m - 1;
#line 551 "zgelsd.f"
	    i__2 = *m - 1;
#line 551 "zgelsd.f"
	    zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 553 "zgelsd.f"
	    itauq = il + ldwork * *m;
#line 554 "zgelsd.f"
	    itaup = itauq + *m;
#line 555 "zgelsd.f"
	    nwork = itaup + *m;
#line 556 "zgelsd.f"
	    ie = 1;
#line 557 "zgelsd.f"
	    nrwork = ie + *m;

/*        Bidiagonalize L in WORK(IL). */
/*        (RWorkspace: need M) */
/*        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB) */

#line 563 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 563 "zgelsd.f"
	    zgebrd_(m, m, &work[il], &ldwork, &s[1], &rwork[ie], &work[itauq],
		     &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L. */
/*        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 570 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 570 "zgelsd.f"
	    zunmbr_("Q", "L", "C", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 576 "zgelsd.f"
	    zlalsd_("U", &smlsiz, m, nrhs, &s[1], &rwork[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1],
		     info, (ftnlen)1);
#line 579 "zgelsd.f"
	    if (*info != 0) {
#line 580 "zgelsd.f"
		goto L10;
#line 581 "zgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of L. */

#line 585 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 585 "zgelsd.f"
	    zunmbr_("P", "L", "N", m, nrhs, m, &work[il], &ldwork, &work[
		    itaup], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Zero out below first M rows of B. */

#line 591 "zgelsd.f"
	    i__1 = *n - *m;
#line 591 "zgelsd.f"
	    zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, 
		    (ftnlen)1);
#line 592 "zgelsd.f"
	    nwork = itau + *m;

/*        Multiply transpose(Q) by B. */
/*        (CWorkspace: need NRHS, prefer NRHS*NB) */

#line 597 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 597 "zgelsd.f"
	    zunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 600 "zgelsd.f"
	} else {

/*        Path 2 - remaining underdetermined cases. */

#line 604 "zgelsd.f"
	    itauq = 1;
#line 605 "zgelsd.f"
	    itaup = itauq + *m;
#line 606 "zgelsd.f"
	    nwork = itaup + *m;
#line 607 "zgelsd.f"
	    ie = 1;
#line 608 "zgelsd.f"
	    nrwork = ie + *m;

/*        Bidiagonalize A. */
/*        (RWorkspace: need M) */
/*        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */

#line 614 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 614 "zgelsd.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors. */
/*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB) */

#line 621 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 621 "zgelsd.f"
	    zunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 626 "zgelsd.f"
	    zlalsd_("L", &smlsiz, m, nrhs, &s[1], &rwork[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1],
		     info, (ftnlen)1);
#line 629 "zgelsd.f"
	    if (*info != 0) {
#line 630 "zgelsd.f"
		goto L10;
#line 631 "zgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of A. */

#line 635 "zgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 635 "zgelsd.f"
	    zunmbr_("P", "L", "N", n, nrhs, m, &a[a_offset], lda, &work[itaup]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

#line 638 "zgelsd.f"
	}
#line 638 "zgelsd.f"
    }

/*     Undo scaling. */

#line 642 "zgelsd.f"
    if (iascl == 1) {
#line 643 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 644 "zgelsd.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 646 "zgelsd.f"
    } else if (iascl == 2) {
#line 647 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 648 "zgelsd.f"
	dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 650 "zgelsd.f"
    }
#line 651 "zgelsd.f"
    if (ibscl == 1) {
#line 652 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 653 "zgelsd.f"
    } else if (ibscl == 2) {
#line 654 "zgelsd.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 655 "zgelsd.f"
    }

#line 657 "zgelsd.f"
L10:
#line 658 "zgelsd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 659 "zgelsd.f"
    iwork[1] = liwork;
#line 660 "zgelsd.f"
    rwork[1] = (doublereal) lrwork;
#line 661 "zgelsd.f"
    return 0;

/*     End of ZGELSD */

} /* zgelsd_ */

