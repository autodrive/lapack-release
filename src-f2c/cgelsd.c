#line 1 "cgelsd.f"
/* cgelsd.f -- translated by f2c (version 20100827).
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

#line 1 "cgelsd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b80 = 0.;

/* > \brief <b> CGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices</b
> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelsd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelsd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelsd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELSD computes the minimum-norm solution to a real linear least */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK must be at least 1. */
/* >          The exact minimum amount of workspace needed depends on M, */
/* >          N and NRHS. As long as LWORK is at least */
/* >              2 * N + N * NRHS */
/* >          if M is greater than or equal to N or */
/* >              2 * M + M * NRHS */
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
/* >          RWORK is REAL array, dimension (MAX(1,LRWORK)) */
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

/* > \ingroup complexGEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* >       California at Berkeley, USA \n */
/* >     Osni Marques, LBNL/NERSC, USA \n */

/*  ===================================================================== */
/* Subroutine */ int cgelsd_(integer *m, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int cgebrd_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), slabad_(
	    doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clalsd_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, integer *,
	     doublecomplex *, doublereal *, integer *, integer *, ftnlen), 
	    clascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), cgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), cunmbr_(char *, char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), cunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer ldwork;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static integer liwork, minwrk, maxwrk;
    static doublereal smlnum;
    static integer lrwork;
    static logical lquery;
    static integer nrwork, smlsiz;


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

#line 276 "cgelsd.f"
    /* Parameter adjustments */
#line 276 "cgelsd.f"
    a_dim1 = *lda;
#line 276 "cgelsd.f"
    a_offset = 1 + a_dim1;
#line 276 "cgelsd.f"
    a -= a_offset;
#line 276 "cgelsd.f"
    b_dim1 = *ldb;
#line 276 "cgelsd.f"
    b_offset = 1 + b_dim1;
#line 276 "cgelsd.f"
    b -= b_offset;
#line 276 "cgelsd.f"
    --s;
#line 276 "cgelsd.f"
    --work;
#line 276 "cgelsd.f"
    --rwork;
#line 276 "cgelsd.f"
    --iwork;
#line 276 "cgelsd.f"

#line 276 "cgelsd.f"
    /* Function Body */
#line 276 "cgelsd.f"
    *info = 0;
#line 277 "cgelsd.f"
    minmn = min(*m,*n);
#line 278 "cgelsd.f"
    maxmn = max(*m,*n);
#line 279 "cgelsd.f"
    lquery = *lwork == -1;
#line 280 "cgelsd.f"
    if (*m < 0) {
#line 281 "cgelsd.f"
	*info = -1;
#line 282 "cgelsd.f"
    } else if (*n < 0) {
#line 283 "cgelsd.f"
	*info = -2;
#line 284 "cgelsd.f"
    } else if (*nrhs < 0) {
#line 285 "cgelsd.f"
	*info = -3;
#line 286 "cgelsd.f"
    } else if (*lda < max(1,*m)) {
#line 287 "cgelsd.f"
	*info = -5;
#line 288 "cgelsd.f"
    } else if (*ldb < max(1,maxmn)) {
#line 289 "cgelsd.f"
	*info = -7;
#line 290 "cgelsd.f"
    }

/*     Compute workspace. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 299 "cgelsd.f"
    if (*info == 0) {
#line 300 "cgelsd.f"
	minwrk = 1;
#line 301 "cgelsd.f"
	maxwrk = 1;
#line 302 "cgelsd.f"
	liwork = 1;
#line 303 "cgelsd.f"
	lrwork = 1;
#line 304 "cgelsd.f"
	if (minmn > 0) {
#line 305 "cgelsd.f"
	    smlsiz = ilaenv_(&c__9, "CGELSD", " ", &c__0, &c__0, &c__0, &c__0,
		     (ftnlen)6, (ftnlen)1);
#line 306 "cgelsd.f"
	    mnthr = ilaenv_(&c__6, "CGELSD", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
/* Computing MAX */
#line 307 "cgelsd.f"
	    i__1 = (integer) (log((doublereal) minmn / (doublereal) (smlsiz + 
		    1)) / log(2.)) + 1;
#line 307 "cgelsd.f"
	    nlvl = max(i__1,0);
#line 309 "cgelsd.f"
	    liwork = minmn * 3 * nlvl + minmn * 11;
#line 310 "cgelsd.f"
	    mm = *m;
#line 311 "cgelsd.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns. */

#line 316 "cgelsd.f"
		mm = *n;
/* Computing MAX */
#line 317 "cgelsd.f"
		i__1 = maxwrk, i__2 = *n * ilaenv_(&c__1, "CGEQRF", " ", m, n,
			 &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 317 "cgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 319 "cgelsd.f"
		i__1 = maxwrk, i__2 = *nrhs * ilaenv_(&c__1, "CUNMQR", "LC", 
			m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 319 "cgelsd.f"
		maxwrk = max(i__1,i__2);
#line 321 "cgelsd.f"
	    }
#line 322 "cgelsd.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined. */

/* Computing MAX */
/* Computing 2nd power */
#line 326 "cgelsd.f"
		i__3 = smlsiz + 1;
#line 326 "cgelsd.f"
		i__1 = i__3 * i__3, i__2 = *n * (*nrhs + 1) + (*nrhs << 1);
#line 326 "cgelsd.f"
		lrwork = *n * 10 + (*n << 1) * smlsiz + (*n << 3) * nlvl + 
			smlsiz * 3 * *nrhs + max(i__1,i__2);
/* Computing MAX */
#line 328 "cgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (mm + *n) * ilaenv_(&c__1, 
			"CGEBRD", " ", &mm, n, &c_n1, &c_n1, (ftnlen)6, (
			ftnlen)1);
#line 328 "cgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 330 "cgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + *nrhs * ilaenv_(&c__1, 
			"CUNMBR", "QLC", &mm, nrhs, n, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 330 "cgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 332 "cgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, 
			"CUNMBR", "PLN", n, nrhs, n, &c_n1, (ftnlen)6, (
			ftnlen)3);
#line 332 "cgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 334 "cgelsd.f"
		i__1 = maxwrk, i__2 = (*n << 1) + *n * *nrhs;
#line 334 "cgelsd.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 335 "cgelsd.f"
		i__1 = (*n << 1) + mm, i__2 = (*n << 1) + *n * *nrhs;
#line 335 "cgelsd.f"
		minwrk = max(i__1,i__2);
#line 336 "cgelsd.f"
	    }
#line 337 "cgelsd.f"
	    if (*n > *m) {
/* Computing MAX */
/* Computing 2nd power */
#line 338 "cgelsd.f"
		i__3 = smlsiz + 1;
#line 338 "cgelsd.f"
		i__1 = i__3 * i__3, i__2 = *n * (*nrhs + 1) + (*nrhs << 1);
#line 338 "cgelsd.f"
		lrwork = *m * 10 + (*m << 1) * smlsiz + (*m << 3) * nlvl + 
			smlsiz * 3 * *nrhs + max(i__1,i__2);
#line 340 "cgelsd.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                           than rows. */

#line 345 "cgelsd.f"
		    maxwrk = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 347 "cgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			    ilaenv_(&c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1, 
			    (ftnlen)6, (ftnlen)1);
#line 347 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 349 "cgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * 
			    ilaenv_(&c__1, "CUNMBR", "QLC", m, nrhs, m, &c_n1,
			     (ftnlen)6, (ftnlen)3);
#line 349 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 351 "cgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			    ilaenv_(&c__1, "CUNMLQ", "LC", n, nrhs, m, &c_n1, 
			    (ftnlen)6, (ftnlen)2);
#line 351 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 353 "cgelsd.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 354 "cgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 354 "cgelsd.f"
			maxwrk = max(i__1,i__2);
#line 355 "cgelsd.f"
		    } else {
/* Computing MAX */
#line 356 "cgelsd.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 356 "cgelsd.f"
			maxwrk = max(i__1,i__2);
#line 357 "cgelsd.f"
		    }
/* Computing MAX */
#line 358 "cgelsd.f"
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *m * *nrhs;
#line 358 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
/*     XXX: Ensure the Path 2a case below is triggered.  The workspace */
/*     calculation should use queries for all routines eventually. */
/* Computing MAX */
/* Computing MAX */
#line 361 "cgelsd.f"
		    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), 
			    i__3 = max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 361 "cgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 2) + *m * *m + max(i__3,i__4)
			    ;
#line 361 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 363 "cgelsd.f"
		} else {

/*                 Path 2 - underdetermined. */

#line 367 "cgelsd.f"
		    maxwrk = (*m << 1) + (*n + *m) * ilaenv_(&c__1, "CGEBRD", 
			    " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 369 "cgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *nrhs * ilaenv_(&c__1, 
			    "CUNMBR", "QLC", m, nrhs, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 369 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 371 "cgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * ilaenv_(&c__1, 
			    "CUNMBR", "PLN", n, nrhs, m, &c_n1, (ftnlen)6, (
			    ftnlen)3);
#line 371 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 373 "cgelsd.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + *m * *nrhs;
#line 373 "cgelsd.f"
		    maxwrk = max(i__1,i__2);
#line 374 "cgelsd.f"
		}
/* Computing MAX */
#line 375 "cgelsd.f"
		i__1 = (*m << 1) + *n, i__2 = (*m << 1) + *m * *nrhs;
#line 375 "cgelsd.f"
		minwrk = max(i__1,i__2);
#line 376 "cgelsd.f"
	    }
#line 377 "cgelsd.f"
	}
#line 378 "cgelsd.f"
	minwrk = min(minwrk,maxwrk);
#line 379 "cgelsd.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 380 "cgelsd.f"
	iwork[1] = liwork;
#line 381 "cgelsd.f"
	rwork[1] = (doublereal) lrwork;

#line 383 "cgelsd.f"
	if (*lwork < minwrk && ! lquery) {
#line 384 "cgelsd.f"
	    *info = -12;
#line 385 "cgelsd.f"
	}
#line 386 "cgelsd.f"
    }

#line 388 "cgelsd.f"
    if (*info != 0) {
#line 389 "cgelsd.f"
	i__1 = -(*info);
#line 389 "cgelsd.f"
	xerbla_("CGELSD", &i__1, (ftnlen)6);
#line 390 "cgelsd.f"
	return 0;
#line 391 "cgelsd.f"
    } else if (lquery) {
#line 392 "cgelsd.f"
	return 0;
#line 393 "cgelsd.f"
    }

/*     Quick return if possible. */

#line 397 "cgelsd.f"
    if (*m == 0 || *n == 0) {
#line 398 "cgelsd.f"
	*rank = 0;
#line 399 "cgelsd.f"
	return 0;
#line 400 "cgelsd.f"
    }

/*     Get machine parameters. */

#line 404 "cgelsd.f"
    eps = slamch_("P", (ftnlen)1);
#line 405 "cgelsd.f"
    sfmin = slamch_("S", (ftnlen)1);
#line 406 "cgelsd.f"
    smlnum = sfmin / eps;
#line 407 "cgelsd.f"
    bignum = 1. / smlnum;
#line 408 "cgelsd.f"
    slabad_(&smlnum, &bignum);

/*     Scale A if max entry outside range [SMLNUM,BIGNUM]. */

#line 412 "cgelsd.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 413 "cgelsd.f"
    iascl = 0;
#line 414 "cgelsd.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 418 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 419 "cgelsd.f"
	iascl = 1;
#line 420 "cgelsd.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 424 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 425 "cgelsd.f"
	iascl = 2;
#line 426 "cgelsd.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 430 "cgelsd.f"
	i__1 = max(*m,*n);
#line 430 "cgelsd.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 431 "cgelsd.f"
	slaset_("F", &minmn, &c__1, &c_b80, &c_b80, &s[1], &c__1, (ftnlen)1);
#line 432 "cgelsd.f"
	*rank = 0;
#line 433 "cgelsd.f"
	goto L10;
#line 434 "cgelsd.f"
    }

/*     Scale B if max entry outside range [SMLNUM,BIGNUM]. */

#line 438 "cgelsd.f"
    bnrm = clange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 439 "cgelsd.f"
    ibscl = 0;
#line 440 "cgelsd.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 444 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 445 "cgelsd.f"
	ibscl = 1;
#line 446 "cgelsd.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 450 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 451 "cgelsd.f"
	ibscl = 2;
#line 452 "cgelsd.f"
    }

/*     If M < N make sure B(M+1:N,:) = 0 */

#line 456 "cgelsd.f"
    if (*m < *n) {
#line 456 "cgelsd.f"
	i__1 = *n - *m;
#line 456 "cgelsd.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, (
		ftnlen)1);
#line 456 "cgelsd.f"
    }

/*     Overdetermined case. */

#line 461 "cgelsd.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined. */

#line 465 "cgelsd.f"
	mm = *m;
#line 466 "cgelsd.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 470 "cgelsd.f"
	    mm = *n;
#line 471 "cgelsd.f"
	    itau = 1;
#line 472 "cgelsd.f"
	    nwork = itau + *n;

/*           Compute A=Q*R. */
/*           (RWorkspace: need N) */
/*           (CWorkspace: need N, prefer N*NB) */

#line 478 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 478 "cgelsd.f"
	    cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);

/*           Multiply B by transpose(Q). */
/*           (RWorkspace: need N) */
/*           (CWorkspace: need NRHS, prefer NRHS*NB) */

#line 485 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 485 "cgelsd.f"
	    cunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R. */

#line 490 "cgelsd.f"
	    if (*n > 1) {
#line 491 "cgelsd.f"
		i__1 = *n - 1;
#line 491 "cgelsd.f"
		i__2 = *n - 1;
#line 491 "cgelsd.f"
		claset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 493 "cgelsd.f"
	    }
#line 494 "cgelsd.f"
	}

#line 496 "cgelsd.f"
	itauq = 1;
#line 497 "cgelsd.f"
	itaup = itauq + *n;
#line 498 "cgelsd.f"
	nwork = itaup + *n;
#line 499 "cgelsd.f"
	ie = 1;
#line 500 "cgelsd.f"
	nrwork = ie + *n;

/*        Bidiagonalize R in A. */
/*        (RWorkspace: need N) */
/*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB) */

#line 506 "cgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 506 "cgelsd.f"
	cgebrd_(&mm, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], &
		work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R. */
/*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB) */

#line 513 "cgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 513 "cgelsd.f"
	cunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 518 "cgelsd.f"
	clalsd_("U", &smlsiz, n, nrhs, &s[1], &rwork[ie], &b[b_offset], ldb, 
		rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1], info, (
		ftnlen)1);
#line 521 "cgelsd.f"
	if (*info != 0) {
#line 522 "cgelsd.f"
	    goto L10;
#line 523 "cgelsd.f"
	}

/*        Multiply B by right bidiagonalizing vectors of R. */

#line 527 "cgelsd.f"
	i__1 = *lwork - nwork + 1;
#line 527 "cgelsd.f"
	cunmbr_("P", "L", "N", n, nrhs, n, &a[a_offset], lda, &work[itaup], &
		b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

#line 530 "cgelsd.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 530 "cgelsd.f"
	i__1 = *m, i__2 = (*m << 1) - 4, i__1 = max(i__1,i__2), i__1 = max(
		i__1,*nrhs), i__2 = *n - *m * 3;
#line 530 "cgelsd.f"
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__1,i__2)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm. */

#line 536 "cgelsd.f"
	    ldwork = *m;
/* Computing MAX */
/* Computing MAX */
#line 537 "cgelsd.f"
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 = 
		    max(i__3,*nrhs), i__4 = *n - *m * 3;
#line 537 "cgelsd.f"
	    i__1 = (*m << 2) + *m * *lda + max(i__3,i__4), i__2 = *m * *lda + 
		    *m + *m * *nrhs;
#line 537 "cgelsd.f"
	    if (*lwork >= max(i__1,i__2)) {
#line 537 "cgelsd.f"
		ldwork = *lda;
#line 537 "cgelsd.f"
	    }
#line 539 "cgelsd.f"
	    itau = 1;
#line 540 "cgelsd.f"
	    nwork = *m + 1;

/*        Compute A=L*Q. */
/*        (CWorkspace: need 2*M, prefer M+M*NB) */

#line 545 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 545 "cgelsd.f"
	    cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1,
		     info);
#line 547 "cgelsd.f"
	    il = nwork;

/*        Copy L to WORK(IL), zeroing out above its diagonal. */

#line 551 "cgelsd.f"
	    clacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 552 "cgelsd.f"
	    i__1 = *m - 1;
#line 552 "cgelsd.f"
	    i__2 = *m - 1;
#line 552 "cgelsd.f"
	    claset_("U", &i__1, &i__2, &c_b1, &c_b1, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 554 "cgelsd.f"
	    itauq = il + ldwork * *m;
#line 555 "cgelsd.f"
	    itaup = itauq + *m;
#line 556 "cgelsd.f"
	    nwork = itaup + *m;
#line 557 "cgelsd.f"
	    ie = 1;
#line 558 "cgelsd.f"
	    nrwork = ie + *m;

/*        Bidiagonalize L in WORK(IL). */
/*        (RWorkspace: need M) */
/*        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB) */

#line 564 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 564 "cgelsd.f"
	    cgebrd_(m, m, &work[il], &ldwork, &s[1], &rwork[ie], &work[itauq],
		     &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L. */
/*        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

#line 571 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 571 "cgelsd.f"
	    cunmbr_("Q", "L", "C", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 577 "cgelsd.f"
	    clalsd_("U", &smlsiz, m, nrhs, &s[1], &rwork[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1],
		     info, (ftnlen)1);
#line 580 "cgelsd.f"
	    if (*info != 0) {
#line 581 "cgelsd.f"
		goto L10;
#line 582 "cgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of L. */

#line 586 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 586 "cgelsd.f"
	    cunmbr_("P", "L", "N", m, nrhs, m, &work[il], &ldwork, &work[
		    itaup], &b[b_offset], ldb, &work[nwork], &i__1, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Zero out below first M rows of B. */

#line 592 "cgelsd.f"
	    i__1 = *n - *m;
#line 592 "cgelsd.f"
	    claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, 
		    (ftnlen)1);
#line 593 "cgelsd.f"
	    nwork = itau + *m;

/*        Multiply transpose(Q) by B. */
/*        (CWorkspace: need NRHS, prefer NRHS*NB) */

#line 598 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 598 "cgelsd.f"
	    cunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 601 "cgelsd.f"
	} else {

/*        Path 2 - remaining underdetermined cases. */

#line 605 "cgelsd.f"
	    itauq = 1;
#line 606 "cgelsd.f"
	    itaup = itauq + *m;
#line 607 "cgelsd.f"
	    nwork = itaup + *m;
#line 608 "cgelsd.f"
	    ie = 1;
#line 609 "cgelsd.f"
	    nrwork = ie + *m;

/*        Bidiagonalize A. */
/*        (RWorkspace: need M) */
/*        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB) */

#line 615 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 615 "cgelsd.f"
	    cgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[nwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors. */
/*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB) */

#line 622 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 622 "cgelsd.f"
	    cunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Solve the bidiagonal least squares problem. */

#line 627 "cgelsd.f"
	    clalsd_("L", &smlsiz, m, nrhs, &s[1], &rwork[ie], &b[b_offset], 
		    ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1],
		     info, (ftnlen)1);
#line 630 "cgelsd.f"
	    if (*info != 0) {
#line 631 "cgelsd.f"
		goto L10;
#line 632 "cgelsd.f"
	    }

/*        Multiply B by right bidiagonalizing vectors of A. */

#line 636 "cgelsd.f"
	    i__1 = *lwork - nwork + 1;
#line 636 "cgelsd.f"
	    cunmbr_("P", "L", "N", n, nrhs, m, &a[a_offset], lda, &work[itaup]
		    , &b[b_offset], ldb, &work[nwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

#line 639 "cgelsd.f"
	}
#line 639 "cgelsd.f"
    }

/*     Undo scaling. */

#line 643 "cgelsd.f"
    if (iascl == 1) {
#line 644 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 645 "cgelsd.f"
	slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 647 "cgelsd.f"
    } else if (iascl == 2) {
#line 648 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 649 "cgelsd.f"
	slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 651 "cgelsd.f"
    }
#line 652 "cgelsd.f"
    if (ibscl == 1) {
#line 653 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 654 "cgelsd.f"
    } else if (ibscl == 2) {
#line 655 "cgelsd.f"
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 656 "cgelsd.f"
    }

#line 658 "cgelsd.f"
L10:
#line 659 "cgelsd.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 660 "cgelsd.f"
    iwork[1] = liwork;
#line 661 "cgelsd.f"
    rwork[1] = (doublereal) lrwork;
#line 662 "cgelsd.f"
    return 0;

/*     End of CGELSD */

} /* cgelsd_ */

