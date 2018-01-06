#line 1 "zgelss.f"
/* zgelss.f -- translated by f2c (version 20100827).
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

#line 1 "zgelss.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b59 = 0.;

/* > \brief <b> ZGELSS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGELSS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelss.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelss.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelss.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELSS computes the minimum norm solution to a complex linear */
/* > least squares problem: */
/* > */
/* > Minimize 2-norm(| b - A*x |). */
/* > */
/* > using the singular value decomposition (SVD) of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix */
/* > X. */
/* > */
/* > The effective rank of A is determined by treating as zero those */
/* > singular values which are less than RCOND times the largest singular */
/* > value. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the first min(m,n) rows of A are overwritten with */
/* >          its right singular vectors, stored rowwise. */
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
/* >          The dimension of the array WORK. LWORK >= 1, and also: */
/* >          LWORK >=  2*min(M,N) + max(M,N,NRHS) */
/* >          For good performance, LWORK should generally be larger. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (5*min(M,N)) */
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

/* > \ingroup complex16GEsolve */

/*  ===================================================================== */
/* Subroutine */ int zgelss_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, bl, ie, il, mm;
    static doublecomplex dum[1];
    static doublereal eps, thr, anrm, bnrm;
    static integer itau, lwork_zgebrd__, lwork_zgelqf__, lwork_zgeqrf__, 
	    lwork_zungbr__, lwork_zunmbr__, iascl, ibscl, lwork_zunmlq__, 
	    chunk, lwork_zunmqr__;
    static doublereal sfmin;
    static integer minmn;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer maxmn, itaup, itauq, mnthr;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer iwork;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zgebrd_();
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zgelqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), zlascl_(char *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen), zgeqrf_(integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, doublecomplex *, integer *, integer *), zdrscl_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static integer ldwork;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zbdsqr_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer minwrk, maxwrk;
    extern /* Subroutine */ int zungbr_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static doublereal smlnum;
    static integer irwork;
    extern /* Subroutine */ int zunmbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int zunmlq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 235 "zgelss.f"
    /* Parameter adjustments */
#line 235 "zgelss.f"
    a_dim1 = *lda;
#line 235 "zgelss.f"
    a_offset = 1 + a_dim1;
#line 235 "zgelss.f"
    a -= a_offset;
#line 235 "zgelss.f"
    b_dim1 = *ldb;
#line 235 "zgelss.f"
    b_offset = 1 + b_dim1;
#line 235 "zgelss.f"
    b -= b_offset;
#line 235 "zgelss.f"
    --s;
#line 235 "zgelss.f"
    --work;
#line 235 "zgelss.f"
    --rwork;
#line 235 "zgelss.f"

#line 235 "zgelss.f"
    /* Function Body */
#line 235 "zgelss.f"
    *info = 0;
#line 236 "zgelss.f"
    minmn = min(*m,*n);
#line 237 "zgelss.f"
    maxmn = max(*m,*n);
#line 238 "zgelss.f"
    lquery = *lwork == -1;
#line 239 "zgelss.f"
    if (*m < 0) {
#line 240 "zgelss.f"
	*info = -1;
#line 241 "zgelss.f"
    } else if (*n < 0) {
#line 242 "zgelss.f"
	*info = -2;
#line 243 "zgelss.f"
    } else if (*nrhs < 0) {
#line 244 "zgelss.f"
	*info = -3;
#line 245 "zgelss.f"
    } else if (*lda < max(1,*m)) {
#line 246 "zgelss.f"
	*info = -5;
#line 247 "zgelss.f"
    } else if (*ldb < max(1,maxmn)) {
#line 248 "zgelss.f"
	*info = -7;
#line 249 "zgelss.f"
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       CWorkspace refers to complex workspace, and RWorkspace refers */
/*       to real workspace. NB refers to the optimal block size for the */
/*       immediately following subroutine, as returned by ILAENV.) */

#line 259 "zgelss.f"
    if (*info == 0) {
#line 260 "zgelss.f"
	minwrk = 1;
#line 261 "zgelss.f"
	maxwrk = 1;
#line 262 "zgelss.f"
	if (minmn > 0) {
#line 263 "zgelss.f"
	    mm = *m;
#line 264 "zgelss.f"
	    mnthr = ilaenv_(&c__6, "ZGELSS", " ", m, n, nrhs, &c_n1, (ftnlen)
		    6, (ftnlen)1);
#line 265 "zgelss.f"
	    if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns */

/*              Compute space needed for ZGEQRF */
#line 271 "zgelss.f"
		zgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 272 "zgelss.f"
		lwork_zgeqrf__ = (integer) dum[0].r;
/*              Compute space needed for ZUNMQR */
#line 274 "zgelss.f"
		zunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, dum, &b[
			b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (ftnlen)
			1);
#line 276 "zgelss.f"
		lwork_zunmqr__ = (integer) dum[0].r;
#line 277 "zgelss.f"
		mm = *n;
/* Computing MAX */
#line 278 "zgelss.f"
		i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "ZGEQRF", 
			" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 278 "zgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 280 "zgelss.f"
		i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "ZUNMQR", 
			"LC", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 280 "zgelss.f"
		maxwrk = max(i__1,i__2);
#line 282 "zgelss.f"
	    }
#line 283 "zgelss.f"
	    if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined */

/*              Compute space needed for ZGEBRD */
#line 288 "zgelss.f"
		zgebrd_(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, 
			&c_n1, info);
#line 290 "zgelss.f"
		lwork_zgebrd__ = (integer) dum[0].r;
/*              Compute space needed for ZUNMBR */
#line 292 "zgelss.f"
		zunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, dum, &
			b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 294 "zgelss.f"
		lwork_zunmbr__ = (integer) dum[0].r;
/*              Compute space needed for ZUNGBR */
#line 296 "zgelss.f"
		zungbr_("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, 
			info, (ftnlen)1);
#line 298 "zgelss.f"
		lwork_zungbr__ = (integer) dum[0].r;
/*              Compute total workspace needed */
/* Computing MAX */
#line 300 "zgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_zgebrd__;
#line 300 "zgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 301 "zgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_zunmbr__;
#line 301 "zgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 302 "zgelss.f"
		i__1 = maxwrk, i__2 = (*n << 1) + lwork_zungbr__;
#line 302 "zgelss.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 303 "zgelss.f"
		i__1 = maxwrk, i__2 = *n * *nrhs;
#line 303 "zgelss.f"
		maxwrk = max(i__1,i__2);
#line 304 "zgelss.f"
		minwrk = (*n << 1) + max(*nrhs,*m);
#line 305 "zgelss.f"
	    }
#line 306 "zgelss.f"
	    if (*n > *m) {
#line 307 "zgelss.f"
		minwrk = (*m << 1) + max(*nrhs,*n);
#line 308 "zgelss.f"
		if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                 than rows */

/*                 Compute space needed for ZGELQF */
#line 314 "zgelss.f"
		    zgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
#line 316 "zgelss.f"
		    lwork_zgelqf__ = (integer) dum[0].r;
/*                 Compute space needed for ZGEBRD */
#line 318 "zgelss.f"
		    zgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 320 "zgelss.f"
		    lwork_zgebrd__ = (integer) dum[0].r;
/*                 Compute space needed for ZUNMBR */
#line 322 "zgelss.f"
		    zunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 324 "zgelss.f"
		    lwork_zunmbr__ = (integer) dum[0].r;
/*                 Compute space needed for ZUNGBR */
#line 326 "zgelss.f"
		    zungbr_("P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 328 "zgelss.f"
		    lwork_zungbr__ = (integer) dum[0].r;
/*                 Compute space needed for ZUNMLQ */
#line 330 "zgelss.f"
		    zunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, dum, &b[
			    b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1);
#line 332 "zgelss.f"
		    lwork_zunmlq__ = (integer) dum[0].r;
/*                 Compute total workspace needed */
#line 334 "zgelss.f"
		    maxwrk = *m + lwork_zgelqf__;
/* Computing MAX */
#line 335 "zgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_zgebrd__;
#line 335 "zgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 336 "zgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_zunmbr__;
#line 336 "zgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 337 "zgelss.f"
		    i__1 = maxwrk, i__2 = *m * 3 + *m * *m + lwork_zungbr__;
#line 337 "zgelss.f"
		    maxwrk = max(i__1,i__2);
#line 338 "zgelss.f"
		    if (*nrhs > 1) {
/* Computing MAX */
#line 339 "zgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
#line 339 "zgelss.f"
			maxwrk = max(i__1,i__2);
#line 340 "zgelss.f"
		    } else {
/* Computing MAX */
#line 341 "zgelss.f"
			i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
#line 341 "zgelss.f"
			maxwrk = max(i__1,i__2);
#line 342 "zgelss.f"
		    }
/* Computing MAX */
#line 343 "zgelss.f"
		    i__1 = maxwrk, i__2 = *m + lwork_zunmlq__;
#line 343 "zgelss.f"
		    maxwrk = max(i__1,i__2);
#line 344 "zgelss.f"
		} else {

/*                 Path 2 - underdetermined */

/*                 Compute space needed for ZGEBRD */
#line 349 "zgelss.f"
		    zgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, 
			    dum, &c_n1, info);
#line 351 "zgelss.f"
		    lwork_zgebrd__ = (integer) dum[0].r;
/*                 Compute space needed for ZUNMBR */
#line 353 "zgelss.f"
		    zunmbr_("Q", "L", "C", m, nrhs, m, &a[a_offset], lda, dum,
			     &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
			    ftnlen)1, (ftnlen)1);
#line 355 "zgelss.f"
		    lwork_zunmbr__ = (integer) dum[0].r;
/*                 Compute space needed for ZUNGBR */
#line 357 "zgelss.f"
		    zungbr_("P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, 
			    info, (ftnlen)1);
#line 359 "zgelss.f"
		    lwork_zungbr__ = (integer) dum[0].r;
#line 360 "zgelss.f"
		    maxwrk = (*m << 1) + lwork_zgebrd__;
/* Computing MAX */
#line 361 "zgelss.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zunmbr__;
#line 361 "zgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 362 "zgelss.f"
		    i__1 = maxwrk, i__2 = (*m << 1) + lwork_zungbr__;
#line 362 "zgelss.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 363 "zgelss.f"
		    i__1 = maxwrk, i__2 = *n * *nrhs;
#line 363 "zgelss.f"
		    maxwrk = max(i__1,i__2);
#line 364 "zgelss.f"
		}
#line 365 "zgelss.f"
	    }
#line 366 "zgelss.f"
	    maxwrk = max(minwrk,maxwrk);
#line 367 "zgelss.f"
	}
#line 368 "zgelss.f"
	work[1].r = (doublereal) maxwrk, work[1].i = 0.;

#line 370 "zgelss.f"
	if (*lwork < minwrk && ! lquery) {
#line 370 "zgelss.f"
	    *info = -12;
#line 370 "zgelss.f"
	}
#line 372 "zgelss.f"
    }

#line 374 "zgelss.f"
    if (*info != 0) {
#line 375 "zgelss.f"
	i__1 = -(*info);
#line 375 "zgelss.f"
	xerbla_("ZGELSS", &i__1, (ftnlen)6);
#line 376 "zgelss.f"
	return 0;
#line 377 "zgelss.f"
    } else if (lquery) {
#line 378 "zgelss.f"
	return 0;
#line 379 "zgelss.f"
    }

/*     Quick return if possible */

#line 383 "zgelss.f"
    if (*m == 0 || *n == 0) {
#line 384 "zgelss.f"
	*rank = 0;
#line 385 "zgelss.f"
	return 0;
#line 386 "zgelss.f"
    }

/*     Get machine parameters */

#line 390 "zgelss.f"
    eps = dlamch_("P", (ftnlen)1);
#line 391 "zgelss.f"
    sfmin = dlamch_("S", (ftnlen)1);
#line 392 "zgelss.f"
    smlnum = sfmin / eps;
#line 393 "zgelss.f"
    bignum = 1. / smlnum;
#line 394 "zgelss.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

#line 398 "zgelss.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 399 "zgelss.f"
    iascl = 0;
#line 400 "zgelss.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 404 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 405 "zgelss.f"
	iascl = 1;
#line 406 "zgelss.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 410 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 411 "zgelss.f"
	iascl = 2;
#line 412 "zgelss.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 416 "zgelss.f"
	i__1 = max(*m,*n);
#line 416 "zgelss.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 417 "zgelss.f"
	dlaset_("F", &minmn, &c__1, &c_b59, &c_b59, &s[1], &minmn, (ftnlen)1);
#line 418 "zgelss.f"
	*rank = 0;
#line 419 "zgelss.f"
	goto L70;
#line 420 "zgelss.f"
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 424 "zgelss.f"
    bnrm = zlange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 425 "zgelss.f"
    ibscl = 0;
#line 426 "zgelss.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 430 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 431 "zgelss.f"
	ibscl = 1;
#line 432 "zgelss.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 436 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 437 "zgelss.f"
	ibscl = 2;
#line 438 "zgelss.f"
    }

/*     Overdetermined case */

#line 442 "zgelss.f"
    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

#line 446 "zgelss.f"
	mm = *m;
#line 447 "zgelss.f"
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

#line 451 "zgelss.f"
	    mm = *n;
#line 452 "zgelss.f"
	    itau = 1;
#line 453 "zgelss.f"
	    iwork = itau + *n;

/*           Compute A=Q*R */
/*           (CWorkspace: need 2*N, prefer N+N*NB) */
/*           (RWorkspace: none) */

#line 459 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 459 "zgelss.f"
	    zgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1,
		     info);

/*           Multiply B by transpose(Q) */
/*           (CWorkspace: need N+NRHS, prefer N+NRHS*NB) */
/*           (RWorkspace: none) */

#line 466 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 466 "zgelss.f"
	    zunmqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

/*           Zero out below R */

#line 471 "zgelss.f"
	    if (*n > 1) {
#line 471 "zgelss.f"
		i__1 = *n - 1;
#line 471 "zgelss.f"
		i__2 = *n - 1;
#line 471 "zgelss.f"
		zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
			(ftnlen)1);
#line 471 "zgelss.f"
	    }
#line 474 "zgelss.f"
	}

#line 476 "zgelss.f"
	ie = 1;
#line 477 "zgelss.f"
	itauq = 1;
#line 478 "zgelss.f"
	itaup = itauq + *n;
#line 479 "zgelss.f"
	iwork = itaup + *n;

/*        Bidiagonalize R in A */
/*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB) */
/*        (RWorkspace: need N) */

#line 485 "zgelss.f"
	i__1 = *lwork - iwork + 1;
#line 485 "zgelss.f"
	zgebrd_(&mm, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], &
		work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R */
/*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB) */
/*        (RWorkspace: none) */

#line 493 "zgelss.f"
	i__1 = *lwork - iwork + 1;
#line 493 "zgelss.f"
	zunmbr_("Q", "L", "C", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], 
		&b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in A */
/*        (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
/*        (RWorkspace: none) */

#line 500 "zgelss.f"
	i__1 = *lwork - iwork + 1;
#line 500 "zgelss.f"
	zungbr_("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &
		i__1, info, (ftnlen)1);
#line 502 "zgelss.f"
	irwork = ie + *n;

/*        Perform bidiagonal QR iteration */
/*          multiply B by transpose of left singular vectors */
/*          compute right singular vectors in A */
/*        (CWorkspace: none) */
/*        (RWorkspace: need BDSPAC) */

#line 510 "zgelss.f"
	zbdsqr_("U", n, n, &c__0, nrhs, &s[1], &rwork[ie], &a[a_offset], lda, 
		dum, &c__1, &b[b_offset], ldb, &rwork[irwork], info, (ftnlen)
		1);
#line 512 "zgelss.f"
	if (*info != 0) {
#line 512 "zgelss.f"
	    goto L70;
#line 512 "zgelss.f"
	}

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 517 "zgelss.f"
	d__1 = *rcond * s[1];
#line 517 "zgelss.f"
	thr = max(d__1,sfmin);
#line 518 "zgelss.f"
	if (*rcond < 0.) {
/* Computing MAX */
#line 518 "zgelss.f"
	    d__1 = eps * s[1];
#line 518 "zgelss.f"
	    thr = max(d__1,sfmin);
#line 518 "zgelss.f"
	}
#line 520 "zgelss.f"
	*rank = 0;
#line 521 "zgelss.f"
	i__1 = *n;
#line 521 "zgelss.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 522 "zgelss.f"
	    if (s[i__] > thr) {
#line 523 "zgelss.f"
		zdrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 524 "zgelss.f"
		++(*rank);
#line 525 "zgelss.f"
	    } else {
#line 526 "zgelss.f"
		zlaset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], ldb,
			 (ftnlen)1);
#line 527 "zgelss.f"
	    }
#line 528 "zgelss.f"
/* L10: */
#line 528 "zgelss.f"
	}

/*        Multiply B by right singular vectors */
/*        (CWorkspace: need N, prefer N*NRHS) */
/*        (RWorkspace: none) */

#line 534 "zgelss.f"
	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 535 "zgelss.f"
	    zgemm_("C", "N", n, nrhs, n, &c_b2, &a[a_offset], lda, &b[
		    b_offset], ldb, &c_b1, &work[1], ldb, (ftnlen)1, (ftnlen)
		    1);
#line 537 "zgelss.f"
	    zlacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1)
		    ;
#line 538 "zgelss.f"
	} else if (*nrhs > 1) {
#line 539 "zgelss.f"
	    chunk = *lwork / *n;
#line 540 "zgelss.f"
	    i__1 = *nrhs;
#line 540 "zgelss.f"
	    i__2 = chunk;
#line 540 "zgelss.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 541 "zgelss.f"
		i__3 = *nrhs - i__ + 1;
#line 541 "zgelss.f"
		bl = min(i__3,chunk);
#line 542 "zgelss.f"
		zgemm_("C", "N", n, &bl, n, &c_b2, &a[a_offset], lda, &b[i__ *
			 b_dim1 + 1], ldb, &c_b1, &work[1], n, (ftnlen)1, (
			ftnlen)1);
#line 544 "zgelss.f"
		zlacpy_("G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (
			ftnlen)1);
#line 545 "zgelss.f"
/* L20: */
#line 545 "zgelss.f"
	    }
#line 546 "zgelss.f"
	} else {
#line 547 "zgelss.f"
	    zgemv_("C", n, n, &c_b2, &a[a_offset], lda, &b[b_offset], &c__1, &
		    c_b1, &work[1], &c__1, (ftnlen)1);
#line 548 "zgelss.f"
	    zcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 549 "zgelss.f"
	}

#line 551 "zgelss.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 551 "zgelss.f"
	i__2 = max(*m,*nrhs), i__1 = *n - (*m << 1);
#line 551 "zgelss.f"
	if (*n >= mnthr && *lwork >= *m * 3 + *m * *m + max(i__2,i__1)) {

/*        Underdetermined case, M much less than N */

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm */

#line 559 "zgelss.f"
	    ldwork = *m;
/* Computing MAX */
#line 560 "zgelss.f"
	    i__2 = max(*m,*nrhs), i__1 = *n - (*m << 1);
#line 560 "zgelss.f"
	    if (*lwork >= *m * 3 + *m * *lda + max(i__2,i__1)) {
#line 560 "zgelss.f"
		ldwork = *lda;
#line 560 "zgelss.f"
	    }
#line 562 "zgelss.f"
	    itau = 1;
#line 563 "zgelss.f"
	    iwork = *m + 1;

/*        Compute A=L*Q */
/*        (CWorkspace: need 2*M, prefer M+M*NB) */
/*        (RWorkspace: none) */

#line 569 "zgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 569 "zgelss.f"
	    zgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
		     info);
#line 571 "zgelss.f"
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

#line 575 "zgelss.f"
	    zlacpy_("L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
		    1);
#line 576 "zgelss.f"
	    i__2 = *m - 1;
#line 576 "zgelss.f"
	    i__1 = *m - 1;
#line 576 "zgelss.f"
	    zlaset_("U", &i__2, &i__1, &c_b1, &c_b1, &work[il + ldwork], &
		    ldwork, (ftnlen)1);
#line 578 "zgelss.f"
	    ie = 1;
#line 579 "zgelss.f"
	    itauq = il + ldwork * *m;
#line 580 "zgelss.f"
	    itaup = itauq + *m;
#line 581 "zgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL) */
/*        (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB) */
/*        (RWorkspace: need M) */

#line 587 "zgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 587 "zgelss.f"
	    zgebrd_(m, m, &work[il], &ldwork, &s[1], &rwork[ie], &work[itauq],
		     &work[itaup], &work[iwork], &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L */
/*        (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB) */
/*        (RWorkspace: none) */

#line 595 "zgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 595 "zgelss.f"
	    zunmbr_("Q", "L", "C", m, nrhs, m, &work[il], &ldwork, &work[
		    itauq], &b[b_offset], ldb, &work[iwork], &i__2, info, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) */
/*        (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB) */
/*        (RWorkspace: none) */

#line 603 "zgelss.f"
	    i__2 = *lwork - iwork + 1;
#line 603 "zgelss.f"
	    zungbr_("P", m, m, m, &work[il], &ldwork, &work[itaup], &work[
		    iwork], &i__2, info, (ftnlen)1);
#line 605 "zgelss.f"
	    irwork = ie + *m;

/*        Perform bidiagonal QR iteration, computing right singular */
/*        vectors of L in WORK(IL) and multiplying B by transpose of */
/*        left singular vectors */
/*        (CWorkspace: need M*M) */
/*        (RWorkspace: need BDSPAC) */

#line 613 "zgelss.f"
	    zbdsqr_("U", m, m, &c__0, nrhs, &s[1], &rwork[ie], &work[il], &
		    ldwork, &a[a_offset], lda, &b[b_offset], ldb, &rwork[
		    irwork], info, (ftnlen)1);
#line 615 "zgelss.f"
	    if (*info != 0) {
#line 615 "zgelss.f"
		goto L70;
#line 615 "zgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 620 "zgelss.f"
	    d__1 = *rcond * s[1];
#line 620 "zgelss.f"
	    thr = max(d__1,sfmin);
#line 621 "zgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 621 "zgelss.f"
		d__1 = eps * s[1];
#line 621 "zgelss.f"
		thr = max(d__1,sfmin);
#line 621 "zgelss.f"
	    }
#line 623 "zgelss.f"
	    *rank = 0;
#line 624 "zgelss.f"
	    i__2 = *m;
#line 624 "zgelss.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 625 "zgelss.f"
		if (s[i__] > thr) {
#line 626 "zgelss.f"
		    zdrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 627 "zgelss.f"
		    ++(*rank);
#line 628 "zgelss.f"
		} else {
#line 629 "zgelss.f"
		    zlaset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], 
			    ldb, (ftnlen)1);
#line 630 "zgelss.f"
		}
#line 631 "zgelss.f"
/* L30: */
#line 631 "zgelss.f"
	    }
#line 632 "zgelss.f"
	    iwork = il + *m * ldwork;

/*        Multiply B by right singular vectors of L in WORK(IL) */
/*        (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS) */
/*        (RWorkspace: none) */

#line 638 "zgelss.f"
	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
#line 639 "zgelss.f"
		zgemm_("C", "N", m, nrhs, m, &c_b2, &work[il], &ldwork, &b[
			b_offset], ldb, &c_b1, &work[iwork], ldb, (ftnlen)1, (
			ftnlen)1);
#line 641 "zgelss.f"
		zlacpy_("G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 642 "zgelss.f"
	    } else if (*nrhs > 1) {
#line 643 "zgelss.f"
		chunk = (*lwork - iwork + 1) / *m;
#line 644 "zgelss.f"
		i__2 = *nrhs;
#line 644 "zgelss.f"
		i__1 = chunk;
#line 644 "zgelss.f"
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__1) {
/* Computing MIN */
#line 645 "zgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 645 "zgelss.f"
		    bl = min(i__3,chunk);
#line 646 "zgelss.f"
		    zgemm_("C", "N", m, &bl, m, &c_b2, &work[il], &ldwork, &b[
			    i__ * b_dim1 + 1], ldb, &c_b1, &work[iwork], m, (
			    ftnlen)1, (ftnlen)1);
#line 648 "zgelss.f"
		    zlacpy_("G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1]
			    , ldb, (ftnlen)1);
#line 650 "zgelss.f"
/* L40: */
#line 650 "zgelss.f"
		}
#line 651 "zgelss.f"
	    } else {
#line 652 "zgelss.f"
		zgemv_("C", m, m, &c_b2, &work[il], &ldwork, &b[b_dim1 + 1], &
			c__1, &c_b1, &work[iwork], &c__1, (ftnlen)1);
#line 654 "zgelss.f"
		zcopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
#line 655 "zgelss.f"
	    }

/*        Zero out below first M rows of B */

#line 659 "zgelss.f"
	    i__1 = *n - *m;
#line 659 "zgelss.f"
	    zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[*m + 1 + b_dim1], ldb, 
		    (ftnlen)1);
#line 660 "zgelss.f"
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B */
/*        (CWorkspace: need M+NRHS, prefer M+NHRS*NB) */
/*        (RWorkspace: none) */

#line 666 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 666 "zgelss.f"
	    zunmlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
		    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
		    ftnlen)1);

#line 669 "zgelss.f"
	} else {

/*        Path 2 - remaining underdetermined cases */

#line 673 "zgelss.f"
	    ie = 1;
#line 674 "zgelss.f"
	    itauq = 1;
#line 675 "zgelss.f"
	    itaup = itauq + *m;
#line 676 "zgelss.f"
	    iwork = itaup + *m;

/*        Bidiagonalize A */
/*        (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB) */
/*        (RWorkspace: need N) */

#line 682 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 682 "zgelss.f"
	    zgebrd_(m, n, &a[a_offset], lda, &s[1], &rwork[ie], &work[itauq], 
		    &work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors */
/*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB) */
/*        (RWorkspace: none) */

#line 690 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 690 "zgelss.f"
	    zunmbr_("Q", "L", "C", m, nrhs, n, &a[a_offset], lda, &work[itauq]
		    , &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1,
		     (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors in A */
/*        (CWorkspace: need 3*M, prefer 2*M+M*NB) */
/*        (RWorkspace: none) */

#line 697 "zgelss.f"
	    i__1 = *lwork - iwork + 1;
#line 697 "zgelss.f"
	    zungbr_("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
		    iwork], &i__1, info, (ftnlen)1);
#line 699 "zgelss.f"
	    irwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of A in A and */
/*           multiplying B by transpose of left singular vectors */
/*        (CWorkspace: none) */
/*        (RWorkspace: need BDSPAC) */

#line 707 "zgelss.f"
	    zbdsqr_("L", m, n, &c__0, nrhs, &s[1], &rwork[ie], &a[a_offset], 
		    lda, dum, &c__1, &b[b_offset], ldb, &rwork[irwork], info, 
		    (ftnlen)1);
#line 709 "zgelss.f"
	    if (*info != 0) {
#line 709 "zgelss.f"
		goto L70;
#line 709 "zgelss.f"
	    }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
#line 714 "zgelss.f"
	    d__1 = *rcond * s[1];
#line 714 "zgelss.f"
	    thr = max(d__1,sfmin);
#line 715 "zgelss.f"
	    if (*rcond < 0.) {
/* Computing MAX */
#line 715 "zgelss.f"
		d__1 = eps * s[1];
#line 715 "zgelss.f"
		thr = max(d__1,sfmin);
#line 715 "zgelss.f"
	    }
#line 717 "zgelss.f"
	    *rank = 0;
#line 718 "zgelss.f"
	    i__1 = *m;
#line 718 "zgelss.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 719 "zgelss.f"
		if (s[i__] > thr) {
#line 720 "zgelss.f"
		    zdrscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
#line 721 "zgelss.f"
		    ++(*rank);
#line 722 "zgelss.f"
		} else {
#line 723 "zgelss.f"
		    zlaset_("F", &c__1, nrhs, &c_b1, &c_b1, &b[i__ + b_dim1], 
			    ldb, (ftnlen)1);
#line 724 "zgelss.f"
		}
#line 725 "zgelss.f"
/* L50: */
#line 725 "zgelss.f"
	    }

/*        Multiply B by right singular vectors of A */
/*        (CWorkspace: need N, prefer N*NRHS) */
/*        (RWorkspace: none) */

#line 731 "zgelss.f"
	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
#line 732 "zgelss.f"
		zgemm_("C", "N", n, nrhs, m, &c_b2, &a[a_offset], lda, &b[
			b_offset], ldb, &c_b1, &work[1], ldb, (ftnlen)1, (
			ftnlen)1);
#line 734 "zgelss.f"
		zlacpy_("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (
			ftnlen)1);
#line 735 "zgelss.f"
	    } else if (*nrhs > 1) {
#line 736 "zgelss.f"
		chunk = *lwork / *n;
#line 737 "zgelss.f"
		i__1 = *nrhs;
#line 737 "zgelss.f"
		i__2 = chunk;
#line 737 "zgelss.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 738 "zgelss.f"
		    i__3 = *nrhs - i__ + 1;
#line 738 "zgelss.f"
		    bl = min(i__3,chunk);
#line 739 "zgelss.f"
		    zgemm_("C", "N", n, &bl, m, &c_b2, &a[a_offset], lda, &b[
			    i__ * b_dim1 + 1], ldb, &c_b1, &work[1], n, (
			    ftnlen)1, (ftnlen)1);
#line 741 "zgelss.f"
		    zlacpy_("F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], 
			    ldb, (ftnlen)1);
#line 742 "zgelss.f"
/* L60: */
#line 742 "zgelss.f"
		}
#line 743 "zgelss.f"
	    } else {
#line 744 "zgelss.f"
		zgemv_("C", m, n, &c_b2, &a[a_offset], lda, &b[b_offset], &
			c__1, &c_b1, &work[1], &c__1, (ftnlen)1);
#line 745 "zgelss.f"
		zcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
#line 746 "zgelss.f"
	    }
#line 747 "zgelss.f"
	}
#line 747 "zgelss.f"
    }

/*     Undo scaling */

#line 751 "zgelss.f"
    if (iascl == 1) {
#line 752 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 753 "zgelss.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 755 "zgelss.f"
    } else if (iascl == 2) {
#line 756 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 757 "zgelss.f"
	dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		minmn, info, (ftnlen)1);
#line 759 "zgelss.f"
    }
#line 760 "zgelss.f"
    if (ibscl == 1) {
#line 761 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 762 "zgelss.f"
    } else if (ibscl == 2) {
#line 763 "zgelss.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 764 "zgelss.f"
    }
#line 765 "zgelss.f"
L70:
#line 766 "zgelss.f"
    work[1].r = (doublereal) maxwrk, work[1].i = 0.;
#line 767 "zgelss.f"
    return 0;

/*     End of ZGELSS */

} /* zgelss_ */

