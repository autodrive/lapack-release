#line 1 "zgelsx.f"
/* zgelsx.f -- translated by f2c (version 20100827).
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

#line 1 "zgelsx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;

/* > \brief <b> ZGELSX solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGELSX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelsx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelsx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelsx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine ZGELSY. */
/* > */
/* > ZGELSX computes the minimum-norm solution to a complex linear least */
/* > squares problem: */
/* >     minimize || A * X - B || */
/* > using a complete orthogonal factorization of A.  A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > */
/* > The routine first computes a QR factorization with column pivoting: */
/* >     A * P = Q * [ R11 R12 ] */
/* >                 [  0  R22 ] */
/* > with R11 defined as the largest leading submatrix whose estimated */
/* > condition number is less than 1/RCOND.  The order of R11, RANK, */
/* > is the effective rank of A. */
/* > */
/* > Then, R22 is considered to be negligible, and R12 is annihilated */
/* > by unitary transformations from the right, arriving at the */
/* > complete orthogonal factorization: */
/* >    A * P = Q * [ T11 0 ] * Z */
/* >                [  0  0 ] */
/* > The minimum-norm solution is then */
/* >    X = P * Z**H [ inv(T11)*Q1**H*B ] */
/* >                 [        0         ] */
/* > where Q1 consists of the first RANK columns of Q. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of */
/* >          columns of matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A has been overwritten by details of its */
/* >          complete orthogonal factorization. */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, the N-by-NRHS solution matrix X. */
/* >          If m >= n and RANK = n, the residual sum-of-squares for */
/* >          the solution in the i-th column is given by the sum of */
/* >          squares of elements N+1:M in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,M,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* >          JPVT is INTEGER array, dimension (N) */
/* >          On entry, if JPVT(i) .ne. 0, the i-th column of A is an */
/* >          initial column, otherwise it is a free column.  Before */
/* >          the QR factorization of A, all initial columns are */
/* >          permuted to the leading positions; only the remaining */
/* >          free columns are moved as a result of column pivoting */
/* >          during the factorization. */
/* >          On exit, if JPVT(i) = k, then the i-th column of A*P */
/* >          was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          RCOND is used to determine the effective rank of A, which */
/* >          is defined as the order of the largest leading triangular */
/* >          submatrix R11 in the QR factorization with pivoting of A, */
/* >          whose estimated condition number < 1/RCOND. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >          The effective rank of A, i.e., the order of the submatrix */
/* >          R11.  This is the same as the order of the submatrix T11 */
/* >          in the complete orthogonal factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension */
/* >                      (min(M,N) + max( N, 2*min(M,N)+NRHS )), */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
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
/* Subroutine */ int zgelsx_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex c1, c2, s1, s2, t1, t2;
    static integer mn;
    static doublereal anrm, bnrm, smin, smax;
    static integer iascl, ibscl, ismin, ismax;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    zlaic1_(integer *, integer *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int zunm2r_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), zgeqpf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *), zlaset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen);
    static doublereal sminpr, smaxpr, smlnum;
    extern /* Subroutine */ int zlatzm_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen), ztzrqf_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *);


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

#line 233 "zgelsx.f"
    /* Parameter adjustments */
#line 233 "zgelsx.f"
    a_dim1 = *lda;
#line 233 "zgelsx.f"
    a_offset = 1 + a_dim1;
#line 233 "zgelsx.f"
    a -= a_offset;
#line 233 "zgelsx.f"
    b_dim1 = *ldb;
#line 233 "zgelsx.f"
    b_offset = 1 + b_dim1;
#line 233 "zgelsx.f"
    b -= b_offset;
#line 233 "zgelsx.f"
    --jpvt;
#line 233 "zgelsx.f"
    --work;
#line 233 "zgelsx.f"
    --rwork;
#line 233 "zgelsx.f"

#line 233 "zgelsx.f"
    /* Function Body */
#line 233 "zgelsx.f"
    mn = min(*m,*n);
#line 234 "zgelsx.f"
    ismin = mn + 1;
#line 235 "zgelsx.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 239 "zgelsx.f"
    *info = 0;
#line 240 "zgelsx.f"
    if (*m < 0) {
#line 241 "zgelsx.f"
	*info = -1;
#line 242 "zgelsx.f"
    } else if (*n < 0) {
#line 243 "zgelsx.f"
	*info = -2;
#line 244 "zgelsx.f"
    } else if (*nrhs < 0) {
#line 245 "zgelsx.f"
	*info = -3;
#line 246 "zgelsx.f"
    } else if (*lda < max(1,*m)) {
#line 247 "zgelsx.f"
	*info = -5;
#line 248 "zgelsx.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 248 "zgelsx.f"
	i__1 = max(1,*m);
#line 248 "zgelsx.f"
	if (*ldb < max(i__1,*n)) {
#line 249 "zgelsx.f"
	    *info = -7;
#line 250 "zgelsx.f"
	}
#line 250 "zgelsx.f"
    }

#line 252 "zgelsx.f"
    if (*info != 0) {
#line 253 "zgelsx.f"
	i__1 = -(*info);
#line 253 "zgelsx.f"
	xerbla_("ZGELSX", &i__1, (ftnlen)6);
#line 254 "zgelsx.f"
	return 0;
#line 255 "zgelsx.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 259 "zgelsx.f"
    i__1 = min(*m,*n);
#line 259 "zgelsx.f"
    if (min(i__1,*nrhs) == 0) {
#line 260 "zgelsx.f"
	*rank = 0;
#line 261 "zgelsx.f"
	return 0;
#line 262 "zgelsx.f"
    }

/*     Get machine parameters */

#line 266 "zgelsx.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 267 "zgelsx.f"
    bignum = 1. / smlnum;
#line 268 "zgelsx.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max elements outside range [SMLNUM,BIGNUM] */

#line 272 "zgelsx.f"
    anrm = zlange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 273 "zgelsx.f"
    iascl = 0;
#line 274 "zgelsx.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 278 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 279 "zgelsx.f"
	iascl = 1;
#line 280 "zgelsx.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 284 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 285 "zgelsx.f"
	iascl = 2;
#line 286 "zgelsx.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 290 "zgelsx.f"
	i__1 = max(*m,*n);
#line 290 "zgelsx.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 291 "zgelsx.f"
	*rank = 0;
#line 292 "zgelsx.f"
	goto L100;
#line 293 "zgelsx.f"
    }

#line 295 "zgelsx.f"
    bnrm = zlange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 296 "zgelsx.f"
    ibscl = 0;
#line 297 "zgelsx.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 301 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 302 "zgelsx.f"
	ibscl = 1;
#line 303 "zgelsx.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 307 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 308 "zgelsx.f"
	ibscl = 2;
#line 309 "zgelsx.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 314 "zgelsx.f"
    zgeqpf_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &
	    rwork[1], info);

/*     complex workspace MN+N. Real workspace 2*N. Details of Householder */
/*     rotations stored in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 322 "zgelsx.f"
    i__1 = ismin;
#line 322 "zgelsx.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 323 "zgelsx.f"
    i__1 = ismax;
#line 323 "zgelsx.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 324 "zgelsx.f"
    smax = z_abs(&a[a_dim1 + 1]);
#line 325 "zgelsx.f"
    smin = smax;
#line 326 "zgelsx.f"
    if (z_abs(&a[a_dim1 + 1]) == 0.) {
#line 327 "zgelsx.f"
	*rank = 0;
#line 328 "zgelsx.f"
	i__1 = max(*m,*n);
#line 328 "zgelsx.f"
	zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 329 "zgelsx.f"
	goto L100;
#line 330 "zgelsx.f"
    } else {
#line 331 "zgelsx.f"
	*rank = 1;
#line 332 "zgelsx.f"
    }

#line 334 "zgelsx.f"
L10:
#line 335 "zgelsx.f"
    if (*rank < mn) {
#line 336 "zgelsx.f"
	i__ = *rank + 1;
#line 337 "zgelsx.f"
	zlaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 339 "zgelsx.f"
	zlaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 342 "zgelsx.f"
	if (smaxpr * *rcond <= sminpr) {
#line 343 "zgelsx.f"
	    i__1 = *rank;
#line 343 "zgelsx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 344 "zgelsx.f"
		i__2 = ismin + i__ - 1;
#line 344 "zgelsx.f"
		i__3 = ismin + i__ - 1;
#line 344 "zgelsx.f"
		z__1.r = s1.r * work[i__3].r - s1.i * work[i__3].i, z__1.i = 
			s1.r * work[i__3].i + s1.i * work[i__3].r;
#line 344 "zgelsx.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 345 "zgelsx.f"
		i__2 = ismax + i__ - 1;
#line 345 "zgelsx.f"
		i__3 = ismax + i__ - 1;
#line 345 "zgelsx.f"
		z__1.r = s2.r * work[i__3].r - s2.i * work[i__3].i, z__1.i = 
			s2.r * work[i__3].i + s2.i * work[i__3].r;
#line 345 "zgelsx.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 346 "zgelsx.f"
/* L20: */
#line 346 "zgelsx.f"
	    }
#line 347 "zgelsx.f"
	    i__1 = ismin + *rank;
#line 347 "zgelsx.f"
	    work[i__1].r = c1.r, work[i__1].i = c1.i;
#line 348 "zgelsx.f"
	    i__1 = ismax + *rank;
#line 348 "zgelsx.f"
	    work[i__1].r = c2.r, work[i__1].i = c2.i;
#line 349 "zgelsx.f"
	    smin = sminpr;
#line 350 "zgelsx.f"
	    smax = smaxpr;
#line 351 "zgelsx.f"
	    ++(*rank);
#line 352 "zgelsx.f"
	    goto L10;
#line 353 "zgelsx.f"
	}
#line 354 "zgelsx.f"
    }

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 362 "zgelsx.f"
    if (*rank < *n) {
#line 362 "zgelsx.f"
	ztzrqf_(rank, n, &a[a_offset], lda, &work[mn + 1], info);
#line 362 "zgelsx.f"
    }

/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */

#line 369 "zgelsx.f"
    zunm2r_("Left", "Conjugate transpose", m, nrhs, &mn, &a[a_offset], lda, &
	    work[1], &b[b_offset], ldb, &work[(mn << 1) + 1], info, (ftnlen)4,
	     (ftnlen)19);

/*     workspace NRHS */

/*      B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 376 "zgelsx.f"
    ztrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b2, &a[
	    a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);

#line 379 "zgelsx.f"
    i__1 = *n;
#line 379 "zgelsx.f"
    for (i__ = *rank + 1; i__ <= i__1; ++i__) {
#line 380 "zgelsx.f"
	i__2 = *nrhs;
#line 380 "zgelsx.f"
	for (j = 1; j <= i__2; ++j) {
#line 381 "zgelsx.f"
	    i__3 = i__ + j * b_dim1;
#line 381 "zgelsx.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 382 "zgelsx.f"
/* L30: */
#line 382 "zgelsx.f"
	}
#line 383 "zgelsx.f"
/* L40: */
#line 383 "zgelsx.f"
    }

/*     B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS) */

#line 387 "zgelsx.f"
    if (*rank < *n) {
#line 388 "zgelsx.f"
	i__1 = *rank;
#line 388 "zgelsx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 389 "zgelsx.f"
	    i__2 = *n - *rank + 1;
#line 389 "zgelsx.f"
	    d_cnjg(&z__1, &work[mn + i__]);
#line 389 "zgelsx.f"
	    zlatzm_("Left", &i__2, nrhs, &a[i__ + (*rank + 1) * a_dim1], lda, 
		    &z__1, &b[i__ + b_dim1], &b[*rank + 1 + b_dim1], ldb, &
		    work[(mn << 1) + 1], (ftnlen)4);
#line 392 "zgelsx.f"
/* L50: */
#line 392 "zgelsx.f"
	}
#line 393 "zgelsx.f"
    }

/*     workspace NRHS */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 399 "zgelsx.f"
    i__1 = *nrhs;
#line 399 "zgelsx.f"
    for (j = 1; j <= i__1; ++j) {
#line 400 "zgelsx.f"
	i__2 = *n;
#line 400 "zgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "zgelsx.f"
	    i__3 = (mn << 1) + i__;
#line 401 "zgelsx.f"
	    work[i__3].r = 1., work[i__3].i = 0.;
#line 402 "zgelsx.f"
/* L60: */
#line 402 "zgelsx.f"
	}
#line 403 "zgelsx.f"
	i__2 = *n;
#line 403 "zgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 404 "zgelsx.f"
	    i__3 = (mn << 1) + i__;
#line 404 "zgelsx.f"
	    if (work[i__3].r == 1. && work[i__3].i == 0.) {
#line 405 "zgelsx.f"
		if (jpvt[i__] != i__) {
#line 406 "zgelsx.f"
		    k = i__;
#line 407 "zgelsx.f"
		    i__3 = k + j * b_dim1;
#line 407 "zgelsx.f"
		    t1.r = b[i__3].r, t1.i = b[i__3].i;
#line 408 "zgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 408 "zgelsx.f"
		    t2.r = b[i__3].r, t2.i = b[i__3].i;
#line 409 "zgelsx.f"
L70:
#line 410 "zgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 410 "zgelsx.f"
		    b[i__3].r = t1.r, b[i__3].i = t1.i;
#line 411 "zgelsx.f"
		    i__3 = (mn << 1) + k;
#line 411 "zgelsx.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 412 "zgelsx.f"
		    t1.r = t2.r, t1.i = t2.i;
#line 413 "zgelsx.f"
		    k = jpvt[k];
#line 414 "zgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 414 "zgelsx.f"
		    t2.r = b[i__3].r, t2.i = b[i__3].i;
#line 415 "zgelsx.f"
		    if (jpvt[k] != i__) {
#line 415 "zgelsx.f"
			goto L70;
#line 415 "zgelsx.f"
		    }
#line 417 "zgelsx.f"
		    i__3 = i__ + j * b_dim1;
#line 417 "zgelsx.f"
		    b[i__3].r = t1.r, b[i__3].i = t1.i;
#line 418 "zgelsx.f"
		    i__3 = (mn << 1) + k;
#line 418 "zgelsx.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 419 "zgelsx.f"
		}
#line 420 "zgelsx.f"
	    }
#line 421 "zgelsx.f"
/* L80: */
#line 421 "zgelsx.f"
	}
#line 422 "zgelsx.f"
/* L90: */
#line 422 "zgelsx.f"
    }

/*     Undo scaling */

#line 426 "zgelsx.f"
    if (iascl == 1) {
#line 427 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 428 "zgelsx.f"
	zlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 430 "zgelsx.f"
    } else if (iascl == 2) {
#line 431 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 432 "zgelsx.f"
	zlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 434 "zgelsx.f"
    }
#line 435 "zgelsx.f"
    if (ibscl == 1) {
#line 436 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 437 "zgelsx.f"
    } else if (ibscl == 2) {
#line 438 "zgelsx.f"
	zlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 439 "zgelsx.f"
    }

#line 441 "zgelsx.f"
L100:

#line 443 "zgelsx.f"
    return 0;

/*     End of ZGELSX */

} /* zgelsx_ */

