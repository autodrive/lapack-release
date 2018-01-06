#line 1 "cgelsx.f"
/* cgelsx.f -- translated by f2c (version 20100827).
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

#line 1 "cgelsx.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;

/* > \brief <b> CGELSX solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGELSX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelsx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelsx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelsx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine CGELSY. */
/* > */
/* > CGELSX computes the minimum-norm solution to a complex linear least */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          RCOND is REAL */
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
/* >          WORK is COMPLEX array, dimension */
/* >                      (min(M,N) + max( N, 2*min(M,N)+NRHS )), */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
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

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgelsx_(integer *m, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    claic1_(integer *, integer *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *), cunm2r_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), slabad_(doublereal *, doublereal *);
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int clascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen), cgeqpf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int clatzm_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static doublereal sminpr;
    extern /* Subroutine */ int ctzrqf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *);
    static doublereal smaxpr, smlnum;


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

#line 233 "cgelsx.f"
    /* Parameter adjustments */
#line 233 "cgelsx.f"
    a_dim1 = *lda;
#line 233 "cgelsx.f"
    a_offset = 1 + a_dim1;
#line 233 "cgelsx.f"
    a -= a_offset;
#line 233 "cgelsx.f"
    b_dim1 = *ldb;
#line 233 "cgelsx.f"
    b_offset = 1 + b_dim1;
#line 233 "cgelsx.f"
    b -= b_offset;
#line 233 "cgelsx.f"
    --jpvt;
#line 233 "cgelsx.f"
    --work;
#line 233 "cgelsx.f"
    --rwork;
#line 233 "cgelsx.f"

#line 233 "cgelsx.f"
    /* Function Body */
#line 233 "cgelsx.f"
    mn = min(*m,*n);
#line 234 "cgelsx.f"
    ismin = mn + 1;
#line 235 "cgelsx.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 239 "cgelsx.f"
    *info = 0;
#line 240 "cgelsx.f"
    if (*m < 0) {
#line 241 "cgelsx.f"
	*info = -1;
#line 242 "cgelsx.f"
    } else if (*n < 0) {
#line 243 "cgelsx.f"
	*info = -2;
#line 244 "cgelsx.f"
    } else if (*nrhs < 0) {
#line 245 "cgelsx.f"
	*info = -3;
#line 246 "cgelsx.f"
    } else if (*lda < max(1,*m)) {
#line 247 "cgelsx.f"
	*info = -5;
#line 248 "cgelsx.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 248 "cgelsx.f"
	i__1 = max(1,*m);
#line 248 "cgelsx.f"
	if (*ldb < max(i__1,*n)) {
#line 249 "cgelsx.f"
	    *info = -7;
#line 250 "cgelsx.f"
	}
#line 250 "cgelsx.f"
    }

#line 252 "cgelsx.f"
    if (*info != 0) {
#line 253 "cgelsx.f"
	i__1 = -(*info);
#line 253 "cgelsx.f"
	xerbla_("CGELSX", &i__1, (ftnlen)6);
#line 254 "cgelsx.f"
	return 0;
#line 255 "cgelsx.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 259 "cgelsx.f"
    i__1 = min(*m,*n);
#line 259 "cgelsx.f"
    if (min(i__1,*nrhs) == 0) {
#line 260 "cgelsx.f"
	*rank = 0;
#line 261 "cgelsx.f"
	return 0;
#line 262 "cgelsx.f"
    }

/*     Get machine parameters */

#line 266 "cgelsx.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 267 "cgelsx.f"
    bignum = 1. / smlnum;
#line 268 "cgelsx.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max elements outside range [SMLNUM,BIGNUM] */

#line 272 "cgelsx.f"
    anrm = clange_("M", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 273 "cgelsx.f"
    iascl = 0;
#line 274 "cgelsx.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 278 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 279 "cgelsx.f"
	iascl = 1;
#line 280 "cgelsx.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 284 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 285 "cgelsx.f"
	iascl = 2;
#line 286 "cgelsx.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 290 "cgelsx.f"
	i__1 = max(*m,*n);
#line 290 "cgelsx.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 291 "cgelsx.f"
	*rank = 0;
#line 292 "cgelsx.f"
	goto L100;
#line 293 "cgelsx.f"
    }

#line 295 "cgelsx.f"
    bnrm = clange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 296 "cgelsx.f"
    ibscl = 0;
#line 297 "cgelsx.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 301 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 302 "cgelsx.f"
	ibscl = 1;
#line 303 "cgelsx.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 307 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 308 "cgelsx.f"
	ibscl = 2;
#line 309 "cgelsx.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 314 "cgelsx.f"
    cgeqpf_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &
	    rwork[1], info);

/*     complex workspace MN+N. Real workspace 2*N. Details of Householder */
/*     rotations stored in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 322 "cgelsx.f"
    i__1 = ismin;
#line 322 "cgelsx.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 323 "cgelsx.f"
    i__1 = ismax;
#line 323 "cgelsx.f"
    work[i__1].r = 1., work[i__1].i = 0.;
#line 324 "cgelsx.f"
    smax = z_abs(&a[a_dim1 + 1]);
#line 325 "cgelsx.f"
    smin = smax;
#line 326 "cgelsx.f"
    if (z_abs(&a[a_dim1 + 1]) == 0.) {
#line 327 "cgelsx.f"
	*rank = 0;
#line 328 "cgelsx.f"
	i__1 = max(*m,*n);
#line 328 "cgelsx.f"
	claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)1);
#line 329 "cgelsx.f"
	goto L100;
#line 330 "cgelsx.f"
    } else {
#line 331 "cgelsx.f"
	*rank = 1;
#line 332 "cgelsx.f"
    }

#line 334 "cgelsx.f"
L10:
#line 335 "cgelsx.f"
    if (*rank < mn) {
#line 336 "cgelsx.f"
	i__ = *rank + 1;
#line 337 "cgelsx.f"
	claic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 339 "cgelsx.f"
	claic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 342 "cgelsx.f"
	if (smaxpr * *rcond <= sminpr) {
#line 343 "cgelsx.f"
	    i__1 = *rank;
#line 343 "cgelsx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 344 "cgelsx.f"
		i__2 = ismin + i__ - 1;
#line 344 "cgelsx.f"
		i__3 = ismin + i__ - 1;
#line 344 "cgelsx.f"
		z__1.r = s1.r * work[i__3].r - s1.i * work[i__3].i, z__1.i = 
			s1.r * work[i__3].i + s1.i * work[i__3].r;
#line 344 "cgelsx.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 345 "cgelsx.f"
		i__2 = ismax + i__ - 1;
#line 345 "cgelsx.f"
		i__3 = ismax + i__ - 1;
#line 345 "cgelsx.f"
		z__1.r = s2.r * work[i__3].r - s2.i * work[i__3].i, z__1.i = 
			s2.r * work[i__3].i + s2.i * work[i__3].r;
#line 345 "cgelsx.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 346 "cgelsx.f"
/* L20: */
#line 346 "cgelsx.f"
	    }
#line 347 "cgelsx.f"
	    i__1 = ismin + *rank;
#line 347 "cgelsx.f"
	    work[i__1].r = c1.r, work[i__1].i = c1.i;
#line 348 "cgelsx.f"
	    i__1 = ismax + *rank;
#line 348 "cgelsx.f"
	    work[i__1].r = c2.r, work[i__1].i = c2.i;
#line 349 "cgelsx.f"
	    smin = sminpr;
#line 350 "cgelsx.f"
	    smax = smaxpr;
#line 351 "cgelsx.f"
	    ++(*rank);
#line 352 "cgelsx.f"
	    goto L10;
#line 353 "cgelsx.f"
	}
#line 354 "cgelsx.f"
    }

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 362 "cgelsx.f"
    if (*rank < *n) {
#line 362 "cgelsx.f"
	ctzrqf_(rank, n, &a[a_offset], lda, &work[mn + 1], info);
#line 362 "cgelsx.f"
    }

/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */

#line 369 "cgelsx.f"
    cunm2r_("Left", "Conjugate transpose", m, nrhs, &mn, &a[a_offset], lda, &
	    work[1], &b[b_offset], ldb, &work[(mn << 1) + 1], info, (ftnlen)4,
	     (ftnlen)19);

/*     workspace NRHS */

/*      B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 376 "cgelsx.f"
    ctrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b2, &a[
	    a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);

#line 379 "cgelsx.f"
    i__1 = *n;
#line 379 "cgelsx.f"
    for (i__ = *rank + 1; i__ <= i__1; ++i__) {
#line 380 "cgelsx.f"
	i__2 = *nrhs;
#line 380 "cgelsx.f"
	for (j = 1; j <= i__2; ++j) {
#line 381 "cgelsx.f"
	    i__3 = i__ + j * b_dim1;
#line 381 "cgelsx.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 382 "cgelsx.f"
/* L30: */
#line 382 "cgelsx.f"
	}
#line 383 "cgelsx.f"
/* L40: */
#line 383 "cgelsx.f"
    }

/*     B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS) */

#line 387 "cgelsx.f"
    if (*rank < *n) {
#line 388 "cgelsx.f"
	i__1 = *rank;
#line 388 "cgelsx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 389 "cgelsx.f"
	    i__2 = *n - *rank + 1;
#line 389 "cgelsx.f"
	    d_cnjg(&z__1, &work[mn + i__]);
#line 389 "cgelsx.f"
	    clatzm_("Left", &i__2, nrhs, &a[i__ + (*rank + 1) * a_dim1], lda, 
		    &z__1, &b[i__ + b_dim1], &b[*rank + 1 + b_dim1], ldb, &
		    work[(mn << 1) + 1], (ftnlen)4);
#line 392 "cgelsx.f"
/* L50: */
#line 392 "cgelsx.f"
	}
#line 393 "cgelsx.f"
    }

/*     workspace NRHS */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 399 "cgelsx.f"
    i__1 = *nrhs;
#line 399 "cgelsx.f"
    for (j = 1; j <= i__1; ++j) {
#line 400 "cgelsx.f"
	i__2 = *n;
#line 400 "cgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "cgelsx.f"
	    i__3 = (mn << 1) + i__;
#line 401 "cgelsx.f"
	    work[i__3].r = 1., work[i__3].i = 0.;
#line 402 "cgelsx.f"
/* L60: */
#line 402 "cgelsx.f"
	}
#line 403 "cgelsx.f"
	i__2 = *n;
#line 403 "cgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 404 "cgelsx.f"
	    i__3 = (mn << 1) + i__;
#line 404 "cgelsx.f"
	    if (work[i__3].r == 1. && work[i__3].i == 0.) {
#line 405 "cgelsx.f"
		if (jpvt[i__] != i__) {
#line 406 "cgelsx.f"
		    k = i__;
#line 407 "cgelsx.f"
		    i__3 = k + j * b_dim1;
#line 407 "cgelsx.f"
		    t1.r = b[i__3].r, t1.i = b[i__3].i;
#line 408 "cgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 408 "cgelsx.f"
		    t2.r = b[i__3].r, t2.i = b[i__3].i;
#line 409 "cgelsx.f"
L70:
#line 410 "cgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 410 "cgelsx.f"
		    b[i__3].r = t1.r, b[i__3].i = t1.i;
#line 411 "cgelsx.f"
		    i__3 = (mn << 1) + k;
#line 411 "cgelsx.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 412 "cgelsx.f"
		    t1.r = t2.r, t1.i = t2.i;
#line 413 "cgelsx.f"
		    k = jpvt[k];
#line 414 "cgelsx.f"
		    i__3 = jpvt[k] + j * b_dim1;
#line 414 "cgelsx.f"
		    t2.r = b[i__3].r, t2.i = b[i__3].i;
#line 415 "cgelsx.f"
		    if (jpvt[k] != i__) {
#line 415 "cgelsx.f"
			goto L70;
#line 415 "cgelsx.f"
		    }
#line 417 "cgelsx.f"
		    i__3 = i__ + j * b_dim1;
#line 417 "cgelsx.f"
		    b[i__3].r = t1.r, b[i__3].i = t1.i;
#line 418 "cgelsx.f"
		    i__3 = (mn << 1) + k;
#line 418 "cgelsx.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 419 "cgelsx.f"
		}
#line 420 "cgelsx.f"
	    }
#line 421 "cgelsx.f"
/* L80: */
#line 421 "cgelsx.f"
	}
#line 422 "cgelsx.f"
/* L90: */
#line 422 "cgelsx.f"
    }

/*     Undo scaling */

#line 426 "cgelsx.f"
    if (iascl == 1) {
#line 427 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 428 "cgelsx.f"
	clascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 430 "cgelsx.f"
    } else if (iascl == 2) {
#line 431 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 432 "cgelsx.f"
	clascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 434 "cgelsx.f"
    }
#line 435 "cgelsx.f"
    if (ibscl == 1) {
#line 436 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 437 "cgelsx.f"
    } else if (ibscl == 2) {
#line 438 "cgelsx.f"
	clascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 439 "cgelsx.f"
    }

#line 441 "cgelsx.f"
L100:

#line 443 "cgelsx.f"
    return 0;

/*     End of CGELSX */

} /* cgelsx_ */

