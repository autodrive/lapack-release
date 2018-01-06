#line 1 "sgelsx.f"
/* sgelsx.f -- translated by f2c (version 20100827).
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

#line 1 "sgelsx.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b13 = 0.;
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b36 = 1.;

/* > \brief <b> SGELSX solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGELSX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelsx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelsx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelsx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, M, N, NRHS, RANK */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine SGELSY. */
/* > */
/* > SGELSX computes the minimum-norm solution to a real linear least */
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
/* > by orthogonal transformations from the right, arriving at the */
/* > complete orthogonal factorization: */
/* >    A * P = Q * [ T11 0 ] * Z */
/* >                [  0  0 ] */
/* > The minimum-norm solution is then */
/* >    X = P * Z**T [ inv(T11)*Q1**T*B ] */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          WORK is REAL array, dimension */
/* >                      (max( min(M,N)+3*N, 2*min(M,N)+NRHS )), */
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

/* > \ingroup realGEsolve */

/*  ===================================================================== */
/* Subroutine */ int sgelsx_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1, c2, s1, s2, t1, t2;
    static integer mn;
    static doublereal anrm, bnrm, smin, smax;
    static integer iascl, ibscl, ismin, ismax;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), slaic1_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), sorm2r_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen), slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), sgeqpf_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *), slaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal sminpr, smaxpr, smlnum;
    extern /* Subroutine */ int slatzm_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, ftnlen), stzrqf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 222 "sgelsx.f"
    /* Parameter adjustments */
#line 222 "sgelsx.f"
    a_dim1 = *lda;
#line 222 "sgelsx.f"
    a_offset = 1 + a_dim1;
#line 222 "sgelsx.f"
    a -= a_offset;
#line 222 "sgelsx.f"
    b_dim1 = *ldb;
#line 222 "sgelsx.f"
    b_offset = 1 + b_dim1;
#line 222 "sgelsx.f"
    b -= b_offset;
#line 222 "sgelsx.f"
    --jpvt;
#line 222 "sgelsx.f"
    --work;
#line 222 "sgelsx.f"

#line 222 "sgelsx.f"
    /* Function Body */
#line 222 "sgelsx.f"
    mn = min(*m,*n);
#line 223 "sgelsx.f"
    ismin = mn + 1;
#line 224 "sgelsx.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 228 "sgelsx.f"
    *info = 0;
#line 229 "sgelsx.f"
    if (*m < 0) {
#line 230 "sgelsx.f"
	*info = -1;
#line 231 "sgelsx.f"
    } else if (*n < 0) {
#line 232 "sgelsx.f"
	*info = -2;
#line 233 "sgelsx.f"
    } else if (*nrhs < 0) {
#line 234 "sgelsx.f"
	*info = -3;
#line 235 "sgelsx.f"
    } else if (*lda < max(1,*m)) {
#line 236 "sgelsx.f"
	*info = -5;
#line 237 "sgelsx.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 237 "sgelsx.f"
	i__1 = max(1,*m);
#line 237 "sgelsx.f"
	if (*ldb < max(i__1,*n)) {
#line 238 "sgelsx.f"
	    *info = -7;
#line 239 "sgelsx.f"
	}
#line 239 "sgelsx.f"
    }

#line 241 "sgelsx.f"
    if (*info != 0) {
#line 242 "sgelsx.f"
	i__1 = -(*info);
#line 242 "sgelsx.f"
	xerbla_("SGELSX", &i__1, (ftnlen)6);
#line 243 "sgelsx.f"
	return 0;
#line 244 "sgelsx.f"
    }

/*     Quick return if possible */

/* Computing MIN */
#line 248 "sgelsx.f"
    i__1 = min(*m,*n);
#line 248 "sgelsx.f"
    if (min(i__1,*nrhs) == 0) {
#line 249 "sgelsx.f"
	*rank = 0;
#line 250 "sgelsx.f"
	return 0;
#line 251 "sgelsx.f"
    }

/*     Get machine parameters */

#line 255 "sgelsx.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 256 "sgelsx.f"
    bignum = 1. / smlnum;
#line 257 "sgelsx.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max elements outside range [SMLNUM,BIGNUM] */

#line 261 "sgelsx.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 262 "sgelsx.f"
    iascl = 0;
#line 263 "sgelsx.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 267 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 268 "sgelsx.f"
	iascl = 1;
#line 269 "sgelsx.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 273 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 274 "sgelsx.f"
	iascl = 2;
#line 275 "sgelsx.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 279 "sgelsx.f"
	i__1 = max(*m,*n);
#line 279 "sgelsx.f"
	slaset_("F", &i__1, nrhs, &c_b13, &c_b13, &b[b_offset], ldb, (ftnlen)
		1);
#line 280 "sgelsx.f"
	*rank = 0;
#line 281 "sgelsx.f"
	goto L100;
#line 282 "sgelsx.f"
    }

#line 284 "sgelsx.f"
    bnrm = slange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 285 "sgelsx.f"
    ibscl = 0;
#line 286 "sgelsx.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 290 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 291 "sgelsx.f"
	ibscl = 1;
#line 292 "sgelsx.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 296 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 297 "sgelsx.f"
	ibscl = 2;
#line 298 "sgelsx.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 303 "sgelsx.f"
    sgeqpf_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], info);

/*     workspace 3*N. Details of Householder rotations stored */
/*     in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 310 "sgelsx.f"
    work[ismin] = 1.;
#line 311 "sgelsx.f"
    work[ismax] = 1.;
#line 312 "sgelsx.f"
    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 313 "sgelsx.f"
    smin = smax;
#line 314 "sgelsx.f"
    if ((d__1 = a[a_dim1 + 1], abs(d__1)) == 0.) {
#line 315 "sgelsx.f"
	*rank = 0;
#line 316 "sgelsx.f"
	i__1 = max(*m,*n);
#line 316 "sgelsx.f"
	slaset_("F", &i__1, nrhs, &c_b13, &c_b13, &b[b_offset], ldb, (ftnlen)
		1);
#line 317 "sgelsx.f"
	goto L100;
#line 318 "sgelsx.f"
    } else {
#line 319 "sgelsx.f"
	*rank = 1;
#line 320 "sgelsx.f"
    }

#line 322 "sgelsx.f"
L10:
#line 323 "sgelsx.f"
    if (*rank < mn) {
#line 324 "sgelsx.f"
	i__ = *rank + 1;
#line 325 "sgelsx.f"
	slaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 327 "sgelsx.f"
	slaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 330 "sgelsx.f"
	if (smaxpr * *rcond <= sminpr) {
#line 331 "sgelsx.f"
	    i__1 = *rank;
#line 331 "sgelsx.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 332 "sgelsx.f"
		work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
#line 333 "sgelsx.f"
		work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
#line 334 "sgelsx.f"
/* L20: */
#line 334 "sgelsx.f"
	    }
#line 335 "sgelsx.f"
	    work[ismin + *rank] = c1;
#line 336 "sgelsx.f"
	    work[ismax + *rank] = c2;
#line 337 "sgelsx.f"
	    smin = sminpr;
#line 338 "sgelsx.f"
	    smax = smaxpr;
#line 339 "sgelsx.f"
	    ++(*rank);
#line 340 "sgelsx.f"
	    goto L10;
#line 341 "sgelsx.f"
	}
#line 342 "sgelsx.f"
    }

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 350 "sgelsx.f"
    if (*rank < *n) {
#line 350 "sgelsx.f"
	stzrqf_(rank, n, &a[a_offset], lda, &work[mn + 1], info);
#line 350 "sgelsx.f"
    }

/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 357 "sgelsx.f"
    sorm2r_("Left", "Transpose", m, nrhs, &mn, &a[a_offset], lda, &work[1], &
	    b[b_offset], ldb, &work[(mn << 1) + 1], info, (ftnlen)4, (ftnlen)
	    9);

/*     workspace NRHS */

/*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 364 "sgelsx.f"
    strsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &
	    a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

#line 367 "sgelsx.f"
    i__1 = *n;
#line 367 "sgelsx.f"
    for (i__ = *rank + 1; i__ <= i__1; ++i__) {
#line 368 "sgelsx.f"
	i__2 = *nrhs;
#line 368 "sgelsx.f"
	for (j = 1; j <= i__2; ++j) {
#line 369 "sgelsx.f"
	    b[i__ + j * b_dim1] = 0.;
#line 370 "sgelsx.f"
/* L30: */
#line 370 "sgelsx.f"
	}
#line 371 "sgelsx.f"
/* L40: */
#line 371 "sgelsx.f"
    }

/*     B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS) */

#line 375 "sgelsx.f"
    if (*rank < *n) {
#line 376 "sgelsx.f"
	i__1 = *rank;
#line 376 "sgelsx.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "sgelsx.f"
	    i__2 = *n - *rank + 1;
#line 377 "sgelsx.f"
	    slatzm_("Left", &i__2, nrhs, &a[i__ + (*rank + 1) * a_dim1], lda, 
		    &work[mn + i__], &b[i__ + b_dim1], &b[*rank + 1 + b_dim1],
		     ldb, &work[(mn << 1) + 1], (ftnlen)4);
#line 380 "sgelsx.f"
/* L50: */
#line 380 "sgelsx.f"
	}
#line 381 "sgelsx.f"
    }

/*     workspace NRHS */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 387 "sgelsx.f"
    i__1 = *nrhs;
#line 387 "sgelsx.f"
    for (j = 1; j <= i__1; ++j) {
#line 388 "sgelsx.f"
	i__2 = *n;
#line 388 "sgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 389 "sgelsx.f"
	    work[(mn << 1) + i__] = 1.;
#line 390 "sgelsx.f"
/* L60: */
#line 390 "sgelsx.f"
	}
#line 391 "sgelsx.f"
	i__2 = *n;
#line 391 "sgelsx.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "sgelsx.f"
	    if (work[(mn << 1) + i__] == 1.) {
#line 393 "sgelsx.f"
		if (jpvt[i__] != i__) {
#line 394 "sgelsx.f"
		    k = i__;
#line 395 "sgelsx.f"
		    t1 = b[k + j * b_dim1];
#line 396 "sgelsx.f"
		    t2 = b[jpvt[k] + j * b_dim1];
#line 397 "sgelsx.f"
L70:
#line 398 "sgelsx.f"
		    b[jpvt[k] + j * b_dim1] = t1;
#line 399 "sgelsx.f"
		    work[(mn << 1) + k] = 0.;
#line 400 "sgelsx.f"
		    t1 = t2;
#line 401 "sgelsx.f"
		    k = jpvt[k];
#line 402 "sgelsx.f"
		    t2 = b[jpvt[k] + j * b_dim1];
#line 403 "sgelsx.f"
		    if (jpvt[k] != i__) {
#line 403 "sgelsx.f"
			goto L70;
#line 403 "sgelsx.f"
		    }
#line 405 "sgelsx.f"
		    b[i__ + j * b_dim1] = t1;
#line 406 "sgelsx.f"
		    work[(mn << 1) + k] = 0.;
#line 407 "sgelsx.f"
		}
#line 408 "sgelsx.f"
	    }
#line 409 "sgelsx.f"
/* L80: */
#line 409 "sgelsx.f"
	}
#line 410 "sgelsx.f"
/* L90: */
#line 410 "sgelsx.f"
    }

/*     Undo scaling */

#line 414 "sgelsx.f"
    if (iascl == 1) {
#line 415 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 416 "sgelsx.f"
	slascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 418 "sgelsx.f"
    } else if (iascl == 2) {
#line 419 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 420 "sgelsx.f"
	slascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 422 "sgelsx.f"
    }
#line 423 "sgelsx.f"
    if (ibscl == 1) {
#line 424 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 425 "sgelsx.f"
    } else if (ibscl == 2) {
#line 426 "sgelsx.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 427 "sgelsx.f"
    }

#line 429 "sgelsx.f"
L100:

#line 431 "sgelsx.f"
    return 0;

/*     End of SGELSX */

} /* sgelsx_ */

