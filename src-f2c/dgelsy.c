#line 1 "dgelsy.f"
/* dgelsy.f -- translated by f2c (version 20100827).
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

#line 1 "dgelsy.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b31 = 0.;
static integer c__2 = 2;
static doublereal c_b54 = 1.;

/* > \brief <b> DGELSY solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGELSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELSY computes the minimum-norm solution to a real linear least */
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
/* > */
/* > This routine is basically identical to the original xGELSX except */
/* > three differences: */
/* >   o The call to the subroutine xGEQPF has been substituted by the */
/* >     the call to the subroutine xGEQP3. This subroutine is a Blas-3 */
/* >     version of the QR factorization with column pivoting. */
/* >   o Matrix B (the right hand side) is updated with Blas-3. */
/* >   o The permutation of matrix B (the right hand side) is faster and */
/* >     more simple. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, the N-by-NRHS solution matrix X. */
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
/* >          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* >          to the front of AP, otherwise column i is a free column. */
/* >          On exit, if JPVT(i) = k, then the i-th column of AP */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          The unblocked strategy requires that: */
/* >             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ), */
/* >          where MN = min( M, N ). */
/* >          The block algorithm requires that: */
/* >             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ), */
/* >          where NB is an upper bound on the blocksize returned */
/* >          by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR, */
/* >          and DORMRZ. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: If INFO = -i, the i-th argument had an illegal value. */
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
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n */
/* >    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgelsy_(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal c1, c2, s1, s2;
    static integer nb, mn, nb1, nb2, nb3, nb4;
    static doublereal anrm, bnrm, smin, smax;
    static integer iascl, ibscl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismin, ismax;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlaic1_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal wsize;
    extern /* Subroutine */ int dgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    static integer lwkmin;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal sminpr, smaxpr, smlnum;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 250 "dgelsy.f"
    /* Parameter adjustments */
#line 250 "dgelsy.f"
    a_dim1 = *lda;
#line 250 "dgelsy.f"
    a_offset = 1 + a_dim1;
#line 250 "dgelsy.f"
    a -= a_offset;
#line 250 "dgelsy.f"
    b_dim1 = *ldb;
#line 250 "dgelsy.f"
    b_offset = 1 + b_dim1;
#line 250 "dgelsy.f"
    b -= b_offset;
#line 250 "dgelsy.f"
    --jpvt;
#line 250 "dgelsy.f"
    --work;
#line 250 "dgelsy.f"

#line 250 "dgelsy.f"
    /* Function Body */
#line 250 "dgelsy.f"
    mn = min(*m,*n);
#line 251 "dgelsy.f"
    ismin = mn + 1;
#line 252 "dgelsy.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 256 "dgelsy.f"
    *info = 0;
#line 257 "dgelsy.f"
    lquery = *lwork == -1;
#line 258 "dgelsy.f"
    if (*m < 0) {
#line 259 "dgelsy.f"
	*info = -1;
#line 260 "dgelsy.f"
    } else if (*n < 0) {
#line 261 "dgelsy.f"
	*info = -2;
#line 262 "dgelsy.f"
    } else if (*nrhs < 0) {
#line 263 "dgelsy.f"
	*info = -3;
#line 264 "dgelsy.f"
    } else if (*lda < max(1,*m)) {
#line 265 "dgelsy.f"
	*info = -5;
#line 266 "dgelsy.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 266 "dgelsy.f"
	i__1 = max(1,*m);
#line 266 "dgelsy.f"
	if (*ldb < max(i__1,*n)) {
#line 267 "dgelsy.f"
	    *info = -7;
#line 268 "dgelsy.f"
	}
#line 268 "dgelsy.f"
    }

/*     Figure out optimal block size */

#line 272 "dgelsy.f"
    if (*info == 0) {
#line 273 "dgelsy.f"
	if (mn == 0 || *nrhs == 0) {
#line 274 "dgelsy.f"
	    lwkmin = 1;
#line 275 "dgelsy.f"
	    lwkopt = 1;
#line 276 "dgelsy.f"
	} else {
#line 277 "dgelsy.f"
	    nb1 = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 278 "dgelsy.f"
	    nb2 = ilaenv_(&c__1, "DGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 279 "dgelsy.f"
	    nb3 = ilaenv_(&c__1, "DORMQR", " ", m, n, nrhs, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 280 "dgelsy.f"
	    nb4 = ilaenv_(&c__1, "DORMRQ", " ", m, n, nrhs, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
/* Computing MAX */
#line 281 "dgelsy.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 281 "dgelsy.f"
	    nb = max(i__1,nb4);
/* Computing MAX */
#line 282 "dgelsy.f"
	    i__1 = mn << 1, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = mn + 
		    *nrhs;
#line 282 "dgelsy.f"
	    lwkmin = mn + max(i__1,i__2);
/* Computing MAX */
#line 283 "dgelsy.f"
	    i__1 = lwkmin, i__2 = mn + (*n << 1) + nb * (*n + 1), i__1 = max(
		    i__1,i__2), i__2 = (mn << 1) + nb * *nrhs;
#line 283 "dgelsy.f"
	    lwkopt = max(i__1,i__2);
#line 285 "dgelsy.f"
	}
#line 286 "dgelsy.f"
	work[1] = (doublereal) lwkopt;

#line 288 "dgelsy.f"
	if (*lwork < lwkmin && ! lquery) {
#line 289 "dgelsy.f"
	    *info = -12;
#line 290 "dgelsy.f"
	}
#line 291 "dgelsy.f"
    }

#line 293 "dgelsy.f"
    if (*info != 0) {
#line 294 "dgelsy.f"
	i__1 = -(*info);
#line 294 "dgelsy.f"
	xerbla_("DGELSY", &i__1, (ftnlen)6);
#line 295 "dgelsy.f"
	return 0;
#line 296 "dgelsy.f"
    } else if (lquery) {
#line 297 "dgelsy.f"
	return 0;
#line 298 "dgelsy.f"
    }

/*     Quick return if possible */

#line 302 "dgelsy.f"
    if (mn == 0 || *nrhs == 0) {
#line 303 "dgelsy.f"
	*rank = 0;
#line 304 "dgelsy.f"
	return 0;
#line 305 "dgelsy.f"
    }

/*     Get machine parameters */

#line 309 "dgelsy.f"
    smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 310 "dgelsy.f"
    bignum = 1. / smlnum;
#line 311 "dgelsy.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max entries outside range [SMLNUM,BIGNUM] */

#line 315 "dgelsy.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 316 "dgelsy.f"
    iascl = 0;
#line 317 "dgelsy.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 321 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 322 "dgelsy.f"
	iascl = 1;
#line 323 "dgelsy.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 327 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 328 "dgelsy.f"
	iascl = 2;
#line 329 "dgelsy.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 333 "dgelsy.f"
	i__1 = max(*m,*n);
#line 333 "dgelsy.f"
	dlaset_("F", &i__1, nrhs, &c_b31, &c_b31, &b[b_offset], ldb, (ftnlen)
		1);
#line 334 "dgelsy.f"
	*rank = 0;
#line 335 "dgelsy.f"
	goto L70;
#line 336 "dgelsy.f"
    }

#line 338 "dgelsy.f"
    bnrm = dlange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 339 "dgelsy.f"
    ibscl = 0;
#line 340 "dgelsy.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 344 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 345 "dgelsy.f"
	ibscl = 1;
#line 346 "dgelsy.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 350 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 351 "dgelsy.f"
	ibscl = 2;
#line 352 "dgelsy.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 357 "dgelsy.f"
    i__1 = *lwork - mn;
#line 357 "dgelsy.f"
    dgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &i__1,
	     info);
#line 359 "dgelsy.f"
    wsize = mn + work[mn + 1];

/*     workspace: MN+2*N+NB*(N+1). */
/*     Details of Householder rotations stored in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 366 "dgelsy.f"
    work[ismin] = 1.;
#line 367 "dgelsy.f"
    work[ismax] = 1.;
#line 368 "dgelsy.f"
    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 369 "dgelsy.f"
    smin = smax;
#line 370 "dgelsy.f"
    if ((d__1 = a[a_dim1 + 1], abs(d__1)) == 0.) {
#line 371 "dgelsy.f"
	*rank = 0;
#line 372 "dgelsy.f"
	i__1 = max(*m,*n);
#line 372 "dgelsy.f"
	dlaset_("F", &i__1, nrhs, &c_b31, &c_b31, &b[b_offset], ldb, (ftnlen)
		1);
#line 373 "dgelsy.f"
	goto L70;
#line 374 "dgelsy.f"
    } else {
#line 375 "dgelsy.f"
	*rank = 1;
#line 376 "dgelsy.f"
    }

#line 378 "dgelsy.f"
L10:
#line 379 "dgelsy.f"
    if (*rank < mn) {
#line 380 "dgelsy.f"
	i__ = *rank + 1;
#line 381 "dgelsy.f"
	dlaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 383 "dgelsy.f"
	dlaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 386 "dgelsy.f"
	if (smaxpr * *rcond <= sminpr) {
#line 387 "dgelsy.f"
	    i__1 = *rank;
#line 387 "dgelsy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "dgelsy.f"
		work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
#line 389 "dgelsy.f"
		work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
#line 390 "dgelsy.f"
/* L20: */
#line 390 "dgelsy.f"
	    }
#line 391 "dgelsy.f"
	    work[ismin + *rank] = c1;
#line 392 "dgelsy.f"
	    work[ismax + *rank] = c2;
#line 393 "dgelsy.f"
	    smin = sminpr;
#line 394 "dgelsy.f"
	    smax = smaxpr;
#line 395 "dgelsy.f"
	    ++(*rank);
#line 396 "dgelsy.f"
	    goto L10;
#line 397 "dgelsy.f"
	}
#line 398 "dgelsy.f"
    }

/*     workspace: 3*MN. */

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 408 "dgelsy.f"
    if (*rank < *n) {
#line 408 "dgelsy.f"
	i__1 = *lwork - (mn << 1);
#line 408 "dgelsy.f"
	dtzrzf_(rank, n, &a[a_offset], lda, &work[mn + 1], &work[(mn << 1) + 
		1], &i__1, info);
#line 408 "dgelsy.f"
    }

/*     workspace: 2*MN. */
/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 417 "dgelsy.f"
    i__1 = *lwork - (mn << 1);
#line 417 "dgelsy.f"
    dormqr_("Left", "Transpose", m, nrhs, &mn, &a[a_offset], lda, &work[1], &
	    b[b_offset], ldb, &work[(mn << 1) + 1], &i__1, info, (ftnlen)4, (
	    ftnlen)9);
/* Computing MAX */
#line 419 "dgelsy.f"
    d__1 = wsize, d__2 = (mn << 1) + work[(mn << 1) + 1];
#line 419 "dgelsy.f"
    wsize = max(d__1,d__2);

/*     workspace: 2*MN+NB*NRHS. */

/*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 425 "dgelsy.f"
    dtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b54, &
	    a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

#line 428 "dgelsy.f"
    i__1 = *nrhs;
#line 428 "dgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 429 "dgelsy.f"
	i__2 = *n;
#line 429 "dgelsy.f"
	for (i__ = *rank + 1; i__ <= i__2; ++i__) {
#line 430 "dgelsy.f"
	    b[i__ + j * b_dim1] = 0.;
#line 431 "dgelsy.f"
/* L30: */
#line 431 "dgelsy.f"
	}
#line 432 "dgelsy.f"
/* L40: */
#line 432 "dgelsy.f"
    }

/*     B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS) */

#line 436 "dgelsy.f"
    if (*rank < *n) {
#line 437 "dgelsy.f"
	i__1 = *n - *rank;
#line 437 "dgelsy.f"
	i__2 = *lwork - (mn << 1);
#line 437 "dgelsy.f"
	dormrz_("Left", "Transpose", n, nrhs, rank, &i__1, &a[a_offset], lda, 
		&work[mn + 1], &b[b_offset], ldb, &work[(mn << 1) + 1], &i__2,
		 info, (ftnlen)4, (ftnlen)9);
#line 440 "dgelsy.f"
    }

/*     workspace: 2*MN+NRHS. */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 446 "dgelsy.f"
    i__1 = *nrhs;
#line 446 "dgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 447 "dgelsy.f"
	i__2 = *n;
#line 447 "dgelsy.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "dgelsy.f"
	    work[jpvt[i__]] = b[i__ + j * b_dim1];
#line 449 "dgelsy.f"
/* L50: */
#line 449 "dgelsy.f"
	}
#line 450 "dgelsy.f"
	dcopy_(n, &work[1], &c__1, &b[j * b_dim1 + 1], &c__1);
#line 451 "dgelsy.f"
/* L60: */
#line 451 "dgelsy.f"
    }

/*     workspace: N. */

/*     Undo scaling */

#line 457 "dgelsy.f"
    if (iascl == 1) {
#line 458 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 459 "dgelsy.f"
	dlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 461 "dgelsy.f"
    } else if (iascl == 2) {
#line 462 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 463 "dgelsy.f"
	dlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 465 "dgelsy.f"
    }
#line 466 "dgelsy.f"
    if (ibscl == 1) {
#line 467 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 468 "dgelsy.f"
    } else if (ibscl == 2) {
#line 469 "dgelsy.f"
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 470 "dgelsy.f"
    }

#line 472 "dgelsy.f"
L70:
#line 473 "dgelsy.f"
    work[1] = (doublereal) lwkopt;

#line 475 "dgelsy.f"
    return 0;

/*     End of DGELSY */

} /* dgelsy_ */

