#line 1 "sgelsy.f"
/* sgelsy.f -- translated by f2c (version 20100827).
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

#line 1 "sgelsy.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b31 = 0.;
static integer c__2 = 2;
static doublereal c_b54 = 1.;

/* > \brief <b> SGELSY solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGELSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelsy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelsy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelsy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
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
/* > SGELSY computes the minimum-norm solution to a real linear least */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* >          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR, */
/* >          and SORMRZ. */
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

/* > \date November 2011 */

/* > \ingroup realGEsolve */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n */
/* >    E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgelsy_(integer *m, integer *n, integer *nrhs, 
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
    static integer iascl, ibscl, ismin, ismax;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal wsize;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), slaic1_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), sgeqp3_(
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), slabad_(
	    doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), slaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static integer lwkmin;
    static doublereal sminpr, smaxpr, smlnum;
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    sormrz_(char *, char *, integer *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), stzrzf_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);


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

#line 250 "sgelsy.f"
    /* Parameter adjustments */
#line 250 "sgelsy.f"
    a_dim1 = *lda;
#line 250 "sgelsy.f"
    a_offset = 1 + a_dim1;
#line 250 "sgelsy.f"
    a -= a_offset;
#line 250 "sgelsy.f"
    b_dim1 = *ldb;
#line 250 "sgelsy.f"
    b_offset = 1 + b_dim1;
#line 250 "sgelsy.f"
    b -= b_offset;
#line 250 "sgelsy.f"
    --jpvt;
#line 250 "sgelsy.f"
    --work;
#line 250 "sgelsy.f"

#line 250 "sgelsy.f"
    /* Function Body */
#line 250 "sgelsy.f"
    mn = min(*m,*n);
#line 251 "sgelsy.f"
    ismin = mn + 1;
#line 252 "sgelsy.f"
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

#line 256 "sgelsy.f"
    *info = 0;
#line 257 "sgelsy.f"
    lquery = *lwork == -1;
#line 258 "sgelsy.f"
    if (*m < 0) {
#line 259 "sgelsy.f"
	*info = -1;
#line 260 "sgelsy.f"
    } else if (*n < 0) {
#line 261 "sgelsy.f"
	*info = -2;
#line 262 "sgelsy.f"
    } else if (*nrhs < 0) {
#line 263 "sgelsy.f"
	*info = -3;
#line 264 "sgelsy.f"
    } else if (*lda < max(1,*m)) {
#line 265 "sgelsy.f"
	*info = -5;
#line 266 "sgelsy.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 266 "sgelsy.f"
	i__1 = max(1,*m);
#line 266 "sgelsy.f"
	if (*ldb < max(i__1,*n)) {
#line 267 "sgelsy.f"
	    *info = -7;
#line 268 "sgelsy.f"
	}
#line 268 "sgelsy.f"
    }

/*     Figure out optimal block size */

#line 272 "sgelsy.f"
    if (*info == 0) {
#line 273 "sgelsy.f"
	if (mn == 0 || *nrhs == 0) {
#line 274 "sgelsy.f"
	    lwkmin = 1;
#line 275 "sgelsy.f"
	    lwkopt = 1;
#line 276 "sgelsy.f"
	} else {
#line 277 "sgelsy.f"
	    nb1 = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 278 "sgelsy.f"
	    nb2 = ilaenv_(&c__1, "SGERQF", " ", m, n, &c_n1, &c_n1, (ftnlen)6,
		     (ftnlen)1);
#line 279 "sgelsy.f"
	    nb3 = ilaenv_(&c__1, "SORMQR", " ", m, n, nrhs, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
#line 280 "sgelsy.f"
	    nb4 = ilaenv_(&c__1, "SORMRQ", " ", m, n, nrhs, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
/* Computing MAX */
#line 281 "sgelsy.f"
	    i__1 = max(nb1,nb2), i__1 = max(i__1,nb3);
#line 281 "sgelsy.f"
	    nb = max(i__1,nb4);
/* Computing MAX */
#line 282 "sgelsy.f"
	    i__1 = mn << 1, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = mn + 
		    *nrhs;
#line 282 "sgelsy.f"
	    lwkmin = mn + max(i__1,i__2);
/* Computing MAX */
#line 283 "sgelsy.f"
	    i__1 = lwkmin, i__2 = mn + (*n << 1) + nb * (*n + 1), i__1 = max(
		    i__1,i__2), i__2 = (mn << 1) + nb * *nrhs;
#line 283 "sgelsy.f"
	    lwkopt = max(i__1,i__2);
#line 285 "sgelsy.f"
	}
#line 286 "sgelsy.f"
	work[1] = (doublereal) lwkopt;

#line 288 "sgelsy.f"
	if (*lwork < lwkmin && ! lquery) {
#line 289 "sgelsy.f"
	    *info = -12;
#line 290 "sgelsy.f"
	}
#line 291 "sgelsy.f"
    }

#line 293 "sgelsy.f"
    if (*info != 0) {
#line 294 "sgelsy.f"
	i__1 = -(*info);
#line 294 "sgelsy.f"
	xerbla_("SGELSY", &i__1, (ftnlen)6);
#line 295 "sgelsy.f"
	return 0;
#line 296 "sgelsy.f"
    } else if (lquery) {
#line 297 "sgelsy.f"
	return 0;
#line 298 "sgelsy.f"
    }

/*     Quick return if possible */

#line 302 "sgelsy.f"
    if (mn == 0 || *nrhs == 0) {
#line 303 "sgelsy.f"
	*rank = 0;
#line 304 "sgelsy.f"
	return 0;
#line 305 "sgelsy.f"
    }

/*     Get machine parameters */

#line 309 "sgelsy.f"
    smlnum = slamch_("S", (ftnlen)1) / slamch_("P", (ftnlen)1);
#line 310 "sgelsy.f"
    bignum = 1. / smlnum;
#line 311 "sgelsy.f"
    slabad_(&smlnum, &bignum);

/*     Scale A, B if max entries outside range [SMLNUM,BIGNUM] */

#line 315 "sgelsy.f"
    anrm = slange_("M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 316 "sgelsy.f"
    iascl = 0;
#line 317 "sgelsy.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 321 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 322 "sgelsy.f"
	iascl = 1;
#line 323 "sgelsy.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 327 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 328 "sgelsy.f"
	iascl = 2;
#line 329 "sgelsy.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 333 "sgelsy.f"
	i__1 = max(*m,*n);
#line 333 "sgelsy.f"
	slaset_("F", &i__1, nrhs, &c_b31, &c_b31, &b[b_offset], ldb, (ftnlen)
		1);
#line 334 "sgelsy.f"
	*rank = 0;
#line 335 "sgelsy.f"
	goto L70;
#line 336 "sgelsy.f"
    }

#line 338 "sgelsy.f"
    bnrm = slange_("M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 339 "sgelsy.f"
    ibscl = 0;
#line 340 "sgelsy.f"
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

#line 344 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 345 "sgelsy.f"
	ibscl = 1;
#line 346 "sgelsy.f"
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

#line 350 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 351 "sgelsy.f"
	ibscl = 2;
#line 352 "sgelsy.f"
    }

/*     Compute QR factorization with column pivoting of A: */
/*        A * P = Q * R */

#line 357 "sgelsy.f"
    i__1 = *lwork - mn;
#line 357 "sgelsy.f"
    sgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &i__1,
	     info);
#line 359 "sgelsy.f"
    wsize = mn + work[mn + 1];

/*     workspace: MN+2*N+NB*(N+1). */
/*     Details of Householder rotations stored in WORK(1:MN). */

/*     Determine RANK using incremental condition estimation */

#line 366 "sgelsy.f"
    work[ismin] = 1.;
#line 367 "sgelsy.f"
    work[ismax] = 1.;
#line 368 "sgelsy.f"
    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 369 "sgelsy.f"
    smin = smax;
#line 370 "sgelsy.f"
    if ((d__1 = a[a_dim1 + 1], abs(d__1)) == 0.) {
#line 371 "sgelsy.f"
	*rank = 0;
#line 372 "sgelsy.f"
	i__1 = max(*m,*n);
#line 372 "sgelsy.f"
	slaset_("F", &i__1, nrhs, &c_b31, &c_b31, &b[b_offset], ldb, (ftnlen)
		1);
#line 373 "sgelsy.f"
	goto L70;
#line 374 "sgelsy.f"
    } else {
#line 375 "sgelsy.f"
	*rank = 1;
#line 376 "sgelsy.f"
    }

#line 378 "sgelsy.f"
L10:
#line 379 "sgelsy.f"
    if (*rank < mn) {
#line 380 "sgelsy.f"
	i__ = *rank + 1;
#line 381 "sgelsy.f"
	slaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 383 "sgelsy.f"
	slaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[
		i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 386 "sgelsy.f"
	if (smaxpr * *rcond <= sminpr) {
#line 387 "sgelsy.f"
	    i__1 = *rank;
#line 387 "sgelsy.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "sgelsy.f"
		work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
#line 389 "sgelsy.f"
		work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
#line 390 "sgelsy.f"
/* L20: */
#line 390 "sgelsy.f"
	    }
#line 391 "sgelsy.f"
	    work[ismin + *rank] = c1;
#line 392 "sgelsy.f"
	    work[ismax + *rank] = c2;
#line 393 "sgelsy.f"
	    smin = sminpr;
#line 394 "sgelsy.f"
	    smax = smaxpr;
#line 395 "sgelsy.f"
	    ++(*rank);
#line 396 "sgelsy.f"
	    goto L10;
#line 397 "sgelsy.f"
	}
#line 398 "sgelsy.f"
    }

/*     workspace: 3*MN. */

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ] */
/*     where R11 = R(1:RANK,1:RANK) */

/*     [R11,R12] = [ T11, 0 ] * Y */

#line 408 "sgelsy.f"
    if (*rank < *n) {
#line 408 "sgelsy.f"
	i__1 = *lwork - (mn << 1);
#line 408 "sgelsy.f"
	stzrzf_(rank, n, &a[a_offset], lda, &work[mn + 1], &work[(mn << 1) + 
		1], &i__1, info);
#line 408 "sgelsy.f"
    }

/*     workspace: 2*MN. */
/*     Details of Householder rotations stored in WORK(MN+1:2*MN) */

/*     B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */

#line 417 "sgelsy.f"
    i__1 = *lwork - (mn << 1);
#line 417 "sgelsy.f"
    sormqr_("Left", "Transpose", m, nrhs, &mn, &a[a_offset], lda, &work[1], &
	    b[b_offset], ldb, &work[(mn << 1) + 1], &i__1, info, (ftnlen)4, (
	    ftnlen)9);
/* Computing MAX */
#line 419 "sgelsy.f"
    d__1 = wsize, d__2 = (mn << 1) + work[(mn << 1) + 1];
#line 419 "sgelsy.f"
    wsize = max(d__1,d__2);

/*     workspace: 2*MN+NB*NRHS. */

/*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */

#line 425 "sgelsy.f"
    strsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b54, &
	    a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

#line 428 "sgelsy.f"
    i__1 = *nrhs;
#line 428 "sgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 429 "sgelsy.f"
	i__2 = *n;
#line 429 "sgelsy.f"
	for (i__ = *rank + 1; i__ <= i__2; ++i__) {
#line 430 "sgelsy.f"
	    b[i__ + j * b_dim1] = 0.;
#line 431 "sgelsy.f"
/* L30: */
#line 431 "sgelsy.f"
	}
#line 432 "sgelsy.f"
/* L40: */
#line 432 "sgelsy.f"
    }

/*     B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS) */

#line 436 "sgelsy.f"
    if (*rank < *n) {
#line 437 "sgelsy.f"
	i__1 = *n - *rank;
#line 437 "sgelsy.f"
	i__2 = *lwork - (mn << 1);
#line 437 "sgelsy.f"
	sormrz_("Left", "Transpose", n, nrhs, rank, &i__1, &a[a_offset], lda, 
		&work[mn + 1], &b[b_offset], ldb, &work[(mn << 1) + 1], &i__2,
		 info, (ftnlen)4, (ftnlen)9);
#line 440 "sgelsy.f"
    }

/*     workspace: 2*MN+NRHS. */

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

#line 446 "sgelsy.f"
    i__1 = *nrhs;
#line 446 "sgelsy.f"
    for (j = 1; j <= i__1; ++j) {
#line 447 "sgelsy.f"
	i__2 = *n;
#line 447 "sgelsy.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "sgelsy.f"
	    work[jpvt[i__]] = b[i__ + j * b_dim1];
#line 449 "sgelsy.f"
/* L50: */
#line 449 "sgelsy.f"
	}
#line 450 "sgelsy.f"
	scopy_(n, &work[1], &c__1, &b[j * b_dim1 + 1], &c__1);
#line 451 "sgelsy.f"
/* L60: */
#line 451 "sgelsy.f"
    }

/*     workspace: N. */

/*     Undo scaling */

#line 457 "sgelsy.f"
    if (iascl == 1) {
#line 458 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 459 "sgelsy.f"
	slascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 461 "sgelsy.f"
    } else if (iascl == 2) {
#line 462 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 463 "sgelsy.f"
	slascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 465 "sgelsy.f"
    }
#line 466 "sgelsy.f"
    if (ibscl == 1) {
#line 467 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 468 "sgelsy.f"
    } else if (ibscl == 2) {
#line 469 "sgelsy.f"
	slascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
		 info, (ftnlen)1);
#line 470 "sgelsy.f"
    }

#line 472 "sgelsy.f"
L70:
#line 473 "sgelsy.f"
    work[1] = (doublereal) lwkopt;

#line 475 "sgelsy.f"
    return 0;

/*     End of SGELSY */

} /* sgelsy_ */

