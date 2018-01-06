#line 1 "stgsyl.f"
/* stgsyl.f -- translated by f2c (version 20100827).
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

#line 1 "stgsyl.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static doublereal c_b14 = 0.;
static integer c__1 = 1;
static doublereal c_b51 = -1.;
static doublereal c_b52 = 1.;

/* > \brief \b STGSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/*                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, */
/*      $                   LWORK, M, N */
/*       REAL               DIF, SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSYL solves the generalized Sylvester equation: */
/* > */
/* >             A * R - L * B = scale * C                 (1) */
/* >             D * R - L * E = scale * F */
/* > */
/* > where R and L are unknown m-by-n matrices, (A, D), (B, E) and */
/* > (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n, */
/* > respectively, with real entries. (A, D) and (B, E) must be in */
/* > generalized (real) Schur canonical form, i.e. A, B are upper quasi */
/* > triangular and D, E are upper triangular. */
/* > */
/* > The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output */
/* > scaling factor chosen to avoid overflow. */
/* > */
/* > In matrix notation (1) is equivalent to solve  Zx = scale b, where */
/* > Z is defined as */
/* > */
/* >            Z = [ kron(In, A)  -kron(B**T, Im) ]         (2) */
/* >                [ kron(In, D)  -kron(E**T, Im) ]. */
/* > */
/* > Here Ik is the identity matrix of size k and X**T is the transpose of */
/* > X. kron(X, Y) is the Kronecker product between the matrices X and Y. */
/* > */
/* > If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b, */
/* > which is equivalent to solve for R and L in */
/* > */
/* >             A**T * R + D**T * L = scale * C           (3) */
/* >             R * B**T + L * E**T = scale * -F */
/* > */
/* > This case (TRANS = 'T') is used to compute an one-norm-based estimate */
/* > of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D) */
/* > and (B,E), using SLACON. */
/* > */
/* > If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate */
/* > of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the */
/* > reciprocal of the smallest singular value of Z. See [1-2] for more */
/* > information. */
/* > */
/* > This is a level 3 BLAS algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N', solve the generalized Sylvester equation (1). */
/* >          = 'T', solve the 'transposed' system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          Specifies what kind of functionality to be performed. */
/* >           =0: solve (1) only. */
/* >           =1: The functionality of 0 and 3. */
/* >           =2: The functionality of 0 and 4. */
/* >           =3: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* >               (look ahead strategy IJOB  = 1 is used). */
/* >           =4: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* >               ( SGECON on sub-systems is used ). */
/* >          Not referenced if TRANS = 'T'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The order of the matrices A and D, and the row dimension of */
/* >          the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices B and E, and the column dimension */
/* >          of the matrices C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA, M) */
/* >          The upper quasi triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB, N) */
/* >          The upper quasi triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC, N) */
/* >          On entry, C contains the right-hand-side of the first matrix */
/* >          equation in (1) or (3). */
/* >          On exit, if IJOB = 0, 1 or 2, C has been overwritten by */
/* >          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R, */
/* >          the solution achieved during the computation of the */
/* >          Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (LDD, M) */
/* >          The upper triangular matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LDD */
/* > \verbatim */
/* >          LDD is INTEGER */
/* >          The leading dimension of the array D. LDD >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (LDE, N) */
/* >          The upper triangular matrix E. */
/* > \endverbatim */
/* > */
/* > \param[in] LDE */
/* > \verbatim */
/* >          LDE is INTEGER */
/* >          The leading dimension of the array E. LDE >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* >          F is REAL array, dimension (LDF, N) */
/* >          On entry, F contains the right-hand-side of the second matrix */
/* >          equation in (1) or (3). */
/* >          On exit, if IJOB = 0, 1 or 2, F has been overwritten by */
/* >          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L, */
/* >          the solution achieved during the computation of the */
/* >          Dif-estimate. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* >          LDF is INTEGER */
/* >          The leading dimension of the array F. LDF >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is REAL */
/* >          On exit DIF is the reciprocal of a lower bound of the */
/* >          reciprocal of the Dif-function, i.e. DIF is an upper bound of */
/* >          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2). */
/* >          IF IJOB = 0 or TRANS = 'T', DIF is not touched. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          On exit SCALE is the scaling factor in (1) or (3). */
/* >          If 0 < SCALE < 1, C and F hold the solutions R and L, resp., */
/* >          to a slightly perturbed system but the input matrices A, B, D */
/* >          and E have not been changed. If SCALE = 0, C and F hold the */
/* >          solutions R and L, respectively, to the homogeneous system */
/* >          with C = F = 0. Normally, SCALE = 1. */
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
/* >          The dimension of the array WORK. LWORK > = 1. */
/* >          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (M+N+6) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >            =0: successful exit */
/* >            <0: If INFO = -i, the i-th argument had an illegal value. */
/* >            >0: (A, D) and (B, E) have common or close eigenvalues. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* >  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* >      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22, */
/* >      No 1, 1996. */
/* > */
/* >  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester */
/* >      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal. */
/* >      Appl., 15(4):1045-1060, 1994 */
/* > */
/* >  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with */
/* >      Condition Estimators for Solving the Generalized Sylvester */
/* >      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7, */
/* >      July 1989, pp 745-751. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	iwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3, 
	    i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p, q, ie, je, mb, nb, is, js, pq;
    static doublereal dsum;
    static integer ppqq;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ifunc;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer linfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer lwmin;
    static doublereal scale2, dscale;
    extern /* Subroutine */ int stgsy2_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static integer iround;
    static logical notran;
    static integer isolve;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*  Replaced various illegal calls to SCOPY by calls to SLASET. */
/*  Sven Hammarling, 1/5/02. */

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

/*     Decode and test input parameters */

#line 349 "stgsyl.f"
    /* Parameter adjustments */
#line 349 "stgsyl.f"
    a_dim1 = *lda;
#line 349 "stgsyl.f"
    a_offset = 1 + a_dim1;
#line 349 "stgsyl.f"
    a -= a_offset;
#line 349 "stgsyl.f"
    b_dim1 = *ldb;
#line 349 "stgsyl.f"
    b_offset = 1 + b_dim1;
#line 349 "stgsyl.f"
    b -= b_offset;
#line 349 "stgsyl.f"
    c_dim1 = *ldc;
#line 349 "stgsyl.f"
    c_offset = 1 + c_dim1;
#line 349 "stgsyl.f"
    c__ -= c_offset;
#line 349 "stgsyl.f"
    d_dim1 = *ldd;
#line 349 "stgsyl.f"
    d_offset = 1 + d_dim1;
#line 349 "stgsyl.f"
    d__ -= d_offset;
#line 349 "stgsyl.f"
    e_dim1 = *lde;
#line 349 "stgsyl.f"
    e_offset = 1 + e_dim1;
#line 349 "stgsyl.f"
    e -= e_offset;
#line 349 "stgsyl.f"
    f_dim1 = *ldf;
#line 349 "stgsyl.f"
    f_offset = 1 + f_dim1;
#line 349 "stgsyl.f"
    f -= f_offset;
#line 349 "stgsyl.f"
    --work;
#line 349 "stgsyl.f"
    --iwork;
#line 349 "stgsyl.f"

#line 349 "stgsyl.f"
    /* Function Body */
#line 349 "stgsyl.f"
    *info = 0;
#line 350 "stgsyl.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 351 "stgsyl.f"
    lquery = *lwork == -1;

#line 353 "stgsyl.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 354 "stgsyl.f"
	*info = -1;
#line 355 "stgsyl.f"
    } else if (notran) {
#line 356 "stgsyl.f"
	if (*ijob < 0 || *ijob > 4) {
#line 357 "stgsyl.f"
	    *info = -2;
#line 358 "stgsyl.f"
	}
#line 359 "stgsyl.f"
    }
#line 360 "stgsyl.f"
    if (*info == 0) {
#line 361 "stgsyl.f"
	if (*m <= 0) {
#line 362 "stgsyl.f"
	    *info = -3;
#line 363 "stgsyl.f"
	} else if (*n <= 0) {
#line 364 "stgsyl.f"
	    *info = -4;
#line 365 "stgsyl.f"
	} else if (*lda < max(1,*m)) {
#line 366 "stgsyl.f"
	    *info = -6;
#line 367 "stgsyl.f"
	} else if (*ldb < max(1,*n)) {
#line 368 "stgsyl.f"
	    *info = -8;
#line 369 "stgsyl.f"
	} else if (*ldc < max(1,*m)) {
#line 370 "stgsyl.f"
	    *info = -10;
#line 371 "stgsyl.f"
	} else if (*ldd < max(1,*m)) {
#line 372 "stgsyl.f"
	    *info = -12;
#line 373 "stgsyl.f"
	} else if (*lde < max(1,*n)) {
#line 374 "stgsyl.f"
	    *info = -14;
#line 375 "stgsyl.f"
	} else if (*ldf < max(1,*m)) {
#line 376 "stgsyl.f"
	    *info = -16;
#line 377 "stgsyl.f"
	}
#line 378 "stgsyl.f"
    }

#line 380 "stgsyl.f"
    if (*info == 0) {
#line 381 "stgsyl.f"
	if (notran) {
#line 382 "stgsyl.f"
	    if (*ijob == 1 || *ijob == 2) {
/* Computing MAX */
#line 383 "stgsyl.f"
		i__1 = 1, i__2 = (*m << 1) * *n;
#line 383 "stgsyl.f"
		lwmin = max(i__1,i__2);
#line 384 "stgsyl.f"
	    } else {
#line 385 "stgsyl.f"
		lwmin = 1;
#line 386 "stgsyl.f"
	    }
#line 387 "stgsyl.f"
	} else {
#line 388 "stgsyl.f"
	    lwmin = 1;
#line 389 "stgsyl.f"
	}
#line 390 "stgsyl.f"
	work[1] = (doublereal) lwmin;

#line 392 "stgsyl.f"
	if (*lwork < lwmin && ! lquery) {
#line 393 "stgsyl.f"
	    *info = -20;
#line 394 "stgsyl.f"
	}
#line 395 "stgsyl.f"
    }

#line 397 "stgsyl.f"
    if (*info != 0) {
#line 398 "stgsyl.f"
	i__1 = -(*info);
#line 398 "stgsyl.f"
	xerbla_("STGSYL", &i__1, (ftnlen)6);
#line 399 "stgsyl.f"
	return 0;
#line 400 "stgsyl.f"
    } else if (lquery) {
#line 401 "stgsyl.f"
	return 0;
#line 402 "stgsyl.f"
    }

/*     Quick return if possible */

#line 406 "stgsyl.f"
    if (*m == 0 || *n == 0) {
#line 407 "stgsyl.f"
	*scale = 1.;
#line 408 "stgsyl.f"
	if (notran) {
#line 409 "stgsyl.f"
	    if (*ijob != 0) {
#line 410 "stgsyl.f"
		*dif = 0.;
#line 411 "stgsyl.f"
	    }
#line 412 "stgsyl.f"
	}
#line 413 "stgsyl.f"
	return 0;
#line 414 "stgsyl.f"
    }

/*     Determine optimal block sizes MB and NB */

#line 418 "stgsyl.f"
    mb = ilaenv_(&c__2, "STGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 419 "stgsyl.f"
    nb = ilaenv_(&c__5, "STGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 421 "stgsyl.f"
    isolve = 1;
#line 422 "stgsyl.f"
    ifunc = 0;
#line 423 "stgsyl.f"
    if (notran) {
#line 424 "stgsyl.f"
	if (*ijob >= 3) {
#line 425 "stgsyl.f"
	    ifunc = *ijob - 2;
#line 426 "stgsyl.f"
	    slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (ftnlen)1)
		    ;
#line 427 "stgsyl.f"
	    slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)1);
#line 428 "stgsyl.f"
	} else if (*ijob >= 1 && notran) {
#line 429 "stgsyl.f"
	    isolve = 2;
#line 430 "stgsyl.f"
	}
#line 431 "stgsyl.f"
    }

#line 433 "stgsyl.f"
    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n) {

#line 436 "stgsyl.f"
	i__1 = isolve;
#line 436 "stgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Use unblocked Level 2 solver */

#line 440 "stgsyl.f"
	    dscale = 0.;
#line 441 "stgsyl.f"
	    dsum = 1.;
#line 442 "stgsyl.f"
	    pq = 0;
#line 443 "stgsyl.f"
	    stgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb,
		     &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], 
		    lde, &f[f_offset], ldf, scale, &dsum, &dscale, &iwork[1], 
		    &pq, info, (ftnlen)1);
#line 446 "stgsyl.f"
	    if (dscale != 0.) {
#line 447 "stgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 448 "stgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 449 "stgsyl.f"
		} else {
#line 450 "stgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 451 "stgsyl.f"
		}
#line 452 "stgsyl.f"
	    }

#line 454 "stgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 455 "stgsyl.f"
		if (notran) {
#line 456 "stgsyl.f"
		    ifunc = *ijob;
#line 457 "stgsyl.f"
		}
#line 458 "stgsyl.f"
		scale2 = *scale;
#line 459 "stgsyl.f"
		slacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 460 "stgsyl.f"
		slacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 461 "stgsyl.f"
		slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (
			ftnlen)1);
#line 462 "stgsyl.f"
		slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)
			1);
#line 463 "stgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 464 "stgsyl.f"
		slacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 465 "stgsyl.f"
		slacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 466 "stgsyl.f"
		*scale = scale2;
#line 467 "stgsyl.f"
	    }
#line 468 "stgsyl.f"
/* L30: */
#line 468 "stgsyl.f"
	}

#line 470 "stgsyl.f"
	return 0;
#line 471 "stgsyl.f"
    }

/*     Determine block structure of A */

#line 475 "stgsyl.f"
    p = 0;
#line 476 "stgsyl.f"
    i__ = 1;
#line 477 "stgsyl.f"
L40:
#line 478 "stgsyl.f"
    if (i__ > *m) {
#line 478 "stgsyl.f"
	goto L50;
#line 478 "stgsyl.f"
    }
#line 480 "stgsyl.f"
    ++p;
#line 481 "stgsyl.f"
    iwork[p] = i__;
#line 482 "stgsyl.f"
    i__ += mb;
#line 483 "stgsyl.f"
    if (i__ >= *m) {
#line 483 "stgsyl.f"
	goto L50;
#line 483 "stgsyl.f"
    }
#line 485 "stgsyl.f"
    if (a[i__ + (i__ - 1) * a_dim1] != 0.) {
#line 485 "stgsyl.f"
	++i__;
#line 485 "stgsyl.f"
    }
#line 487 "stgsyl.f"
    goto L40;
#line 488 "stgsyl.f"
L50:

#line 490 "stgsyl.f"
    iwork[p + 1] = *m + 1;
#line 491 "stgsyl.f"
    if (iwork[p] == iwork[p + 1]) {
#line 491 "stgsyl.f"
	--p;
#line 491 "stgsyl.f"
    }

/*     Determine block structure of B */

#line 496 "stgsyl.f"
    q = p + 1;
#line 497 "stgsyl.f"
    j = 1;
#line 498 "stgsyl.f"
L60:
#line 499 "stgsyl.f"
    if (j > *n) {
#line 499 "stgsyl.f"
	goto L70;
#line 499 "stgsyl.f"
    }
#line 501 "stgsyl.f"
    ++q;
#line 502 "stgsyl.f"
    iwork[q] = j;
#line 503 "stgsyl.f"
    j += nb;
#line 504 "stgsyl.f"
    if (j >= *n) {
#line 504 "stgsyl.f"
	goto L70;
#line 504 "stgsyl.f"
    }
#line 506 "stgsyl.f"
    if (b[j + (j - 1) * b_dim1] != 0.) {
#line 506 "stgsyl.f"
	++j;
#line 506 "stgsyl.f"
    }
#line 508 "stgsyl.f"
    goto L60;
#line 509 "stgsyl.f"
L70:

#line 511 "stgsyl.f"
    iwork[q + 1] = *n + 1;
#line 512 "stgsyl.f"
    if (iwork[q] == iwork[q + 1]) {
#line 512 "stgsyl.f"
	--q;
#line 512 "stgsyl.f"
    }

#line 515 "stgsyl.f"
    if (notran) {

#line 517 "stgsyl.f"
	i__1 = isolve;
#line 517 "stgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Solve (I, J)-subsystem */
/*               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*           for I = P, P - 1,..., 1; J = 1, 2,..., Q */

#line 524 "stgsyl.f"
	    dscale = 0.;
#line 525 "stgsyl.f"
	    dsum = 1.;
#line 526 "stgsyl.f"
	    pq = 0;
#line 527 "stgsyl.f"
	    *scale = 1.;
#line 528 "stgsyl.f"
	    i__2 = q;
#line 528 "stgsyl.f"
	    for (j = p + 2; j <= i__2; ++j) {
#line 529 "stgsyl.f"
		js = iwork[j];
#line 530 "stgsyl.f"
		je = iwork[j + 1] - 1;
#line 531 "stgsyl.f"
		nb = je - js + 1;
#line 532 "stgsyl.f"
		for (i__ = p; i__ >= 1; --i__) {
#line 533 "stgsyl.f"
		    is = iwork[i__];
#line 534 "stgsyl.f"
		    ie = iwork[i__ + 1] - 1;
#line 535 "stgsyl.f"
		    mb = ie - is + 1;
#line 536 "stgsyl.f"
		    ppqq = 0;
#line 537 "stgsyl.f"
		    stgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], 
			    lda, &b[js + js * b_dim1], ldb, &c__[is + js * 
			    c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js 
			    + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &
			    scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, &
			    linfo, (ftnlen)1);
#line 542 "stgsyl.f"
		    if (linfo > 0) {
#line 542 "stgsyl.f"
			*info = linfo;
#line 542 "stgsyl.f"
		    }

#line 545 "stgsyl.f"
		    pq += ppqq;
#line 546 "stgsyl.f"
		    if (scaloc != 1.) {
#line 547 "stgsyl.f"
			i__3 = js - 1;
#line 547 "stgsyl.f"
			for (k = 1; k <= i__3; ++k) {
#line 548 "stgsyl.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 549 "stgsyl.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 550 "stgsyl.f"
/* L80: */
#line 550 "stgsyl.f"
			}
#line 551 "stgsyl.f"
			i__3 = je;
#line 551 "stgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 552 "stgsyl.f"
			    i__4 = is - 1;
#line 552 "stgsyl.f"
			    sscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], &
				    c__1);
#line 553 "stgsyl.f"
			    i__4 = is - 1;
#line 553 "stgsyl.f"
			    sscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 554 "stgsyl.f"
/* L90: */
#line 554 "stgsyl.f"
			}
#line 555 "stgsyl.f"
			i__3 = je;
#line 555 "stgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 556 "stgsyl.f"
			    i__4 = *m - ie;
#line 556 "stgsyl.f"
			    sscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], 
				    &c__1);
#line 557 "stgsyl.f"
			    i__4 = *m - ie;
#line 557 "stgsyl.f"
			    sscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], &
				    c__1);
#line 558 "stgsyl.f"
/* L100: */
#line 558 "stgsyl.f"
			}
#line 559 "stgsyl.f"
			i__3 = *n;
#line 559 "stgsyl.f"
			for (k = je + 1; k <= i__3; ++k) {
#line 560 "stgsyl.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 561 "stgsyl.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 562 "stgsyl.f"
/* L110: */
#line 562 "stgsyl.f"
			}
#line 563 "stgsyl.f"
			*scale *= scaloc;
#line 564 "stgsyl.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 569 "stgsyl.f"
		    if (i__ > 1) {
#line 570 "stgsyl.f"
			i__3 = is - 1;
#line 570 "stgsyl.f"
			sgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &a[is * 
				a_dim1 + 1], lda, &c__[is + js * c_dim1], ldc,
				 &c_b52, &c__[js * c_dim1 + 1], ldc, (ftnlen)
				1, (ftnlen)1);
#line 573 "stgsyl.f"
			i__3 = is - 1;
#line 573 "stgsyl.f"
			sgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &d__[is * 
				d_dim1 + 1], ldd, &c__[is + js * c_dim1], ldc,
				 &c_b52, &f[js * f_dim1 + 1], ldf, (ftnlen)1, 
				(ftnlen)1);
#line 576 "stgsyl.f"
		    }
#line 577 "stgsyl.f"
		    if (j < q) {
#line 578 "stgsyl.f"
			i__3 = *n - je;
#line 578 "stgsyl.f"
			sgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js *
				 f_dim1], ldf, &b[js + (je + 1) * b_dim1], 
				ldb, &c_b52, &c__[is + (je + 1) * c_dim1], 
				ldc, (ftnlen)1, (ftnlen)1);
#line 581 "stgsyl.f"
			i__3 = *n - je;
#line 581 "stgsyl.f"
			sgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js *
				 f_dim1], ldf, &e[js + (je + 1) * e_dim1], 
				lde, &c_b52, &f[is + (je + 1) * f_dim1], ldf, 
				(ftnlen)1, (ftnlen)1);
#line 584 "stgsyl.f"
		    }
#line 585 "stgsyl.f"
/* L120: */
#line 585 "stgsyl.f"
		}
#line 586 "stgsyl.f"
/* L130: */
#line 586 "stgsyl.f"
	    }
#line 587 "stgsyl.f"
	    if (dscale != 0.) {
#line 588 "stgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 589 "stgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 590 "stgsyl.f"
		} else {
#line 591 "stgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 592 "stgsyl.f"
		}
#line 593 "stgsyl.f"
	    }
#line 594 "stgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 595 "stgsyl.f"
		if (notran) {
#line 596 "stgsyl.f"
		    ifunc = *ijob;
#line 597 "stgsyl.f"
		}
#line 598 "stgsyl.f"
		scale2 = *scale;
#line 599 "stgsyl.f"
		slacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 600 "stgsyl.f"
		slacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 601 "stgsyl.f"
		slaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (
			ftnlen)1);
#line 602 "stgsyl.f"
		slaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)
			1);
#line 603 "stgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 604 "stgsyl.f"
		slacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 605 "stgsyl.f"
		slacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 606 "stgsyl.f"
		*scale = scale2;
#line 607 "stgsyl.f"
	    }
#line 608 "stgsyl.f"
/* L150: */
#line 608 "stgsyl.f"
	}

#line 610 "stgsyl.f"
    } else {

/*        Solve transposed (I, J)-subsystem */
/*             A(I, I)**T * R(I, J)  + D(I, I)**T * L(I, J)  =  C(I, J) */
/*             R(I, J)  * B(J, J)**T + L(I, J)  * E(J, J)**T = -F(I, J) */
/*        for I = 1,2,..., P; J = Q, Q-1,..., 1 */

#line 617 "stgsyl.f"
	*scale = 1.;
#line 618 "stgsyl.f"
	i__1 = p;
#line 618 "stgsyl.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 619 "stgsyl.f"
	    is = iwork[i__];
#line 620 "stgsyl.f"
	    ie = iwork[i__ + 1] - 1;
#line 621 "stgsyl.f"
	    mb = ie - is + 1;
#line 622 "stgsyl.f"
	    i__2 = p + 2;
#line 622 "stgsyl.f"
	    for (j = q; j >= i__2; --j) {
#line 623 "stgsyl.f"
		js = iwork[j];
#line 624 "stgsyl.f"
		je = iwork[j + 1] - 1;
#line 625 "stgsyl.f"
		nb = je - js + 1;
#line 626 "stgsyl.f"
		stgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &
			b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc,
			 &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], 
			lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum, &
			dscale, &iwork[q + 2], &ppqq, &linfo, (ftnlen)1);
#line 631 "stgsyl.f"
		if (linfo > 0) {
#line 631 "stgsyl.f"
		    *info = linfo;
#line 631 "stgsyl.f"
		}
#line 633 "stgsyl.f"
		if (scaloc != 1.) {
#line 634 "stgsyl.f"
		    i__3 = js - 1;
#line 634 "stgsyl.f"
		    for (k = 1; k <= i__3; ++k) {
#line 635 "stgsyl.f"
			sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 636 "stgsyl.f"
			sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 637 "stgsyl.f"
/* L160: */
#line 637 "stgsyl.f"
		    }
#line 638 "stgsyl.f"
		    i__3 = je;
#line 638 "stgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 639 "stgsyl.f"
			i__4 = is - 1;
#line 639 "stgsyl.f"
			sscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 640 "stgsyl.f"
			i__4 = is - 1;
#line 640 "stgsyl.f"
			sscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 641 "stgsyl.f"
/* L170: */
#line 641 "stgsyl.f"
		    }
#line 642 "stgsyl.f"
		    i__3 = je;
#line 642 "stgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 643 "stgsyl.f"
			i__4 = *m - ie;
#line 643 "stgsyl.f"
			sscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], &
				c__1);
#line 644 "stgsyl.f"
			i__4 = *m - ie;
#line 644 "stgsyl.f"
			sscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], &c__1)
				;
#line 645 "stgsyl.f"
/* L180: */
#line 645 "stgsyl.f"
		    }
#line 646 "stgsyl.f"
		    i__3 = *n;
#line 646 "stgsyl.f"
		    for (k = je + 1; k <= i__3; ++k) {
#line 647 "stgsyl.f"
			sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 648 "stgsyl.f"
			sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 649 "stgsyl.f"
/* L190: */
#line 649 "stgsyl.f"
		    }
#line 650 "stgsyl.f"
		    *scale *= scaloc;
#line 651 "stgsyl.f"
		}

/*              Substitute R(I, J) and L(I, J) into remaining equation. */

#line 655 "stgsyl.f"
		if (j > p + 2) {
#line 656 "stgsyl.f"
		    i__3 = js - 1;
#line 656 "stgsyl.f"
		    sgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &c__[is + js * 
			    c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &c_b52, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 659 "stgsyl.f"
		    i__3 = js - 1;
#line 659 "stgsyl.f"
		    sgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &f[is + js * 
			    f_dim1], ldf, &e[js * e_dim1 + 1], lde, &c_b52, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 662 "stgsyl.f"
		}
#line 663 "stgsyl.f"
		if (i__ < p) {
#line 664 "stgsyl.f"
		    i__3 = *m - ie;
#line 664 "stgsyl.f"
		    sgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &a[is + (ie + 1)
			     * a_dim1], lda, &c__[is + js * c_dim1], ldc, &
			    c_b52, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 667 "stgsyl.f"
		    i__3 = *m - ie;
#line 667 "stgsyl.f"
		    sgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &d__[is + (ie + 
			    1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, &
			    c_b52, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 670 "stgsyl.f"
		}
#line 671 "stgsyl.f"
/* L200: */
#line 671 "stgsyl.f"
	    }
#line 672 "stgsyl.f"
/* L210: */
#line 672 "stgsyl.f"
	}

#line 674 "stgsyl.f"
    }

#line 676 "stgsyl.f"
    work[1] = (doublereal) lwmin;

#line 678 "stgsyl.f"
    return 0;

/*     End of STGSYL */

} /* stgsyl_ */

