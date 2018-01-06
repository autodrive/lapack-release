#line 1 "dtgsyl.f"
/* dtgsyl.f -- translated by f2c (version 20100827).
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

#line 1 "dtgsyl.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static doublereal c_b14 = 0.;
static integer c__1 = 1;
static doublereal c_b51 = -1.;
static doublereal c_b52 = 1.;

/* > \brief \b DTGSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/*                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, */
/*      $                   LWORK, M, N */
/*       DOUBLE PRECISION   DIF, SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGSYL solves the generalized Sylvester equation: */
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
/* > If TRANS = 'T', DTGSYL solves the transposed system Z**T*y = scale*b, */
/* > which is equivalent to solve for R and L in */
/* > */
/* >             A**T * R + D**T * L = scale * C           (3) */
/* >             R * B**T + L * E**T = scale * -F */
/* > */
/* > This case (TRANS = 'T') is used to compute an one-norm-based estimate */
/* > of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D) */
/* > and (B,E), using DLACON. */
/* > */
/* > If IJOB >= 1, DTGSYL computes a Frobenius norm-based estimate */
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
/* >               ( DGECON on sub-systems is used ). */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, M) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
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
/* >          C is DOUBLE PRECISION array, dimension (LDC, N) */
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
/* >          D is DOUBLE PRECISION array, dimension (LDD, M) */
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
/* >          E is DOUBLE PRECISION array, dimension (LDE, N) */
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
/* >          F is DOUBLE PRECISION array, dimension (LDF, N) */
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
/* >          DIF is DOUBLE PRECISION */
/* >          On exit DIF is the reciprocal of a lower bound of the */
/* >          reciprocal of the Dif-function, i.e. DIF is an upper bound of */
/* >          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2). */
/* >          IF IJOB = 0 or TRANS = 'T', DIF is not touched. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup doubleSYcomputational */

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
/* Subroutine */ int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ifunc, linfo, lwmin;
    static doublereal scale2;
    extern /* Subroutine */ int dtgsy2_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static doublereal dscale, scaloc;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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
/*  Replaced various illegal calls to DCOPY by calls to DLASET. */
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

#line 349 "dtgsyl.f"
    /* Parameter adjustments */
#line 349 "dtgsyl.f"
    a_dim1 = *lda;
#line 349 "dtgsyl.f"
    a_offset = 1 + a_dim1;
#line 349 "dtgsyl.f"
    a -= a_offset;
#line 349 "dtgsyl.f"
    b_dim1 = *ldb;
#line 349 "dtgsyl.f"
    b_offset = 1 + b_dim1;
#line 349 "dtgsyl.f"
    b -= b_offset;
#line 349 "dtgsyl.f"
    c_dim1 = *ldc;
#line 349 "dtgsyl.f"
    c_offset = 1 + c_dim1;
#line 349 "dtgsyl.f"
    c__ -= c_offset;
#line 349 "dtgsyl.f"
    d_dim1 = *ldd;
#line 349 "dtgsyl.f"
    d_offset = 1 + d_dim1;
#line 349 "dtgsyl.f"
    d__ -= d_offset;
#line 349 "dtgsyl.f"
    e_dim1 = *lde;
#line 349 "dtgsyl.f"
    e_offset = 1 + e_dim1;
#line 349 "dtgsyl.f"
    e -= e_offset;
#line 349 "dtgsyl.f"
    f_dim1 = *ldf;
#line 349 "dtgsyl.f"
    f_offset = 1 + f_dim1;
#line 349 "dtgsyl.f"
    f -= f_offset;
#line 349 "dtgsyl.f"
    --work;
#line 349 "dtgsyl.f"
    --iwork;
#line 349 "dtgsyl.f"

#line 349 "dtgsyl.f"
    /* Function Body */
#line 349 "dtgsyl.f"
    *info = 0;
#line 350 "dtgsyl.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 351 "dtgsyl.f"
    lquery = *lwork == -1;

#line 353 "dtgsyl.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 354 "dtgsyl.f"
	*info = -1;
#line 355 "dtgsyl.f"
    } else if (notran) {
#line 356 "dtgsyl.f"
	if (*ijob < 0 || *ijob > 4) {
#line 357 "dtgsyl.f"
	    *info = -2;
#line 358 "dtgsyl.f"
	}
#line 359 "dtgsyl.f"
    }
#line 360 "dtgsyl.f"
    if (*info == 0) {
#line 361 "dtgsyl.f"
	if (*m <= 0) {
#line 362 "dtgsyl.f"
	    *info = -3;
#line 363 "dtgsyl.f"
	} else if (*n <= 0) {
#line 364 "dtgsyl.f"
	    *info = -4;
#line 365 "dtgsyl.f"
	} else if (*lda < max(1,*m)) {
#line 366 "dtgsyl.f"
	    *info = -6;
#line 367 "dtgsyl.f"
	} else if (*ldb < max(1,*n)) {
#line 368 "dtgsyl.f"
	    *info = -8;
#line 369 "dtgsyl.f"
	} else if (*ldc < max(1,*m)) {
#line 370 "dtgsyl.f"
	    *info = -10;
#line 371 "dtgsyl.f"
	} else if (*ldd < max(1,*m)) {
#line 372 "dtgsyl.f"
	    *info = -12;
#line 373 "dtgsyl.f"
	} else if (*lde < max(1,*n)) {
#line 374 "dtgsyl.f"
	    *info = -14;
#line 375 "dtgsyl.f"
	} else if (*ldf < max(1,*m)) {
#line 376 "dtgsyl.f"
	    *info = -16;
#line 377 "dtgsyl.f"
	}
#line 378 "dtgsyl.f"
    }

#line 380 "dtgsyl.f"
    if (*info == 0) {
#line 381 "dtgsyl.f"
	if (notran) {
#line 382 "dtgsyl.f"
	    if (*ijob == 1 || *ijob == 2) {
/* Computing MAX */
#line 383 "dtgsyl.f"
		i__1 = 1, i__2 = (*m << 1) * *n;
#line 383 "dtgsyl.f"
		lwmin = max(i__1,i__2);
#line 384 "dtgsyl.f"
	    } else {
#line 385 "dtgsyl.f"
		lwmin = 1;
#line 386 "dtgsyl.f"
	    }
#line 387 "dtgsyl.f"
	} else {
#line 388 "dtgsyl.f"
	    lwmin = 1;
#line 389 "dtgsyl.f"
	}
#line 390 "dtgsyl.f"
	work[1] = (doublereal) lwmin;

#line 392 "dtgsyl.f"
	if (*lwork < lwmin && ! lquery) {
#line 393 "dtgsyl.f"
	    *info = -20;
#line 394 "dtgsyl.f"
	}
#line 395 "dtgsyl.f"
    }

#line 397 "dtgsyl.f"
    if (*info != 0) {
#line 398 "dtgsyl.f"
	i__1 = -(*info);
#line 398 "dtgsyl.f"
	xerbla_("DTGSYL", &i__1, (ftnlen)6);
#line 399 "dtgsyl.f"
	return 0;
#line 400 "dtgsyl.f"
    } else if (lquery) {
#line 401 "dtgsyl.f"
	return 0;
#line 402 "dtgsyl.f"
    }

/*     Quick return if possible */

#line 406 "dtgsyl.f"
    if (*m == 0 || *n == 0) {
#line 407 "dtgsyl.f"
	*scale = 1.;
#line 408 "dtgsyl.f"
	if (notran) {
#line 409 "dtgsyl.f"
	    if (*ijob != 0) {
#line 410 "dtgsyl.f"
		*dif = 0.;
#line 411 "dtgsyl.f"
	    }
#line 412 "dtgsyl.f"
	}
#line 413 "dtgsyl.f"
	return 0;
#line 414 "dtgsyl.f"
    }

/*     Determine optimal block sizes MB and NB */

#line 418 "dtgsyl.f"
    mb = ilaenv_(&c__2, "DTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 419 "dtgsyl.f"
    nb = ilaenv_(&c__5, "DTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 421 "dtgsyl.f"
    isolve = 1;
#line 422 "dtgsyl.f"
    ifunc = 0;
#line 423 "dtgsyl.f"
    if (notran) {
#line 424 "dtgsyl.f"
	if (*ijob >= 3) {
#line 425 "dtgsyl.f"
	    ifunc = *ijob - 2;
#line 426 "dtgsyl.f"
	    dlaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (ftnlen)1)
		    ;
#line 427 "dtgsyl.f"
	    dlaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)1);
#line 428 "dtgsyl.f"
	} else if (*ijob >= 1) {
#line 429 "dtgsyl.f"
	    isolve = 2;
#line 430 "dtgsyl.f"
	}
#line 431 "dtgsyl.f"
    }

#line 433 "dtgsyl.f"
    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n) {

#line 436 "dtgsyl.f"
	i__1 = isolve;
#line 436 "dtgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Use unblocked Level 2 solver */

#line 440 "dtgsyl.f"
	    dscale = 0.;
#line 441 "dtgsyl.f"
	    dsum = 1.;
#line 442 "dtgsyl.f"
	    pq = 0;
#line 443 "dtgsyl.f"
	    dtgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb,
		     &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], 
		    lde, &f[f_offset], ldf, scale, &dsum, &dscale, &iwork[1], 
		    &pq, info, (ftnlen)1);
#line 446 "dtgsyl.f"
	    if (dscale != 0.) {
#line 447 "dtgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 448 "dtgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 449 "dtgsyl.f"
		} else {
#line 450 "dtgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 451 "dtgsyl.f"
		}
#line 452 "dtgsyl.f"
	    }

#line 454 "dtgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 455 "dtgsyl.f"
		if (notran) {
#line 456 "dtgsyl.f"
		    ifunc = *ijob;
#line 457 "dtgsyl.f"
		}
#line 458 "dtgsyl.f"
		scale2 = *scale;
#line 459 "dtgsyl.f"
		dlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 460 "dtgsyl.f"
		dlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 461 "dtgsyl.f"
		dlaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (
			ftnlen)1);
#line 462 "dtgsyl.f"
		dlaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)
			1);
#line 463 "dtgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 464 "dtgsyl.f"
		dlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 465 "dtgsyl.f"
		dlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 466 "dtgsyl.f"
		*scale = scale2;
#line 467 "dtgsyl.f"
	    }
#line 468 "dtgsyl.f"
/* L30: */
#line 468 "dtgsyl.f"
	}

#line 470 "dtgsyl.f"
	return 0;
#line 471 "dtgsyl.f"
    }

/*     Determine block structure of A */

#line 475 "dtgsyl.f"
    p = 0;
#line 476 "dtgsyl.f"
    i__ = 1;
#line 477 "dtgsyl.f"
L40:
#line 478 "dtgsyl.f"
    if (i__ > *m) {
#line 478 "dtgsyl.f"
	goto L50;
#line 478 "dtgsyl.f"
    }
#line 480 "dtgsyl.f"
    ++p;
#line 481 "dtgsyl.f"
    iwork[p] = i__;
#line 482 "dtgsyl.f"
    i__ += mb;
#line 483 "dtgsyl.f"
    if (i__ >= *m) {
#line 483 "dtgsyl.f"
	goto L50;
#line 483 "dtgsyl.f"
    }
#line 485 "dtgsyl.f"
    if (a[i__ + (i__ - 1) * a_dim1] != 0.) {
#line 485 "dtgsyl.f"
	++i__;
#line 485 "dtgsyl.f"
    }
#line 487 "dtgsyl.f"
    goto L40;
#line 488 "dtgsyl.f"
L50:

#line 490 "dtgsyl.f"
    iwork[p + 1] = *m + 1;
#line 491 "dtgsyl.f"
    if (iwork[p] == iwork[p + 1]) {
#line 491 "dtgsyl.f"
	--p;
#line 491 "dtgsyl.f"
    }

/*     Determine block structure of B */

#line 496 "dtgsyl.f"
    q = p + 1;
#line 497 "dtgsyl.f"
    j = 1;
#line 498 "dtgsyl.f"
L60:
#line 499 "dtgsyl.f"
    if (j > *n) {
#line 499 "dtgsyl.f"
	goto L70;
#line 499 "dtgsyl.f"
    }
#line 501 "dtgsyl.f"
    ++q;
#line 502 "dtgsyl.f"
    iwork[q] = j;
#line 503 "dtgsyl.f"
    j += nb;
#line 504 "dtgsyl.f"
    if (j >= *n) {
#line 504 "dtgsyl.f"
	goto L70;
#line 504 "dtgsyl.f"
    }
#line 506 "dtgsyl.f"
    if (b[j + (j - 1) * b_dim1] != 0.) {
#line 506 "dtgsyl.f"
	++j;
#line 506 "dtgsyl.f"
    }
#line 508 "dtgsyl.f"
    goto L60;
#line 509 "dtgsyl.f"
L70:

#line 511 "dtgsyl.f"
    iwork[q + 1] = *n + 1;
#line 512 "dtgsyl.f"
    if (iwork[q] == iwork[q + 1]) {
#line 512 "dtgsyl.f"
	--q;
#line 512 "dtgsyl.f"
    }

#line 515 "dtgsyl.f"
    if (notran) {

#line 517 "dtgsyl.f"
	i__1 = isolve;
#line 517 "dtgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Solve (I, J)-subsystem */
/*               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*           for I = P, P - 1,..., 1; J = 1, 2,..., Q */

#line 524 "dtgsyl.f"
	    dscale = 0.;
#line 525 "dtgsyl.f"
	    dsum = 1.;
#line 526 "dtgsyl.f"
	    pq = 0;
#line 527 "dtgsyl.f"
	    *scale = 1.;
#line 528 "dtgsyl.f"
	    i__2 = q;
#line 528 "dtgsyl.f"
	    for (j = p + 2; j <= i__2; ++j) {
#line 529 "dtgsyl.f"
		js = iwork[j];
#line 530 "dtgsyl.f"
		je = iwork[j + 1] - 1;
#line 531 "dtgsyl.f"
		nb = je - js + 1;
#line 532 "dtgsyl.f"
		for (i__ = p; i__ >= 1; --i__) {
#line 533 "dtgsyl.f"
		    is = iwork[i__];
#line 534 "dtgsyl.f"
		    ie = iwork[i__ + 1] - 1;
#line 535 "dtgsyl.f"
		    mb = ie - is + 1;
#line 536 "dtgsyl.f"
		    ppqq = 0;
#line 537 "dtgsyl.f"
		    dtgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], 
			    lda, &b[js + js * b_dim1], ldb, &c__[is + js * 
			    c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js 
			    + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &
			    scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, &
			    linfo, (ftnlen)1);
#line 542 "dtgsyl.f"
		    if (linfo > 0) {
#line 542 "dtgsyl.f"
			*info = linfo;
#line 542 "dtgsyl.f"
		    }

#line 545 "dtgsyl.f"
		    pq += ppqq;
#line 546 "dtgsyl.f"
		    if (scaloc != 1.) {
#line 547 "dtgsyl.f"
			i__3 = js - 1;
#line 547 "dtgsyl.f"
			for (k = 1; k <= i__3; ++k) {
#line 548 "dtgsyl.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 549 "dtgsyl.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 550 "dtgsyl.f"
/* L80: */
#line 550 "dtgsyl.f"
			}
#line 551 "dtgsyl.f"
			i__3 = je;
#line 551 "dtgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 552 "dtgsyl.f"
			    i__4 = is - 1;
#line 552 "dtgsyl.f"
			    dscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], &
				    c__1);
#line 553 "dtgsyl.f"
			    i__4 = is - 1;
#line 553 "dtgsyl.f"
			    dscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 554 "dtgsyl.f"
/* L90: */
#line 554 "dtgsyl.f"
			}
#line 555 "dtgsyl.f"
			i__3 = je;
#line 555 "dtgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 556 "dtgsyl.f"
			    i__4 = *m - ie;
#line 556 "dtgsyl.f"
			    dscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], 
				    &c__1);
#line 557 "dtgsyl.f"
			    i__4 = *m - ie;
#line 557 "dtgsyl.f"
			    dscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], &
				    c__1);
#line 558 "dtgsyl.f"
/* L100: */
#line 558 "dtgsyl.f"
			}
#line 559 "dtgsyl.f"
			i__3 = *n;
#line 559 "dtgsyl.f"
			for (k = je + 1; k <= i__3; ++k) {
#line 560 "dtgsyl.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 561 "dtgsyl.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 562 "dtgsyl.f"
/* L110: */
#line 562 "dtgsyl.f"
			}
#line 563 "dtgsyl.f"
			*scale *= scaloc;
#line 564 "dtgsyl.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 569 "dtgsyl.f"
		    if (i__ > 1) {
#line 570 "dtgsyl.f"
			i__3 = is - 1;
#line 570 "dtgsyl.f"
			dgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &a[is * 
				a_dim1 + 1], lda, &c__[is + js * c_dim1], ldc,
				 &c_b52, &c__[js * c_dim1 + 1], ldc, (ftnlen)
				1, (ftnlen)1);
#line 573 "dtgsyl.f"
			i__3 = is - 1;
#line 573 "dtgsyl.f"
			dgemm_("N", "N", &i__3, &nb, &mb, &c_b51, &d__[is * 
				d_dim1 + 1], ldd, &c__[is + js * c_dim1], ldc,
				 &c_b52, &f[js * f_dim1 + 1], ldf, (ftnlen)1, 
				(ftnlen)1);
#line 576 "dtgsyl.f"
		    }
#line 577 "dtgsyl.f"
		    if (j < q) {
#line 578 "dtgsyl.f"
			i__3 = *n - je;
#line 578 "dtgsyl.f"
			dgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js *
				 f_dim1], ldf, &b[js + (je + 1) * b_dim1], 
				ldb, &c_b52, &c__[is + (je + 1) * c_dim1], 
				ldc, (ftnlen)1, (ftnlen)1);
#line 581 "dtgsyl.f"
			i__3 = *n - je;
#line 581 "dtgsyl.f"
			dgemm_("N", "N", &mb, &i__3, &nb, &c_b52, &f[is + js *
				 f_dim1], ldf, &e[js + (je + 1) * e_dim1], 
				lde, &c_b52, &f[is + (je + 1) * f_dim1], ldf, 
				(ftnlen)1, (ftnlen)1);
#line 584 "dtgsyl.f"
		    }
#line 585 "dtgsyl.f"
/* L120: */
#line 585 "dtgsyl.f"
		}
#line 586 "dtgsyl.f"
/* L130: */
#line 586 "dtgsyl.f"
	    }
#line 587 "dtgsyl.f"
	    if (dscale != 0.) {
#line 588 "dtgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 589 "dtgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 590 "dtgsyl.f"
		} else {
#line 591 "dtgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 592 "dtgsyl.f"
		}
#line 593 "dtgsyl.f"
	    }
#line 594 "dtgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 595 "dtgsyl.f"
		if (notran) {
#line 596 "dtgsyl.f"
		    ifunc = *ijob;
#line 597 "dtgsyl.f"
		}
#line 598 "dtgsyl.f"
		scale2 = *scale;
#line 599 "dtgsyl.f"
		dlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 600 "dtgsyl.f"
		dlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 601 "dtgsyl.f"
		dlaset_("F", m, n, &c_b14, &c_b14, &c__[c_offset], ldc, (
			ftnlen)1);
#line 602 "dtgsyl.f"
		dlaset_("F", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)
			1);
#line 603 "dtgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 604 "dtgsyl.f"
		dlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 605 "dtgsyl.f"
		dlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 606 "dtgsyl.f"
		*scale = scale2;
#line 607 "dtgsyl.f"
	    }
#line 608 "dtgsyl.f"
/* L150: */
#line 608 "dtgsyl.f"
	}

#line 610 "dtgsyl.f"
    } else {

/*        Solve transposed (I, J)-subsystem */
/*             A(I, I)**T * R(I, J)  + D(I, I)**T * L(I, J)  =  C(I, J) */
/*             R(I, J)  * B(J, J)**T + L(I, J)  * E(J, J)**T = -F(I, J) */
/*        for I = 1,2,..., P; J = Q, Q-1,..., 1 */

#line 617 "dtgsyl.f"
	*scale = 1.;
#line 618 "dtgsyl.f"
	i__1 = p;
#line 618 "dtgsyl.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 619 "dtgsyl.f"
	    is = iwork[i__];
#line 620 "dtgsyl.f"
	    ie = iwork[i__ + 1] - 1;
#line 621 "dtgsyl.f"
	    mb = ie - is + 1;
#line 622 "dtgsyl.f"
	    i__2 = p + 2;
#line 622 "dtgsyl.f"
	    for (j = q; j >= i__2; --j) {
#line 623 "dtgsyl.f"
		js = iwork[j];
#line 624 "dtgsyl.f"
		je = iwork[j + 1] - 1;
#line 625 "dtgsyl.f"
		nb = je - js + 1;
#line 626 "dtgsyl.f"
		dtgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &
			b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc,
			 &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], 
			lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum, &
			dscale, &iwork[q + 2], &ppqq, &linfo, (ftnlen)1);
#line 631 "dtgsyl.f"
		if (linfo > 0) {
#line 631 "dtgsyl.f"
		    *info = linfo;
#line 631 "dtgsyl.f"
		}
#line 633 "dtgsyl.f"
		if (scaloc != 1.) {
#line 634 "dtgsyl.f"
		    i__3 = js - 1;
#line 634 "dtgsyl.f"
		    for (k = 1; k <= i__3; ++k) {
#line 635 "dtgsyl.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 636 "dtgsyl.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 637 "dtgsyl.f"
/* L160: */
#line 637 "dtgsyl.f"
		    }
#line 638 "dtgsyl.f"
		    i__3 = je;
#line 638 "dtgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 639 "dtgsyl.f"
			i__4 = is - 1;
#line 639 "dtgsyl.f"
			dscal_(&i__4, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 640 "dtgsyl.f"
			i__4 = is - 1;
#line 640 "dtgsyl.f"
			dscal_(&i__4, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 641 "dtgsyl.f"
/* L170: */
#line 641 "dtgsyl.f"
		    }
#line 642 "dtgsyl.f"
		    i__3 = je;
#line 642 "dtgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 643 "dtgsyl.f"
			i__4 = *m - ie;
#line 643 "dtgsyl.f"
			dscal_(&i__4, &scaloc, &c__[ie + 1 + k * c_dim1], &
				c__1);
#line 644 "dtgsyl.f"
			i__4 = *m - ie;
#line 644 "dtgsyl.f"
			dscal_(&i__4, &scaloc, &f[ie + 1 + k * f_dim1], &c__1)
				;
#line 645 "dtgsyl.f"
/* L180: */
#line 645 "dtgsyl.f"
		    }
#line 646 "dtgsyl.f"
		    i__3 = *n;
#line 646 "dtgsyl.f"
		    for (k = je + 1; k <= i__3; ++k) {
#line 647 "dtgsyl.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 648 "dtgsyl.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 649 "dtgsyl.f"
/* L190: */
#line 649 "dtgsyl.f"
		    }
#line 650 "dtgsyl.f"
		    *scale *= scaloc;
#line 651 "dtgsyl.f"
		}

/*              Substitute R(I, J) and L(I, J) into remaining equation. */

#line 655 "dtgsyl.f"
		if (j > p + 2) {
#line 656 "dtgsyl.f"
		    i__3 = js - 1;
#line 656 "dtgsyl.f"
		    dgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &c__[is + js * 
			    c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &c_b52, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 659 "dtgsyl.f"
		    i__3 = js - 1;
#line 659 "dtgsyl.f"
		    dgemm_("N", "T", &mb, &i__3, &nb, &c_b52, &f[is + js * 
			    f_dim1], ldf, &e[js * e_dim1 + 1], lde, &c_b52, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 662 "dtgsyl.f"
		}
#line 663 "dtgsyl.f"
		if (i__ < p) {
#line 664 "dtgsyl.f"
		    i__3 = *m - ie;
#line 664 "dtgsyl.f"
		    dgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &a[is + (ie + 1)
			     * a_dim1], lda, &c__[is + js * c_dim1], ldc, &
			    c_b52, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 667 "dtgsyl.f"
		    i__3 = *m - ie;
#line 667 "dtgsyl.f"
		    dgemm_("T", "N", &i__3, &nb, &mb, &c_b51, &d__[is + (ie + 
			    1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, &
			    c_b52, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 670 "dtgsyl.f"
		}
#line 671 "dtgsyl.f"
/* L200: */
#line 671 "dtgsyl.f"
	    }
#line 672 "dtgsyl.f"
/* L210: */
#line 672 "dtgsyl.f"
	}

#line 674 "dtgsyl.f"
    }

#line 676 "dtgsyl.f"
    work[1] = (doublereal) lwmin;

#line 678 "dtgsyl.f"
    return 0;

/*     End of DTGSYL */

} /* dtgsyl_ */

