#line 1 "ztgsyl.f"
/* ztgsyl.f -- translated by f2c (version 20100827).
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

#line 1 "ztgsyl.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static integer c__1 = 1;
static doublecomplex c_b44 = {-1.,0.};
static doublecomplex c_b45 = {1.,0.};

/* > \brief \b ZTGSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTGSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
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
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSYL solves the generalized Sylvester equation: */
/* > */
/* >             A * R - L * B = scale * C            (1) */
/* >             D * R - L * E = scale * F */
/* > */
/* > where R and L are unknown m-by-n matrices, (A, D), (B, E) and */
/* > (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n, */
/* > respectively, with complex entries. A, B, D and E are upper */
/* > triangular (i.e., (A,D) and (B,E) in generalized Schur form). */
/* > */
/* > The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 */
/* > is an output scaling factor chosen to avoid overflow. */
/* > */
/* > In matrix notation (1) is equivalent to solve Zx = scale*b, where Z */
/* > is defined as */
/* > */
/* >        Z = [ kron(In, A)  -kron(B**H, Im) ]        (2) */
/* >            [ kron(In, D)  -kron(E**H, Im) ], */
/* > */
/* > Here Ix is the identity matrix of size x and X**H is the conjugate */
/* > transpose of X. Kron(X, Y) is the Kronecker product between the */
/* > matrices X and Y. */
/* > */
/* > If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b */
/* > is solved for, which is equivalent to solve for R and L in */
/* > */
/* >             A**H * R + D**H * L = scale * C           (3) */
/* >             R * B**H + L * E**H = scale * -F */
/* > */
/* > This case (TRANS = 'C') is used to compute an one-norm-based estimate */
/* > of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D) */
/* > and (B,E), using ZLACON. */
/* > */
/* > If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of */
/* > Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the */
/* > reciprocal of the smallest singular value of Z. */
/* > */
/* > This is a level-3 BLAS algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': solve the generalized sylvester equation (1). */
/* >          = 'C': solve the "conjugate transposed" system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          Specifies what kind of functionality to be performed. */
/* >          =0: solve (1) only. */
/* >          =1: The functionality of 0 and 3. */
/* >          =2: The functionality of 0 and 4. */
/* >          =3: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* >              (look ahead strategy is used). */
/* >          =4: Only an estimate of Dif[(A,D), (B,E)] is computed. */
/* >              (ZGECON on sub-systems is used). */
/* >          Not referenced if TRANS = 'C'. */
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
/* >          A is COMPLEX*16 array, dimension (LDA, M) */
/* >          The upper triangular matrix A. */
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
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
/* >          The upper triangular matrix B. */
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
/* >          C is COMPLEX*16 array, dimension (LDC, N) */
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
/* >          D is COMPLEX*16 array, dimension (LDD, M) */
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
/* >          E is COMPLEX*16 array, dimension (LDE, N) */
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
/* >          F is COMPLEX*16 array, dimension (LDF, N) */
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
/* >          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2). */
/* >          IF IJOB = 0 or TRANS = 'C', DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          On exit SCALE is the scaling factor in (1) or (3). */
/* >          If 0 < SCALE < 1, C and F hold the solutions R and L, resp., */
/* >          to a slightly perturbed system but the input matrices A, B, */
/* >          D and E have not been changed. If SCALE = 0, R and L will */
/* >          hold the solutions to the homogenious system with C = F = 0. */
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
/* >          IWORK is INTEGER array, dimension (M+N+2) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >            =0: successful exit */
/* >            <0: If INFO = -i, the i-th argument had an illegal value. */
/* >            >0: (A, D) and (B, E) have common or very close */
/* >                eigenvalues. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16SYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* >      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22, */
/* >      No 1, 1996. */
/* > \n */
/* >  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester */
/* >      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal. */
/* >      Appl., 15(4):1045-1060, 1994. */
/* > \n */
/* >  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with */
/* >      Condition Estimators for Solving the Generalized Sylvester */
/* >      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7, */
/* >      July 1989, pp 745-751. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *dif, doublecomplex *work, integer *
	lwork, integer *iwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3, 
	    i__4;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p, q, ie, je, mb, nb, is, js, pq;
    static doublereal dsum;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ifunc, linfo;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer lwmin;
    static doublereal scale2, dscale;
    extern /* Subroutine */ int ztgsy2_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iround;
    static logical notran;
    static integer isolve;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
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
/*  Replaced various illegal calls to CCOPY by calls to CLASET. */
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

#line 347 "ztgsyl.f"
    /* Parameter adjustments */
#line 347 "ztgsyl.f"
    a_dim1 = *lda;
#line 347 "ztgsyl.f"
    a_offset = 1 + a_dim1;
#line 347 "ztgsyl.f"
    a -= a_offset;
#line 347 "ztgsyl.f"
    b_dim1 = *ldb;
#line 347 "ztgsyl.f"
    b_offset = 1 + b_dim1;
#line 347 "ztgsyl.f"
    b -= b_offset;
#line 347 "ztgsyl.f"
    c_dim1 = *ldc;
#line 347 "ztgsyl.f"
    c_offset = 1 + c_dim1;
#line 347 "ztgsyl.f"
    c__ -= c_offset;
#line 347 "ztgsyl.f"
    d_dim1 = *ldd;
#line 347 "ztgsyl.f"
    d_offset = 1 + d_dim1;
#line 347 "ztgsyl.f"
    d__ -= d_offset;
#line 347 "ztgsyl.f"
    e_dim1 = *lde;
#line 347 "ztgsyl.f"
    e_offset = 1 + e_dim1;
#line 347 "ztgsyl.f"
    e -= e_offset;
#line 347 "ztgsyl.f"
    f_dim1 = *ldf;
#line 347 "ztgsyl.f"
    f_offset = 1 + f_dim1;
#line 347 "ztgsyl.f"
    f -= f_offset;
#line 347 "ztgsyl.f"
    --work;
#line 347 "ztgsyl.f"
    --iwork;
#line 347 "ztgsyl.f"

#line 347 "ztgsyl.f"
    /* Function Body */
#line 347 "ztgsyl.f"
    *info = 0;
#line 348 "ztgsyl.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 349 "ztgsyl.f"
    lquery = *lwork == -1;

#line 351 "ztgsyl.f"
    if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 352 "ztgsyl.f"
	*info = -1;
#line 353 "ztgsyl.f"
    } else if (notran) {
#line 354 "ztgsyl.f"
	if (*ijob < 0 || *ijob > 4) {
#line 355 "ztgsyl.f"
	    *info = -2;
#line 356 "ztgsyl.f"
	}
#line 357 "ztgsyl.f"
    }
#line 358 "ztgsyl.f"
    if (*info == 0) {
#line 359 "ztgsyl.f"
	if (*m <= 0) {
#line 360 "ztgsyl.f"
	    *info = -3;
#line 361 "ztgsyl.f"
	} else if (*n <= 0) {
#line 362 "ztgsyl.f"
	    *info = -4;
#line 363 "ztgsyl.f"
	} else if (*lda < max(1,*m)) {
#line 364 "ztgsyl.f"
	    *info = -6;
#line 365 "ztgsyl.f"
	} else if (*ldb < max(1,*n)) {
#line 366 "ztgsyl.f"
	    *info = -8;
#line 367 "ztgsyl.f"
	} else if (*ldc < max(1,*m)) {
#line 368 "ztgsyl.f"
	    *info = -10;
#line 369 "ztgsyl.f"
	} else if (*ldd < max(1,*m)) {
#line 370 "ztgsyl.f"
	    *info = -12;
#line 371 "ztgsyl.f"
	} else if (*lde < max(1,*n)) {
#line 372 "ztgsyl.f"
	    *info = -14;
#line 373 "ztgsyl.f"
	} else if (*ldf < max(1,*m)) {
#line 374 "ztgsyl.f"
	    *info = -16;
#line 375 "ztgsyl.f"
	}
#line 376 "ztgsyl.f"
    }

#line 378 "ztgsyl.f"
    if (*info == 0) {
#line 379 "ztgsyl.f"
	if (notran) {
#line 380 "ztgsyl.f"
	    if (*ijob == 1 || *ijob == 2) {
/* Computing MAX */
#line 381 "ztgsyl.f"
		i__1 = 1, i__2 = (*m << 1) * *n;
#line 381 "ztgsyl.f"
		lwmin = max(i__1,i__2);
#line 382 "ztgsyl.f"
	    } else {
#line 383 "ztgsyl.f"
		lwmin = 1;
#line 384 "ztgsyl.f"
	    }
#line 385 "ztgsyl.f"
	} else {
#line 386 "ztgsyl.f"
	    lwmin = 1;
#line 387 "ztgsyl.f"
	}
#line 388 "ztgsyl.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 390 "ztgsyl.f"
	if (*lwork < lwmin && ! lquery) {
#line 391 "ztgsyl.f"
	    *info = -20;
#line 392 "ztgsyl.f"
	}
#line 393 "ztgsyl.f"
    }

#line 395 "ztgsyl.f"
    if (*info != 0) {
#line 396 "ztgsyl.f"
	i__1 = -(*info);
#line 396 "ztgsyl.f"
	xerbla_("ZTGSYL", &i__1, (ftnlen)6);
#line 397 "ztgsyl.f"
	return 0;
#line 398 "ztgsyl.f"
    } else if (lquery) {
#line 399 "ztgsyl.f"
	return 0;
#line 400 "ztgsyl.f"
    }

/*     Quick return if possible */

#line 404 "ztgsyl.f"
    if (*m == 0 || *n == 0) {
#line 405 "ztgsyl.f"
	*scale = 1.;
#line 406 "ztgsyl.f"
	if (notran) {
#line 407 "ztgsyl.f"
	    if (*ijob != 0) {
#line 408 "ztgsyl.f"
		*dif = 0.;
#line 409 "ztgsyl.f"
	    }
#line 410 "ztgsyl.f"
	}
#line 411 "ztgsyl.f"
	return 0;
#line 412 "ztgsyl.f"
    }

/*     Determine  optimal block sizes MB and NB */

#line 416 "ztgsyl.f"
    mb = ilaenv_(&c__2, "ZTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 417 "ztgsyl.f"
    nb = ilaenv_(&c__5, "ZTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 419 "ztgsyl.f"
    isolve = 1;
#line 420 "ztgsyl.f"
    ifunc = 0;
#line 421 "ztgsyl.f"
    if (notran) {
#line 422 "ztgsyl.f"
	if (*ijob >= 3) {
#line 423 "ztgsyl.f"
	    ifunc = *ijob - 2;
#line 424 "ztgsyl.f"
	    zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)1);
#line 425 "ztgsyl.f"
	    zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1);
#line 426 "ztgsyl.f"
	} else if (*ijob >= 1 && notran) {
#line 427 "ztgsyl.f"
	    isolve = 2;
#line 428 "ztgsyl.f"
	}
#line 429 "ztgsyl.f"
    }

#line 431 "ztgsyl.f"
    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n) {

/*        Use unblocked Level 2 solver */

#line 436 "ztgsyl.f"
	i__1 = isolve;
#line 436 "ztgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

#line 438 "ztgsyl.f"
	    *scale = 1.;
#line 439 "ztgsyl.f"
	    dscale = 0.;
#line 440 "ztgsyl.f"
	    dsum = 1.;
#line 441 "ztgsyl.f"
	    pq = *m * *n;
#line 442 "ztgsyl.f"
	    ztgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb,
		     &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], 
		    lde, &f[f_offset], ldf, scale, &dsum, &dscale, info, (
		    ftnlen)1);
#line 445 "ztgsyl.f"
	    if (dscale != 0.) {
#line 446 "ztgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 447 "ztgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 448 "ztgsyl.f"
		} else {
#line 449 "ztgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 450 "ztgsyl.f"
		}
#line 451 "ztgsyl.f"
	    }
#line 452 "ztgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 453 "ztgsyl.f"
		if (notran) {
#line 454 "ztgsyl.f"
		    ifunc = *ijob;
#line 455 "ztgsyl.f"
		}
#line 456 "ztgsyl.f"
		scale2 = *scale;
#line 457 "ztgsyl.f"
		zlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 458 "ztgsyl.f"
		zlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 459 "ztgsyl.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)
			1);
#line 460 "ztgsyl.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1)
			;
#line 461 "ztgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 462 "ztgsyl.f"
		zlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 463 "ztgsyl.f"
		zlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 464 "ztgsyl.f"
		*scale = scale2;
#line 465 "ztgsyl.f"
	    }
#line 466 "ztgsyl.f"
/* L30: */
#line 466 "ztgsyl.f"
	}

#line 468 "ztgsyl.f"
	return 0;

#line 470 "ztgsyl.f"
    }

/*     Determine block structure of A */

#line 474 "ztgsyl.f"
    p = 0;
#line 475 "ztgsyl.f"
    i__ = 1;
#line 476 "ztgsyl.f"
L40:
#line 477 "ztgsyl.f"
    if (i__ > *m) {
#line 477 "ztgsyl.f"
	goto L50;
#line 477 "ztgsyl.f"
    }
#line 479 "ztgsyl.f"
    ++p;
#line 480 "ztgsyl.f"
    iwork[p] = i__;
#line 481 "ztgsyl.f"
    i__ += mb;
#line 482 "ztgsyl.f"
    if (i__ >= *m) {
#line 482 "ztgsyl.f"
	goto L50;
#line 482 "ztgsyl.f"
    }
#line 484 "ztgsyl.f"
    goto L40;
#line 485 "ztgsyl.f"
L50:
#line 486 "ztgsyl.f"
    iwork[p + 1] = *m + 1;
#line 487 "ztgsyl.f"
    if (iwork[p] == iwork[p + 1]) {
#line 487 "ztgsyl.f"
	--p;
#line 487 "ztgsyl.f"
    }

/*     Determine block structure of B */

#line 492 "ztgsyl.f"
    q = p + 1;
#line 493 "ztgsyl.f"
    j = 1;
#line 494 "ztgsyl.f"
L60:
#line 495 "ztgsyl.f"
    if (j > *n) {
#line 495 "ztgsyl.f"
	goto L70;
#line 495 "ztgsyl.f"
    }

#line 498 "ztgsyl.f"
    ++q;
#line 499 "ztgsyl.f"
    iwork[q] = j;
#line 500 "ztgsyl.f"
    j += nb;
#line 501 "ztgsyl.f"
    if (j >= *n) {
#line 501 "ztgsyl.f"
	goto L70;
#line 501 "ztgsyl.f"
    }
#line 503 "ztgsyl.f"
    goto L60;

#line 505 "ztgsyl.f"
L70:
#line 506 "ztgsyl.f"
    iwork[q + 1] = *n + 1;
#line 507 "ztgsyl.f"
    if (iwork[q] == iwork[q + 1]) {
#line 507 "ztgsyl.f"
	--q;
#line 507 "ztgsyl.f"
    }

#line 510 "ztgsyl.f"
    if (notran) {
#line 511 "ztgsyl.f"
	i__1 = isolve;
#line 511 "ztgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Solve (I, J) - subsystem */
/*               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*           for I = P, P - 1, ..., 1; J = 1, 2, ..., Q */

#line 518 "ztgsyl.f"
	    pq = 0;
#line 519 "ztgsyl.f"
	    *scale = 1.;
#line 520 "ztgsyl.f"
	    dscale = 0.;
#line 521 "ztgsyl.f"
	    dsum = 1.;
#line 522 "ztgsyl.f"
	    i__2 = q;
#line 522 "ztgsyl.f"
	    for (j = p + 2; j <= i__2; ++j) {
#line 523 "ztgsyl.f"
		js = iwork[j];
#line 524 "ztgsyl.f"
		je = iwork[j + 1] - 1;
#line 525 "ztgsyl.f"
		nb = je - js + 1;
#line 526 "ztgsyl.f"
		for (i__ = p; i__ >= 1; --i__) {
#line 527 "ztgsyl.f"
		    is = iwork[i__];
#line 528 "ztgsyl.f"
		    ie = iwork[i__ + 1] - 1;
#line 529 "ztgsyl.f"
		    mb = ie - is + 1;
#line 530 "ztgsyl.f"
		    ztgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], 
			    lda, &b[js + js * b_dim1], ldb, &c__[is + js * 
			    c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js 
			    + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &
			    scaloc, &dsum, &dscale, &linfo, (ftnlen)1);
#line 535 "ztgsyl.f"
		    if (linfo > 0) {
#line 535 "ztgsyl.f"
			*info = linfo;
#line 535 "ztgsyl.f"
		    }
#line 537 "ztgsyl.f"
		    pq += mb * nb;
#line 538 "ztgsyl.f"
		    if (scaloc != 1.) {
#line 539 "ztgsyl.f"
			i__3 = js - 1;
#line 539 "ztgsyl.f"
			for (k = 1; k <= i__3; ++k) {
#line 540 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 540 "ztgsyl.f"
			    zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 542 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 542 "ztgsyl.f"
			    zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 544 "ztgsyl.f"
/* L80: */
#line 544 "ztgsyl.f"
			}
#line 545 "ztgsyl.f"
			i__3 = je;
#line 545 "ztgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 546 "ztgsyl.f"
			    i__4 = is - 1;
#line 546 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 546 "ztgsyl.f"
			    zscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 548 "ztgsyl.f"
			    i__4 = is - 1;
#line 548 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 548 "ztgsyl.f"
			    zscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 550 "ztgsyl.f"
/* L90: */
#line 550 "ztgsyl.f"
			}
#line 551 "ztgsyl.f"
			i__3 = je;
#line 551 "ztgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 552 "ztgsyl.f"
			    i__4 = *m - ie;
#line 552 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 552 "ztgsyl.f"
			    zscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &
				    c__1);
#line 554 "ztgsyl.f"
			    i__4 = *m - ie;
#line 554 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 554 "ztgsyl.f"
			    zscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &
				    c__1);
#line 556 "ztgsyl.f"
/* L100: */
#line 556 "ztgsyl.f"
			}
#line 557 "ztgsyl.f"
			i__3 = *n;
#line 557 "ztgsyl.f"
			for (k = je + 1; k <= i__3; ++k) {
#line 558 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 558 "ztgsyl.f"
			    zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 560 "ztgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 560 "ztgsyl.f"
			    zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 562 "ztgsyl.f"
/* L110: */
#line 562 "ztgsyl.f"
			}
#line 563 "ztgsyl.f"
			*scale *= scaloc;
#line 564 "ztgsyl.f"
		    }

/*                 Substitute R(I,J) and L(I,J) into remaining equation. */

#line 568 "ztgsyl.f"
		    if (i__ > 1) {
#line 569 "ztgsyl.f"
			i__3 = is - 1;
#line 569 "ztgsyl.f"
			zgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &a[is * 
				a_dim1 + 1], lda, &c__[is + js * c_dim1], ldc,
				 &c_b45, &c__[js * c_dim1 + 1], ldc, (ftnlen)
				1, (ftnlen)1);
#line 573 "ztgsyl.f"
			i__3 = is - 1;
#line 573 "ztgsyl.f"
			zgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &d__[is * 
				d_dim1 + 1], ldd, &c__[is + js * c_dim1], ldc,
				 &c_b45, &f[js * f_dim1 + 1], ldf, (ftnlen)1, 
				(ftnlen)1);
#line 577 "ztgsyl.f"
		    }
#line 578 "ztgsyl.f"
		    if (j < q) {
#line 579 "ztgsyl.f"
			i__3 = *n - je;
#line 579 "ztgsyl.f"
			zgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js *
				 f_dim1], ldf, &b[js + (je + 1) * b_dim1], 
				ldb, &c_b45, &c__[is + (je + 1) * c_dim1], 
				ldc, (ftnlen)1, (ftnlen)1);
#line 584 "ztgsyl.f"
			i__3 = *n - je;
#line 584 "ztgsyl.f"
			zgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js *
				 f_dim1], ldf, &e[js + (je + 1) * e_dim1], 
				lde, &c_b45, &f[is + (je + 1) * f_dim1], ldf, 
				(ftnlen)1, (ftnlen)1);
#line 589 "ztgsyl.f"
		    }
#line 590 "ztgsyl.f"
/* L120: */
#line 590 "ztgsyl.f"
		}
#line 591 "ztgsyl.f"
/* L130: */
#line 591 "ztgsyl.f"
	    }
#line 592 "ztgsyl.f"
	    if (dscale != 0.) {
#line 593 "ztgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 594 "ztgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 595 "ztgsyl.f"
		} else {
#line 596 "ztgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 597 "ztgsyl.f"
		}
#line 598 "ztgsyl.f"
	    }
#line 599 "ztgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 600 "ztgsyl.f"
		if (notran) {
#line 601 "ztgsyl.f"
		    ifunc = *ijob;
#line 602 "ztgsyl.f"
		}
#line 603 "ztgsyl.f"
		scale2 = *scale;
#line 604 "ztgsyl.f"
		zlacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 605 "ztgsyl.f"
		zlacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 606 "ztgsyl.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)
			1);
#line 607 "ztgsyl.f"
		zlaset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1)
			;
#line 608 "ztgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 609 "ztgsyl.f"
		zlacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 610 "ztgsyl.f"
		zlacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 611 "ztgsyl.f"
		*scale = scale2;
#line 612 "ztgsyl.f"
	    }
#line 613 "ztgsyl.f"
/* L150: */
#line 613 "ztgsyl.f"
	}
#line 614 "ztgsyl.f"
    } else {

/*        Solve transposed (I, J)-subsystem */
/*            A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J) */
/*            R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J) */
/*        for I = 1,2,..., P; J = Q, Q-1,..., 1 */

#line 621 "ztgsyl.f"
	*scale = 1.;
#line 622 "ztgsyl.f"
	i__1 = p;
#line 622 "ztgsyl.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 623 "ztgsyl.f"
	    is = iwork[i__];
#line 624 "ztgsyl.f"
	    ie = iwork[i__ + 1] - 1;
#line 625 "ztgsyl.f"
	    mb = ie - is + 1;
#line 626 "ztgsyl.f"
	    i__2 = p + 2;
#line 626 "ztgsyl.f"
	    for (j = q; j >= i__2; --j) {
#line 627 "ztgsyl.f"
		js = iwork[j];
#line 628 "ztgsyl.f"
		je = iwork[j + 1] - 1;
#line 629 "ztgsyl.f"
		nb = je - js + 1;
#line 630 "ztgsyl.f"
		ztgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &
			b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc,
			 &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], 
			lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum, &
			dscale, &linfo, (ftnlen)1);
#line 635 "ztgsyl.f"
		if (linfo > 0) {
#line 635 "ztgsyl.f"
		    *info = linfo;
#line 635 "ztgsyl.f"
		}
#line 637 "ztgsyl.f"
		if (scaloc != 1.) {
#line 638 "ztgsyl.f"
		    i__3 = js - 1;
#line 638 "ztgsyl.f"
		    for (k = 1; k <= i__3; ++k) {
#line 639 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 639 "ztgsyl.f"
			zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 641 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 641 "ztgsyl.f"
			zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 643 "ztgsyl.f"
/* L160: */
#line 643 "ztgsyl.f"
		    }
#line 644 "ztgsyl.f"
		    i__3 = je;
#line 644 "ztgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 645 "ztgsyl.f"
			i__4 = is - 1;
#line 645 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 645 "ztgsyl.f"
			zscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 647 "ztgsyl.f"
			i__4 = is - 1;
#line 647 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 647 "ztgsyl.f"
			zscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 649 "ztgsyl.f"
/* L170: */
#line 649 "ztgsyl.f"
		    }
#line 650 "ztgsyl.f"
		    i__3 = je;
#line 650 "ztgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 651 "ztgsyl.f"
			i__4 = *m - ie;
#line 651 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 651 "ztgsyl.f"
			zscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &c__1)
				;
#line 653 "ztgsyl.f"
			i__4 = *m - ie;
#line 653 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 653 "ztgsyl.f"
			zscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &c__1);
#line 655 "ztgsyl.f"
/* L180: */
#line 655 "ztgsyl.f"
		    }
#line 656 "ztgsyl.f"
		    i__3 = *n;
#line 656 "ztgsyl.f"
		    for (k = je + 1; k <= i__3; ++k) {
#line 657 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 657 "ztgsyl.f"
			zscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 659 "ztgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 659 "ztgsyl.f"
			zscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 661 "ztgsyl.f"
/* L190: */
#line 661 "ztgsyl.f"
		    }
#line 662 "ztgsyl.f"
		    *scale *= scaloc;
#line 663 "ztgsyl.f"
		}

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 667 "ztgsyl.f"
		if (j > p + 2) {
#line 668 "ztgsyl.f"
		    i__3 = js - 1;
#line 668 "ztgsyl.f"
		    zgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &c__[is + js * 
			    c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &c_b45, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 672 "ztgsyl.f"
		    i__3 = js - 1;
#line 672 "ztgsyl.f"
		    zgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &f[is + js * 
			    f_dim1], ldf, &e[js * e_dim1 + 1], lde, &c_b45, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 676 "ztgsyl.f"
		}
#line 677 "ztgsyl.f"
		if (i__ < p) {
#line 678 "ztgsyl.f"
		    i__3 = *m - ie;
#line 678 "ztgsyl.f"
		    zgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &a[is + (ie + 1)
			     * a_dim1], lda, &c__[is + js * c_dim1], ldc, &
			    c_b45, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 682 "ztgsyl.f"
		    i__3 = *m - ie;
#line 682 "ztgsyl.f"
		    zgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &d__[is + (ie + 
			    1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, &
			    c_b45, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 686 "ztgsyl.f"
		}
#line 687 "ztgsyl.f"
/* L200: */
#line 687 "ztgsyl.f"
	    }
#line 688 "ztgsyl.f"
/* L210: */
#line 688 "ztgsyl.f"
	}
#line 689 "ztgsyl.f"
    }

#line 691 "ztgsyl.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 693 "ztgsyl.f"
    return 0;

/*     End of ZTGSYL */

} /* ztgsyl_ */

