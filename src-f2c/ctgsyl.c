#line 1 "ctgsyl.f"
/* ctgsyl.f -- translated by f2c (version 20100827).
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

#line 1 "ctgsyl.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__5 = 5;
static integer c__1 = 1;
static doublecomplex c_b44 = {-1.,0.};
static doublecomplex c_b45 = {1.,0.};

/* > \brief \b CTGSYL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGSYL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgsyl.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgsyl.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgsyl.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
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
/*       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGSYL solves the generalized Sylvester equation: */
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
/* > and (B,E), using CLACON. */
/* > */
/* > If IJOB >= 1, CTGSYL computes a Frobenius norm-based estimate of */
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
/* >              (CGECON on sub-systems is used). */
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
/* >          A is COMPLEX array, dimension (LDA, M) */
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
/* >          B is COMPLEX array, dimension (LDB, N) */
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
/* >          C is COMPLEX array, dimension (LDC, N) */
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
/* >          D is COMPLEX array, dimension (LDD, M) */
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
/* >          E is COMPLEX array, dimension (LDE, N) */
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
/* >          F is COMPLEX array, dimension (LDF, N) */
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
/* >          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2). */
/* >          IF IJOB = 0 or TRANS = 'C', DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is REAL */
/* >          On exit SCALE is the scaling factor in (1) or (3). */
/* >          If 0 < SCALE < 1, C and F hold the solutions R and L, resp., */
/* >          to a slightly perturbed system but the input matrices A, B, */
/* >          D and E have not been changed. If SCALE = 0, R and L will */
/* >          hold the solutions to the homogenious system with C = F = 0. */
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

/* > \date December 2016 */

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int ctgsyl_(char *trans, integer *ijob, integer *m, integer *
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cgemm_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ifunc, linfo, lwmin;
    static doublereal scale2;
    extern /* Subroutine */ int ctgsy2_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal dscale, scaloc;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer iround;
    static logical notran;
    static integer isolve;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 347 "ctgsyl.f"
    /* Parameter adjustments */
#line 347 "ctgsyl.f"
    a_dim1 = *lda;
#line 347 "ctgsyl.f"
    a_offset = 1 + a_dim1;
#line 347 "ctgsyl.f"
    a -= a_offset;
#line 347 "ctgsyl.f"
    b_dim1 = *ldb;
#line 347 "ctgsyl.f"
    b_offset = 1 + b_dim1;
#line 347 "ctgsyl.f"
    b -= b_offset;
#line 347 "ctgsyl.f"
    c_dim1 = *ldc;
#line 347 "ctgsyl.f"
    c_offset = 1 + c_dim1;
#line 347 "ctgsyl.f"
    c__ -= c_offset;
#line 347 "ctgsyl.f"
    d_dim1 = *ldd;
#line 347 "ctgsyl.f"
    d_offset = 1 + d_dim1;
#line 347 "ctgsyl.f"
    d__ -= d_offset;
#line 347 "ctgsyl.f"
    e_dim1 = *lde;
#line 347 "ctgsyl.f"
    e_offset = 1 + e_dim1;
#line 347 "ctgsyl.f"
    e -= e_offset;
#line 347 "ctgsyl.f"
    f_dim1 = *ldf;
#line 347 "ctgsyl.f"
    f_offset = 1 + f_dim1;
#line 347 "ctgsyl.f"
    f -= f_offset;
#line 347 "ctgsyl.f"
    --work;
#line 347 "ctgsyl.f"
    --iwork;
#line 347 "ctgsyl.f"

#line 347 "ctgsyl.f"
    /* Function Body */
#line 347 "ctgsyl.f"
    *info = 0;
#line 348 "ctgsyl.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 349 "ctgsyl.f"
    lquery = *lwork == -1;

#line 351 "ctgsyl.f"
    if (! notran && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 352 "ctgsyl.f"
	*info = -1;
#line 353 "ctgsyl.f"
    } else if (notran) {
#line 354 "ctgsyl.f"
	if (*ijob < 0 || *ijob > 4) {
#line 355 "ctgsyl.f"
	    *info = -2;
#line 356 "ctgsyl.f"
	}
#line 357 "ctgsyl.f"
    }
#line 358 "ctgsyl.f"
    if (*info == 0) {
#line 359 "ctgsyl.f"
	if (*m <= 0) {
#line 360 "ctgsyl.f"
	    *info = -3;
#line 361 "ctgsyl.f"
	} else if (*n <= 0) {
#line 362 "ctgsyl.f"
	    *info = -4;
#line 363 "ctgsyl.f"
	} else if (*lda < max(1,*m)) {
#line 364 "ctgsyl.f"
	    *info = -6;
#line 365 "ctgsyl.f"
	} else if (*ldb < max(1,*n)) {
#line 366 "ctgsyl.f"
	    *info = -8;
#line 367 "ctgsyl.f"
	} else if (*ldc < max(1,*m)) {
#line 368 "ctgsyl.f"
	    *info = -10;
#line 369 "ctgsyl.f"
	} else if (*ldd < max(1,*m)) {
#line 370 "ctgsyl.f"
	    *info = -12;
#line 371 "ctgsyl.f"
	} else if (*lde < max(1,*n)) {
#line 372 "ctgsyl.f"
	    *info = -14;
#line 373 "ctgsyl.f"
	} else if (*ldf < max(1,*m)) {
#line 374 "ctgsyl.f"
	    *info = -16;
#line 375 "ctgsyl.f"
	}
#line 376 "ctgsyl.f"
    }

#line 378 "ctgsyl.f"
    if (*info == 0) {
#line 379 "ctgsyl.f"
	if (notran) {
#line 380 "ctgsyl.f"
	    if (*ijob == 1 || *ijob == 2) {
/* Computing MAX */
#line 381 "ctgsyl.f"
		i__1 = 1, i__2 = (*m << 1) * *n;
#line 381 "ctgsyl.f"
		lwmin = max(i__1,i__2);
#line 382 "ctgsyl.f"
	    } else {
#line 383 "ctgsyl.f"
		lwmin = 1;
#line 384 "ctgsyl.f"
	    }
#line 385 "ctgsyl.f"
	} else {
#line 386 "ctgsyl.f"
	    lwmin = 1;
#line 387 "ctgsyl.f"
	}
#line 388 "ctgsyl.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 390 "ctgsyl.f"
	if (*lwork < lwmin && ! lquery) {
#line 391 "ctgsyl.f"
	    *info = -20;
#line 392 "ctgsyl.f"
	}
#line 393 "ctgsyl.f"
    }

#line 395 "ctgsyl.f"
    if (*info != 0) {
#line 396 "ctgsyl.f"
	i__1 = -(*info);
#line 396 "ctgsyl.f"
	xerbla_("CTGSYL", &i__1, (ftnlen)6);
#line 397 "ctgsyl.f"
	return 0;
#line 398 "ctgsyl.f"
    } else if (lquery) {
#line 399 "ctgsyl.f"
	return 0;
#line 400 "ctgsyl.f"
    }

/*     Quick return if possible */

#line 404 "ctgsyl.f"
    if (*m == 0 || *n == 0) {
#line 405 "ctgsyl.f"
	*scale = 1.;
#line 406 "ctgsyl.f"
	if (notran) {
#line 407 "ctgsyl.f"
	    if (*ijob != 0) {
#line 408 "ctgsyl.f"
		*dif = 0.;
#line 409 "ctgsyl.f"
	    }
#line 410 "ctgsyl.f"
	}
#line 411 "ctgsyl.f"
	return 0;
#line 412 "ctgsyl.f"
    }

/*     Determine  optimal block sizes MB and NB */

#line 416 "ctgsyl.f"
    mb = ilaenv_(&c__2, "CTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 417 "ctgsyl.f"
    nb = ilaenv_(&c__5, "CTGSYL", trans, m, n, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 419 "ctgsyl.f"
    isolve = 1;
#line 420 "ctgsyl.f"
    ifunc = 0;
#line 421 "ctgsyl.f"
    if (notran) {
#line 422 "ctgsyl.f"
	if (*ijob >= 3) {
#line 423 "ctgsyl.f"
	    ifunc = *ijob - 2;
#line 424 "ctgsyl.f"
	    claset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)1);
#line 425 "ctgsyl.f"
	    claset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1);
#line 426 "ctgsyl.f"
	} else if (*ijob >= 1 && notran) {
#line 427 "ctgsyl.f"
	    isolve = 2;
#line 428 "ctgsyl.f"
	}
#line 429 "ctgsyl.f"
    }

#line 431 "ctgsyl.f"
    if (mb <= 1 && nb <= 1 || mb >= *m && nb >= *n) {

/*        Use unblocked Level 2 solver */

#line 436 "ctgsyl.f"
	i__1 = isolve;
#line 436 "ctgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

#line 438 "ctgsyl.f"
	    *scale = 1.;
#line 439 "ctgsyl.f"
	    dscale = 0.;
#line 440 "ctgsyl.f"
	    dsum = 1.;
#line 441 "ctgsyl.f"
	    pq = *m * *n;
#line 442 "ctgsyl.f"
	    ctgsy2_(trans, &ifunc, m, n, &a[a_offset], lda, &b[b_offset], ldb,
		     &c__[c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], 
		    lde, &f[f_offset], ldf, scale, &dsum, &dscale, info, (
		    ftnlen)1);
#line 445 "ctgsyl.f"
	    if (dscale != 0.) {
#line 446 "ctgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 447 "ctgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 448 "ctgsyl.f"
		} else {
#line 449 "ctgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 450 "ctgsyl.f"
		}
#line 451 "ctgsyl.f"
	    }
#line 452 "ctgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 453 "ctgsyl.f"
		if (notran) {
#line 454 "ctgsyl.f"
		    ifunc = *ijob;
#line 455 "ctgsyl.f"
		}
#line 456 "ctgsyl.f"
		scale2 = *scale;
#line 457 "ctgsyl.f"
		clacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 458 "ctgsyl.f"
		clacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 459 "ctgsyl.f"
		claset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)
			1);
#line 460 "ctgsyl.f"
		claset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1)
			;
#line 461 "ctgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 462 "ctgsyl.f"
		clacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 463 "ctgsyl.f"
		clacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 464 "ctgsyl.f"
		*scale = scale2;
#line 465 "ctgsyl.f"
	    }
#line 466 "ctgsyl.f"
/* L30: */
#line 466 "ctgsyl.f"
	}

#line 468 "ctgsyl.f"
	return 0;

#line 470 "ctgsyl.f"
    }

/*     Determine block structure of A */

#line 474 "ctgsyl.f"
    p = 0;
#line 475 "ctgsyl.f"
    i__ = 1;
#line 476 "ctgsyl.f"
L40:
#line 477 "ctgsyl.f"
    if (i__ > *m) {
#line 477 "ctgsyl.f"
	goto L50;
#line 477 "ctgsyl.f"
    }
#line 479 "ctgsyl.f"
    ++p;
#line 480 "ctgsyl.f"
    iwork[p] = i__;
#line 481 "ctgsyl.f"
    i__ += mb;
#line 482 "ctgsyl.f"
    if (i__ >= *m) {
#line 482 "ctgsyl.f"
	goto L50;
#line 482 "ctgsyl.f"
    }
#line 484 "ctgsyl.f"
    goto L40;
#line 485 "ctgsyl.f"
L50:
#line 486 "ctgsyl.f"
    iwork[p + 1] = *m + 1;
#line 487 "ctgsyl.f"
    if (iwork[p] == iwork[p + 1]) {
#line 487 "ctgsyl.f"
	--p;
#line 487 "ctgsyl.f"
    }

/*     Determine block structure of B */

#line 492 "ctgsyl.f"
    q = p + 1;
#line 493 "ctgsyl.f"
    j = 1;
#line 494 "ctgsyl.f"
L60:
#line 495 "ctgsyl.f"
    if (j > *n) {
#line 495 "ctgsyl.f"
	goto L70;
#line 495 "ctgsyl.f"
    }

#line 498 "ctgsyl.f"
    ++q;
#line 499 "ctgsyl.f"
    iwork[q] = j;
#line 500 "ctgsyl.f"
    j += nb;
#line 501 "ctgsyl.f"
    if (j >= *n) {
#line 501 "ctgsyl.f"
	goto L70;
#line 501 "ctgsyl.f"
    }
#line 503 "ctgsyl.f"
    goto L60;

#line 505 "ctgsyl.f"
L70:
#line 506 "ctgsyl.f"
    iwork[q + 1] = *n + 1;
#line 507 "ctgsyl.f"
    if (iwork[q] == iwork[q + 1]) {
#line 507 "ctgsyl.f"
	--q;
#line 507 "ctgsyl.f"
    }

#line 510 "ctgsyl.f"
    if (notran) {
#line 511 "ctgsyl.f"
	i__1 = isolve;
#line 511 "ctgsyl.f"
	for (iround = 1; iround <= i__1; ++iround) {

/*           Solve (I, J) - subsystem */
/*               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*           for I = P, P - 1, ..., 1; J = 1, 2, ..., Q */

#line 518 "ctgsyl.f"
	    pq = 0;
#line 519 "ctgsyl.f"
	    *scale = 1.;
#line 520 "ctgsyl.f"
	    dscale = 0.;
#line 521 "ctgsyl.f"
	    dsum = 1.;
#line 522 "ctgsyl.f"
	    i__2 = q;
#line 522 "ctgsyl.f"
	    for (j = p + 2; j <= i__2; ++j) {
#line 523 "ctgsyl.f"
		js = iwork[j];
#line 524 "ctgsyl.f"
		je = iwork[j + 1] - 1;
#line 525 "ctgsyl.f"
		nb = je - js + 1;
#line 526 "ctgsyl.f"
		for (i__ = p; i__ >= 1; --i__) {
#line 527 "ctgsyl.f"
		    is = iwork[i__];
#line 528 "ctgsyl.f"
		    ie = iwork[i__ + 1] - 1;
#line 529 "ctgsyl.f"
		    mb = ie - is + 1;
#line 530 "ctgsyl.f"
		    ctgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], 
			    lda, &b[js + js * b_dim1], ldb, &c__[is + js * 
			    c_dim1], ldc, &d__[is + is * d_dim1], ldd, &e[js 
			    + js * e_dim1], lde, &f[is + js * f_dim1], ldf, &
			    scaloc, &dsum, &dscale, &linfo, (ftnlen)1);
#line 535 "ctgsyl.f"
		    if (linfo > 0) {
#line 535 "ctgsyl.f"
			*info = linfo;
#line 535 "ctgsyl.f"
		    }
#line 537 "ctgsyl.f"
		    pq += mb * nb;
#line 538 "ctgsyl.f"
		    if (scaloc != 1.) {
#line 539 "ctgsyl.f"
			i__3 = js - 1;
#line 539 "ctgsyl.f"
			for (k = 1; k <= i__3; ++k) {
#line 540 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 540 "ctgsyl.f"
			    cscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 542 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 542 "ctgsyl.f"
			    cscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 544 "ctgsyl.f"
/* L80: */
#line 544 "ctgsyl.f"
			}
#line 545 "ctgsyl.f"
			i__3 = je;
#line 545 "ctgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 546 "ctgsyl.f"
			    i__4 = is - 1;
#line 546 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 546 "ctgsyl.f"
			    cscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 548 "ctgsyl.f"
			    i__4 = is - 1;
#line 548 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 548 "ctgsyl.f"
			    cscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 550 "ctgsyl.f"
/* L90: */
#line 550 "ctgsyl.f"
			}
#line 551 "ctgsyl.f"
			i__3 = je;
#line 551 "ctgsyl.f"
			for (k = js; k <= i__3; ++k) {
#line 552 "ctgsyl.f"
			    i__4 = *m - ie;
#line 552 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 552 "ctgsyl.f"
			    cscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &
				    c__1);
#line 554 "ctgsyl.f"
			    i__4 = *m - ie;
#line 554 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 554 "ctgsyl.f"
			    cscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &
				    c__1);
#line 556 "ctgsyl.f"
/* L100: */
#line 556 "ctgsyl.f"
			}
#line 557 "ctgsyl.f"
			i__3 = *n;
#line 557 "ctgsyl.f"
			for (k = je + 1; k <= i__3; ++k) {
#line 558 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 558 "ctgsyl.f"
			    cscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 560 "ctgsyl.f"
			    z__1.r = scaloc, z__1.i = 0.;
#line 560 "ctgsyl.f"
			    cscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 562 "ctgsyl.f"
/* L110: */
#line 562 "ctgsyl.f"
			}
#line 563 "ctgsyl.f"
			*scale *= scaloc;
#line 564 "ctgsyl.f"
		    }

/*                 Substitute R(I,J) and L(I,J) into remaining equation. */

#line 568 "ctgsyl.f"
		    if (i__ > 1) {
#line 569 "ctgsyl.f"
			i__3 = is - 1;
#line 569 "ctgsyl.f"
			cgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &a[is * 
				a_dim1 + 1], lda, &c__[is + js * c_dim1], ldc,
				 &c_b45, &c__[js * c_dim1 + 1], ldc, (ftnlen)
				1, (ftnlen)1);
#line 573 "ctgsyl.f"
			i__3 = is - 1;
#line 573 "ctgsyl.f"
			cgemm_("N", "N", &i__3, &nb, &mb, &c_b44, &d__[is * 
				d_dim1 + 1], ldd, &c__[is + js * c_dim1], ldc,
				 &c_b45, &f[js * f_dim1 + 1], ldf, (ftnlen)1, 
				(ftnlen)1);
#line 577 "ctgsyl.f"
		    }
#line 578 "ctgsyl.f"
		    if (j < q) {
#line 579 "ctgsyl.f"
			i__3 = *n - je;
#line 579 "ctgsyl.f"
			cgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js *
				 f_dim1], ldf, &b[js + (je + 1) * b_dim1], 
				ldb, &c_b45, &c__[is + (je + 1) * c_dim1], 
				ldc, (ftnlen)1, (ftnlen)1);
#line 583 "ctgsyl.f"
			i__3 = *n - je;
#line 583 "ctgsyl.f"
			cgemm_("N", "N", &mb, &i__3, &nb, &c_b45, &f[is + js *
				 f_dim1], ldf, &e[js + (je + 1) * e_dim1], 
				lde, &c_b45, &f[is + (je + 1) * f_dim1], ldf, 
				(ftnlen)1, (ftnlen)1);
#line 587 "ctgsyl.f"
		    }
#line 588 "ctgsyl.f"
/* L120: */
#line 588 "ctgsyl.f"
		}
#line 589 "ctgsyl.f"
/* L130: */
#line 589 "ctgsyl.f"
	    }
#line 590 "ctgsyl.f"
	    if (dscale != 0.) {
#line 591 "ctgsyl.f"
		if (*ijob == 1 || *ijob == 3) {
#line 592 "ctgsyl.f"
		    *dif = sqrt((doublereal) ((*m << 1) * *n)) / (dscale * 
			    sqrt(dsum));
#line 593 "ctgsyl.f"
		} else {
#line 594 "ctgsyl.f"
		    *dif = sqrt((doublereal) pq) / (dscale * sqrt(dsum));
#line 595 "ctgsyl.f"
		}
#line 596 "ctgsyl.f"
	    }
#line 597 "ctgsyl.f"
	    if (isolve == 2 && iround == 1) {
#line 598 "ctgsyl.f"
		if (notran) {
#line 599 "ctgsyl.f"
		    ifunc = *ijob;
#line 600 "ctgsyl.f"
		}
#line 601 "ctgsyl.f"
		scale2 = *scale;
#line 602 "ctgsyl.f"
		clacpy_("F", m, n, &c__[c_offset], ldc, &work[1], m, (ftnlen)
			1);
#line 603 "ctgsyl.f"
		clacpy_("F", m, n, &f[f_offset], ldf, &work[*m * *n + 1], m, (
			ftnlen)1);
#line 604 "ctgsyl.f"
		claset_("F", m, n, &c_b1, &c_b1, &c__[c_offset], ldc, (ftnlen)
			1);
#line 605 "ctgsyl.f"
		claset_("F", m, n, &c_b1, &c_b1, &f[f_offset], ldf, (ftnlen)1)
			;
#line 606 "ctgsyl.f"
	    } else if (isolve == 2 && iround == 2) {
#line 607 "ctgsyl.f"
		clacpy_("F", m, n, &work[1], m, &c__[c_offset], ldc, (ftnlen)
			1);
#line 608 "ctgsyl.f"
		clacpy_("F", m, n, &work[*m * *n + 1], m, &f[f_offset], ldf, (
			ftnlen)1);
#line 609 "ctgsyl.f"
		*scale = scale2;
#line 610 "ctgsyl.f"
	    }
#line 611 "ctgsyl.f"
/* L150: */
#line 611 "ctgsyl.f"
	}
#line 612 "ctgsyl.f"
    } else {

/*        Solve transposed (I, J)-subsystem */
/*            A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J) */
/*            R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J) */
/*        for I = 1,2,..., P; J = Q, Q-1,..., 1 */

#line 619 "ctgsyl.f"
	*scale = 1.;
#line 620 "ctgsyl.f"
	i__1 = p;
#line 620 "ctgsyl.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 621 "ctgsyl.f"
	    is = iwork[i__];
#line 622 "ctgsyl.f"
	    ie = iwork[i__ + 1] - 1;
#line 623 "ctgsyl.f"
	    mb = ie - is + 1;
#line 624 "ctgsyl.f"
	    i__2 = p + 2;
#line 624 "ctgsyl.f"
	    for (j = q; j >= i__2; --j) {
#line 625 "ctgsyl.f"
		js = iwork[j];
#line 626 "ctgsyl.f"
		je = iwork[j + 1] - 1;
#line 627 "ctgsyl.f"
		nb = je - js + 1;
#line 628 "ctgsyl.f"
		ctgsy2_(trans, &ifunc, &mb, &nb, &a[is + is * a_dim1], lda, &
			b[js + js * b_dim1], ldb, &c__[is + js * c_dim1], ldc,
			 &d__[is + is * d_dim1], ldd, &e[js + js * e_dim1], 
			lde, &f[is + js * f_dim1], ldf, &scaloc, &dsum, &
			dscale, &linfo, (ftnlen)1);
#line 633 "ctgsyl.f"
		if (linfo > 0) {
#line 633 "ctgsyl.f"
		    *info = linfo;
#line 633 "ctgsyl.f"
		}
#line 635 "ctgsyl.f"
		if (scaloc != 1.) {
#line 636 "ctgsyl.f"
		    i__3 = js - 1;
#line 636 "ctgsyl.f"
		    for (k = 1; k <= i__3; ++k) {
#line 637 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 637 "ctgsyl.f"
			cscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 639 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 639 "ctgsyl.f"
			cscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 641 "ctgsyl.f"
/* L160: */
#line 641 "ctgsyl.f"
		    }
#line 642 "ctgsyl.f"
		    i__3 = je;
#line 642 "ctgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 643 "ctgsyl.f"
			i__4 = is - 1;
#line 643 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 643 "ctgsyl.f"
			cscal_(&i__4, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 645 "ctgsyl.f"
			i__4 = is - 1;
#line 645 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 645 "ctgsyl.f"
			cscal_(&i__4, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 647 "ctgsyl.f"
/* L170: */
#line 647 "ctgsyl.f"
		    }
#line 648 "ctgsyl.f"
		    i__3 = je;
#line 648 "ctgsyl.f"
		    for (k = js; k <= i__3; ++k) {
#line 649 "ctgsyl.f"
			i__4 = *m - ie;
#line 649 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 649 "ctgsyl.f"
			cscal_(&i__4, &z__1, &c__[ie + 1 + k * c_dim1], &c__1)
				;
#line 651 "ctgsyl.f"
			i__4 = *m - ie;
#line 651 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 651 "ctgsyl.f"
			cscal_(&i__4, &z__1, &f[ie + 1 + k * f_dim1], &c__1);
#line 653 "ctgsyl.f"
/* L180: */
#line 653 "ctgsyl.f"
		    }
#line 654 "ctgsyl.f"
		    i__3 = *n;
#line 654 "ctgsyl.f"
		    for (k = je + 1; k <= i__3; ++k) {
#line 655 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 655 "ctgsyl.f"
			cscal_(m, &z__1, &c__[k * c_dim1 + 1], &c__1);
#line 657 "ctgsyl.f"
			z__1.r = scaloc, z__1.i = 0.;
#line 657 "ctgsyl.f"
			cscal_(m, &z__1, &f[k * f_dim1 + 1], &c__1);
#line 659 "ctgsyl.f"
/* L190: */
#line 659 "ctgsyl.f"
		    }
#line 660 "ctgsyl.f"
		    *scale *= scaloc;
#line 661 "ctgsyl.f"
		}

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 665 "ctgsyl.f"
		if (j > p + 2) {
#line 666 "ctgsyl.f"
		    i__3 = js - 1;
#line 666 "ctgsyl.f"
		    cgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &c__[is + js * 
			    c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &c_b45, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 670 "ctgsyl.f"
		    i__3 = js - 1;
#line 670 "ctgsyl.f"
		    cgemm_("N", "C", &mb, &i__3, &nb, &c_b45, &f[is + js * 
			    f_dim1], ldf, &e[js * e_dim1 + 1], lde, &c_b45, &
			    f[is + f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 674 "ctgsyl.f"
		}
#line 675 "ctgsyl.f"
		if (i__ < p) {
#line 676 "ctgsyl.f"
		    i__3 = *m - ie;
#line 676 "ctgsyl.f"
		    cgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &a[is + (ie + 1)
			     * a_dim1], lda, &c__[is + js * c_dim1], ldc, &
			    c_b45, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 680 "ctgsyl.f"
		    i__3 = *m - ie;
#line 680 "ctgsyl.f"
		    cgemm_("C", "N", &i__3, &nb, &mb, &c_b44, &d__[is + (ie + 
			    1) * d_dim1], ldd, &f[is + js * f_dim1], ldf, &
			    c_b45, &c__[ie + 1 + js * c_dim1], ldc, (ftnlen)1,
			     (ftnlen)1);
#line 684 "ctgsyl.f"
		}
#line 685 "ctgsyl.f"
/* L200: */
#line 685 "ctgsyl.f"
	    }
#line 686 "ctgsyl.f"
/* L210: */
#line 686 "ctgsyl.f"
	}
#line 687 "ctgsyl.f"
    }

#line 689 "ctgsyl.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 691 "ctgsyl.f"
    return 0;

/*     End of CTGSYL */

} /* ctgsyl_ */

