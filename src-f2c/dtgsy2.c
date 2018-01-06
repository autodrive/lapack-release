#line 1 "dtgsy2.f"
/* dtgsy2.f -- translated by f2c (version 20100827).
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

#line 1 "dtgsy2.f"
/* Table of constant values */

static integer c__8 = 8;
static integer c__1 = 1;
static doublereal c_b27 = -1.;
static doublereal c_b42 = 1.;
static doublereal c_b56 = 0.;

/* > \brief \b DTGSY2 solves the generalized Sylvester equation (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGSY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsy2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsy2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsy2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/*                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, */
/*                          IWORK, PQ, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N, */
/*      $                   PQ */
/*       DOUBLE PRECISION   RDSCAL, RDSUM, SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGSY2 solves the generalized Sylvester equation: */
/* > */
/* >             A * R - L * B = scale * C                (1) */
/* >             D * R - L * E = scale * F, */
/* > */
/* > using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices, */
/* > (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M, */
/* > N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E) */
/* > must be in generalized Schur canonical form, i.e. A, B are upper */
/* > quasi triangular and D, E are upper triangular. The solution (R, L) */
/* > overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor */
/* > chosen to avoid overflow. */
/* > */
/* > In matrix notation solving equation (1) corresponds to solve */
/* > Z*x = scale*b, where Z is defined as */
/* > */
/* >        Z = [ kron(In, A)  -kron(B**T, Im) ]             (2) */
/* >            [ kron(In, D)  -kron(E**T, Im) ], */
/* > */
/* > Ik is the identity matrix of size k and X**T is the transpose of X. */
/* > kron(X, Y) is the Kronecker product between the matrices X and Y. */
/* > In the process of solving (1), we solve a number of such systems */
/* > where Dim(In), Dim(In) = 1 or 2. */
/* > */
/* > If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y, */
/* > which is equivalent to solve for R and L in */
/* > */
/* >             A**T * R  + D**T * L   = scale * C           (3) */
/* >             R  * B**T + L  * E**T  = scale * -F */
/* > */
/* > This case is used to compute an estimate of Dif[(A, D), (B, E)] = */
/* > sigma_min(Z) using reverse communicaton with DLACON. */
/* > */
/* > DTGSY2 also (IJOB >= 1) contributes to the computation in DTGSYL */
/* > of an upper bound on the separation between to matrix pairs. Then */
/* > the input (A, D), (B, E) are sub-pencils of the matrix pair in */
/* > DTGSYL. See DTGSYL for details. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N', solve the generalized Sylvester equation (1). */
/* >          = 'T': solve the 'transposed' system (3). */
/* > \endverbatim */
/* > */
/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          Specifies what kind of functionality to be performed. */
/* >          = 0: solve (1) only. */
/* >          = 1: A contribution from this subsystem to a Frobenius */
/* >               norm-based estimate of the separation between two matrix */
/* >               pairs is computed. (look ahead strategy is used). */
/* >          = 2: A contribution from this subsystem to a Frobenius */
/* >               norm-based estimate of the separation between two matrix */
/* >               pairs is computed. (DGECON on sub-systems is used.) */
/* >          Not referenced if TRANS = 'T'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          On entry, M specifies the order of A and D, and the row */
/* >          dimension of C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          On entry, N specifies the order of B and E, and the column */
/* >          dimension of C, F, R and L. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA, M) */
/* >          On entry, A contains an upper quasi triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the matrix A. LDA >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, B contains an upper quasi triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the matrix B. LDB >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC, N) */
/* >          On entry, C contains the right-hand-side of the first matrix */
/* >          equation in (1). */
/* >          On exit, if IJOB = 0, C has been overwritten by the */
/* >          solution R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the matrix C. LDC >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (LDD, M) */
/* >          On entry, D contains an upper triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDD */
/* > \verbatim */
/* >          LDD is INTEGER */
/* >          The leading dimension of the matrix D. LDD >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (LDE, N) */
/* >          On entry, E contains an upper triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDE */
/* > \verbatim */
/* >          LDE is INTEGER */
/* >          The leading dimension of the matrix E. LDE >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* >          F is DOUBLE PRECISION array, dimension (LDF, N) */
/* >          On entry, F contains the right-hand-side of the second matrix */
/* >          equation in (1). */
/* >          On exit, if IJOB = 0, F has been overwritten by the */
/* >          solution L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* >          LDF is INTEGER */
/* >          The leading dimension of the matrix F. LDF >= max(1, M). */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* >          SCALE is DOUBLE PRECISION */
/* >          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions */
/* >          R and L (C and F on entry) will hold the solutions to a */
/* >          slightly perturbed system but the input matrices A, B, D and */
/* >          E have not been changed. If SCALE = 0, R and L will hold the */
/* >          solutions to the homogeneous system with C = F = 0. Normally, */
/* >          SCALE = 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSUM */
/* > \verbatim */
/* >          RDSUM is DOUBLE PRECISION */
/* >          On entry, the sum of squares of computed contributions to */
/* >          the Dif-estimate under computation by DTGSYL, where the */
/* >          scaling factor RDSCAL (see below) has been factored out. */
/* >          On exit, the corresponding sum of squares updated with the */
/* >          contributions from the current sub-system. */
/* >          If TRANS = 'T' RDSUM is not touched. */
/* >          NOTE: RDSUM only makes sense when DTGSY2 is called by DTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* >          RDSCAL is DOUBLE PRECISION */
/* >          On entry, scaling factor used to prevent overflow in RDSUM. */
/* >          On exit, RDSCAL is updated w.r.t. the current contributions */
/* >          in RDSUM. */
/* >          If TRANS = 'T', RDSCAL is not touched. */
/* >          NOTE: RDSCAL only makes sense when DTGSY2 is called by */
/* >                DTGSYL. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (M+N+2) */
/* > \endverbatim */
/* > */
/* > \param[out] PQ */
/* > \verbatim */
/* >          PQ is INTEGER */
/* >          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and */
/* >          8-by-8) solved by this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          On exit, if INFO is set to */
/* >            =0: Successful exit */
/* >            <0: If INFO = -i, the i-th argument had an illegal value. */
/* >            >0: The matrix pairs (A, D) and (B, E) have common or very */
/* >                close eigenvalues. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleSYauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int dtgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer 
	*pq, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, p, q;
    static doublereal z__[64]	/* was [8][8] */;
    static integer ie, je, mb, nb, ii, jj, is, js;
    static doublereal rhs[8];
    static integer isp1, jsp1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr, zdim, ipiv[8], jpiv[8];
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    , dgesc2_(integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *), dgetc2_(integer *, 
	    doublereal *, integer *, integer *, integer *, integer *), 
	    dlatdf_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical notran;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*  Replaced various illegal calls to DCOPY by calls to DLASET. */
/*  Sven Hammarling, 27/5/02. */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test input parameters */

#line 329 "dtgsy2.f"
    /* Parameter adjustments */
#line 329 "dtgsy2.f"
    a_dim1 = *lda;
#line 329 "dtgsy2.f"
    a_offset = 1 + a_dim1;
#line 329 "dtgsy2.f"
    a -= a_offset;
#line 329 "dtgsy2.f"
    b_dim1 = *ldb;
#line 329 "dtgsy2.f"
    b_offset = 1 + b_dim1;
#line 329 "dtgsy2.f"
    b -= b_offset;
#line 329 "dtgsy2.f"
    c_dim1 = *ldc;
#line 329 "dtgsy2.f"
    c_offset = 1 + c_dim1;
#line 329 "dtgsy2.f"
    c__ -= c_offset;
#line 329 "dtgsy2.f"
    d_dim1 = *ldd;
#line 329 "dtgsy2.f"
    d_offset = 1 + d_dim1;
#line 329 "dtgsy2.f"
    d__ -= d_offset;
#line 329 "dtgsy2.f"
    e_dim1 = *lde;
#line 329 "dtgsy2.f"
    e_offset = 1 + e_dim1;
#line 329 "dtgsy2.f"
    e -= e_offset;
#line 329 "dtgsy2.f"
    f_dim1 = *ldf;
#line 329 "dtgsy2.f"
    f_offset = 1 + f_dim1;
#line 329 "dtgsy2.f"
    f -= f_offset;
#line 329 "dtgsy2.f"
    --iwork;
#line 329 "dtgsy2.f"

#line 329 "dtgsy2.f"
    /* Function Body */
#line 329 "dtgsy2.f"
    *info = 0;
#line 330 "dtgsy2.f"
    ierr = 0;
#line 331 "dtgsy2.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 332 "dtgsy2.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 333 "dtgsy2.f"
	*info = -1;
#line 334 "dtgsy2.f"
    } else if (notran) {
#line 335 "dtgsy2.f"
	if (*ijob < 0 || *ijob > 2) {
#line 336 "dtgsy2.f"
	    *info = -2;
#line 337 "dtgsy2.f"
	}
#line 338 "dtgsy2.f"
    }
#line 339 "dtgsy2.f"
    if (*info == 0) {
#line 340 "dtgsy2.f"
	if (*m <= 0) {
#line 341 "dtgsy2.f"
	    *info = -3;
#line 342 "dtgsy2.f"
	} else if (*n <= 0) {
#line 343 "dtgsy2.f"
	    *info = -4;
#line 344 "dtgsy2.f"
	} else if (*lda < max(1,*m)) {
#line 345 "dtgsy2.f"
	    *info = -6;
#line 346 "dtgsy2.f"
	} else if (*ldb < max(1,*n)) {
#line 347 "dtgsy2.f"
	    *info = -8;
#line 348 "dtgsy2.f"
	} else if (*ldc < max(1,*m)) {
#line 349 "dtgsy2.f"
	    *info = -10;
#line 350 "dtgsy2.f"
	} else if (*ldd < max(1,*m)) {
#line 351 "dtgsy2.f"
	    *info = -12;
#line 352 "dtgsy2.f"
	} else if (*lde < max(1,*n)) {
#line 353 "dtgsy2.f"
	    *info = -14;
#line 354 "dtgsy2.f"
	} else if (*ldf < max(1,*m)) {
#line 355 "dtgsy2.f"
	    *info = -16;
#line 356 "dtgsy2.f"
	}
#line 357 "dtgsy2.f"
    }
#line 358 "dtgsy2.f"
    if (*info != 0) {
#line 359 "dtgsy2.f"
	i__1 = -(*info);
#line 359 "dtgsy2.f"
	xerbla_("DTGSY2", &i__1, (ftnlen)6);
#line 360 "dtgsy2.f"
	return 0;
#line 361 "dtgsy2.f"
    }

/*     Determine block structure of A */

#line 365 "dtgsy2.f"
    *pq = 0;
#line 366 "dtgsy2.f"
    p = 0;
#line 367 "dtgsy2.f"
    i__ = 1;
#line 368 "dtgsy2.f"
L10:
#line 369 "dtgsy2.f"
    if (i__ > *m) {
#line 369 "dtgsy2.f"
	goto L20;
#line 369 "dtgsy2.f"
    }
#line 371 "dtgsy2.f"
    ++p;
#line 372 "dtgsy2.f"
    iwork[p] = i__;
#line 373 "dtgsy2.f"
    if (i__ == *m) {
#line 373 "dtgsy2.f"
	goto L20;
#line 373 "dtgsy2.f"
    }
#line 375 "dtgsy2.f"
    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 376 "dtgsy2.f"
	i__ += 2;
#line 377 "dtgsy2.f"
    } else {
#line 378 "dtgsy2.f"
	++i__;
#line 379 "dtgsy2.f"
    }
#line 380 "dtgsy2.f"
    goto L10;
#line 381 "dtgsy2.f"
L20:
#line 382 "dtgsy2.f"
    iwork[p + 1] = *m + 1;

/*     Determine block structure of B */

#line 386 "dtgsy2.f"
    q = p + 1;
#line 387 "dtgsy2.f"
    j = 1;
#line 388 "dtgsy2.f"
L30:
#line 389 "dtgsy2.f"
    if (j > *n) {
#line 389 "dtgsy2.f"
	goto L40;
#line 389 "dtgsy2.f"
    }
#line 391 "dtgsy2.f"
    ++q;
#line 392 "dtgsy2.f"
    iwork[q] = j;
#line 393 "dtgsy2.f"
    if (j == *n) {
#line 393 "dtgsy2.f"
	goto L40;
#line 393 "dtgsy2.f"
    }
#line 395 "dtgsy2.f"
    if (b[j + 1 + j * b_dim1] != 0.) {
#line 396 "dtgsy2.f"
	j += 2;
#line 397 "dtgsy2.f"
    } else {
#line 398 "dtgsy2.f"
	++j;
#line 399 "dtgsy2.f"
    }
#line 400 "dtgsy2.f"
    goto L30;
#line 401 "dtgsy2.f"
L40:
#line 402 "dtgsy2.f"
    iwork[q + 1] = *n + 1;
#line 403 "dtgsy2.f"
    *pq = p * (q - p - 1);

#line 405 "dtgsy2.f"
    if (notran) {

/*        Solve (I, J) - subsystem */
/*           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*        for I = P, P - 1, ..., 1; J = 1, 2, ..., Q */

#line 412 "dtgsy2.f"
	*scale = 1.;
#line 413 "dtgsy2.f"
	scaloc = 1.;
#line 414 "dtgsy2.f"
	i__1 = q;
#line 414 "dtgsy2.f"
	for (j = p + 2; j <= i__1; ++j) {
#line 415 "dtgsy2.f"
	    js = iwork[j];
#line 416 "dtgsy2.f"
	    jsp1 = js + 1;
#line 417 "dtgsy2.f"
	    je = iwork[j + 1] - 1;
#line 418 "dtgsy2.f"
	    nb = je - js + 1;
#line 419 "dtgsy2.f"
	    for (i__ = p; i__ >= 1; --i__) {

#line 421 "dtgsy2.f"
		is = iwork[i__];
#line 422 "dtgsy2.f"
		isp1 = is + 1;
#line 423 "dtgsy2.f"
		ie = iwork[i__ + 1] - 1;
#line 424 "dtgsy2.f"
		mb = ie - is + 1;
#line 425 "dtgsy2.f"
		zdim = mb * nb << 1;

#line 427 "dtgsy2.f"
		if (mb == 1 && nb == 1) {

/*                 Build a 2-by-2 system Z * x = RHS */

#line 431 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 432 "dtgsy2.f"
		    z__[1] = d__[is + is * d_dim1];
#line 433 "dtgsy2.f"
		    z__[8] = -b[js + js * b_dim1];
#line 434 "dtgsy2.f"
		    z__[9] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 438 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 439 "dtgsy2.f"
		    rhs[1] = f[is + js * f_dim1];

/*                 Solve Z * x = RHS */

#line 443 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 444 "dtgsy2.f"
		    if (ierr > 0) {
#line 444 "dtgsy2.f"
			*info = ierr;
#line 444 "dtgsy2.f"
		    }

#line 447 "dtgsy2.f"
		    if (*ijob == 0) {
#line 448 "dtgsy2.f"
			dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 450 "dtgsy2.f"
			if (scaloc != 1.) {
#line 451 "dtgsy2.f"
			    i__2 = *n;
#line 451 "dtgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 452 "dtgsy2.f"
				dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 453 "dtgsy2.f"
				dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 454 "dtgsy2.f"
/* L50: */
#line 454 "dtgsy2.f"
			    }
#line 455 "dtgsy2.f"
			    *scale *= scaloc;
#line 456 "dtgsy2.f"
			}
#line 457 "dtgsy2.f"
		    } else {
#line 458 "dtgsy2.f"
			dlatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 460 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 464 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 465 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[1];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 470 "dtgsy2.f"
		    if (i__ > 1) {
#line 471 "dtgsy2.f"
			alpha = -rhs[0];
#line 472 "dtgsy2.f"
			i__2 = is - 1;
#line 472 "dtgsy2.f"
			daxpy_(&i__2, &alpha, &a[is * a_dim1 + 1], &c__1, &
				c__[js * c_dim1 + 1], &c__1);
#line 474 "dtgsy2.f"
			i__2 = is - 1;
#line 474 "dtgsy2.f"
			daxpy_(&i__2, &alpha, &d__[is * d_dim1 + 1], &c__1, &
				f[js * f_dim1 + 1], &c__1);
#line 476 "dtgsy2.f"
		    }
#line 477 "dtgsy2.f"
		    if (j < q) {
#line 478 "dtgsy2.f"
			i__2 = *n - je;
#line 478 "dtgsy2.f"
			daxpy_(&i__2, &rhs[1], &b[js + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 480 "dtgsy2.f"
			i__2 = *n - je;
#line 480 "dtgsy2.f"
			daxpy_(&i__2, &rhs[1], &e[js + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 482 "dtgsy2.f"
		    }

#line 484 "dtgsy2.f"
		} else if (mb == 1 && nb == 2) {

/*                 Build a 4-by-4 system Z * x = RHS */

#line 488 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 489 "dtgsy2.f"
		    z__[1] = 0.;
#line 490 "dtgsy2.f"
		    z__[2] = d__[is + is * d_dim1];
#line 491 "dtgsy2.f"
		    z__[3] = 0.;

#line 493 "dtgsy2.f"
		    z__[8] = 0.;
#line 494 "dtgsy2.f"
		    z__[9] = a[is + is * a_dim1];
#line 495 "dtgsy2.f"
		    z__[10] = 0.;
#line 496 "dtgsy2.f"
		    z__[11] = d__[is + is * d_dim1];

#line 498 "dtgsy2.f"
		    z__[16] = -b[js + js * b_dim1];
#line 499 "dtgsy2.f"
		    z__[17] = -b[js + jsp1 * b_dim1];
#line 500 "dtgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 501 "dtgsy2.f"
		    z__[19] = -e[js + jsp1 * e_dim1];

#line 503 "dtgsy2.f"
		    z__[24] = -b[jsp1 + js * b_dim1];
#line 504 "dtgsy2.f"
		    z__[25] = -b[jsp1 + jsp1 * b_dim1];
#line 505 "dtgsy2.f"
		    z__[26] = 0.;
#line 506 "dtgsy2.f"
		    z__[27] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 510 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 511 "dtgsy2.f"
		    rhs[1] = c__[is + jsp1 * c_dim1];
#line 512 "dtgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 513 "dtgsy2.f"
		    rhs[3] = f[is + jsp1 * f_dim1];

/*                 Solve Z * x = RHS */

#line 517 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 518 "dtgsy2.f"
		    if (ierr > 0) {
#line 518 "dtgsy2.f"
			*info = ierr;
#line 518 "dtgsy2.f"
		    }

#line 521 "dtgsy2.f"
		    if (*ijob == 0) {
#line 522 "dtgsy2.f"
			dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 524 "dtgsy2.f"
			if (scaloc != 1.) {
#line 525 "dtgsy2.f"
			    i__2 = *n;
#line 525 "dtgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 526 "dtgsy2.f"
				dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 527 "dtgsy2.f"
				dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 528 "dtgsy2.f"
/* L60: */
#line 528 "dtgsy2.f"
			    }
#line 529 "dtgsy2.f"
			    *scale *= scaloc;
#line 530 "dtgsy2.f"
			}
#line 531 "dtgsy2.f"
		    } else {
#line 532 "dtgsy2.f"
			dlatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 534 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 538 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 539 "dtgsy2.f"
		    c__[is + jsp1 * c_dim1] = rhs[1];
#line 540 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 541 "dtgsy2.f"
		    f[is + jsp1 * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 546 "dtgsy2.f"
		    if (i__ > 1) {
#line 547 "dtgsy2.f"
			i__2 = is - 1;
#line 547 "dtgsy2.f"
			dger_(&i__2, &nb, &c_b27, &a[is * a_dim1 + 1], &c__1, 
				rhs, &c__1, &c__[js * c_dim1 + 1], ldc);
#line 549 "dtgsy2.f"
			i__2 = is - 1;
#line 549 "dtgsy2.f"
			dger_(&i__2, &nb, &c_b27, &d__[is * d_dim1 + 1], &
				c__1, rhs, &c__1, &f[js * f_dim1 + 1], ldf);
#line 551 "dtgsy2.f"
		    }
#line 552 "dtgsy2.f"
		    if (j < q) {
#line 553 "dtgsy2.f"
			i__2 = *n - je;
#line 553 "dtgsy2.f"
			daxpy_(&i__2, &rhs[2], &b[js + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 555 "dtgsy2.f"
			i__2 = *n - je;
#line 555 "dtgsy2.f"
			daxpy_(&i__2, &rhs[2], &e[js + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 557 "dtgsy2.f"
			i__2 = *n - je;
#line 557 "dtgsy2.f"
			daxpy_(&i__2, &rhs[3], &b[jsp1 + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 559 "dtgsy2.f"
			i__2 = *n - je;
#line 559 "dtgsy2.f"
			daxpy_(&i__2, &rhs[3], &e[jsp1 + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 561 "dtgsy2.f"
		    }

#line 563 "dtgsy2.f"
		} else if (mb == 2 && nb == 1) {

/*                 Build a 4-by-4 system Z * x = RHS */

#line 567 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 568 "dtgsy2.f"
		    z__[1] = a[isp1 + is * a_dim1];
#line 569 "dtgsy2.f"
		    z__[2] = d__[is + is * d_dim1];
#line 570 "dtgsy2.f"
		    z__[3] = 0.;

#line 572 "dtgsy2.f"
		    z__[8] = a[is + isp1 * a_dim1];
#line 573 "dtgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 574 "dtgsy2.f"
		    z__[10] = d__[is + isp1 * d_dim1];
#line 575 "dtgsy2.f"
		    z__[11] = d__[isp1 + isp1 * d_dim1];

#line 577 "dtgsy2.f"
		    z__[16] = -b[js + js * b_dim1];
#line 578 "dtgsy2.f"
		    z__[17] = 0.;
#line 579 "dtgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 580 "dtgsy2.f"
		    z__[19] = 0.;

#line 582 "dtgsy2.f"
		    z__[24] = 0.;
#line 583 "dtgsy2.f"
		    z__[25] = -b[js + js * b_dim1];
#line 584 "dtgsy2.f"
		    z__[26] = 0.;
#line 585 "dtgsy2.f"
		    z__[27] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 589 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 590 "dtgsy2.f"
		    rhs[1] = c__[isp1 + js * c_dim1];
#line 591 "dtgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 592 "dtgsy2.f"
		    rhs[3] = f[isp1 + js * f_dim1];

/*                 Solve Z * x = RHS */

#line 596 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 597 "dtgsy2.f"
		    if (ierr > 0) {
#line 597 "dtgsy2.f"
			*info = ierr;
#line 597 "dtgsy2.f"
		    }
#line 599 "dtgsy2.f"
		    if (*ijob == 0) {
#line 600 "dtgsy2.f"
			dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 602 "dtgsy2.f"
			if (scaloc != 1.) {
#line 603 "dtgsy2.f"
			    i__2 = *n;
#line 603 "dtgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 604 "dtgsy2.f"
				dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 605 "dtgsy2.f"
				dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 606 "dtgsy2.f"
/* L70: */
#line 606 "dtgsy2.f"
			    }
#line 607 "dtgsy2.f"
			    *scale *= scaloc;
#line 608 "dtgsy2.f"
			}
#line 609 "dtgsy2.f"
		    } else {
#line 610 "dtgsy2.f"
			dlatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 612 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 616 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 617 "dtgsy2.f"
		    c__[isp1 + js * c_dim1] = rhs[1];
#line 618 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 619 "dtgsy2.f"
		    f[isp1 + js * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 624 "dtgsy2.f"
		    if (i__ > 1) {
#line 625 "dtgsy2.f"
			i__2 = is - 1;
#line 625 "dtgsy2.f"
			dgemv_("N", &i__2, &mb, &c_b27, &a[is * a_dim1 + 1], 
				lda, rhs, &c__1, &c_b42, &c__[js * c_dim1 + 1]
				, &c__1, (ftnlen)1);
#line 627 "dtgsy2.f"
			i__2 = is - 1;
#line 627 "dtgsy2.f"
			dgemv_("N", &i__2, &mb, &c_b27, &d__[is * d_dim1 + 1],
				 ldd, rhs, &c__1, &c_b42, &f[js * f_dim1 + 1],
				 &c__1, (ftnlen)1);
#line 629 "dtgsy2.f"
		    }
#line 630 "dtgsy2.f"
		    if (j < q) {
#line 631 "dtgsy2.f"
			i__2 = *n - je;
#line 631 "dtgsy2.f"
			dger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &b[js + (je 
				+ 1) * b_dim1], ldb, &c__[is + (je + 1) * 
				c_dim1], ldc);
#line 633 "dtgsy2.f"
			i__2 = *n - je;
#line 633 "dtgsy2.f"
			dger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &e[js + (je 
				+ 1) * e_dim1], lde, &f[is + (je + 1) * 
				f_dim1], ldf);
#line 635 "dtgsy2.f"
		    }

#line 637 "dtgsy2.f"
		} else if (mb == 2 && nb == 2) {

/*                 Build an 8-by-8 system Z * x = RHS */

#line 641 "dtgsy2.f"
		    dlaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8, (
			    ftnlen)1);

#line 643 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 644 "dtgsy2.f"
		    z__[1] = a[isp1 + is * a_dim1];
#line 645 "dtgsy2.f"
		    z__[4] = d__[is + is * d_dim1];

#line 647 "dtgsy2.f"
		    z__[8] = a[is + isp1 * a_dim1];
#line 648 "dtgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 649 "dtgsy2.f"
		    z__[12] = d__[is + isp1 * d_dim1];
#line 650 "dtgsy2.f"
		    z__[13] = d__[isp1 + isp1 * d_dim1];

#line 652 "dtgsy2.f"
		    z__[18] = a[is + is * a_dim1];
#line 653 "dtgsy2.f"
		    z__[19] = a[isp1 + is * a_dim1];
#line 654 "dtgsy2.f"
		    z__[22] = d__[is + is * d_dim1];

#line 656 "dtgsy2.f"
		    z__[26] = a[is + isp1 * a_dim1];
#line 657 "dtgsy2.f"
		    z__[27] = a[isp1 + isp1 * a_dim1];
#line 658 "dtgsy2.f"
		    z__[30] = d__[is + isp1 * d_dim1];
#line 659 "dtgsy2.f"
		    z__[31] = d__[isp1 + isp1 * d_dim1];

#line 661 "dtgsy2.f"
		    z__[32] = -b[js + js * b_dim1];
#line 662 "dtgsy2.f"
		    z__[34] = -b[js + jsp1 * b_dim1];
#line 663 "dtgsy2.f"
		    z__[36] = -e[js + js * e_dim1];
#line 664 "dtgsy2.f"
		    z__[38] = -e[js + jsp1 * e_dim1];

#line 666 "dtgsy2.f"
		    z__[41] = -b[js + js * b_dim1];
#line 667 "dtgsy2.f"
		    z__[43] = -b[js + jsp1 * b_dim1];
#line 668 "dtgsy2.f"
		    z__[45] = -e[js + js * e_dim1];
#line 669 "dtgsy2.f"
		    z__[47] = -e[js + jsp1 * e_dim1];

#line 671 "dtgsy2.f"
		    z__[48] = -b[jsp1 + js * b_dim1];
#line 672 "dtgsy2.f"
		    z__[50] = -b[jsp1 + jsp1 * b_dim1];
#line 673 "dtgsy2.f"
		    z__[54] = -e[jsp1 + jsp1 * e_dim1];

#line 675 "dtgsy2.f"
		    z__[57] = -b[jsp1 + js * b_dim1];
#line 676 "dtgsy2.f"
		    z__[59] = -b[jsp1 + jsp1 * b_dim1];
#line 677 "dtgsy2.f"
		    z__[63] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 681 "dtgsy2.f"
		    k = 1;
#line 682 "dtgsy2.f"
		    ii = mb * nb + 1;
#line 683 "dtgsy2.f"
		    i__2 = nb - 1;
#line 683 "dtgsy2.f"
		    for (jj = 0; jj <= i__2; ++jj) {
#line 684 "dtgsy2.f"
			dcopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &
				rhs[k - 1], &c__1);
#line 685 "dtgsy2.f"
			dcopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[
				ii - 1], &c__1);
#line 686 "dtgsy2.f"
			k += mb;
#line 687 "dtgsy2.f"
			ii += mb;
#line 688 "dtgsy2.f"
/* L80: */
#line 688 "dtgsy2.f"
		    }

/*                 Solve Z * x = RHS */

#line 692 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 693 "dtgsy2.f"
		    if (ierr > 0) {
#line 693 "dtgsy2.f"
			*info = ierr;
#line 693 "dtgsy2.f"
		    }
#line 695 "dtgsy2.f"
		    if (*ijob == 0) {
#line 696 "dtgsy2.f"
			dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 698 "dtgsy2.f"
			if (scaloc != 1.) {
#line 699 "dtgsy2.f"
			    i__2 = *n;
#line 699 "dtgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 700 "dtgsy2.f"
				dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 701 "dtgsy2.f"
				dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 702 "dtgsy2.f"
/* L90: */
#line 702 "dtgsy2.f"
			    }
#line 703 "dtgsy2.f"
			    *scale *= scaloc;
#line 704 "dtgsy2.f"
			}
#line 705 "dtgsy2.f"
		    } else {
#line 706 "dtgsy2.f"
			dlatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 708 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 712 "dtgsy2.f"
		    k = 1;
#line 713 "dtgsy2.f"
		    ii = mb * nb + 1;
#line 714 "dtgsy2.f"
		    i__2 = nb - 1;
#line 714 "dtgsy2.f"
		    for (jj = 0; jj <= i__2; ++jj) {
#line 715 "dtgsy2.f"
			dcopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
				c_dim1], &c__1);
#line 716 "dtgsy2.f"
			dcopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
				f_dim1], &c__1);
#line 717 "dtgsy2.f"
			k += mb;
#line 718 "dtgsy2.f"
			ii += mb;
#line 719 "dtgsy2.f"
/* L100: */
#line 719 "dtgsy2.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 724 "dtgsy2.f"
		    if (i__ > 1) {
#line 725 "dtgsy2.f"
			i__2 = is - 1;
#line 725 "dtgsy2.f"
			dgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &a[is * 
				a_dim1 + 1], lda, rhs, &mb, &c_b42, &c__[js * 
				c_dim1 + 1], ldc, (ftnlen)1, (ftnlen)1);
#line 728 "dtgsy2.f"
			i__2 = is - 1;
#line 728 "dtgsy2.f"
			dgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &d__[is * 
				d_dim1 + 1], ldd, rhs, &mb, &c_b42, &f[js * 
				f_dim1 + 1], ldf, (ftnlen)1, (ftnlen)1);
#line 731 "dtgsy2.f"
		    }
#line 732 "dtgsy2.f"
		    if (j < q) {
#line 733 "dtgsy2.f"
			k = mb * nb + 1;
#line 734 "dtgsy2.f"
			i__2 = *n - je;
#line 734 "dtgsy2.f"
			dgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1],
				 &mb, &b[js + (je + 1) * b_dim1], ldb, &c_b42,
				 &c__[is + (je + 1) * c_dim1], ldc, (ftnlen)1,
				 (ftnlen)1);
#line 737 "dtgsy2.f"
			i__2 = *n - je;
#line 737 "dtgsy2.f"
			dgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1],
				 &mb, &e[js + (je + 1) * e_dim1], lde, &c_b42,
				 &f[is + (je + 1) * f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 740 "dtgsy2.f"
		    }

#line 742 "dtgsy2.f"
		}

#line 744 "dtgsy2.f"
/* L110: */
#line 744 "dtgsy2.f"
	    }
#line 745 "dtgsy2.f"
/* L120: */
#line 745 "dtgsy2.f"
	}
#line 746 "dtgsy2.f"
    } else {

/*        Solve (I, J) - subsystem */
/*             A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J)  =  C(I, J) */
/*             R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J) */
/*        for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1 */

#line 753 "dtgsy2.f"
	*scale = 1.;
#line 754 "dtgsy2.f"
	scaloc = 1.;
#line 755 "dtgsy2.f"
	i__1 = p;
#line 755 "dtgsy2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 757 "dtgsy2.f"
	    is = iwork[i__];
#line 758 "dtgsy2.f"
	    isp1 = is + 1;
#line 759 "dtgsy2.f"
	    ie = iwork[i__ + 1] - 1;
#line 760 "dtgsy2.f"
	    mb = ie - is + 1;
#line 761 "dtgsy2.f"
	    i__2 = p + 2;
#line 761 "dtgsy2.f"
	    for (j = q; j >= i__2; --j) {

#line 763 "dtgsy2.f"
		js = iwork[j];
#line 764 "dtgsy2.f"
		jsp1 = js + 1;
#line 765 "dtgsy2.f"
		je = iwork[j + 1] - 1;
#line 766 "dtgsy2.f"
		nb = je - js + 1;
#line 767 "dtgsy2.f"
		zdim = mb * nb << 1;
#line 768 "dtgsy2.f"
		if (mb == 1 && nb == 1) {

/*                 Build a 2-by-2 system Z**T * x = RHS */

#line 772 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 773 "dtgsy2.f"
		    z__[1] = -b[js + js * b_dim1];
#line 774 "dtgsy2.f"
		    z__[8] = d__[is + is * d_dim1];
#line 775 "dtgsy2.f"
		    z__[9] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 779 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 780 "dtgsy2.f"
		    rhs[1] = f[is + js * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 784 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 785 "dtgsy2.f"
		    if (ierr > 0) {
#line 785 "dtgsy2.f"
			*info = ierr;
#line 785 "dtgsy2.f"
		    }

#line 788 "dtgsy2.f"
		    dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 789 "dtgsy2.f"
		    if (scaloc != 1.) {
#line 790 "dtgsy2.f"
			i__3 = *n;
#line 790 "dtgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 791 "dtgsy2.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 792 "dtgsy2.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 793 "dtgsy2.f"
/* L130: */
#line 793 "dtgsy2.f"
			}
#line 794 "dtgsy2.f"
			*scale *= scaloc;
#line 795 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 799 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 800 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[1];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 805 "dtgsy2.f"
		    if (j > p + 2) {
#line 806 "dtgsy2.f"
			alpha = rhs[0];
#line 807 "dtgsy2.f"
			i__3 = js - 1;
#line 807 "dtgsy2.f"
			daxpy_(&i__3, &alpha, &b[js * b_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 809 "dtgsy2.f"
			alpha = rhs[1];
#line 810 "dtgsy2.f"
			i__3 = js - 1;
#line 810 "dtgsy2.f"
			daxpy_(&i__3, &alpha, &e[js * e_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 812 "dtgsy2.f"
		    }
#line 813 "dtgsy2.f"
		    if (i__ < p) {
#line 814 "dtgsy2.f"
			alpha = -rhs[0];
#line 815 "dtgsy2.f"
			i__3 = *m - ie;
#line 815 "dtgsy2.f"
			daxpy_(&i__3, &alpha, &a[is + (ie + 1) * a_dim1], lda,
				 &c__[ie + 1 + js * c_dim1], &c__1);
#line 817 "dtgsy2.f"
			alpha = -rhs[1];
#line 818 "dtgsy2.f"
			i__3 = *m - ie;
#line 818 "dtgsy2.f"
			daxpy_(&i__3, &alpha, &d__[is + (ie + 1) * d_dim1], 
				ldd, &c__[ie + 1 + js * c_dim1], &c__1);
#line 820 "dtgsy2.f"
		    }

#line 822 "dtgsy2.f"
		} else if (mb == 1 && nb == 2) {

/*                 Build a 4-by-4 system Z**T * x = RHS */

#line 826 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 827 "dtgsy2.f"
		    z__[1] = 0.;
#line 828 "dtgsy2.f"
		    z__[2] = -b[js + js * b_dim1];
#line 829 "dtgsy2.f"
		    z__[3] = -b[jsp1 + js * b_dim1];

#line 831 "dtgsy2.f"
		    z__[8] = 0.;
#line 832 "dtgsy2.f"
		    z__[9] = a[is + is * a_dim1];
#line 833 "dtgsy2.f"
		    z__[10] = -b[js + jsp1 * b_dim1];
#line 834 "dtgsy2.f"
		    z__[11] = -b[jsp1 + jsp1 * b_dim1];

#line 836 "dtgsy2.f"
		    z__[16] = d__[is + is * d_dim1];
#line 837 "dtgsy2.f"
		    z__[17] = 0.;
#line 838 "dtgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 839 "dtgsy2.f"
		    z__[19] = 0.;

#line 841 "dtgsy2.f"
		    z__[24] = 0.;
#line 842 "dtgsy2.f"
		    z__[25] = d__[is + is * d_dim1];
#line 843 "dtgsy2.f"
		    z__[26] = -e[js + jsp1 * e_dim1];
#line 844 "dtgsy2.f"
		    z__[27] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 848 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 849 "dtgsy2.f"
		    rhs[1] = c__[is + jsp1 * c_dim1];
#line 850 "dtgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 851 "dtgsy2.f"
		    rhs[3] = f[is + jsp1 * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 855 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 856 "dtgsy2.f"
		    if (ierr > 0) {
#line 856 "dtgsy2.f"
			*info = ierr;
#line 856 "dtgsy2.f"
		    }
#line 858 "dtgsy2.f"
		    dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 859 "dtgsy2.f"
		    if (scaloc != 1.) {
#line 860 "dtgsy2.f"
			i__3 = *n;
#line 860 "dtgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 861 "dtgsy2.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 862 "dtgsy2.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 863 "dtgsy2.f"
/* L140: */
#line 863 "dtgsy2.f"
			}
#line 864 "dtgsy2.f"
			*scale *= scaloc;
#line 865 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 869 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 870 "dtgsy2.f"
		    c__[is + jsp1 * c_dim1] = rhs[1];
#line 871 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 872 "dtgsy2.f"
		    f[is + jsp1 * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 877 "dtgsy2.f"
		    if (j > p + 2) {
#line 878 "dtgsy2.f"
			i__3 = js - 1;
#line 878 "dtgsy2.f"
			daxpy_(&i__3, rhs, &b[js * b_dim1 + 1], &c__1, &f[is 
				+ f_dim1], ldf);
#line 880 "dtgsy2.f"
			i__3 = js - 1;
#line 880 "dtgsy2.f"
			daxpy_(&i__3, &rhs[1], &b[jsp1 * b_dim1 + 1], &c__1, &
				f[is + f_dim1], ldf);
#line 882 "dtgsy2.f"
			i__3 = js - 1;
#line 882 "dtgsy2.f"
			daxpy_(&i__3, &rhs[2], &e[js * e_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 884 "dtgsy2.f"
			i__3 = js - 1;
#line 884 "dtgsy2.f"
			daxpy_(&i__3, &rhs[3], &e[jsp1 * e_dim1 + 1], &c__1, &
				f[is + f_dim1], ldf);
#line 886 "dtgsy2.f"
		    }
#line 887 "dtgsy2.f"
		    if (i__ < p) {
#line 888 "dtgsy2.f"
			i__3 = *m - ie;
#line 888 "dtgsy2.f"
			dger_(&i__3, &nb, &c_b27, &a[is + (ie + 1) * a_dim1], 
				lda, rhs, &c__1, &c__[ie + 1 + js * c_dim1], 
				ldc);
#line 890 "dtgsy2.f"
			i__3 = *m - ie;
#line 890 "dtgsy2.f"
			dger_(&i__3, &nb, &c_b27, &d__[is + (ie + 1) * d_dim1]
				, ldd, &rhs[2], &c__1, &c__[ie + 1 + js * 
				c_dim1], ldc);
#line 892 "dtgsy2.f"
		    }

#line 894 "dtgsy2.f"
		} else if (mb == 2 && nb == 1) {

/*                 Build a 4-by-4 system Z**T * x = RHS */

#line 898 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 899 "dtgsy2.f"
		    z__[1] = a[is + isp1 * a_dim1];
#line 900 "dtgsy2.f"
		    z__[2] = -b[js + js * b_dim1];
#line 901 "dtgsy2.f"
		    z__[3] = 0.;

#line 903 "dtgsy2.f"
		    z__[8] = a[isp1 + is * a_dim1];
#line 904 "dtgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 905 "dtgsy2.f"
		    z__[10] = 0.;
#line 906 "dtgsy2.f"
		    z__[11] = -b[js + js * b_dim1];

#line 908 "dtgsy2.f"
		    z__[16] = d__[is + is * d_dim1];
#line 909 "dtgsy2.f"
		    z__[17] = d__[is + isp1 * d_dim1];
#line 910 "dtgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 911 "dtgsy2.f"
		    z__[19] = 0.;

#line 913 "dtgsy2.f"
		    z__[24] = 0.;
#line 914 "dtgsy2.f"
		    z__[25] = d__[isp1 + isp1 * d_dim1];
#line 915 "dtgsy2.f"
		    z__[26] = 0.;
#line 916 "dtgsy2.f"
		    z__[27] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 920 "dtgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 921 "dtgsy2.f"
		    rhs[1] = c__[isp1 + js * c_dim1];
#line 922 "dtgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 923 "dtgsy2.f"
		    rhs[3] = f[isp1 + js * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 927 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 928 "dtgsy2.f"
		    if (ierr > 0) {
#line 928 "dtgsy2.f"
			*info = ierr;
#line 928 "dtgsy2.f"
		    }

#line 931 "dtgsy2.f"
		    dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 932 "dtgsy2.f"
		    if (scaloc != 1.) {
#line 933 "dtgsy2.f"
			i__3 = *n;
#line 933 "dtgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 934 "dtgsy2.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 935 "dtgsy2.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 936 "dtgsy2.f"
/* L150: */
#line 936 "dtgsy2.f"
			}
#line 937 "dtgsy2.f"
			*scale *= scaloc;
#line 938 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 942 "dtgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 943 "dtgsy2.f"
		    c__[isp1 + js * c_dim1] = rhs[1];
#line 944 "dtgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 945 "dtgsy2.f"
		    f[isp1 + js * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 950 "dtgsy2.f"
		    if (j > p + 2) {
#line 951 "dtgsy2.f"
			i__3 = js - 1;
#line 951 "dtgsy2.f"
			dger_(&mb, &i__3, &c_b42, rhs, &c__1, &b[js * b_dim1 
				+ 1], &c__1, &f[is + f_dim1], ldf);
#line 953 "dtgsy2.f"
			i__3 = js - 1;
#line 953 "dtgsy2.f"
			dger_(&mb, &i__3, &c_b42, &rhs[2], &c__1, &e[js * 
				e_dim1 + 1], &c__1, &f[is + f_dim1], ldf);
#line 955 "dtgsy2.f"
		    }
#line 956 "dtgsy2.f"
		    if (i__ < p) {
#line 957 "dtgsy2.f"
			i__3 = *m - ie;
#line 957 "dtgsy2.f"
			dgemv_("T", &mb, &i__3, &c_b27, &a[is + (ie + 1) * 
				a_dim1], lda, rhs, &c__1, &c_b42, &c__[ie + 1 
				+ js * c_dim1], &c__1, (ftnlen)1);
#line 960 "dtgsy2.f"
			i__3 = *m - ie;
#line 960 "dtgsy2.f"
			dgemv_("T", &mb, &i__3, &c_b27, &d__[is + (ie + 1) * 
				d_dim1], ldd, &rhs[2], &c__1, &c_b42, &c__[ie 
				+ 1 + js * c_dim1], &c__1, (ftnlen)1);
#line 963 "dtgsy2.f"
		    }

#line 965 "dtgsy2.f"
		} else if (mb == 2 && nb == 2) {

/*                 Build an 8-by-8 system Z**T * x = RHS */

#line 969 "dtgsy2.f"
		    dlaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8, (
			    ftnlen)1);

#line 971 "dtgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 972 "dtgsy2.f"
		    z__[1] = a[is + isp1 * a_dim1];
#line 973 "dtgsy2.f"
		    z__[4] = -b[js + js * b_dim1];
#line 974 "dtgsy2.f"
		    z__[6] = -b[jsp1 + js * b_dim1];

#line 976 "dtgsy2.f"
		    z__[8] = a[isp1 + is * a_dim1];
#line 977 "dtgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 978 "dtgsy2.f"
		    z__[13] = -b[js + js * b_dim1];
#line 979 "dtgsy2.f"
		    z__[15] = -b[jsp1 + js * b_dim1];

#line 981 "dtgsy2.f"
		    z__[18] = a[is + is * a_dim1];
#line 982 "dtgsy2.f"
		    z__[19] = a[is + isp1 * a_dim1];
#line 983 "dtgsy2.f"
		    z__[20] = -b[js + jsp1 * b_dim1];
#line 984 "dtgsy2.f"
		    z__[22] = -b[jsp1 + jsp1 * b_dim1];

#line 986 "dtgsy2.f"
		    z__[26] = a[isp1 + is * a_dim1];
#line 987 "dtgsy2.f"
		    z__[27] = a[isp1 + isp1 * a_dim1];
#line 988 "dtgsy2.f"
		    z__[29] = -b[js + jsp1 * b_dim1];
#line 989 "dtgsy2.f"
		    z__[31] = -b[jsp1 + jsp1 * b_dim1];

#line 991 "dtgsy2.f"
		    z__[32] = d__[is + is * d_dim1];
#line 992 "dtgsy2.f"
		    z__[33] = d__[is + isp1 * d_dim1];
#line 993 "dtgsy2.f"
		    z__[36] = -e[js + js * e_dim1];

#line 995 "dtgsy2.f"
		    z__[41] = d__[isp1 + isp1 * d_dim1];
#line 996 "dtgsy2.f"
		    z__[45] = -e[js + js * e_dim1];

#line 998 "dtgsy2.f"
		    z__[50] = d__[is + is * d_dim1];
#line 999 "dtgsy2.f"
		    z__[51] = d__[is + isp1 * d_dim1];
#line 1000 "dtgsy2.f"
		    z__[52] = -e[js + jsp1 * e_dim1];
#line 1001 "dtgsy2.f"
		    z__[54] = -e[jsp1 + jsp1 * e_dim1];

#line 1003 "dtgsy2.f"
		    z__[59] = d__[isp1 + isp1 * d_dim1];
#line 1004 "dtgsy2.f"
		    z__[61] = -e[js + jsp1 * e_dim1];
#line 1005 "dtgsy2.f"
		    z__[63] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 1009 "dtgsy2.f"
		    k = 1;
#line 1010 "dtgsy2.f"
		    ii = mb * nb + 1;
#line 1011 "dtgsy2.f"
		    i__3 = nb - 1;
#line 1011 "dtgsy2.f"
		    for (jj = 0; jj <= i__3; ++jj) {
#line 1012 "dtgsy2.f"
			dcopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &
				rhs[k - 1], &c__1);
#line 1013 "dtgsy2.f"
			dcopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[
				ii - 1], &c__1);
#line 1014 "dtgsy2.f"
			k += mb;
#line 1015 "dtgsy2.f"
			ii += mb;
#line 1016 "dtgsy2.f"
/* L160: */
#line 1016 "dtgsy2.f"
		    }


/*                 Solve Z**T * x = RHS */

#line 1021 "dtgsy2.f"
		    dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 1022 "dtgsy2.f"
		    if (ierr > 0) {
#line 1022 "dtgsy2.f"
			*info = ierr;
#line 1022 "dtgsy2.f"
		    }

#line 1025 "dtgsy2.f"
		    dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 1026 "dtgsy2.f"
		    if (scaloc != 1.) {
#line 1027 "dtgsy2.f"
			i__3 = *n;
#line 1027 "dtgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 1028 "dtgsy2.f"
			    dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 1029 "dtgsy2.f"
			    dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 1030 "dtgsy2.f"
/* L170: */
#line 1030 "dtgsy2.f"
			}
#line 1031 "dtgsy2.f"
			*scale *= scaloc;
#line 1032 "dtgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 1036 "dtgsy2.f"
		    k = 1;
#line 1037 "dtgsy2.f"
		    ii = mb * nb + 1;
#line 1038 "dtgsy2.f"
		    i__3 = nb - 1;
#line 1038 "dtgsy2.f"
		    for (jj = 0; jj <= i__3; ++jj) {
#line 1039 "dtgsy2.f"
			dcopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
				c_dim1], &c__1);
#line 1040 "dtgsy2.f"
			dcopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
				f_dim1], &c__1);
#line 1041 "dtgsy2.f"
			k += mb;
#line 1042 "dtgsy2.f"
			ii += mb;
#line 1043 "dtgsy2.f"
/* L180: */
#line 1043 "dtgsy2.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 1048 "dtgsy2.f"
		    if (j > p + 2) {
#line 1049 "dtgsy2.f"
			i__3 = js - 1;
#line 1049 "dtgsy2.f"
			dgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &c__[is + 
				js * c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &
				c_b42, &f[is + f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 1052 "dtgsy2.f"
			i__3 = js - 1;
#line 1052 "dtgsy2.f"
			dgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &f[is + js *
				 f_dim1], ldf, &e[js * e_dim1 + 1], lde, &
				c_b42, &f[is + f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 1055 "dtgsy2.f"
		    }
#line 1056 "dtgsy2.f"
		    if (i__ < p) {
#line 1057 "dtgsy2.f"
			i__3 = *m - ie;
#line 1057 "dtgsy2.f"
			dgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &a[is + (ie 
				+ 1) * a_dim1], lda, &c__[is + js * c_dim1], 
				ldc, &c_b42, &c__[ie + 1 + js * c_dim1], ldc, 
				(ftnlen)1, (ftnlen)1);
#line 1060 "dtgsy2.f"
			i__3 = *m - ie;
#line 1060 "dtgsy2.f"
			dgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &d__[is + (
				ie + 1) * d_dim1], ldd, &f[is + js * f_dim1], 
				ldf, &c_b42, &c__[ie + 1 + js * c_dim1], ldc, 
				(ftnlen)1, (ftnlen)1);
#line 1063 "dtgsy2.f"
		    }

#line 1065 "dtgsy2.f"
		}

#line 1067 "dtgsy2.f"
/* L190: */
#line 1067 "dtgsy2.f"
	    }
#line 1068 "dtgsy2.f"
/* L200: */
#line 1068 "dtgsy2.f"
	}

#line 1070 "dtgsy2.f"
    }
#line 1071 "dtgsy2.f"
    return 0;

/*     End of DTGSY2 */

} /* dtgsy2_ */

