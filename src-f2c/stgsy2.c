#line 1 "stgsy2.f"
/* stgsy2.f -- translated by f2c (version 20100827).
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

#line 1 "stgsy2.f"
/* Table of constant values */

static integer c__8 = 8;
static integer c__1 = 1;
static doublereal c_b27 = -1.;
static doublereal c_b42 = 1.;
static doublereal c_b56 = 0.;

/* > \brief \b STGSY2 solves the generalized Sylvester equation (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGSY2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsy2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsy2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsy2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, */
/*                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL, */
/*                          IWORK, PQ, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N, */
/*      $                   PQ */
/*       REAL               RDSCAL, RDSUM, SCALE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ), */
/*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSY2 solves the generalized Sylvester equation: */
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
/* > sigma_min(Z) using reverse communicaton with SLACON. */
/* > */
/* > STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL */
/* > of an upper bound on the separation between to matrix pairs. Then */
/* > the input (A, D), (B, E) are sub-pencils of the matrix pair in */
/* > STGSYL. See STGSYL for details. */
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
/* >               pairs is computed. (SGECON on sub-systems is used.) */
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
/* >          A is REAL array, dimension (LDA, M) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          C is REAL array, dimension (LDC, N) */
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
/* >          D is REAL array, dimension (LDD, M) */
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
/* >          E is REAL array, dimension (LDE, N) */
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
/* >          F is REAL array, dimension (LDF, N) */
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
/* >          SCALE is REAL */
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
/* >          RDSUM is REAL */
/* >          On entry, the sum of squares of computed contributions to */
/* >          the Dif-estimate under computation by STGSYL, where the */
/* >          scaling factor RDSCAL (see below) has been factored out. */
/* >          On exit, the corresponding sum of squares updated with the */
/* >          contributions from the current sub-system. */
/* >          If TRANS = 'T' RDSUM is not touched. */
/* >          NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RDSCAL */
/* > \verbatim */
/* >          RDSCAL is REAL */
/* >          On entry, scaling factor used to prevent overflow in RDSUM. */
/* >          On exit, RDSCAL is updated w.r.t. the current contributions */
/* >          in RDSUM. */
/* >          If TRANS = 'T', RDSCAL is not touched. */
/* >          NOTE: RDSCAL only makes sense when STGSY2 is called by */
/* >                STGSYL. */
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

/* > \date November 2015 */

/* > \ingroup realSYauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/*  ===================================================================== */
/* Subroutine */ int stgsy2_(char *trans, integer *ijob, integer *m, integer *
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
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr, zdim, ipiv[8], jpiv[8];
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), sgemv_(
	    char *, integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), scopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *), saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), sgesc2_(integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *), sgetc2_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int slatdf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static logical notran;


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*  Replaced various illegal calls to SCOPY by calls to SLASET. */
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

#line 329 "stgsy2.f"
    /* Parameter adjustments */
#line 329 "stgsy2.f"
    a_dim1 = *lda;
#line 329 "stgsy2.f"
    a_offset = 1 + a_dim1;
#line 329 "stgsy2.f"
    a -= a_offset;
#line 329 "stgsy2.f"
    b_dim1 = *ldb;
#line 329 "stgsy2.f"
    b_offset = 1 + b_dim1;
#line 329 "stgsy2.f"
    b -= b_offset;
#line 329 "stgsy2.f"
    c_dim1 = *ldc;
#line 329 "stgsy2.f"
    c_offset = 1 + c_dim1;
#line 329 "stgsy2.f"
    c__ -= c_offset;
#line 329 "stgsy2.f"
    d_dim1 = *ldd;
#line 329 "stgsy2.f"
    d_offset = 1 + d_dim1;
#line 329 "stgsy2.f"
    d__ -= d_offset;
#line 329 "stgsy2.f"
    e_dim1 = *lde;
#line 329 "stgsy2.f"
    e_offset = 1 + e_dim1;
#line 329 "stgsy2.f"
    e -= e_offset;
#line 329 "stgsy2.f"
    f_dim1 = *ldf;
#line 329 "stgsy2.f"
    f_offset = 1 + f_dim1;
#line 329 "stgsy2.f"
    f -= f_offset;
#line 329 "stgsy2.f"
    --iwork;
#line 329 "stgsy2.f"

#line 329 "stgsy2.f"
    /* Function Body */
#line 329 "stgsy2.f"
    *info = 0;
#line 330 "stgsy2.f"
    ierr = 0;
#line 331 "stgsy2.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 332 "stgsy2.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 333 "stgsy2.f"
	*info = -1;
#line 334 "stgsy2.f"
    } else if (notran) {
#line 335 "stgsy2.f"
	if (*ijob < 0 || *ijob > 2) {
#line 336 "stgsy2.f"
	    *info = -2;
#line 337 "stgsy2.f"
	}
#line 338 "stgsy2.f"
    }
#line 339 "stgsy2.f"
    if (*info == 0) {
#line 340 "stgsy2.f"
	if (*m <= 0) {
#line 341 "stgsy2.f"
	    *info = -3;
#line 342 "stgsy2.f"
	} else if (*n <= 0) {
#line 343 "stgsy2.f"
	    *info = -4;
#line 344 "stgsy2.f"
	} else if (*lda < max(1,*m)) {
#line 345 "stgsy2.f"
	    *info = -6;
#line 346 "stgsy2.f"
	} else if (*ldb < max(1,*n)) {
#line 347 "stgsy2.f"
	    *info = -8;
#line 348 "stgsy2.f"
	} else if (*ldc < max(1,*m)) {
#line 349 "stgsy2.f"
	    *info = -10;
#line 350 "stgsy2.f"
	} else if (*ldd < max(1,*m)) {
#line 351 "stgsy2.f"
	    *info = -12;
#line 352 "stgsy2.f"
	} else if (*lde < max(1,*n)) {
#line 353 "stgsy2.f"
	    *info = -14;
#line 354 "stgsy2.f"
	} else if (*ldf < max(1,*m)) {
#line 355 "stgsy2.f"
	    *info = -16;
#line 356 "stgsy2.f"
	}
#line 357 "stgsy2.f"
    }
#line 358 "stgsy2.f"
    if (*info != 0) {
#line 359 "stgsy2.f"
	i__1 = -(*info);
#line 359 "stgsy2.f"
	xerbla_("STGSY2", &i__1, (ftnlen)6);
#line 360 "stgsy2.f"
	return 0;
#line 361 "stgsy2.f"
    }

/*     Determine block structure of A */

#line 365 "stgsy2.f"
    *pq = 0;
#line 366 "stgsy2.f"
    p = 0;
#line 367 "stgsy2.f"
    i__ = 1;
#line 368 "stgsy2.f"
L10:
#line 369 "stgsy2.f"
    if (i__ > *m) {
#line 369 "stgsy2.f"
	goto L20;
#line 369 "stgsy2.f"
    }
#line 371 "stgsy2.f"
    ++p;
#line 372 "stgsy2.f"
    iwork[p] = i__;
#line 373 "stgsy2.f"
    if (i__ == *m) {
#line 373 "stgsy2.f"
	goto L20;
#line 373 "stgsy2.f"
    }
#line 375 "stgsy2.f"
    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 376 "stgsy2.f"
	i__ += 2;
#line 377 "stgsy2.f"
    } else {
#line 378 "stgsy2.f"
	++i__;
#line 379 "stgsy2.f"
    }
#line 380 "stgsy2.f"
    goto L10;
#line 381 "stgsy2.f"
L20:
#line 382 "stgsy2.f"
    iwork[p + 1] = *m + 1;

/*     Determine block structure of B */

#line 386 "stgsy2.f"
    q = p + 1;
#line 387 "stgsy2.f"
    j = 1;
#line 388 "stgsy2.f"
L30:
#line 389 "stgsy2.f"
    if (j > *n) {
#line 389 "stgsy2.f"
	goto L40;
#line 389 "stgsy2.f"
    }
#line 391 "stgsy2.f"
    ++q;
#line 392 "stgsy2.f"
    iwork[q] = j;
#line 393 "stgsy2.f"
    if (j == *n) {
#line 393 "stgsy2.f"
	goto L40;
#line 393 "stgsy2.f"
    }
#line 395 "stgsy2.f"
    if (b[j + 1 + j * b_dim1] != 0.) {
#line 396 "stgsy2.f"
	j += 2;
#line 397 "stgsy2.f"
    } else {
#line 398 "stgsy2.f"
	++j;
#line 399 "stgsy2.f"
    }
#line 400 "stgsy2.f"
    goto L30;
#line 401 "stgsy2.f"
L40:
#line 402 "stgsy2.f"
    iwork[q + 1] = *n + 1;
#line 403 "stgsy2.f"
    *pq = p * (q - p - 1);

#line 405 "stgsy2.f"
    if (notran) {

/*        Solve (I, J) - subsystem */
/*           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J) */
/*           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J) */
/*        for I = P, P - 1, ..., 1; J = 1, 2, ..., Q */

#line 412 "stgsy2.f"
	*scale = 1.;
#line 413 "stgsy2.f"
	scaloc = 1.;
#line 414 "stgsy2.f"
	i__1 = q;
#line 414 "stgsy2.f"
	for (j = p + 2; j <= i__1; ++j) {
#line 415 "stgsy2.f"
	    js = iwork[j];
#line 416 "stgsy2.f"
	    jsp1 = js + 1;
#line 417 "stgsy2.f"
	    je = iwork[j + 1] - 1;
#line 418 "stgsy2.f"
	    nb = je - js + 1;
#line 419 "stgsy2.f"
	    for (i__ = p; i__ >= 1; --i__) {

#line 421 "stgsy2.f"
		is = iwork[i__];
#line 422 "stgsy2.f"
		isp1 = is + 1;
#line 423 "stgsy2.f"
		ie = iwork[i__ + 1] - 1;
#line 424 "stgsy2.f"
		mb = ie - is + 1;
#line 425 "stgsy2.f"
		zdim = mb * nb << 1;

#line 427 "stgsy2.f"
		if (mb == 1 && nb == 1) {

/*                 Build a 2-by-2 system Z * x = RHS */

#line 431 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 432 "stgsy2.f"
		    z__[1] = d__[is + is * d_dim1];
#line 433 "stgsy2.f"
		    z__[8] = -b[js + js * b_dim1];
#line 434 "stgsy2.f"
		    z__[9] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 438 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 439 "stgsy2.f"
		    rhs[1] = f[is + js * f_dim1];

/*                 Solve Z * x = RHS */

#line 443 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 444 "stgsy2.f"
		    if (ierr > 0) {
#line 444 "stgsy2.f"
			*info = ierr;
#line 444 "stgsy2.f"
		    }

#line 447 "stgsy2.f"
		    if (*ijob == 0) {
#line 448 "stgsy2.f"
			sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 450 "stgsy2.f"
			if (scaloc != 1.) {
#line 451 "stgsy2.f"
			    i__2 = *n;
#line 451 "stgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 452 "stgsy2.f"
				sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 453 "stgsy2.f"
				sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 454 "stgsy2.f"
/* L50: */
#line 454 "stgsy2.f"
			    }
#line 455 "stgsy2.f"
			    *scale *= scaloc;
#line 456 "stgsy2.f"
			}
#line 457 "stgsy2.f"
		    } else {
#line 458 "stgsy2.f"
			slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 460 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 464 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 465 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[1];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 470 "stgsy2.f"
		    if (i__ > 1) {
#line 471 "stgsy2.f"
			alpha = -rhs[0];
#line 472 "stgsy2.f"
			i__2 = is - 1;
#line 472 "stgsy2.f"
			saxpy_(&i__2, &alpha, &a[is * a_dim1 + 1], &c__1, &
				c__[js * c_dim1 + 1], &c__1);
#line 474 "stgsy2.f"
			i__2 = is - 1;
#line 474 "stgsy2.f"
			saxpy_(&i__2, &alpha, &d__[is * d_dim1 + 1], &c__1, &
				f[js * f_dim1 + 1], &c__1);
#line 476 "stgsy2.f"
		    }
#line 477 "stgsy2.f"
		    if (j < q) {
#line 478 "stgsy2.f"
			i__2 = *n - je;
#line 478 "stgsy2.f"
			saxpy_(&i__2, &rhs[1], &b[js + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 480 "stgsy2.f"
			i__2 = *n - je;
#line 480 "stgsy2.f"
			saxpy_(&i__2, &rhs[1], &e[js + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 482 "stgsy2.f"
		    }

#line 484 "stgsy2.f"
		} else if (mb == 1 && nb == 2) {

/*                 Build a 4-by-4 system Z * x = RHS */

#line 488 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 489 "stgsy2.f"
		    z__[1] = 0.;
#line 490 "stgsy2.f"
		    z__[2] = d__[is + is * d_dim1];
#line 491 "stgsy2.f"
		    z__[3] = 0.;

#line 493 "stgsy2.f"
		    z__[8] = 0.;
#line 494 "stgsy2.f"
		    z__[9] = a[is + is * a_dim1];
#line 495 "stgsy2.f"
		    z__[10] = 0.;
#line 496 "stgsy2.f"
		    z__[11] = d__[is + is * d_dim1];

#line 498 "stgsy2.f"
		    z__[16] = -b[js + js * b_dim1];
#line 499 "stgsy2.f"
		    z__[17] = -b[js + jsp1 * b_dim1];
#line 500 "stgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 501 "stgsy2.f"
		    z__[19] = -e[js + jsp1 * e_dim1];

#line 503 "stgsy2.f"
		    z__[24] = -b[jsp1 + js * b_dim1];
#line 504 "stgsy2.f"
		    z__[25] = -b[jsp1 + jsp1 * b_dim1];
#line 505 "stgsy2.f"
		    z__[26] = 0.;
#line 506 "stgsy2.f"
		    z__[27] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 510 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 511 "stgsy2.f"
		    rhs[1] = c__[is + jsp1 * c_dim1];
#line 512 "stgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 513 "stgsy2.f"
		    rhs[3] = f[is + jsp1 * f_dim1];

/*                 Solve Z * x = RHS */

#line 517 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 518 "stgsy2.f"
		    if (ierr > 0) {
#line 518 "stgsy2.f"
			*info = ierr;
#line 518 "stgsy2.f"
		    }

#line 521 "stgsy2.f"
		    if (*ijob == 0) {
#line 522 "stgsy2.f"
			sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 524 "stgsy2.f"
			if (scaloc != 1.) {
#line 525 "stgsy2.f"
			    i__2 = *n;
#line 525 "stgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 526 "stgsy2.f"
				sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 527 "stgsy2.f"
				sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 528 "stgsy2.f"
/* L60: */
#line 528 "stgsy2.f"
			    }
#line 529 "stgsy2.f"
			    *scale *= scaloc;
#line 530 "stgsy2.f"
			}
#line 531 "stgsy2.f"
		    } else {
#line 532 "stgsy2.f"
			slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 534 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 538 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 539 "stgsy2.f"
		    c__[is + jsp1 * c_dim1] = rhs[1];
#line 540 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 541 "stgsy2.f"
		    f[is + jsp1 * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 546 "stgsy2.f"
		    if (i__ > 1) {
#line 547 "stgsy2.f"
			i__2 = is - 1;
#line 547 "stgsy2.f"
			sger_(&i__2, &nb, &c_b27, &a[is * a_dim1 + 1], &c__1, 
				rhs, &c__1, &c__[js * c_dim1 + 1], ldc);
#line 549 "stgsy2.f"
			i__2 = is - 1;
#line 549 "stgsy2.f"
			sger_(&i__2, &nb, &c_b27, &d__[is * d_dim1 + 1], &
				c__1, rhs, &c__1, &f[js * f_dim1 + 1], ldf);
#line 551 "stgsy2.f"
		    }
#line 552 "stgsy2.f"
		    if (j < q) {
#line 553 "stgsy2.f"
			i__2 = *n - je;
#line 553 "stgsy2.f"
			saxpy_(&i__2, &rhs[2], &b[js + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 555 "stgsy2.f"
			i__2 = *n - je;
#line 555 "stgsy2.f"
			saxpy_(&i__2, &rhs[2], &e[js + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 557 "stgsy2.f"
			i__2 = *n - je;
#line 557 "stgsy2.f"
			saxpy_(&i__2, &rhs[3], &b[jsp1 + (je + 1) * b_dim1], 
				ldb, &c__[is + (je + 1) * c_dim1], ldc);
#line 559 "stgsy2.f"
			i__2 = *n - je;
#line 559 "stgsy2.f"
			saxpy_(&i__2, &rhs[3], &e[jsp1 + (je + 1) * e_dim1], 
				lde, &f[is + (je + 1) * f_dim1], ldf);
#line 561 "stgsy2.f"
		    }

#line 563 "stgsy2.f"
		} else if (mb == 2 && nb == 1) {

/*                 Build a 4-by-4 system Z * x = RHS */

#line 567 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 568 "stgsy2.f"
		    z__[1] = a[isp1 + is * a_dim1];
#line 569 "stgsy2.f"
		    z__[2] = d__[is + is * d_dim1];
#line 570 "stgsy2.f"
		    z__[3] = 0.;

#line 572 "stgsy2.f"
		    z__[8] = a[is + isp1 * a_dim1];
#line 573 "stgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 574 "stgsy2.f"
		    z__[10] = d__[is + isp1 * d_dim1];
#line 575 "stgsy2.f"
		    z__[11] = d__[isp1 + isp1 * d_dim1];

#line 577 "stgsy2.f"
		    z__[16] = -b[js + js * b_dim1];
#line 578 "stgsy2.f"
		    z__[17] = 0.;
#line 579 "stgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 580 "stgsy2.f"
		    z__[19] = 0.;

#line 582 "stgsy2.f"
		    z__[24] = 0.;
#line 583 "stgsy2.f"
		    z__[25] = -b[js + js * b_dim1];
#line 584 "stgsy2.f"
		    z__[26] = 0.;
#line 585 "stgsy2.f"
		    z__[27] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 589 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 590 "stgsy2.f"
		    rhs[1] = c__[isp1 + js * c_dim1];
#line 591 "stgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 592 "stgsy2.f"
		    rhs[3] = f[isp1 + js * f_dim1];

/*                 Solve Z * x = RHS */

#line 596 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 597 "stgsy2.f"
		    if (ierr > 0) {
#line 597 "stgsy2.f"
			*info = ierr;
#line 597 "stgsy2.f"
		    }
#line 599 "stgsy2.f"
		    if (*ijob == 0) {
#line 600 "stgsy2.f"
			sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 602 "stgsy2.f"
			if (scaloc != 1.) {
#line 603 "stgsy2.f"
			    i__2 = *n;
#line 603 "stgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 604 "stgsy2.f"
				sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 605 "stgsy2.f"
				sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 606 "stgsy2.f"
/* L70: */
#line 606 "stgsy2.f"
			    }
#line 607 "stgsy2.f"
			    *scale *= scaloc;
#line 608 "stgsy2.f"
			}
#line 609 "stgsy2.f"
		    } else {
#line 610 "stgsy2.f"
			slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 612 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 616 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 617 "stgsy2.f"
		    c__[isp1 + js * c_dim1] = rhs[1];
#line 618 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 619 "stgsy2.f"
		    f[isp1 + js * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 624 "stgsy2.f"
		    if (i__ > 1) {
#line 625 "stgsy2.f"
			i__2 = is - 1;
#line 625 "stgsy2.f"
			sgemv_("N", &i__2, &mb, &c_b27, &a[is * a_dim1 + 1], 
				lda, rhs, &c__1, &c_b42, &c__[js * c_dim1 + 1]
				, &c__1, (ftnlen)1);
#line 627 "stgsy2.f"
			i__2 = is - 1;
#line 627 "stgsy2.f"
			sgemv_("N", &i__2, &mb, &c_b27, &d__[is * d_dim1 + 1],
				 ldd, rhs, &c__1, &c_b42, &f[js * f_dim1 + 1],
				 &c__1, (ftnlen)1);
#line 629 "stgsy2.f"
		    }
#line 630 "stgsy2.f"
		    if (j < q) {
#line 631 "stgsy2.f"
			i__2 = *n - je;
#line 631 "stgsy2.f"
			sger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &b[js + (je 
				+ 1) * b_dim1], ldb, &c__[is + (je + 1) * 
				c_dim1], ldc);
#line 633 "stgsy2.f"
			i__2 = *n - je;
#line 633 "stgsy2.f"
			sger_(&mb, &i__2, &c_b42, &rhs[2], &c__1, &e[js + (je 
				+ 1) * e_dim1], lde, &f[is + (je + 1) * 
				f_dim1], ldf);
#line 635 "stgsy2.f"
		    }

#line 637 "stgsy2.f"
		} else if (mb == 2 && nb == 2) {

/*                 Build an 8-by-8 system Z * x = RHS */

#line 641 "stgsy2.f"
		    slaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8, (
			    ftnlen)1);

#line 643 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 644 "stgsy2.f"
		    z__[1] = a[isp1 + is * a_dim1];
#line 645 "stgsy2.f"
		    z__[4] = d__[is + is * d_dim1];

#line 647 "stgsy2.f"
		    z__[8] = a[is + isp1 * a_dim1];
#line 648 "stgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 649 "stgsy2.f"
		    z__[12] = d__[is + isp1 * d_dim1];
#line 650 "stgsy2.f"
		    z__[13] = d__[isp1 + isp1 * d_dim1];

#line 652 "stgsy2.f"
		    z__[18] = a[is + is * a_dim1];
#line 653 "stgsy2.f"
		    z__[19] = a[isp1 + is * a_dim1];
#line 654 "stgsy2.f"
		    z__[22] = d__[is + is * d_dim1];

#line 656 "stgsy2.f"
		    z__[26] = a[is + isp1 * a_dim1];
#line 657 "stgsy2.f"
		    z__[27] = a[isp1 + isp1 * a_dim1];
#line 658 "stgsy2.f"
		    z__[30] = d__[is + isp1 * d_dim1];
#line 659 "stgsy2.f"
		    z__[31] = d__[isp1 + isp1 * d_dim1];

#line 661 "stgsy2.f"
		    z__[32] = -b[js + js * b_dim1];
#line 662 "stgsy2.f"
		    z__[34] = -b[js + jsp1 * b_dim1];
#line 663 "stgsy2.f"
		    z__[36] = -e[js + js * e_dim1];
#line 664 "stgsy2.f"
		    z__[38] = -e[js + jsp1 * e_dim1];

#line 666 "stgsy2.f"
		    z__[41] = -b[js + js * b_dim1];
#line 667 "stgsy2.f"
		    z__[43] = -b[js + jsp1 * b_dim1];
#line 668 "stgsy2.f"
		    z__[45] = -e[js + js * e_dim1];
#line 669 "stgsy2.f"
		    z__[47] = -e[js + jsp1 * e_dim1];

#line 671 "stgsy2.f"
		    z__[48] = -b[jsp1 + js * b_dim1];
#line 672 "stgsy2.f"
		    z__[50] = -b[jsp1 + jsp1 * b_dim1];
#line 673 "stgsy2.f"
		    z__[54] = -e[jsp1 + jsp1 * e_dim1];

#line 675 "stgsy2.f"
		    z__[57] = -b[jsp1 + js * b_dim1];
#line 676 "stgsy2.f"
		    z__[59] = -b[jsp1 + jsp1 * b_dim1];
#line 677 "stgsy2.f"
		    z__[63] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 681 "stgsy2.f"
		    k = 1;
#line 682 "stgsy2.f"
		    ii = mb * nb + 1;
#line 683 "stgsy2.f"
		    i__2 = nb - 1;
#line 683 "stgsy2.f"
		    for (jj = 0; jj <= i__2; ++jj) {
#line 684 "stgsy2.f"
			scopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &
				rhs[k - 1], &c__1);
#line 685 "stgsy2.f"
			scopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[
				ii - 1], &c__1);
#line 686 "stgsy2.f"
			k += mb;
#line 687 "stgsy2.f"
			ii += mb;
#line 688 "stgsy2.f"
/* L80: */
#line 688 "stgsy2.f"
		    }

/*                 Solve Z * x = RHS */

#line 692 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 693 "stgsy2.f"
		    if (ierr > 0) {
#line 693 "stgsy2.f"
			*info = ierr;
#line 693 "stgsy2.f"
		    }
#line 695 "stgsy2.f"
		    if (*ijob == 0) {
#line 696 "stgsy2.f"
			sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 698 "stgsy2.f"
			if (scaloc != 1.) {
#line 699 "stgsy2.f"
			    i__2 = *n;
#line 699 "stgsy2.f"
			    for (k = 1; k <= i__2; ++k) {
#line 700 "stgsy2.f"
				sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &
					c__1);
#line 701 "stgsy2.f"
				sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 702 "stgsy2.f"
/* L90: */
#line 702 "stgsy2.f"
			    }
#line 703 "stgsy2.f"
			    *scale *= scaloc;
#line 704 "stgsy2.f"
			}
#line 705 "stgsy2.f"
		    } else {
#line 706 "stgsy2.f"
			slatdf_(ijob, &zdim, z__, &c__8, rhs, rdsum, rdscal, 
				ipiv, jpiv);
#line 708 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 712 "stgsy2.f"
		    k = 1;
#line 713 "stgsy2.f"
		    ii = mb * nb + 1;
#line 714 "stgsy2.f"
		    i__2 = nb - 1;
#line 714 "stgsy2.f"
		    for (jj = 0; jj <= i__2; ++jj) {
#line 715 "stgsy2.f"
			scopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
				c_dim1], &c__1);
#line 716 "stgsy2.f"
			scopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
				f_dim1], &c__1);
#line 717 "stgsy2.f"
			k += mb;
#line 718 "stgsy2.f"
			ii += mb;
#line 719 "stgsy2.f"
/* L100: */
#line 719 "stgsy2.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 724 "stgsy2.f"
		    if (i__ > 1) {
#line 725 "stgsy2.f"
			i__2 = is - 1;
#line 725 "stgsy2.f"
			sgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &a[is * 
				a_dim1 + 1], lda, rhs, &mb, &c_b42, &c__[js * 
				c_dim1 + 1], ldc, (ftnlen)1, (ftnlen)1);
#line 728 "stgsy2.f"
			i__2 = is - 1;
#line 728 "stgsy2.f"
			sgemm_("N", "N", &i__2, &nb, &mb, &c_b27, &d__[is * 
				d_dim1 + 1], ldd, rhs, &mb, &c_b42, &f[js * 
				f_dim1 + 1], ldf, (ftnlen)1, (ftnlen)1);
#line 731 "stgsy2.f"
		    }
#line 732 "stgsy2.f"
		    if (j < q) {
#line 733 "stgsy2.f"
			k = mb * nb + 1;
#line 734 "stgsy2.f"
			i__2 = *n - je;
#line 734 "stgsy2.f"
			sgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1],
				 &mb, &b[js + (je + 1) * b_dim1], ldb, &c_b42,
				 &c__[is + (je + 1) * c_dim1], ldc, (ftnlen)1,
				 (ftnlen)1);
#line 737 "stgsy2.f"
			i__2 = *n - je;
#line 737 "stgsy2.f"
			sgemm_("N", "N", &mb, &i__2, &nb, &c_b42, &rhs[k - 1],
				 &mb, &e[js + (je + 1) * e_dim1], lde, &c_b42,
				 &f[is + (je + 1) * f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 740 "stgsy2.f"
		    }

#line 742 "stgsy2.f"
		}

#line 744 "stgsy2.f"
/* L110: */
#line 744 "stgsy2.f"
	    }
#line 745 "stgsy2.f"
/* L120: */
#line 745 "stgsy2.f"
	}
#line 746 "stgsy2.f"
    } else {

/*        Solve (I, J) - subsystem */
/*             A(I, I)**T * R(I, J) + D(I, I)**T * L(J, J)  =  C(I, J) */
/*             R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J) */
/*        for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1 */

#line 753 "stgsy2.f"
	*scale = 1.;
#line 754 "stgsy2.f"
	scaloc = 1.;
#line 755 "stgsy2.f"
	i__1 = p;
#line 755 "stgsy2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 757 "stgsy2.f"
	    is = iwork[i__];
#line 758 "stgsy2.f"
	    isp1 = is + 1;
#line 759 "stgsy2.f"
	    ie = iwork[i__ + 1] - 1;
#line 760 "stgsy2.f"
	    mb = ie - is + 1;
#line 761 "stgsy2.f"
	    i__2 = p + 2;
#line 761 "stgsy2.f"
	    for (j = q; j >= i__2; --j) {

#line 763 "stgsy2.f"
		js = iwork[j];
#line 764 "stgsy2.f"
		jsp1 = js + 1;
#line 765 "stgsy2.f"
		je = iwork[j + 1] - 1;
#line 766 "stgsy2.f"
		nb = je - js + 1;
#line 767 "stgsy2.f"
		zdim = mb * nb << 1;
#line 768 "stgsy2.f"
		if (mb == 1 && nb == 1) {

/*                 Build a 2-by-2 system Z**T * x = RHS */

#line 772 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 773 "stgsy2.f"
		    z__[1] = -b[js + js * b_dim1];
#line 774 "stgsy2.f"
		    z__[8] = d__[is + is * d_dim1];
#line 775 "stgsy2.f"
		    z__[9] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 779 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 780 "stgsy2.f"
		    rhs[1] = f[is + js * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 784 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 785 "stgsy2.f"
		    if (ierr > 0) {
#line 785 "stgsy2.f"
			*info = ierr;
#line 785 "stgsy2.f"
		    }

#line 788 "stgsy2.f"
		    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 789 "stgsy2.f"
		    if (scaloc != 1.) {
#line 790 "stgsy2.f"
			i__3 = *n;
#line 790 "stgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 791 "stgsy2.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 792 "stgsy2.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 793 "stgsy2.f"
/* L130: */
#line 793 "stgsy2.f"
			}
#line 794 "stgsy2.f"
			*scale *= scaloc;
#line 795 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 799 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 800 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[1];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 805 "stgsy2.f"
		    if (j > p + 2) {
#line 806 "stgsy2.f"
			alpha = rhs[0];
#line 807 "stgsy2.f"
			i__3 = js - 1;
#line 807 "stgsy2.f"
			saxpy_(&i__3, &alpha, &b[js * b_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 809 "stgsy2.f"
			alpha = rhs[1];
#line 810 "stgsy2.f"
			i__3 = js - 1;
#line 810 "stgsy2.f"
			saxpy_(&i__3, &alpha, &e[js * e_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 812 "stgsy2.f"
		    }
#line 813 "stgsy2.f"
		    if (i__ < p) {
#line 814 "stgsy2.f"
			alpha = -rhs[0];
#line 815 "stgsy2.f"
			i__3 = *m - ie;
#line 815 "stgsy2.f"
			saxpy_(&i__3, &alpha, &a[is + (ie + 1) * a_dim1], lda,
				 &c__[ie + 1 + js * c_dim1], &c__1);
#line 817 "stgsy2.f"
			alpha = -rhs[1];
#line 818 "stgsy2.f"
			i__3 = *m - ie;
#line 818 "stgsy2.f"
			saxpy_(&i__3, &alpha, &d__[is + (ie + 1) * d_dim1], 
				ldd, &c__[ie + 1 + js * c_dim1], &c__1);
#line 820 "stgsy2.f"
		    }

#line 822 "stgsy2.f"
		} else if (mb == 1 && nb == 2) {

/*                 Build a 4-by-4 system Z**T * x = RHS */

#line 826 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 827 "stgsy2.f"
		    z__[1] = 0.;
#line 828 "stgsy2.f"
		    z__[2] = -b[js + js * b_dim1];
#line 829 "stgsy2.f"
		    z__[3] = -b[jsp1 + js * b_dim1];

#line 831 "stgsy2.f"
		    z__[8] = 0.;
#line 832 "stgsy2.f"
		    z__[9] = a[is + is * a_dim1];
#line 833 "stgsy2.f"
		    z__[10] = -b[js + jsp1 * b_dim1];
#line 834 "stgsy2.f"
		    z__[11] = -b[jsp1 + jsp1 * b_dim1];

#line 836 "stgsy2.f"
		    z__[16] = d__[is + is * d_dim1];
#line 837 "stgsy2.f"
		    z__[17] = 0.;
#line 838 "stgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 839 "stgsy2.f"
		    z__[19] = 0.;

#line 841 "stgsy2.f"
		    z__[24] = 0.;
#line 842 "stgsy2.f"
		    z__[25] = d__[is + is * d_dim1];
#line 843 "stgsy2.f"
		    z__[26] = -e[js + jsp1 * e_dim1];
#line 844 "stgsy2.f"
		    z__[27] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 848 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 849 "stgsy2.f"
		    rhs[1] = c__[is + jsp1 * c_dim1];
#line 850 "stgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 851 "stgsy2.f"
		    rhs[3] = f[is + jsp1 * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 855 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 856 "stgsy2.f"
		    if (ierr > 0) {
#line 856 "stgsy2.f"
			*info = ierr;
#line 856 "stgsy2.f"
		    }
#line 858 "stgsy2.f"
		    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 859 "stgsy2.f"
		    if (scaloc != 1.) {
#line 860 "stgsy2.f"
			i__3 = *n;
#line 860 "stgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 861 "stgsy2.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 862 "stgsy2.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 863 "stgsy2.f"
/* L140: */
#line 863 "stgsy2.f"
			}
#line 864 "stgsy2.f"
			*scale *= scaloc;
#line 865 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 869 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 870 "stgsy2.f"
		    c__[is + jsp1 * c_dim1] = rhs[1];
#line 871 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 872 "stgsy2.f"
		    f[is + jsp1 * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 877 "stgsy2.f"
		    if (j > p + 2) {
#line 878 "stgsy2.f"
			i__3 = js - 1;
#line 878 "stgsy2.f"
			saxpy_(&i__3, rhs, &b[js * b_dim1 + 1], &c__1, &f[is 
				+ f_dim1], ldf);
#line 880 "stgsy2.f"
			i__3 = js - 1;
#line 880 "stgsy2.f"
			saxpy_(&i__3, &rhs[1], &b[jsp1 * b_dim1 + 1], &c__1, &
				f[is + f_dim1], ldf);
#line 882 "stgsy2.f"
			i__3 = js - 1;
#line 882 "stgsy2.f"
			saxpy_(&i__3, &rhs[2], &e[js * e_dim1 + 1], &c__1, &f[
				is + f_dim1], ldf);
#line 884 "stgsy2.f"
			i__3 = js - 1;
#line 884 "stgsy2.f"
			saxpy_(&i__3, &rhs[3], &e[jsp1 * e_dim1 + 1], &c__1, &
				f[is + f_dim1], ldf);
#line 886 "stgsy2.f"
		    }
#line 887 "stgsy2.f"
		    if (i__ < p) {
#line 888 "stgsy2.f"
			i__3 = *m - ie;
#line 888 "stgsy2.f"
			sger_(&i__3, &nb, &c_b27, &a[is + (ie + 1) * a_dim1], 
				lda, rhs, &c__1, &c__[ie + 1 + js * c_dim1], 
				ldc);
#line 890 "stgsy2.f"
			i__3 = *m - ie;
#line 890 "stgsy2.f"
			sger_(&i__3, &nb, &c_b27, &d__[is + (ie + 1) * d_dim1]
				, ldd, &rhs[2], &c__1, &c__[ie + 1 + js * 
				c_dim1], ldc);
#line 892 "stgsy2.f"
		    }

#line 894 "stgsy2.f"
		} else if (mb == 2 && nb == 1) {

/*                 Build a 4-by-4 system Z**T * x = RHS */

#line 898 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 899 "stgsy2.f"
		    z__[1] = a[is + isp1 * a_dim1];
#line 900 "stgsy2.f"
		    z__[2] = -b[js + js * b_dim1];
#line 901 "stgsy2.f"
		    z__[3] = 0.;

#line 903 "stgsy2.f"
		    z__[8] = a[isp1 + is * a_dim1];
#line 904 "stgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 905 "stgsy2.f"
		    z__[10] = 0.;
#line 906 "stgsy2.f"
		    z__[11] = -b[js + js * b_dim1];

#line 908 "stgsy2.f"
		    z__[16] = d__[is + is * d_dim1];
#line 909 "stgsy2.f"
		    z__[17] = d__[is + isp1 * d_dim1];
#line 910 "stgsy2.f"
		    z__[18] = -e[js + js * e_dim1];
#line 911 "stgsy2.f"
		    z__[19] = 0.;

#line 913 "stgsy2.f"
		    z__[24] = 0.;
#line 914 "stgsy2.f"
		    z__[25] = d__[isp1 + isp1 * d_dim1];
#line 915 "stgsy2.f"
		    z__[26] = 0.;
#line 916 "stgsy2.f"
		    z__[27] = -e[js + js * e_dim1];

/*                 Set up right hand side(s) */

#line 920 "stgsy2.f"
		    rhs[0] = c__[is + js * c_dim1];
#line 921 "stgsy2.f"
		    rhs[1] = c__[isp1 + js * c_dim1];
#line 922 "stgsy2.f"
		    rhs[2] = f[is + js * f_dim1];
#line 923 "stgsy2.f"
		    rhs[3] = f[isp1 + js * f_dim1];

/*                 Solve Z**T * x = RHS */

#line 927 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 928 "stgsy2.f"
		    if (ierr > 0) {
#line 928 "stgsy2.f"
			*info = ierr;
#line 928 "stgsy2.f"
		    }

#line 931 "stgsy2.f"
		    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 932 "stgsy2.f"
		    if (scaloc != 1.) {
#line 933 "stgsy2.f"
			i__3 = *n;
#line 933 "stgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 934 "stgsy2.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 935 "stgsy2.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 936 "stgsy2.f"
/* L150: */
#line 936 "stgsy2.f"
			}
#line 937 "stgsy2.f"
			*scale *= scaloc;
#line 938 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 942 "stgsy2.f"
		    c__[is + js * c_dim1] = rhs[0];
#line 943 "stgsy2.f"
		    c__[isp1 + js * c_dim1] = rhs[1];
#line 944 "stgsy2.f"
		    f[is + js * f_dim1] = rhs[2];
#line 945 "stgsy2.f"
		    f[isp1 + js * f_dim1] = rhs[3];

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 950 "stgsy2.f"
		    if (j > p + 2) {
#line 951 "stgsy2.f"
			i__3 = js - 1;
#line 951 "stgsy2.f"
			sger_(&mb, &i__3, &c_b42, rhs, &c__1, &b[js * b_dim1 
				+ 1], &c__1, &f[is + f_dim1], ldf);
#line 953 "stgsy2.f"
			i__3 = js - 1;
#line 953 "stgsy2.f"
			sger_(&mb, &i__3, &c_b42, &rhs[2], &c__1, &e[js * 
				e_dim1 + 1], &c__1, &f[is + f_dim1], ldf);
#line 955 "stgsy2.f"
		    }
#line 956 "stgsy2.f"
		    if (i__ < p) {
#line 957 "stgsy2.f"
			i__3 = *m - ie;
#line 957 "stgsy2.f"
			sgemv_("T", &mb, &i__3, &c_b27, &a[is + (ie + 1) * 
				a_dim1], lda, rhs, &c__1, &c_b42, &c__[ie + 1 
				+ js * c_dim1], &c__1, (ftnlen)1);
#line 960 "stgsy2.f"
			i__3 = *m - ie;
#line 960 "stgsy2.f"
			sgemv_("T", &mb, &i__3, &c_b27, &d__[is + (ie + 1) * 
				d_dim1], ldd, &rhs[2], &c__1, &c_b42, &c__[ie 
				+ 1 + js * c_dim1], &c__1, (ftnlen)1);
#line 963 "stgsy2.f"
		    }

#line 965 "stgsy2.f"
		} else if (mb == 2 && nb == 2) {

/*                 Build an 8-by-8 system Z**T * x = RHS */

#line 969 "stgsy2.f"
		    slaset_("F", &c__8, &c__8, &c_b56, &c_b56, z__, &c__8, (
			    ftnlen)1);

#line 971 "stgsy2.f"
		    z__[0] = a[is + is * a_dim1];
#line 972 "stgsy2.f"
		    z__[1] = a[is + isp1 * a_dim1];
#line 973 "stgsy2.f"
		    z__[4] = -b[js + js * b_dim1];
#line 974 "stgsy2.f"
		    z__[6] = -b[jsp1 + js * b_dim1];

#line 976 "stgsy2.f"
		    z__[8] = a[isp1 + is * a_dim1];
#line 977 "stgsy2.f"
		    z__[9] = a[isp1 + isp1 * a_dim1];
#line 978 "stgsy2.f"
		    z__[13] = -b[js + js * b_dim1];
#line 979 "stgsy2.f"
		    z__[15] = -b[jsp1 + js * b_dim1];

#line 981 "stgsy2.f"
		    z__[18] = a[is + is * a_dim1];
#line 982 "stgsy2.f"
		    z__[19] = a[is + isp1 * a_dim1];
#line 983 "stgsy2.f"
		    z__[20] = -b[js + jsp1 * b_dim1];
#line 984 "stgsy2.f"
		    z__[22] = -b[jsp1 + jsp1 * b_dim1];

#line 986 "stgsy2.f"
		    z__[26] = a[isp1 + is * a_dim1];
#line 987 "stgsy2.f"
		    z__[27] = a[isp1 + isp1 * a_dim1];
#line 988 "stgsy2.f"
		    z__[29] = -b[js + jsp1 * b_dim1];
#line 989 "stgsy2.f"
		    z__[31] = -b[jsp1 + jsp1 * b_dim1];

#line 991 "stgsy2.f"
		    z__[32] = d__[is + is * d_dim1];
#line 992 "stgsy2.f"
		    z__[33] = d__[is + isp1 * d_dim1];
#line 993 "stgsy2.f"
		    z__[36] = -e[js + js * e_dim1];

#line 995 "stgsy2.f"
		    z__[41] = d__[isp1 + isp1 * d_dim1];
#line 996 "stgsy2.f"
		    z__[45] = -e[js + js * e_dim1];

#line 998 "stgsy2.f"
		    z__[50] = d__[is + is * d_dim1];
#line 999 "stgsy2.f"
		    z__[51] = d__[is + isp1 * d_dim1];
#line 1000 "stgsy2.f"
		    z__[52] = -e[js + jsp1 * e_dim1];
#line 1001 "stgsy2.f"
		    z__[54] = -e[jsp1 + jsp1 * e_dim1];

#line 1003 "stgsy2.f"
		    z__[59] = d__[isp1 + isp1 * d_dim1];
#line 1004 "stgsy2.f"
		    z__[61] = -e[js + jsp1 * e_dim1];
#line 1005 "stgsy2.f"
		    z__[63] = -e[jsp1 + jsp1 * e_dim1];

/*                 Set up right hand side(s) */

#line 1009 "stgsy2.f"
		    k = 1;
#line 1010 "stgsy2.f"
		    ii = mb * nb + 1;
#line 1011 "stgsy2.f"
		    i__3 = nb - 1;
#line 1011 "stgsy2.f"
		    for (jj = 0; jj <= i__3; ++jj) {
#line 1012 "stgsy2.f"
			scopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &
				rhs[k - 1], &c__1);
#line 1013 "stgsy2.f"
			scopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[
				ii - 1], &c__1);
#line 1014 "stgsy2.f"
			k += mb;
#line 1015 "stgsy2.f"
			ii += mb;
#line 1016 "stgsy2.f"
/* L160: */
#line 1016 "stgsy2.f"
		    }


/*                 Solve Z**T * x = RHS */

#line 1021 "stgsy2.f"
		    sgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 1022 "stgsy2.f"
		    if (ierr > 0) {
#line 1022 "stgsy2.f"
			*info = ierr;
#line 1022 "stgsy2.f"
		    }

#line 1025 "stgsy2.f"
		    sgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 1026 "stgsy2.f"
		    if (scaloc != 1.) {
#line 1027 "stgsy2.f"
			i__3 = *n;
#line 1027 "stgsy2.f"
			for (k = 1; k <= i__3; ++k) {
#line 1028 "stgsy2.f"
			    sscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 1029 "stgsy2.f"
			    sscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 1030 "stgsy2.f"
/* L170: */
#line 1030 "stgsy2.f"
			}
#line 1031 "stgsy2.f"
			*scale *= scaloc;
#line 1032 "stgsy2.f"
		    }

/*                 Unpack solution vector(s) */

#line 1036 "stgsy2.f"
		    k = 1;
#line 1037 "stgsy2.f"
		    ii = mb * nb + 1;
#line 1038 "stgsy2.f"
		    i__3 = nb - 1;
#line 1038 "stgsy2.f"
		    for (jj = 0; jj <= i__3; ++jj) {
#line 1039 "stgsy2.f"
			scopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
				c_dim1], &c__1);
#line 1040 "stgsy2.f"
			scopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
				f_dim1], &c__1);
#line 1041 "stgsy2.f"
			k += mb;
#line 1042 "stgsy2.f"
			ii += mb;
#line 1043 "stgsy2.f"
/* L180: */
#line 1043 "stgsy2.f"
		    }

/*                 Substitute R(I, J) and L(I, J) into remaining */
/*                 equation. */

#line 1048 "stgsy2.f"
		    if (j > p + 2) {
#line 1049 "stgsy2.f"
			i__3 = js - 1;
#line 1049 "stgsy2.f"
			sgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &c__[is + 
				js * c_dim1], ldc, &b[js * b_dim1 + 1], ldb, &
				c_b42, &f[is + f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 1052 "stgsy2.f"
			i__3 = js - 1;
#line 1052 "stgsy2.f"
			sgemm_("N", "T", &mb, &i__3, &nb, &c_b42, &f[is + js *
				 f_dim1], ldf, &e[js * e_dim1 + 1], lde, &
				c_b42, &f[is + f_dim1], ldf, (ftnlen)1, (
				ftnlen)1);
#line 1055 "stgsy2.f"
		    }
#line 1056 "stgsy2.f"
		    if (i__ < p) {
#line 1057 "stgsy2.f"
			i__3 = *m - ie;
#line 1057 "stgsy2.f"
			sgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &a[is + (ie 
				+ 1) * a_dim1], lda, &c__[is + js * c_dim1], 
				ldc, &c_b42, &c__[ie + 1 + js * c_dim1], ldc, 
				(ftnlen)1, (ftnlen)1);
#line 1060 "stgsy2.f"
			i__3 = *m - ie;
#line 1060 "stgsy2.f"
			sgemm_("T", "N", &i__3, &nb, &mb, &c_b27, &d__[is + (
				ie + 1) * d_dim1], ldd, &f[is + js * f_dim1], 
				ldf, &c_b42, &c__[ie + 1 + js * c_dim1], ldc, 
				(ftnlen)1, (ftnlen)1);
#line 1063 "stgsy2.f"
		    }

#line 1065 "stgsy2.f"
		}

#line 1067 "stgsy2.f"
/* L190: */
#line 1067 "stgsy2.f"
	    }
#line 1068 "stgsy2.f"
/* L200: */
#line 1068 "stgsy2.f"
	}

#line 1070 "stgsy2.f"
    }
#line 1071 "stgsy2.f"
    return 0;

/*     End of STGSY2 */

} /* stgsy2_ */

