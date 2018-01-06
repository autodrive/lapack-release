#line 1 "ztfsm.f"
/* ztfsm.f -- translated by f2c (version 20100827).
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

#line 1 "ztfsm.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZTFSM solves a matrix equation (one operand is a triangular matrix in RFP format). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTFSM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztfsm.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztfsm.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztfsm.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, */
/*                         B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO */
/*       INTEGER            LDB, M, N */
/*       COMPLEX*16         ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( 0: * ), B( 0: LDB-1, 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for A in RFP Format. */
/* > */
/* > ZTFSM  solves the matrix equation */
/* > */
/* >    op( A )*X = alpha*B  or  X*op( A ) = alpha*B */
/* > */
/* > where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**H. */
/* > */
/* > A is in Rectangular Full Packed (RFP) Format. */
/* > */
/* > The matrix X is overwritten on B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  The Normal Form of RFP A is stored; */
/* >          = 'C':  The Conjugate-transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry, SIDE specifies whether op( A ) appears on the left */
/* >           or right of X as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */
/* > */
/* >              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the RFP matrix A came from */
/* >           an upper or lower triangular matrix as follows: */
/* >           UPLO = 'U' or 'u' RFP A came from an upper triangular matrix */
/* >           UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS  specifies the form of op( A ) to be used */
/* >           in the matrix multiplication as follows: */
/* > */
/* >              TRANS  = 'N' or 'n'   op( A ) = A. */
/* > */
/* >              TRANS  = 'C' or 'c'   op( A ) = conjg( A' ). */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not RFP A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of B. M must be at */
/* >           least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of B.  N must be */
/* >           at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >           NT = N*(N+1)/2. On entry, the matrix A in RFP Format. */
/* >           RFP Format is described by TRANSR, UPLO and N as follows: */
/* >           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even; */
/* >           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If */
/* >           TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A as */
/* >           defined when TRANSR = 'N'. The contents of RFP A are defined */
/* >           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT */
/* >           elements of upper packed A either in normal or */
/* >           conjugate-transpose Format. If UPLO = 'L' the RFP A contains */
/* >           the NT elements of lower packed A either in normal or */
/* >           conjugate-transpose Format. The LDA of RFP A is (N+1)/2 when */
/* >           TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is */
/* >           even and is N when is odd. */
/* >           See the Note below for more details. Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain  the  right-hand  side  matrix  B,  and  on exit  is */
/* >           overwritten by the solution matrix  X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           max( 1, m ). */
/* >           Unchanged on exit. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  We first consider Standard Packed Format when N is even. */
/* >  We give an example where N = 6. */
/* > */
/* >      AP is Upper             AP is Lower */
/* > */
/* >   00 01 02 03 04 05       00 */
/* >      11 12 13 14 15       10 11 */
/* >         22 23 24 25       20 21 22 */
/* >            33 34 35       30 31 32 33 */
/* >               44 45       40 41 42 43 44 */
/* >                  55       50 51 52 53 54 55 */
/* > */
/* > */
/* >  Let TRANSR = 'N'. RFP holds AP as follows: */
/* >  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last */
/* >  three columns of AP upper. The lower triangle A(4:6,0:2) consists of */
/* >  conjugate-transpose of the first three columns of AP upper. */
/* >  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* >  three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* >  conjugate-transpose of the last three columns of AP lower. */
/* >  To denote conjugate we place -- above the element. This covers the */
/* >  case N even and TRANSR = 'N'. */
/* > */
/* >         RFP A                   RFP A */
/* > */
/* >                                -- -- -- */
/* >        03 04 05                33 43 53 */
/* >                                   -- -- */
/* >        13 14 15                00 44 54 */
/* >                                      -- */
/* >        23 24 25                10 11 55 */
/* > */
/* >        33 34 35                20 21 22 */
/* >        -- */
/* >        00 44 45                30 31 32 */
/* >        -- -- */
/* >        01 11 55                40 41 42 */
/* >        -- -- -- */
/* >        02 12 22                50 51 52 */
/* > */
/* >  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* >  transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* >           RFP A                   RFP A */
/* > */
/* >     -- -- -- --                -- -- -- -- -- -- */
/* >     03 13 23 33 00 01 02    33 00 10 20 30 40 50 */
/* >     -- -- -- -- --                -- -- -- -- -- */
/* >     04 14 24 34 44 11 12    43 44 11 21 31 41 51 */
/* >     -- -- -- -- -- --                -- -- -- -- */
/* >     05 15 25 35 45 55 22    53 54 55 22 32 42 52 */
/* > */
/* > */
/* >  We next  consider Standard Packed Format when N is odd. */
/* >  We give an example where N = 5. */
/* > */
/* >     AP is Upper                 AP is Lower */
/* > */
/* >   00 01 02 03 04              00 */
/* >      11 12 13 14              10 11 */
/* >         22 23 24              20 21 22 */
/* >            33 34              30 31 32 33 */
/* >               44              40 41 42 43 44 */
/* > */
/* > */
/* >  Let TRANSR = 'N'. RFP holds AP as follows: */
/* >  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last */
/* >  three columns of AP upper. The lower triangle A(3:4,0:1) consists of */
/* >  conjugate-transpose of the first two   columns of AP upper. */
/* >  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* >  three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* >  conjugate-transpose of the last two   columns of AP lower. */
/* >  To denote conjugate we place -- above the element. This covers the */
/* >  case N odd  and TRANSR = 'N'. */
/* > */
/* >         RFP A                   RFP A */
/* > */
/* >                                   -- -- */
/* >        02 03 04                00 33 43 */
/* >                                      -- */
/* >        12 13 14                10 11 44 */
/* > */
/* >        22 23 24                20 21 22 */
/* >        -- */
/* >        00 33 34                30 31 32 */
/* >        -- -- */
/* >        01 11 44                40 41 42 */
/* > */
/* >  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* >  transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* >           RFP A                   RFP A */
/* > */
/* >     -- -- --                   -- -- -- -- -- -- */
/* >     02 12 22 00 01             00 10 20 30 40 50 */
/* >     -- -- -- --                   -- -- -- -- -- */
/* >     03 13 23 33 11             33 11 21 31 41 51 */
/* >     -- -- -- -- --                   -- -- -- -- */
/* >     04 14 24 34 44             43 44 22 32 42 52 */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztfsm_(char *transr, char *side, char *uplo, char *trans,
	 char *diag, integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *a, doublecomplex *b, integer *ldb, ftnlen transr_len, 
	ftnlen side_len, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, k, m1, m2, n1, n2, info;
    static logical normaltransr, lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical misodd, nisodd, notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*     .. */
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

/*     Test the input parameters. */

#line 341 "ztfsm.f"
    /* Parameter adjustments */
#line 341 "ztfsm.f"
    b_dim1 = *ldb - 1 - 0 + 1;
#line 341 "ztfsm.f"
    b_offset = 0 + b_dim1 * 0;
#line 341 "ztfsm.f"
    b -= b_offset;
#line 341 "ztfsm.f"

#line 341 "ztfsm.f"
    /* Function Body */
#line 341 "ztfsm.f"
    info = 0;
#line 342 "ztfsm.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 343 "ztfsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 344 "ztfsm.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 345 "ztfsm.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 346 "ztfsm.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 347 "ztfsm.f"
	info = -1;
#line 348 "ztfsm.f"
    } else if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 349 "ztfsm.f"
	info = -2;
#line 350 "ztfsm.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 351 "ztfsm.f"
	info = -3;
#line 352 "ztfsm.f"
    } else if (! notrans && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 353 "ztfsm.f"
	info = -4;
#line 354 "ztfsm.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 356 "ztfsm.f"
	info = -5;
#line 357 "ztfsm.f"
    } else if (*m < 0) {
#line 358 "ztfsm.f"
	info = -6;
#line 359 "ztfsm.f"
    } else if (*n < 0) {
#line 360 "ztfsm.f"
	info = -7;
#line 361 "ztfsm.f"
    } else if (*ldb < max(1,*m)) {
#line 362 "ztfsm.f"
	info = -11;
#line 363 "ztfsm.f"
    }
#line 364 "ztfsm.f"
    if (info != 0) {
#line 365 "ztfsm.f"
	i__1 = -info;
#line 365 "ztfsm.f"
	xerbla_("ZTFSM ", &i__1, (ftnlen)6);
#line 366 "ztfsm.f"
	return 0;
#line 367 "ztfsm.f"
    }

/*     Quick return when ( (N.EQ.0).OR.(M.EQ.0) ) */

#line 371 "ztfsm.f"
    if (*m == 0 || *n == 0) {
#line 371 "ztfsm.f"
	return 0;
#line 371 "ztfsm.f"
    }

/*     Quick return when ALPHA.EQ.(0D+0,0D+0) */

#line 376 "ztfsm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 377 "ztfsm.f"
	i__1 = *n - 1;
#line 377 "ztfsm.f"
	for (j = 0; j <= i__1; ++j) {
#line 378 "ztfsm.f"
	    i__2 = *m - 1;
#line 378 "ztfsm.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {
#line 379 "ztfsm.f"
		i__3 = i__ + j * b_dim1;
#line 379 "ztfsm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 380 "ztfsm.f"
/* L10: */
#line 380 "ztfsm.f"
	    }
#line 381 "ztfsm.f"
/* L20: */
#line 381 "ztfsm.f"
	}
#line 382 "ztfsm.f"
	return 0;
#line 383 "ztfsm.f"
    }

#line 385 "ztfsm.f"
    if (lside) {

/*        SIDE = 'L' */

/*        A is M-by-M. */
/*        If M is odd, set NISODD = .TRUE., and M1 and M2. */
/*        If M is even, NISODD = .FALSE., and M. */

#line 393 "ztfsm.f"
	if (*m % 2 == 0) {
#line 394 "ztfsm.f"
	    misodd = FALSE_;
#line 395 "ztfsm.f"
	    k = *m / 2;
#line 396 "ztfsm.f"
	} else {
#line 397 "ztfsm.f"
	    misodd = TRUE_;
#line 398 "ztfsm.f"
	    if (lower) {
#line 399 "ztfsm.f"
		m2 = *m / 2;
#line 400 "ztfsm.f"
		m1 = *m - m2;
#line 401 "ztfsm.f"
	    } else {
#line 402 "ztfsm.f"
		m1 = *m / 2;
#line 403 "ztfsm.f"
		m2 = *m - m1;
#line 404 "ztfsm.f"
	    }
#line 405 "ztfsm.f"
	}

#line 407 "ztfsm.f"
	if (misodd) {

/*           SIDE = 'L' and N is odd */

#line 411 "ztfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is odd, and TRANSR = 'N' */

#line 415 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 419 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 424 "ztfsm.f"
			if (*m == 1) {
#line 425 "ztfsm.f"
			    ztrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 427 "ztfsm.f"
			} else {
#line 428 "ztfsm.f"
			    ztrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 430 "ztfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 430 "ztfsm.f"
			    zgemm_("N", "N", &m2, n, &m1, &z__1, &a[m1], m, &
				    b[b_offset], ldb, alpha, &b[m1], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 432 "ztfsm.f"
			    ztrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 434 "ztfsm.f"
			}

#line 436 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 441 "ztfsm.f"
			if (*m == 1) {
#line 442 "ztfsm.f"
			    ztrsm_("L", "L", "C", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 444 "ztfsm.f"
			} else {
#line 445 "ztfsm.f"
			    ztrsm_("L", "U", "N", diag, &m2, n, alpha, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 447 "ztfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 447 "ztfsm.f"
			    zgemm_("C", "N", &m1, n, &m2, &z__1, &a[m1], m, &
				    b[m1], ldb, alpha, &b[b_offset], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 449 "ztfsm.f"
			    ztrsm_("L", "L", "C", diag, &m1, n, &c_b1, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 451 "ztfsm.f"
			}

#line 453 "ztfsm.f"
		    }

#line 455 "ztfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 459 "ztfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 464 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &m1, n, alpha, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 466 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 466 "ztfsm.f"
			zgemm_("C", "N", &m2, n, &m1, &z__1, a, m, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 468 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);

#line 471 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 476 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &m2, n, alpha, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);
#line 478 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 478 "ztfsm.f"
			zgemm_("N", "N", &m1, n, &m2, &z__1, a, m, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 480 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &m1, n, &c_b1, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 483 "ztfsm.f"
		    }

#line 485 "ztfsm.f"
		}

#line 487 "ztfsm.f"
	    } else {

/*              SIDE = 'L', N is odd, and TRANSR = 'C' */

#line 491 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 495 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 500 "ztfsm.f"
			if (*m == 1) {
#line 501 "ztfsm.f"
			    ztrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 503 "ztfsm.f"
			} else {
#line 504 "ztfsm.f"
			    ztrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 506 "ztfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 506 "ztfsm.f"
			    zgemm_("C", "N", &m2, n, &m1, &z__1, &a[m1 * m1], 
				    &m1, &b[b_offset], ldb, alpha, &b[m1], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 509 "ztfsm.f"
			    ztrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 511 "ztfsm.f"
			}

#line 513 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 518 "ztfsm.f"
			if (*m == 1) {
#line 519 "ztfsm.f"
			    ztrsm_("L", "U", "N", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 521 "ztfsm.f"
			} else {
#line 522 "ztfsm.f"
			    ztrsm_("L", "L", "C", diag, &m2, n, alpha, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 524 "ztfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 524 "ztfsm.f"
			    zgemm_("N", "N", &m1, n, &m2, &z__1, &a[m1 * m1], 
				    &m1, &b[m1], ldb, alpha, &b[b_offset], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 527 "ztfsm.f"
			    ztrsm_("L", "U", "N", diag, &m1, n, &c_b1, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 529 "ztfsm.f"
			}

#line 531 "ztfsm.f"
		    }

#line 533 "ztfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 537 "ztfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 542 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &m1, n, alpha, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 544 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 544 "ztfsm.f"
			zgemm_("N", "N", &m2, n, &m1, &z__1, a, &m2, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 546 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 549 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 554 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &m2, n, alpha, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 556 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 556 "ztfsm.f"
			zgemm_("C", "N", &m1, n, &m2, &z__1, a, &m2, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 558 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &m1, n, &c_b1, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 561 "ztfsm.f"
		    }

#line 563 "ztfsm.f"
		}

#line 565 "ztfsm.f"
	    }

#line 567 "ztfsm.f"
	} else {

/*           SIDE = 'L' and N is even */

#line 571 "ztfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is even, and TRANSR = 'N' */

#line 575 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 579 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 584 "ztfsm.f"
			i__1 = *m + 1;
#line 584 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &k, n, alpha, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 586 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 586 "ztfsm.f"
			i__1 = *m + 1;
#line 586 "ztfsm.f"
			zgemm_("N", "N", &k, n, &k, &z__1, &a[k + 1], &i__1, &
				b[b_offset], ldb, alpha, &b[k], ldb, (ftnlen)
				1, (ftnlen)1);
#line 588 "ztfsm.f"
			i__1 = *m + 1;
#line 588 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &k, n, &c_b1, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 591 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 596 "ztfsm.f"
			i__1 = *m + 1;
#line 596 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &k, n, alpha, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 598 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 598 "ztfsm.f"
			i__1 = *m + 1;
#line 598 "ztfsm.f"
			zgemm_("C", "N", &k, n, &k, &z__1, &a[k + 1], &i__1, &
				b[k], ldb, alpha, &b[b_offset], ldb, (ftnlen)
				1, (ftnlen)1);
#line 600 "ztfsm.f"
			i__1 = *m + 1;
#line 600 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 603 "ztfsm.f"
		    }

#line 605 "ztfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 609 "ztfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 614 "ztfsm.f"
			i__1 = *m + 1;
#line 614 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &k, n, alpha, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 616 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 616 "ztfsm.f"
			i__1 = *m + 1;
#line 616 "ztfsm.f"
			zgemm_("C", "N", &k, n, &k, &z__1, a, &i__1, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 618 "ztfsm.f"
			i__1 = *m + 1;
#line 618 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &k, n, &c_b1, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 621 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'C' */
#line 625 "ztfsm.f"
			i__1 = *m + 1;
#line 625 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &k, n, alpha, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 627 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 627 "ztfsm.f"
			i__1 = *m + 1;
#line 627 "ztfsm.f"
			zgemm_("N", "N", &k, n, &k, &z__1, a, &i__1, &b[k], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 629 "ztfsm.f"
			i__1 = *m + 1;
#line 629 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 632 "ztfsm.f"
		    }

#line 634 "ztfsm.f"
		}

#line 636 "ztfsm.f"
	    } else {

/*              SIDE = 'L', N is even, and TRANSR = 'C' */

#line 640 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L' */

#line 644 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 649 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &k, n, alpha, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 651 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 651 "ztfsm.f"
			zgemm_("C", "N", &k, n, &k, &z__1, &a[k * (k + 1)], &
				k, &b[b_offset], ldb, alpha, &b[k], ldb, (
				ftnlen)1, (ftnlen)1);
#line 654 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &k, n, &c_b1, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 657 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 662 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &k, n, alpha, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 664 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 664 "ztfsm.f"
			zgemm_("N", "N", &k, n, &k, &z__1, &a[k * (k + 1)], &
				k, &b[k], ldb, alpha, &b[b_offset], ldb, (
				ftnlen)1, (ftnlen)1);
#line 667 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 670 "ztfsm.f"
		    }

#line 672 "ztfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U' */

#line 676 "ztfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 681 "ztfsm.f"
			ztrsm_("L", "U", "C", diag, &k, n, alpha, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 683 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 683 "ztfsm.f"
			zgemm_("N", "N", &k, n, &k, &z__1, a, &k, &b[b_offset]
				, ldb, alpha, &b[k], ldb, (ftnlen)1, (ftnlen)
				1);
#line 685 "ztfsm.f"
			ztrsm_("L", "L", "N", diag, &k, n, &c_b1, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 688 "ztfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 693 "ztfsm.f"
			ztrsm_("L", "L", "C", diag, &k, n, alpha, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 695 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 695 "ztfsm.f"
			zgemm_("C", "N", &k, n, &k, &z__1, a, &k, &b[k], ldb, 
				alpha, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1);
#line 697 "ztfsm.f"
			ztrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 700 "ztfsm.f"
		    }

#line 702 "ztfsm.f"
		}

#line 704 "ztfsm.f"
	    }

#line 706 "ztfsm.f"
	}

#line 708 "ztfsm.f"
    } else {

/*        SIDE = 'R' */

/*        A is N-by-N. */
/*        If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*        If N is even, NISODD = .FALSE., and K. */

#line 716 "ztfsm.f"
	if (*n % 2 == 0) {
#line 717 "ztfsm.f"
	    nisodd = FALSE_;
#line 718 "ztfsm.f"
	    k = *n / 2;
#line 719 "ztfsm.f"
	} else {
#line 720 "ztfsm.f"
	    nisodd = TRUE_;
#line 721 "ztfsm.f"
	    if (lower) {
#line 722 "ztfsm.f"
		n2 = *n / 2;
#line 723 "ztfsm.f"
		n1 = *n - n2;
#line 724 "ztfsm.f"
	    } else {
#line 725 "ztfsm.f"
		n1 = *n / 2;
#line 726 "ztfsm.f"
		n2 = *n - n1;
#line 727 "ztfsm.f"
	    }
#line 728 "ztfsm.f"
	}

#line 730 "ztfsm.f"
	if (nisodd) {

/*           SIDE = 'R' and N is odd */

#line 734 "ztfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is odd, and TRANSR = 'N' */

#line 738 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 742 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 747 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &n2, alpha, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 749 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 749 "ztfsm.f"
			zgemm_("N", "N", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, &a[n1], n, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 752 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &n1, &c_b1, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 755 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 760 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &n1, alpha, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 762 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 762 "ztfsm.f"
			zgemm_("N", "C", m, &n2, &n1, &z__1, b, ldb, &a[n1], 
				n, alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 765 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 768 "ztfsm.f"
		    }

#line 770 "ztfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 774 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 779 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &n1, alpha, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 781 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 781 "ztfsm.f"
			zgemm_("N", "N", m, &n2, &n1, &z__1, b, ldb, a, n, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 784 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 787 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 792 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &n2, alpha, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 794 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 794 "ztfsm.f"
			zgemm_("N", "C", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, a, n, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 796 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &n1, &c_b1, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 799 "ztfsm.f"
		    }

#line 801 "ztfsm.f"
		}

#line 803 "ztfsm.f"
	    } else {

/*              SIDE = 'R', N is odd, and TRANSR = 'C' */

#line 807 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 811 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 816 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &n2, alpha, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 818 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 818 "ztfsm.f"
			zgemm_("N", "C", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, &a[n1 * n1], &n1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 821 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &n1, &c_b1, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 824 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 829 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &n1, alpha, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 831 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 831 "ztfsm.f"
			zgemm_("N", "N", m, &n2, &n1, &z__1, b, ldb, &a[n1 * 
				n1], &n1, alpha, &b[n1 * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 834 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 837 "ztfsm.f"
		    }

#line 839 "ztfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 843 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 848 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &n1, alpha, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 850 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 850 "ztfsm.f"
			zgemm_("N", "C", m, &n2, &n1, &z__1, b, ldb, a, &n2, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 853 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 856 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 861 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &n2, alpha, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 863 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 863 "ztfsm.f"
			zgemm_("N", "N", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, a, &n2, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 866 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &n1, &c_b1, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 869 "ztfsm.f"
		    }

#line 871 "ztfsm.f"
		}

#line 873 "ztfsm.f"
	    }

#line 875 "ztfsm.f"
	} else {

/*           SIDE = 'R' and N is even */

#line 879 "ztfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is even, and TRANSR = 'N' */

#line 883 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 887 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 892 "ztfsm.f"
			i__1 = *n + 1;
#line 892 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &k, alpha, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 894 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 894 "ztfsm.f"
			i__1 = *n + 1;
#line 894 "ztfsm.f"
			zgemm_("N", "N", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, &a[k + 1], &i__1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 897 "ztfsm.f"
			i__1 = *n + 1;
#line 897 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 900 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 905 "ztfsm.f"
			i__1 = *n + 1;
#line 905 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &k, alpha, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 907 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 907 "ztfsm.f"
			i__1 = *n + 1;
#line 907 "ztfsm.f"
			zgemm_("N", "C", m, &k, &k, &z__1, b, ldb, &a[k + 1], 
				&i__1, alpha, &b[k * b_dim1], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 910 "ztfsm.f"
			i__1 = *n + 1;
#line 910 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &k, &c_b1, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 913 "ztfsm.f"
		    }

#line 915 "ztfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 919 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 924 "ztfsm.f"
			i__1 = *n + 1;
#line 924 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &k, alpha, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 926 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 926 "ztfsm.f"
			i__1 = *n + 1;
#line 926 "ztfsm.f"
			zgemm_("N", "N", m, &k, &k, &z__1, b, ldb, a, &i__1, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 929 "ztfsm.f"
			i__1 = *n + 1;
#line 929 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &k, &c_b1, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 932 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 937 "ztfsm.f"
			i__1 = *n + 1;
#line 937 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &k, alpha, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 939 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 939 "ztfsm.f"
			i__1 = *n + 1;
#line 939 "ztfsm.f"
			zgemm_("N", "C", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, a, &i__1, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 942 "ztfsm.f"
			i__1 = *n + 1;
#line 942 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 945 "ztfsm.f"
		    }

#line 947 "ztfsm.f"
		}

#line 949 "ztfsm.f"
	    } else {

/*              SIDE = 'R', N is even, and TRANSR = 'C' */

#line 953 "ztfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L' */

#line 957 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 962 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &k, alpha, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 964 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 964 "ztfsm.f"
			zgemm_("N", "C", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, &a[(k + 1) * k], &k, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 967 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 970 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 975 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &k, alpha, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 977 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 977 "ztfsm.f"
			zgemm_("N", "N", m, &k, &k, &z__1, b, ldb, &a[(k + 1) 
				* k], &k, alpha, &b[k * b_dim1], ldb, (ftnlen)
				1, (ftnlen)1);
#line 980 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &k, &c_b1, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 983 "ztfsm.f"
		    }

#line 985 "ztfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U' */

#line 989 "ztfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 994 "ztfsm.f"
			ztrsm_("R", "U", "N", diag, m, &k, alpha, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 996 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 996 "ztfsm.f"
			zgemm_("N", "C", m, &k, &k, &z__1, b, ldb, a, &k, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 998 "ztfsm.f"
			ztrsm_("R", "L", "C", diag, m, &k, &c_b1, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);

#line 1001 "ztfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 1006 "ztfsm.f"
			ztrsm_("R", "L", "N", diag, m, &k, alpha, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);
#line 1008 "ztfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 1008 "ztfsm.f"
			zgemm_("N", "N", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, a, &k, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 1010 "ztfsm.f"
			ztrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 1013 "ztfsm.f"
		    }

#line 1015 "ztfsm.f"
		}

#line 1017 "ztfsm.f"
	    }

#line 1019 "ztfsm.f"
	}
#line 1020 "ztfsm.f"
    }

#line 1022 "ztfsm.f"
    return 0;

/*     End of ZTFSM */

} /* ztfsm_ */

