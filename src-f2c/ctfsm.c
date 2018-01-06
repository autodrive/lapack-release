#line 1 "ctfsm.f"
/* ctfsm.f -- translated by f2c (version 20100827).
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

#line 1 "ctfsm.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CTFSM solves a matrix equation (one operand is a triangular matrix in RFP format). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTFSM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctfsm.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctfsm.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctfsm.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, */
/*                         B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO */
/*       INTEGER            LDB, M, N */
/*       COMPLEX            ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( 0: * ), B( 0: LDB-1, 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for A in RFP Format. */
/* > */
/* > CTFSM solves the matrix equation */
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
/* >          ALPHA is COMPLEX */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (N*(N+1)/2) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int ctfsm_(char *transr, char *side, char *uplo, char *trans,
	 char *diag, integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *a, doublecomplex *b, integer *ldb, ftnlen transr_len, 
	ftnlen side_len, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, k, m1, m2, n1, n2, info;
    static logical normaltransr;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
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

#line 341 "ctfsm.f"
    /* Parameter adjustments */
#line 341 "ctfsm.f"
    b_dim1 = *ldb - 1 - 0 + 1;
#line 341 "ctfsm.f"
    b_offset = 0 + b_dim1 * 0;
#line 341 "ctfsm.f"
    b -= b_offset;
#line 341 "ctfsm.f"

#line 341 "ctfsm.f"
    /* Function Body */
#line 341 "ctfsm.f"
    info = 0;
#line 342 "ctfsm.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 343 "ctfsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 344 "ctfsm.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 345 "ctfsm.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 346 "ctfsm.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 347 "ctfsm.f"
	info = -1;
#line 348 "ctfsm.f"
    } else if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 349 "ctfsm.f"
	info = -2;
#line 350 "ctfsm.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 351 "ctfsm.f"
	info = -3;
#line 352 "ctfsm.f"
    } else if (! notrans && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 353 "ctfsm.f"
	info = -4;
#line 354 "ctfsm.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 356 "ctfsm.f"
	info = -5;
#line 357 "ctfsm.f"
    } else if (*m < 0) {
#line 358 "ctfsm.f"
	info = -6;
#line 359 "ctfsm.f"
    } else if (*n < 0) {
#line 360 "ctfsm.f"
	info = -7;
#line 361 "ctfsm.f"
    } else if (*ldb < max(1,*m)) {
#line 362 "ctfsm.f"
	info = -11;
#line 363 "ctfsm.f"
    }
#line 364 "ctfsm.f"
    if (info != 0) {
#line 365 "ctfsm.f"
	i__1 = -info;
#line 365 "ctfsm.f"
	xerbla_("CTFSM ", &i__1, (ftnlen)6);
#line 366 "ctfsm.f"
	return 0;
#line 367 "ctfsm.f"
    }

/*     Quick return when ( (N.EQ.0).OR.(M.EQ.0) ) */

#line 371 "ctfsm.f"
    if (*m == 0 || *n == 0) {
#line 371 "ctfsm.f"
	return 0;
#line 371 "ctfsm.f"
    }

/*     Quick return when ALPHA.EQ.(0E+0,0E+0) */

#line 376 "ctfsm.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 377 "ctfsm.f"
	i__1 = *n - 1;
#line 377 "ctfsm.f"
	for (j = 0; j <= i__1; ++j) {
#line 378 "ctfsm.f"
	    i__2 = *m - 1;
#line 378 "ctfsm.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {
#line 379 "ctfsm.f"
		i__3 = i__ + j * b_dim1;
#line 379 "ctfsm.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 380 "ctfsm.f"
/* L10: */
#line 380 "ctfsm.f"
	    }
#line 381 "ctfsm.f"
/* L20: */
#line 381 "ctfsm.f"
	}
#line 382 "ctfsm.f"
	return 0;
#line 383 "ctfsm.f"
    }

#line 385 "ctfsm.f"
    if (lside) {

/*        SIDE = 'L' */

/*        A is M-by-M. */
/*        If M is odd, set NISODD = .TRUE., and M1 and M2. */
/*        If M is even, NISODD = .FALSE., and M. */

#line 393 "ctfsm.f"
	if (*m % 2 == 0) {
#line 394 "ctfsm.f"
	    misodd = FALSE_;
#line 395 "ctfsm.f"
	    k = *m / 2;
#line 396 "ctfsm.f"
	} else {
#line 397 "ctfsm.f"
	    misodd = TRUE_;
#line 398 "ctfsm.f"
	    if (lower) {
#line 399 "ctfsm.f"
		m2 = *m / 2;
#line 400 "ctfsm.f"
		m1 = *m - m2;
#line 401 "ctfsm.f"
	    } else {
#line 402 "ctfsm.f"
		m1 = *m / 2;
#line 403 "ctfsm.f"
		m2 = *m - m1;
#line 404 "ctfsm.f"
	    }
#line 405 "ctfsm.f"
	}

#line 407 "ctfsm.f"
	if (misodd) {

/*           SIDE = 'L' and N is odd */

#line 411 "ctfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is odd, and TRANSR = 'N' */

#line 415 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 419 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 424 "ctfsm.f"
			if (*m == 1) {
#line 425 "ctfsm.f"
			    ctrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 427 "ctfsm.f"
			} else {
#line 428 "ctfsm.f"
			    ctrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 430 "ctfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 430 "ctfsm.f"
			    cgemm_("N", "N", &m2, n, &m1, &z__1, &a[m1], m, &
				    b[b_offset], ldb, alpha, &b[m1], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 432 "ctfsm.f"
			    ctrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 434 "ctfsm.f"
			}

#line 436 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 441 "ctfsm.f"
			if (*m == 1) {
#line 442 "ctfsm.f"
			    ctrsm_("L", "L", "C", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 444 "ctfsm.f"
			} else {
#line 445 "ctfsm.f"
			    ctrsm_("L", "U", "N", diag, &m2, n, alpha, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 447 "ctfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 447 "ctfsm.f"
			    cgemm_("C", "N", &m1, n, &m2, &z__1, &a[m1], m, &
				    b[m1], ldb, alpha, &b[b_offset], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 449 "ctfsm.f"
			    ctrsm_("L", "L", "C", diag, &m1, n, &c_b1, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 451 "ctfsm.f"
			}

#line 453 "ctfsm.f"
		    }

#line 455 "ctfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 459 "ctfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 464 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &m1, n, alpha, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 466 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 466 "ctfsm.f"
			cgemm_("C", "N", &m2, n, &m1, &z__1, a, m, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 468 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);

#line 471 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 476 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &m2, n, alpha, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);
#line 478 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 478 "ctfsm.f"
			cgemm_("N", "N", &m1, n, &m2, &z__1, a, m, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 480 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &m1, n, &c_b1, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 483 "ctfsm.f"
		    }

#line 485 "ctfsm.f"
		}

#line 487 "ctfsm.f"
	    } else {

/*              SIDE = 'L', N is odd, and TRANSR = 'C' */

#line 491 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 495 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 500 "ctfsm.f"
			if (*m == 1) {
#line 501 "ctfsm.f"
			    ctrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 503 "ctfsm.f"
			} else {
#line 504 "ctfsm.f"
			    ctrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 506 "ctfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 506 "ctfsm.f"
			    cgemm_("C", "N", &m2, n, &m1, &z__1, &a[m1 * m1], 
				    &m1, &b[b_offset], ldb, alpha, &b[m1], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 509 "ctfsm.f"
			    ctrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 511 "ctfsm.f"
			}

#line 513 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 518 "ctfsm.f"
			if (*m == 1) {
#line 519 "ctfsm.f"
			    ctrsm_("L", "U", "N", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 521 "ctfsm.f"
			} else {
#line 522 "ctfsm.f"
			    ctrsm_("L", "L", "C", diag, &m2, n, alpha, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 524 "ctfsm.f"
			    z__1.r = -1., z__1.i = -0.;
#line 524 "ctfsm.f"
			    cgemm_("N", "N", &m1, n, &m2, &z__1, &a[m1 * m1], 
				    &m1, &b[m1], ldb, alpha, &b[b_offset], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 527 "ctfsm.f"
			    ctrsm_("L", "U", "N", diag, &m1, n, &c_b1, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 529 "ctfsm.f"
			}

#line 531 "ctfsm.f"
		    }

#line 533 "ctfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 537 "ctfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 542 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &m1, n, alpha, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 544 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 544 "ctfsm.f"
			cgemm_("N", "N", &m2, n, &m1, &z__1, a, &m2, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 546 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 549 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 554 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &m2, n, alpha, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 556 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 556 "ctfsm.f"
			cgemm_("C", "N", &m1, n, &m2, &z__1, a, &m2, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 558 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &m1, n, &c_b1, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 561 "ctfsm.f"
		    }

#line 563 "ctfsm.f"
		}

#line 565 "ctfsm.f"
	    }

#line 567 "ctfsm.f"
	} else {

/*           SIDE = 'L' and N is even */

#line 571 "ctfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is even, and TRANSR = 'N' */

#line 575 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 579 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 584 "ctfsm.f"
			i__1 = *m + 1;
#line 584 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &k, n, alpha, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 586 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 586 "ctfsm.f"
			i__1 = *m + 1;
#line 586 "ctfsm.f"
			cgemm_("N", "N", &k, n, &k, &z__1, &a[k + 1], &i__1, &
				b[b_offset], ldb, alpha, &b[k], ldb, (ftnlen)
				1, (ftnlen)1);
#line 588 "ctfsm.f"
			i__1 = *m + 1;
#line 588 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &k, n, &c_b1, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 591 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 596 "ctfsm.f"
			i__1 = *m + 1;
#line 596 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &k, n, alpha, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 598 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 598 "ctfsm.f"
			i__1 = *m + 1;
#line 598 "ctfsm.f"
			cgemm_("C", "N", &k, n, &k, &z__1, &a[k + 1], &i__1, &
				b[k], ldb, alpha, &b[b_offset], ldb, (ftnlen)
				1, (ftnlen)1);
#line 600 "ctfsm.f"
			i__1 = *m + 1;
#line 600 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 603 "ctfsm.f"
		    }

#line 605 "ctfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 609 "ctfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 614 "ctfsm.f"
			i__1 = *m + 1;
#line 614 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &k, n, alpha, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 616 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 616 "ctfsm.f"
			i__1 = *m + 1;
#line 616 "ctfsm.f"
			cgemm_("C", "N", &k, n, &k, &z__1, a, &i__1, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 618 "ctfsm.f"
			i__1 = *m + 1;
#line 618 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &k, n, &c_b1, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 621 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'C' */
#line 625 "ctfsm.f"
			i__1 = *m + 1;
#line 625 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &k, n, alpha, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 627 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 627 "ctfsm.f"
			i__1 = *m + 1;
#line 627 "ctfsm.f"
			cgemm_("N", "N", &k, n, &k, &z__1, a, &i__1, &b[k], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 629 "ctfsm.f"
			i__1 = *m + 1;
#line 629 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 632 "ctfsm.f"
		    }

#line 634 "ctfsm.f"
		}

#line 636 "ctfsm.f"
	    } else {

/*              SIDE = 'L', N is even, and TRANSR = 'C' */

#line 640 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L' */

#line 644 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 649 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &k, n, alpha, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 651 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 651 "ctfsm.f"
			cgemm_("C", "N", &k, n, &k, &z__1, &a[k * (k + 1)], &
				k, &b[b_offset], ldb, alpha, &b[k], ldb, (
				ftnlen)1, (ftnlen)1);
#line 654 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &k, n, &c_b1, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 657 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 662 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &k, n, alpha, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 664 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 664 "ctfsm.f"
			cgemm_("N", "N", &k, n, &k, &z__1, &a[k * (k + 1)], &
				k, &b[k], ldb, alpha, &b[b_offset], ldb, (
				ftnlen)1, (ftnlen)1);
#line 667 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 670 "ctfsm.f"
		    }

#line 672 "ctfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U' */

#line 676 "ctfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 681 "ctfsm.f"
			ctrsm_("L", "U", "C", diag, &k, n, alpha, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 683 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 683 "ctfsm.f"
			cgemm_("N", "N", &k, n, &k, &z__1, a, &k, &b[b_offset]
				, ldb, alpha, &b[k], ldb, (ftnlen)1, (ftnlen)
				1);
#line 685 "ctfsm.f"
			ctrsm_("L", "L", "N", diag, &k, n, &c_b1, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 688 "ctfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 693 "ctfsm.f"
			ctrsm_("L", "L", "C", diag, &k, n, alpha, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 695 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 695 "ctfsm.f"
			cgemm_("C", "N", &k, n, &k, &z__1, a, &k, &b[k], ldb, 
				alpha, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1);
#line 697 "ctfsm.f"
			ctrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 700 "ctfsm.f"
		    }

#line 702 "ctfsm.f"
		}

#line 704 "ctfsm.f"
	    }

#line 706 "ctfsm.f"
	}

#line 708 "ctfsm.f"
    } else {

/*        SIDE = 'R' */

/*        A is N-by-N. */
/*        If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*        If N is even, NISODD = .FALSE., and K. */

#line 716 "ctfsm.f"
	if (*n % 2 == 0) {
#line 717 "ctfsm.f"
	    nisodd = FALSE_;
#line 718 "ctfsm.f"
	    k = *n / 2;
#line 719 "ctfsm.f"
	} else {
#line 720 "ctfsm.f"
	    nisodd = TRUE_;
#line 721 "ctfsm.f"
	    if (lower) {
#line 722 "ctfsm.f"
		n2 = *n / 2;
#line 723 "ctfsm.f"
		n1 = *n - n2;
#line 724 "ctfsm.f"
	    } else {
#line 725 "ctfsm.f"
		n1 = *n / 2;
#line 726 "ctfsm.f"
		n2 = *n - n1;
#line 727 "ctfsm.f"
	    }
#line 728 "ctfsm.f"
	}

#line 730 "ctfsm.f"
	if (nisodd) {

/*           SIDE = 'R' and N is odd */

#line 734 "ctfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is odd, and TRANSR = 'N' */

#line 738 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 742 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 747 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &n2, alpha, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 749 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 749 "ctfsm.f"
			cgemm_("N", "N", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, &a[n1], n, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 752 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &n1, &c_b1, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 755 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 760 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &n1, alpha, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 762 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 762 "ctfsm.f"
			cgemm_("N", "C", m, &n2, &n1, &z__1, b, ldb, &a[n1], 
				n, alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 765 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 768 "ctfsm.f"
		    }

#line 770 "ctfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 774 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 779 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &n1, alpha, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 781 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 781 "ctfsm.f"
			cgemm_("N", "N", m, &n2, &n1, &z__1, b, ldb, a, n, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 784 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 787 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 792 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &n2, alpha, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 794 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 794 "ctfsm.f"
			cgemm_("N", "C", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, a, n, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 796 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &n1, &c_b1, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 799 "ctfsm.f"
		    }

#line 801 "ctfsm.f"
		}

#line 803 "ctfsm.f"
	    } else {

/*              SIDE = 'R', N is odd, and TRANSR = 'C' */

#line 807 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 811 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 816 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &n2, alpha, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 818 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 818 "ctfsm.f"
			cgemm_("N", "C", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, &a[n1 * n1], &n1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 821 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &n1, &c_b1, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 824 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
/*                    TRANS = 'C' */

#line 829 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &n1, alpha, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 831 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 831 "ctfsm.f"
			cgemm_("N", "N", m, &n2, &n1, &z__1, b, ldb, &a[n1 * 
				n1], &n1, alpha, &b[n1 * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 834 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 837 "ctfsm.f"
		    }

#line 839 "ctfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 843 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 848 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &n1, alpha, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 850 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 850 "ctfsm.f"
			cgemm_("N", "C", m, &n2, &n1, &z__1, b, ldb, a, &n2, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 853 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 856 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
/*                    TRANS = 'C' */

#line 861 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &n2, alpha, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 863 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 863 "ctfsm.f"
			cgemm_("N", "N", m, &n1, &n2, &z__1, &b[n1 * b_dim1], 
				ldb, a, &n2, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 866 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &n1, &c_b1, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 869 "ctfsm.f"
		    }

#line 871 "ctfsm.f"
		}

#line 873 "ctfsm.f"
	    }

#line 875 "ctfsm.f"
	} else {

/*           SIDE = 'R' and N is even */

#line 879 "ctfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is even, and TRANSR = 'N' */

#line 883 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 887 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 892 "ctfsm.f"
			i__1 = *n + 1;
#line 892 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &k, alpha, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 894 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 894 "ctfsm.f"
			i__1 = *n + 1;
#line 894 "ctfsm.f"
			cgemm_("N", "N", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, &a[k + 1], &i__1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 897 "ctfsm.f"
			i__1 = *n + 1;
#line 897 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 900 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 905 "ctfsm.f"
			i__1 = *n + 1;
#line 905 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &k, alpha, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 907 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 907 "ctfsm.f"
			i__1 = *n + 1;
#line 907 "ctfsm.f"
			cgemm_("N", "C", m, &k, &k, &z__1, b, ldb, &a[k + 1], 
				&i__1, alpha, &b[k * b_dim1], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 910 "ctfsm.f"
			i__1 = *n + 1;
#line 910 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &k, &c_b1, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 913 "ctfsm.f"
		    }

#line 915 "ctfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 919 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 924 "ctfsm.f"
			i__1 = *n + 1;
#line 924 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &k, alpha, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 926 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 926 "ctfsm.f"
			i__1 = *n + 1;
#line 926 "ctfsm.f"
			cgemm_("N", "N", m, &k, &k, &z__1, b, ldb, a, &i__1, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 929 "ctfsm.f"
			i__1 = *n + 1;
#line 929 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &k, &c_b1, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 932 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 937 "ctfsm.f"
			i__1 = *n + 1;
#line 937 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &k, alpha, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 939 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 939 "ctfsm.f"
			i__1 = *n + 1;
#line 939 "ctfsm.f"
			cgemm_("N", "C", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, a, &i__1, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 942 "ctfsm.f"
			i__1 = *n + 1;
#line 942 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 945 "ctfsm.f"
		    }

#line 947 "ctfsm.f"
		}

#line 949 "ctfsm.f"
	    } else {

/*              SIDE = 'R', N is even, and TRANSR = 'C' */

#line 953 "ctfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L' */

#line 957 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 962 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &k, alpha, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 964 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 964 "ctfsm.f"
			cgemm_("N", "C", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, &a[(k + 1) * k], &k, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 967 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 970 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L', */
/*                    and TRANS = 'C' */

#line 975 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &k, alpha, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 977 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 977 "ctfsm.f"
			cgemm_("N", "N", m, &k, &k, &z__1, b, ldb, &a[(k + 1) 
				* k], &k, alpha, &b[k * b_dim1], ldb, (ftnlen)
				1, (ftnlen)1);
#line 980 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &k, &c_b1, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 983 "ctfsm.f"
		    }

#line 985 "ctfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U' */

#line 989 "ctfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 994 "ctfsm.f"
			ctrsm_("R", "U", "N", diag, m, &k, alpha, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 996 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 996 "ctfsm.f"
			cgemm_("N", "C", m, &k, &k, &z__1, b, ldb, a, &k, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 998 "ctfsm.f"
			ctrsm_("R", "L", "C", diag, m, &k, &c_b1, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);

#line 1001 "ctfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U', */
/*                    and TRANS = 'C' */

#line 1006 "ctfsm.f"
			ctrsm_("R", "L", "N", diag, m, &k, alpha, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);
#line 1008 "ctfsm.f"
			z__1.r = -1., z__1.i = -0.;
#line 1008 "ctfsm.f"
			cgemm_("N", "N", m, &k, &k, &z__1, &b[k * b_dim1], 
				ldb, a, &k, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 1010 "ctfsm.f"
			ctrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 1013 "ctfsm.f"
		    }

#line 1015 "ctfsm.f"
		}

#line 1017 "ctfsm.f"
	    }

#line 1019 "ctfsm.f"
	}
#line 1020 "ctfsm.f"
    }

#line 1022 "ctfsm.f"
    return 0;

/*     End of CTFSM */

} /* ctfsm_ */

