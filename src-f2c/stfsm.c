#line 1 "stfsm.f"
/* stfsm.f -- translated by f2c (version 20100827).
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

#line 1 "stfsm.f"
/* Table of constant values */

static doublereal c_b23 = -1.;
static doublereal c_b27 = 1.;

/* > \brief \b STFSM solves a matrix equation (one operand is a triangular matrix in RFP format). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STFSM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfsm.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfsm.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfsm.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, */
/*                         B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO */
/*       INTEGER            LDB, M, N */
/*       REAL               ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( 0: * ), B( 0: LDB-1, 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for A in RFP Format. */
/* > */
/* > STFSM  solves the matrix equation */
/* > */
/* >    op( A )*X = alpha*B  or  X*op( A ) = alpha*B */
/* > */
/* > where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T. */
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
/* >          = 'T':  The Transpose Form of RFP A is stored. */
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
/* >              TRANS  = 'T' or 't'   op( A ) = A'. */
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
/* >          ALPHA is REAL */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (NT) */
/* >           NT = N*(N+1)/2. On entry, the matrix A in RFP Format. */
/* >           RFP Format is described by TRANSR, UPLO and N as follows: */
/* >           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even; */
/* >           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If */
/* >           TRANSR = 'T' then RFP is the transpose of RFP A as */
/* >           defined when TRANSR = 'N'. The contents of RFP A are defined */
/* >           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT */
/* >           elements of upper packed A either in normal or */
/* >           transpose Format. If UPLO = 'L' the RFP A contains */
/* >           the NT elements of lower packed A either in normal or */
/* >           transpose Format. The LDA of RFP A is (N+1)/2 when */
/* >           TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is */
/* >           even and is N when is odd. */
/* >           See the Note below for more details. Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
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

/* > \date June 2017 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  We first consider Rectangular Full Packed (RFP) Format when N is */
/* >  even. We give an example where N = 6. */
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
/* >  the transpose of the first three columns of AP upper. */
/* >  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* >  three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* >  the transpose of the last three columns of AP lower. */
/* >  This covers the case N even and TRANSR = 'N'. */
/* > */
/* >         RFP A                   RFP A */
/* > */
/* >        03 04 05                33 43 53 */
/* >        13 14 15                00 44 54 */
/* >        23 24 25                10 11 55 */
/* >        33 34 35                20 21 22 */
/* >        00 44 45                30 31 32 */
/* >        01 11 55                40 41 42 */
/* >        02 12 22                50 51 52 */
/* > */
/* >  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* >  transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* >           RFP A                   RFP A */
/* > */
/* >     03 13 23 33 00 01 02    33 00 10 20 30 40 50 */
/* >     04 14 24 34 44 11 12    43 44 11 21 31 41 51 */
/* >     05 15 25 35 45 55 22    53 54 55 22 32 42 52 */
/* > */
/* > */
/* >  We then consider Rectangular Full Packed (RFP) Format when N is */
/* >  odd. We give an example where N = 5. */
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
/* >  the transpose of the first two columns of AP upper. */
/* >  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* >  three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* >  the transpose of the last two columns of AP lower. */
/* >  This covers the case N odd and TRANSR = 'N'. */
/* > */
/* >         RFP A                   RFP A */
/* > */
/* >        02 03 04                00 33 43 */
/* >        12 13 14                10 11 44 */
/* >        22 23 24                20 21 22 */
/* >        00 33 34                30 31 32 */
/* >        01 11 44                40 41 42 */
/* > */
/* >  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* >  transpose of RFP A above. One therefore gets: */
/* > */
/* >           RFP A                   RFP A */
/* > */
/* >     02 12 22 00 01             00 10 20 30 40 50 */
/* >     03 13 23 33 11             33 11 21 31 41 51 */
/* >     04 14 24 34 44             43 44 22 32 42 52 */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int stfsm_(char *transr, char *side, char *uplo, char *trans,
	 char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a,
	 doublereal *b, integer *ldb, ftnlen transr_len, ftnlen side_len, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m1, m2, n1, n2, info;
    static logical normaltransr, lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    static logical misodd, nisodd, notrans;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 320 "stfsm.f"
    /* Parameter adjustments */
#line 320 "stfsm.f"
    b_dim1 = *ldb - 1 - 0 + 1;
#line 320 "stfsm.f"
    b_offset = 0 + b_dim1 * 0;
#line 320 "stfsm.f"
    b -= b_offset;
#line 320 "stfsm.f"

#line 320 "stfsm.f"
    /* Function Body */
#line 320 "stfsm.f"
    info = 0;
#line 321 "stfsm.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 322 "stfsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 323 "stfsm.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 324 "stfsm.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 325 "stfsm.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 326 "stfsm.f"
	info = -1;
#line 327 "stfsm.f"
    } else if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 328 "stfsm.f"
	info = -2;
#line 329 "stfsm.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 330 "stfsm.f"
	info = -3;
#line 331 "stfsm.f"
    } else if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 332 "stfsm.f"
	info = -4;
#line 333 "stfsm.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 335 "stfsm.f"
	info = -5;
#line 336 "stfsm.f"
    } else if (*m < 0) {
#line 337 "stfsm.f"
	info = -6;
#line 338 "stfsm.f"
    } else if (*n < 0) {
#line 339 "stfsm.f"
	info = -7;
#line 340 "stfsm.f"
    } else if (*ldb < max(1,*m)) {
#line 341 "stfsm.f"
	info = -11;
#line 342 "stfsm.f"
    }
#line 343 "stfsm.f"
    if (info != 0) {
#line 344 "stfsm.f"
	i__1 = -info;
#line 344 "stfsm.f"
	xerbla_("STFSM ", &i__1, (ftnlen)6);
#line 345 "stfsm.f"
	return 0;
#line 346 "stfsm.f"
    }

/*     Quick return when ( (N.EQ.0).OR.(M.EQ.0) ) */

#line 350 "stfsm.f"
    if (*m == 0 || *n == 0) {
#line 350 "stfsm.f"
	return 0;
#line 350 "stfsm.f"
    }

/*     Quick return when ALPHA.EQ.(0D+0) */

#line 355 "stfsm.f"
    if (*alpha == 0.) {
#line 356 "stfsm.f"
	i__1 = *n - 1;
#line 356 "stfsm.f"
	for (j = 0; j <= i__1; ++j) {
#line 357 "stfsm.f"
	    i__2 = *m - 1;
#line 357 "stfsm.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {
#line 358 "stfsm.f"
		b[i__ + j * b_dim1] = 0.;
#line 359 "stfsm.f"
/* L10: */
#line 359 "stfsm.f"
	    }
#line 360 "stfsm.f"
/* L20: */
#line 360 "stfsm.f"
	}
#line 361 "stfsm.f"
	return 0;
#line 362 "stfsm.f"
    }

#line 364 "stfsm.f"
    if (lside) {

/*        SIDE = 'L' */

/*        A is M-by-M. */
/*        If M is odd, set NISODD = .TRUE., and M1 and M2. */
/*        If M is even, NISODD = .FALSE., and M. */

#line 372 "stfsm.f"
	if (*m % 2 == 0) {
#line 373 "stfsm.f"
	    misodd = FALSE_;
#line 374 "stfsm.f"
	    k = *m / 2;
#line 375 "stfsm.f"
	} else {
#line 376 "stfsm.f"
	    misodd = TRUE_;
#line 377 "stfsm.f"
	    if (lower) {
#line 378 "stfsm.f"
		m2 = *m / 2;
#line 379 "stfsm.f"
		m1 = *m - m2;
#line 380 "stfsm.f"
	    } else {
#line 381 "stfsm.f"
		m1 = *m / 2;
#line 382 "stfsm.f"
		m2 = *m - m1;
#line 383 "stfsm.f"
	    }
#line 384 "stfsm.f"
	}

#line 386 "stfsm.f"
	if (misodd) {

/*           SIDE = 'L' and N is odd */

#line 390 "stfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is odd, and TRANSR = 'N' */

#line 394 "stfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 398 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 403 "stfsm.f"
			if (*m == 1) {
#line 404 "stfsm.f"
			    strsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 406 "stfsm.f"
			} else {
#line 407 "stfsm.f"
			    strsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 409 "stfsm.f"
			    sgemm_("N", "N", &m2, n, &m1, &c_b23, &a[m1], m, &
				    b[b_offset], ldb, alpha, &b[m1], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 411 "stfsm.f"
			    strsm_("L", "U", "T", diag, &m2, n, &c_b27, &a[*m]
				    , m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 413 "stfsm.f"
			}

#line 415 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 420 "stfsm.f"
			if (*m == 1) {
#line 421 "stfsm.f"
			    strsm_("L", "L", "T", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 423 "stfsm.f"
			} else {
#line 424 "stfsm.f"
			    strsm_("L", "U", "N", diag, &m2, n, alpha, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 426 "stfsm.f"
			    sgemm_("T", "N", &m1, n, &m2, &c_b23, &a[m1], m, &
				    b[m1], ldb, alpha, &b[b_offset], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 428 "stfsm.f"
			    strsm_("L", "L", "T", diag, &m1, n, &c_b27, a, m, 
				    &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 430 "stfsm.f"
			}

#line 432 "stfsm.f"
		    }

#line 434 "stfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 438 "stfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 443 "stfsm.f"
			strsm_("L", "L", "N", diag, &m1, n, alpha, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 445 "stfsm.f"
			sgemm_("T", "N", &m2, n, &m1, &c_b23, a, m, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 447 "stfsm.f"
			strsm_("L", "U", "T", diag, &m2, n, &c_b27, &a[m1], m,
				 &b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 450 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 455 "stfsm.f"
			strsm_("L", "U", "N", diag, &m2, n, alpha, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);
#line 457 "stfsm.f"
			sgemm_("N", "N", &m1, n, &m2, &c_b23, a, m, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 459 "stfsm.f"
			strsm_("L", "L", "T", diag, &m1, n, &c_b27, &a[m2], m,
				 &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 462 "stfsm.f"
		    }

#line 464 "stfsm.f"
		}

#line 466 "stfsm.f"
	    } else {

/*              SIDE = 'L', N is odd, and TRANSR = 'T' */

#line 470 "stfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 474 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 479 "stfsm.f"
			if (*m == 1) {
#line 480 "stfsm.f"
			    strsm_("L", "U", "T", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 482 "stfsm.f"
			} else {
#line 483 "stfsm.f"
			    strsm_("L", "U", "T", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 485 "stfsm.f"
			    sgemm_("T", "N", &m2, n, &m1, &c_b23, &a[m1 * m1],
				     &m1, &b[b_offset], ldb, alpha, &b[m1], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 488 "stfsm.f"
			    strsm_("L", "L", "N", diag, &m2, n, &c_b27, &a[1],
				     &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 490 "stfsm.f"
			}

#line 492 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 497 "stfsm.f"
			if (*m == 1) {
#line 498 "stfsm.f"
			    strsm_("L", "U", "N", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 500 "stfsm.f"
			} else {
#line 501 "stfsm.f"
			    strsm_("L", "L", "T", diag, &m2, n, alpha, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 503 "stfsm.f"
			    sgemm_("N", "N", &m1, n, &m2, &c_b23, &a[m1 * m1],
				     &m1, &b[m1], ldb, alpha, &b[b_offset], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 506 "stfsm.f"
			    strsm_("L", "U", "N", diag, &m1, n, &c_b27, a, &
				    m1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				    1, (ftnlen)1, (ftnlen)1);
#line 508 "stfsm.f"
			}

#line 510 "stfsm.f"
		    }

#line 512 "stfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 516 "stfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 521 "stfsm.f"
			strsm_("L", "U", "T", diag, &m1, n, alpha, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 523 "stfsm.f"
			sgemm_("N", "N", &m2, n, &m1, &c_b23, a, &m2, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 525 "stfsm.f"
			strsm_("L", "L", "N", diag, &m2, n, &c_b27, &a[m1 * 
				m2], &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 528 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 533 "stfsm.f"
			strsm_("L", "L", "T", diag, &m2, n, alpha, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 535 "stfsm.f"
			sgemm_("T", "N", &m1, n, &m2, &c_b23, a, &m2, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 537 "stfsm.f"
			strsm_("L", "U", "N", diag, &m1, n, &c_b27, &a[m2 * 
				m2], &m2, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 540 "stfsm.f"
		    }

#line 542 "stfsm.f"
		}

#line 544 "stfsm.f"
	    }

#line 546 "stfsm.f"
	} else {

/*           SIDE = 'L' and N is even */

#line 550 "stfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is even, and TRANSR = 'N' */

#line 554 "stfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 558 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 563 "stfsm.f"
			i__1 = *m + 1;
#line 563 "stfsm.f"
			strsm_("L", "L", "N", diag, &k, n, alpha, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 565 "stfsm.f"
			i__1 = *m + 1;
#line 565 "stfsm.f"
			sgemm_("N", "N", &k, n, &k, &c_b23, &a[k + 1], &i__1, 
				&b[b_offset], ldb, alpha, &b[k], ldb, (ftnlen)
				1, (ftnlen)1);
#line 567 "stfsm.f"
			i__1 = *m + 1;
#line 567 "stfsm.f"
			strsm_("L", "U", "T", diag, &k, n, &c_b27, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 570 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 575 "stfsm.f"
			i__1 = *m + 1;
#line 575 "stfsm.f"
			strsm_("L", "U", "N", diag, &k, n, alpha, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 577 "stfsm.f"
			i__1 = *m + 1;
#line 577 "stfsm.f"
			sgemm_("T", "N", &k, n, &k, &c_b23, &a[k + 1], &i__1, 
				&b[k], ldb, alpha, &b[b_offset], ldb, (ftnlen)
				1, (ftnlen)1);
#line 579 "stfsm.f"
			i__1 = *m + 1;
#line 579 "stfsm.f"
			strsm_("L", "L", "T", diag, &k, n, &c_b27, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 582 "stfsm.f"
		    }

#line 584 "stfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 588 "stfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 593 "stfsm.f"
			i__1 = *m + 1;
#line 593 "stfsm.f"
			strsm_("L", "L", "N", diag, &k, n, alpha, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 595 "stfsm.f"
			i__1 = *m + 1;
#line 595 "stfsm.f"
			sgemm_("T", "N", &k, n, &k, &c_b23, a, &i__1, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 597 "stfsm.f"
			i__1 = *m + 1;
#line 597 "stfsm.f"
			strsm_("L", "U", "T", diag, &k, n, &c_b27, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 600 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'T' */
#line 604 "stfsm.f"
			i__1 = *m + 1;
#line 604 "stfsm.f"
			strsm_("L", "U", "N", diag, &k, n, alpha, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 606 "stfsm.f"
			i__1 = *m + 1;
#line 606 "stfsm.f"
			sgemm_("N", "N", &k, n, &k, &c_b23, a, &i__1, &b[k], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 608 "stfsm.f"
			i__1 = *m + 1;
#line 608 "stfsm.f"
			strsm_("L", "L", "T", diag, &k, n, &c_b27, &a[k + 1], 
				&i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 611 "stfsm.f"
		    }

#line 613 "stfsm.f"
		}

#line 615 "stfsm.f"
	    } else {

/*              SIDE = 'L', N is even, and TRANSR = 'T' */

#line 619 "stfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L' */

#line 623 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 628 "stfsm.f"
			strsm_("L", "U", "T", diag, &k, n, alpha, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 630 "stfsm.f"
			sgemm_("T", "N", &k, n, &k, &c_b23, &a[k * (k + 1)], &
				k, &b[b_offset], ldb, alpha, &b[k], ldb, (
				ftnlen)1, (ftnlen)1);
#line 633 "stfsm.f"
			strsm_("L", "L", "N", diag, &k, n, &c_b27, a, &k, &b[
				k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 636 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 641 "stfsm.f"
			strsm_("L", "L", "T", diag, &k, n, alpha, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 643 "stfsm.f"
			sgemm_("N", "N", &k, n, &k, &c_b23, &a[k * (k + 1)], &
				k, &b[k], ldb, alpha, &b[b_offset], ldb, (
				ftnlen)1, (ftnlen)1);
#line 646 "stfsm.f"
			strsm_("L", "U", "N", diag, &k, n, &c_b27, &a[k], &k, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 649 "stfsm.f"
		    }

#line 651 "stfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U' */

#line 655 "stfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 660 "stfsm.f"
			strsm_("L", "U", "T", diag, &k, n, alpha, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 662 "stfsm.f"
			sgemm_("N", "N", &k, n, &k, &c_b23, a, &k, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 664 "stfsm.f"
			strsm_("L", "L", "N", diag, &k, n, &c_b27, &a[k * k], 
				&k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 667 "stfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 672 "stfsm.f"
			strsm_("L", "L", "T", diag, &k, n, alpha, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 674 "stfsm.f"
			sgemm_("T", "N", &k, n, &k, &c_b23, a, &k, &b[k], ldb,
				 alpha, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1);
#line 676 "stfsm.f"
			strsm_("L", "U", "N", diag, &k, n, &c_b27, &a[k * (k 
				+ 1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 679 "stfsm.f"
		    }

#line 681 "stfsm.f"
		}

#line 683 "stfsm.f"
	    }

#line 685 "stfsm.f"
	}

#line 687 "stfsm.f"
    } else {

/*        SIDE = 'R' */

/*        A is N-by-N. */
/*        If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*        If N is even, NISODD = .FALSE., and K. */

#line 695 "stfsm.f"
	if (*n % 2 == 0) {
#line 696 "stfsm.f"
	    nisodd = FALSE_;
#line 697 "stfsm.f"
	    k = *n / 2;
#line 698 "stfsm.f"
	} else {
#line 699 "stfsm.f"
	    nisodd = TRUE_;
#line 700 "stfsm.f"
	    if (lower) {
#line 701 "stfsm.f"
		n2 = *n / 2;
#line 702 "stfsm.f"
		n1 = *n - n2;
#line 703 "stfsm.f"
	    } else {
#line 704 "stfsm.f"
		n1 = *n / 2;
#line 705 "stfsm.f"
		n2 = *n - n1;
#line 706 "stfsm.f"
	    }
#line 707 "stfsm.f"
	}

#line 709 "stfsm.f"
	if (nisodd) {

/*           SIDE = 'R' and N is odd */

#line 713 "stfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is odd, and TRANSR = 'N' */

#line 717 "stfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 721 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 726 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &n2, alpha, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 728 "stfsm.f"
			sgemm_("N", "N", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, &a[n1], n, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 731 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &n1, &c_b27, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 734 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 739 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &n1, alpha, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 741 "stfsm.f"
			sgemm_("N", "T", m, &n2, &n1, &c_b23, b, ldb, &a[n1], 
				n, alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 744 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &n2, &c_b27, &a[*n], n,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 747 "stfsm.f"
		    }

#line 749 "stfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 753 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 758 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &n1, alpha, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 760 "stfsm.f"
			sgemm_("N", "N", m, &n2, &n1, &c_b23, b, ldb, a, n, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 763 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &n2, &c_b27, &a[n1], n,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 766 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 771 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &n2, alpha, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 773 "stfsm.f"
			sgemm_("N", "T", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, a, n, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 775 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &n1, &c_b27, &a[n2], n,
				 b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 778 "stfsm.f"
		    }

#line 780 "stfsm.f"
		}

#line 782 "stfsm.f"
	    } else {

/*              SIDE = 'R', N is odd, and TRANSR = 'T' */

#line 786 "stfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 790 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 795 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &n2, alpha, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 797 "stfsm.f"
			sgemm_("N", "T", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, &a[n1 * n1], &n1, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 800 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &n1, &c_b27, a, &n1, b,
				 ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 803 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 808 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &n1, alpha, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 810 "stfsm.f"
			sgemm_("N", "N", m, &n2, &n1, &c_b23, b, ldb, &a[n1 * 
				n1], &n1, alpha, &b[n1 * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 813 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &n2, &c_b27, &a[1], &
				n1, &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 816 "stfsm.f"
		    }

#line 818 "stfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 822 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 827 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &n1, alpha, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 829 "stfsm.f"
			sgemm_("N", "T", m, &n2, &n1, &c_b23, b, ldb, a, &n2, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 832 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &n2, &c_b27, &a[n1 * 
				n2], &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 835 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 840 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &n2, alpha, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 842 "stfsm.f"
			sgemm_("N", "N", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, a, &n2, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 845 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &n1, &c_b27, &a[n2 * 
				n2], &n2, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 848 "stfsm.f"
		    }

#line 850 "stfsm.f"
		}

#line 852 "stfsm.f"
	    }

#line 854 "stfsm.f"
	} else {

/*           SIDE = 'R' and N is even */

#line 858 "stfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is even, and TRANSR = 'N' */

#line 862 "stfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 866 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 871 "stfsm.f"
			i__1 = *n + 1;
#line 871 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &k, alpha, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 873 "stfsm.f"
			i__1 = *n + 1;
#line 873 "stfsm.f"
			sgemm_("N", "N", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, &a[k + 1], &i__1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 876 "stfsm.f"
			i__1 = *n + 1;
#line 876 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &k, &c_b27, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 879 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 884 "stfsm.f"
			i__1 = *n + 1;
#line 884 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &k, alpha, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 886 "stfsm.f"
			i__1 = *n + 1;
#line 886 "stfsm.f"
			sgemm_("N", "T", m, &k, &k, &c_b23, b, ldb, &a[k + 1],
				 &i__1, alpha, &b[k * b_dim1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 889 "stfsm.f"
			i__1 = *n + 1;
#line 889 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &k, &c_b27, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 892 "stfsm.f"
		    }

#line 894 "stfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 898 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 903 "stfsm.f"
			i__1 = *n + 1;
#line 903 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &k, alpha, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 905 "stfsm.f"
			i__1 = *n + 1;
#line 905 "stfsm.f"
			sgemm_("N", "N", m, &k, &k, &c_b23, b, ldb, a, &i__1, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 908 "stfsm.f"
			i__1 = *n + 1;
#line 908 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &k, &c_b27, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 911 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 916 "stfsm.f"
			i__1 = *n + 1;
#line 916 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &k, alpha, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 918 "stfsm.f"
			i__1 = *n + 1;
#line 918 "stfsm.f"
			sgemm_("N", "T", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, a, &i__1, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 921 "stfsm.f"
			i__1 = *n + 1;
#line 921 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &k, &c_b27, &a[k + 1], 
				&i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 924 "stfsm.f"
		    }

#line 926 "stfsm.f"
		}

#line 928 "stfsm.f"
	    } else {

/*              SIDE = 'R', N is even, and TRANSR = 'T' */

#line 932 "stfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L' */

#line 936 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 941 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &k, alpha, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 943 "stfsm.f"
			sgemm_("N", "T", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, &a[(k + 1) * k], &k, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 946 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &k, &c_b27, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 949 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 954 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &k, alpha, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 956 "stfsm.f"
			sgemm_("N", "N", m, &k, &k, &c_b23, b, ldb, &a[(k + 1)
				 * k], &k, alpha, &b[k * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 959 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &k, &c_b27, a, &k, &b[
				k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 962 "stfsm.f"
		    }

#line 964 "stfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U' */

#line 968 "stfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 973 "stfsm.f"
			strsm_("R", "U", "N", diag, m, &k, alpha, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 975 "stfsm.f"
			sgemm_("N", "T", m, &k, &k, &c_b23, b, ldb, a, &k, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 977 "stfsm.f"
			strsm_("R", "L", "T", diag, m, &k, &c_b27, &a[k * k], 
				&k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 980 "stfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 985 "stfsm.f"
			strsm_("R", "L", "N", diag, m, &k, alpha, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);
#line 987 "stfsm.f"
			sgemm_("N", "N", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, a, &k, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 989 "stfsm.f"
			strsm_("R", "U", "T", diag, m, &k, &c_b27, &a[(k + 1) 
				* k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 992 "stfsm.f"
		    }

#line 994 "stfsm.f"
		}

#line 996 "stfsm.f"
	    }

#line 998 "stfsm.f"
	}
#line 999 "stfsm.f"
    }

#line 1001 "stfsm.f"
    return 0;

/*     End of STFSM */

} /* stfsm_ */

