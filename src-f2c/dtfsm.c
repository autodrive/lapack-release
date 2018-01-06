#line 1 "dtfsm.f"
/* dtfsm.f -- translated by f2c (version 20100827).
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

#line 1 "dtfsm.f"
/* Table of constant values */

static doublereal c_b23 = -1.;
static doublereal c_b27 = 1.;

/* > \brief \b DTFSM solves a matrix equation (one operand is a triangular matrix in RFP format). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTFSM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtfsm.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtfsm.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtfsm.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, */
/*                         B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO */
/*       INTEGER            LDB, M, N */
/*       DOUBLE PRECISION   ALPHA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( 0: * ), B( 0: LDB-1, 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for A in RFP Format. */
/* > */
/* > DTFSM  solves the matrix equation */
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
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (NT) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dtfsm_(char *transr, char *side, char *uplo, char *trans,
	 char *diag, integer *m, integer *n, doublereal *alpha, doublereal *a,
	 doublereal *b, integer *ldb, ftnlen transr_len, ftnlen side_len, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m1, m2, n1, n2, info;
    static logical normaltransr;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    static logical misodd, nisodd, notrans;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 320 "dtfsm.f"
    /* Parameter adjustments */
#line 320 "dtfsm.f"
    b_dim1 = *ldb - 1 - 0 + 1;
#line 320 "dtfsm.f"
    b_offset = 0 + b_dim1 * 0;
#line 320 "dtfsm.f"
    b -= b_offset;
#line 320 "dtfsm.f"

#line 320 "dtfsm.f"
    /* Function Body */
#line 320 "dtfsm.f"
    info = 0;
#line 321 "dtfsm.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 322 "dtfsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 323 "dtfsm.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 324 "dtfsm.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 325 "dtfsm.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 326 "dtfsm.f"
	info = -1;
#line 327 "dtfsm.f"
    } else if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 328 "dtfsm.f"
	info = -2;
#line 329 "dtfsm.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 330 "dtfsm.f"
	info = -3;
#line 331 "dtfsm.f"
    } else if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 332 "dtfsm.f"
	info = -4;
#line 333 "dtfsm.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 335 "dtfsm.f"
	info = -5;
#line 336 "dtfsm.f"
    } else if (*m < 0) {
#line 337 "dtfsm.f"
	info = -6;
#line 338 "dtfsm.f"
    } else if (*n < 0) {
#line 339 "dtfsm.f"
	info = -7;
#line 340 "dtfsm.f"
    } else if (*ldb < max(1,*m)) {
#line 341 "dtfsm.f"
	info = -11;
#line 342 "dtfsm.f"
    }
#line 343 "dtfsm.f"
    if (info != 0) {
#line 344 "dtfsm.f"
	i__1 = -info;
#line 344 "dtfsm.f"
	xerbla_("DTFSM ", &i__1, (ftnlen)6);
#line 345 "dtfsm.f"
	return 0;
#line 346 "dtfsm.f"
    }

/*     Quick return when ( (N.EQ.0).OR.(M.EQ.0) ) */

#line 350 "dtfsm.f"
    if (*m == 0 || *n == 0) {
#line 350 "dtfsm.f"
	return 0;
#line 350 "dtfsm.f"
    }

/*     Quick return when ALPHA.EQ.(0D+0) */

#line 355 "dtfsm.f"
    if (*alpha == 0.) {
#line 356 "dtfsm.f"
	i__1 = *n - 1;
#line 356 "dtfsm.f"
	for (j = 0; j <= i__1; ++j) {
#line 357 "dtfsm.f"
	    i__2 = *m - 1;
#line 357 "dtfsm.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {
#line 358 "dtfsm.f"
		b[i__ + j * b_dim1] = 0.;
#line 359 "dtfsm.f"
/* L10: */
#line 359 "dtfsm.f"
	    }
#line 360 "dtfsm.f"
/* L20: */
#line 360 "dtfsm.f"
	}
#line 361 "dtfsm.f"
	return 0;
#line 362 "dtfsm.f"
    }

#line 364 "dtfsm.f"
    if (lside) {

/*        SIDE = 'L' */

/*        A is M-by-M. */
/*        If M is odd, set NISODD = .TRUE., and M1 and M2. */
/*        If M is even, NISODD = .FALSE., and M. */

#line 372 "dtfsm.f"
	if (*m % 2 == 0) {
#line 373 "dtfsm.f"
	    misodd = FALSE_;
#line 374 "dtfsm.f"
	    k = *m / 2;
#line 375 "dtfsm.f"
	} else {
#line 376 "dtfsm.f"
	    misodd = TRUE_;
#line 377 "dtfsm.f"
	    if (lower) {
#line 378 "dtfsm.f"
		m2 = *m / 2;
#line 379 "dtfsm.f"
		m1 = *m - m2;
#line 380 "dtfsm.f"
	    } else {
#line 381 "dtfsm.f"
		m1 = *m / 2;
#line 382 "dtfsm.f"
		m2 = *m - m1;
#line 383 "dtfsm.f"
	    }
#line 384 "dtfsm.f"
	}


#line 387 "dtfsm.f"
	if (misodd) {

/*           SIDE = 'L' and N is odd */

#line 391 "dtfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is odd, and TRANSR = 'N' */

#line 395 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 399 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 404 "dtfsm.f"
			if (*m == 1) {
#line 405 "dtfsm.f"
			    dtrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 407 "dtfsm.f"
			} else {
#line 408 "dtfsm.f"
			    dtrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 410 "dtfsm.f"
			    dgemm_("N", "N", &m2, n, &m1, &c_b23, &a[m1], m, &
				    b[b_offset], ldb, alpha, &b[m1], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 412 "dtfsm.f"
			    dtrsm_("L", "U", "T", diag, &m2, n, &c_b27, &a[*m]
				    , m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 414 "dtfsm.f"
			}

#line 416 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 421 "dtfsm.f"
			if (*m == 1) {
#line 422 "dtfsm.f"
			    dtrsm_("L", "L", "T", diag, &m1, n, alpha, a, m, &
				    b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 424 "dtfsm.f"
			} else {
#line 425 "dtfsm.f"
			    dtrsm_("L", "U", "N", diag, &m2, n, alpha, &a[*m],
				     m, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 427 "dtfsm.f"
			    dgemm_("T", "N", &m1, n, &m2, &c_b23, &a[m1], m, &
				    b[m1], ldb, alpha, &b[b_offset], ldb, (
				    ftnlen)1, (ftnlen)1);
#line 429 "dtfsm.f"
			    dtrsm_("L", "L", "T", diag, &m1, n, &c_b27, a, m, 
				    &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 431 "dtfsm.f"
			}

#line 433 "dtfsm.f"
		    }

#line 435 "dtfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 439 "dtfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 444 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &m1, n, alpha, &a[m2], m, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 446 "dtfsm.f"
			dgemm_("T", "N", &m2, n, &m1, &c_b23, a, m, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 448 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &m2, n, &c_b27, &a[m1], m,
				 &b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 451 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 456 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &m2, n, alpha, &a[m1], m, 
				&b[m1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1);
#line 458 "dtfsm.f"
			dgemm_("N", "N", &m1, n, &m2, &c_b23, a, m, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 460 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &m1, n, &c_b27, &a[m2], m,
				 &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 463 "dtfsm.f"
		    }

#line 465 "dtfsm.f"
		}

#line 467 "dtfsm.f"
	    } else {

/*              SIDE = 'L', N is odd, and TRANSR = 'T' */

#line 471 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 475 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 480 "dtfsm.f"
			if (*m == 1) {
#line 481 "dtfsm.f"
			    dtrsm_("L", "U", "T", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 483 "dtfsm.f"
			} else {
#line 484 "dtfsm.f"
			    dtrsm_("L", "U", "T", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 486 "dtfsm.f"
			    dgemm_("T", "N", &m2, n, &m1, &c_b23, &a[m1 * m1],
				     &m1, &b[b_offset], ldb, alpha, &b[m1], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 489 "dtfsm.f"
			    dtrsm_("L", "L", "N", diag, &m2, n, &c_b27, &a[1],
				     &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 491 "dtfsm.f"
			}

#line 493 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 498 "dtfsm.f"
			if (*m == 1) {
#line 499 "dtfsm.f"
			    dtrsm_("L", "U", "N", diag, &m1, n, alpha, a, &m1,
				     &b[b_offset], ldb, (ftnlen)1, (ftnlen)1, 
				    (ftnlen)1, (ftnlen)1);
#line 501 "dtfsm.f"
			} else {
#line 502 "dtfsm.f"
			    dtrsm_("L", "L", "T", diag, &m2, n, alpha, &a[1], 
				    &m1, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				    ftnlen)1, (ftnlen)1);
#line 504 "dtfsm.f"
			    dgemm_("N", "N", &m1, n, &m2, &c_b23, &a[m1 * m1],
				     &m1, &b[m1], ldb, alpha, &b[b_offset], 
				    ldb, (ftnlen)1, (ftnlen)1);
#line 507 "dtfsm.f"
			    dtrsm_("L", "U", "N", diag, &m1, n, &c_b27, a, &
				    m1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				    1, (ftnlen)1, (ftnlen)1);
#line 509 "dtfsm.f"
			}

#line 511 "dtfsm.f"
		    }

#line 513 "dtfsm.f"
		} else {

/*                 SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 517 "dtfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 522 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &m1, n, alpha, &a[m2 * m2]
				, &m2, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 524 "dtfsm.f"
			dgemm_("N", "N", &m2, n, &m1, &c_b23, a, &m2, &b[
				b_offset], ldb, alpha, &b[m1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 526 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &m2, n, &c_b27, &a[m1 * 
				m2], &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 529 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 534 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &m2, n, alpha, &a[m1 * m2]
				, &m2, &b[m1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 536 "dtfsm.f"
			dgemm_("T", "N", &m1, n, &m2, &c_b23, a, &m2, &b[m1], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 538 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &m1, n, &c_b27, &a[m2 * 
				m2], &m2, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 541 "dtfsm.f"
		    }

#line 543 "dtfsm.f"
		}

#line 545 "dtfsm.f"
	    }

#line 547 "dtfsm.f"
	} else {

/*           SIDE = 'L' and N is even */

#line 551 "dtfsm.f"
	    if (normaltransr) {

/*              SIDE = 'L', N is even, and TRANSR = 'N' */

#line 555 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 559 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 564 "dtfsm.f"
			i__1 = *m + 1;
#line 564 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &k, n, alpha, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 566 "dtfsm.f"
			i__1 = *m + 1;
#line 566 "dtfsm.f"
			dgemm_("N", "N", &k, n, &k, &c_b23, &a[k + 1], &i__1, 
				&b[b_offset], ldb, alpha, &b[k], ldb, (ftnlen)
				1, (ftnlen)1);
#line 568 "dtfsm.f"
			i__1 = *m + 1;
#line 568 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &k, n, &c_b27, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 571 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 576 "dtfsm.f"
			i__1 = *m + 1;
#line 576 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &k, n, alpha, a, &i__1, &
				b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 578 "dtfsm.f"
			i__1 = *m + 1;
#line 578 "dtfsm.f"
			dgemm_("T", "N", &k, n, &k, &c_b23, &a[k + 1], &i__1, 
				&b[k], ldb, alpha, &b[b_offset], ldb, (ftnlen)
				1, (ftnlen)1);
#line 580 "dtfsm.f"
			i__1 = *m + 1;
#line 580 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &k, n, &c_b27, &a[1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 583 "dtfsm.f"
		    }

#line 585 "dtfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 589 "dtfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 594 "dtfsm.f"
			i__1 = *m + 1;
#line 594 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &k, n, alpha, &a[k + 1], &
				i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);
#line 596 "dtfsm.f"
			i__1 = *m + 1;
#line 596 "dtfsm.f"
			dgemm_("T", "N", &k, n, &k, &c_b23, a, &i__1, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 598 "dtfsm.f"
			i__1 = *m + 1;
#line 598 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &k, n, &c_b27, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 601 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'T' */
#line 605 "dtfsm.f"
			i__1 = *m + 1;
#line 605 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &k, n, alpha, &a[k], &
				i__1, &b[k], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 607 "dtfsm.f"
			i__1 = *m + 1;
#line 607 "dtfsm.f"
			dgemm_("N", "N", &k, n, &k, &c_b23, a, &i__1, &b[k], 
				ldb, alpha, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1);
#line 609 "dtfsm.f"
			i__1 = *m + 1;
#line 609 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &k, n, &c_b27, &a[k + 1], 
				&i__1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 612 "dtfsm.f"
		    }

#line 614 "dtfsm.f"
		}

#line 616 "dtfsm.f"
	    } else {

/*              SIDE = 'L', N is even, and TRANSR = 'T' */

#line 620 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L' */

#line 624 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 629 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &k, n, alpha, &a[k], &k, &
				b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 631 "dtfsm.f"
			dgemm_("T", "N", &k, n, &k, &c_b23, &a[k * (k + 1)], &
				k, &b[b_offset], ldb, alpha, &b[k], ldb, (
				ftnlen)1, (ftnlen)1);
#line 634 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &k, n, &c_b27, a, &k, &b[
				k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 637 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 642 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &k, n, alpha, a, &k, &b[k]
				, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 644 "dtfsm.f"
			dgemm_("N", "N", &k, n, &k, &c_b23, &a[k * (k + 1)], &
				k, &b[k], ldb, alpha, &b[b_offset], ldb, (
				ftnlen)1, (ftnlen)1);
#line 647 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &k, n, &c_b27, &a[k], &k, 
				&b[b_offset], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 650 "dtfsm.f"
		    }

#line 652 "dtfsm.f"
		} else {

/*                 SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U' */

#line 656 "dtfsm.f"
		    if (! notrans) {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 661 "dtfsm.f"
			dtrsm_("L", "U", "T", diag, &k, n, alpha, &a[k * (k + 
				1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 663 "dtfsm.f"
			dgemm_("N", "N", &k, n, &k, &c_b23, a, &k, &b[
				b_offset], ldb, alpha, &b[k], ldb, (ftnlen)1, 
				(ftnlen)1);
#line 665 "dtfsm.f"
			dtrsm_("L", "L", "N", diag, &k, n, &c_b27, &a[k * k], 
				&k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 668 "dtfsm.f"
		    } else {

/*                    SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 673 "dtfsm.f"
			dtrsm_("L", "L", "T", diag, &k, n, alpha, &a[k * k], &
				k, &b[k], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 675 "dtfsm.f"
			dgemm_("T", "N", &k, n, &k, &c_b23, a, &k, &b[k], ldb,
				 alpha, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
				1);
#line 677 "dtfsm.f"
			dtrsm_("L", "U", "N", diag, &k, n, &c_b27, &a[k * (k 
				+ 1)], &k, &b[b_offset], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 680 "dtfsm.f"
		    }

#line 682 "dtfsm.f"
		}

#line 684 "dtfsm.f"
	    }

#line 686 "dtfsm.f"
	}

#line 688 "dtfsm.f"
    } else {

/*        SIDE = 'R' */

/*        A is N-by-N. */
/*        If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*        If N is even, NISODD = .FALSE., and K. */

#line 696 "dtfsm.f"
	if (*n % 2 == 0) {
#line 697 "dtfsm.f"
	    nisodd = FALSE_;
#line 698 "dtfsm.f"
	    k = *n / 2;
#line 699 "dtfsm.f"
	} else {
#line 700 "dtfsm.f"
	    nisodd = TRUE_;
#line 701 "dtfsm.f"
	    if (lower) {
#line 702 "dtfsm.f"
		n2 = *n / 2;
#line 703 "dtfsm.f"
		n1 = *n - n2;
#line 704 "dtfsm.f"
	    } else {
#line 705 "dtfsm.f"
		n1 = *n / 2;
#line 706 "dtfsm.f"
		n2 = *n - n1;
#line 707 "dtfsm.f"
	    }
#line 708 "dtfsm.f"
	}

#line 710 "dtfsm.f"
	if (nisodd) {

/*           SIDE = 'R' and N is odd */

#line 714 "dtfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is odd, and TRANSR = 'N' */

#line 718 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 722 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 727 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &n2, alpha, &a[*n], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 729 "dtfsm.f"
			dgemm_("N", "N", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, &a[n1], n, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 732 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &n1, &c_b27, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);

#line 735 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 740 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &n1, alpha, a, n, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 742 "dtfsm.f"
			dgemm_("N", "T", m, &n2, &n1, &c_b23, b, ldb, &a[n1], 
				n, alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 745 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &n2, &c_b27, &a[*n], n,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 748 "dtfsm.f"
		    }

#line 750 "dtfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 754 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 759 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &n1, alpha, &a[n2], n, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 761 "dtfsm.f"
			dgemm_("N", "N", m, &n2, &n1, &c_b23, b, ldb, a, n, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 764 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &n2, &c_b27, &a[n1], n,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 767 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 772 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &n2, alpha, &a[n1], n, 
				&b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 774 "dtfsm.f"
			dgemm_("N", "T", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, a, n, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 776 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &n1, &c_b27, &a[n2], n,
				 b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 779 "dtfsm.f"
		    }

#line 781 "dtfsm.f"
		}

#line 783 "dtfsm.f"
	    } else {

/*              SIDE = 'R', N is odd, and TRANSR = 'T' */

#line 787 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 791 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'N' */

#line 796 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &n2, alpha, &a[1], &n1,
				 &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 798 "dtfsm.f"
			dgemm_("N", "T", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, &a[n1 * n1], &n1, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 801 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &n1, &c_b27, a, &n1, b,
				 ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 804 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and */
/*                    TRANS = 'T' */

#line 809 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &n1, alpha, a, &n1, b, 
				ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
				1);
#line 811 "dtfsm.f"
			dgemm_("N", "N", m, &n2, &n1, &c_b23, b, ldb, &a[n1 * 
				n1], &n1, alpha, &b[n1 * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 814 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &n2, &c_b27, &a[1], &
				n1, &b[n1 * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 817 "dtfsm.f"
		    }

#line 819 "dtfsm.f"
		} else {

/*                 SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 823 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'N' */

#line 828 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &n1, alpha, &a[n2 * n2]
				, &n2, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 830 "dtfsm.f"
			dgemm_("N", "T", m, &n2, &n1, &c_b23, b, ldb, a, &n2, 
				alpha, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 833 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &n2, &c_b27, &a[n1 * 
				n2], &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 836 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and */
/*                    TRANS = 'T' */

#line 841 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &n2, alpha, &a[n1 * n2]
				, &n2, &b[n1 * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 843 "dtfsm.f"
			dgemm_("N", "N", m, &n1, &n2, &c_b23, &b[n1 * b_dim1],
				 ldb, a, &n2, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 846 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &n1, &c_b27, &a[n2 * 
				n2], &n2, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 849 "dtfsm.f"
		    }

#line 851 "dtfsm.f"
		}

#line 853 "dtfsm.f"
	    }

#line 855 "dtfsm.f"
	} else {

/*           SIDE = 'R' and N is even */

#line 859 "dtfsm.f"
	    if (normaltransr) {

/*              SIDE = 'R', N is even, and TRANSR = 'N' */

#line 863 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L' */

#line 867 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 872 "dtfsm.f"
			i__1 = *n + 1;
#line 872 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &k, alpha, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 874 "dtfsm.f"
			i__1 = *n + 1;
#line 874 "dtfsm.f"
			dgemm_("N", "N", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, &a[k + 1], &i__1, alpha, b, ldb, (ftnlen)
				1, (ftnlen)1);
#line 877 "dtfsm.f"
			i__1 = *n + 1;
#line 877 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &k, &c_b27, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);

#line 880 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 885 "dtfsm.f"
			i__1 = *n + 1;
#line 885 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &k, alpha, &a[1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 887 "dtfsm.f"
			i__1 = *n + 1;
#line 887 "dtfsm.f"
			dgemm_("N", "T", m, &k, &k, &c_b23, b, ldb, &a[k + 1],
				 &i__1, alpha, &b[k * b_dim1], ldb, (ftnlen)1,
				 (ftnlen)1);
#line 890 "dtfsm.f"
			i__1 = *n + 1;
#line 890 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &k, &c_b27, a, &i__1, &
				b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 893 "dtfsm.f"
		    }

#line 895 "dtfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U' */

#line 899 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 904 "dtfsm.f"
			i__1 = *n + 1;
#line 904 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &k, alpha, &a[k + 1], &
				i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1);
#line 906 "dtfsm.f"
			i__1 = *n + 1;
#line 906 "dtfsm.f"
			dgemm_("N", "N", m, &k, &k, &c_b23, b, ldb, a, &i__1, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 909 "dtfsm.f"
			i__1 = *n + 1;
#line 909 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &k, &c_b27, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);

#line 912 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 917 "dtfsm.f"
			i__1 = *n + 1;
#line 917 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &k, alpha, &a[k], &
				i__1, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)
				1, (ftnlen)1, (ftnlen)1);
#line 919 "dtfsm.f"
			i__1 = *n + 1;
#line 919 "dtfsm.f"
			dgemm_("N", "T", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, a, &i__1, alpha, b, ldb, (ftnlen)1, (
				ftnlen)1);
#line 922 "dtfsm.f"
			i__1 = *n + 1;
#line 922 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &k, &c_b27, &a[k + 1], 
				&i__1, b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);

#line 925 "dtfsm.f"
		    }

#line 927 "dtfsm.f"
		}

#line 929 "dtfsm.f"
	    } else {

/*              SIDE = 'R', N is even, and TRANSR = 'T' */

#line 933 "dtfsm.f"
		if (lower) {

/*                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L' */

#line 937 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'N' */

#line 942 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &k, alpha, a, &k, &b[k 
				* b_dim1], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)
				1, (ftnlen)1);
#line 944 "dtfsm.f"
			dgemm_("N", "T", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, &a[(k + 1) * k], &k, alpha, b, ldb, (
				ftnlen)1, (ftnlen)1);
#line 947 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &k, &c_b27, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);

#line 950 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L', */
/*                    and TRANS = 'T' */

#line 955 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &k, alpha, &a[k], &k, 
				b, ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
				ftnlen)1);
#line 957 "dtfsm.f"
			dgemm_("N", "N", m, &k, &k, &c_b23, b, ldb, &a[(k + 1)
				 * k], &k, alpha, &b[k * b_dim1], ldb, (
				ftnlen)1, (ftnlen)1);
#line 960 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &k, &c_b27, a, &k, &b[
				k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 963 "dtfsm.f"
		    }

#line 965 "dtfsm.f"
		} else {

/*                 SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U' */

#line 969 "dtfsm.f"
		    if (notrans) {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'N' */

#line 974 "dtfsm.f"
			dtrsm_("R", "U", "N", diag, m, &k, alpha, &a[(k + 1) *
				 k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);
#line 976 "dtfsm.f"
			dgemm_("N", "T", m, &k, &k, &c_b23, b, ldb, a, &k, 
				alpha, &b[k * b_dim1], ldb, (ftnlen)1, (
				ftnlen)1);
#line 978 "dtfsm.f"
			dtrsm_("R", "L", "T", diag, m, &k, &c_b27, &a[k * k], 
				&k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1);

#line 981 "dtfsm.f"
		    } else {

/*                    SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U', */
/*                    and TRANS = 'T' */

#line 986 "dtfsm.f"
			dtrsm_("R", "L", "N", diag, m, &k, alpha, &a[k * k], &
				k, &b[k * b_dim1], ldb, (ftnlen)1, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1);
#line 988 "dtfsm.f"
			dgemm_("N", "N", m, &k, &k, &c_b23, &b[k * b_dim1], 
				ldb, a, &k, alpha, b, ldb, (ftnlen)1, (ftnlen)
				1);
#line 990 "dtfsm.f"
			dtrsm_("R", "U", "T", diag, m, &k, &c_b27, &a[(k + 1) 
				* k], &k, b, ldb, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

#line 993 "dtfsm.f"
		    }

#line 995 "dtfsm.f"
		}

#line 997 "dtfsm.f"
	    }

#line 999 "dtfsm.f"
	}
#line 1000 "dtfsm.f"
    }

#line 1002 "dtfsm.f"
    return 0;

/*     End of DTFSM */

} /* dtfsm_ */

