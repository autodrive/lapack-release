#line 1 "ztftri.f"
/* ztftri.f -- translated by f2c (version 20100827).
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

#line 1 "ztftri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZTFTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTFTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztftri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztftri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztftri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTFTRI( TRANSR, UPLO, DIAG, N, A, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO, DIAG */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTFTRI computes the inverse of a triangular matrix A stored in RFP */
/* > format. */
/* > */
/* > This is a Level 3 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  The Normal TRANSR of RFP A is stored; */
/* >          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( N*(N+1)/2 ); */
/* >          On entry, the triangular matrix A in RFP format. RFP format */
/* >          is described by TRANSR, UPLO, and N as follows: If TRANSR = */
/* >          'N' then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is */
/* >          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is */
/* >          the Conjugate-transpose of RFP A as defined when */
/* >          TRANSR = 'N'. The contents of RFP A are defined by UPLO as */
/* >          follows: If UPLO = 'U' the RFP A contains the nt elements of */
/* >          upper packed A; If UPLO = 'L' the RFP A contains the nt */
/* >          elements of lower packed A. The LDA of RFP A is (N+1)/2 when */
/* >          TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is */
/* >          even and N is odd. See the Note below for more details. */
/* > */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same storage format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular */
/* >               matrix is singular and its inverse can not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

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
/* Subroutine */ int ztftri_(char *transr, char *uplo, char *diag, integer *n,
	 doublecomplex *a, integer *info, ftnlen transr_len, ftnlen uplo_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Local variables */
    static integer k, n1, n2;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int ztrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical nisodd;
    extern /* Subroutine */ int ztrtri_(char *, char *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 261 "ztftri.f"
    *info = 0;
#line 262 "ztftri.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 263 "ztftri.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 264 "ztftri.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 265 "ztftri.f"
	*info = -1;
#line 266 "ztftri.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 267 "ztftri.f"
	*info = -2;
#line 268 "ztftri.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 270 "ztftri.f"
	*info = -3;
#line 271 "ztftri.f"
    } else if (*n < 0) {
#line 272 "ztftri.f"
	*info = -4;
#line 273 "ztftri.f"
    }
#line 274 "ztftri.f"
    if (*info != 0) {
#line 275 "ztftri.f"
	i__1 = -(*info);
#line 275 "ztftri.f"
	xerbla_("ZTFTRI", &i__1, (ftnlen)6);
#line 276 "ztftri.f"
	return 0;
#line 277 "ztftri.f"
    }

/*     Quick return if possible */

#line 281 "ztftri.f"
    if (*n == 0) {
#line 281 "ztftri.f"
	return 0;
#line 281 "ztftri.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

#line 287 "ztftri.f"
    if (*n % 2 == 0) {
#line 288 "ztftri.f"
	k = *n / 2;
#line 289 "ztftri.f"
	nisodd = FALSE_;
#line 290 "ztftri.f"
    } else {
#line 291 "ztftri.f"
	nisodd = TRUE_;
#line 292 "ztftri.f"
    }

/*     Set N1 and N2 depending on LOWER */

#line 296 "ztftri.f"
    if (lower) {
#line 297 "ztftri.f"
	n2 = *n / 2;
#line 298 "ztftri.f"
	n1 = *n - n2;
#line 299 "ztftri.f"
    } else {
#line 300 "ztftri.f"
	n1 = *n / 2;
#line 301 "ztftri.f"
	n2 = *n - n1;
#line 302 "ztftri.f"
    }


/*     start execution: there are eight cases */

#line 307 "ztftri.f"
    if (nisodd) {

/*        N is odd */

#line 311 "ztftri.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 315 "ztftri.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1) */

#line 321 "ztftri.f"
		ztrtri_("L", diag, &n1, a, n, info, (ftnlen)1, (ftnlen)1);
#line 322 "ztftri.f"
		if (*info > 0) {
#line 322 "ztftri.f"
		    return 0;
#line 322 "ztftri.f"
		}
#line 324 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 324 "ztftri.f"
		ztrmm_("R", "L", "N", diag, &n2, &n1, &z__1, a, n, &a[n1], n, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 326 "ztftri.f"
		ztrtri_("U", diag, &n2, &a[*n], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 327 "ztftri.f"
		if (*info > 0) {
#line 327 "ztftri.f"
		    *info += n1;
#line 327 "ztftri.f"
		}
#line 329 "ztftri.f"
		if (*info > 0) {
#line 329 "ztftri.f"
		    return 0;
#line 329 "ztftri.f"
		}
#line 331 "ztftri.f"
		ztrmm_("L", "U", "C", diag, &n2, &n1, &c_b1, &a[*n], n, &a[n1]
			, n, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 334 "ztftri.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 340 "ztftri.f"
		ztrtri_("L", diag, &n1, &a[n2], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 341 "ztftri.f"
		if (*info > 0) {
#line 341 "ztftri.f"
		    return 0;
#line 341 "ztftri.f"
		}
#line 343 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 343 "ztftri.f"
		ztrmm_("L", "L", "C", diag, &n1, &n2, &z__1, &a[n2], n, a, n, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 345 "ztftri.f"
		ztrtri_("U", diag, &n2, &a[n1], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 346 "ztftri.f"
		if (*info > 0) {
#line 346 "ztftri.f"
		    *info += n1;
#line 346 "ztftri.f"
		}
#line 348 "ztftri.f"
		if (*info > 0) {
#line 348 "ztftri.f"
		    return 0;
#line 348 "ztftri.f"
		}
#line 350 "ztftri.f"
		ztrmm_("R", "U", "N", diag, &n1, &n2, &c_b1, &a[n1], n, a, n, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 353 "ztftri.f"
	    }

#line 355 "ztftri.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 359 "ztftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1) */

#line 364 "ztftri.f"
		ztrtri_("U", diag, &n1, a, &n1, info, (ftnlen)1, (ftnlen)1);
#line 365 "ztftri.f"
		if (*info > 0) {
#line 365 "ztftri.f"
		    return 0;
#line 365 "ztftri.f"
		}
#line 367 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 367 "ztftri.f"
		ztrmm_("L", "U", "N", diag, &n1, &n2, &z__1, a, &n1, &a[n1 * 
			n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 369 "ztftri.f"
		ztrtri_("L", diag, &n2, &a[1], &n1, info, (ftnlen)1, (ftnlen)
			1);
#line 370 "ztftri.f"
		if (*info > 0) {
#line 370 "ztftri.f"
		    *info += n1;
#line 370 "ztftri.f"
		}
#line 372 "ztftri.f"
		if (*info > 0) {
#line 372 "ztftri.f"
		    return 0;
#line 372 "ztftri.f"
		}
#line 374 "ztftri.f"
		ztrmm_("R", "L", "C", diag, &n1, &n2, &c_b1, &a[1], &n1, &a[
			n1 * n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);

#line 377 "ztftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0) */

#line 382 "ztftri.f"
		ztrtri_("U", diag, &n1, &a[n2 * n2], &n2, info, (ftnlen)1, (
			ftnlen)1);
#line 383 "ztftri.f"
		if (*info > 0) {
#line 383 "ztftri.f"
		    return 0;
#line 383 "ztftri.f"
		}
#line 385 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 385 "ztftri.f"
		ztrmm_("R", "U", "C", diag, &n2, &n1, &z__1, &a[n2 * n2], &n2,
			 a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 387 "ztftri.f"
		ztrtri_("L", diag, &n2, &a[n1 * n2], &n2, info, (ftnlen)1, (
			ftnlen)1);
#line 388 "ztftri.f"
		if (*info > 0) {
#line 388 "ztftri.f"
		    *info += n1;
#line 388 "ztftri.f"
		}
#line 390 "ztftri.f"
		if (*info > 0) {
#line 390 "ztftri.f"
		    return 0;
#line 390 "ztftri.f"
		}
#line 392 "ztftri.f"
		ztrmm_("L", "L", "N", diag, &n2, &n1, &c_b1, &a[n1 * n2], &n2,
			 a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 394 "ztftri.f"
	    }

#line 396 "ztftri.f"
	}

#line 398 "ztftri.f"
    } else {

/*        N is even */

#line 402 "ztftri.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 406 "ztftri.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 412 "ztftri.f"
		i__1 = *n + 1;
#line 412 "ztftri.f"
		ztrtri_("L", diag, &k, &a[1], &i__1, info, (ftnlen)1, (ftnlen)
			1);
#line 413 "ztftri.f"
		if (*info > 0) {
#line 413 "ztftri.f"
		    return 0;
#line 413 "ztftri.f"
		}
#line 415 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 415 "ztftri.f"
		i__1 = *n + 1;
#line 415 "ztftri.f"
		i__2 = *n + 1;
#line 415 "ztftri.f"
		ztrmm_("R", "L", "N", diag, &k, &k, &z__1, &a[1], &i__1, &a[k 
			+ 1], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 417 "ztftri.f"
		i__1 = *n + 1;
#line 417 "ztftri.f"
		ztrtri_("U", diag, &k, a, &i__1, info, (ftnlen)1, (ftnlen)1);
#line 418 "ztftri.f"
		if (*info > 0) {
#line 418 "ztftri.f"
		    *info += k;
#line 418 "ztftri.f"
		}
#line 420 "ztftri.f"
		if (*info > 0) {
#line 420 "ztftri.f"
		    return 0;
#line 420 "ztftri.f"
		}
#line 422 "ztftri.f"
		i__1 = *n + 1;
#line 422 "ztftri.f"
		i__2 = *n + 1;
#line 422 "ztftri.f"
		ztrmm_("L", "U", "C", diag, &k, &k, &c_b1, a, &i__1, &a[k + 1]
			, &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 425 "ztftri.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 431 "ztftri.f"
		i__1 = *n + 1;
#line 431 "ztftri.f"
		ztrtri_("L", diag, &k, &a[k + 1], &i__1, info, (ftnlen)1, (
			ftnlen)1);
#line 432 "ztftri.f"
		if (*info > 0) {
#line 432 "ztftri.f"
		    return 0;
#line 432 "ztftri.f"
		}
#line 434 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 434 "ztftri.f"
		i__1 = *n + 1;
#line 434 "ztftri.f"
		i__2 = *n + 1;
#line 434 "ztftri.f"
		ztrmm_("L", "L", "C", diag, &k, &k, &z__1, &a[k + 1], &i__1, 
			a, &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 436 "ztftri.f"
		i__1 = *n + 1;
#line 436 "ztftri.f"
		ztrtri_("U", diag, &k, &a[k], &i__1, info, (ftnlen)1, (ftnlen)
			1);
#line 437 "ztftri.f"
		if (*info > 0) {
#line 437 "ztftri.f"
		    *info += k;
#line 437 "ztftri.f"
		}
#line 439 "ztftri.f"
		if (*info > 0) {
#line 439 "ztftri.f"
		    return 0;
#line 439 "ztftri.f"
		}
#line 441 "ztftri.f"
		i__1 = *n + 1;
#line 441 "ztftri.f"
		i__2 = *n + 1;
#line 441 "ztftri.f"
		ztrmm_("R", "U", "N", diag, &k, &k, &c_b1, &a[k], &i__1, a, &
			i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 443 "ztftri.f"
	    }
#line 444 "ztftri.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 448 "ztftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 454 "ztftri.f"
		ztrtri_("U", diag, &k, &a[k], &k, info, (ftnlen)1, (ftnlen)1);
#line 455 "ztftri.f"
		if (*info > 0) {
#line 455 "ztftri.f"
		    return 0;
#line 455 "ztftri.f"
		}
#line 457 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 457 "ztftri.f"
		ztrmm_("L", "U", "N", diag, &k, &k, &z__1, &a[k], &k, &a[k * (
			k + 1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 459 "ztftri.f"
		ztrtri_("L", diag, &k, a, &k, info, (ftnlen)1, (ftnlen)1);
#line 460 "ztftri.f"
		if (*info > 0) {
#line 460 "ztftri.f"
		    *info += k;
#line 460 "ztftri.f"
		}
#line 462 "ztftri.f"
		if (*info > 0) {
#line 462 "ztftri.f"
		    return 0;
#line 462 "ztftri.f"
		}
#line 464 "ztftri.f"
		ztrmm_("R", "L", "C", diag, &k, &k, &c_b1, a, &k, &a[k * (k + 
			1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 466 "ztftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 472 "ztftri.f"
		ztrtri_("U", diag, &k, &a[k * (k + 1)], &k, info, (ftnlen)1, (
			ftnlen)1);
#line 473 "ztftri.f"
		if (*info > 0) {
#line 473 "ztftri.f"
		    return 0;
#line 473 "ztftri.f"
		}
#line 475 "ztftri.f"
		z__1.r = -1., z__1.i = -0.;
#line 475 "ztftri.f"
		ztrmm_("R", "U", "C", diag, &k, &k, &z__1, &a[k * (k + 1)], &
			k, a, &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 477 "ztftri.f"
		ztrtri_("L", diag, &k, &a[k * k], &k, info, (ftnlen)1, (
			ftnlen)1);
#line 478 "ztftri.f"
		if (*info > 0) {
#line 478 "ztftri.f"
		    *info += k;
#line 478 "ztftri.f"
		}
#line 480 "ztftri.f"
		if (*info > 0) {
#line 480 "ztftri.f"
		    return 0;
#line 480 "ztftri.f"
		}
#line 482 "ztftri.f"
		ztrmm_("L", "L", "N", diag, &k, &k, &c_b1, &a[k * k], &k, a, &
			k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 484 "ztftri.f"
	    }
#line 485 "ztftri.f"
	}
#line 486 "ztftri.f"
    }

#line 488 "ztftri.f"
    return 0;

/*     End of ZTFTRI */

} /* ztftri_ */

