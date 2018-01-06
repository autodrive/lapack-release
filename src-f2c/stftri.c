#line 1 "stftri.f"
/* stftri.f -- translated by f2c (version 20100827).
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

#line 1 "stftri.f"
/* Table of constant values */

static doublereal c_b13 = -1.;
static doublereal c_b18 = 1.;

/* > \brief \b STFTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STFTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stftri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stftri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stftri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STFTRI( TRANSR, UPLO, DIAG, N, A, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO, DIAG */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STFTRI computes the inverse of a triangular matrix A stored in RFP */
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
/* >          = 'T':  The Transpose TRANSR of RFP A is stored. */
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
/* >          A is REAL array, dimension (NT); */
/* >          NT=N*(N+1)/2. On entry, the triangular factor of a Hermitian */
/* >          Positive Definite matrix A in RFP format. RFP format is */
/* >          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N' */
/* >          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is */
/* >          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is */
/* >          the transpose of RFP A as defined when */
/* >          TRANSR = 'N'. The contents of RFP A are defined by UPLO as */
/* >          follows: If UPLO = 'U' the RFP A contains the nt elements of */
/* >          upper packed A; If UPLO = 'L' the RFP A contains the nt */
/* >          elements of lower packed A. The LDA of RFP A is (N+1)/2 when */
/* >          TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is */
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

/* > \date November 2011 */

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
/* > */
/*  ===================================================================== */
/* Subroutine */ int stftri_(char *transr, char *uplo, char *diag, integer *n,
	 doublereal *a, integer *info, ftnlen transr_len, ftnlen uplo_len, 
	ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, n1, n2;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen);
    static logical nisodd;
    extern /* Subroutine */ int strtri_(char *, char *, integer *, doublereal 
	    *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 241 "stftri.f"
    *info = 0;
#line 242 "stftri.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 243 "stftri.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 244 "stftri.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 245 "stftri.f"
	*info = -1;
#line 246 "stftri.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 247 "stftri.f"
	*info = -2;
#line 248 "stftri.f"
    } else if (! lsame_(diag, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "U", (ftnlen)1, (ftnlen)1)) {
#line 250 "stftri.f"
	*info = -3;
#line 251 "stftri.f"
    } else if (*n < 0) {
#line 252 "stftri.f"
	*info = -4;
#line 253 "stftri.f"
    }
#line 254 "stftri.f"
    if (*info != 0) {
#line 255 "stftri.f"
	i__1 = -(*info);
#line 255 "stftri.f"
	xerbla_("STFTRI", &i__1, (ftnlen)6);
#line 256 "stftri.f"
	return 0;
#line 257 "stftri.f"
    }

/*     Quick return if possible */

#line 261 "stftri.f"
    if (*n == 0) {
#line 261 "stftri.f"
	return 0;
#line 261 "stftri.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

#line 267 "stftri.f"
    if (*n % 2 == 0) {
#line 268 "stftri.f"
	k = *n / 2;
#line 269 "stftri.f"
	nisodd = FALSE_;
#line 270 "stftri.f"
    } else {
#line 271 "stftri.f"
	nisodd = TRUE_;
#line 272 "stftri.f"
    }

/*     Set N1 and N2 depending on LOWER */

#line 276 "stftri.f"
    if (lower) {
#line 277 "stftri.f"
	n2 = *n / 2;
#line 278 "stftri.f"
	n1 = *n - n2;
#line 279 "stftri.f"
    } else {
#line 280 "stftri.f"
	n1 = *n / 2;
#line 281 "stftri.f"
	n2 = *n - n1;
#line 282 "stftri.f"
    }


/*     start execution: there are eight cases */

#line 287 "stftri.f"
    if (nisodd) {

/*        N is odd */

#line 291 "stftri.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 295 "stftri.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1) */

#line 301 "stftri.f"
		strtri_("L", diag, &n1, a, n, info, (ftnlen)1, (ftnlen)1);
#line 302 "stftri.f"
		if (*info > 0) {
#line 302 "stftri.f"
		    return 0;
#line 302 "stftri.f"
		}
#line 304 "stftri.f"
		strmm_("R", "L", "N", diag, &n2, &n1, &c_b13, a, n, &a[n1], n,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 306 "stftri.f"
		strtri_("U", diag, &n2, &a[*n], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 307 "stftri.f"
		if (*info > 0) {
#line 307 "stftri.f"
		    *info += n1;
#line 307 "stftri.f"
		}
#line 309 "stftri.f"
		if (*info > 0) {
#line 309 "stftri.f"
		    return 0;
#line 309 "stftri.f"
		}
#line 311 "stftri.f"
		strmm_("L", "U", "T", diag, &n2, &n1, &c_b18, &a[*n], n, &a[
			n1], n, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 314 "stftri.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 320 "stftri.f"
		strtri_("L", diag, &n1, &a[n2], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 321 "stftri.f"
		if (*info > 0) {
#line 321 "stftri.f"
		    return 0;
#line 321 "stftri.f"
		}
#line 323 "stftri.f"
		strmm_("L", "L", "T", diag, &n1, &n2, &c_b13, &a[n2], n, a, n,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 325 "stftri.f"
		strtri_("U", diag, &n2, &a[n1], n, info, (ftnlen)1, (ftnlen)1)
			;
#line 326 "stftri.f"
		if (*info > 0) {
#line 326 "stftri.f"
		    *info += n1;
#line 326 "stftri.f"
		}
#line 328 "stftri.f"
		if (*info > 0) {
#line 328 "stftri.f"
		    return 0;
#line 328 "stftri.f"
		}
#line 330 "stftri.f"
		strmm_("R", "U", "N", diag, &n1, &n2, &c_b18, &a[n1], n, a, n,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 333 "stftri.f"
	    }

#line 335 "stftri.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 339 "stftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1) */

#line 344 "stftri.f"
		strtri_("U", diag, &n1, a, &n1, info, (ftnlen)1, (ftnlen)1);
#line 345 "stftri.f"
		if (*info > 0) {
#line 345 "stftri.f"
		    return 0;
#line 345 "stftri.f"
		}
#line 347 "stftri.f"
		strmm_("L", "U", "N", diag, &n1, &n2, &c_b13, a, &n1, &a[n1 * 
			n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 349 "stftri.f"
		strtri_("L", diag, &n2, &a[1], &n1, info, (ftnlen)1, (ftnlen)
			1);
#line 350 "stftri.f"
		if (*info > 0) {
#line 350 "stftri.f"
		    *info += n1;
#line 350 "stftri.f"
		}
#line 352 "stftri.f"
		if (*info > 0) {
#line 352 "stftri.f"
		    return 0;
#line 352 "stftri.f"
		}
#line 354 "stftri.f"
		strmm_("R", "L", "T", diag, &n1, &n2, &c_b18, &a[1], &n1, &a[
			n1 * n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);

#line 357 "stftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0) */

#line 362 "stftri.f"
		strtri_("U", diag, &n1, &a[n2 * n2], &n2, info, (ftnlen)1, (
			ftnlen)1);
#line 363 "stftri.f"
		if (*info > 0) {
#line 363 "stftri.f"
		    return 0;
#line 363 "stftri.f"
		}
#line 365 "stftri.f"
		strmm_("R", "U", "T", diag, &n2, &n1, &c_b13, &a[n2 * n2], &
			n2, a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 367 "stftri.f"
		strtri_("L", diag, &n2, &a[n1 * n2], &n2, info, (ftnlen)1, (
			ftnlen)1);
#line 368 "stftri.f"
		if (*info > 0) {
#line 368 "stftri.f"
		    *info += n1;
#line 368 "stftri.f"
		}
#line 370 "stftri.f"
		if (*info > 0) {
#line 370 "stftri.f"
		    return 0;
#line 370 "stftri.f"
		}
#line 372 "stftri.f"
		strmm_("L", "L", "N", diag, &n2, &n1, &c_b18, &a[n1 * n2], &
			n2, a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 374 "stftri.f"
	    }

#line 376 "stftri.f"
	}

#line 378 "stftri.f"
    } else {

/*        N is even */

#line 382 "stftri.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 386 "stftri.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 392 "stftri.f"
		i__1 = *n + 1;
#line 392 "stftri.f"
		strtri_("L", diag, &k, &a[1], &i__1, info, (ftnlen)1, (ftnlen)
			1);
#line 393 "stftri.f"
		if (*info > 0) {
#line 393 "stftri.f"
		    return 0;
#line 393 "stftri.f"
		}
#line 395 "stftri.f"
		i__1 = *n + 1;
#line 395 "stftri.f"
		i__2 = *n + 1;
#line 395 "stftri.f"
		strmm_("R", "L", "N", diag, &k, &k, &c_b13, &a[1], &i__1, &a[
			k + 1], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 397 "stftri.f"
		i__1 = *n + 1;
#line 397 "stftri.f"
		strtri_("U", diag, &k, a, &i__1, info, (ftnlen)1, (ftnlen)1);
#line 398 "stftri.f"
		if (*info > 0) {
#line 398 "stftri.f"
		    *info += k;
#line 398 "stftri.f"
		}
#line 400 "stftri.f"
		if (*info > 0) {
#line 400 "stftri.f"
		    return 0;
#line 400 "stftri.f"
		}
#line 402 "stftri.f"
		i__1 = *n + 1;
#line 402 "stftri.f"
		i__2 = *n + 1;
#line 402 "stftri.f"
		strmm_("L", "U", "T", diag, &k, &k, &c_b18, a, &i__1, &a[k + 
			1], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1)
			;

#line 405 "stftri.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 411 "stftri.f"
		i__1 = *n + 1;
#line 411 "stftri.f"
		strtri_("L", diag, &k, &a[k + 1], &i__1, info, (ftnlen)1, (
			ftnlen)1);
#line 412 "stftri.f"
		if (*info > 0) {
#line 412 "stftri.f"
		    return 0;
#line 412 "stftri.f"
		}
#line 414 "stftri.f"
		i__1 = *n + 1;
#line 414 "stftri.f"
		i__2 = *n + 1;
#line 414 "stftri.f"
		strmm_("L", "L", "T", diag, &k, &k, &c_b13, &a[k + 1], &i__1, 
			a, &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 416 "stftri.f"
		i__1 = *n + 1;
#line 416 "stftri.f"
		strtri_("U", diag, &k, &a[k], &i__1, info, (ftnlen)1, (ftnlen)
			1);
#line 417 "stftri.f"
		if (*info > 0) {
#line 417 "stftri.f"
		    *info += k;
#line 417 "stftri.f"
		}
#line 419 "stftri.f"
		if (*info > 0) {
#line 419 "stftri.f"
		    return 0;
#line 419 "stftri.f"
		}
#line 421 "stftri.f"
		i__1 = *n + 1;
#line 421 "stftri.f"
		i__2 = *n + 1;
#line 421 "stftri.f"
		strmm_("R", "U", "N", diag, &k, &k, &c_b18, &a[k], &i__1, a, &
			i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 423 "stftri.f"
	    }
#line 424 "stftri.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 428 "stftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 434 "stftri.f"
		strtri_("U", diag, &k, &a[k], &k, info, (ftnlen)1, (ftnlen)1);
#line 435 "stftri.f"
		if (*info > 0) {
#line 435 "stftri.f"
		    return 0;
#line 435 "stftri.f"
		}
#line 437 "stftri.f"
		strmm_("L", "U", "N", diag, &k, &k, &c_b13, &a[k], &k, &a[k * 
			(k + 1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 439 "stftri.f"
		strtri_("L", diag, &k, a, &k, info, (ftnlen)1, (ftnlen)1);
#line 440 "stftri.f"
		if (*info > 0) {
#line 440 "stftri.f"
		    *info += k;
#line 440 "stftri.f"
		}
#line 442 "stftri.f"
		if (*info > 0) {
#line 442 "stftri.f"
		    return 0;
#line 442 "stftri.f"
		}
#line 444 "stftri.f"
		strmm_("R", "L", "T", diag, &k, &k, &c_b18, a, &k, &a[k * (k 
			+ 1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1)
			;
#line 446 "stftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 452 "stftri.f"
		strtri_("U", diag, &k, &a[k * (k + 1)], &k, info, (ftnlen)1, (
			ftnlen)1);
#line 453 "stftri.f"
		if (*info > 0) {
#line 453 "stftri.f"
		    return 0;
#line 453 "stftri.f"
		}
#line 455 "stftri.f"
		strmm_("R", "U", "T", diag, &k, &k, &c_b13, &a[k * (k + 1)], &
			k, a, &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 457 "stftri.f"
		strtri_("L", diag, &k, &a[k * k], &k, info, (ftnlen)1, (
			ftnlen)1);
#line 458 "stftri.f"
		if (*info > 0) {
#line 458 "stftri.f"
		    *info += k;
#line 458 "stftri.f"
		}
#line 460 "stftri.f"
		if (*info > 0) {
#line 460 "stftri.f"
		    return 0;
#line 460 "stftri.f"
		}
#line 462 "stftri.f"
		strmm_("L", "L", "N", diag, &k, &k, &c_b18, &a[k * k], &k, a, 
			&k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 464 "stftri.f"
	    }
#line 465 "stftri.f"
	}
#line 466 "stftri.f"
    }

#line 468 "stftri.f"
    return 0;

/*     End of STFTRI */

} /* stftri_ */

