#line 1 "spftrf.f"
/* spftrf.f -- translated by f2c (version 20100827).
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

#line 1 "spftrf.f"
/* Table of constant values */

static doublereal c_b12 = 1.;
static doublereal c_b15 = -1.;

/* > \brief \b SPFTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPFTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spftrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spftrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spftrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPFTRF( TRANSR, UPLO, N, A, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            N, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( 0: * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPFTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the block version of the algorithm, calling Level 3 BLAS. */
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
/* >          = 'U':  Upper triangle of RFP A is stored; */
/* >          = 'L':  Lower triangle of RFP A is stored. */
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
/* >          A is REAL array, dimension ( N*(N+1)/2 ); */
/* >          On entry, the symmetric matrix A in RFP format. RFP format is */
/* >          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N' */
/* >          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is */
/* >          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is */
/* >          the transpose of RFP A as defined when */
/* >          TRANSR = 'N'. The contents of RFP A are defined by UPLO as */
/* >          follows: If UPLO = 'U' the RFP A contains the NT elements of */
/* >          upper packed A. If UPLO = 'L' the RFP A contains the elements */
/* >          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR = */
/* >          'T'. When TRANSR is 'N' the LDA is N+1 when N is even and N */
/* >          is odd. See the Note below for more details. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization RFP A = U**T*U or RFP A = L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the factorization could not be */
/* >                completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

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
/* Subroutine */ int spftrf_(char *transr, char *uplo, integer *n, doublereal 
	*a, integer *info, ftnlen transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, n1, n2;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), ssyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     xerbla_(char *, integer *, ftnlen);
    static logical nisodd;
    extern /* Subroutine */ int spotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

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

#line 237 "spftrf.f"
    *info = 0;
#line 238 "spftrf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 239 "spftrf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 240 "spftrf.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 241 "spftrf.f"
	*info = -1;
#line 242 "spftrf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 243 "spftrf.f"
	*info = -2;
#line 244 "spftrf.f"
    } else if (*n < 0) {
#line 245 "spftrf.f"
	*info = -3;
#line 246 "spftrf.f"
    }
#line 247 "spftrf.f"
    if (*info != 0) {
#line 248 "spftrf.f"
	i__1 = -(*info);
#line 248 "spftrf.f"
	xerbla_("SPFTRF", &i__1, (ftnlen)6);
#line 249 "spftrf.f"
	return 0;
#line 250 "spftrf.f"
    }

/*     Quick return if possible */

#line 254 "spftrf.f"
    if (*n == 0) {
#line 254 "spftrf.f"
	return 0;
#line 254 "spftrf.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

#line 260 "spftrf.f"
    if (*n % 2 == 0) {
#line 261 "spftrf.f"
	k = *n / 2;
#line 262 "spftrf.f"
	nisodd = FALSE_;
#line 263 "spftrf.f"
    } else {
#line 264 "spftrf.f"
	nisodd = TRUE_;
#line 265 "spftrf.f"
    }

/*     Set N1 and N2 depending on LOWER */

#line 269 "spftrf.f"
    if (lower) {
#line 270 "spftrf.f"
	n2 = *n / 2;
#line 271 "spftrf.f"
	n1 = *n - n2;
#line 272 "spftrf.f"
    } else {
#line 273 "spftrf.f"
	n1 = *n / 2;
#line 274 "spftrf.f"
	n2 = *n - n1;
#line 275 "spftrf.f"
    }

/*     start execution: there are eight cases */

#line 279 "spftrf.f"
    if (nisodd) {

/*        N is odd */

#line 283 "spftrf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 287 "spftrf.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1) */

#line 293 "spftrf.f"
		spotrf_("L", &n1, a, n, info, (ftnlen)1);
#line 294 "spftrf.f"
		if (*info > 0) {
#line 294 "spftrf.f"
		    return 0;
#line 294 "spftrf.f"
		}
#line 296 "spftrf.f"
		strsm_("R", "L", "T", "N", &n2, &n1, &c_b12, a, n, &a[n1], n, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 298 "spftrf.f"
		ssyrk_("U", "N", &n2, &n1, &c_b15, &a[n1], n, &c_b12, &a[*n], 
			n, (ftnlen)1, (ftnlen)1);
#line 300 "spftrf.f"
		spotrf_("U", &n2, &a[*n], n, info, (ftnlen)1);
#line 301 "spftrf.f"
		if (*info > 0) {
#line 301 "spftrf.f"
		    *info += n1;
#line 301 "spftrf.f"
		}

#line 304 "spftrf.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 310 "spftrf.f"
		spotrf_("L", &n1, &a[n2], n, info, (ftnlen)1);
#line 311 "spftrf.f"
		if (*info > 0) {
#line 311 "spftrf.f"
		    return 0;
#line 311 "spftrf.f"
		}
#line 313 "spftrf.f"
		strsm_("L", "L", "N", "N", &n1, &n2, &c_b12, &a[n2], n, a, n, 
			(ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 315 "spftrf.f"
		ssyrk_("U", "T", &n2, &n1, &c_b15, a, n, &c_b12, &a[n1], n, (
			ftnlen)1, (ftnlen)1);
#line 317 "spftrf.f"
		spotrf_("U", &n2, &a[n1], n, info, (ftnlen)1);
#line 318 "spftrf.f"
		if (*info > 0) {
#line 318 "spftrf.f"
		    *info += n1;
#line 318 "spftrf.f"
		}

#line 321 "spftrf.f"
	    }

#line 323 "spftrf.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 327 "spftrf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 333 "spftrf.f"
		spotrf_("U", &n1, a, &n1, info, (ftnlen)1);
#line 334 "spftrf.f"
		if (*info > 0) {
#line 334 "spftrf.f"
		    return 0;
#line 334 "spftrf.f"
		}
#line 336 "spftrf.f"
		strsm_("L", "U", "T", "N", &n1, &n2, &c_b12, a, &n1, &a[n1 * 
			n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 338 "spftrf.f"
		ssyrk_("L", "T", &n2, &n1, &c_b15, &a[n1 * n1], &n1, &c_b12, &
			a[1], &n1, (ftnlen)1, (ftnlen)1);
#line 340 "spftrf.f"
		spotrf_("L", &n2, &a[1], &n1, info, (ftnlen)1);
#line 341 "spftrf.f"
		if (*info > 0) {
#line 341 "spftrf.f"
		    *info += n1;
#line 341 "spftrf.f"
		}

#line 344 "spftrf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 350 "spftrf.f"
		spotrf_("U", &n1, &a[n2 * n2], &n2, info, (ftnlen)1);
#line 351 "spftrf.f"
		if (*info > 0) {
#line 351 "spftrf.f"
		    return 0;
#line 351 "spftrf.f"
		}
#line 353 "spftrf.f"
		strsm_("R", "U", "N", "N", &n2, &n1, &c_b12, &a[n2 * n2], &n2,
			 a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 355 "spftrf.f"
		ssyrk_("L", "N", &n2, &n1, &c_b15, a, &n2, &c_b12, &a[n1 * n2]
			, &n2, (ftnlen)1, (ftnlen)1);
#line 357 "spftrf.f"
		spotrf_("L", &n2, &a[n1 * n2], &n2, info, (ftnlen)1);
#line 358 "spftrf.f"
		if (*info > 0) {
#line 358 "spftrf.f"
		    *info += n1;
#line 358 "spftrf.f"
		}

#line 361 "spftrf.f"
	    }

#line 363 "spftrf.f"
	}

#line 365 "spftrf.f"
    } else {

/*        N is even */

#line 369 "spftrf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 373 "spftrf.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 379 "spftrf.f"
		i__1 = *n + 1;
#line 379 "spftrf.f"
		spotrf_("L", &k, &a[1], &i__1, info, (ftnlen)1);
#line 380 "spftrf.f"
		if (*info > 0) {
#line 380 "spftrf.f"
		    return 0;
#line 380 "spftrf.f"
		}
#line 382 "spftrf.f"
		i__1 = *n + 1;
#line 382 "spftrf.f"
		i__2 = *n + 1;
#line 382 "spftrf.f"
		strsm_("R", "L", "T", "N", &k, &k, &c_b12, &a[1], &i__1, &a[k 
			+ 1], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 384 "spftrf.f"
		i__1 = *n + 1;
#line 384 "spftrf.f"
		i__2 = *n + 1;
#line 384 "spftrf.f"
		ssyrk_("U", "N", &k, &k, &c_b15, &a[k + 1], &i__1, &c_b12, a, 
			&i__2, (ftnlen)1, (ftnlen)1);
#line 386 "spftrf.f"
		i__1 = *n + 1;
#line 386 "spftrf.f"
		spotrf_("U", &k, a, &i__1, info, (ftnlen)1);
#line 387 "spftrf.f"
		if (*info > 0) {
#line 387 "spftrf.f"
		    *info += k;
#line 387 "spftrf.f"
		}

#line 390 "spftrf.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 396 "spftrf.f"
		i__1 = *n + 1;
#line 396 "spftrf.f"
		spotrf_("L", &k, &a[k + 1], &i__1, info, (ftnlen)1);
#line 397 "spftrf.f"
		if (*info > 0) {
#line 397 "spftrf.f"
		    return 0;
#line 397 "spftrf.f"
		}
#line 399 "spftrf.f"
		i__1 = *n + 1;
#line 399 "spftrf.f"
		i__2 = *n + 1;
#line 399 "spftrf.f"
		strsm_("L", "L", "N", "N", &k, &k, &c_b12, &a[k + 1], &i__1, 
			a, &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 401 "spftrf.f"
		i__1 = *n + 1;
#line 401 "spftrf.f"
		i__2 = *n + 1;
#line 401 "spftrf.f"
		ssyrk_("U", "T", &k, &k, &c_b15, a, &i__1, &c_b12, &a[k], &
			i__2, (ftnlen)1, (ftnlen)1);
#line 403 "spftrf.f"
		i__1 = *n + 1;
#line 403 "spftrf.f"
		spotrf_("U", &k, &a[k], &i__1, info, (ftnlen)1);
#line 404 "spftrf.f"
		if (*info > 0) {
#line 404 "spftrf.f"
		    *info += k;
#line 404 "spftrf.f"
		}

#line 407 "spftrf.f"
	    }

#line 409 "spftrf.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 413 "spftrf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 419 "spftrf.f"
		spotrf_("U", &k, &a[k], &k, info, (ftnlen)1);
#line 420 "spftrf.f"
		if (*info > 0) {
#line 420 "spftrf.f"
		    return 0;
#line 420 "spftrf.f"
		}
#line 422 "spftrf.f"
		strsm_("L", "U", "T", "N", &k, &k, &c_b12, &a[k], &n1, &a[k * 
			(k + 1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 424 "spftrf.f"
		ssyrk_("L", "T", &k, &k, &c_b15, &a[k * (k + 1)], &k, &c_b12, 
			a, &k, (ftnlen)1, (ftnlen)1);
#line 426 "spftrf.f"
		spotrf_("L", &k, a, &k, info, (ftnlen)1);
#line 427 "spftrf.f"
		if (*info > 0) {
#line 427 "spftrf.f"
		    *info += k;
#line 427 "spftrf.f"
		}

#line 430 "spftrf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 436 "spftrf.f"
		spotrf_("U", &k, &a[k * (k + 1)], &k, info, (ftnlen)1);
#line 437 "spftrf.f"
		if (*info > 0) {
#line 437 "spftrf.f"
		    return 0;
#line 437 "spftrf.f"
		}
#line 439 "spftrf.f"
		strsm_("R", "U", "N", "N", &k, &k, &c_b12, &a[k * (k + 1)], &
			k, a, &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 441 "spftrf.f"
		ssyrk_("L", "N", &k, &k, &c_b15, a, &k, &c_b12, &a[k * k], &k,
			 (ftnlen)1, (ftnlen)1);
#line 443 "spftrf.f"
		spotrf_("L", &k, &a[k * k], &k, info, (ftnlen)1);
#line 444 "spftrf.f"
		if (*info > 0) {
#line 444 "spftrf.f"
		    *info += k;
#line 444 "spftrf.f"
		}

#line 447 "spftrf.f"
	    }

#line 449 "spftrf.f"
	}

#line 451 "spftrf.f"
    }

#line 453 "spftrf.f"
    return 0;

/*     End of SPFTRF */

} /* spftrf_ */

