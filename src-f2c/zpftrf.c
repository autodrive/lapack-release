#line 1 "zpftrf.f"
/* zpftrf.f -- translated by f2c (version 20100827).
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

#line 1 "zpftrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublereal c_b15 = -1.;
static doublereal c_b16 = 1.;

/* > \brief \b ZPFTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPFTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpftrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpftrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpftrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPFTRF( TRANSR, UPLO, N, A, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            N, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( 0: * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPFTRF computes the Cholesky factorization of a complex Hermitian */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**H * U,  if UPLO = 'U', or */
/* >    A = L  * L**H,  if UPLO = 'L', */
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
/* >          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored. */
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
/* >          A is COMPLEX*16 array, dimension ( N*(N+1)/2 ); */
/* >          On entry, the Hermitian matrix A in RFP format. RFP format is */
/* >          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N' */
/* >          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is */
/* >          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'C' then RFP is */
/* >          the Conjugate-transpose of RFP A as defined when */
/* >          TRANSR = 'N'. The contents of RFP A are defined by UPLO as */
/* >          follows: If UPLO = 'U' the RFP A contains the nt elements of */
/* >          upper packed A. If UPLO = 'L' the RFP A contains the elements */
/* >          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR = */
/* >          'C'. When TRANSR is 'N' the LDA is N+1 when N is even and N */
/* >          is odd. See the Note below for more details. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization RFP A = U**H*U or RFP A = L*L**H. */
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
/* > */
/* >  Further Notes on RFP Format: */
/* >  ============================ */
/* > */
/* >  We first consider Standard Packed Format when N is even. */
/* >  We give an example where N = 6. */
/* > */
/* >     AP is Upper             AP is Lower */
/* > */
/* >   00 01 02 03 04 05       00 */
/* >      11 12 13 14 15       10 11 */
/* >         22 23 24 25       20 21 22 */
/* >            33 34 35       30 31 32 33 */
/* >               44 45       40 41 42 43 44 */
/* >                  55       50 51 52 53 54 55 */
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
/* >           RFP A                   RFP A */
/* > */
/* >     -- -- -- --                -- -- -- -- -- -- */
/* >     03 13 23 33 00 01 02    33 00 10 20 30 40 50 */
/* >     -- -- -- -- --                -- -- -- -- -- */
/* >     04 14 24 34 44 11 12    43 44 11 21 31 41 51 */
/* >     -- -- -- -- -- --                -- -- -- -- */
/* >     05 15 25 35 45 55 22    53 54 55 22 32 42 52 */
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
/* >           RFP A                   RFP A */
/* > */
/* >     -- -- --                   -- -- -- -- -- -- */
/* >     02 12 22 00 01             00 10 20 30 40 50 */
/* >     -- -- -- --                   -- -- -- -- -- */
/* >     03 13 23 33 11             33 11 21 31 41 51 */
/* >     -- -- -- -- --                   -- -- -- -- */
/* >     04 14 24 34 44             43 44 22 32 42 52 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpftrf_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *info, ftnlen transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, n1, n2;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zherk_(char *, char *, integer *, integer *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical nisodd;
    extern /* Subroutine */ int zpotrf_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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

#line 251 "zpftrf.f"
    *info = 0;
#line 252 "zpftrf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 253 "zpftrf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 254 "zpftrf.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 255 "zpftrf.f"
	*info = -1;
#line 256 "zpftrf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 257 "zpftrf.f"
	*info = -2;
#line 258 "zpftrf.f"
    } else if (*n < 0) {
#line 259 "zpftrf.f"
	*info = -3;
#line 260 "zpftrf.f"
    }
#line 261 "zpftrf.f"
    if (*info != 0) {
#line 262 "zpftrf.f"
	i__1 = -(*info);
#line 262 "zpftrf.f"
	xerbla_("ZPFTRF", &i__1, (ftnlen)6);
#line 263 "zpftrf.f"
	return 0;
#line 264 "zpftrf.f"
    }

/*     Quick return if possible */

#line 268 "zpftrf.f"
    if (*n == 0) {
#line 268 "zpftrf.f"
	return 0;
#line 268 "zpftrf.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

#line 274 "zpftrf.f"
    if (*n % 2 == 0) {
#line 275 "zpftrf.f"
	k = *n / 2;
#line 276 "zpftrf.f"
	nisodd = FALSE_;
#line 277 "zpftrf.f"
    } else {
#line 278 "zpftrf.f"
	nisodd = TRUE_;
#line 279 "zpftrf.f"
    }

/*     Set N1 and N2 depending on LOWER */

#line 283 "zpftrf.f"
    if (lower) {
#line 284 "zpftrf.f"
	n2 = *n / 2;
#line 285 "zpftrf.f"
	n1 = *n - n2;
#line 286 "zpftrf.f"
    } else {
#line 287 "zpftrf.f"
	n1 = *n / 2;
#line 288 "zpftrf.f"
	n2 = *n - n1;
#line 289 "zpftrf.f"
    }

/*     start execution: there are eight cases */

#line 293 "zpftrf.f"
    if (nisodd) {

/*        N is odd */

#line 297 "zpftrf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 301 "zpftrf.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1) */

#line 307 "zpftrf.f"
		zpotrf_("L", &n1, a, n, info, (ftnlen)1);
#line 308 "zpftrf.f"
		if (*info > 0) {
#line 308 "zpftrf.f"
		    return 0;
#line 308 "zpftrf.f"
		}
#line 310 "zpftrf.f"
		ztrsm_("R", "L", "C", "N", &n2, &n1, &c_b1, a, n, &a[n1], n, (
			ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 312 "zpftrf.f"
		zherk_("U", "N", &n2, &n1, &c_b15, &a[n1], n, &c_b16, &a[*n], 
			n, (ftnlen)1, (ftnlen)1);
#line 314 "zpftrf.f"
		zpotrf_("U", &n2, &a[*n], n, info, (ftnlen)1);
#line 315 "zpftrf.f"
		if (*info > 0) {
#line 315 "zpftrf.f"
		    *info += n1;
#line 315 "zpftrf.f"
		}

#line 318 "zpftrf.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 324 "zpftrf.f"
		zpotrf_("L", &n1, &a[n2], n, info, (ftnlen)1);
#line 325 "zpftrf.f"
		if (*info > 0) {
#line 325 "zpftrf.f"
		    return 0;
#line 325 "zpftrf.f"
		}
#line 327 "zpftrf.f"
		ztrsm_("L", "L", "N", "N", &n1, &n2, &c_b1, &a[n2], n, a, n, (
			ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 329 "zpftrf.f"
		zherk_("U", "C", &n2, &n1, &c_b15, a, n, &c_b16, &a[n1], n, (
			ftnlen)1, (ftnlen)1);
#line 331 "zpftrf.f"
		zpotrf_("U", &n2, &a[n1], n, info, (ftnlen)1);
#line 332 "zpftrf.f"
		if (*info > 0) {
#line 332 "zpftrf.f"
		    *info += n1;
#line 332 "zpftrf.f"
		}

#line 335 "zpftrf.f"
	    }

#line 337 "zpftrf.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 341 "zpftrf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 347 "zpftrf.f"
		zpotrf_("U", &n1, a, &n1, info, (ftnlen)1);
#line 348 "zpftrf.f"
		if (*info > 0) {
#line 348 "zpftrf.f"
		    return 0;
#line 348 "zpftrf.f"
		}
#line 350 "zpftrf.f"
		ztrsm_("L", "U", "C", "N", &n1, &n2, &c_b1, a, &n1, &a[n1 * 
			n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 352 "zpftrf.f"
		zherk_("L", "C", &n2, &n1, &c_b15, &a[n1 * n1], &n1, &c_b16, &
			a[1], &n1, (ftnlen)1, (ftnlen)1);
#line 354 "zpftrf.f"
		zpotrf_("L", &n2, &a[1], &n1, info, (ftnlen)1);
#line 355 "zpftrf.f"
		if (*info > 0) {
#line 355 "zpftrf.f"
		    *info += n1;
#line 355 "zpftrf.f"
		}

#line 358 "zpftrf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 364 "zpftrf.f"
		zpotrf_("U", &n1, &a[n2 * n2], &n2, info, (ftnlen)1);
#line 365 "zpftrf.f"
		if (*info > 0) {
#line 365 "zpftrf.f"
		    return 0;
#line 365 "zpftrf.f"
		}
#line 367 "zpftrf.f"
		ztrsm_("R", "U", "N", "N", &n2, &n1, &c_b1, &a[n2 * n2], &n2, 
			a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 369 "zpftrf.f"
		zherk_("L", "N", &n2, &n1, &c_b15, a, &n2, &c_b16, &a[n1 * n2]
			, &n2, (ftnlen)1, (ftnlen)1);
#line 371 "zpftrf.f"
		zpotrf_("L", &n2, &a[n1 * n2], &n2, info, (ftnlen)1);
#line 372 "zpftrf.f"
		if (*info > 0) {
#line 372 "zpftrf.f"
		    *info += n1;
#line 372 "zpftrf.f"
		}

#line 375 "zpftrf.f"
	    }

#line 377 "zpftrf.f"
	}

#line 379 "zpftrf.f"
    } else {

/*        N is even */

#line 383 "zpftrf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 387 "zpftrf.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 393 "zpftrf.f"
		i__1 = *n + 1;
#line 393 "zpftrf.f"
		zpotrf_("L", &k, &a[1], &i__1, info, (ftnlen)1);
#line 394 "zpftrf.f"
		if (*info > 0) {
#line 394 "zpftrf.f"
		    return 0;
#line 394 "zpftrf.f"
		}
#line 396 "zpftrf.f"
		i__1 = *n + 1;
#line 396 "zpftrf.f"
		i__2 = *n + 1;
#line 396 "zpftrf.f"
		ztrsm_("R", "L", "C", "N", &k, &k, &c_b1, &a[1], &i__1, &a[k 
			+ 1], &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 398 "zpftrf.f"
		i__1 = *n + 1;
#line 398 "zpftrf.f"
		i__2 = *n + 1;
#line 398 "zpftrf.f"
		zherk_("U", "N", &k, &k, &c_b15, &a[k + 1], &i__1, &c_b16, a, 
			&i__2, (ftnlen)1, (ftnlen)1);
#line 400 "zpftrf.f"
		i__1 = *n + 1;
#line 400 "zpftrf.f"
		zpotrf_("U", &k, a, &i__1, info, (ftnlen)1);
#line 401 "zpftrf.f"
		if (*info > 0) {
#line 401 "zpftrf.f"
		    *info += k;
#line 401 "zpftrf.f"
		}

#line 404 "zpftrf.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 410 "zpftrf.f"
		i__1 = *n + 1;
#line 410 "zpftrf.f"
		zpotrf_("L", &k, &a[k + 1], &i__1, info, (ftnlen)1);
#line 411 "zpftrf.f"
		if (*info > 0) {
#line 411 "zpftrf.f"
		    return 0;
#line 411 "zpftrf.f"
		}
#line 413 "zpftrf.f"
		i__1 = *n + 1;
#line 413 "zpftrf.f"
		i__2 = *n + 1;
#line 413 "zpftrf.f"
		ztrsm_("L", "L", "N", "N", &k, &k, &c_b1, &a[k + 1], &i__1, a,
			 &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 415 "zpftrf.f"
		i__1 = *n + 1;
#line 415 "zpftrf.f"
		i__2 = *n + 1;
#line 415 "zpftrf.f"
		zherk_("U", "C", &k, &k, &c_b15, a, &i__1, &c_b16, &a[k], &
			i__2, (ftnlen)1, (ftnlen)1);
#line 417 "zpftrf.f"
		i__1 = *n + 1;
#line 417 "zpftrf.f"
		zpotrf_("U", &k, &a[k], &i__1, info, (ftnlen)1);
#line 418 "zpftrf.f"
		if (*info > 0) {
#line 418 "zpftrf.f"
		    *info += k;
#line 418 "zpftrf.f"
		}

#line 421 "zpftrf.f"
	    }

#line 423 "zpftrf.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 427 "zpftrf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 433 "zpftrf.f"
		zpotrf_("U", &k, &a[k], &k, info, (ftnlen)1);
#line 434 "zpftrf.f"
		if (*info > 0) {
#line 434 "zpftrf.f"
		    return 0;
#line 434 "zpftrf.f"
		}
#line 436 "zpftrf.f"
		ztrsm_("L", "U", "C", "N", &k, &k, &c_b1, &a[k], &n1, &a[k * (
			k + 1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 438 "zpftrf.f"
		zherk_("L", "C", &k, &k, &c_b15, &a[k * (k + 1)], &k, &c_b16, 
			a, &k, (ftnlen)1, (ftnlen)1);
#line 440 "zpftrf.f"
		zpotrf_("L", &k, a, &k, info, (ftnlen)1);
#line 441 "zpftrf.f"
		if (*info > 0) {
#line 441 "zpftrf.f"
		    *info += k;
#line 441 "zpftrf.f"
		}

#line 444 "zpftrf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 450 "zpftrf.f"
		zpotrf_("U", &k, &a[k * (k + 1)], &k, info, (ftnlen)1);
#line 451 "zpftrf.f"
		if (*info > 0) {
#line 451 "zpftrf.f"
		    return 0;
#line 451 "zpftrf.f"
		}
#line 453 "zpftrf.f"
		ztrsm_("R", "U", "N", "N", &k, &k, &c_b1, &a[k * (k + 1)], &k,
			 a, &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 455 "zpftrf.f"
		zherk_("L", "N", &k, &k, &c_b15, a, &k, &c_b16, &a[k * k], &k,
			 (ftnlen)1, (ftnlen)1);
#line 457 "zpftrf.f"
		zpotrf_("L", &k, &a[k * k], &k, info, (ftnlen)1);
#line 458 "zpftrf.f"
		if (*info > 0) {
#line 458 "zpftrf.f"
		    *info += k;
#line 458 "zpftrf.f"
		}

#line 461 "zpftrf.f"
	    }

#line 463 "zpftrf.f"
	}

#line 465 "zpftrf.f"
    }

#line 467 "zpftrf.f"
    return 0;

/*     End of ZPFTRF */

} /* zpftrf_ */

