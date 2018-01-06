#line 1 "cpftri.f"
/* cpftri.f -- translated by f2c (version 20100827).
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

#line 1 "cpftri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublereal c_b12 = 1.;

/* > \brief \b CPFTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPFTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpftri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpftri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpftri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPFTRI( TRANSR, UPLO, N, A, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. Array Arguments .. */
/*       COMPLEX            A( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPFTRI computes the inverse of a complex Hermitian positive definite */
/* > matrix A using the Cholesky factorization A = U**H*U or A = L*L**H */
/* > computed by CPFTRF. */
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
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
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
/* >          A is COMPLEX array, dimension ( N*(N+1)/2 ); */
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
/* >          On exit, the Hermitian inverse of the original matrix, in the */
/* >          same storage format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the (i,i) element of the factor U or L is */
/* >                zero, and the inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
/* Subroutine */ int cpftri_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *info, ftnlen transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, n1, n2;
    static logical normaltransr;
    extern /* Subroutine */ int cherk_(char *, char *, integer *, integer *, 
	    doublereal *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nisodd;
    extern /* Subroutine */ int clauum_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen), ctftri_(char *, char *, char *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
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

#line 252 "cpftri.f"
    *info = 0;
#line 253 "cpftri.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 254 "cpftri.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 255 "cpftri.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 256 "cpftri.f"
	*info = -1;
#line 257 "cpftri.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 258 "cpftri.f"
	*info = -2;
#line 259 "cpftri.f"
    } else if (*n < 0) {
#line 260 "cpftri.f"
	*info = -3;
#line 261 "cpftri.f"
    }
#line 262 "cpftri.f"
    if (*info != 0) {
#line 263 "cpftri.f"
	i__1 = -(*info);
#line 263 "cpftri.f"
	xerbla_("CPFTRI", &i__1, (ftnlen)6);
#line 264 "cpftri.f"
	return 0;
#line 265 "cpftri.f"
    }

/*     Quick return if possible */

#line 269 "cpftri.f"
    if (*n == 0) {
#line 269 "cpftri.f"
	return 0;
#line 269 "cpftri.f"
    }

/*     Invert the triangular Cholesky factor U or L. */

#line 274 "cpftri.f"
    ctftri_(transr, uplo, "N", n, a, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 275 "cpftri.f"
    if (*info > 0) {
#line 275 "cpftri.f"
	return 0;
#line 275 "cpftri.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

#line 281 "cpftri.f"
    if (*n % 2 == 0) {
#line 282 "cpftri.f"
	k = *n / 2;
#line 283 "cpftri.f"
	nisodd = FALSE_;
#line 284 "cpftri.f"
    } else {
#line 285 "cpftri.f"
	nisodd = TRUE_;
#line 286 "cpftri.f"
    }

/*     Set N1 and N2 depending on LOWER */

#line 290 "cpftri.f"
    if (lower) {
#line 291 "cpftri.f"
	n2 = *n / 2;
#line 292 "cpftri.f"
	n1 = *n - n2;
#line 293 "cpftri.f"
    } else {
#line 294 "cpftri.f"
	n1 = *n / 2;
#line 295 "cpftri.f"
	n2 = *n - n1;
#line 296 "cpftri.f"
    }

/*     Start execution of triangular matrix multiply: inv(U)*inv(U)^C or */
/*     inv(L)^C*inv(L). There are eight cases. */

#line 301 "cpftri.f"
    if (nisodd) {

/*        N is odd */

#line 305 "cpftri.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 309 "cpftri.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) ) */
/*              T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0) */
/*              T1 -> a(0), T2 -> a(n), S -> a(N1) */

#line 315 "cpftri.f"
		clauum_("L", &n1, a, n, info, (ftnlen)1);
#line 316 "cpftri.f"
		cherk_("L", "C", &n1, &n2, &c_b12, &a[n1], n, &c_b12, a, n, (
			ftnlen)1, (ftnlen)1);
#line 318 "cpftri.f"
		ctrmm_("L", "U", "N", "N", &n2, &n1, &c_b1, &a[*n], n, &a[n1],
			 n, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 320 "cpftri.f"
		clauum_("U", &n2, &a[*n], n, info, (ftnlen)1);

#line 322 "cpftri.f"
	    } else {

/*              SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1) */
/*              T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0) */
/*              T1 -> a(N2), T2 -> a(N1), S -> a(0) */

#line 328 "cpftri.f"
		clauum_("L", &n1, &a[n2], n, info, (ftnlen)1);
#line 329 "cpftri.f"
		cherk_("L", "N", &n1, &n2, &c_b12, a, n, &c_b12, &a[n2], n, (
			ftnlen)1, (ftnlen)1);
#line 331 "cpftri.f"
		ctrmm_("R", "U", "C", "N", &n1, &n2, &c_b1, &a[n1], n, a, n, (
			ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 333 "cpftri.f"
		clauum_("U", &n2, &a[n1], n, info, (ftnlen)1);

#line 335 "cpftri.f"
	    }

#line 337 "cpftri.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 341 "cpftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE, and N is odd */
/*              T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1) */

#line 346 "cpftri.f"
		clauum_("U", &n1, a, &n1, info, (ftnlen)1);
#line 347 "cpftri.f"
		cherk_("U", "N", &n1, &n2, &c_b12, &a[n1 * n1], &n1, &c_b12, 
			a, &n1, (ftnlen)1, (ftnlen)1);
#line 349 "cpftri.f"
		ctrmm_("R", "L", "N", "N", &n1, &n2, &c_b1, &a[1], &n1, &a[n1 
			* n1], &n1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 351 "cpftri.f"
		clauum_("L", &n2, &a[1], &n1, info, (ftnlen)1);

#line 353 "cpftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE, and N is odd */
/*              T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0) */

#line 358 "cpftri.f"
		clauum_("U", &n1, &a[n2 * n2], &n2, info, (ftnlen)1);
#line 359 "cpftri.f"
		cherk_("U", "C", &n1, &n2, &c_b12, a, &n2, &c_b12, &a[n2 * n2]
			, &n2, (ftnlen)1, (ftnlen)1);
#line 361 "cpftri.f"
		ctrmm_("L", "L", "C", "N", &n2, &n1, &c_b1, &a[n1 * n2], &n2, 
			a, &n2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 363 "cpftri.f"
		clauum_("L", &n2, &a[n1 * n2], &n2, info, (ftnlen)1);

#line 365 "cpftri.f"
	    }

#line 367 "cpftri.f"
	}

#line 369 "cpftri.f"
    } else {

/*        N is even */

#line 373 "cpftri.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 377 "cpftri.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 383 "cpftri.f"
		i__1 = *n + 1;
#line 383 "cpftri.f"
		clauum_("L", &k, &a[1], &i__1, info, (ftnlen)1);
#line 384 "cpftri.f"
		i__1 = *n + 1;
#line 384 "cpftri.f"
		i__2 = *n + 1;
#line 384 "cpftri.f"
		cherk_("L", "C", &k, &k, &c_b12, &a[k + 1], &i__1, &c_b12, &a[
			1], &i__2, (ftnlen)1, (ftnlen)1);
#line 386 "cpftri.f"
		i__1 = *n + 1;
#line 386 "cpftri.f"
		i__2 = *n + 1;
#line 386 "cpftri.f"
		ctrmm_("L", "U", "N", "N", &k, &k, &c_b1, a, &i__1, &a[k + 1],
			 &i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 388 "cpftri.f"
		i__1 = *n + 1;
#line 388 "cpftri.f"
		clauum_("U", &k, a, &i__1, info, (ftnlen)1);

#line 390 "cpftri.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 396 "cpftri.f"
		i__1 = *n + 1;
#line 396 "cpftri.f"
		clauum_("L", &k, &a[k + 1], &i__1, info, (ftnlen)1);
#line 397 "cpftri.f"
		i__1 = *n + 1;
#line 397 "cpftri.f"
		i__2 = *n + 1;
#line 397 "cpftri.f"
		cherk_("L", "N", &k, &k, &c_b12, a, &i__1, &c_b12, &a[k + 1], 
			&i__2, (ftnlen)1, (ftnlen)1);
#line 399 "cpftri.f"
		i__1 = *n + 1;
#line 399 "cpftri.f"
		i__2 = *n + 1;
#line 399 "cpftri.f"
		ctrmm_("R", "U", "C", "N", &k, &k, &c_b1, &a[k], &i__1, a, &
			i__2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 401 "cpftri.f"
		i__1 = *n + 1;
#line 401 "cpftri.f"
		clauum_("U", &k, &a[k], &i__1, info, (ftnlen)1);

#line 403 "cpftri.f"
	    }

#line 405 "cpftri.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 409 "cpftri.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE, and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1), */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 415 "cpftri.f"
		clauum_("U", &k, &a[k], &k, info, (ftnlen)1);
#line 416 "cpftri.f"
		cherk_("U", "N", &k, &k, &c_b12, &a[k * (k + 1)], &k, &c_b12, 
			&a[k], &k, (ftnlen)1, (ftnlen)1);
#line 418 "cpftri.f"
		ctrmm_("R", "L", "N", "N", &k, &k, &c_b1, a, &k, &a[k * (k + 
			1)], &k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 420 "cpftri.f"
		clauum_("L", &k, a, &k, info, (ftnlen)1);

#line 422 "cpftri.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE, and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0), */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 428 "cpftri.f"
		clauum_("U", &k, &a[k * (k + 1)], &k, info, (ftnlen)1);
#line 429 "cpftri.f"
		cherk_("U", "C", &k, &k, &c_b12, a, &k, &c_b12, &a[k * (k + 1)
			], &k, (ftnlen)1, (ftnlen)1);
#line 431 "cpftri.f"
		ctrmm_("L", "L", "C", "N", &k, &k, &c_b1, &a[k * k], &k, a, &
			k, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 433 "cpftri.f"
		clauum_("L", &k, &a[k * k], &k, info, (ftnlen)1);

#line 435 "cpftri.f"
	    }

#line 437 "cpftri.f"
	}

#line 439 "cpftri.f"
    }

#line 441 "cpftri.f"
    return 0;

/*     End of CPFTRI */

} /* cpftri_ */

