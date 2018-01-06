#line 1 "strttf.f"
/* strttf.f -- translated by f2c (version 20100827).
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

#line 1 "strttf.f"
/* > \brief \b STRTTF copies a triangular matrix from the standard full format (TR) to the rectangular full pa
cked format (TF). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRTTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strttf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strttf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strttf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( 0: LDA-1, 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRTTF copies a triangular matrix A from standard full format (TR) */
/* > to rectangular full packed format (TF) . */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF in Normal form is wanted; */
/* >          = 'T':  ARF in Transpose form is wanted. */
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
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N). */
/* >          On entry, the triangular matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of the array A contains */
/* >          the upper triangular matrix, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of the array A contains */
/* >          the lower triangular matrix, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the matrix A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ARF */
/* > \verbatim */
/* >          ARF is REAL array, dimension (NT). */
/* >          NT=N*(N+1)/2. On exit, the triangular matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

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
/* Subroutine */ int strttf_(char *transr, char *uplo, integer *n, doublereal 
	*a, integer *lda, doublereal *arf, integer *info, ftnlen transr_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, n1, n2, ij, nt, nx2, np1x2;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nisodd;


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

#line 231 "strttf.f"
    /* Parameter adjustments */
#line 231 "strttf.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 231 "strttf.f"
    a_offset = 0 + a_dim1 * 0;
#line 231 "strttf.f"
    a -= a_offset;
#line 231 "strttf.f"

#line 231 "strttf.f"
    /* Function Body */
#line 231 "strttf.f"
    *info = 0;
#line 232 "strttf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 233 "strttf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 234 "strttf.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 235 "strttf.f"
	*info = -1;
#line 236 "strttf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 237 "strttf.f"
	*info = -2;
#line 238 "strttf.f"
    } else if (*n < 0) {
#line 239 "strttf.f"
	*info = -3;
#line 240 "strttf.f"
    } else if (*lda < max(1,*n)) {
#line 241 "strttf.f"
	*info = -5;
#line 242 "strttf.f"
    }
#line 243 "strttf.f"
    if (*info != 0) {
#line 244 "strttf.f"
	i__1 = -(*info);
#line 244 "strttf.f"
	xerbla_("STRTTF", &i__1, (ftnlen)6);
#line 245 "strttf.f"
	return 0;
#line 246 "strttf.f"
    }

/*     Quick return if possible */

#line 250 "strttf.f"
    if (*n <= 1) {
#line 251 "strttf.f"
	if (*n == 1) {
#line 252 "strttf.f"
	    arf[0] = a[0];
#line 253 "strttf.f"
	}
#line 254 "strttf.f"
	return 0;
#line 255 "strttf.f"
    }

/*     Size of array ARF(0:nt-1) */

#line 259 "strttf.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 263 "strttf.f"
    if (lower) {
#line 264 "strttf.f"
	n2 = *n / 2;
#line 265 "strttf.f"
	n1 = *n - n2;
#line 266 "strttf.f"
    } else {
#line 267 "strttf.f"
	n1 = *n / 2;
#line 268 "strttf.f"
	n2 = *n - n1;
#line 269 "strttf.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 275 "strttf.f"
    if (*n % 2 == 0) {
#line 276 "strttf.f"
	k = *n / 2;
#line 277 "strttf.f"
	nisodd = FALSE_;
#line 278 "strttf.f"
	if (! lower) {
#line 278 "strttf.f"
	    np1x2 = *n + *n + 2;
#line 278 "strttf.f"
	}
#line 280 "strttf.f"
    } else {
#line 281 "strttf.f"
	nisodd = TRUE_;
#line 282 "strttf.f"
	if (! lower) {
#line 282 "strttf.f"
	    nx2 = *n + *n;
#line 282 "strttf.f"
	}
#line 284 "strttf.f"
    }

#line 286 "strttf.f"
    if (nisodd) {

/*        N is odd */

#line 290 "strttf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 294 "strttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 298 "strttf.f"
		ij = 0;
#line 299 "strttf.f"
		i__1 = n2;
#line 299 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 300 "strttf.f"
		    i__2 = n2 + j;
#line 300 "strttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 301 "strttf.f"
			arf[ij] = a[n2 + j + i__ * a_dim1];
#line 302 "strttf.f"
			++ij;
#line 303 "strttf.f"
		    }
#line 304 "strttf.f"
		    i__2 = *n - 1;
#line 304 "strttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 305 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 306 "strttf.f"
			++ij;
#line 307 "strttf.f"
		    }
#line 308 "strttf.f"
		}

#line 310 "strttf.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 314 "strttf.f"
		ij = nt - *n;
#line 315 "strttf.f"
		i__1 = n1;
#line 315 "strttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 316 "strttf.f"
		    i__2 = j;
#line 316 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 317 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 318 "strttf.f"
			++ij;
#line 319 "strttf.f"
		    }
#line 320 "strttf.f"
		    i__2 = n1 - 1;
#line 320 "strttf.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 321 "strttf.f"
			arf[ij] = a[j - n1 + l * a_dim1];
#line 322 "strttf.f"
			++ij;
#line 323 "strttf.f"
		    }
#line 324 "strttf.f"
		    ij -= nx2;
#line 325 "strttf.f"
		}

#line 327 "strttf.f"
	    }

#line 329 "strttf.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 333 "strttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 337 "strttf.f"
		ij = 0;
#line 338 "strttf.f"
		i__1 = n2 - 1;
#line 338 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 339 "strttf.f"
		    i__2 = j;
#line 339 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 340 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 341 "strttf.f"
			++ij;
#line 342 "strttf.f"
		    }
#line 343 "strttf.f"
		    i__2 = *n - 1;
#line 343 "strttf.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 344 "strttf.f"
			arf[ij] = a[i__ + (n1 + j) * a_dim1];
#line 345 "strttf.f"
			++ij;
#line 346 "strttf.f"
		    }
#line 347 "strttf.f"
		}
#line 348 "strttf.f"
		i__1 = *n - 1;
#line 348 "strttf.f"
		for (j = n2; j <= i__1; ++j) {
#line 349 "strttf.f"
		    i__2 = n1 - 1;
#line 349 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 350 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 351 "strttf.f"
			++ij;
#line 352 "strttf.f"
		    }
#line 353 "strttf.f"
		}

#line 355 "strttf.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 359 "strttf.f"
		ij = 0;
#line 360 "strttf.f"
		i__1 = n1;
#line 360 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 361 "strttf.f"
		    i__2 = *n - 1;
#line 361 "strttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 362 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 363 "strttf.f"
			++ij;
#line 364 "strttf.f"
		    }
#line 365 "strttf.f"
		}
#line 366 "strttf.f"
		i__1 = n1 - 1;
#line 366 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 367 "strttf.f"
		    i__2 = j;
#line 367 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 368 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 369 "strttf.f"
			++ij;
#line 370 "strttf.f"
		    }
#line 371 "strttf.f"
		    i__2 = *n - 1;
#line 371 "strttf.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 372 "strttf.f"
			arf[ij] = a[n2 + j + l * a_dim1];
#line 373 "strttf.f"
			++ij;
#line 374 "strttf.f"
		    }
#line 375 "strttf.f"
		}

#line 377 "strttf.f"
	    }

#line 379 "strttf.f"
	}

#line 381 "strttf.f"
    } else {

/*        N is even */

#line 385 "strttf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 389 "strttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 393 "strttf.f"
		ij = 0;
#line 394 "strttf.f"
		i__1 = k - 1;
#line 394 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 395 "strttf.f"
		    i__2 = k + j;
#line 395 "strttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 396 "strttf.f"
			arf[ij] = a[k + j + i__ * a_dim1];
#line 397 "strttf.f"
			++ij;
#line 398 "strttf.f"
		    }
#line 399 "strttf.f"
		    i__2 = *n - 1;
#line 399 "strttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 400 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 401 "strttf.f"
			++ij;
#line 402 "strttf.f"
		    }
#line 403 "strttf.f"
		}

#line 405 "strttf.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 409 "strttf.f"
		ij = nt - *n - 1;
#line 410 "strttf.f"
		i__1 = k;
#line 410 "strttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 411 "strttf.f"
		    i__2 = j;
#line 411 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 412 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 413 "strttf.f"
			++ij;
#line 414 "strttf.f"
		    }
#line 415 "strttf.f"
		    i__2 = k - 1;
#line 415 "strttf.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 416 "strttf.f"
			arf[ij] = a[j - k + l * a_dim1];
#line 417 "strttf.f"
			++ij;
#line 418 "strttf.f"
		    }
#line 419 "strttf.f"
		    ij -= np1x2;
#line 420 "strttf.f"
		}

#line 422 "strttf.f"
	    }

#line 424 "strttf.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 428 "strttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 432 "strttf.f"
		ij = 0;
#line 433 "strttf.f"
		j = k;
#line 434 "strttf.f"
		i__1 = *n - 1;
#line 434 "strttf.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 435 "strttf.f"
		    arf[ij] = a[i__ + j * a_dim1];
#line 436 "strttf.f"
		    ++ij;
#line 437 "strttf.f"
		}
#line 438 "strttf.f"
		i__1 = k - 2;
#line 438 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 439 "strttf.f"
		    i__2 = j;
#line 439 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 440 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 441 "strttf.f"
			++ij;
#line 442 "strttf.f"
		    }
#line 443 "strttf.f"
		    i__2 = *n - 1;
#line 443 "strttf.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 444 "strttf.f"
			arf[ij] = a[i__ + (k + 1 + j) * a_dim1];
#line 445 "strttf.f"
			++ij;
#line 446 "strttf.f"
		    }
#line 447 "strttf.f"
		}
#line 448 "strttf.f"
		i__1 = *n - 1;
#line 448 "strttf.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 449 "strttf.f"
		    i__2 = k - 1;
#line 449 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 450 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 451 "strttf.f"
			++ij;
#line 452 "strttf.f"
		    }
#line 453 "strttf.f"
		}

#line 455 "strttf.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 459 "strttf.f"
		ij = 0;
#line 460 "strttf.f"
		i__1 = k;
#line 460 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 461 "strttf.f"
		    i__2 = *n - 1;
#line 461 "strttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 462 "strttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 463 "strttf.f"
			++ij;
#line 464 "strttf.f"
		    }
#line 465 "strttf.f"
		}
#line 466 "strttf.f"
		i__1 = k - 2;
#line 466 "strttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 467 "strttf.f"
		    i__2 = j;
#line 467 "strttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 468 "strttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 469 "strttf.f"
			++ij;
#line 470 "strttf.f"
		    }
#line 471 "strttf.f"
		    i__2 = *n - 1;
#line 471 "strttf.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 472 "strttf.f"
			arf[ij] = a[k + 1 + j + l * a_dim1];
#line 473 "strttf.f"
			++ij;
#line 474 "strttf.f"
		    }
#line 475 "strttf.f"
		}
/*              Note that here, on exit of the loop, J = K-1 */
#line 477 "strttf.f"
		i__1 = j;
#line 477 "strttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 478 "strttf.f"
		    arf[ij] = a[i__ + j * a_dim1];
#line 479 "strttf.f"
		    ++ij;
#line 480 "strttf.f"
		}

#line 482 "strttf.f"
	    }

#line 484 "strttf.f"
	}

#line 486 "strttf.f"
    }

#line 488 "strttf.f"
    return 0;

/*     End of STRTTF */

} /* strttf_ */

