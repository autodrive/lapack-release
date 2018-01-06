#line 1 "stfttr.f"
/* stfttr.f -- translated by f2c (version 20100827).
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

#line 1 "stfttr.f"
/* > \brief \b STFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard 
full format (TR). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STFTTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfttr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfttr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfttr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO ) */

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
/* > STFTTR copies a triangular matrix A from rectangular full packed */
/* > format (TF) to standard full format (TR). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF is in Normal format; */
/* >          = 'T':  ARF is in Transpose format. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices ARF and A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ARF */
/* > \verbatim */
/* >          ARF is REAL array, dimension (N*(N+1)/2). */
/* >          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') */
/* >          matrix A in RFP format. See the "Notes" below for more */
/* >          details. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On exit, the triangular matrix A.  If UPLO = 'U', the */
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
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
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

/*  ===================================================================== */
/* Subroutine */ int stfttr_(char *transr, char *uplo, integer *n, doublereal 
	*arf, doublereal *a, integer *lda, integer *info, ftnlen transr_len, 
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

#line 234 "stfttr.f"
    /* Parameter adjustments */
#line 234 "stfttr.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 234 "stfttr.f"
    a_offset = 0 + a_dim1 * 0;
#line 234 "stfttr.f"
    a -= a_offset;
#line 234 "stfttr.f"

#line 234 "stfttr.f"
    /* Function Body */
#line 234 "stfttr.f"
    *info = 0;
#line 235 "stfttr.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 236 "stfttr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 237 "stfttr.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 238 "stfttr.f"
	*info = -1;
#line 239 "stfttr.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 240 "stfttr.f"
	*info = -2;
#line 241 "stfttr.f"
    } else if (*n < 0) {
#line 242 "stfttr.f"
	*info = -3;
#line 243 "stfttr.f"
    } else if (*lda < max(1,*n)) {
#line 244 "stfttr.f"
	*info = -6;
#line 245 "stfttr.f"
    }
#line 246 "stfttr.f"
    if (*info != 0) {
#line 247 "stfttr.f"
	i__1 = -(*info);
#line 247 "stfttr.f"
	xerbla_("STFTTR", &i__1, (ftnlen)6);
#line 248 "stfttr.f"
	return 0;
#line 249 "stfttr.f"
    }

/*     Quick return if possible */

#line 253 "stfttr.f"
    if (*n <= 1) {
#line 254 "stfttr.f"
	if (*n == 1) {
#line 255 "stfttr.f"
	    a[0] = arf[0];
#line 256 "stfttr.f"
	}
#line 257 "stfttr.f"
	return 0;
#line 258 "stfttr.f"
    }

/*     Size of array ARF(0:nt-1) */

#line 262 "stfttr.f"
    nt = *n * (*n + 1) / 2;

/*     set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 266 "stfttr.f"
    if (lower) {
#line 267 "stfttr.f"
	n2 = *n / 2;
#line 268 "stfttr.f"
	n1 = *n - n2;
#line 269 "stfttr.f"
    } else {
#line 270 "stfttr.f"
	n1 = *n / 2;
#line 271 "stfttr.f"
	n2 = *n - n1;
#line 272 "stfttr.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 278 "stfttr.f"
    if (*n % 2 == 0) {
#line 279 "stfttr.f"
	k = *n / 2;
#line 280 "stfttr.f"
	nisodd = FALSE_;
#line 281 "stfttr.f"
	if (! lower) {
#line 281 "stfttr.f"
	    np1x2 = *n + *n + 2;
#line 281 "stfttr.f"
	}
#line 283 "stfttr.f"
    } else {
#line 284 "stfttr.f"
	nisodd = TRUE_;
#line 285 "stfttr.f"
	if (! lower) {
#line 285 "stfttr.f"
	    nx2 = *n + *n;
#line 285 "stfttr.f"
	}
#line 287 "stfttr.f"
    }

#line 289 "stfttr.f"
    if (nisodd) {

/*        N is odd */

#line 293 "stfttr.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 297 "stfttr.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 301 "stfttr.f"
		ij = 0;
#line 302 "stfttr.f"
		i__1 = n2;
#line 302 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 303 "stfttr.f"
		    i__2 = n2 + j;
#line 303 "stfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 304 "stfttr.f"
			a[n2 + j + i__ * a_dim1] = arf[ij];
#line 305 "stfttr.f"
			++ij;
#line 306 "stfttr.f"
		    }
#line 307 "stfttr.f"
		    i__2 = *n - 1;
#line 307 "stfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 308 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 309 "stfttr.f"
			++ij;
#line 310 "stfttr.f"
		    }
#line 311 "stfttr.f"
		}

#line 313 "stfttr.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 317 "stfttr.f"
		ij = nt - *n;
#line 318 "stfttr.f"
		i__1 = n1;
#line 318 "stfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 319 "stfttr.f"
		    i__2 = j;
#line 319 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 320 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 321 "stfttr.f"
			++ij;
#line 322 "stfttr.f"
		    }
#line 323 "stfttr.f"
		    i__2 = n1 - 1;
#line 323 "stfttr.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 324 "stfttr.f"
			a[j - n1 + l * a_dim1] = arf[ij];
#line 325 "stfttr.f"
			++ij;
#line 326 "stfttr.f"
		    }
#line 327 "stfttr.f"
		    ij -= nx2;
#line 328 "stfttr.f"
		}

#line 330 "stfttr.f"
	    }

#line 332 "stfttr.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 336 "stfttr.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 340 "stfttr.f"
		ij = 0;
#line 341 "stfttr.f"
		i__1 = n2 - 1;
#line 341 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 342 "stfttr.f"
		    i__2 = j;
#line 342 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 343 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 344 "stfttr.f"
			++ij;
#line 345 "stfttr.f"
		    }
#line 346 "stfttr.f"
		    i__2 = *n - 1;
#line 346 "stfttr.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 347 "stfttr.f"
			a[i__ + (n1 + j) * a_dim1] = arf[ij];
#line 348 "stfttr.f"
			++ij;
#line 349 "stfttr.f"
		    }
#line 350 "stfttr.f"
		}
#line 351 "stfttr.f"
		i__1 = *n - 1;
#line 351 "stfttr.f"
		for (j = n2; j <= i__1; ++j) {
#line 352 "stfttr.f"
		    i__2 = n1 - 1;
#line 352 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 353 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 354 "stfttr.f"
			++ij;
#line 355 "stfttr.f"
		    }
#line 356 "stfttr.f"
		}

#line 358 "stfttr.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 362 "stfttr.f"
		ij = 0;
#line 363 "stfttr.f"
		i__1 = n1;
#line 363 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 364 "stfttr.f"
		    i__2 = *n - 1;
#line 364 "stfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 365 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 366 "stfttr.f"
			++ij;
#line 367 "stfttr.f"
		    }
#line 368 "stfttr.f"
		}
#line 369 "stfttr.f"
		i__1 = n1 - 1;
#line 369 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 370 "stfttr.f"
		    i__2 = j;
#line 370 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 371 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 372 "stfttr.f"
			++ij;
#line 373 "stfttr.f"
		    }
#line 374 "stfttr.f"
		    i__2 = *n - 1;
#line 374 "stfttr.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 375 "stfttr.f"
			a[n2 + j + l * a_dim1] = arf[ij];
#line 376 "stfttr.f"
			++ij;
#line 377 "stfttr.f"
		    }
#line 378 "stfttr.f"
		}

#line 380 "stfttr.f"
	    }

#line 382 "stfttr.f"
	}

#line 384 "stfttr.f"
    } else {

/*        N is even */

#line 388 "stfttr.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 392 "stfttr.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 396 "stfttr.f"
		ij = 0;
#line 397 "stfttr.f"
		i__1 = k - 1;
#line 397 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 398 "stfttr.f"
		    i__2 = k + j;
#line 398 "stfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 399 "stfttr.f"
			a[k + j + i__ * a_dim1] = arf[ij];
#line 400 "stfttr.f"
			++ij;
#line 401 "stfttr.f"
		    }
#line 402 "stfttr.f"
		    i__2 = *n - 1;
#line 402 "stfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 403 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 404 "stfttr.f"
			++ij;
#line 405 "stfttr.f"
		    }
#line 406 "stfttr.f"
		}

#line 408 "stfttr.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 412 "stfttr.f"
		ij = nt - *n - 1;
#line 413 "stfttr.f"
		i__1 = k;
#line 413 "stfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 414 "stfttr.f"
		    i__2 = j;
#line 414 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 415 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 416 "stfttr.f"
			++ij;
#line 417 "stfttr.f"
		    }
#line 418 "stfttr.f"
		    i__2 = k - 1;
#line 418 "stfttr.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 419 "stfttr.f"
			a[j - k + l * a_dim1] = arf[ij];
#line 420 "stfttr.f"
			++ij;
#line 421 "stfttr.f"
		    }
#line 422 "stfttr.f"
		    ij -= np1x2;
#line 423 "stfttr.f"
		}

#line 425 "stfttr.f"
	    }

#line 427 "stfttr.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 431 "stfttr.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 435 "stfttr.f"
		ij = 0;
#line 436 "stfttr.f"
		j = k;
#line 437 "stfttr.f"
		i__1 = *n - 1;
#line 437 "stfttr.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 438 "stfttr.f"
		    a[i__ + j * a_dim1] = arf[ij];
#line 439 "stfttr.f"
		    ++ij;
#line 440 "stfttr.f"
		}
#line 441 "stfttr.f"
		i__1 = k - 2;
#line 441 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 442 "stfttr.f"
		    i__2 = j;
#line 442 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 443 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 444 "stfttr.f"
			++ij;
#line 445 "stfttr.f"
		    }
#line 446 "stfttr.f"
		    i__2 = *n - 1;
#line 446 "stfttr.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 447 "stfttr.f"
			a[i__ + (k + 1 + j) * a_dim1] = arf[ij];
#line 448 "stfttr.f"
			++ij;
#line 449 "stfttr.f"
		    }
#line 450 "stfttr.f"
		}
#line 451 "stfttr.f"
		i__1 = *n - 1;
#line 451 "stfttr.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 452 "stfttr.f"
		    i__2 = k - 1;
#line 452 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 453 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 454 "stfttr.f"
			++ij;
#line 455 "stfttr.f"
		    }
#line 456 "stfttr.f"
		}

#line 458 "stfttr.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 462 "stfttr.f"
		ij = 0;
#line 463 "stfttr.f"
		i__1 = k;
#line 463 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 464 "stfttr.f"
		    i__2 = *n - 1;
#line 464 "stfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 465 "stfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 466 "stfttr.f"
			++ij;
#line 467 "stfttr.f"
		    }
#line 468 "stfttr.f"
		}
#line 469 "stfttr.f"
		i__1 = k - 2;
#line 469 "stfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 470 "stfttr.f"
		    i__2 = j;
#line 470 "stfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 471 "stfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 472 "stfttr.f"
			++ij;
#line 473 "stfttr.f"
		    }
#line 474 "stfttr.f"
		    i__2 = *n - 1;
#line 474 "stfttr.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 475 "stfttr.f"
			a[k + 1 + j + l * a_dim1] = arf[ij];
#line 476 "stfttr.f"
			++ij;
#line 477 "stfttr.f"
		    }
#line 478 "stfttr.f"
		}
/*              Note that here, on exit of the loop, J = K-1 */
#line 480 "stfttr.f"
		i__1 = j;
#line 480 "stfttr.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 481 "stfttr.f"
		    a[i__ + j * a_dim1] = arf[ij];
#line 482 "stfttr.f"
		    ++ij;
#line 483 "stfttr.f"
		}

#line 485 "stfttr.f"
	    }

#line 487 "stfttr.f"
	}

#line 489 "stfttr.f"
    }

#line 491 "stfttr.f"
    return 0;

/*     End of STFTTR */

} /* stfttr_ */

