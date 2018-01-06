#line 1 "dtfttr.f"
/* dtfttr.f -- translated by f2c (version 20100827).
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

#line 1 "dtfttr.f"
/* > \brief \b DTFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard 
full format (TR). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTFTTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtfttr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtfttr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtfttr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( 0: LDA-1, 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTFTTR copies a triangular matrix A from rectangular full packed */
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
/* >          ARF is DOUBLE PRECISION array, dimension (N*(N+1)/2). */
/* >          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') */
/* >          matrix A in RFP format. See the "Notes" below for more */
/* >          details. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* Subroutine */ int dtfttr_(char *transr, char *uplo, integer *n, doublereal 
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

#line 234 "dtfttr.f"
    /* Parameter adjustments */
#line 234 "dtfttr.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 234 "dtfttr.f"
    a_offset = 0 + a_dim1 * 0;
#line 234 "dtfttr.f"
    a -= a_offset;
#line 234 "dtfttr.f"

#line 234 "dtfttr.f"
    /* Function Body */
#line 234 "dtfttr.f"
    *info = 0;
#line 235 "dtfttr.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 236 "dtfttr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 237 "dtfttr.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 238 "dtfttr.f"
	*info = -1;
#line 239 "dtfttr.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 240 "dtfttr.f"
	*info = -2;
#line 241 "dtfttr.f"
    } else if (*n < 0) {
#line 242 "dtfttr.f"
	*info = -3;
#line 243 "dtfttr.f"
    } else if (*lda < max(1,*n)) {
#line 244 "dtfttr.f"
	*info = -6;
#line 245 "dtfttr.f"
    }
#line 246 "dtfttr.f"
    if (*info != 0) {
#line 247 "dtfttr.f"
	i__1 = -(*info);
#line 247 "dtfttr.f"
	xerbla_("DTFTTR", &i__1, (ftnlen)6);
#line 248 "dtfttr.f"
	return 0;
#line 249 "dtfttr.f"
    }

/*     Quick return if possible */

#line 253 "dtfttr.f"
    if (*n <= 1) {
#line 254 "dtfttr.f"
	if (*n == 1) {
#line 255 "dtfttr.f"
	    a[0] = arf[0];
#line 256 "dtfttr.f"
	}
#line 257 "dtfttr.f"
	return 0;
#line 258 "dtfttr.f"
    }

/*     Size of array ARF(0:nt-1) */

#line 262 "dtfttr.f"
    nt = *n * (*n + 1) / 2;

/*     set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 266 "dtfttr.f"
    if (lower) {
#line 267 "dtfttr.f"
	n2 = *n / 2;
#line 268 "dtfttr.f"
	n1 = *n - n2;
#line 269 "dtfttr.f"
    } else {
#line 270 "dtfttr.f"
	n1 = *n / 2;
#line 271 "dtfttr.f"
	n2 = *n - n1;
#line 272 "dtfttr.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 278 "dtfttr.f"
    if (*n % 2 == 0) {
#line 279 "dtfttr.f"
	k = *n / 2;
#line 280 "dtfttr.f"
	nisodd = FALSE_;
#line 281 "dtfttr.f"
	if (! lower) {
#line 281 "dtfttr.f"
	    np1x2 = *n + *n + 2;
#line 281 "dtfttr.f"
	}
#line 283 "dtfttr.f"
    } else {
#line 284 "dtfttr.f"
	nisodd = TRUE_;
#line 285 "dtfttr.f"
	if (! lower) {
#line 285 "dtfttr.f"
	    nx2 = *n + *n;
#line 285 "dtfttr.f"
	}
#line 287 "dtfttr.f"
    }

#line 289 "dtfttr.f"
    if (nisodd) {

/*        N is odd */

#line 293 "dtfttr.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 297 "dtfttr.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 301 "dtfttr.f"
		ij = 0;
#line 302 "dtfttr.f"
		i__1 = n2;
#line 302 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 303 "dtfttr.f"
		    i__2 = n2 + j;
#line 303 "dtfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 304 "dtfttr.f"
			a[n2 + j + i__ * a_dim1] = arf[ij];
#line 305 "dtfttr.f"
			++ij;
#line 306 "dtfttr.f"
		    }
#line 307 "dtfttr.f"
		    i__2 = *n - 1;
#line 307 "dtfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 308 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 309 "dtfttr.f"
			++ij;
#line 310 "dtfttr.f"
		    }
#line 311 "dtfttr.f"
		}

#line 313 "dtfttr.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 317 "dtfttr.f"
		ij = nt - *n;
#line 318 "dtfttr.f"
		i__1 = n1;
#line 318 "dtfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 319 "dtfttr.f"
		    i__2 = j;
#line 319 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 320 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 321 "dtfttr.f"
			++ij;
#line 322 "dtfttr.f"
		    }
#line 323 "dtfttr.f"
		    i__2 = n1 - 1;
#line 323 "dtfttr.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 324 "dtfttr.f"
			a[j - n1 + l * a_dim1] = arf[ij];
#line 325 "dtfttr.f"
			++ij;
#line 326 "dtfttr.f"
		    }
#line 327 "dtfttr.f"
		    ij -= nx2;
#line 328 "dtfttr.f"
		}

#line 330 "dtfttr.f"
	    }

#line 332 "dtfttr.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 336 "dtfttr.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 340 "dtfttr.f"
		ij = 0;
#line 341 "dtfttr.f"
		i__1 = n2 - 1;
#line 341 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 342 "dtfttr.f"
		    i__2 = j;
#line 342 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 343 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 344 "dtfttr.f"
			++ij;
#line 345 "dtfttr.f"
		    }
#line 346 "dtfttr.f"
		    i__2 = *n - 1;
#line 346 "dtfttr.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 347 "dtfttr.f"
			a[i__ + (n1 + j) * a_dim1] = arf[ij];
#line 348 "dtfttr.f"
			++ij;
#line 349 "dtfttr.f"
		    }
#line 350 "dtfttr.f"
		}
#line 351 "dtfttr.f"
		i__1 = *n - 1;
#line 351 "dtfttr.f"
		for (j = n2; j <= i__1; ++j) {
#line 352 "dtfttr.f"
		    i__2 = n1 - 1;
#line 352 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 353 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 354 "dtfttr.f"
			++ij;
#line 355 "dtfttr.f"
		    }
#line 356 "dtfttr.f"
		}

#line 358 "dtfttr.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 362 "dtfttr.f"
		ij = 0;
#line 363 "dtfttr.f"
		i__1 = n1;
#line 363 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 364 "dtfttr.f"
		    i__2 = *n - 1;
#line 364 "dtfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 365 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 366 "dtfttr.f"
			++ij;
#line 367 "dtfttr.f"
		    }
#line 368 "dtfttr.f"
		}
#line 369 "dtfttr.f"
		i__1 = n1 - 1;
#line 369 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 370 "dtfttr.f"
		    i__2 = j;
#line 370 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 371 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 372 "dtfttr.f"
			++ij;
#line 373 "dtfttr.f"
		    }
#line 374 "dtfttr.f"
		    i__2 = *n - 1;
#line 374 "dtfttr.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 375 "dtfttr.f"
			a[n2 + j + l * a_dim1] = arf[ij];
#line 376 "dtfttr.f"
			++ij;
#line 377 "dtfttr.f"
		    }
#line 378 "dtfttr.f"
		}

#line 380 "dtfttr.f"
	    }

#line 382 "dtfttr.f"
	}

#line 384 "dtfttr.f"
    } else {

/*        N is even */

#line 388 "dtfttr.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 392 "dtfttr.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 396 "dtfttr.f"
		ij = 0;
#line 397 "dtfttr.f"
		i__1 = k - 1;
#line 397 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 398 "dtfttr.f"
		    i__2 = k + j;
#line 398 "dtfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 399 "dtfttr.f"
			a[k + j + i__ * a_dim1] = arf[ij];
#line 400 "dtfttr.f"
			++ij;
#line 401 "dtfttr.f"
		    }
#line 402 "dtfttr.f"
		    i__2 = *n - 1;
#line 402 "dtfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 403 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 404 "dtfttr.f"
			++ij;
#line 405 "dtfttr.f"
		    }
#line 406 "dtfttr.f"
		}

#line 408 "dtfttr.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 412 "dtfttr.f"
		ij = nt - *n - 1;
#line 413 "dtfttr.f"
		i__1 = k;
#line 413 "dtfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 414 "dtfttr.f"
		    i__2 = j;
#line 414 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 415 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 416 "dtfttr.f"
			++ij;
#line 417 "dtfttr.f"
		    }
#line 418 "dtfttr.f"
		    i__2 = k - 1;
#line 418 "dtfttr.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 419 "dtfttr.f"
			a[j - k + l * a_dim1] = arf[ij];
#line 420 "dtfttr.f"
			++ij;
#line 421 "dtfttr.f"
		    }
#line 422 "dtfttr.f"
		    ij -= np1x2;
#line 423 "dtfttr.f"
		}

#line 425 "dtfttr.f"
	    }

#line 427 "dtfttr.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 431 "dtfttr.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 435 "dtfttr.f"
		ij = 0;
#line 436 "dtfttr.f"
		j = k;
#line 437 "dtfttr.f"
		i__1 = *n - 1;
#line 437 "dtfttr.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 438 "dtfttr.f"
		    a[i__ + j * a_dim1] = arf[ij];
#line 439 "dtfttr.f"
		    ++ij;
#line 440 "dtfttr.f"
		}
#line 441 "dtfttr.f"
		i__1 = k - 2;
#line 441 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 442 "dtfttr.f"
		    i__2 = j;
#line 442 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 443 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 444 "dtfttr.f"
			++ij;
#line 445 "dtfttr.f"
		    }
#line 446 "dtfttr.f"
		    i__2 = *n - 1;
#line 446 "dtfttr.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 447 "dtfttr.f"
			a[i__ + (k + 1 + j) * a_dim1] = arf[ij];
#line 448 "dtfttr.f"
			++ij;
#line 449 "dtfttr.f"
		    }
#line 450 "dtfttr.f"
		}
#line 451 "dtfttr.f"
		i__1 = *n - 1;
#line 451 "dtfttr.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 452 "dtfttr.f"
		    i__2 = k - 1;
#line 452 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 453 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 454 "dtfttr.f"
			++ij;
#line 455 "dtfttr.f"
		    }
#line 456 "dtfttr.f"
		}

#line 458 "dtfttr.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 462 "dtfttr.f"
		ij = 0;
#line 463 "dtfttr.f"
		i__1 = k;
#line 463 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 464 "dtfttr.f"
		    i__2 = *n - 1;
#line 464 "dtfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 465 "dtfttr.f"
			a[j + i__ * a_dim1] = arf[ij];
#line 466 "dtfttr.f"
			++ij;
#line 467 "dtfttr.f"
		    }
#line 468 "dtfttr.f"
		}
#line 469 "dtfttr.f"
		i__1 = k - 2;
#line 469 "dtfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 470 "dtfttr.f"
		    i__2 = j;
#line 470 "dtfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 471 "dtfttr.f"
			a[i__ + j * a_dim1] = arf[ij];
#line 472 "dtfttr.f"
			++ij;
#line 473 "dtfttr.f"
		    }
#line 474 "dtfttr.f"
		    i__2 = *n - 1;
#line 474 "dtfttr.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 475 "dtfttr.f"
			a[k + 1 + j + l * a_dim1] = arf[ij];
#line 476 "dtfttr.f"
			++ij;
#line 477 "dtfttr.f"
		    }
#line 478 "dtfttr.f"
		}
/*              Note that here, on exit of the loop, J = K-1 */
#line 480 "dtfttr.f"
		i__1 = j;
#line 480 "dtfttr.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 481 "dtfttr.f"
		    a[i__ + j * a_dim1] = arf[ij];
#line 482 "dtfttr.f"
		    ++ij;
#line 483 "dtfttr.f"
		}

#line 485 "dtfttr.f"
	    }

#line 487 "dtfttr.f"
	}

#line 489 "dtfttr.f"
    }

#line 491 "dtfttr.f"
    return 0;

/*     End of DTFTTR */

} /* dtfttr_ */

