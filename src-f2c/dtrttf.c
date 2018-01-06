#line 1 "dtrttf.f"
/* dtrttf.f -- translated by f2c (version 20100827).
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

#line 1 "dtrttf.f"
/* > \brief \b DTRTTF copies a triangular matrix from the standard full format (TR) to the rectangular full pa
cked format (TF). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRTTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrttf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrttf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrttf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO ) */

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
/* > DTRTTF copies a triangular matrix A from standard full format (TR) */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N). */
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
/* >          ARF is DOUBLE PRECISION array, dimension (NT). */
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
/* Subroutine */ int dtrttf_(char *transr, char *uplo, integer *n, doublereal 
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

#line 231 "dtrttf.f"
    /* Parameter adjustments */
#line 231 "dtrttf.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 231 "dtrttf.f"
    a_offset = 0 + a_dim1 * 0;
#line 231 "dtrttf.f"
    a -= a_offset;
#line 231 "dtrttf.f"

#line 231 "dtrttf.f"
    /* Function Body */
#line 231 "dtrttf.f"
    *info = 0;
#line 232 "dtrttf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 233 "dtrttf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 234 "dtrttf.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 235 "dtrttf.f"
	*info = -1;
#line 236 "dtrttf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 237 "dtrttf.f"
	*info = -2;
#line 238 "dtrttf.f"
    } else if (*n < 0) {
#line 239 "dtrttf.f"
	*info = -3;
#line 240 "dtrttf.f"
    } else if (*lda < max(1,*n)) {
#line 241 "dtrttf.f"
	*info = -5;
#line 242 "dtrttf.f"
    }
#line 243 "dtrttf.f"
    if (*info != 0) {
#line 244 "dtrttf.f"
	i__1 = -(*info);
#line 244 "dtrttf.f"
	xerbla_("DTRTTF", &i__1, (ftnlen)6);
#line 245 "dtrttf.f"
	return 0;
#line 246 "dtrttf.f"
    }

/*     Quick return if possible */

#line 250 "dtrttf.f"
    if (*n <= 1) {
#line 251 "dtrttf.f"
	if (*n == 1) {
#line 252 "dtrttf.f"
	    arf[0] = a[0];
#line 253 "dtrttf.f"
	}
#line 254 "dtrttf.f"
	return 0;
#line 255 "dtrttf.f"
    }

/*     Size of array ARF(0:nt-1) */

#line 259 "dtrttf.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 263 "dtrttf.f"
    if (lower) {
#line 264 "dtrttf.f"
	n2 = *n / 2;
#line 265 "dtrttf.f"
	n1 = *n - n2;
#line 266 "dtrttf.f"
    } else {
#line 267 "dtrttf.f"
	n1 = *n / 2;
#line 268 "dtrttf.f"
	n2 = *n - n1;
#line 269 "dtrttf.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 275 "dtrttf.f"
    if (*n % 2 == 0) {
#line 276 "dtrttf.f"
	k = *n / 2;
#line 277 "dtrttf.f"
	nisodd = FALSE_;
#line 278 "dtrttf.f"
	if (! lower) {
#line 278 "dtrttf.f"
	    np1x2 = *n + *n + 2;
#line 278 "dtrttf.f"
	}
#line 280 "dtrttf.f"
    } else {
#line 281 "dtrttf.f"
	nisodd = TRUE_;
#line 282 "dtrttf.f"
	if (! lower) {
#line 282 "dtrttf.f"
	    nx2 = *n + *n;
#line 282 "dtrttf.f"
	}
#line 284 "dtrttf.f"
    }

#line 286 "dtrttf.f"
    if (nisodd) {

/*        N is odd */

#line 290 "dtrttf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 294 "dtrttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 298 "dtrttf.f"
		ij = 0;
#line 299 "dtrttf.f"
		i__1 = n2;
#line 299 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 300 "dtrttf.f"
		    i__2 = n2 + j;
#line 300 "dtrttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 301 "dtrttf.f"
			arf[ij] = a[n2 + j + i__ * a_dim1];
#line 302 "dtrttf.f"
			++ij;
#line 303 "dtrttf.f"
		    }
#line 304 "dtrttf.f"
		    i__2 = *n - 1;
#line 304 "dtrttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 305 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 306 "dtrttf.f"
			++ij;
#line 307 "dtrttf.f"
		    }
#line 308 "dtrttf.f"
		}

#line 310 "dtrttf.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 314 "dtrttf.f"
		ij = nt - *n;
#line 315 "dtrttf.f"
		i__1 = n1;
#line 315 "dtrttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 316 "dtrttf.f"
		    i__2 = j;
#line 316 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 317 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 318 "dtrttf.f"
			++ij;
#line 319 "dtrttf.f"
		    }
#line 320 "dtrttf.f"
		    i__2 = n1 - 1;
#line 320 "dtrttf.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 321 "dtrttf.f"
			arf[ij] = a[j - n1 + l * a_dim1];
#line 322 "dtrttf.f"
			++ij;
#line 323 "dtrttf.f"
		    }
#line 324 "dtrttf.f"
		    ij -= nx2;
#line 325 "dtrttf.f"
		}

#line 327 "dtrttf.f"
	    }

#line 329 "dtrttf.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 333 "dtrttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 337 "dtrttf.f"
		ij = 0;
#line 338 "dtrttf.f"
		i__1 = n2 - 1;
#line 338 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 339 "dtrttf.f"
		    i__2 = j;
#line 339 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 340 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 341 "dtrttf.f"
			++ij;
#line 342 "dtrttf.f"
		    }
#line 343 "dtrttf.f"
		    i__2 = *n - 1;
#line 343 "dtrttf.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 344 "dtrttf.f"
			arf[ij] = a[i__ + (n1 + j) * a_dim1];
#line 345 "dtrttf.f"
			++ij;
#line 346 "dtrttf.f"
		    }
#line 347 "dtrttf.f"
		}
#line 348 "dtrttf.f"
		i__1 = *n - 1;
#line 348 "dtrttf.f"
		for (j = n2; j <= i__1; ++j) {
#line 349 "dtrttf.f"
		    i__2 = n1 - 1;
#line 349 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 350 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 351 "dtrttf.f"
			++ij;
#line 352 "dtrttf.f"
		    }
#line 353 "dtrttf.f"
		}

#line 355 "dtrttf.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 359 "dtrttf.f"
		ij = 0;
#line 360 "dtrttf.f"
		i__1 = n1;
#line 360 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 361 "dtrttf.f"
		    i__2 = *n - 1;
#line 361 "dtrttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 362 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 363 "dtrttf.f"
			++ij;
#line 364 "dtrttf.f"
		    }
#line 365 "dtrttf.f"
		}
#line 366 "dtrttf.f"
		i__1 = n1 - 1;
#line 366 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 367 "dtrttf.f"
		    i__2 = j;
#line 367 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 368 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 369 "dtrttf.f"
			++ij;
#line 370 "dtrttf.f"
		    }
#line 371 "dtrttf.f"
		    i__2 = *n - 1;
#line 371 "dtrttf.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 372 "dtrttf.f"
			arf[ij] = a[n2 + j + l * a_dim1];
#line 373 "dtrttf.f"
			++ij;
#line 374 "dtrttf.f"
		    }
#line 375 "dtrttf.f"
		}

#line 377 "dtrttf.f"
	    }

#line 379 "dtrttf.f"
	}

#line 381 "dtrttf.f"
    } else {

/*        N is even */

#line 385 "dtrttf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 389 "dtrttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 393 "dtrttf.f"
		ij = 0;
#line 394 "dtrttf.f"
		i__1 = k - 1;
#line 394 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 395 "dtrttf.f"
		    i__2 = k + j;
#line 395 "dtrttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 396 "dtrttf.f"
			arf[ij] = a[k + j + i__ * a_dim1];
#line 397 "dtrttf.f"
			++ij;
#line 398 "dtrttf.f"
		    }
#line 399 "dtrttf.f"
		    i__2 = *n - 1;
#line 399 "dtrttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 400 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 401 "dtrttf.f"
			++ij;
#line 402 "dtrttf.f"
		    }
#line 403 "dtrttf.f"
		}

#line 405 "dtrttf.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 409 "dtrttf.f"
		ij = nt - *n - 1;
#line 410 "dtrttf.f"
		i__1 = k;
#line 410 "dtrttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 411 "dtrttf.f"
		    i__2 = j;
#line 411 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 412 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 413 "dtrttf.f"
			++ij;
#line 414 "dtrttf.f"
		    }
#line 415 "dtrttf.f"
		    i__2 = k - 1;
#line 415 "dtrttf.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 416 "dtrttf.f"
			arf[ij] = a[j - k + l * a_dim1];
#line 417 "dtrttf.f"
			++ij;
#line 418 "dtrttf.f"
		    }
#line 419 "dtrttf.f"
		    ij -= np1x2;
#line 420 "dtrttf.f"
		}

#line 422 "dtrttf.f"
	    }

#line 424 "dtrttf.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 428 "dtrttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 432 "dtrttf.f"
		ij = 0;
#line 433 "dtrttf.f"
		j = k;
#line 434 "dtrttf.f"
		i__1 = *n - 1;
#line 434 "dtrttf.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 435 "dtrttf.f"
		    arf[ij] = a[i__ + j * a_dim1];
#line 436 "dtrttf.f"
		    ++ij;
#line 437 "dtrttf.f"
		}
#line 438 "dtrttf.f"
		i__1 = k - 2;
#line 438 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 439 "dtrttf.f"
		    i__2 = j;
#line 439 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 440 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 441 "dtrttf.f"
			++ij;
#line 442 "dtrttf.f"
		    }
#line 443 "dtrttf.f"
		    i__2 = *n - 1;
#line 443 "dtrttf.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 444 "dtrttf.f"
			arf[ij] = a[i__ + (k + 1 + j) * a_dim1];
#line 445 "dtrttf.f"
			++ij;
#line 446 "dtrttf.f"
		    }
#line 447 "dtrttf.f"
		}
#line 448 "dtrttf.f"
		i__1 = *n - 1;
#line 448 "dtrttf.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 449 "dtrttf.f"
		    i__2 = k - 1;
#line 449 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 450 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 451 "dtrttf.f"
			++ij;
#line 452 "dtrttf.f"
		    }
#line 453 "dtrttf.f"
		}

#line 455 "dtrttf.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 459 "dtrttf.f"
		ij = 0;
#line 460 "dtrttf.f"
		i__1 = k;
#line 460 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 461 "dtrttf.f"
		    i__2 = *n - 1;
#line 461 "dtrttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 462 "dtrttf.f"
			arf[ij] = a[j + i__ * a_dim1];
#line 463 "dtrttf.f"
			++ij;
#line 464 "dtrttf.f"
		    }
#line 465 "dtrttf.f"
		}
#line 466 "dtrttf.f"
		i__1 = k - 2;
#line 466 "dtrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 467 "dtrttf.f"
		    i__2 = j;
#line 467 "dtrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 468 "dtrttf.f"
			arf[ij] = a[i__ + j * a_dim1];
#line 469 "dtrttf.f"
			++ij;
#line 470 "dtrttf.f"
		    }
#line 471 "dtrttf.f"
		    i__2 = *n - 1;
#line 471 "dtrttf.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 472 "dtrttf.f"
			arf[ij] = a[k + 1 + j + l * a_dim1];
#line 473 "dtrttf.f"
			++ij;
#line 474 "dtrttf.f"
		    }
#line 475 "dtrttf.f"
		}
/*              Note that here, on exit of the loop, J = K-1 */
#line 477 "dtrttf.f"
		i__1 = j;
#line 477 "dtrttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 478 "dtrttf.f"
		    arf[ij] = a[i__ + j * a_dim1];
#line 479 "dtrttf.f"
		    ++ij;
#line 480 "dtrttf.f"
		}

#line 482 "dtrttf.f"
	    }

#line 484 "dtrttf.f"
	}

#line 486 "dtrttf.f"
    }

#line 488 "dtrttf.f"
    return 0;

/*     End of DTRTTF */

} /* dtrttf_ */

