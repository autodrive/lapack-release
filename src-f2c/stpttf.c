#line 1 "stpttf.f"
/* stpttf.f -- translated by f2c (version 20100827).
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

#line 1 "stpttf.f"
/* > \brief \b STPTTF copies a triangular matrix from the standard packed format (TP) to the rectangular full 
packed format (TF). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPTTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpttf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpttf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpttf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPTTF( TRANSR, UPLO, N, AP, ARF, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( 0: * ), ARF( 0: * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPTTF copies a triangular matrix A from standard packed format (TP) */
/* > to rectangular full packed format (TF). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF in Normal format is wanted; */
/* >          = 'T':  ARF in Conjugate-transpose format is wanted. */
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
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A, packed */
/* >          columnwise in a linear array. The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[out] ARF */
/* > \verbatim */
/* >          ARF is REAL array, dimension ( N*(N+1)/2 ), */
/* >          On exit, the upper or lower triangular matrix A stored in */
/* >          RFP format. For a further discussion see Notes below. */
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
/* > */
/*  ===================================================================== */
/* Subroutine */ int stpttf_(char *transr, char *uplo, integer *n, doublereal 
	*ap, doublereal *arf, integer *info, ftnlen transr_len, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, n1, n2, ij, jp, js, nt, lda, ijp;
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

#line 225 "stpttf.f"
    *info = 0;
#line 226 "stpttf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 227 "stpttf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 228 "stpttf.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 229 "stpttf.f"
	*info = -1;
#line 230 "stpttf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 231 "stpttf.f"
	*info = -2;
#line 232 "stpttf.f"
    } else if (*n < 0) {
#line 233 "stpttf.f"
	*info = -3;
#line 234 "stpttf.f"
    }
#line 235 "stpttf.f"
    if (*info != 0) {
#line 236 "stpttf.f"
	i__1 = -(*info);
#line 236 "stpttf.f"
	xerbla_("STPTTF", &i__1, (ftnlen)6);
#line 237 "stpttf.f"
	return 0;
#line 238 "stpttf.f"
    }

/*     Quick return if possible */

#line 242 "stpttf.f"
    if (*n == 0) {
#line 242 "stpttf.f"
	return 0;
#line 242 "stpttf.f"
    }

#line 245 "stpttf.f"
    if (*n == 1) {
#line 246 "stpttf.f"
	if (normaltransr) {
#line 247 "stpttf.f"
	    arf[0] = ap[0];
#line 248 "stpttf.f"
	} else {
#line 249 "stpttf.f"
	    arf[0] = ap[0];
#line 250 "stpttf.f"
	}
#line 251 "stpttf.f"
	return 0;
#line 252 "stpttf.f"
    }

/*     Size of array ARF(0:NT-1) */

#line 256 "stpttf.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER */

#line 260 "stpttf.f"
    if (lower) {
#line 261 "stpttf.f"
	n2 = *n / 2;
#line 262 "stpttf.f"
	n1 = *n - n2;
#line 263 "stpttf.f"
    } else {
#line 264 "stpttf.f"
	n1 = *n / 2;
#line 265 "stpttf.f"
	n2 = *n - n1;
#line 266 "stpttf.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

/*     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe) */
/*     where noe = 0 if n is even, noe = 1 if n is odd */

#line 274 "stpttf.f"
    if (*n % 2 == 0) {
#line 275 "stpttf.f"
	k = *n / 2;
#line 276 "stpttf.f"
	nisodd = FALSE_;
#line 277 "stpttf.f"
	lda = *n + 1;
#line 278 "stpttf.f"
    } else {
#line 279 "stpttf.f"
	nisodd = TRUE_;
#line 280 "stpttf.f"
	lda = *n;
#line 281 "stpttf.f"
    }

/*     ARF^C has lda rows and n+1-noe cols */

#line 285 "stpttf.f"
    if (! normaltransr) {
#line 285 "stpttf.f"
	lda = (*n + 1) / 2;
#line 285 "stpttf.f"
    }

/*     start execution: there are eight cases */

#line 290 "stpttf.f"
    if (nisodd) {

/*        N is odd */

#line 294 "stpttf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 298 "stpttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 302 "stpttf.f"
		ijp = 0;
#line 303 "stpttf.f"
		jp = 0;
#line 304 "stpttf.f"
		i__1 = n2;
#line 304 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 305 "stpttf.f"
		    i__2 = *n - 1;
#line 305 "stpttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 306 "stpttf.f"
			ij = i__ + jp;
#line 307 "stpttf.f"
			arf[ij] = ap[ijp];
#line 308 "stpttf.f"
			++ijp;
#line 309 "stpttf.f"
		    }
#line 310 "stpttf.f"
		    jp += lda;
#line 311 "stpttf.f"
		}
#line 312 "stpttf.f"
		i__1 = n2 - 1;
#line 312 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 313 "stpttf.f"
		    i__2 = n2;
#line 313 "stpttf.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 314 "stpttf.f"
			ij = i__ + j * lda;
#line 315 "stpttf.f"
			arf[ij] = ap[ijp];
#line 316 "stpttf.f"
			++ijp;
#line 317 "stpttf.f"
		    }
#line 318 "stpttf.f"
		}

#line 320 "stpttf.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 324 "stpttf.f"
		ijp = 0;
#line 325 "stpttf.f"
		i__1 = n1 - 1;
#line 325 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 326 "stpttf.f"
		    ij = n2 + j;
#line 327 "stpttf.f"
		    i__2 = j;
#line 327 "stpttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 328 "stpttf.f"
			arf[ij] = ap[ijp];
#line 329 "stpttf.f"
			++ijp;
#line 330 "stpttf.f"
			ij += lda;
#line 331 "stpttf.f"
		    }
#line 332 "stpttf.f"
		}
#line 333 "stpttf.f"
		js = 0;
#line 334 "stpttf.f"
		i__1 = *n - 1;
#line 334 "stpttf.f"
		for (j = n1; j <= i__1; ++j) {
#line 335 "stpttf.f"
		    ij = js;
#line 336 "stpttf.f"
		    i__2 = js + j;
#line 336 "stpttf.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 337 "stpttf.f"
			arf[ij] = ap[ijp];
#line 338 "stpttf.f"
			++ijp;
#line 339 "stpttf.f"
		    }
#line 340 "stpttf.f"
		    js += lda;
#line 341 "stpttf.f"
		}

#line 343 "stpttf.f"
	    }

#line 345 "stpttf.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 349 "stpttf.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 353 "stpttf.f"
		ijp = 0;
#line 354 "stpttf.f"
		i__1 = n2;
#line 354 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 355 "stpttf.f"
		    i__2 = *n * lda - 1;
#line 355 "stpttf.f"
		    i__3 = lda;
#line 355 "stpttf.f"
		    for (ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= 
			    i__2; ij += i__3) {
#line 356 "stpttf.f"
			arf[ij] = ap[ijp];
#line 357 "stpttf.f"
			++ijp;
#line 358 "stpttf.f"
		    }
#line 359 "stpttf.f"
		}
#line 360 "stpttf.f"
		js = 1;
#line 361 "stpttf.f"
		i__1 = n2 - 1;
#line 361 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 362 "stpttf.f"
		    i__3 = js + n2 - j - 1;
#line 362 "stpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 363 "stpttf.f"
			arf[ij] = ap[ijp];
#line 364 "stpttf.f"
			++ijp;
#line 365 "stpttf.f"
		    }
#line 366 "stpttf.f"
		    js = js + lda + 1;
#line 367 "stpttf.f"
		}

#line 369 "stpttf.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 373 "stpttf.f"
		ijp = 0;
#line 374 "stpttf.f"
		js = n2 * lda;
#line 375 "stpttf.f"
		i__1 = n1 - 1;
#line 375 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 376 "stpttf.f"
		    i__3 = js + j;
#line 376 "stpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 377 "stpttf.f"
			arf[ij] = ap[ijp];
#line 378 "stpttf.f"
			++ijp;
#line 379 "stpttf.f"
		    }
#line 380 "stpttf.f"
		    js += lda;
#line 381 "stpttf.f"
		}
#line 382 "stpttf.f"
		i__1 = n1;
#line 382 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 383 "stpttf.f"
		    i__3 = i__ + (n1 + i__) * lda;
#line 383 "stpttf.f"
		    i__2 = lda;
#line 383 "stpttf.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 384 "stpttf.f"
			arf[ij] = ap[ijp];
#line 385 "stpttf.f"
			++ijp;
#line 386 "stpttf.f"
		    }
#line 387 "stpttf.f"
		}

#line 389 "stpttf.f"
	    }

#line 391 "stpttf.f"
	}

#line 393 "stpttf.f"
    } else {

/*        N is even */

#line 397 "stpttf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 401 "stpttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 405 "stpttf.f"
		ijp = 0;
#line 406 "stpttf.f"
		jp = 0;
#line 407 "stpttf.f"
		i__1 = k - 1;
#line 407 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 408 "stpttf.f"
		    i__2 = *n - 1;
#line 408 "stpttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 409 "stpttf.f"
			ij = i__ + 1 + jp;
#line 410 "stpttf.f"
			arf[ij] = ap[ijp];
#line 411 "stpttf.f"
			++ijp;
#line 412 "stpttf.f"
		    }
#line 413 "stpttf.f"
		    jp += lda;
#line 414 "stpttf.f"
		}
#line 415 "stpttf.f"
		i__1 = k - 1;
#line 415 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 416 "stpttf.f"
		    i__2 = k - 1;
#line 416 "stpttf.f"
		    for (j = i__; j <= i__2; ++j) {
#line 417 "stpttf.f"
			ij = i__ + j * lda;
#line 418 "stpttf.f"
			arf[ij] = ap[ijp];
#line 419 "stpttf.f"
			++ijp;
#line 420 "stpttf.f"
		    }
#line 421 "stpttf.f"
		}

#line 423 "stpttf.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 427 "stpttf.f"
		ijp = 0;
#line 428 "stpttf.f"
		i__1 = k - 1;
#line 428 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 429 "stpttf.f"
		    ij = k + 1 + j;
#line 430 "stpttf.f"
		    i__2 = j;
#line 430 "stpttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 431 "stpttf.f"
			arf[ij] = ap[ijp];
#line 432 "stpttf.f"
			++ijp;
#line 433 "stpttf.f"
			ij += lda;
#line 434 "stpttf.f"
		    }
#line 435 "stpttf.f"
		}
#line 436 "stpttf.f"
		js = 0;
#line 437 "stpttf.f"
		i__1 = *n - 1;
#line 437 "stpttf.f"
		for (j = k; j <= i__1; ++j) {
#line 438 "stpttf.f"
		    ij = js;
#line 439 "stpttf.f"
		    i__2 = js + j;
#line 439 "stpttf.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 440 "stpttf.f"
			arf[ij] = ap[ijp];
#line 441 "stpttf.f"
			++ijp;
#line 442 "stpttf.f"
		    }
#line 443 "stpttf.f"
		    js += lda;
#line 444 "stpttf.f"
		}

#line 446 "stpttf.f"
	    }

#line 448 "stpttf.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 452 "stpttf.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 456 "stpttf.f"
		ijp = 0;
#line 457 "stpttf.f"
		i__1 = k - 1;
#line 457 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 458 "stpttf.f"
		    i__2 = (*n + 1) * lda - 1;
#line 458 "stpttf.f"
		    i__3 = lda;
#line 458 "stpttf.f"
		    for (ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : 
			    ij <= i__2; ij += i__3) {
#line 459 "stpttf.f"
			arf[ij] = ap[ijp];
#line 460 "stpttf.f"
			++ijp;
#line 461 "stpttf.f"
		    }
#line 462 "stpttf.f"
		}
#line 463 "stpttf.f"
		js = 0;
#line 464 "stpttf.f"
		i__1 = k - 1;
#line 464 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 465 "stpttf.f"
		    i__3 = js + k - j - 1;
#line 465 "stpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 466 "stpttf.f"
			arf[ij] = ap[ijp];
#line 467 "stpttf.f"
			++ijp;
#line 468 "stpttf.f"
		    }
#line 469 "stpttf.f"
		    js = js + lda + 1;
#line 470 "stpttf.f"
		}

#line 472 "stpttf.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 476 "stpttf.f"
		ijp = 0;
#line 477 "stpttf.f"
		js = (k + 1) * lda;
#line 478 "stpttf.f"
		i__1 = k - 1;
#line 478 "stpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 479 "stpttf.f"
		    i__3 = js + j;
#line 479 "stpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 480 "stpttf.f"
			arf[ij] = ap[ijp];
#line 481 "stpttf.f"
			++ijp;
#line 482 "stpttf.f"
		    }
#line 483 "stpttf.f"
		    js += lda;
#line 484 "stpttf.f"
		}
#line 485 "stpttf.f"
		i__1 = k - 1;
#line 485 "stpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 486 "stpttf.f"
		    i__3 = i__ + (k + i__) * lda;
#line 486 "stpttf.f"
		    i__2 = lda;
#line 486 "stpttf.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 487 "stpttf.f"
			arf[ij] = ap[ijp];
#line 488 "stpttf.f"
			++ijp;
#line 489 "stpttf.f"
		    }
#line 490 "stpttf.f"
		}

#line 492 "stpttf.f"
	    }

#line 494 "stpttf.f"
	}

#line 496 "stpttf.f"
    }

#line 498 "stpttf.f"
    return 0;

/*     End of STPTTF */

} /* stpttf_ */

