#line 1 "stfttp.f"
/* stfttp.f -- translated by f2c (version 20100827).
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

#line 1 "stfttp.f"
/* > \brief \b STFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard 
packed format (TP). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STFTTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfttp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfttp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfttp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STFTTP( TRANSR, UPLO, N, ARF, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STFTTP copies a triangular matrix A from rectangular full packed */
/* > format (TF) to standard packed format (TP). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF is in Normal format; */
/* >          = 'T':  ARF is in Transpose format; */
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
/* >          The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ARF */
/* > \verbatim */
/* >          ARF is REAL array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A stored in */
/* >          RFP format. For a further discussion see Notes below. */
/* > \endverbatim */
/* > */
/* > \param[out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension ( N*(N+1)/2 ), */
/* >          On exit, the upper or lower triangular matrix A, packed */
/* >          columnwise in a linear array. The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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
/* Subroutine */ int stfttp_(char *transr, char *uplo, integer *n, doublereal 
	*arf, doublereal *ap, integer *info, ftnlen transr_len, ftnlen 
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 224 "stfttp.f"
    *info = 0;
#line 225 "stfttp.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 226 "stfttp.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 227 "stfttp.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 228 "stfttp.f"
	*info = -1;
#line 229 "stfttp.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "stfttp.f"
	*info = -2;
#line 231 "stfttp.f"
    } else if (*n < 0) {
#line 232 "stfttp.f"
	*info = -3;
#line 233 "stfttp.f"
    }
#line 234 "stfttp.f"
    if (*info != 0) {
#line 235 "stfttp.f"
	i__1 = -(*info);
#line 235 "stfttp.f"
	xerbla_("STFTTP", &i__1, (ftnlen)6);
#line 236 "stfttp.f"
	return 0;
#line 237 "stfttp.f"
    }

/*     Quick return if possible */

#line 241 "stfttp.f"
    if (*n == 0) {
#line 241 "stfttp.f"
	return 0;
#line 241 "stfttp.f"
    }

#line 244 "stfttp.f"
    if (*n == 1) {
#line 245 "stfttp.f"
	if (normaltransr) {
#line 246 "stfttp.f"
	    ap[0] = arf[0];
#line 247 "stfttp.f"
	} else {
#line 248 "stfttp.f"
	    ap[0] = arf[0];
#line 249 "stfttp.f"
	}
#line 250 "stfttp.f"
	return 0;
#line 251 "stfttp.f"
    }

/*     Size of array ARF(0:NT-1) */

#line 255 "stfttp.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER */

#line 259 "stfttp.f"
    if (lower) {
#line 260 "stfttp.f"
	n2 = *n / 2;
#line 261 "stfttp.f"
	n1 = *n - n2;
#line 262 "stfttp.f"
    } else {
#line 263 "stfttp.f"
	n1 = *n / 2;
#line 264 "stfttp.f"
	n2 = *n - n1;
#line 265 "stfttp.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

/*     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe) */
/*     where noe = 0 if n is even, noe = 1 if n is odd */

#line 273 "stfttp.f"
    if (*n % 2 == 0) {
#line 274 "stfttp.f"
	k = *n / 2;
#line 275 "stfttp.f"
	nisodd = FALSE_;
#line 276 "stfttp.f"
	lda = *n + 1;
#line 277 "stfttp.f"
    } else {
#line 278 "stfttp.f"
	nisodd = TRUE_;
#line 279 "stfttp.f"
	lda = *n;
#line 280 "stfttp.f"
    }

/*     ARF^C has lda rows and n+1-noe cols */

#line 284 "stfttp.f"
    if (! normaltransr) {
#line 284 "stfttp.f"
	lda = (*n + 1) / 2;
#line 284 "stfttp.f"
    }

/*     start execution: there are eight cases */

#line 289 "stfttp.f"
    if (nisodd) {

/*        N is odd */

#line 293 "stfttp.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 297 "stfttp.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n */

#line 303 "stfttp.f"
		ijp = 0;
#line 304 "stfttp.f"
		jp = 0;
#line 305 "stfttp.f"
		i__1 = n2;
#line 305 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 306 "stfttp.f"
		    i__2 = *n - 1;
#line 306 "stfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 307 "stfttp.f"
			ij = i__ + jp;
#line 308 "stfttp.f"
			ap[ijp] = arf[ij];
#line 309 "stfttp.f"
			++ijp;
#line 310 "stfttp.f"
		    }
#line 311 "stfttp.f"
		    jp += lda;
#line 312 "stfttp.f"
		}
#line 313 "stfttp.f"
		i__1 = n2 - 1;
#line 313 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 314 "stfttp.f"
		    i__2 = n2;
#line 314 "stfttp.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 315 "stfttp.f"
			ij = i__ + j * lda;
#line 316 "stfttp.f"
			ap[ijp] = arf[ij];
#line 317 "stfttp.f"
			++ijp;
#line 318 "stfttp.f"
		    }
#line 319 "stfttp.f"
		}

#line 321 "stfttp.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 327 "stfttp.f"
		ijp = 0;
#line 328 "stfttp.f"
		i__1 = n1 - 1;
#line 328 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 329 "stfttp.f"
		    ij = n2 + j;
#line 330 "stfttp.f"
		    i__2 = j;
#line 330 "stfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 331 "stfttp.f"
			ap[ijp] = arf[ij];
#line 332 "stfttp.f"
			++ijp;
#line 333 "stfttp.f"
			ij += lda;
#line 334 "stfttp.f"
		    }
#line 335 "stfttp.f"
		}
#line 336 "stfttp.f"
		js = 0;
#line 337 "stfttp.f"
		i__1 = *n - 1;
#line 337 "stfttp.f"
		for (j = n1; j <= i__1; ++j) {
#line 338 "stfttp.f"
		    ij = js;
#line 339 "stfttp.f"
		    i__2 = js + j;
#line 339 "stfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 340 "stfttp.f"
			ap[ijp] = arf[ij];
#line 341 "stfttp.f"
			++ijp;
#line 342 "stfttp.f"
		    }
#line 343 "stfttp.f"
		    js += lda;
#line 344 "stfttp.f"
		}

#line 346 "stfttp.f"
	    }

#line 348 "stfttp.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 352 "stfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 358 "stfttp.f"
		ijp = 0;
#line 359 "stfttp.f"
		i__1 = n2;
#line 359 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 360 "stfttp.f"
		    i__2 = *n * lda - 1;
#line 360 "stfttp.f"
		    i__3 = lda;
#line 360 "stfttp.f"
		    for (ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= 
			    i__2; ij += i__3) {
#line 361 "stfttp.f"
			ap[ijp] = arf[ij];
#line 362 "stfttp.f"
			++ijp;
#line 363 "stfttp.f"
		    }
#line 364 "stfttp.f"
		}
#line 365 "stfttp.f"
		js = 1;
#line 366 "stfttp.f"
		i__1 = n2 - 1;
#line 366 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 367 "stfttp.f"
		    i__3 = js + n2 - j - 1;
#line 367 "stfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 368 "stfttp.f"
			ap[ijp] = arf[ij];
#line 369 "stfttp.f"
			++ijp;
#line 370 "stfttp.f"
		    }
#line 371 "stfttp.f"
		    js = js + lda + 1;
#line 372 "stfttp.f"
		}

#line 374 "stfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 380 "stfttp.f"
		ijp = 0;
#line 381 "stfttp.f"
		js = n2 * lda;
#line 382 "stfttp.f"
		i__1 = n1 - 1;
#line 382 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 383 "stfttp.f"
		    i__3 = js + j;
#line 383 "stfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 384 "stfttp.f"
			ap[ijp] = arf[ij];
#line 385 "stfttp.f"
			++ijp;
#line 386 "stfttp.f"
		    }
#line 387 "stfttp.f"
		    js += lda;
#line 388 "stfttp.f"
		}
#line 389 "stfttp.f"
		i__1 = n1;
#line 389 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 390 "stfttp.f"
		    i__3 = i__ + (n1 + i__) * lda;
#line 390 "stfttp.f"
		    i__2 = lda;
#line 390 "stfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 391 "stfttp.f"
			ap[ijp] = arf[ij];
#line 392 "stfttp.f"
			++ijp;
#line 393 "stfttp.f"
		    }
#line 394 "stfttp.f"
		}

#line 396 "stfttp.f"
	    }

#line 398 "stfttp.f"
	}

#line 400 "stfttp.f"
    } else {

/*        N is even */

#line 404 "stfttp.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 408 "stfttp.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 414 "stfttp.f"
		ijp = 0;
#line 415 "stfttp.f"
		jp = 0;
#line 416 "stfttp.f"
		i__1 = k - 1;
#line 416 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 417 "stfttp.f"
		    i__2 = *n - 1;
#line 417 "stfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 418 "stfttp.f"
			ij = i__ + 1 + jp;
#line 419 "stfttp.f"
			ap[ijp] = arf[ij];
#line 420 "stfttp.f"
			++ijp;
#line 421 "stfttp.f"
		    }
#line 422 "stfttp.f"
		    jp += lda;
#line 423 "stfttp.f"
		}
#line 424 "stfttp.f"
		i__1 = k - 1;
#line 424 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 425 "stfttp.f"
		    i__2 = k - 1;
#line 425 "stfttp.f"
		    for (j = i__; j <= i__2; ++j) {
#line 426 "stfttp.f"
			ij = i__ + j * lda;
#line 427 "stfttp.f"
			ap[ijp] = arf[ij];
#line 428 "stfttp.f"
			++ijp;
#line 429 "stfttp.f"
		    }
#line 430 "stfttp.f"
		}

#line 432 "stfttp.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 438 "stfttp.f"
		ijp = 0;
#line 439 "stfttp.f"
		i__1 = k - 1;
#line 439 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 440 "stfttp.f"
		    ij = k + 1 + j;
#line 441 "stfttp.f"
		    i__2 = j;
#line 441 "stfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 442 "stfttp.f"
			ap[ijp] = arf[ij];
#line 443 "stfttp.f"
			++ijp;
#line 444 "stfttp.f"
			ij += lda;
#line 445 "stfttp.f"
		    }
#line 446 "stfttp.f"
		}
#line 447 "stfttp.f"
		js = 0;
#line 448 "stfttp.f"
		i__1 = *n - 1;
#line 448 "stfttp.f"
		for (j = k; j <= i__1; ++j) {
#line 449 "stfttp.f"
		    ij = js;
#line 450 "stfttp.f"
		    i__2 = js + j;
#line 450 "stfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 451 "stfttp.f"
			ap[ijp] = arf[ij];
#line 452 "stfttp.f"
			++ijp;
#line 453 "stfttp.f"
		    }
#line 454 "stfttp.f"
		    js += lda;
#line 455 "stfttp.f"
		}

#line 457 "stfttp.f"
	    }

#line 459 "stfttp.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 463 "stfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 469 "stfttp.f"
		ijp = 0;
#line 470 "stfttp.f"
		i__1 = k - 1;
#line 470 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 471 "stfttp.f"
		    i__2 = (*n + 1) * lda - 1;
#line 471 "stfttp.f"
		    i__3 = lda;
#line 471 "stfttp.f"
		    for (ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : 
			    ij <= i__2; ij += i__3) {
#line 472 "stfttp.f"
			ap[ijp] = arf[ij];
#line 473 "stfttp.f"
			++ijp;
#line 474 "stfttp.f"
		    }
#line 475 "stfttp.f"
		}
#line 476 "stfttp.f"
		js = 0;
#line 477 "stfttp.f"
		i__1 = k - 1;
#line 477 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 478 "stfttp.f"
		    i__3 = js + k - j - 1;
#line 478 "stfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 479 "stfttp.f"
			ap[ijp] = arf[ij];
#line 480 "stfttp.f"
			++ijp;
#line 481 "stfttp.f"
		    }
#line 482 "stfttp.f"
		    js = js + lda + 1;
#line 483 "stfttp.f"
		}

#line 485 "stfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 491 "stfttp.f"
		ijp = 0;
#line 492 "stfttp.f"
		js = (k + 1) * lda;
#line 493 "stfttp.f"
		i__1 = k - 1;
#line 493 "stfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 494 "stfttp.f"
		    i__3 = js + j;
#line 494 "stfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 495 "stfttp.f"
			ap[ijp] = arf[ij];
#line 496 "stfttp.f"
			++ijp;
#line 497 "stfttp.f"
		    }
#line 498 "stfttp.f"
		    js += lda;
#line 499 "stfttp.f"
		}
#line 500 "stfttp.f"
		i__1 = k - 1;
#line 500 "stfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 501 "stfttp.f"
		    i__3 = i__ + (k + i__) * lda;
#line 501 "stfttp.f"
		    i__2 = lda;
#line 501 "stfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 502 "stfttp.f"
			ap[ijp] = arf[ij];
#line 503 "stfttp.f"
			++ijp;
#line 504 "stfttp.f"
		    }
#line 505 "stfttp.f"
		}

#line 507 "stfttp.f"
	    }

#line 509 "stfttp.f"
	}

#line 511 "stfttp.f"
    }

#line 513 "stfttp.f"
    return 0;

/*     End of STFTTP */

} /* stfttp_ */

