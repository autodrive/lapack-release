#line 1 "dtfttp.f"
/* dtfttp.f -- translated by f2c (version 20100827).
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

#line 1 "dtfttp.f"
/* > \brief \b DTFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard 
packed format (TP). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTFTTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtfttp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtfttp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtfttp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTFTTP( TRANSR, UPLO, N, ARF, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTFTTP copies a triangular matrix A from rectangular full packed */
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
/* >          ARF is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A stored in */
/* >          RFP format. For a further discussion see Notes below. */
/* > \endverbatim */
/* > */
/* > \param[out] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ), */
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
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtfttp_(char *transr, char *uplo, integer *n, doublereal 
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

#line 224 "dtfttp.f"
    *info = 0;
#line 225 "dtfttp.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 226 "dtfttp.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 227 "dtfttp.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 228 "dtfttp.f"
	*info = -1;
#line 229 "dtfttp.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 230 "dtfttp.f"
	*info = -2;
#line 231 "dtfttp.f"
    } else if (*n < 0) {
#line 232 "dtfttp.f"
	*info = -3;
#line 233 "dtfttp.f"
    }
#line 234 "dtfttp.f"
    if (*info != 0) {
#line 235 "dtfttp.f"
	i__1 = -(*info);
#line 235 "dtfttp.f"
	xerbla_("DTFTTP", &i__1, (ftnlen)6);
#line 236 "dtfttp.f"
	return 0;
#line 237 "dtfttp.f"
    }

/*     Quick return if possible */

#line 241 "dtfttp.f"
    if (*n == 0) {
#line 241 "dtfttp.f"
	return 0;
#line 241 "dtfttp.f"
    }

#line 244 "dtfttp.f"
    if (*n == 1) {
#line 245 "dtfttp.f"
	if (normaltransr) {
#line 246 "dtfttp.f"
	    ap[0] = arf[0];
#line 247 "dtfttp.f"
	} else {
#line 248 "dtfttp.f"
	    ap[0] = arf[0];
#line 249 "dtfttp.f"
	}
#line 250 "dtfttp.f"
	return 0;
#line 251 "dtfttp.f"
    }

/*     Size of array ARF(0:NT-1) */

#line 255 "dtfttp.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER */

#line 259 "dtfttp.f"
    if (lower) {
#line 260 "dtfttp.f"
	n2 = *n / 2;
#line 261 "dtfttp.f"
	n1 = *n - n2;
#line 262 "dtfttp.f"
    } else {
#line 263 "dtfttp.f"
	n1 = *n / 2;
#line 264 "dtfttp.f"
	n2 = *n - n1;
#line 265 "dtfttp.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

/*     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe) */
/*     where noe = 0 if n is even, noe = 1 if n is odd */

#line 273 "dtfttp.f"
    if (*n % 2 == 0) {
#line 274 "dtfttp.f"
	k = *n / 2;
#line 275 "dtfttp.f"
	nisodd = FALSE_;
#line 276 "dtfttp.f"
	lda = *n + 1;
#line 277 "dtfttp.f"
    } else {
#line 278 "dtfttp.f"
	nisodd = TRUE_;
#line 279 "dtfttp.f"
	lda = *n;
#line 280 "dtfttp.f"
    }

/*     ARF^C has lda rows and n+1-noe cols */

#line 284 "dtfttp.f"
    if (! normaltransr) {
#line 284 "dtfttp.f"
	lda = (*n + 1) / 2;
#line 284 "dtfttp.f"
    }

/*     start execution: there are eight cases */

#line 289 "dtfttp.f"
    if (nisodd) {

/*        N is odd */

#line 293 "dtfttp.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 297 "dtfttp.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n */

#line 303 "dtfttp.f"
		ijp = 0;
#line 304 "dtfttp.f"
		jp = 0;
#line 305 "dtfttp.f"
		i__1 = n2;
#line 305 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 306 "dtfttp.f"
		    i__2 = *n - 1;
#line 306 "dtfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 307 "dtfttp.f"
			ij = i__ + jp;
#line 308 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 309 "dtfttp.f"
			++ijp;
#line 310 "dtfttp.f"
		    }
#line 311 "dtfttp.f"
		    jp += lda;
#line 312 "dtfttp.f"
		}
#line 313 "dtfttp.f"
		i__1 = n2 - 1;
#line 313 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 314 "dtfttp.f"
		    i__2 = n2;
#line 314 "dtfttp.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 315 "dtfttp.f"
			ij = i__ + j * lda;
#line 316 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 317 "dtfttp.f"
			++ijp;
#line 318 "dtfttp.f"
		    }
#line 319 "dtfttp.f"
		}

#line 321 "dtfttp.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 327 "dtfttp.f"
		ijp = 0;
#line 328 "dtfttp.f"
		i__1 = n1 - 1;
#line 328 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 329 "dtfttp.f"
		    ij = n2 + j;
#line 330 "dtfttp.f"
		    i__2 = j;
#line 330 "dtfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 331 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 332 "dtfttp.f"
			++ijp;
#line 333 "dtfttp.f"
			ij += lda;
#line 334 "dtfttp.f"
		    }
#line 335 "dtfttp.f"
		}
#line 336 "dtfttp.f"
		js = 0;
#line 337 "dtfttp.f"
		i__1 = *n - 1;
#line 337 "dtfttp.f"
		for (j = n1; j <= i__1; ++j) {
#line 338 "dtfttp.f"
		    ij = js;
#line 339 "dtfttp.f"
		    i__2 = js + j;
#line 339 "dtfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 340 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 341 "dtfttp.f"
			++ijp;
#line 342 "dtfttp.f"
		    }
#line 343 "dtfttp.f"
		    js += lda;
#line 344 "dtfttp.f"
		}

#line 346 "dtfttp.f"
	    }

#line 348 "dtfttp.f"
	} else {

/*           N is odd and TRANSR = 'T' */

#line 352 "dtfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 358 "dtfttp.f"
		ijp = 0;
#line 359 "dtfttp.f"
		i__1 = n2;
#line 359 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 360 "dtfttp.f"
		    i__2 = *n * lda - 1;
#line 360 "dtfttp.f"
		    i__3 = lda;
#line 360 "dtfttp.f"
		    for (ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= 
			    i__2; ij += i__3) {
#line 361 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 362 "dtfttp.f"
			++ijp;
#line 363 "dtfttp.f"
		    }
#line 364 "dtfttp.f"
		}
#line 365 "dtfttp.f"
		js = 1;
#line 366 "dtfttp.f"
		i__1 = n2 - 1;
#line 366 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 367 "dtfttp.f"
		    i__3 = js + n2 - j - 1;
#line 367 "dtfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 368 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 369 "dtfttp.f"
			++ijp;
#line 370 "dtfttp.f"
		    }
#line 371 "dtfttp.f"
		    js = js + lda + 1;
#line 372 "dtfttp.f"
		}

#line 374 "dtfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 380 "dtfttp.f"
		ijp = 0;
#line 381 "dtfttp.f"
		js = n2 * lda;
#line 382 "dtfttp.f"
		i__1 = n1 - 1;
#line 382 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 383 "dtfttp.f"
		    i__3 = js + j;
#line 383 "dtfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 384 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 385 "dtfttp.f"
			++ijp;
#line 386 "dtfttp.f"
		    }
#line 387 "dtfttp.f"
		    js += lda;
#line 388 "dtfttp.f"
		}
#line 389 "dtfttp.f"
		i__1 = n1;
#line 389 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 390 "dtfttp.f"
		    i__3 = i__ + (n1 + i__) * lda;
#line 390 "dtfttp.f"
		    i__2 = lda;
#line 390 "dtfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 391 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 392 "dtfttp.f"
			++ijp;
#line 393 "dtfttp.f"
		    }
#line 394 "dtfttp.f"
		}

#line 396 "dtfttp.f"
	    }

#line 398 "dtfttp.f"
	}

#line 400 "dtfttp.f"
    } else {

/*        N is even */

#line 404 "dtfttp.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 408 "dtfttp.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 414 "dtfttp.f"
		ijp = 0;
#line 415 "dtfttp.f"
		jp = 0;
#line 416 "dtfttp.f"
		i__1 = k - 1;
#line 416 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 417 "dtfttp.f"
		    i__2 = *n - 1;
#line 417 "dtfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 418 "dtfttp.f"
			ij = i__ + 1 + jp;
#line 419 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 420 "dtfttp.f"
			++ijp;
#line 421 "dtfttp.f"
		    }
#line 422 "dtfttp.f"
		    jp += lda;
#line 423 "dtfttp.f"
		}
#line 424 "dtfttp.f"
		i__1 = k - 1;
#line 424 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 425 "dtfttp.f"
		    i__2 = k - 1;
#line 425 "dtfttp.f"
		    for (j = i__; j <= i__2; ++j) {
#line 426 "dtfttp.f"
			ij = i__ + j * lda;
#line 427 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 428 "dtfttp.f"
			++ijp;
#line 429 "dtfttp.f"
		    }
#line 430 "dtfttp.f"
		}

#line 432 "dtfttp.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 438 "dtfttp.f"
		ijp = 0;
#line 439 "dtfttp.f"
		i__1 = k - 1;
#line 439 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 440 "dtfttp.f"
		    ij = k + 1 + j;
#line 441 "dtfttp.f"
		    i__2 = j;
#line 441 "dtfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 442 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 443 "dtfttp.f"
			++ijp;
#line 444 "dtfttp.f"
			ij += lda;
#line 445 "dtfttp.f"
		    }
#line 446 "dtfttp.f"
		}
#line 447 "dtfttp.f"
		js = 0;
#line 448 "dtfttp.f"
		i__1 = *n - 1;
#line 448 "dtfttp.f"
		for (j = k; j <= i__1; ++j) {
#line 449 "dtfttp.f"
		    ij = js;
#line 450 "dtfttp.f"
		    i__2 = js + j;
#line 450 "dtfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 451 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 452 "dtfttp.f"
			++ijp;
#line 453 "dtfttp.f"
		    }
#line 454 "dtfttp.f"
		    js += lda;
#line 455 "dtfttp.f"
		}

#line 457 "dtfttp.f"
	    }

#line 459 "dtfttp.f"
	} else {

/*           N is even and TRANSR = 'T' */

#line 463 "dtfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 469 "dtfttp.f"
		ijp = 0;
#line 470 "dtfttp.f"
		i__1 = k - 1;
#line 470 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 471 "dtfttp.f"
		    i__2 = (*n + 1) * lda - 1;
#line 471 "dtfttp.f"
		    i__3 = lda;
#line 471 "dtfttp.f"
		    for (ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : 
			    ij <= i__2; ij += i__3) {
#line 472 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 473 "dtfttp.f"
			++ijp;
#line 474 "dtfttp.f"
		    }
#line 475 "dtfttp.f"
		}
#line 476 "dtfttp.f"
		js = 0;
#line 477 "dtfttp.f"
		i__1 = k - 1;
#line 477 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 478 "dtfttp.f"
		    i__3 = js + k - j - 1;
#line 478 "dtfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 479 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 480 "dtfttp.f"
			++ijp;
#line 481 "dtfttp.f"
		    }
#line 482 "dtfttp.f"
		    js = js + lda + 1;
#line 483 "dtfttp.f"
		}

#line 485 "dtfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 491 "dtfttp.f"
		ijp = 0;
#line 492 "dtfttp.f"
		js = (k + 1) * lda;
#line 493 "dtfttp.f"
		i__1 = k - 1;
#line 493 "dtfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 494 "dtfttp.f"
		    i__3 = js + j;
#line 494 "dtfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 495 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 496 "dtfttp.f"
			++ijp;
#line 497 "dtfttp.f"
		    }
#line 498 "dtfttp.f"
		    js += lda;
#line 499 "dtfttp.f"
		}
#line 500 "dtfttp.f"
		i__1 = k - 1;
#line 500 "dtfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 501 "dtfttp.f"
		    i__3 = i__ + (k + i__) * lda;
#line 501 "dtfttp.f"
		    i__2 = lda;
#line 501 "dtfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 502 "dtfttp.f"
			ap[ijp] = arf[ij];
#line 503 "dtfttp.f"
			++ijp;
#line 504 "dtfttp.f"
		    }
#line 505 "dtfttp.f"
		}

#line 507 "dtfttp.f"
	    }

#line 509 "dtfttp.f"
	}

#line 511 "dtfttp.f"
    }

#line 513 "dtfttp.f"
    return 0;

/*     End of DTFTTP */

} /* dtfttp_ */

