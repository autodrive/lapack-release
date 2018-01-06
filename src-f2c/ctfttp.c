#line 1 "ctfttp.f"
/* ctfttp.f -- translated by f2c (version 20100827).
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

#line 1 "ctfttp.f"
/* > \brief \b CTFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard 
packed format (TP). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTFTTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctfttp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctfttp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctfttp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTFTTP( TRANSR, UPLO, N, ARF, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTFTTP copies a triangular matrix A from rectangular full packed */
/* > format (TF) to standard packed format (TP). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF is in Normal format; */
/* >          = 'C':  ARF is in Conjugate-transpose format; */
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
/* >          ARF is COMPLEX array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A stored in */
/* >          RFP format. For a further discussion see Notes below. */
/* > \endverbatim */
/* > */
/* > \param[out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension ( N*(N+1)/2 ), */
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
/* Subroutine */ int ctfttp_(char *transr, char *uplo, integer *n, 
	doublecomplex *arf, doublecomplex *ap, integer *info, ftnlen 
	transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 250 "ctfttp.f"
    *info = 0;
#line 251 "ctfttp.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 252 "ctfttp.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 253 "ctfttp.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "ctfttp.f"
	*info = -1;
#line 255 "ctfttp.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "ctfttp.f"
	*info = -2;
#line 257 "ctfttp.f"
    } else if (*n < 0) {
#line 258 "ctfttp.f"
	*info = -3;
#line 259 "ctfttp.f"
    }
#line 260 "ctfttp.f"
    if (*info != 0) {
#line 261 "ctfttp.f"
	i__1 = -(*info);
#line 261 "ctfttp.f"
	xerbla_("CTFTTP", &i__1, (ftnlen)6);
#line 262 "ctfttp.f"
	return 0;
#line 263 "ctfttp.f"
    }

/*     Quick return if possible */

#line 267 "ctfttp.f"
    if (*n == 0) {
#line 267 "ctfttp.f"
	return 0;
#line 267 "ctfttp.f"
    }

#line 270 "ctfttp.f"
    if (*n == 1) {
#line 271 "ctfttp.f"
	if (normaltransr) {
#line 272 "ctfttp.f"
	    ap[0].r = arf[0].r, ap[0].i = arf[0].i;
#line 273 "ctfttp.f"
	} else {
#line 274 "ctfttp.f"
	    d_cnjg(&z__1, arf);
#line 274 "ctfttp.f"
	    ap[0].r = z__1.r, ap[0].i = z__1.i;
#line 275 "ctfttp.f"
	}
#line 276 "ctfttp.f"
	return 0;
#line 277 "ctfttp.f"
    }

/*     Size of array ARF(0:NT-1) */

#line 281 "ctfttp.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER */

#line 285 "ctfttp.f"
    if (lower) {
#line 286 "ctfttp.f"
	n2 = *n / 2;
#line 287 "ctfttp.f"
	n1 = *n - n2;
#line 288 "ctfttp.f"
    } else {
#line 289 "ctfttp.f"
	n1 = *n / 2;
#line 290 "ctfttp.f"
	n2 = *n - n1;
#line 291 "ctfttp.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

/*     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe) */
/*     where noe = 0 if n is even, noe = 1 if n is odd */

#line 299 "ctfttp.f"
    if (*n % 2 == 0) {
#line 300 "ctfttp.f"
	k = *n / 2;
#line 301 "ctfttp.f"
	nisodd = FALSE_;
#line 302 "ctfttp.f"
	lda = *n + 1;
#line 303 "ctfttp.f"
    } else {
#line 304 "ctfttp.f"
	nisodd = TRUE_;
#line 305 "ctfttp.f"
	lda = *n;
#line 306 "ctfttp.f"
    }

/*     ARF^C has lda rows and n+1-noe cols */

#line 310 "ctfttp.f"
    if (! normaltransr) {
#line 310 "ctfttp.f"
	lda = (*n + 1) / 2;
#line 310 "ctfttp.f"
    }

/*     start execution: there are eight cases */

#line 315 "ctfttp.f"
    if (nisodd) {

/*        N is odd */

#line 319 "ctfttp.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 323 "ctfttp.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n */

#line 329 "ctfttp.f"
		ijp = 0;
#line 330 "ctfttp.f"
		jp = 0;
#line 331 "ctfttp.f"
		i__1 = n2;
#line 331 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 332 "ctfttp.f"
		    i__2 = *n - 1;
#line 332 "ctfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 333 "ctfttp.f"
			ij = i__ + jp;
#line 334 "ctfttp.f"
			i__3 = ijp;
#line 334 "ctfttp.f"
			i__4 = ij;
#line 334 "ctfttp.f"
			ap[i__3].r = arf[i__4].r, ap[i__3].i = arf[i__4].i;
#line 335 "ctfttp.f"
			++ijp;
#line 336 "ctfttp.f"
		    }
#line 337 "ctfttp.f"
		    jp += lda;
#line 338 "ctfttp.f"
		}
#line 339 "ctfttp.f"
		i__1 = n2 - 1;
#line 339 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 340 "ctfttp.f"
		    i__2 = n2;
#line 340 "ctfttp.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 341 "ctfttp.f"
			ij = i__ + j * lda;
#line 342 "ctfttp.f"
			i__3 = ijp;
#line 342 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 342 "ctfttp.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 343 "ctfttp.f"
			++ijp;
#line 344 "ctfttp.f"
		    }
#line 345 "ctfttp.f"
		}

#line 347 "ctfttp.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 353 "ctfttp.f"
		ijp = 0;
#line 354 "ctfttp.f"
		i__1 = n1 - 1;
#line 354 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 355 "ctfttp.f"
		    ij = n2 + j;
#line 356 "ctfttp.f"
		    i__2 = j;
#line 356 "ctfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 357 "ctfttp.f"
			i__3 = ijp;
#line 357 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 357 "ctfttp.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 358 "ctfttp.f"
			++ijp;
#line 359 "ctfttp.f"
			ij += lda;
#line 360 "ctfttp.f"
		    }
#line 361 "ctfttp.f"
		}
#line 362 "ctfttp.f"
		js = 0;
#line 363 "ctfttp.f"
		i__1 = *n - 1;
#line 363 "ctfttp.f"
		for (j = n1; j <= i__1; ++j) {
#line 364 "ctfttp.f"
		    ij = js;
#line 365 "ctfttp.f"
		    i__2 = js + j;
#line 365 "ctfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 366 "ctfttp.f"
			i__3 = ijp;
#line 366 "ctfttp.f"
			i__4 = ij;
#line 366 "ctfttp.f"
			ap[i__3].r = arf[i__4].r, ap[i__3].i = arf[i__4].i;
#line 367 "ctfttp.f"
			++ijp;
#line 368 "ctfttp.f"
		    }
#line 369 "ctfttp.f"
		    js += lda;
#line 370 "ctfttp.f"
		}

#line 372 "ctfttp.f"
	    }

#line 374 "ctfttp.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 378 "ctfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 384 "ctfttp.f"
		ijp = 0;
#line 385 "ctfttp.f"
		i__1 = n2;
#line 385 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 386 "ctfttp.f"
		    i__2 = *n * lda - 1;
#line 386 "ctfttp.f"
		    i__3 = lda;
#line 386 "ctfttp.f"
		    for (ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= 
			    i__2; ij += i__3) {
#line 387 "ctfttp.f"
			i__4 = ijp;
#line 387 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 387 "ctfttp.f"
			ap[i__4].r = z__1.r, ap[i__4].i = z__1.i;
#line 388 "ctfttp.f"
			++ijp;
#line 389 "ctfttp.f"
		    }
#line 390 "ctfttp.f"
		}
#line 391 "ctfttp.f"
		js = 1;
#line 392 "ctfttp.f"
		i__1 = n2 - 1;
#line 392 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 393 "ctfttp.f"
		    i__3 = js + n2 - j - 1;
#line 393 "ctfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 394 "ctfttp.f"
			i__2 = ijp;
#line 394 "ctfttp.f"
			i__4 = ij;
#line 394 "ctfttp.f"
			ap[i__2].r = arf[i__4].r, ap[i__2].i = arf[i__4].i;
#line 395 "ctfttp.f"
			++ijp;
#line 396 "ctfttp.f"
		    }
#line 397 "ctfttp.f"
		    js = js + lda + 1;
#line 398 "ctfttp.f"
		}

#line 400 "ctfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 406 "ctfttp.f"
		ijp = 0;
#line 407 "ctfttp.f"
		js = n2 * lda;
#line 408 "ctfttp.f"
		i__1 = n1 - 1;
#line 408 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 409 "ctfttp.f"
		    i__3 = js + j;
#line 409 "ctfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 410 "ctfttp.f"
			i__2 = ijp;
#line 410 "ctfttp.f"
			i__4 = ij;
#line 410 "ctfttp.f"
			ap[i__2].r = arf[i__4].r, ap[i__2].i = arf[i__4].i;
#line 411 "ctfttp.f"
			++ijp;
#line 412 "ctfttp.f"
		    }
#line 413 "ctfttp.f"
		    js += lda;
#line 414 "ctfttp.f"
		}
#line 415 "ctfttp.f"
		i__1 = n1;
#line 415 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 416 "ctfttp.f"
		    i__3 = i__ + (n1 + i__) * lda;
#line 416 "ctfttp.f"
		    i__2 = lda;
#line 416 "ctfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 417 "ctfttp.f"
			i__4 = ijp;
#line 417 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 417 "ctfttp.f"
			ap[i__4].r = z__1.r, ap[i__4].i = z__1.i;
#line 418 "ctfttp.f"
			++ijp;
#line 419 "ctfttp.f"
		    }
#line 420 "ctfttp.f"
		}

#line 422 "ctfttp.f"
	    }

#line 424 "ctfttp.f"
	}

#line 426 "ctfttp.f"
    } else {

/*        N is even */

#line 430 "ctfttp.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 434 "ctfttp.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 440 "ctfttp.f"
		ijp = 0;
#line 441 "ctfttp.f"
		jp = 0;
#line 442 "ctfttp.f"
		i__1 = k - 1;
#line 442 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 443 "ctfttp.f"
		    i__2 = *n - 1;
#line 443 "ctfttp.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 444 "ctfttp.f"
			ij = i__ + 1 + jp;
#line 445 "ctfttp.f"
			i__3 = ijp;
#line 445 "ctfttp.f"
			i__4 = ij;
#line 445 "ctfttp.f"
			ap[i__3].r = arf[i__4].r, ap[i__3].i = arf[i__4].i;
#line 446 "ctfttp.f"
			++ijp;
#line 447 "ctfttp.f"
		    }
#line 448 "ctfttp.f"
		    jp += lda;
#line 449 "ctfttp.f"
		}
#line 450 "ctfttp.f"
		i__1 = k - 1;
#line 450 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 451 "ctfttp.f"
		    i__2 = k - 1;
#line 451 "ctfttp.f"
		    for (j = i__; j <= i__2; ++j) {
#line 452 "ctfttp.f"
			ij = i__ + j * lda;
#line 453 "ctfttp.f"
			i__3 = ijp;
#line 453 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 453 "ctfttp.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 454 "ctfttp.f"
			++ijp;
#line 455 "ctfttp.f"
		    }
#line 456 "ctfttp.f"
		}

#line 458 "ctfttp.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 464 "ctfttp.f"
		ijp = 0;
#line 465 "ctfttp.f"
		i__1 = k - 1;
#line 465 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 466 "ctfttp.f"
		    ij = k + 1 + j;
#line 467 "ctfttp.f"
		    i__2 = j;
#line 467 "ctfttp.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 468 "ctfttp.f"
			i__3 = ijp;
#line 468 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 468 "ctfttp.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 469 "ctfttp.f"
			++ijp;
#line 470 "ctfttp.f"
			ij += lda;
#line 471 "ctfttp.f"
		    }
#line 472 "ctfttp.f"
		}
#line 473 "ctfttp.f"
		js = 0;
#line 474 "ctfttp.f"
		i__1 = *n - 1;
#line 474 "ctfttp.f"
		for (j = k; j <= i__1; ++j) {
#line 475 "ctfttp.f"
		    ij = js;
#line 476 "ctfttp.f"
		    i__2 = js + j;
#line 476 "ctfttp.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 477 "ctfttp.f"
			i__3 = ijp;
#line 477 "ctfttp.f"
			i__4 = ij;
#line 477 "ctfttp.f"
			ap[i__3].r = arf[i__4].r, ap[i__3].i = arf[i__4].i;
#line 478 "ctfttp.f"
			++ijp;
#line 479 "ctfttp.f"
		    }
#line 480 "ctfttp.f"
		    js += lda;
#line 481 "ctfttp.f"
		}

#line 483 "ctfttp.f"
	    }

#line 485 "ctfttp.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 489 "ctfttp.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 495 "ctfttp.f"
		ijp = 0;
#line 496 "ctfttp.f"
		i__1 = k - 1;
#line 496 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 497 "ctfttp.f"
		    i__2 = (*n + 1) * lda - 1;
#line 497 "ctfttp.f"
		    i__3 = lda;
#line 497 "ctfttp.f"
		    for (ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : 
			    ij <= i__2; ij += i__3) {
#line 498 "ctfttp.f"
			i__4 = ijp;
#line 498 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 498 "ctfttp.f"
			ap[i__4].r = z__1.r, ap[i__4].i = z__1.i;
#line 499 "ctfttp.f"
			++ijp;
#line 500 "ctfttp.f"
		    }
#line 501 "ctfttp.f"
		}
#line 502 "ctfttp.f"
		js = 0;
#line 503 "ctfttp.f"
		i__1 = k - 1;
#line 503 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 504 "ctfttp.f"
		    i__3 = js + k - j - 1;
#line 504 "ctfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 505 "ctfttp.f"
			i__2 = ijp;
#line 505 "ctfttp.f"
			i__4 = ij;
#line 505 "ctfttp.f"
			ap[i__2].r = arf[i__4].r, ap[i__2].i = arf[i__4].i;
#line 506 "ctfttp.f"
			++ijp;
#line 507 "ctfttp.f"
		    }
#line 508 "ctfttp.f"
		    js = js + lda + 1;
#line 509 "ctfttp.f"
		}

#line 511 "ctfttp.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 517 "ctfttp.f"
		ijp = 0;
#line 518 "ctfttp.f"
		js = (k + 1) * lda;
#line 519 "ctfttp.f"
		i__1 = k - 1;
#line 519 "ctfttp.f"
		for (j = 0; j <= i__1; ++j) {
#line 520 "ctfttp.f"
		    i__3 = js + j;
#line 520 "ctfttp.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 521 "ctfttp.f"
			i__2 = ijp;
#line 521 "ctfttp.f"
			i__4 = ij;
#line 521 "ctfttp.f"
			ap[i__2].r = arf[i__4].r, ap[i__2].i = arf[i__4].i;
#line 522 "ctfttp.f"
			++ijp;
#line 523 "ctfttp.f"
		    }
#line 524 "ctfttp.f"
		    js += lda;
#line 525 "ctfttp.f"
		}
#line 526 "ctfttp.f"
		i__1 = k - 1;
#line 526 "ctfttp.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 527 "ctfttp.f"
		    i__3 = i__ + (k + i__) * lda;
#line 527 "ctfttp.f"
		    i__2 = lda;
#line 527 "ctfttp.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 528 "ctfttp.f"
			i__4 = ijp;
#line 528 "ctfttp.f"
			d_cnjg(&z__1, &arf[ij]);
#line 528 "ctfttp.f"
			ap[i__4].r = z__1.r, ap[i__4].i = z__1.i;
#line 529 "ctfttp.f"
			++ijp;
#line 530 "ctfttp.f"
		    }
#line 531 "ctfttp.f"
		}

#line 533 "ctfttp.f"
	    }

#line 535 "ctfttp.f"
	}

#line 537 "ctfttp.f"
    }

#line 539 "ctfttp.f"
    return 0;

/*     End of CTFTTP */

} /* ctfttp_ */

