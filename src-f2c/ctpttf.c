#line 1 "ctpttf.f"
/* ctpttf.f -- translated by f2c (version 20100827).
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

#line 1 "ctpttf.f"
/* > \brief \b CTPTTF copies a triangular matrix from the standard packed format (TP) to the rectangular full 
packed format (TF). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPTTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpttf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpttf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpttf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPTTF( TRANSR, UPLO, N, AP, ARF, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( 0: * ), ARF( 0: * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPTTF copies a triangular matrix A from standard packed format (TP) */
/* > to rectangular full packed format (TF). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF in Normal format is wanted; */
/* >          = 'C':  ARF in Conjugate-transpose format is wanted. */
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
/* >          AP is COMPLEX array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A, packed */
/* >          columnwise in a linear array. The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[out] ARF */
/* > \verbatim */
/* >          ARF is COMPLEX array, dimension ( N*(N+1)/2 ), */
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

/* > \date December 2016 */

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
/* Subroutine */ int ctpttf_(char *transr, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *arf, integer *info, ftnlen 
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 246 "ctpttf.f"
    *info = 0;
#line 247 "ctpttf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 248 "ctpttf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 249 "ctpttf.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 250 "ctpttf.f"
	*info = -1;
#line 251 "ctpttf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 252 "ctpttf.f"
	*info = -2;
#line 253 "ctpttf.f"
    } else if (*n < 0) {
#line 254 "ctpttf.f"
	*info = -3;
#line 255 "ctpttf.f"
    }
#line 256 "ctpttf.f"
    if (*info != 0) {
#line 257 "ctpttf.f"
	i__1 = -(*info);
#line 257 "ctpttf.f"
	xerbla_("CTPTTF", &i__1, (ftnlen)6);
#line 258 "ctpttf.f"
	return 0;
#line 259 "ctpttf.f"
    }

/*     Quick return if possible */

#line 263 "ctpttf.f"
    if (*n == 0) {
#line 263 "ctpttf.f"
	return 0;
#line 263 "ctpttf.f"
    }

#line 266 "ctpttf.f"
    if (*n == 1) {
#line 267 "ctpttf.f"
	if (normaltransr) {
#line 268 "ctpttf.f"
	    arf[0].r = ap[0].r, arf[0].i = ap[0].i;
#line 269 "ctpttf.f"
	} else {
#line 270 "ctpttf.f"
	    d_cnjg(&z__1, ap);
#line 270 "ctpttf.f"
	    arf[0].r = z__1.r, arf[0].i = z__1.i;
#line 271 "ctpttf.f"
	}
#line 272 "ctpttf.f"
	return 0;
#line 273 "ctpttf.f"
    }

/*     Size of array ARF(0:NT-1) */

#line 277 "ctpttf.f"
    nt = *n * (*n + 1) / 2;

/*     Set N1 and N2 depending on LOWER */

#line 281 "ctpttf.f"
    if (lower) {
#line 282 "ctpttf.f"
	n2 = *n / 2;
#line 283 "ctpttf.f"
	n1 = *n - n2;
#line 284 "ctpttf.f"
    } else {
#line 285 "ctpttf.f"
	n1 = *n / 2;
#line 286 "ctpttf.f"
	n2 = *n - n1;
#line 287 "ctpttf.f"
    }

/*     If N is odd, set NISODD = .TRUE. */
/*     If N is even, set K = N/2 and NISODD = .FALSE. */

/*     set lda of ARF^C; ARF^C is (0:(N+1)/2-1,0:N-noe) */
/*     where noe = 0 if n is even, noe = 1 if n is odd */

#line 295 "ctpttf.f"
    if (*n % 2 == 0) {
#line 296 "ctpttf.f"
	k = *n / 2;
#line 297 "ctpttf.f"
	nisodd = FALSE_;
#line 298 "ctpttf.f"
	lda = *n + 1;
#line 299 "ctpttf.f"
    } else {
#line 300 "ctpttf.f"
	nisodd = TRUE_;
#line 301 "ctpttf.f"
	lda = *n;
#line 302 "ctpttf.f"
    }

/*     ARF^C has lda rows and n+1-noe cols */

#line 306 "ctpttf.f"
    if (! normaltransr) {
#line 306 "ctpttf.f"
	lda = (*n + 1) / 2;
#line 306 "ctpttf.f"
    }

/*     start execution: there are eight cases */

#line 311 "ctpttf.f"
    if (nisodd) {

/*        N is odd */

#line 315 "ctpttf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 319 "ctpttf.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda = n */

#line 325 "ctpttf.f"
		ijp = 0;
#line 326 "ctpttf.f"
		jp = 0;
#line 327 "ctpttf.f"
		i__1 = n2;
#line 327 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 328 "ctpttf.f"
		    i__2 = *n - 1;
#line 328 "ctpttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 329 "ctpttf.f"
			ij = i__ + jp;
#line 330 "ctpttf.f"
			i__3 = ij;
#line 330 "ctpttf.f"
			i__4 = ijp;
#line 330 "ctpttf.f"
			arf[i__3].r = ap[i__4].r, arf[i__3].i = ap[i__4].i;
#line 331 "ctpttf.f"
			++ijp;
#line 332 "ctpttf.f"
		    }
#line 333 "ctpttf.f"
		    jp += lda;
#line 334 "ctpttf.f"
		}
#line 335 "ctpttf.f"
		i__1 = n2 - 1;
#line 335 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 336 "ctpttf.f"
		    i__2 = n2;
#line 336 "ctpttf.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 337 "ctpttf.f"
			ij = i__ + j * lda;
#line 338 "ctpttf.f"
			i__3 = ij;
#line 338 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 338 "ctpttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 339 "ctpttf.f"
			++ijp;
#line 340 "ctpttf.f"
		    }
#line 341 "ctpttf.f"
		}

#line 343 "ctpttf.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0) */

#line 349 "ctpttf.f"
		ijp = 0;
#line 350 "ctpttf.f"
		i__1 = n1 - 1;
#line 350 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 351 "ctpttf.f"
		    ij = n2 + j;
#line 352 "ctpttf.f"
		    i__2 = j;
#line 352 "ctpttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 353 "ctpttf.f"
			i__3 = ij;
#line 353 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 353 "ctpttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 354 "ctpttf.f"
			++ijp;
#line 355 "ctpttf.f"
			ij += lda;
#line 356 "ctpttf.f"
		    }
#line 357 "ctpttf.f"
		}
#line 358 "ctpttf.f"
		js = 0;
#line 359 "ctpttf.f"
		i__1 = *n - 1;
#line 359 "ctpttf.f"
		for (j = n1; j <= i__1; ++j) {
#line 360 "ctpttf.f"
		    ij = js;
#line 361 "ctpttf.f"
		    i__2 = js + j;
#line 361 "ctpttf.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 362 "ctpttf.f"
			i__3 = ij;
#line 362 "ctpttf.f"
			i__4 = ijp;
#line 362 "ctpttf.f"
			arf[i__3].r = ap[i__4].r, arf[i__3].i = ap[i__4].i;
#line 363 "ctpttf.f"
			++ijp;
#line 364 "ctpttf.f"
		    }
#line 365 "ctpttf.f"
		    js += lda;
#line 366 "ctpttf.f"
		}

#line 368 "ctpttf.f"
	    }

#line 370 "ctpttf.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 374 "ctpttf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1 */

#line 380 "ctpttf.f"
		ijp = 0;
#line 381 "ctpttf.f"
		i__1 = n2;
#line 381 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 382 "ctpttf.f"
		    i__2 = *n * lda - 1;
#line 382 "ctpttf.f"
		    i__3 = lda;
#line 382 "ctpttf.f"
		    for (ij = i__ * (lda + 1); i__3 < 0 ? ij >= i__2 : ij <= 
			    i__2; ij += i__3) {
#line 383 "ctpttf.f"
			i__4 = ij;
#line 383 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 383 "ctpttf.f"
			arf[i__4].r = z__1.r, arf[i__4].i = z__1.i;
#line 384 "ctpttf.f"
			++ijp;
#line 385 "ctpttf.f"
		    }
#line 386 "ctpttf.f"
		}
#line 387 "ctpttf.f"
		js = 1;
#line 388 "ctpttf.f"
		i__1 = n2 - 1;
#line 388 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 389 "ctpttf.f"
		    i__3 = js + n2 - j - 1;
#line 389 "ctpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 390 "ctpttf.f"
			i__2 = ij;
#line 390 "ctpttf.f"
			i__4 = ijp;
#line 390 "ctpttf.f"
			arf[i__2].r = ap[i__4].r, arf[i__2].i = ap[i__4].i;
#line 391 "ctpttf.f"
			++ijp;
#line 392 "ctpttf.f"
		    }
#line 393 "ctpttf.f"
		    js = js + lda + 1;
#line 394 "ctpttf.f"
		}

#line 396 "ctpttf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2 */

#line 402 "ctpttf.f"
		ijp = 0;
#line 403 "ctpttf.f"
		js = n2 * lda;
#line 404 "ctpttf.f"
		i__1 = n1 - 1;
#line 404 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 405 "ctpttf.f"
		    i__3 = js + j;
#line 405 "ctpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 406 "ctpttf.f"
			i__2 = ij;
#line 406 "ctpttf.f"
			i__4 = ijp;
#line 406 "ctpttf.f"
			arf[i__2].r = ap[i__4].r, arf[i__2].i = ap[i__4].i;
#line 407 "ctpttf.f"
			++ijp;
#line 408 "ctpttf.f"
		    }
#line 409 "ctpttf.f"
		    js += lda;
#line 410 "ctpttf.f"
		}
#line 411 "ctpttf.f"
		i__1 = n1;
#line 411 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 412 "ctpttf.f"
		    i__3 = i__ + (n1 + i__) * lda;
#line 412 "ctpttf.f"
		    i__2 = lda;
#line 412 "ctpttf.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 413 "ctpttf.f"
			i__4 = ij;
#line 413 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 413 "ctpttf.f"
			arf[i__4].r = z__1.r, arf[i__4].i = z__1.i;
#line 414 "ctpttf.f"
			++ijp;
#line 415 "ctpttf.f"
		    }
#line 416 "ctpttf.f"
		}

#line 418 "ctpttf.f"
	    }

#line 420 "ctpttf.f"
	}

#line 422 "ctpttf.f"
    } else {

/*        N is even */

#line 426 "ctpttf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 430 "ctpttf.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1) */

#line 436 "ctpttf.f"
		ijp = 0;
#line 437 "ctpttf.f"
		jp = 0;
#line 438 "ctpttf.f"
		i__1 = k - 1;
#line 438 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 439 "ctpttf.f"
		    i__2 = *n - 1;
#line 439 "ctpttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 440 "ctpttf.f"
			ij = i__ + 1 + jp;
#line 441 "ctpttf.f"
			i__3 = ij;
#line 441 "ctpttf.f"
			i__4 = ijp;
#line 441 "ctpttf.f"
			arf[i__3].r = ap[i__4].r, arf[i__3].i = ap[i__4].i;
#line 442 "ctpttf.f"
			++ijp;
#line 443 "ctpttf.f"
		    }
#line 444 "ctpttf.f"
		    jp += lda;
#line 445 "ctpttf.f"
		}
#line 446 "ctpttf.f"
		i__1 = k - 1;
#line 446 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 447 "ctpttf.f"
		    i__2 = k - 1;
#line 447 "ctpttf.f"
		    for (j = i__; j <= i__2; ++j) {
#line 448 "ctpttf.f"
			ij = i__ + j * lda;
#line 449 "ctpttf.f"
			i__3 = ij;
#line 449 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 449 "ctpttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 450 "ctpttf.f"
			++ijp;
#line 451 "ctpttf.f"
		    }
#line 452 "ctpttf.f"
		}

#line 454 "ctpttf.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0) */

#line 460 "ctpttf.f"
		ijp = 0;
#line 461 "ctpttf.f"
		i__1 = k - 1;
#line 461 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 462 "ctpttf.f"
		    ij = k + 1 + j;
#line 463 "ctpttf.f"
		    i__2 = j;
#line 463 "ctpttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 464 "ctpttf.f"
			i__3 = ij;
#line 464 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 464 "ctpttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 465 "ctpttf.f"
			++ijp;
#line 466 "ctpttf.f"
			ij += lda;
#line 467 "ctpttf.f"
		    }
#line 468 "ctpttf.f"
		}
#line 469 "ctpttf.f"
		js = 0;
#line 470 "ctpttf.f"
		i__1 = *n - 1;
#line 470 "ctpttf.f"
		for (j = k; j <= i__1; ++j) {
#line 471 "ctpttf.f"
		    ij = js;
#line 472 "ctpttf.f"
		    i__2 = js + j;
#line 472 "ctpttf.f"
		    for (ij = js; ij <= i__2; ++ij) {
#line 473 "ctpttf.f"
			i__3 = ij;
#line 473 "ctpttf.f"
			i__4 = ijp;
#line 473 "ctpttf.f"
			arf[i__3].r = ap[i__4].r, arf[i__3].i = ap[i__4].i;
#line 474 "ctpttf.f"
			++ijp;
#line 475 "ctpttf.f"
		    }
#line 476 "ctpttf.f"
		    js += lda;
#line 477 "ctpttf.f"
		}

#line 479 "ctpttf.f"
	    }

#line 481 "ctpttf.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 485 "ctpttf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1) */
/*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k */

#line 491 "ctpttf.f"
		ijp = 0;
#line 492 "ctpttf.f"
		i__1 = k - 1;
#line 492 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 493 "ctpttf.f"
		    i__2 = (*n + 1) * lda - 1;
#line 493 "ctpttf.f"
		    i__3 = lda;
#line 493 "ctpttf.f"
		    for (ij = i__ + (i__ + 1) * lda; i__3 < 0 ? ij >= i__2 : 
			    ij <= i__2; ij += i__3) {
#line 494 "ctpttf.f"
			i__4 = ij;
#line 494 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 494 "ctpttf.f"
			arf[i__4].r = z__1.r, arf[i__4].i = z__1.i;
#line 495 "ctpttf.f"
			++ijp;
#line 496 "ctpttf.f"
		    }
#line 497 "ctpttf.f"
		}
#line 498 "ctpttf.f"
		js = 0;
#line 499 "ctpttf.f"
		i__1 = k - 1;
#line 499 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 500 "ctpttf.f"
		    i__3 = js + k - j - 1;
#line 500 "ctpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 501 "ctpttf.f"
			i__2 = ij;
#line 501 "ctpttf.f"
			i__4 = ijp;
#line 501 "ctpttf.f"
			arf[i__2].r = ap[i__4].r, arf[i__2].i = ap[i__4].i;
#line 502 "ctpttf.f"
			++ijp;
#line 503 "ctpttf.f"
		    }
#line 504 "ctpttf.f"
		    js = js + lda + 1;
#line 505 "ctpttf.f"
		}

#line 507 "ctpttf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper) */
/*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0) */
/*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k */

#line 513 "ctpttf.f"
		ijp = 0;
#line 514 "ctpttf.f"
		js = (k + 1) * lda;
#line 515 "ctpttf.f"
		i__1 = k - 1;
#line 515 "ctpttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 516 "ctpttf.f"
		    i__3 = js + j;
#line 516 "ctpttf.f"
		    for (ij = js; ij <= i__3; ++ij) {
#line 517 "ctpttf.f"
			i__2 = ij;
#line 517 "ctpttf.f"
			i__4 = ijp;
#line 517 "ctpttf.f"
			arf[i__2].r = ap[i__4].r, arf[i__2].i = ap[i__4].i;
#line 518 "ctpttf.f"
			++ijp;
#line 519 "ctpttf.f"
		    }
#line 520 "ctpttf.f"
		    js += lda;
#line 521 "ctpttf.f"
		}
#line 522 "ctpttf.f"
		i__1 = k - 1;
#line 522 "ctpttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 523 "ctpttf.f"
		    i__3 = i__ + (k + i__) * lda;
#line 523 "ctpttf.f"
		    i__2 = lda;
#line 523 "ctpttf.f"
		    for (ij = i__; i__2 < 0 ? ij >= i__3 : ij <= i__3; ij += 
			    i__2) {
#line 524 "ctpttf.f"
			i__4 = ij;
#line 524 "ctpttf.f"
			d_cnjg(&z__1, &ap[ijp]);
#line 524 "ctpttf.f"
			arf[i__4].r = z__1.r, arf[i__4].i = z__1.i;
#line 525 "ctpttf.f"
			++ijp;
#line 526 "ctpttf.f"
		    }
#line 527 "ctpttf.f"
		}

#line 529 "ctpttf.f"
	    }

#line 531 "ctpttf.f"
	}

#line 533 "ctpttf.f"
    }

#line 535 "ctpttf.f"
    return 0;

/*     End of CTPTTF */

} /* ctpttf_ */

