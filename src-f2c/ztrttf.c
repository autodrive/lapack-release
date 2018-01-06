#line 1 "ztrttf.f"
/* ztrttf.f -- translated by f2c (version 20100827).
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

#line 1 "ztrttf.f"
/* > \brief \b ZTRTTF copies a triangular matrix from the standard full format (TR) to the rectangular full pa
cked format (TF). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTRTTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrttf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrttf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrttf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( 0: LDA-1, 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRTTF copies a triangular matrix A from standard full format (TR) */
/* > to rectangular full packed format (TF) . */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  ARF in Normal mode is wanted; */
/* >          = 'C':  ARF in Conjugate Transpose mode is wanted; */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
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
/* >          The leading dimension of the matrix A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ARF */
/* > \verbatim */
/* >          ARF is COMPLEX*16 array, dimension ( N*(N+1)/2 ), */
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

/* > \ingroup complex16OTHERcomputational */

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
/* Subroutine */ int ztrttf_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *arf, integer *info, 
	ftnlen transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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

#line 254 "ztrttf.f"
    /* Parameter adjustments */
#line 254 "ztrttf.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 254 "ztrttf.f"
    a_offset = 0 + a_dim1 * 0;
#line 254 "ztrttf.f"
    a -= a_offset;
#line 254 "ztrttf.f"

#line 254 "ztrttf.f"
    /* Function Body */
#line 254 "ztrttf.f"
    *info = 0;
#line 255 "ztrttf.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 256 "ztrttf.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 257 "ztrttf.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 258 "ztrttf.f"
	*info = -1;
#line 259 "ztrttf.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 260 "ztrttf.f"
	*info = -2;
#line 261 "ztrttf.f"
    } else if (*n < 0) {
#line 262 "ztrttf.f"
	*info = -3;
#line 263 "ztrttf.f"
    } else if (*lda < max(1,*n)) {
#line 264 "ztrttf.f"
	*info = -5;
#line 265 "ztrttf.f"
    }
#line 266 "ztrttf.f"
    if (*info != 0) {
#line 267 "ztrttf.f"
	i__1 = -(*info);
#line 267 "ztrttf.f"
	xerbla_("ZTRTTF", &i__1, (ftnlen)6);
#line 268 "ztrttf.f"
	return 0;
#line 269 "ztrttf.f"
    }

/*     Quick return if possible */

#line 273 "ztrttf.f"
    if (*n <= 1) {
#line 274 "ztrttf.f"
	if (*n == 1) {
#line 275 "ztrttf.f"
	    if (normaltransr) {
#line 276 "ztrttf.f"
		arf[0].r = a[0].r, arf[0].i = a[0].i;
#line 277 "ztrttf.f"
	    } else {
#line 278 "ztrttf.f"
		d_cnjg(&z__1, a);
#line 278 "ztrttf.f"
		arf[0].r = z__1.r, arf[0].i = z__1.i;
#line 279 "ztrttf.f"
	    }
#line 280 "ztrttf.f"
	}
#line 281 "ztrttf.f"
	return 0;
#line 282 "ztrttf.f"
    }

/*     Size of array ARF(1:2,0:nt-1) */

#line 286 "ztrttf.f"
    nt = *n * (*n + 1) / 2;

/*     set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 290 "ztrttf.f"
    if (lower) {
#line 291 "ztrttf.f"
	n2 = *n / 2;
#line 292 "ztrttf.f"
	n1 = *n - n2;
#line 293 "ztrttf.f"
    } else {
#line 294 "ztrttf.f"
	n1 = *n / 2;
#line 295 "ztrttf.f"
	n2 = *n - n1;
#line 296 "ztrttf.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 302 "ztrttf.f"
    if (*n % 2 == 0) {
#line 303 "ztrttf.f"
	k = *n / 2;
#line 304 "ztrttf.f"
	nisodd = FALSE_;
#line 305 "ztrttf.f"
	if (! lower) {
#line 305 "ztrttf.f"
	    np1x2 = *n + *n + 2;
#line 305 "ztrttf.f"
	}
#line 307 "ztrttf.f"
    } else {
#line 308 "ztrttf.f"
	nisodd = TRUE_;
#line 309 "ztrttf.f"
	if (! lower) {
#line 309 "ztrttf.f"
	    nx2 = *n + *n;
#line 309 "ztrttf.f"
	}
#line 311 "ztrttf.f"
    }

#line 313 "ztrttf.f"
    if (nisodd) {

/*        N is odd */

#line 317 "ztrttf.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 321 "ztrttf.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n */

#line 327 "ztrttf.f"
		ij = 0;
#line 328 "ztrttf.f"
		i__1 = n2;
#line 328 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 329 "ztrttf.f"
		    i__2 = n2 + j;
#line 329 "ztrttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 330 "ztrttf.f"
			i__3 = ij;
#line 330 "ztrttf.f"
			d_cnjg(&z__1, &a[n2 + j + i__ * a_dim1]);
#line 330 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 331 "ztrttf.f"
			++ij;
#line 332 "ztrttf.f"
		    }
#line 333 "ztrttf.f"
		    i__2 = *n - 1;
#line 333 "ztrttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 334 "ztrttf.f"
			i__3 = ij;
#line 334 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 334 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 335 "ztrttf.f"
			++ij;
#line 336 "ztrttf.f"
		    }
#line 337 "ztrttf.f"
		}

#line 339 "ztrttf.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n */

#line 345 "ztrttf.f"
		ij = nt - *n;
#line 346 "ztrttf.f"
		i__1 = n1;
#line 346 "ztrttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 347 "ztrttf.f"
		    i__2 = j;
#line 347 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 348 "ztrttf.f"
			i__3 = ij;
#line 348 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 348 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 349 "ztrttf.f"
			++ij;
#line 350 "ztrttf.f"
		    }
#line 351 "ztrttf.f"
		    i__2 = n1 - 1;
#line 351 "ztrttf.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 352 "ztrttf.f"
			i__3 = ij;
#line 352 "ztrttf.f"
			d_cnjg(&z__1, &a[j - n1 + l * a_dim1]);
#line 352 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 353 "ztrttf.f"
			++ij;
#line 354 "ztrttf.f"
		    }
#line 355 "ztrttf.f"
		    ij -= nx2;
#line 356 "ztrttf.f"
		}

#line 358 "ztrttf.f"
	    }

#line 360 "ztrttf.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 364 "ztrttf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1 */

#line 370 "ztrttf.f"
		ij = 0;
#line 371 "ztrttf.f"
		i__1 = n2 - 1;
#line 371 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 372 "ztrttf.f"
		    i__2 = j;
#line 372 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 373 "ztrttf.f"
			i__3 = ij;
#line 373 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 373 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 374 "ztrttf.f"
			++ij;
#line 375 "ztrttf.f"
		    }
#line 376 "ztrttf.f"
		    i__2 = *n - 1;
#line 376 "ztrttf.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 377 "ztrttf.f"
			i__3 = ij;
#line 377 "ztrttf.f"
			i__4 = i__ + (n1 + j) * a_dim1;
#line 377 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 378 "ztrttf.f"
			++ij;
#line 379 "ztrttf.f"
		    }
#line 380 "ztrttf.f"
		}
#line 381 "ztrttf.f"
		i__1 = *n - 1;
#line 381 "ztrttf.f"
		for (j = n2; j <= i__1; ++j) {
#line 382 "ztrttf.f"
		    i__2 = n1 - 1;
#line 382 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 383 "ztrttf.f"
			i__3 = ij;
#line 383 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 383 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 384 "ztrttf.f"
			++ij;
#line 385 "ztrttf.f"
		    }
#line 386 "ztrttf.f"
		}

#line 388 "ztrttf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda=n2 */

#line 394 "ztrttf.f"
		ij = 0;
#line 395 "ztrttf.f"
		i__1 = n1;
#line 395 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 396 "ztrttf.f"
		    i__2 = *n - 1;
#line 396 "ztrttf.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 397 "ztrttf.f"
			i__3 = ij;
#line 397 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 397 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 398 "ztrttf.f"
			++ij;
#line 399 "ztrttf.f"
		    }
#line 400 "ztrttf.f"
		}
#line 401 "ztrttf.f"
		i__1 = n1 - 1;
#line 401 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 402 "ztrttf.f"
		    i__2 = j;
#line 402 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 403 "ztrttf.f"
			i__3 = ij;
#line 403 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 403 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 404 "ztrttf.f"
			++ij;
#line 405 "ztrttf.f"
		    }
#line 406 "ztrttf.f"
		    i__2 = *n - 1;
#line 406 "ztrttf.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 407 "ztrttf.f"
			i__3 = ij;
#line 407 "ztrttf.f"
			d_cnjg(&z__1, &a[n2 + j + l * a_dim1]);
#line 407 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 408 "ztrttf.f"
			++ij;
#line 409 "ztrttf.f"
		    }
#line 410 "ztrttf.f"
		}

#line 412 "ztrttf.f"
	    }

#line 414 "ztrttf.f"
	}

#line 416 "ztrttf.f"
    } else {

/*        N is even */

#line 420 "ztrttf.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 424 "ztrttf.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1 */

#line 430 "ztrttf.f"
		ij = 0;
#line 431 "ztrttf.f"
		i__1 = k - 1;
#line 431 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 432 "ztrttf.f"
		    i__2 = k + j;
#line 432 "ztrttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 433 "ztrttf.f"
			i__3 = ij;
#line 433 "ztrttf.f"
			d_cnjg(&z__1, &a[k + j + i__ * a_dim1]);
#line 433 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 434 "ztrttf.f"
			++ij;
#line 435 "ztrttf.f"
		    }
#line 436 "ztrttf.f"
		    i__2 = *n - 1;
#line 436 "ztrttf.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 437 "ztrttf.f"
			i__3 = ij;
#line 437 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 437 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 438 "ztrttf.f"
			++ij;
#line 439 "ztrttf.f"
		    }
#line 440 "ztrttf.f"
		}

#line 442 "ztrttf.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1 */

#line 448 "ztrttf.f"
		ij = nt - *n - 1;
#line 449 "ztrttf.f"
		i__1 = k;
#line 449 "ztrttf.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 450 "ztrttf.f"
		    i__2 = j;
#line 450 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 451 "ztrttf.f"
			i__3 = ij;
#line 451 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 451 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 452 "ztrttf.f"
			++ij;
#line 453 "ztrttf.f"
		    }
#line 454 "ztrttf.f"
		    i__2 = k - 1;
#line 454 "ztrttf.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 455 "ztrttf.f"
			i__3 = ij;
#line 455 "ztrttf.f"
			d_cnjg(&z__1, &a[j - k + l * a_dim1]);
#line 455 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 456 "ztrttf.f"
			++ij;
#line 457 "ztrttf.f"
		    }
#line 458 "ztrttf.f"
		    ij -= np1x2;
#line 459 "ztrttf.f"
		}

#line 461 "ztrttf.f"
	    }

#line 463 "ztrttf.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 467 "ztrttf.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B) */
/*              T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) : */
/*              T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k */

#line 473 "ztrttf.f"
		ij = 0;
#line 474 "ztrttf.f"
		j = k;
#line 475 "ztrttf.f"
		i__1 = *n - 1;
#line 475 "ztrttf.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 476 "ztrttf.f"
		    i__2 = ij;
#line 476 "ztrttf.f"
		    i__3 = i__ + j * a_dim1;
#line 476 "ztrttf.f"
		    arf[i__2].r = a[i__3].r, arf[i__2].i = a[i__3].i;
#line 477 "ztrttf.f"
		    ++ij;
#line 478 "ztrttf.f"
		}
#line 479 "ztrttf.f"
		i__1 = k - 2;
#line 479 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 480 "ztrttf.f"
		    i__2 = j;
#line 480 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 481 "ztrttf.f"
			i__3 = ij;
#line 481 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 481 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 482 "ztrttf.f"
			++ij;
#line 483 "ztrttf.f"
		    }
#line 484 "ztrttf.f"
		    i__2 = *n - 1;
#line 484 "ztrttf.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 485 "ztrttf.f"
			i__3 = ij;
#line 485 "ztrttf.f"
			i__4 = i__ + (k + 1 + j) * a_dim1;
#line 485 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 486 "ztrttf.f"
			++ij;
#line 487 "ztrttf.f"
		    }
#line 488 "ztrttf.f"
		}
#line 489 "ztrttf.f"
		i__1 = *n - 1;
#line 489 "ztrttf.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 490 "ztrttf.f"
		    i__2 = k - 1;
#line 490 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 491 "ztrttf.f"
			i__3 = ij;
#line 491 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 491 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 492 "ztrttf.f"
			++ij;
#line 493 "ztrttf.f"
		    }
#line 494 "ztrttf.f"
		}

#line 496 "ztrttf.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B) */
/*              T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0) */
/*              T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k */

#line 502 "ztrttf.f"
		ij = 0;
#line 503 "ztrttf.f"
		i__1 = k;
#line 503 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 504 "ztrttf.f"
		    i__2 = *n - 1;
#line 504 "ztrttf.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 505 "ztrttf.f"
			i__3 = ij;
#line 505 "ztrttf.f"
			d_cnjg(&z__1, &a[j + i__ * a_dim1]);
#line 505 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 506 "ztrttf.f"
			++ij;
#line 507 "ztrttf.f"
		    }
#line 508 "ztrttf.f"
		}
#line 509 "ztrttf.f"
		i__1 = k - 2;
#line 509 "ztrttf.f"
		for (j = 0; j <= i__1; ++j) {
#line 510 "ztrttf.f"
		    i__2 = j;
#line 510 "ztrttf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 511 "ztrttf.f"
			i__3 = ij;
#line 511 "ztrttf.f"
			i__4 = i__ + j * a_dim1;
#line 511 "ztrttf.f"
			arf[i__3].r = a[i__4].r, arf[i__3].i = a[i__4].i;
#line 512 "ztrttf.f"
			++ij;
#line 513 "ztrttf.f"
		    }
#line 514 "ztrttf.f"
		    i__2 = *n - 1;
#line 514 "ztrttf.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 515 "ztrttf.f"
			i__3 = ij;
#line 515 "ztrttf.f"
			d_cnjg(&z__1, &a[k + 1 + j + l * a_dim1]);
#line 515 "ztrttf.f"
			arf[i__3].r = z__1.r, arf[i__3].i = z__1.i;
#line 516 "ztrttf.f"
			++ij;
#line 517 "ztrttf.f"
		    }
#line 518 "ztrttf.f"
		}

/*              Note that here J = K-1 */

#line 522 "ztrttf.f"
		i__1 = j;
#line 522 "ztrttf.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 523 "ztrttf.f"
		    i__2 = ij;
#line 523 "ztrttf.f"
		    i__3 = i__ + j * a_dim1;
#line 523 "ztrttf.f"
		    arf[i__2].r = a[i__3].r, arf[i__2].i = a[i__3].i;
#line 524 "ztrttf.f"
		    ++ij;
#line 525 "ztrttf.f"
		}

#line 527 "ztrttf.f"
	    }

#line 529 "ztrttf.f"
	}

#line 531 "ztrttf.f"
    }

#line 533 "ztrttf.f"
    return 0;

/*     End of ZTRTTF */

} /* ztrttf_ */

