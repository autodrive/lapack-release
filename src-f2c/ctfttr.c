#line 1 "ctfttr.f"
/* ctfttr.f -- translated by f2c (version 20100827).
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

#line 1 "ctfttr.f"
/* > \brief \b CTFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard 
full format (TR). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTFTTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctfttr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctfttr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctfttr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANSR, UPLO */
/*       INTEGER            INFO, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( 0: LDA-1, 0: * ), ARF( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTFTTR copies a triangular matrix A from rectangular full packed */
/* > format (TF) to standard full format (TR). */
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
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ARF */
/* > \verbatim */
/* >          ARF is COMPLEX array, dimension ( N*(N+1)/2 ), */
/* >          On entry, the upper or lower triangular matrix A stored in */
/* >          RFP format. For a further discussion see Notes below. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension ( LDA, N ) */
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
/* Subroutine */ int ctfttr_(char *transr, char *uplo, integer *n, 
	doublecomplex *arf, doublecomplex *a, integer *lda, integer *info, 
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 255 "ctfttr.f"
    /* Parameter adjustments */
#line 255 "ctfttr.f"
    a_dim1 = *lda - 1 - 0 + 1;
#line 255 "ctfttr.f"
    a_offset = 0 + a_dim1 * 0;
#line 255 "ctfttr.f"
    a -= a_offset;
#line 255 "ctfttr.f"

#line 255 "ctfttr.f"
    /* Function Body */
#line 255 "ctfttr.f"
    *info = 0;
#line 256 "ctfttr.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 257 "ctfttr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 258 "ctfttr.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 259 "ctfttr.f"
	*info = -1;
#line 260 "ctfttr.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 261 "ctfttr.f"
	*info = -2;
#line 262 "ctfttr.f"
    } else if (*n < 0) {
#line 263 "ctfttr.f"
	*info = -3;
#line 264 "ctfttr.f"
    } else if (*lda < max(1,*n)) {
#line 265 "ctfttr.f"
	*info = -6;
#line 266 "ctfttr.f"
    }
#line 267 "ctfttr.f"
    if (*info != 0) {
#line 268 "ctfttr.f"
	i__1 = -(*info);
#line 268 "ctfttr.f"
	xerbla_("CTFTTR", &i__1, (ftnlen)6);
#line 269 "ctfttr.f"
	return 0;
#line 270 "ctfttr.f"
    }

/*     Quick return if possible */

#line 274 "ctfttr.f"
    if (*n <= 1) {
#line 275 "ctfttr.f"
	if (*n == 1) {
#line 276 "ctfttr.f"
	    if (normaltransr) {
#line 277 "ctfttr.f"
		a[0].r = arf[0].r, a[0].i = arf[0].i;
#line 278 "ctfttr.f"
	    } else {
#line 279 "ctfttr.f"
		d_cnjg(&z__1, arf);
#line 279 "ctfttr.f"
		a[0].r = z__1.r, a[0].i = z__1.i;
#line 280 "ctfttr.f"
	    }
#line 281 "ctfttr.f"
	}
#line 282 "ctfttr.f"
	return 0;
#line 283 "ctfttr.f"
    }

/*     Size of array ARF(1:2,0:nt-1) */

#line 287 "ctfttr.f"
    nt = *n * (*n + 1) / 2;

/*     set N1 and N2 depending on LOWER: for N even N1=N2=K */

#line 291 "ctfttr.f"
    if (lower) {
#line 292 "ctfttr.f"
	n2 = *n / 2;
#line 293 "ctfttr.f"
	n1 = *n - n2;
#line 294 "ctfttr.f"
    } else {
#line 295 "ctfttr.f"
	n1 = *n / 2;
#line 296 "ctfttr.f"
	n2 = *n - n1;
#line 297 "ctfttr.f"
    }

/*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2. */
/*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is */
/*     N--by--(N+1)/2. */

#line 303 "ctfttr.f"
    if (*n % 2 == 0) {
#line 304 "ctfttr.f"
	k = *n / 2;
#line 305 "ctfttr.f"
	nisodd = FALSE_;
#line 306 "ctfttr.f"
	if (! lower) {
#line 306 "ctfttr.f"
	    np1x2 = *n + *n + 2;
#line 306 "ctfttr.f"
	}
#line 308 "ctfttr.f"
    } else {
#line 309 "ctfttr.f"
	nisodd = TRUE_;
#line 310 "ctfttr.f"
	if (! lower) {
#line 310 "ctfttr.f"
	    nx2 = *n + *n;
#line 310 "ctfttr.f"
	}
#line 312 "ctfttr.f"
    }

#line 314 "ctfttr.f"
    if (nisodd) {

/*        N is odd */

#line 318 "ctfttr.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 322 "ctfttr.f"
	    if (lower) {

/*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) ) */
/*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0) */
/*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n */

#line 328 "ctfttr.f"
		ij = 0;
#line 329 "ctfttr.f"
		i__1 = n2;
#line 329 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 330 "ctfttr.f"
		    i__2 = n2 + j;
#line 330 "ctfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 331 "ctfttr.f"
			i__3 = n2 + j + i__ * a_dim1;
#line 331 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 331 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 332 "ctfttr.f"
			++ij;
#line 333 "ctfttr.f"
		    }
#line 334 "ctfttr.f"
		    i__2 = *n - 1;
#line 334 "ctfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 335 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 335 "ctfttr.f"
			i__4 = ij;
#line 335 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 336 "ctfttr.f"
			++ij;
#line 337 "ctfttr.f"
		    }
#line 338 "ctfttr.f"
		}

#line 340 "ctfttr.f"
	    } else {

/*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1) */
/*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0) */
/*             T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n */

#line 346 "ctfttr.f"
		ij = nt - *n;
#line 347 "ctfttr.f"
		i__1 = n1;
#line 347 "ctfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 348 "ctfttr.f"
		    i__2 = j;
#line 348 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 349 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 349 "ctfttr.f"
			i__4 = ij;
#line 349 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 350 "ctfttr.f"
			++ij;
#line 351 "ctfttr.f"
		    }
#line 352 "ctfttr.f"
		    i__2 = n1 - 1;
#line 352 "ctfttr.f"
		    for (l = j - n1; l <= i__2; ++l) {
#line 353 "ctfttr.f"
			i__3 = j - n1 + l * a_dim1;
#line 353 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 353 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 354 "ctfttr.f"
			++ij;
#line 355 "ctfttr.f"
		    }
#line 356 "ctfttr.f"
		    ij -= nx2;
#line 357 "ctfttr.f"
		}

#line 359 "ctfttr.f"
	    }

#line 361 "ctfttr.f"
	} else {

/*           N is odd and TRANSR = 'C' */

#line 365 "ctfttr.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is odd */
/*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1) */
/*              T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1 */

#line 371 "ctfttr.f"
		ij = 0;
#line 372 "ctfttr.f"
		i__1 = n2 - 1;
#line 372 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 373 "ctfttr.f"
		    i__2 = j;
#line 373 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 374 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 374 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 374 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 375 "ctfttr.f"
			++ij;
#line 376 "ctfttr.f"
		    }
#line 377 "ctfttr.f"
		    i__2 = *n - 1;
#line 377 "ctfttr.f"
		    for (i__ = n1 + j; i__ <= i__2; ++i__) {
#line 378 "ctfttr.f"
			i__3 = i__ + (n1 + j) * a_dim1;
#line 378 "ctfttr.f"
			i__4 = ij;
#line 378 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 379 "ctfttr.f"
			++ij;
#line 380 "ctfttr.f"
		    }
#line 381 "ctfttr.f"
		}
#line 382 "ctfttr.f"
		i__1 = *n - 1;
#line 382 "ctfttr.f"
		for (j = n2; j <= i__1; ++j) {
#line 383 "ctfttr.f"
		    i__2 = n1 - 1;
#line 383 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 384 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 384 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 384 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 385 "ctfttr.f"
			++ij;
#line 386 "ctfttr.f"
		    }
#line 387 "ctfttr.f"
		}

#line 389 "ctfttr.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is odd */
/*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0) */
/*              T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda = n2 */

#line 395 "ctfttr.f"
		ij = 0;
#line 396 "ctfttr.f"
		i__1 = n1;
#line 396 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 397 "ctfttr.f"
		    i__2 = *n - 1;
#line 397 "ctfttr.f"
		    for (i__ = n1; i__ <= i__2; ++i__) {
#line 398 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 398 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 398 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 399 "ctfttr.f"
			++ij;
#line 400 "ctfttr.f"
		    }
#line 401 "ctfttr.f"
		}
#line 402 "ctfttr.f"
		i__1 = n1 - 1;
#line 402 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 403 "ctfttr.f"
		    i__2 = j;
#line 403 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 404 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 404 "ctfttr.f"
			i__4 = ij;
#line 404 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 405 "ctfttr.f"
			++ij;
#line 406 "ctfttr.f"
		    }
#line 407 "ctfttr.f"
		    i__2 = *n - 1;
#line 407 "ctfttr.f"
		    for (l = n2 + j; l <= i__2; ++l) {
#line 408 "ctfttr.f"
			i__3 = n2 + j + l * a_dim1;
#line 408 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 408 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 409 "ctfttr.f"
			++ij;
#line 410 "ctfttr.f"
		    }
#line 411 "ctfttr.f"
		}

#line 413 "ctfttr.f"
	    }

#line 415 "ctfttr.f"
	}

#line 417 "ctfttr.f"
    } else {

/*        N is even */

#line 421 "ctfttr.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 425 "ctfttr.f"
	    if (lower) {

/*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0) */
/*              T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1 */

#line 431 "ctfttr.f"
		ij = 0;
#line 432 "ctfttr.f"
		i__1 = k - 1;
#line 432 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 433 "ctfttr.f"
		    i__2 = k + j;
#line 433 "ctfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 434 "ctfttr.f"
			i__3 = k + j + i__ * a_dim1;
#line 434 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 434 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 435 "ctfttr.f"
			++ij;
#line 436 "ctfttr.f"
		    }
#line 437 "ctfttr.f"
		    i__2 = *n - 1;
#line 437 "ctfttr.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 438 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 438 "ctfttr.f"
			i__4 = ij;
#line 438 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 439 "ctfttr.f"
			++ij;
#line 440 "ctfttr.f"
		    }
#line 441 "ctfttr.f"
		}

#line 443 "ctfttr.f"
	    } else {

/*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) ) */
/*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0) */
/*              T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1 */

#line 449 "ctfttr.f"
		ij = nt - *n - 1;
#line 450 "ctfttr.f"
		i__1 = k;
#line 450 "ctfttr.f"
		for (j = *n - 1; j >= i__1; --j) {
#line 451 "ctfttr.f"
		    i__2 = j;
#line 451 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 452 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 452 "ctfttr.f"
			i__4 = ij;
#line 452 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 453 "ctfttr.f"
			++ij;
#line 454 "ctfttr.f"
		    }
#line 455 "ctfttr.f"
		    i__2 = k - 1;
#line 455 "ctfttr.f"
		    for (l = j - k; l <= i__2; ++l) {
#line 456 "ctfttr.f"
			i__3 = j - k + l * a_dim1;
#line 456 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 456 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 457 "ctfttr.f"
			++ij;
#line 458 "ctfttr.f"
		    }
#line 459 "ctfttr.f"
		    ij -= np1x2;
#line 460 "ctfttr.f"
		}

#line 462 "ctfttr.f"
	    }

#line 464 "ctfttr.f"
	} else {

/*           N is even and TRANSR = 'C' */

#line 468 "ctfttr.f"
	    if (lower) {

/*              SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B) */
/*              T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) : */
/*              T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k */

#line 474 "ctfttr.f"
		ij = 0;
#line 475 "ctfttr.f"
		j = k;
#line 476 "ctfttr.f"
		i__1 = *n - 1;
#line 476 "ctfttr.f"
		for (i__ = k; i__ <= i__1; ++i__) {
#line 477 "ctfttr.f"
		    i__2 = i__ + j * a_dim1;
#line 477 "ctfttr.f"
		    i__3 = ij;
#line 477 "ctfttr.f"
		    a[i__2].r = arf[i__3].r, a[i__2].i = arf[i__3].i;
#line 478 "ctfttr.f"
		    ++ij;
#line 479 "ctfttr.f"
		}
#line 480 "ctfttr.f"
		i__1 = k - 2;
#line 480 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 481 "ctfttr.f"
		    i__2 = j;
#line 481 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 482 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 482 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 482 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 483 "ctfttr.f"
			++ij;
#line 484 "ctfttr.f"
		    }
#line 485 "ctfttr.f"
		    i__2 = *n - 1;
#line 485 "ctfttr.f"
		    for (i__ = k + 1 + j; i__ <= i__2; ++i__) {
#line 486 "ctfttr.f"
			i__3 = i__ + (k + 1 + j) * a_dim1;
#line 486 "ctfttr.f"
			i__4 = ij;
#line 486 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 487 "ctfttr.f"
			++ij;
#line 488 "ctfttr.f"
		    }
#line 489 "ctfttr.f"
		}
#line 490 "ctfttr.f"
		i__1 = *n - 1;
#line 490 "ctfttr.f"
		for (j = k - 1; j <= i__1; ++j) {
#line 491 "ctfttr.f"
		    i__2 = k - 1;
#line 491 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 492 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 492 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 492 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 493 "ctfttr.f"
			++ij;
#line 494 "ctfttr.f"
		    }
#line 495 "ctfttr.f"
		}

#line 497 "ctfttr.f"
	    } else {

/*              SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B) */
/*              T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0) */
/*              T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k */

#line 503 "ctfttr.f"
		ij = 0;
#line 504 "ctfttr.f"
		i__1 = k;
#line 504 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 505 "ctfttr.f"
		    i__2 = *n - 1;
#line 505 "ctfttr.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 506 "ctfttr.f"
			i__3 = j + i__ * a_dim1;
#line 506 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 506 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 507 "ctfttr.f"
			++ij;
#line 508 "ctfttr.f"
		    }
#line 509 "ctfttr.f"
		}
#line 510 "ctfttr.f"
		i__1 = k - 2;
#line 510 "ctfttr.f"
		for (j = 0; j <= i__1; ++j) {
#line 511 "ctfttr.f"
		    i__2 = j;
#line 511 "ctfttr.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 512 "ctfttr.f"
			i__3 = i__ + j * a_dim1;
#line 512 "ctfttr.f"
			i__4 = ij;
#line 512 "ctfttr.f"
			a[i__3].r = arf[i__4].r, a[i__3].i = arf[i__4].i;
#line 513 "ctfttr.f"
			++ij;
#line 514 "ctfttr.f"
		    }
#line 515 "ctfttr.f"
		    i__2 = *n - 1;
#line 515 "ctfttr.f"
		    for (l = k + 1 + j; l <= i__2; ++l) {
#line 516 "ctfttr.f"
			i__3 = k + 1 + j + l * a_dim1;
#line 516 "ctfttr.f"
			d_cnjg(&z__1, &arf[ij]);
#line 516 "ctfttr.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 517 "ctfttr.f"
			++ij;
#line 518 "ctfttr.f"
		    }
#line 519 "ctfttr.f"
		}

/*              Note that here J = K-1 */

#line 523 "ctfttr.f"
		i__1 = j;
#line 523 "ctfttr.f"
		for (i__ = 0; i__ <= i__1; ++i__) {
#line 524 "ctfttr.f"
		    i__2 = i__ + j * a_dim1;
#line 524 "ctfttr.f"
		    i__3 = ij;
#line 524 "ctfttr.f"
		    a[i__2].r = arf[i__3].r, a[i__2].i = arf[i__3].i;
#line 525 "ctfttr.f"
		    ++ij;
#line 526 "ctfttr.f"
		}

#line 528 "ctfttr.f"
	    }

#line 530 "ctfttr.f"
	}

#line 532 "ctfttr.f"
    }

#line 534 "ctfttr.f"
    return 0;

/*     End of CTFTTR */

} /* ctfttr_ */

