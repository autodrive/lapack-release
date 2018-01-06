#line 1 "zlanhf.f"
/* zlanhf.f -- translated by f2c (version 20100827).
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

#line 1 "zlanhf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a Hermitian matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHF( NORM, TRANSR, UPLO, N, A, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, TRANSR, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( 0: * ) */
/*       COMPLEX*16         A( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANHF  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex Hermitian matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \return ZLANHF */
/* > \verbatim */
/* > */
/* >    ZLANHF = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
/* >             ( */
/* >             ( norm1(A),         NORM = '1', 'O' or 'o' */
/* >             ( */
/* >             ( normI(A),         NORM = 'I' or 'i' */
/* >             ( */
/* >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/* > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/* > normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/* > squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER */
/* >            Specifies the value to be returned in ZLANHF as described */
/* >            above. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER */
/* >            Specifies whether the RFP format of A is normal or */
/* >            conjugate-transposed format. */
/* >            = 'N':  RFP format is Normal */
/* >            = 'C':  RFP format is Conjugate-transposed */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER */
/* >            On entry, UPLO specifies whether the RFP matrix A came from */
/* >            an upper or lower triangular matrix as follows: */
/* > */
/* >            UPLO = 'U' or 'u' RFP A came from an upper triangular */
/* >            matrix */
/* > */
/* >            UPLO = 'L' or 'l' RFP A came from a  lower triangular */
/* >            matrix */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >            The order of the matrix A.  N >= 0.  When N = 0, ZLANHF is */
/* >            set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( N*(N+1)/2 ); */
/* >            On entry, the matrix A in RFP Format. */
/* >            RFP Format is described by TRANSR, UPLO and N as follows: */
/* >            If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even; */
/* >            K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If */
/* >            TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A */
/* >            as defined when TRANSR = 'N'. The contents of RFP A are */
/* >            defined by UPLO as follows: If UPLO = 'U' the RFP A */
/* >            contains the ( N*(N+1)/2 ) elements of upper packed A */
/* >            either in normal or conjugate-transpose Format. If */
/* >            UPLO = 'L' the RFP A contains the ( N*(N+1) /2 ) elements */
/* >            of lower packed A either in normal or conjugate-transpose */
/* >            Format. The LDA of RFP A is (N+1)/2 when TRANSR = 'C'. When */
/* >            TRANSR is 'N' the LDA is N+1 when N is even and is N when */
/* >            is odd. See the Note below for more details. */
/* >            Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK), */
/* >            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise, */
/* >            WORK is not referenced. */
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
doublereal zlanhf_(char *norm, char *transr, char *uplo, integer *n, 
	doublecomplex *a, doublereal *work, ftnlen norm_len, ftnlen 
	transr_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal s;
    static integer n1;
    static doublereal aa;
    static integer lda, ifm, noe, ilu;
    static doublereal temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);


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

#line 285 "zlanhf.f"
    if (*n == 0) {
#line 286 "zlanhf.f"
	ret_val = 0.;
#line 287 "zlanhf.f"
	return ret_val;
#line 288 "zlanhf.f"
    } else if (*n == 1) {
#line 289 "zlanhf.f"
	ret_val = z_abs(a);
#line 290 "zlanhf.f"
	return ret_val;
#line 291 "zlanhf.f"
    }

/*     set noe = 1 if n is odd. if n is even set noe=0 */

#line 295 "zlanhf.f"
    noe = 1;
#line 296 "zlanhf.f"
    if (*n % 2 == 0) {
#line 296 "zlanhf.f"
	noe = 0;
#line 296 "zlanhf.f"
    }

/*     set ifm = 0 when form='C' or 'c' and 1 otherwise */

#line 301 "zlanhf.f"
    ifm = 1;
#line 302 "zlanhf.f"
    if (lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 302 "zlanhf.f"
	ifm = 0;
#line 302 "zlanhf.f"
    }

/*     set ilu = 0 when uplo='U or 'u' and 1 otherwise */

#line 307 "zlanhf.f"
    ilu = 1;
#line 308 "zlanhf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 308 "zlanhf.f"
	ilu = 0;
#line 308 "zlanhf.f"
    }

/*     set lda = (n+1)/2 when ifm = 0 */
/*     set lda = n when ifm = 1 and noe = 1 */
/*     set lda = n+1 when ifm = 1 and noe = 0 */

#line 315 "zlanhf.f"
    if (ifm == 1) {
#line 316 "zlanhf.f"
	if (noe == 1) {
#line 317 "zlanhf.f"
	    lda = *n;
#line 318 "zlanhf.f"
	} else {
/*           noe=0 */
#line 320 "zlanhf.f"
	    lda = *n + 1;
#line 321 "zlanhf.f"
	}
#line 322 "zlanhf.f"
    } else {
/*        ifm=0 */
#line 324 "zlanhf.f"
	lda = (*n + 1) / 2;
#line 325 "zlanhf.f"
    }

#line 327 "zlanhf.f"
    if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*       Find max(abs(A(i,j))). */

#line 331 "zlanhf.f"
	k = (*n + 1) / 2;
#line 332 "zlanhf.f"
	value = 0.;
#line 333 "zlanhf.f"
	if (noe == 1) {
/*           n is odd & n = k + k - 1 */
#line 335 "zlanhf.f"
	    if (ifm == 1) {
/*              A is n by k */
#line 337 "zlanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 339 "zlanhf.f"
		    j = 0;
/*                 -> L(0,0) */
#line 341 "zlanhf.f"
		    i__1 = j + j * lda;
#line 341 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 342 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 342 "zlanhf.f"
			value = temp;
#line 342 "zlanhf.f"
		    }
#line 344 "zlanhf.f"
		    i__1 = *n - 1;
#line 344 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 346 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 346 "zlanhf.f"
			    value = temp;
#line 346 "zlanhf.f"
			}
#line 348 "zlanhf.f"
		    }
#line 349 "zlanhf.f"
		    i__1 = k - 1;
#line 349 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 350 "zlanhf.f"
			i__2 = j - 2;
#line 350 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 351 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 352 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 352 "zlanhf.f"
				value = temp;
#line 352 "zlanhf.f"
			    }
#line 354 "zlanhf.f"
			}
#line 355 "zlanhf.f"
			i__ = j - 1;
/*                    L(k+j,k+j) */
#line 357 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 357 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 358 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 358 "zlanhf.f"
			    value = temp;
#line 358 "zlanhf.f"
			}
#line 360 "zlanhf.f"
			i__ = j;
/*                    -> L(j,j) */
#line 362 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 362 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 363 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 363 "zlanhf.f"
			    value = temp;
#line 363 "zlanhf.f"
			}
#line 365 "zlanhf.f"
			i__2 = *n - 1;
#line 365 "zlanhf.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 366 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 367 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 367 "zlanhf.f"
				value = temp;
#line 367 "zlanhf.f"
			    }
#line 369 "zlanhf.f"
			}
#line 370 "zlanhf.f"
		    }
#line 371 "zlanhf.f"
		} else {
/*                 uplo = 'U' */
#line 373 "zlanhf.f"
		    i__1 = k - 2;
#line 373 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 374 "zlanhf.f"
			i__2 = k + j - 2;
#line 374 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 375 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 376 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 376 "zlanhf.f"
				value = temp;
#line 376 "zlanhf.f"
			    }
#line 378 "zlanhf.f"
			}
#line 379 "zlanhf.f"
			i__ = k + j - 1;
/*                    -> U(i,i) */
#line 381 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 381 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 382 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 382 "zlanhf.f"
			    value = temp;
#line 382 "zlanhf.f"
			}
#line 384 "zlanhf.f"
			++i__;
/*                    =k+j; i -> U(j,j) */
#line 386 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 386 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 387 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 387 "zlanhf.f"
			    value = temp;
#line 387 "zlanhf.f"
			}
#line 389 "zlanhf.f"
			i__2 = *n - 1;
#line 389 "zlanhf.f"
			for (i__ = k + j + 1; i__ <= i__2; ++i__) {
#line 390 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 391 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 391 "zlanhf.f"
				value = temp;
#line 391 "zlanhf.f"
			    }
#line 393 "zlanhf.f"
			}
#line 394 "zlanhf.f"
		    }
#line 395 "zlanhf.f"
		    i__1 = *n - 2;
#line 395 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 396 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 397 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 397 "zlanhf.f"
			    value = temp;
#line 397 "zlanhf.f"
			}
/*                    j=k-1 */
#line 400 "zlanhf.f"
		    }
/*                 i=n-1 -> U(n-1,n-1) */
#line 402 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 402 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 403 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 403 "zlanhf.f"
			value = temp;
#line 403 "zlanhf.f"
		    }
#line 405 "zlanhf.f"
		}
#line 406 "zlanhf.f"
	    } else {
/*              xpose case; A is k by n */
#line 408 "zlanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 410 "zlanhf.f"
		    i__1 = k - 2;
#line 410 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 411 "zlanhf.f"
			i__2 = j - 1;
#line 411 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 412 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 413 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 413 "zlanhf.f"
				value = temp;
#line 413 "zlanhf.f"
			    }
#line 415 "zlanhf.f"
			}
#line 416 "zlanhf.f"
			i__ = j;
/*                    L(i,i) */
#line 418 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 418 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 419 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 419 "zlanhf.f"
			    value = temp;
#line 419 "zlanhf.f"
			}
#line 421 "zlanhf.f"
			i__ = j + 1;
/*                    L(j+k,j+k) */
#line 423 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 423 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 424 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 424 "zlanhf.f"
			    value = temp;
#line 424 "zlanhf.f"
			}
#line 426 "zlanhf.f"
			i__2 = k - 1;
#line 426 "zlanhf.f"
			for (i__ = j + 2; i__ <= i__2; ++i__) {
#line 427 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 428 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 428 "zlanhf.f"
				value = temp;
#line 428 "zlanhf.f"
			    }
#line 430 "zlanhf.f"
			}
#line 431 "zlanhf.f"
		    }
#line 432 "zlanhf.f"
		    j = k - 1;
#line 433 "zlanhf.f"
		    i__1 = k - 2;
#line 433 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 434 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 435 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 435 "zlanhf.f"
			    value = temp;
#line 435 "zlanhf.f"
			}
#line 437 "zlanhf.f"
		    }
#line 438 "zlanhf.f"
		    i__ = k - 1;
/*                 -> L(i,i) is at A(i,j) */
#line 440 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 440 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 441 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 441 "zlanhf.f"
			value = temp;
#line 441 "zlanhf.f"
		    }
#line 443 "zlanhf.f"
		    i__1 = *n - 1;
#line 443 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 444 "zlanhf.f"
			i__2 = k - 1;
#line 444 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 445 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 446 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 446 "zlanhf.f"
				value = temp;
#line 446 "zlanhf.f"
			    }
#line 448 "zlanhf.f"
			}
#line 449 "zlanhf.f"
		    }
#line 450 "zlanhf.f"
		} else {
/*                 uplo = 'U' */
#line 452 "zlanhf.f"
		    i__1 = k - 2;
#line 452 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 453 "zlanhf.f"
			i__2 = k - 1;
#line 453 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 454 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 455 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 455 "zlanhf.f"
				value = temp;
#line 455 "zlanhf.f"
			    }
#line 457 "zlanhf.f"
			}
#line 458 "zlanhf.f"
		    }
#line 459 "zlanhf.f"
		    j = k - 1;
/*                 -> U(j,j) is at A(0,j) */
#line 461 "zlanhf.f"
		    i__1 = j * lda;
#line 461 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 462 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 462 "zlanhf.f"
			value = temp;
#line 462 "zlanhf.f"
		    }
#line 464 "zlanhf.f"
		    i__1 = k - 1;
#line 464 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 465 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 466 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 466 "zlanhf.f"
			    value = temp;
#line 466 "zlanhf.f"
			}
#line 468 "zlanhf.f"
		    }
#line 469 "zlanhf.f"
		    i__1 = *n - 1;
#line 469 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 470 "zlanhf.f"
			i__2 = j - k - 1;
#line 470 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 471 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 472 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 472 "zlanhf.f"
				value = temp;
#line 472 "zlanhf.f"
			    }
#line 474 "zlanhf.f"
			}
#line 475 "zlanhf.f"
			i__ = j - k;
/*                    -> U(i,i) at A(i,j) */
#line 477 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 477 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 478 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 478 "zlanhf.f"
			    value = temp;
#line 478 "zlanhf.f"
			}
#line 480 "zlanhf.f"
			i__ = j - k + 1;
/*                    U(j,j) */
#line 482 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 482 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 483 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 483 "zlanhf.f"
			    value = temp;
#line 483 "zlanhf.f"
			}
#line 485 "zlanhf.f"
			i__2 = k - 1;
#line 485 "zlanhf.f"
			for (i__ = j - k + 2; i__ <= i__2; ++i__) {
#line 486 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 487 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 487 "zlanhf.f"
				value = temp;
#line 487 "zlanhf.f"
			    }
#line 489 "zlanhf.f"
			}
#line 490 "zlanhf.f"
		    }
#line 491 "zlanhf.f"
		}
#line 492 "zlanhf.f"
	    }
#line 493 "zlanhf.f"
	} else {
/*           n is even & k = n/2 */
#line 495 "zlanhf.f"
	    if (ifm == 1) {
/*              A is n+1 by k */
#line 497 "zlanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 499 "zlanhf.f"
		    j = 0;
/*                 -> L(k,k) & j=1 -> L(0,0) */
#line 501 "zlanhf.f"
		    i__1 = j + j * lda;
#line 501 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 502 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 502 "zlanhf.f"
			value = temp;
#line 502 "zlanhf.f"
		    }
#line 504 "zlanhf.f"
		    i__1 = j + 1 + j * lda;
#line 504 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 505 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 505 "zlanhf.f"
			value = temp;
#line 505 "zlanhf.f"
		    }
#line 507 "zlanhf.f"
		    i__1 = *n;
#line 507 "zlanhf.f"
		    for (i__ = 2; i__ <= i__1; ++i__) {
#line 508 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 509 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 509 "zlanhf.f"
			    value = temp;
#line 509 "zlanhf.f"
			}
#line 511 "zlanhf.f"
		    }
#line 512 "zlanhf.f"
		    i__1 = k - 1;
#line 512 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 513 "zlanhf.f"
			i__2 = j - 1;
#line 513 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 514 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 515 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 515 "zlanhf.f"
				value = temp;
#line 515 "zlanhf.f"
			    }
#line 517 "zlanhf.f"
			}
#line 518 "zlanhf.f"
			i__ = j;
/*                    L(k+j,k+j) */
#line 520 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 520 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 521 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 521 "zlanhf.f"
			    value = temp;
#line 521 "zlanhf.f"
			}
#line 523 "zlanhf.f"
			i__ = j + 1;
/*                    -> L(j,j) */
#line 525 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 525 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 526 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 526 "zlanhf.f"
			    value = temp;
#line 526 "zlanhf.f"
			}
#line 528 "zlanhf.f"
			i__2 = *n;
#line 528 "zlanhf.f"
			for (i__ = j + 2; i__ <= i__2; ++i__) {
#line 529 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 530 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 530 "zlanhf.f"
				value = temp;
#line 530 "zlanhf.f"
			    }
#line 532 "zlanhf.f"
			}
#line 533 "zlanhf.f"
		    }
#line 534 "zlanhf.f"
		} else {
/*                 uplo = 'U' */
#line 536 "zlanhf.f"
		    i__1 = k - 2;
#line 536 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 537 "zlanhf.f"
			i__2 = k + j - 1;
#line 537 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 538 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 539 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 539 "zlanhf.f"
				value = temp;
#line 539 "zlanhf.f"
			    }
#line 541 "zlanhf.f"
			}
#line 542 "zlanhf.f"
			i__ = k + j;
/*                    -> U(i,i) */
#line 544 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 544 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 545 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 545 "zlanhf.f"
			    value = temp;
#line 545 "zlanhf.f"
			}
#line 547 "zlanhf.f"
			++i__;
/*                    =k+j+1; i -> U(j,j) */
#line 549 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 549 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 550 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 550 "zlanhf.f"
			    value = temp;
#line 550 "zlanhf.f"
			}
#line 552 "zlanhf.f"
			i__2 = *n;
#line 552 "zlanhf.f"
			for (i__ = k + j + 2; i__ <= i__2; ++i__) {
#line 553 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 554 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 554 "zlanhf.f"
				value = temp;
#line 554 "zlanhf.f"
			    }
#line 556 "zlanhf.f"
			}
#line 557 "zlanhf.f"
		    }
#line 558 "zlanhf.f"
		    i__1 = *n - 2;
#line 558 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 559 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 560 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 560 "zlanhf.f"
			    value = temp;
#line 560 "zlanhf.f"
			}
/*                    j=k-1 */
#line 563 "zlanhf.f"
		    }
/*                 i=n-1 -> U(n-1,n-1) */
#line 565 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 565 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 566 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 566 "zlanhf.f"
			value = temp;
#line 566 "zlanhf.f"
		    }
#line 568 "zlanhf.f"
		    i__ = *n;
/*                 -> U(k-1,k-1) */
#line 570 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 570 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 571 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 571 "zlanhf.f"
			value = temp;
#line 571 "zlanhf.f"
		    }
#line 573 "zlanhf.f"
		}
#line 574 "zlanhf.f"
	    } else {
/*              xpose case; A is k by n+1 */
#line 576 "zlanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 578 "zlanhf.f"
		    j = 0;
/*                 -> L(k,k) at A(0,0) */
#line 580 "zlanhf.f"
		    i__1 = j + j * lda;
#line 580 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 581 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 581 "zlanhf.f"
			value = temp;
#line 581 "zlanhf.f"
		    }
#line 583 "zlanhf.f"
		    i__1 = k - 1;
#line 583 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 584 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 585 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 585 "zlanhf.f"
			    value = temp;
#line 585 "zlanhf.f"
			}
#line 587 "zlanhf.f"
		    }
#line 588 "zlanhf.f"
		    i__1 = k - 1;
#line 588 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 589 "zlanhf.f"
			i__2 = j - 2;
#line 589 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 590 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 591 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 591 "zlanhf.f"
				value = temp;
#line 591 "zlanhf.f"
			    }
#line 593 "zlanhf.f"
			}
#line 594 "zlanhf.f"
			i__ = j - 1;
/*                    L(i,i) */
#line 596 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 596 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 597 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 597 "zlanhf.f"
			    value = temp;
#line 597 "zlanhf.f"
			}
#line 599 "zlanhf.f"
			i__ = j;
/*                    L(j+k,j+k) */
#line 601 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 601 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 602 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 602 "zlanhf.f"
			    value = temp;
#line 602 "zlanhf.f"
			}
#line 604 "zlanhf.f"
			i__2 = k - 1;
#line 604 "zlanhf.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 605 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 606 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 606 "zlanhf.f"
				value = temp;
#line 606 "zlanhf.f"
			    }
#line 608 "zlanhf.f"
			}
#line 609 "zlanhf.f"
		    }
#line 610 "zlanhf.f"
		    j = k;
#line 611 "zlanhf.f"
		    i__1 = k - 2;
#line 611 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 612 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 613 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 613 "zlanhf.f"
			    value = temp;
#line 613 "zlanhf.f"
			}
#line 615 "zlanhf.f"
		    }
#line 616 "zlanhf.f"
		    i__ = k - 1;
/*                 -> L(i,i) is at A(i,j) */
#line 618 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 618 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 619 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 619 "zlanhf.f"
			value = temp;
#line 619 "zlanhf.f"
		    }
#line 621 "zlanhf.f"
		    i__1 = *n;
#line 621 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 622 "zlanhf.f"
			i__2 = k - 1;
#line 622 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 623 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 624 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 624 "zlanhf.f"
				value = temp;
#line 624 "zlanhf.f"
			    }
#line 626 "zlanhf.f"
			}
#line 627 "zlanhf.f"
		    }
#line 628 "zlanhf.f"
		} else {
/*                 uplo = 'U' */
#line 630 "zlanhf.f"
		    i__1 = k - 1;
#line 630 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 631 "zlanhf.f"
			i__2 = k - 1;
#line 631 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 632 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 633 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 633 "zlanhf.f"
				value = temp;
#line 633 "zlanhf.f"
			    }
#line 635 "zlanhf.f"
			}
#line 636 "zlanhf.f"
		    }
#line 637 "zlanhf.f"
		    j = k;
/*                 -> U(j,j) is at A(0,j) */
#line 639 "zlanhf.f"
		    i__1 = j * lda;
#line 639 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 640 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 640 "zlanhf.f"
			value = temp;
#line 640 "zlanhf.f"
		    }
#line 642 "zlanhf.f"
		    i__1 = k - 1;
#line 642 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 643 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 644 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 644 "zlanhf.f"
			    value = temp;
#line 644 "zlanhf.f"
			}
#line 646 "zlanhf.f"
		    }
#line 647 "zlanhf.f"
		    i__1 = *n - 1;
#line 647 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 648 "zlanhf.f"
			i__2 = j - k - 2;
#line 648 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 649 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 650 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 650 "zlanhf.f"
				value = temp;
#line 650 "zlanhf.f"
			    }
#line 652 "zlanhf.f"
			}
#line 653 "zlanhf.f"
			i__ = j - k - 1;
/*                    -> U(i,i) at A(i,j) */
#line 655 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 655 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 656 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 656 "zlanhf.f"
			    value = temp;
#line 656 "zlanhf.f"
			}
#line 658 "zlanhf.f"
			i__ = j - k;
/*                    U(j,j) */
#line 660 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 660 "zlanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 661 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 661 "zlanhf.f"
			    value = temp;
#line 661 "zlanhf.f"
			}
#line 663 "zlanhf.f"
			i__2 = k - 1;
#line 663 "zlanhf.f"
			for (i__ = j - k + 1; i__ <= i__2; ++i__) {
#line 664 "zlanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 665 "zlanhf.f"
			    if (value < temp || disnan_(&temp)) {
#line 665 "zlanhf.f"
				value = temp;
#line 665 "zlanhf.f"
			    }
#line 667 "zlanhf.f"
			}
#line 668 "zlanhf.f"
		    }
#line 669 "zlanhf.f"
		    j = *n;
#line 670 "zlanhf.f"
		    i__1 = k - 2;
#line 670 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 671 "zlanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 672 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 672 "zlanhf.f"
			    value = temp;
#line 672 "zlanhf.f"
			}
#line 674 "zlanhf.f"
		    }
#line 675 "zlanhf.f"
		    i__ = k - 1;
/*                 U(k,k) at A(i,j) */
#line 677 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 677 "zlanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 678 "zlanhf.f"
		    if (value < temp || disnan_(&temp)) {
#line 678 "zlanhf.f"
			value = temp;
#line 678 "zlanhf.f"
		    }
#line 680 "zlanhf.f"
		}
#line 681 "zlanhf.f"
	    }
#line 682 "zlanhf.f"
	}
#line 683 "zlanhf.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*       Find normI(A) ( = norm1(A), since A is Hermitian). */

#line 688 "zlanhf.f"
	if (ifm == 1) {
/*           A is 'N' */
#line 690 "zlanhf.f"
	    k = *n / 2;
#line 691 "zlanhf.f"
	    if (noe == 1) {
/*              n is odd & A is n by (n+1)/2 */
#line 693 "zlanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 695 "zlanhf.f"
		    i__1 = k - 1;
#line 695 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 696 "zlanhf.f"
			work[i__] = 0.;
#line 697 "zlanhf.f"
		    }
#line 698 "zlanhf.f"
		    i__1 = k;
#line 698 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 699 "zlanhf.f"
			s = 0.;
#line 700 "zlanhf.f"
			i__2 = k + j - 1;
#line 700 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 701 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(i,j+k) */
#line 703 "zlanhf.f"
			    s += aa;
#line 704 "zlanhf.f"
			    work[i__] += aa;
#line 705 "zlanhf.f"
			}
#line 706 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 706 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 708 "zlanhf.f"
			work[j + k] = s + aa;
#line 709 "zlanhf.f"
			if (i__ == k + k) {
#line 709 "zlanhf.f"
			    goto L10;
#line 709 "zlanhf.f"
			}
#line 711 "zlanhf.f"
			++i__;
#line 712 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 712 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j,j) */
#line 714 "zlanhf.f"
			work[j] += aa;
#line 715 "zlanhf.f"
			s = 0.;
#line 716 "zlanhf.f"
			i__2 = k - 1;
#line 716 "zlanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 717 "zlanhf.f"
			    ++i__;
#line 718 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 720 "zlanhf.f"
			    s += aa;
#line 721 "zlanhf.f"
			    work[l] += aa;
#line 722 "zlanhf.f"
			}
#line 723 "zlanhf.f"
			work[j] += s;
#line 724 "zlanhf.f"
		    }
#line 725 "zlanhf.f"
L10:
#line 726 "zlanhf.f"
		    value = work[0];
#line 727 "zlanhf.f"
		    i__1 = *n - 1;
#line 727 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 728 "zlanhf.f"
			temp = work[i__];
#line 729 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 729 "zlanhf.f"
			    value = temp;
#line 729 "zlanhf.f"
			}
#line 731 "zlanhf.f"
		    }
#line 732 "zlanhf.f"
		} else {
/*                 ilu = 1 & uplo = 'L' */
#line 734 "zlanhf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 736 "zlanhf.f"
		    i__1 = *n - 1;
#line 736 "zlanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 737 "zlanhf.f"
			work[i__] = 0.;
#line 738 "zlanhf.f"
		    }
#line 739 "zlanhf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 740 "zlanhf.f"
			s = 0.;
#line 741 "zlanhf.f"
			i__1 = j - 2;
#line 741 "zlanhf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 742 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(j+k,i+k) */
#line 744 "zlanhf.f"
			    s += aa;
#line 745 "zlanhf.f"
			    work[i__ + k] += aa;
#line 746 "zlanhf.f"
			}
#line 747 "zlanhf.f"
			if (j > 0) {
#line 748 "zlanhf.f"
			    i__1 = i__ + j * lda;
#line 748 "zlanhf.f"
			    aa = (d__1 = a[i__1].r, abs(d__1));
/*                       -> A(j+k,j+k) */
#line 750 "zlanhf.f"
			    s += aa;
#line 751 "zlanhf.f"
			    work[i__ + k] += s;
/*                       i=j */
#line 753 "zlanhf.f"
			    ++i__;
#line 754 "zlanhf.f"
			}
#line 755 "zlanhf.f"
			i__1 = i__ + j * lda;
#line 755 "zlanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j,j) */
#line 757 "zlanhf.f"
			work[j] = aa;
#line 758 "zlanhf.f"
			s = 0.;
#line 759 "zlanhf.f"
			i__1 = *n - 1;
#line 759 "zlanhf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 760 "zlanhf.f"
			    ++i__;
#line 761 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 763 "zlanhf.f"
			    s += aa;
#line 764 "zlanhf.f"
			    work[l] += aa;
#line 765 "zlanhf.f"
			}
#line 766 "zlanhf.f"
			work[j] += s;
#line 767 "zlanhf.f"
		    }
#line 768 "zlanhf.f"
		    value = work[0];
#line 769 "zlanhf.f"
		    i__1 = *n - 1;
#line 769 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 770 "zlanhf.f"
			temp = work[i__];
#line 771 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 771 "zlanhf.f"
			    value = temp;
#line 771 "zlanhf.f"
			}
#line 773 "zlanhf.f"
		    }
#line 774 "zlanhf.f"
		}
#line 775 "zlanhf.f"
	    } else {
/*              n is even & A is n+1 by k = n/2 */
#line 777 "zlanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 779 "zlanhf.f"
		    i__1 = k - 1;
#line 779 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 780 "zlanhf.f"
			work[i__] = 0.;
#line 781 "zlanhf.f"
		    }
#line 782 "zlanhf.f"
		    i__1 = k - 1;
#line 782 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 783 "zlanhf.f"
			s = 0.;
#line 784 "zlanhf.f"
			i__2 = k + j - 1;
#line 784 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 785 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(i,j+k) */
#line 787 "zlanhf.f"
			    s += aa;
#line 788 "zlanhf.f"
			    work[i__] += aa;
#line 789 "zlanhf.f"
			}
#line 790 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 790 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 792 "zlanhf.f"
			work[j + k] = s + aa;
#line 793 "zlanhf.f"
			++i__;
#line 794 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 794 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j,j) */
#line 796 "zlanhf.f"
			work[j] += aa;
#line 797 "zlanhf.f"
			s = 0.;
#line 798 "zlanhf.f"
			i__2 = k - 1;
#line 798 "zlanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 799 "zlanhf.f"
			    ++i__;
#line 800 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 802 "zlanhf.f"
			    s += aa;
#line 803 "zlanhf.f"
			    work[l] += aa;
#line 804 "zlanhf.f"
			}
#line 805 "zlanhf.f"
			work[j] += s;
#line 806 "zlanhf.f"
		    }
#line 807 "zlanhf.f"
		    value = work[0];
#line 808 "zlanhf.f"
		    i__1 = *n - 1;
#line 808 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 809 "zlanhf.f"
			temp = work[i__];
#line 810 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 810 "zlanhf.f"
			    value = temp;
#line 810 "zlanhf.f"
			}
#line 812 "zlanhf.f"
		    }
#line 813 "zlanhf.f"
		} else {
/*                 ilu = 1 & uplo = 'L' */
#line 815 "zlanhf.f"
		    i__1 = *n - 1;
#line 815 "zlanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 816 "zlanhf.f"
			work[i__] = 0.;
#line 817 "zlanhf.f"
		    }
#line 818 "zlanhf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 819 "zlanhf.f"
			s = 0.;
#line 820 "zlanhf.f"
			i__1 = j - 1;
#line 820 "zlanhf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 821 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(j+k,i+k) */
#line 823 "zlanhf.f"
			    s += aa;
#line 824 "zlanhf.f"
			    work[i__ + k] += aa;
#line 825 "zlanhf.f"
			}
#line 826 "zlanhf.f"
			i__1 = i__ + j * lda;
#line 826 "zlanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 828 "zlanhf.f"
			s += aa;
#line 829 "zlanhf.f"
			work[i__ + k] += s;
/*                    i=j */
#line 831 "zlanhf.f"
			++i__;
#line 832 "zlanhf.f"
			i__1 = i__ + j * lda;
#line 832 "zlanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j,j) */
#line 834 "zlanhf.f"
			work[j] = aa;
#line 835 "zlanhf.f"
			s = 0.;
#line 836 "zlanhf.f"
			i__1 = *n - 1;
#line 836 "zlanhf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 837 "zlanhf.f"
			    ++i__;
#line 838 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 840 "zlanhf.f"
			    s += aa;
#line 841 "zlanhf.f"
			    work[l] += aa;
#line 842 "zlanhf.f"
			}
#line 843 "zlanhf.f"
			work[j] += s;
#line 844 "zlanhf.f"
		    }
#line 845 "zlanhf.f"
		    value = work[0];
#line 846 "zlanhf.f"
		    i__1 = *n - 1;
#line 846 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 847 "zlanhf.f"
			temp = work[i__];
#line 848 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 848 "zlanhf.f"
			    value = temp;
#line 848 "zlanhf.f"
			}
#line 850 "zlanhf.f"
		    }
#line 851 "zlanhf.f"
		}
#line 852 "zlanhf.f"
	    }
#line 853 "zlanhf.f"
	} else {
/*           ifm=0 */
#line 855 "zlanhf.f"
	    k = *n / 2;
#line 856 "zlanhf.f"
	    if (noe == 1) {
/*              n is odd & A is (n+1)/2 by n */
#line 858 "zlanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 860 "zlanhf.f"
		    n1 = k;
/*                 n/2 */
#line 862 "zlanhf.f"
		    ++k;
/*                 k is the row size and lda */
#line 864 "zlanhf.f"
		    i__1 = *n - 1;
#line 864 "zlanhf.f"
		    for (i__ = n1; i__ <= i__1; ++i__) {
#line 865 "zlanhf.f"
			work[i__] = 0.;
#line 866 "zlanhf.f"
		    }
#line 867 "zlanhf.f"
		    i__1 = n1 - 1;
#line 867 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 868 "zlanhf.f"
			s = 0.;
#line 869 "zlanhf.f"
			i__2 = k - 1;
#line 869 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 870 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,n1+i) */
#line 872 "zlanhf.f"
			    work[i__ + n1] += aa;
#line 873 "zlanhf.f"
			    s += aa;
#line 874 "zlanhf.f"
			}
#line 875 "zlanhf.f"
			work[j] = s;
#line 876 "zlanhf.f"
		    }
/*                 j=n1=k-1 is special */
#line 878 "zlanhf.f"
		    i__1 = j * lda;
#line 878 "zlanhf.f"
		    s = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 880 "zlanhf.f"
		    i__1 = k - 1;
#line 880 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 881 "zlanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k-1,i+n1) */
#line 883 "zlanhf.f"
			work[i__ + n1] += aa;
#line 884 "zlanhf.f"
			s += aa;
#line 885 "zlanhf.f"
		    }
#line 886 "zlanhf.f"
		    work[j] += s;
#line 887 "zlanhf.f"
		    i__1 = *n - 1;
#line 887 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 888 "zlanhf.f"
			s = 0.;
#line 889 "zlanhf.f"
			i__2 = j - k - 1;
#line 889 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 890 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(i,j-k) */
#line 892 "zlanhf.f"
			    work[i__] += aa;
#line 893 "zlanhf.f"
			    s += aa;
#line 894 "zlanhf.f"
			}
/*                    i=j-k */
#line 896 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 896 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j-k,j-k) */
#line 898 "zlanhf.f"
			s += aa;
#line 899 "zlanhf.f"
			work[j - k] += s;
#line 900 "zlanhf.f"
			++i__;
#line 901 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 901 "zlanhf.f"
			s = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j,j) */
#line 903 "zlanhf.f"
			i__2 = *n - 1;
#line 903 "zlanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 904 "zlanhf.f"
			    ++i__;
#line 905 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,l) */
#line 907 "zlanhf.f"
			    work[l] += aa;
#line 908 "zlanhf.f"
			    s += aa;
#line 909 "zlanhf.f"
			}
#line 910 "zlanhf.f"
			work[j] += s;
#line 911 "zlanhf.f"
		    }
#line 912 "zlanhf.f"
		    value = work[0];
#line 913 "zlanhf.f"
		    i__1 = *n - 1;
#line 913 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 914 "zlanhf.f"
			temp = work[i__];
#line 915 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 915 "zlanhf.f"
			    value = temp;
#line 915 "zlanhf.f"
			}
#line 917 "zlanhf.f"
		    }
#line 918 "zlanhf.f"
		} else {
/*                 ilu=1 & uplo = 'L' */
#line 920 "zlanhf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 922 "zlanhf.f"
		    i__1 = *n - 1;
#line 922 "zlanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 923 "zlanhf.f"
			work[i__] = 0.;
#line 924 "zlanhf.f"
		    }
#line 925 "zlanhf.f"
		    i__1 = k - 2;
#line 925 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
/*                    process */
#line 927 "zlanhf.f"
			s = 0.;
#line 928 "zlanhf.f"
			i__2 = j - 1;
#line 928 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 929 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i) */
#line 931 "zlanhf.f"
			    work[i__] += aa;
#line 932 "zlanhf.f"
			    s += aa;
#line 933 "zlanhf.f"
			}
#line 934 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 934 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    i=j so process of A(j,j) */
#line 936 "zlanhf.f"
			s += aa;
#line 937 "zlanhf.f"
			work[j] = s;
/*                    is initialised here */
#line 939 "zlanhf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 941 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 941 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
#line 942 "zlanhf.f"
			s = aa;
#line 943 "zlanhf.f"
			i__2 = *n - 1;
#line 943 "zlanhf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 944 "zlanhf.f"
			    ++i__;
#line 945 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(l,k+j) */
#line 947 "zlanhf.f"
			    s += aa;
#line 948 "zlanhf.f"
			    work[l] += aa;
#line 949 "zlanhf.f"
			}
#line 950 "zlanhf.f"
			work[k + j] += s;
#line 951 "zlanhf.f"
		    }
/*                 j=k-1 is special :process col A(k-1,0:k-1) */
#line 953 "zlanhf.f"
		    s = 0.;
#line 954 "zlanhf.f"
		    i__1 = k - 2;
#line 954 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 955 "zlanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,i) */
#line 957 "zlanhf.f"
			work[i__] += aa;
#line 958 "zlanhf.f"
			s += aa;
#line 959 "zlanhf.f"
		    }
/*                 i=k-1 */
#line 961 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 961 "zlanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 963 "zlanhf.f"
		    s += aa;
#line 964 "zlanhf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 966 "zlanhf.f"
		    i__1 = *n - 1;
#line 966 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
/*                    process col j of A = A(j,0:k-1) */
#line 968 "zlanhf.f"
			s = 0.;
#line 969 "zlanhf.f"
			i__2 = k - 1;
#line 969 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 970 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i) */
#line 972 "zlanhf.f"
			    work[i__] += aa;
#line 973 "zlanhf.f"
			    s += aa;
#line 974 "zlanhf.f"
			}
#line 975 "zlanhf.f"
			work[j] += s;
#line 976 "zlanhf.f"
		    }
#line 977 "zlanhf.f"
		    value = work[0];
#line 978 "zlanhf.f"
		    i__1 = *n - 1;
#line 978 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 979 "zlanhf.f"
			temp = work[i__];
#line 980 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 980 "zlanhf.f"
			    value = temp;
#line 980 "zlanhf.f"
			}
#line 982 "zlanhf.f"
		    }
#line 983 "zlanhf.f"
		}
#line 984 "zlanhf.f"
	    } else {
/*              n is even & A is k=n/2 by n+1 */
#line 986 "zlanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 988 "zlanhf.f"
		    i__1 = *n - 1;
#line 988 "zlanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 989 "zlanhf.f"
			work[i__] = 0.;
#line 990 "zlanhf.f"
		    }
#line 991 "zlanhf.f"
		    i__1 = k - 1;
#line 991 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 992 "zlanhf.f"
			s = 0.;
#line 993 "zlanhf.f"
			i__2 = k - 1;
#line 993 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 994 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i+k) */
#line 996 "zlanhf.f"
			    work[i__ + k] += aa;
#line 997 "zlanhf.f"
			    s += aa;
#line 998 "zlanhf.f"
			}
#line 999 "zlanhf.f"
			work[j] = s;
#line 1000 "zlanhf.f"
		    }
/*                 j=k */
#line 1002 "zlanhf.f"
		    i__1 = j * lda;
#line 1002 "zlanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k,k) */
#line 1004 "zlanhf.f"
		    s = aa;
#line 1005 "zlanhf.f"
		    i__1 = k - 1;
#line 1005 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1006 "zlanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,k+i) */
#line 1008 "zlanhf.f"
			work[i__ + k] += aa;
#line 1009 "zlanhf.f"
			s += aa;
#line 1010 "zlanhf.f"
		    }
#line 1011 "zlanhf.f"
		    work[j] += s;
#line 1012 "zlanhf.f"
		    i__1 = *n - 1;
#line 1012 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1013 "zlanhf.f"
			s = 0.;
#line 1014 "zlanhf.f"
			i__2 = j - 2 - k;
#line 1014 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1015 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(i,j-k-1) */
#line 1017 "zlanhf.f"
			    work[i__] += aa;
#line 1018 "zlanhf.f"
			    s += aa;
#line 1019 "zlanhf.f"
			}
/*                    i=j-1-k */
#line 1021 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 1021 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j-k-1,j-k-1) */
#line 1023 "zlanhf.f"
			s += aa;
#line 1024 "zlanhf.f"
			work[j - k - 1] += s;
#line 1025 "zlanhf.f"
			++i__;
#line 1026 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 1026 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j,j) */
#line 1028 "zlanhf.f"
			s = aa;
#line 1029 "zlanhf.f"
			i__2 = *n - 1;
#line 1029 "zlanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 1030 "zlanhf.f"
			    ++i__;
#line 1031 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,l) */
#line 1033 "zlanhf.f"
			    work[l] += aa;
#line 1034 "zlanhf.f"
			    s += aa;
#line 1035 "zlanhf.f"
			}
#line 1036 "zlanhf.f"
			work[j] += s;
#line 1037 "zlanhf.f"
		    }
/*                 j=n */
#line 1039 "zlanhf.f"
		    s = 0.;
#line 1040 "zlanhf.f"
		    i__1 = k - 2;
#line 1040 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1041 "zlanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(i,k-1) */
#line 1043 "zlanhf.f"
			work[i__] += aa;
#line 1044 "zlanhf.f"
			s += aa;
#line 1045 "zlanhf.f"
		    }
/*                 i=k-1 */
#line 1047 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 1047 "zlanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 1049 "zlanhf.f"
		    s += aa;
#line 1050 "zlanhf.f"
		    work[i__] += s;
#line 1051 "zlanhf.f"
		    value = work[0];
#line 1052 "zlanhf.f"
		    i__1 = *n - 1;
#line 1052 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1053 "zlanhf.f"
			temp = work[i__];
#line 1054 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 1054 "zlanhf.f"
			    value = temp;
#line 1054 "zlanhf.f"
			}
#line 1056 "zlanhf.f"
		    }
#line 1057 "zlanhf.f"
		} else {
/*                 ilu=1 & uplo = 'L' */
#line 1059 "zlanhf.f"
		    i__1 = *n - 1;
#line 1059 "zlanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 1060 "zlanhf.f"
			work[i__] = 0.;
#line 1061 "zlanhf.f"
		    }
/*                 j=0 is special :process col A(k:n-1,k) */
#line 1063 "zlanhf.f"
		    s = (d__1 = a[0].r, abs(d__1));
/*                 A(k,k) */
#line 1065 "zlanhf.f"
		    i__1 = k - 1;
#line 1065 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1066 "zlanhf.f"
			aa = z_abs(&a[i__]);
/*                    A(k+i,k) */
#line 1068 "zlanhf.f"
			work[i__ + k] += aa;
#line 1069 "zlanhf.f"
			s += aa;
#line 1070 "zlanhf.f"
		    }
#line 1071 "zlanhf.f"
		    work[k] += s;
#line 1072 "zlanhf.f"
		    i__1 = k - 1;
#line 1072 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
/*                    process */
#line 1074 "zlanhf.f"
			s = 0.;
#line 1075 "zlanhf.f"
			i__2 = j - 2;
#line 1075 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1076 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j-1,i) */
#line 1078 "zlanhf.f"
			    work[i__] += aa;
#line 1079 "zlanhf.f"
			    s += aa;
#line 1080 "zlanhf.f"
			}
#line 1081 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 1081 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    i=j-1 so process of A(j-1,j-1) */
#line 1083 "zlanhf.f"
			s += aa;
#line 1084 "zlanhf.f"
			work[j - 1] = s;
/*                    is initialised here */
#line 1086 "zlanhf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 1088 "zlanhf.f"
			i__2 = i__ + j * lda;
#line 1088 "zlanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
#line 1089 "zlanhf.f"
			s = aa;
#line 1090 "zlanhf.f"
			i__2 = *n - 1;
#line 1090 "zlanhf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 1091 "zlanhf.f"
			    ++i__;
#line 1092 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(l,k+j) */
#line 1094 "zlanhf.f"
			    s += aa;
#line 1095 "zlanhf.f"
			    work[l] += aa;
#line 1096 "zlanhf.f"
			}
#line 1097 "zlanhf.f"
			work[k + j] += s;
#line 1098 "zlanhf.f"
		    }
/*                 j=k is special :process col A(k,0:k-1) */
#line 1100 "zlanhf.f"
		    s = 0.;
#line 1101 "zlanhf.f"
		    i__1 = k - 2;
#line 1101 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1102 "zlanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,i) */
#line 1104 "zlanhf.f"
			work[i__] += aa;
#line 1105 "zlanhf.f"
			s += aa;
#line 1106 "zlanhf.f"
		    }

/*                 i=k-1 */
#line 1109 "zlanhf.f"
		    i__1 = i__ + j * lda;
#line 1109 "zlanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 1111 "zlanhf.f"
		    s += aa;
#line 1112 "zlanhf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 1114 "zlanhf.f"
		    i__1 = *n;
#line 1114 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {

/*                    process col j-1 of A = A(j-1,0:k-1) */
#line 1117 "zlanhf.f"
			s = 0.;
#line 1118 "zlanhf.f"
			i__2 = k - 1;
#line 1118 "zlanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1119 "zlanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j-1,i) */
#line 1121 "zlanhf.f"
			    work[i__] += aa;
#line 1122 "zlanhf.f"
			    s += aa;
#line 1123 "zlanhf.f"
			}
#line 1124 "zlanhf.f"
			work[j - 1] += s;
#line 1125 "zlanhf.f"
		    }
#line 1126 "zlanhf.f"
		    value = work[0];
#line 1127 "zlanhf.f"
		    i__1 = *n - 1;
#line 1127 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1128 "zlanhf.f"
			temp = work[i__];
#line 1129 "zlanhf.f"
			if (value < temp || disnan_(&temp)) {
#line 1129 "zlanhf.f"
			    value = temp;
#line 1129 "zlanhf.f"
			}
#line 1131 "zlanhf.f"
		    }
#line 1132 "zlanhf.f"
		}
#line 1133 "zlanhf.f"
	    }
#line 1134 "zlanhf.f"
	}
#line 1135 "zlanhf.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*       Find normF(A). */

#line 1139 "zlanhf.f"
	k = (*n + 1) / 2;
#line 1140 "zlanhf.f"
	scale = 0.;
#line 1141 "zlanhf.f"
	s = 1.;
#line 1142 "zlanhf.f"
	if (noe == 1) {
/*           n is odd */
#line 1144 "zlanhf.f"
	    if (ifm == 1) {
/*              A is normal & A is n by k */
#line 1146 "zlanhf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 1148 "zlanhf.f"
		    i__1 = k - 3;
#line 1148 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1149 "zlanhf.f"
			i__2 = k - j - 2;
#line 1149 "zlanhf.f"
			zlassq_(&i__2, &a[k + j + 1 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k,0) */
#line 1151 "zlanhf.f"
		    }
#line 1152 "zlanhf.f"
		    i__1 = k - 1;
#line 1152 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1153 "zlanhf.f"
			i__2 = k + j - 1;
#line 1153 "zlanhf.f"
			zlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 1155 "zlanhf.f"
		    }
#line 1156 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1158 "zlanhf.f"
		    l = k - 1;
/*                 -> U(k,k) at A(k-1,0) */
#line 1160 "zlanhf.f"
		    i__1 = k - 2;
#line 1160 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1161 "zlanhf.f"
			i__2 = l;
#line 1161 "zlanhf.f"
			aa = a[i__2].r;
/*                    U(k+i,k+i) */
#line 1163 "zlanhf.f"
			if (aa != 0.) {
#line 1164 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1165 "zlanhf.f"
				d__1 = scale / aa;
#line 1165 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1166 "zlanhf.f"
				scale = aa;
#line 1167 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1168 "zlanhf.f"
				d__1 = aa / scale;
#line 1168 "zlanhf.f"
				s += d__1 * d__1;
#line 1169 "zlanhf.f"
			    }
#line 1170 "zlanhf.f"
			}
#line 1171 "zlanhf.f"
			i__2 = l + 1;
#line 1171 "zlanhf.f"
			aa = a[i__2].r;
/*                    U(i,i) */
#line 1173 "zlanhf.f"
			if (aa != 0.) {
#line 1174 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1175 "zlanhf.f"
				d__1 = scale / aa;
#line 1175 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1176 "zlanhf.f"
				scale = aa;
#line 1177 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1178 "zlanhf.f"
				d__1 = aa / scale;
#line 1178 "zlanhf.f"
				s += d__1 * d__1;
#line 1179 "zlanhf.f"
			    }
#line 1180 "zlanhf.f"
			}
#line 1181 "zlanhf.f"
			l = l + lda + 1;
#line 1182 "zlanhf.f"
		    }
#line 1183 "zlanhf.f"
		    i__1 = l;
#line 1183 "zlanhf.f"
		    aa = a[i__1].r;
/*                 U(n-1,n-1) */
#line 1185 "zlanhf.f"
		    if (aa != 0.) {
#line 1186 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1187 "zlanhf.f"
			    d__1 = scale / aa;
#line 1187 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1188 "zlanhf.f"
			    scale = aa;
#line 1189 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1190 "zlanhf.f"
			    d__1 = aa / scale;
#line 1190 "zlanhf.f"
			    s += d__1 * d__1;
#line 1191 "zlanhf.f"
			}
#line 1192 "zlanhf.f"
		    }
#line 1193 "zlanhf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 1195 "zlanhf.f"
		    i__1 = k - 1;
#line 1195 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1196 "zlanhf.f"
			i__2 = *n - j - 1;
#line 1196 "zlanhf.f"
			zlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(0,0) */
#line 1198 "zlanhf.f"
		    }
#line 1199 "zlanhf.f"
		    i__1 = k - 2;
#line 1199 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1200 "zlanhf.f"
			zlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 1202 "zlanhf.f"
		    }
#line 1203 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1205 "zlanhf.f"
		    aa = a[0].r;
/*                 L(0,0) at A(0,0) */
#line 1207 "zlanhf.f"
		    if (aa != 0.) {
#line 1208 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1209 "zlanhf.f"
			    d__1 = scale / aa;
#line 1209 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1210 "zlanhf.f"
			    scale = aa;
#line 1211 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1212 "zlanhf.f"
			    d__1 = aa / scale;
#line 1212 "zlanhf.f"
			    s += d__1 * d__1;
#line 1213 "zlanhf.f"
			}
#line 1214 "zlanhf.f"
		    }
#line 1215 "zlanhf.f"
		    l = lda;
/*                 -> L(k,k) at A(0,1) */
#line 1217 "zlanhf.f"
		    i__1 = k - 1;
#line 1217 "zlanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1218 "zlanhf.f"
			i__2 = l;
#line 1218 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(k-1+i,k-1+i) */
#line 1220 "zlanhf.f"
			if (aa != 0.) {
#line 1221 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1222 "zlanhf.f"
				d__1 = scale / aa;
#line 1222 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1223 "zlanhf.f"
				scale = aa;
#line 1224 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1225 "zlanhf.f"
				d__1 = aa / scale;
#line 1225 "zlanhf.f"
				s += d__1 * d__1;
#line 1226 "zlanhf.f"
			    }
#line 1227 "zlanhf.f"
			}
#line 1228 "zlanhf.f"
			i__2 = l + 1;
#line 1228 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1230 "zlanhf.f"
			if (aa != 0.) {
#line 1231 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1232 "zlanhf.f"
				d__1 = scale / aa;
#line 1232 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1233 "zlanhf.f"
				scale = aa;
#line 1234 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1235 "zlanhf.f"
				d__1 = aa / scale;
#line 1235 "zlanhf.f"
				s += d__1 * d__1;
#line 1236 "zlanhf.f"
			    }
#line 1237 "zlanhf.f"
			}
#line 1238 "zlanhf.f"
			l = l + lda + 1;
#line 1239 "zlanhf.f"
		    }
#line 1240 "zlanhf.f"
		}
#line 1241 "zlanhf.f"
	    } else {
/*              A is xpose & A is k by n */
#line 1243 "zlanhf.f"
		if (ilu == 0) {
/*                 A**H is upper */
#line 1245 "zlanhf.f"
		    i__1 = k - 2;
#line 1245 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1246 "zlanhf.f"
			zlassq_(&j, &a[(k + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k) */
#line 1248 "zlanhf.f"
		    }
#line 1249 "zlanhf.f"
		    i__1 = k - 2;
#line 1249 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1250 "zlanhf.f"
			zlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,0) */
#line 1252 "zlanhf.f"
		    }
#line 1253 "zlanhf.f"
		    i__1 = k - 2;
#line 1253 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1254 "zlanhf.f"
			i__2 = k - j - 1;
#line 1254 "zlanhf.f"
			zlassq_(&i__2, &a[j + 1 + (j + k - 1) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k-1) */
#line 1257 "zlanhf.f"
		    }
#line 1258 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1260 "zlanhf.f"
		    l = k * lda - lda;
/*                 -> U(k-1,k-1) at A(0,k-1) */
#line 1262 "zlanhf.f"
		    i__1 = l;
#line 1262 "zlanhf.f"
		    aa = a[i__1].r;
/*                 U(k-1,k-1) */
#line 1264 "zlanhf.f"
		    if (aa != 0.) {
#line 1265 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1266 "zlanhf.f"
			    d__1 = scale / aa;
#line 1266 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1267 "zlanhf.f"
			    scale = aa;
#line 1268 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1269 "zlanhf.f"
			    d__1 = aa / scale;
#line 1269 "zlanhf.f"
			    s += d__1 * d__1;
#line 1270 "zlanhf.f"
			}
#line 1271 "zlanhf.f"
		    }
#line 1272 "zlanhf.f"
		    l += lda;
/*                 -> U(0,0) at A(0,k) */
#line 1274 "zlanhf.f"
		    i__1 = *n - 1;
#line 1274 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 1275 "zlanhf.f"
			i__2 = l;
#line 1275 "zlanhf.f"
			aa = a[i__2].r;
/*                    -> U(j-k,j-k) */
#line 1277 "zlanhf.f"
			if (aa != 0.) {
#line 1278 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1279 "zlanhf.f"
				d__1 = scale / aa;
#line 1279 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1280 "zlanhf.f"
				scale = aa;
#line 1281 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1282 "zlanhf.f"
				d__1 = aa / scale;
#line 1282 "zlanhf.f"
				s += d__1 * d__1;
#line 1283 "zlanhf.f"
			    }
#line 1284 "zlanhf.f"
			}
#line 1285 "zlanhf.f"
			i__2 = l + 1;
#line 1285 "zlanhf.f"
			aa = a[i__2].r;
/*                    -> U(j,j) */
#line 1287 "zlanhf.f"
			if (aa != 0.) {
#line 1288 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1289 "zlanhf.f"
				d__1 = scale / aa;
#line 1289 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1290 "zlanhf.f"
				scale = aa;
#line 1291 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1292 "zlanhf.f"
				d__1 = aa / scale;
#line 1292 "zlanhf.f"
				s += d__1 * d__1;
#line 1293 "zlanhf.f"
			    }
#line 1294 "zlanhf.f"
			}
#line 1295 "zlanhf.f"
			l = l + lda + 1;
#line 1296 "zlanhf.f"
		    }
#line 1297 "zlanhf.f"
		} else {
/*                 A**H is lower */
#line 1299 "zlanhf.f"
		    i__1 = k - 1;
#line 1299 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1300 "zlanhf.f"
			zlassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 1302 "zlanhf.f"
		    }
#line 1303 "zlanhf.f"
		    i__1 = *n - 1;
#line 1303 "zlanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 1304 "zlanhf.f"
			zlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,k) */
#line 1306 "zlanhf.f"
		    }
#line 1307 "zlanhf.f"
		    i__1 = k - 3;
#line 1307 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1308 "zlanhf.f"
			i__2 = k - j - 2;
#line 1308 "zlanhf.f"
			zlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(1,0) */
#line 1310 "zlanhf.f"
		    }
#line 1311 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1313 "zlanhf.f"
		    l = 0;
/*                 -> L(0,0) at A(0,0) */
#line 1315 "zlanhf.f"
		    i__1 = k - 2;
#line 1315 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1316 "zlanhf.f"
			i__2 = l;
#line 1316 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1318 "zlanhf.f"
			if (aa != 0.) {
#line 1319 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1320 "zlanhf.f"
				d__1 = scale / aa;
#line 1320 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1321 "zlanhf.f"
				scale = aa;
#line 1322 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1323 "zlanhf.f"
				d__1 = aa / scale;
#line 1323 "zlanhf.f"
				s += d__1 * d__1;
#line 1324 "zlanhf.f"
			    }
#line 1325 "zlanhf.f"
			}
#line 1326 "zlanhf.f"
			i__2 = l + 1;
#line 1326 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(k+i,k+i) */
#line 1328 "zlanhf.f"
			if (aa != 0.) {
#line 1329 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1330 "zlanhf.f"
				d__1 = scale / aa;
#line 1330 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1331 "zlanhf.f"
				scale = aa;
#line 1332 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1333 "zlanhf.f"
				d__1 = aa / scale;
#line 1333 "zlanhf.f"
				s += d__1 * d__1;
#line 1334 "zlanhf.f"
			    }
#line 1335 "zlanhf.f"
			}
#line 1336 "zlanhf.f"
			l = l + lda + 1;
#line 1337 "zlanhf.f"
		    }
/*                 L-> k-1 + (k-1)*lda or L(k-1,k-1) at A(k-1,k-1) */
#line 1339 "zlanhf.f"
		    i__1 = l;
#line 1339 "zlanhf.f"
		    aa = a[i__1].r;
/*                 L(k-1,k-1) at A(k-1,k-1) */
#line 1341 "zlanhf.f"
		    if (aa != 0.) {
#line 1342 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1343 "zlanhf.f"
			    d__1 = scale / aa;
#line 1343 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1344 "zlanhf.f"
			    scale = aa;
#line 1345 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1346 "zlanhf.f"
			    d__1 = aa / scale;
#line 1346 "zlanhf.f"
			    s += d__1 * d__1;
#line 1347 "zlanhf.f"
			}
#line 1348 "zlanhf.f"
		    }
#line 1349 "zlanhf.f"
		}
#line 1350 "zlanhf.f"
	    }
#line 1351 "zlanhf.f"
	} else {
/*           n is even */
#line 1353 "zlanhf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 1355 "zlanhf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 1357 "zlanhf.f"
		    i__1 = k - 2;
#line 1357 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1358 "zlanhf.f"
			i__2 = k - j - 1;
#line 1358 "zlanhf.f"
			zlassq_(&i__2, &a[k + j + 2 + j * lda], &c__1, &scale,
				 &s);
/*                 L at A(k+1,0) */
#line 1360 "zlanhf.f"
		    }
#line 1361 "zlanhf.f"
		    i__1 = k - 1;
#line 1361 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1362 "zlanhf.f"
			i__2 = k + j;
#line 1362 "zlanhf.f"
			zlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                 trap U at A(0,0) */
#line 1364 "zlanhf.f"
		    }
#line 1365 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1367 "zlanhf.f"
		    l = k;
/*                 -> U(k,k) at A(k,0) */
#line 1369 "zlanhf.f"
		    i__1 = k - 1;
#line 1369 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1370 "zlanhf.f"
			i__2 = l;
#line 1370 "zlanhf.f"
			aa = a[i__2].r;
/*                    U(k+i,k+i) */
#line 1372 "zlanhf.f"
			if (aa != 0.) {
#line 1373 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1374 "zlanhf.f"
				d__1 = scale / aa;
#line 1374 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1375 "zlanhf.f"
				scale = aa;
#line 1376 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1377 "zlanhf.f"
				d__1 = aa / scale;
#line 1377 "zlanhf.f"
				s += d__1 * d__1;
#line 1378 "zlanhf.f"
			    }
#line 1379 "zlanhf.f"
			}
#line 1380 "zlanhf.f"
			i__2 = l + 1;
#line 1380 "zlanhf.f"
			aa = a[i__2].r;
/*                    U(i,i) */
#line 1382 "zlanhf.f"
			if (aa != 0.) {
#line 1383 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1384 "zlanhf.f"
				d__1 = scale / aa;
#line 1384 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1385 "zlanhf.f"
				scale = aa;
#line 1386 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1387 "zlanhf.f"
				d__1 = aa / scale;
#line 1387 "zlanhf.f"
				s += d__1 * d__1;
#line 1388 "zlanhf.f"
			    }
#line 1389 "zlanhf.f"
			}
#line 1390 "zlanhf.f"
			l = l + lda + 1;
#line 1391 "zlanhf.f"
		    }
#line 1392 "zlanhf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 1394 "zlanhf.f"
		    i__1 = k - 1;
#line 1394 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1395 "zlanhf.f"
			i__2 = *n - j - 1;
#line 1395 "zlanhf.f"
			zlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(1,0) */
#line 1397 "zlanhf.f"
		    }
#line 1398 "zlanhf.f"
		    i__1 = k - 1;
#line 1398 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1399 "zlanhf.f"
			zlassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 1401 "zlanhf.f"
		    }
#line 1402 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1404 "zlanhf.f"
		    l = 0;
/*                 -> L(k,k) at A(0,0) */
#line 1406 "zlanhf.f"
		    i__1 = k - 1;
#line 1406 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1407 "zlanhf.f"
			i__2 = l;
#line 1407 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(k-1+i,k-1+i) */
#line 1409 "zlanhf.f"
			if (aa != 0.) {
#line 1410 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1411 "zlanhf.f"
				d__1 = scale / aa;
#line 1411 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1412 "zlanhf.f"
				scale = aa;
#line 1413 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1414 "zlanhf.f"
				d__1 = aa / scale;
#line 1414 "zlanhf.f"
				s += d__1 * d__1;
#line 1415 "zlanhf.f"
			    }
#line 1416 "zlanhf.f"
			}
#line 1417 "zlanhf.f"
			i__2 = l + 1;
#line 1417 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1419 "zlanhf.f"
			if (aa != 0.) {
#line 1420 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1421 "zlanhf.f"
				d__1 = scale / aa;
#line 1421 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1422 "zlanhf.f"
				scale = aa;
#line 1423 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1424 "zlanhf.f"
				d__1 = aa / scale;
#line 1424 "zlanhf.f"
				s += d__1 * d__1;
#line 1425 "zlanhf.f"
			    }
#line 1426 "zlanhf.f"
			}
#line 1427 "zlanhf.f"
			l = l + lda + 1;
#line 1428 "zlanhf.f"
		    }
#line 1429 "zlanhf.f"
		}
#line 1430 "zlanhf.f"
	    } else {
/*              A is xpose */
#line 1432 "zlanhf.f"
		if (ilu == 0) {
/*                 A**H is upper */
#line 1434 "zlanhf.f"
		    i__1 = k - 1;
#line 1434 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1435 "zlanhf.f"
			zlassq_(&j, &a[(k + 1 + j) * lda], &c__1, &scale, &s);
/*                 U at A(0,k+1) */
#line 1437 "zlanhf.f"
		    }
#line 1438 "zlanhf.f"
		    i__1 = k - 1;
#line 1438 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1439 "zlanhf.f"
			zlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                 k by k rect. at A(0,0) */
#line 1441 "zlanhf.f"
		    }
#line 1442 "zlanhf.f"
		    i__1 = k - 2;
#line 1442 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1443 "zlanhf.f"
			i__2 = k - j - 1;
#line 1443 "zlanhf.f"
			zlassq_(&i__2, &a[j + 1 + (j + k) * lda], &c__1, &
				scale, &s);
/*                 L at A(0,k) */
#line 1446 "zlanhf.f"
		    }
#line 1447 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1449 "zlanhf.f"
		    l = k * lda;
/*                 -> U(k,k) at A(0,k) */
#line 1451 "zlanhf.f"
		    i__1 = l;
#line 1451 "zlanhf.f"
		    aa = a[i__1].r;
/*                 U(k,k) */
#line 1453 "zlanhf.f"
		    if (aa != 0.) {
#line 1454 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1455 "zlanhf.f"
			    d__1 = scale / aa;
#line 1455 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1456 "zlanhf.f"
			    scale = aa;
#line 1457 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1458 "zlanhf.f"
			    d__1 = aa / scale;
#line 1458 "zlanhf.f"
			    s += d__1 * d__1;
#line 1459 "zlanhf.f"
			}
#line 1460 "zlanhf.f"
		    }
#line 1461 "zlanhf.f"
		    l += lda;
/*                 -> U(0,0) at A(0,k+1) */
#line 1463 "zlanhf.f"
		    i__1 = *n - 1;
#line 1463 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1464 "zlanhf.f"
			i__2 = l;
#line 1464 "zlanhf.f"
			aa = a[i__2].r;
/*                    -> U(j-k-1,j-k-1) */
#line 1466 "zlanhf.f"
			if (aa != 0.) {
#line 1467 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1468 "zlanhf.f"
				d__1 = scale / aa;
#line 1468 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1469 "zlanhf.f"
				scale = aa;
#line 1470 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1471 "zlanhf.f"
				d__1 = aa / scale;
#line 1471 "zlanhf.f"
				s += d__1 * d__1;
#line 1472 "zlanhf.f"
			    }
#line 1473 "zlanhf.f"
			}
#line 1474 "zlanhf.f"
			i__2 = l + 1;
#line 1474 "zlanhf.f"
			aa = a[i__2].r;
/*                    -> U(j,j) */
#line 1476 "zlanhf.f"
			if (aa != 0.) {
#line 1477 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1478 "zlanhf.f"
				d__1 = scale / aa;
#line 1478 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1479 "zlanhf.f"
				scale = aa;
#line 1480 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1481 "zlanhf.f"
				d__1 = aa / scale;
#line 1481 "zlanhf.f"
				s += d__1 * d__1;
#line 1482 "zlanhf.f"
			    }
#line 1483 "zlanhf.f"
			}
#line 1484 "zlanhf.f"
			l = l + lda + 1;
#line 1485 "zlanhf.f"
		    }
/*                 L=k-1+n*lda */
/*                 -> U(k-1,k-1) at A(k-1,n) */
#line 1488 "zlanhf.f"
		    i__1 = l;
#line 1488 "zlanhf.f"
		    aa = a[i__1].r;
/*                 U(k,k) */
#line 1490 "zlanhf.f"
		    if (aa != 0.) {
#line 1491 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1492 "zlanhf.f"
			    d__1 = scale / aa;
#line 1492 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1493 "zlanhf.f"
			    scale = aa;
#line 1494 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1495 "zlanhf.f"
			    d__1 = aa / scale;
#line 1495 "zlanhf.f"
			    s += d__1 * d__1;
#line 1496 "zlanhf.f"
			}
#line 1497 "zlanhf.f"
		    }
#line 1498 "zlanhf.f"
		} else {
/*                 A**H is lower */
#line 1500 "zlanhf.f"
		    i__1 = k - 1;
#line 1500 "zlanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1501 "zlanhf.f"
			zlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                 U at A(0,1) */
#line 1503 "zlanhf.f"
		    }
#line 1504 "zlanhf.f"
		    i__1 = *n;
#line 1504 "zlanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1505 "zlanhf.f"
			zlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                 k by k rect. at A(0,k+1) */
#line 1507 "zlanhf.f"
		    }
#line 1508 "zlanhf.f"
		    i__1 = k - 2;
#line 1508 "zlanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1509 "zlanhf.f"
			i__2 = k - j - 1;
#line 1509 "zlanhf.f"
			zlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                 L at A(0,0) */
#line 1511 "zlanhf.f"
		    }
#line 1512 "zlanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1514 "zlanhf.f"
		    l = 0;
/*                 -> L(k,k) at A(0,0) */
#line 1516 "zlanhf.f"
		    i__1 = l;
#line 1516 "zlanhf.f"
		    aa = a[i__1].r;
/*                 L(k,k) at A(0,0) */
#line 1518 "zlanhf.f"
		    if (aa != 0.) {
#line 1519 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1520 "zlanhf.f"
			    d__1 = scale / aa;
#line 1520 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1521 "zlanhf.f"
			    scale = aa;
#line 1522 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1523 "zlanhf.f"
			    d__1 = aa / scale;
#line 1523 "zlanhf.f"
			    s += d__1 * d__1;
#line 1524 "zlanhf.f"
			}
#line 1525 "zlanhf.f"
		    }
#line 1526 "zlanhf.f"
		    l = lda;
/*                 -> L(0,0) at A(0,1) */
#line 1528 "zlanhf.f"
		    i__1 = k - 2;
#line 1528 "zlanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1529 "zlanhf.f"
			i__2 = l;
#line 1529 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1531 "zlanhf.f"
			if (aa != 0.) {
#line 1532 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1533 "zlanhf.f"
				d__1 = scale / aa;
#line 1533 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1534 "zlanhf.f"
				scale = aa;
#line 1535 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1536 "zlanhf.f"
				d__1 = aa / scale;
#line 1536 "zlanhf.f"
				s += d__1 * d__1;
#line 1537 "zlanhf.f"
			    }
#line 1538 "zlanhf.f"
			}
#line 1539 "zlanhf.f"
			i__2 = l + 1;
#line 1539 "zlanhf.f"
			aa = a[i__2].r;
/*                    L(k+i+1,k+i+1) */
#line 1541 "zlanhf.f"
			if (aa != 0.) {
#line 1542 "zlanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1543 "zlanhf.f"
				d__1 = scale / aa;
#line 1543 "zlanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1544 "zlanhf.f"
				scale = aa;
#line 1545 "zlanhf.f"
			    } else {
/* Computing 2nd power */
#line 1546 "zlanhf.f"
				d__1 = aa / scale;
#line 1546 "zlanhf.f"
				s += d__1 * d__1;
#line 1547 "zlanhf.f"
			    }
#line 1548 "zlanhf.f"
			}
#line 1549 "zlanhf.f"
			l = l + lda + 1;
#line 1550 "zlanhf.f"
		    }
/*                 L-> k - 1 + k*lda or L(k-1,k-1) at A(k-1,k) */
#line 1552 "zlanhf.f"
		    i__1 = l;
#line 1552 "zlanhf.f"
		    aa = a[i__1].r;
/*                 L(k-1,k-1) at A(k-1,k) */
#line 1554 "zlanhf.f"
		    if (aa != 0.) {
#line 1555 "zlanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1556 "zlanhf.f"
			    d__1 = scale / aa;
#line 1556 "zlanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1557 "zlanhf.f"
			    scale = aa;
#line 1558 "zlanhf.f"
			} else {
/* Computing 2nd power */
#line 1559 "zlanhf.f"
			    d__1 = aa / scale;
#line 1559 "zlanhf.f"
			    s += d__1 * d__1;
#line 1560 "zlanhf.f"
			}
#line 1561 "zlanhf.f"
		    }
#line 1562 "zlanhf.f"
		}
#line 1563 "zlanhf.f"
	    }
#line 1564 "zlanhf.f"
	}
#line 1565 "zlanhf.f"
	value = scale * sqrt(s);
#line 1566 "zlanhf.f"
    }

#line 1568 "zlanhf.f"
    ret_val = value;
#line 1569 "zlanhf.f"
    return ret_val;

/*     End of ZLANHF */

} /* zlanhf_ */

