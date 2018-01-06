#line 1 "clanhf.f"
/* clanhf.f -- translated by f2c (version 20100827).
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

#line 1 "clanhf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANHF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a Hermitian matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANHF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION CLANHF( NORM, TRANSR, UPLO, N, A, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, TRANSR, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               WORK( 0: * ) */
/*       COMPLEX            A( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANHF  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex Hermitian matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \return CLANHF */
/* > \verbatim */
/* > */
/* >    CLANHF = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >            Specifies the value to be returned in CLANHF as described */
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
/* >            The order of the matrix A.  N >= 0.  When N = 0, CLANHF is */
/* >            set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension ( N*(N+1)/2 ); */
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
/* >          WORK is REAL array, dimension (LWORK), */
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
doublereal clanhf_(char *norm, char *transr, char *uplo, integer *n, 
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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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

#line 285 "clanhf.f"
    if (*n == 0) {
#line 286 "clanhf.f"
	ret_val = 0.;
#line 287 "clanhf.f"
	return ret_val;
#line 288 "clanhf.f"
    } else if (*n == 1) {
#line 289 "clanhf.f"
	ret_val = z_abs(a);
#line 290 "clanhf.f"
	return ret_val;
#line 291 "clanhf.f"
    }

/*     set noe = 1 if n is odd. if n is even set noe=0 */

#line 295 "clanhf.f"
    noe = 1;
#line 296 "clanhf.f"
    if (*n % 2 == 0) {
#line 296 "clanhf.f"
	noe = 0;
#line 296 "clanhf.f"
    }

/*     set ifm = 0 when form='C' or 'c' and 1 otherwise */

#line 301 "clanhf.f"
    ifm = 1;
#line 302 "clanhf.f"
    if (lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 302 "clanhf.f"
	ifm = 0;
#line 302 "clanhf.f"
    }

/*     set ilu = 0 when uplo='U or 'u' and 1 otherwise */

#line 307 "clanhf.f"
    ilu = 1;
#line 308 "clanhf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 308 "clanhf.f"
	ilu = 0;
#line 308 "clanhf.f"
    }

/*     set lda = (n+1)/2 when ifm = 0 */
/*     set lda = n when ifm = 1 and noe = 1 */
/*     set lda = n+1 when ifm = 1 and noe = 0 */

#line 315 "clanhf.f"
    if (ifm == 1) {
#line 316 "clanhf.f"
	if (noe == 1) {
#line 317 "clanhf.f"
	    lda = *n;
#line 318 "clanhf.f"
	} else {
/*           noe=0 */
#line 320 "clanhf.f"
	    lda = *n + 1;
#line 321 "clanhf.f"
	}
#line 322 "clanhf.f"
    } else {
/*        ifm=0 */
#line 324 "clanhf.f"
	lda = (*n + 1) / 2;
#line 325 "clanhf.f"
    }

#line 327 "clanhf.f"
    if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*       Find max(abs(A(i,j))). */

#line 331 "clanhf.f"
	k = (*n + 1) / 2;
#line 332 "clanhf.f"
	value = 0.;
#line 333 "clanhf.f"
	if (noe == 1) {
/*           n is odd & n = k + k - 1 */
#line 335 "clanhf.f"
	    if (ifm == 1) {
/*              A is n by k */
#line 337 "clanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 339 "clanhf.f"
		    j = 0;
/*                 -> L(0,0) */
#line 341 "clanhf.f"
		    i__1 = j + j * lda;
#line 341 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 342 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 342 "clanhf.f"
			value = temp;
#line 342 "clanhf.f"
		    }
#line 344 "clanhf.f"
		    i__1 = *n - 1;
#line 344 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 346 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 346 "clanhf.f"
			    value = temp;
#line 346 "clanhf.f"
			}
#line 348 "clanhf.f"
		    }
#line 349 "clanhf.f"
		    i__1 = k - 1;
#line 349 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 350 "clanhf.f"
			i__2 = j - 2;
#line 350 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 351 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 352 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 352 "clanhf.f"
				value = temp;
#line 352 "clanhf.f"
			    }
#line 354 "clanhf.f"
			}
#line 355 "clanhf.f"
			i__ = j - 1;
/*                    L(k+j,k+j) */
#line 357 "clanhf.f"
			i__2 = i__ + j * lda;
#line 357 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 358 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 358 "clanhf.f"
			    value = temp;
#line 358 "clanhf.f"
			}
#line 360 "clanhf.f"
			i__ = j;
/*                    -> L(j,j) */
#line 362 "clanhf.f"
			i__2 = i__ + j * lda;
#line 362 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 363 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 363 "clanhf.f"
			    value = temp;
#line 363 "clanhf.f"
			}
#line 365 "clanhf.f"
			i__2 = *n - 1;
#line 365 "clanhf.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 366 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 367 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 367 "clanhf.f"
				value = temp;
#line 367 "clanhf.f"
			    }
#line 369 "clanhf.f"
			}
#line 370 "clanhf.f"
		    }
#line 371 "clanhf.f"
		} else {
/*                 uplo = 'U' */
#line 373 "clanhf.f"
		    i__1 = k - 2;
#line 373 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 374 "clanhf.f"
			i__2 = k + j - 2;
#line 374 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 375 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 376 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 376 "clanhf.f"
				value = temp;
#line 376 "clanhf.f"
			    }
#line 378 "clanhf.f"
			}
#line 379 "clanhf.f"
			i__ = k + j - 1;
/*                    -> U(i,i) */
#line 381 "clanhf.f"
			i__2 = i__ + j * lda;
#line 381 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 382 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 382 "clanhf.f"
			    value = temp;
#line 382 "clanhf.f"
			}
#line 384 "clanhf.f"
			++i__;
/*                    =k+j; i -> U(j,j) */
#line 386 "clanhf.f"
			i__2 = i__ + j * lda;
#line 386 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 387 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 387 "clanhf.f"
			    value = temp;
#line 387 "clanhf.f"
			}
#line 389 "clanhf.f"
			i__2 = *n - 1;
#line 389 "clanhf.f"
			for (i__ = k + j + 1; i__ <= i__2; ++i__) {
#line 390 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 391 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 391 "clanhf.f"
				value = temp;
#line 391 "clanhf.f"
			    }
#line 393 "clanhf.f"
			}
#line 394 "clanhf.f"
		    }
#line 395 "clanhf.f"
		    i__1 = *n - 2;
#line 395 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 396 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 397 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 397 "clanhf.f"
			    value = temp;
#line 397 "clanhf.f"
			}
/*                    j=k-1 */
#line 400 "clanhf.f"
		    }
/*                 i=n-1 -> U(n-1,n-1) */
#line 402 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 402 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 403 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 403 "clanhf.f"
			value = temp;
#line 403 "clanhf.f"
		    }
#line 405 "clanhf.f"
		}
#line 406 "clanhf.f"
	    } else {
/*              xpose case; A is k by n */
#line 408 "clanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 410 "clanhf.f"
		    i__1 = k - 2;
#line 410 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 411 "clanhf.f"
			i__2 = j - 1;
#line 411 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 412 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 413 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 413 "clanhf.f"
				value = temp;
#line 413 "clanhf.f"
			    }
#line 415 "clanhf.f"
			}
#line 416 "clanhf.f"
			i__ = j;
/*                    L(i,i) */
#line 418 "clanhf.f"
			i__2 = i__ + j * lda;
#line 418 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 419 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 419 "clanhf.f"
			    value = temp;
#line 419 "clanhf.f"
			}
#line 421 "clanhf.f"
			i__ = j + 1;
/*                    L(j+k,j+k) */
#line 423 "clanhf.f"
			i__2 = i__ + j * lda;
#line 423 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 424 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 424 "clanhf.f"
			    value = temp;
#line 424 "clanhf.f"
			}
#line 426 "clanhf.f"
			i__2 = k - 1;
#line 426 "clanhf.f"
			for (i__ = j + 2; i__ <= i__2; ++i__) {
#line 427 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 428 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 428 "clanhf.f"
				value = temp;
#line 428 "clanhf.f"
			    }
#line 430 "clanhf.f"
			}
#line 431 "clanhf.f"
		    }
#line 432 "clanhf.f"
		    j = k - 1;
#line 433 "clanhf.f"
		    i__1 = k - 2;
#line 433 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 434 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 435 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 435 "clanhf.f"
			    value = temp;
#line 435 "clanhf.f"
			}
#line 437 "clanhf.f"
		    }
#line 438 "clanhf.f"
		    i__ = k - 1;
/*                 -> L(i,i) is at A(i,j) */
#line 440 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 440 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 441 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 441 "clanhf.f"
			value = temp;
#line 441 "clanhf.f"
		    }
#line 443 "clanhf.f"
		    i__1 = *n - 1;
#line 443 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 444 "clanhf.f"
			i__2 = k - 1;
#line 444 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 445 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 446 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 446 "clanhf.f"
				value = temp;
#line 446 "clanhf.f"
			    }
#line 448 "clanhf.f"
			}
#line 449 "clanhf.f"
		    }
#line 450 "clanhf.f"
		} else {
/*                 uplo = 'U' */
#line 452 "clanhf.f"
		    i__1 = k - 2;
#line 452 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 453 "clanhf.f"
			i__2 = k - 1;
#line 453 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 454 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 455 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 455 "clanhf.f"
				value = temp;
#line 455 "clanhf.f"
			    }
#line 457 "clanhf.f"
			}
#line 458 "clanhf.f"
		    }
#line 459 "clanhf.f"
		    j = k - 1;
/*                 -> U(j,j) is at A(0,j) */
#line 461 "clanhf.f"
		    i__1 = j * lda;
#line 461 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 462 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 462 "clanhf.f"
			value = temp;
#line 462 "clanhf.f"
		    }
#line 464 "clanhf.f"
		    i__1 = k - 1;
#line 464 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 465 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 466 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 466 "clanhf.f"
			    value = temp;
#line 466 "clanhf.f"
			}
#line 468 "clanhf.f"
		    }
#line 469 "clanhf.f"
		    i__1 = *n - 1;
#line 469 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 470 "clanhf.f"
			i__2 = j - k - 1;
#line 470 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 471 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 472 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 472 "clanhf.f"
				value = temp;
#line 472 "clanhf.f"
			    }
#line 474 "clanhf.f"
			}
#line 475 "clanhf.f"
			i__ = j - k;
/*                    -> U(i,i) at A(i,j) */
#line 477 "clanhf.f"
			i__2 = i__ + j * lda;
#line 477 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 478 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 478 "clanhf.f"
			    value = temp;
#line 478 "clanhf.f"
			}
#line 480 "clanhf.f"
			i__ = j - k + 1;
/*                    U(j,j) */
#line 482 "clanhf.f"
			i__2 = i__ + j * lda;
#line 482 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 483 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 483 "clanhf.f"
			    value = temp;
#line 483 "clanhf.f"
			}
#line 485 "clanhf.f"
			i__2 = k - 1;
#line 485 "clanhf.f"
			for (i__ = j - k + 2; i__ <= i__2; ++i__) {
#line 486 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 487 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 487 "clanhf.f"
				value = temp;
#line 487 "clanhf.f"
			    }
#line 489 "clanhf.f"
			}
#line 490 "clanhf.f"
		    }
#line 491 "clanhf.f"
		}
#line 492 "clanhf.f"
	    }
#line 493 "clanhf.f"
	} else {
/*           n is even & k = n/2 */
#line 495 "clanhf.f"
	    if (ifm == 1) {
/*              A is n+1 by k */
#line 497 "clanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 499 "clanhf.f"
		    j = 0;
/*                 -> L(k,k) & j=1 -> L(0,0) */
#line 501 "clanhf.f"
		    i__1 = j + j * lda;
#line 501 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 502 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 502 "clanhf.f"
			value = temp;
#line 502 "clanhf.f"
		    }
#line 504 "clanhf.f"
		    i__1 = j + 1 + j * lda;
#line 504 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 505 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 505 "clanhf.f"
			value = temp;
#line 505 "clanhf.f"
		    }
#line 507 "clanhf.f"
		    i__1 = *n;
#line 507 "clanhf.f"
		    for (i__ = 2; i__ <= i__1; ++i__) {
#line 508 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 509 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 509 "clanhf.f"
			    value = temp;
#line 509 "clanhf.f"
			}
#line 511 "clanhf.f"
		    }
#line 512 "clanhf.f"
		    i__1 = k - 1;
#line 512 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 513 "clanhf.f"
			i__2 = j - 1;
#line 513 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 514 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 515 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 515 "clanhf.f"
				value = temp;
#line 515 "clanhf.f"
			    }
#line 517 "clanhf.f"
			}
#line 518 "clanhf.f"
			i__ = j;
/*                    L(k+j,k+j) */
#line 520 "clanhf.f"
			i__2 = i__ + j * lda;
#line 520 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 521 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 521 "clanhf.f"
			    value = temp;
#line 521 "clanhf.f"
			}
#line 523 "clanhf.f"
			i__ = j + 1;
/*                    -> L(j,j) */
#line 525 "clanhf.f"
			i__2 = i__ + j * lda;
#line 525 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 526 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 526 "clanhf.f"
			    value = temp;
#line 526 "clanhf.f"
			}
#line 528 "clanhf.f"
			i__2 = *n;
#line 528 "clanhf.f"
			for (i__ = j + 2; i__ <= i__2; ++i__) {
#line 529 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 530 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 530 "clanhf.f"
				value = temp;
#line 530 "clanhf.f"
			    }
#line 532 "clanhf.f"
			}
#line 533 "clanhf.f"
		    }
#line 534 "clanhf.f"
		} else {
/*                 uplo = 'U' */
#line 536 "clanhf.f"
		    i__1 = k - 2;
#line 536 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 537 "clanhf.f"
			i__2 = k + j - 1;
#line 537 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 538 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 539 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 539 "clanhf.f"
				value = temp;
#line 539 "clanhf.f"
			    }
#line 541 "clanhf.f"
			}
#line 542 "clanhf.f"
			i__ = k + j;
/*                    -> U(i,i) */
#line 544 "clanhf.f"
			i__2 = i__ + j * lda;
#line 544 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 545 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 545 "clanhf.f"
			    value = temp;
#line 545 "clanhf.f"
			}
#line 547 "clanhf.f"
			++i__;
/*                    =k+j+1; i -> U(j,j) */
#line 549 "clanhf.f"
			i__2 = i__ + j * lda;
#line 549 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 550 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 550 "clanhf.f"
			    value = temp;
#line 550 "clanhf.f"
			}
#line 552 "clanhf.f"
			i__2 = *n;
#line 552 "clanhf.f"
			for (i__ = k + j + 2; i__ <= i__2; ++i__) {
#line 553 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 554 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 554 "clanhf.f"
				value = temp;
#line 554 "clanhf.f"
			    }
#line 556 "clanhf.f"
			}
#line 557 "clanhf.f"
		    }
#line 558 "clanhf.f"
		    i__1 = *n - 2;
#line 558 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 559 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 560 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 560 "clanhf.f"
			    value = temp;
#line 560 "clanhf.f"
			}
/*                 j=k-1 */
#line 563 "clanhf.f"
		    }
/*                 i=n-1 -> U(n-1,n-1) */
#line 565 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 565 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 566 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 566 "clanhf.f"
			value = temp;
#line 566 "clanhf.f"
		    }
#line 568 "clanhf.f"
		    i__ = *n;
/*                 -> U(k-1,k-1) */
#line 570 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 570 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 571 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 571 "clanhf.f"
			value = temp;
#line 571 "clanhf.f"
		    }
#line 573 "clanhf.f"
		}
#line 574 "clanhf.f"
	    } else {
/*              xpose case; A is k by n+1 */
#line 576 "clanhf.f"
		if (ilu == 1) {
/*                 uplo ='L' */
#line 578 "clanhf.f"
		    j = 0;
/*                 -> L(k,k) at A(0,0) */
#line 580 "clanhf.f"
		    i__1 = j + j * lda;
#line 580 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 581 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 581 "clanhf.f"
			value = temp;
#line 581 "clanhf.f"
		    }
#line 583 "clanhf.f"
		    i__1 = k - 1;
#line 583 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 584 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 585 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 585 "clanhf.f"
			    value = temp;
#line 585 "clanhf.f"
			}
#line 587 "clanhf.f"
		    }
#line 588 "clanhf.f"
		    i__1 = k - 1;
#line 588 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 589 "clanhf.f"
			i__2 = j - 2;
#line 589 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 590 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 591 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 591 "clanhf.f"
				value = temp;
#line 591 "clanhf.f"
			    }
#line 593 "clanhf.f"
			}
#line 594 "clanhf.f"
			i__ = j - 1;
/*                    L(i,i) */
#line 596 "clanhf.f"
			i__2 = i__ + j * lda;
#line 596 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 597 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 597 "clanhf.f"
			    value = temp;
#line 597 "clanhf.f"
			}
#line 599 "clanhf.f"
			i__ = j;
/*                    L(j+k,j+k) */
#line 601 "clanhf.f"
			i__2 = i__ + j * lda;
#line 601 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 602 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 602 "clanhf.f"
			    value = temp;
#line 602 "clanhf.f"
			}
#line 604 "clanhf.f"
			i__2 = k - 1;
#line 604 "clanhf.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 605 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 606 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 606 "clanhf.f"
				value = temp;
#line 606 "clanhf.f"
			    }
#line 608 "clanhf.f"
			}
#line 609 "clanhf.f"
		    }
#line 610 "clanhf.f"
		    j = k;
#line 611 "clanhf.f"
		    i__1 = k - 2;
#line 611 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 612 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 613 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 613 "clanhf.f"
			    value = temp;
#line 613 "clanhf.f"
			}
#line 615 "clanhf.f"
		    }
#line 616 "clanhf.f"
		    i__ = k - 1;
/*                 -> L(i,i) is at A(i,j) */
#line 618 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 618 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 619 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 619 "clanhf.f"
			value = temp;
#line 619 "clanhf.f"
		    }
#line 621 "clanhf.f"
		    i__1 = *n;
#line 621 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 622 "clanhf.f"
			i__2 = k - 1;
#line 622 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 623 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 624 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 624 "clanhf.f"
				value = temp;
#line 624 "clanhf.f"
			    }
#line 626 "clanhf.f"
			}
#line 627 "clanhf.f"
		    }
#line 628 "clanhf.f"
		} else {
/*                 uplo = 'U' */
#line 630 "clanhf.f"
		    i__1 = k - 1;
#line 630 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 631 "clanhf.f"
			i__2 = k - 1;
#line 631 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 632 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 633 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 633 "clanhf.f"
				value = temp;
#line 633 "clanhf.f"
			    }
#line 635 "clanhf.f"
			}
#line 636 "clanhf.f"
		    }
#line 637 "clanhf.f"
		    j = k;
/*                 -> U(j,j) is at A(0,j) */
#line 639 "clanhf.f"
		    i__1 = j * lda;
#line 639 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 640 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 640 "clanhf.f"
			value = temp;
#line 640 "clanhf.f"
		    }
#line 642 "clanhf.f"
		    i__1 = k - 1;
#line 642 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 643 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 644 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 644 "clanhf.f"
			    value = temp;
#line 644 "clanhf.f"
			}
#line 646 "clanhf.f"
		    }
#line 647 "clanhf.f"
		    i__1 = *n - 1;
#line 647 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 648 "clanhf.f"
			i__2 = j - k - 2;
#line 648 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 649 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 650 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 650 "clanhf.f"
				value = temp;
#line 650 "clanhf.f"
			    }
#line 652 "clanhf.f"
			}
#line 653 "clanhf.f"
			i__ = j - k - 1;
/*                    -> U(i,i) at A(i,j) */
#line 655 "clanhf.f"
			i__2 = i__ + j * lda;
#line 655 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 656 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 656 "clanhf.f"
			    value = temp;
#line 656 "clanhf.f"
			}
#line 658 "clanhf.f"
			i__ = j - k;
/*                    U(j,j) */
#line 660 "clanhf.f"
			i__2 = i__ + j * lda;
#line 660 "clanhf.f"
			temp = (d__1 = a[i__2].r, abs(d__1));
#line 661 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 661 "clanhf.f"
			    value = temp;
#line 661 "clanhf.f"
			}
#line 663 "clanhf.f"
			i__2 = k - 1;
#line 663 "clanhf.f"
			for (i__ = j - k + 1; i__ <= i__2; ++i__) {
#line 664 "clanhf.f"
			    temp = z_abs(&a[i__ + j * lda]);
#line 665 "clanhf.f"
			    if (value < temp || sisnan_(&temp)) {
#line 665 "clanhf.f"
				value = temp;
#line 665 "clanhf.f"
			    }
#line 667 "clanhf.f"
			}
#line 668 "clanhf.f"
		    }
#line 669 "clanhf.f"
		    j = *n;
#line 670 "clanhf.f"
		    i__1 = k - 2;
#line 670 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 671 "clanhf.f"
			temp = z_abs(&a[i__ + j * lda]);
#line 672 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 672 "clanhf.f"
			    value = temp;
#line 672 "clanhf.f"
			}
#line 674 "clanhf.f"
		    }
#line 675 "clanhf.f"
		    i__ = k - 1;
/*                 U(k,k) at A(i,j) */
#line 677 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 677 "clanhf.f"
		    temp = (d__1 = a[i__1].r, abs(d__1));
#line 678 "clanhf.f"
		    if (value < temp || sisnan_(&temp)) {
#line 678 "clanhf.f"
			value = temp;
#line 678 "clanhf.f"
		    }
#line 680 "clanhf.f"
		}
#line 681 "clanhf.f"
	    }
#line 682 "clanhf.f"
	}
#line 683 "clanhf.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*       Find normI(A) ( = norm1(A), since A is Hermitian). */

#line 688 "clanhf.f"
	if (ifm == 1) {
/*           A is 'N' */
#line 690 "clanhf.f"
	    k = *n / 2;
#line 691 "clanhf.f"
	    if (noe == 1) {
/*              n is odd & A is n by (n+1)/2 */
#line 693 "clanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 695 "clanhf.f"
		    i__1 = k - 1;
#line 695 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 696 "clanhf.f"
			work[i__] = 0.;
#line 697 "clanhf.f"
		    }
#line 698 "clanhf.f"
		    i__1 = k;
#line 698 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 699 "clanhf.f"
			s = 0.;
#line 700 "clanhf.f"
			i__2 = k + j - 1;
#line 700 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 701 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(i,j+k) */
#line 703 "clanhf.f"
			    s += aa;
#line 704 "clanhf.f"
			    work[i__] += aa;
#line 705 "clanhf.f"
			}
#line 706 "clanhf.f"
			i__2 = i__ + j * lda;
#line 706 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 708 "clanhf.f"
			work[j + k] = s + aa;
#line 709 "clanhf.f"
			if (i__ == k + k) {
#line 709 "clanhf.f"
			    goto L10;
#line 709 "clanhf.f"
			}
#line 711 "clanhf.f"
			++i__;
#line 712 "clanhf.f"
			i__2 = i__ + j * lda;
#line 712 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j,j) */
#line 714 "clanhf.f"
			work[j] += aa;
#line 715 "clanhf.f"
			s = 0.;
#line 716 "clanhf.f"
			i__2 = k - 1;
#line 716 "clanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 717 "clanhf.f"
			    ++i__;
#line 718 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 720 "clanhf.f"
			    s += aa;
#line 721 "clanhf.f"
			    work[l] += aa;
#line 722 "clanhf.f"
			}
#line 723 "clanhf.f"
			work[j] += s;
#line 724 "clanhf.f"
		    }
#line 725 "clanhf.f"
L10:
#line 726 "clanhf.f"
		    value = work[0];
#line 727 "clanhf.f"
		    i__1 = *n - 1;
#line 727 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 728 "clanhf.f"
			temp = work[i__];
#line 729 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 729 "clanhf.f"
			    value = temp;
#line 729 "clanhf.f"
			}
#line 731 "clanhf.f"
		    }
#line 732 "clanhf.f"
		} else {
/*                 ilu = 1 & uplo = 'L' */
#line 734 "clanhf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 736 "clanhf.f"
		    i__1 = *n - 1;
#line 736 "clanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 737 "clanhf.f"
			work[i__] = 0.;
#line 738 "clanhf.f"
		    }
#line 739 "clanhf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 740 "clanhf.f"
			s = 0.;
#line 741 "clanhf.f"
			i__1 = j - 2;
#line 741 "clanhf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 742 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(j+k,i+k) */
#line 744 "clanhf.f"
			    s += aa;
#line 745 "clanhf.f"
			    work[i__ + k] += aa;
#line 746 "clanhf.f"
			}
#line 747 "clanhf.f"
			if (j > 0) {
#line 748 "clanhf.f"
			    i__1 = i__ + j * lda;
#line 748 "clanhf.f"
			    aa = (d__1 = a[i__1].r, abs(d__1));
/*                       -> A(j+k,j+k) */
#line 750 "clanhf.f"
			    s += aa;
#line 751 "clanhf.f"
			    work[i__ + k] += s;
/*                       i=j */
#line 753 "clanhf.f"
			    ++i__;
#line 754 "clanhf.f"
			}
#line 755 "clanhf.f"
			i__1 = i__ + j * lda;
#line 755 "clanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j,j) */
#line 757 "clanhf.f"
			work[j] = aa;
#line 758 "clanhf.f"
			s = 0.;
#line 759 "clanhf.f"
			i__1 = *n - 1;
#line 759 "clanhf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 760 "clanhf.f"
			    ++i__;
#line 761 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 763 "clanhf.f"
			    s += aa;
#line 764 "clanhf.f"
			    work[l] += aa;
#line 765 "clanhf.f"
			}
#line 766 "clanhf.f"
			work[j] += s;
#line 767 "clanhf.f"
		    }
#line 768 "clanhf.f"
		    value = work[0];
#line 769 "clanhf.f"
		    i__1 = *n - 1;
#line 769 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 770 "clanhf.f"
			temp = work[i__];
#line 771 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 771 "clanhf.f"
			    value = temp;
#line 771 "clanhf.f"
			}
#line 773 "clanhf.f"
		    }
#line 774 "clanhf.f"
		}
#line 775 "clanhf.f"
	    } else {
/*              n is even & A is n+1 by k = n/2 */
#line 777 "clanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 779 "clanhf.f"
		    i__1 = k - 1;
#line 779 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 780 "clanhf.f"
			work[i__] = 0.;
#line 781 "clanhf.f"
		    }
#line 782 "clanhf.f"
		    i__1 = k - 1;
#line 782 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 783 "clanhf.f"
			s = 0.;
#line 784 "clanhf.f"
			i__2 = k + j - 1;
#line 784 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 785 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(i,j+k) */
#line 787 "clanhf.f"
			    s += aa;
#line 788 "clanhf.f"
			    work[i__] += aa;
#line 789 "clanhf.f"
			}
#line 790 "clanhf.f"
			i__2 = i__ + j * lda;
#line 790 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 792 "clanhf.f"
			work[j + k] = s + aa;
#line 793 "clanhf.f"
			++i__;
#line 794 "clanhf.f"
			i__2 = i__ + j * lda;
#line 794 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    -> A(j,j) */
#line 796 "clanhf.f"
			work[j] += aa;
#line 797 "clanhf.f"
			s = 0.;
#line 798 "clanhf.f"
			i__2 = k - 1;
#line 798 "clanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 799 "clanhf.f"
			    ++i__;
#line 800 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 802 "clanhf.f"
			    s += aa;
#line 803 "clanhf.f"
			    work[l] += aa;
#line 804 "clanhf.f"
			}
#line 805 "clanhf.f"
			work[j] += s;
#line 806 "clanhf.f"
		    }
#line 807 "clanhf.f"
		    value = work[0];
#line 808 "clanhf.f"
		    i__1 = *n - 1;
#line 808 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 809 "clanhf.f"
			temp = work[i__];
#line 810 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 810 "clanhf.f"
			    value = temp;
#line 810 "clanhf.f"
			}
#line 812 "clanhf.f"
		    }
#line 813 "clanhf.f"
		} else {
/*                 ilu = 1 & uplo = 'L' */
#line 815 "clanhf.f"
		    i__1 = *n - 1;
#line 815 "clanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 816 "clanhf.f"
			work[i__] = 0.;
#line 817 "clanhf.f"
		    }
#line 818 "clanhf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 819 "clanhf.f"
			s = 0.;
#line 820 "clanhf.f"
			i__1 = j - 1;
#line 820 "clanhf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 821 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(j+k,i+k) */
#line 823 "clanhf.f"
			    s += aa;
#line 824 "clanhf.f"
			    work[i__ + k] += aa;
#line 825 "clanhf.f"
			}
#line 826 "clanhf.f"
			i__1 = i__ + j * lda;
#line 826 "clanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j+k,j+k) */
#line 828 "clanhf.f"
			s += aa;
#line 829 "clanhf.f"
			work[i__ + k] += s;
/*                    i=j */
#line 831 "clanhf.f"
			++i__;
#line 832 "clanhf.f"
			i__1 = i__ + j * lda;
#line 832 "clanhf.f"
			aa = (d__1 = a[i__1].r, abs(d__1));
/*                    -> A(j,j) */
#line 834 "clanhf.f"
			work[j] = aa;
#line 835 "clanhf.f"
			s = 0.;
#line 836 "clanhf.f"
			i__1 = *n - 1;
#line 836 "clanhf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 837 "clanhf.f"
			    ++i__;
#line 838 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       -> A(l,j) */
#line 840 "clanhf.f"
			    s += aa;
#line 841 "clanhf.f"
			    work[l] += aa;
#line 842 "clanhf.f"
			}
#line 843 "clanhf.f"
			work[j] += s;
#line 844 "clanhf.f"
		    }
#line 845 "clanhf.f"
		    value = work[0];
#line 846 "clanhf.f"
		    i__1 = *n - 1;
#line 846 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 847 "clanhf.f"
			temp = work[i__];
#line 848 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 848 "clanhf.f"
			    value = temp;
#line 848 "clanhf.f"
			}
#line 850 "clanhf.f"
		    }
#line 851 "clanhf.f"
		}
#line 852 "clanhf.f"
	    }
#line 853 "clanhf.f"
	} else {
/*           ifm=0 */
#line 855 "clanhf.f"
	    k = *n / 2;
#line 856 "clanhf.f"
	    if (noe == 1) {
/*              n is odd & A is (n+1)/2 by n */
#line 858 "clanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 860 "clanhf.f"
		    n1 = k;
/*                 n/2 */
#line 862 "clanhf.f"
		    ++k;
/*                 k is the row size and lda */
#line 864 "clanhf.f"
		    i__1 = *n - 1;
#line 864 "clanhf.f"
		    for (i__ = n1; i__ <= i__1; ++i__) {
#line 865 "clanhf.f"
			work[i__] = 0.;
#line 866 "clanhf.f"
		    }
#line 867 "clanhf.f"
		    i__1 = n1 - 1;
#line 867 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 868 "clanhf.f"
			s = 0.;
#line 869 "clanhf.f"
			i__2 = k - 1;
#line 869 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 870 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,n1+i) */
#line 872 "clanhf.f"
			    work[i__ + n1] += aa;
#line 873 "clanhf.f"
			    s += aa;
#line 874 "clanhf.f"
			}
#line 875 "clanhf.f"
			work[j] = s;
#line 876 "clanhf.f"
		    }
/*                 j=n1=k-1 is special */
#line 878 "clanhf.f"
		    i__1 = j * lda;
#line 878 "clanhf.f"
		    s = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 880 "clanhf.f"
		    i__1 = k - 1;
#line 880 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 881 "clanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k-1,i+n1) */
#line 883 "clanhf.f"
			work[i__ + n1] += aa;
#line 884 "clanhf.f"
			s += aa;
#line 885 "clanhf.f"
		    }
#line 886 "clanhf.f"
		    work[j] += s;
#line 887 "clanhf.f"
		    i__1 = *n - 1;
#line 887 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 888 "clanhf.f"
			s = 0.;
#line 889 "clanhf.f"
			i__2 = j - k - 1;
#line 889 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 890 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(i,j-k) */
#line 892 "clanhf.f"
			    work[i__] += aa;
#line 893 "clanhf.f"
			    s += aa;
#line 894 "clanhf.f"
			}
/*                    i=j-k */
#line 896 "clanhf.f"
			i__2 = i__ + j * lda;
#line 896 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j-k,j-k) */
#line 898 "clanhf.f"
			s += aa;
#line 899 "clanhf.f"
			work[j - k] += s;
#line 900 "clanhf.f"
			++i__;
#line 901 "clanhf.f"
			i__2 = i__ + j * lda;
#line 901 "clanhf.f"
			s = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j,j) */
#line 903 "clanhf.f"
			i__2 = *n - 1;
#line 903 "clanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 904 "clanhf.f"
			    ++i__;
#line 905 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,l) */
#line 907 "clanhf.f"
			    work[l] += aa;
#line 908 "clanhf.f"
			    s += aa;
#line 909 "clanhf.f"
			}
#line 910 "clanhf.f"
			work[j] += s;
#line 911 "clanhf.f"
		    }
#line 912 "clanhf.f"
		    value = work[0];
#line 913 "clanhf.f"
		    i__1 = *n - 1;
#line 913 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 914 "clanhf.f"
			temp = work[i__];
#line 915 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 915 "clanhf.f"
			    value = temp;
#line 915 "clanhf.f"
			}
#line 917 "clanhf.f"
		    }
#line 918 "clanhf.f"
		} else {
/*                 ilu=1 & uplo = 'L' */
#line 920 "clanhf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 922 "clanhf.f"
		    i__1 = *n - 1;
#line 922 "clanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 923 "clanhf.f"
			work[i__] = 0.;
#line 924 "clanhf.f"
		    }
#line 925 "clanhf.f"
		    i__1 = k - 2;
#line 925 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
/*                    process */
#line 927 "clanhf.f"
			s = 0.;
#line 928 "clanhf.f"
			i__2 = j - 1;
#line 928 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 929 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i) */
#line 931 "clanhf.f"
			    work[i__] += aa;
#line 932 "clanhf.f"
			    s += aa;
#line 933 "clanhf.f"
			}
#line 934 "clanhf.f"
			i__2 = i__ + j * lda;
#line 934 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    i=j so process of A(j,j) */
#line 936 "clanhf.f"
			s += aa;
#line 937 "clanhf.f"
			work[j] = s;
/*                    is initialised here */
#line 939 "clanhf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 941 "clanhf.f"
			i__2 = i__ + j * lda;
#line 941 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
#line 942 "clanhf.f"
			s = aa;
#line 943 "clanhf.f"
			i__2 = *n - 1;
#line 943 "clanhf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 944 "clanhf.f"
			    ++i__;
#line 945 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(l,k+j) */
#line 947 "clanhf.f"
			    s += aa;
#line 948 "clanhf.f"
			    work[l] += aa;
#line 949 "clanhf.f"
			}
#line 950 "clanhf.f"
			work[k + j] += s;
#line 951 "clanhf.f"
		    }
/*                 j=k-1 is special :process col A(k-1,0:k-1) */
#line 953 "clanhf.f"
		    s = 0.;
#line 954 "clanhf.f"
		    i__1 = k - 2;
#line 954 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 955 "clanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,i) */
#line 957 "clanhf.f"
			work[i__] += aa;
#line 958 "clanhf.f"
			s += aa;
#line 959 "clanhf.f"
		    }
/*                 i=k-1 */
#line 961 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 961 "clanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 963 "clanhf.f"
		    s += aa;
#line 964 "clanhf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 966 "clanhf.f"
		    i__1 = *n - 1;
#line 966 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
/*                    process col j of A = A(j,0:k-1) */
#line 968 "clanhf.f"
			s = 0.;
#line 969 "clanhf.f"
			i__2 = k - 1;
#line 969 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 970 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i) */
#line 972 "clanhf.f"
			    work[i__] += aa;
#line 973 "clanhf.f"
			    s += aa;
#line 974 "clanhf.f"
			}
#line 975 "clanhf.f"
			work[j] += s;
#line 976 "clanhf.f"
		    }
#line 977 "clanhf.f"
		    value = work[0];
#line 978 "clanhf.f"
		    i__1 = *n - 1;
#line 978 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 979 "clanhf.f"
			temp = work[i__];
#line 980 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 980 "clanhf.f"
			    value = temp;
#line 980 "clanhf.f"
			}
#line 982 "clanhf.f"
		    }
#line 983 "clanhf.f"
		}
#line 984 "clanhf.f"
	    } else {
/*              n is even & A is k=n/2 by n+1 */
#line 986 "clanhf.f"
		if (ilu == 0) {
/*                 uplo = 'U' */
#line 988 "clanhf.f"
		    i__1 = *n - 1;
#line 988 "clanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 989 "clanhf.f"
			work[i__] = 0.;
#line 990 "clanhf.f"
		    }
#line 991 "clanhf.f"
		    i__1 = k - 1;
#line 991 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 992 "clanhf.f"
			s = 0.;
#line 993 "clanhf.f"
			i__2 = k - 1;
#line 993 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 994 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,i+k) */
#line 996 "clanhf.f"
			    work[i__ + k] += aa;
#line 997 "clanhf.f"
			    s += aa;
#line 998 "clanhf.f"
			}
#line 999 "clanhf.f"
			work[j] = s;
#line 1000 "clanhf.f"
		    }
/*                 j=k */
#line 1002 "clanhf.f"
		    i__1 = j * lda;
#line 1002 "clanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k,k) */
#line 1004 "clanhf.f"
		    s = aa;
#line 1005 "clanhf.f"
		    i__1 = k - 1;
#line 1005 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1006 "clanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,k+i) */
#line 1008 "clanhf.f"
			work[i__ + k] += aa;
#line 1009 "clanhf.f"
			s += aa;
#line 1010 "clanhf.f"
		    }
#line 1011 "clanhf.f"
		    work[j] += s;
#line 1012 "clanhf.f"
		    i__1 = *n - 1;
#line 1012 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1013 "clanhf.f"
			s = 0.;
#line 1014 "clanhf.f"
			i__2 = j - 2 - k;
#line 1014 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1015 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(i,j-k-1) */
#line 1017 "clanhf.f"
			    work[i__] += aa;
#line 1018 "clanhf.f"
			    s += aa;
#line 1019 "clanhf.f"
			}
/*                    i=j-1-k */
#line 1021 "clanhf.f"
			i__2 = i__ + j * lda;
#line 1021 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j-k-1,j-k-1) */
#line 1023 "clanhf.f"
			s += aa;
#line 1024 "clanhf.f"
			work[j - k - 1] += s;
#line 1025 "clanhf.f"
			++i__;
#line 1026 "clanhf.f"
			i__2 = i__ + j * lda;
#line 1026 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    A(j,j) */
#line 1028 "clanhf.f"
			s = aa;
#line 1029 "clanhf.f"
			i__2 = *n - 1;
#line 1029 "clanhf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 1030 "clanhf.f"
			    ++i__;
#line 1031 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j,l) */
#line 1033 "clanhf.f"
			    work[l] += aa;
#line 1034 "clanhf.f"
			    s += aa;
#line 1035 "clanhf.f"
			}
#line 1036 "clanhf.f"
			work[j] += s;
#line 1037 "clanhf.f"
		    }
/*                 j=n */
#line 1039 "clanhf.f"
		    s = 0.;
#line 1040 "clanhf.f"
		    i__1 = k - 2;
#line 1040 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1041 "clanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(i,k-1) */
#line 1043 "clanhf.f"
			work[i__] += aa;
#line 1044 "clanhf.f"
			s += aa;
#line 1045 "clanhf.f"
		    }
/*                 i=k-1 */
#line 1047 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 1047 "clanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 1049 "clanhf.f"
		    s += aa;
#line 1050 "clanhf.f"
		    work[i__] += s;
#line 1051 "clanhf.f"
		    value = work[0];
#line 1052 "clanhf.f"
		    i__1 = *n - 1;
#line 1052 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1053 "clanhf.f"
			temp = work[i__];
#line 1054 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 1054 "clanhf.f"
			    value = temp;
#line 1054 "clanhf.f"
			}
#line 1056 "clanhf.f"
		    }
#line 1057 "clanhf.f"
		} else {
/*                 ilu=1 & uplo = 'L' */
#line 1059 "clanhf.f"
		    i__1 = *n - 1;
#line 1059 "clanhf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 1060 "clanhf.f"
			work[i__] = 0.;
#line 1061 "clanhf.f"
		    }
/*                 j=0 is special :process col A(k:n-1,k) */
#line 1063 "clanhf.f"
		    s = (d__1 = a[0].r, abs(d__1));
/*                 A(k,k) */
#line 1065 "clanhf.f"
		    i__1 = k - 1;
#line 1065 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1066 "clanhf.f"
			aa = z_abs(&a[i__]);
/*                    A(k+i,k) */
#line 1068 "clanhf.f"
			work[i__ + k] += aa;
#line 1069 "clanhf.f"
			s += aa;
#line 1070 "clanhf.f"
		    }
#line 1071 "clanhf.f"
		    work[k] += s;
#line 1072 "clanhf.f"
		    i__1 = k - 1;
#line 1072 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
/*                    process */
#line 1074 "clanhf.f"
			s = 0.;
#line 1075 "clanhf.f"
			i__2 = j - 2;
#line 1075 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1076 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j-1,i) */
#line 1078 "clanhf.f"
			    work[i__] += aa;
#line 1079 "clanhf.f"
			    s += aa;
#line 1080 "clanhf.f"
			}
#line 1081 "clanhf.f"
			i__2 = i__ + j * lda;
#line 1081 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
/*                    i=j-1 so process of A(j-1,j-1) */
#line 1083 "clanhf.f"
			s += aa;
#line 1084 "clanhf.f"
			work[j - 1] = s;
/*                    is initialised here */
#line 1086 "clanhf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 1088 "clanhf.f"
			i__2 = i__ + j * lda;
#line 1088 "clanhf.f"
			aa = (d__1 = a[i__2].r, abs(d__1));
#line 1089 "clanhf.f"
			s = aa;
#line 1090 "clanhf.f"
			i__2 = *n - 1;
#line 1090 "clanhf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 1091 "clanhf.f"
			    ++i__;
#line 1092 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(l,k+j) */
#line 1094 "clanhf.f"
			    s += aa;
#line 1095 "clanhf.f"
			    work[l] += aa;
#line 1096 "clanhf.f"
			}
#line 1097 "clanhf.f"
			work[k + j] += s;
#line 1098 "clanhf.f"
		    }
/*                 j=k is special :process col A(k,0:k-1) */
#line 1100 "clanhf.f"
		    s = 0.;
#line 1101 "clanhf.f"
		    i__1 = k - 2;
#line 1101 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1102 "clanhf.f"
			aa = z_abs(&a[i__ + j * lda]);
/*                    A(k,i) */
#line 1104 "clanhf.f"
			work[i__] += aa;
#line 1105 "clanhf.f"
			s += aa;
#line 1106 "clanhf.f"
		    }

/*                 i=k-1 */
#line 1109 "clanhf.f"
		    i__1 = i__ + j * lda;
#line 1109 "clanhf.f"
		    aa = (d__1 = a[i__1].r, abs(d__1));
/*                 A(k-1,k-1) */
#line 1111 "clanhf.f"
		    s += aa;
#line 1112 "clanhf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 1114 "clanhf.f"
		    i__1 = *n;
#line 1114 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {

/*                    process col j-1 of A = A(j-1,0:k-1) */
#line 1117 "clanhf.f"
			s = 0.;
#line 1118 "clanhf.f"
			i__2 = k - 1;
#line 1118 "clanhf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 1119 "clanhf.f"
			    aa = z_abs(&a[i__ + j * lda]);
/*                       A(j-1,i) */
#line 1121 "clanhf.f"
			    work[i__] += aa;
#line 1122 "clanhf.f"
			    s += aa;
#line 1123 "clanhf.f"
			}
#line 1124 "clanhf.f"
			work[j - 1] += s;
#line 1125 "clanhf.f"
		    }
#line 1126 "clanhf.f"
		    value = work[0];
#line 1127 "clanhf.f"
		    i__1 = *n - 1;
#line 1127 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1128 "clanhf.f"
			temp = work[i__];
#line 1129 "clanhf.f"
			if (value < temp || sisnan_(&temp)) {
#line 1129 "clanhf.f"
			    value = temp;
#line 1129 "clanhf.f"
			}
#line 1131 "clanhf.f"
		    }
#line 1132 "clanhf.f"
		}
#line 1133 "clanhf.f"
	    }
#line 1134 "clanhf.f"
	}
#line 1135 "clanhf.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*       Find normF(A). */

#line 1139 "clanhf.f"
	k = (*n + 1) / 2;
#line 1140 "clanhf.f"
	scale = 0.;
#line 1141 "clanhf.f"
	s = 1.;
#line 1142 "clanhf.f"
	if (noe == 1) {
/*           n is odd */
#line 1144 "clanhf.f"
	    if (ifm == 1) {
/*              A is normal & A is n by k */
#line 1146 "clanhf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 1148 "clanhf.f"
		    i__1 = k - 3;
#line 1148 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1149 "clanhf.f"
			i__2 = k - j - 2;
#line 1149 "clanhf.f"
			classq_(&i__2, &a[k + j + 1 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k,0) */
#line 1151 "clanhf.f"
		    }
#line 1152 "clanhf.f"
		    i__1 = k - 1;
#line 1152 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1153 "clanhf.f"
			i__2 = k + j - 1;
#line 1153 "clanhf.f"
			classq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 1155 "clanhf.f"
		    }
#line 1156 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1158 "clanhf.f"
		    l = k - 1;
/*                 -> U(k,k) at A(k-1,0) */
#line 1160 "clanhf.f"
		    i__1 = k - 2;
#line 1160 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1161 "clanhf.f"
			i__2 = l;
#line 1161 "clanhf.f"
			aa = a[i__2].r;
/*                    U(k+i,k+i) */
#line 1163 "clanhf.f"
			if (aa != 0.) {
#line 1164 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1165 "clanhf.f"
				d__1 = scale / aa;
#line 1165 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1166 "clanhf.f"
				scale = aa;
#line 1167 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1168 "clanhf.f"
				d__1 = aa / scale;
#line 1168 "clanhf.f"
				s += d__1 * d__1;
#line 1169 "clanhf.f"
			    }
#line 1170 "clanhf.f"
			}
#line 1171 "clanhf.f"
			i__2 = l + 1;
#line 1171 "clanhf.f"
			aa = a[i__2].r;
/*                    U(i,i) */
#line 1173 "clanhf.f"
			if (aa != 0.) {
#line 1174 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1175 "clanhf.f"
				d__1 = scale / aa;
#line 1175 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1176 "clanhf.f"
				scale = aa;
#line 1177 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1178 "clanhf.f"
				d__1 = aa / scale;
#line 1178 "clanhf.f"
				s += d__1 * d__1;
#line 1179 "clanhf.f"
			    }
#line 1180 "clanhf.f"
			}
#line 1181 "clanhf.f"
			l = l + lda + 1;
#line 1182 "clanhf.f"
		    }
#line 1183 "clanhf.f"
		    i__1 = l;
#line 1183 "clanhf.f"
		    aa = a[i__1].r;
/*                 U(n-1,n-1) */
#line 1185 "clanhf.f"
		    if (aa != 0.) {
#line 1186 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1187 "clanhf.f"
			    d__1 = scale / aa;
#line 1187 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1188 "clanhf.f"
			    scale = aa;
#line 1189 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1190 "clanhf.f"
			    d__1 = aa / scale;
#line 1190 "clanhf.f"
			    s += d__1 * d__1;
#line 1191 "clanhf.f"
			}
#line 1192 "clanhf.f"
		    }
#line 1193 "clanhf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 1195 "clanhf.f"
		    i__1 = k - 1;
#line 1195 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1196 "clanhf.f"
			i__2 = *n - j - 1;
#line 1196 "clanhf.f"
			classq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(0,0) */
#line 1198 "clanhf.f"
		    }
#line 1199 "clanhf.f"
		    i__1 = k - 2;
#line 1199 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1200 "clanhf.f"
			classq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 1202 "clanhf.f"
		    }
#line 1203 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1205 "clanhf.f"
		    aa = a[0].r;
/*                 L(0,0) at A(0,0) */
#line 1207 "clanhf.f"
		    if (aa != 0.) {
#line 1208 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1209 "clanhf.f"
			    d__1 = scale / aa;
#line 1209 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1210 "clanhf.f"
			    scale = aa;
#line 1211 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1212 "clanhf.f"
			    d__1 = aa / scale;
#line 1212 "clanhf.f"
			    s += d__1 * d__1;
#line 1213 "clanhf.f"
			}
#line 1214 "clanhf.f"
		    }
#line 1215 "clanhf.f"
		    l = lda;
/*                 -> L(k,k) at A(0,1) */
#line 1217 "clanhf.f"
		    i__1 = k - 1;
#line 1217 "clanhf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1218 "clanhf.f"
			i__2 = l;
#line 1218 "clanhf.f"
			aa = a[i__2].r;
/*                    L(k-1+i,k-1+i) */
#line 1220 "clanhf.f"
			if (aa != 0.) {
#line 1221 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1222 "clanhf.f"
				d__1 = scale / aa;
#line 1222 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1223 "clanhf.f"
				scale = aa;
#line 1224 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1225 "clanhf.f"
				d__1 = aa / scale;
#line 1225 "clanhf.f"
				s += d__1 * d__1;
#line 1226 "clanhf.f"
			    }
#line 1227 "clanhf.f"
			}
#line 1228 "clanhf.f"
			i__2 = l + 1;
#line 1228 "clanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1230 "clanhf.f"
			if (aa != 0.) {
#line 1231 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1232 "clanhf.f"
				d__1 = scale / aa;
#line 1232 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1233 "clanhf.f"
				scale = aa;
#line 1234 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1235 "clanhf.f"
				d__1 = aa / scale;
#line 1235 "clanhf.f"
				s += d__1 * d__1;
#line 1236 "clanhf.f"
			    }
#line 1237 "clanhf.f"
			}
#line 1238 "clanhf.f"
			l = l + lda + 1;
#line 1239 "clanhf.f"
		    }
#line 1240 "clanhf.f"
		}
#line 1241 "clanhf.f"
	    } else {
/*              A is xpose & A is k by n */
#line 1243 "clanhf.f"
		if (ilu == 0) {
/*                 A**H is upper */
#line 1245 "clanhf.f"
		    i__1 = k - 2;
#line 1245 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1246 "clanhf.f"
			classq_(&j, &a[(k + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k) */
#line 1248 "clanhf.f"
		    }
#line 1249 "clanhf.f"
		    i__1 = k - 2;
#line 1249 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1250 "clanhf.f"
			classq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,0) */
#line 1252 "clanhf.f"
		    }
#line 1253 "clanhf.f"
		    i__1 = k - 2;
#line 1253 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1254 "clanhf.f"
			i__2 = k - j - 1;
#line 1254 "clanhf.f"
			classq_(&i__2, &a[j + 1 + (j + k - 1) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k-1) */
#line 1257 "clanhf.f"
		    }
#line 1258 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1260 "clanhf.f"
		    l = k * lda - lda;
/*                 -> U(k-1,k-1) at A(0,k-1) */
#line 1262 "clanhf.f"
		    i__1 = l;
#line 1262 "clanhf.f"
		    aa = a[i__1].r;
/*                 U(k-1,k-1) */
#line 1264 "clanhf.f"
		    if (aa != 0.) {
#line 1265 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1266 "clanhf.f"
			    d__1 = scale / aa;
#line 1266 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1267 "clanhf.f"
			    scale = aa;
#line 1268 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1269 "clanhf.f"
			    d__1 = aa / scale;
#line 1269 "clanhf.f"
			    s += d__1 * d__1;
#line 1270 "clanhf.f"
			}
#line 1271 "clanhf.f"
		    }
#line 1272 "clanhf.f"
		    l += lda;
/*                 -> U(0,0) at A(0,k) */
#line 1274 "clanhf.f"
		    i__1 = *n - 1;
#line 1274 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 1275 "clanhf.f"
			i__2 = l;
#line 1275 "clanhf.f"
			aa = a[i__2].r;
/*                    -> U(j-k,j-k) */
#line 1277 "clanhf.f"
			if (aa != 0.) {
#line 1278 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1279 "clanhf.f"
				d__1 = scale / aa;
#line 1279 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1280 "clanhf.f"
				scale = aa;
#line 1281 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1282 "clanhf.f"
				d__1 = aa / scale;
#line 1282 "clanhf.f"
				s += d__1 * d__1;
#line 1283 "clanhf.f"
			    }
#line 1284 "clanhf.f"
			}
#line 1285 "clanhf.f"
			i__2 = l + 1;
#line 1285 "clanhf.f"
			aa = a[i__2].r;
/*                    -> U(j,j) */
#line 1287 "clanhf.f"
			if (aa != 0.) {
#line 1288 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1289 "clanhf.f"
				d__1 = scale / aa;
#line 1289 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1290 "clanhf.f"
				scale = aa;
#line 1291 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1292 "clanhf.f"
				d__1 = aa / scale;
#line 1292 "clanhf.f"
				s += d__1 * d__1;
#line 1293 "clanhf.f"
			    }
#line 1294 "clanhf.f"
			}
#line 1295 "clanhf.f"
			l = l + lda + 1;
#line 1296 "clanhf.f"
		    }
#line 1297 "clanhf.f"
		} else {
/*                 A**H is lower */
#line 1299 "clanhf.f"
		    i__1 = k - 1;
#line 1299 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1300 "clanhf.f"
			classq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 1302 "clanhf.f"
		    }
#line 1303 "clanhf.f"
		    i__1 = *n - 1;
#line 1303 "clanhf.f"
		    for (j = k; j <= i__1; ++j) {
#line 1304 "clanhf.f"
			classq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,k) */
#line 1306 "clanhf.f"
		    }
#line 1307 "clanhf.f"
		    i__1 = k - 3;
#line 1307 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1308 "clanhf.f"
			i__2 = k - j - 2;
#line 1308 "clanhf.f"
			classq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(1,0) */
#line 1310 "clanhf.f"
		    }
#line 1311 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1313 "clanhf.f"
		    l = 0;
/*                 -> L(0,0) at A(0,0) */
#line 1315 "clanhf.f"
		    i__1 = k - 2;
#line 1315 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1316 "clanhf.f"
			i__2 = l;
#line 1316 "clanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1318 "clanhf.f"
			if (aa != 0.) {
#line 1319 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1320 "clanhf.f"
				d__1 = scale / aa;
#line 1320 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1321 "clanhf.f"
				scale = aa;
#line 1322 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1323 "clanhf.f"
				d__1 = aa / scale;
#line 1323 "clanhf.f"
				s += d__1 * d__1;
#line 1324 "clanhf.f"
			    }
#line 1325 "clanhf.f"
			}
#line 1326 "clanhf.f"
			i__2 = l + 1;
#line 1326 "clanhf.f"
			aa = a[i__2].r;
/*                    L(k+i,k+i) */
#line 1328 "clanhf.f"
			if (aa != 0.) {
#line 1329 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1330 "clanhf.f"
				d__1 = scale / aa;
#line 1330 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1331 "clanhf.f"
				scale = aa;
#line 1332 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1333 "clanhf.f"
				d__1 = aa / scale;
#line 1333 "clanhf.f"
				s += d__1 * d__1;
#line 1334 "clanhf.f"
			    }
#line 1335 "clanhf.f"
			}
#line 1336 "clanhf.f"
			l = l + lda + 1;
#line 1337 "clanhf.f"
		    }
/*                 L-> k-1 + (k-1)*lda or L(k-1,k-1) at A(k-1,k-1) */
#line 1339 "clanhf.f"
		    i__1 = l;
#line 1339 "clanhf.f"
		    aa = a[i__1].r;
/*                 L(k-1,k-1) at A(k-1,k-1) */
#line 1341 "clanhf.f"
		    if (aa != 0.) {
#line 1342 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1343 "clanhf.f"
			    d__1 = scale / aa;
#line 1343 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1344 "clanhf.f"
			    scale = aa;
#line 1345 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1346 "clanhf.f"
			    d__1 = aa / scale;
#line 1346 "clanhf.f"
			    s += d__1 * d__1;
#line 1347 "clanhf.f"
			}
#line 1348 "clanhf.f"
		    }
#line 1349 "clanhf.f"
		}
#line 1350 "clanhf.f"
	    }
#line 1351 "clanhf.f"
	} else {
/*           n is even */
#line 1353 "clanhf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 1355 "clanhf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 1357 "clanhf.f"
		    i__1 = k - 2;
#line 1357 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1358 "clanhf.f"
			i__2 = k - j - 1;
#line 1358 "clanhf.f"
			classq_(&i__2, &a[k + j + 2 + j * lda], &c__1, &scale,
				 &s);
/*                 L at A(k+1,0) */
#line 1360 "clanhf.f"
		    }
#line 1361 "clanhf.f"
		    i__1 = k - 1;
#line 1361 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1362 "clanhf.f"
			i__2 = k + j;
#line 1362 "clanhf.f"
			classq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                 trap U at A(0,0) */
#line 1364 "clanhf.f"
		    }
#line 1365 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1367 "clanhf.f"
		    l = k;
/*                 -> U(k,k) at A(k,0) */
#line 1369 "clanhf.f"
		    i__1 = k - 1;
#line 1369 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1370 "clanhf.f"
			i__2 = l;
#line 1370 "clanhf.f"
			aa = a[i__2].r;
/*                    U(k+i,k+i) */
#line 1372 "clanhf.f"
			if (aa != 0.) {
#line 1373 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1374 "clanhf.f"
				d__1 = scale / aa;
#line 1374 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1375 "clanhf.f"
				scale = aa;
#line 1376 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1377 "clanhf.f"
				d__1 = aa / scale;
#line 1377 "clanhf.f"
				s += d__1 * d__1;
#line 1378 "clanhf.f"
			    }
#line 1379 "clanhf.f"
			}
#line 1380 "clanhf.f"
			i__2 = l + 1;
#line 1380 "clanhf.f"
			aa = a[i__2].r;
/*                    U(i,i) */
#line 1382 "clanhf.f"
			if (aa != 0.) {
#line 1383 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1384 "clanhf.f"
				d__1 = scale / aa;
#line 1384 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1385 "clanhf.f"
				scale = aa;
#line 1386 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1387 "clanhf.f"
				d__1 = aa / scale;
#line 1387 "clanhf.f"
				s += d__1 * d__1;
#line 1388 "clanhf.f"
			    }
#line 1389 "clanhf.f"
			}
#line 1390 "clanhf.f"
			l = l + lda + 1;
#line 1391 "clanhf.f"
		    }
#line 1392 "clanhf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 1394 "clanhf.f"
		    i__1 = k - 1;
#line 1394 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1395 "clanhf.f"
			i__2 = *n - j - 1;
#line 1395 "clanhf.f"
			classq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(1,0) */
#line 1397 "clanhf.f"
		    }
#line 1398 "clanhf.f"
		    i__1 = k - 1;
#line 1398 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1399 "clanhf.f"
			classq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 1401 "clanhf.f"
		    }
#line 1402 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1404 "clanhf.f"
		    l = 0;
/*                 -> L(k,k) at A(0,0) */
#line 1406 "clanhf.f"
		    i__1 = k - 1;
#line 1406 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1407 "clanhf.f"
			i__2 = l;
#line 1407 "clanhf.f"
			aa = a[i__2].r;
/*                    L(k-1+i,k-1+i) */
#line 1409 "clanhf.f"
			if (aa != 0.) {
#line 1410 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1411 "clanhf.f"
				d__1 = scale / aa;
#line 1411 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1412 "clanhf.f"
				scale = aa;
#line 1413 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1414 "clanhf.f"
				d__1 = aa / scale;
#line 1414 "clanhf.f"
				s += d__1 * d__1;
#line 1415 "clanhf.f"
			    }
#line 1416 "clanhf.f"
			}
#line 1417 "clanhf.f"
			i__2 = l + 1;
#line 1417 "clanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1419 "clanhf.f"
			if (aa != 0.) {
#line 1420 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1421 "clanhf.f"
				d__1 = scale / aa;
#line 1421 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1422 "clanhf.f"
				scale = aa;
#line 1423 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1424 "clanhf.f"
				d__1 = aa / scale;
#line 1424 "clanhf.f"
				s += d__1 * d__1;
#line 1425 "clanhf.f"
			    }
#line 1426 "clanhf.f"
			}
#line 1427 "clanhf.f"
			l = l + lda + 1;
#line 1428 "clanhf.f"
		    }
#line 1429 "clanhf.f"
		}
#line 1430 "clanhf.f"
	    } else {
/*              A is xpose */
#line 1432 "clanhf.f"
		if (ilu == 0) {
/*                 A**H is upper */
#line 1434 "clanhf.f"
		    i__1 = k - 1;
#line 1434 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1435 "clanhf.f"
			classq_(&j, &a[(k + 1 + j) * lda], &c__1, &scale, &s);
/*                 U at A(0,k+1) */
#line 1437 "clanhf.f"
		    }
#line 1438 "clanhf.f"
		    i__1 = k - 1;
#line 1438 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1439 "clanhf.f"
			classq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                 k by k rect. at A(0,0) */
#line 1441 "clanhf.f"
		    }
#line 1442 "clanhf.f"
		    i__1 = k - 2;
#line 1442 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1443 "clanhf.f"
			i__2 = k - j - 1;
#line 1443 "clanhf.f"
			classq_(&i__2, &a[j + 1 + (j + k) * lda], &c__1, &
				scale, &s);
/*                 L at A(0,k) */
#line 1446 "clanhf.f"
		    }
#line 1447 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1449 "clanhf.f"
		    l = k * lda;
/*                 -> U(k,k) at A(0,k) */
#line 1451 "clanhf.f"
		    i__1 = l;
#line 1451 "clanhf.f"
		    aa = a[i__1].r;
/*                 U(k,k) */
#line 1453 "clanhf.f"
		    if (aa != 0.) {
#line 1454 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1455 "clanhf.f"
			    d__1 = scale / aa;
#line 1455 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1456 "clanhf.f"
			    scale = aa;
#line 1457 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1458 "clanhf.f"
			    d__1 = aa / scale;
#line 1458 "clanhf.f"
			    s += d__1 * d__1;
#line 1459 "clanhf.f"
			}
#line 1460 "clanhf.f"
		    }
#line 1461 "clanhf.f"
		    l += lda;
/*                 -> U(0,0) at A(0,k+1) */
#line 1463 "clanhf.f"
		    i__1 = *n - 1;
#line 1463 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1464 "clanhf.f"
			i__2 = l;
#line 1464 "clanhf.f"
			aa = a[i__2].r;
/*                    -> U(j-k-1,j-k-1) */
#line 1466 "clanhf.f"
			if (aa != 0.) {
#line 1467 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1468 "clanhf.f"
				d__1 = scale / aa;
#line 1468 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1469 "clanhf.f"
				scale = aa;
#line 1470 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1471 "clanhf.f"
				d__1 = aa / scale;
#line 1471 "clanhf.f"
				s += d__1 * d__1;
#line 1472 "clanhf.f"
			    }
#line 1473 "clanhf.f"
			}
#line 1474 "clanhf.f"
			i__2 = l + 1;
#line 1474 "clanhf.f"
			aa = a[i__2].r;
/*                    -> U(j,j) */
#line 1476 "clanhf.f"
			if (aa != 0.) {
#line 1477 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1478 "clanhf.f"
				d__1 = scale / aa;
#line 1478 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1479 "clanhf.f"
				scale = aa;
#line 1480 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1481 "clanhf.f"
				d__1 = aa / scale;
#line 1481 "clanhf.f"
				s += d__1 * d__1;
#line 1482 "clanhf.f"
			    }
#line 1483 "clanhf.f"
			}
#line 1484 "clanhf.f"
			l = l + lda + 1;
#line 1485 "clanhf.f"
		    }
/*                 L=k-1+n*lda */
/*                 -> U(k-1,k-1) at A(k-1,n) */
#line 1488 "clanhf.f"
		    i__1 = l;
#line 1488 "clanhf.f"
		    aa = a[i__1].r;
/*                 U(k,k) */
#line 1490 "clanhf.f"
		    if (aa != 0.) {
#line 1491 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1492 "clanhf.f"
			    d__1 = scale / aa;
#line 1492 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1493 "clanhf.f"
			    scale = aa;
#line 1494 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1495 "clanhf.f"
			    d__1 = aa / scale;
#line 1495 "clanhf.f"
			    s += d__1 * d__1;
#line 1496 "clanhf.f"
			}
#line 1497 "clanhf.f"
		    }
#line 1498 "clanhf.f"
		} else {
/*                 A**H is lower */
#line 1500 "clanhf.f"
		    i__1 = k - 1;
#line 1500 "clanhf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 1501 "clanhf.f"
			classq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                 U at A(0,1) */
#line 1503 "clanhf.f"
		    }
#line 1504 "clanhf.f"
		    i__1 = *n;
#line 1504 "clanhf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 1505 "clanhf.f"
			classq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                 k by k rect. at A(0,k+1) */
#line 1507 "clanhf.f"
		    }
#line 1508 "clanhf.f"
		    i__1 = k - 2;
#line 1508 "clanhf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 1509 "clanhf.f"
			i__2 = k - j - 1;
#line 1509 "clanhf.f"
			classq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                 L at A(0,0) */
#line 1511 "clanhf.f"
		    }
#line 1512 "clanhf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 1514 "clanhf.f"
		    l = 0;
/*                 -> L(k,k) at A(0,0) */
#line 1516 "clanhf.f"
		    i__1 = l;
#line 1516 "clanhf.f"
		    aa = a[i__1].r;
/*                 L(k,k) at A(0,0) */
#line 1518 "clanhf.f"
		    if (aa != 0.) {
#line 1519 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1520 "clanhf.f"
			    d__1 = scale / aa;
#line 1520 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1521 "clanhf.f"
			    scale = aa;
#line 1522 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1523 "clanhf.f"
			    d__1 = aa / scale;
#line 1523 "clanhf.f"
			    s += d__1 * d__1;
#line 1524 "clanhf.f"
			}
#line 1525 "clanhf.f"
		    }
#line 1526 "clanhf.f"
		    l = lda;
/*                 -> L(0,0) at A(0,1) */
#line 1528 "clanhf.f"
		    i__1 = k - 2;
#line 1528 "clanhf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1529 "clanhf.f"
			i__2 = l;
#line 1529 "clanhf.f"
			aa = a[i__2].r;
/*                    L(i,i) */
#line 1531 "clanhf.f"
			if (aa != 0.) {
#line 1532 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1533 "clanhf.f"
				d__1 = scale / aa;
#line 1533 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1534 "clanhf.f"
				scale = aa;
#line 1535 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1536 "clanhf.f"
				d__1 = aa / scale;
#line 1536 "clanhf.f"
				s += d__1 * d__1;
#line 1537 "clanhf.f"
			    }
#line 1538 "clanhf.f"
			}
#line 1539 "clanhf.f"
			i__2 = l + 1;
#line 1539 "clanhf.f"
			aa = a[i__2].r;
/*                    L(k+i+1,k+i+1) */
#line 1541 "clanhf.f"
			if (aa != 0.) {
#line 1542 "clanhf.f"
			    if (scale < aa) {
/* Computing 2nd power */
#line 1543 "clanhf.f"
				d__1 = scale / aa;
#line 1543 "clanhf.f"
				s = s * (d__1 * d__1) + 1.;
#line 1544 "clanhf.f"
				scale = aa;
#line 1545 "clanhf.f"
			    } else {
/* Computing 2nd power */
#line 1546 "clanhf.f"
				d__1 = aa / scale;
#line 1546 "clanhf.f"
				s += d__1 * d__1;
#line 1547 "clanhf.f"
			    }
#line 1548 "clanhf.f"
			}
#line 1549 "clanhf.f"
			l = l + lda + 1;
#line 1550 "clanhf.f"
		    }
/*                 L-> k - 1 + k*lda or L(k-1,k-1) at A(k-1,k) */
#line 1552 "clanhf.f"
		    i__1 = l;
#line 1552 "clanhf.f"
		    aa = a[i__1].r;
/*                 L(k-1,k-1) at A(k-1,k) */
#line 1554 "clanhf.f"
		    if (aa != 0.) {
#line 1555 "clanhf.f"
			if (scale < aa) {
/* Computing 2nd power */
#line 1556 "clanhf.f"
			    d__1 = scale / aa;
#line 1556 "clanhf.f"
			    s = s * (d__1 * d__1) + 1.;
#line 1557 "clanhf.f"
			    scale = aa;
#line 1558 "clanhf.f"
			} else {
/* Computing 2nd power */
#line 1559 "clanhf.f"
			    d__1 = aa / scale;
#line 1559 "clanhf.f"
			    s += d__1 * d__1;
#line 1560 "clanhf.f"
			}
#line 1561 "clanhf.f"
		    }
#line 1562 "clanhf.f"
		}
#line 1563 "clanhf.f"
	    }
#line 1564 "clanhf.f"
	}
#line 1565 "clanhf.f"
	value = scale * sqrt(s);
#line 1566 "clanhf.f"
    }

#line 1568 "clanhf.f"
    ret_val = value;
#line 1569 "clanhf.f"
    return ret_val;

/*     End of CLANHF */

} /* clanhf_ */

