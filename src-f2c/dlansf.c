#line 1 "dlansf.f"
/* dlansf.f -- translated by f2c (version 20100827).
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

#line 1 "dlansf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANSF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a symmetric matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANSF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANSF( NORM, TRANSR, UPLO, N, A, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, TRANSR, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( 0: * ), WORK( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSF returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > real symmetric matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \return DLANSF */
/* > \verbatim */
/* > */
/* >    DLANSF = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          NORM is CHARACTER*1 */
/* >          Specifies the value to be returned in DLANSF as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          Specifies whether the RFP format of A is normal or */
/* >          transposed format. */
/* >          = 'N':  RFP format is Normal; */
/* >          = 'T':  RFP format is Transpose. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the RFP matrix A came from */
/* >           an upper or lower triangular matrix as follows: */
/* >           = 'U': RFP A came from an upper triangular matrix; */
/* >           = 'L': RFP A came from a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. N >= 0. When N = 0, DLANSF is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( N*(N+1)/2 ); */
/* >          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') */
/* >          part of the symmetric matrix A stored in RFP format. See the */
/* >          "Notes" below for more details. */
/* >          Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise, */
/* >          WORK is not referenced. */
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
doublereal dlansf_(char *norm, char *transr, char *uplo, integer *n, 
	doublereal *a, doublereal *work, ftnlen norm_len, ftnlen transr_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

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
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


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

#line 247 "dlansf.f"
    if (*n == 0) {
#line 248 "dlansf.f"
	ret_val = 0.;
#line 249 "dlansf.f"
	return ret_val;
#line 250 "dlansf.f"
    } else if (*n == 1) {
#line 251 "dlansf.f"
	ret_val = abs(a[0]);
#line 252 "dlansf.f"
	return ret_val;
#line 253 "dlansf.f"
    }

/*     set noe = 1 if n is odd. if n is even set noe=0 */

#line 257 "dlansf.f"
    noe = 1;
#line 258 "dlansf.f"
    if (*n % 2 == 0) {
#line 258 "dlansf.f"
	noe = 0;
#line 258 "dlansf.f"
    }

/*     set ifm = 0 when form='T or 't' and 1 otherwise */

#line 263 "dlansf.f"
    ifm = 1;
#line 264 "dlansf.f"
    if (lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 264 "dlansf.f"
	ifm = 0;
#line 264 "dlansf.f"
    }

/*     set ilu = 0 when uplo='U or 'u' and 1 otherwise */

#line 269 "dlansf.f"
    ilu = 1;
#line 270 "dlansf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 270 "dlansf.f"
	ilu = 0;
#line 270 "dlansf.f"
    }

/*     set lda = (n+1)/2 when ifm = 0 */
/*     set lda = n when ifm = 1 and noe = 1 */
/*     set lda = n+1 when ifm = 1 and noe = 0 */

#line 277 "dlansf.f"
    if (ifm == 1) {
#line 278 "dlansf.f"
	if (noe == 1) {
#line 279 "dlansf.f"
	    lda = *n;
#line 280 "dlansf.f"
	} else {
/*           noe=0 */
#line 282 "dlansf.f"
	    lda = *n + 1;
#line 283 "dlansf.f"
	}
#line 284 "dlansf.f"
    } else {
/*        ifm=0 */
#line 286 "dlansf.f"
	lda = (*n + 1) / 2;
#line 287 "dlansf.f"
    }

#line 289 "dlansf.f"
    if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*       Find max(abs(A(i,j))). */

#line 293 "dlansf.f"
	k = (*n + 1) / 2;
#line 294 "dlansf.f"
	value = 0.;
#line 295 "dlansf.f"
	if (noe == 1) {
/*           n is odd */
#line 297 "dlansf.f"
	    if (ifm == 1) {
/*           A is n by k */
#line 299 "dlansf.f"
		i__1 = k - 1;
#line 299 "dlansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 300 "dlansf.f"
		    i__2 = *n - 1;
#line 300 "dlansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 301 "dlansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 302 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 302 "dlansf.f"
			    value = temp;
#line 302 "dlansf.f"
			}
#line 304 "dlansf.f"
		    }
#line 305 "dlansf.f"
		}
#line 306 "dlansf.f"
	    } else {
/*              xpose case; A is k by n */
#line 308 "dlansf.f"
		i__1 = *n - 1;
#line 308 "dlansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 309 "dlansf.f"
		    i__2 = k - 1;
#line 309 "dlansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 310 "dlansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 311 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 311 "dlansf.f"
			    value = temp;
#line 311 "dlansf.f"
			}
#line 313 "dlansf.f"
		    }
#line 314 "dlansf.f"
		}
#line 315 "dlansf.f"
	    }
#line 316 "dlansf.f"
	} else {
/*           n is even */
#line 318 "dlansf.f"
	    if (ifm == 1) {
/*              A is n+1 by k */
#line 320 "dlansf.f"
		i__1 = k - 1;
#line 320 "dlansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 321 "dlansf.f"
		    i__2 = *n;
#line 321 "dlansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 322 "dlansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 323 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 323 "dlansf.f"
			    value = temp;
#line 323 "dlansf.f"
			}
#line 325 "dlansf.f"
		    }
#line 326 "dlansf.f"
		}
#line 327 "dlansf.f"
	    } else {
/*              xpose case; A is k by n+1 */
#line 329 "dlansf.f"
		i__1 = *n;
#line 329 "dlansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 330 "dlansf.f"
		    i__2 = k - 1;
#line 330 "dlansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 331 "dlansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 332 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 332 "dlansf.f"
			    value = temp;
#line 332 "dlansf.f"
			}
#line 334 "dlansf.f"
		    }
#line 335 "dlansf.f"
		}
#line 336 "dlansf.f"
	    }
#line 337 "dlansf.f"
	}
#line 338 "dlansf.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 343 "dlansf.f"
	if (ifm == 1) {
#line 344 "dlansf.f"
	    k = *n / 2;
#line 345 "dlansf.f"
	    if (noe == 1) {
/*              n is odd */
#line 347 "dlansf.f"
		if (ilu == 0) {
#line 348 "dlansf.f"
		    i__1 = k - 1;
#line 348 "dlansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 349 "dlansf.f"
			work[i__] = 0.;
#line 350 "dlansf.f"
		    }
#line 351 "dlansf.f"
		    i__1 = k;
#line 351 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 352 "dlansf.f"
			s = 0.;
#line 353 "dlansf.f"
			i__2 = k + j - 1;
#line 353 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 354 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(i,j+k) */
#line 356 "dlansf.f"
			    s += aa;
#line 357 "dlansf.f"
			    work[i__] += aa;
#line 358 "dlansf.f"
			}
#line 359 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 361 "dlansf.f"
			work[j + k] = s + aa;
#line 362 "dlansf.f"
			if (i__ == k + k) {
#line 362 "dlansf.f"
			    goto L10;
#line 362 "dlansf.f"
			}
#line 364 "dlansf.f"
			++i__;
#line 365 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 367 "dlansf.f"
			work[j] += aa;
#line 368 "dlansf.f"
			s = 0.;
#line 369 "dlansf.f"
			i__2 = k - 1;
#line 369 "dlansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 370 "dlansf.f"
			    ++i__;
#line 371 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 373 "dlansf.f"
			    s += aa;
#line 374 "dlansf.f"
			    work[l] += aa;
#line 375 "dlansf.f"
			}
#line 376 "dlansf.f"
			work[j] += s;
#line 377 "dlansf.f"
		    }
#line 378 "dlansf.f"
L10:
#line 379 "dlansf.f"
		    value = work[0];
#line 380 "dlansf.f"
		    i__1 = *n - 1;
#line 380 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 381 "dlansf.f"
			temp = work[i__];
#line 382 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 382 "dlansf.f"
			    value = temp;
#line 382 "dlansf.f"
			}
#line 384 "dlansf.f"
		    }
#line 385 "dlansf.f"
		} else {
/*                 ilu = 1 */
#line 387 "dlansf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 389 "dlansf.f"
		    i__1 = *n - 1;
#line 389 "dlansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 390 "dlansf.f"
			work[i__] = 0.;
#line 391 "dlansf.f"
		    }
#line 392 "dlansf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 393 "dlansf.f"
			s = 0.;
#line 394 "dlansf.f"
			i__1 = j - 2;
#line 394 "dlansf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 395 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,i+k) */
#line 397 "dlansf.f"
			    s += aa;
#line 398 "dlansf.f"
			    work[i__ + k] += aa;
#line 399 "dlansf.f"
			}
#line 400 "dlansf.f"
			if (j > 0) {
#line 401 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,j+k) */
#line 403 "dlansf.f"
			    s += aa;
#line 404 "dlansf.f"
			    work[i__ + k] += s;
/*                       i=j */
#line 406 "dlansf.f"
			    ++i__;
#line 407 "dlansf.f"
			}
#line 408 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 410 "dlansf.f"
			work[j] = aa;
#line 411 "dlansf.f"
			s = 0.;
#line 412 "dlansf.f"
			i__1 = *n - 1;
#line 412 "dlansf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 413 "dlansf.f"
			    ++i__;
#line 414 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 416 "dlansf.f"
			    s += aa;
#line 417 "dlansf.f"
			    work[l] += aa;
#line 418 "dlansf.f"
			}
#line 419 "dlansf.f"
			work[j] += s;
#line 420 "dlansf.f"
		    }
#line 421 "dlansf.f"
		    value = work[0];
#line 422 "dlansf.f"
		    i__1 = *n - 1;
#line 422 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 423 "dlansf.f"
			temp = work[i__];
#line 424 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 424 "dlansf.f"
			    value = temp;
#line 424 "dlansf.f"
			}
#line 426 "dlansf.f"
		    }
#line 427 "dlansf.f"
		}
#line 428 "dlansf.f"
	    } else {
/*              n is even */
#line 430 "dlansf.f"
		if (ilu == 0) {
#line 431 "dlansf.f"
		    i__1 = k - 1;
#line 431 "dlansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 432 "dlansf.f"
			work[i__] = 0.;
#line 433 "dlansf.f"
		    }
#line 434 "dlansf.f"
		    i__1 = k - 1;
#line 434 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 435 "dlansf.f"
			s = 0.;
#line 436 "dlansf.f"
			i__2 = k + j - 1;
#line 436 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 437 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(i,j+k) */
#line 439 "dlansf.f"
			    s += aa;
#line 440 "dlansf.f"
			    work[i__] += aa;
#line 441 "dlansf.f"
			}
#line 442 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 444 "dlansf.f"
			work[j + k] = s + aa;
#line 445 "dlansf.f"
			++i__;
#line 446 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 448 "dlansf.f"
			work[j] += aa;
#line 449 "dlansf.f"
			s = 0.;
#line 450 "dlansf.f"
			i__2 = k - 1;
#line 450 "dlansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 451 "dlansf.f"
			    ++i__;
#line 452 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 454 "dlansf.f"
			    s += aa;
#line 455 "dlansf.f"
			    work[l] += aa;
#line 456 "dlansf.f"
			}
#line 457 "dlansf.f"
			work[j] += s;
#line 458 "dlansf.f"
		    }
#line 459 "dlansf.f"
		    value = work[0];
#line 460 "dlansf.f"
		    i__1 = *n - 1;
#line 460 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 461 "dlansf.f"
			temp = work[i__];
#line 462 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 462 "dlansf.f"
			    value = temp;
#line 462 "dlansf.f"
			}
#line 464 "dlansf.f"
		    }
#line 465 "dlansf.f"
		} else {
/*                 ilu = 1 */
#line 467 "dlansf.f"
		    i__1 = *n - 1;
#line 467 "dlansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 468 "dlansf.f"
			work[i__] = 0.;
#line 469 "dlansf.f"
		    }
#line 470 "dlansf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 471 "dlansf.f"
			s = 0.;
#line 472 "dlansf.f"
			i__1 = j - 1;
#line 472 "dlansf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 473 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,i+k) */
#line 475 "dlansf.f"
			    s += aa;
#line 476 "dlansf.f"
			    work[i__ + k] += aa;
#line 477 "dlansf.f"
			}
#line 478 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 480 "dlansf.f"
			s += aa;
#line 481 "dlansf.f"
			work[i__ + k] += s;
/*                    i=j */
#line 483 "dlansf.f"
			++i__;
#line 484 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 486 "dlansf.f"
			work[j] = aa;
#line 487 "dlansf.f"
			s = 0.;
#line 488 "dlansf.f"
			i__1 = *n - 1;
#line 488 "dlansf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 489 "dlansf.f"
			    ++i__;
#line 490 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 492 "dlansf.f"
			    s += aa;
#line 493 "dlansf.f"
			    work[l] += aa;
#line 494 "dlansf.f"
			}
#line 495 "dlansf.f"
			work[j] += s;
#line 496 "dlansf.f"
		    }
#line 497 "dlansf.f"
		    value = work[0];
#line 498 "dlansf.f"
		    i__1 = *n - 1;
#line 498 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 499 "dlansf.f"
			temp = work[i__];
#line 500 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 500 "dlansf.f"
			    value = temp;
#line 500 "dlansf.f"
			}
#line 502 "dlansf.f"
		    }
#line 503 "dlansf.f"
		}
#line 504 "dlansf.f"
	    }
#line 505 "dlansf.f"
	} else {
/*           ifm=0 */
#line 507 "dlansf.f"
	    k = *n / 2;
#line 508 "dlansf.f"
	    if (noe == 1) {
/*              n is odd */
#line 510 "dlansf.f"
		if (ilu == 0) {
#line 511 "dlansf.f"
		    n1 = k;
/*                 n/2 */
#line 513 "dlansf.f"
		    ++k;
/*                 k is the row size and lda */
#line 515 "dlansf.f"
		    i__1 = *n - 1;
#line 515 "dlansf.f"
		    for (i__ = n1; i__ <= i__1; ++i__) {
#line 516 "dlansf.f"
			work[i__] = 0.;
#line 517 "dlansf.f"
		    }
#line 518 "dlansf.f"
		    i__1 = n1 - 1;
#line 518 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 519 "dlansf.f"
			s = 0.;
#line 520 "dlansf.f"
			i__2 = k - 1;
#line 520 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 521 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,n1+i) */
#line 523 "dlansf.f"
			    work[i__ + n1] += aa;
#line 524 "dlansf.f"
			    s += aa;
#line 525 "dlansf.f"
			}
#line 526 "dlansf.f"
			work[j] = s;
#line 527 "dlansf.f"
		    }
/*                 j=n1=k-1 is special */
#line 529 "dlansf.f"
		    s = (d__1 = a[j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 531 "dlansf.f"
		    i__1 = k - 1;
#line 531 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 532 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k-1,i+n1) */
#line 534 "dlansf.f"
			work[i__ + n1] += aa;
#line 535 "dlansf.f"
			s += aa;
#line 536 "dlansf.f"
		    }
#line 537 "dlansf.f"
		    work[j] += s;
#line 538 "dlansf.f"
		    i__1 = *n - 1;
#line 538 "dlansf.f"
		    for (j = k; j <= i__1; ++j) {
#line 539 "dlansf.f"
			s = 0.;
#line 540 "dlansf.f"
			i__2 = j - k - 1;
#line 540 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 541 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(i,j-k) */
#line 543 "dlansf.f"
			    work[i__] += aa;
#line 544 "dlansf.f"
			    s += aa;
#line 545 "dlansf.f"
			}
/*                    i=j-k */
#line 547 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j-k,j-k) */
#line 549 "dlansf.f"
			s += aa;
#line 550 "dlansf.f"
			work[j - k] += s;
#line 551 "dlansf.f"
			++i__;
#line 552 "dlansf.f"
			s = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j,j) */
#line 554 "dlansf.f"
			i__2 = *n - 1;
#line 554 "dlansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 555 "dlansf.f"
			    ++i__;
#line 556 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,l) */
#line 558 "dlansf.f"
			    work[l] += aa;
#line 559 "dlansf.f"
			    s += aa;
#line 560 "dlansf.f"
			}
#line 561 "dlansf.f"
			work[j] += s;
#line 562 "dlansf.f"
		    }
#line 563 "dlansf.f"
		    value = work[0];
#line 564 "dlansf.f"
		    i__1 = *n - 1;
#line 564 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 565 "dlansf.f"
			temp = work[i__];
#line 566 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 566 "dlansf.f"
			    value = temp;
#line 566 "dlansf.f"
			}
#line 568 "dlansf.f"
		    }
#line 569 "dlansf.f"
		} else {
/*                 ilu=1 */
#line 571 "dlansf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 573 "dlansf.f"
		    i__1 = *n - 1;
#line 573 "dlansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 574 "dlansf.f"
			work[i__] = 0.;
#line 575 "dlansf.f"
		    }
#line 576 "dlansf.f"
		    i__1 = k - 2;
#line 576 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
/*                    process */
#line 578 "dlansf.f"
			s = 0.;
#line 579 "dlansf.f"
			i__2 = j - 1;
#line 579 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 580 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i) */
#line 582 "dlansf.f"
			    work[i__] += aa;
#line 583 "dlansf.f"
			    s += aa;
#line 584 "dlansf.f"
			}
#line 585 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    i=j so process of A(j,j) */
#line 587 "dlansf.f"
			s += aa;
#line 588 "dlansf.f"
			work[j] = s;
/*                    is initialised here */
#line 590 "dlansf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 592 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
#line 593 "dlansf.f"
			s = aa;
#line 594 "dlansf.f"
			i__2 = *n - 1;
#line 594 "dlansf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 595 "dlansf.f"
			    ++i__;
#line 596 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(l,k+j) */
#line 598 "dlansf.f"
			    s += aa;
#line 599 "dlansf.f"
			    work[l] += aa;
#line 600 "dlansf.f"
			}
#line 601 "dlansf.f"
			work[k + j] += s;
#line 602 "dlansf.f"
		    }
/*                 j=k-1 is special :process col A(k-1,0:k-1) */
#line 604 "dlansf.f"
		    s = 0.;
#line 605 "dlansf.f"
		    i__1 = k - 2;
#line 605 "dlansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 606 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,i) */
#line 608 "dlansf.f"
			work[i__] += aa;
#line 609 "dlansf.f"
			s += aa;
#line 610 "dlansf.f"
		    }
/*                 i=k-1 */
#line 612 "dlansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 614 "dlansf.f"
		    s += aa;
#line 615 "dlansf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 617 "dlansf.f"
		    i__1 = *n - 1;
#line 617 "dlansf.f"
		    for (j = k; j <= i__1; ++j) {
/*                    process col j of A = A(j,0:k-1) */
#line 619 "dlansf.f"
			s = 0.;
#line 620 "dlansf.f"
			i__2 = k - 1;
#line 620 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 621 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i) */
#line 623 "dlansf.f"
			    work[i__] += aa;
#line 624 "dlansf.f"
			    s += aa;
#line 625 "dlansf.f"
			}
#line 626 "dlansf.f"
			work[j] += s;
#line 627 "dlansf.f"
		    }
#line 628 "dlansf.f"
		    value = work[0];
#line 629 "dlansf.f"
		    i__1 = *n - 1;
#line 629 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 630 "dlansf.f"
			temp = work[i__];
#line 631 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 631 "dlansf.f"
			    value = temp;
#line 631 "dlansf.f"
			}
#line 633 "dlansf.f"
		    }
#line 634 "dlansf.f"
		}
#line 635 "dlansf.f"
	    } else {
/*              n is even */
#line 637 "dlansf.f"
		if (ilu == 0) {
#line 638 "dlansf.f"
		    i__1 = *n - 1;
#line 638 "dlansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 639 "dlansf.f"
			work[i__] = 0.;
#line 640 "dlansf.f"
		    }
#line 641 "dlansf.f"
		    i__1 = k - 1;
#line 641 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 642 "dlansf.f"
			s = 0.;
#line 643 "dlansf.f"
			i__2 = k - 1;
#line 643 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 644 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i+k) */
#line 646 "dlansf.f"
			    work[i__ + k] += aa;
#line 647 "dlansf.f"
			    s += aa;
#line 648 "dlansf.f"
			}
#line 649 "dlansf.f"
			work[j] = s;
#line 650 "dlansf.f"
		    }
/*                 j=k */
#line 652 "dlansf.f"
		    aa = (d__1 = a[j * lda], abs(d__1));
/*                 A(k,k) */
#line 654 "dlansf.f"
		    s = aa;
#line 655 "dlansf.f"
		    i__1 = k - 1;
#line 655 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 656 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,k+i) */
#line 658 "dlansf.f"
			work[i__ + k] += aa;
#line 659 "dlansf.f"
			s += aa;
#line 660 "dlansf.f"
		    }
#line 661 "dlansf.f"
		    work[j] += s;
#line 662 "dlansf.f"
		    i__1 = *n - 1;
#line 662 "dlansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 663 "dlansf.f"
			s = 0.;
#line 664 "dlansf.f"
			i__2 = j - 2 - k;
#line 664 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 665 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(i,j-k-1) */
#line 667 "dlansf.f"
			    work[i__] += aa;
#line 668 "dlansf.f"
			    s += aa;
#line 669 "dlansf.f"
			}
/*                     i=j-1-k */
#line 671 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j-k-1,j-k-1) */
#line 673 "dlansf.f"
			s += aa;
#line 674 "dlansf.f"
			work[j - k - 1] += s;
#line 675 "dlansf.f"
			++i__;
#line 676 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j,j) */
#line 678 "dlansf.f"
			s = aa;
#line 679 "dlansf.f"
			i__2 = *n - 1;
#line 679 "dlansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 680 "dlansf.f"
			    ++i__;
#line 681 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,l) */
#line 683 "dlansf.f"
			    work[l] += aa;
#line 684 "dlansf.f"
			    s += aa;
#line 685 "dlansf.f"
			}
#line 686 "dlansf.f"
			work[j] += s;
#line 687 "dlansf.f"
		    }
/*                 j=n */
#line 689 "dlansf.f"
		    s = 0.;
#line 690 "dlansf.f"
		    i__1 = k - 2;
#line 690 "dlansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 691 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(i,k-1) */
#line 693 "dlansf.f"
			work[i__] += aa;
#line 694 "dlansf.f"
			s += aa;
#line 695 "dlansf.f"
		    }
/*                 i=k-1 */
#line 697 "dlansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 699 "dlansf.f"
		    s += aa;
#line 700 "dlansf.f"
		    work[i__] += s;
#line 701 "dlansf.f"
		    value = work[0];
#line 702 "dlansf.f"
		    i__1 = *n - 1;
#line 702 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 703 "dlansf.f"
			temp = work[i__];
#line 704 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 704 "dlansf.f"
			    value = temp;
#line 704 "dlansf.f"
			}
#line 706 "dlansf.f"
		    }
#line 707 "dlansf.f"
		} else {
/*                 ilu=1 */
#line 709 "dlansf.f"
		    i__1 = *n - 1;
#line 709 "dlansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 710 "dlansf.f"
			work[i__] = 0.;
#line 711 "dlansf.f"
		    }
/*                 j=0 is special :process col A(k:n-1,k) */
#line 713 "dlansf.f"
		    s = abs(a[0]);
/*                 A(k,k) */
#line 715 "dlansf.f"
		    i__1 = k - 1;
#line 715 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 716 "dlansf.f"
			aa = (d__1 = a[i__], abs(d__1));
/*                    A(k+i,k) */
#line 718 "dlansf.f"
			work[i__ + k] += aa;
#line 719 "dlansf.f"
			s += aa;
#line 720 "dlansf.f"
		    }
#line 721 "dlansf.f"
		    work[k] += s;
#line 722 "dlansf.f"
		    i__1 = k - 1;
#line 722 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
/*                    process */
#line 724 "dlansf.f"
			s = 0.;
#line 725 "dlansf.f"
			i__2 = j - 2;
#line 725 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 726 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j-1,i) */
#line 728 "dlansf.f"
			    work[i__] += aa;
#line 729 "dlansf.f"
			    s += aa;
#line 730 "dlansf.f"
			}
#line 731 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    i=j-1 so process of A(j-1,j-1) */
#line 733 "dlansf.f"
			s += aa;
#line 734 "dlansf.f"
			work[j - 1] = s;
/*                    is initialised here */
#line 736 "dlansf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 738 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
#line 739 "dlansf.f"
			s = aa;
#line 740 "dlansf.f"
			i__2 = *n - 1;
#line 740 "dlansf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 741 "dlansf.f"
			    ++i__;
#line 742 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(l,k+j) */
#line 744 "dlansf.f"
			    s += aa;
#line 745 "dlansf.f"
			    work[l] += aa;
#line 746 "dlansf.f"
			}
#line 747 "dlansf.f"
			work[k + j] += s;
#line 748 "dlansf.f"
		    }
/*                 j=k is special :process col A(k,0:k-1) */
#line 750 "dlansf.f"
		    s = 0.;
#line 751 "dlansf.f"
		    i__1 = k - 2;
#line 751 "dlansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 752 "dlansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,i) */
#line 754 "dlansf.f"
			work[i__] += aa;
#line 755 "dlansf.f"
			s += aa;
#line 756 "dlansf.f"
		    }
/*                 i=k-1 */
#line 758 "dlansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 760 "dlansf.f"
		    s += aa;
#line 761 "dlansf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 763 "dlansf.f"
		    i__1 = *n;
#line 763 "dlansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
/*                    process col j-1 of A = A(j-1,0:k-1) */
#line 765 "dlansf.f"
			s = 0.;
#line 766 "dlansf.f"
			i__2 = k - 1;
#line 766 "dlansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 767 "dlansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j-1,i) */
#line 769 "dlansf.f"
			    work[i__] += aa;
#line 770 "dlansf.f"
			    s += aa;
#line 771 "dlansf.f"
			}
#line 772 "dlansf.f"
			work[j - 1] += s;
#line 773 "dlansf.f"
		    }
#line 774 "dlansf.f"
		    value = work[0];
#line 775 "dlansf.f"
		    i__1 = *n - 1;
#line 775 "dlansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 776 "dlansf.f"
			temp = work[i__];
#line 777 "dlansf.f"
			if (value < temp || disnan_(&temp)) {
#line 777 "dlansf.f"
			    value = temp;
#line 777 "dlansf.f"
			}
#line 779 "dlansf.f"
		    }
#line 780 "dlansf.f"
		}
#line 781 "dlansf.f"
	    }
#line 782 "dlansf.f"
	}
#line 783 "dlansf.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*       Find normF(A). */

#line 787 "dlansf.f"
	k = (*n + 1) / 2;
#line 788 "dlansf.f"
	scale = 0.;
#line 789 "dlansf.f"
	s = 1.;
#line 790 "dlansf.f"
	if (noe == 1) {
/*           n is odd */
#line 792 "dlansf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 794 "dlansf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 796 "dlansf.f"
		    i__1 = k - 3;
#line 796 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 797 "dlansf.f"
			i__2 = k - j - 2;
#line 797 "dlansf.f"
			dlassq_(&i__2, &a[k + j + 1 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k,0) */
#line 799 "dlansf.f"
		    }
#line 800 "dlansf.f"
		    i__1 = k - 1;
#line 800 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 801 "dlansf.f"
			i__2 = k + j - 1;
#line 801 "dlansf.f"
			dlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 803 "dlansf.f"
		    }
#line 804 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 806 "dlansf.f"
		    i__1 = k - 1;
#line 806 "dlansf.f"
		    i__2 = lda + 1;
#line 806 "dlansf.f"
		    dlassq_(&i__1, &a[k], &i__2, &scale, &s);
/*                 tri L at A(k,0) */
#line 808 "dlansf.f"
		    i__1 = lda + 1;
#line 808 "dlansf.f"
		    dlassq_(&k, &a[k - 1], &i__1, &scale, &s);
/*                 tri U at A(k-1,0) */
#line 810 "dlansf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 812 "dlansf.f"
		    i__1 = k - 1;
#line 812 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 813 "dlansf.f"
			i__2 = *n - j - 1;
#line 813 "dlansf.f"
			dlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(0,0) */
#line 815 "dlansf.f"
		    }
#line 816 "dlansf.f"
		    i__1 = k - 2;
#line 816 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 817 "dlansf.f"
			dlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 819 "dlansf.f"
		    }
#line 820 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 822 "dlansf.f"
		    i__1 = lda + 1;
#line 822 "dlansf.f"
		    dlassq_(&k, a, &i__1, &scale, &s);
/*                 tri L at A(0,0) */
#line 824 "dlansf.f"
		    i__1 = k - 1;
#line 824 "dlansf.f"
		    i__2 = lda + 1;
#line 824 "dlansf.f"
		    dlassq_(&i__1, &a[lda], &i__2, &scale, &s);
/*                 tri U at A(0,1) */
#line 826 "dlansf.f"
		}
#line 827 "dlansf.f"
	    } else {
/*              A is xpose */
#line 829 "dlansf.f"
		if (ilu == 0) {
/*                 A**T is upper */
#line 831 "dlansf.f"
		    i__1 = k - 2;
#line 831 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 832 "dlansf.f"
			dlassq_(&j, &a[(k + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k) */
#line 834 "dlansf.f"
		    }
#line 835 "dlansf.f"
		    i__1 = k - 2;
#line 835 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 836 "dlansf.f"
			dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,0) */
#line 838 "dlansf.f"
		    }
#line 839 "dlansf.f"
		    i__1 = k - 2;
#line 839 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 840 "dlansf.f"
			i__2 = k - j - 1;
#line 840 "dlansf.f"
			dlassq_(&i__2, &a[j + 1 + (j + k - 1) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k-1) */
#line 843 "dlansf.f"
		    }
#line 844 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 846 "dlansf.f"
		    i__1 = k - 1;
#line 846 "dlansf.f"
		    i__2 = lda + 1;
#line 846 "dlansf.f"
		    dlassq_(&i__1, &a[k * lda], &i__2, &scale, &s);
/*                 tri U at A(0,k) */
#line 848 "dlansf.f"
		    i__1 = lda + 1;
#line 848 "dlansf.f"
		    dlassq_(&k, &a[(k - 1) * lda], &i__1, &scale, &s);
/*                 tri L at A(0,k-1) */
#line 850 "dlansf.f"
		} else {
/*                 A**T is lower */
#line 852 "dlansf.f"
		    i__1 = k - 1;
#line 852 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 853 "dlansf.f"
			dlassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 855 "dlansf.f"
		    }
#line 856 "dlansf.f"
		    i__1 = *n - 1;
#line 856 "dlansf.f"
		    for (j = k; j <= i__1; ++j) {
#line 857 "dlansf.f"
			dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,k) */
#line 859 "dlansf.f"
		    }
#line 860 "dlansf.f"
		    i__1 = k - 3;
#line 860 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 861 "dlansf.f"
			i__2 = k - j - 2;
#line 861 "dlansf.f"
			dlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(1,0) */
#line 863 "dlansf.f"
		    }
#line 864 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 866 "dlansf.f"
		    i__1 = lda + 1;
#line 866 "dlansf.f"
		    dlassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 868 "dlansf.f"
		    i__1 = k - 1;
#line 868 "dlansf.f"
		    i__2 = lda + 1;
#line 868 "dlansf.f"
		    dlassq_(&i__1, &a[1], &i__2, &scale, &s);
/*                 tri L at A(1,0) */
#line 870 "dlansf.f"
		}
#line 871 "dlansf.f"
	    }
#line 872 "dlansf.f"
	} else {
/*           n is even */
#line 874 "dlansf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 876 "dlansf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 878 "dlansf.f"
		    i__1 = k - 2;
#line 878 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 879 "dlansf.f"
			i__2 = k - j - 1;
#line 879 "dlansf.f"
			dlassq_(&i__2, &a[k + j + 2 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k+1,0) */
#line 881 "dlansf.f"
		    }
#line 882 "dlansf.f"
		    i__1 = k - 1;
#line 882 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 883 "dlansf.f"
			i__2 = k + j;
#line 883 "dlansf.f"
			dlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 885 "dlansf.f"
		    }
#line 886 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 888 "dlansf.f"
		    i__1 = lda + 1;
#line 888 "dlansf.f"
		    dlassq_(&k, &a[k + 1], &i__1, &scale, &s);
/*                 tri L at A(k+1,0) */
#line 890 "dlansf.f"
		    i__1 = lda + 1;
#line 890 "dlansf.f"
		    dlassq_(&k, &a[k], &i__1, &scale, &s);
/*                 tri U at A(k,0) */
#line 892 "dlansf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 894 "dlansf.f"
		    i__1 = k - 1;
#line 894 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 895 "dlansf.f"
			i__2 = *n - j - 1;
#line 895 "dlansf.f"
			dlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(1,0) */
#line 897 "dlansf.f"
		    }
#line 898 "dlansf.f"
		    i__1 = k - 1;
#line 898 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 899 "dlansf.f"
			dlassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 901 "dlansf.f"
		    }
#line 902 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 904 "dlansf.f"
		    i__1 = lda + 1;
#line 904 "dlansf.f"
		    dlassq_(&k, &a[1], &i__1, &scale, &s);
/*                 tri L at A(1,0) */
#line 906 "dlansf.f"
		    i__1 = lda + 1;
#line 906 "dlansf.f"
		    dlassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 908 "dlansf.f"
		}
#line 909 "dlansf.f"
	    } else {
/*              A is xpose */
#line 911 "dlansf.f"
		if (ilu == 0) {
/*                 A**T is upper */
#line 913 "dlansf.f"
		    i__1 = k - 1;
#line 913 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 914 "dlansf.f"
			dlassq_(&j, &a[(k + 1 + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k+1) */
#line 916 "dlansf.f"
		    }
#line 917 "dlansf.f"
		    i__1 = k - 1;
#line 917 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 918 "dlansf.f"
			dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k rect. at A(0,0) */
#line 920 "dlansf.f"
		    }
#line 921 "dlansf.f"
		    i__1 = k - 2;
#line 921 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 922 "dlansf.f"
			i__2 = k - j - 1;
#line 922 "dlansf.f"
			dlassq_(&i__2, &a[j + 1 + (j + k) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k) */
#line 925 "dlansf.f"
		    }
#line 926 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 928 "dlansf.f"
		    i__1 = lda + 1;
#line 928 "dlansf.f"
		    dlassq_(&k, &a[(k + 1) * lda], &i__1, &scale, &s);
/*                 tri U at A(0,k+1) */
#line 930 "dlansf.f"
		    i__1 = lda + 1;
#line 930 "dlansf.f"
		    dlassq_(&k, &a[k * lda], &i__1, &scale, &s);
/*                 tri L at A(0,k) */
#line 932 "dlansf.f"
		} else {
/*                 A**T is lower */
#line 934 "dlansf.f"
		    i__1 = k - 1;
#line 934 "dlansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 935 "dlansf.f"
			dlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 937 "dlansf.f"
		    }
#line 938 "dlansf.f"
		    i__1 = *n;
#line 938 "dlansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 939 "dlansf.f"
			dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k rect. at A(0,k+1) */
#line 941 "dlansf.f"
		    }
#line 942 "dlansf.f"
		    i__1 = k - 2;
#line 942 "dlansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 943 "dlansf.f"
			i__2 = k - j - 1;
#line 943 "dlansf.f"
			dlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(0,0) */
#line 945 "dlansf.f"
		    }
#line 946 "dlansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 948 "dlansf.f"
		    i__1 = lda + 1;
#line 948 "dlansf.f"
		    dlassq_(&k, &a[lda], &i__1, &scale, &s);
/*                 tri L at A(0,1) */
#line 950 "dlansf.f"
		    i__1 = lda + 1;
#line 950 "dlansf.f"
		    dlassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 952 "dlansf.f"
		}
#line 953 "dlansf.f"
	    }
#line 954 "dlansf.f"
	}
#line 955 "dlansf.f"
	value = scale * sqrt(s);
#line 956 "dlansf.f"
    }

#line 958 "dlansf.f"
    ret_val = value;
#line 959 "dlansf.f"
    return ret_val;

/*     End of DLANSF */

} /* dlansf_ */

