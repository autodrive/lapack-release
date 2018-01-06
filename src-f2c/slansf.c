#line 1 "slansf.f"
/* slansf.f -- translated by f2c (version 20100827).
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

#line 1 "slansf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANSF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANSF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slansf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slansf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slansf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, TRANSR, UPLO */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( 0: * ), WORK( 0: * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANSF returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > real symmetric matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \return SLANSF */
/* > \verbatim */
/* > */
/* >    SLANSF = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANSF as described */
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
/* >          The order of the matrix A. N >= 0. When N = 0, SLANSF is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension ( N*(N+1)/2 ); */
/* >          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') */
/* >          part of the symmetric matrix A stored in RFP format. See the */
/* >          "Notes" below for more details. */
/* >          Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
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

/*  ===================================================================== */
doublereal slansf_(char *norm, char *transr, char *uplo, integer *n, 
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
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

/*     .. */
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

#line 248 "slansf.f"
    if (*n == 0) {
#line 249 "slansf.f"
	ret_val = 0.;
#line 250 "slansf.f"
	return ret_val;
#line 251 "slansf.f"
    } else if (*n == 1) {
#line 252 "slansf.f"
	ret_val = abs(a[0]);
#line 253 "slansf.f"
	return ret_val;
#line 254 "slansf.f"
    }

/*     set noe = 1 if n is odd. if n is even set noe=0 */

#line 258 "slansf.f"
    noe = 1;
#line 259 "slansf.f"
    if (*n % 2 == 0) {
#line 259 "slansf.f"
	noe = 0;
#line 259 "slansf.f"
    }

/*     set ifm = 0 when form='T or 't' and 1 otherwise */

#line 264 "slansf.f"
    ifm = 1;
#line 265 "slansf.f"
    if (lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 265 "slansf.f"
	ifm = 0;
#line 265 "slansf.f"
    }

/*     set ilu = 0 when uplo='U or 'u' and 1 otherwise */

#line 270 "slansf.f"
    ilu = 1;
#line 271 "slansf.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 271 "slansf.f"
	ilu = 0;
#line 271 "slansf.f"
    }

/*     set lda = (n+1)/2 when ifm = 0 */
/*     set lda = n when ifm = 1 and noe = 1 */
/*     set lda = n+1 when ifm = 1 and noe = 0 */

#line 278 "slansf.f"
    if (ifm == 1) {
#line 279 "slansf.f"
	if (noe == 1) {
#line 280 "slansf.f"
	    lda = *n;
#line 281 "slansf.f"
	} else {
/*           noe=0 */
#line 283 "slansf.f"
	    lda = *n + 1;
#line 284 "slansf.f"
	}
#line 285 "slansf.f"
    } else {
/*        ifm=0 */
#line 287 "slansf.f"
	lda = (*n + 1) / 2;
#line 288 "slansf.f"
    }

#line 290 "slansf.f"
    if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*       Find max(abs(A(i,j))). */

#line 294 "slansf.f"
	k = (*n + 1) / 2;
#line 295 "slansf.f"
	value = 0.;
#line 296 "slansf.f"
	if (noe == 1) {
/*           n is odd */
#line 298 "slansf.f"
	    if (ifm == 1) {
/*           A is n by k */
#line 300 "slansf.f"
		i__1 = k - 1;
#line 300 "slansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 301 "slansf.f"
		    i__2 = *n - 1;
#line 301 "slansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 302 "slansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 303 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 303 "slansf.f"
			    value = temp;
#line 303 "slansf.f"
			}
#line 305 "slansf.f"
		    }
#line 306 "slansf.f"
		}
#line 307 "slansf.f"
	    } else {
/*              xpose case; A is k by n */
#line 309 "slansf.f"
		i__1 = *n - 1;
#line 309 "slansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 310 "slansf.f"
		    i__2 = k - 1;
#line 310 "slansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 311 "slansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 312 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 312 "slansf.f"
			    value = temp;
#line 312 "slansf.f"
			}
#line 314 "slansf.f"
		    }
#line 315 "slansf.f"
		}
#line 316 "slansf.f"
	    }
#line 317 "slansf.f"
	} else {
/*           n is even */
#line 319 "slansf.f"
	    if (ifm == 1) {
/*              A is n+1 by k */
#line 321 "slansf.f"
		i__1 = k - 1;
#line 321 "slansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 322 "slansf.f"
		    i__2 = *n;
#line 322 "slansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 323 "slansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 324 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 324 "slansf.f"
			    value = temp;
#line 324 "slansf.f"
			}
#line 326 "slansf.f"
		    }
#line 327 "slansf.f"
		}
#line 328 "slansf.f"
	    } else {
/*              xpose case; A is k by n+1 */
#line 330 "slansf.f"
		i__1 = *n;
#line 330 "slansf.f"
		for (j = 0; j <= i__1; ++j) {
#line 331 "slansf.f"
		    i__2 = k - 1;
#line 331 "slansf.f"
		    for (i__ = 0; i__ <= i__2; ++i__) {
#line 332 "slansf.f"
			temp = (d__1 = a[i__ + j * lda], abs(d__1));
#line 333 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 333 "slansf.f"
			    value = temp;
#line 333 "slansf.f"
			}
#line 335 "slansf.f"
		    }
#line 336 "slansf.f"
		}
#line 337 "slansf.f"
	    }
#line 338 "slansf.f"
	}
#line 339 "slansf.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1) || lsame_(norm, "O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

#line 344 "slansf.f"
	if (ifm == 1) {
#line 345 "slansf.f"
	    k = *n / 2;
#line 346 "slansf.f"
	    if (noe == 1) {
/*              n is odd */
#line 348 "slansf.f"
		if (ilu == 0) {
#line 349 "slansf.f"
		    i__1 = k - 1;
#line 349 "slansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 350 "slansf.f"
			work[i__] = 0.;
#line 351 "slansf.f"
		    }
#line 352 "slansf.f"
		    i__1 = k;
#line 352 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 353 "slansf.f"
			s = 0.;
#line 354 "slansf.f"
			i__2 = k + j - 1;
#line 354 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 355 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(i,j+k) */
#line 357 "slansf.f"
			    s += aa;
#line 358 "slansf.f"
			    work[i__] += aa;
#line 359 "slansf.f"
			}
#line 360 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 362 "slansf.f"
			work[j + k] = s + aa;
#line 363 "slansf.f"
			if (i__ == k + k) {
#line 363 "slansf.f"
			    goto L10;
#line 363 "slansf.f"
			}
#line 365 "slansf.f"
			++i__;
#line 366 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 368 "slansf.f"
			work[j] += aa;
#line 369 "slansf.f"
			s = 0.;
#line 370 "slansf.f"
			i__2 = k - 1;
#line 370 "slansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 371 "slansf.f"
			    ++i__;
#line 372 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 374 "slansf.f"
			    s += aa;
#line 375 "slansf.f"
			    work[l] += aa;
#line 376 "slansf.f"
			}
#line 377 "slansf.f"
			work[j] += s;
#line 378 "slansf.f"
		    }
#line 379 "slansf.f"
L10:
#line 380 "slansf.f"
		    value = work[0];
#line 381 "slansf.f"
		    i__1 = *n - 1;
#line 381 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "slansf.f"
			temp = work[i__];
#line 383 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 383 "slansf.f"
			    value = temp;
#line 383 "slansf.f"
			}
#line 385 "slansf.f"
		    }
#line 386 "slansf.f"
		} else {
/*                 ilu = 1 */
#line 388 "slansf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 390 "slansf.f"
		    i__1 = *n - 1;
#line 390 "slansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 391 "slansf.f"
			work[i__] = 0.;
#line 392 "slansf.f"
		    }
#line 393 "slansf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 394 "slansf.f"
			s = 0.;
#line 395 "slansf.f"
			i__1 = j - 2;
#line 395 "slansf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 396 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,i+k) */
#line 398 "slansf.f"
			    s += aa;
#line 399 "slansf.f"
			    work[i__ + k] += aa;
#line 400 "slansf.f"
			}
#line 401 "slansf.f"
			if (j > 0) {
#line 402 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,j+k) */
#line 404 "slansf.f"
			    s += aa;
#line 405 "slansf.f"
			    work[i__ + k] += s;
/*                       i=j */
#line 407 "slansf.f"
			    ++i__;
#line 408 "slansf.f"
			}
#line 409 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 411 "slansf.f"
			work[j] = aa;
#line 412 "slansf.f"
			s = 0.;
#line 413 "slansf.f"
			i__1 = *n - 1;
#line 413 "slansf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 414 "slansf.f"
			    ++i__;
#line 415 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 417 "slansf.f"
			    s += aa;
#line 418 "slansf.f"
			    work[l] += aa;
#line 419 "slansf.f"
			}
#line 420 "slansf.f"
			work[j] += s;
#line 421 "slansf.f"
		    }
#line 422 "slansf.f"
		    value = work[0];
#line 423 "slansf.f"
		    i__1 = *n - 1;
#line 423 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 424 "slansf.f"
			temp = work[i__];
#line 425 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 425 "slansf.f"
			    value = temp;
#line 425 "slansf.f"
			}
#line 427 "slansf.f"
		    }
#line 428 "slansf.f"
		}
#line 429 "slansf.f"
	    } else {
/*              n is even */
#line 431 "slansf.f"
		if (ilu == 0) {
#line 432 "slansf.f"
		    i__1 = k - 1;
#line 432 "slansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 433 "slansf.f"
			work[i__] = 0.;
#line 434 "slansf.f"
		    }
#line 435 "slansf.f"
		    i__1 = k - 1;
#line 435 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 436 "slansf.f"
			s = 0.;
#line 437 "slansf.f"
			i__2 = k + j - 1;
#line 437 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 438 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(i,j+k) */
#line 440 "slansf.f"
			    s += aa;
#line 441 "slansf.f"
			    work[i__] += aa;
#line 442 "slansf.f"
			}
#line 443 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 445 "slansf.f"
			work[j + k] = s + aa;
#line 446 "slansf.f"
			++i__;
#line 447 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 449 "slansf.f"
			work[j] += aa;
#line 450 "slansf.f"
			s = 0.;
#line 451 "slansf.f"
			i__2 = k - 1;
#line 451 "slansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 452 "slansf.f"
			    ++i__;
#line 453 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 455 "slansf.f"
			    s += aa;
#line 456 "slansf.f"
			    work[l] += aa;
#line 457 "slansf.f"
			}
#line 458 "slansf.f"
			work[j] += s;
#line 459 "slansf.f"
		    }
#line 460 "slansf.f"
		    value = work[0];
#line 461 "slansf.f"
		    i__1 = *n - 1;
#line 461 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 462 "slansf.f"
			temp = work[i__];
#line 463 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 463 "slansf.f"
			    value = temp;
#line 463 "slansf.f"
			}
#line 465 "slansf.f"
		    }
#line 466 "slansf.f"
		} else {
/*                 ilu = 1 */
#line 468 "slansf.f"
		    i__1 = *n - 1;
#line 468 "slansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 469 "slansf.f"
			work[i__] = 0.;
#line 470 "slansf.f"
		    }
#line 471 "slansf.f"
		    for (j = k - 1; j >= 0; --j) {
#line 472 "slansf.f"
			s = 0.;
#line 473 "slansf.f"
			i__1 = j - 1;
#line 473 "slansf.f"
			for (i__ = 0; i__ <= i__1; ++i__) {
#line 474 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(j+k,i+k) */
#line 476 "slansf.f"
			    s += aa;
#line 477 "slansf.f"
			    work[i__ + k] += aa;
#line 478 "slansf.f"
			}
#line 479 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j+k,j+k) */
#line 481 "slansf.f"
			s += aa;
#line 482 "slansf.f"
			work[i__ + k] += s;
/*                    i=j */
#line 484 "slansf.f"
			++i__;
#line 485 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    -> A(j,j) */
#line 487 "slansf.f"
			work[j] = aa;
#line 488 "slansf.f"
			s = 0.;
#line 489 "slansf.f"
			i__1 = *n - 1;
#line 489 "slansf.f"
			for (l = j + 1; l <= i__1; ++l) {
#line 490 "slansf.f"
			    ++i__;
#line 491 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       -> A(l,j) */
#line 493 "slansf.f"
			    s += aa;
#line 494 "slansf.f"
			    work[l] += aa;
#line 495 "slansf.f"
			}
#line 496 "slansf.f"
			work[j] += s;
#line 497 "slansf.f"
		    }
#line 498 "slansf.f"
		    value = work[0];
#line 499 "slansf.f"
		    i__1 = *n - 1;
#line 499 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 500 "slansf.f"
			temp = work[i__];
#line 501 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 501 "slansf.f"
			    value = temp;
#line 501 "slansf.f"
			}
#line 503 "slansf.f"
		    }
#line 504 "slansf.f"
		}
#line 505 "slansf.f"
	    }
#line 506 "slansf.f"
	} else {
/*           ifm=0 */
#line 508 "slansf.f"
	    k = *n / 2;
#line 509 "slansf.f"
	    if (noe == 1) {
/*              n is odd */
#line 511 "slansf.f"
		if (ilu == 0) {
#line 512 "slansf.f"
		    n1 = k;
/*                 n/2 */
#line 514 "slansf.f"
		    ++k;
/*                 k is the row size and lda */
#line 516 "slansf.f"
		    i__1 = *n - 1;
#line 516 "slansf.f"
		    for (i__ = n1; i__ <= i__1; ++i__) {
#line 517 "slansf.f"
			work[i__] = 0.;
#line 518 "slansf.f"
		    }
#line 519 "slansf.f"
		    i__1 = n1 - 1;
#line 519 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 520 "slansf.f"
			s = 0.;
#line 521 "slansf.f"
			i__2 = k - 1;
#line 521 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 522 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,n1+i) */
#line 524 "slansf.f"
			    work[i__ + n1] += aa;
#line 525 "slansf.f"
			    s += aa;
#line 526 "slansf.f"
			}
#line 527 "slansf.f"
			work[j] = s;
#line 528 "slansf.f"
		    }
/*                 j=n1=k-1 is special */
#line 530 "slansf.f"
		    s = (d__1 = a[j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 532 "slansf.f"
		    i__1 = k - 1;
#line 532 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 533 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k-1,i+n1) */
#line 535 "slansf.f"
			work[i__ + n1] += aa;
#line 536 "slansf.f"
			s += aa;
#line 537 "slansf.f"
		    }
#line 538 "slansf.f"
		    work[j] += s;
#line 539 "slansf.f"
		    i__1 = *n - 1;
#line 539 "slansf.f"
		    for (j = k; j <= i__1; ++j) {
#line 540 "slansf.f"
			s = 0.;
#line 541 "slansf.f"
			i__2 = j - k - 1;
#line 541 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 542 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(i,j-k) */
#line 544 "slansf.f"
			    work[i__] += aa;
#line 545 "slansf.f"
			    s += aa;
#line 546 "slansf.f"
			}
/*                    i=j-k */
#line 548 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j-k,j-k) */
#line 550 "slansf.f"
			s += aa;
#line 551 "slansf.f"
			work[j - k] += s;
#line 552 "slansf.f"
			++i__;
#line 553 "slansf.f"
			s = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j,j) */
#line 555 "slansf.f"
			i__2 = *n - 1;
#line 555 "slansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 556 "slansf.f"
			    ++i__;
#line 557 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,l) */
#line 559 "slansf.f"
			    work[l] += aa;
#line 560 "slansf.f"
			    s += aa;
#line 561 "slansf.f"
			}
#line 562 "slansf.f"
			work[j] += s;
#line 563 "slansf.f"
		    }
#line 564 "slansf.f"
		    value = work[0];
#line 565 "slansf.f"
		    i__1 = *n - 1;
#line 565 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 566 "slansf.f"
			temp = work[i__];
#line 567 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 567 "slansf.f"
			    value = temp;
#line 567 "slansf.f"
			}
#line 569 "slansf.f"
		    }
#line 570 "slansf.f"
		} else {
/*                 ilu=1 */
#line 572 "slansf.f"
		    ++k;
/*                 k=(n+1)/2 for n odd and ilu=1 */
#line 574 "slansf.f"
		    i__1 = *n - 1;
#line 574 "slansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 575 "slansf.f"
			work[i__] = 0.;
#line 576 "slansf.f"
		    }
#line 577 "slansf.f"
		    i__1 = k - 2;
#line 577 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
/*                    process */
#line 579 "slansf.f"
			s = 0.;
#line 580 "slansf.f"
			i__2 = j - 1;
#line 580 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 581 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i) */
#line 583 "slansf.f"
			    work[i__] += aa;
#line 584 "slansf.f"
			    s += aa;
#line 585 "slansf.f"
			}
#line 586 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    i=j so process of A(j,j) */
#line 588 "slansf.f"
			s += aa;
#line 589 "slansf.f"
			work[j] = s;
/*                    is initialised here */
#line 591 "slansf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 593 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
#line 594 "slansf.f"
			s = aa;
#line 595 "slansf.f"
			i__2 = *n - 1;
#line 595 "slansf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 596 "slansf.f"
			    ++i__;
#line 597 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(l,k+j) */
#line 599 "slansf.f"
			    s += aa;
#line 600 "slansf.f"
			    work[l] += aa;
#line 601 "slansf.f"
			}
#line 602 "slansf.f"
			work[k + j] += s;
#line 603 "slansf.f"
		    }
/*                 j=k-1 is special :process col A(k-1,0:k-1) */
#line 605 "slansf.f"
		    s = 0.;
#line 606 "slansf.f"
		    i__1 = k - 2;
#line 606 "slansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 607 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,i) */
#line 609 "slansf.f"
			work[i__] += aa;
#line 610 "slansf.f"
			s += aa;
#line 611 "slansf.f"
		    }
/*                 i=k-1 */
#line 613 "slansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 615 "slansf.f"
		    s += aa;
#line 616 "slansf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 618 "slansf.f"
		    i__1 = *n - 1;
#line 618 "slansf.f"
		    for (j = k; j <= i__1; ++j) {
/*                    process col j of A = A(j,0:k-1) */
#line 620 "slansf.f"
			s = 0.;
#line 621 "slansf.f"
			i__2 = k - 1;
#line 621 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 622 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i) */
#line 624 "slansf.f"
			    work[i__] += aa;
#line 625 "slansf.f"
			    s += aa;
#line 626 "slansf.f"
			}
#line 627 "slansf.f"
			work[j] += s;
#line 628 "slansf.f"
		    }
#line 629 "slansf.f"
		    value = work[0];
#line 630 "slansf.f"
		    i__1 = *n - 1;
#line 630 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 631 "slansf.f"
			temp = work[i__];
#line 632 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 632 "slansf.f"
			    value = temp;
#line 632 "slansf.f"
			}
#line 634 "slansf.f"
		    }
#line 635 "slansf.f"
		}
#line 636 "slansf.f"
	    } else {
/*              n is even */
#line 638 "slansf.f"
		if (ilu == 0) {
#line 639 "slansf.f"
		    i__1 = *n - 1;
#line 639 "slansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 640 "slansf.f"
			work[i__] = 0.;
#line 641 "slansf.f"
		    }
#line 642 "slansf.f"
		    i__1 = k - 1;
#line 642 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 643 "slansf.f"
			s = 0.;
#line 644 "slansf.f"
			i__2 = k - 1;
#line 644 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 645 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,i+k) */
#line 647 "slansf.f"
			    work[i__ + k] += aa;
#line 648 "slansf.f"
			    s += aa;
#line 649 "slansf.f"
			}
#line 650 "slansf.f"
			work[j] = s;
#line 651 "slansf.f"
		    }
/*                 j=k */
#line 653 "slansf.f"
		    aa = (d__1 = a[j * lda], abs(d__1));
/*                 A(k,k) */
#line 655 "slansf.f"
		    s = aa;
#line 656 "slansf.f"
		    i__1 = k - 1;
#line 656 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 657 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,k+i) */
#line 659 "slansf.f"
			work[i__ + k] += aa;
#line 660 "slansf.f"
			s += aa;
#line 661 "slansf.f"
		    }
#line 662 "slansf.f"
		    work[j] += s;
#line 663 "slansf.f"
		    i__1 = *n - 1;
#line 663 "slansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 664 "slansf.f"
			s = 0.;
#line 665 "slansf.f"
			i__2 = j - 2 - k;
#line 665 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 666 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(i,j-k-1) */
#line 668 "slansf.f"
			    work[i__] += aa;
#line 669 "slansf.f"
			    s += aa;
#line 670 "slansf.f"
			}
/*                     i=j-1-k */
#line 672 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j-k-1,j-k-1) */
#line 674 "slansf.f"
			s += aa;
#line 675 "slansf.f"
			work[j - k - 1] += s;
#line 676 "slansf.f"
			++i__;
#line 677 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(j,j) */
#line 679 "slansf.f"
			s = aa;
#line 680 "slansf.f"
			i__2 = *n - 1;
#line 680 "slansf.f"
			for (l = j + 1; l <= i__2; ++l) {
#line 681 "slansf.f"
			    ++i__;
#line 682 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j,l) */
#line 684 "slansf.f"
			    work[l] += aa;
#line 685 "slansf.f"
			    s += aa;
#line 686 "slansf.f"
			}
#line 687 "slansf.f"
			work[j] += s;
#line 688 "slansf.f"
		    }
/*                 j=n */
#line 690 "slansf.f"
		    s = 0.;
#line 691 "slansf.f"
		    i__1 = k - 2;
#line 691 "slansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 692 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(i,k-1) */
#line 694 "slansf.f"
			work[i__] += aa;
#line 695 "slansf.f"
			s += aa;
#line 696 "slansf.f"
		    }
/*                 i=k-1 */
#line 698 "slansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 700 "slansf.f"
		    s += aa;
#line 701 "slansf.f"
		    work[i__] += s;
#line 702 "slansf.f"
		    value = work[0];
#line 703 "slansf.f"
		    i__1 = *n - 1;
#line 703 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 704 "slansf.f"
			temp = work[i__];
#line 705 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 705 "slansf.f"
			    value = temp;
#line 705 "slansf.f"
			}
#line 707 "slansf.f"
		    }
#line 708 "slansf.f"
		} else {
/*                 ilu=1 */
#line 710 "slansf.f"
		    i__1 = *n - 1;
#line 710 "slansf.f"
		    for (i__ = k; i__ <= i__1; ++i__) {
#line 711 "slansf.f"
			work[i__] = 0.;
#line 712 "slansf.f"
		    }
/*                 j=0 is special :process col A(k:n-1,k) */
#line 714 "slansf.f"
		    s = abs(a[0]);
/*                 A(k,k) */
#line 716 "slansf.f"
		    i__1 = k - 1;
#line 716 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 717 "slansf.f"
			aa = (d__1 = a[i__], abs(d__1));
/*                    A(k+i,k) */
#line 719 "slansf.f"
			work[i__ + k] += aa;
#line 720 "slansf.f"
			s += aa;
#line 721 "slansf.f"
		    }
#line 722 "slansf.f"
		    work[k] += s;
#line 723 "slansf.f"
		    i__1 = k - 1;
#line 723 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
/*                    process */
#line 725 "slansf.f"
			s = 0.;
#line 726 "slansf.f"
			i__2 = j - 2;
#line 726 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 727 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j-1,i) */
#line 729 "slansf.f"
			    work[i__] += aa;
#line 730 "slansf.f"
			    s += aa;
#line 731 "slansf.f"
			}
#line 732 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    i=j-1 so process of A(j-1,j-1) */
#line 734 "slansf.f"
			s += aa;
#line 735 "slansf.f"
			work[j - 1] = s;
/*                    is initialised here */
#line 737 "slansf.f"
			++i__;
/*                    i=j process A(j+k,j+k) */
#line 739 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
#line 740 "slansf.f"
			s = aa;
#line 741 "slansf.f"
			i__2 = *n - 1;
#line 741 "slansf.f"
			for (l = k + j + 1; l <= i__2; ++l) {
#line 742 "slansf.f"
			    ++i__;
#line 743 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(l,k+j) */
#line 745 "slansf.f"
			    s += aa;
#line 746 "slansf.f"
			    work[l] += aa;
#line 747 "slansf.f"
			}
#line 748 "slansf.f"
			work[k + j] += s;
#line 749 "slansf.f"
		    }
/*                 j=k is special :process col A(k,0:k-1) */
#line 751 "slansf.f"
		    s = 0.;
#line 752 "slansf.f"
		    i__1 = k - 2;
#line 752 "slansf.f"
		    for (i__ = 0; i__ <= i__1; ++i__) {
#line 753 "slansf.f"
			aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                    A(k,i) */
#line 755 "slansf.f"
			work[i__] += aa;
#line 756 "slansf.f"
			s += aa;
#line 757 "slansf.f"
		    }
/*                 i=k-1 */
#line 759 "slansf.f"
		    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                 A(k-1,k-1) */
#line 761 "slansf.f"
		    s += aa;
#line 762 "slansf.f"
		    work[i__] = s;
/*                 done with col j=k+1 */
#line 764 "slansf.f"
		    i__1 = *n;
#line 764 "slansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
/*                    process col j-1 of A = A(j-1,0:k-1) */
#line 766 "slansf.f"
			s = 0.;
#line 767 "slansf.f"
			i__2 = k - 1;
#line 767 "slansf.f"
			for (i__ = 0; i__ <= i__2; ++i__) {
#line 768 "slansf.f"
			    aa = (d__1 = a[i__ + j * lda], abs(d__1));
/*                       A(j-1,i) */
#line 770 "slansf.f"
			    work[i__] += aa;
#line 771 "slansf.f"
			    s += aa;
#line 772 "slansf.f"
			}
#line 773 "slansf.f"
			work[j - 1] += s;
#line 774 "slansf.f"
		    }
#line 775 "slansf.f"
		    value = work[0];
#line 776 "slansf.f"
		    i__1 = *n - 1;
#line 776 "slansf.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 777 "slansf.f"
			temp = work[i__];
#line 778 "slansf.f"
			if (value < temp || sisnan_(&temp)) {
#line 778 "slansf.f"
			    value = temp;
#line 778 "slansf.f"
			}
#line 780 "slansf.f"
		    }
#line 781 "slansf.f"
		}
#line 782 "slansf.f"
	    }
#line 783 "slansf.f"
	}
#line 784 "slansf.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*       Find normF(A). */

#line 788 "slansf.f"
	k = (*n + 1) / 2;
#line 789 "slansf.f"
	scale = 0.;
#line 790 "slansf.f"
	s = 1.;
#line 791 "slansf.f"
	if (noe == 1) {
/*           n is odd */
#line 793 "slansf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 795 "slansf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 797 "slansf.f"
		    i__1 = k - 3;
#line 797 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 798 "slansf.f"
			i__2 = k - j - 2;
#line 798 "slansf.f"
			slassq_(&i__2, &a[k + j + 1 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k,0) */
#line 800 "slansf.f"
		    }
#line 801 "slansf.f"
		    i__1 = k - 1;
#line 801 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 802 "slansf.f"
			i__2 = k + j - 1;
#line 802 "slansf.f"
			slassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 804 "slansf.f"
		    }
#line 805 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 807 "slansf.f"
		    i__1 = k - 1;
#line 807 "slansf.f"
		    i__2 = lda + 1;
#line 807 "slansf.f"
		    slassq_(&i__1, &a[k], &i__2, &scale, &s);
/*                 tri L at A(k,0) */
#line 809 "slansf.f"
		    i__1 = lda + 1;
#line 809 "slansf.f"
		    slassq_(&k, &a[k - 1], &i__1, &scale, &s);
/*                 tri U at A(k-1,0) */
#line 811 "slansf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 813 "slansf.f"
		    i__1 = k - 1;
#line 813 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 814 "slansf.f"
			i__2 = *n - j - 1;
#line 814 "slansf.f"
			slassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(0,0) */
#line 816 "slansf.f"
		    }
#line 817 "slansf.f"
		    i__1 = k - 2;
#line 817 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 818 "slansf.f"
			slassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 820 "slansf.f"
		    }
#line 821 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 823 "slansf.f"
		    i__1 = lda + 1;
#line 823 "slansf.f"
		    slassq_(&k, a, &i__1, &scale, &s);
/*                 tri L at A(0,0) */
#line 825 "slansf.f"
		    i__1 = k - 1;
#line 825 "slansf.f"
		    i__2 = lda + 1;
#line 825 "slansf.f"
		    slassq_(&i__1, &a[lda], &i__2, &scale, &s);
/*                 tri U at A(0,1) */
#line 827 "slansf.f"
		}
#line 828 "slansf.f"
	    } else {
/*              A is xpose */
#line 830 "slansf.f"
		if (ilu == 0) {
/*                 A**T is upper */
#line 832 "slansf.f"
		    i__1 = k - 2;
#line 832 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 833 "slansf.f"
			slassq_(&j, &a[(k + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k) */
#line 835 "slansf.f"
		    }
#line 836 "slansf.f"
		    i__1 = k - 2;
#line 836 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 837 "slansf.f"
			slassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,0) */
#line 839 "slansf.f"
		    }
#line 840 "slansf.f"
		    i__1 = k - 2;
#line 840 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 841 "slansf.f"
			i__2 = k - j - 1;
#line 841 "slansf.f"
			slassq_(&i__2, &a[j + 1 + (j + k - 1) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k-1) */
#line 844 "slansf.f"
		    }
#line 845 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 847 "slansf.f"
		    i__1 = k - 1;
#line 847 "slansf.f"
		    i__2 = lda + 1;
#line 847 "slansf.f"
		    slassq_(&i__1, &a[k * lda], &i__2, &scale, &s);
/*                 tri U at A(0,k) */
#line 849 "slansf.f"
		    i__1 = lda + 1;
#line 849 "slansf.f"
		    slassq_(&k, &a[(k - 1) * lda], &i__1, &scale, &s);
/*                 tri L at A(0,k-1) */
#line 851 "slansf.f"
		} else {
/*                 A**T is lower */
#line 853 "slansf.f"
		    i__1 = k - 1;
#line 853 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 854 "slansf.f"
			slassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 856 "slansf.f"
		    }
#line 857 "slansf.f"
		    i__1 = *n - 1;
#line 857 "slansf.f"
		    for (j = k; j <= i__1; ++j) {
#line 858 "slansf.f"
			slassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k-1 rect. at A(0,k) */
#line 860 "slansf.f"
		    }
#line 861 "slansf.f"
		    i__1 = k - 3;
#line 861 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 862 "slansf.f"
			i__2 = k - j - 2;
#line 862 "slansf.f"
			slassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(1,0) */
#line 864 "slansf.f"
		    }
#line 865 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 867 "slansf.f"
		    i__1 = lda + 1;
#line 867 "slansf.f"
		    slassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 869 "slansf.f"
		    i__1 = k - 1;
#line 869 "slansf.f"
		    i__2 = lda + 1;
#line 869 "slansf.f"
		    slassq_(&i__1, &a[1], &i__2, &scale, &s);
/*                 tri L at A(1,0) */
#line 871 "slansf.f"
		}
#line 872 "slansf.f"
	    }
#line 873 "slansf.f"
	} else {
/*           n is even */
#line 875 "slansf.f"
	    if (ifm == 1) {
/*              A is normal */
#line 877 "slansf.f"
		if (ilu == 0) {
/*                 A is upper */
#line 879 "slansf.f"
		    i__1 = k - 2;
#line 879 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 880 "slansf.f"
			i__2 = k - j - 1;
#line 880 "slansf.f"
			slassq_(&i__2, &a[k + j + 2 + j * lda], &c__1, &scale,
				 &s);
/*                    L at A(k+1,0) */
#line 882 "slansf.f"
		    }
#line 883 "slansf.f"
		    i__1 = k - 1;
#line 883 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 884 "slansf.f"
			i__2 = k + j;
#line 884 "slansf.f"
			slassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
/*                    trap U at A(0,0) */
#line 886 "slansf.f"
		    }
#line 887 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 889 "slansf.f"
		    i__1 = lda + 1;
#line 889 "slansf.f"
		    slassq_(&k, &a[k + 1], &i__1, &scale, &s);
/*                 tri L at A(k+1,0) */
#line 891 "slansf.f"
		    i__1 = lda + 1;
#line 891 "slansf.f"
		    slassq_(&k, &a[k], &i__1, &scale, &s);
/*                 tri U at A(k,0) */
#line 893 "slansf.f"
		} else {
/*                 ilu=1 & A is lower */
#line 895 "slansf.f"
		    i__1 = k - 1;
#line 895 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 896 "slansf.f"
			i__2 = *n - j - 1;
#line 896 "slansf.f"
			slassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s)
				;
/*                    trap L at A(1,0) */
#line 898 "slansf.f"
		    }
#line 899 "slansf.f"
		    i__1 = k - 1;
#line 899 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 900 "slansf.f"
			slassq_(&j, &a[j * lda], &c__1, &scale, &s);
/*                    U at A(0,0) */
#line 902 "slansf.f"
		    }
#line 903 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 905 "slansf.f"
		    i__1 = lda + 1;
#line 905 "slansf.f"
		    slassq_(&k, &a[1], &i__1, &scale, &s);
/*                 tri L at A(1,0) */
#line 907 "slansf.f"
		    i__1 = lda + 1;
#line 907 "slansf.f"
		    slassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 909 "slansf.f"
		}
#line 910 "slansf.f"
	    } else {
/*              A is xpose */
#line 912 "slansf.f"
		if (ilu == 0) {
/*                 A**T is upper */
#line 914 "slansf.f"
		    i__1 = k - 1;
#line 914 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 915 "slansf.f"
			slassq_(&j, &a[(k + 1 + j) * lda], &c__1, &scale, &s);
/*                    U at A(0,k+1) */
#line 917 "slansf.f"
		    }
#line 918 "slansf.f"
		    i__1 = k - 1;
#line 918 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 919 "slansf.f"
			slassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k rect. at A(0,0) */
#line 921 "slansf.f"
		    }
#line 922 "slansf.f"
		    i__1 = k - 2;
#line 922 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 923 "slansf.f"
			i__2 = k - j - 1;
#line 923 "slansf.f"
			slassq_(&i__2, &a[j + 1 + (j + k) * lda], &c__1, &
				scale, &s);
/*                    L at A(0,k) */
#line 926 "slansf.f"
		    }
#line 927 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 929 "slansf.f"
		    i__1 = lda + 1;
#line 929 "slansf.f"
		    slassq_(&k, &a[(k + 1) * lda], &i__1, &scale, &s);
/*                 tri U at A(0,k+1) */
#line 931 "slansf.f"
		    i__1 = lda + 1;
#line 931 "slansf.f"
		    slassq_(&k, &a[k * lda], &i__1, &scale, &s);
/*                 tri L at A(0,k) */
#line 933 "slansf.f"
		} else {
/*                 A**T is lower */
#line 935 "slansf.f"
		    i__1 = k - 1;
#line 935 "slansf.f"
		    for (j = 1; j <= i__1; ++j) {
#line 936 "slansf.f"
			slassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
/*                    U at A(0,1) */
#line 938 "slansf.f"
		    }
#line 939 "slansf.f"
		    i__1 = *n;
#line 939 "slansf.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 940 "slansf.f"
			slassq_(&k, &a[j * lda], &c__1, &scale, &s);
/*                    k by k rect. at A(0,k+1) */
#line 942 "slansf.f"
		    }
#line 943 "slansf.f"
		    i__1 = k - 2;
#line 943 "slansf.f"
		    for (j = 0; j <= i__1; ++j) {
#line 944 "slansf.f"
			i__2 = k - j - 1;
#line 944 "slansf.f"
			slassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s)
				;
/*                    L at A(0,0) */
#line 946 "slansf.f"
		    }
#line 947 "slansf.f"
		    s += s;
/*                 double s for the off diagonal elements */
#line 949 "slansf.f"
		    i__1 = lda + 1;
#line 949 "slansf.f"
		    slassq_(&k, &a[lda], &i__1, &scale, &s);
/*                 tri L at A(0,1) */
#line 951 "slansf.f"
		    i__1 = lda + 1;
#line 951 "slansf.f"
		    slassq_(&k, a, &i__1, &scale, &s);
/*                 tri U at A(0,0) */
#line 953 "slansf.f"
		}
#line 954 "slansf.f"
	    }
#line 955 "slansf.f"
	}
#line 956 "slansf.f"
	value = scale * sqrt(s);
#line 957 "slansf.f"
    }

#line 959 "slansf.f"
    ret_val = value;
#line 960 "slansf.f"
    return ret_val;

/*     End of SLANSF */

} /* slansf_ */

