#line 1 "ssyconvf.f"
/* ssyconvf.f -- translated by f2c (version 20100827).
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

#line 1 "ssyconvf.f"
/* > \brief \b SSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONVF( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > SSYCONVF converts the factorization output format used in */
/* > SSYTRF provided on entry in parameter A into the factorization */
/* > output format used in SSYTRF_RK (or SSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also coverts in place details of */
/* > the intechanges stored in IPIV from the format used in SSYTRF into */
/* > the format used in SSYTRF_RK (or SSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > SSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in SSYTRF_RK */
/* > (or SSYTRF_BK) provided on entry in parametes A and E into */
/* > the factorization output format used in SSYTRF that is stored */
/* > on exit in parameter A. It also coverts in place details of */
/* > the intechanges stored in IPIV from the format used in SSYTRF_RK */
/* > (or SSYTRF_BK) into the format used in SSYTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are */
/* >          stored as an upper or lower triangular matrix A. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* >          WAY is CHARACTER*1 */
/* >          = 'C': Convert */
/* >          = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          SSYTRF_RK or SSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          SSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, just a workspace. */
/* > */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          2) If WAY = 'R': */
/* > */
/* >          On entry, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; */
/* >          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* >          On exit, is not changed */
/* > \endverbatim */
/* . */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* > */
/* >          1) If WAY ='C': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF_RK */
/* >          ( or SSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF_RK */
/* >          ( or SSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in SSYTRF. */
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

/* > \ingroup singleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int ssyconvf_(char *uplo, char *way, integer *n, doublereal *
	a, integer *lda, doublereal *e, integer *ipiv, integer *info, ftnlen 
	uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical convert;


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
/*     .. External Functions .. */

/*     .. External Subroutines .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 241 "ssyconvf.f"
    /* Parameter adjustments */
#line 241 "ssyconvf.f"
    a_dim1 = *lda;
#line 241 "ssyconvf.f"
    a_offset = 1 + a_dim1;
#line 241 "ssyconvf.f"
    a -= a_offset;
#line 241 "ssyconvf.f"
    --e;
#line 241 "ssyconvf.f"
    --ipiv;
#line 241 "ssyconvf.f"

#line 241 "ssyconvf.f"
    /* Function Body */
#line 241 "ssyconvf.f"
    *info = 0;
#line 242 "ssyconvf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 243 "ssyconvf.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 244 "ssyconvf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 245 "ssyconvf.f"
	*info = -1;
#line 246 "ssyconvf.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 247 "ssyconvf.f"
	*info = -2;
#line 248 "ssyconvf.f"
    } else if (*n < 0) {
#line 249 "ssyconvf.f"
	*info = -3;
#line 250 "ssyconvf.f"
    } else if (*lda < max(1,*n)) {
#line 251 "ssyconvf.f"
	*info = -5;
#line 253 "ssyconvf.f"
    }
#line 254 "ssyconvf.f"
    if (*info != 0) {
#line 255 "ssyconvf.f"
	i__1 = -(*info);
#line 255 "ssyconvf.f"
	xerbla_("SSYCONVF", &i__1, (ftnlen)8);
#line 256 "ssyconvf.f"
	return 0;
#line 257 "ssyconvf.f"
    }

/*     Quick return if possible */

#line 261 "ssyconvf.f"
    if (*n == 0) {
#line 261 "ssyconvf.f"
	return 0;
#line 261 "ssyconvf.f"
    }

#line 264 "ssyconvf.f"
    if (upper) {

/*        Begin A is UPPER */

#line 268 "ssyconvf.f"
	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 278 "ssyconvf.f"
	    i__ = *n;
#line 279 "ssyconvf.f"
	    e[1] = 0.;
#line 280 "ssyconvf.f"
	    while(i__ > 1) {
#line 281 "ssyconvf.f"
		if (ipiv[i__] < 0) {
#line 282 "ssyconvf.f"
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
#line 283 "ssyconvf.f"
		    e[i__ - 1] = 0.;
#line 284 "ssyconvf.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 285 "ssyconvf.f"
		    --i__;
#line 286 "ssyconvf.f"
		} else {
#line 287 "ssyconvf.f"
		    e[i__] = 0.;
#line 288 "ssyconvf.f"
		}
#line 289 "ssyconvf.f"
		--i__;
#line 290 "ssyconvf.f"
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

#line 297 "ssyconvf.f"
	    i__ = *n;
#line 298 "ssyconvf.f"
	    while(i__ >= 1) {
#line 299 "ssyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 305 "ssyconvf.f"
		    ip = ipiv[i__];
#line 306 "ssyconvf.f"
		    if (i__ < *n) {
#line 307 "ssyconvf.f"
			if (ip != i__) {
#line 308 "ssyconvf.f"
			    i__1 = *n - i__;
#line 308 "ssyconvf.f"
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 310 "ssyconvf.f"
			}
#line 311 "ssyconvf.f"
		    }

#line 313 "ssyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

#line 319 "ssyconvf.f"
		    ip = -ipiv[i__];
#line 320 "ssyconvf.f"
		    if (i__ < *n) {
#line 321 "ssyconvf.f"
			if (ip != i__ - 1) {
#line 322 "ssyconvf.f"
			    i__1 = *n - i__;
#line 322 "ssyconvf.f"
			    sswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
#line 324 "ssyconvf.f"
			}
#line 325 "ssyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

#line 332 "ssyconvf.f"
		    ipiv[i__] = i__;

#line 334 "ssyconvf.f"
		    --i__;

#line 336 "ssyconvf.f"
		}
#line 337 "ssyconvf.f"
		--i__;
#line 338 "ssyconvf.f"
	    }

#line 340 "ssyconvf.f"
	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

#line 350 "ssyconvf.f"
	    i__ = 1;
#line 351 "ssyconvf.f"
	    while(i__ <= *n) {
#line 352 "ssyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 358 "ssyconvf.f"
		    ip = ipiv[i__];
#line 359 "ssyconvf.f"
		    if (i__ < *n) {
#line 360 "ssyconvf.f"
			if (ip != i__) {
#line 361 "ssyconvf.f"
			    i__1 = *n - i__;
#line 361 "ssyconvf.f"
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 363 "ssyconvf.f"
			}
#line 364 "ssyconvf.f"
		    }

#line 366 "ssyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

#line 372 "ssyconvf.f"
		    ++i__;
#line 373 "ssyconvf.f"
		    ip = -ipiv[i__];
#line 374 "ssyconvf.f"
		    if (i__ < *n) {
#line 375 "ssyconvf.f"
			if (ip != i__ - 1) {
#line 376 "ssyconvf.f"
			    i__1 = *n - i__;
#line 376 "ssyconvf.f"
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
#line 378 "ssyconvf.f"
			}
#line 379 "ssyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

#line 386 "ssyconvf.f"
		    ipiv[i__] = ipiv[i__ - 1];

#line 388 "ssyconvf.f"
		}
#line 389 "ssyconvf.f"
		++i__;
#line 390 "ssyconvf.f"
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

#line 396 "ssyconvf.f"
	    i__ = *n;
#line 397 "ssyconvf.f"
	    while(i__ > 1) {
#line 398 "ssyconvf.f"
		if (ipiv[i__] < 0) {
#line 399 "ssyconvf.f"
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
#line 400 "ssyconvf.f"
		    --i__;
#line 401 "ssyconvf.f"
		}
#line 402 "ssyconvf.f"
		--i__;
#line 403 "ssyconvf.f"
	    }

/*        End A is UPPER */

#line 407 "ssyconvf.f"
	}

#line 409 "ssyconvf.f"
    } else {

/*        Begin A is LOWER */

#line 413 "ssyconvf.f"
	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 422 "ssyconvf.f"
	    i__ = 1;
#line 423 "ssyconvf.f"
	    e[*n] = 0.;
#line 424 "ssyconvf.f"
	    while(i__ <= *n) {
#line 425 "ssyconvf.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 426 "ssyconvf.f"
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 427 "ssyconvf.f"
		    e[i__ + 1] = 0.;
#line 428 "ssyconvf.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 429 "ssyconvf.f"
		    ++i__;
#line 430 "ssyconvf.f"
		} else {
#line 431 "ssyconvf.f"
		    e[i__] = 0.;
#line 432 "ssyconvf.f"
		}
#line 433 "ssyconvf.f"
		++i__;
#line 434 "ssyconvf.f"
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

#line 441 "ssyconvf.f"
	    i__ = 1;
#line 442 "ssyconvf.f"
	    while(i__ <= *n) {
#line 443 "ssyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 449 "ssyconvf.f"
		    ip = ipiv[i__];
#line 450 "ssyconvf.f"
		    if (i__ > 1) {
#line 451 "ssyconvf.f"
			if (ip != i__) {
#line 452 "ssyconvf.f"
			    i__1 = i__ - 1;
#line 452 "ssyconvf.f"
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 454 "ssyconvf.f"
			}
#line 455 "ssyconvf.f"
		    }

#line 457 "ssyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

#line 463 "ssyconvf.f"
		    ip = -ipiv[i__];
#line 464 "ssyconvf.f"
		    if (i__ > 1) {
#line 465 "ssyconvf.f"
			if (ip != i__ + 1) {
#line 466 "ssyconvf.f"
			    i__1 = i__ - 1;
#line 466 "ssyconvf.f"
			    sswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 468 "ssyconvf.f"
			}
#line 469 "ssyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

#line 476 "ssyconvf.f"
		    ipiv[i__] = i__;

#line 478 "ssyconvf.f"
		    ++i__;

#line 480 "ssyconvf.f"
		}
#line 481 "ssyconvf.f"
		++i__;
#line 482 "ssyconvf.f"
	    }

#line 484 "ssyconvf.f"
	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

#line 494 "ssyconvf.f"
	    i__ = *n;
#line 495 "ssyconvf.f"
	    while(i__ >= 1) {
#line 496 "ssyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 502 "ssyconvf.f"
		    ip = ipiv[i__];
#line 503 "ssyconvf.f"
		    if (i__ > 1) {
#line 504 "ssyconvf.f"
			if (ip != i__) {
#line 505 "ssyconvf.f"
			    i__1 = i__ - 1;
#line 505 "ssyconvf.f"
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 507 "ssyconvf.f"
			}
#line 508 "ssyconvf.f"
		    }

#line 510 "ssyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

#line 516 "ssyconvf.f"
		    --i__;
#line 517 "ssyconvf.f"
		    ip = -ipiv[i__];
#line 518 "ssyconvf.f"
		    if (i__ > 1) {
#line 519 "ssyconvf.f"
			if (ip != i__ + 1) {
#line 520 "ssyconvf.f"
			    i__1 = i__ - 1;
#line 520 "ssyconvf.f"
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
#line 522 "ssyconvf.f"
			}
#line 523 "ssyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

#line 530 "ssyconvf.f"
		    ipiv[i__] = ipiv[i__ + 1];

#line 532 "ssyconvf.f"
		}
#line 533 "ssyconvf.f"
		--i__;
#line 534 "ssyconvf.f"
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

#line 540 "ssyconvf.f"
	    i__ = 1;
#line 541 "ssyconvf.f"
	    while(i__ <= *n - 1) {
#line 542 "ssyconvf.f"
		if (ipiv[i__] < 0) {
#line 543 "ssyconvf.f"
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 544 "ssyconvf.f"
		    ++i__;
#line 545 "ssyconvf.f"
		}
#line 546 "ssyconvf.f"
		++i__;
#line 547 "ssyconvf.f"
	    }

#line 549 "ssyconvf.f"
	}

/*        End A is LOWER */

#line 553 "ssyconvf.f"
    }
#line 555 "ssyconvf.f"
    return 0;

/*     End of SSYCONVF */

} /* ssyconvf_ */

