#line 1 "csyconvf.f"
/* csyconvf.f -- translated by f2c (version 20100827).
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

#line 1 "csyconvf.f"
/* > \brief \b CSYCONVF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONVF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
f.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONVF( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > CSYCONVF converts the factorization output format used in */
/* > CSYTRF provided on entry in parameter A into the factorization */
/* > output format used in CSYTRF_RK (or CSYTRF_BK) that is stored */
/* > on exit in parameters A and E. It also coverts in place details of */
/* > the intechanges stored in IPIV from the format used in CSYTRF into */
/* > the format used in CSYTRF_RK (or CSYTRF_BK). */
/* > */
/* > If parameter WAY = 'R': */
/* > CSYCONVF performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in CSYTRF_RK */
/* > (or CSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in CSYTRF that is stored */
/* > on exit in parameter A. It also coverts in place details of */
/* > the intechanges stored in IPIV from the format used in CSYTRF_RK */
/* > (or CSYTRF_BK) into the format used in CSYTRF. */
/* > */
/* > CSYCONVF can also convert in Hermitian matrix case, i.e. between */
/* > formats used in CHETRF and CHETRF_RK (or CHETRF_BK). */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          CSYTRF: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF_RK or CSYTRF_BK: */
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
/* >          CSYTRF_RK or CSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          CSYTRF: */
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
/* >          E is COMPLEX array, dimension (N) */
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
/* >          structure of D in the format used in CSYTRF. */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF_RK */
/* >          ( or CSYTRF_BK). */
/* > */
/* >          1) If WAY ='R': */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF_RK */
/* >          ( or CSYTRF_BK). */
/* >          On exit, details of the interchanges and the block */
/* >          structure of D in the format used in CSYTRF. */
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

/* > \date November 2017 */

/* > \ingroup complexSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2017,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* > \endverbatim */
/*  ===================================================================== */
/* Subroutine */ int csyconvf_(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical convert;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 244 "csyconvf.f"
    /* Parameter adjustments */
#line 244 "csyconvf.f"
    a_dim1 = *lda;
#line 244 "csyconvf.f"
    a_offset = 1 + a_dim1;
#line 244 "csyconvf.f"
    a -= a_offset;
#line 244 "csyconvf.f"
    --e;
#line 244 "csyconvf.f"
    --ipiv;
#line 244 "csyconvf.f"

#line 244 "csyconvf.f"
    /* Function Body */
#line 244 "csyconvf.f"
    *info = 0;
#line 245 "csyconvf.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 246 "csyconvf.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 247 "csyconvf.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 248 "csyconvf.f"
	*info = -1;
#line 249 "csyconvf.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 250 "csyconvf.f"
	*info = -2;
#line 251 "csyconvf.f"
    } else if (*n < 0) {
#line 252 "csyconvf.f"
	*info = -3;
#line 253 "csyconvf.f"
    } else if (*lda < max(1,*n)) {
#line 254 "csyconvf.f"
	*info = -5;
#line 256 "csyconvf.f"
    }
#line 257 "csyconvf.f"
    if (*info != 0) {
#line 258 "csyconvf.f"
	i__1 = -(*info);
#line 258 "csyconvf.f"
	xerbla_("CSYCONVF", &i__1, (ftnlen)8);
#line 259 "csyconvf.f"
	return 0;
#line 260 "csyconvf.f"
    }

/*     Quick return if possible */

#line 264 "csyconvf.f"
    if (*n == 0) {
#line 264 "csyconvf.f"
	return 0;
#line 264 "csyconvf.f"
    }

#line 267 "csyconvf.f"
    if (upper) {

/*        Begin A is UPPER */

#line 271 "csyconvf.f"
	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 281 "csyconvf.f"
	    i__ = *n;
#line 282 "csyconvf.f"
	    e[1].r = 0., e[1].i = 0.;
#line 283 "csyconvf.f"
	    while(i__ > 1) {
#line 284 "csyconvf.f"
		if (ipiv[i__] < 0) {
#line 285 "csyconvf.f"
		    i__1 = i__;
#line 285 "csyconvf.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 285 "csyconvf.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 286 "csyconvf.f"
		    i__1 = i__ - 1;
#line 286 "csyconvf.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 287 "csyconvf.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 287 "csyconvf.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 288 "csyconvf.f"
		    --i__;
#line 289 "csyconvf.f"
		} else {
#line 290 "csyconvf.f"
		    i__1 = i__;
#line 290 "csyconvf.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 291 "csyconvf.f"
		}
#line 292 "csyconvf.f"
		--i__;
#line 293 "csyconvf.f"
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

#line 300 "csyconvf.f"
	    i__ = *n;
#line 301 "csyconvf.f"
	    while(i__ >= 1) {
#line 302 "csyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 308 "csyconvf.f"
		    ip = ipiv[i__];
#line 309 "csyconvf.f"
		    if (i__ < *n) {
#line 310 "csyconvf.f"
			if (ip != i__) {
#line 311 "csyconvf.f"
			    i__1 = *n - i__;
#line 311 "csyconvf.f"
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 313 "csyconvf.f"
			}
#line 314 "csyconvf.f"
		    }

#line 316 "csyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

#line 322 "csyconvf.f"
		    ip = -ipiv[i__];
#line 323 "csyconvf.f"
		    if (i__ < *n) {
#line 324 "csyconvf.f"
			if (ip != i__ - 1) {
#line 325 "csyconvf.f"
			    i__1 = *n - i__;
#line 325 "csyconvf.f"
			    cswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip + (i__ + 1) * a_dim1], lda);
#line 327 "csyconvf.f"
			}
#line 328 "csyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

#line 335 "csyconvf.f"
		    ipiv[i__] = i__;

#line 337 "csyconvf.f"
		    --i__;

#line 339 "csyconvf.f"
		}
#line 340 "csyconvf.f"
		--i__;
#line 341 "csyconvf.f"
	    }

#line 343 "csyconvf.f"
	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

#line 353 "csyconvf.f"
	    i__ = 1;
#line 354 "csyconvf.f"
	    while(i__ <= *n) {
#line 355 "csyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 361 "csyconvf.f"
		    ip = ipiv[i__];
#line 362 "csyconvf.f"
		    if (i__ < *n) {
#line 363 "csyconvf.f"
			if (ip != i__) {
#line 364 "csyconvf.f"
			    i__1 = *n - i__;
#line 364 "csyconvf.f"
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 366 "csyconvf.f"
			}
#line 367 "csyconvf.f"
		    }

#line 369 "csyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i) in A(1:i,N-i:N) */

#line 375 "csyconvf.f"
		    ++i__;
#line 376 "csyconvf.f"
		    ip = -ipiv[i__];
#line 377 "csyconvf.f"
		    if (i__ < *n) {
#line 378 "csyconvf.f"
			if (ip != i__ - 1) {
#line 379 "csyconvf.f"
			    i__1 = *n - i__;
#line 379 "csyconvf.f"
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
#line 381 "csyconvf.f"
			}
#line 382 "csyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i-1 and IPIV(i-1), */
/*                 so this should be recorded in two consecutive entries */
/*                 in IPIV format for *SYTRF */

#line 389 "csyconvf.f"
		    ipiv[i__] = ipiv[i__ - 1];

#line 391 "csyconvf.f"
		}
#line 392 "csyconvf.f"
		++i__;
#line 393 "csyconvf.f"
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

#line 399 "csyconvf.f"
	    i__ = *n;
#line 400 "csyconvf.f"
	    while(i__ > 1) {
#line 401 "csyconvf.f"
		if (ipiv[i__] < 0) {
#line 402 "csyconvf.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 402 "csyconvf.f"
		    i__2 = i__;
#line 402 "csyconvf.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 403 "csyconvf.f"
		    --i__;
#line 404 "csyconvf.f"
		}
#line 405 "csyconvf.f"
		--i__;
#line 406 "csyconvf.f"
	    }

/*        End A is UPPER */

#line 410 "csyconvf.f"
	}

#line 412 "csyconvf.f"
    } else {

/*        Begin A is LOWER */

#line 416 "csyconvf.f"
	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 425 "csyconvf.f"
	    i__ = 1;
#line 426 "csyconvf.f"
	    i__1 = *n;
#line 426 "csyconvf.f"
	    e[i__1].r = 0., e[i__1].i = 0.;
#line 427 "csyconvf.f"
	    while(i__ <= *n) {
#line 428 "csyconvf.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 429 "csyconvf.f"
		    i__1 = i__;
#line 429 "csyconvf.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 429 "csyconvf.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 430 "csyconvf.f"
		    i__1 = i__ + 1;
#line 430 "csyconvf.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 431 "csyconvf.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 431 "csyconvf.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 432 "csyconvf.f"
		    ++i__;
#line 433 "csyconvf.f"
		} else {
#line 434 "csyconvf.f"
		    i__1 = i__;
#line 434 "csyconvf.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 435 "csyconvf.f"
		}
#line 436 "csyconvf.f"
		++i__;
#line 437 "csyconvf.f"
	    }

/*           Convert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of lower part of A */
/*           in factorization order where k increases from 1 to N */

#line 444 "csyconvf.f"
	    i__ = 1;
#line 445 "csyconvf.f"
	    while(i__ <= *n) {
#line 446 "csyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 452 "csyconvf.f"
		    ip = ipiv[i__];
#line 453 "csyconvf.f"
		    if (i__ > 1) {
#line 454 "csyconvf.f"
			if (ip != i__) {
#line 455 "csyconvf.f"
			    i__1 = i__ - 1;
#line 455 "csyconvf.f"
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 457 "csyconvf.f"
			}
#line 458 "csyconvf.f"
		    }

#line 460 "csyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

#line 466 "csyconvf.f"
		    ip = -ipiv[i__];
#line 467 "csyconvf.f"
		    if (i__ > 1) {
#line 468 "csyconvf.f"
			if (ip != i__ + 1) {
#line 469 "csyconvf.f"
			    i__1 = i__ - 1;
#line 469 "csyconvf.f"
			    cswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 471 "csyconvf.f"
			}
#line 472 "csyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is no interchnge of rows i and and IPIV(i), */
/*                 so this should be reflected in IPIV format for */
/*                 *SYTRF_RK ( or *SYTRF_BK) */

#line 479 "csyconvf.f"
		    ipiv[i__] = i__;

#line 481 "csyconvf.f"
		    ++i__;

#line 483 "csyconvf.f"
		}
#line 484 "csyconvf.f"
		++i__;
#line 485 "csyconvf.f"
	    }

#line 487 "csyconvf.f"
	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS and IPIV */

/*           Apply permutaions to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

#line 497 "csyconvf.f"
	    i__ = *n;
#line 498 "csyconvf.f"
	    while(i__ >= 1) {
#line 499 "csyconvf.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 505 "csyconvf.f"
		    ip = ipiv[i__];
#line 506 "csyconvf.f"
		    if (i__ > 1) {
#line 507 "csyconvf.f"
			if (ip != i__) {
#line 508 "csyconvf.f"
			    i__1 = i__ - 1;
#line 508 "csyconvf.f"
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 510 "csyconvf.f"
			}
#line 511 "csyconvf.f"
		    }

#line 513 "csyconvf.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i) in A(i:N,1:i-1) */

#line 519 "csyconvf.f"
		    --i__;
#line 520 "csyconvf.f"
		    ip = -ipiv[i__];
#line 521 "csyconvf.f"
		    if (i__ > 1) {
#line 522 "csyconvf.f"
			if (ip != i__ + 1) {
#line 523 "csyconvf.f"
			    i__1 = i__ - 1;
#line 523 "csyconvf.f"
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
#line 525 "csyconvf.f"
			}
#line 526 "csyconvf.f"
		    }

/*                 Convert IPIV */
/*                 There is one interchange of rows i+1 and IPIV(i+1), */
/*                 so this should be recorded in consecutive entries */
/*                 in IPIV format for *SYTRF */

#line 533 "csyconvf.f"
		    ipiv[i__] = ipiv[i__ + 1];

#line 535 "csyconvf.f"
		}
#line 536 "csyconvf.f"
		--i__;
#line 537 "csyconvf.f"
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

#line 543 "csyconvf.f"
	    i__ = 1;
#line 544 "csyconvf.f"
	    while(i__ <= *n - 1) {
#line 545 "csyconvf.f"
		if (ipiv[i__] < 0) {
#line 546 "csyconvf.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 546 "csyconvf.f"
		    i__2 = i__;
#line 546 "csyconvf.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 547 "csyconvf.f"
		    ++i__;
#line 548 "csyconvf.f"
		}
#line 549 "csyconvf.f"
		++i__;
#line 550 "csyconvf.f"
	    }

#line 552 "csyconvf.f"
	}

/*        End A is LOWER */

#line 556 "csyconvf.f"
    }
#line 558 "csyconvf.f"
    return 0;

/*     End of CSYCONVF */

} /* csyconvf_ */

