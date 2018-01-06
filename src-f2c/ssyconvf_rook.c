#line 1 "ssyconvf_rook.f"
/* ssyconvf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "ssyconvf_rook.f"
/* > \brief \b SSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

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
/* > SSYCONVF_ROOK converts the factorization output format used in */
/* > SSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in SSYTRF_RK (or SSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for SSYTRF_ROOK and */
/* > SSYTRF_RK (or SSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > SSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in SSYTRF_RK */
/* > (or SSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in SSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for SSYTRF_ROOK and */
/* > SSYTRF_RK (or SSYTRF_BK) is the same and is not converted. */
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
/* >          SSYTRF_ROOK: */
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
/* >          SSYTRF_ROOK: */
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
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by SSYTRF_ROOK, if WAY ='C'; */
/* >          2) by SSYTRF_RK (or SSYTRF_BK), if WAY ='R'. */
/* >          The IPIV format is the same for all these routines. */
/* > */
/* >          On exit, is not changed. */
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

/* > \ingroup singleSYcomputational */

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
/* Subroutine */ int ssyconvf_rook__(char *uplo, char *way, integer *n, 
	doublereal *a, integer *lda, doublereal *e, integer *ipiv, integer *
	info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
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

#line 232 "ssyconvf_rook.f"
    /* Parameter adjustments */
#line 232 "ssyconvf_rook.f"
    a_dim1 = *lda;
#line 232 "ssyconvf_rook.f"
    a_offset = 1 + a_dim1;
#line 232 "ssyconvf_rook.f"
    a -= a_offset;
#line 232 "ssyconvf_rook.f"
    --e;
#line 232 "ssyconvf_rook.f"
    --ipiv;
#line 232 "ssyconvf_rook.f"

#line 232 "ssyconvf_rook.f"
    /* Function Body */
#line 232 "ssyconvf_rook.f"
    *info = 0;
#line 233 "ssyconvf_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 234 "ssyconvf_rook.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 235 "ssyconvf_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 236 "ssyconvf_rook.f"
	*info = -1;
#line 237 "ssyconvf_rook.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 238 "ssyconvf_rook.f"
	*info = -2;
#line 239 "ssyconvf_rook.f"
    } else if (*n < 0) {
#line 240 "ssyconvf_rook.f"
	*info = -3;
#line 241 "ssyconvf_rook.f"
    } else if (*lda < max(1,*n)) {
#line 242 "ssyconvf_rook.f"
	*info = -5;
#line 244 "ssyconvf_rook.f"
    }
#line 245 "ssyconvf_rook.f"
    if (*info != 0) {
#line 246 "ssyconvf_rook.f"
	i__1 = -(*info);
#line 246 "ssyconvf_rook.f"
	xerbla_("SSYCONVF_ROOK", &i__1, (ftnlen)13);
#line 247 "ssyconvf_rook.f"
	return 0;
#line 248 "ssyconvf_rook.f"
    }

/*     Quick return if possible */

#line 252 "ssyconvf_rook.f"
    if (*n == 0) {
#line 252 "ssyconvf_rook.f"
	return 0;
#line 252 "ssyconvf_rook.f"
    }

#line 255 "ssyconvf_rook.f"
    if (upper) {

/*        Begin A is UPPER */

#line 259 "ssyconvf_rook.f"
	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 269 "ssyconvf_rook.f"
	    i__ = *n;
#line 270 "ssyconvf_rook.f"
	    e[1] = 0.;
#line 271 "ssyconvf_rook.f"
	    while(i__ > 1) {
#line 272 "ssyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 273 "ssyconvf_rook.f"
		    e[i__] = a[i__ - 1 + i__ * a_dim1];
#line 274 "ssyconvf_rook.f"
		    e[i__ - 1] = 0.;
#line 275 "ssyconvf_rook.f"
		    a[i__ - 1 + i__ * a_dim1] = 0.;
#line 276 "ssyconvf_rook.f"
		    --i__;
#line 277 "ssyconvf_rook.f"
		} else {
#line 278 "ssyconvf_rook.f"
		    e[i__] = 0.;
#line 279 "ssyconvf_rook.f"
		}
#line 280 "ssyconvf_rook.f"
		--i__;
#line 281 "ssyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

#line 288 "ssyconvf_rook.f"
	    i__ = *n;
#line 289 "ssyconvf_rook.f"
	    while(i__ >= 1) {
#line 290 "ssyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 296 "ssyconvf_rook.f"
		    ip = ipiv[i__];
#line 297 "ssyconvf_rook.f"
		    if (i__ < *n) {
#line 298 "ssyconvf_rook.f"
			if (ip != i__) {
#line 299 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 299 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 301 "ssyconvf_rook.f"
			}
#line 302 "ssyconvf_rook.f"
		    }

#line 304 "ssyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

#line 311 "ssyconvf_rook.f"
		    ip = -ipiv[i__];
#line 312 "ssyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 313 "ssyconvf_rook.f"
		    if (i__ < *n) {
#line 314 "ssyconvf_rook.f"
			if (ip != i__) {
#line 315 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 315 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 317 "ssyconvf_rook.f"
			}
#line 318 "ssyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 319 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 319 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
#line 321 "ssyconvf_rook.f"
			}
#line 322 "ssyconvf_rook.f"
		    }
#line 323 "ssyconvf_rook.f"
		    --i__;

#line 325 "ssyconvf_rook.f"
		}
#line 326 "ssyconvf_rook.f"
		--i__;
#line 327 "ssyconvf_rook.f"
	    }

#line 329 "ssyconvf_rook.f"
	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

#line 339 "ssyconvf_rook.f"
	    i__ = 1;
#line 340 "ssyconvf_rook.f"
	    while(i__ <= *n) {
#line 341 "ssyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 347 "ssyconvf_rook.f"
		    ip = ipiv[i__];
#line 348 "ssyconvf_rook.f"
		    if (i__ < *n) {
#line 349 "ssyconvf_rook.f"
			if (ip != i__) {
#line 350 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 350 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 352 "ssyconvf_rook.f"
			}
#line 353 "ssyconvf_rook.f"
		    }

#line 355 "ssyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

#line 362 "ssyconvf_rook.f"
		    ++i__;
#line 363 "ssyconvf_rook.f"
		    ip = -ipiv[i__];
#line 364 "ssyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 365 "ssyconvf_rook.f"
		    if (i__ < *n) {
#line 366 "ssyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 367 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 367 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
#line 369 "ssyconvf_rook.f"
			}
#line 370 "ssyconvf_rook.f"
			if (ip != i__) {
#line 371 "ssyconvf_rook.f"
			    i__1 = *n - i__;
#line 371 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 373 "ssyconvf_rook.f"
			}
#line 374 "ssyconvf_rook.f"
		    }

#line 376 "ssyconvf_rook.f"
		}
#line 377 "ssyconvf_rook.f"
		++i__;
#line 378 "ssyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

#line 384 "ssyconvf_rook.f"
	    i__ = *n;
#line 385 "ssyconvf_rook.f"
	    while(i__ > 1) {
#line 386 "ssyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 387 "ssyconvf_rook.f"
		    a[i__ - 1 + i__ * a_dim1] = e[i__];
#line 388 "ssyconvf_rook.f"
		    --i__;
#line 389 "ssyconvf_rook.f"
		}
#line 390 "ssyconvf_rook.f"
		--i__;
#line 391 "ssyconvf_rook.f"
	    }

/*        End A is UPPER */

#line 395 "ssyconvf_rook.f"
	}

#line 397 "ssyconvf_rook.f"
    } else {

/*        Begin A is LOWER */

#line 401 "ssyconvf_rook.f"
	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 410 "ssyconvf_rook.f"
	    i__ = 1;
#line 411 "ssyconvf_rook.f"
	    e[*n] = 0.;
#line 412 "ssyconvf_rook.f"
	    while(i__ <= *n) {
#line 413 "ssyconvf_rook.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 414 "ssyconvf_rook.f"
		    e[i__] = a[i__ + 1 + i__ * a_dim1];
#line 415 "ssyconvf_rook.f"
		    e[i__ + 1] = 0.;
#line 416 "ssyconvf_rook.f"
		    a[i__ + 1 + i__ * a_dim1] = 0.;
#line 417 "ssyconvf_rook.f"
		    ++i__;
#line 418 "ssyconvf_rook.f"
		} else {
#line 419 "ssyconvf_rook.f"
		    e[i__] = 0.;
#line 420 "ssyconvf_rook.f"
		}
#line 421 "ssyconvf_rook.f"
		++i__;
#line 422 "ssyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

#line 429 "ssyconvf_rook.f"
	    i__ = 1;
#line 430 "ssyconvf_rook.f"
	    while(i__ <= *n) {
#line 431 "ssyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 437 "ssyconvf_rook.f"
		    ip = ipiv[i__];
#line 438 "ssyconvf_rook.f"
		    if (i__ > 1) {
#line 439 "ssyconvf_rook.f"
			if (ip != i__) {
#line 440 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 440 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 442 "ssyconvf_rook.f"
			}
#line 443 "ssyconvf_rook.f"
		    }

#line 445 "ssyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

#line 452 "ssyconvf_rook.f"
		    ip = -ipiv[i__];
#line 453 "ssyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 454 "ssyconvf_rook.f"
		    if (i__ > 1) {
#line 455 "ssyconvf_rook.f"
			if (ip != i__) {
#line 456 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 456 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 458 "ssyconvf_rook.f"
			}
#line 459 "ssyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 460 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 460 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
#line 462 "ssyconvf_rook.f"
			}
#line 463 "ssyconvf_rook.f"
		    }
#line 464 "ssyconvf_rook.f"
		    ++i__;

#line 466 "ssyconvf_rook.f"
		}
#line 467 "ssyconvf_rook.f"
		++i__;
#line 468 "ssyconvf_rook.f"
	    }

#line 470 "ssyconvf_rook.f"
	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

#line 480 "ssyconvf_rook.f"
	    i__ = *n;
#line 481 "ssyconvf_rook.f"
	    while(i__ >= 1) {
#line 482 "ssyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 488 "ssyconvf_rook.f"
		    ip = ipiv[i__];
#line 489 "ssyconvf_rook.f"
		    if (i__ > 1) {
#line 490 "ssyconvf_rook.f"
			if (ip != i__) {
#line 491 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 491 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 493 "ssyconvf_rook.f"
			}
#line 494 "ssyconvf_rook.f"
		    }

#line 496 "ssyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

#line 503 "ssyconvf_rook.f"
		    --i__;
#line 504 "ssyconvf_rook.f"
		    ip = -ipiv[i__];
#line 505 "ssyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 506 "ssyconvf_rook.f"
		    if (i__ > 1) {
#line 507 "ssyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 508 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 508 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
#line 510 "ssyconvf_rook.f"
			}
#line 511 "ssyconvf_rook.f"
			if (ip != i__) {
#line 512 "ssyconvf_rook.f"
			    i__1 = i__ - 1;
#line 512 "ssyconvf_rook.f"
			    sswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 514 "ssyconvf_rook.f"
			}
#line 515 "ssyconvf_rook.f"
		    }

#line 517 "ssyconvf_rook.f"
		}
#line 518 "ssyconvf_rook.f"
		--i__;
#line 519 "ssyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

#line 525 "ssyconvf_rook.f"
	    i__ = 1;
#line 526 "ssyconvf_rook.f"
	    while(i__ <= *n - 1) {
#line 527 "ssyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 528 "ssyconvf_rook.f"
		    a[i__ + 1 + i__ * a_dim1] = e[i__];
#line 529 "ssyconvf_rook.f"
		    ++i__;
#line 530 "ssyconvf_rook.f"
		}
#line 531 "ssyconvf_rook.f"
		++i__;
#line 532 "ssyconvf_rook.f"
	    }

#line 534 "ssyconvf_rook.f"
	}

/*        End A is LOWER */

#line 538 "ssyconvf_rook.f"
    }
#line 540 "ssyconvf_rook.f"
    return 0;

/*     End of SSYCONVF_ROOK */

} /* ssyconvf_rook__ */

