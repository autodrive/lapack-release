#line 1 "zsyconvf_rook.f"
/* zsyconvf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "zsyconvf_rook.f"
/* > \brief \b ZSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYCONVF_ROOK( UPLO, WAY, N, A, LDA, E, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, WAY */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > If parameter WAY = 'C': */
/* > ZSYCONVF_ROOK converts the factorization output format used in */
/* > ZSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in ZSYTRF_RK (or ZSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for ZSYTRF_ROOK and */
/* > ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > ZSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in ZSYTRF_RK */
/* > (or ZSYTRF_BK) provided on entry in parameters A and E into */
/* > the factorization output format used in ZSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for ZSYTRF_ROOK and */
/* > ZSYTRF_RK (or ZSYTRF_BK) is the same and is not converted. */
/* > */
/* > ZSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between */
/* > formats used in ZHETRF_ROOK and ZHETRF_RK (or ZHETRF_BK). */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* > */
/* >          1) If WAY ='C': */
/* > */
/* >          On entry, contains factorization details in format used in */
/* >          ZSYTRF_ROOK: */
/* >            a) all elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A and on superdiagonal */
/* >               (or subdiagonal) of A, and */
/* >            b) If UPLO = 'U': multipliers used to obtain factor U */
/* >               in the superdiagonal part of A. */
/* >               If UPLO = 'L': multipliers used to obtain factor L */
/* >               in the superdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF_RK or ZSYTRF_BK: */
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
/* >          ZSYTRF_RK or ZSYTRF_BK: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* >          On exit, contains factorization details in format used in */
/* >          ZSYTRF_ROOK: */
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
/* >          E is COMPLEX*16 array, dimension (N) */
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
/* >          1) by ZSYTRF_ROOK, if WAY ='C'; */
/* >          2) by ZSYTRF_RK (or ZSYTRF_BK), if WAY ='R'. */
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

/* > \ingroup complex16SYcomputational */

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
/* Subroutine */ int zsyconvf_rook__(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
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

#line 235 "zsyconvf_rook.f"
    /* Parameter adjustments */
#line 235 "zsyconvf_rook.f"
    a_dim1 = *lda;
#line 235 "zsyconvf_rook.f"
    a_offset = 1 + a_dim1;
#line 235 "zsyconvf_rook.f"
    a -= a_offset;
#line 235 "zsyconvf_rook.f"
    --e;
#line 235 "zsyconvf_rook.f"
    --ipiv;
#line 235 "zsyconvf_rook.f"

#line 235 "zsyconvf_rook.f"
    /* Function Body */
#line 235 "zsyconvf_rook.f"
    *info = 0;
#line 236 "zsyconvf_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 237 "zsyconvf_rook.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 238 "zsyconvf_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 239 "zsyconvf_rook.f"
	*info = -1;
#line 240 "zsyconvf_rook.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 241 "zsyconvf_rook.f"
	*info = -2;
#line 242 "zsyconvf_rook.f"
    } else if (*n < 0) {
#line 243 "zsyconvf_rook.f"
	*info = -3;
#line 244 "zsyconvf_rook.f"
    } else if (*lda < max(1,*n)) {
#line 245 "zsyconvf_rook.f"
	*info = -5;
#line 247 "zsyconvf_rook.f"
    }
#line 248 "zsyconvf_rook.f"
    if (*info != 0) {
#line 249 "zsyconvf_rook.f"
	i__1 = -(*info);
#line 249 "zsyconvf_rook.f"
	xerbla_("ZSYCONVF_ROOK", &i__1, (ftnlen)13);
#line 250 "zsyconvf_rook.f"
	return 0;
#line 251 "zsyconvf_rook.f"
    }

/*     Quick return if possible */

#line 255 "zsyconvf_rook.f"
    if (*n == 0) {
#line 255 "zsyconvf_rook.f"
	return 0;
#line 255 "zsyconvf_rook.f"
    }

#line 258 "zsyconvf_rook.f"
    if (upper) {

/*        Begin A is UPPER */

#line 262 "zsyconvf_rook.f"
	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 272 "zsyconvf_rook.f"
	    i__ = *n;
#line 273 "zsyconvf_rook.f"
	    e[1].r = 0., e[1].i = 0.;
#line 274 "zsyconvf_rook.f"
	    while(i__ > 1) {
#line 275 "zsyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 276 "zsyconvf_rook.f"
		    i__1 = i__;
#line 276 "zsyconvf_rook.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 276 "zsyconvf_rook.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 277 "zsyconvf_rook.f"
		    i__1 = i__ - 1;
#line 277 "zsyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 278 "zsyconvf_rook.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 278 "zsyconvf_rook.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 279 "zsyconvf_rook.f"
		    --i__;
#line 280 "zsyconvf_rook.f"
		} else {
#line 281 "zsyconvf_rook.f"
		    i__1 = i__;
#line 281 "zsyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 282 "zsyconvf_rook.f"
		}
#line 283 "zsyconvf_rook.f"
		--i__;
#line 284 "zsyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

#line 291 "zsyconvf_rook.f"
	    i__ = *n;
#line 292 "zsyconvf_rook.f"
	    while(i__ >= 1) {
#line 293 "zsyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 299 "zsyconvf_rook.f"
		    ip = ipiv[i__];
#line 300 "zsyconvf_rook.f"
		    if (i__ < *n) {
#line 301 "zsyconvf_rook.f"
			if (ip != i__) {
#line 302 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 302 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 304 "zsyconvf_rook.f"
			}
#line 305 "zsyconvf_rook.f"
		    }

#line 307 "zsyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

#line 314 "zsyconvf_rook.f"
		    ip = -ipiv[i__];
#line 315 "zsyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 316 "zsyconvf_rook.f"
		    if (i__ < *n) {
#line 317 "zsyconvf_rook.f"
			if (ip != i__) {
#line 318 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 318 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 320 "zsyconvf_rook.f"
			}
#line 321 "zsyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 322 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 322 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
#line 324 "zsyconvf_rook.f"
			}
#line 325 "zsyconvf_rook.f"
		    }
#line 326 "zsyconvf_rook.f"
		    --i__;

#line 328 "zsyconvf_rook.f"
		}
#line 329 "zsyconvf_rook.f"
		--i__;
#line 330 "zsyconvf_rook.f"
	    }

#line 332 "zsyconvf_rook.f"
	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

#line 342 "zsyconvf_rook.f"
	    i__ = 1;
#line 343 "zsyconvf_rook.f"
	    while(i__ <= *n) {
#line 344 "zsyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 350 "zsyconvf_rook.f"
		    ip = ipiv[i__];
#line 351 "zsyconvf_rook.f"
		    if (i__ < *n) {
#line 352 "zsyconvf_rook.f"
			if (ip != i__) {
#line 353 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 353 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 355 "zsyconvf_rook.f"
			}
#line 356 "zsyconvf_rook.f"
		    }

#line 358 "zsyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

#line 365 "zsyconvf_rook.f"
		    ++i__;
#line 366 "zsyconvf_rook.f"
		    ip = -ipiv[i__];
#line 367 "zsyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 368 "zsyconvf_rook.f"
		    if (i__ < *n) {
#line 369 "zsyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 370 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 370 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
#line 372 "zsyconvf_rook.f"
			}
#line 373 "zsyconvf_rook.f"
			if (ip != i__) {
#line 374 "zsyconvf_rook.f"
			    i__1 = *n - i__;
#line 374 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 376 "zsyconvf_rook.f"
			}
#line 377 "zsyconvf_rook.f"
		    }

#line 379 "zsyconvf_rook.f"
		}
#line 380 "zsyconvf_rook.f"
		++i__;
#line 381 "zsyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

#line 387 "zsyconvf_rook.f"
	    i__ = *n;
#line 388 "zsyconvf_rook.f"
	    while(i__ > 1) {
#line 389 "zsyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 390 "zsyconvf_rook.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 390 "zsyconvf_rook.f"
		    i__2 = i__;
#line 390 "zsyconvf_rook.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 391 "zsyconvf_rook.f"
		    --i__;
#line 392 "zsyconvf_rook.f"
		}
#line 393 "zsyconvf_rook.f"
		--i__;
#line 394 "zsyconvf_rook.f"
	    }

/*        End A is UPPER */

#line 398 "zsyconvf_rook.f"
	}

#line 400 "zsyconvf_rook.f"
    } else {

/*        Begin A is LOWER */

#line 404 "zsyconvf_rook.f"
	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 413 "zsyconvf_rook.f"
	    i__ = 1;
#line 414 "zsyconvf_rook.f"
	    i__1 = *n;
#line 414 "zsyconvf_rook.f"
	    e[i__1].r = 0., e[i__1].i = 0.;
#line 415 "zsyconvf_rook.f"
	    while(i__ <= *n) {
#line 416 "zsyconvf_rook.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 417 "zsyconvf_rook.f"
		    i__1 = i__;
#line 417 "zsyconvf_rook.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 417 "zsyconvf_rook.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 418 "zsyconvf_rook.f"
		    i__1 = i__ + 1;
#line 418 "zsyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 419 "zsyconvf_rook.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 419 "zsyconvf_rook.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 420 "zsyconvf_rook.f"
		    ++i__;
#line 421 "zsyconvf_rook.f"
		} else {
#line 422 "zsyconvf_rook.f"
		    i__1 = i__;
#line 422 "zsyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 423 "zsyconvf_rook.f"
		}
#line 424 "zsyconvf_rook.f"
		++i__;
#line 425 "zsyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

#line 432 "zsyconvf_rook.f"
	    i__ = 1;
#line 433 "zsyconvf_rook.f"
	    while(i__ <= *n) {
#line 434 "zsyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 440 "zsyconvf_rook.f"
		    ip = ipiv[i__];
#line 441 "zsyconvf_rook.f"
		    if (i__ > 1) {
#line 442 "zsyconvf_rook.f"
			if (ip != i__) {
#line 443 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 443 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 445 "zsyconvf_rook.f"
			}
#line 446 "zsyconvf_rook.f"
		    }

#line 448 "zsyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

#line 455 "zsyconvf_rook.f"
		    ip = -ipiv[i__];
#line 456 "zsyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 457 "zsyconvf_rook.f"
		    if (i__ > 1) {
#line 458 "zsyconvf_rook.f"
			if (ip != i__) {
#line 459 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 459 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 461 "zsyconvf_rook.f"
			}
#line 462 "zsyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 463 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 463 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
#line 465 "zsyconvf_rook.f"
			}
#line 466 "zsyconvf_rook.f"
		    }
#line 467 "zsyconvf_rook.f"
		    ++i__;

#line 469 "zsyconvf_rook.f"
		}
#line 470 "zsyconvf_rook.f"
		++i__;
#line 471 "zsyconvf_rook.f"
	    }

#line 473 "zsyconvf_rook.f"
	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

#line 483 "zsyconvf_rook.f"
	    i__ = *n;
#line 484 "zsyconvf_rook.f"
	    while(i__ >= 1) {
#line 485 "zsyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 491 "zsyconvf_rook.f"
		    ip = ipiv[i__];
#line 492 "zsyconvf_rook.f"
		    if (i__ > 1) {
#line 493 "zsyconvf_rook.f"
			if (ip != i__) {
#line 494 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 494 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 496 "zsyconvf_rook.f"
			}
#line 497 "zsyconvf_rook.f"
		    }

#line 499 "zsyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

#line 506 "zsyconvf_rook.f"
		    --i__;
#line 507 "zsyconvf_rook.f"
		    ip = -ipiv[i__];
#line 508 "zsyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 509 "zsyconvf_rook.f"
		    if (i__ > 1) {
#line 510 "zsyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 511 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 511 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
#line 513 "zsyconvf_rook.f"
			}
#line 514 "zsyconvf_rook.f"
			if (ip != i__) {
#line 515 "zsyconvf_rook.f"
			    i__1 = i__ - 1;
#line 515 "zsyconvf_rook.f"
			    zswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 517 "zsyconvf_rook.f"
			}
#line 518 "zsyconvf_rook.f"
		    }

#line 520 "zsyconvf_rook.f"
		}
#line 521 "zsyconvf_rook.f"
		--i__;
#line 522 "zsyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

#line 528 "zsyconvf_rook.f"
	    i__ = 1;
#line 529 "zsyconvf_rook.f"
	    while(i__ <= *n - 1) {
#line 530 "zsyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 531 "zsyconvf_rook.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 531 "zsyconvf_rook.f"
		    i__2 = i__;
#line 531 "zsyconvf_rook.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 532 "zsyconvf_rook.f"
		    ++i__;
#line 533 "zsyconvf_rook.f"
		}
#line 534 "zsyconvf_rook.f"
		++i__;
#line 535 "zsyconvf_rook.f"
	    }

#line 537 "zsyconvf_rook.f"
	}

/*        End A is LOWER */

#line 541 "zsyconvf_rook.f"
    }
#line 543 "zsyconvf_rook.f"
    return 0;

/*     End of ZSYCONVF_ROOK */

} /* zsyconvf_rook__ */

