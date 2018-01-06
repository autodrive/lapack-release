#line 1 "csyconvf_rook.f"
/* csyconvf_rook.f -- translated by f2c (version 20100827).
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

#line 1 "csyconvf_rook.f"
/* > \brief \b CSYCONVF_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYCONVF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyconv
f_rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYCONVF_ROOK( UPLO, WAY, N, A, LDA, IPIV, E, INFO ) */

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
/* > CSYCONVF_ROOK converts the factorization output format used in */
/* > CSYTRF_ROOK provided on entry in parameter A into the factorization */
/* > output format used in CSYTRF_RK (or CSYTRF_BK) that is stored */
/* > on exit in parameters A and E. IPIV format for CSYTRF_ROOK and */
/* > CSYTRF_RK (or CSYTRF_BK) is the same and is not converted. */
/* > */
/* > If parameter WAY = 'R': */
/* > CSYCONVF_ROOK performs the conversion in reverse direction, i.e. */
/* > converts the factorization output format used in CSYTRF_RK */
/* > (or CSYTRF_BK) provided on entry in parametes A and E into */
/* > the factorization output format used in CSYTRF_ROOK that is stored */
/* > on exit in parameter A. IPIV format for CSYTRF_ROOK and */
/* > CSYTRF_RK (or CSYTRF_BK) is the same and is not converted. */
/* > */
/* > CSYCONVF_ROOK can also convert in Hermitian matrix case, i.e. between */
/* > formats used in CHETRF_ROOK and CHETRF_RK (or CHETRF_BK). */
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
/* >          CSYTRF_ROOK: */
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
/* >          CSYTRF_ROOK: */
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
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          On entry, details of the interchanges and the block */
/* >          structure of D as determined: */
/* >          1) by CSYTRF_ROOK, if WAY ='C'; */
/* >          2) by CSYTRF_RK (or CSYTRF_BK), if WAY ='R'. */
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

/* > \date December 2016 */

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int csyconvf_rook__(char *uplo, char *way, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, 
	integer *info, ftnlen uplo_len, ftnlen way_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, ip, ip2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 235 "csyconvf_rook.f"
    /* Parameter adjustments */
#line 235 "csyconvf_rook.f"
    a_dim1 = *lda;
#line 235 "csyconvf_rook.f"
    a_offset = 1 + a_dim1;
#line 235 "csyconvf_rook.f"
    a -= a_offset;
#line 235 "csyconvf_rook.f"
    --e;
#line 235 "csyconvf_rook.f"
    --ipiv;
#line 235 "csyconvf_rook.f"

#line 235 "csyconvf_rook.f"
    /* Function Body */
#line 235 "csyconvf_rook.f"
    *info = 0;
#line 236 "csyconvf_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 237 "csyconvf_rook.f"
    convert = lsame_(way, "C", (ftnlen)1, (ftnlen)1);
#line 238 "csyconvf_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 239 "csyconvf_rook.f"
	*info = -1;
#line 240 "csyconvf_rook.f"
    } else if (! convert && ! lsame_(way, "R", (ftnlen)1, (ftnlen)1)) {
#line 241 "csyconvf_rook.f"
	*info = -2;
#line 242 "csyconvf_rook.f"
    } else if (*n < 0) {
#line 243 "csyconvf_rook.f"
	*info = -3;
#line 244 "csyconvf_rook.f"
    } else if (*lda < max(1,*n)) {
#line 245 "csyconvf_rook.f"
	*info = -5;
#line 247 "csyconvf_rook.f"
    }
#line 248 "csyconvf_rook.f"
    if (*info != 0) {
#line 249 "csyconvf_rook.f"
	i__1 = -(*info);
#line 249 "csyconvf_rook.f"
	xerbla_("CSYCONVF_ROOK", &i__1, (ftnlen)13);
#line 250 "csyconvf_rook.f"
	return 0;
#line 251 "csyconvf_rook.f"
    }

/*     Quick return if possible */

#line 255 "csyconvf_rook.f"
    if (*n == 0) {
#line 255 "csyconvf_rook.f"
	return 0;
#line 255 "csyconvf_rook.f"
    }

#line 258 "csyconvf_rook.f"
    if (upper) {

/*        Begin A is UPPER */

#line 262 "csyconvf_rook.f"
	if (convert) {

/*           Convert A (A is upper) */


/*           Convert VALUE */

/*           Assign superdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 272 "csyconvf_rook.f"
	    i__ = *n;
#line 273 "csyconvf_rook.f"
	    e[1].r = 0., e[1].i = 0.;
#line 274 "csyconvf_rook.f"
	    while(i__ > 1) {
#line 275 "csyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 276 "csyconvf_rook.f"
		    i__1 = i__;
#line 276 "csyconvf_rook.f"
		    i__2 = i__ - 1 + i__ * a_dim1;
#line 276 "csyconvf_rook.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 277 "csyconvf_rook.f"
		    i__1 = i__ - 1;
#line 277 "csyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 278 "csyconvf_rook.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 278 "csyconvf_rook.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 279 "csyconvf_rook.f"
		    --i__;
#line 280 "csyconvf_rook.f"
		} else {
#line 281 "csyconvf_rook.f"
		    i__1 = i__;
#line 281 "csyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 282 "csyconvf_rook.f"
		}
#line 283 "csyconvf_rook.f"
		--i__;
#line 284 "csyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in factorization order where i decreases from N to 1 */

#line 291 "csyconvf_rook.f"
	    i__ = *n;
#line 292 "csyconvf_rook.f"
	    while(i__ >= 1) {
#line 293 "csyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 299 "csyconvf_rook.f"
		    ip = ipiv[i__];
#line 300 "csyconvf_rook.f"
		    if (i__ < *n) {
#line 301 "csyconvf_rook.f"
			if (ip != i__) {
#line 302 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 302 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 304 "csyconvf_rook.f"
			}
#line 305 "csyconvf_rook.f"
		    }

#line 307 "csyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1) */
/*                 in A(1:i,N-i:N) */

#line 314 "csyconvf_rook.f"
		    ip = -ipiv[i__];
#line 315 "csyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 316 "csyconvf_rook.f"
		    if (i__ < *n) {
#line 317 "csyconvf_rook.f"
			if (ip != i__) {
#line 318 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 318 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &
				    a[ip + (i__ + 1) * a_dim1], lda);
#line 320 "csyconvf_rook.f"
			}
#line 321 "csyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 322 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 322 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], 
				    lda, &a[ip2 + (i__ + 1) * a_dim1], lda);
#line 324 "csyconvf_rook.f"
			}
#line 325 "csyconvf_rook.f"
		    }
#line 326 "csyconvf_rook.f"
		    --i__;

#line 328 "csyconvf_rook.f"
		}
#line 329 "csyconvf_rook.f"
		--i__;
#line 330 "csyconvf_rook.f"
	    }

#line 332 "csyconvf_rook.f"
	} else {

/*           Revert A (A is upper) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of upper part of A */
/*           in reverse factorization order where i increases from 1 to N */

#line 342 "csyconvf_rook.f"
	    i__ = 1;
#line 343 "csyconvf_rook.f"
	    while(i__ <= *n) {
#line 344 "csyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(1:i,N-i:N) */

#line 350 "csyconvf_rook.f"
		    ip = ipiv[i__];
#line 351 "csyconvf_rook.f"
		    if (i__ < *n) {
#line 352 "csyconvf_rook.f"
			if (ip != i__) {
#line 353 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 353 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 355 "csyconvf_rook.f"
			}
#line 356 "csyconvf_rook.f"
		    }

#line 358 "csyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i) */
/*                 in A(1:i,N-i:N) */

#line 365 "csyconvf_rook.f"
		    ++i__;
#line 366 "csyconvf_rook.f"
		    ip = -ipiv[i__];
#line 367 "csyconvf_rook.f"
		    ip2 = -ipiv[i__ - 1];
#line 368 "csyconvf_rook.f"
		    if (i__ < *n) {
#line 369 "csyconvf_rook.f"
			if (ip2 != i__ - 1) {
#line 370 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 370 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip2 + (i__ + 1) * a_dim1], lda, &
				    a[i__ - 1 + (i__ + 1) * a_dim1], lda);
#line 372 "csyconvf_rook.f"
			}
#line 373 "csyconvf_rook.f"
			if (ip != i__) {
#line 374 "csyconvf_rook.f"
			    i__1 = *n - i__;
#line 374 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip + (i__ + 1) * a_dim1], lda, &
				    a[i__ + (i__ + 1) * a_dim1], lda);
#line 376 "csyconvf_rook.f"
			}
#line 377 "csyconvf_rook.f"
		    }

#line 379 "csyconvf_rook.f"
		}
#line 380 "csyconvf_rook.f"
		++i__;
#line 381 "csyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign superdiagonal entries of D from array E to */
/*           superdiagonal entries of A. */

#line 387 "csyconvf_rook.f"
	    i__ = *n;
#line 388 "csyconvf_rook.f"
	    while(i__ > 1) {
#line 389 "csyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 390 "csyconvf_rook.f"
		    i__1 = i__ - 1 + i__ * a_dim1;
#line 390 "csyconvf_rook.f"
		    i__2 = i__;
#line 390 "csyconvf_rook.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 391 "csyconvf_rook.f"
		    --i__;
#line 392 "csyconvf_rook.f"
		}
#line 393 "csyconvf_rook.f"
		--i__;
#line 394 "csyconvf_rook.f"
	    }

/*        End A is UPPER */

#line 398 "csyconvf_rook.f"
	}

#line 400 "csyconvf_rook.f"
    } else {

/*        Begin A is LOWER */

#line 404 "csyconvf_rook.f"
	if (convert) {

/*           Convert A (A is lower) */


/*           Convert VALUE */
/*           Assign subdiagonal entries of D to array E and zero out */
/*           corresponding entries in input storage A */

#line 413 "csyconvf_rook.f"
	    i__ = 1;
#line 414 "csyconvf_rook.f"
	    i__1 = *n;
#line 414 "csyconvf_rook.f"
	    e[i__1].r = 0., e[i__1].i = 0.;
#line 415 "csyconvf_rook.f"
	    while(i__ <= *n) {
#line 416 "csyconvf_rook.f"
		if (i__ < *n && ipiv[i__] < 0) {
#line 417 "csyconvf_rook.f"
		    i__1 = i__;
#line 417 "csyconvf_rook.f"
		    i__2 = i__ + 1 + i__ * a_dim1;
#line 417 "csyconvf_rook.f"
		    e[i__1].r = a[i__2].r, e[i__1].i = a[i__2].i;
#line 418 "csyconvf_rook.f"
		    i__1 = i__ + 1;
#line 418 "csyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 419 "csyconvf_rook.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 419 "csyconvf_rook.f"
		    a[i__1].r = 0., a[i__1].i = 0.;
#line 420 "csyconvf_rook.f"
		    ++i__;
#line 421 "csyconvf_rook.f"
		} else {
#line 422 "csyconvf_rook.f"
		    i__1 = i__;
#line 422 "csyconvf_rook.f"
		    e[i__1].r = 0., e[i__1].i = 0.;
#line 423 "csyconvf_rook.f"
		}
#line 424 "csyconvf_rook.f"
		++i__;
#line 425 "csyconvf_rook.f"
	    }

/*           Convert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in factorization order where i increases from 1 to N */

#line 432 "csyconvf_rook.f"
	    i__ = 1;
#line 433 "csyconvf_rook.f"
	    while(i__ <= *n) {
#line 434 "csyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 440 "csyconvf_rook.f"
		    ip = ipiv[i__];
#line 441 "csyconvf_rook.f"
		    if (i__ > 1) {
#line 442 "csyconvf_rook.f"
			if (ip != i__) {
#line 443 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 443 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 445 "csyconvf_rook.f"
			}
#line 446 "csyconvf_rook.f"
		    }

#line 448 "csyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1) */
/*                 in A(i:N,1:i-1) */

#line 455 "csyconvf_rook.f"
		    ip = -ipiv[i__];
#line 456 "csyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 457 "csyconvf_rook.f"
		    if (i__ > 1) {
#line 458 "csyconvf_rook.f"
			if (ip != i__) {
#line 459 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 459 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ + a_dim1], lda, &a[ip + 
				    a_dim1], lda);
#line 461 "csyconvf_rook.f"
			}
#line 462 "csyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 463 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 463 "csyconvf_rook.f"
			    cswap_(&i__1, &a[i__ + 1 + a_dim1], lda, &a[ip2 + 
				    a_dim1], lda);
#line 465 "csyconvf_rook.f"
			}
#line 466 "csyconvf_rook.f"
		    }
#line 467 "csyconvf_rook.f"
		    ++i__;

#line 469 "csyconvf_rook.f"
		}
#line 470 "csyconvf_rook.f"
		++i__;
#line 471 "csyconvf_rook.f"
	    }

#line 473 "csyconvf_rook.f"
	} else {

/*           Revert A (A is lower) */


/*           Revert PERMUTATIONS */

/*           Apply permutaions to submatrices of lower part of A */
/*           in reverse factorization order where i decreases from N to 1 */

#line 483 "csyconvf_rook.f"
	    i__ = *n;
#line 484 "csyconvf_rook.f"
	    while(i__ >= 1) {
#line 485 "csyconvf_rook.f"
		if (ipiv[i__] > 0) {

/*                 1-by-1 pivot interchange */

/*                 Swap rows i and IPIV(i) in A(i:N,1:i-1) */

#line 491 "csyconvf_rook.f"
		    ip = ipiv[i__];
#line 492 "csyconvf_rook.f"
		    if (i__ > 1) {
#line 493 "csyconvf_rook.f"
			if (ip != i__) {
#line 494 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 494 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 496 "csyconvf_rook.f"
			}
#line 497 "csyconvf_rook.f"
		    }

#line 499 "csyconvf_rook.f"
		} else {

/*                 2-by-2 pivot interchange */

/*                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i) */
/*                 in A(i:N,1:i-1) */

#line 506 "csyconvf_rook.f"
		    --i__;
#line 507 "csyconvf_rook.f"
		    ip = -ipiv[i__];
#line 508 "csyconvf_rook.f"
		    ip2 = -ipiv[i__ + 1];
#line 509 "csyconvf_rook.f"
		    if (i__ > 1) {
#line 510 "csyconvf_rook.f"
			if (ip2 != i__ + 1) {
#line 511 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 511 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip2 + a_dim1], lda, &a[i__ + 1 + 
				    a_dim1], lda);
#line 513 "csyconvf_rook.f"
			}
#line 514 "csyconvf_rook.f"
			if (ip != i__) {
#line 515 "csyconvf_rook.f"
			    i__1 = i__ - 1;
#line 515 "csyconvf_rook.f"
			    cswap_(&i__1, &a[ip + a_dim1], lda, &a[i__ + 
				    a_dim1], lda);
#line 517 "csyconvf_rook.f"
			}
#line 518 "csyconvf_rook.f"
		    }

#line 520 "csyconvf_rook.f"
		}
#line 521 "csyconvf_rook.f"
		--i__;
#line 522 "csyconvf_rook.f"
	    }

/*           Revert VALUE */
/*           Assign subdiagonal entries of D from array E to */
/*           subgiagonal entries of A. */

#line 528 "csyconvf_rook.f"
	    i__ = 1;
#line 529 "csyconvf_rook.f"
	    while(i__ <= *n - 1) {
#line 530 "csyconvf_rook.f"
		if (ipiv[i__] < 0) {
#line 531 "csyconvf_rook.f"
		    i__1 = i__ + 1 + i__ * a_dim1;
#line 531 "csyconvf_rook.f"
		    i__2 = i__;
#line 531 "csyconvf_rook.f"
		    a[i__1].r = e[i__2].r, a[i__1].i = e[i__2].i;
#line 532 "csyconvf_rook.f"
		    ++i__;
#line 533 "csyconvf_rook.f"
		}
#line 534 "csyconvf_rook.f"
		++i__;
#line 535 "csyconvf_rook.f"
	    }

#line 537 "csyconvf_rook.f"
	}

/*        End A is LOWER */

#line 541 "csyconvf_rook.f"
    }
#line 543 "csyconvf_rook.f"
    return 0;

/*     End of CSYCONVF_ROOK */

} /* csyconvf_rook__ */

