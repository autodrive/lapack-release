#line 1 "ssytri_rook.f"
/* ssytri_rook.f -- translated by f2c (version 20100827).
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

#line 1 "ssytri_rook.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b SSYTRI_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRI_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRI_ROOK computes the inverse of a real symmetric */
/* > matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* > computed by SSYTRF_ROOK. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
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
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by SSYTRF_ROOK. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
/* >          matrix.  If UPLO = 'U', the upper triangular part of the */
/* >          inverse is formed and the part of A below the diagonal is not */
/* >          referenced; if UPLO = 'L' the lower triangular part of the */
/* >          inverse is formed and the part of A above the diagonal is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by SSYTRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its */
/* >               inverse could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup realSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >   April 2012, Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int ssytri_rook__(char *uplo, integer *n, doublereal *a, 
	integer *lda, integer *ipiv, doublereal *work, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer k;
    static doublereal t, ak;
    static integer kp;
    static doublereal akp1, temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), ssymv_(char *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


/*  -- LAPACK computational routine (version 3.4.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 172 "ssytri_rook.f"
    /* Parameter adjustments */
#line 172 "ssytri_rook.f"
    a_dim1 = *lda;
#line 172 "ssytri_rook.f"
    a_offset = 1 + a_dim1;
#line 172 "ssytri_rook.f"
    a -= a_offset;
#line 172 "ssytri_rook.f"
    --ipiv;
#line 172 "ssytri_rook.f"
    --work;
#line 172 "ssytri_rook.f"

#line 172 "ssytri_rook.f"
    /* Function Body */
#line 172 "ssytri_rook.f"
    *info = 0;
#line 173 "ssytri_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 174 "ssytri_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 175 "ssytri_rook.f"
	*info = -1;
#line 176 "ssytri_rook.f"
    } else if (*n < 0) {
#line 177 "ssytri_rook.f"
	*info = -2;
#line 178 "ssytri_rook.f"
    } else if (*lda < max(1,*n)) {
#line 179 "ssytri_rook.f"
	*info = -4;
#line 180 "ssytri_rook.f"
    }
#line 181 "ssytri_rook.f"
    if (*info != 0) {
#line 182 "ssytri_rook.f"
	i__1 = -(*info);
#line 182 "ssytri_rook.f"
	xerbla_("SSYTRI_ROOK", &i__1, (ftnlen)11);
#line 183 "ssytri_rook.f"
	return 0;
#line 184 "ssytri_rook.f"
    }

/*     Quick return if possible */

#line 188 "ssytri_rook.f"
    if (*n == 0) {
#line 188 "ssytri_rook.f"
	return 0;
#line 188 "ssytri_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 193 "ssytri_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 197 "ssytri_rook.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 198 "ssytri_rook.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 198 "ssytri_rook.f"
		return 0;
#line 198 "ssytri_rook.f"
	    }
#line 200 "ssytri_rook.f"
/* L10: */
#line 200 "ssytri_rook.f"
	}
#line 201 "ssytri_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 205 "ssytri_rook.f"
	i__1 = *n;
#line 205 "ssytri_rook.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 206 "ssytri_rook.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 206 "ssytri_rook.f"
		return 0;
#line 206 "ssytri_rook.f"
	    }
#line 208 "ssytri_rook.f"
/* L20: */
#line 208 "ssytri_rook.f"
	}
#line 209 "ssytri_rook.f"
    }
#line 210 "ssytri_rook.f"
    *info = 0;

#line 212 "ssytri_rook.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 219 "ssytri_rook.f"
	k = 1;
#line 220 "ssytri_rook.f"
L30:

/*        If K > N, exit from loop. */

#line 224 "ssytri_rook.f"
	if (k > *n) {
#line 224 "ssytri_rook.f"
	    goto L40;
#line 224 "ssytri_rook.f"
	}

#line 227 "ssytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 233 "ssytri_rook.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 237 "ssytri_rook.f"
	    if (k > 1) {
#line 238 "ssytri_rook.f"
		i__1 = k - 1;
#line 238 "ssytri_rook.f"
		scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 239 "ssytri_rook.f"
		i__1 = k - 1;
#line 239 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 241 "ssytri_rook.f"
		i__1 = k - 1;
#line 241 "ssytri_rook.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 243 "ssytri_rook.f"
	    }
#line 244 "ssytri_rook.f"
	    kstep = 1;
#line 245 "ssytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 251 "ssytri_rook.f"
	    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
#line 252 "ssytri_rook.f"
	    ak = a[k + k * a_dim1] / t;
#line 253 "ssytri_rook.f"
	    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 254 "ssytri_rook.f"
	    akkp1 = a[k + (k + 1) * a_dim1] / t;
#line 255 "ssytri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 256 "ssytri_rook.f"
	    a[k + k * a_dim1] = akp1 / d__;
#line 257 "ssytri_rook.f"
	    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
#line 258 "ssytri_rook.f"
	    a[k + (k + 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 262 "ssytri_rook.f"
	    if (k > 1) {
#line 263 "ssytri_rook.f"
		i__1 = k - 1;
#line 263 "ssytri_rook.f"
		scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 264 "ssytri_rook.f"
		i__1 = k - 1;
#line 264 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 266 "ssytri_rook.f"
		i__1 = k - 1;
#line 266 "ssytri_rook.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 268 "ssytri_rook.f"
		i__1 = k - 1;
#line 268 "ssytri_rook.f"
		a[k + (k + 1) * a_dim1] -= sdot_(&i__1, &a[k * a_dim1 + 1], &
			c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
#line 270 "ssytri_rook.f"
		i__1 = k - 1;
#line 270 "ssytri_rook.f"
		scopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 271 "ssytri_rook.f"
		i__1 = k - 1;
#line 271 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[(k + 1) * a_dim1 + 1], &c__1, (
			ftnlen)1);
#line 273 "ssytri_rook.f"
		i__1 = k - 1;
#line 273 "ssytri_rook.f"
		a[k + 1 + (k + 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &
			a[(k + 1) * a_dim1 + 1], &c__1);
#line 275 "ssytri_rook.f"
	    }
#line 276 "ssytri_rook.f"
	    kstep = 2;
#line 277 "ssytri_rook.f"
	}

#line 279 "ssytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 284 "ssytri_rook.f"
	    kp = ipiv[k];
#line 285 "ssytri_rook.f"
	    if (kp != k) {
#line 286 "ssytri_rook.f"
		if (kp > 1) {
#line 286 "ssytri_rook.f"
		    i__1 = kp - 1;
#line 286 "ssytri_rook.f"
		    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 286 "ssytri_rook.f"
		}
#line 288 "ssytri_rook.f"
		i__1 = k - kp - 1;
#line 288 "ssytri_rook.f"
		sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 289 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 290 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 291 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 292 "ssytri_rook.f"
	    }
#line 293 "ssytri_rook.f"
	} else {

/*           Interchange rows and columns K and K+1 with -IPIV(K) and */
/*           -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1) */

#line 298 "ssytri_rook.f"
	    kp = -ipiv[k];
#line 299 "ssytri_rook.f"
	    if (kp != k) {
#line 300 "ssytri_rook.f"
		if (kp > 1) {
#line 300 "ssytri_rook.f"
		    i__1 = kp - 1;
#line 300 "ssytri_rook.f"
		    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 300 "ssytri_rook.f"
		}
#line 302 "ssytri_rook.f"
		i__1 = k - kp - 1;
#line 302 "ssytri_rook.f"
		sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);

#line 304 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 305 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 306 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 307 "ssytri_rook.f"
		temp = a[k + (k + 1) * a_dim1];
#line 308 "ssytri_rook.f"
		a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
#line 309 "ssytri_rook.f"
		a[kp + (k + 1) * a_dim1] = temp;
#line 310 "ssytri_rook.f"
	    }

#line 312 "ssytri_rook.f"
	    ++k;
#line 313 "ssytri_rook.f"
	    kp = -ipiv[k];
#line 314 "ssytri_rook.f"
	    if (kp != k) {
#line 315 "ssytri_rook.f"
		if (kp > 1) {
#line 315 "ssytri_rook.f"
		    i__1 = kp - 1;
#line 315 "ssytri_rook.f"
		    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 315 "ssytri_rook.f"
		}
#line 317 "ssytri_rook.f"
		i__1 = k - kp - 1;
#line 317 "ssytri_rook.f"
		sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 318 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 319 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 320 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 321 "ssytri_rook.f"
	    }
#line 322 "ssytri_rook.f"
	}

#line 324 "ssytri_rook.f"
	++k;
#line 325 "ssytri_rook.f"
	goto L30;
#line 326 "ssytri_rook.f"
L40:

#line 328 "ssytri_rook.f"
	;
#line 328 "ssytri_rook.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 335 "ssytri_rook.f"
	k = *n;
#line 336 "ssytri_rook.f"
L50:

/*        If K < 1, exit from loop. */

#line 340 "ssytri_rook.f"
	if (k < 1) {
#line 340 "ssytri_rook.f"
	    goto L60;
#line 340 "ssytri_rook.f"
	}

#line 343 "ssytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 349 "ssytri_rook.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 353 "ssytri_rook.f"
	    if (k < *n) {
#line 354 "ssytri_rook.f"
		i__1 = *n - k;
#line 354 "ssytri_rook.f"
		scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 355 "ssytri_rook.f"
		i__1 = *n - k;
#line 355 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 357 "ssytri_rook.f"
		i__1 = *n - k;
#line 357 "ssytri_rook.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 359 "ssytri_rook.f"
	    }
#line 360 "ssytri_rook.f"
	    kstep = 1;
#line 361 "ssytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 367 "ssytri_rook.f"
	    t = (d__1 = a[k + (k - 1) * a_dim1], abs(d__1));
#line 368 "ssytri_rook.f"
	    ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 369 "ssytri_rook.f"
	    akp1 = a[k + k * a_dim1] / t;
#line 370 "ssytri_rook.f"
	    akkp1 = a[k + (k - 1) * a_dim1] / t;
#line 371 "ssytri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 372 "ssytri_rook.f"
	    a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
#line 373 "ssytri_rook.f"
	    a[k + k * a_dim1] = ak / d__;
#line 374 "ssytri_rook.f"
	    a[k + (k - 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 378 "ssytri_rook.f"
	    if (k < *n) {
#line 379 "ssytri_rook.f"
		i__1 = *n - k;
#line 379 "ssytri_rook.f"
		scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 380 "ssytri_rook.f"
		i__1 = *n - k;
#line 380 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 382 "ssytri_rook.f"
		i__1 = *n - k;
#line 382 "ssytri_rook.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 384 "ssytri_rook.f"
		i__1 = *n - k;
#line 384 "ssytri_rook.f"
		a[k + (k - 1) * a_dim1] -= sdot_(&i__1, &a[k + 1 + k * a_dim1]
			, &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 387 "ssytri_rook.f"
		i__1 = *n - k;
#line 387 "ssytri_rook.f"
		scopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 388 "ssytri_rook.f"
		i__1 = *n - k;
#line 388 "ssytri_rook.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + (k - 1) * a_dim1]
			, &c__1, (ftnlen)1);
#line 390 "ssytri_rook.f"
		i__1 = *n - k;
#line 390 "ssytri_rook.f"
		a[k - 1 + (k - 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &
			a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 392 "ssytri_rook.f"
	    }
#line 393 "ssytri_rook.f"
	    kstep = 2;
#line 394 "ssytri_rook.f"
	}

#line 396 "ssytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 401 "ssytri_rook.f"
	    kp = ipiv[k];
#line 402 "ssytri_rook.f"
	    if (kp != k) {
#line 403 "ssytri_rook.f"
		if (kp < *n) {
#line 403 "ssytri_rook.f"
		    i__1 = *n - kp;
#line 403 "ssytri_rook.f"
		    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 403 "ssytri_rook.f"
		}
#line 405 "ssytri_rook.f"
		i__1 = kp - k - 1;
#line 405 "ssytri_rook.f"
		sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 406 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 407 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 408 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 409 "ssytri_rook.f"
	    }
#line 410 "ssytri_rook.f"
	} else {

/*           Interchange rows and columns K and K-1 with -IPIV(K) and */
/*           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */

#line 415 "ssytri_rook.f"
	    kp = -ipiv[k];
#line 416 "ssytri_rook.f"
	    if (kp != k) {
#line 417 "ssytri_rook.f"
		if (kp < *n) {
#line 417 "ssytri_rook.f"
		    i__1 = *n - kp;
#line 417 "ssytri_rook.f"
		    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 417 "ssytri_rook.f"
		}
#line 419 "ssytri_rook.f"
		i__1 = kp - k - 1;
#line 419 "ssytri_rook.f"
		sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);

#line 421 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 422 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 423 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 424 "ssytri_rook.f"
		temp = a[k + (k - 1) * a_dim1];
#line 425 "ssytri_rook.f"
		a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
#line 426 "ssytri_rook.f"
		a[kp + (k - 1) * a_dim1] = temp;
#line 427 "ssytri_rook.f"
	    }

#line 429 "ssytri_rook.f"
	    --k;
#line 430 "ssytri_rook.f"
	    kp = -ipiv[k];
#line 431 "ssytri_rook.f"
	    if (kp != k) {
#line 432 "ssytri_rook.f"
		if (kp < *n) {
#line 432 "ssytri_rook.f"
		    i__1 = *n - kp;
#line 432 "ssytri_rook.f"
		    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 432 "ssytri_rook.f"
		}
#line 434 "ssytri_rook.f"
		i__1 = kp - k - 1;
#line 434 "ssytri_rook.f"
		sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 435 "ssytri_rook.f"
		temp = a[k + k * a_dim1];
#line 436 "ssytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 437 "ssytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 438 "ssytri_rook.f"
	    }
#line 439 "ssytri_rook.f"
	}

#line 441 "ssytri_rook.f"
	--k;
#line 442 "ssytri_rook.f"
	goto L50;
#line 443 "ssytri_rook.f"
L60:
#line 444 "ssytri_rook.f"
	;
#line 444 "ssytri_rook.f"
    }

#line 446 "ssytri_rook.f"
    return 0;

/*     End of SSYTRI_ROOK */

} /* ssytri_rook__ */

