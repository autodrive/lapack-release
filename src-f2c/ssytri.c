#line 1 "ssytri.f"
/* ssytri.f -- translated by f2c (version 20100827).
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

#line 1 "ssytri.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b SSYTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

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
/* > SSYTRI computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > SSYTRF. */
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
/* >          used to obtain the factor U or L as computed by SSYTRF. */
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
/* >          as determined by SSYTRF. */
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

/* > \date November 2011 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssytri_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *info, ftnlen uplo_len)
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


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 157 "ssytri.f"
    /* Parameter adjustments */
#line 157 "ssytri.f"
    a_dim1 = *lda;
#line 157 "ssytri.f"
    a_offset = 1 + a_dim1;
#line 157 "ssytri.f"
    a -= a_offset;
#line 157 "ssytri.f"
    --ipiv;
#line 157 "ssytri.f"
    --work;
#line 157 "ssytri.f"

#line 157 "ssytri.f"
    /* Function Body */
#line 157 "ssytri.f"
    *info = 0;
#line 158 "ssytri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 159 "ssytri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 160 "ssytri.f"
	*info = -1;
#line 161 "ssytri.f"
    } else if (*n < 0) {
#line 162 "ssytri.f"
	*info = -2;
#line 163 "ssytri.f"
    } else if (*lda < max(1,*n)) {
#line 164 "ssytri.f"
	*info = -4;
#line 165 "ssytri.f"
    }
#line 166 "ssytri.f"
    if (*info != 0) {
#line 167 "ssytri.f"
	i__1 = -(*info);
#line 167 "ssytri.f"
	xerbla_("SSYTRI", &i__1, (ftnlen)6);
#line 168 "ssytri.f"
	return 0;
#line 169 "ssytri.f"
    }

/*     Quick return if possible */

#line 173 "ssytri.f"
    if (*n == 0) {
#line 173 "ssytri.f"
	return 0;
#line 173 "ssytri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 178 "ssytri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 182 "ssytri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 183 "ssytri.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 183 "ssytri.f"
		return 0;
#line 183 "ssytri.f"
	    }
#line 185 "ssytri.f"
/* L10: */
#line 185 "ssytri.f"
	}
#line 186 "ssytri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 190 "ssytri.f"
	i__1 = *n;
#line 190 "ssytri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 191 "ssytri.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 191 "ssytri.f"
		return 0;
#line 191 "ssytri.f"
	    }
#line 193 "ssytri.f"
/* L20: */
#line 193 "ssytri.f"
	}
#line 194 "ssytri.f"
    }
#line 195 "ssytri.f"
    *info = 0;

#line 197 "ssytri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 204 "ssytri.f"
	k = 1;
#line 205 "ssytri.f"
L30:

/*        If K > N, exit from loop. */

#line 209 "ssytri.f"
	if (k > *n) {
#line 209 "ssytri.f"
	    goto L40;
#line 209 "ssytri.f"
	}

#line 212 "ssytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 218 "ssytri.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 222 "ssytri.f"
	    if (k > 1) {
#line 223 "ssytri.f"
		i__1 = k - 1;
#line 223 "ssytri.f"
		scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 224 "ssytri.f"
		i__1 = k - 1;
#line 224 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 226 "ssytri.f"
		i__1 = k - 1;
#line 226 "ssytri.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 228 "ssytri.f"
	    }
#line 229 "ssytri.f"
	    kstep = 1;
#line 230 "ssytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 236 "ssytri.f"
	    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
#line 237 "ssytri.f"
	    ak = a[k + k * a_dim1] / t;
#line 238 "ssytri.f"
	    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 239 "ssytri.f"
	    akkp1 = a[k + (k + 1) * a_dim1] / t;
#line 240 "ssytri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 241 "ssytri.f"
	    a[k + k * a_dim1] = akp1 / d__;
#line 242 "ssytri.f"
	    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
#line 243 "ssytri.f"
	    a[k + (k + 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 247 "ssytri.f"
	    if (k > 1) {
#line 248 "ssytri.f"
		i__1 = k - 1;
#line 248 "ssytri.f"
		scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 249 "ssytri.f"
		i__1 = k - 1;
#line 249 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 251 "ssytri.f"
		i__1 = k - 1;
#line 251 "ssytri.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 253 "ssytri.f"
		i__1 = k - 1;
#line 253 "ssytri.f"
		a[k + (k + 1) * a_dim1] -= sdot_(&i__1, &a[k * a_dim1 + 1], &
			c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
#line 255 "ssytri.f"
		i__1 = k - 1;
#line 255 "ssytri.f"
		scopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 256 "ssytri.f"
		i__1 = k - 1;
#line 256 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[(k + 1) * a_dim1 + 1], &c__1, (
			ftnlen)1);
#line 258 "ssytri.f"
		i__1 = k - 1;
#line 258 "ssytri.f"
		a[k + 1 + (k + 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &
			a[(k + 1) * a_dim1 + 1], &c__1);
#line 260 "ssytri.f"
	    }
#line 261 "ssytri.f"
	    kstep = 2;
#line 262 "ssytri.f"
	}

#line 264 "ssytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 265 "ssytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 270 "ssytri.f"
	    i__1 = kp - 1;
#line 270 "ssytri.f"
	    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
		    c__1);
#line 271 "ssytri.f"
	    i__1 = k - kp - 1;
#line 271 "ssytri.f"
	    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * 
		    a_dim1], lda);
#line 272 "ssytri.f"
	    temp = a[k + k * a_dim1];
#line 273 "ssytri.f"
	    a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 274 "ssytri.f"
	    a[kp + kp * a_dim1] = temp;
#line 275 "ssytri.f"
	    if (kstep == 2) {
#line 276 "ssytri.f"
		temp = a[k + (k + 1) * a_dim1];
#line 277 "ssytri.f"
		a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
#line 278 "ssytri.f"
		a[kp + (k + 1) * a_dim1] = temp;
#line 279 "ssytri.f"
	    }
#line 280 "ssytri.f"
	}

#line 282 "ssytri.f"
	k += kstep;
#line 283 "ssytri.f"
	goto L30;
#line 284 "ssytri.f"
L40:

#line 286 "ssytri.f"
	;
#line 286 "ssytri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 293 "ssytri.f"
	k = *n;
#line 294 "ssytri.f"
L50:

/*        If K < 1, exit from loop. */

#line 298 "ssytri.f"
	if (k < 1) {
#line 298 "ssytri.f"
	    goto L60;
#line 298 "ssytri.f"
	}

#line 301 "ssytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 307 "ssytri.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 311 "ssytri.f"
	    if (k < *n) {
#line 312 "ssytri.f"
		i__1 = *n - k;
#line 312 "ssytri.f"
		scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 313 "ssytri.f"
		i__1 = *n - k;
#line 313 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 315 "ssytri.f"
		i__1 = *n - k;
#line 315 "ssytri.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 317 "ssytri.f"
	    }
#line 318 "ssytri.f"
	    kstep = 1;
#line 319 "ssytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 325 "ssytri.f"
	    t = (d__1 = a[k + (k - 1) * a_dim1], abs(d__1));
#line 326 "ssytri.f"
	    ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 327 "ssytri.f"
	    akp1 = a[k + k * a_dim1] / t;
#line 328 "ssytri.f"
	    akkp1 = a[k + (k - 1) * a_dim1] / t;
#line 329 "ssytri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 330 "ssytri.f"
	    a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
#line 331 "ssytri.f"
	    a[k + k * a_dim1] = ak / d__;
#line 332 "ssytri.f"
	    a[k + (k - 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 336 "ssytri.f"
	    if (k < *n) {
#line 337 "ssytri.f"
		i__1 = *n - k;
#line 337 "ssytri.f"
		scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 338 "ssytri.f"
		i__1 = *n - k;
#line 338 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 340 "ssytri.f"
		i__1 = *n - k;
#line 340 "ssytri.f"
		a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 342 "ssytri.f"
		i__1 = *n - k;
#line 342 "ssytri.f"
		a[k + (k - 1) * a_dim1] -= sdot_(&i__1, &a[k + 1 + k * a_dim1]
			, &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 345 "ssytri.f"
		i__1 = *n - k;
#line 345 "ssytri.f"
		scopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 346 "ssytri.f"
		i__1 = *n - k;
#line 346 "ssytri.f"
		ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + (k - 1) * a_dim1]
			, &c__1, (ftnlen)1);
#line 348 "ssytri.f"
		i__1 = *n - k;
#line 348 "ssytri.f"
		a[k - 1 + (k - 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &
			a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 350 "ssytri.f"
	    }
#line 351 "ssytri.f"
	    kstep = 2;
#line 352 "ssytri.f"
	}

#line 354 "ssytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 355 "ssytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 360 "ssytri.f"
	    if (kp < *n) {
#line 360 "ssytri.f"
		i__1 = *n - kp;
#line 360 "ssytri.f"
		sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp *
			 a_dim1], &c__1);
#line 360 "ssytri.f"
	    }
#line 362 "ssytri.f"
	    i__1 = kp - k - 1;
#line 362 "ssytri.f"
	    sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * 
		    a_dim1], lda);
#line 363 "ssytri.f"
	    temp = a[k + k * a_dim1];
#line 364 "ssytri.f"
	    a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 365 "ssytri.f"
	    a[kp + kp * a_dim1] = temp;
#line 366 "ssytri.f"
	    if (kstep == 2) {
#line 367 "ssytri.f"
		temp = a[k + (k - 1) * a_dim1];
#line 368 "ssytri.f"
		a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
#line 369 "ssytri.f"
		a[kp + (k - 1) * a_dim1] = temp;
#line 370 "ssytri.f"
	    }
#line 371 "ssytri.f"
	}

#line 373 "ssytri.f"
	k -= kstep;
#line 374 "ssytri.f"
	goto L50;
#line 375 "ssytri.f"
L60:
#line 376 "ssytri.f"
	;
#line 376 "ssytri.f"
    }

#line 378 "ssytri.f"
    return 0;

/*     End of SSYTRI */

} /* ssytri_ */

