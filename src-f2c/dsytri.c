#line 1 "dsytri.f"
/* dsytri.f -- translated by f2c (version 20100827).
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

#line 1 "dsytri.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b DSYTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRI computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > DSYTRF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by DSYTRF. */
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
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \date December 2016 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsytri_(char *uplo, integer *n, doublereal *a, integer *
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
    static doublereal akp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     Test the input parameters. */

#line 157 "dsytri.f"
    /* Parameter adjustments */
#line 157 "dsytri.f"
    a_dim1 = *lda;
#line 157 "dsytri.f"
    a_offset = 1 + a_dim1;
#line 157 "dsytri.f"
    a -= a_offset;
#line 157 "dsytri.f"
    --ipiv;
#line 157 "dsytri.f"
    --work;
#line 157 "dsytri.f"

#line 157 "dsytri.f"
    /* Function Body */
#line 157 "dsytri.f"
    *info = 0;
#line 158 "dsytri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 159 "dsytri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 160 "dsytri.f"
	*info = -1;
#line 161 "dsytri.f"
    } else if (*n < 0) {
#line 162 "dsytri.f"
	*info = -2;
#line 163 "dsytri.f"
    } else if (*lda < max(1,*n)) {
#line 164 "dsytri.f"
	*info = -4;
#line 165 "dsytri.f"
    }
#line 166 "dsytri.f"
    if (*info != 0) {
#line 167 "dsytri.f"
	i__1 = -(*info);
#line 167 "dsytri.f"
	xerbla_("DSYTRI", &i__1, (ftnlen)6);
#line 168 "dsytri.f"
	return 0;
#line 169 "dsytri.f"
    }

/*     Quick return if possible */

#line 173 "dsytri.f"
    if (*n == 0) {
#line 173 "dsytri.f"
	return 0;
#line 173 "dsytri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 178 "dsytri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 182 "dsytri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 183 "dsytri.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 183 "dsytri.f"
		return 0;
#line 183 "dsytri.f"
	    }
#line 185 "dsytri.f"
/* L10: */
#line 185 "dsytri.f"
	}
#line 186 "dsytri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 190 "dsytri.f"
	i__1 = *n;
#line 190 "dsytri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 191 "dsytri.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 191 "dsytri.f"
		return 0;
#line 191 "dsytri.f"
	    }
#line 193 "dsytri.f"
/* L20: */
#line 193 "dsytri.f"
	}
#line 194 "dsytri.f"
    }
#line 195 "dsytri.f"
    *info = 0;

#line 197 "dsytri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 204 "dsytri.f"
	k = 1;
#line 205 "dsytri.f"
L30:

/*        If K > N, exit from loop. */

#line 209 "dsytri.f"
	if (k > *n) {
#line 209 "dsytri.f"
	    goto L40;
#line 209 "dsytri.f"
	}

#line 212 "dsytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 218 "dsytri.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 222 "dsytri.f"
	    if (k > 1) {
#line 223 "dsytri.f"
		i__1 = k - 1;
#line 223 "dsytri.f"
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 224 "dsytri.f"
		i__1 = k - 1;
#line 224 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 226 "dsytri.f"
		i__1 = k - 1;
#line 226 "dsytri.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 228 "dsytri.f"
	    }
#line 229 "dsytri.f"
	    kstep = 1;
#line 230 "dsytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 236 "dsytri.f"
	    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
#line 237 "dsytri.f"
	    ak = a[k + k * a_dim1] / t;
#line 238 "dsytri.f"
	    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 239 "dsytri.f"
	    akkp1 = a[k + (k + 1) * a_dim1] / t;
#line 240 "dsytri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 241 "dsytri.f"
	    a[k + k * a_dim1] = akp1 / d__;
#line 242 "dsytri.f"
	    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
#line 243 "dsytri.f"
	    a[k + (k + 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 247 "dsytri.f"
	    if (k > 1) {
#line 248 "dsytri.f"
		i__1 = k - 1;
#line 248 "dsytri.f"
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 249 "dsytri.f"
		i__1 = k - 1;
#line 249 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 251 "dsytri.f"
		i__1 = k - 1;
#line 251 "dsytri.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 253 "dsytri.f"
		i__1 = k - 1;
#line 253 "dsytri.f"
		a[k + (k + 1) * a_dim1] -= ddot_(&i__1, &a[k * a_dim1 + 1], &
			c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
#line 255 "dsytri.f"
		i__1 = k - 1;
#line 255 "dsytri.f"
		dcopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 256 "dsytri.f"
		i__1 = k - 1;
#line 256 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[(k + 1) * a_dim1 + 1], &c__1, (
			ftnlen)1);
#line 258 "dsytri.f"
		i__1 = k - 1;
#line 258 "dsytri.f"
		a[k + 1 + (k + 1) * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &
			a[(k + 1) * a_dim1 + 1], &c__1);
#line 260 "dsytri.f"
	    }
#line 261 "dsytri.f"
	    kstep = 2;
#line 262 "dsytri.f"
	}

#line 264 "dsytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 265 "dsytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 270 "dsytri.f"
	    i__1 = kp - 1;
#line 270 "dsytri.f"
	    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
		    c__1);
#line 271 "dsytri.f"
	    i__1 = k - kp - 1;
#line 271 "dsytri.f"
	    dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * 
		    a_dim1], lda);
#line 272 "dsytri.f"
	    temp = a[k + k * a_dim1];
#line 273 "dsytri.f"
	    a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 274 "dsytri.f"
	    a[kp + kp * a_dim1] = temp;
#line 275 "dsytri.f"
	    if (kstep == 2) {
#line 276 "dsytri.f"
		temp = a[k + (k + 1) * a_dim1];
#line 277 "dsytri.f"
		a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
#line 278 "dsytri.f"
		a[kp + (k + 1) * a_dim1] = temp;
#line 279 "dsytri.f"
	    }
#line 280 "dsytri.f"
	}

#line 282 "dsytri.f"
	k += kstep;
#line 283 "dsytri.f"
	goto L30;
#line 284 "dsytri.f"
L40:

#line 286 "dsytri.f"
	;
#line 286 "dsytri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 293 "dsytri.f"
	k = *n;
#line 294 "dsytri.f"
L50:

/*        If K < 1, exit from loop. */

#line 298 "dsytri.f"
	if (k < 1) {
#line 298 "dsytri.f"
	    goto L60;
#line 298 "dsytri.f"
	}

#line 301 "dsytri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 307 "dsytri.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 311 "dsytri.f"
	    if (k < *n) {
#line 312 "dsytri.f"
		i__1 = *n - k;
#line 312 "dsytri.f"
		dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 313 "dsytri.f"
		i__1 = *n - k;
#line 313 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 315 "dsytri.f"
		i__1 = *n - k;
#line 315 "dsytri.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 317 "dsytri.f"
	    }
#line 318 "dsytri.f"
	    kstep = 1;
#line 319 "dsytri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 325 "dsytri.f"
	    t = (d__1 = a[k + (k - 1) * a_dim1], abs(d__1));
#line 326 "dsytri.f"
	    ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 327 "dsytri.f"
	    akp1 = a[k + k * a_dim1] / t;
#line 328 "dsytri.f"
	    akkp1 = a[k + (k - 1) * a_dim1] / t;
#line 329 "dsytri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 330 "dsytri.f"
	    a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
#line 331 "dsytri.f"
	    a[k + k * a_dim1] = ak / d__;
#line 332 "dsytri.f"
	    a[k + (k - 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 336 "dsytri.f"
	    if (k < *n) {
#line 337 "dsytri.f"
		i__1 = *n - k;
#line 337 "dsytri.f"
		dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 338 "dsytri.f"
		i__1 = *n - k;
#line 338 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 340 "dsytri.f"
		i__1 = *n - k;
#line 340 "dsytri.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 342 "dsytri.f"
		i__1 = *n - k;
#line 342 "dsytri.f"
		a[k + (k - 1) * a_dim1] -= ddot_(&i__1, &a[k + 1 + k * a_dim1]
			, &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 345 "dsytri.f"
		i__1 = *n - k;
#line 345 "dsytri.f"
		dcopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 346 "dsytri.f"
		i__1 = *n - k;
#line 346 "dsytri.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + (k - 1) * a_dim1]
			, &c__1, (ftnlen)1);
#line 348 "dsytri.f"
		i__1 = *n - k;
#line 348 "dsytri.f"
		a[k - 1 + (k - 1) * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &
			a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 350 "dsytri.f"
	    }
#line 351 "dsytri.f"
	    kstep = 2;
#line 352 "dsytri.f"
	}

#line 354 "dsytri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 355 "dsytri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 360 "dsytri.f"
	    if (kp < *n) {
#line 360 "dsytri.f"
		i__1 = *n - kp;
#line 360 "dsytri.f"
		dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp *
			 a_dim1], &c__1);
#line 360 "dsytri.f"
	    }
#line 362 "dsytri.f"
	    i__1 = kp - k - 1;
#line 362 "dsytri.f"
	    dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * 
		    a_dim1], lda);
#line 363 "dsytri.f"
	    temp = a[k + k * a_dim1];
#line 364 "dsytri.f"
	    a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 365 "dsytri.f"
	    a[kp + kp * a_dim1] = temp;
#line 366 "dsytri.f"
	    if (kstep == 2) {
#line 367 "dsytri.f"
		temp = a[k + (k - 1) * a_dim1];
#line 368 "dsytri.f"
		a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
#line 369 "dsytri.f"
		a[kp + (k - 1) * a_dim1] = temp;
#line 370 "dsytri.f"
	    }
#line 371 "dsytri.f"
	}

#line 373 "dsytri.f"
	k -= kstep;
#line 374 "dsytri.f"
	goto L50;
#line 375 "dsytri.f"
L60:
#line 376 "dsytri.f"
	;
#line 376 "dsytri.f"
    }

#line 378 "dsytri.f"
    return 0;

/*     End of DSYTRI */

} /* dsytri_ */

