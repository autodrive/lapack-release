#line 1 "ssptri.f"
/* ssptri.f -- translated by f2c (version 20100827).
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

#line 1 "ssptri.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b SSPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPTRI( UPLO, N, AP, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPTRI computes the inverse of a real symmetric indefinite matrix */
/* > A in packed storage using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by SSPTRF. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by SSPTRF, */
/* >          stored as a packed triangular matrix. */
/* > */
/* >          On exit, if INFO = 0, the (symmetric) inverse of the original */
/* >          matrix, stored as a packed triangular matrix. The j-th column */
/* >          of inv(A) is stored in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', */
/* >             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by SSPTRF. */
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

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssptri_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer kc, kp, kx, kpc, npp;
    static doublereal akp1, temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kstep;
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), sspmv_(char *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer kcnext;


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

#line 152 "ssptri.f"
    /* Parameter adjustments */
#line 152 "ssptri.f"
    --work;
#line 152 "ssptri.f"
    --ipiv;
#line 152 "ssptri.f"
    --ap;
#line 152 "ssptri.f"

#line 152 "ssptri.f"
    /* Function Body */
#line 152 "ssptri.f"
    *info = 0;
#line 153 "ssptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 154 "ssptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 155 "ssptri.f"
	*info = -1;
#line 156 "ssptri.f"
    } else if (*n < 0) {
#line 157 "ssptri.f"
	*info = -2;
#line 158 "ssptri.f"
    }
#line 159 "ssptri.f"
    if (*info != 0) {
#line 160 "ssptri.f"
	i__1 = -(*info);
#line 160 "ssptri.f"
	xerbla_("SSPTRI", &i__1, (ftnlen)6);
#line 161 "ssptri.f"
	return 0;
#line 162 "ssptri.f"
    }

/*     Quick return if possible */

#line 166 "ssptri.f"
    if (*n == 0) {
#line 166 "ssptri.f"
	return 0;
#line 166 "ssptri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 171 "ssptri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 175 "ssptri.f"
	kp = *n * (*n + 1) / 2;
#line 176 "ssptri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 177 "ssptri.f"
	    if (ipiv[*info] > 0 && ap[kp] == 0.) {
#line 177 "ssptri.f"
		return 0;
#line 177 "ssptri.f"
	    }
#line 179 "ssptri.f"
	    kp -= *info;
#line 180 "ssptri.f"
/* L10: */
#line 180 "ssptri.f"
	}
#line 181 "ssptri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 185 "ssptri.f"
	kp = 1;
#line 186 "ssptri.f"
	i__1 = *n;
#line 186 "ssptri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 187 "ssptri.f"
	    if (ipiv[*info] > 0 && ap[kp] == 0.) {
#line 187 "ssptri.f"
		return 0;
#line 187 "ssptri.f"
	    }
#line 189 "ssptri.f"
	    kp = kp + *n - *info + 1;
#line 190 "ssptri.f"
/* L20: */
#line 190 "ssptri.f"
	}
#line 191 "ssptri.f"
    }
#line 192 "ssptri.f"
    *info = 0;

#line 194 "ssptri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 201 "ssptri.f"
	k = 1;
#line 202 "ssptri.f"
	kc = 1;
#line 203 "ssptri.f"
L30:

/*        If K > N, exit from loop. */

#line 207 "ssptri.f"
	if (k > *n) {
#line 207 "ssptri.f"
	    goto L50;
#line 207 "ssptri.f"
	}

#line 210 "ssptri.f"
	kcnext = kc + k;
#line 211 "ssptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 217 "ssptri.f"
	    ap[kc + k - 1] = 1. / ap[kc + k - 1];

/*           Compute column K of the inverse. */

#line 221 "ssptri.f"
	    if (k > 1) {
#line 222 "ssptri.f"
		i__1 = k - 1;
#line 222 "ssptri.f"
		scopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 223 "ssptri.f"
		i__1 = k - 1;
#line 223 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kc], &c__1, (ftnlen)1);
#line 225 "ssptri.f"
		i__1 = k - 1;
#line 225 "ssptri.f"
		ap[kc + k - 1] -= sdot_(&i__1, &work[1], &c__1, &ap[kc], &
			c__1);
#line 227 "ssptri.f"
	    }
#line 228 "ssptri.f"
	    kstep = 1;
#line 229 "ssptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 235 "ssptri.f"
	    t = (d__1 = ap[kcnext + k - 1], abs(d__1));
#line 236 "ssptri.f"
	    ak = ap[kc + k - 1] / t;
#line 237 "ssptri.f"
	    akp1 = ap[kcnext + k] / t;
#line 238 "ssptri.f"
	    akkp1 = ap[kcnext + k - 1] / t;
#line 239 "ssptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 240 "ssptri.f"
	    ap[kc + k - 1] = akp1 / d__;
#line 241 "ssptri.f"
	    ap[kcnext + k] = ak / d__;
#line 242 "ssptri.f"
	    ap[kcnext + k - 1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 246 "ssptri.f"
	    if (k > 1) {
#line 247 "ssptri.f"
		i__1 = k - 1;
#line 247 "ssptri.f"
		scopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 248 "ssptri.f"
		i__1 = k - 1;
#line 248 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kc], &c__1, (ftnlen)1);
#line 250 "ssptri.f"
		i__1 = k - 1;
#line 250 "ssptri.f"
		ap[kc + k - 1] -= sdot_(&i__1, &work[1], &c__1, &ap[kc], &
			c__1);
#line 252 "ssptri.f"
		i__1 = k - 1;
#line 252 "ssptri.f"
		ap[kcnext + k - 1] -= sdot_(&i__1, &ap[kc], &c__1, &ap[kcnext]
			, &c__1);
#line 255 "ssptri.f"
		i__1 = k - 1;
#line 255 "ssptri.f"
		scopy_(&i__1, &ap[kcnext], &c__1, &work[1], &c__1);
#line 256 "ssptri.f"
		i__1 = k - 1;
#line 256 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kcnext], &c__1, (ftnlen)1);
#line 258 "ssptri.f"
		i__1 = k - 1;
#line 258 "ssptri.f"
		ap[kcnext + k] -= sdot_(&i__1, &work[1], &c__1, &ap[kcnext], &
			c__1);
#line 260 "ssptri.f"
	    }
#line 261 "ssptri.f"
	    kstep = 2;
#line 262 "ssptri.f"
	    kcnext = kcnext + k + 1;
#line 263 "ssptri.f"
	}

#line 265 "ssptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 266 "ssptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 271 "ssptri.f"
	    kpc = (kp - 1) * kp / 2 + 1;
#line 272 "ssptri.f"
	    i__1 = kp - 1;
#line 272 "ssptri.f"
	    sswap_(&i__1, &ap[kc], &c__1, &ap[kpc], &c__1);
#line 273 "ssptri.f"
	    kx = kpc + kp - 1;
#line 274 "ssptri.f"
	    i__1 = k - 1;
#line 274 "ssptri.f"
	    for (j = kp + 1; j <= i__1; ++j) {
#line 275 "ssptri.f"
		kx = kx + j - 1;
#line 276 "ssptri.f"
		temp = ap[kc + j - 1];
#line 277 "ssptri.f"
		ap[kc + j - 1] = ap[kx];
#line 278 "ssptri.f"
		ap[kx] = temp;
#line 279 "ssptri.f"
/* L40: */
#line 279 "ssptri.f"
	    }
#line 280 "ssptri.f"
	    temp = ap[kc + k - 1];
#line 281 "ssptri.f"
	    ap[kc + k - 1] = ap[kpc + kp - 1];
#line 282 "ssptri.f"
	    ap[kpc + kp - 1] = temp;
#line 283 "ssptri.f"
	    if (kstep == 2) {
#line 284 "ssptri.f"
		temp = ap[kc + k + k - 1];
#line 285 "ssptri.f"
		ap[kc + k + k - 1] = ap[kc + k + kp - 1];
#line 286 "ssptri.f"
		ap[kc + k + kp - 1] = temp;
#line 287 "ssptri.f"
	    }
#line 288 "ssptri.f"
	}

#line 290 "ssptri.f"
	k += kstep;
#line 291 "ssptri.f"
	kc = kcnext;
#line 292 "ssptri.f"
	goto L30;
#line 293 "ssptri.f"
L50:

#line 295 "ssptri.f"
	;
#line 295 "ssptri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 302 "ssptri.f"
	npp = *n * (*n + 1) / 2;
#line 303 "ssptri.f"
	k = *n;
#line 304 "ssptri.f"
	kc = npp;
#line 305 "ssptri.f"
L60:

/*        If K < 1, exit from loop. */

#line 309 "ssptri.f"
	if (k < 1) {
#line 309 "ssptri.f"
	    goto L80;
#line 309 "ssptri.f"
	}

#line 312 "ssptri.f"
	kcnext = kc - (*n - k + 2);
#line 313 "ssptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 319 "ssptri.f"
	    ap[kc] = 1. / ap[kc];

/*           Compute column K of the inverse. */

#line 323 "ssptri.f"
	    if (k < *n) {
#line 324 "ssptri.f"
		i__1 = *n - k;
#line 324 "ssptri.f"
		scopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 325 "ssptri.f"
		i__1 = *n - k;
#line 325 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[kc + *n - k + 1], &work[1], &
			c__1, &c_b13, &ap[kc + 1], &c__1, (ftnlen)1);
#line 327 "ssptri.f"
		i__1 = *n - k;
#line 327 "ssptri.f"
		ap[kc] -= sdot_(&i__1, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 328 "ssptri.f"
	    }
#line 329 "ssptri.f"
	    kstep = 1;
#line 330 "ssptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 336 "ssptri.f"
	    t = (d__1 = ap[kcnext + 1], abs(d__1));
#line 337 "ssptri.f"
	    ak = ap[kcnext] / t;
#line 338 "ssptri.f"
	    akp1 = ap[kc] / t;
#line 339 "ssptri.f"
	    akkp1 = ap[kcnext + 1] / t;
#line 340 "ssptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 341 "ssptri.f"
	    ap[kcnext] = akp1 / d__;
#line 342 "ssptri.f"
	    ap[kc] = ak / d__;
#line 343 "ssptri.f"
	    ap[kcnext + 1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 347 "ssptri.f"
	    if (k < *n) {
#line 348 "ssptri.f"
		i__1 = *n - k;
#line 348 "ssptri.f"
		scopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 349 "ssptri.f"
		i__1 = *n - k;
#line 349 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[kc + (*n - k + 1)], &work[1], 
			&c__1, &c_b13, &ap[kc + 1], &c__1, (ftnlen)1);
#line 351 "ssptri.f"
		i__1 = *n - k;
#line 351 "ssptri.f"
		ap[kc] -= sdot_(&i__1, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 352 "ssptri.f"
		i__1 = *n - k;
#line 352 "ssptri.f"
		ap[kcnext + 1] -= sdot_(&i__1, &ap[kc + 1], &c__1, &ap[kcnext 
			+ 2], &c__1);
#line 355 "ssptri.f"
		i__1 = *n - k;
#line 355 "ssptri.f"
		scopy_(&i__1, &ap[kcnext + 2], &c__1, &work[1], &c__1);
#line 356 "ssptri.f"
		i__1 = *n - k;
#line 356 "ssptri.f"
		sspmv_(uplo, &i__1, &c_b11, &ap[kc + (*n - k + 1)], &work[1], 
			&c__1, &c_b13, &ap[kcnext + 2], &c__1, (ftnlen)1);
#line 358 "ssptri.f"
		i__1 = *n - k;
#line 358 "ssptri.f"
		ap[kcnext] -= sdot_(&i__1, &work[1], &c__1, &ap[kcnext + 2], &
			c__1);
#line 360 "ssptri.f"
	    }
#line 361 "ssptri.f"
	    kstep = 2;
#line 362 "ssptri.f"
	    kcnext -= *n - k + 3;
#line 363 "ssptri.f"
	}

#line 365 "ssptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 366 "ssptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 371 "ssptri.f"
	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
#line 372 "ssptri.f"
	    if (kp < *n) {
#line 372 "ssptri.f"
		i__1 = *n - kp;
#line 372 "ssptri.f"
		sswap_(&i__1, &ap[kc + kp - k + 1], &c__1, &ap[kpc + 1], &
			c__1);
#line 372 "ssptri.f"
	    }
#line 374 "ssptri.f"
	    kx = kc + kp - k;
#line 375 "ssptri.f"
	    i__1 = kp - 1;
#line 375 "ssptri.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 376 "ssptri.f"
		kx = kx + *n - j + 1;
#line 377 "ssptri.f"
		temp = ap[kc + j - k];
#line 378 "ssptri.f"
		ap[kc + j - k] = ap[kx];
#line 379 "ssptri.f"
		ap[kx] = temp;
#line 380 "ssptri.f"
/* L70: */
#line 380 "ssptri.f"
	    }
#line 381 "ssptri.f"
	    temp = ap[kc];
#line 382 "ssptri.f"
	    ap[kc] = ap[kpc];
#line 383 "ssptri.f"
	    ap[kpc] = temp;
#line 384 "ssptri.f"
	    if (kstep == 2) {
#line 385 "ssptri.f"
		temp = ap[kc - *n + k - 1];
#line 386 "ssptri.f"
		ap[kc - *n + k - 1] = ap[kc - *n + kp - 1];
#line 387 "ssptri.f"
		ap[kc - *n + kp - 1] = temp;
#line 388 "ssptri.f"
	    }
#line 389 "ssptri.f"
	}

#line 391 "ssptri.f"
	k -= kstep;
#line 392 "ssptri.f"
	kc = kcnext;
#line 393 "ssptri.f"
	goto L60;
#line 394 "ssptri.f"
L80:
#line 395 "ssptri.f"
	;
#line 395 "ssptri.f"
    }

#line 397 "ssptri.f"
    return 0;

/*     End of SSPTRI */

} /* ssptri_ */

