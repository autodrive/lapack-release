#line 1 "dsptri.f"
/* dsptri.f -- translated by f2c (version 20100827).
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

#line 1 "dsptri.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b DSPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   AP( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPTRI computes the inverse of a real symmetric indefinite matrix */
/* > A in packed storage using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSPTRF. */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          On entry, the block diagonal matrix D and the multipliers */
/* >          used to obtain the factor U or L as computed by DSPTRF, */
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
/* >          as determined by DSPTRF. */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsptri_(char *uplo, integer *n, doublereal *ap, integer *
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
    static doublereal akp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer kstep;
    extern /* Subroutine */ int dspmv_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer kcnext;


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

#line 152 "dsptri.f"
    /* Parameter adjustments */
#line 152 "dsptri.f"
    --work;
#line 152 "dsptri.f"
    --ipiv;
#line 152 "dsptri.f"
    --ap;
#line 152 "dsptri.f"

#line 152 "dsptri.f"
    /* Function Body */
#line 152 "dsptri.f"
    *info = 0;
#line 153 "dsptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 154 "dsptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 155 "dsptri.f"
	*info = -1;
#line 156 "dsptri.f"
    } else if (*n < 0) {
#line 157 "dsptri.f"
	*info = -2;
#line 158 "dsptri.f"
    }
#line 159 "dsptri.f"
    if (*info != 0) {
#line 160 "dsptri.f"
	i__1 = -(*info);
#line 160 "dsptri.f"
	xerbla_("DSPTRI", &i__1, (ftnlen)6);
#line 161 "dsptri.f"
	return 0;
#line 162 "dsptri.f"
    }

/*     Quick return if possible */

#line 166 "dsptri.f"
    if (*n == 0) {
#line 166 "dsptri.f"
	return 0;
#line 166 "dsptri.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 171 "dsptri.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 175 "dsptri.f"
	kp = *n * (*n + 1) / 2;
#line 176 "dsptri.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 177 "dsptri.f"
	    if (ipiv[*info] > 0 && ap[kp] == 0.) {
#line 177 "dsptri.f"
		return 0;
#line 177 "dsptri.f"
	    }
#line 179 "dsptri.f"
	    kp -= *info;
#line 180 "dsptri.f"
/* L10: */
#line 180 "dsptri.f"
	}
#line 181 "dsptri.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 185 "dsptri.f"
	kp = 1;
#line 186 "dsptri.f"
	i__1 = *n;
#line 186 "dsptri.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 187 "dsptri.f"
	    if (ipiv[*info] > 0 && ap[kp] == 0.) {
#line 187 "dsptri.f"
		return 0;
#line 187 "dsptri.f"
	    }
#line 189 "dsptri.f"
	    kp = kp + *n - *info + 1;
#line 190 "dsptri.f"
/* L20: */
#line 190 "dsptri.f"
	}
#line 191 "dsptri.f"
    }
#line 192 "dsptri.f"
    *info = 0;

#line 194 "dsptri.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 201 "dsptri.f"
	k = 1;
#line 202 "dsptri.f"
	kc = 1;
#line 203 "dsptri.f"
L30:

/*        If K > N, exit from loop. */

#line 207 "dsptri.f"
	if (k > *n) {
#line 207 "dsptri.f"
	    goto L50;
#line 207 "dsptri.f"
	}

#line 210 "dsptri.f"
	kcnext = kc + k;
#line 211 "dsptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 217 "dsptri.f"
	    ap[kc + k - 1] = 1. / ap[kc + k - 1];

/*           Compute column K of the inverse. */

#line 221 "dsptri.f"
	    if (k > 1) {
#line 222 "dsptri.f"
		i__1 = k - 1;
#line 222 "dsptri.f"
		dcopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 223 "dsptri.f"
		i__1 = k - 1;
#line 223 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kc], &c__1, (ftnlen)1);
#line 225 "dsptri.f"
		i__1 = k - 1;
#line 225 "dsptri.f"
		ap[kc + k - 1] -= ddot_(&i__1, &work[1], &c__1, &ap[kc], &
			c__1);
#line 227 "dsptri.f"
	    }
#line 228 "dsptri.f"
	    kstep = 1;
#line 229 "dsptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 235 "dsptri.f"
	    t = (d__1 = ap[kcnext + k - 1], abs(d__1));
#line 236 "dsptri.f"
	    ak = ap[kc + k - 1] / t;
#line 237 "dsptri.f"
	    akp1 = ap[kcnext + k] / t;
#line 238 "dsptri.f"
	    akkp1 = ap[kcnext + k - 1] / t;
#line 239 "dsptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 240 "dsptri.f"
	    ap[kc + k - 1] = akp1 / d__;
#line 241 "dsptri.f"
	    ap[kcnext + k] = ak / d__;
#line 242 "dsptri.f"
	    ap[kcnext + k - 1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 246 "dsptri.f"
	    if (k > 1) {
#line 247 "dsptri.f"
		i__1 = k - 1;
#line 247 "dsptri.f"
		dcopy_(&i__1, &ap[kc], &c__1, &work[1], &c__1);
#line 248 "dsptri.f"
		i__1 = k - 1;
#line 248 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kc], &c__1, (ftnlen)1);
#line 250 "dsptri.f"
		i__1 = k - 1;
#line 250 "dsptri.f"
		ap[kc + k - 1] -= ddot_(&i__1, &work[1], &c__1, &ap[kc], &
			c__1);
#line 252 "dsptri.f"
		i__1 = k - 1;
#line 252 "dsptri.f"
		ap[kcnext + k - 1] -= ddot_(&i__1, &ap[kc], &c__1, &ap[kcnext]
			, &c__1);
#line 255 "dsptri.f"
		i__1 = k - 1;
#line 255 "dsptri.f"
		dcopy_(&i__1, &ap[kcnext], &c__1, &work[1], &c__1);
#line 256 "dsptri.f"
		i__1 = k - 1;
#line 256 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[1], &work[1], &c__1, &c_b13, &
			ap[kcnext], &c__1, (ftnlen)1);
#line 258 "dsptri.f"
		i__1 = k - 1;
#line 258 "dsptri.f"
		ap[kcnext + k] -= ddot_(&i__1, &work[1], &c__1, &ap[kcnext], &
			c__1);
#line 260 "dsptri.f"
	    }
#line 261 "dsptri.f"
	    kstep = 2;
#line 262 "dsptri.f"
	    kcnext = kcnext + k + 1;
#line 263 "dsptri.f"
	}

#line 265 "dsptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 266 "dsptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 271 "dsptri.f"
	    kpc = (kp - 1) * kp / 2 + 1;
#line 272 "dsptri.f"
	    i__1 = kp - 1;
#line 272 "dsptri.f"
	    dswap_(&i__1, &ap[kc], &c__1, &ap[kpc], &c__1);
#line 273 "dsptri.f"
	    kx = kpc + kp - 1;
#line 274 "dsptri.f"
	    i__1 = k - 1;
#line 274 "dsptri.f"
	    for (j = kp + 1; j <= i__1; ++j) {
#line 275 "dsptri.f"
		kx = kx + j - 1;
#line 276 "dsptri.f"
		temp = ap[kc + j - 1];
#line 277 "dsptri.f"
		ap[kc + j - 1] = ap[kx];
#line 278 "dsptri.f"
		ap[kx] = temp;
#line 279 "dsptri.f"
/* L40: */
#line 279 "dsptri.f"
	    }
#line 280 "dsptri.f"
	    temp = ap[kc + k - 1];
#line 281 "dsptri.f"
	    ap[kc + k - 1] = ap[kpc + kp - 1];
#line 282 "dsptri.f"
	    ap[kpc + kp - 1] = temp;
#line 283 "dsptri.f"
	    if (kstep == 2) {
#line 284 "dsptri.f"
		temp = ap[kc + k + k - 1];
#line 285 "dsptri.f"
		ap[kc + k + k - 1] = ap[kc + k + kp - 1];
#line 286 "dsptri.f"
		ap[kc + k + kp - 1] = temp;
#line 287 "dsptri.f"
	    }
#line 288 "dsptri.f"
	}

#line 290 "dsptri.f"
	k += kstep;
#line 291 "dsptri.f"
	kc = kcnext;
#line 292 "dsptri.f"
	goto L30;
#line 293 "dsptri.f"
L50:

#line 295 "dsptri.f"
	;
#line 295 "dsptri.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 302 "dsptri.f"
	npp = *n * (*n + 1) / 2;
#line 303 "dsptri.f"
	k = *n;
#line 304 "dsptri.f"
	kc = npp;
#line 305 "dsptri.f"
L60:

/*        If K < 1, exit from loop. */

#line 309 "dsptri.f"
	if (k < 1) {
#line 309 "dsptri.f"
	    goto L80;
#line 309 "dsptri.f"
	}

#line 312 "dsptri.f"
	kcnext = kc - (*n - k + 2);
#line 313 "dsptri.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 319 "dsptri.f"
	    ap[kc] = 1. / ap[kc];

/*           Compute column K of the inverse. */

#line 323 "dsptri.f"
	    if (k < *n) {
#line 324 "dsptri.f"
		i__1 = *n - k;
#line 324 "dsptri.f"
		dcopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 325 "dsptri.f"
		i__1 = *n - k;
#line 325 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[kc + *n - k + 1], &work[1], &
			c__1, &c_b13, &ap[kc + 1], &c__1, (ftnlen)1);
#line 327 "dsptri.f"
		i__1 = *n - k;
#line 327 "dsptri.f"
		ap[kc] -= ddot_(&i__1, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 328 "dsptri.f"
	    }
#line 329 "dsptri.f"
	    kstep = 1;
#line 330 "dsptri.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 336 "dsptri.f"
	    t = (d__1 = ap[kcnext + 1], abs(d__1));
#line 337 "dsptri.f"
	    ak = ap[kcnext] / t;
#line 338 "dsptri.f"
	    akp1 = ap[kc] / t;
#line 339 "dsptri.f"
	    akkp1 = ap[kcnext + 1] / t;
#line 340 "dsptri.f"
	    d__ = t * (ak * akp1 - 1.);
#line 341 "dsptri.f"
	    ap[kcnext] = akp1 / d__;
#line 342 "dsptri.f"
	    ap[kc] = ak / d__;
#line 343 "dsptri.f"
	    ap[kcnext + 1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 347 "dsptri.f"
	    if (k < *n) {
#line 348 "dsptri.f"
		i__1 = *n - k;
#line 348 "dsptri.f"
		dcopy_(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
#line 349 "dsptri.f"
		i__1 = *n - k;
#line 349 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[kc + (*n - k + 1)], &work[1], 
			&c__1, &c_b13, &ap[kc + 1], &c__1, (ftnlen)1);
#line 351 "dsptri.f"
		i__1 = *n - k;
#line 351 "dsptri.f"
		ap[kc] -= ddot_(&i__1, &work[1], &c__1, &ap[kc + 1], &c__1);
#line 352 "dsptri.f"
		i__1 = *n - k;
#line 352 "dsptri.f"
		ap[kcnext + 1] -= ddot_(&i__1, &ap[kc + 1], &c__1, &ap[kcnext 
			+ 2], &c__1);
#line 355 "dsptri.f"
		i__1 = *n - k;
#line 355 "dsptri.f"
		dcopy_(&i__1, &ap[kcnext + 2], &c__1, &work[1], &c__1);
#line 356 "dsptri.f"
		i__1 = *n - k;
#line 356 "dsptri.f"
		dspmv_(uplo, &i__1, &c_b11, &ap[kc + (*n - k + 1)], &work[1], 
			&c__1, &c_b13, &ap[kcnext + 2], &c__1, (ftnlen)1);
#line 358 "dsptri.f"
		i__1 = *n - k;
#line 358 "dsptri.f"
		ap[kcnext] -= ddot_(&i__1, &work[1], &c__1, &ap[kcnext + 2], &
			c__1);
#line 360 "dsptri.f"
	    }
#line 361 "dsptri.f"
	    kstep = 2;
#line 362 "dsptri.f"
	    kcnext -= *n - k + 3;
#line 363 "dsptri.f"
	}

#line 365 "dsptri.f"
	kp = (i__1 = ipiv[k], abs(i__1));
#line 366 "dsptri.f"
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 371 "dsptri.f"
	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
#line 372 "dsptri.f"
	    if (kp < *n) {
#line 372 "dsptri.f"
		i__1 = *n - kp;
#line 372 "dsptri.f"
		dswap_(&i__1, &ap[kc + kp - k + 1], &c__1, &ap[kpc + 1], &
			c__1);
#line 372 "dsptri.f"
	    }
#line 374 "dsptri.f"
	    kx = kc + kp - k;
#line 375 "dsptri.f"
	    i__1 = kp - 1;
#line 375 "dsptri.f"
	    for (j = k + 1; j <= i__1; ++j) {
#line 376 "dsptri.f"
		kx = kx + *n - j + 1;
#line 377 "dsptri.f"
		temp = ap[kc + j - k];
#line 378 "dsptri.f"
		ap[kc + j - k] = ap[kx];
#line 379 "dsptri.f"
		ap[kx] = temp;
#line 380 "dsptri.f"
/* L70: */
#line 380 "dsptri.f"
	    }
#line 381 "dsptri.f"
	    temp = ap[kc];
#line 382 "dsptri.f"
	    ap[kc] = ap[kpc];
#line 383 "dsptri.f"
	    ap[kpc] = temp;
#line 384 "dsptri.f"
	    if (kstep == 2) {
#line 385 "dsptri.f"
		temp = ap[kc - *n + k - 1];
#line 386 "dsptri.f"
		ap[kc - *n + k - 1] = ap[kc - *n + kp - 1];
#line 387 "dsptri.f"
		ap[kc - *n + kp - 1] = temp;
#line 388 "dsptri.f"
	    }
#line 389 "dsptri.f"
	}

#line 391 "dsptri.f"
	k -= kstep;
#line 392 "dsptri.f"
	kc = kcnext;
#line 393 "dsptri.f"
	goto L60;
#line 394 "dsptri.f"
L80:
#line 395 "dsptri.f"
	;
#line 395 "dsptri.f"
    }

#line 397 "dsptri.f"
    return 0;

/*     End of DSPTRI */

} /* dsptri_ */

