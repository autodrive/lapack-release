#line 1 "dsytri_rook.f"
/* dsytri_rook.f -- translated by f2c (version 20100827).
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

#line 1 "dsytri_rook.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;

/* > \brief \b DSYTRI_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRI_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */

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
/* > DSYTRI_ROOK computes the inverse of a real symmetric */
/* > matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* > computed by DSYTRF_ROOK. */
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
/* >          used to obtain the factor U or L as computed by DSYTRF_ROOK. */
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
/* >          as determined by DSYTRF_ROOK. */
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

/* > \date April 2012 */

/* > \ingroup doubleSYcomputational */

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
/* Subroutine */ int dsytri_rook__(char *uplo, integer *n, doublereal *a, 
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

#line 172 "dsytri_rook.f"
    /* Parameter adjustments */
#line 172 "dsytri_rook.f"
    a_dim1 = *lda;
#line 172 "dsytri_rook.f"
    a_offset = 1 + a_dim1;
#line 172 "dsytri_rook.f"
    a -= a_offset;
#line 172 "dsytri_rook.f"
    --ipiv;
#line 172 "dsytri_rook.f"
    --work;
#line 172 "dsytri_rook.f"

#line 172 "dsytri_rook.f"
    /* Function Body */
#line 172 "dsytri_rook.f"
    *info = 0;
#line 173 "dsytri_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 174 "dsytri_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 175 "dsytri_rook.f"
	*info = -1;
#line 176 "dsytri_rook.f"
    } else if (*n < 0) {
#line 177 "dsytri_rook.f"
	*info = -2;
#line 178 "dsytri_rook.f"
    } else if (*lda < max(1,*n)) {
#line 179 "dsytri_rook.f"
	*info = -4;
#line 180 "dsytri_rook.f"
    }
#line 181 "dsytri_rook.f"
    if (*info != 0) {
#line 182 "dsytri_rook.f"
	i__1 = -(*info);
#line 182 "dsytri_rook.f"
	xerbla_("DSYTRI_ROOK", &i__1, (ftnlen)11);
#line 183 "dsytri_rook.f"
	return 0;
#line 184 "dsytri_rook.f"
    }

/*     Quick return if possible */

#line 188 "dsytri_rook.f"
    if (*n == 0) {
#line 188 "dsytri_rook.f"
	return 0;
#line 188 "dsytri_rook.f"
    }

/*     Check that the diagonal matrix D is nonsingular. */

#line 193 "dsytri_rook.f"
    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

#line 197 "dsytri_rook.f"
	for (*info = *n; *info >= 1; --(*info)) {
#line 198 "dsytri_rook.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 198 "dsytri_rook.f"
		return 0;
#line 198 "dsytri_rook.f"
	    }
#line 200 "dsytri_rook.f"
/* L10: */
#line 200 "dsytri_rook.f"
	}
#line 201 "dsytri_rook.f"
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

#line 205 "dsytri_rook.f"
	i__1 = *n;
#line 205 "dsytri_rook.f"
	for (*info = 1; *info <= i__1; ++(*info)) {
#line 206 "dsytri_rook.f"
	    if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
#line 206 "dsytri_rook.f"
		return 0;
#line 206 "dsytri_rook.f"
	    }
#line 208 "dsytri_rook.f"
/* L20: */
#line 208 "dsytri_rook.f"
	}
#line 209 "dsytri_rook.f"
    }
#line 210 "dsytri_rook.f"
    *info = 0;

#line 212 "dsytri_rook.f"
    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 219 "dsytri_rook.f"
	k = 1;
#line 220 "dsytri_rook.f"
L30:

/*        If K > N, exit from loop. */

#line 224 "dsytri_rook.f"
	if (k > *n) {
#line 224 "dsytri_rook.f"
	    goto L40;
#line 224 "dsytri_rook.f"
	}

#line 227 "dsytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 233 "dsytri_rook.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 237 "dsytri_rook.f"
	    if (k > 1) {
#line 238 "dsytri_rook.f"
		i__1 = k - 1;
#line 238 "dsytri_rook.f"
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 239 "dsytri_rook.f"
		i__1 = k - 1;
#line 239 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 241 "dsytri_rook.f"
		i__1 = k - 1;
#line 241 "dsytri_rook.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 243 "dsytri_rook.f"
	    }
#line 244 "dsytri_rook.f"
	    kstep = 1;
#line 245 "dsytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 251 "dsytri_rook.f"
	    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
#line 252 "dsytri_rook.f"
	    ak = a[k + k * a_dim1] / t;
#line 253 "dsytri_rook.f"
	    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
#line 254 "dsytri_rook.f"
	    akkp1 = a[k + (k + 1) * a_dim1] / t;
#line 255 "dsytri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 256 "dsytri_rook.f"
	    a[k + k * a_dim1] = akp1 / d__;
#line 257 "dsytri_rook.f"
	    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
#line 258 "dsytri_rook.f"
	    a[k + (k + 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K and K+1 of the inverse. */

#line 262 "dsytri_rook.f"
	    if (k > 1) {
#line 263 "dsytri_rook.f"
		i__1 = k - 1;
#line 263 "dsytri_rook.f"
		dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
#line 264 "dsytri_rook.f"
		i__1 = k - 1;
#line 264 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
#line 266 "dsytri_rook.f"
		i__1 = k - 1;
#line 266 "dsytri_rook.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * 
			a_dim1 + 1], &c__1);
#line 268 "dsytri_rook.f"
		i__1 = k - 1;
#line 268 "dsytri_rook.f"
		a[k + (k + 1) * a_dim1] -= ddot_(&i__1, &a[k * a_dim1 + 1], &
			c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
#line 270 "dsytri_rook.f"
		i__1 = k - 1;
#line 270 "dsytri_rook.f"
		dcopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &
			c__1);
#line 271 "dsytri_rook.f"
		i__1 = k - 1;
#line 271 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &
			c__1, &c_b13, &a[(k + 1) * a_dim1 + 1], &c__1, (
			ftnlen)1);
#line 273 "dsytri_rook.f"
		i__1 = k - 1;
#line 273 "dsytri_rook.f"
		a[k + 1 + (k + 1) * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &
			a[(k + 1) * a_dim1 + 1], &c__1);
#line 275 "dsytri_rook.f"
	    }
#line 276 "dsytri_rook.f"
	    kstep = 2;
#line 277 "dsytri_rook.f"
	}

#line 279 "dsytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the leading */
/*           submatrix A(1:k+1,1:k+1) */

#line 284 "dsytri_rook.f"
	    kp = ipiv[k];
#line 285 "dsytri_rook.f"
	    if (kp != k) {
#line 286 "dsytri_rook.f"
		if (kp > 1) {
#line 286 "dsytri_rook.f"
		    i__1 = kp - 1;
#line 286 "dsytri_rook.f"
		    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 286 "dsytri_rook.f"
		}
#line 288 "dsytri_rook.f"
		i__1 = k - kp - 1;
#line 288 "dsytri_rook.f"
		dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 289 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 290 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 291 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 292 "dsytri_rook.f"
	    }
#line 293 "dsytri_rook.f"
	} else {

/*           Interchange rows and columns K and K+1 with -IPIV(K) and */
/*           -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1) */

#line 298 "dsytri_rook.f"
	    kp = -ipiv[k];
#line 299 "dsytri_rook.f"
	    if (kp != k) {
#line 300 "dsytri_rook.f"
		if (kp > 1) {
#line 300 "dsytri_rook.f"
		    i__1 = kp - 1;
#line 300 "dsytri_rook.f"
		    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 300 "dsytri_rook.f"
		}
#line 302 "dsytri_rook.f"
		i__1 = k - kp - 1;
#line 302 "dsytri_rook.f"
		dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);

#line 304 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 305 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 306 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 307 "dsytri_rook.f"
		temp = a[k + (k + 1) * a_dim1];
#line 308 "dsytri_rook.f"
		a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
#line 309 "dsytri_rook.f"
		a[kp + (k + 1) * a_dim1] = temp;
#line 310 "dsytri_rook.f"
	    }

#line 312 "dsytri_rook.f"
	    ++k;
#line 313 "dsytri_rook.f"
	    kp = -ipiv[k];
#line 314 "dsytri_rook.f"
	    if (kp != k) {
#line 315 "dsytri_rook.f"
		if (kp > 1) {
#line 315 "dsytri_rook.f"
		    i__1 = kp - 1;
#line 315 "dsytri_rook.f"
		    dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 
			    1], &c__1);
#line 315 "dsytri_rook.f"
		}
#line 317 "dsytri_rook.f"
		i__1 = k - kp - 1;
#line 317 "dsytri_rook.f"
		dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1)
			 * a_dim1], lda);
#line 318 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 319 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 320 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 321 "dsytri_rook.f"
	    }
#line 322 "dsytri_rook.f"
	}

#line 324 "dsytri_rook.f"
	++k;
#line 325 "dsytri_rook.f"
	goto L30;
#line 326 "dsytri_rook.f"
L40:

#line 328 "dsytri_rook.f"
	;
#line 328 "dsytri_rook.f"
    } else {

/*        Compute inv(A) from the factorization A = L*D*L**T. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 335 "dsytri_rook.f"
	k = *n;
#line 336 "dsytri_rook.f"
L50:

/*        If K < 1, exit from loop. */

#line 340 "dsytri_rook.f"
	if (k < 1) {
#line 340 "dsytri_rook.f"
	    goto L60;
#line 340 "dsytri_rook.f"
	}

#line 343 "dsytri_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Invert the diagonal block. */

#line 349 "dsytri_rook.f"
	    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];

/*           Compute column K of the inverse. */

#line 353 "dsytri_rook.f"
	    if (k < *n) {
#line 354 "dsytri_rook.f"
		i__1 = *n - k;
#line 354 "dsytri_rook.f"
		dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 355 "dsytri_rook.f"
		i__1 = *n - k;
#line 355 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 357 "dsytri_rook.f"
		i__1 = *n - k;
#line 357 "dsytri_rook.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 359 "dsytri_rook.f"
	    }
#line 360 "dsytri_rook.f"
	    kstep = 1;
#line 361 "dsytri_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Invert the diagonal block. */

#line 367 "dsytri_rook.f"
	    t = (d__1 = a[k + (k - 1) * a_dim1], abs(d__1));
#line 368 "dsytri_rook.f"
	    ak = a[k - 1 + (k - 1) * a_dim1] / t;
#line 369 "dsytri_rook.f"
	    akp1 = a[k + k * a_dim1] / t;
#line 370 "dsytri_rook.f"
	    akkp1 = a[k + (k - 1) * a_dim1] / t;
#line 371 "dsytri_rook.f"
	    d__ = t * (ak * akp1 - 1.);
#line 372 "dsytri_rook.f"
	    a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
#line 373 "dsytri_rook.f"
	    a[k + k * a_dim1] = ak / d__;
#line 374 "dsytri_rook.f"
	    a[k + (k - 1) * a_dim1] = -akkp1 / d__;

/*           Compute columns K-1 and K of the inverse. */

#line 378 "dsytri_rook.f"
	    if (k < *n) {
#line 379 "dsytri_rook.f"
		i__1 = *n - k;
#line 379 "dsytri_rook.f"
		dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
#line 380 "dsytri_rook.f"
		i__1 = *n - k;
#line 380 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], &
			c__1, (ftnlen)1);
#line 382 "dsytri_rook.f"
		i__1 = *n - k;
#line 382 "dsytri_rook.f"
		a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + 
			k * a_dim1], &c__1);
#line 384 "dsytri_rook.f"
		i__1 = *n - k;
#line 384 "dsytri_rook.f"
		a[k + (k - 1) * a_dim1] -= ddot_(&i__1, &a[k + 1 + k * a_dim1]
			, &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 387 "dsytri_rook.f"
		i__1 = *n - k;
#line 387 "dsytri_rook.f"
		dcopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &
			c__1);
#line 388 "dsytri_rook.f"
		i__1 = *n - k;
#line 388 "dsytri_rook.f"
		dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda,
			 &work[1], &c__1, &c_b13, &a[k + 1 + (k - 1) * a_dim1]
			, &c__1, (ftnlen)1);
#line 390 "dsytri_rook.f"
		i__1 = *n - k;
#line 390 "dsytri_rook.f"
		a[k - 1 + (k - 1) * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &
			a[k + 1 + (k - 1) * a_dim1], &c__1);
#line 392 "dsytri_rook.f"
	    }
#line 393 "dsytri_rook.f"
	    kstep = 2;
#line 394 "dsytri_rook.f"
	}

#line 396 "dsytri_rook.f"
	if (kstep == 1) {

/*           Interchange rows and columns K and IPIV(K) in the trailing */
/*           submatrix A(k-1:n,k-1:n) */

#line 401 "dsytri_rook.f"
	    kp = ipiv[k];
#line 402 "dsytri_rook.f"
	    if (kp != k) {
#line 403 "dsytri_rook.f"
		if (kp < *n) {
#line 403 "dsytri_rook.f"
		    i__1 = *n - kp;
#line 403 "dsytri_rook.f"
		    dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 403 "dsytri_rook.f"
		}
#line 405 "dsytri_rook.f"
		i__1 = kp - k - 1;
#line 405 "dsytri_rook.f"
		dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 406 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 407 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 408 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 409 "dsytri_rook.f"
	    }
#line 410 "dsytri_rook.f"
	} else {

/*           Interchange rows and columns K and K-1 with -IPIV(K) and */
/*           -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */

#line 415 "dsytri_rook.f"
	    kp = -ipiv[k];
#line 416 "dsytri_rook.f"
	    if (kp != k) {
#line 417 "dsytri_rook.f"
		if (kp < *n) {
#line 417 "dsytri_rook.f"
		    i__1 = *n - kp;
#line 417 "dsytri_rook.f"
		    dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 417 "dsytri_rook.f"
		}
#line 419 "dsytri_rook.f"
		i__1 = kp - k - 1;
#line 419 "dsytri_rook.f"
		dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);

#line 421 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 422 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 423 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 424 "dsytri_rook.f"
		temp = a[k + (k - 1) * a_dim1];
#line 425 "dsytri_rook.f"
		a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
#line 426 "dsytri_rook.f"
		a[kp + (k - 1) * a_dim1] = temp;
#line 427 "dsytri_rook.f"
	    }

#line 429 "dsytri_rook.f"
	    --k;
#line 430 "dsytri_rook.f"
	    kp = -ipiv[k];
#line 431 "dsytri_rook.f"
	    if (kp != k) {
#line 432 "dsytri_rook.f"
		if (kp < *n) {
#line 432 "dsytri_rook.f"
		    i__1 = *n - kp;
#line 432 "dsytri_rook.f"
		    dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + 
			    kp * a_dim1], &c__1);
#line 432 "dsytri_rook.f"
		}
#line 434 "dsytri_rook.f"
		i__1 = kp - k - 1;
#line 434 "dsytri_rook.f"
		dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) *
			 a_dim1], lda);
#line 435 "dsytri_rook.f"
		temp = a[k + k * a_dim1];
#line 436 "dsytri_rook.f"
		a[k + k * a_dim1] = a[kp + kp * a_dim1];
#line 437 "dsytri_rook.f"
		a[kp + kp * a_dim1] = temp;
#line 438 "dsytri_rook.f"
	    }
#line 439 "dsytri_rook.f"
	}

#line 441 "dsytri_rook.f"
	--k;
#line 442 "dsytri_rook.f"
	goto L50;
#line 443 "dsytri_rook.f"
L60:
#line 444 "dsytri_rook.f"
	;
#line 444 "dsytri_rook.f"
    }

#line 446 "dsytri_rook.f"
    return 0;

/*     End of DSYTRI_ROOK */

} /* dsytri_rook__ */

