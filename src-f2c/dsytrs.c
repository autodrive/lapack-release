#line 1 "dsytrs.f"
/* dsytrs.f -- translated by f2c (version 20100827).
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

#line 1 "dsytrs.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b DSYTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSYTRF. */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSYTRF. */
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
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
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

/* > \date November 2011 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsytrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j, k;
    static doublereal ak, bk;
    static integer kp;
    static doublereal akm1, bkm1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akm1k;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 160 "dsytrs.f"
    /* Parameter adjustments */
#line 160 "dsytrs.f"
    a_dim1 = *lda;
#line 160 "dsytrs.f"
    a_offset = 1 + a_dim1;
#line 160 "dsytrs.f"
    a -= a_offset;
#line 160 "dsytrs.f"
    --ipiv;
#line 160 "dsytrs.f"
    b_dim1 = *ldb;
#line 160 "dsytrs.f"
    b_offset = 1 + b_dim1;
#line 160 "dsytrs.f"
    b -= b_offset;
#line 160 "dsytrs.f"

#line 160 "dsytrs.f"
    /* Function Body */
#line 160 "dsytrs.f"
    *info = 0;
#line 161 "dsytrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "dsytrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "dsytrs.f"
	*info = -1;
#line 164 "dsytrs.f"
    } else if (*n < 0) {
#line 165 "dsytrs.f"
	*info = -2;
#line 166 "dsytrs.f"
    } else if (*nrhs < 0) {
#line 167 "dsytrs.f"
	*info = -3;
#line 168 "dsytrs.f"
    } else if (*lda < max(1,*n)) {
#line 169 "dsytrs.f"
	*info = -5;
#line 170 "dsytrs.f"
    } else if (*ldb < max(1,*n)) {
#line 171 "dsytrs.f"
	*info = -8;
#line 172 "dsytrs.f"
    }
#line 173 "dsytrs.f"
    if (*info != 0) {
#line 174 "dsytrs.f"
	i__1 = -(*info);
#line 174 "dsytrs.f"
	xerbla_("DSYTRS", &i__1, (ftnlen)6);
#line 175 "dsytrs.f"
	return 0;
#line 176 "dsytrs.f"
    }

/*     Quick return if possible */

#line 180 "dsytrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 180 "dsytrs.f"
	return 0;
#line 180 "dsytrs.f"
    }

#line 183 "dsytrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 192 "dsytrs.f"
	k = *n;
#line 193 "dsytrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 197 "dsytrs.f"
	if (k < 1) {
#line 197 "dsytrs.f"
	    goto L30;
#line 197 "dsytrs.f"
	}

#line 200 "dsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 206 "dsytrs.f"
	    kp = ipiv[k];
#line 207 "dsytrs.f"
	    if (kp != k) {
#line 207 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 207 "dsytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 213 "dsytrs.f"
	    i__1 = k - 1;
#line 213 "dsytrs.f"
	    dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 218 "dsytrs.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 218 "dsytrs.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 219 "dsytrs.f"
	    --k;
#line 220 "dsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 226 "dsytrs.f"
	    kp = -ipiv[k];
#line 227 "dsytrs.f"
	    if (kp != k - 1) {
#line 227 "dsytrs.f"
		dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 227 "dsytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 233 "dsytrs.f"
	    i__1 = k - 2;
#line 233 "dsytrs.f"
	    dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 235 "dsytrs.f"
	    i__1 = k - 2;
#line 235 "dsytrs.f"
	    dger_(&i__1, nrhs, &c_b7, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 
		    1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 240 "dsytrs.f"
	    akm1k = a[k - 1 + k * a_dim1];
#line 241 "dsytrs.f"
	    akm1 = a[k - 1 + (k - 1) * a_dim1] / akm1k;
#line 242 "dsytrs.f"
	    ak = a[k + k * a_dim1] / akm1k;
#line 243 "dsytrs.f"
	    denom = akm1 * ak - 1.;
#line 244 "dsytrs.f"
	    i__1 = *nrhs;
#line 244 "dsytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "dsytrs.f"
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
#line 246 "dsytrs.f"
		bk = b[k + j * b_dim1] / akm1k;
#line 247 "dsytrs.f"
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 248 "dsytrs.f"
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 249 "dsytrs.f"
/* L20: */
#line 249 "dsytrs.f"
	    }
#line 250 "dsytrs.f"
	    k += -2;
#line 251 "dsytrs.f"
	}

#line 253 "dsytrs.f"
	goto L10;
#line 254 "dsytrs.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 261 "dsytrs.f"
	k = 1;
#line 262 "dsytrs.f"
L40:

/*        If K > N, exit from loop. */

#line 266 "dsytrs.f"
	if (k > *n) {
#line 266 "dsytrs.f"
	    goto L50;
#line 266 "dsytrs.f"
	}

#line 269 "dsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 276 "dsytrs.f"
	    i__1 = k - 1;
#line 276 "dsytrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)
		    9);

/*           Interchange rows K and IPIV(K). */

#line 281 "dsytrs.f"
	    kp = ipiv[k];
#line 282 "dsytrs.f"
	    if (kp != k) {
#line 282 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "dsytrs.f"
	    }
#line 284 "dsytrs.f"
	    ++k;
#line 285 "dsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 292 "dsytrs.f"
	    i__1 = k - 1;
#line 292 "dsytrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)
		    9);
#line 294 "dsytrs.f"
	    i__1 = k - 1;
#line 294 "dsytrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[(k 
		    + 1) * a_dim1 + 1], &c__1, &c_b19, &b[k + 1 + b_dim1], 
		    ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 299 "dsytrs.f"
	    kp = -ipiv[k];
#line 300 "dsytrs.f"
	    if (kp != k) {
#line 300 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 300 "dsytrs.f"
	    }
#line 302 "dsytrs.f"
	    k += 2;
#line 303 "dsytrs.f"
	}

#line 305 "dsytrs.f"
	goto L40;
#line 306 "dsytrs.f"
L50:

#line 308 "dsytrs.f"
	;
#line 308 "dsytrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 317 "dsytrs.f"
	k = 1;
#line 318 "dsytrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "dsytrs.f"
	if (k > *n) {
#line 322 "dsytrs.f"
	    goto L80;
#line 322 "dsytrs.f"
	}

#line 325 "dsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "dsytrs.f"
	    kp = ipiv[k];
#line 332 "dsytrs.f"
	    if (kp != k) {
#line 332 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "dsytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "dsytrs.f"
	    if (k < *n) {
#line 338 "dsytrs.f"
		i__1 = *n - k;
#line 338 "dsytrs.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 1 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "dsytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "dsytrs.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 344 "dsytrs.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 345 "dsytrs.f"
	    ++k;
#line 346 "dsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 352 "dsytrs.f"
	    kp = -ipiv[k];
#line 353 "dsytrs.f"
	    if (kp != k + 1) {
#line 353 "dsytrs.f"
		dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 353 "dsytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 359 "dsytrs.f"
	    if (k < *n - 1) {
#line 360 "dsytrs.f"
		i__1 = *n - k - 1;
#line 360 "dsytrs.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 2 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 362 "dsytrs.f"
		i__1 = *n - k - 1;
#line 362 "dsytrs.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 2 + (k + 1) * a_dim1], &c__1,
			 &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 364 "dsytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 368 "dsytrs.f"
	    akm1k = a[k + 1 + k * a_dim1];
#line 369 "dsytrs.f"
	    akm1 = a[k + k * a_dim1] / akm1k;
#line 370 "dsytrs.f"
	    ak = a[k + 1 + (k + 1) * a_dim1] / akm1k;
#line 371 "dsytrs.f"
	    denom = akm1 * ak - 1.;
#line 372 "dsytrs.f"
	    i__1 = *nrhs;
#line 372 "dsytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 373 "dsytrs.f"
		bkm1 = b[k + j * b_dim1] / akm1k;
#line 374 "dsytrs.f"
		bk = b[k + 1 + j * b_dim1] / akm1k;
#line 375 "dsytrs.f"
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 376 "dsytrs.f"
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 377 "dsytrs.f"
/* L70: */
#line 377 "dsytrs.f"
	    }
#line 378 "dsytrs.f"
	    k += 2;
#line 379 "dsytrs.f"
	}

#line 381 "dsytrs.f"
	goto L60;
#line 382 "dsytrs.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 389 "dsytrs.f"
	k = *n;
#line 390 "dsytrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 394 "dsytrs.f"
	if (k < 1) {
#line 394 "dsytrs.f"
	    goto L100;
#line 394 "dsytrs.f"
	}

#line 397 "dsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 404 "dsytrs.f"
	    if (k < *n) {
#line 404 "dsytrs.f"
		i__1 = *n - k;
#line 404 "dsytrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 404 "dsytrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 410 "dsytrs.f"
	    kp = ipiv[k];
#line 411 "dsytrs.f"
	    if (kp != k) {
#line 411 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 411 "dsytrs.f"
	    }
#line 413 "dsytrs.f"
	    --k;
#line 414 "dsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 421 "dsytrs.f"
	    if (k < *n) {
#line 422 "dsytrs.f"
		i__1 = *n - k;
#line 422 "dsytrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 424 "dsytrs.f"
		i__1 = *n - k;
#line 424 "dsytrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b19, &b[
			k - 1 + b_dim1], ldb, (ftnlen)9);
#line 427 "dsytrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 431 "dsytrs.f"
	    kp = -ipiv[k];
#line 432 "dsytrs.f"
	    if (kp != k) {
#line 432 "dsytrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 432 "dsytrs.f"
	    }
#line 434 "dsytrs.f"
	    k += -2;
#line 435 "dsytrs.f"
	}

#line 437 "dsytrs.f"
	goto L90;
#line 438 "dsytrs.f"
L100:
#line 439 "dsytrs.f"
	;
#line 439 "dsytrs.f"
    }

#line 441 "dsytrs.f"
    return 0;

/*     End of DSYTRS */

} /* dsytrs_ */

