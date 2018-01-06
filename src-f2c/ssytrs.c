#line 1 "ssytrs.f"
/* ssytrs.f -- translated by f2c (version 20100827).
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

#line 1 "ssytrs.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b SSYTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by SSYTRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSYTRF. */
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
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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

/* > \date December 2016 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssytrs_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

#line 160 "ssytrs.f"
    /* Parameter adjustments */
#line 160 "ssytrs.f"
    a_dim1 = *lda;
#line 160 "ssytrs.f"
    a_offset = 1 + a_dim1;
#line 160 "ssytrs.f"
    a -= a_offset;
#line 160 "ssytrs.f"
    --ipiv;
#line 160 "ssytrs.f"
    b_dim1 = *ldb;
#line 160 "ssytrs.f"
    b_offset = 1 + b_dim1;
#line 160 "ssytrs.f"
    b -= b_offset;
#line 160 "ssytrs.f"

#line 160 "ssytrs.f"
    /* Function Body */
#line 160 "ssytrs.f"
    *info = 0;
#line 161 "ssytrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "ssytrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "ssytrs.f"
	*info = -1;
#line 164 "ssytrs.f"
    } else if (*n < 0) {
#line 165 "ssytrs.f"
	*info = -2;
#line 166 "ssytrs.f"
    } else if (*nrhs < 0) {
#line 167 "ssytrs.f"
	*info = -3;
#line 168 "ssytrs.f"
    } else if (*lda < max(1,*n)) {
#line 169 "ssytrs.f"
	*info = -5;
#line 170 "ssytrs.f"
    } else if (*ldb < max(1,*n)) {
#line 171 "ssytrs.f"
	*info = -8;
#line 172 "ssytrs.f"
    }
#line 173 "ssytrs.f"
    if (*info != 0) {
#line 174 "ssytrs.f"
	i__1 = -(*info);
#line 174 "ssytrs.f"
	xerbla_("SSYTRS", &i__1, (ftnlen)6);
#line 175 "ssytrs.f"
	return 0;
#line 176 "ssytrs.f"
    }

/*     Quick return if possible */

#line 180 "ssytrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 180 "ssytrs.f"
	return 0;
#line 180 "ssytrs.f"
    }

#line 183 "ssytrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 192 "ssytrs.f"
	k = *n;
#line 193 "ssytrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 197 "ssytrs.f"
	if (k < 1) {
#line 197 "ssytrs.f"
	    goto L30;
#line 197 "ssytrs.f"
	}

#line 200 "ssytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 206 "ssytrs.f"
	    kp = ipiv[k];
#line 207 "ssytrs.f"
	    if (kp != k) {
#line 207 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 207 "ssytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 213 "ssytrs.f"
	    i__1 = k - 1;
#line 213 "ssytrs.f"
	    sger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 218 "ssytrs.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 218 "ssytrs.f"
	    sscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 219 "ssytrs.f"
	    --k;
#line 220 "ssytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 226 "ssytrs.f"
	    kp = -ipiv[k];
#line 227 "ssytrs.f"
	    if (kp != k - 1) {
#line 227 "ssytrs.f"
		sswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 227 "ssytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 233 "ssytrs.f"
	    i__1 = k - 2;
#line 233 "ssytrs.f"
	    sger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 235 "ssytrs.f"
	    i__1 = k - 2;
#line 235 "ssytrs.f"
	    sger_(&i__1, nrhs, &c_b7, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 
		    1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 240 "ssytrs.f"
	    akm1k = a[k - 1 + k * a_dim1];
#line 241 "ssytrs.f"
	    akm1 = a[k - 1 + (k - 1) * a_dim1] / akm1k;
#line 242 "ssytrs.f"
	    ak = a[k + k * a_dim1] / akm1k;
#line 243 "ssytrs.f"
	    denom = akm1 * ak - 1.;
#line 244 "ssytrs.f"
	    i__1 = *nrhs;
#line 244 "ssytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "ssytrs.f"
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
#line 246 "ssytrs.f"
		bk = b[k + j * b_dim1] / akm1k;
#line 247 "ssytrs.f"
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 248 "ssytrs.f"
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 249 "ssytrs.f"
/* L20: */
#line 249 "ssytrs.f"
	    }
#line 250 "ssytrs.f"
	    k += -2;
#line 251 "ssytrs.f"
	}

#line 253 "ssytrs.f"
	goto L10;
#line 254 "ssytrs.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 261 "ssytrs.f"
	k = 1;
#line 262 "ssytrs.f"
L40:

/*        If K > N, exit from loop. */

#line 266 "ssytrs.f"
	if (k > *n) {
#line 266 "ssytrs.f"
	    goto L50;
#line 266 "ssytrs.f"
	}

#line 269 "ssytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 276 "ssytrs.f"
	    i__1 = k - 1;
#line 276 "ssytrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)
		    9);

/*           Interchange rows K and IPIV(K). */

#line 281 "ssytrs.f"
	    kp = ipiv[k];
#line 282 "ssytrs.f"
	    if (kp != k) {
#line 282 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "ssytrs.f"
	    }
#line 284 "ssytrs.f"
	    ++k;
#line 285 "ssytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 292 "ssytrs.f"
	    i__1 = k - 1;
#line 292 "ssytrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)
		    9);
#line 294 "ssytrs.f"
	    i__1 = k - 1;
#line 294 "ssytrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[(k 
		    + 1) * a_dim1 + 1], &c__1, &c_b19, &b[k + 1 + b_dim1], 
		    ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 299 "ssytrs.f"
	    kp = -ipiv[k];
#line 300 "ssytrs.f"
	    if (kp != k) {
#line 300 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 300 "ssytrs.f"
	    }
#line 302 "ssytrs.f"
	    k += 2;
#line 303 "ssytrs.f"
	}

#line 305 "ssytrs.f"
	goto L40;
#line 306 "ssytrs.f"
L50:

#line 308 "ssytrs.f"
	;
#line 308 "ssytrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 317 "ssytrs.f"
	k = 1;
#line 318 "ssytrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "ssytrs.f"
	if (k > *n) {
#line 322 "ssytrs.f"
	    goto L80;
#line 322 "ssytrs.f"
	}

#line 325 "ssytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "ssytrs.f"
	    kp = ipiv[k];
#line 332 "ssytrs.f"
	    if (kp != k) {
#line 332 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "ssytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "ssytrs.f"
	    if (k < *n) {
#line 338 "ssytrs.f"
		i__1 = *n - k;
#line 338 "ssytrs.f"
		sger_(&i__1, nrhs, &c_b7, &a[k + 1 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "ssytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "ssytrs.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 344 "ssytrs.f"
	    sscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 345 "ssytrs.f"
	    ++k;
#line 346 "ssytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 352 "ssytrs.f"
	    kp = -ipiv[k];
#line 353 "ssytrs.f"
	    if (kp != k + 1) {
#line 353 "ssytrs.f"
		sswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 353 "ssytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 359 "ssytrs.f"
	    if (k < *n - 1) {
#line 360 "ssytrs.f"
		i__1 = *n - k - 1;
#line 360 "ssytrs.f"
		sger_(&i__1, nrhs, &c_b7, &a[k + 2 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 362 "ssytrs.f"
		i__1 = *n - k - 1;
#line 362 "ssytrs.f"
		sger_(&i__1, nrhs, &c_b7, &a[k + 2 + (k + 1) * a_dim1], &c__1,
			 &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 364 "ssytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 368 "ssytrs.f"
	    akm1k = a[k + 1 + k * a_dim1];
#line 369 "ssytrs.f"
	    akm1 = a[k + k * a_dim1] / akm1k;
#line 370 "ssytrs.f"
	    ak = a[k + 1 + (k + 1) * a_dim1] / akm1k;
#line 371 "ssytrs.f"
	    denom = akm1 * ak - 1.;
#line 372 "ssytrs.f"
	    i__1 = *nrhs;
#line 372 "ssytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 373 "ssytrs.f"
		bkm1 = b[k + j * b_dim1] / akm1k;
#line 374 "ssytrs.f"
		bk = b[k + 1 + j * b_dim1] / akm1k;
#line 375 "ssytrs.f"
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 376 "ssytrs.f"
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 377 "ssytrs.f"
/* L70: */
#line 377 "ssytrs.f"
	    }
#line 378 "ssytrs.f"
	    k += 2;
#line 379 "ssytrs.f"
	}

#line 381 "ssytrs.f"
	goto L60;
#line 382 "ssytrs.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 389 "ssytrs.f"
	k = *n;
#line 390 "ssytrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 394 "ssytrs.f"
	if (k < 1) {
#line 394 "ssytrs.f"
	    goto L100;
#line 394 "ssytrs.f"
	}

#line 397 "ssytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 404 "ssytrs.f"
	    if (k < *n) {
#line 404 "ssytrs.f"
		i__1 = *n - k;
#line 404 "ssytrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 404 "ssytrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 410 "ssytrs.f"
	    kp = ipiv[k];
#line 411 "ssytrs.f"
	    if (kp != k) {
#line 411 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 411 "ssytrs.f"
	    }
#line 413 "ssytrs.f"
	    --k;
#line 414 "ssytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 421 "ssytrs.f"
	    if (k < *n) {
#line 422 "ssytrs.f"
		i__1 = *n - k;
#line 422 "ssytrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 424 "ssytrs.f"
		i__1 = *n - k;
#line 424 "ssytrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b19, &b[
			k - 1 + b_dim1], ldb, (ftnlen)9);
#line 427 "ssytrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 431 "ssytrs.f"
	    kp = -ipiv[k];
#line 432 "ssytrs.f"
	    if (kp != k) {
#line 432 "ssytrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 432 "ssytrs.f"
	    }
#line 434 "ssytrs.f"
	    k += -2;
#line 435 "ssytrs.f"
	}

#line 437 "ssytrs.f"
	goto L90;
#line 438 "ssytrs.f"
L100:
#line 439 "ssytrs.f"
	;
#line 439 "ssytrs.f"
    }

#line 441 "ssytrs.f"
    return 0;

/*     End of SSYTRS */

} /* ssytrs_ */

