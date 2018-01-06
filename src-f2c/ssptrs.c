#line 1 "ssptrs.f"
/* ssptrs.f -- translated by f2c (version 20100827).
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

#line 1 "ssptrs.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b SSPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A stored in packed format using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by SSPTRF. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by SSPTRF. */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
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
/* Subroutine */ int ssptrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j, k;
    static doublereal ak, bk;
    static integer kc, kp;
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

#line 155 "ssptrs.f"
    /* Parameter adjustments */
#line 155 "ssptrs.f"
    --ap;
#line 155 "ssptrs.f"
    --ipiv;
#line 155 "ssptrs.f"
    b_dim1 = *ldb;
#line 155 "ssptrs.f"
    b_offset = 1 + b_dim1;
#line 155 "ssptrs.f"
    b -= b_offset;
#line 155 "ssptrs.f"

#line 155 "ssptrs.f"
    /* Function Body */
#line 155 "ssptrs.f"
    *info = 0;
#line 156 "ssptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "ssptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "ssptrs.f"
	*info = -1;
#line 159 "ssptrs.f"
    } else if (*n < 0) {
#line 160 "ssptrs.f"
	*info = -2;
#line 161 "ssptrs.f"
    } else if (*nrhs < 0) {
#line 162 "ssptrs.f"
	*info = -3;
#line 163 "ssptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 164 "ssptrs.f"
	*info = -7;
#line 165 "ssptrs.f"
    }
#line 166 "ssptrs.f"
    if (*info != 0) {
#line 167 "ssptrs.f"
	i__1 = -(*info);
#line 167 "ssptrs.f"
	xerbla_("SSPTRS", &i__1, (ftnlen)6);
#line 168 "ssptrs.f"
	return 0;
#line 169 "ssptrs.f"
    }

/*     Quick return if possible */

#line 173 "ssptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 173 "ssptrs.f"
	return 0;
#line 173 "ssptrs.f"
    }

#line 176 "ssptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 185 "ssptrs.f"
	k = *n;
#line 186 "ssptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 187 "ssptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 191 "ssptrs.f"
	if (k < 1) {
#line 191 "ssptrs.f"
	    goto L30;
#line 191 "ssptrs.f"
	}

#line 194 "ssptrs.f"
	kc -= k;
#line 195 "ssptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 201 "ssptrs.f"
	    kp = ipiv[k];
#line 202 "ssptrs.f"
	    if (kp != k) {
#line 202 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 202 "ssptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 208 "ssptrs.f"
	    i__1 = k - 1;
#line 208 "ssptrs.f"
	    sger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 213 "ssptrs.f"
	    d__1 = 1. / ap[kc + k - 1];
#line 213 "ssptrs.f"
	    sscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 214 "ssptrs.f"
	    --k;
#line 215 "ssptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 221 "ssptrs.f"
	    kp = -ipiv[k];
#line 222 "ssptrs.f"
	    if (kp != k - 1) {
#line 222 "ssptrs.f"
		sswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 222 "ssptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 228 "ssptrs.f"
	    i__1 = k - 2;
#line 228 "ssptrs.f"
	    sger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);
#line 230 "ssptrs.f"
	    i__1 = k - 2;
#line 230 "ssptrs.f"
	    sger_(&i__1, nrhs, &c_b7, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 235 "ssptrs.f"
	    akm1k = ap[kc + k - 2];
#line 236 "ssptrs.f"
	    akm1 = ap[kc - 1] / akm1k;
#line 237 "ssptrs.f"
	    ak = ap[kc + k - 1] / akm1k;
#line 238 "ssptrs.f"
	    denom = akm1 * ak - 1.;
#line 239 "ssptrs.f"
	    i__1 = *nrhs;
#line 239 "ssptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 240 "ssptrs.f"
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
#line 241 "ssptrs.f"
		bk = b[k + j * b_dim1] / akm1k;
#line 242 "ssptrs.f"
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 243 "ssptrs.f"
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 244 "ssptrs.f"
/* L20: */
#line 244 "ssptrs.f"
	    }
#line 245 "ssptrs.f"
	    kc = kc - k + 1;
#line 246 "ssptrs.f"
	    k += -2;
#line 247 "ssptrs.f"
	}

#line 249 "ssptrs.f"
	goto L10;
#line 250 "ssptrs.f"
L30:

/*        Next solve U**T*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 257 "ssptrs.f"
	k = 1;
#line 258 "ssptrs.f"
	kc = 1;
#line 259 "ssptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 263 "ssptrs.f"
	if (k > *n) {
#line 263 "ssptrs.f"
	    goto L50;
#line 263 "ssptrs.f"
	}

#line 266 "ssptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 273 "ssptrs.f"
	    i__1 = k - 1;
#line 273 "ssptrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and IPIV(K). */

#line 278 "ssptrs.f"
	    kp = ipiv[k];
#line 279 "ssptrs.f"
	    if (kp != k) {
#line 279 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 279 "ssptrs.f"
	    }
#line 281 "ssptrs.f"
	    kc += k;
#line 282 "ssptrs.f"
	    ++k;
#line 283 "ssptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 290 "ssptrs.f"
	    i__1 = k - 1;
#line 290 "ssptrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)9);
#line 292 "ssptrs.f"
	    i__1 = k - 1;
#line 292 "ssptrs.f"
	    sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc 
		    + k], &c__1, &c_b19, &b[k + 1 + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 297 "ssptrs.f"
	    kp = -ipiv[k];
#line 298 "ssptrs.f"
	    if (kp != k) {
#line 298 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 298 "ssptrs.f"
	    }
#line 300 "ssptrs.f"
	    kc = kc + (k << 1) + 1;
#line 301 "ssptrs.f"
	    k += 2;
#line 302 "ssptrs.f"
	}

#line 304 "ssptrs.f"
	goto L40;
#line 305 "ssptrs.f"
L50:

#line 307 "ssptrs.f"
	;
#line 307 "ssptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 316 "ssptrs.f"
	k = 1;
#line 317 "ssptrs.f"
	kc = 1;
#line 318 "ssptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "ssptrs.f"
	if (k > *n) {
#line 322 "ssptrs.f"
	    goto L80;
#line 322 "ssptrs.f"
	}

#line 325 "ssptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "ssptrs.f"
	    kp = ipiv[k];
#line 332 "ssptrs.f"
	    if (kp != k) {
#line 332 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "ssptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "ssptrs.f"
	    if (k < *n) {
#line 338 "ssptrs.f"
		i__1 = *n - k;
#line 338 "ssptrs.f"
		sger_(&i__1, nrhs, &c_b7, &ap[kc + 1], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "ssptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "ssptrs.f"
	    d__1 = 1. / ap[kc];
#line 344 "ssptrs.f"
	    sscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 345 "ssptrs.f"
	    kc = kc + *n - k + 1;
#line 346 "ssptrs.f"
	    ++k;
#line 347 "ssptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 353 "ssptrs.f"
	    kp = -ipiv[k];
#line 354 "ssptrs.f"
	    if (kp != k + 1) {
#line 354 "ssptrs.f"
		sswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 354 "ssptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 360 "ssptrs.f"
	    if (k < *n - 1) {
#line 361 "ssptrs.f"
		i__1 = *n - k - 1;
#line 361 "ssptrs.f"
		sger_(&i__1, nrhs, &c_b7, &ap[kc + 2], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 2 + b_dim1], ldb);
#line 363 "ssptrs.f"
		i__1 = *n - k - 1;
#line 363 "ssptrs.f"
		sger_(&i__1, nrhs, &c_b7, &ap[kc + *n - k + 2], &c__1, &b[k + 
			1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 365 "ssptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 369 "ssptrs.f"
	    akm1k = ap[kc + 1];
#line 370 "ssptrs.f"
	    akm1 = ap[kc] / akm1k;
#line 371 "ssptrs.f"
	    ak = ap[kc + *n - k + 1] / akm1k;
#line 372 "ssptrs.f"
	    denom = akm1 * ak - 1.;
#line 373 "ssptrs.f"
	    i__1 = *nrhs;
#line 373 "ssptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 374 "ssptrs.f"
		bkm1 = b[k + j * b_dim1] / akm1k;
#line 375 "ssptrs.f"
		bk = b[k + 1 + j * b_dim1] / akm1k;
#line 376 "ssptrs.f"
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 377 "ssptrs.f"
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 378 "ssptrs.f"
/* L70: */
#line 378 "ssptrs.f"
	    }
#line 379 "ssptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 380 "ssptrs.f"
	    k += 2;
#line 381 "ssptrs.f"
	}

#line 383 "ssptrs.f"
	goto L60;
#line 384 "ssptrs.f"
L80:

/*        Next solve L**T*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 391 "ssptrs.f"
	k = *n;
#line 392 "ssptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 393 "ssptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 397 "ssptrs.f"
	if (k < 1) {
#line 397 "ssptrs.f"
	    goto L100;
#line 397 "ssptrs.f"
	}

#line 400 "ssptrs.f"
	kc -= *n - k + 1;
#line 401 "ssptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 408 "ssptrs.f"
	    if (k < *n) {
#line 408 "ssptrs.f"
		i__1 = *n - k;
#line 408 "ssptrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, 
			(ftnlen)9);
#line 408 "ssptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 414 "ssptrs.f"
	    kp = ipiv[k];
#line 415 "ssptrs.f"
	    if (kp != k) {
#line 415 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 415 "ssptrs.f"
	    }
#line 417 "ssptrs.f"
	    --k;
#line 418 "ssptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 425 "ssptrs.f"
	    if (k < *n) {
#line 426 "ssptrs.f"
		i__1 = *n - k;
#line 426 "ssptrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, 
			(ftnlen)9);
#line 428 "ssptrs.f"
		i__1 = *n - k;
#line 428 "ssptrs.f"
		sgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc - (*n - k)], &c__1, &c_b19, &b[k - 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 431 "ssptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 435 "ssptrs.f"
	    kp = -ipiv[k];
#line 436 "ssptrs.f"
	    if (kp != k) {
#line 436 "ssptrs.f"
		sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 436 "ssptrs.f"
	    }
#line 438 "ssptrs.f"
	    kc -= *n - k + 2;
#line 439 "ssptrs.f"
	    k += -2;
#line 440 "ssptrs.f"
	}

#line 442 "ssptrs.f"
	goto L90;
#line 443 "ssptrs.f"
L100:
#line 444 "ssptrs.f"
	;
#line 444 "ssptrs.f"
    }

#line 446 "ssptrs.f"
    return 0;

/*     End of SSPTRS */

} /* ssptrs_ */

