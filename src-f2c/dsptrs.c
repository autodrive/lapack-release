#line 1 "dsptrs.f"
/* dsptrs.f -- translated by f2c (version 20100827).
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

#line 1 "dsptrs.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b DSPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A stored in packed format using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by DSPTRF. */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSPTRF. */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
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
/* Subroutine */ int dsptrs_(char *uplo, integer *n, integer *nrhs, 
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

#line 155 "dsptrs.f"
    /* Parameter adjustments */
#line 155 "dsptrs.f"
    --ap;
#line 155 "dsptrs.f"
    --ipiv;
#line 155 "dsptrs.f"
    b_dim1 = *ldb;
#line 155 "dsptrs.f"
    b_offset = 1 + b_dim1;
#line 155 "dsptrs.f"
    b -= b_offset;
#line 155 "dsptrs.f"

#line 155 "dsptrs.f"
    /* Function Body */
#line 155 "dsptrs.f"
    *info = 0;
#line 156 "dsptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "dsptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "dsptrs.f"
	*info = -1;
#line 159 "dsptrs.f"
    } else if (*n < 0) {
#line 160 "dsptrs.f"
	*info = -2;
#line 161 "dsptrs.f"
    } else if (*nrhs < 0) {
#line 162 "dsptrs.f"
	*info = -3;
#line 163 "dsptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 164 "dsptrs.f"
	*info = -7;
#line 165 "dsptrs.f"
    }
#line 166 "dsptrs.f"
    if (*info != 0) {
#line 167 "dsptrs.f"
	i__1 = -(*info);
#line 167 "dsptrs.f"
	xerbla_("DSPTRS", &i__1, (ftnlen)6);
#line 168 "dsptrs.f"
	return 0;
#line 169 "dsptrs.f"
    }

/*     Quick return if possible */

#line 173 "dsptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 173 "dsptrs.f"
	return 0;
#line 173 "dsptrs.f"
    }

#line 176 "dsptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 185 "dsptrs.f"
	k = *n;
#line 186 "dsptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 187 "dsptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 191 "dsptrs.f"
	if (k < 1) {
#line 191 "dsptrs.f"
	    goto L30;
#line 191 "dsptrs.f"
	}

#line 194 "dsptrs.f"
	kc -= k;
#line 195 "dsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 201 "dsptrs.f"
	    kp = ipiv[k];
#line 202 "dsptrs.f"
	    if (kp != k) {
#line 202 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 202 "dsptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 208 "dsptrs.f"
	    i__1 = k - 1;
#line 208 "dsptrs.f"
	    dger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 213 "dsptrs.f"
	    d__1 = 1. / ap[kc + k - 1];
#line 213 "dsptrs.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 214 "dsptrs.f"
	    --k;
#line 215 "dsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 221 "dsptrs.f"
	    kp = -ipiv[k];
#line 222 "dsptrs.f"
	    if (kp != k - 1) {
#line 222 "dsptrs.f"
		dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 222 "dsptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 228 "dsptrs.f"
	    i__1 = k - 2;
#line 228 "dsptrs.f"
	    dger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);
#line 230 "dsptrs.f"
	    i__1 = k - 2;
#line 230 "dsptrs.f"
	    dger_(&i__1, nrhs, &c_b7, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 235 "dsptrs.f"
	    akm1k = ap[kc + k - 2];
#line 236 "dsptrs.f"
	    akm1 = ap[kc - 1] / akm1k;
#line 237 "dsptrs.f"
	    ak = ap[kc + k - 1] / akm1k;
#line 238 "dsptrs.f"
	    denom = akm1 * ak - 1.;
#line 239 "dsptrs.f"
	    i__1 = *nrhs;
#line 239 "dsptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 240 "dsptrs.f"
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
#line 241 "dsptrs.f"
		bk = b[k + j * b_dim1] / akm1k;
#line 242 "dsptrs.f"
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 243 "dsptrs.f"
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 244 "dsptrs.f"
/* L20: */
#line 244 "dsptrs.f"
	    }
#line 245 "dsptrs.f"
	    kc = kc - k + 1;
#line 246 "dsptrs.f"
	    k += -2;
#line 247 "dsptrs.f"
	}

#line 249 "dsptrs.f"
	goto L10;
#line 250 "dsptrs.f"
L30:

/*        Next solve U**T*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 257 "dsptrs.f"
	k = 1;
#line 258 "dsptrs.f"
	kc = 1;
#line 259 "dsptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 263 "dsptrs.f"
	if (k > *n) {
#line 263 "dsptrs.f"
	    goto L50;
#line 263 "dsptrs.f"
	}

#line 266 "dsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 273 "dsptrs.f"
	    i__1 = k - 1;
#line 273 "dsptrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and IPIV(K). */

#line 278 "dsptrs.f"
	    kp = ipiv[k];
#line 279 "dsptrs.f"
	    if (kp != k) {
#line 279 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 279 "dsptrs.f"
	    }
#line 281 "dsptrs.f"
	    kc += k;
#line 282 "dsptrs.f"
	    ++k;
#line 283 "dsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 290 "dsptrs.f"
	    i__1 = k - 1;
#line 290 "dsptrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)9);
#line 292 "dsptrs.f"
	    i__1 = k - 1;
#line 292 "dsptrs.f"
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc 
		    + k], &c__1, &c_b19, &b[k + 1 + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 297 "dsptrs.f"
	    kp = -ipiv[k];
#line 298 "dsptrs.f"
	    if (kp != k) {
#line 298 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 298 "dsptrs.f"
	    }
#line 300 "dsptrs.f"
	    kc = kc + (k << 1) + 1;
#line 301 "dsptrs.f"
	    k += 2;
#line 302 "dsptrs.f"
	}

#line 304 "dsptrs.f"
	goto L40;
#line 305 "dsptrs.f"
L50:

#line 307 "dsptrs.f"
	;
#line 307 "dsptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 316 "dsptrs.f"
	k = 1;
#line 317 "dsptrs.f"
	kc = 1;
#line 318 "dsptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "dsptrs.f"
	if (k > *n) {
#line 322 "dsptrs.f"
	    goto L80;
#line 322 "dsptrs.f"
	}

#line 325 "dsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "dsptrs.f"
	    kp = ipiv[k];
#line 332 "dsptrs.f"
	    if (kp != k) {
#line 332 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "dsptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "dsptrs.f"
	    if (k < *n) {
#line 338 "dsptrs.f"
		i__1 = *n - k;
#line 338 "dsptrs.f"
		dger_(&i__1, nrhs, &c_b7, &ap[kc + 1], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "dsptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "dsptrs.f"
	    d__1 = 1. / ap[kc];
#line 344 "dsptrs.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 345 "dsptrs.f"
	    kc = kc + *n - k + 1;
#line 346 "dsptrs.f"
	    ++k;
#line 347 "dsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 353 "dsptrs.f"
	    kp = -ipiv[k];
#line 354 "dsptrs.f"
	    if (kp != k + 1) {
#line 354 "dsptrs.f"
		dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 354 "dsptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 360 "dsptrs.f"
	    if (k < *n - 1) {
#line 361 "dsptrs.f"
		i__1 = *n - k - 1;
#line 361 "dsptrs.f"
		dger_(&i__1, nrhs, &c_b7, &ap[kc + 2], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 2 + b_dim1], ldb);
#line 363 "dsptrs.f"
		i__1 = *n - k - 1;
#line 363 "dsptrs.f"
		dger_(&i__1, nrhs, &c_b7, &ap[kc + *n - k + 2], &c__1, &b[k + 
			1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 365 "dsptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 369 "dsptrs.f"
	    akm1k = ap[kc + 1];
#line 370 "dsptrs.f"
	    akm1 = ap[kc] / akm1k;
#line 371 "dsptrs.f"
	    ak = ap[kc + *n - k + 1] / akm1k;
#line 372 "dsptrs.f"
	    denom = akm1 * ak - 1.;
#line 373 "dsptrs.f"
	    i__1 = *nrhs;
#line 373 "dsptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 374 "dsptrs.f"
		bkm1 = b[k + j * b_dim1] / akm1k;
#line 375 "dsptrs.f"
		bk = b[k + 1 + j * b_dim1] / akm1k;
#line 376 "dsptrs.f"
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 377 "dsptrs.f"
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 378 "dsptrs.f"
/* L70: */
#line 378 "dsptrs.f"
	    }
#line 379 "dsptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 380 "dsptrs.f"
	    k += 2;
#line 381 "dsptrs.f"
	}

#line 383 "dsptrs.f"
	goto L60;
#line 384 "dsptrs.f"
L80:

/*        Next solve L**T*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 391 "dsptrs.f"
	k = *n;
#line 392 "dsptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 393 "dsptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 397 "dsptrs.f"
	if (k < 1) {
#line 397 "dsptrs.f"
	    goto L100;
#line 397 "dsptrs.f"
	}

#line 400 "dsptrs.f"
	kc -= *n - k + 1;
#line 401 "dsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 408 "dsptrs.f"
	    if (k < *n) {
#line 408 "dsptrs.f"
		i__1 = *n - k;
#line 408 "dsptrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, 
			(ftnlen)9);
#line 408 "dsptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 414 "dsptrs.f"
	    kp = ipiv[k];
#line 415 "dsptrs.f"
	    if (kp != k) {
#line 415 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 415 "dsptrs.f"
	    }
#line 417 "dsptrs.f"
	    --k;
#line 418 "dsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 425 "dsptrs.f"
	    if (k < *n) {
#line 426 "dsptrs.f"
		i__1 = *n - k;
#line 426 "dsptrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, 
			(ftnlen)9);
#line 428 "dsptrs.f"
		i__1 = *n - k;
#line 428 "dsptrs.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc - (*n - k)], &c__1, &c_b19, &b[k - 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 431 "dsptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 435 "dsptrs.f"
	    kp = -ipiv[k];
#line 436 "dsptrs.f"
	    if (kp != k) {
#line 436 "dsptrs.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 436 "dsptrs.f"
	    }
#line 438 "dsptrs.f"
	    kc -= *n - k + 2;
#line 439 "dsptrs.f"
	    k += -2;
#line 440 "dsptrs.f"
	}

#line 442 "dsptrs.f"
	goto L90;
#line 443 "dsptrs.f"
L100:
#line 444 "dsptrs.f"
	;
#line 444 "dsptrs.f"
    }

#line 446 "dsptrs.f"
    return 0;

/*     End of DSPTRS */

} /* dsptrs_ */

