#line 1 "dsytrs_rook.f"
/* dsytrs_rook.f -- translated by f2c (version 20100827).
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

#line 1 "dsytrs_rook.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b DSYTRS_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRS_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > DSYTRS_ROOK solves a system of linear equations A*X = B with */
/* > a real symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSYTRF_ROOK. */
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
/* >          obtain the factor U or L as computed by DSYTRF_ROOK. */
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
/* Subroutine */ int dsytrs_rook__(char *uplo, integer *n, integer *nrhs, 
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

#line 176 "dsytrs_rook.f"
    /* Parameter adjustments */
#line 176 "dsytrs_rook.f"
    a_dim1 = *lda;
#line 176 "dsytrs_rook.f"
    a_offset = 1 + a_dim1;
#line 176 "dsytrs_rook.f"
    a -= a_offset;
#line 176 "dsytrs_rook.f"
    --ipiv;
#line 176 "dsytrs_rook.f"
    b_dim1 = *ldb;
#line 176 "dsytrs_rook.f"
    b_offset = 1 + b_dim1;
#line 176 "dsytrs_rook.f"
    b -= b_offset;
#line 176 "dsytrs_rook.f"

#line 176 "dsytrs_rook.f"
    /* Function Body */
#line 176 "dsytrs_rook.f"
    *info = 0;
#line 177 "dsytrs_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 178 "dsytrs_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 179 "dsytrs_rook.f"
	*info = -1;
#line 180 "dsytrs_rook.f"
    } else if (*n < 0) {
#line 181 "dsytrs_rook.f"
	*info = -2;
#line 182 "dsytrs_rook.f"
    } else if (*nrhs < 0) {
#line 183 "dsytrs_rook.f"
	*info = -3;
#line 184 "dsytrs_rook.f"
    } else if (*lda < max(1,*n)) {
#line 185 "dsytrs_rook.f"
	*info = -5;
#line 186 "dsytrs_rook.f"
    } else if (*ldb < max(1,*n)) {
#line 187 "dsytrs_rook.f"
	*info = -8;
#line 188 "dsytrs_rook.f"
    }
#line 189 "dsytrs_rook.f"
    if (*info != 0) {
#line 190 "dsytrs_rook.f"
	i__1 = -(*info);
#line 190 "dsytrs_rook.f"
	xerbla_("DSYTRS_ROOK", &i__1, (ftnlen)11);
#line 191 "dsytrs_rook.f"
	return 0;
#line 192 "dsytrs_rook.f"
    }

/*     Quick return if possible */

#line 196 "dsytrs_rook.f"
    if (*n == 0 || *nrhs == 0) {
#line 196 "dsytrs_rook.f"
	return 0;
#line 196 "dsytrs_rook.f"
    }

#line 199 "dsytrs_rook.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 208 "dsytrs_rook.f"
	k = *n;
#line 209 "dsytrs_rook.f"
L10:

/*        If K < 1, exit from loop. */

#line 213 "dsytrs_rook.f"
	if (k < 1) {
#line 213 "dsytrs_rook.f"
	    goto L30;
#line 213 "dsytrs_rook.f"
	}

#line 216 "dsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 222 "dsytrs_rook.f"
	    kp = ipiv[k];
#line 223 "dsytrs_rook.f"
	    if (kp != k) {
#line 223 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 223 "dsytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 229 "dsytrs_rook.f"
	    i__1 = k - 1;
#line 229 "dsytrs_rook.f"
	    dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 234 "dsytrs_rook.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 234 "dsytrs_rook.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 235 "dsytrs_rook.f"
	    --k;
#line 236 "dsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 242 "dsytrs_rook.f"
	    kp = -ipiv[k];
#line 243 "dsytrs_rook.f"
	    if (kp != k) {
#line 243 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 243 "dsytrs_rook.f"
	    }

#line 246 "dsytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 247 "dsytrs_rook.f"
	    if (kp != k - 1) {
#line 247 "dsytrs_rook.f"
		dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "dsytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 253 "dsytrs_rook.f"
	    if (k > 2) {
#line 254 "dsytrs_rook.f"
		i__1 = k - 2;
#line 254 "dsytrs_rook.f"
		dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + 
			b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 256 "dsytrs_rook.f"
		i__1 = k - 2;
#line 256 "dsytrs_rook.f"
		dger_(&i__1, nrhs, &c_b7, &a[(k - 1) * a_dim1 + 1], &c__1, &b[
			k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 258 "dsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 262 "dsytrs_rook.f"
	    akm1k = a[k - 1 + k * a_dim1];
#line 263 "dsytrs_rook.f"
	    akm1 = a[k - 1 + (k - 1) * a_dim1] / akm1k;
#line 264 "dsytrs_rook.f"
	    ak = a[k + k * a_dim1] / akm1k;
#line 265 "dsytrs_rook.f"
	    denom = akm1 * ak - 1.;
#line 266 "dsytrs_rook.f"
	    i__1 = *nrhs;
#line 266 "dsytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 267 "dsytrs_rook.f"
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
#line 268 "dsytrs_rook.f"
		bk = b[k + j * b_dim1] / akm1k;
#line 269 "dsytrs_rook.f"
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 270 "dsytrs_rook.f"
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 271 "dsytrs_rook.f"
/* L20: */
#line 271 "dsytrs_rook.f"
	    }
#line 272 "dsytrs_rook.f"
	    k += -2;
#line 273 "dsytrs_rook.f"
	}

#line 275 "dsytrs_rook.f"
	goto L10;
#line 276 "dsytrs_rook.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 283 "dsytrs_rook.f"
	k = 1;
#line 284 "dsytrs_rook.f"
L40:

/*        If K > N, exit from loop. */

#line 288 "dsytrs_rook.f"
	if (k > *n) {
#line 288 "dsytrs_rook.f"
	    goto L50;
#line 288 "dsytrs_rook.f"
	}

#line 291 "dsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 298 "dsytrs_rook.f"
	    if (k > 1) {
#line 298 "dsytrs_rook.f"
		i__1 = k - 1;
#line 298 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 298 "dsytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 304 "dsytrs_rook.f"
	    kp = ipiv[k];
#line 305 "dsytrs_rook.f"
	    if (kp != k) {
#line 305 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 305 "dsytrs_rook.f"
	    }
#line 307 "dsytrs_rook.f"
	    ++k;
#line 308 "dsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 315 "dsytrs_rook.f"
	    if (k > 1) {
#line 316 "dsytrs_rook.f"
		i__1 = k - 1;
#line 316 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 318 "dsytrs_rook.f"
		i__1 = k - 1;
#line 318 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[
			(k + 1) * a_dim1 + 1], &c__1, &c_b19, &b[k + 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 320 "dsytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1). */

#line 324 "dsytrs_rook.f"
	    kp = -ipiv[k];
#line 325 "dsytrs_rook.f"
	    if (kp != k) {
#line 325 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 325 "dsytrs_rook.f"
	    }

#line 328 "dsytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 329 "dsytrs_rook.f"
	    if (kp != k + 1) {
#line 329 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 329 "dsytrs_rook.f"
	    }

#line 332 "dsytrs_rook.f"
	    k += 2;
#line 333 "dsytrs_rook.f"
	}

#line 335 "dsytrs_rook.f"
	goto L40;
#line 336 "dsytrs_rook.f"
L50:

#line 338 "dsytrs_rook.f"
	;
#line 338 "dsytrs_rook.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 347 "dsytrs_rook.f"
	k = 1;
#line 348 "dsytrs_rook.f"
L60:

/*        If K > N, exit from loop. */

#line 352 "dsytrs_rook.f"
	if (k > *n) {
#line 352 "dsytrs_rook.f"
	    goto L80;
#line 352 "dsytrs_rook.f"
	}

#line 355 "dsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 361 "dsytrs_rook.f"
	    kp = ipiv[k];
#line 362 "dsytrs_rook.f"
	    if (kp != k) {
#line 362 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 362 "dsytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 368 "dsytrs_rook.f"
	    if (k < *n) {
#line 368 "dsytrs_rook.f"
		i__1 = *n - k;
#line 368 "dsytrs_rook.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 1 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 368 "dsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 374 "dsytrs_rook.f"
	    d__1 = 1. / a[k + k * a_dim1];
#line 374 "dsytrs_rook.f"
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
#line 375 "dsytrs_rook.f"
	    ++k;
#line 376 "dsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1) */

#line 382 "dsytrs_rook.f"
	    kp = -ipiv[k];
#line 383 "dsytrs_rook.f"
	    if (kp != k) {
#line 383 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 383 "dsytrs_rook.f"
	    }

#line 386 "dsytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 387 "dsytrs_rook.f"
	    if (kp != k + 1) {
#line 387 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 387 "dsytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 393 "dsytrs_rook.f"
	    if (k < *n - 1) {
#line 394 "dsytrs_rook.f"
		i__1 = *n - k - 1;
#line 394 "dsytrs_rook.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 2 + k * a_dim1], &c__1, &b[k 
			+ b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 396 "dsytrs_rook.f"
		i__1 = *n - k - 1;
#line 396 "dsytrs_rook.f"
		dger_(&i__1, nrhs, &c_b7, &a[k + 2 + (k + 1) * a_dim1], &c__1,
			 &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 398 "dsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 402 "dsytrs_rook.f"
	    akm1k = a[k + 1 + k * a_dim1];
#line 403 "dsytrs_rook.f"
	    akm1 = a[k + k * a_dim1] / akm1k;
#line 404 "dsytrs_rook.f"
	    ak = a[k + 1 + (k + 1) * a_dim1] / akm1k;
#line 405 "dsytrs_rook.f"
	    denom = akm1 * ak - 1.;
#line 406 "dsytrs_rook.f"
	    i__1 = *nrhs;
#line 406 "dsytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 407 "dsytrs_rook.f"
		bkm1 = b[k + j * b_dim1] / akm1k;
#line 408 "dsytrs_rook.f"
		bk = b[k + 1 + j * b_dim1] / akm1k;
#line 409 "dsytrs_rook.f"
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 410 "dsytrs_rook.f"
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 411 "dsytrs_rook.f"
/* L70: */
#line 411 "dsytrs_rook.f"
	    }
#line 412 "dsytrs_rook.f"
	    k += 2;
#line 413 "dsytrs_rook.f"
	}

#line 415 "dsytrs_rook.f"
	goto L60;
#line 416 "dsytrs_rook.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 423 "dsytrs_rook.f"
	k = *n;
#line 424 "dsytrs_rook.f"
L90:

/*        If K < 1, exit from loop. */

#line 428 "dsytrs_rook.f"
	if (k < 1) {
#line 428 "dsytrs_rook.f"
	    goto L100;
#line 428 "dsytrs_rook.f"
	}

#line 431 "dsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 438 "dsytrs_rook.f"
	    if (k < *n) {
#line 438 "dsytrs_rook.f"
		i__1 = *n - k;
#line 438 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 438 "dsytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 444 "dsytrs_rook.f"
	    kp = ipiv[k];
#line 445 "dsytrs_rook.f"
	    if (kp != k) {
#line 445 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 445 "dsytrs_rook.f"
	    }
#line 447 "dsytrs_rook.f"
	    --k;
#line 448 "dsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 455 "dsytrs_rook.f"
	    if (k < *n) {
#line 456 "dsytrs_rook.f"
		i__1 = *n - k;
#line 456 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 458 "dsytrs_rook.f"
		i__1 = *n - k;
#line 458 "dsytrs_rook.f"
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b19, &b[
			k - 1 + b_dim1], ldb, (ftnlen)9);
#line 461 "dsytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 465 "dsytrs_rook.f"
	    kp = -ipiv[k];
#line 466 "dsytrs_rook.f"
	    if (kp != k) {
#line 466 "dsytrs_rook.f"
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 466 "dsytrs_rook.f"
	    }

#line 469 "dsytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 470 "dsytrs_rook.f"
	    if (kp != k - 1) {
#line 470 "dsytrs_rook.f"
		dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 470 "dsytrs_rook.f"
	    }

#line 473 "dsytrs_rook.f"
	    k += -2;
#line 474 "dsytrs_rook.f"
	}

#line 476 "dsytrs_rook.f"
	goto L90;
#line 477 "dsytrs_rook.f"
L100:
#line 478 "dsytrs_rook.f"
	;
#line 478 "dsytrs_rook.f"
    }

#line 480 "dsytrs_rook.f"
    return 0;

/*     End of DSYTRS_ROOK */

} /* dsytrs_rook__ */

