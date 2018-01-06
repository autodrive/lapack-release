#line 1 "zsytrs.f"
/* zsytrs.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by ZSYTRF. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSYTRF. */
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
/* >          as determined by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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

/* > \ingroup complex16SYcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsytrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), xerbla_(char *, integer *,
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

#line 160 "zsytrs.f"
    /* Parameter adjustments */
#line 160 "zsytrs.f"
    a_dim1 = *lda;
#line 160 "zsytrs.f"
    a_offset = 1 + a_dim1;
#line 160 "zsytrs.f"
    a -= a_offset;
#line 160 "zsytrs.f"
    --ipiv;
#line 160 "zsytrs.f"
    b_dim1 = *ldb;
#line 160 "zsytrs.f"
    b_offset = 1 + b_dim1;
#line 160 "zsytrs.f"
    b -= b_offset;
#line 160 "zsytrs.f"

#line 160 "zsytrs.f"
    /* Function Body */
#line 160 "zsytrs.f"
    *info = 0;
#line 161 "zsytrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "zsytrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "zsytrs.f"
	*info = -1;
#line 164 "zsytrs.f"
    } else if (*n < 0) {
#line 165 "zsytrs.f"
	*info = -2;
#line 166 "zsytrs.f"
    } else if (*nrhs < 0) {
#line 167 "zsytrs.f"
	*info = -3;
#line 168 "zsytrs.f"
    } else if (*lda < max(1,*n)) {
#line 169 "zsytrs.f"
	*info = -5;
#line 170 "zsytrs.f"
    } else if (*ldb < max(1,*n)) {
#line 171 "zsytrs.f"
	*info = -8;
#line 172 "zsytrs.f"
    }
#line 173 "zsytrs.f"
    if (*info != 0) {
#line 174 "zsytrs.f"
	i__1 = -(*info);
#line 174 "zsytrs.f"
	xerbla_("ZSYTRS", &i__1, (ftnlen)6);
#line 175 "zsytrs.f"
	return 0;
#line 176 "zsytrs.f"
    }

/*     Quick return if possible */

#line 180 "zsytrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 180 "zsytrs.f"
	return 0;
#line 180 "zsytrs.f"
    }

#line 183 "zsytrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 192 "zsytrs.f"
	k = *n;
#line 193 "zsytrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 197 "zsytrs.f"
	if (k < 1) {
#line 197 "zsytrs.f"
	    goto L30;
#line 197 "zsytrs.f"
	}

#line 200 "zsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 206 "zsytrs.f"
	    kp = ipiv[k];
#line 207 "zsytrs.f"
	    if (kp != k) {
#line 207 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 207 "zsytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 213 "zsytrs.f"
	    i__1 = k - 1;
#line 213 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 213 "zsytrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 218 "zsytrs.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 218 "zsytrs.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 219 "zsytrs.f"
	    --k;
#line 220 "zsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 226 "zsytrs.f"
	    kp = -ipiv[k];
#line 227 "zsytrs.f"
	    if (kp != k - 1) {
#line 227 "zsytrs.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 227 "zsytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 233 "zsytrs.f"
	    i__1 = k - 2;
#line 233 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 233 "zsytrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 235 "zsytrs.f"
	    i__1 = k - 2;
#line 235 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 235 "zsytrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k 
		    - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 240 "zsytrs.f"
	    i__1 = k - 1 + k * a_dim1;
#line 240 "zsytrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 241 "zsytrs.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 241 "zsytrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 242 "zsytrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 242 "zsytrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 243 "zsytrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 243 "zsytrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 243 "zsytrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 244 "zsytrs.f"
	    i__1 = *nrhs;
#line 244 "zsytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "zsytrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 245 "zsytrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 246 "zsytrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 246 "zsytrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 247 "zsytrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 247 "zsytrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 247 "zsytrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 247 "zsytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 247 "zsytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 248 "zsytrs.f"
		i__2 = k + j * b_dim1;
#line 248 "zsytrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 248 "zsytrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 248 "zsytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 248 "zsytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 249 "zsytrs.f"
/* L20: */
#line 249 "zsytrs.f"
	    }
#line 250 "zsytrs.f"
	    k += -2;
#line 251 "zsytrs.f"
	}

#line 253 "zsytrs.f"
	goto L10;
#line 254 "zsytrs.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 261 "zsytrs.f"
	k = 1;
#line 262 "zsytrs.f"
L40:

/*        If K > N, exit from loop. */

#line 266 "zsytrs.f"
	if (k > *n) {
#line 266 "zsytrs.f"
	    goto L50;
#line 266 "zsytrs.f"
	}

#line 269 "zsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 276 "zsytrs.f"
	    i__1 = k - 1;
#line 276 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 276 "zsytrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9)
		    ;

/*           Interchange rows K and IPIV(K). */

#line 281 "zsytrs.f"
	    kp = ipiv[k];
#line 282 "zsytrs.f"
	    if (kp != k) {
#line 282 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "zsytrs.f"
	    }
#line 284 "zsytrs.f"
	    ++k;
#line 285 "zsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 292 "zsytrs.f"
	    i__1 = k - 1;
#line 292 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 292 "zsytrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9)
		    ;
#line 294 "zsytrs.f"
	    i__1 = k - 1;
#line 294 "zsytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 294 "zsytrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[(k 
		    + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb,
		     (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 299 "zsytrs.f"
	    kp = -ipiv[k];
#line 300 "zsytrs.f"
	    if (kp != k) {
#line 300 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 300 "zsytrs.f"
	    }
#line 302 "zsytrs.f"
	    k += 2;
#line 303 "zsytrs.f"
	}

#line 305 "zsytrs.f"
	goto L40;
#line 306 "zsytrs.f"
L50:

#line 308 "zsytrs.f"
	;
#line 308 "zsytrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 317 "zsytrs.f"
	k = 1;
#line 318 "zsytrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "zsytrs.f"
	if (k > *n) {
#line 322 "zsytrs.f"
	    goto L80;
#line 322 "zsytrs.f"
	}

#line 325 "zsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "zsytrs.f"
	    kp = ipiv[k];
#line 332 "zsytrs.f"
	    if (kp != k) {
#line 332 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "zsytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "zsytrs.f"
	    if (k < *n) {
#line 338 "zsytrs.f"
		i__1 = *n - k;
#line 338 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 338 "zsytrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "zsytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "zsytrs.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 344 "zsytrs.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 345 "zsytrs.f"
	    ++k;
#line 346 "zsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 352 "zsytrs.f"
	    kp = -ipiv[k];
#line 353 "zsytrs.f"
	    if (kp != k + 1) {
#line 353 "zsytrs.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 353 "zsytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 359 "zsytrs.f"
	    if (k < *n - 1) {
#line 360 "zsytrs.f"
		i__1 = *n - k - 1;
#line 360 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 360 "zsytrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 362 "zsytrs.f"
		i__1 = *n - k - 1;
#line 362 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 362 "zsytrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 364 "zsytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 368 "zsytrs.f"
	    i__1 = k + 1 + k * a_dim1;
#line 368 "zsytrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 369 "zsytrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 369 "zsytrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 370 "zsytrs.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 370 "zsytrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 371 "zsytrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 371 "zsytrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 371 "zsytrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 372 "zsytrs.f"
	    i__1 = *nrhs;
#line 372 "zsytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 373 "zsytrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 373 "zsytrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 374 "zsytrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 374 "zsytrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 375 "zsytrs.f"
		i__2 = k + j * b_dim1;
#line 375 "zsytrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 375 "zsytrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 375 "zsytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 375 "zsytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 376 "zsytrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 376 "zsytrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 376 "zsytrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 376 "zsytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 376 "zsytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 377 "zsytrs.f"
/* L70: */
#line 377 "zsytrs.f"
	    }
#line 378 "zsytrs.f"
	    k += 2;
#line 379 "zsytrs.f"
	}

#line 381 "zsytrs.f"
	goto L60;
#line 382 "zsytrs.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 389 "zsytrs.f"
	k = *n;
#line 390 "zsytrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 394 "zsytrs.f"
	if (k < 1) {
#line 394 "zsytrs.f"
	    goto L100;
#line 394 "zsytrs.f"
	}

#line 397 "zsytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 404 "zsytrs.f"
	    if (k < *n) {
#line 404 "zsytrs.f"
		i__1 = *n - k;
#line 404 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 404 "zsytrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 404 "zsytrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 410 "zsytrs.f"
	    kp = ipiv[k];
#line 411 "zsytrs.f"
	    if (kp != k) {
#line 411 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 411 "zsytrs.f"
	    }
#line 413 "zsytrs.f"
	    --k;
#line 414 "zsytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 421 "zsytrs.f"
	    if (k < *n) {
#line 422 "zsytrs.f"
		i__1 = *n - k;
#line 422 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 422 "zsytrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 424 "zsytrs.f"
		i__1 = *n - k;
#line 424 "zsytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 424 "zsytrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)9);
#line 427 "zsytrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 431 "zsytrs.f"
	    kp = -ipiv[k];
#line 432 "zsytrs.f"
	    if (kp != k) {
#line 432 "zsytrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 432 "zsytrs.f"
	    }
#line 434 "zsytrs.f"
	    k += -2;
#line 435 "zsytrs.f"
	}

#line 437 "zsytrs.f"
	goto L90;
#line 438 "zsytrs.f"
L100:
#line 439 "zsytrs.f"
	;
#line 439 "zsytrs.f"
    }

#line 441 "zsytrs.f"
    return 0;

/*     End of ZSYTRS */

} /* zsytrs_ */

