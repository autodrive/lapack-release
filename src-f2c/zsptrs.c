#line 1 "zsptrs.f"
/* zsptrs.f -- translated by f2c (version 20100827).
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

#line 1 "zsptrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSPTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A stored in packed format using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by ZSPTRF. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZSPTRF. */
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
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex ak, bk;
    static integer kc, kp;
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

#line 155 "zsptrs.f"
    /* Parameter adjustments */
#line 155 "zsptrs.f"
    --ap;
#line 155 "zsptrs.f"
    --ipiv;
#line 155 "zsptrs.f"
    b_dim1 = *ldb;
#line 155 "zsptrs.f"
    b_offset = 1 + b_dim1;
#line 155 "zsptrs.f"
    b -= b_offset;
#line 155 "zsptrs.f"

#line 155 "zsptrs.f"
    /* Function Body */
#line 155 "zsptrs.f"
    *info = 0;
#line 156 "zsptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "zsptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "zsptrs.f"
	*info = -1;
#line 159 "zsptrs.f"
    } else if (*n < 0) {
#line 160 "zsptrs.f"
	*info = -2;
#line 161 "zsptrs.f"
    } else if (*nrhs < 0) {
#line 162 "zsptrs.f"
	*info = -3;
#line 163 "zsptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 164 "zsptrs.f"
	*info = -7;
#line 165 "zsptrs.f"
    }
#line 166 "zsptrs.f"
    if (*info != 0) {
#line 167 "zsptrs.f"
	i__1 = -(*info);
#line 167 "zsptrs.f"
	xerbla_("ZSPTRS", &i__1, (ftnlen)6);
#line 168 "zsptrs.f"
	return 0;
#line 169 "zsptrs.f"
    }

/*     Quick return if possible */

#line 173 "zsptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 173 "zsptrs.f"
	return 0;
#line 173 "zsptrs.f"
    }

#line 176 "zsptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 185 "zsptrs.f"
	k = *n;
#line 186 "zsptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 187 "zsptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 191 "zsptrs.f"
	if (k < 1) {
#line 191 "zsptrs.f"
	    goto L30;
#line 191 "zsptrs.f"
	}

#line 194 "zsptrs.f"
	kc -= k;
#line 195 "zsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 201 "zsptrs.f"
	    kp = ipiv[k];
#line 202 "zsptrs.f"
	    if (kp != k) {
#line 202 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 202 "zsptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 208 "zsptrs.f"
	    i__1 = k - 1;
#line 208 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 208 "zsptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 213 "zsptrs.f"
	    z_div(&z__1, &c_b1, &ap[kc + k - 1]);
#line 213 "zsptrs.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 214 "zsptrs.f"
	    --k;
#line 215 "zsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 221 "zsptrs.f"
	    kp = -ipiv[k];
#line 222 "zsptrs.f"
	    if (kp != k - 1) {
#line 222 "zsptrs.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 222 "zsptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 228 "zsptrs.f"
	    i__1 = k - 2;
#line 228 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 228 "zsptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);
#line 230 "zsptrs.f"
	    i__1 = k - 2;
#line 230 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "zsptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 235 "zsptrs.f"
	    i__1 = kc + k - 2;
#line 235 "zsptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 236 "zsptrs.f"
	    z_div(&z__1, &ap[kc - 1], &akm1k);
#line 236 "zsptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 237 "zsptrs.f"
	    z_div(&z__1, &ap[kc + k - 1], &akm1k);
#line 237 "zsptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 238 "zsptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 238 "zsptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 238 "zsptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 239 "zsptrs.f"
	    i__1 = *nrhs;
#line 239 "zsptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 240 "zsptrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 240 "zsptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 241 "zsptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 241 "zsptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 242 "zsptrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 242 "zsptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 242 "zsptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 242 "zsptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 242 "zsptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 243 "zsptrs.f"
		i__2 = k + j * b_dim1;
#line 243 "zsptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 243 "zsptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 243 "zsptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 243 "zsptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 244 "zsptrs.f"
/* L20: */
#line 244 "zsptrs.f"
	    }
#line 245 "zsptrs.f"
	    kc = kc - k + 1;
#line 246 "zsptrs.f"
	    k += -2;
#line 247 "zsptrs.f"
	}

#line 249 "zsptrs.f"
	goto L10;
#line 250 "zsptrs.f"
L30:

/*        Next solve U**T*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 257 "zsptrs.f"
	k = 1;
#line 258 "zsptrs.f"
	kc = 1;
#line 259 "zsptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 263 "zsptrs.f"
	if (k > *n) {
#line 263 "zsptrs.f"
	    goto L50;
#line 263 "zsptrs.f"
	}

#line 266 "zsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 273 "zsptrs.f"
	    i__1 = k - 1;
#line 273 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 273 "zsptrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and IPIV(K). */

#line 278 "zsptrs.f"
	    kp = ipiv[k];
#line 279 "zsptrs.f"
	    if (kp != k) {
#line 279 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 279 "zsptrs.f"
	    }
#line 281 "zsptrs.f"
	    kc += k;
#line 282 "zsptrs.f"
	    ++k;
#line 283 "zsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 290 "zsptrs.f"
	    i__1 = k - 1;
#line 290 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 290 "zsptrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9);
#line 292 "zsptrs.f"
	    i__1 = k - 1;
#line 292 "zsptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 292 "zsptrs.f"
	    zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc 
		    + k], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 297 "zsptrs.f"
	    kp = -ipiv[k];
#line 298 "zsptrs.f"
	    if (kp != k) {
#line 298 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 298 "zsptrs.f"
	    }
#line 300 "zsptrs.f"
	    kc = kc + (k << 1) + 1;
#line 301 "zsptrs.f"
	    k += 2;
#line 302 "zsptrs.f"
	}

#line 304 "zsptrs.f"
	goto L40;
#line 305 "zsptrs.f"
L50:

#line 307 "zsptrs.f"
	;
#line 307 "zsptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 316 "zsptrs.f"
	k = 1;
#line 317 "zsptrs.f"
	kc = 1;
#line 318 "zsptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "zsptrs.f"
	if (k > *n) {
#line 322 "zsptrs.f"
	    goto L80;
#line 322 "zsptrs.f"
	}

#line 325 "zsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "zsptrs.f"
	    kp = ipiv[k];
#line 332 "zsptrs.f"
	    if (kp != k) {
#line 332 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "zsptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "zsptrs.f"
	    if (k < *n) {
#line 338 "zsptrs.f"
		i__1 = *n - k;
#line 338 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 338 "zsptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + 1], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "zsptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "zsptrs.f"
	    z_div(&z__1, &c_b1, &ap[kc]);
#line 344 "zsptrs.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 345 "zsptrs.f"
	    kc = kc + *n - k + 1;
#line 346 "zsptrs.f"
	    ++k;
#line 347 "zsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 353 "zsptrs.f"
	    kp = -ipiv[k];
#line 354 "zsptrs.f"
	    if (kp != k + 1) {
#line 354 "zsptrs.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 354 "zsptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 360 "zsptrs.f"
	    if (k < *n - 1) {
#line 361 "zsptrs.f"
		i__1 = *n - k - 1;
#line 361 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 361 "zsptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + 2], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 2 + b_dim1], ldb);
#line 363 "zsptrs.f"
		i__1 = *n - k - 1;
#line 363 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 363 "zsptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + *n - k + 2], &c__1, &b[k 
			+ 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 365 "zsptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 369 "zsptrs.f"
	    i__1 = kc + 1;
#line 369 "zsptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 370 "zsptrs.f"
	    z_div(&z__1, &ap[kc], &akm1k);
#line 370 "zsptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 371 "zsptrs.f"
	    z_div(&z__1, &ap[kc + *n - k + 1], &akm1k);
#line 371 "zsptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 372 "zsptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 372 "zsptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 372 "zsptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 373 "zsptrs.f"
	    i__1 = *nrhs;
#line 373 "zsptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 374 "zsptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 374 "zsptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 375 "zsptrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 375 "zsptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 376 "zsptrs.f"
		i__2 = k + j * b_dim1;
#line 376 "zsptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 376 "zsptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 376 "zsptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 376 "zsptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 377 "zsptrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 377 "zsptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 377 "zsptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 377 "zsptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 377 "zsptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 378 "zsptrs.f"
/* L70: */
#line 378 "zsptrs.f"
	    }
#line 379 "zsptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 380 "zsptrs.f"
	    k += 2;
#line 381 "zsptrs.f"
	}

#line 383 "zsptrs.f"
	goto L60;
#line 384 "zsptrs.f"
L80:

/*        Next solve L**T*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 391 "zsptrs.f"
	k = *n;
#line 392 "zsptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 393 "zsptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 397 "zsptrs.f"
	if (k < 1) {
#line 397 "zsptrs.f"
	    goto L100;
#line 397 "zsptrs.f"
	}

#line 400 "zsptrs.f"
	kc -= *n - k + 1;
#line 401 "zsptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 408 "zsptrs.f"
	    if (k < *n) {
#line 408 "zsptrs.f"
		i__1 = *n - k;
#line 408 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 408 "zsptrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 408 "zsptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 414 "zsptrs.f"
	    kp = ipiv[k];
#line 415 "zsptrs.f"
	    if (kp != k) {
#line 415 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 415 "zsptrs.f"
	    }
#line 417 "zsptrs.f"
	    --k;
#line 418 "zsptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 425 "zsptrs.f"
	    if (k < *n) {
#line 426 "zsptrs.f"
		i__1 = *n - k;
#line 426 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 426 "zsptrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 428 "zsptrs.f"
		i__1 = *n - k;
#line 428 "zsptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 428 "zsptrs.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc - (*n - k)], &c__1, &c_b1, &b[k - 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 431 "zsptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 435 "zsptrs.f"
	    kp = -ipiv[k];
#line 436 "zsptrs.f"
	    if (kp != k) {
#line 436 "zsptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 436 "zsptrs.f"
	    }
#line 438 "zsptrs.f"
	    kc -= *n - k + 2;
#line 439 "zsptrs.f"
	    k += -2;
#line 440 "zsptrs.f"
	}

#line 442 "zsptrs.f"
	goto L90;
#line 443 "zsptrs.f"
L100:
#line 444 "zsptrs.f"
	;
#line 444 "zsptrs.f"
    }

#line 446 "zsptrs.f"
    return 0;

/*     End of ZSPTRS */

} /* zsptrs_ */

