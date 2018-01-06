#line 1 "zhetrs.f"
/* zhetrs.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHETRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > ZHETRS solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHETRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**H; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**H. */
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
/* >          obtain the factor U or L as computed by ZHETRF. */
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
/* >          as determined by ZHETRF. */
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

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhetrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), xerbla_(char *, integer *,
	     ftnlen), zdscal_(integer *, doublereal *, doublecomplex *, 
	    integer *), zlacgv_(integer *, doublecomplex *, integer *);


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

#line 161 "zhetrs.f"
    /* Parameter adjustments */
#line 161 "zhetrs.f"
    a_dim1 = *lda;
#line 161 "zhetrs.f"
    a_offset = 1 + a_dim1;
#line 161 "zhetrs.f"
    a -= a_offset;
#line 161 "zhetrs.f"
    --ipiv;
#line 161 "zhetrs.f"
    b_dim1 = *ldb;
#line 161 "zhetrs.f"
    b_offset = 1 + b_dim1;
#line 161 "zhetrs.f"
    b -= b_offset;
#line 161 "zhetrs.f"

#line 161 "zhetrs.f"
    /* Function Body */
#line 161 "zhetrs.f"
    *info = 0;
#line 162 "zhetrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "zhetrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "zhetrs.f"
	*info = -1;
#line 165 "zhetrs.f"
    } else if (*n < 0) {
#line 166 "zhetrs.f"
	*info = -2;
#line 167 "zhetrs.f"
    } else if (*nrhs < 0) {
#line 168 "zhetrs.f"
	*info = -3;
#line 169 "zhetrs.f"
    } else if (*lda < max(1,*n)) {
#line 170 "zhetrs.f"
	*info = -5;
#line 171 "zhetrs.f"
    } else if (*ldb < max(1,*n)) {
#line 172 "zhetrs.f"
	*info = -8;
#line 173 "zhetrs.f"
    }
#line 174 "zhetrs.f"
    if (*info != 0) {
#line 175 "zhetrs.f"
	i__1 = -(*info);
#line 175 "zhetrs.f"
	xerbla_("ZHETRS", &i__1, (ftnlen)6);
#line 176 "zhetrs.f"
	return 0;
#line 177 "zhetrs.f"
    }

/*     Quick return if possible */

#line 181 "zhetrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 181 "zhetrs.f"
	return 0;
#line 181 "zhetrs.f"
    }

#line 184 "zhetrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 193 "zhetrs.f"
	k = *n;
#line 194 "zhetrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 198 "zhetrs.f"
	if (k < 1) {
#line 198 "zhetrs.f"
	    goto L30;
#line 198 "zhetrs.f"
	}

#line 201 "zhetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 207 "zhetrs.f"
	    kp = ipiv[k];
#line 208 "zhetrs.f"
	    if (kp != k) {
#line 208 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 208 "zhetrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 214 "zhetrs.f"
	    i__1 = k - 1;
#line 214 "zhetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 214 "zhetrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 219 "zhetrs.f"
	    i__1 = k + k * a_dim1;
#line 219 "zhetrs.f"
	    s = 1. / a[i__1].r;
#line 220 "zhetrs.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 221 "zhetrs.f"
	    --k;
#line 222 "zhetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 228 "zhetrs.f"
	    kp = -ipiv[k];
#line 229 "zhetrs.f"
	    if (kp != k - 1) {
#line 229 "zhetrs.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 229 "zhetrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 235 "zhetrs.f"
	    i__1 = k - 2;
#line 235 "zhetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 235 "zhetrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 237 "zhetrs.f"
	    i__1 = k - 2;
#line 237 "zhetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 237 "zhetrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k 
		    - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 242 "zhetrs.f"
	    i__1 = k - 1 + k * a_dim1;
#line 242 "zhetrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 243 "zhetrs.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 243 "zhetrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 244 "zhetrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 244 "zhetrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 244 "zhetrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 245 "zhetrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 245 "zhetrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 245 "zhetrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 246 "zhetrs.f"
	    i__1 = *nrhs;
#line 246 "zhetrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 247 "zhetrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 247 "zhetrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 248 "zhetrs.f"
		d_cnjg(&z__2, &akm1k);
#line 248 "zhetrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 248 "zhetrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 249 "zhetrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 249 "zhetrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 249 "zhetrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 249 "zhetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 249 "zhetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 250 "zhetrs.f"
		i__2 = k + j * b_dim1;
#line 250 "zhetrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 250 "zhetrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 250 "zhetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 250 "zhetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 251 "zhetrs.f"
/* L20: */
#line 251 "zhetrs.f"
	    }
#line 252 "zhetrs.f"
	    k += -2;
#line 253 "zhetrs.f"
	}

#line 255 "zhetrs.f"
	goto L10;
#line 256 "zhetrs.f"
L30:

/*        Next solve U**H *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 263 "zhetrs.f"
	k = 1;
#line 264 "zhetrs.f"
L40:

/*        If K > N, exit from loop. */

#line 268 "zhetrs.f"
	if (k > *n) {
#line 268 "zhetrs.f"
	    goto L50;
#line 268 "zhetrs.f"
	}

#line 271 "zhetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**H(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 278 "zhetrs.f"
	    if (k > 1) {
#line 279 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 280 "zhetrs.f"
		i__1 = k - 1;
#line 280 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 280 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 282 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 283 "zhetrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 287 "zhetrs.f"
	    kp = ipiv[k];
#line 288 "zhetrs.f"
	    if (kp != k) {
#line 288 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 288 "zhetrs.f"
	    }
#line 290 "zhetrs.f"
	    ++k;
#line 291 "zhetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 298 "zhetrs.f"
	    if (k > 1) {
#line 299 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 300 "zhetrs.f"
		i__1 = k - 1;
#line 300 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 300 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 302 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 304 "zhetrs.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 305 "zhetrs.f"
		i__1 = k - 1;
#line 305 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 305 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			1 + b_dim1], ldb, (ftnlen)19);
#line 307 "zhetrs.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 308 "zhetrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 312 "zhetrs.f"
	    kp = -ipiv[k];
#line 313 "zhetrs.f"
	    if (kp != k) {
#line 313 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 313 "zhetrs.f"
	    }
#line 315 "zhetrs.f"
	    k += 2;
#line 316 "zhetrs.f"
	}

#line 318 "zhetrs.f"
	goto L40;
#line 319 "zhetrs.f"
L50:

#line 321 "zhetrs.f"
	;
#line 321 "zhetrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 330 "zhetrs.f"
	k = 1;
#line 331 "zhetrs.f"
L60:

/*        If K > N, exit from loop. */

#line 335 "zhetrs.f"
	if (k > *n) {
#line 335 "zhetrs.f"
	    goto L80;
#line 335 "zhetrs.f"
	}

#line 338 "zhetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 344 "zhetrs.f"
	    kp = ipiv[k];
#line 345 "zhetrs.f"
	    if (kp != k) {
#line 345 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 345 "zhetrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 351 "zhetrs.f"
	    if (k < *n) {
#line 351 "zhetrs.f"
		i__1 = *n - k;
#line 351 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 351 "zhetrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 351 "zhetrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 357 "zhetrs.f"
	    i__1 = k + k * a_dim1;
#line 357 "zhetrs.f"
	    s = 1. / a[i__1].r;
#line 358 "zhetrs.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 359 "zhetrs.f"
	    ++k;
#line 360 "zhetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 366 "zhetrs.f"
	    kp = -ipiv[k];
#line 367 "zhetrs.f"
	    if (kp != k + 1) {
#line 367 "zhetrs.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 367 "zhetrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 373 "zhetrs.f"
	    if (k < *n - 1) {
#line 374 "zhetrs.f"
		i__1 = *n - k - 1;
#line 374 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 374 "zhetrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 376 "zhetrs.f"
		i__1 = *n - k - 1;
#line 376 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 376 "zhetrs.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 378 "zhetrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 382 "zhetrs.f"
	    i__1 = k + 1 + k * a_dim1;
#line 382 "zhetrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 383 "zhetrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 383 "zhetrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 383 "zhetrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 384 "zhetrs.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 384 "zhetrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 385 "zhetrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 385 "zhetrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 385 "zhetrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 386 "zhetrs.f"
	    i__1 = *nrhs;
#line 386 "zhetrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 387 "zhetrs.f"
		d_cnjg(&z__2, &akm1k);
#line 387 "zhetrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 387 "zhetrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 388 "zhetrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 388 "zhetrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 389 "zhetrs.f"
		i__2 = k + j * b_dim1;
#line 389 "zhetrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 389 "zhetrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 389 "zhetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 389 "zhetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 390 "zhetrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 390 "zhetrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 390 "zhetrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 390 "zhetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 390 "zhetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 391 "zhetrs.f"
/* L70: */
#line 391 "zhetrs.f"
	    }
#line 392 "zhetrs.f"
	    k += 2;
#line 393 "zhetrs.f"
	}

#line 395 "zhetrs.f"
	goto L60;
#line 396 "zhetrs.f"
L80:

/*        Next solve L**H *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 403 "zhetrs.f"
	k = *n;
#line 404 "zhetrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 408 "zhetrs.f"
	if (k < 1) {
#line 408 "zhetrs.f"
	    goto L100;
#line 408 "zhetrs.f"
	}

#line 411 "zhetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**H(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 418 "zhetrs.f"
	    if (k < *n) {
#line 419 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 420 "zhetrs.f"
		i__1 = *n - k;
#line 420 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 420 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 423 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 424 "zhetrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 428 "zhetrs.f"
	    kp = ipiv[k];
#line 429 "zhetrs.f"
	    if (kp != k) {
#line 429 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 429 "zhetrs.f"
	    }
#line 431 "zhetrs.f"
	    --k;
#line 432 "zhetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 439 "zhetrs.f"
	    if (k < *n) {
#line 440 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 441 "zhetrs.f"
		i__1 = *n - k;
#line 441 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 441 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 444 "zhetrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 446 "zhetrs.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 447 "zhetrs.f"
		i__1 = *n - k;
#line 447 "zhetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 447 "zhetrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &
			c_b1, &b[k - 1 + b_dim1], ldb, (ftnlen)19);
#line 450 "zhetrs.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 451 "zhetrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 455 "zhetrs.f"
	    kp = -ipiv[k];
#line 456 "zhetrs.f"
	    if (kp != k) {
#line 456 "zhetrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 456 "zhetrs.f"
	    }
#line 458 "zhetrs.f"
	    k += -2;
#line 459 "zhetrs.f"
	}

#line 461 "zhetrs.f"
	goto L90;
#line 462 "zhetrs.f"
L100:
#line 463 "zhetrs.f"
	;
#line 463 "zhetrs.f"
    }

#line 465 "zhetrs.f"
    return 0;

/*     End of ZHETRS */

} /* zhetrs_ */

