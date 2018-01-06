#line 1 "zsytrs_rook.f"
/* zsytrs_rook.f -- translated by f2c (version 20100827).
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

#line 1 "zsytrs_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZSYTRS_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYTRS_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > ZSYTRS_ROOK solves a system of linear equations A*X = B with */
/* > a complex symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by ZSYTRF_ROOK. */
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
/* >          obtain the factor U or L as computed by ZSYTRF_ROOK. */
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
/* >          as determined by ZSYTRF_ROOK. */
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

/* > \date December 2016 */

/* > \ingroup complex16SYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >   December 2016, Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zsytrs_rook__(char *uplo, integer *n, integer *nrhs, 
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

#line 176 "zsytrs_rook.f"
    /* Parameter adjustments */
#line 176 "zsytrs_rook.f"
    a_dim1 = *lda;
#line 176 "zsytrs_rook.f"
    a_offset = 1 + a_dim1;
#line 176 "zsytrs_rook.f"
    a -= a_offset;
#line 176 "zsytrs_rook.f"
    --ipiv;
#line 176 "zsytrs_rook.f"
    b_dim1 = *ldb;
#line 176 "zsytrs_rook.f"
    b_offset = 1 + b_dim1;
#line 176 "zsytrs_rook.f"
    b -= b_offset;
#line 176 "zsytrs_rook.f"

#line 176 "zsytrs_rook.f"
    /* Function Body */
#line 176 "zsytrs_rook.f"
    *info = 0;
#line 177 "zsytrs_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 178 "zsytrs_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 179 "zsytrs_rook.f"
	*info = -1;
#line 180 "zsytrs_rook.f"
    } else if (*n < 0) {
#line 181 "zsytrs_rook.f"
	*info = -2;
#line 182 "zsytrs_rook.f"
    } else if (*nrhs < 0) {
#line 183 "zsytrs_rook.f"
	*info = -3;
#line 184 "zsytrs_rook.f"
    } else if (*lda < max(1,*n)) {
#line 185 "zsytrs_rook.f"
	*info = -5;
#line 186 "zsytrs_rook.f"
    } else if (*ldb < max(1,*n)) {
#line 187 "zsytrs_rook.f"
	*info = -8;
#line 188 "zsytrs_rook.f"
    }
#line 189 "zsytrs_rook.f"
    if (*info != 0) {
#line 190 "zsytrs_rook.f"
	i__1 = -(*info);
#line 190 "zsytrs_rook.f"
	xerbla_("ZSYTRS_ROOK", &i__1, (ftnlen)11);
#line 191 "zsytrs_rook.f"
	return 0;
#line 192 "zsytrs_rook.f"
    }

/*     Quick return if possible */

#line 196 "zsytrs_rook.f"
    if (*n == 0 || *nrhs == 0) {
#line 196 "zsytrs_rook.f"
	return 0;
#line 196 "zsytrs_rook.f"
    }

#line 199 "zsytrs_rook.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 208 "zsytrs_rook.f"
	k = *n;
#line 209 "zsytrs_rook.f"
L10:

/*        If K < 1, exit from loop. */

#line 213 "zsytrs_rook.f"
	if (k < 1) {
#line 213 "zsytrs_rook.f"
	    goto L30;
#line 213 "zsytrs_rook.f"
	}

#line 216 "zsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 222 "zsytrs_rook.f"
	    kp = ipiv[k];
#line 223 "zsytrs_rook.f"
	    if (kp != k) {
#line 223 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 223 "zsytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 229 "zsytrs_rook.f"
	    i__1 = k - 1;
#line 229 "zsytrs_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 229 "zsytrs_rook.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 234 "zsytrs_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 234 "zsytrs_rook.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 235 "zsytrs_rook.f"
	    --k;
#line 236 "zsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 242 "zsytrs_rook.f"
	    kp = -ipiv[k];
#line 243 "zsytrs_rook.f"
	    if (kp != k) {
#line 243 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 243 "zsytrs_rook.f"
	    }

#line 246 "zsytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 247 "zsytrs_rook.f"
	    if (kp != k - 1) {
#line 247 "zsytrs_rook.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "zsytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 253 "zsytrs_rook.f"
	    if (k > 2) {
#line 254 "zsytrs_rook.f"
		i__1 = k - 2;
#line 254 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 254 "zsytrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
			b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 256 "zsytrs_rook.f"
		i__1 = k - 2;
#line 256 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 256 "zsytrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &
			b[k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 258 "zsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 262 "zsytrs_rook.f"
	    i__1 = k - 1 + k * a_dim1;
#line 262 "zsytrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 263 "zsytrs_rook.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 263 "zsytrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 264 "zsytrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 264 "zsytrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 265 "zsytrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 265 "zsytrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 265 "zsytrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 266 "zsytrs_rook.f"
	    i__1 = *nrhs;
#line 266 "zsytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 267 "zsytrs_rook.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 267 "zsytrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 268 "zsytrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 268 "zsytrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 269 "zsytrs_rook.f"
		i__2 = k - 1 + j * b_dim1;
#line 269 "zsytrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 269 "zsytrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 269 "zsytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 269 "zsytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 270 "zsytrs_rook.f"
		i__2 = k + j * b_dim1;
#line 270 "zsytrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 270 "zsytrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 270 "zsytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 270 "zsytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 271 "zsytrs_rook.f"
/* L20: */
#line 271 "zsytrs_rook.f"
	    }
#line 272 "zsytrs_rook.f"
	    k += -2;
#line 273 "zsytrs_rook.f"
	}

#line 275 "zsytrs_rook.f"
	goto L10;
#line 276 "zsytrs_rook.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 283 "zsytrs_rook.f"
	k = 1;
#line 284 "zsytrs_rook.f"
L40:

/*        If K > N, exit from loop. */

#line 288 "zsytrs_rook.f"
	if (k > *n) {
#line 288 "zsytrs_rook.f"
	    goto L50;
#line 288 "zsytrs_rook.f"
	}

#line 291 "zsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 298 "zsytrs_rook.f"
	    if (k > 1) {
#line 298 "zsytrs_rook.f"
		i__1 = k - 1;
#line 298 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 298 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 298 "zsytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 304 "zsytrs_rook.f"
	    kp = ipiv[k];
#line 305 "zsytrs_rook.f"
	    if (kp != k) {
#line 305 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 305 "zsytrs_rook.f"
	    }
#line 307 "zsytrs_rook.f"
	    ++k;
#line 308 "zsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 315 "zsytrs_rook.f"
	    if (k > 1) {
#line 316 "zsytrs_rook.f"
		i__1 = k - 1;
#line 316 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 316 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 318 "zsytrs_rook.f"
		i__1 = k - 1;
#line 318 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 318 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 320 "zsytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1). */

#line 324 "zsytrs_rook.f"
	    kp = -ipiv[k];
#line 325 "zsytrs_rook.f"
	    if (kp != k) {
#line 325 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 325 "zsytrs_rook.f"
	    }

#line 328 "zsytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 329 "zsytrs_rook.f"
	    if (kp != k + 1) {
#line 329 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 329 "zsytrs_rook.f"
	    }

#line 332 "zsytrs_rook.f"
	    k += 2;
#line 333 "zsytrs_rook.f"
	}

#line 335 "zsytrs_rook.f"
	goto L40;
#line 336 "zsytrs_rook.f"
L50:

#line 338 "zsytrs_rook.f"
	;
#line 338 "zsytrs_rook.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 347 "zsytrs_rook.f"
	k = 1;
#line 348 "zsytrs_rook.f"
L60:

/*        If K > N, exit from loop. */

#line 352 "zsytrs_rook.f"
	if (k > *n) {
#line 352 "zsytrs_rook.f"
	    goto L80;
#line 352 "zsytrs_rook.f"
	}

#line 355 "zsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 361 "zsytrs_rook.f"
	    kp = ipiv[k];
#line 362 "zsytrs_rook.f"
	    if (kp != k) {
#line 362 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 362 "zsytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 368 "zsytrs_rook.f"
	    if (k < *n) {
#line 368 "zsytrs_rook.f"
		i__1 = *n - k;
#line 368 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 368 "zsytrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 368 "zsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 374 "zsytrs_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 374 "zsytrs_rook.f"
	    zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 375 "zsytrs_rook.f"
	    ++k;
#line 376 "zsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1) */

#line 382 "zsytrs_rook.f"
	    kp = -ipiv[k];
#line 383 "zsytrs_rook.f"
	    if (kp != k) {
#line 383 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 383 "zsytrs_rook.f"
	    }

#line 386 "zsytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 387 "zsytrs_rook.f"
	    if (kp != k + 1) {
#line 387 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 387 "zsytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 393 "zsytrs_rook.f"
	    if (k < *n - 1) {
#line 394 "zsytrs_rook.f"
		i__1 = *n - k - 1;
#line 394 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 394 "zsytrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 396 "zsytrs_rook.f"
		i__1 = *n - k - 1;
#line 396 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 396 "zsytrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 398 "zsytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 402 "zsytrs_rook.f"
	    i__1 = k + 1 + k * a_dim1;
#line 402 "zsytrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 403 "zsytrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 403 "zsytrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 404 "zsytrs_rook.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 404 "zsytrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 405 "zsytrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 405 "zsytrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 405 "zsytrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 406 "zsytrs_rook.f"
	    i__1 = *nrhs;
#line 406 "zsytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 407 "zsytrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 407 "zsytrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 408 "zsytrs_rook.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 408 "zsytrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 409 "zsytrs_rook.f"
		i__2 = k + j * b_dim1;
#line 409 "zsytrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 409 "zsytrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 409 "zsytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 409 "zsytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 410 "zsytrs_rook.f"
		i__2 = k + 1 + j * b_dim1;
#line 410 "zsytrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 410 "zsytrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 410 "zsytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 410 "zsytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 411 "zsytrs_rook.f"
/* L70: */
#line 411 "zsytrs_rook.f"
	    }
#line 412 "zsytrs_rook.f"
	    k += 2;
#line 413 "zsytrs_rook.f"
	}

#line 415 "zsytrs_rook.f"
	goto L60;
#line 416 "zsytrs_rook.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 423 "zsytrs_rook.f"
	k = *n;
#line 424 "zsytrs_rook.f"
L90:

/*        If K < 1, exit from loop. */

#line 428 "zsytrs_rook.f"
	if (k < 1) {
#line 428 "zsytrs_rook.f"
	    goto L100;
#line 428 "zsytrs_rook.f"
	}

#line 431 "zsytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 438 "zsytrs_rook.f"
	    if (k < *n) {
#line 438 "zsytrs_rook.f"
		i__1 = *n - k;
#line 438 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 438 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 438 "zsytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 444 "zsytrs_rook.f"
	    kp = ipiv[k];
#line 445 "zsytrs_rook.f"
	    if (kp != k) {
#line 445 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 445 "zsytrs_rook.f"
	    }
#line 447 "zsytrs_rook.f"
	    --k;
#line 448 "zsytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 455 "zsytrs_rook.f"
	    if (k < *n) {
#line 456 "zsytrs_rook.f"
		i__1 = *n - k;
#line 456 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 456 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 458 "zsytrs_rook.f"
		i__1 = *n - k;
#line 458 "zsytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 458 "zsytrs_rook.f"
		zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)9);
#line 461 "zsytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 465 "zsytrs_rook.f"
	    kp = -ipiv[k];
#line 466 "zsytrs_rook.f"
	    if (kp != k) {
#line 466 "zsytrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 466 "zsytrs_rook.f"
	    }

#line 469 "zsytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 470 "zsytrs_rook.f"
	    if (kp != k - 1) {
#line 470 "zsytrs_rook.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 470 "zsytrs_rook.f"
	    }

#line 473 "zsytrs_rook.f"
	    k += -2;
#line 474 "zsytrs_rook.f"
	}

#line 476 "zsytrs_rook.f"
	goto L90;
#line 477 "zsytrs_rook.f"
L100:
#line 478 "zsytrs_rook.f"
	;
#line 478 "zsytrs_rook.f"
    }

#line 480 "zsytrs_rook.f"
    return 0;

/*     End of ZSYTRS_ROOK */

} /* zsytrs_rook__ */

