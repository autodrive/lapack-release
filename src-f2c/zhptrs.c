#line 1 "zhptrs.f"
/* zhptrs.f -- translated by f2c (version 20100827).
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

#line 1 "zhptrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

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
/* > ZHPTRS solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A stored in packed format using the factorization */
/* > A = U*D*U**H or A = L*D*L**H computed by ZHPTRF. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by ZHPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHPTRF. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static doublecomplex ak, bk;
    static integer kc, kp;
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

#line 156 "zhptrs.f"
    /* Parameter adjustments */
#line 156 "zhptrs.f"
    --ap;
#line 156 "zhptrs.f"
    --ipiv;
#line 156 "zhptrs.f"
    b_dim1 = *ldb;
#line 156 "zhptrs.f"
    b_offset = 1 + b_dim1;
#line 156 "zhptrs.f"
    b -= b_offset;
#line 156 "zhptrs.f"

#line 156 "zhptrs.f"
    /* Function Body */
#line 156 "zhptrs.f"
    *info = 0;
#line 157 "zhptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 158 "zhptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 159 "zhptrs.f"
	*info = -1;
#line 160 "zhptrs.f"
    } else if (*n < 0) {
#line 161 "zhptrs.f"
	*info = -2;
#line 162 "zhptrs.f"
    } else if (*nrhs < 0) {
#line 163 "zhptrs.f"
	*info = -3;
#line 164 "zhptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 165 "zhptrs.f"
	*info = -7;
#line 166 "zhptrs.f"
    }
#line 167 "zhptrs.f"
    if (*info != 0) {
#line 168 "zhptrs.f"
	i__1 = -(*info);
#line 168 "zhptrs.f"
	xerbla_("ZHPTRS", &i__1, (ftnlen)6);
#line 169 "zhptrs.f"
	return 0;
#line 170 "zhptrs.f"
    }

/*     Quick return if possible */

#line 174 "zhptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 174 "zhptrs.f"
	return 0;
#line 174 "zhptrs.f"
    }

#line 177 "zhptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 186 "zhptrs.f"
	k = *n;
#line 187 "zhptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 188 "zhptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 192 "zhptrs.f"
	if (k < 1) {
#line 192 "zhptrs.f"
	    goto L30;
#line 192 "zhptrs.f"
	}

#line 195 "zhptrs.f"
	kc -= k;
#line 196 "zhptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 202 "zhptrs.f"
	    kp = ipiv[k];
#line 203 "zhptrs.f"
	    if (kp != k) {
#line 203 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 203 "zhptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 209 "zhptrs.f"
	    i__1 = k - 1;
#line 209 "zhptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 209 "zhptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 214 "zhptrs.f"
	    i__1 = kc + k - 1;
#line 214 "zhptrs.f"
	    s = 1. / ap[i__1].r;
#line 215 "zhptrs.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 216 "zhptrs.f"
	    --k;
#line 217 "zhptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 223 "zhptrs.f"
	    kp = -ipiv[k];
#line 224 "zhptrs.f"
	    if (kp != k - 1) {
#line 224 "zhptrs.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 224 "zhptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 230 "zhptrs.f"
	    i__1 = k - 2;
#line 230 "zhptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "zhptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);
#line 232 "zhptrs.f"
	    i__1 = k - 2;
#line 232 "zhptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 232 "zhptrs.f"
	    zgeru_(&i__1, nrhs, &z__1, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 237 "zhptrs.f"
	    i__1 = kc + k - 2;
#line 237 "zhptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 238 "zhptrs.f"
	    z_div(&z__1, &ap[kc - 1], &akm1k);
#line 238 "zhptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 239 "zhptrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 239 "zhptrs.f"
	    z_div(&z__1, &ap[kc + k - 1], &z__2);
#line 239 "zhptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 240 "zhptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 240 "zhptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 240 "zhptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 241 "zhptrs.f"
	    i__1 = *nrhs;
#line 241 "zhptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 242 "zhptrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 242 "zhptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 243 "zhptrs.f"
		d_cnjg(&z__2, &akm1k);
#line 243 "zhptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 243 "zhptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 244 "zhptrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 244 "zhptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 244 "zhptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 244 "zhptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 244 "zhptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 245 "zhptrs.f"
		i__2 = k + j * b_dim1;
#line 245 "zhptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 245 "zhptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 245 "zhptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 245 "zhptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 246 "zhptrs.f"
/* L20: */
#line 246 "zhptrs.f"
	    }
#line 247 "zhptrs.f"
	    kc = kc - k + 1;
#line 248 "zhptrs.f"
	    k += -2;
#line 249 "zhptrs.f"
	}

#line 251 "zhptrs.f"
	goto L10;
#line 252 "zhptrs.f"
L30:

/*        Next solve U**H *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 259 "zhptrs.f"
	k = 1;
#line 260 "zhptrs.f"
	kc = 1;
#line 261 "zhptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 265 "zhptrs.f"
	if (k > *n) {
#line 265 "zhptrs.f"
	    goto L50;
#line 265 "zhptrs.f"
	}

#line 268 "zhptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**H(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 275 "zhptrs.f"
	    if (k > 1) {
#line 276 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 277 "zhptrs.f"
		i__1 = k - 1;
#line 277 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 277 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)19);
#line 279 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 280 "zhptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 284 "zhptrs.f"
	    kp = ipiv[k];
#line 285 "zhptrs.f"
	    if (kp != k) {
#line 285 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 285 "zhptrs.f"
	    }
#line 287 "zhptrs.f"
	    kc += k;
#line 288 "zhptrs.f"
	    ++k;
#line 289 "zhptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 296 "zhptrs.f"
	    if (k > 1) {
#line 297 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 298 "zhptrs.f"
		i__1 = k - 1;
#line 298 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 298 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)19);
#line 300 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 302 "zhptrs.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 303 "zhptrs.f"
		i__1 = k - 1;
#line 303 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 303 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc + k], &c__1, &c_b1, &b[k + 1 + b_dim1], 
			ldb, (ftnlen)19);
#line 305 "zhptrs.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 306 "zhptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 310 "zhptrs.f"
	    kp = -ipiv[k];
#line 311 "zhptrs.f"
	    if (kp != k) {
#line 311 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 311 "zhptrs.f"
	    }
#line 313 "zhptrs.f"
	    kc = kc + (k << 1) + 1;
#line 314 "zhptrs.f"
	    k += 2;
#line 315 "zhptrs.f"
	}

#line 317 "zhptrs.f"
	goto L40;
#line 318 "zhptrs.f"
L50:

#line 320 "zhptrs.f"
	;
#line 320 "zhptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 329 "zhptrs.f"
	k = 1;
#line 330 "zhptrs.f"
	kc = 1;
#line 331 "zhptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 335 "zhptrs.f"
	if (k > *n) {
#line 335 "zhptrs.f"
	    goto L80;
#line 335 "zhptrs.f"
	}

#line 338 "zhptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 344 "zhptrs.f"
	    kp = ipiv[k];
#line 345 "zhptrs.f"
	    if (kp != k) {
#line 345 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 345 "zhptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 351 "zhptrs.f"
	    if (k < *n) {
#line 351 "zhptrs.f"
		i__1 = *n - k;
#line 351 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 351 "zhptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + 1], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 1 + b_dim1], ldb);
#line 351 "zhptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 357 "zhptrs.f"
	    i__1 = kc;
#line 357 "zhptrs.f"
	    s = 1. / ap[i__1].r;
#line 358 "zhptrs.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 359 "zhptrs.f"
	    kc = kc + *n - k + 1;
#line 360 "zhptrs.f"
	    ++k;
#line 361 "zhptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 367 "zhptrs.f"
	    kp = -ipiv[k];
#line 368 "zhptrs.f"
	    if (kp != k + 1) {
#line 368 "zhptrs.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 368 "zhptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 374 "zhptrs.f"
	    if (k < *n - 1) {
#line 375 "zhptrs.f"
		i__1 = *n - k - 1;
#line 375 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 375 "zhptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + 2], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 2 + b_dim1], ldb);
#line 377 "zhptrs.f"
		i__1 = *n - k - 1;
#line 377 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 377 "zhptrs.f"
		zgeru_(&i__1, nrhs, &z__1, &ap[kc + *n - k + 2], &c__1, &b[k 
			+ 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 379 "zhptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 383 "zhptrs.f"
	    i__1 = kc + 1;
#line 383 "zhptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 384 "zhptrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 384 "zhptrs.f"
	    z_div(&z__1, &ap[kc], &z__2);
#line 384 "zhptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 385 "zhptrs.f"
	    z_div(&z__1, &ap[kc + *n - k + 1], &akm1k);
#line 385 "zhptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 386 "zhptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 386 "zhptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 386 "zhptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 387 "zhptrs.f"
	    i__1 = *nrhs;
#line 387 "zhptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 388 "zhptrs.f"
		d_cnjg(&z__2, &akm1k);
#line 388 "zhptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 388 "zhptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 389 "zhptrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 389 "zhptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 390 "zhptrs.f"
		i__2 = k + j * b_dim1;
#line 390 "zhptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 390 "zhptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 390 "zhptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 390 "zhptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 391 "zhptrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 391 "zhptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 391 "zhptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 391 "zhptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 391 "zhptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 392 "zhptrs.f"
/* L70: */
#line 392 "zhptrs.f"
	    }
#line 393 "zhptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 394 "zhptrs.f"
	    k += 2;
#line 395 "zhptrs.f"
	}

#line 397 "zhptrs.f"
	goto L60;
#line 398 "zhptrs.f"
L80:

/*        Next solve L**H *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 405 "zhptrs.f"
	k = *n;
#line 406 "zhptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 407 "zhptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 411 "zhptrs.f"
	if (k < 1) {
#line 411 "zhptrs.f"
	    goto L100;
#line 411 "zhptrs.f"
	}

#line 414 "zhptrs.f"
	kc -= *n - k + 1;
#line 415 "zhptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**H(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 422 "zhptrs.f"
	    if (k < *n) {
#line 423 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 424 "zhptrs.f"
		i__1 = *n - k;
#line 424 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 424 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 427 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 428 "zhptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 432 "zhptrs.f"
	    kp = ipiv[k];
#line 433 "zhptrs.f"
	    if (kp != k) {
#line 433 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 433 "zhptrs.f"
	    }
#line 435 "zhptrs.f"
	    --k;
#line 436 "zhptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 443 "zhptrs.f"
	    if (k < *n) {
#line 444 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 445 "zhptrs.f"
		i__1 = *n - k;
#line 445 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 445 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 448 "zhptrs.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 450 "zhptrs.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 451 "zhptrs.f"
		i__1 = *n - k;
#line 451 "zhptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 451 "zhptrs.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc - (*n - k)], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)19);
#line 454 "zhptrs.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 455 "zhptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 459 "zhptrs.f"
	    kp = -ipiv[k];
#line 460 "zhptrs.f"
	    if (kp != k) {
#line 460 "zhptrs.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 460 "zhptrs.f"
	    }
#line 462 "zhptrs.f"
	    kc -= *n - k + 2;
#line 463 "zhptrs.f"
	    k += -2;
#line 464 "zhptrs.f"
	}

#line 466 "zhptrs.f"
	goto L90;
#line 467 "zhptrs.f"
L100:
#line 468 "zhptrs.f"
	;
#line 468 "zhptrs.f"
    }

#line 470 "zhptrs.f"
    return 0;

/*     End of ZHPTRS */

} /* zhptrs_ */

