#line 1 "csptrs.f"
/* csptrs.f -- translated by f2c (version 20100827).
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

#line 1 "csptrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            AP( * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A stored in packed format using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by CSPTRF. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by CSPTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int csptrs_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    cgeru_(integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 155 "csptrs.f"
    /* Parameter adjustments */
#line 155 "csptrs.f"
    --ap;
#line 155 "csptrs.f"
    --ipiv;
#line 155 "csptrs.f"
    b_dim1 = *ldb;
#line 155 "csptrs.f"
    b_offset = 1 + b_dim1;
#line 155 "csptrs.f"
    b -= b_offset;
#line 155 "csptrs.f"

#line 155 "csptrs.f"
    /* Function Body */
#line 155 "csptrs.f"
    *info = 0;
#line 156 "csptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "csptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "csptrs.f"
	*info = -1;
#line 159 "csptrs.f"
    } else if (*n < 0) {
#line 160 "csptrs.f"
	*info = -2;
#line 161 "csptrs.f"
    } else if (*nrhs < 0) {
#line 162 "csptrs.f"
	*info = -3;
#line 163 "csptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 164 "csptrs.f"
	*info = -7;
#line 165 "csptrs.f"
    }
#line 166 "csptrs.f"
    if (*info != 0) {
#line 167 "csptrs.f"
	i__1 = -(*info);
#line 167 "csptrs.f"
	xerbla_("CSPTRS", &i__1, (ftnlen)6);
#line 168 "csptrs.f"
	return 0;
#line 169 "csptrs.f"
    }

/*     Quick return if possible */

#line 173 "csptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 173 "csptrs.f"
	return 0;
#line 173 "csptrs.f"
    }

#line 176 "csptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 185 "csptrs.f"
	k = *n;
#line 186 "csptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 187 "csptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 191 "csptrs.f"
	if (k < 1) {
#line 191 "csptrs.f"
	    goto L30;
#line 191 "csptrs.f"
	}

#line 194 "csptrs.f"
	kc -= k;
#line 195 "csptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 201 "csptrs.f"
	    kp = ipiv[k];
#line 202 "csptrs.f"
	    if (kp != k) {
#line 202 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 202 "csptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 208 "csptrs.f"
	    i__1 = k - 1;
#line 208 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 208 "csptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 213 "csptrs.f"
	    z_div(&z__1, &c_b1, &ap[kc + k - 1]);
#line 213 "csptrs.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 214 "csptrs.f"
	    --k;
#line 215 "csptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 221 "csptrs.f"
	    kp = -ipiv[k];
#line 222 "csptrs.f"
	    if (kp != k - 1) {
#line 222 "csptrs.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 222 "csptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 228 "csptrs.f"
	    i__1 = k - 2;
#line 228 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 228 "csptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);
#line 230 "csptrs.f"
	    i__1 = k - 2;
#line 230 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "csptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 235 "csptrs.f"
	    i__1 = kc + k - 2;
#line 235 "csptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 236 "csptrs.f"
	    z_div(&z__1, &ap[kc - 1], &akm1k);
#line 236 "csptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 237 "csptrs.f"
	    z_div(&z__1, &ap[kc + k - 1], &akm1k);
#line 237 "csptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 238 "csptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 238 "csptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 238 "csptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 239 "csptrs.f"
	    i__1 = *nrhs;
#line 239 "csptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 240 "csptrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 240 "csptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 241 "csptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 241 "csptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 242 "csptrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 242 "csptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 242 "csptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 242 "csptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 242 "csptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 243 "csptrs.f"
		i__2 = k + j * b_dim1;
#line 243 "csptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 243 "csptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 243 "csptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 243 "csptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 244 "csptrs.f"
/* L20: */
#line 244 "csptrs.f"
	    }
#line 245 "csptrs.f"
	    kc = kc - k + 1;
#line 246 "csptrs.f"
	    k += -2;
#line 247 "csptrs.f"
	}

#line 249 "csptrs.f"
	goto L10;
#line 250 "csptrs.f"
L30:

/*        Next solve U**T*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 257 "csptrs.f"
	k = 1;
#line 258 "csptrs.f"
	kc = 1;
#line 259 "csptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 263 "csptrs.f"
	if (k > *n) {
#line 263 "csptrs.f"
	    goto L50;
#line 263 "csptrs.f"
	}

#line 266 "csptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 273 "csptrs.f"
	    i__1 = k - 1;
#line 273 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 273 "csptrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and IPIV(K). */

#line 278 "csptrs.f"
	    kp = ipiv[k];
#line 279 "csptrs.f"
	    if (kp != k) {
#line 279 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 279 "csptrs.f"
	    }
#line 281 "csptrs.f"
	    kc += k;
#line 282 "csptrs.f"
	    ++k;
#line 283 "csptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 290 "csptrs.f"
	    i__1 = k - 1;
#line 290 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 290 "csptrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9);
#line 292 "csptrs.f"
	    i__1 = k - 1;
#line 292 "csptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 292 "csptrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &ap[kc 
		    + k], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb, (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 297 "csptrs.f"
	    kp = -ipiv[k];
#line 298 "csptrs.f"
	    if (kp != k) {
#line 298 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 298 "csptrs.f"
	    }
#line 300 "csptrs.f"
	    kc = kc + (k << 1) + 1;
#line 301 "csptrs.f"
	    k += 2;
#line 302 "csptrs.f"
	}

#line 304 "csptrs.f"
	goto L40;
#line 305 "csptrs.f"
L50:

#line 307 "csptrs.f"
	;
#line 307 "csptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 316 "csptrs.f"
	k = 1;
#line 317 "csptrs.f"
	kc = 1;
#line 318 "csptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "csptrs.f"
	if (k > *n) {
#line 322 "csptrs.f"
	    goto L80;
#line 322 "csptrs.f"
	}

#line 325 "csptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "csptrs.f"
	    kp = ipiv[k];
#line 332 "csptrs.f"
	    if (kp != k) {
#line 332 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "csptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "csptrs.f"
	    if (k < *n) {
#line 338 "csptrs.f"
		i__1 = *n - k;
#line 338 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 338 "csptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + 1], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "csptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "csptrs.f"
	    z_div(&z__1, &c_b1, &ap[kc]);
#line 344 "csptrs.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 345 "csptrs.f"
	    kc = kc + *n - k + 1;
#line 346 "csptrs.f"
	    ++k;
#line 347 "csptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 353 "csptrs.f"
	    kp = -ipiv[k];
#line 354 "csptrs.f"
	    if (kp != k + 1) {
#line 354 "csptrs.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 354 "csptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 360 "csptrs.f"
	    if (k < *n - 1) {
#line 361 "csptrs.f"
		i__1 = *n - k - 1;
#line 361 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 361 "csptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + 2], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 2 + b_dim1], ldb);
#line 363 "csptrs.f"
		i__1 = *n - k - 1;
#line 363 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 363 "csptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + *n - k + 2], &c__1, &b[k 
			+ 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 365 "csptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 369 "csptrs.f"
	    i__1 = kc + 1;
#line 369 "csptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 370 "csptrs.f"
	    z_div(&z__1, &ap[kc], &akm1k);
#line 370 "csptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 371 "csptrs.f"
	    z_div(&z__1, &ap[kc + *n - k + 1], &akm1k);
#line 371 "csptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 372 "csptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 372 "csptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 372 "csptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 373 "csptrs.f"
	    i__1 = *nrhs;
#line 373 "csptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 374 "csptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 374 "csptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 375 "csptrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 375 "csptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 376 "csptrs.f"
		i__2 = k + j * b_dim1;
#line 376 "csptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 376 "csptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 376 "csptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 376 "csptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 377 "csptrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 377 "csptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 377 "csptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 377 "csptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 377 "csptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 378 "csptrs.f"
/* L70: */
#line 378 "csptrs.f"
	    }
#line 379 "csptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 380 "csptrs.f"
	    k += 2;
#line 381 "csptrs.f"
	}

#line 383 "csptrs.f"
	goto L60;
#line 384 "csptrs.f"
L80:

/*        Next solve L**T*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 391 "csptrs.f"
	k = *n;
#line 392 "csptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 393 "csptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 397 "csptrs.f"
	if (k < 1) {
#line 397 "csptrs.f"
	    goto L100;
#line 397 "csptrs.f"
	}

#line 400 "csptrs.f"
	kc -= *n - k + 1;
#line 401 "csptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 408 "csptrs.f"
	    if (k < *n) {
#line 408 "csptrs.f"
		i__1 = *n - k;
#line 408 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 408 "csptrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 408 "csptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 414 "csptrs.f"
	    kp = ipiv[k];
#line 415 "csptrs.f"
	    if (kp != k) {
#line 415 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 415 "csptrs.f"
	    }
#line 417 "csptrs.f"
	    --k;
#line 418 "csptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 425 "csptrs.f"
	    if (k < *n) {
#line 426 "csptrs.f"
		i__1 = *n - k;
#line 426 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 426 "csptrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 428 "csptrs.f"
		i__1 = *n - k;
#line 428 "csptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 428 "csptrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &ap[kc - (*n - k)], &c__1, &c_b1, &b[k - 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 431 "csptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 435 "csptrs.f"
	    kp = -ipiv[k];
#line 436 "csptrs.f"
	    if (kp != k) {
#line 436 "csptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 436 "csptrs.f"
	    }
#line 438 "csptrs.f"
	    kc -= *n - k + 2;
#line 439 "csptrs.f"
	    k += -2;
#line 440 "csptrs.f"
	}

#line 442 "csptrs.f"
	goto L90;
#line 443 "csptrs.f"
L100:
#line 444 "csptrs.f"
	;
#line 444 "csptrs.f"
    }

#line 446 "csptrs.f"
    return 0;

/*     End of CSPTRS */

} /* csptrs_ */

