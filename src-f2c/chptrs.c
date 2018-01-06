#line 1 "chptrs.f"
/* chptrs.f -- translated by f2c (version 20100827).
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

#line 1 "chptrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

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
/* > CHPTRS solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A stored in packed format using the factorization */
/* > A = U*D*U**H or A = L*D*L**H computed by CHPTRF. */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CHPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by CHPTRF. */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int chptrs_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    cgeru_(integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , cswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , csscal_(integer *, doublereal *, doublecomplex *, integer *), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 156 "chptrs.f"
    /* Parameter adjustments */
#line 156 "chptrs.f"
    --ap;
#line 156 "chptrs.f"
    --ipiv;
#line 156 "chptrs.f"
    b_dim1 = *ldb;
#line 156 "chptrs.f"
    b_offset = 1 + b_dim1;
#line 156 "chptrs.f"
    b -= b_offset;
#line 156 "chptrs.f"

#line 156 "chptrs.f"
    /* Function Body */
#line 156 "chptrs.f"
    *info = 0;
#line 157 "chptrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 158 "chptrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 159 "chptrs.f"
	*info = -1;
#line 160 "chptrs.f"
    } else if (*n < 0) {
#line 161 "chptrs.f"
	*info = -2;
#line 162 "chptrs.f"
    } else if (*nrhs < 0) {
#line 163 "chptrs.f"
	*info = -3;
#line 164 "chptrs.f"
    } else if (*ldb < max(1,*n)) {
#line 165 "chptrs.f"
	*info = -7;
#line 166 "chptrs.f"
    }
#line 167 "chptrs.f"
    if (*info != 0) {
#line 168 "chptrs.f"
	i__1 = -(*info);
#line 168 "chptrs.f"
	xerbla_("CHPTRS", &i__1, (ftnlen)6);
#line 169 "chptrs.f"
	return 0;
#line 170 "chptrs.f"
    }

/*     Quick return if possible */

#line 174 "chptrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 174 "chptrs.f"
	return 0;
#line 174 "chptrs.f"
    }

#line 177 "chptrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 186 "chptrs.f"
	k = *n;
#line 187 "chptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 188 "chptrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 192 "chptrs.f"
	if (k < 1) {
#line 192 "chptrs.f"
	    goto L30;
#line 192 "chptrs.f"
	}

#line 195 "chptrs.f"
	kc -= k;
#line 196 "chptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 202 "chptrs.f"
	    kp = ipiv[k];
#line 203 "chptrs.f"
	    if (kp != k) {
#line 203 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 203 "chptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 209 "chptrs.f"
	    i__1 = k - 1;
#line 209 "chptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 209 "chptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 214 "chptrs.f"
	    i__1 = kc + k - 1;
#line 214 "chptrs.f"
	    s = 1. / ap[i__1].r;
#line 215 "chptrs.f"
	    csscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 216 "chptrs.f"
	    --k;
#line 217 "chptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 223 "chptrs.f"
	    kp = -ipiv[k];
#line 224 "chptrs.f"
	    if (kp != k - 1) {
#line 224 "chptrs.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 224 "chptrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 230 "chptrs.f"
	    i__1 = k - 2;
#line 230 "chptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "chptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc], &c__1, &b[k + b_dim1], ldb, &
		    b[b_dim1 + 1], ldb);
#line 232 "chptrs.f"
	    i__1 = k - 2;
#line 232 "chptrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 232 "chptrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 237 "chptrs.f"
	    i__1 = kc + k - 2;
#line 237 "chptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 238 "chptrs.f"
	    z_div(&z__1, &ap[kc - 1], &akm1k);
#line 238 "chptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 239 "chptrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 239 "chptrs.f"
	    z_div(&z__1, &ap[kc + k - 1], &z__2);
#line 239 "chptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 240 "chptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 240 "chptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 240 "chptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 241 "chptrs.f"
	    i__1 = *nrhs;
#line 241 "chptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 242 "chptrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 242 "chptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 243 "chptrs.f"
		d_cnjg(&z__2, &akm1k);
#line 243 "chptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 243 "chptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 244 "chptrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 244 "chptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 244 "chptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 244 "chptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 244 "chptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 245 "chptrs.f"
		i__2 = k + j * b_dim1;
#line 245 "chptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 245 "chptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 245 "chptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 245 "chptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 246 "chptrs.f"
/* L20: */
#line 246 "chptrs.f"
	    }
#line 247 "chptrs.f"
	    kc = kc - k + 1;
#line 248 "chptrs.f"
	    k += -2;
#line 249 "chptrs.f"
	}

#line 251 "chptrs.f"
	goto L10;
#line 252 "chptrs.f"
L30:

/*        Next solve U**H *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 259 "chptrs.f"
	k = 1;
#line 260 "chptrs.f"
	kc = 1;
#line 261 "chptrs.f"
L40:

/*        If K > N, exit from loop. */

#line 265 "chptrs.f"
	if (k > *n) {
#line 265 "chptrs.f"
	    goto L50;
#line 265 "chptrs.f"
	}

#line 268 "chptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**H(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 275 "chptrs.f"
	    if (k > 1) {
#line 276 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 277 "chptrs.f"
		i__1 = k - 1;
#line 277 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 277 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)19);
#line 279 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 280 "chptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 284 "chptrs.f"
	    kp = ipiv[k];
#line 285 "chptrs.f"
	    if (kp != k) {
#line 285 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 285 "chptrs.f"
	    }
#line 287 "chptrs.f"
	    kc += k;
#line 288 "chptrs.f"
	    ++k;
#line 289 "chptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 296 "chptrs.f"
	    if (k > 1) {
#line 297 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 298 "chptrs.f"
		i__1 = k - 1;
#line 298 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 298 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)19);
#line 300 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);

#line 302 "chptrs.f"
		clacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 303 "chptrs.f"
		i__1 = k - 1;
#line 303 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 303 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &ap[kc + k], &c__1, &c_b1, &b[k + 1 + b_dim1], 
			ldb, (ftnlen)19);
#line 305 "chptrs.f"
		clacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 306 "chptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 310 "chptrs.f"
	    kp = -ipiv[k];
#line 311 "chptrs.f"
	    if (kp != k) {
#line 311 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 311 "chptrs.f"
	    }
#line 313 "chptrs.f"
	    kc = kc + (k << 1) + 1;
#line 314 "chptrs.f"
	    k += 2;
#line 315 "chptrs.f"
	}

#line 317 "chptrs.f"
	goto L40;
#line 318 "chptrs.f"
L50:

#line 320 "chptrs.f"
	;
#line 320 "chptrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 329 "chptrs.f"
	k = 1;
#line 330 "chptrs.f"
	kc = 1;
#line 331 "chptrs.f"
L60:

/*        If K > N, exit from loop. */

#line 335 "chptrs.f"
	if (k > *n) {
#line 335 "chptrs.f"
	    goto L80;
#line 335 "chptrs.f"
	}

#line 338 "chptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 344 "chptrs.f"
	    kp = ipiv[k];
#line 345 "chptrs.f"
	    if (kp != k) {
#line 345 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 345 "chptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 351 "chptrs.f"
	    if (k < *n) {
#line 351 "chptrs.f"
		i__1 = *n - k;
#line 351 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 351 "chptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + 1], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 1 + b_dim1], ldb);
#line 351 "chptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 357 "chptrs.f"
	    i__1 = kc;
#line 357 "chptrs.f"
	    s = 1. / ap[i__1].r;
#line 358 "chptrs.f"
	    csscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 359 "chptrs.f"
	    kc = kc + *n - k + 1;
#line 360 "chptrs.f"
	    ++k;
#line 361 "chptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 367 "chptrs.f"
	    kp = -ipiv[k];
#line 368 "chptrs.f"
	    if (kp != k + 1) {
#line 368 "chptrs.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 368 "chptrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 374 "chptrs.f"
	    if (k < *n - 1) {
#line 375 "chptrs.f"
		i__1 = *n - k - 1;
#line 375 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 375 "chptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + 2], &c__1, &b[k + b_dim1],
			 ldb, &b[k + 2 + b_dim1], ldb);
#line 377 "chptrs.f"
		i__1 = *n - k - 1;
#line 377 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 377 "chptrs.f"
		cgeru_(&i__1, nrhs, &z__1, &ap[kc + *n - k + 2], &c__1, &b[k 
			+ 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 379 "chptrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 383 "chptrs.f"
	    i__1 = kc + 1;
#line 383 "chptrs.f"
	    akm1k.r = ap[i__1].r, akm1k.i = ap[i__1].i;
#line 384 "chptrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 384 "chptrs.f"
	    z_div(&z__1, &ap[kc], &z__2);
#line 384 "chptrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 385 "chptrs.f"
	    z_div(&z__1, &ap[kc + *n - k + 1], &akm1k);
#line 385 "chptrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 386 "chptrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 386 "chptrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 386 "chptrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 387 "chptrs.f"
	    i__1 = *nrhs;
#line 387 "chptrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 388 "chptrs.f"
		d_cnjg(&z__2, &akm1k);
#line 388 "chptrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 388 "chptrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 389 "chptrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 389 "chptrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 390 "chptrs.f"
		i__2 = k + j * b_dim1;
#line 390 "chptrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 390 "chptrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 390 "chptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 390 "chptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 391 "chptrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 391 "chptrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 391 "chptrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 391 "chptrs.f"
		z_div(&z__1, &z__2, &denom);
#line 391 "chptrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 392 "chptrs.f"
/* L70: */
#line 392 "chptrs.f"
	    }
#line 393 "chptrs.f"
	    kc = kc + (*n - k << 1) + 1;
#line 394 "chptrs.f"
	    k += 2;
#line 395 "chptrs.f"
	}

#line 397 "chptrs.f"
	goto L60;
#line 398 "chptrs.f"
L80:

/*        Next solve L**H *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 405 "chptrs.f"
	k = *n;
#line 406 "chptrs.f"
	kc = *n * (*n + 1) / 2 + 1;
#line 407 "chptrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 411 "chptrs.f"
	if (k < 1) {
#line 411 "chptrs.f"
	    goto L100;
#line 411 "chptrs.f"
	}

#line 414 "chptrs.f"
	kc -= *n - k + 1;
#line 415 "chptrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**H(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 422 "chptrs.f"
	    if (k < *n) {
#line 423 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 424 "chptrs.f"
		i__1 = *n - k;
#line 424 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 424 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 427 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 428 "chptrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 432 "chptrs.f"
	    kp = ipiv[k];
#line 433 "chptrs.f"
	    if (kp != k) {
#line 433 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 433 "chptrs.f"
	    }
#line 435 "chptrs.f"
	    --k;
#line 436 "chptrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 443 "chptrs.f"
	    if (k < *n) {
#line 444 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 445 "chptrs.f"
		i__1 = *n - k;
#line 445 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 445 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 448 "chptrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);

#line 450 "chptrs.f"
		clacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 451 "chptrs.f"
		i__1 = *n - k;
#line 451 "chptrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 451 "chptrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &ap[kc - (*n - k)], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)19);
#line 454 "chptrs.f"
		clacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 455 "chptrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 459 "chptrs.f"
	    kp = -ipiv[k];
#line 460 "chptrs.f"
	    if (kp != k) {
#line 460 "chptrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 460 "chptrs.f"
	    }
#line 462 "chptrs.f"
	    kc -= *n - k + 2;
#line 463 "chptrs.f"
	    k += -2;
#line 464 "chptrs.f"
	}

#line 466 "chptrs.f"
	goto L90;
#line 467 "chptrs.f"
L100:
#line 468 "chptrs.f"
	;
#line 468 "chptrs.f"
    }

#line 470 "chptrs.f"
    return 0;

/*     End of CHPTRS */

} /* chptrs_ */

