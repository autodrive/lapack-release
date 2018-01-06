#line 1 "csytrs.f"
/* csytrs.f -- translated by f2c (version 20100827).
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

#line 1 "csytrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSYTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSYTRF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSYTRF. */
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
/* >          as determined by CSYTRF. */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csytrs_(char *uplo, integer *n, integer *nrhs, 
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

#line 160 "csytrs.f"
    /* Parameter adjustments */
#line 160 "csytrs.f"
    a_dim1 = *lda;
#line 160 "csytrs.f"
    a_offset = 1 + a_dim1;
#line 160 "csytrs.f"
    a -= a_offset;
#line 160 "csytrs.f"
    --ipiv;
#line 160 "csytrs.f"
    b_dim1 = *ldb;
#line 160 "csytrs.f"
    b_offset = 1 + b_dim1;
#line 160 "csytrs.f"
    b -= b_offset;
#line 160 "csytrs.f"

#line 160 "csytrs.f"
    /* Function Body */
#line 160 "csytrs.f"
    *info = 0;
#line 161 "csytrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "csytrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "csytrs.f"
	*info = -1;
#line 164 "csytrs.f"
    } else if (*n < 0) {
#line 165 "csytrs.f"
	*info = -2;
#line 166 "csytrs.f"
    } else if (*nrhs < 0) {
#line 167 "csytrs.f"
	*info = -3;
#line 168 "csytrs.f"
    } else if (*lda < max(1,*n)) {
#line 169 "csytrs.f"
	*info = -5;
#line 170 "csytrs.f"
    } else if (*ldb < max(1,*n)) {
#line 171 "csytrs.f"
	*info = -8;
#line 172 "csytrs.f"
    }
#line 173 "csytrs.f"
    if (*info != 0) {
#line 174 "csytrs.f"
	i__1 = -(*info);
#line 174 "csytrs.f"
	xerbla_("CSYTRS", &i__1, (ftnlen)6);
#line 175 "csytrs.f"
	return 0;
#line 176 "csytrs.f"
    }

/*     Quick return if possible */

#line 180 "csytrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 180 "csytrs.f"
	return 0;
#line 180 "csytrs.f"
    }

#line 183 "csytrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 192 "csytrs.f"
	k = *n;
#line 193 "csytrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 197 "csytrs.f"
	if (k < 1) {
#line 197 "csytrs.f"
	    goto L30;
#line 197 "csytrs.f"
	}

#line 200 "csytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 206 "csytrs.f"
	    kp = ipiv[k];
#line 207 "csytrs.f"
	    if (kp != k) {
#line 207 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 207 "csytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 213 "csytrs.f"
	    i__1 = k - 1;
#line 213 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 213 "csytrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 218 "csytrs.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 218 "csytrs.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 219 "csytrs.f"
	    --k;
#line 220 "csytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 226 "csytrs.f"
	    kp = -ipiv[k];
#line 227 "csytrs.f"
	    if (kp != k - 1) {
#line 227 "csytrs.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 227 "csytrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 233 "csytrs.f"
	    i__1 = k - 2;
#line 233 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 233 "csytrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 235 "csytrs.f"
	    i__1 = k - 2;
#line 235 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 235 "csytrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k 
		    - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 240 "csytrs.f"
	    i__1 = k - 1 + k * a_dim1;
#line 240 "csytrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 241 "csytrs.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 241 "csytrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 242 "csytrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 242 "csytrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 243 "csytrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 243 "csytrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 243 "csytrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 244 "csytrs.f"
	    i__1 = *nrhs;
#line 244 "csytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 245 "csytrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 245 "csytrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 246 "csytrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 246 "csytrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 247 "csytrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 247 "csytrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 247 "csytrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 247 "csytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 247 "csytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 248 "csytrs.f"
		i__2 = k + j * b_dim1;
#line 248 "csytrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 248 "csytrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 248 "csytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 248 "csytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 249 "csytrs.f"
/* L20: */
#line 249 "csytrs.f"
	    }
#line 250 "csytrs.f"
	    k += -2;
#line 251 "csytrs.f"
	}

#line 253 "csytrs.f"
	goto L10;
#line 254 "csytrs.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 261 "csytrs.f"
	k = 1;
#line 262 "csytrs.f"
L40:

/*        If K > N, exit from loop. */

#line 266 "csytrs.f"
	if (k > *n) {
#line 266 "csytrs.f"
	    goto L50;
#line 266 "csytrs.f"
	}

#line 269 "csytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 276 "csytrs.f"
	    i__1 = k - 1;
#line 276 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 276 "csytrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9)
		    ;

/*           Interchange rows K and IPIV(K). */

#line 281 "csytrs.f"
	    kp = ipiv[k];
#line 282 "csytrs.f"
	    if (kp != k) {
#line 282 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "csytrs.f"
	    }
#line 284 "csytrs.f"
	    ++k;
#line 285 "csytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 292 "csytrs.f"
	    i__1 = k - 1;
#line 292 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 292 "csytrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * 
		    a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (ftnlen)9)
		    ;
#line 294 "csytrs.f"
	    i__1 = k - 1;
#line 294 "csytrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 294 "csytrs.f"
	    cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[(k 
		    + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb,
		     (ftnlen)9);

/*           Interchange rows K and -IPIV(K). */

#line 299 "csytrs.f"
	    kp = -ipiv[k];
#line 300 "csytrs.f"
	    if (kp != k) {
#line 300 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 300 "csytrs.f"
	    }
#line 302 "csytrs.f"
	    k += 2;
#line 303 "csytrs.f"
	}

#line 305 "csytrs.f"
	goto L40;
#line 306 "csytrs.f"
L50:

#line 308 "csytrs.f"
	;
#line 308 "csytrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 317 "csytrs.f"
	k = 1;
#line 318 "csytrs.f"
L60:

/*        If K > N, exit from loop. */

#line 322 "csytrs.f"
	if (k > *n) {
#line 322 "csytrs.f"
	    goto L80;
#line 322 "csytrs.f"
	}

#line 325 "csytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 331 "csytrs.f"
	    kp = ipiv[k];
#line 332 "csytrs.f"
	    if (kp != k) {
#line 332 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 332 "csytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 338 "csytrs.f"
	    if (k < *n) {
#line 338 "csytrs.f"
		i__1 = *n - k;
#line 338 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 338 "csytrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 338 "csytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 344 "csytrs.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 344 "csytrs.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 345 "csytrs.f"
	    ++k;
#line 346 "csytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 352 "csytrs.f"
	    kp = -ipiv[k];
#line 353 "csytrs.f"
	    if (kp != k + 1) {
#line 353 "csytrs.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 353 "csytrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 359 "csytrs.f"
	    if (k < *n - 1) {
#line 360 "csytrs.f"
		i__1 = *n - k - 1;
#line 360 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 360 "csytrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 362 "csytrs.f"
		i__1 = *n - k - 1;
#line 362 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 362 "csytrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 364 "csytrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 368 "csytrs.f"
	    i__1 = k + 1 + k * a_dim1;
#line 368 "csytrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 369 "csytrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 369 "csytrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 370 "csytrs.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 370 "csytrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 371 "csytrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 371 "csytrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 371 "csytrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 372 "csytrs.f"
	    i__1 = *nrhs;
#line 372 "csytrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 373 "csytrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 373 "csytrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 374 "csytrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 374 "csytrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 375 "csytrs.f"
		i__2 = k + j * b_dim1;
#line 375 "csytrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 375 "csytrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 375 "csytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 375 "csytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 376 "csytrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 376 "csytrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 376 "csytrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 376 "csytrs.f"
		z_div(&z__1, &z__2, &denom);
#line 376 "csytrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 377 "csytrs.f"
/* L70: */
#line 377 "csytrs.f"
	    }
#line 378 "csytrs.f"
	    k += 2;
#line 379 "csytrs.f"
	}

#line 381 "csytrs.f"
	goto L60;
#line 382 "csytrs.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 389 "csytrs.f"
	k = *n;
#line 390 "csytrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 394 "csytrs.f"
	if (k < 1) {
#line 394 "csytrs.f"
	    goto L100;
#line 394 "csytrs.f"
	}

#line 397 "csytrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 404 "csytrs.f"
	    if (k < *n) {
#line 404 "csytrs.f"
		i__1 = *n - k;
#line 404 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 404 "csytrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 404 "csytrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 410 "csytrs.f"
	    kp = ipiv[k];
#line 411 "csytrs.f"
	    if (kp != k) {
#line 411 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 411 "csytrs.f"
	    }
#line 413 "csytrs.f"
	    --k;
#line 414 "csytrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 421 "csytrs.f"
	    if (k < *n) {
#line 422 "csytrs.f"
		i__1 = *n - k;
#line 422 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 422 "csytrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 424 "csytrs.f"
		i__1 = *n - k;
#line 424 "csytrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 424 "csytrs.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)9);
#line 427 "csytrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 431 "csytrs.f"
	    kp = -ipiv[k];
#line 432 "csytrs.f"
	    if (kp != k) {
#line 432 "csytrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 432 "csytrs.f"
	    }
#line 434 "csytrs.f"
	    k += -2;
#line 435 "csytrs.f"
	}

#line 437 "csytrs.f"
	goto L90;
#line 438 "csytrs.f"
L100:
#line 439 "csytrs.f"
	;
#line 439 "csytrs.f"
    }

#line 441 "csytrs.f"
    return 0;

/*     End of CSYTRS */

} /* csytrs_ */

