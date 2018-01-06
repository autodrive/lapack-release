#line 1 "chetrs.f"
/* chetrs.f -- translated by f2c (version 20100827).
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

#line 1 "chetrs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHETRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > CHETRS solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by CHETRF. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CHETRF. */
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
/* >          as determined by CHETRF. */
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

/* > \date November 2011 */

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
/* Subroutine */ int chetrs_(char *uplo, integer *n, integer *nrhs, 
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

#line 161 "chetrs.f"
    /* Parameter adjustments */
#line 161 "chetrs.f"
    a_dim1 = *lda;
#line 161 "chetrs.f"
    a_offset = 1 + a_dim1;
#line 161 "chetrs.f"
    a -= a_offset;
#line 161 "chetrs.f"
    --ipiv;
#line 161 "chetrs.f"
    b_dim1 = *ldb;
#line 161 "chetrs.f"
    b_offset = 1 + b_dim1;
#line 161 "chetrs.f"
    b -= b_offset;
#line 161 "chetrs.f"

#line 161 "chetrs.f"
    /* Function Body */
#line 161 "chetrs.f"
    *info = 0;
#line 162 "chetrs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 163 "chetrs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 164 "chetrs.f"
	*info = -1;
#line 165 "chetrs.f"
    } else if (*n < 0) {
#line 166 "chetrs.f"
	*info = -2;
#line 167 "chetrs.f"
    } else if (*nrhs < 0) {
#line 168 "chetrs.f"
	*info = -3;
#line 169 "chetrs.f"
    } else if (*lda < max(1,*n)) {
#line 170 "chetrs.f"
	*info = -5;
#line 171 "chetrs.f"
    } else if (*ldb < max(1,*n)) {
#line 172 "chetrs.f"
	*info = -8;
#line 173 "chetrs.f"
    }
#line 174 "chetrs.f"
    if (*info != 0) {
#line 175 "chetrs.f"
	i__1 = -(*info);
#line 175 "chetrs.f"
	xerbla_("CHETRS", &i__1, (ftnlen)6);
#line 176 "chetrs.f"
	return 0;
#line 177 "chetrs.f"
    }

/*     Quick return if possible */

#line 181 "chetrs.f"
    if (*n == 0 || *nrhs == 0) {
#line 181 "chetrs.f"
	return 0;
#line 181 "chetrs.f"
    }

#line 184 "chetrs.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 193 "chetrs.f"
	k = *n;
#line 194 "chetrs.f"
L10:

/*        If K < 1, exit from loop. */

#line 198 "chetrs.f"
	if (k < 1) {
#line 198 "chetrs.f"
	    goto L30;
#line 198 "chetrs.f"
	}

#line 201 "chetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 207 "chetrs.f"
	    kp = ipiv[k];
#line 208 "chetrs.f"
	    if (kp != k) {
#line 208 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 208 "chetrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 214 "chetrs.f"
	    i__1 = k - 1;
#line 214 "chetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 214 "chetrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 219 "chetrs.f"
	    i__1 = k + k * a_dim1;
#line 219 "chetrs.f"
	    s = 1. / a[i__1].r;
#line 220 "chetrs.f"
	    csscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 221 "chetrs.f"
	    --k;
#line 222 "chetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

#line 228 "chetrs.f"
	    kp = -ipiv[k];
#line 229 "chetrs.f"
	    if (kp != k - 1) {
#line 229 "chetrs.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 229 "chetrs.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 235 "chetrs.f"
	    i__1 = k - 2;
#line 235 "chetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 235 "chetrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 237 "chetrs.f"
	    i__1 = k - 2;
#line 237 "chetrs.f"
	    z__1.r = -1., z__1.i = -0.;
#line 237 "chetrs.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k 
		    - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 242 "chetrs.f"
	    i__1 = k - 1 + k * a_dim1;
#line 242 "chetrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 243 "chetrs.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 243 "chetrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 244 "chetrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 244 "chetrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 244 "chetrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 245 "chetrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 245 "chetrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 245 "chetrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 246 "chetrs.f"
	    i__1 = *nrhs;
#line 246 "chetrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 247 "chetrs.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 247 "chetrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 248 "chetrs.f"
		d_cnjg(&z__2, &akm1k);
#line 248 "chetrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 248 "chetrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 249 "chetrs.f"
		i__2 = k - 1 + j * b_dim1;
#line 249 "chetrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 249 "chetrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 249 "chetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 249 "chetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 250 "chetrs.f"
		i__2 = k + j * b_dim1;
#line 250 "chetrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 250 "chetrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 250 "chetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 250 "chetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 251 "chetrs.f"
/* L20: */
#line 251 "chetrs.f"
	    }
#line 252 "chetrs.f"
	    k += -2;
#line 253 "chetrs.f"
	}

#line 255 "chetrs.f"
	goto L10;
#line 256 "chetrs.f"
L30:

/*        Next solve U**H *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 263 "chetrs.f"
	k = 1;
#line 264 "chetrs.f"
L40:

/*        If K > N, exit from loop. */

#line 268 "chetrs.f"
	if (k > *n) {
#line 268 "chetrs.f"
	    goto L50;
#line 268 "chetrs.f"
	}

#line 271 "chetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**H(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 278 "chetrs.f"
	    if (k > 1) {
#line 279 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 280 "chetrs.f"
		i__1 = k - 1;
#line 280 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 280 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 282 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 283 "chetrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 287 "chetrs.f"
	    kp = ipiv[k];
#line 288 "chetrs.f"
	    if (kp != k) {
#line 288 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 288 "chetrs.f"
	    }
#line 290 "chetrs.f"
	    ++k;
#line 291 "chetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 298 "chetrs.f"
	    if (k > 1) {
#line 299 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 300 "chetrs.f"
		i__1 = k - 1;
#line 300 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 300 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 302 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);

#line 304 "chetrs.f"
		clacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 305 "chetrs.f"
		i__1 = k - 1;
#line 305 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 305 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			1 + b_dim1], ldb, (ftnlen)19);
#line 307 "chetrs.f"
		clacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 308 "chetrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 312 "chetrs.f"
	    kp = -ipiv[k];
#line 313 "chetrs.f"
	    if (kp != k) {
#line 313 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 313 "chetrs.f"
	    }
#line 315 "chetrs.f"
	    k += 2;
#line 316 "chetrs.f"
	}

#line 318 "chetrs.f"
	goto L40;
#line 319 "chetrs.f"
L50:

#line 321 "chetrs.f"
	;
#line 321 "chetrs.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 330 "chetrs.f"
	k = 1;
#line 331 "chetrs.f"
L60:

/*        If K > N, exit from loop. */

#line 335 "chetrs.f"
	if (k > *n) {
#line 335 "chetrs.f"
	    goto L80;
#line 335 "chetrs.f"
	}

#line 338 "chetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 344 "chetrs.f"
	    kp = ipiv[k];
#line 345 "chetrs.f"
	    if (kp != k) {
#line 345 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 345 "chetrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 351 "chetrs.f"
	    if (k < *n) {
#line 351 "chetrs.f"
		i__1 = *n - k;
#line 351 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 351 "chetrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 351 "chetrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 357 "chetrs.f"
	    i__1 = k + k * a_dim1;
#line 357 "chetrs.f"
	    s = 1. / a[i__1].r;
#line 358 "chetrs.f"
	    csscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 359 "chetrs.f"
	    ++k;
#line 360 "chetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

#line 366 "chetrs.f"
	    kp = -ipiv[k];
#line 367 "chetrs.f"
	    if (kp != k + 1) {
#line 367 "chetrs.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 367 "chetrs.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 373 "chetrs.f"
	    if (k < *n - 1) {
#line 374 "chetrs.f"
		i__1 = *n - k - 1;
#line 374 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 374 "chetrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 376 "chetrs.f"
		i__1 = *n - k - 1;
#line 376 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 376 "chetrs.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 378 "chetrs.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 382 "chetrs.f"
	    i__1 = k + 1 + k * a_dim1;
#line 382 "chetrs.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 383 "chetrs.f"
	    d_cnjg(&z__2, &akm1k);
#line 383 "chetrs.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 383 "chetrs.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 384 "chetrs.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 384 "chetrs.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 385 "chetrs.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 385 "chetrs.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 385 "chetrs.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 386 "chetrs.f"
	    i__1 = *nrhs;
#line 386 "chetrs.f"
	    for (j = 1; j <= i__1; ++j) {
#line 387 "chetrs.f"
		d_cnjg(&z__2, &akm1k);
#line 387 "chetrs.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 387 "chetrs.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 388 "chetrs.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 388 "chetrs.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 389 "chetrs.f"
		i__2 = k + j * b_dim1;
#line 389 "chetrs.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 389 "chetrs.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 389 "chetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 389 "chetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 390 "chetrs.f"
		i__2 = k + 1 + j * b_dim1;
#line 390 "chetrs.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 390 "chetrs.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 390 "chetrs.f"
		z_div(&z__1, &z__2, &denom);
#line 390 "chetrs.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 391 "chetrs.f"
/* L70: */
#line 391 "chetrs.f"
	    }
#line 392 "chetrs.f"
	    k += 2;
#line 393 "chetrs.f"
	}

#line 395 "chetrs.f"
	goto L60;
#line 396 "chetrs.f"
L80:

/*        Next solve L**H *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 403 "chetrs.f"
	k = *n;
#line 404 "chetrs.f"
L90:

/*        If K < 1, exit from loop. */

#line 408 "chetrs.f"
	if (k < 1) {
#line 408 "chetrs.f"
	    goto L100;
#line 408 "chetrs.f"
	}

#line 411 "chetrs.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**H(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 418 "chetrs.f"
	    if (k < *n) {
#line 419 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 420 "chetrs.f"
		i__1 = *n - k;
#line 420 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 420 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 423 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 424 "chetrs.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 428 "chetrs.f"
	    kp = ipiv[k];
#line 429 "chetrs.f"
	    if (kp != k) {
#line 429 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 429 "chetrs.f"
	    }
#line 431 "chetrs.f"
	    --k;
#line 432 "chetrs.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 439 "chetrs.f"
	    if (k < *n) {
#line 440 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);
#line 441 "chetrs.f"
		i__1 = *n - k;
#line 441 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 441 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 444 "chetrs.f"
		clacgv_(nrhs, &b[k + b_dim1], ldb);

#line 446 "chetrs.f"
		clacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 447 "chetrs.f"
		i__1 = *n - k;
#line 447 "chetrs.f"
		z__1.r = -1., z__1.i = -0.;
#line 447 "chetrs.f"
		cgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &
			c_b1, &b[k - 1 + b_dim1], ldb, (ftnlen)19);
#line 450 "chetrs.f"
		clacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 451 "chetrs.f"
	    }

/*           Interchange rows K and -IPIV(K). */

#line 455 "chetrs.f"
	    kp = -ipiv[k];
#line 456 "chetrs.f"
	    if (kp != k) {
#line 456 "chetrs.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 456 "chetrs.f"
	    }
#line 458 "chetrs.f"
	    k += -2;
#line 459 "chetrs.f"
	}

#line 461 "chetrs.f"
	goto L90;
#line 462 "chetrs.f"
L100:
#line 463 "chetrs.f"
	;
#line 463 "chetrs.f"
    }

#line 465 "chetrs.f"
    return 0;

/*     End of CHETRS */

} /* chetrs_ */

