#line 1 "zhetrs_rook.f"
/* zhetrs_rook.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrs_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHETRS_ROOK computes the solution to a system of linear equations A * X = B for HE matrices usi
ng factorization obtained with one of the bounded diagonal pivoting methods (max 2 interchanges) */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRS_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > ZHETRS_ROOK solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHETRF_ROOK. */
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
/* >          obtain the factor U or L as computed by ZHETRF_ROOK. */
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
/* >          as determined by ZHETRF_ROOK. */
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

/* > \date November 2013 */

/* > \ingroup complex16HEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  November 2013,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int zhetrs_rook__(char *uplo, integer *n, integer *nrhs, 
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


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

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

#line 177 "zhetrs_rook.f"
    /* Parameter adjustments */
#line 177 "zhetrs_rook.f"
    a_dim1 = *lda;
#line 177 "zhetrs_rook.f"
    a_offset = 1 + a_dim1;
#line 177 "zhetrs_rook.f"
    a -= a_offset;
#line 177 "zhetrs_rook.f"
    --ipiv;
#line 177 "zhetrs_rook.f"
    b_dim1 = *ldb;
#line 177 "zhetrs_rook.f"
    b_offset = 1 + b_dim1;
#line 177 "zhetrs_rook.f"
    b -= b_offset;
#line 177 "zhetrs_rook.f"

#line 177 "zhetrs_rook.f"
    /* Function Body */
#line 177 "zhetrs_rook.f"
    *info = 0;
#line 178 "zhetrs_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 179 "zhetrs_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 180 "zhetrs_rook.f"
	*info = -1;
#line 181 "zhetrs_rook.f"
    } else if (*n < 0) {
#line 182 "zhetrs_rook.f"
	*info = -2;
#line 183 "zhetrs_rook.f"
    } else if (*nrhs < 0) {
#line 184 "zhetrs_rook.f"
	*info = -3;
#line 185 "zhetrs_rook.f"
    } else if (*lda < max(1,*n)) {
#line 186 "zhetrs_rook.f"
	*info = -5;
#line 187 "zhetrs_rook.f"
    } else if (*ldb < max(1,*n)) {
#line 188 "zhetrs_rook.f"
	*info = -8;
#line 189 "zhetrs_rook.f"
    }
#line 190 "zhetrs_rook.f"
    if (*info != 0) {
#line 191 "zhetrs_rook.f"
	i__1 = -(*info);
#line 191 "zhetrs_rook.f"
	xerbla_("ZHETRS_ROOK", &i__1, (ftnlen)11);
#line 192 "zhetrs_rook.f"
	return 0;
#line 193 "zhetrs_rook.f"
    }

/*     Quick return if possible */

#line 197 "zhetrs_rook.f"
    if (*n == 0 || *nrhs == 0) {
#line 197 "zhetrs_rook.f"
	return 0;
#line 197 "zhetrs_rook.f"
    }

#line 200 "zhetrs_rook.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 209 "zhetrs_rook.f"
	k = *n;
#line 210 "zhetrs_rook.f"
L10:

/*        If K < 1, exit from loop. */

#line 214 "zhetrs_rook.f"
	if (k < 1) {
#line 214 "zhetrs_rook.f"
	    goto L30;
#line 214 "zhetrs_rook.f"
	}

#line 217 "zhetrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 223 "zhetrs_rook.f"
	    kp = ipiv[k];
#line 224 "zhetrs_rook.f"
	    if (kp != k) {
#line 224 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 224 "zhetrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 230 "zhetrs_rook.f"
	    i__1 = k - 1;
#line 230 "zhetrs_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 230 "zhetrs_rook.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 235 "zhetrs_rook.f"
	    i__1 = k + k * a_dim1;
#line 235 "zhetrs_rook.f"
	    s = 1. / a[i__1].r;
#line 236 "zhetrs_rook.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 237 "zhetrs_rook.f"
	    --k;
#line 238 "zhetrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1) */

#line 244 "zhetrs_rook.f"
	    kp = -ipiv[k];
#line 245 "zhetrs_rook.f"
	    if (kp != k) {
#line 245 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 245 "zhetrs_rook.f"
	    }

#line 248 "zhetrs_rook.f"
	    kp = -ipiv[k - 1];
#line 249 "zhetrs_rook.f"
	    if (kp != k - 1) {
#line 249 "zhetrs_rook.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 249 "zhetrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 255 "zhetrs_rook.f"
	    i__1 = k - 2;
#line 255 "zhetrs_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 255 "zhetrs_rook.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 257 "zhetrs_rook.f"
	    i__1 = k - 2;
#line 257 "zhetrs_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 257 "zhetrs_rook.f"
	    zgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k 
		    - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 262 "zhetrs_rook.f"
	    i__1 = k - 1 + k * a_dim1;
#line 262 "zhetrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 263 "zhetrs_rook.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 263 "zhetrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 264 "zhetrs_rook.f"
	    d_cnjg(&z__2, &akm1k);
#line 264 "zhetrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 264 "zhetrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 265 "zhetrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 265 "zhetrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 265 "zhetrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 266 "zhetrs_rook.f"
	    i__1 = *nrhs;
#line 266 "zhetrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 267 "zhetrs_rook.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 267 "zhetrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 268 "zhetrs_rook.f"
		d_cnjg(&z__2, &akm1k);
#line 268 "zhetrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 268 "zhetrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 269 "zhetrs_rook.f"
		i__2 = k - 1 + j * b_dim1;
#line 269 "zhetrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 269 "zhetrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 269 "zhetrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 269 "zhetrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 270 "zhetrs_rook.f"
		i__2 = k + j * b_dim1;
#line 270 "zhetrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 270 "zhetrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 270 "zhetrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 270 "zhetrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 271 "zhetrs_rook.f"
/* L20: */
#line 271 "zhetrs_rook.f"
	    }
#line 272 "zhetrs_rook.f"
	    k += -2;
#line 273 "zhetrs_rook.f"
	}

#line 275 "zhetrs_rook.f"
	goto L10;
#line 276 "zhetrs_rook.f"
L30:

/*        Next solve U**H *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 283 "zhetrs_rook.f"
	k = 1;
#line 284 "zhetrs_rook.f"
L40:

/*        If K > N, exit from loop. */

#line 288 "zhetrs_rook.f"
	if (k > *n) {
#line 288 "zhetrs_rook.f"
	    goto L50;
#line 288 "zhetrs_rook.f"
	}

#line 291 "zhetrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**H(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 298 "zhetrs_rook.f"
	    if (k > 1) {
#line 299 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 300 "zhetrs_rook.f"
		i__1 = k - 1;
#line 300 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 300 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 302 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 303 "zhetrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 307 "zhetrs_rook.f"
	    kp = ipiv[k];
#line 308 "zhetrs_rook.f"
	    if (kp != k) {
#line 308 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 308 "zhetrs_rook.f"
	    }
#line 310 "zhetrs_rook.f"
	    ++k;
#line 311 "zhetrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**H(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 318 "zhetrs_rook.f"
	    if (k > 1) {
#line 319 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 320 "zhetrs_rook.f"
		i__1 = k - 1;
#line 320 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 320 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)19);
#line 322 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 324 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 325 "zhetrs_rook.f"
		i__1 = k - 1;
#line 325 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 325 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[b_offset]
			, ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 
			1 + b_dim1], ldb, (ftnlen)19);
#line 327 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + 1 + b_dim1], ldb);
#line 328 "zhetrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1) */

#line 332 "zhetrs_rook.f"
	    kp = -ipiv[k];
#line 333 "zhetrs_rook.f"
	    if (kp != k) {
#line 333 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 333 "zhetrs_rook.f"
	    }

#line 336 "zhetrs_rook.f"
	    kp = -ipiv[k + 1];
#line 337 "zhetrs_rook.f"
	    if (kp != k + 1) {
#line 337 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 337 "zhetrs_rook.f"
	    }

#line 340 "zhetrs_rook.f"
	    k += 2;
#line 341 "zhetrs_rook.f"
	}

#line 343 "zhetrs_rook.f"
	goto L40;
#line 344 "zhetrs_rook.f"
L50:

#line 346 "zhetrs_rook.f"
	;
#line 346 "zhetrs_rook.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 355 "zhetrs_rook.f"
	k = 1;
#line 356 "zhetrs_rook.f"
L60:

/*        If K > N, exit from loop. */

#line 360 "zhetrs_rook.f"
	if (k > *n) {
#line 360 "zhetrs_rook.f"
	    goto L80;
#line 360 "zhetrs_rook.f"
	}

#line 363 "zhetrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 369 "zhetrs_rook.f"
	    kp = ipiv[k];
#line 370 "zhetrs_rook.f"
	    if (kp != k) {
#line 370 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 370 "zhetrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 376 "zhetrs_rook.f"
	    if (k < *n) {
#line 376 "zhetrs_rook.f"
		i__1 = *n - k;
#line 376 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 376 "zhetrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 376 "zhetrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 382 "zhetrs_rook.f"
	    i__1 = k + k * a_dim1;
#line 382 "zhetrs_rook.f"
	    s = 1. / a[i__1].r;
#line 383 "zhetrs_rook.f"
	    zdscal_(nrhs, &s, &b[k + b_dim1], ldb);
#line 384 "zhetrs_rook.f"
	    ++k;
#line 385 "zhetrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K), then K+1 and -IPIV(K+1) */

#line 391 "zhetrs_rook.f"
	    kp = -ipiv[k];
#line 392 "zhetrs_rook.f"
	    if (kp != k) {
#line 392 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 392 "zhetrs_rook.f"
	    }

#line 395 "zhetrs_rook.f"
	    kp = -ipiv[k + 1];
#line 396 "zhetrs_rook.f"
	    if (kp != k + 1) {
#line 396 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 396 "zhetrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 402 "zhetrs_rook.f"
	    if (k < *n - 1) {
#line 403 "zhetrs_rook.f"
		i__1 = *n - k - 1;
#line 403 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 403 "zhetrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 405 "zhetrs_rook.f"
		i__1 = *n - k - 1;
#line 405 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 405 "zhetrs_rook.f"
		zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 407 "zhetrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 411 "zhetrs_rook.f"
	    i__1 = k + 1 + k * a_dim1;
#line 411 "zhetrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 412 "zhetrs_rook.f"
	    d_cnjg(&z__2, &akm1k);
#line 412 "zhetrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &z__2);
#line 412 "zhetrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 413 "zhetrs_rook.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 413 "zhetrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 414 "zhetrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 414 "zhetrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 414 "zhetrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 415 "zhetrs_rook.f"
	    i__1 = *nrhs;
#line 415 "zhetrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 416 "zhetrs_rook.f"
		d_cnjg(&z__2, &akm1k);
#line 416 "zhetrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &z__2);
#line 416 "zhetrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 417 "zhetrs_rook.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 417 "zhetrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 418 "zhetrs_rook.f"
		i__2 = k + j * b_dim1;
#line 418 "zhetrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 418 "zhetrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 418 "zhetrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 418 "zhetrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 419 "zhetrs_rook.f"
		i__2 = k + 1 + j * b_dim1;
#line 419 "zhetrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 419 "zhetrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 419 "zhetrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 419 "zhetrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 420 "zhetrs_rook.f"
/* L70: */
#line 420 "zhetrs_rook.f"
	    }
#line 421 "zhetrs_rook.f"
	    k += 2;
#line 422 "zhetrs_rook.f"
	}

#line 424 "zhetrs_rook.f"
	goto L60;
#line 425 "zhetrs_rook.f"
L80:

/*        Next solve L**H *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 432 "zhetrs_rook.f"
	k = *n;
#line 433 "zhetrs_rook.f"
L90:

/*        If K < 1, exit from loop. */

#line 437 "zhetrs_rook.f"
	if (k < 1) {
#line 437 "zhetrs_rook.f"
	    goto L100;
#line 437 "zhetrs_rook.f"
	}

#line 440 "zhetrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**H(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 447 "zhetrs_rook.f"
	    if (k < *n) {
#line 448 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 449 "zhetrs_rook.f"
		i__1 = *n - k;
#line 449 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 449 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 452 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 453 "zhetrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 457 "zhetrs_rook.f"
	    kp = ipiv[k];
#line 458 "zhetrs_rook.f"
	    if (kp != k) {
#line 458 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 458 "zhetrs_rook.f"
	    }
#line 460 "zhetrs_rook.f"
	    --k;
#line 461 "zhetrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**H(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 468 "zhetrs_rook.f"
	    if (k < *n) {
#line 469 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);
#line 470 "zhetrs_rook.f"
		i__1 = *n - k;
#line 470 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 470 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &
			b[k + b_dim1], ldb, (ftnlen)19);
#line 473 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k + b_dim1], ldb);

#line 475 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 476 "zhetrs_rook.f"
		i__1 = *n - k;
#line 476 "zhetrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 476 "zhetrs_rook.f"
		zgemv_("Conjugate transpose", &i__1, nrhs, &z__1, &b[k + 1 + 
			b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &
			c_b1, &b[k - 1 + b_dim1], ldb, (ftnlen)19);
#line 479 "zhetrs_rook.f"
		zlacgv_(nrhs, &b[k - 1 + b_dim1], ldb);
#line 480 "zhetrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K), then K-1 and -IPIV(K-1) */

#line 484 "zhetrs_rook.f"
	    kp = -ipiv[k];
#line 485 "zhetrs_rook.f"
	    if (kp != k) {
#line 485 "zhetrs_rook.f"
		zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 485 "zhetrs_rook.f"
	    }

#line 488 "zhetrs_rook.f"
	    kp = -ipiv[k - 1];
#line 489 "zhetrs_rook.f"
	    if (kp != k - 1) {
#line 489 "zhetrs_rook.f"
		zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 489 "zhetrs_rook.f"
	    }

#line 492 "zhetrs_rook.f"
	    k += -2;
#line 493 "zhetrs_rook.f"
	}

#line 495 "zhetrs_rook.f"
	goto L90;
#line 496 "zhetrs_rook.f"
L100:
#line 497 "zhetrs_rook.f"
	;
#line 497 "zhetrs_rook.f"
    }

#line 499 "zhetrs_rook.f"
    return 0;

/*     End of ZHETRS_ROOK */

} /* zhetrs_rook__ */

