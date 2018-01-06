#line 1 "csytrs_rook.f"
/* csytrs_rook.f -- translated by f2c (version 20100827).
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

#line 1 "csytrs_rook.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSYTRS_ROOK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRS_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs_
rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs_
rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs_
rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */

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
/* > CSYTRS_ROOK solves a system of linear equations A*X = B with */
/* > a complex symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSYTRF_ROOK. */
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
/* >          obtain the factor U or L as computed by CSYTRF_ROOK. */
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
/* >          as determined by CSYTRF_ROOK. */
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

/* > \ingroup complexSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >   November 2011, Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int csytrs_rook__(char *uplo, integer *n, integer *nrhs, 
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

#line 176 "csytrs_rook.f"
    /* Parameter adjustments */
#line 176 "csytrs_rook.f"
    a_dim1 = *lda;
#line 176 "csytrs_rook.f"
    a_offset = 1 + a_dim1;
#line 176 "csytrs_rook.f"
    a -= a_offset;
#line 176 "csytrs_rook.f"
    --ipiv;
#line 176 "csytrs_rook.f"
    b_dim1 = *ldb;
#line 176 "csytrs_rook.f"
    b_offset = 1 + b_dim1;
#line 176 "csytrs_rook.f"
    b -= b_offset;
#line 176 "csytrs_rook.f"

#line 176 "csytrs_rook.f"
    /* Function Body */
#line 176 "csytrs_rook.f"
    *info = 0;
#line 177 "csytrs_rook.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 178 "csytrs_rook.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 179 "csytrs_rook.f"
	*info = -1;
#line 180 "csytrs_rook.f"
    } else if (*n < 0) {
#line 181 "csytrs_rook.f"
	*info = -2;
#line 182 "csytrs_rook.f"
    } else if (*nrhs < 0) {
#line 183 "csytrs_rook.f"
	*info = -3;
#line 184 "csytrs_rook.f"
    } else if (*lda < max(1,*n)) {
#line 185 "csytrs_rook.f"
	*info = -5;
#line 186 "csytrs_rook.f"
    } else if (*ldb < max(1,*n)) {
#line 187 "csytrs_rook.f"
	*info = -8;
#line 188 "csytrs_rook.f"
    }
#line 189 "csytrs_rook.f"
    if (*info != 0) {
#line 190 "csytrs_rook.f"
	i__1 = -(*info);
#line 190 "csytrs_rook.f"
	xerbla_("CSYTRS_ROOK", &i__1, (ftnlen)11);
#line 191 "csytrs_rook.f"
	return 0;
#line 192 "csytrs_rook.f"
    }

/*     Quick return if possible */

#line 196 "csytrs_rook.f"
    if (*n == 0 || *nrhs == 0) {
#line 196 "csytrs_rook.f"
	return 0;
#line 196 "csytrs_rook.f"
    }

#line 199 "csytrs_rook.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 208 "csytrs_rook.f"
	k = *n;
#line 209 "csytrs_rook.f"
L10:

/*        If K < 1, exit from loop. */

#line 213 "csytrs_rook.f"
	if (k < 1) {
#line 213 "csytrs_rook.f"
	    goto L30;
#line 213 "csytrs_rook.f"
	}

#line 216 "csytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 222 "csytrs_rook.f"
	    kp = ipiv[k];
#line 223 "csytrs_rook.f"
	    if (kp != k) {
#line 223 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 223 "csytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 229 "csytrs_rook.f"
	    i__1 = k - 1;
#line 229 "csytrs_rook.f"
	    z__1.r = -1., z__1.i = -0.;
#line 229 "csytrs_rook.f"
	    cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

#line 234 "csytrs_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 234 "csytrs_rook.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 235 "csytrs_rook.f"
	    --k;
#line 236 "csytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 242 "csytrs_rook.f"
	    kp = -ipiv[k];
#line 243 "csytrs_rook.f"
	    if (kp != k) {
#line 243 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 243 "csytrs_rook.f"
	    }

#line 246 "csytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 247 "csytrs_rook.f"
	    if (kp != k - 1) {
#line 247 "csytrs_rook.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 247 "csytrs_rook.f"
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 253 "csytrs_rook.f"
	    if (k > 2) {
#line 254 "csytrs_rook.f"
		i__1 = k - 2;
#line 254 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 254 "csytrs_rook.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + 
			b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 256 "csytrs_rook.f"
		i__1 = k - 2;
#line 256 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 256 "csytrs_rook.f"
		cgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &
			b[k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
#line 258 "csytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 262 "csytrs_rook.f"
	    i__1 = k - 1 + k * a_dim1;
#line 262 "csytrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 263 "csytrs_rook.f"
	    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
#line 263 "csytrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 264 "csytrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 264 "csytrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 265 "csytrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 265 "csytrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 265 "csytrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 266 "csytrs_rook.f"
	    i__1 = *nrhs;
#line 266 "csytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 267 "csytrs_rook.f"
		z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
#line 267 "csytrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 268 "csytrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 268 "csytrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 269 "csytrs_rook.f"
		i__2 = k - 1 + j * b_dim1;
#line 269 "csytrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 269 "csytrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 269 "csytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 269 "csytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 270 "csytrs_rook.f"
		i__2 = k + j * b_dim1;
#line 270 "csytrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 270 "csytrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 270 "csytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 270 "csytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 271 "csytrs_rook.f"
/* L20: */
#line 271 "csytrs_rook.f"
	    }
#line 272 "csytrs_rook.f"
	    k += -2;
#line 273 "csytrs_rook.f"
	}

#line 275 "csytrs_rook.f"
	goto L10;
#line 276 "csytrs_rook.f"
L30:

/*        Next solve U**T *X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 283 "csytrs_rook.f"
	k = 1;
#line 284 "csytrs_rook.f"
L40:

/*        If K > N, exit from loop. */

#line 288 "csytrs_rook.f"
	if (k > *n) {
#line 288 "csytrs_rook.f"
	    goto L50;
#line 288 "csytrs_rook.f"
	}

#line 291 "csytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

#line 298 "csytrs_rook.f"
	    if (k > 1) {
#line 298 "csytrs_rook.f"
		i__1 = k - 1;
#line 298 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 298 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 298 "csytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 304 "csytrs_rook.f"
	    kp = ipiv[k];
#line 305 "csytrs_rook.f"
	    if (kp != k) {
#line 305 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 305 "csytrs_rook.f"
	    }
#line 307 "csytrs_rook.f"
	    ++k;
#line 308 "csytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 315 "csytrs_rook.f"
	    if (k > 1) {
#line 316 "csytrs_rook.f"
		i__1 = k - 1;
#line 316 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 316 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb, (
			ftnlen)9);
#line 318 "csytrs_rook.f"
		i__1 = k - 1;
#line 318 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 318 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[
			(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + 
			b_dim1], ldb, (ftnlen)9);
#line 320 "csytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1). */

#line 324 "csytrs_rook.f"
	    kp = -ipiv[k];
#line 325 "csytrs_rook.f"
	    if (kp != k) {
#line 325 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 325 "csytrs_rook.f"
	    }

#line 328 "csytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 329 "csytrs_rook.f"
	    if (kp != k + 1) {
#line 329 "csytrs_rook.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 329 "csytrs_rook.f"
	    }

#line 332 "csytrs_rook.f"
	    k += 2;
#line 333 "csytrs_rook.f"
	}

#line 335 "csytrs_rook.f"
	goto L40;
#line 336 "csytrs_rook.f"
L50:

#line 338 "csytrs_rook.f"
	;
#line 338 "csytrs_rook.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 347 "csytrs_rook.f"
	k = 1;
#line 348 "csytrs_rook.f"
L60:

/*        If K > N, exit from loop. */

#line 352 "csytrs_rook.f"
	if (k > *n) {
#line 352 "csytrs_rook.f"
	    goto L80;
#line 352 "csytrs_rook.f"
	}

#line 355 "csytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

#line 361 "csytrs_rook.f"
	    kp = ipiv[k];
#line 362 "csytrs_rook.f"
	    if (kp != k) {
#line 362 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 362 "csytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 368 "csytrs_rook.f"
	    if (k < *n) {
#line 368 "csytrs_rook.f"
		i__1 = *n - k;
#line 368 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 368 "csytrs_rook.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
#line 368 "csytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 374 "csytrs_rook.f"
	    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
#line 374 "csytrs_rook.f"
	    cscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
#line 375 "csytrs_rook.f"
	    ++k;
#line 376 "csytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1) */

#line 382 "csytrs_rook.f"
	    kp = -ipiv[k];
#line 383 "csytrs_rook.f"
	    if (kp != k) {
#line 383 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 383 "csytrs_rook.f"
	    }

#line 386 "csytrs_rook.f"
	    kp = -ipiv[k + 1];
#line 387 "csytrs_rook.f"
	    if (kp != k + 1) {
#line 387 "csytrs_rook.f"
		cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 387 "csytrs_rook.f"
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

#line 393 "csytrs_rook.f"
	    if (k < *n - 1) {
#line 394 "csytrs_rook.f"
		i__1 = *n - k - 1;
#line 394 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 394 "csytrs_rook.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[
			k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
#line 396 "csytrs_rook.f"
		i__1 = *n - k - 1;
#line 396 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 396 "csytrs_rook.f"
		cgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], &
			c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], 
			ldb);
#line 398 "csytrs_rook.f"
	    }

/*           Multiply by the inverse of the diagonal block. */

#line 402 "csytrs_rook.f"
	    i__1 = k + 1 + k * a_dim1;
#line 402 "csytrs_rook.f"
	    akm1k.r = a[i__1].r, akm1k.i = a[i__1].i;
#line 403 "csytrs_rook.f"
	    z_div(&z__1, &a[k + k * a_dim1], &akm1k);
#line 403 "csytrs_rook.f"
	    akm1.r = z__1.r, akm1.i = z__1.i;
#line 404 "csytrs_rook.f"
	    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
#line 404 "csytrs_rook.f"
	    ak.r = z__1.r, ak.i = z__1.i;
#line 405 "csytrs_rook.f"
	    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * ak.i + 
		    akm1.i * ak.r;
#line 405 "csytrs_rook.f"
	    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 405 "csytrs_rook.f"
	    denom.r = z__1.r, denom.i = z__1.i;
#line 406 "csytrs_rook.f"
	    i__1 = *nrhs;
#line 406 "csytrs_rook.f"
	    for (j = 1; j <= i__1; ++j) {
#line 407 "csytrs_rook.f"
		z_div(&z__1, &b[k + j * b_dim1], &akm1k);
#line 407 "csytrs_rook.f"
		bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 408 "csytrs_rook.f"
		z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
#line 408 "csytrs_rook.f"
		bk.r = z__1.r, bk.i = z__1.i;
#line 409 "csytrs_rook.f"
		i__2 = k + j * b_dim1;
#line 409 "csytrs_rook.f"
		z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			bkm1.i + ak.i * bkm1.r;
#line 409 "csytrs_rook.f"
		z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 409 "csytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 409 "csytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 410 "csytrs_rook.f"
		i__2 = k + 1 + j * b_dim1;
#line 410 "csytrs_rook.f"
		z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			bk.i + akm1.i * bk.r;
#line 410 "csytrs_rook.f"
		z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 410 "csytrs_rook.f"
		z_div(&z__1, &z__2, &denom);
#line 410 "csytrs_rook.f"
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 411 "csytrs_rook.f"
/* L70: */
#line 411 "csytrs_rook.f"
	    }
#line 412 "csytrs_rook.f"
	    k += 2;
#line 413 "csytrs_rook.f"
	}

#line 415 "csytrs_rook.f"
	goto L60;
#line 416 "csytrs_rook.f"
L80:

/*        Next solve L**T *X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

#line 423 "csytrs_rook.f"
	k = *n;
#line 424 "csytrs_rook.f"
L90:

/*        If K < 1, exit from loop. */

#line 428 "csytrs_rook.f"
	if (k < 1) {
#line 428 "csytrs_rook.f"
	    goto L100;
#line 428 "csytrs_rook.f"
	}

#line 431 "csytrs_rook.f"
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

#line 438 "csytrs_rook.f"
	    if (k < *n) {
#line 438 "csytrs_rook.f"
		i__1 = *n - k;
#line 438 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 438 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 438 "csytrs_rook.f"
	    }

/*           Interchange rows K and IPIV(K). */

#line 444 "csytrs_rook.f"
	    kp = ipiv[k];
#line 445 "csytrs_rook.f"
	    if (kp != k) {
#line 445 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 445 "csytrs_rook.f"
	    }
#line 447 "csytrs_rook.f"
	    --k;
#line 448 "csytrs_rook.f"
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

#line 455 "csytrs_rook.f"
	    if (k < *n) {
#line 456 "csytrs_rook.f"
		i__1 = *n - k;
#line 456 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 456 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + 
			b_dim1], ldb, (ftnlen)9);
#line 458 "csytrs_rook.f"
		i__1 = *n - k;
#line 458 "csytrs_rook.f"
		z__1.r = -1., z__1.i = -0.;
#line 458 "csytrs_rook.f"
		cgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], 
			ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k 
			- 1 + b_dim1], ldb, (ftnlen)9);
#line 461 "csytrs_rook.f"
	    }

/*           Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1) */

#line 465 "csytrs_rook.f"
	    kp = -ipiv[k];
#line 466 "csytrs_rook.f"
	    if (kp != k) {
#line 466 "csytrs_rook.f"
		cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 466 "csytrs_rook.f"
	    }

#line 469 "csytrs_rook.f"
	    kp = -ipiv[k - 1];
#line 470 "csytrs_rook.f"
	    if (kp != k - 1) {
#line 470 "csytrs_rook.f"
		cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 470 "csytrs_rook.f"
	    }

#line 473 "csytrs_rook.f"
	    k += -2;
#line 474 "csytrs_rook.f"
	}

#line 476 "csytrs_rook.f"
	goto L90;
#line 477 "csytrs_rook.f"
L100:
#line 478 "csytrs_rook.f"
	;
#line 478 "csytrs_rook.f"
    }

#line 480 "csytrs_rook.f"
    return 0;

/*     End of CSYTRS_ROOK */

} /* csytrs_rook__ */

