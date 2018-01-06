#line 1 "zhetrs2.f"
/* zhetrs2.f -- translated by f2c (version 20100827).
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

#line 1 "zhetrs2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b ZHETRS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                           WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16       A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRS2 solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHETRF and converted by ZSYCONV. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
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
/* Subroutine */ int zhetrs2_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), zdscal_(integer *, doublereal 
	    *, doublecomplex *, integer *), zsyconv_(char *, char *, integer *
	    , doublecomplex *, integer *, integer *, doublecomplex *, integer 
	    *, ftnlen, ftnlen);


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

#line 168 "zhetrs2.f"
    /* Parameter adjustments */
#line 168 "zhetrs2.f"
    a_dim1 = *lda;
#line 168 "zhetrs2.f"
    a_offset = 1 + a_dim1;
#line 168 "zhetrs2.f"
    a -= a_offset;
#line 168 "zhetrs2.f"
    --ipiv;
#line 168 "zhetrs2.f"
    b_dim1 = *ldb;
#line 168 "zhetrs2.f"
    b_offset = 1 + b_dim1;
#line 168 "zhetrs2.f"
    b -= b_offset;
#line 168 "zhetrs2.f"
    --work;
#line 168 "zhetrs2.f"

#line 168 "zhetrs2.f"
    /* Function Body */
#line 168 "zhetrs2.f"
    *info = 0;
#line 169 "zhetrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "zhetrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 171 "zhetrs2.f"
	*info = -1;
#line 172 "zhetrs2.f"
    } else if (*n < 0) {
#line 173 "zhetrs2.f"
	*info = -2;
#line 174 "zhetrs2.f"
    } else if (*nrhs < 0) {
#line 175 "zhetrs2.f"
	*info = -3;
#line 176 "zhetrs2.f"
    } else if (*lda < max(1,*n)) {
#line 177 "zhetrs2.f"
	*info = -5;
#line 178 "zhetrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 179 "zhetrs2.f"
	*info = -8;
#line 180 "zhetrs2.f"
    }
#line 181 "zhetrs2.f"
    if (*info != 0) {
#line 182 "zhetrs2.f"
	i__1 = -(*info);
#line 182 "zhetrs2.f"
	xerbla_("ZHETRS2", &i__1, (ftnlen)7);
#line 183 "zhetrs2.f"
	return 0;
#line 184 "zhetrs2.f"
    }

/*     Quick return if possible */

#line 188 "zhetrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 188 "zhetrs2.f"
	return 0;
#line 188 "zhetrs2.f"
    }

/*     Convert A */

#line 193 "zhetrs2.f"
    zsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 195 "zhetrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*       P**T * B */
#line 200 "zhetrs2.f"
	k = *n;
#line 201 "zhetrs2.f"
	while(k >= 1) {
#line 202 "zhetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 205 "zhetrs2.f"
		kp = ipiv[k];
#line 206 "zhetrs2.f"
		if (kp != k) {
#line 206 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 206 "zhetrs2.f"
		}
#line 208 "zhetrs2.f"
		--k;
#line 209 "zhetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 212 "zhetrs2.f"
		kp = -ipiv[k];
#line 213 "zhetrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 213 "zhetrs2.f"
		    zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 213 "zhetrs2.f"
		}
#line 215 "zhetrs2.f"
		k += -2;
#line 216 "zhetrs2.f"
	    }
#line 217 "zhetrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 221 "zhetrs2.f"
	ztrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 225 "zhetrs2.f"
	i__ = *n;
#line 226 "zhetrs2.f"
	while(i__ >= 1) {
#line 227 "zhetrs2.f"
	    if (ipiv[i__] > 0) {
#line 228 "zhetrs2.f"
		i__1 = i__ + i__ * a_dim1;
#line 228 "zhetrs2.f"
		s = 1. / a[i__1].r;
#line 229 "zhetrs2.f"
		zdscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 230 "zhetrs2.f"
	    } else if (i__ > 1) {
#line 231 "zhetrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 232 "zhetrs2.f"
		    i__1 = i__;
#line 232 "zhetrs2.f"
		    akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 233 "zhetrs2.f"
		    z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 233 "zhetrs2.f"
		    akm1.r = z__1.r, akm1.i = z__1.i;
#line 234 "zhetrs2.f"
		    d_cnjg(&z__2, &akm1k);
#line 234 "zhetrs2.f"
		    z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 234 "zhetrs2.f"
		    ak.r = z__1.r, ak.i = z__1.i;
#line 235 "zhetrs2.f"
		    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			    ak.i + akm1.i * ak.r;
#line 235 "zhetrs2.f"
		    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 235 "zhetrs2.f"
		    denom.r = z__1.r, denom.i = z__1.i;
#line 236 "zhetrs2.f"
		    i__1 = *nrhs;
#line 236 "zhetrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 237 "zhetrs2.f"
			z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 237 "zhetrs2.f"
			bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 238 "zhetrs2.f"
			d_cnjg(&z__2, &akm1k);
#line 238 "zhetrs2.f"
			z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 238 "zhetrs2.f"
			bk.r = z__1.r, bk.i = z__1.i;
#line 239 "zhetrs2.f"
			i__2 = i__ - 1 + j * b_dim1;
#line 239 "zhetrs2.f"
			z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r 
				* bkm1.i + ak.i * bkm1.r;
#line 239 "zhetrs2.f"
			z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 239 "zhetrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 239 "zhetrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 240 "zhetrs2.f"
			i__2 = i__ + j * b_dim1;
#line 240 "zhetrs2.f"
			z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = 
				akm1.r * bk.i + akm1.i * bk.r;
#line 240 "zhetrs2.f"
			z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 240 "zhetrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 240 "zhetrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 241 "zhetrs2.f"
/* L15: */
#line 241 "zhetrs2.f"
		    }
#line 242 "zhetrs2.f"
		    --i__;
#line 243 "zhetrs2.f"
		}
#line 244 "zhetrs2.f"
	    }
#line 245 "zhetrs2.f"
	    --i__;
#line 246 "zhetrs2.f"
	}

/*      Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ] */

#line 250 "zhetrs2.f"
	ztrsm_("L", "U", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ] */

#line 254 "zhetrs2.f"
	k = 1;
#line 255 "zhetrs2.f"
	while(k <= *n) {
#line 256 "zhetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 259 "zhetrs2.f"
		kp = ipiv[k];
#line 260 "zhetrs2.f"
		if (kp != k) {
#line 260 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 260 "zhetrs2.f"
		}
#line 262 "zhetrs2.f"
		++k;
#line 263 "zhetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 266 "zhetrs2.f"
		kp = -ipiv[k];
#line 267 "zhetrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 267 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 267 "zhetrs2.f"
		}
#line 269 "zhetrs2.f"
		k += 2;
#line 270 "zhetrs2.f"
	    }
#line 271 "zhetrs2.f"
	}

#line 273 "zhetrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*       P**T * B */
#line 278 "zhetrs2.f"
	k = 1;
#line 279 "zhetrs2.f"
	while(k <= *n) {
#line 280 "zhetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 283 "zhetrs2.f"
		kp = ipiv[k];
#line 284 "zhetrs2.f"
		if (kp != k) {
#line 284 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 284 "zhetrs2.f"
		}
#line 286 "zhetrs2.f"
		++k;
#line 287 "zhetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 290 "zhetrs2.f"
		kp = -ipiv[k + 1];
#line 291 "zhetrs2.f"
		if (kp == -ipiv[k]) {
#line 291 "zhetrs2.f"
		    zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 291 "zhetrs2.f"
		}
#line 293 "zhetrs2.f"
		k += 2;
#line 294 "zhetrs2.f"
	    }
#line 295 "zhetrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 299 "zhetrs2.f"
	ztrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 303 "zhetrs2.f"
	i__ = 1;
#line 304 "zhetrs2.f"
	while(i__ <= *n) {
#line 305 "zhetrs2.f"
	    if (ipiv[i__] > 0) {
#line 306 "zhetrs2.f"
		i__1 = i__ + i__ * a_dim1;
#line 306 "zhetrs2.f"
		s = 1. / a[i__1].r;
#line 307 "zhetrs2.f"
		zdscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 308 "zhetrs2.f"
	    } else {
#line 309 "zhetrs2.f"
		i__1 = i__;
#line 309 "zhetrs2.f"
		akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 310 "zhetrs2.f"
		d_cnjg(&z__2, &akm1k);
#line 310 "zhetrs2.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 310 "zhetrs2.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 311 "zhetrs2.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 311 "zhetrs2.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 312 "zhetrs2.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 312 "zhetrs2.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 312 "zhetrs2.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 313 "zhetrs2.f"
		i__1 = *nrhs;
#line 313 "zhetrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "zhetrs2.f"
		    d_cnjg(&z__2, &akm1k);
#line 314 "zhetrs2.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 314 "zhetrs2.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 315 "zhetrs2.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 315 "zhetrs2.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 316 "zhetrs2.f"
		    i__2 = i__ + j * b_dim1;
#line 316 "zhetrs2.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 316 "zhetrs2.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 316 "zhetrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 316 "zhetrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 317 "zhetrs2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 317 "zhetrs2.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 317 "zhetrs2.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 317 "zhetrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 317 "zhetrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 318 "zhetrs2.f"
/* L25: */
#line 318 "zhetrs2.f"
		}
#line 319 "zhetrs2.f"
		++i__;
#line 320 "zhetrs2.f"
	    }
#line 321 "zhetrs2.f"
	    ++i__;
#line 322 "zhetrs2.f"
	}

/*  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ] */

#line 326 "zhetrs2.f"
	ztrsm_("L", "L", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ] */

#line 330 "zhetrs2.f"
	k = *n;
#line 331 "zhetrs2.f"
	while(k >= 1) {
#line 332 "zhetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 335 "zhetrs2.f"
		kp = ipiv[k];
#line 336 "zhetrs2.f"
		if (kp != k) {
#line 336 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 336 "zhetrs2.f"
		}
#line 338 "zhetrs2.f"
		--k;
#line 339 "zhetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 342 "zhetrs2.f"
		kp = -ipiv[k];
#line 343 "zhetrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 343 "zhetrs2.f"
		    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 343 "zhetrs2.f"
		}
#line 345 "zhetrs2.f"
		k += -2;
#line 346 "zhetrs2.f"
	    }
#line 347 "zhetrs2.f"
	}

#line 349 "zhetrs2.f"
    }

/*     Revert A */

#line 353 "zhetrs2.f"
    zsyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 355 "zhetrs2.f"
    return 0;

/*     End of ZHETRS2 */

} /* zhetrs2_ */

