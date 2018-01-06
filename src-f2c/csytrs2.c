#line 1 "csytrs2.f"
/* csytrs2.f -- translated by f2c (version 20100827).
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

#line 1 "csytrs2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CSYTRS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYTRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                           WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRS2 solves a system of linear equations A*X = B with a COMPLEX */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSYTRF and converted by CSYCONV. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by CSYTRF. */
/* >          Note that A is input / output. This might be counter-intuitive, */
/* >          and one may think that A is input only. A is input / output. This */
/* >          is because, at the start of the subroutine, we permute A in a */
/* >          "better" form and then we permute A back to its original form at */
/* >          the end. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
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

/* > \date November 2015 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csytrs2_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, akm1k;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex denom;
    static integer iinfo;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), csyconv_(
	    char *, char *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 172 "csytrs2.f"
    /* Parameter adjustments */
#line 172 "csytrs2.f"
    a_dim1 = *lda;
#line 172 "csytrs2.f"
    a_offset = 1 + a_dim1;
#line 172 "csytrs2.f"
    a -= a_offset;
#line 172 "csytrs2.f"
    --ipiv;
#line 172 "csytrs2.f"
    b_dim1 = *ldb;
#line 172 "csytrs2.f"
    b_offset = 1 + b_dim1;
#line 172 "csytrs2.f"
    b -= b_offset;
#line 172 "csytrs2.f"
    --work;
#line 172 "csytrs2.f"

#line 172 "csytrs2.f"
    /* Function Body */
#line 172 "csytrs2.f"
    *info = 0;
#line 173 "csytrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 174 "csytrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 175 "csytrs2.f"
	*info = -1;
#line 176 "csytrs2.f"
    } else if (*n < 0) {
#line 177 "csytrs2.f"
	*info = -2;
#line 178 "csytrs2.f"
    } else if (*nrhs < 0) {
#line 179 "csytrs2.f"
	*info = -3;
#line 180 "csytrs2.f"
    } else if (*lda < max(1,*n)) {
#line 181 "csytrs2.f"
	*info = -5;
#line 182 "csytrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 183 "csytrs2.f"
	*info = -8;
#line 184 "csytrs2.f"
    }
#line 185 "csytrs2.f"
    if (*info != 0) {
#line 186 "csytrs2.f"
	i__1 = -(*info);
#line 186 "csytrs2.f"
	xerbla_("CSYTRS2", &i__1, (ftnlen)7);
#line 187 "csytrs2.f"
	return 0;
#line 188 "csytrs2.f"
    }

/*     Quick return if possible */

#line 192 "csytrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 192 "csytrs2.f"
	return 0;
#line 192 "csytrs2.f"
    }

/*     Convert A */

#line 197 "csytrs2.f"
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 199 "csytrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*       P**T * B */
#line 204 "csytrs2.f"
	k = *n;
#line 205 "csytrs2.f"
	while(k >= 1) {
#line 206 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 209 "csytrs2.f"
		kp = ipiv[k];
#line 210 "csytrs2.f"
		if (kp != k) {
#line 210 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 210 "csytrs2.f"
		}
#line 212 "csytrs2.f"
		--k;
#line 213 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 216 "csytrs2.f"
		kp = -ipiv[k];
#line 217 "csytrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 217 "csytrs2.f"
		    cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 217 "csytrs2.f"
		}
#line 219 "csytrs2.f"
		k += -2;
#line 220 "csytrs2.f"
	    }
#line 221 "csytrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 225 "csytrs2.f"
	ctrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 229 "csytrs2.f"
	i__ = *n;
#line 230 "csytrs2.f"
	while(i__ >= 1) {
#line 231 "csytrs2.f"
	    if (ipiv[i__] > 0) {
#line 232 "csytrs2.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 232 "csytrs2.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 233 "csytrs2.f"
	    } else if (i__ > 1) {
#line 234 "csytrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 235 "csytrs2.f"
		    i__1 = i__;
#line 235 "csytrs2.f"
		    akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 236 "csytrs2.f"
		    z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 236 "csytrs2.f"
		    akm1.r = z__1.r, akm1.i = z__1.i;
#line 237 "csytrs2.f"
		    z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 237 "csytrs2.f"
		    ak.r = z__1.r, ak.i = z__1.i;
#line 238 "csytrs2.f"
		    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			    ak.i + akm1.i * ak.r;
#line 238 "csytrs2.f"
		    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 238 "csytrs2.f"
		    denom.r = z__1.r, denom.i = z__1.i;
#line 239 "csytrs2.f"
		    i__1 = *nrhs;
#line 239 "csytrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 240 "csytrs2.f"
			z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 240 "csytrs2.f"
			bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 241 "csytrs2.f"
			z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 241 "csytrs2.f"
			bk.r = z__1.r, bk.i = z__1.i;
#line 242 "csytrs2.f"
			i__2 = i__ - 1 + j * b_dim1;
#line 242 "csytrs2.f"
			z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r 
				* bkm1.i + ak.i * bkm1.r;
#line 242 "csytrs2.f"
			z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 242 "csytrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 242 "csytrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 243 "csytrs2.f"
			i__2 = i__ + j * b_dim1;
#line 243 "csytrs2.f"
			z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = 
				akm1.r * bk.i + akm1.i * bk.r;
#line 243 "csytrs2.f"
			z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 243 "csytrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 243 "csytrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 244 "csytrs2.f"
/* L15: */
#line 244 "csytrs2.f"
		    }
#line 245 "csytrs2.f"
		    --i__;
#line 246 "csytrs2.f"
		}
#line 247 "csytrs2.f"
	    }
#line 248 "csytrs2.f"
	    --i__;
#line 249 "csytrs2.f"
	}

/*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 253 "csytrs2.f"
	ctrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

#line 257 "csytrs2.f"
	k = 1;
#line 258 "csytrs2.f"
	while(k <= *n) {
#line 259 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 262 "csytrs2.f"
		kp = ipiv[k];
#line 263 "csytrs2.f"
		if (kp != k) {
#line 263 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 263 "csytrs2.f"
		}
#line 265 "csytrs2.f"
		++k;
#line 266 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 269 "csytrs2.f"
		kp = -ipiv[k];
#line 270 "csytrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 270 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 270 "csytrs2.f"
		}
#line 272 "csytrs2.f"
		k += 2;
#line 273 "csytrs2.f"
	    }
#line 274 "csytrs2.f"
	}

#line 276 "csytrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*       P**T * B */
#line 281 "csytrs2.f"
	k = 1;
#line 282 "csytrs2.f"
	while(k <= *n) {
#line 283 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 286 "csytrs2.f"
		kp = ipiv[k];
#line 287 "csytrs2.f"
		if (kp != k) {
#line 287 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 287 "csytrs2.f"
		}
#line 289 "csytrs2.f"
		++k;
#line 290 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 293 "csytrs2.f"
		kp = -ipiv[k + 1];
#line 294 "csytrs2.f"
		if (kp == -ipiv[k]) {
#line 294 "csytrs2.f"
		    cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 294 "csytrs2.f"
		}
#line 296 "csytrs2.f"
		k += 2;
#line 297 "csytrs2.f"
	    }
#line 298 "csytrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 302 "csytrs2.f"
	ctrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 306 "csytrs2.f"
	i__ = 1;
#line 307 "csytrs2.f"
	while(i__ <= *n) {
#line 308 "csytrs2.f"
	    if (ipiv[i__] > 0) {
#line 309 "csytrs2.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 309 "csytrs2.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 310 "csytrs2.f"
	    } else {
#line 311 "csytrs2.f"
		i__1 = i__;
#line 311 "csytrs2.f"
		akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 312 "csytrs2.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 312 "csytrs2.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 313 "csytrs2.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 313 "csytrs2.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 314 "csytrs2.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 314 "csytrs2.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 314 "csytrs2.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 315 "csytrs2.f"
		i__1 = *nrhs;
#line 315 "csytrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 316 "csytrs2.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 316 "csytrs2.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 317 "csytrs2.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 317 "csytrs2.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 318 "csytrs2.f"
		    i__2 = i__ + j * b_dim1;
#line 318 "csytrs2.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 318 "csytrs2.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 318 "csytrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 318 "csytrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 319 "csytrs2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 319 "csytrs2.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 319 "csytrs2.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 319 "csytrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 319 "csytrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 320 "csytrs2.f"
/* L25: */
#line 320 "csytrs2.f"
		}
#line 321 "csytrs2.f"
		++i__;
#line 322 "csytrs2.f"
	    }
#line 323 "csytrs2.f"
	    ++i__;
#line 324 "csytrs2.f"
	}

/*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 328 "csytrs2.f"
	ctrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

#line 332 "csytrs2.f"
	k = *n;
#line 333 "csytrs2.f"
	while(k >= 1) {
#line 334 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 337 "csytrs2.f"
		kp = ipiv[k];
#line 338 "csytrs2.f"
		if (kp != k) {
#line 338 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 338 "csytrs2.f"
		}
#line 340 "csytrs2.f"
		--k;
#line 341 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 344 "csytrs2.f"
		kp = -ipiv[k];
#line 345 "csytrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 345 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 345 "csytrs2.f"
		}
#line 347 "csytrs2.f"
		k += -2;
#line 348 "csytrs2.f"
	    }
#line 349 "csytrs2.f"
	}

#line 351 "csytrs2.f"
    }

/*     Revert A */

#line 355 "csytrs2.f"
    csyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 357 "csytrs2.f"
    return 0;

/*     End of CSYTRS2 */

} /* csytrs2_ */

