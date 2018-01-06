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

/* > \date November 2011 */

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

#line 167 "csytrs2.f"
    /* Parameter adjustments */
#line 167 "csytrs2.f"
    a_dim1 = *lda;
#line 167 "csytrs2.f"
    a_offset = 1 + a_dim1;
#line 167 "csytrs2.f"
    a -= a_offset;
#line 167 "csytrs2.f"
    --ipiv;
#line 167 "csytrs2.f"
    b_dim1 = *ldb;
#line 167 "csytrs2.f"
    b_offset = 1 + b_dim1;
#line 167 "csytrs2.f"
    b -= b_offset;
#line 167 "csytrs2.f"
    --work;
#line 167 "csytrs2.f"

#line 167 "csytrs2.f"
    /* Function Body */
#line 167 "csytrs2.f"
    *info = 0;
#line 168 "csytrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 169 "csytrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 170 "csytrs2.f"
	*info = -1;
#line 171 "csytrs2.f"
    } else if (*n < 0) {
#line 172 "csytrs2.f"
	*info = -2;
#line 173 "csytrs2.f"
    } else if (*nrhs < 0) {
#line 174 "csytrs2.f"
	*info = -3;
#line 175 "csytrs2.f"
    } else if (*lda < max(1,*n)) {
#line 176 "csytrs2.f"
	*info = -5;
#line 177 "csytrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 178 "csytrs2.f"
	*info = -8;
#line 179 "csytrs2.f"
    }
#line 180 "csytrs2.f"
    if (*info != 0) {
#line 181 "csytrs2.f"
	i__1 = -(*info);
#line 181 "csytrs2.f"
	xerbla_("CSYTRS2", &i__1, (ftnlen)7);
#line 182 "csytrs2.f"
	return 0;
#line 183 "csytrs2.f"
    }

/*     Quick return if possible */

#line 187 "csytrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 187 "csytrs2.f"
	return 0;
#line 187 "csytrs2.f"
    }

/*     Convert A */

#line 192 "csytrs2.f"
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 194 "csytrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*       P**T * B */
#line 199 "csytrs2.f"
	k = *n;
#line 200 "csytrs2.f"
	while(k >= 1) {
#line 201 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 204 "csytrs2.f"
		kp = ipiv[k];
#line 205 "csytrs2.f"
		if (kp != k) {
#line 205 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 205 "csytrs2.f"
		}
#line 207 "csytrs2.f"
		--k;
#line 208 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 211 "csytrs2.f"
		kp = -ipiv[k];
#line 212 "csytrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 212 "csytrs2.f"
		    cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 212 "csytrs2.f"
		}
#line 214 "csytrs2.f"
		k += -2;
#line 215 "csytrs2.f"
	    }
#line 216 "csytrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 220 "csytrs2.f"
	ctrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 224 "csytrs2.f"
	i__ = *n;
#line 225 "csytrs2.f"
	while(i__ >= 1) {
#line 226 "csytrs2.f"
	    if (ipiv[i__] > 0) {
#line 227 "csytrs2.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 227 "csytrs2.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 228 "csytrs2.f"
	    } else if (i__ > 1) {
#line 229 "csytrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 230 "csytrs2.f"
		    i__1 = i__;
#line 230 "csytrs2.f"
		    akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 231 "csytrs2.f"
		    z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 231 "csytrs2.f"
		    akm1.r = z__1.r, akm1.i = z__1.i;
#line 232 "csytrs2.f"
		    z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 232 "csytrs2.f"
		    ak.r = z__1.r, ak.i = z__1.i;
#line 233 "csytrs2.f"
		    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			    ak.i + akm1.i * ak.r;
#line 233 "csytrs2.f"
		    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 233 "csytrs2.f"
		    denom.r = z__1.r, denom.i = z__1.i;
#line 234 "csytrs2.f"
		    i__1 = *nrhs;
#line 234 "csytrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 235 "csytrs2.f"
			z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 235 "csytrs2.f"
			bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 236 "csytrs2.f"
			z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 236 "csytrs2.f"
			bk.r = z__1.r, bk.i = z__1.i;
#line 237 "csytrs2.f"
			i__2 = i__ - 1 + j * b_dim1;
#line 237 "csytrs2.f"
			z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r 
				* bkm1.i + ak.i * bkm1.r;
#line 237 "csytrs2.f"
			z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 237 "csytrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 237 "csytrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 238 "csytrs2.f"
			i__2 = i__ + j * b_dim1;
#line 238 "csytrs2.f"
			z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = 
				akm1.r * bk.i + akm1.i * bk.r;
#line 238 "csytrs2.f"
			z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 238 "csytrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 238 "csytrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 239 "csytrs2.f"
/* L15: */
#line 239 "csytrs2.f"
		    }
#line 240 "csytrs2.f"
		    --i__;
#line 241 "csytrs2.f"
		}
#line 242 "csytrs2.f"
	    }
#line 243 "csytrs2.f"
	    --i__;
#line 244 "csytrs2.f"
	}

/*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 248 "csytrs2.f"
	ctrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

#line 252 "csytrs2.f"
	k = 1;
#line 253 "csytrs2.f"
	while(k <= *n) {
#line 254 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 257 "csytrs2.f"
		kp = ipiv[k];
#line 258 "csytrs2.f"
		if (kp != k) {
#line 258 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 258 "csytrs2.f"
		}
#line 260 "csytrs2.f"
		++k;
#line 261 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 264 "csytrs2.f"
		kp = -ipiv[k];
#line 265 "csytrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 265 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 265 "csytrs2.f"
		}
#line 267 "csytrs2.f"
		k += 2;
#line 268 "csytrs2.f"
	    }
#line 269 "csytrs2.f"
	}

#line 271 "csytrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*       P**T * B */
#line 276 "csytrs2.f"
	k = 1;
#line 277 "csytrs2.f"
	while(k <= *n) {
#line 278 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 281 "csytrs2.f"
		kp = ipiv[k];
#line 282 "csytrs2.f"
		if (kp != k) {
#line 282 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "csytrs2.f"
		}
#line 284 "csytrs2.f"
		++k;
#line 285 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 288 "csytrs2.f"
		kp = -ipiv[k + 1];
#line 289 "csytrs2.f"
		if (kp == -ipiv[k]) {
#line 289 "csytrs2.f"
		    cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 289 "csytrs2.f"
		}
#line 291 "csytrs2.f"
		k += 2;
#line 292 "csytrs2.f"
	    }
#line 293 "csytrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 297 "csytrs2.f"
	ctrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 301 "csytrs2.f"
	i__ = 1;
#line 302 "csytrs2.f"
	while(i__ <= *n) {
#line 303 "csytrs2.f"
	    if (ipiv[i__] > 0) {
#line 304 "csytrs2.f"
		z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
#line 304 "csytrs2.f"
		cscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
#line 305 "csytrs2.f"
	    } else {
#line 306 "csytrs2.f"
		i__1 = i__;
#line 306 "csytrs2.f"
		akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 307 "csytrs2.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
#line 307 "csytrs2.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 308 "csytrs2.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 308 "csytrs2.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 309 "csytrs2.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 309 "csytrs2.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 309 "csytrs2.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 310 "csytrs2.f"
		i__1 = *nrhs;
#line 310 "csytrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 311 "csytrs2.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
#line 311 "csytrs2.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 312 "csytrs2.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 312 "csytrs2.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 313 "csytrs2.f"
		    i__2 = i__ + j * b_dim1;
#line 313 "csytrs2.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 313 "csytrs2.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 313 "csytrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 313 "csytrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 314 "csytrs2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 314 "csytrs2.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 314 "csytrs2.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 314 "csytrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 314 "csytrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 315 "csytrs2.f"
/* L25: */
#line 315 "csytrs2.f"
		}
#line 316 "csytrs2.f"
		++i__;
#line 317 "csytrs2.f"
	    }
#line 318 "csytrs2.f"
	    ++i__;
#line 319 "csytrs2.f"
	}

/*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 323 "csytrs2.f"
	ctrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

#line 327 "csytrs2.f"
	k = *n;
#line 328 "csytrs2.f"
	while(k >= 1) {
#line 329 "csytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 332 "csytrs2.f"
		kp = ipiv[k];
#line 333 "csytrs2.f"
		if (kp != k) {
#line 333 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 333 "csytrs2.f"
		}
#line 335 "csytrs2.f"
		--k;
#line 336 "csytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 339 "csytrs2.f"
		kp = -ipiv[k];
#line 340 "csytrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 340 "csytrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 340 "csytrs2.f"
		}
#line 342 "csytrs2.f"
		k += -2;
#line 343 "csytrs2.f"
	    }
#line 344 "csytrs2.f"
	}

#line 346 "csytrs2.f"
    }

/*     Revert A */

#line 350 "csytrs2.f"
    csyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 352 "csytrs2.f"
    return 0;

/*     End of CSYTRS2 */

} /* csytrs2_ */

