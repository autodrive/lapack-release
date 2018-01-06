#line 1 "ssytrs2.f"
/* ssytrs2.f -- translated by f2c (version 20100827).
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

#line 1 "ssytrs2.f"
/* Table of constant values */

static doublereal c_b10 = 1.;

/* > \brief \b SSYTRS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYTRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                           WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRS2 solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by SSYTRF and converted by SSYCONV. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by SSYTRF. */
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
/* >          as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssytrs2_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, doublereal *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ak, bk;
    static integer kp;
    static doublereal akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    static integer iinfo;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen), ssyconv_(char *, char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);


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

#line 167 "ssytrs2.f"
    /* Parameter adjustments */
#line 167 "ssytrs2.f"
    a_dim1 = *lda;
#line 167 "ssytrs2.f"
    a_offset = 1 + a_dim1;
#line 167 "ssytrs2.f"
    a -= a_offset;
#line 167 "ssytrs2.f"
    --ipiv;
#line 167 "ssytrs2.f"
    b_dim1 = *ldb;
#line 167 "ssytrs2.f"
    b_offset = 1 + b_dim1;
#line 167 "ssytrs2.f"
    b -= b_offset;
#line 167 "ssytrs2.f"
    --work;
#line 167 "ssytrs2.f"

#line 167 "ssytrs2.f"
    /* Function Body */
#line 167 "ssytrs2.f"
    *info = 0;
#line 168 "ssytrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 169 "ssytrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 170 "ssytrs2.f"
	*info = -1;
#line 171 "ssytrs2.f"
    } else if (*n < 0) {
#line 172 "ssytrs2.f"
	*info = -2;
#line 173 "ssytrs2.f"
    } else if (*nrhs < 0) {
#line 174 "ssytrs2.f"
	*info = -3;
#line 175 "ssytrs2.f"
    } else if (*lda < max(1,*n)) {
#line 176 "ssytrs2.f"
	*info = -5;
#line 177 "ssytrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 178 "ssytrs2.f"
	*info = -8;
#line 179 "ssytrs2.f"
    }
#line 180 "ssytrs2.f"
    if (*info != 0) {
#line 181 "ssytrs2.f"
	i__1 = -(*info);
#line 181 "ssytrs2.f"
	xerbla_("SSYTRS2", &i__1, (ftnlen)7);
#line 182 "ssytrs2.f"
	return 0;
#line 183 "ssytrs2.f"
    }

/*     Quick return if possible */

#line 187 "ssytrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 187 "ssytrs2.f"
	return 0;
#line 187 "ssytrs2.f"
    }

/*     Convert A */

#line 192 "ssytrs2.f"
    ssyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 194 "ssytrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*       P**T * B */
#line 199 "ssytrs2.f"
	k = *n;
#line 200 "ssytrs2.f"
	while(k >= 1) {
#line 201 "ssytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 204 "ssytrs2.f"
		kp = ipiv[k];
#line 205 "ssytrs2.f"
		if (kp != k) {
#line 205 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 205 "ssytrs2.f"
		}
#line 207 "ssytrs2.f"
		--k;
#line 208 "ssytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 211 "ssytrs2.f"
		kp = -ipiv[k];
#line 212 "ssytrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 212 "ssytrs2.f"
		    sswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 212 "ssytrs2.f"
		}
#line 214 "ssytrs2.f"
		k += -2;
#line 215 "ssytrs2.f"
	    }
#line 216 "ssytrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 220 "ssytrs2.f"
	strsm_("L", "U", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 224 "ssytrs2.f"
	i__ = *n;
#line 225 "ssytrs2.f"
	while(i__ >= 1) {
#line 226 "ssytrs2.f"
	    if (ipiv[i__] > 0) {
#line 227 "ssytrs2.f"
		d__1 = 1. / a[i__ + i__ * a_dim1];
#line 227 "ssytrs2.f"
		sscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
#line 228 "ssytrs2.f"
	    } else if (i__ > 1) {
#line 229 "ssytrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 230 "ssytrs2.f"
		    akm1k = work[i__];
#line 231 "ssytrs2.f"
		    akm1 = a[i__ - 1 + (i__ - 1) * a_dim1] / akm1k;
#line 232 "ssytrs2.f"
		    ak = a[i__ + i__ * a_dim1] / akm1k;
#line 233 "ssytrs2.f"
		    denom = akm1 * ak - 1.;
#line 234 "ssytrs2.f"
		    i__1 = *nrhs;
#line 234 "ssytrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 235 "ssytrs2.f"
			bkm1 = b[i__ - 1 + j * b_dim1] / akm1k;
#line 236 "ssytrs2.f"
			bk = b[i__ + j * b_dim1] / akm1k;
#line 237 "ssytrs2.f"
			b[i__ - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 238 "ssytrs2.f"
			b[i__ + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 239 "ssytrs2.f"
/* L15: */
#line 239 "ssytrs2.f"
		    }
#line 240 "ssytrs2.f"
		    --i__;
#line 241 "ssytrs2.f"
		}
#line 242 "ssytrs2.f"
	    }
#line 243 "ssytrs2.f"
	    --i__;
#line 244 "ssytrs2.f"
	}

/*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 248 "ssytrs2.f"
	strsm_("L", "U", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

#line 252 "ssytrs2.f"
	k = 1;
#line 253 "ssytrs2.f"
	while(k <= *n) {
#line 254 "ssytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 257 "ssytrs2.f"
		kp = ipiv[k];
#line 258 "ssytrs2.f"
		if (kp != k) {
#line 258 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 258 "ssytrs2.f"
		}
#line 260 "ssytrs2.f"
		++k;
#line 261 "ssytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 264 "ssytrs2.f"
		kp = -ipiv[k];
#line 265 "ssytrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 265 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 265 "ssytrs2.f"
		}
#line 267 "ssytrs2.f"
		k += 2;
#line 268 "ssytrs2.f"
	    }
#line 269 "ssytrs2.f"
	}

#line 271 "ssytrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*       P**T * B */
#line 276 "ssytrs2.f"
	k = 1;
#line 277 "ssytrs2.f"
	while(k <= *n) {
#line 278 "ssytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 281 "ssytrs2.f"
		kp = ipiv[k];
#line 282 "ssytrs2.f"
		if (kp != k) {
#line 282 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "ssytrs2.f"
		}
#line 284 "ssytrs2.f"
		++k;
#line 285 "ssytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 288 "ssytrs2.f"
		kp = -ipiv[k + 1];
#line 289 "ssytrs2.f"
		if (kp == -ipiv[k]) {
#line 289 "ssytrs2.f"
		    sswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 289 "ssytrs2.f"
		}
#line 291 "ssytrs2.f"
		k += 2;
#line 292 "ssytrs2.f"
	    }
#line 293 "ssytrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 297 "ssytrs2.f"
	strsm_("L", "L", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 301 "ssytrs2.f"
	i__ = 1;
#line 302 "ssytrs2.f"
	while(i__ <= *n) {
#line 303 "ssytrs2.f"
	    if (ipiv[i__] > 0) {
#line 304 "ssytrs2.f"
		d__1 = 1. / a[i__ + i__ * a_dim1];
#line 304 "ssytrs2.f"
		sscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
#line 305 "ssytrs2.f"
	    } else {
#line 306 "ssytrs2.f"
		akm1k = work[i__];
#line 307 "ssytrs2.f"
		akm1 = a[i__ + i__ * a_dim1] / akm1k;
#line 308 "ssytrs2.f"
		ak = a[i__ + 1 + (i__ + 1) * a_dim1] / akm1k;
#line 309 "ssytrs2.f"
		denom = akm1 * ak - 1.;
#line 310 "ssytrs2.f"
		i__1 = *nrhs;
#line 310 "ssytrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 311 "ssytrs2.f"
		    bkm1 = b[i__ + j * b_dim1] / akm1k;
#line 312 "ssytrs2.f"
		    bk = b[i__ + 1 + j * b_dim1] / akm1k;
#line 313 "ssytrs2.f"
		    b[i__ + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 314 "ssytrs2.f"
		    b[i__ + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 315 "ssytrs2.f"
/* L25: */
#line 315 "ssytrs2.f"
		}
#line 316 "ssytrs2.f"
		++i__;
#line 317 "ssytrs2.f"
	    }
#line 318 "ssytrs2.f"
	    ++i__;
#line 319 "ssytrs2.f"
	}

/*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 323 "ssytrs2.f"
	strsm_("L", "L", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

#line 327 "ssytrs2.f"
	k = *n;
#line 328 "ssytrs2.f"
	while(k >= 1) {
#line 329 "ssytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 332 "ssytrs2.f"
		kp = ipiv[k];
#line 333 "ssytrs2.f"
		if (kp != k) {
#line 333 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 333 "ssytrs2.f"
		}
#line 335 "ssytrs2.f"
		--k;
#line 336 "ssytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 339 "ssytrs2.f"
		kp = -ipiv[k];
#line 340 "ssytrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 340 "ssytrs2.f"
		    sswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 340 "ssytrs2.f"
		}
#line 342 "ssytrs2.f"
		k += -2;
#line 343 "ssytrs2.f"
	    }
#line 344 "ssytrs2.f"
	}

#line 346 "ssytrs2.f"
    }

/*     Revert A */

#line 350 "ssytrs2.f"
    ssyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 352 "ssytrs2.f"
    return 0;

/*     End of SSYTRS2 */

} /* ssytrs2_ */

