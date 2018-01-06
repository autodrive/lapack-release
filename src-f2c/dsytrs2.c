#line 1 "dsytrs2.f"
/* dsytrs2.f -- translated by f2c (version 20100827).
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

#line 1 "dsytrs2.f"
/* Table of constant values */

static doublereal c_b10 = 1.;

/* > \brief \b DSYTRS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYTRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/*                           WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRS2 solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSYTRF and converted by DSYCONV. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSYTRF. */
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
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsytrs2_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    static integer iinfo;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dsyconv_(
	    char *, char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);


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

#line 167 "dsytrs2.f"
    /* Parameter adjustments */
#line 167 "dsytrs2.f"
    a_dim1 = *lda;
#line 167 "dsytrs2.f"
    a_offset = 1 + a_dim1;
#line 167 "dsytrs2.f"
    a -= a_offset;
#line 167 "dsytrs2.f"
    --ipiv;
#line 167 "dsytrs2.f"
    b_dim1 = *ldb;
#line 167 "dsytrs2.f"
    b_offset = 1 + b_dim1;
#line 167 "dsytrs2.f"
    b -= b_offset;
#line 167 "dsytrs2.f"
    --work;
#line 167 "dsytrs2.f"

#line 167 "dsytrs2.f"
    /* Function Body */
#line 167 "dsytrs2.f"
    *info = 0;
#line 168 "dsytrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 169 "dsytrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 170 "dsytrs2.f"
	*info = -1;
#line 171 "dsytrs2.f"
    } else if (*n < 0) {
#line 172 "dsytrs2.f"
	*info = -2;
#line 173 "dsytrs2.f"
    } else if (*nrhs < 0) {
#line 174 "dsytrs2.f"
	*info = -3;
#line 175 "dsytrs2.f"
    } else if (*lda < max(1,*n)) {
#line 176 "dsytrs2.f"
	*info = -5;
#line 177 "dsytrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 178 "dsytrs2.f"
	*info = -8;
#line 179 "dsytrs2.f"
    }
#line 180 "dsytrs2.f"
    if (*info != 0) {
#line 181 "dsytrs2.f"
	i__1 = -(*info);
#line 181 "dsytrs2.f"
	xerbla_("DSYTRS2", &i__1, (ftnlen)7);
#line 182 "dsytrs2.f"
	return 0;
#line 183 "dsytrs2.f"
    }

/*     Quick return if possible */

#line 187 "dsytrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 187 "dsytrs2.f"
	return 0;
#line 187 "dsytrs2.f"
    }

/*     Convert A */

#line 192 "dsytrs2.f"
    dsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 194 "dsytrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*       P**T * B */
#line 199 "dsytrs2.f"
	k = *n;
#line 200 "dsytrs2.f"
	while(k >= 1) {
#line 201 "dsytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 204 "dsytrs2.f"
		kp = ipiv[k];
#line 205 "dsytrs2.f"
		if (kp != k) {
#line 205 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 205 "dsytrs2.f"
		}
#line 207 "dsytrs2.f"
		--k;
#line 208 "dsytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 211 "dsytrs2.f"
		kp = -ipiv[k];
#line 212 "dsytrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 212 "dsytrs2.f"
		    dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 212 "dsytrs2.f"
		}
#line 214 "dsytrs2.f"
		k += -2;
#line 215 "dsytrs2.f"
	    }
#line 216 "dsytrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 220 "dsytrs2.f"
	dtrsm_("L", "U", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 224 "dsytrs2.f"
	i__ = *n;
#line 225 "dsytrs2.f"
	while(i__ >= 1) {
#line 226 "dsytrs2.f"
	    if (ipiv[i__] > 0) {
#line 227 "dsytrs2.f"
		d__1 = 1. / a[i__ + i__ * a_dim1];
#line 227 "dsytrs2.f"
		dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
#line 228 "dsytrs2.f"
	    } else if (i__ > 1) {
#line 229 "dsytrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 230 "dsytrs2.f"
		    akm1k = work[i__];
#line 231 "dsytrs2.f"
		    akm1 = a[i__ - 1 + (i__ - 1) * a_dim1] / akm1k;
#line 232 "dsytrs2.f"
		    ak = a[i__ + i__ * a_dim1] / akm1k;
#line 233 "dsytrs2.f"
		    denom = akm1 * ak - 1.;
#line 234 "dsytrs2.f"
		    i__1 = *nrhs;
#line 234 "dsytrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 235 "dsytrs2.f"
			bkm1 = b[i__ - 1 + j * b_dim1] / akm1k;
#line 236 "dsytrs2.f"
			bk = b[i__ + j * b_dim1] / akm1k;
#line 237 "dsytrs2.f"
			b[i__ - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 238 "dsytrs2.f"
			b[i__ + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 239 "dsytrs2.f"
/* L15: */
#line 239 "dsytrs2.f"
		    }
#line 240 "dsytrs2.f"
		    --i__;
#line 241 "dsytrs2.f"
		}
#line 242 "dsytrs2.f"
	    }
#line 243 "dsytrs2.f"
	    --i__;
#line 244 "dsytrs2.f"
	}

/*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ] */

#line 248 "dsytrs2.f"
	dtrsm_("L", "U", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ] */

#line 252 "dsytrs2.f"
	k = 1;
#line 253 "dsytrs2.f"
	while(k <= *n) {
#line 254 "dsytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 257 "dsytrs2.f"
		kp = ipiv[k];
#line 258 "dsytrs2.f"
		if (kp != k) {
#line 258 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 258 "dsytrs2.f"
		}
#line 260 "dsytrs2.f"
		++k;
#line 261 "dsytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 264 "dsytrs2.f"
		kp = -ipiv[k];
#line 265 "dsytrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 265 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 265 "dsytrs2.f"
		}
#line 267 "dsytrs2.f"
		k += 2;
#line 268 "dsytrs2.f"
	    }
#line 269 "dsytrs2.f"
	}

#line 271 "dsytrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*       P**T * B */
#line 276 "dsytrs2.f"
	k = 1;
#line 277 "dsytrs2.f"
	while(k <= *n) {
#line 278 "dsytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 281 "dsytrs2.f"
		kp = ipiv[k];
#line 282 "dsytrs2.f"
		if (kp != k) {
#line 282 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 282 "dsytrs2.f"
		}
#line 284 "dsytrs2.f"
		++k;
#line 285 "dsytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 288 "dsytrs2.f"
		kp = -ipiv[k + 1];
#line 289 "dsytrs2.f"
		if (kp == -ipiv[k]) {
#line 289 "dsytrs2.f"
		    dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 289 "dsytrs2.f"
		}
#line 291 "dsytrs2.f"
		k += 2;
#line 292 "dsytrs2.f"
	    }
#line 293 "dsytrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 297 "dsytrs2.f"
	dtrsm_("L", "L", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 301 "dsytrs2.f"
	i__ = 1;
#line 302 "dsytrs2.f"
	while(i__ <= *n) {
#line 303 "dsytrs2.f"
	    if (ipiv[i__] > 0) {
#line 304 "dsytrs2.f"
		d__1 = 1. / a[i__ + i__ * a_dim1];
#line 304 "dsytrs2.f"
		dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
#line 305 "dsytrs2.f"
	    } else {
#line 306 "dsytrs2.f"
		akm1k = work[i__];
#line 307 "dsytrs2.f"
		akm1 = a[i__ + i__ * a_dim1] / akm1k;
#line 308 "dsytrs2.f"
		ak = a[i__ + 1 + (i__ + 1) * a_dim1] / akm1k;
#line 309 "dsytrs2.f"
		denom = akm1 * ak - 1.;
#line 310 "dsytrs2.f"
		i__1 = *nrhs;
#line 310 "dsytrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 311 "dsytrs2.f"
		    bkm1 = b[i__ + j * b_dim1] / akm1k;
#line 312 "dsytrs2.f"
		    bk = b[i__ + 1 + j * b_dim1] / akm1k;
#line 313 "dsytrs2.f"
		    b[i__ + j * b_dim1] = (ak * bkm1 - bk) / denom;
#line 314 "dsytrs2.f"
		    b[i__ + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
#line 315 "dsytrs2.f"
/* L25: */
#line 315 "dsytrs2.f"
		}
#line 316 "dsytrs2.f"
		++i__;
#line 317 "dsytrs2.f"
	    }
#line 318 "dsytrs2.f"
	    ++i__;
#line 319 "dsytrs2.f"
	}

/*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ] */

#line 323 "dsytrs2.f"
	dtrsm_("L", "L", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ] */

#line 327 "dsytrs2.f"
	k = *n;
#line 328 "dsytrs2.f"
	while(k >= 1) {
#line 329 "dsytrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 332 "dsytrs2.f"
		kp = ipiv[k];
#line 333 "dsytrs2.f"
		if (kp != k) {
#line 333 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 333 "dsytrs2.f"
		}
#line 335 "dsytrs2.f"
		--k;
#line 336 "dsytrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 339 "dsytrs2.f"
		kp = -ipiv[k];
#line 340 "dsytrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 340 "dsytrs2.f"
		    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 340 "dsytrs2.f"
		}
#line 342 "dsytrs2.f"
		k += -2;
#line 343 "dsytrs2.f"
	    }
#line 344 "dsytrs2.f"
	}

#line 346 "dsytrs2.f"
    }

/*     Revert A */

#line 350 "dsytrs2.f"
    dsyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 352 "dsytrs2.f"
    return 0;

/*     End of DSYTRS2 */

} /* dsytrs2_ */

