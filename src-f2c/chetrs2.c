#line 1 "chetrs2.f"
/* chetrs2.f -- translated by f2c (version 20100827).
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

#line 1 "chetrs2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};

/* > \brief \b CHETRS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHETRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrs2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrs2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrs2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
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
/* > CHETRS2 solves a system of linear equations A*X = B with a complex */
/* > Hermitian matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by CHETRF and converted by CSYCONV. */
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

/* > \date December 2016 */

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
/* Subroutine */ int chetrs2_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ctrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublecomplex *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen), 
	    csyconv_(char *, char *, integer *, doublecomplex *, integer *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen);


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

#line 168 "chetrs2.f"
    /* Parameter adjustments */
#line 168 "chetrs2.f"
    a_dim1 = *lda;
#line 168 "chetrs2.f"
    a_offset = 1 + a_dim1;
#line 168 "chetrs2.f"
    a -= a_offset;
#line 168 "chetrs2.f"
    --ipiv;
#line 168 "chetrs2.f"
    b_dim1 = *ldb;
#line 168 "chetrs2.f"
    b_offset = 1 + b_dim1;
#line 168 "chetrs2.f"
    b -= b_offset;
#line 168 "chetrs2.f"
    --work;
#line 168 "chetrs2.f"

#line 168 "chetrs2.f"
    /* Function Body */
#line 168 "chetrs2.f"
    *info = 0;
#line 169 "chetrs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "chetrs2.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 171 "chetrs2.f"
	*info = -1;
#line 172 "chetrs2.f"
    } else if (*n < 0) {
#line 173 "chetrs2.f"
	*info = -2;
#line 174 "chetrs2.f"
    } else if (*nrhs < 0) {
#line 175 "chetrs2.f"
	*info = -3;
#line 176 "chetrs2.f"
    } else if (*lda < max(1,*n)) {
#line 177 "chetrs2.f"
	*info = -5;
#line 178 "chetrs2.f"
    } else if (*ldb < max(1,*n)) {
#line 179 "chetrs2.f"
	*info = -8;
#line 180 "chetrs2.f"
    }
#line 181 "chetrs2.f"
    if (*info != 0) {
#line 182 "chetrs2.f"
	i__1 = -(*info);
#line 182 "chetrs2.f"
	xerbla_("CHETRS2", &i__1, (ftnlen)7);
#line 183 "chetrs2.f"
	return 0;
#line 184 "chetrs2.f"
    }

/*     Quick return if possible */

#line 188 "chetrs2.f"
    if (*n == 0 || *nrhs == 0) {
#line 188 "chetrs2.f"
	return 0;
#line 188 "chetrs2.f"
    }

/*     Convert A */

#line 193 "chetrs2.f"
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 195 "chetrs2.f"
    if (upper) {

/*        Solve A*X = B, where A = U*D*U**H. */

/*       P**T * B */
#line 200 "chetrs2.f"
	k = *n;
#line 201 "chetrs2.f"
	while(k >= 1) {
#line 202 "chetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 205 "chetrs2.f"
		kp = ipiv[k];
#line 206 "chetrs2.f"
		if (kp != k) {
#line 206 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 206 "chetrs2.f"
		}
#line 208 "chetrs2.f"
		--k;
#line 209 "chetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 212 "chetrs2.f"
		kp = -ipiv[k];
#line 213 "chetrs2.f"
		if (kp == -ipiv[k - 1]) {
#line 213 "chetrs2.f"
		    cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 213 "chetrs2.f"
		}
#line 215 "chetrs2.f"
		k += -2;
#line 216 "chetrs2.f"
	    }
#line 217 "chetrs2.f"
	}

/*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ] */

#line 221 "chetrs2.f"
	ctrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (U \P**T * B) ] */

#line 225 "chetrs2.f"
	i__ = *n;
#line 226 "chetrs2.f"
	while(i__ >= 1) {
#line 227 "chetrs2.f"
	    if (ipiv[i__] > 0) {
#line 228 "chetrs2.f"
		i__1 = i__ + i__ * a_dim1;
#line 228 "chetrs2.f"
		s = 1. / a[i__1].r;
#line 229 "chetrs2.f"
		csscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 230 "chetrs2.f"
	    } else if (i__ > 1) {
#line 231 "chetrs2.f"
		if (ipiv[i__ - 1] == ipiv[i__]) {
#line 232 "chetrs2.f"
		    i__1 = i__;
#line 232 "chetrs2.f"
		    akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 233 "chetrs2.f"
		    z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
#line 233 "chetrs2.f"
		    akm1.r = z__1.r, akm1.i = z__1.i;
#line 234 "chetrs2.f"
		    d_cnjg(&z__2, &akm1k);
#line 234 "chetrs2.f"
		    z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 234 "chetrs2.f"
		    ak.r = z__1.r, ak.i = z__1.i;
#line 235 "chetrs2.f"
		    z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			    ak.i + akm1.i * ak.r;
#line 235 "chetrs2.f"
		    z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 235 "chetrs2.f"
		    denom.r = z__1.r, denom.i = z__1.i;
#line 236 "chetrs2.f"
		    i__1 = *nrhs;
#line 236 "chetrs2.f"
		    for (j = 1; j <= i__1; ++j) {
#line 237 "chetrs2.f"
			z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
#line 237 "chetrs2.f"
			bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 238 "chetrs2.f"
			d_cnjg(&z__2, &akm1k);
#line 238 "chetrs2.f"
			z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 238 "chetrs2.f"
			bk.r = z__1.r, bk.i = z__1.i;
#line 239 "chetrs2.f"
			i__2 = i__ - 1 + j * b_dim1;
#line 239 "chetrs2.f"
			z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r 
				* bkm1.i + ak.i * bkm1.r;
#line 239 "chetrs2.f"
			z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 239 "chetrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 239 "chetrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 240 "chetrs2.f"
			i__2 = i__ + j * b_dim1;
#line 240 "chetrs2.f"
			z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = 
				akm1.r * bk.i + akm1.i * bk.r;
#line 240 "chetrs2.f"
			z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 240 "chetrs2.f"
			z_div(&z__1, &z__2, &denom);
#line 240 "chetrs2.f"
			b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 241 "chetrs2.f"
/* L15: */
#line 241 "chetrs2.f"
		    }
#line 242 "chetrs2.f"
		    --i__;
#line 243 "chetrs2.f"
		}
#line 244 "chetrs2.f"
	    }
#line 245 "chetrs2.f"
	    --i__;
#line 246 "chetrs2.f"
	}

/*      Compute (U**H \ B) -> B   [ U**H \ (D \ (U \P**T * B) ) ] */

#line 250 "chetrs2.f"
	ctrsm_("L", "U", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (U**H \ (D \ (U \P**T * B) )) ] */

#line 254 "chetrs2.f"
	k = 1;
#line 255 "chetrs2.f"
	while(k <= *n) {
#line 256 "chetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 259 "chetrs2.f"
		kp = ipiv[k];
#line 260 "chetrs2.f"
		if (kp != k) {
#line 260 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 260 "chetrs2.f"
		}
#line 262 "chetrs2.f"
		++k;
#line 263 "chetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 266 "chetrs2.f"
		kp = -ipiv[k];
#line 267 "chetrs2.f"
		if (k < *n && kp == -ipiv[k + 1]) {
#line 267 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 267 "chetrs2.f"
		}
#line 269 "chetrs2.f"
		k += 2;
#line 270 "chetrs2.f"
	    }
#line 271 "chetrs2.f"
	}

#line 273 "chetrs2.f"
    } else {

/*        Solve A*X = B, where A = L*D*L**H. */

/*       P**T * B */
#line 278 "chetrs2.f"
	k = 1;
#line 279 "chetrs2.f"
	while(k <= *n) {
#line 280 "chetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 283 "chetrs2.f"
		kp = ipiv[k];
#line 284 "chetrs2.f"
		if (kp != k) {
#line 284 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 284 "chetrs2.f"
		}
#line 286 "chetrs2.f"
		++k;
#line 287 "chetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K and -IPIV(K+1). */
#line 290 "chetrs2.f"
		kp = -ipiv[k + 1];
#line 291 "chetrs2.f"
		if (kp == -ipiv[k]) {
#line 291 "chetrs2.f"
		    cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], 
			    ldb);
#line 291 "chetrs2.f"
		}
#line 293 "chetrs2.f"
		k += 2;
#line 294 "chetrs2.f"
	    }
#line 295 "chetrs2.f"
	}

/*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ] */

#line 299 "chetrs2.f"
	ctrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*  Compute D \ B -> B   [ D \ (L \P**T * B) ] */

#line 303 "chetrs2.f"
	i__ = 1;
#line 304 "chetrs2.f"
	while(i__ <= *n) {
#line 305 "chetrs2.f"
	    if (ipiv[i__] > 0) {
#line 306 "chetrs2.f"
		i__1 = i__ + i__ * a_dim1;
#line 306 "chetrs2.f"
		s = 1. / a[i__1].r;
#line 307 "chetrs2.f"
		csscal_(nrhs, &s, &b[i__ + b_dim1], ldb);
#line 308 "chetrs2.f"
	    } else {
#line 309 "chetrs2.f"
		i__1 = i__;
#line 309 "chetrs2.f"
		akm1k.r = work[i__1].r, akm1k.i = work[i__1].i;
#line 310 "chetrs2.f"
		d_cnjg(&z__2, &akm1k);
#line 310 "chetrs2.f"
		z_div(&z__1, &a[i__ + i__ * a_dim1], &z__2);
#line 310 "chetrs2.f"
		akm1.r = z__1.r, akm1.i = z__1.i;
#line 311 "chetrs2.f"
		z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
#line 311 "chetrs2.f"
		ak.r = z__1.r, ak.i = z__1.i;
#line 312 "chetrs2.f"
		z__2.r = akm1.r * ak.r - akm1.i * ak.i, z__2.i = akm1.r * 
			ak.i + akm1.i * ak.r;
#line 312 "chetrs2.f"
		z__1.r = z__2.r - 1., z__1.i = z__2.i - 0.;
#line 312 "chetrs2.f"
		denom.r = z__1.r, denom.i = z__1.i;
#line 313 "chetrs2.f"
		i__1 = *nrhs;
#line 313 "chetrs2.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "chetrs2.f"
		    d_cnjg(&z__2, &akm1k);
#line 314 "chetrs2.f"
		    z_div(&z__1, &b[i__ + j * b_dim1], &z__2);
#line 314 "chetrs2.f"
		    bkm1.r = z__1.r, bkm1.i = z__1.i;
#line 315 "chetrs2.f"
		    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
#line 315 "chetrs2.f"
		    bk.r = z__1.r, bk.i = z__1.i;
#line 316 "chetrs2.f"
		    i__2 = i__ + j * b_dim1;
#line 316 "chetrs2.f"
		    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * 
			    bkm1.i + ak.i * bkm1.r;
#line 316 "chetrs2.f"
		    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
#line 316 "chetrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 316 "chetrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 317 "chetrs2.f"
		    i__2 = i__ + 1 + j * b_dim1;
#line 317 "chetrs2.f"
		    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * 
			    bk.i + akm1.i * bk.r;
#line 317 "chetrs2.f"
		    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
#line 317 "chetrs2.f"
		    z_div(&z__1, &z__2, &denom);
#line 317 "chetrs2.f"
		    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
#line 318 "chetrs2.f"
/* L25: */
#line 318 "chetrs2.f"
		}
#line 319 "chetrs2.f"
		++i__;
#line 320 "chetrs2.f"
	    }
#line 321 "chetrs2.f"
	    ++i__;
#line 322 "chetrs2.f"
	}

/*  Compute (L**H \ B) -> B   [ L**H \ (D \ (L \P**T * B) ) ] */

#line 326 "chetrs2.f"
	ctrsm_("L", "L", "C", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[
		b_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*       P * B  [ P * (L**H \ (D \ (L \P**T * B) )) ] */

#line 330 "chetrs2.f"
	k = *n;
#line 331 "chetrs2.f"
	while(k >= 1) {
#line 332 "chetrs2.f"
	    if (ipiv[k] > 0) {
/*           1 x 1 diagonal block */
/*           Interchange rows K and IPIV(K). */
#line 335 "chetrs2.f"
		kp = ipiv[k];
#line 336 "chetrs2.f"
		if (kp != k) {
#line 336 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 336 "chetrs2.f"
		}
#line 338 "chetrs2.f"
		--k;
#line 339 "chetrs2.f"
	    } else {
/*           2 x 2 diagonal block */
/*           Interchange rows K-1 and -IPIV(K). */
#line 342 "chetrs2.f"
		kp = -ipiv[k];
#line 343 "chetrs2.f"
		if (k > 1 && kp == -ipiv[k - 1]) {
#line 343 "chetrs2.f"
		    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
#line 343 "chetrs2.f"
		}
#line 345 "chetrs2.f"
		k += -2;
#line 346 "chetrs2.f"
	    }
#line 347 "chetrs2.f"
	}

#line 349 "chetrs2.f"
    }

/*     Revert A */

#line 353 "chetrs2.f"
    csyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (
	    ftnlen)1, (ftnlen)1);

#line 355 "chetrs2.f"
    return 0;

/*     End of CHETRS2 */

} /* chetrs2_ */

