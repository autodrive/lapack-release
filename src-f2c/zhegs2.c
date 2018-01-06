#line 1 "zhegs2.f"
/* zhegs2.f -- translated by f2c (version 20100827).
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

#line 1 "zhegs2.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHEGS2 reduces a Hermitian definite generalized eigenproblem to standard form, using the factor
ization results obtained from cpotrf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegs2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegs2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegs2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEGS2 reduces a complex Hermitian-definite generalized */
/* > eigenproblem to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H *A*L. */
/* > */
/* > B must have been previously factorized as U**H *U or L*L**H by ZPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); */
/* >          = 2 or 3: compute U*A*U**H or L**H *A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          Hermitian matrix A is stored, and how B has been factorized. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          n by n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the transformed matrix, stored in the */
/* >          same format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          as returned by ZPOTRF. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhegs2_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    static integer k;
    static doublecomplex ct;
    static doublereal akk, bkk;
    extern /* Subroutine */ int zher2_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztrmv_(
	    char *, char *, char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), ztrsv_(char *
	    , char *, char *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), xerbla_(char 
	    *, integer *, ftnlen), zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zlacgv_(integer *, doublecomplex *, 
	    integer *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 172 "zhegs2.f"
    /* Parameter adjustments */
#line 172 "zhegs2.f"
    a_dim1 = *lda;
#line 172 "zhegs2.f"
    a_offset = 1 + a_dim1;
#line 172 "zhegs2.f"
    a -= a_offset;
#line 172 "zhegs2.f"
    b_dim1 = *ldb;
#line 172 "zhegs2.f"
    b_offset = 1 + b_dim1;
#line 172 "zhegs2.f"
    b -= b_offset;
#line 172 "zhegs2.f"

#line 172 "zhegs2.f"
    /* Function Body */
#line 172 "zhegs2.f"
    *info = 0;
#line 173 "zhegs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 174 "zhegs2.f"
    if (*itype < 1 || *itype > 3) {
#line 175 "zhegs2.f"
	*info = -1;
#line 176 "zhegs2.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 177 "zhegs2.f"
	*info = -2;
#line 178 "zhegs2.f"
    } else if (*n < 0) {
#line 179 "zhegs2.f"
	*info = -3;
#line 180 "zhegs2.f"
    } else if (*lda < max(1,*n)) {
#line 181 "zhegs2.f"
	*info = -5;
#line 182 "zhegs2.f"
    } else if (*ldb < max(1,*n)) {
#line 183 "zhegs2.f"
	*info = -7;
#line 184 "zhegs2.f"
    }
#line 185 "zhegs2.f"
    if (*info != 0) {
#line 186 "zhegs2.f"
	i__1 = -(*info);
#line 186 "zhegs2.f"
	xerbla_("ZHEGS2", &i__1, (ftnlen)6);
#line 187 "zhegs2.f"
	return 0;
#line 188 "zhegs2.f"
    }

#line 190 "zhegs2.f"
    if (*itype == 1) {
#line 191 "zhegs2.f"
	if (upper) {

/*           Compute inv(U**H)*A*inv(U) */

#line 195 "zhegs2.f"
	    i__1 = *n;
#line 195 "zhegs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

#line 199 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 199 "zhegs2.f"
		akk = a[i__2].r;
#line 200 "zhegs2.f"
		i__2 = k + k * b_dim1;
#line 200 "zhegs2.f"
		bkk = b[i__2].r;
/* Computing 2nd power */
#line 201 "zhegs2.f"
		d__1 = bkk;
#line 201 "zhegs2.f"
		akk /= d__1 * d__1;
#line 202 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 202 "zhegs2.f"
		a[i__2].r = akk, a[i__2].i = 0.;
#line 203 "zhegs2.f"
		if (k < *n) {
#line 204 "zhegs2.f"
		    i__2 = *n - k;
#line 204 "zhegs2.f"
		    d__1 = 1. / bkk;
#line 204 "zhegs2.f"
		    zdscal_(&i__2, &d__1, &a[k + (k + 1) * a_dim1], lda);
#line 205 "zhegs2.f"
		    d__1 = akk * -.5;
#line 205 "zhegs2.f"
		    ct.r = d__1, ct.i = 0.;
#line 206 "zhegs2.f"
		    i__2 = *n - k;
#line 206 "zhegs2.f"
		    zlacgv_(&i__2, &a[k + (k + 1) * a_dim1], lda);
#line 207 "zhegs2.f"
		    i__2 = *n - k;
#line 207 "zhegs2.f"
		    zlacgv_(&i__2, &b[k + (k + 1) * b_dim1], ldb);
#line 208 "zhegs2.f"
		    i__2 = *n - k;
#line 208 "zhegs2.f"
		    zaxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
#line 210 "zhegs2.f"
		    i__2 = *n - k;
#line 210 "zhegs2.f"
		    z__1.r = -1., z__1.i = -0.;
#line 210 "zhegs2.f"
		    zher2_(uplo, &i__2, &z__1, &a[k + (k + 1) * a_dim1], lda, 
			    &b[k + (k + 1) * b_dim1], ldb, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
#line 212 "zhegs2.f"
		    i__2 = *n - k;
#line 212 "zhegs2.f"
		    zaxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
#line 214 "zhegs2.f"
		    i__2 = *n - k;
#line 214 "zhegs2.f"
		    zlacgv_(&i__2, &b[k + (k + 1) * b_dim1], ldb);
#line 215 "zhegs2.f"
		    i__2 = *n - k;
#line 215 "zhegs2.f"
		    ztrsv_(uplo, "Conjugate transpose", "Non-unit", &i__2, &b[
			    k + 1 + (k + 1) * b_dim1], ldb, &a[k + (k + 1) * 
			    a_dim1], lda, (ftnlen)1, (ftnlen)19, (ftnlen)8);
#line 218 "zhegs2.f"
		    i__2 = *n - k;
#line 218 "zhegs2.f"
		    zlacgv_(&i__2, &a[k + (k + 1) * a_dim1], lda);
#line 219 "zhegs2.f"
		}
#line 220 "zhegs2.f"
/* L10: */
#line 220 "zhegs2.f"
	    }
#line 221 "zhegs2.f"
	} else {

/*           Compute inv(L)*A*inv(L**H) */

#line 225 "zhegs2.f"
	    i__1 = *n;
#line 225 "zhegs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

#line 229 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 229 "zhegs2.f"
		akk = a[i__2].r;
#line 230 "zhegs2.f"
		i__2 = k + k * b_dim1;
#line 230 "zhegs2.f"
		bkk = b[i__2].r;
/* Computing 2nd power */
#line 231 "zhegs2.f"
		d__1 = bkk;
#line 231 "zhegs2.f"
		akk /= d__1 * d__1;
#line 232 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 232 "zhegs2.f"
		a[i__2].r = akk, a[i__2].i = 0.;
#line 233 "zhegs2.f"
		if (k < *n) {
#line 234 "zhegs2.f"
		    i__2 = *n - k;
#line 234 "zhegs2.f"
		    d__1 = 1. / bkk;
#line 234 "zhegs2.f"
		    zdscal_(&i__2, &d__1, &a[k + 1 + k * a_dim1], &c__1);
#line 235 "zhegs2.f"
		    d__1 = akk * -.5;
#line 235 "zhegs2.f"
		    ct.r = d__1, ct.i = 0.;
#line 236 "zhegs2.f"
		    i__2 = *n - k;
#line 236 "zhegs2.f"
		    zaxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
#line 237 "zhegs2.f"
		    i__2 = *n - k;
#line 237 "zhegs2.f"
		    z__1.r = -1., z__1.i = -0.;
#line 237 "zhegs2.f"
		    zher2_(uplo, &i__2, &z__1, &a[k + 1 + k * a_dim1], &c__1, 
			    &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
#line 239 "zhegs2.f"
		    i__2 = *n - k;
#line 239 "zhegs2.f"
		    zaxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
#line 240 "zhegs2.f"
		    i__2 = *n - k;
#line 240 "zhegs2.f"
		    ztrsv_(uplo, "No transpose", "Non-unit", &i__2, &b[k + 1 
			    + (k + 1) * b_dim1], ldb, &a[k + 1 + k * a_dim1], 
			    &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 242 "zhegs2.f"
		}
#line 243 "zhegs2.f"
/* L20: */
#line 243 "zhegs2.f"
	    }
#line 244 "zhegs2.f"
	}
#line 245 "zhegs2.f"
    } else {
#line 246 "zhegs2.f"
	if (upper) {

/*           Compute U*A*U**H */

#line 250 "zhegs2.f"
	    i__1 = *n;
#line 250 "zhegs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

#line 254 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 254 "zhegs2.f"
		akk = a[i__2].r;
#line 255 "zhegs2.f"
		i__2 = k + k * b_dim1;
#line 255 "zhegs2.f"
		bkk = b[i__2].r;
#line 256 "zhegs2.f"
		i__2 = k - 1;
#line 256 "zhegs2.f"
		ztrmv_(uplo, "No transpose", "Non-unit", &i__2, &b[b_offset], 
			ldb, &a[k * a_dim1 + 1], &c__1, (ftnlen)1, (ftnlen)12,
			 (ftnlen)8);
#line 258 "zhegs2.f"
		d__1 = akk * .5;
#line 258 "zhegs2.f"
		ct.r = d__1, ct.i = 0.;
#line 259 "zhegs2.f"
		i__2 = k - 1;
#line 259 "zhegs2.f"
		zaxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
#line 260 "zhegs2.f"
		i__2 = k - 1;
#line 260 "zhegs2.f"
		zher2_(uplo, &i__2, &c_b1, &a[k * a_dim1 + 1], &c__1, &b[k * 
			b_dim1 + 1], &c__1, &a[a_offset], lda, (ftnlen)1);
#line 262 "zhegs2.f"
		i__2 = k - 1;
#line 262 "zhegs2.f"
		zaxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
#line 263 "zhegs2.f"
		i__2 = k - 1;
#line 263 "zhegs2.f"
		zdscal_(&i__2, &bkk, &a[k * a_dim1 + 1], &c__1);
#line 264 "zhegs2.f"
		i__2 = k + k * a_dim1;
/* Computing 2nd power */
#line 264 "zhegs2.f"
		d__2 = bkk;
#line 264 "zhegs2.f"
		d__1 = akk * (d__2 * d__2);
#line 264 "zhegs2.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 265 "zhegs2.f"
/* L30: */
#line 265 "zhegs2.f"
	    }
#line 266 "zhegs2.f"
	} else {

/*           Compute L**H *A*L */

#line 270 "zhegs2.f"
	    i__1 = *n;
#line 270 "zhegs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

#line 274 "zhegs2.f"
		i__2 = k + k * a_dim1;
#line 274 "zhegs2.f"
		akk = a[i__2].r;
#line 275 "zhegs2.f"
		i__2 = k + k * b_dim1;
#line 275 "zhegs2.f"
		bkk = b[i__2].r;
#line 276 "zhegs2.f"
		i__2 = k - 1;
#line 276 "zhegs2.f"
		zlacgv_(&i__2, &a[k + a_dim1], lda);
#line 277 "zhegs2.f"
		i__2 = k - 1;
#line 277 "zhegs2.f"
		ztrmv_(uplo, "Conjugate transpose", "Non-unit", &i__2, &b[
			b_offset], ldb, &a[k + a_dim1], lda, (ftnlen)1, (
			ftnlen)19, (ftnlen)8);
#line 279 "zhegs2.f"
		d__1 = akk * .5;
#line 279 "zhegs2.f"
		ct.r = d__1, ct.i = 0.;
#line 280 "zhegs2.f"
		i__2 = k - 1;
#line 280 "zhegs2.f"
		zlacgv_(&i__2, &b[k + b_dim1], ldb);
#line 281 "zhegs2.f"
		i__2 = k - 1;
#line 281 "zhegs2.f"
		zaxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
#line 282 "zhegs2.f"
		i__2 = k - 1;
#line 282 "zhegs2.f"
		zher2_(uplo, &i__2, &c_b1, &a[k + a_dim1], lda, &b[k + b_dim1]
			, ldb, &a[a_offset], lda, (ftnlen)1);
#line 284 "zhegs2.f"
		i__2 = k - 1;
#line 284 "zhegs2.f"
		zaxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
#line 285 "zhegs2.f"
		i__2 = k - 1;
#line 285 "zhegs2.f"
		zlacgv_(&i__2, &b[k + b_dim1], ldb);
#line 286 "zhegs2.f"
		i__2 = k - 1;
#line 286 "zhegs2.f"
		zdscal_(&i__2, &bkk, &a[k + a_dim1], lda);
#line 287 "zhegs2.f"
		i__2 = k - 1;
#line 287 "zhegs2.f"
		zlacgv_(&i__2, &a[k + a_dim1], lda);
#line 288 "zhegs2.f"
		i__2 = k + k * a_dim1;
/* Computing 2nd power */
#line 288 "zhegs2.f"
		d__2 = bkk;
#line 288 "zhegs2.f"
		d__1 = akk * (d__2 * d__2);
#line 288 "zhegs2.f"
		a[i__2].r = d__1, a[i__2].i = 0.;
#line 289 "zhegs2.f"
/* L40: */
#line 289 "zhegs2.f"
	    }
#line 290 "zhegs2.f"
	}
#line 291 "zhegs2.f"
    }
#line 292 "zhegs2.f"
    return 0;

/*     End of ZHEGS2 */

} /* zhegs2_ */

