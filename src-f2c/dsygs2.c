#line 1 "dsygs2.f"
/* dsygs2.f -- translated by f2c (version 20100827).
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

#line 1 "dsygs2.f"
/* Table of constant values */

static doublereal c_b6 = -1.;
static integer c__1 = 1;
static doublereal c_b27 = 1.;

/* > \brief \b DSYGS2 reduces a symmetric definite generalized eigenproblem to standard form, using the factor
ization results obtained from spotrf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYGS2 reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T *A*L. */
/* > */
/* > B must have been previously factorized as U**T *U or L*L**T by DPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); */
/* >          = 2 or 3: compute U*A*U**T or L**T *A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored, and how B has been factorized. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
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
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          as returned by DPOTRF. */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsygs2_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal ct, akk, bkk;
    extern /* Subroutine */ int dsyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dtrsv_(char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

#line 168 "dsygs2.f"
    /* Parameter adjustments */
#line 168 "dsygs2.f"
    a_dim1 = *lda;
#line 168 "dsygs2.f"
    a_offset = 1 + a_dim1;
#line 168 "dsygs2.f"
    a -= a_offset;
#line 168 "dsygs2.f"
    b_dim1 = *ldb;
#line 168 "dsygs2.f"
    b_offset = 1 + b_dim1;
#line 168 "dsygs2.f"
    b -= b_offset;
#line 168 "dsygs2.f"

#line 168 "dsygs2.f"
    /* Function Body */
#line 168 "dsygs2.f"
    *info = 0;
#line 169 "dsygs2.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "dsygs2.f"
    if (*itype < 1 || *itype > 3) {
#line 171 "dsygs2.f"
	*info = -1;
#line 172 "dsygs2.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "dsygs2.f"
	*info = -2;
#line 174 "dsygs2.f"
    } else if (*n < 0) {
#line 175 "dsygs2.f"
	*info = -3;
#line 176 "dsygs2.f"
    } else if (*lda < max(1,*n)) {
#line 177 "dsygs2.f"
	*info = -5;
#line 178 "dsygs2.f"
    } else if (*ldb < max(1,*n)) {
#line 179 "dsygs2.f"
	*info = -7;
#line 180 "dsygs2.f"
    }
#line 181 "dsygs2.f"
    if (*info != 0) {
#line 182 "dsygs2.f"
	i__1 = -(*info);
#line 182 "dsygs2.f"
	xerbla_("DSYGS2", &i__1, (ftnlen)6);
#line 183 "dsygs2.f"
	return 0;
#line 184 "dsygs2.f"
    }

#line 186 "dsygs2.f"
    if (*itype == 1) {
#line 187 "dsygs2.f"
	if (upper) {

/*           Compute inv(U**T)*A*inv(U) */

#line 191 "dsygs2.f"
	    i__1 = *n;
#line 191 "dsygs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

#line 195 "dsygs2.f"
		akk = a[k + k * a_dim1];
#line 196 "dsygs2.f"
		bkk = b[k + k * b_dim1];
/* Computing 2nd power */
#line 197 "dsygs2.f"
		d__1 = bkk;
#line 197 "dsygs2.f"
		akk /= d__1 * d__1;
#line 198 "dsygs2.f"
		a[k + k * a_dim1] = akk;
#line 199 "dsygs2.f"
		if (k < *n) {
#line 200 "dsygs2.f"
		    i__2 = *n - k;
#line 200 "dsygs2.f"
		    d__1 = 1. / bkk;
#line 200 "dsygs2.f"
		    dscal_(&i__2, &d__1, &a[k + (k + 1) * a_dim1], lda);
#line 201 "dsygs2.f"
		    ct = akk * -.5;
#line 202 "dsygs2.f"
		    i__2 = *n - k;
#line 202 "dsygs2.f"
		    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
#line 204 "dsygs2.f"
		    i__2 = *n - k;
#line 204 "dsygs2.f"
		    dsyr2_(uplo, &i__2, &c_b6, &a[k + (k + 1) * a_dim1], lda, 
			    &b[k + (k + 1) * b_dim1], ldb, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
#line 206 "dsygs2.f"
		    i__2 = *n - k;
#line 206 "dsygs2.f"
		    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
#line 208 "dsygs2.f"
		    i__2 = *n - k;
#line 208 "dsygs2.f"
		    dtrsv_(uplo, "Transpose", "Non-unit", &i__2, &b[k + 1 + (
			    k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1], 
			    lda, (ftnlen)1, (ftnlen)9, (ftnlen)8);
#line 210 "dsygs2.f"
		}
#line 211 "dsygs2.f"
/* L10: */
#line 211 "dsygs2.f"
	    }
#line 212 "dsygs2.f"
	} else {

/*           Compute inv(L)*A*inv(L**T) */

#line 216 "dsygs2.f"
	    i__1 = *n;
#line 216 "dsygs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

#line 220 "dsygs2.f"
		akk = a[k + k * a_dim1];
#line 221 "dsygs2.f"
		bkk = b[k + k * b_dim1];
/* Computing 2nd power */
#line 222 "dsygs2.f"
		d__1 = bkk;
#line 222 "dsygs2.f"
		akk /= d__1 * d__1;
#line 223 "dsygs2.f"
		a[k + k * a_dim1] = akk;
#line 224 "dsygs2.f"
		if (k < *n) {
#line 225 "dsygs2.f"
		    i__2 = *n - k;
#line 225 "dsygs2.f"
		    d__1 = 1. / bkk;
#line 225 "dsygs2.f"
		    dscal_(&i__2, &d__1, &a[k + 1 + k * a_dim1], &c__1);
#line 226 "dsygs2.f"
		    ct = akk * -.5;
#line 227 "dsygs2.f"
		    i__2 = *n - k;
#line 227 "dsygs2.f"
		    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
#line 228 "dsygs2.f"
		    i__2 = *n - k;
#line 228 "dsygs2.f"
		    dsyr2_(uplo, &i__2, &c_b6, &a[k + 1 + k * a_dim1], &c__1, 
			    &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
#line 230 "dsygs2.f"
		    i__2 = *n - k;
#line 230 "dsygs2.f"
		    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
#line 231 "dsygs2.f"
		    i__2 = *n - k;
#line 231 "dsygs2.f"
		    dtrsv_(uplo, "No transpose", "Non-unit", &i__2, &b[k + 1 
			    + (k + 1) * b_dim1], ldb, &a[k + 1 + k * a_dim1], 
			    &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 233 "dsygs2.f"
		}
#line 234 "dsygs2.f"
/* L20: */
#line 234 "dsygs2.f"
	    }
#line 235 "dsygs2.f"
	}
#line 236 "dsygs2.f"
    } else {
#line 237 "dsygs2.f"
	if (upper) {

/*           Compute U*A*U**T */

#line 241 "dsygs2.f"
	    i__1 = *n;
#line 241 "dsygs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

#line 245 "dsygs2.f"
		akk = a[k + k * a_dim1];
#line 246 "dsygs2.f"
		bkk = b[k + k * b_dim1];
#line 247 "dsygs2.f"
		i__2 = k - 1;
#line 247 "dsygs2.f"
		dtrmv_(uplo, "No transpose", "Non-unit", &i__2, &b[b_offset], 
			ldb, &a[k * a_dim1 + 1], &c__1, (ftnlen)1, (ftnlen)12,
			 (ftnlen)8);
#line 249 "dsygs2.f"
		ct = akk * .5;
#line 250 "dsygs2.f"
		i__2 = k - 1;
#line 250 "dsygs2.f"
		daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
#line 251 "dsygs2.f"
		i__2 = k - 1;
#line 251 "dsygs2.f"
		dsyr2_(uplo, &i__2, &c_b27, &a[k * a_dim1 + 1], &c__1, &b[k * 
			b_dim1 + 1], &c__1, &a[a_offset], lda, (ftnlen)1);
#line 253 "dsygs2.f"
		i__2 = k - 1;
#line 253 "dsygs2.f"
		daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
#line 254 "dsygs2.f"
		i__2 = k - 1;
#line 254 "dsygs2.f"
		dscal_(&i__2, &bkk, &a[k * a_dim1 + 1], &c__1);
/* Computing 2nd power */
#line 255 "dsygs2.f"
		d__1 = bkk;
#line 255 "dsygs2.f"
		a[k + k * a_dim1] = akk * (d__1 * d__1);
#line 256 "dsygs2.f"
/* L30: */
#line 256 "dsygs2.f"
	    }
#line 257 "dsygs2.f"
	} else {

/*           Compute L**T *A*L */

#line 261 "dsygs2.f"
	    i__1 = *n;
#line 261 "dsygs2.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

#line 265 "dsygs2.f"
		akk = a[k + k * a_dim1];
#line 266 "dsygs2.f"
		bkk = b[k + k * b_dim1];
#line 267 "dsygs2.f"
		i__2 = k - 1;
#line 267 "dsygs2.f"
		dtrmv_(uplo, "Transpose", "Non-unit", &i__2, &b[b_offset], 
			ldb, &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)9, (
			ftnlen)8);
#line 269 "dsygs2.f"
		ct = akk * .5;
#line 270 "dsygs2.f"
		i__2 = k - 1;
#line 270 "dsygs2.f"
		daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
#line 271 "dsygs2.f"
		i__2 = k - 1;
#line 271 "dsygs2.f"
		dsyr2_(uplo, &i__2, &c_b27, &a[k + a_dim1], lda, &b[k + 
			b_dim1], ldb, &a[a_offset], lda, (ftnlen)1);
#line 273 "dsygs2.f"
		i__2 = k - 1;
#line 273 "dsygs2.f"
		daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
#line 274 "dsygs2.f"
		i__2 = k - 1;
#line 274 "dsygs2.f"
		dscal_(&i__2, &bkk, &a[k + a_dim1], lda);
/* Computing 2nd power */
#line 275 "dsygs2.f"
		d__1 = bkk;
#line 275 "dsygs2.f"
		a[k + k * a_dim1] = akk * (d__1 * d__1);
#line 276 "dsygs2.f"
/* L40: */
#line 276 "dsygs2.f"
	    }
#line 277 "dsygs2.f"
	}
#line 278 "dsygs2.f"
    }
#line 279 "dsygs2.f"
    return 0;

/*     End of DSYGS2 */

} /* dsygs2_ */

