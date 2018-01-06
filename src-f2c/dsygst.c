#line 1 "dsygst.f"
/* dsygst.f -- translated by f2c (version 20100827).
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

#line 1 "dsygst.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 1.;
static doublereal c_b16 = -.5;
static doublereal c_b19 = -1.;
static doublereal c_b52 = .5;

/* > \brief \b DSYGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

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
/* > DSYGST reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L. */
/* > */
/* > B must have been previously factorized as U**T*U or L*L**T by DPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); */
/* >          = 2 or 3: compute U*A*U**T or L**T*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored and B is factored as */
/* >                  U**T*U; */
/* >          = 'L':  Lower triangle of A is stored and B is factored as */
/* >                  L*L**T. */
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
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsygst_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, kb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsymm_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsygs2_(
	    integer *, char *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, ftnlen), dsyr2k_(char *, char *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen)
	    , xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 168 "dsygst.f"
    /* Parameter adjustments */
#line 168 "dsygst.f"
    a_dim1 = *lda;
#line 168 "dsygst.f"
    a_offset = 1 + a_dim1;
#line 168 "dsygst.f"
    a -= a_offset;
#line 168 "dsygst.f"
    b_dim1 = *ldb;
#line 168 "dsygst.f"
    b_offset = 1 + b_dim1;
#line 168 "dsygst.f"
    b -= b_offset;
#line 168 "dsygst.f"

#line 168 "dsygst.f"
    /* Function Body */
#line 168 "dsygst.f"
    *info = 0;
#line 169 "dsygst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "dsygst.f"
    if (*itype < 1 || *itype > 3) {
#line 171 "dsygst.f"
	*info = -1;
#line 172 "dsygst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "dsygst.f"
	*info = -2;
#line 174 "dsygst.f"
    } else if (*n < 0) {
#line 175 "dsygst.f"
	*info = -3;
#line 176 "dsygst.f"
    } else if (*lda < max(1,*n)) {
#line 177 "dsygst.f"
	*info = -5;
#line 178 "dsygst.f"
    } else if (*ldb < max(1,*n)) {
#line 179 "dsygst.f"
	*info = -7;
#line 180 "dsygst.f"
    }
#line 181 "dsygst.f"
    if (*info != 0) {
#line 182 "dsygst.f"
	i__1 = -(*info);
#line 182 "dsygst.f"
	xerbla_("DSYGST", &i__1, (ftnlen)6);
#line 183 "dsygst.f"
	return 0;
#line 184 "dsygst.f"
    }

/*     Quick return if possible */

#line 188 "dsygst.f"
    if (*n == 0) {
#line 188 "dsygst.f"
	return 0;
#line 188 "dsygst.f"
    }

/*     Determine the block size for this environment. */

#line 193 "dsygst.f"
    nb = ilaenv_(&c__1, "DSYGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 195 "dsygst.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 199 "dsygst.f"
	dsygs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);
#line 200 "dsygst.f"
    } else {

/*        Use blocked code */

#line 204 "dsygst.f"
	if (*itype == 1) {
#line 205 "dsygst.f"
	    if (upper) {

/*              Compute inv(U**T)*A*inv(U) */

#line 209 "dsygst.f"
		i__1 = *n;
#line 209 "dsygst.f"
		i__2 = nb;
#line 209 "dsygst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 210 "dsygst.f"
		    i__3 = *n - k + 1;
#line 210 "dsygst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n) */

#line 214 "dsygst.f"
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 216 "dsygst.f"
		    if (k + kb <= *n) {
#line 217 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 217 "dsygst.f"
			dtrsm_("Left", uplo, "Transpose", "Non-unit", &kb, &
				i__3, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				(k + kb) * a_dim1], lda, (ftnlen)4, (ftnlen)1,
				 (ftnlen)9, (ftnlen)8);
#line 220 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 220 "dsygst.f"
			dsymm_("Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 223 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 223 "dsygst.f"
			dsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b19, &a[k + 
				(k + kb) * a_dim1], lda, &b[k + (k + kb) * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)9);
#line 226 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 226 "dsygst.f"
			dsymm_("Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 229 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 229 "dsygst.f"
			dtrsm_("Right", uplo, "No transpose", "Non-unit", &kb,
				 &i__3, &c_b14, &b[k + kb + (k + kb) * b_dim1]
				, ldb, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 233 "dsygst.f"
		    }
#line 234 "dsygst.f"
/* L10: */
#line 234 "dsygst.f"
		}
#line 235 "dsygst.f"
	    } else {

/*              Compute inv(L)*A*inv(L**T) */

#line 239 "dsygst.f"
		i__2 = *n;
#line 239 "dsygst.f"
		i__1 = nb;
#line 239 "dsygst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 240 "dsygst.f"
		    i__3 = *n - k + 1;
#line 240 "dsygst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n) */

#line 244 "dsygst.f"
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 246 "dsygst.f"
		    if (k + kb <= *n) {
#line 247 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 247 "dsygst.f"
			dtrsm_("Right", uplo, "Transpose", "Non-unit", &i__3, 
				&kb, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				kb + k * a_dim1], lda, (ftnlen)5, (ftnlen)1, (
				ftnlen)9, (ftnlen)8);
#line 250 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 250 "dsygst.f"
			dsymm_("Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
#line 253 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 253 "dsygst.f"
			dsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b19, &a[
				k + kb + k * a_dim1], lda, &b[k + kb + k * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)12);
#line 256 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 256 "dsygst.f"
			dsymm_("Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
#line 259 "dsygst.f"
			i__3 = *n - k - kb + 1;
#line 259 "dsygst.f"
			dtrsm_("Left", uplo, "No transpose", "Non-unit", &
				i__3, &kb, &c_b14, &b[k + kb + (k + kb) * 
				b_dim1], ldb, &a[k + kb + k * a_dim1], lda, (
				ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 263 "dsygst.f"
		    }
#line 264 "dsygst.f"
/* L20: */
#line 264 "dsygst.f"
		}
#line 265 "dsygst.f"
	    }
#line 266 "dsygst.f"
	} else {
#line 267 "dsygst.f"
	    if (upper) {

/*              Compute U*A*U**T */

#line 271 "dsygst.f"
		i__1 = *n;
#line 271 "dsygst.f"
		i__2 = nb;
#line 271 "dsygst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 272 "dsygst.f"
		    i__3 = *n - k + 1;
#line 272 "dsygst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */

#line 276 "dsygst.f"
		    i__3 = k - 1;
#line 276 "dsygst.f"
		    dtrmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
			    kb, &c_b14, &b[b_offset], ldb, &a[k * a_dim1 + 1],
			     lda, (ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8)
			    ;
#line 278 "dsygst.f"
		    i__3 = k - 1;
#line 278 "dsygst.f"
		    dsymm_("Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 280 "dsygst.f"
		    i__3 = k - 1;
#line 280 "dsygst.f"
		    dsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b14, &a[k * 
			    a_dim1 + 1], lda, &b[k * b_dim1 + 1], ldb, &c_b14,
			     &a[a_offset], lda, (ftnlen)1, (ftnlen)12);
#line 283 "dsygst.f"
		    i__3 = k - 1;
#line 283 "dsygst.f"
		    dsymm_("Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 285 "dsygst.f"
		    i__3 = k - 1;
#line 285 "dsygst.f"
		    dtrmm_("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
			     &c_b14, &b[k + k * b_dim1], ldb, &a[k * a_dim1 + 
			    1], lda, (ftnlen)5, (ftnlen)1, (ftnlen)9, (ftnlen)
			    8);
#line 288 "dsygst.f"
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 290 "dsygst.f"
/* L30: */
#line 290 "dsygst.f"
		}
#line 291 "dsygst.f"
	    } else {

/*              Compute L**T*A*L */

#line 295 "dsygst.f"
		i__2 = *n;
#line 295 "dsygst.f"
		i__1 = nb;
#line 295 "dsygst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 296 "dsygst.f"
		    i__3 = *n - k + 1;
#line 296 "dsygst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */

#line 300 "dsygst.f"
		    i__3 = k - 1;
#line 300 "dsygst.f"
		    dtrmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
			    i__3, &c_b14, &b[b_offset], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 302 "dsygst.f"
		    i__3 = k - 1;
#line 302 "dsygst.f"
		    dsymm_("Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
#line 304 "dsygst.f"
		    i__3 = k - 1;
#line 304 "dsygst.f"
		    dsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b14, &a[k + 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[
			    a_offset], lda, (ftnlen)1, (ftnlen)9);
#line 307 "dsygst.f"
		    i__3 = k - 1;
#line 307 "dsygst.f"
		    dsymm_("Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
#line 309 "dsygst.f"
		    i__3 = k - 1;
#line 309 "dsygst.f"
		    dtrmm_("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
			    &c_b14, &b[k + k * b_dim1], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)4, (ftnlen)1, (ftnlen)9, (ftnlen)8);
#line 311 "dsygst.f"
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 313 "dsygst.f"
/* L40: */
#line 313 "dsygst.f"
		}
#line 314 "dsygst.f"
	    }
#line 315 "dsygst.f"
	}
#line 316 "dsygst.f"
    }
#line 317 "dsygst.f"
    return 0;

/*     End of DSYGST */

} /* dsygst_ */

