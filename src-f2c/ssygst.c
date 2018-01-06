#line 1 "ssygst.f"
/* ssygst.f -- translated by f2c (version 20100827).
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

#line 1 "ssygst.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 1.;
static doublereal c_b16 = -.5;
static doublereal c_b19 = -1.;
static doublereal c_b52 = .5;

/* > \brief \b SSYGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSYGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYGST reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L. */
/* > */
/* > B must have been previously factorized as U**T*U or L*L**T by SPOTRF. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,N) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          as returned by SPOTRF. */
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

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssygst_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, kb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), ssymm_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), strsm_(char *, char *, char *, char *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), ssygs2_(
	    integer *, char *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, ftnlen), ssyr2k_(char *, char *, integer 
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

#line 168 "ssygst.f"
    /* Parameter adjustments */
#line 168 "ssygst.f"
    a_dim1 = *lda;
#line 168 "ssygst.f"
    a_offset = 1 + a_dim1;
#line 168 "ssygst.f"
    a -= a_offset;
#line 168 "ssygst.f"
    b_dim1 = *ldb;
#line 168 "ssygst.f"
    b_offset = 1 + b_dim1;
#line 168 "ssygst.f"
    b -= b_offset;
#line 168 "ssygst.f"

#line 168 "ssygst.f"
    /* Function Body */
#line 168 "ssygst.f"
    *info = 0;
#line 169 "ssygst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 170 "ssygst.f"
    if (*itype < 1 || *itype > 3) {
#line 171 "ssygst.f"
	*info = -1;
#line 172 "ssygst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 173 "ssygst.f"
	*info = -2;
#line 174 "ssygst.f"
    } else if (*n < 0) {
#line 175 "ssygst.f"
	*info = -3;
#line 176 "ssygst.f"
    } else if (*lda < max(1,*n)) {
#line 177 "ssygst.f"
	*info = -5;
#line 178 "ssygst.f"
    } else if (*ldb < max(1,*n)) {
#line 179 "ssygst.f"
	*info = -7;
#line 180 "ssygst.f"
    }
#line 181 "ssygst.f"
    if (*info != 0) {
#line 182 "ssygst.f"
	i__1 = -(*info);
#line 182 "ssygst.f"
	xerbla_("SSYGST", &i__1, (ftnlen)6);
#line 183 "ssygst.f"
	return 0;
#line 184 "ssygst.f"
    }

/*     Quick return if possible */

#line 188 "ssygst.f"
    if (*n == 0) {
#line 188 "ssygst.f"
	return 0;
#line 188 "ssygst.f"
    }

/*     Determine the block size for this environment. */

#line 193 "ssygst.f"
    nb = ilaenv_(&c__1, "SSYGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 195 "ssygst.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 199 "ssygst.f"
	ssygs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);
#line 200 "ssygst.f"
    } else {

/*        Use blocked code */

#line 204 "ssygst.f"
	if (*itype == 1) {
#line 205 "ssygst.f"
	    if (upper) {

/*              Compute inv(U**T)*A*inv(U) */

#line 209 "ssygst.f"
		i__1 = *n;
#line 209 "ssygst.f"
		i__2 = nb;
#line 209 "ssygst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 210 "ssygst.f"
		    i__3 = *n - k + 1;
#line 210 "ssygst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n) */

#line 214 "ssygst.f"
		    ssygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 216 "ssygst.f"
		    if (k + kb <= *n) {
#line 217 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 217 "ssygst.f"
			strsm_("Left", uplo, "Transpose", "Non-unit", &kb, &
				i__3, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				(k + kb) * a_dim1], lda, (ftnlen)4, (ftnlen)1,
				 (ftnlen)9, (ftnlen)8);
#line 220 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 220 "ssygst.f"
			ssymm_("Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 223 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 223 "ssygst.f"
			ssyr2k_(uplo, "Transpose", &i__3, &kb, &c_b19, &a[k + 
				(k + kb) * a_dim1], lda, &b[k + (k + kb) * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)9);
#line 226 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 226 "ssygst.f"
			ssymm_("Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 229 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 229 "ssygst.f"
			strsm_("Right", uplo, "No transpose", "Non-unit", &kb,
				 &i__3, &c_b14, &b[k + kb + (k + kb) * b_dim1]
				, ldb, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 233 "ssygst.f"
		    }
#line 234 "ssygst.f"
/* L10: */
#line 234 "ssygst.f"
		}
#line 235 "ssygst.f"
	    } else {

/*              Compute inv(L)*A*inv(L**T) */

#line 239 "ssygst.f"
		i__2 = *n;
#line 239 "ssygst.f"
		i__1 = nb;
#line 239 "ssygst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 240 "ssygst.f"
		    i__3 = *n - k + 1;
#line 240 "ssygst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n) */

#line 244 "ssygst.f"
		    ssygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 246 "ssygst.f"
		    if (k + kb <= *n) {
#line 247 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 247 "ssygst.f"
			strsm_("Right", uplo, "Transpose", "Non-unit", &i__3, 
				&kb, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				kb + k * a_dim1], lda, (ftnlen)5, (ftnlen)1, (
				ftnlen)9, (ftnlen)8);
#line 250 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 250 "ssygst.f"
			ssymm_("Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
#line 253 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 253 "ssygst.f"
			ssyr2k_(uplo, "No transpose", &i__3, &kb, &c_b19, &a[
				k + kb + k * a_dim1], lda, &b[k + kb + k * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)12);
#line 256 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 256 "ssygst.f"
			ssymm_("Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
#line 259 "ssygst.f"
			i__3 = *n - k - kb + 1;
#line 259 "ssygst.f"
			strsm_("Left", uplo, "No transpose", "Non-unit", &
				i__3, &kb, &c_b14, &b[k + kb + (k + kb) * 
				b_dim1], ldb, &a[k + kb + k * a_dim1], lda, (
				ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 263 "ssygst.f"
		    }
#line 264 "ssygst.f"
/* L20: */
#line 264 "ssygst.f"
		}
#line 265 "ssygst.f"
	    }
#line 266 "ssygst.f"
	} else {
#line 267 "ssygst.f"
	    if (upper) {

/*              Compute U*A*U**T */

#line 271 "ssygst.f"
		i__1 = *n;
#line 271 "ssygst.f"
		i__2 = nb;
#line 271 "ssygst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 272 "ssygst.f"
		    i__3 = *n - k + 1;
#line 272 "ssygst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */

#line 276 "ssygst.f"
		    i__3 = k - 1;
#line 276 "ssygst.f"
		    strmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
			    kb, &c_b14, &b[b_offset], ldb, &a[k * a_dim1 + 1],
			     lda, (ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8)
			    ;
#line 278 "ssygst.f"
		    i__3 = k - 1;
#line 278 "ssygst.f"
		    ssymm_("Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 280 "ssygst.f"
		    i__3 = k - 1;
#line 280 "ssygst.f"
		    ssyr2k_(uplo, "No transpose", &i__3, &kb, &c_b14, &a[k * 
			    a_dim1 + 1], lda, &b[k * b_dim1 + 1], ldb, &c_b14,
			     &a[a_offset], lda, (ftnlen)1, (ftnlen)12);
#line 283 "ssygst.f"
		    i__3 = k - 1;
#line 283 "ssygst.f"
		    ssymm_("Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 285 "ssygst.f"
		    i__3 = k - 1;
#line 285 "ssygst.f"
		    strmm_("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
			     &c_b14, &b[k + k * b_dim1], ldb, &a[k * a_dim1 + 
			    1], lda, (ftnlen)5, (ftnlen)1, (ftnlen)9, (ftnlen)
			    8);
#line 288 "ssygst.f"
		    ssygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 290 "ssygst.f"
/* L30: */
#line 290 "ssygst.f"
		}
#line 291 "ssygst.f"
	    } else {

/*              Compute L**T*A*L */

#line 295 "ssygst.f"
		i__2 = *n;
#line 295 "ssygst.f"
		i__1 = nb;
#line 295 "ssygst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 296 "ssygst.f"
		    i__3 = *n - k + 1;
#line 296 "ssygst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */

#line 300 "ssygst.f"
		    i__3 = k - 1;
#line 300 "ssygst.f"
		    strmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
			    i__3, &c_b14, &b[b_offset], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 302 "ssygst.f"
		    i__3 = k - 1;
#line 302 "ssygst.f"
		    ssymm_("Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
#line 304 "ssygst.f"
		    i__3 = k - 1;
#line 304 "ssygst.f"
		    ssyr2k_(uplo, "Transpose", &i__3, &kb, &c_b14, &a[k + 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[
			    a_offset], lda, (ftnlen)1, (ftnlen)9);
#line 307 "ssygst.f"
		    i__3 = k - 1;
#line 307 "ssygst.f"
		    ssymm_("Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
#line 309 "ssygst.f"
		    i__3 = k - 1;
#line 309 "ssygst.f"
		    strmm_("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
			    &c_b14, &b[k + k * b_dim1], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)4, (ftnlen)1, (ftnlen)9, (ftnlen)8);
#line 311 "ssygst.f"
		    ssygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 313 "ssygst.f"
/* L40: */
#line 313 "ssygst.f"
		}
#line 314 "ssygst.f"
	    }
#line 315 "ssygst.f"
	}
#line 316 "ssygst.f"
    }
#line 317 "ssygst.f"
    return 0;

/*     End of SSYGST */

} /* ssygst_ */

