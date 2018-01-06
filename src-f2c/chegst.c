#line 1 "chegst.f"
/* chegst.f -- translated by f2c (version 20100827).
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

#line 1 "chegst.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {.5,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 1.;

/* > \brief \b CHEGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHEGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGST reduces a complex Hermitian-definite generalized */
/* > eigenproblem to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L. */
/* > */
/* > B must have been previously factorized as U**H*U or L*L**H by CPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); */
/* >          = 2 or 3: compute U*A*U**H or L**H*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored and B is factored as */
/* >                  U**H*U; */
/* >          = 'L':  Lower triangle of A is stored and B is factored as */
/* >                  L*L**H. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
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
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          as returned by CPOTRF. */
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

/* > \ingroup complexHEcomputational */

/*  ===================================================================== */
/* Subroutine */ int chegst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer k, kb, nb;
    extern /* Subroutine */ int chemm_(char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    ctrsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int chegs2_(integer *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen), cher2k_(char *, char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 171 "chegst.f"
    /* Parameter adjustments */
#line 171 "chegst.f"
    a_dim1 = *lda;
#line 171 "chegst.f"
    a_offset = 1 + a_dim1;
#line 171 "chegst.f"
    a -= a_offset;
#line 171 "chegst.f"
    b_dim1 = *ldb;
#line 171 "chegst.f"
    b_offset = 1 + b_dim1;
#line 171 "chegst.f"
    b -= b_offset;
#line 171 "chegst.f"

#line 171 "chegst.f"
    /* Function Body */
#line 171 "chegst.f"
    *info = 0;
#line 172 "chegst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 173 "chegst.f"
    if (*itype < 1 || *itype > 3) {
#line 174 "chegst.f"
	*info = -1;
#line 175 "chegst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 176 "chegst.f"
	*info = -2;
#line 177 "chegst.f"
    } else if (*n < 0) {
#line 178 "chegst.f"
	*info = -3;
#line 179 "chegst.f"
    } else if (*lda < max(1,*n)) {
#line 180 "chegst.f"
	*info = -5;
#line 181 "chegst.f"
    } else if (*ldb < max(1,*n)) {
#line 182 "chegst.f"
	*info = -7;
#line 183 "chegst.f"
    }
#line 184 "chegst.f"
    if (*info != 0) {
#line 185 "chegst.f"
	i__1 = -(*info);
#line 185 "chegst.f"
	xerbla_("CHEGST", &i__1, (ftnlen)6);
#line 186 "chegst.f"
	return 0;
#line 187 "chegst.f"
    }

/*     Quick return if possible */

#line 191 "chegst.f"
    if (*n == 0) {
#line 191 "chegst.f"
	return 0;
#line 191 "chegst.f"
    }

/*     Determine the block size for this environment. */

#line 196 "chegst.f"
    nb = ilaenv_(&c__1, "CHEGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 198 "chegst.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

#line 202 "chegst.f"
	chegs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);
#line 203 "chegst.f"
    } else {

/*        Use blocked code */

#line 207 "chegst.f"
	if (*itype == 1) {
#line 208 "chegst.f"
	    if (upper) {

/*              Compute inv(U**H)*A*inv(U) */

#line 212 "chegst.f"
		i__1 = *n;
#line 212 "chegst.f"
		i__2 = nb;
#line 212 "chegst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 213 "chegst.f"
		    i__3 = *n - k + 1;
#line 213 "chegst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n) */

#line 217 "chegst.f"
		    chegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 219 "chegst.f"
		    if (k + kb <= *n) {
#line 220 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 220 "chegst.f"
			ctrsm_("Left", uplo, "Conjugate transpose", "Non-unit"
				, &kb, &i__3, &c_b1, &b[k + k * b_dim1], ldb, 
				&a[k + (k + kb) * a_dim1], lda, (ftnlen)4, (
				ftnlen)1, (ftnlen)19, (ftnlen)8);
#line 223 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 223 "chegst.f"
			z__1.r = -.5, z__1.i = -0.;
#line 223 "chegst.f"
			chemm_("Left", uplo, &kb, &i__3, &z__1, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b1, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 226 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 226 "chegst.f"
			z__1.r = -1., z__1.i = -0.;
#line 226 "chegst.f"
			cher2k_(uplo, "Conjugate transpose", &i__3, &kb, &
				z__1, &a[k + (k + kb) * a_dim1], lda, &b[k + (
				k + kb) * b_dim1], ldb, &c_b18, &a[k + kb + (
				k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)19)
				;
#line 230 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 230 "chegst.f"
			z__1.r = -.5, z__1.i = -0.;
#line 230 "chegst.f"
			chemm_("Left", uplo, &kb, &i__3, &z__1, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b1, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
#line 233 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 233 "chegst.f"
			ctrsm_("Right", uplo, "No transpose", "Non-unit", &kb,
				 &i__3, &c_b1, &b[k + kb + (k + kb) * b_dim1],
				 ldb, &a[k + (k + kb) * a_dim1], lda, (ftnlen)
				5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 237 "chegst.f"
		    }
#line 238 "chegst.f"
/* L10: */
#line 238 "chegst.f"
		}
#line 239 "chegst.f"
	    } else {

/*              Compute inv(L)*A*inv(L**H) */

#line 243 "chegst.f"
		i__2 = *n;
#line 243 "chegst.f"
		i__1 = nb;
#line 243 "chegst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 244 "chegst.f"
		    i__3 = *n - k + 1;
#line 244 "chegst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n) */

#line 248 "chegst.f"
		    chegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 250 "chegst.f"
		    if (k + kb <= *n) {
#line 251 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 251 "chegst.f"
			ctrsm_("Right", uplo, "Conjugate transpose", "Non-un"\
				"it", &i__3, &kb, &c_b1, &b[k + k * b_dim1], 
				ldb, &a[k + kb + k * a_dim1], lda, (ftnlen)5, 
				(ftnlen)1, (ftnlen)19, (ftnlen)8);
#line 254 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 254 "chegst.f"
			z__1.r = -.5, z__1.i = -0.;
#line 254 "chegst.f"
			chemm_("Right", uplo, &i__3, &kb, &z__1, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b1, &a[k + kb + k * a_dim1], lda, (ftnlen)5,
				 (ftnlen)1);
#line 257 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 257 "chegst.f"
			z__1.r = -1., z__1.i = -0.;
#line 257 "chegst.f"
			cher2k_(uplo, "No transpose", &i__3, &kb, &z__1, &a[k 
				+ kb + k * a_dim1], lda, &b[k + kb + k * 
				b_dim1], ldb, &c_b18, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)12);
#line 261 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 261 "chegst.f"
			z__1.r = -.5, z__1.i = -0.;
#line 261 "chegst.f"
			chemm_("Right", uplo, &i__3, &kb, &z__1, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b1, &a[k + kb + k * a_dim1], lda, (ftnlen)5,
				 (ftnlen)1);
#line 264 "chegst.f"
			i__3 = *n - k - kb + 1;
#line 264 "chegst.f"
			ctrsm_("Left", uplo, "No transpose", "Non-unit", &
				i__3, &kb, &c_b1, &b[k + kb + (k + kb) * 
				b_dim1], ldb, &a[k + kb + k * a_dim1], lda, (
				ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 268 "chegst.f"
		    }
#line 269 "chegst.f"
/* L20: */
#line 269 "chegst.f"
		}
#line 270 "chegst.f"
	    }
#line 271 "chegst.f"
	} else {
#line 272 "chegst.f"
	    if (upper) {

/*              Compute U*A*U**H */

#line 276 "chegst.f"
		i__1 = *n;
#line 276 "chegst.f"
		i__2 = nb;
#line 276 "chegst.f"
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
#line 277 "chegst.f"
		    i__3 = *n - k + 1;
#line 277 "chegst.f"
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */

#line 281 "chegst.f"
		    i__3 = k - 1;
#line 281 "chegst.f"
		    ctrmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
			    kb, &c_b1, &b[b_offset], ldb, &a[k * a_dim1 + 1], 
			    lda, (ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 283 "chegst.f"
		    i__3 = k - 1;
#line 283 "chegst.f"
		    chemm_("Right", uplo, &i__3, &kb, &c_b2, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b1, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 286 "chegst.f"
		    i__3 = k - 1;
#line 286 "chegst.f"
		    cher2k_(uplo, "No transpose", &i__3, &kb, &c_b1, &a[k * 
			    a_dim1 + 1], lda, &b[k * b_dim1 + 1], ldb, &c_b18,
			     &a[a_offset], lda, (ftnlen)1, (ftnlen)12);
#line 289 "chegst.f"
		    i__3 = k - 1;
#line 289 "chegst.f"
		    chemm_("Right", uplo, &i__3, &kb, &c_b2, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b1, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
#line 292 "chegst.f"
		    i__3 = k - 1;
#line 292 "chegst.f"
		    ctrmm_("Right", uplo, "Conjugate transpose", "Non-unit", &
			    i__3, &kb, &c_b1, &b[k + k * b_dim1], ldb, &a[k * 
			    a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1, (ftnlen)
			    19, (ftnlen)8);
#line 295 "chegst.f"
		    chegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 297 "chegst.f"
/* L30: */
#line 297 "chegst.f"
		}
#line 298 "chegst.f"
	    } else {

/*              Compute L**H*A*L */

#line 302 "chegst.f"
		i__2 = *n;
#line 302 "chegst.f"
		i__1 = nb;
#line 302 "chegst.f"
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
#line 303 "chegst.f"
		    i__3 = *n - k + 1;
#line 303 "chegst.f"
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */

#line 307 "chegst.f"
		    i__3 = k - 1;
#line 307 "chegst.f"
		    ctrmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
			    i__3, &c_b1, &b[b_offset], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
#line 309 "chegst.f"
		    i__3 = k - 1;
#line 309 "chegst.f"
		    chemm_("Left", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1]
			    , lda, &b[k + b_dim1], ldb, &c_b1, &a[k + a_dim1],
			     lda, (ftnlen)4, (ftnlen)1);
#line 312 "chegst.f"
		    i__3 = k - 1;
#line 312 "chegst.f"
		    cher2k_(uplo, "Conjugate transpose", &i__3, &kb, &c_b1, &
			    a[k + a_dim1], lda, &b[k + b_dim1], ldb, &c_b18, &
			    a[a_offset], lda, (ftnlen)1, (ftnlen)19);
#line 315 "chegst.f"
		    i__3 = k - 1;
#line 315 "chegst.f"
		    chemm_("Left", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1]
			    , lda, &b[k + b_dim1], ldb, &c_b1, &a[k + a_dim1],
			     lda, (ftnlen)4, (ftnlen)1);
#line 318 "chegst.f"
		    i__3 = k - 1;
#line 318 "chegst.f"
		    ctrmm_("Left", uplo, "Conjugate transpose", "Non-unit", &
			    kb, &i__3, &c_b1, &b[k + k * b_dim1], ldb, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1, (ftnlen)19, (
			    ftnlen)8);
#line 321 "chegst.f"
		    chegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
#line 323 "chegst.f"
/* L40: */
#line 323 "chegst.f"
		}
#line 324 "chegst.f"
	    }
#line 325 "chegst.f"
	}
#line 326 "chegst.f"
    }
#line 327 "chegst.f"
    return 0;

/*     End of CHEGST */

} /* chegst_ */

