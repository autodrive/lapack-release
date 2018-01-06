#line 1 "dsfrk.f"
/* dsfrk.f -- translated by f2c (version 20100827).
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

#line 1 "dsfrk.f"
/* > \brief \b DSFRK performs a symmetric rank-k operation for matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSFRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsfrk.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsfrk.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsfrk.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/*                         C ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            K, LDA, N */
/*       CHARACTER          TRANS, TRANSR, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > DSFRK performs one of the symmetric rank--k operations */
/* > */
/* >    C := alpha*A*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*A + beta*C, */
/* > */
/* > where alpha and beta are real scalars, C is an n--by--n symmetric */
/* > matrix and A is an n--by--k matrix in the first case and a k--by--n */
/* > matrix in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  The Normal Form of RFP A is stored; */
/* >          = 'T':  The Transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array C is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of C */
/* >                                  is to be referenced. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix C. N must be */
/* >           at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with TRANS = 'N' or 'n', K specifies the number */
/* >           of  columns of the matrix A, and on entry with TRANS = 'T' */
/* >           or 't', K specifies the number of rows of the matrix A. K */
/* >           must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,ka) */
/* >           where KA */
/* >           is K  when TRANS = 'N' or 'n', and is N otherwise. Before */
/* >           entry with TRANS = 'N' or 'n', the leading N--by--K part of */
/* >           the array A must contain the matrix A, otherwise the leading */
/* >           K--by--N part of the array A must contain the matrix A. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/* >           be at least  max( 1, k ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION */
/* >           On entry, BETA specifies the scalar beta. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (NT) */
/* >           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP */
/* >           Format. RFP Format is described by TRANSR, UPLO and N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsfrk_(char *transr, char *uplo, char *trans, integer *n,
	 integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *beta, doublereal *c__, ftnlen transr_len, ftnlen uplo_len,
	 ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer j, n1, n2, nk, info;
    static logical normaltransr;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical lower;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical nisodd, notrans;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. */
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

/*     Test the input parameters. */

#line 208 "dsfrk.f"
    /* Parameter adjustments */
#line 208 "dsfrk.f"
    a_dim1 = *lda;
#line 208 "dsfrk.f"
    a_offset = 1 + a_dim1;
#line 208 "dsfrk.f"
    a -= a_offset;
#line 208 "dsfrk.f"
    --c__;
#line 208 "dsfrk.f"

#line 208 "dsfrk.f"
    /* Function Body */
#line 208 "dsfrk.f"
    info = 0;
#line 209 "dsfrk.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 210 "dsfrk.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 211 "dsfrk.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 213 "dsfrk.f"
    if (notrans) {
#line 214 "dsfrk.f"
	nrowa = *n;
#line 215 "dsfrk.f"
    } else {
#line 216 "dsfrk.f"
	nrowa = *k;
#line 217 "dsfrk.f"
    }

#line 219 "dsfrk.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 220 "dsfrk.f"
	info = -1;
#line 221 "dsfrk.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 222 "dsfrk.f"
	info = -2;
#line 223 "dsfrk.f"
    } else if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 224 "dsfrk.f"
	info = -3;
#line 225 "dsfrk.f"
    } else if (*n < 0) {
#line 226 "dsfrk.f"
	info = -4;
#line 227 "dsfrk.f"
    } else if (*k < 0) {
#line 228 "dsfrk.f"
	info = -5;
#line 229 "dsfrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 230 "dsfrk.f"
	info = -8;
#line 231 "dsfrk.f"
    }
#line 232 "dsfrk.f"
    if (info != 0) {
#line 233 "dsfrk.f"
	i__1 = -info;
#line 233 "dsfrk.f"
	xerbla_("DSFRK ", &i__1, (ftnlen)6);
#line 234 "dsfrk.f"
	return 0;
#line 235 "dsfrk.f"
    }

/*     Quick return if possible. */

/*     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
/*     done (it is in DSYRK for example) and left in the general case. */

#line 242 "dsfrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 242 "dsfrk.f"
	return 0;
#line 242 "dsfrk.f"
    }

#line 245 "dsfrk.f"
    if (*alpha == 0. && *beta == 0.) {
#line 246 "dsfrk.f"
	i__1 = *n * (*n + 1) / 2;
#line 246 "dsfrk.f"
	for (j = 1; j <= i__1; ++j) {
#line 247 "dsfrk.f"
	    c__[j] = 0.;
#line 248 "dsfrk.f"
	}
#line 249 "dsfrk.f"
	return 0;
#line 250 "dsfrk.f"
    }

/*     C is N-by-N. */
/*     If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*     If N is even, NISODD = .FALSE., and NK. */

#line 256 "dsfrk.f"
    if (*n % 2 == 0) {
#line 257 "dsfrk.f"
	nisodd = FALSE_;
#line 258 "dsfrk.f"
	nk = *n / 2;
#line 259 "dsfrk.f"
    } else {
#line 260 "dsfrk.f"
	nisodd = TRUE_;
#line 261 "dsfrk.f"
	if (lower) {
#line 262 "dsfrk.f"
	    n2 = *n / 2;
#line 263 "dsfrk.f"
	    n1 = *n - n2;
#line 264 "dsfrk.f"
	} else {
#line 265 "dsfrk.f"
	    n1 = *n / 2;
#line 266 "dsfrk.f"
	    n2 = *n - n1;
#line 267 "dsfrk.f"
	}
#line 268 "dsfrk.f"
    }

#line 270 "dsfrk.f"
    if (nisodd) {

/*        N is odd */

#line 274 "dsfrk.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 278 "dsfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 282 "dsfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 286 "dsfrk.f"
		    dsyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 288 "dsfrk.f"
		    dsyrk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1);
#line 290 "dsfrk.f"
		    dgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1], n, (
			    ftnlen)1, (ftnlen)1);

#line 293 "dsfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */

#line 297 "dsfrk.f"
		    dsyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 299 "dsfrk.f"
		    dsyrk_("U", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1)
			    ;
#line 301 "dsfrk.f"
		    dgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1]
			    , n, (ftnlen)1, (ftnlen)1);

#line 304 "dsfrk.f"
		}

#line 306 "dsfrk.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 310 "dsfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 314 "dsfrk.f"
		    dsyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 316 "dsfrk.f"
		    dsyrk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 318 "dsfrk.f"
		    dgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n2 + a_dim1], lda, beta, &c__[1], n, (ftnlen)1,
			     (ftnlen)1);

#line 321 "dsfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */

#line 325 "dsfrk.f"
		    dsyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 327 "dsfrk.f"
		    dsyrk_("U", "T", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 329 "dsfrk.f"
		    dgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n2 * a_dim1 + 1], lda, beta, &c__[1], n, (
			    ftnlen)1, (ftnlen)1);

#line 332 "dsfrk.f"
		}

#line 334 "dsfrk.f"
	    }

#line 336 "dsfrk.f"
	} else {

/*           N is odd, and TRANSR = 'T' */

#line 340 "dsfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 344 "dsfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */

#line 348 "dsfrk.f"
		    dsyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 350 "dsfrk.f"
		    dsyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 352 "dsfrk.f"
		    dgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n1 + 1 + a_dim1], lda, beta, &c__[n1 * n1 + 1],
			     &n1, (ftnlen)1, (ftnlen)1);

#line 356 "dsfrk.f"
		} else {

/*                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */

#line 360 "dsfrk.f"
		    dsyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 362 "dsfrk.f"
		    dsyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 364 "dsfrk.f"
		    dgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[n1 * 
			    n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 368 "dsfrk.f"
		}

#line 370 "dsfrk.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 374 "dsfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */

#line 378 "dsfrk.f"
		    dsyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 380 "dsfrk.f"
		    dsyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (ftnlen)
			    1);
#line 382 "dsfrk.f"
		    dgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[1], &n2, (
			    ftnlen)1, (ftnlen)1);

#line 385 "dsfrk.f"
		} else {

/*                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */

#line 389 "dsfrk.f"
		    dsyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 391 "dsfrk.f"
		    dsyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (
			    ftnlen)1);
#line 393 "dsfrk.f"
		    dgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &
			    n2, (ftnlen)1, (ftnlen)1);

#line 396 "dsfrk.f"
		}

#line 398 "dsfrk.f"
	    }

#line 400 "dsfrk.f"
	}

#line 402 "dsfrk.f"
    } else {

/*        N is even */

#line 406 "dsfrk.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 410 "dsfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 414 "dsfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 418 "dsfrk.f"
		    i__1 = *n + 1;
#line 418 "dsfrk.f"
		    dsyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 420 "dsfrk.f"
		    i__1 = *n + 1;
#line 420 "dsfrk.f"
		    dsyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 422 "dsfrk.f"
		    i__1 = *n + 1;
#line 422 "dsfrk.f"
		    dgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 426 "dsfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */

#line 430 "dsfrk.f"
		    i__1 = *n + 1;
#line 430 "dsfrk.f"
		    dsyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 432 "dsfrk.f"
		    i__1 = *n + 1;
#line 432 "dsfrk.f"
		    dsyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 434 "dsfrk.f"
		    i__1 = *n + 1;
#line 434 "dsfrk.f"
		    dgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2]
			    , &i__1, (ftnlen)1, (ftnlen)1);

#line 438 "dsfrk.f"
		}

#line 440 "dsfrk.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 444 "dsfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 448 "dsfrk.f"
		    i__1 = *n + 1;
#line 448 "dsfrk.f"
		    dsyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 450 "dsfrk.f"
		    i__1 = *n + 1;
#line 450 "dsfrk.f"
		    dsyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 452 "dsfrk.f"
		    i__1 = *n + 1;
#line 452 "dsfrk.f"
		    dgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[nk + 1 + a_dim1], lda, beta, &c__[1], &i__1, (
			    ftnlen)1, (ftnlen)1);

#line 456 "dsfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */

#line 460 "dsfrk.f"
		    i__1 = *n + 1;
#line 460 "dsfrk.f"
		    dsyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 462 "dsfrk.f"
		    i__1 = *n + 1;
#line 462 "dsfrk.f"
		    dsyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk + 1], &i__1, (ftnlen)1, (
			    ftnlen)1);
#line 464 "dsfrk.f"
		    i__1 = *n + 1;
#line 464 "dsfrk.f"
		    dgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 468 "dsfrk.f"
		}

#line 470 "dsfrk.f"
	    }

#line 472 "dsfrk.f"
	} else {

/*           N is even, and TRANSR = 'T' */

#line 476 "dsfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 480 "dsfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */

#line 484 "dsfrk.f"
		    dsyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 486 "dsfrk.f"
		    dsyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 488 "dsfrk.f"
		    dgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[nk + 1 + a_dim1], lda, beta, &c__[(nk + 1) * 
			    nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 492 "dsfrk.f"
		} else {

/*                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */

#line 496 "dsfrk.f"
		    dsyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 498 "dsfrk.f"
		    dsyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 500 "dsfrk.f"
		    dgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[(nk + 
			    1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 504 "dsfrk.f"
		}

#line 506 "dsfrk.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 510 "dsfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */

#line 514 "dsfrk.f"
		    dsyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 516 "dsfrk.f"
		    dsyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 518 "dsfrk.f"
		    dgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[1], &nk, (
			    ftnlen)1, (ftnlen)1);

#line 521 "dsfrk.f"
		} else {

/*                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */

#line 525 "dsfrk.f"
		    dsyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 527 "dsfrk.f"
		    dsyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (
			    ftnlen)1);
#line 529 "dsfrk.f"
		    dgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &
			    nk, (ftnlen)1, (ftnlen)1);

#line 532 "dsfrk.f"
		}

#line 534 "dsfrk.f"
	    }

#line 536 "dsfrk.f"
	}

#line 538 "dsfrk.f"
    }

#line 540 "dsfrk.f"
    return 0;

/*     End of DSFRK */

} /* dsfrk_ */

