#line 1 "ssfrk.f"
/* ssfrk.f -- translated by f2c (version 20100827).
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

#line 1 "ssfrk.f"
/* > \brief \b SSFRK performs a symmetric rank-k operation for matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSFRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssfrk.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssfrk.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssfrk.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/*                         C ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            K, LDA, N */
/*       CHARACTER          TRANS, TRANSR, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > SSFRK performs one of the symmetric rank--k operations */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION (LDA,ka) */
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
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (NT) */
/* >           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP */
/* >           Format. RFP Format is described by TRANSR, UPLO and N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssfrk_(char *transr, char *uplo, char *trans, integer *n,
	 integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *beta, doublereal *c__, ftnlen transr_len, ftnlen uplo_len,
	 ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer j, n1, n2, nk, info;
    static logical normaltransr;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer nrowa;
    static logical lower;
    extern /* Subroutine */ int ssyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical nisodd, notrans;


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

/*     Test the input parameters. */

#line 207 "ssfrk.f"
    /* Parameter adjustments */
#line 207 "ssfrk.f"
    a_dim1 = *lda;
#line 207 "ssfrk.f"
    a_offset = 1 + a_dim1;
#line 207 "ssfrk.f"
    a -= a_offset;
#line 207 "ssfrk.f"
    --c__;
#line 207 "ssfrk.f"

#line 207 "ssfrk.f"
    /* Function Body */
#line 207 "ssfrk.f"
    info = 0;
#line 208 "ssfrk.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 209 "ssfrk.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 210 "ssfrk.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 212 "ssfrk.f"
    if (notrans) {
#line 213 "ssfrk.f"
	nrowa = *n;
#line 214 "ssfrk.f"
    } else {
#line 215 "ssfrk.f"
	nrowa = *k;
#line 216 "ssfrk.f"
    }

#line 218 "ssfrk.f"
    if (! normaltransr && ! lsame_(transr, "T", (ftnlen)1, (ftnlen)1)) {
#line 219 "ssfrk.f"
	info = -1;
#line 220 "ssfrk.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 221 "ssfrk.f"
	info = -2;
#line 222 "ssfrk.f"
    } else if (! notrans && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 223 "ssfrk.f"
	info = -3;
#line 224 "ssfrk.f"
    } else if (*n < 0) {
#line 225 "ssfrk.f"
	info = -4;
#line 226 "ssfrk.f"
    } else if (*k < 0) {
#line 227 "ssfrk.f"
	info = -5;
#line 228 "ssfrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 229 "ssfrk.f"
	info = -8;
#line 230 "ssfrk.f"
    }
#line 231 "ssfrk.f"
    if (info != 0) {
#line 232 "ssfrk.f"
	i__1 = -info;
#line 232 "ssfrk.f"
	xerbla_("SSFRK ", &i__1, (ftnlen)6);
#line 233 "ssfrk.f"
	return 0;
#line 234 "ssfrk.f"
    }

/*     Quick return if possible. */

/*     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
/*     done (it is in SSYRK for example) and left in the general case. */

#line 241 "ssfrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 241 "ssfrk.f"
	return 0;
#line 241 "ssfrk.f"
    }

#line 244 "ssfrk.f"
    if (*alpha == 0. && *beta == 0.) {
#line 245 "ssfrk.f"
	i__1 = *n * (*n + 1) / 2;
#line 245 "ssfrk.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "ssfrk.f"
	    c__[j] = 0.;
#line 247 "ssfrk.f"
	}
#line 248 "ssfrk.f"
	return 0;
#line 249 "ssfrk.f"
    }

/*     C is N-by-N. */
/*     If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*     If N is even, NISODD = .FALSE., and NK. */

#line 255 "ssfrk.f"
    if (*n % 2 == 0) {
#line 256 "ssfrk.f"
	nisodd = FALSE_;
#line 257 "ssfrk.f"
	nk = *n / 2;
#line 258 "ssfrk.f"
    } else {
#line 259 "ssfrk.f"
	nisodd = TRUE_;
#line 260 "ssfrk.f"
	if (lower) {
#line 261 "ssfrk.f"
	    n2 = *n / 2;
#line 262 "ssfrk.f"
	    n1 = *n - n2;
#line 263 "ssfrk.f"
	} else {
#line 264 "ssfrk.f"
	    n1 = *n / 2;
#line 265 "ssfrk.f"
	    n2 = *n - n1;
#line 266 "ssfrk.f"
	}
#line 267 "ssfrk.f"
    }

#line 269 "ssfrk.f"
    if (nisodd) {

/*        N is odd */

#line 273 "ssfrk.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 277 "ssfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 281 "ssfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 285 "ssfrk.f"
		    ssyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 287 "ssfrk.f"
		    ssyrk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1);
#line 289 "ssfrk.f"
		    sgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1], n, (
			    ftnlen)1, (ftnlen)1);

#line 292 "ssfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */

#line 296 "ssfrk.f"
		    ssyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 298 "ssfrk.f"
		    ssyrk_("U", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1)
			    ;
#line 300 "ssfrk.f"
		    sgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1]
			    , n, (ftnlen)1, (ftnlen)1);

#line 303 "ssfrk.f"
		}

#line 305 "ssfrk.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 309 "ssfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 313 "ssfrk.f"
		    ssyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 315 "ssfrk.f"
		    ssyrk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 317 "ssfrk.f"
		    sgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n2 + a_dim1], lda, beta, &c__[1], n, (ftnlen)1,
			     (ftnlen)1);

#line 320 "ssfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */

#line 324 "ssfrk.f"
		    ssyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 326 "ssfrk.f"
		    ssyrk_("U", "T", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 328 "ssfrk.f"
		    sgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n2 * a_dim1 + 1], lda, beta, &c__[1], n, (
			    ftnlen)1, (ftnlen)1);

#line 331 "ssfrk.f"
		}

#line 333 "ssfrk.f"
	    }

#line 335 "ssfrk.f"
	} else {

/*           N is odd, and TRANSR = 'T' */

#line 339 "ssfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'T', and UPLO = 'L' */

#line 343 "ssfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */

#line 347 "ssfrk.f"
		    ssyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 349 "ssfrk.f"
		    ssyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 351 "ssfrk.f"
		    sgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[n1 + 1 + a_dim1], lda, beta, &c__[n1 * n1 + 1],
			     &n1, (ftnlen)1, (ftnlen)1);

#line 355 "ssfrk.f"
		} else {

/*                 N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */

#line 359 "ssfrk.f"
		    ssyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 361 "ssfrk.f"
		    ssyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 363 "ssfrk.f"
		    sgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[n1 * 
			    n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 367 "ssfrk.f"
		}

#line 369 "ssfrk.f"
	    } else {

/*              N is odd, TRANSR = 'T', and UPLO = 'U' */

#line 373 "ssfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */

#line 377 "ssfrk.f"
		    ssyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 379 "ssfrk.f"
		    ssyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (ftnlen)
			    1);
#line 381 "ssfrk.f"
		    sgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[1], &n2, (
			    ftnlen)1, (ftnlen)1);

#line 384 "ssfrk.f"
		} else {

/*                 N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */

#line 388 "ssfrk.f"
		    ssyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 390 "ssfrk.f"
		    ssyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (
			    ftnlen)1);
#line 392 "ssfrk.f"
		    sgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &
			    n2, (ftnlen)1, (ftnlen)1);

#line 395 "ssfrk.f"
		}

#line 397 "ssfrk.f"
	    }

#line 399 "ssfrk.f"
	}

#line 401 "ssfrk.f"
    } else {

/*        N is even */

#line 405 "ssfrk.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 409 "ssfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 413 "ssfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 417 "ssfrk.f"
		    i__1 = *n + 1;
#line 417 "ssfrk.f"
		    ssyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 419 "ssfrk.f"
		    i__1 = *n + 1;
#line 419 "ssfrk.f"
		    ssyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 421 "ssfrk.f"
		    i__1 = *n + 1;
#line 421 "ssfrk.f"
		    sgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 425 "ssfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */

#line 429 "ssfrk.f"
		    i__1 = *n + 1;
#line 429 "ssfrk.f"
		    ssyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 431 "ssfrk.f"
		    i__1 = *n + 1;
#line 431 "ssfrk.f"
		    ssyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 433 "ssfrk.f"
		    i__1 = *n + 1;
#line 433 "ssfrk.f"
		    sgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2]
			    , &i__1, (ftnlen)1, (ftnlen)1);

#line 437 "ssfrk.f"
		}

#line 439 "ssfrk.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 443 "ssfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 447 "ssfrk.f"
		    i__1 = *n + 1;
#line 447 "ssfrk.f"
		    ssyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 449 "ssfrk.f"
		    i__1 = *n + 1;
#line 449 "ssfrk.f"
		    ssyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 451 "ssfrk.f"
		    i__1 = *n + 1;
#line 451 "ssfrk.f"
		    sgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[nk + 1 + a_dim1], lda, beta, &c__[1], &i__1, (
			    ftnlen)1, (ftnlen)1);

#line 455 "ssfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */

#line 459 "ssfrk.f"
		    i__1 = *n + 1;
#line 459 "ssfrk.f"
		    ssyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 461 "ssfrk.f"
		    i__1 = *n + 1;
#line 461 "ssfrk.f"
		    ssyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk + 1], &i__1, (ftnlen)1, (
			    ftnlen)1);
#line 463 "ssfrk.f"
		    i__1 = *n + 1;
#line 463 "ssfrk.f"
		    sgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 467 "ssfrk.f"
		}

#line 469 "ssfrk.f"
	    }

#line 471 "ssfrk.f"
	} else {

/*           N is even, and TRANSR = 'T' */

#line 475 "ssfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'T', and UPLO = 'L' */

#line 479 "ssfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */

#line 483 "ssfrk.f"
		    ssyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 485 "ssfrk.f"
		    ssyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 487 "ssfrk.f"
		    sgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[nk + 1 + a_dim1], lda, beta, &c__[(nk + 1) * 
			    nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 491 "ssfrk.f"
		} else {

/*                 N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */

#line 495 "ssfrk.f"
		    ssyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 497 "ssfrk.f"
		    ssyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 499 "ssfrk.f"
		    sgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, 
			    &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[(nk + 
			    1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 503 "ssfrk.f"
		}

#line 505 "ssfrk.f"
	    } else {

/*              N is even, TRANSR = 'T', and UPLO = 'U' */

#line 509 "ssfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */

#line 513 "ssfrk.f"
		    ssyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 515 "ssfrk.f"
		    ssyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 517 "ssfrk.f"
		    sgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], 
			    lda, &a[a_dim1 + 1], lda, beta, &c__[1], &nk, (
			    ftnlen)1, (ftnlen)1);

#line 520 "ssfrk.f"
		} else {

/*                 N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */

#line 524 "ssfrk.f"
		    ssyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 526 "ssfrk.f"
		    ssyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (
			    ftnlen)1);
#line 528 "ssfrk.f"
		    sgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 
			    + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &
			    nk, (ftnlen)1, (ftnlen)1);

#line 531 "ssfrk.f"
		}

#line 533 "ssfrk.f"
	    }

#line 535 "ssfrk.f"
	}

#line 537 "ssfrk.f"
    }

#line 539 "ssfrk.f"
    return 0;

/*     End of SSFRK */

} /* ssfrk_ */

