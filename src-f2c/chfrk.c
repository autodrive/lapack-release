#line 1 "chfrk.f"
/* chfrk.f -- translated by f2c (version 20100827).
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

#line 1 "chfrk.f"
/* > \brief \b CHFRK performs a Hermitian rank-k operation for matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHFRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chfrk.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chfrk.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chfrk.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/*                         C ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            K, LDA, N */
/*       CHARACTER          TRANS, TRANSR, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > CHFRK performs one of the Hermitian rank--k operations */
/* > */
/* >    C := alpha*A*A**H + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**H*A + beta*C, */
/* > */
/* > where alpha and beta are real scalars, C is an n--by--n Hermitian */
/* > matrix and A is an n--by--k matrix in the first case and a k--by--n */
/* > matrix in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSR */
/* > \verbatim */
/* >          TRANSR is CHARACTER*1 */
/* >          = 'N':  The Normal Form of RFP A is stored; */
/* >          = 'C':  The Conjugate-transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of the  array  C  is to be  referenced  as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry,  TRANS  specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the order of the matrix C.  N must be */
/* >           at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/* >           of  columns   of  the   matrix   A,   and  on   entry   with */
/* >           TRANS = 'C' or 'c',  K  specifies  the number of rows of the */
/* >           matrix A.  K must be at least zero. */
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
/* >          A is COMPLEX array, dimension (LDA,ka) */
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
/* >          C is COMPLEX array, dimension (N*(N+1)/2) */
/* >           On entry, the matrix A in RFP Format. RFP Format is */
/* >           described by TRANSR, UPLO and N. Note that the imaginary */
/* >           parts of the diagonal elements need not be set, they are */
/* >           assumed to be zero, and on exit they are set to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int chfrk_(char *transr, char *uplo, char *trans, integer *n,
	 integer *k, doublereal *alpha, doublecomplex *a, integer *lda, 
	doublereal *beta, doublecomplex *c__, ftnlen transr_len, ftnlen 
	uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublecomplex z__1;

    /* Local variables */
    static integer j, n1, n2, nk, info;
    static doublecomplex cbeta;
    static logical normaltransr;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cherk_(char *, char *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, doublereal *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical lower;
    static doublecomplex calpha;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 214 "chfrk.f"
    /* Parameter adjustments */
#line 214 "chfrk.f"
    a_dim1 = *lda;
#line 214 "chfrk.f"
    a_offset = 1 + a_dim1;
#line 214 "chfrk.f"
    a -= a_offset;
#line 214 "chfrk.f"
    --c__;
#line 214 "chfrk.f"

#line 214 "chfrk.f"
    /* Function Body */
#line 214 "chfrk.f"
    info = 0;
#line 215 "chfrk.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 216 "chfrk.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 217 "chfrk.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 219 "chfrk.f"
    if (notrans) {
#line 220 "chfrk.f"
	nrowa = *n;
#line 221 "chfrk.f"
    } else {
#line 222 "chfrk.f"
	nrowa = *k;
#line 223 "chfrk.f"
    }

#line 225 "chfrk.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 226 "chfrk.f"
	info = -1;
#line 227 "chfrk.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 228 "chfrk.f"
	info = -2;
#line 229 "chfrk.f"
    } else if (! notrans && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 230 "chfrk.f"
	info = -3;
#line 231 "chfrk.f"
    } else if (*n < 0) {
#line 232 "chfrk.f"
	info = -4;
#line 233 "chfrk.f"
    } else if (*k < 0) {
#line 234 "chfrk.f"
	info = -5;
#line 235 "chfrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 236 "chfrk.f"
	info = -8;
#line 237 "chfrk.f"
    }
#line 238 "chfrk.f"
    if (info != 0) {
#line 239 "chfrk.f"
	i__1 = -info;
#line 239 "chfrk.f"
	xerbla_("CHFRK ", &i__1, (ftnlen)6);
#line 240 "chfrk.f"
	return 0;
#line 241 "chfrk.f"
    }

/*     Quick return if possible. */

/*     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
/*     done (it is in CHERK for example) and left in the general case. */

#line 248 "chfrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 248 "chfrk.f"
	return 0;
#line 248 "chfrk.f"
    }

#line 251 "chfrk.f"
    if (*alpha == 0. && *beta == 0.) {
#line 252 "chfrk.f"
	i__1 = *n * (*n + 1) / 2;
#line 252 "chfrk.f"
	for (j = 1; j <= i__1; ++j) {
#line 253 "chfrk.f"
	    i__2 = j;
#line 253 "chfrk.f"
	    c__[i__2].r = 0., c__[i__2].i = 0.;
#line 254 "chfrk.f"
	}
#line 255 "chfrk.f"
	return 0;
#line 256 "chfrk.f"
    }

#line 258 "chfrk.f"
    z__1.r = *alpha, z__1.i = 0.;
#line 258 "chfrk.f"
    calpha.r = z__1.r, calpha.i = z__1.i;
#line 259 "chfrk.f"
    z__1.r = *beta, z__1.i = 0.;
#line 259 "chfrk.f"
    cbeta.r = z__1.r, cbeta.i = z__1.i;

/*     C is N-by-N. */
/*     If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*     If N is even, NISODD = .FALSE., and NK. */

#line 265 "chfrk.f"
    if (*n % 2 == 0) {
#line 266 "chfrk.f"
	nisodd = FALSE_;
#line 267 "chfrk.f"
	nk = *n / 2;
#line 268 "chfrk.f"
    } else {
#line 269 "chfrk.f"
	nisodd = TRUE_;
#line 270 "chfrk.f"
	if (lower) {
#line 271 "chfrk.f"
	    n2 = *n / 2;
#line 272 "chfrk.f"
	    n1 = *n - n2;
#line 273 "chfrk.f"
	} else {
#line 274 "chfrk.f"
	    n1 = *n / 2;
#line 275 "chfrk.f"
	    n2 = *n - n1;
#line 276 "chfrk.f"
	}
#line 277 "chfrk.f"
    }

#line 279 "chfrk.f"
    if (nisodd) {

/*        N is odd */

#line 283 "chfrk.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 287 "chfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 291 "chfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 295 "chfrk.f"
		    cherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 297 "chfrk.f"
		    cherk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1);
#line 299 "chfrk.f"
		    cgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[n1 + 1], 
			    n, (ftnlen)1, (ftnlen)1);

#line 302 "chfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */

#line 306 "chfrk.f"
		    cherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 308 "chfrk.f"
		    cherk_("U", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1)
			    ;
#line 310 "chfrk.f"
		    cgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);

#line 313 "chfrk.f"
		}

#line 315 "chfrk.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 319 "chfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 323 "chfrk.f"
		    cherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 325 "chfrk.f"
		    cherk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 327 "chfrk.f"
		    cgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n2 + a_dim1], lda, &cbeta, &c__[1], n, (
			    ftnlen)1, (ftnlen)1);

#line 330 "chfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */

#line 334 "chfrk.f"
		    cherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 336 "chfrk.f"
		    cherk_("U", "C", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 338 "chfrk.f"
		    cgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n2 * a_dim1 + 1], lda, &cbeta, &c__[1], n,
			     (ftnlen)1, (ftnlen)1);

#line 341 "chfrk.f"
		}

#line 343 "chfrk.f"
	    }

#line 345 "chfrk.f"
	} else {

/*           N is odd, and TRANSR = 'C' */

#line 349 "chfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 353 "chfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */

#line 357 "chfrk.f"
		    cherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 359 "chfrk.f"
		    cherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 361 "chfrk.f"
		    cgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n1 + 1 + a_dim1], lda, &cbeta, &c__[n1 * 
			    n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 365 "chfrk.f"
		} else {

/*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */

#line 369 "chfrk.f"
		    cherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 371 "chfrk.f"
		    cherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 373 "chfrk.f"
		    cgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(n1 + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    n1 * n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 377 "chfrk.f"
		}

#line 379 "chfrk.f"
	    } else {

/*              N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 383 "chfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */

#line 387 "chfrk.f"
		    cherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 389 "chfrk.f"
		    cherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (ftnlen)
			    1);
#line 391 "chfrk.f"
		    cgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[1], &n2, 
			    (ftnlen)1, (ftnlen)1);

#line 394 "chfrk.f"
		} else {

/*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */

#line 398 "chfrk.f"
		    cherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 400 "chfrk.f"
		    cherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (
			    ftnlen)1);
#line 402 "chfrk.f"
		    cgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[1], &n2, (ftnlen)1, (ftnlen)1);

#line 405 "chfrk.f"
		}

#line 407 "chfrk.f"
	    }

#line 409 "chfrk.f"
	}

#line 411 "chfrk.f"
    } else {

/*        N is even */

#line 415 "chfrk.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 419 "chfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 423 "chfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 427 "chfrk.f"
		    i__1 = *n + 1;
#line 427 "chfrk.f"
		    cherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 429 "chfrk.f"
		    i__1 = *n + 1;
#line 429 "chfrk.f"
		    cherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 431 "chfrk.f"
		    i__1 = *n + 1;
#line 431 "chfrk.f"
		    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[nk + 2], 
			    &i__1, (ftnlen)1, (ftnlen)1);

#line 435 "chfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */

#line 439 "chfrk.f"
		    i__1 = *n + 1;
#line 439 "chfrk.f"
		    cherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 441 "chfrk.f"
		    i__1 = *n + 1;
#line 441 "chfrk.f"
		    cherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 443 "chfrk.f"
		    i__1 = *n + 1;
#line 443 "chfrk.f"
		    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);

#line 447 "chfrk.f"
		}

#line 449 "chfrk.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 453 "chfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 457 "chfrk.f"
		    i__1 = *n + 1;
#line 457 "chfrk.f"
		    cherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 459 "chfrk.f"
		    i__1 = *n + 1;
#line 459 "chfrk.f"
		    cherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 461 "chfrk.f"
		    i__1 = *n + 1;
#line 461 "chfrk.f"
		    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[nk + 1 + a_dim1], lda, &cbeta, &c__[1], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 465 "chfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */

#line 469 "chfrk.f"
		    i__1 = *n + 1;
#line 469 "chfrk.f"
		    cherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 471 "chfrk.f"
		    i__1 = *n + 1;
#line 471 "chfrk.f"
		    cherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk + 1], &i__1, (ftnlen)1, (
			    ftnlen)1);
#line 473 "chfrk.f"
		    i__1 = *n + 1;
#line 473 "chfrk.f"
		    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    1], &i__1, (ftnlen)1, (ftnlen)1);

#line 477 "chfrk.f"
		}

#line 479 "chfrk.f"
	    }

#line 481 "chfrk.f"
	} else {

/*           N is even, and TRANSR = 'C' */

#line 485 "chfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'C', and UPLO = 'L' */

#line 489 "chfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */

#line 493 "chfrk.f"
		    cherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 495 "chfrk.f"
		    cherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 497 "chfrk.f"
		    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[nk + 1 + a_dim1], lda, &cbeta, &c__[(nk + 
			    1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 501 "chfrk.f"
		} else {

/*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */

#line 505 "chfrk.f"
		    cherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 507 "chfrk.f"
		    cherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 509 "chfrk.f"
		    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    (nk + 1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 513 "chfrk.f"
		}

#line 515 "chfrk.f"
	    } else {

/*              N is even, TRANSR = 'C', and UPLO = 'U' */

#line 519 "chfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */

#line 523 "chfrk.f"
		    cherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 525 "chfrk.f"
		    cherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 527 "chfrk.f"
		    cgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[1], &nk, 
			    (ftnlen)1, (ftnlen)1);

#line 530 "chfrk.f"
		} else {

/*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */

#line 534 "chfrk.f"
		    cherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 536 "chfrk.f"
		    cherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (
			    ftnlen)1);
#line 538 "chfrk.f"
		    cgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[1], &nk, (ftnlen)1, (ftnlen)1);

#line 541 "chfrk.f"
		}

#line 543 "chfrk.f"
	    }

#line 545 "chfrk.f"
	}

#line 547 "chfrk.f"
    }

#line 549 "chfrk.f"
    return 0;

/*     End of CHFRK */

} /* chfrk_ */

