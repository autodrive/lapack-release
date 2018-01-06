#line 1 "zhfrk.f"
/* zhfrk.f -- translated by f2c (version 20100827).
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

#line 1 "zhfrk.f"
/* > \brief \b ZHFRK performs a Hermitian rank-k operation for matrix in RFP format. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHFRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhfrk.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhfrk.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhfrk.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/*                         C ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            K, LDA, N */
/*       CHARACTER          TRANS, TRANSR, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), C( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > ZHFRK performs one of the Hermitian rank--k operations */
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
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION (LDA,ka) */
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
/* >          C is COMPLEX*16 array, dimension (N*(N+1)/2) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhfrk_(char *transr, char *uplo, char *trans, integer *n,
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zherk_(char *, char *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, doublereal *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    static integer nrowa;
    static logical lower;
    static doublecomplex calpha;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

#line 213 "zhfrk.f"
    /* Parameter adjustments */
#line 213 "zhfrk.f"
    a_dim1 = *lda;
#line 213 "zhfrk.f"
    a_offset = 1 + a_dim1;
#line 213 "zhfrk.f"
    a -= a_offset;
#line 213 "zhfrk.f"
    --c__;
#line 213 "zhfrk.f"

#line 213 "zhfrk.f"
    /* Function Body */
#line 213 "zhfrk.f"
    info = 0;
#line 214 "zhfrk.f"
    normaltransr = lsame_(transr, "N", (ftnlen)1, (ftnlen)1);
#line 215 "zhfrk.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 216 "zhfrk.f"
    notrans = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

#line 218 "zhfrk.f"
    if (notrans) {
#line 219 "zhfrk.f"
	nrowa = *n;
#line 220 "zhfrk.f"
    } else {
#line 221 "zhfrk.f"
	nrowa = *k;
#line 222 "zhfrk.f"
    }

#line 224 "zhfrk.f"
    if (! normaltransr && ! lsame_(transr, "C", (ftnlen)1, (ftnlen)1)) {
#line 225 "zhfrk.f"
	info = -1;
#line 226 "zhfrk.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 227 "zhfrk.f"
	info = -2;
#line 228 "zhfrk.f"
    } else if (! notrans && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 229 "zhfrk.f"
	info = -3;
#line 230 "zhfrk.f"
    } else if (*n < 0) {
#line 231 "zhfrk.f"
	info = -4;
#line 232 "zhfrk.f"
    } else if (*k < 0) {
#line 233 "zhfrk.f"
	info = -5;
#line 234 "zhfrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 235 "zhfrk.f"
	info = -8;
#line 236 "zhfrk.f"
    }
#line 237 "zhfrk.f"
    if (info != 0) {
#line 238 "zhfrk.f"
	i__1 = -info;
#line 238 "zhfrk.f"
	xerbla_("ZHFRK ", &i__1, (ftnlen)6);
#line 239 "zhfrk.f"
	return 0;
#line 240 "zhfrk.f"
    }

/*     Quick return if possible. */

/*     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
/*     done (it is in ZHERK for example) and left in the general case. */

#line 247 "zhfrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 247 "zhfrk.f"
	return 0;
#line 247 "zhfrk.f"
    }

#line 250 "zhfrk.f"
    if (*alpha == 0. && *beta == 0.) {
#line 251 "zhfrk.f"
	i__1 = *n * (*n + 1) / 2;
#line 251 "zhfrk.f"
	for (j = 1; j <= i__1; ++j) {
#line 252 "zhfrk.f"
	    i__2 = j;
#line 252 "zhfrk.f"
	    c__[i__2].r = 0., c__[i__2].i = 0.;
#line 253 "zhfrk.f"
	}
#line 254 "zhfrk.f"
	return 0;
#line 255 "zhfrk.f"
    }

#line 257 "zhfrk.f"
    z__1.r = *alpha, z__1.i = 0.;
#line 257 "zhfrk.f"
    calpha.r = z__1.r, calpha.i = z__1.i;
#line 258 "zhfrk.f"
    z__1.r = *beta, z__1.i = 0.;
#line 258 "zhfrk.f"
    cbeta.r = z__1.r, cbeta.i = z__1.i;

/*     C is N-by-N. */
/*     If N is odd, set NISODD = .TRUE., and N1 and N2. */
/*     If N is even, NISODD = .FALSE., and NK. */

#line 264 "zhfrk.f"
    if (*n % 2 == 0) {
#line 265 "zhfrk.f"
	nisodd = FALSE_;
#line 266 "zhfrk.f"
	nk = *n / 2;
#line 267 "zhfrk.f"
    } else {
#line 268 "zhfrk.f"
	nisodd = TRUE_;
#line 269 "zhfrk.f"
	if (lower) {
#line 270 "zhfrk.f"
	    n2 = *n / 2;
#line 271 "zhfrk.f"
	    n1 = *n - n2;
#line 272 "zhfrk.f"
	} else {
#line 273 "zhfrk.f"
	    n1 = *n / 2;
#line 274 "zhfrk.f"
	    n2 = *n - n1;
#line 275 "zhfrk.f"
	}
#line 276 "zhfrk.f"
    }

#line 278 "zhfrk.f"
    if (nisodd) {

/*        N is odd */

#line 282 "zhfrk.f"
	if (normaltransr) {

/*           N is odd and TRANSR = 'N' */

#line 286 "zhfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'N', and UPLO = 'L' */

#line 290 "zhfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 294 "zhfrk.f"
		    zherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 296 "zhfrk.f"
		    zherk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1);
#line 298 "zhfrk.f"
		    zgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[n1 + 1], 
			    n, (ftnlen)1, (ftnlen)1);

#line 301 "zhfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */

#line 305 "zhfrk.f"
		    zherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], n, (ftnlen)1, (ftnlen)1);
#line 307 "zhfrk.f"
		    zherk_("U", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[*n + 1], n, (ftnlen)1, (ftnlen)1)
			    ;
#line 309 "zhfrk.f"
		    zgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);

#line 312 "zhfrk.f"
		}

#line 314 "zhfrk.f"
	    } else {

/*              N is odd, TRANSR = 'N', and UPLO = 'U' */

#line 318 "zhfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 322 "zhfrk.f"
		    zherk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 324 "zhfrk.f"
		    zherk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 326 "zhfrk.f"
		    zgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n2 + a_dim1], lda, &cbeta, &c__[1], n, (
			    ftnlen)1, (ftnlen)1);

#line 329 "zhfrk.f"
		} else {

/*                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */

#line 333 "zhfrk.f"
		    zherk_("L", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 + 1], n, (ftnlen)1, (ftnlen)1);
#line 335 "zhfrk.f"
		    zherk_("U", "C", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, 
			    beta, &c__[n1 + 1], n, (ftnlen)1, (ftnlen)1);
#line 337 "zhfrk.f"
		    zgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n2 * a_dim1 + 1], lda, &cbeta, &c__[1], n,
			     (ftnlen)1, (ftnlen)1);

#line 340 "zhfrk.f"
		}

#line 342 "zhfrk.f"
	    }

#line 344 "zhfrk.f"
	} else {

/*           N is odd, and TRANSR = 'C' */

#line 348 "zhfrk.f"
	    if (lower) {

/*              N is odd, TRANSR = 'C', and UPLO = 'L' */

#line 352 "zhfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */

#line 356 "zhfrk.f"
		    zherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 358 "zhfrk.f"
		    zherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 360 "zhfrk.f"
		    zgemm_("N", "C", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[n1 + 1 + a_dim1], lda, &cbeta, &c__[n1 * 
			    n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 364 "zhfrk.f"
		} else {

/*                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */

#line 368 "zhfrk.f"
		    zherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[1], &n1, (ftnlen)1, (ftnlen)1);
#line 370 "zhfrk.f"
		    zherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[2], &n1, (ftnlen)1, (ftnlen)1);
#line 372 "zhfrk.f"
		    zgemm_("C", "N", &n1, &n2, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(n1 + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    n1 * n1 + 1], &n1, (ftnlen)1, (ftnlen)1);

#line 376 "zhfrk.f"
		}

#line 378 "zhfrk.f"
	    } else {

/*              N is odd, TRANSR = 'C', and UPLO = 'U' */

#line 382 "zhfrk.f"
		if (notrans) {

/*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */

#line 386 "zhfrk.f"
		    zherk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 388 "zhfrk.f"
		    zherk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, 
			    beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (ftnlen)
			    1);
#line 390 "zhfrk.f"
		    zgemm_("N", "C", &n2, &n1, k, &calpha, &a[n1 + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[1], &n2, 
			    (ftnlen)1, (ftnlen)1);

#line 393 "zhfrk.f"
		} else {

/*                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */

#line 397 "zhfrk.f"
		    zherk_("U", "C", &n1, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[n2 * n2 + 1], &n2, (ftnlen)1, (ftnlen)1);
#line 399 "zhfrk.f"
		    zherk_("L", "C", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1],
			     lda, beta, &c__[n1 * n2 + 1], &n2, (ftnlen)1, (
			    ftnlen)1);
#line 401 "zhfrk.f"
		    zgemm_("C", "N", &n2, &n1, k, &calpha, &a[(n1 + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[1], &n2, (ftnlen)1, (ftnlen)1);

#line 404 "zhfrk.f"
		}

#line 406 "zhfrk.f"
	    }

#line 408 "zhfrk.f"
	}

#line 410 "zhfrk.f"
    } else {

/*        N is even */

#line 414 "zhfrk.f"
	if (normaltransr) {

/*           N is even and TRANSR = 'N' */

#line 418 "zhfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'N', and UPLO = 'L' */

#line 422 "zhfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */

#line 426 "zhfrk.f"
		    i__1 = *n + 1;
#line 426 "zhfrk.f"
		    zherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 428 "zhfrk.f"
		    i__1 = *n + 1;
#line 428 "zhfrk.f"
		    zherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 430 "zhfrk.f"
		    i__1 = *n + 1;
#line 430 "zhfrk.f"
		    zgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[nk + 2], 
			    &i__1, (ftnlen)1, (ftnlen)1);

#line 434 "zhfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C' */

#line 438 "zhfrk.f"
		    i__1 = *n + 1;
#line 438 "zhfrk.f"
		    zherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[2], &i__1, (ftnlen)1, (ftnlen)1);
#line 440 "zhfrk.f"
		    i__1 = *n + 1;
#line 440 "zhfrk.f"
		    zherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &i__1, (ftnlen)1, (ftnlen)1);
#line 442 "zhfrk.f"
		    i__1 = *n + 1;
#line 442 "zhfrk.f"
		    zgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);

#line 446 "zhfrk.f"
		}

#line 448 "zhfrk.f"
	    } else {

/*              N is even, TRANSR = 'N', and UPLO = 'U' */

#line 452 "zhfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */

#line 456 "zhfrk.f"
		    i__1 = *n + 1;
#line 456 "zhfrk.f"
		    zherk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 458 "zhfrk.f"
		    i__1 = *n + 1;
#line 458 "zhfrk.f"
		    zherk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 460 "zhfrk.f"
		    i__1 = *n + 1;
#line 460 "zhfrk.f"
		    zgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[nk + 1 + a_dim1], lda, &cbeta, &c__[1], &
			    i__1, (ftnlen)1, (ftnlen)1);

#line 464 "zhfrk.f"
		} else {

/*                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C' */

#line 468 "zhfrk.f"
		    i__1 = *n + 1;
#line 468 "zhfrk.f"
		    zherk_("L", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 2], &i__1, (ftnlen)1, (ftnlen)1);
#line 470 "zhfrk.f"
		    i__1 = *n + 1;
#line 470 "zhfrk.f"
		    zherk_("U", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk + 1], &i__1, (ftnlen)1, (
			    ftnlen)1);
#line 472 "zhfrk.f"
		    i__1 = *n + 1;
#line 472 "zhfrk.f"
		    zgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    1], &i__1, (ftnlen)1, (ftnlen)1);

#line 476 "zhfrk.f"
		}

#line 478 "zhfrk.f"
	    }

#line 480 "zhfrk.f"
	} else {

/*           N is even, and TRANSR = 'C' */

#line 484 "zhfrk.f"
	    if (lower) {

/*              N is even, TRANSR = 'C', and UPLO = 'L' */

#line 488 "zhfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N' */

#line 492 "zhfrk.f"
		    zherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 494 "zhfrk.f"
		    zherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 496 "zhfrk.f"
		    zgemm_("N", "C", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[nk + 1 + a_dim1], lda, &cbeta, &c__[(nk + 
			    1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 500 "zhfrk.f"
		} else {

/*                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C' */

#line 504 "zhfrk.f"
		    zherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk + 1], &nk, (ftnlen)1, (ftnlen)1);
#line 506 "zhfrk.f"
		    zherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[1], &nk, (ftnlen)1, (ftnlen)1);
#line 508 "zhfrk.f"
		    zgemm_("C", "N", &nk, &nk, k, &calpha, &a[a_dim1 + 1], 
			    lda, &a[(nk + 1) * a_dim1 + 1], lda, &cbeta, &c__[
			    (nk + 1) * nk + 1], &nk, (ftnlen)1, (ftnlen)1);

#line 512 "zhfrk.f"
		}

#line 514 "zhfrk.f"
	    } else {

/*              N is even, TRANSR = 'C', and UPLO = 'U' */

#line 518 "zhfrk.f"
		if (notrans) {

/*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N' */

#line 522 "zhfrk.f"
		    zherk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 524 "zhfrk.f"
		    zherk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, 
			    beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 526 "zhfrk.f"
		    zgemm_("N", "C", &nk, &nk, k, &calpha, &a[nk + 1 + a_dim1]
			    , lda, &a[a_dim1 + 1], lda, &cbeta, &c__[1], &nk, 
			    (ftnlen)1, (ftnlen)1);

#line 529 "zhfrk.f"
		} else {

/*                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C' */

#line 533 "zhfrk.f"
		    zherk_("U", "C", &nk, k, alpha, &a[a_dim1 + 1], lda, beta,
			     &c__[nk * (nk + 1) + 1], &nk, (ftnlen)1, (ftnlen)
			    1);
#line 535 "zhfrk.f"
		    zherk_("L", "C", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1],
			     lda, beta, &c__[nk * nk + 1], &nk, (ftnlen)1, (
			    ftnlen)1);
#line 537 "zhfrk.f"
		    zgemm_("C", "N", &nk, &nk, k, &calpha, &a[(nk + 1) * 
			    a_dim1 + 1], lda, &a[a_dim1 + 1], lda, &cbeta, &
			    c__[1], &nk, (ftnlen)1, (ftnlen)1);

#line 540 "zhfrk.f"
		}

#line 542 "zhfrk.f"
	    }

#line 544 "zhfrk.f"
	}

#line 546 "zhfrk.f"
    }

#line 548 "zhfrk.f"
    return 0;

/*     End of ZHFRK */

} /* zhfrk_ */

