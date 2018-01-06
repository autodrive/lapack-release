#line 1 "slagv2.f"
/* slagv2.f -- translated by f2c (version 20100827).
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

#line 1 "slagv2.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* > \brief \b SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 matrix pencil (A,B) where 
B is upper triangular. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAGV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagv2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagv2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagv2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL, */
/*                          CSR, SNR ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDB */
/*       REAL               CSL, CSR, SNL, SNR */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ), */
/*      $                   B( LDB, * ), BETA( 2 ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 */
/* > matrix pencil (A,B) where B is upper triangular. This routine */
/* > computes orthogonal (rotation) matrices given by CSL, SNL and CSR, */
/* > SNR such that */
/* > */
/* > 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0 */
/* >    types), then */
/* > */
/* >    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/* >    [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */
/* > */
/* >    [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ] */
/* >    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ], */
/* > */
/* > 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues, */
/* >    then */
/* > */
/* >    [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/* >    [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */
/* > */
/* >    [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ] */
/* >    [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ] */
/* > */
/* >    where b11 >= b22 > 0. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA, 2) */
/* >          On entry, the 2 x 2 matrix A. */
/* >          On exit, A is overwritten by the ``A-part'' of the */
/* >          generalized Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          THe leading dimension of the array A.  LDA >= 2. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB, 2) */
/* >          On entry, the upper triangular 2 x 2 matrix B. */
/* >          On exit, B is overwritten by the ``B-part'' of the */
/* >          generalized Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          THe leading dimension of the array B.  LDB >= 2. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* >          ALPHAR is REAL array, dimension (2) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is REAL array, dimension (2) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is REAL array, dimension (2) */
/* >          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the */
/* >          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may */
/* >          be zero. */
/* > \endverbatim */
/* > */
/* > \param[out] CSL */
/* > \verbatim */
/* >          CSL is REAL */
/* >          The cosine of the left rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SNL */
/* > \verbatim */
/* >          SNL is REAL */
/* >          The sine of the left rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] CSR */
/* > \verbatim */
/* >          CSR is REAL */
/* >          The cosine of the right rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SNR */
/* > \verbatim */
/* >          SNR is REAL */
/* >          The sine of the right rotation matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */

/*  ===================================================================== */
/* Subroutine */ int slagv2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal r__, t, h1, h2, h3, wi, qq, rr, wr1, wr2, ulp;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), slag2_(
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal anorm, bnorm, scale1, scale2;
    extern /* Subroutine */ int slasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal slapy2_(doublereal *, doublereal *);
    static doublereal ascale, bscale;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int slartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 197 "slagv2.f"
    /* Parameter adjustments */
#line 197 "slagv2.f"
    a_dim1 = *lda;
#line 197 "slagv2.f"
    a_offset = 1 + a_dim1;
#line 197 "slagv2.f"
    a -= a_offset;
#line 197 "slagv2.f"
    b_dim1 = *ldb;
#line 197 "slagv2.f"
    b_offset = 1 + b_dim1;
#line 197 "slagv2.f"
    b -= b_offset;
#line 197 "slagv2.f"
    --alphar;
#line 197 "slagv2.f"
    --alphai;
#line 197 "slagv2.f"
    --beta;
#line 197 "slagv2.f"

#line 197 "slagv2.f"
    /* Function Body */
#line 197 "slagv2.f"
    safmin = slamch_("S", (ftnlen)1);
#line 198 "slagv2.f"
    ulp = slamch_("P", (ftnlen)1);

/*     Scale A */

/* Computing MAX */
#line 202 "slagv2.f"
    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[a_dim1 + 2], abs(
	    d__2)), d__6 = (d__3 = a[(a_dim1 << 1) + 1], abs(d__3)) + (d__4 = 
	    a[(a_dim1 << 1) + 2], abs(d__4)), d__5 = max(d__5,d__6);
#line 202 "slagv2.f"
    anorm = max(d__5,safmin);
#line 204 "slagv2.f"
    ascale = 1. / anorm;
#line 205 "slagv2.f"
    a[a_dim1 + 1] = ascale * a[a_dim1 + 1];
#line 206 "slagv2.f"
    a[(a_dim1 << 1) + 1] = ascale * a[(a_dim1 << 1) + 1];
#line 207 "slagv2.f"
    a[a_dim1 + 2] = ascale * a[a_dim1 + 2];
#line 208 "slagv2.f"
    a[(a_dim1 << 1) + 2] = ascale * a[(a_dim1 << 1) + 2];

/*     Scale B */

/* Computing MAX */
#line 212 "slagv2.f"
    d__4 = (d__3 = b[b_dim1 + 1], abs(d__3)), d__5 = (d__1 = b[(b_dim1 << 1) 
	    + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1) + 2], abs(d__2)), d__4 
	    = max(d__4,d__5);
#line 212 "slagv2.f"
    bnorm = max(d__4,safmin);
#line 214 "slagv2.f"
    bscale = 1. / bnorm;
#line 215 "slagv2.f"
    b[b_dim1 + 1] = bscale * b[b_dim1 + 1];
#line 216 "slagv2.f"
    b[(b_dim1 << 1) + 1] = bscale * b[(b_dim1 << 1) + 1];
#line 217 "slagv2.f"
    b[(b_dim1 << 1) + 2] = bscale * b[(b_dim1 << 1) + 2];

/*     Check if A can be deflated */

#line 221 "slagv2.f"
    if ((d__1 = a[a_dim1 + 2], abs(d__1)) <= ulp) {
#line 222 "slagv2.f"
	*csl = 1.;
#line 223 "slagv2.f"
	*snl = 0.;
#line 224 "slagv2.f"
	*csr = 1.;
#line 225 "slagv2.f"
	*snr = 0.;
#line 226 "slagv2.f"
	a[a_dim1 + 2] = 0.;
#line 227 "slagv2.f"
	b[b_dim1 + 2] = 0.;
#line 228 "slagv2.f"
	wi = 0.;

/*     Check if B is singular */

#line 232 "slagv2.f"
    } else if ((d__1 = b[b_dim1 + 1], abs(d__1)) <= ulp) {
#line 233 "slagv2.f"
	slartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);
#line 234 "slagv2.f"
	*csr = 1.;
#line 235 "slagv2.f"
	*snr = 0.;
#line 236 "slagv2.f"
	srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 237 "slagv2.f"
	srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);
#line 238 "slagv2.f"
	a[a_dim1 + 2] = 0.;
#line 239 "slagv2.f"
	b[b_dim1 + 1] = 0.;
#line 240 "slagv2.f"
	b[b_dim1 + 2] = 0.;
#line 241 "slagv2.f"
	wi = 0.;

#line 243 "slagv2.f"
    } else if ((d__1 = b[(b_dim1 << 1) + 2], abs(d__1)) <= ulp) {
#line 244 "slagv2.f"
	slartg_(&a[(a_dim1 << 1) + 2], &a[a_dim1 + 2], csr, snr, &t);
#line 245 "slagv2.f"
	*snr = -(*snr);
#line 246 "slagv2.f"
	srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr,
		 snr);
#line 247 "slagv2.f"
	srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csr,
		 snr);
#line 248 "slagv2.f"
	*csl = 1.;
#line 249 "slagv2.f"
	*snl = 0.;
#line 250 "slagv2.f"
	a[a_dim1 + 2] = 0.;
#line 251 "slagv2.f"
	b[b_dim1 + 2] = 0.;
#line 252 "slagv2.f"
	b[(b_dim1 << 1) + 2] = 0.;
#line 253 "slagv2.f"
	wi = 0.;

#line 255 "slagv2.f"
    } else {

/*        B is nonsingular, first compute the eigenvalues of (A,B) */

#line 259 "slagv2.f"
	slag2_(&a[a_offset], lda, &b[b_offset], ldb, &safmin, &scale1, &
		scale2, &wr1, &wr2, &wi);

#line 262 "slagv2.f"
	if (wi == 0.) {

/*           two real eigenvalues, compute s*A-w*B */

#line 266 "slagv2.f"
	    h1 = scale1 * a[a_dim1 + 1] - wr1 * b[b_dim1 + 1];
#line 267 "slagv2.f"
	    h2 = scale1 * a[(a_dim1 << 1) + 1] - wr1 * b[(b_dim1 << 1) + 1];
#line 268 "slagv2.f"
	    h3 = scale1 * a[(a_dim1 << 1) + 2] - wr1 * b[(b_dim1 << 1) + 2];

#line 270 "slagv2.f"
	    rr = slapy2_(&h1, &h2);
#line 271 "slagv2.f"
	    d__1 = scale1 * a[a_dim1 + 2];
#line 271 "slagv2.f"
	    qq = slapy2_(&d__1, &h3);

#line 273 "slagv2.f"
	    if (rr > qq) {

/*              find right rotation matrix to zero 1,1 element of */
/*              (sA - wB) */

#line 278 "slagv2.f"
		slartg_(&h2, &h1, csr, snr, &t);

#line 280 "slagv2.f"
	    } else {

/*              find right rotation matrix to zero 2,1 element of */
/*              (sA - wB) */

#line 285 "slagv2.f"
		d__1 = scale1 * a[a_dim1 + 2];
#line 285 "slagv2.f"
		slartg_(&h3, &d__1, csr, snr, &t);

#line 287 "slagv2.f"
	    }

#line 289 "slagv2.f"
	    *snr = -(*snr);
#line 290 "slagv2.f"
	    srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
#line 291 "slagv2.f"
	    srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csr, snr);

/*           compute inf norms of A and B */

/* Computing MAX */
#line 295 "slagv2.f"
	    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[(a_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = a[a_dim1 + 2], abs(d__3)
		    ) + (d__4 = a[(a_dim1 << 1) + 2], abs(d__4));
#line 295 "slagv2.f"
	    h1 = max(d__5,d__6);
/* Computing MAX */
#line 297 "slagv2.f"
	    d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = b[b_dim1 + 2], abs(d__3)
		    ) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
#line 297 "slagv2.f"
	    h2 = max(d__5,d__6);

#line 300 "slagv2.f"
	    if (scale1 * h1 >= abs(wr1) * h2) {

/*              find left rotation matrix Q to zero out B(2,1) */

#line 304 "slagv2.f"
		slartg_(&b[b_dim1 + 1], &b[b_dim1 + 2], csl, snl, &r__);

#line 306 "slagv2.f"
	    } else {

/*              find left rotation matrix Q to zero out A(2,1) */

#line 310 "slagv2.f"
		slartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);

#line 312 "slagv2.f"
	    }

#line 314 "slagv2.f"
	    srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 315 "slagv2.f"
	    srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);

#line 317 "slagv2.f"
	    a[a_dim1 + 2] = 0.;
#line 318 "slagv2.f"
	    b[b_dim1 + 2] = 0.;

#line 320 "slagv2.f"
	} else {

/*           a pair of complex conjugate eigenvalues */
/*           first compute the SVD of the matrix B */

#line 325 "slagv2.f"
	    slasv2_(&b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], &b[(b_dim1 << 1) + 
		    2], &r__, &t, snr, csr, snl, csl);

/*           Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and */
/*           Z is right rotation matrix computed from SLASV2 */

#line 331 "slagv2.f"
	    srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 332 "slagv2.f"
	    srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);
#line 333 "slagv2.f"
	    srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
#line 334 "slagv2.f"
	    srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csr, snr);

#line 336 "slagv2.f"
	    b[b_dim1 + 2] = 0.;
#line 337 "slagv2.f"
	    b[(b_dim1 << 1) + 1] = 0.;

#line 339 "slagv2.f"
	}

#line 341 "slagv2.f"
    }

/*     Unscaling */

#line 345 "slagv2.f"
    a[a_dim1 + 1] = anorm * a[a_dim1 + 1];
#line 346 "slagv2.f"
    a[a_dim1 + 2] = anorm * a[a_dim1 + 2];
#line 347 "slagv2.f"
    a[(a_dim1 << 1) + 1] = anorm * a[(a_dim1 << 1) + 1];
#line 348 "slagv2.f"
    a[(a_dim1 << 1) + 2] = anorm * a[(a_dim1 << 1) + 2];
#line 349 "slagv2.f"
    b[b_dim1 + 1] = bnorm * b[b_dim1 + 1];
#line 350 "slagv2.f"
    b[b_dim1 + 2] = bnorm * b[b_dim1 + 2];
#line 351 "slagv2.f"
    b[(b_dim1 << 1) + 1] = bnorm * b[(b_dim1 << 1) + 1];
#line 352 "slagv2.f"
    b[(b_dim1 << 1) + 2] = bnorm * b[(b_dim1 << 1) + 2];

#line 354 "slagv2.f"
    if (wi == 0.) {
#line 355 "slagv2.f"
	alphar[1] = a[a_dim1 + 1];
#line 356 "slagv2.f"
	alphar[2] = a[(a_dim1 << 1) + 2];
#line 357 "slagv2.f"
	alphai[1] = 0.;
#line 358 "slagv2.f"
	alphai[2] = 0.;
#line 359 "slagv2.f"
	beta[1] = b[b_dim1 + 1];
#line 360 "slagv2.f"
	beta[2] = b[(b_dim1 << 1) + 2];
#line 361 "slagv2.f"
    } else {
#line 362 "slagv2.f"
	alphar[1] = anorm * wr1 / scale1 / bnorm;
#line 363 "slagv2.f"
	alphai[1] = anorm * wi / scale1 / bnorm;
#line 364 "slagv2.f"
	alphar[2] = alphar[1];
#line 365 "slagv2.f"
	alphai[2] = -alphai[1];
#line 366 "slagv2.f"
	beta[1] = 1.;
#line 367 "slagv2.f"
	beta[2] = 1.;
#line 368 "slagv2.f"
    }

#line 370 "slagv2.f"
    return 0;

/*     End of SLAGV2 */

} /* slagv2_ */

