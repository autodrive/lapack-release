#line 1 "sorbdb4.f"
/* sorbdb4.f -- translated by f2c (version 20100827).
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

#line 1 "sorbdb4.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = -1.;

/* > \brief \b SORBDB4 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORBDB4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb4
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb4
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb4
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, */
/*                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, */
/*                           INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               PHI(*), THETA(*) */
/*       REAL               PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), */
/*      $                   WORK(*), X11(LDX11,*), X21(LDX21,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny */
/* > matrix X with orthonomal columns: */
/* > */
/* >                            [ B11 ] */
/* >      [ X11 ]   [ P1 |    ] [  0  ] */
/* >      [-----] = [---------] [-----] Q1**T . */
/* >      [ X21 ]   [    | P2 ] [ B21 ] */
/* >                            [  0  ] */
/* > */
/* > X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P, */
/* > M-P, or Q. Routines SORBDB1, SORBDB2, and SORBDB3 handle cases in */
/* > which M-Q is not the minimum dimension. */
/* > */
/* > The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P), */
/* > and (M-Q)-by-(M-Q), respectively. They are represented implicitly by */
/* > Householder vectors. */
/* > */
/* > B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented */
/* > implicitly by angles THETA, PHI. */
/* > */
/* >\endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           The number of rows X11 plus the number of rows in X21. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >           The number of rows in X11. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >           The number of columns in X11 and X21. 0 <= Q <= M and */
/* >           M-Q <= min(P,M-P,Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* >          X11 is REAL array, dimension (LDX11,Q) */
/* >           On entry, the top block of the matrix X to be reduced. On */
/* >           exit, the columns of tril(X11) specify reflectors for P1 and */
/* >           the rows of triu(X11,1) specify reflectors for Q1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* >          LDX11 is INTEGER */
/* >           The leading dimension of X11. LDX11 >= P. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* >          X21 is REAL array, dimension (LDX21,Q) */
/* >           On entry, the bottom block of the matrix X to be reduced. On */
/* >           exit, the columns of tril(X21) specify reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* >          LDX21 is INTEGER */
/* >           The leading dimension of X21. LDX21 >= M-P. */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* >          THETA is REAL array, dimension (Q) */
/* >           The entries of the bidiagonal blocks B11, B21 are defined by */
/* >           THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* >          PHI is REAL array, dimension (Q-1) */
/* >           The entries of the bidiagonal blocks B11, B21 are defined by */
/* >           THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* >          TAUP1 is REAL array, dimension (P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is REAL array, dimension (M-P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is REAL array, dimension (Q) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] PHANTOM */
/* > \verbatim */
/* >          PHANTOM is REAL array, dimension (M) */
/* >           The routine computes an M-by-1 column vector Y that is */
/* >           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and */
/* >           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and */
/* >           Y(P+1:M), respectively. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK. LWORK >= M-Q. */
/* > */
/* >           If LWORK = -1, then a workspace query is assumed; the routine */
/* >           only calculates the optimal size of the WORK array, returns */
/* >           this value as the first entry of the WORK array, and no error */
/* >           message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit. */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* > */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date July 2012 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The upper-bidiagonal blocks B11, B21 are represented implicitly by */
/* >  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry */
/* >  in each bidiagonal band is a product of a sine or cosine of a THETA */
/* >  with a sine or cosine of a PHI. See [1] or SORCSD for details. */
/* > */
/* >  P1, P2, and Q1 are represented as products of elementary reflectors. */
/* >  See SORCSD2BY1 for details on generating P1, P2, and Q1 using SORGQR */
/* >  and SORGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sorbdb4_(integer *m, integer *p, integer *q, doublereal *
	x11, integer *ldx11, doublereal *x21, integer *ldx21, doublereal *
	theta, doublereal *phi, doublereal *taup1, doublereal *taup2, 
	doublereal *tauq1, doublereal *phantom, doublereal *work, integer *
	lwork, integer *info)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan2(doublereal, doublereal), cos(doublereal), sin(doublereal), 
	    sqrt(doublereal);

    /* Local variables */
    static integer lworkmin, lworkopt;
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer childinfo;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer ilarf, llarf;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), slarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical lquery;
    static integer iorbdb5, lorbdb5;
    extern /* Subroutine */ int sorbdb5_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), slarfgp_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     July 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ==================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Function .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */

#line 257 "sorbdb4.f"
    /* Parameter adjustments */
#line 257 "sorbdb4.f"
    x11_dim1 = *ldx11;
#line 257 "sorbdb4.f"
    x11_offset = 1 + x11_dim1;
#line 257 "sorbdb4.f"
    x11 -= x11_offset;
#line 257 "sorbdb4.f"
    x21_dim1 = *ldx21;
#line 257 "sorbdb4.f"
    x21_offset = 1 + x21_dim1;
#line 257 "sorbdb4.f"
    x21 -= x21_offset;
#line 257 "sorbdb4.f"
    --theta;
#line 257 "sorbdb4.f"
    --phi;
#line 257 "sorbdb4.f"
    --taup1;
#line 257 "sorbdb4.f"
    --taup2;
#line 257 "sorbdb4.f"
    --tauq1;
#line 257 "sorbdb4.f"
    --phantom;
#line 257 "sorbdb4.f"
    --work;
#line 257 "sorbdb4.f"

#line 257 "sorbdb4.f"
    /* Function Body */
#line 257 "sorbdb4.f"
    *info = 0;
#line 258 "sorbdb4.f"
    lquery = *lwork == -1;

#line 260 "sorbdb4.f"
    if (*m < 0) {
#line 261 "sorbdb4.f"
	*info = -1;
#line 262 "sorbdb4.f"
    } else if (*p < *m - *q || *m - *p < *m - *q) {
#line 263 "sorbdb4.f"
	*info = -2;
#line 264 "sorbdb4.f"
    } else if (*q < *m - *q || *q > *m) {
#line 265 "sorbdb4.f"
	*info = -3;
#line 266 "sorbdb4.f"
    } else if (*ldx11 < max(1,*p)) {
#line 267 "sorbdb4.f"
	*info = -5;
#line 268 "sorbdb4.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 268 "sorbdb4.f"
	i__1 = 1, i__2 = *m - *p;
#line 268 "sorbdb4.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 269 "sorbdb4.f"
	    *info = -7;
#line 270 "sorbdb4.f"
	}
#line 270 "sorbdb4.f"
    }

/*     Compute workspace */

#line 274 "sorbdb4.f"
    if (*info == 0) {
#line 275 "sorbdb4.f"
	ilarf = 2;
/* Computing MAX */
#line 276 "sorbdb4.f"
	i__1 = *q - 1, i__2 = *p - 1, i__1 = max(i__1,i__2), i__2 = *m - *p - 
		1;
#line 276 "sorbdb4.f"
	llarf = max(i__1,i__2);
#line 277 "sorbdb4.f"
	iorbdb5 = 2;
#line 278 "sorbdb4.f"
	lorbdb5 = *q;
#line 279 "sorbdb4.f"
	lworkopt = ilarf + llarf - 1;
/* Computing MAX */
#line 280 "sorbdb4.f"
	i__1 = lworkopt, i__2 = iorbdb5 + lorbdb5 - 1;
#line 280 "sorbdb4.f"
	lworkopt = max(i__1,i__2);
#line 281 "sorbdb4.f"
	lworkmin = lworkopt;
#line 282 "sorbdb4.f"
	work[1] = (doublereal) lworkopt;
#line 283 "sorbdb4.f"
	if (*lwork < lworkmin && ! lquery) {
#line 284 "sorbdb4.f"
	    *info = -14;
#line 285 "sorbdb4.f"
	}
#line 286 "sorbdb4.f"
    }
#line 287 "sorbdb4.f"
    if (*info != 0) {
#line 288 "sorbdb4.f"
	i__1 = -(*info);
#line 288 "sorbdb4.f"
	xerbla_("SORBDB4", &i__1, (ftnlen)7);
#line 289 "sorbdb4.f"
	return 0;
#line 290 "sorbdb4.f"
    } else if (lquery) {
#line 291 "sorbdb4.f"
	return 0;
#line 292 "sorbdb4.f"
    }

/*     Reduce columns 1, ..., M-Q of X11 and X21 */

#line 296 "sorbdb4.f"
    i__1 = *m - *q;
#line 296 "sorbdb4.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 298 "sorbdb4.f"
	if (i__ == 1) {
#line 299 "sorbdb4.f"
	    i__2 = *m;
#line 299 "sorbdb4.f"
	    for (j = 1; j <= i__2; ++j) {
#line 300 "sorbdb4.f"
		phantom[j] = 0.;
#line 301 "sorbdb4.f"
	    }
#line 302 "sorbdb4.f"
	    i__2 = *m - *p;
#line 302 "sorbdb4.f"
	    sorbdb5_(p, &i__2, q, &phantom[1], &c__1, &phantom[*p + 1], &c__1,
		     &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &work[
		    iorbdb5], &lorbdb5, &childinfo);
#line 305 "sorbdb4.f"
	    sscal_(p, &c_b5, &phantom[1], &c__1);
#line 306 "sorbdb4.f"
	    slarfgp_(p, &phantom[1], &phantom[2], &c__1, &taup1[1]);
#line 307 "sorbdb4.f"
	    i__2 = *m - *p;
#line 307 "sorbdb4.f"
	    slarfgp_(&i__2, &phantom[*p + 1], &phantom[*p + 2], &c__1, &taup2[
		    1]);
#line 308 "sorbdb4.f"
	    theta[i__] = atan2(phantom[1], phantom[*p + 1]);
#line 309 "sorbdb4.f"
	    c__ = cos(theta[i__]);
#line 310 "sorbdb4.f"
	    s = sin(theta[i__]);
#line 311 "sorbdb4.f"
	    phantom[1] = 1.;
#line 312 "sorbdb4.f"
	    phantom[*p + 1] = 1.;
#line 313 "sorbdb4.f"
	    slarf_("L", p, q, &phantom[1], &c__1, &taup1[1], &x11[x11_offset],
		     ldx11, &work[ilarf], (ftnlen)1);
#line 315 "sorbdb4.f"
	    i__2 = *m - *p;
#line 315 "sorbdb4.f"
	    slarf_("L", &i__2, q, &phantom[*p + 1], &c__1, &taup2[1], &x21[
		    x21_offset], ldx21, &work[ilarf], (ftnlen)1);
#line 317 "sorbdb4.f"
	} else {
#line 318 "sorbdb4.f"
	    i__2 = *p - i__ + 1;
#line 318 "sorbdb4.f"
	    i__3 = *m - *p - i__ + 1;
#line 318 "sorbdb4.f"
	    i__4 = *q - i__ + 1;
#line 318 "sorbdb4.f"
	    sorbdb5_(&i__2, &i__3, &i__4, &x11[i__ + (i__ - 1) * x11_dim1], &
		    c__1, &x21[i__ + (i__ - 1) * x21_dim1], &c__1, &x11[i__ + 
		    i__ * x11_dim1], ldx11, &x21[i__ + i__ * x21_dim1], ldx21,
		     &work[iorbdb5], &lorbdb5, &childinfo);
#line 321 "sorbdb4.f"
	    i__2 = *p - i__ + 1;
#line 321 "sorbdb4.f"
	    sscal_(&i__2, &c_b5, &x11[i__ + (i__ - 1) * x11_dim1], &c__1);
#line 322 "sorbdb4.f"
	    i__2 = *p - i__ + 1;
#line 322 "sorbdb4.f"
	    slarfgp_(&i__2, &x11[i__ + (i__ - 1) * x11_dim1], &x11[i__ + 1 + (
		    i__ - 1) * x11_dim1], &c__1, &taup1[i__]);
#line 323 "sorbdb4.f"
	    i__2 = *m - *p - i__ + 1;
#line 323 "sorbdb4.f"
	    slarfgp_(&i__2, &x21[i__ + (i__ - 1) * x21_dim1], &x21[i__ + 1 + (
		    i__ - 1) * x21_dim1], &c__1, &taup2[i__]);
#line 325 "sorbdb4.f"
	    theta[i__] = atan2(x11[i__ + (i__ - 1) * x11_dim1], x21[i__ + (
		    i__ - 1) * x21_dim1]);
#line 326 "sorbdb4.f"
	    c__ = cos(theta[i__]);
#line 327 "sorbdb4.f"
	    s = sin(theta[i__]);
#line 328 "sorbdb4.f"
	    x11[i__ + (i__ - 1) * x11_dim1] = 1.;
#line 329 "sorbdb4.f"
	    x21[i__ + (i__ - 1) * x21_dim1] = 1.;
#line 330 "sorbdb4.f"
	    i__2 = *p - i__ + 1;
#line 330 "sorbdb4.f"
	    i__3 = *q - i__ + 1;
#line 330 "sorbdb4.f"
	    slarf_("L", &i__2, &i__3, &x11[i__ + (i__ - 1) * x11_dim1], &c__1,
		     &taup1[i__], &x11[i__ + i__ * x11_dim1], ldx11, &work[
		    ilarf], (ftnlen)1);
#line 332 "sorbdb4.f"
	    i__2 = *m - *p - i__ + 1;
#line 332 "sorbdb4.f"
	    i__3 = *q - i__ + 1;
#line 332 "sorbdb4.f"
	    slarf_("L", &i__2, &i__3, &x21[i__ + (i__ - 1) * x21_dim1], &c__1,
		     &taup2[i__], &x21[i__ + i__ * x21_dim1], ldx21, &work[
		    ilarf], (ftnlen)1);
#line 334 "sorbdb4.f"
	}

#line 336 "sorbdb4.f"
	i__2 = *q - i__ + 1;
#line 336 "sorbdb4.f"
	d__1 = -c__;
#line 336 "sorbdb4.f"
	srot_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11, &x21[i__ + i__ * 
		x21_dim1], ldx21, &s, &d__1);
#line 337 "sorbdb4.f"
	i__2 = *q - i__ + 1;
#line 337 "sorbdb4.f"
	slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 1) * 
		x21_dim1], ldx21, &tauq1[i__]);
#line 338 "sorbdb4.f"
	c__ = x21[i__ + i__ * x21_dim1];
#line 339 "sorbdb4.f"
	x21[i__ + i__ * x21_dim1] = 1.;
#line 340 "sorbdb4.f"
	i__2 = *p - i__;
#line 340 "sorbdb4.f"
	i__3 = *q - i__ + 1;
#line 340 "sorbdb4.f"
	slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf], (
		ftnlen)1);
#line 342 "sorbdb4.f"
	i__2 = *m - *p - i__;
#line 342 "sorbdb4.f"
	i__3 = *q - i__ + 1;
#line 342 "sorbdb4.f"
	slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[ilarf], (
		ftnlen)1);
#line 344 "sorbdb4.f"
	if (i__ < *m - *q) {
#line 345 "sorbdb4.f"
	    i__2 = *p - i__;
/* Computing 2nd power */
#line 345 "sorbdb4.f"
	    d__1 = snrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 345 "sorbdb4.f"
	    i__3 = *m - *p - i__;
/* Computing 2nd power */
#line 345 "sorbdb4.f"
	    d__2 = snrm2_(&i__3, &x21[i__ + 1 + i__ * x21_dim1], &c__1);
#line 345 "sorbdb4.f"
	    s = sqrt(d__1 * d__1 + d__2 * d__2);
#line 347 "sorbdb4.f"
	    phi[i__] = atan2(s, c__);
#line 348 "sorbdb4.f"
	}

#line 350 "sorbdb4.f"
    }

/*     Reduce the bottom-right portion of X11 to [ I 0 ] */

#line 354 "sorbdb4.f"
    i__1 = *p;
#line 354 "sorbdb4.f"
    for (i__ = *m - *q + 1; i__ <= i__1; ++i__) {
#line 355 "sorbdb4.f"
	i__2 = *q - i__ + 1;
#line 355 "sorbdb4.f"
	slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) * 
		x11_dim1], ldx11, &tauq1[i__]);
#line 356 "sorbdb4.f"
	x11[i__ + i__ * x11_dim1] = 1.;
#line 357 "sorbdb4.f"
	i__2 = *p - i__;
#line 357 "sorbdb4.f"
	i__3 = *q - i__ + 1;
#line 357 "sorbdb4.f"
	slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf], (
		ftnlen)1);
#line 359 "sorbdb4.f"
	i__2 = *q - *p;
#line 359 "sorbdb4.f"
	i__3 = *q - i__ + 1;
#line 359 "sorbdb4.f"
	slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x21[*m - *q + 1 + i__ * x21_dim1], ldx21, &work[ilarf],
		 (ftnlen)1);
#line 361 "sorbdb4.f"
    }

/*     Reduce the bottom-right portion of X21 to [ 0 I ] */

#line 365 "sorbdb4.f"
    i__1 = *q;
#line 365 "sorbdb4.f"
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
#line 366 "sorbdb4.f"
	i__2 = *q - i__ + 1;
#line 366 "sorbdb4.f"
	slarfgp_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], &x21[*m - *
		q + i__ - *p + (i__ + 1) * x21_dim1], ldx21, &tauq1[i__]);
#line 368 "sorbdb4.f"
	x21[*m - *q + i__ - *p + i__ * x21_dim1] = 1.;
#line 369 "sorbdb4.f"
	i__2 = *q - i__;
#line 369 "sorbdb4.f"
	i__3 = *q - i__ + 1;
#line 369 "sorbdb4.f"
	slarf_("R", &i__2, &i__3, &x21[*m - *q + i__ - *p + i__ * x21_dim1], 
		ldx21, &tauq1[i__], &x21[*m - *q + i__ - *p + 1 + i__ * 
		x21_dim1], ldx21, &work[ilarf], (ftnlen)1);
#line 371 "sorbdb4.f"
    }

#line 373 "sorbdb4.f"
    return 0;

/*     End of SORBDB4 */

} /* sorbdb4_ */

