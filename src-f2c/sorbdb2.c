#line 1 "sorbdb2.f"
/* sorbdb2.f -- translated by f2c (version 20100827).
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

#line 1 "sorbdb2.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = -1.;

/* > \brief \b SORBDB2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORBDB2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb2
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb2
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb2
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, */
/*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               PHI(*), THETA(*) */
/*       REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), */
/*      $                   X11(LDX11,*), X21(LDX21,*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny */
/* > matrix X with orthonomal columns: */
/* > */
/* >                            [ B11 ] */
/* >      [ X11 ]   [ P1 |    ] [  0  ] */
/* >      [-----] = [---------] [-----] Q1**T . */
/* >      [ X21 ]   [    | P2 ] [ B21 ] */
/* >                            [  0  ] */
/* > */
/* > X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P, */
/* > Q, or M-Q. Routines SORBDB1, SORBDB3, and SORBDB4 handle cases in */
/* > which P is not the minimum dimension. */
/* > */
/* > The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P), */
/* > and (M-Q)-by-(M-Q), respectively. They are represented implicitly by */
/* > Householder vectors. */
/* > */
/* > B11 and B12 are P-by-P bidiagonal matrices represented implicitly by */
/* > angles THETA, PHI. */
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
/* >           The number of rows in X11. 0 <= P <= min(M-P,Q,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >           The number of columns in X11 and X21. 0 <= Q <= M. */
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
/* Subroutine */ int sorbdb2_(integer *m, integer *p, integer *q, doublereal *
	x11, integer *ldx11, doublereal *x21, integer *ldx21, doublereal *
	theta, doublereal *phi, doublereal *taup1, doublereal *taup2, 
	doublereal *tauq1, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), cos(doublereal), 
	    sin(doublereal);

    /* Local variables */
    static integer lworkmin, lworkopt;
    static doublereal c__;
    static integer i__;
    static doublereal s;
    static integer childinfo;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *, doublereal *,
	     integer *);
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


/*  -- LAPACK computational routine (version 3.5.0) -- */
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

#line 244 "sorbdb2.f"
    /* Parameter adjustments */
#line 244 "sorbdb2.f"
    x11_dim1 = *ldx11;
#line 244 "sorbdb2.f"
    x11_offset = 1 + x11_dim1;
#line 244 "sorbdb2.f"
    x11 -= x11_offset;
#line 244 "sorbdb2.f"
    x21_dim1 = *ldx21;
#line 244 "sorbdb2.f"
    x21_offset = 1 + x21_dim1;
#line 244 "sorbdb2.f"
    x21 -= x21_offset;
#line 244 "sorbdb2.f"
    --theta;
#line 244 "sorbdb2.f"
    --phi;
#line 244 "sorbdb2.f"
    --taup1;
#line 244 "sorbdb2.f"
    --taup2;
#line 244 "sorbdb2.f"
    --tauq1;
#line 244 "sorbdb2.f"
    --work;
#line 244 "sorbdb2.f"

#line 244 "sorbdb2.f"
    /* Function Body */
#line 244 "sorbdb2.f"
    *info = 0;
#line 245 "sorbdb2.f"
    lquery = *lwork == -1;

#line 247 "sorbdb2.f"
    if (*m < 0) {
#line 248 "sorbdb2.f"
	*info = -1;
#line 249 "sorbdb2.f"
    } else if (*p < 0 || *p > *m - *p) {
#line 250 "sorbdb2.f"
	*info = -2;
#line 251 "sorbdb2.f"
    } else if (*q < 0 || *q < *p || *m - *q < *p) {
#line 252 "sorbdb2.f"
	*info = -3;
#line 253 "sorbdb2.f"
    } else if (*ldx11 < max(1,*p)) {
#line 254 "sorbdb2.f"
	*info = -5;
#line 255 "sorbdb2.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 255 "sorbdb2.f"
	i__1 = 1, i__2 = *m - *p;
#line 255 "sorbdb2.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 256 "sorbdb2.f"
	    *info = -7;
#line 257 "sorbdb2.f"
	}
#line 257 "sorbdb2.f"
    }

/*     Compute workspace */

#line 261 "sorbdb2.f"
    if (*info == 0) {
#line 262 "sorbdb2.f"
	ilarf = 2;
/* Computing MAX */
#line 263 "sorbdb2.f"
	i__1 = *p - 1, i__2 = *m - *p, i__1 = max(i__1,i__2), i__2 = *q - 1;
#line 263 "sorbdb2.f"
	llarf = max(i__1,i__2);
#line 264 "sorbdb2.f"
	iorbdb5 = 2;
#line 265 "sorbdb2.f"
	lorbdb5 = *q - 1;
/* Computing MAX */
#line 266 "sorbdb2.f"
	i__1 = ilarf + llarf - 1, i__2 = iorbdb5 + lorbdb5 - 1;
#line 266 "sorbdb2.f"
	lworkopt = max(i__1,i__2);
#line 267 "sorbdb2.f"
	lworkmin = lworkopt;
#line 268 "sorbdb2.f"
	work[1] = (doublereal) lworkopt;
#line 269 "sorbdb2.f"
	if (*lwork < lworkmin && ! lquery) {
#line 270 "sorbdb2.f"
	    *info = -14;
#line 271 "sorbdb2.f"
	}
#line 272 "sorbdb2.f"
    }
#line 273 "sorbdb2.f"
    if (*info != 0) {
#line 274 "sorbdb2.f"
	i__1 = -(*info);
#line 274 "sorbdb2.f"
	xerbla_("SORBDB2", &i__1, (ftnlen)7);
#line 275 "sorbdb2.f"
	return 0;
#line 276 "sorbdb2.f"
    } else if (lquery) {
#line 277 "sorbdb2.f"
	return 0;
#line 278 "sorbdb2.f"
    }

/*     Reduce rows 1, ..., P of X11 and X21 */

#line 282 "sorbdb2.f"
    i__1 = *p;
#line 282 "sorbdb2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 284 "sorbdb2.f"
	if (i__ > 1) {
#line 285 "sorbdb2.f"
	    i__2 = *q - i__ + 1;
#line 285 "sorbdb2.f"
	    srot_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11, &x21[i__ - 1 + 
		    i__ * x21_dim1], ldx21, &c__, &s);
#line 286 "sorbdb2.f"
	}
#line 287 "sorbdb2.f"
	i__2 = *q - i__ + 1;
#line 287 "sorbdb2.f"
	slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) * 
		x11_dim1], ldx11, &tauq1[i__]);
#line 288 "sorbdb2.f"
	c__ = x11[i__ + i__ * x11_dim1];
#line 289 "sorbdb2.f"
	x11[i__ + i__ * x11_dim1] = 1.;
#line 290 "sorbdb2.f"
	i__2 = *p - i__;
#line 290 "sorbdb2.f"
	i__3 = *q - i__ + 1;
#line 290 "sorbdb2.f"
	slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf], (
		ftnlen)1);
#line 292 "sorbdb2.f"
	i__2 = *m - *p - i__ + 1;
#line 292 "sorbdb2.f"
	i__3 = *q - i__ + 1;
#line 292 "sorbdb2.f"
	slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x21[i__ + i__ * x21_dim1], ldx21, &work[ilarf], (
		ftnlen)1);
#line 294 "sorbdb2.f"
	i__2 = *p - i__;
/* Computing 2nd power */
#line 294 "sorbdb2.f"
	d__1 = snrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &c__1, &x11[i__ 
		+ 1 + i__ * x11_dim1], &c__1);
#line 294 "sorbdb2.f"
	i__3 = *m - *p - i__ + 1;
/* Computing 2nd power */
#line 294 "sorbdb2.f"
	d__2 = snrm2_(&i__3, &x21[i__ + i__ * x21_dim1], &c__1, &x21[i__ + 
		i__ * x21_dim1], &c__1);
#line 294 "sorbdb2.f"
	s = sqrt(d__1 * d__1 + d__2 * d__2);
#line 296 "sorbdb2.f"
	theta[i__] = atan2(s, c__);

#line 298 "sorbdb2.f"
	i__2 = *p - i__;
#line 298 "sorbdb2.f"
	i__3 = *m - *p - i__ + 1;
#line 298 "sorbdb2.f"
	i__4 = *q - i__;
#line 298 "sorbdb2.f"
	sorbdb5_(&i__2, &i__3, &i__4, &x11[i__ + 1 + i__ * x11_dim1], &c__1, &
		x21[i__ + i__ * x21_dim1], &c__1, &x11[i__ + 1 + (i__ + 1) * 
		x11_dim1], ldx11, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &
		work[iorbdb5], &lorbdb5, &childinfo);
#line 301 "sorbdb2.f"
	i__2 = *p - i__;
#line 301 "sorbdb2.f"
	sscal_(&i__2, &c_b11, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 302 "sorbdb2.f"
	i__2 = *m - *p - i__ + 1;
#line 302 "sorbdb2.f"
	slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + i__ * 
		x21_dim1], &c__1, &taup2[i__]);
#line 303 "sorbdb2.f"
	if (i__ < *p) {
#line 304 "sorbdb2.f"
	    i__2 = *p - i__;
#line 304 "sorbdb2.f"
	    slarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ + 2 + 
		    i__ * x11_dim1], &c__1, &taup1[i__]);
#line 305 "sorbdb2.f"
	    phi[i__] = atan2(x11[i__ + 1 + i__ * x11_dim1], x21[i__ + i__ * 
		    x21_dim1]);
#line 306 "sorbdb2.f"
	    c__ = cos(phi[i__]);
#line 307 "sorbdb2.f"
	    s = sin(phi[i__]);
#line 308 "sorbdb2.f"
	    x11[i__ + 1 + i__ * x11_dim1] = 1.;
#line 309 "sorbdb2.f"
	    i__2 = *p - i__;
#line 309 "sorbdb2.f"
	    i__3 = *q - i__;
#line 309 "sorbdb2.f"
	    slarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &c__1, &
		    taup1[i__], &x11[i__ + 1 + (i__ + 1) * x11_dim1], ldx11, &
		    work[ilarf], (ftnlen)1);
#line 311 "sorbdb2.f"
	}
#line 312 "sorbdb2.f"
	x21[i__ + i__ * x21_dim1] = 1.;
#line 313 "sorbdb2.f"
	i__2 = *m - *p - i__ + 1;
#line 313 "sorbdb2.f"
	i__3 = *q - i__;
#line 313 "sorbdb2.f"
	slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &taup2[
		i__], &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[ilarf], (
		ftnlen)1);

#line 316 "sorbdb2.f"
    }

/*     Reduce the bottom-right portion of X21 to the identity matrix */

#line 320 "sorbdb2.f"
    i__1 = *q;
#line 320 "sorbdb2.f"
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
#line 321 "sorbdb2.f"
	i__2 = *m - *p - i__ + 1;
#line 321 "sorbdb2.f"
	slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + i__ * 
		x21_dim1], &c__1, &taup2[i__]);
#line 322 "sorbdb2.f"
	x21[i__ + i__ * x21_dim1] = 1.;
#line 323 "sorbdb2.f"
	i__2 = *m - *p - i__ + 1;
#line 323 "sorbdb2.f"
	i__3 = *q - i__;
#line 323 "sorbdb2.f"
	slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &taup2[
		i__], &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[ilarf], (
		ftnlen)1);
#line 325 "sorbdb2.f"
    }

#line 327 "sorbdb2.f"
    return 0;

/*     End of SORBDB2 */

} /* sorbdb2_ */

