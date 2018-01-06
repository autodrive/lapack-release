#line 1 "cunbdb4.f"
/* cunbdb4.f -- translated by f2c (version 20100827).
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

#line 1 "cunbdb4.f"
/* Table of constant values */

static doublecomplex c_b1 = {-1.,0.};
static integer c__1 = 1;

/* > \brief \b CUNBDB4 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNBDB4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb4
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb4
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb4
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, */
/*                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, */
/*                           INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               PHI(*), THETA(*) */
/*       COMPLEX            PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), */
/*      $                   WORK(*), X11(LDX11,*), X21(LDX21,*) */
/*       .. */


/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > CUNBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny */
/* > matrix X with orthonomal columns: */
/* > */
/* >                            [ B11 ] */
/* >      [ X11 ]   [ P1 |    ] [  0  ] */
/* >      [-----] = [---------] [-----] Q1**T . */
/* >      [ X21 ]   [    | P2 ] [ B21 ] */
/* >                            [  0  ] */
/* > */
/* > X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P, */
/* > M-P, or Q. Routines CUNBDB1, CUNBDB2, and CUNBDB3 handle cases in */
/* > which M-Q is not the minimum dimension. */
/* > */
/* > The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P), */
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
/* >          X11 is COMPLEX array, dimension (LDX11,Q) */
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
/* >          X21 is COMPLEX array, dimension (LDX21,Q) */
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
/* >          TAUP1 is COMPLEX array, dimension (P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is COMPLEX array, dimension (M-P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is COMPLEX array, dimension (Q) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] PHANTOM */
/* > \verbatim */
/* >          PHANTOM is COMPLEX array, dimension (M) */
/* >           The routine computes an M-by-1 column vector Y that is */
/* >           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and */
/* >           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and */
/* >           Y(P+1:M), respectively. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (LWORK) */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > \verbatim */
/* > */
/* >  The upper-bidiagonal blocks B11, B21 are represented implicitly by */
/* >  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry */
/* >  in each bidiagonal band is a product of a sine or cosine of a THETA */
/* >  with a sine or cosine of a PHI. See [1] or CUNCSD for details. */
/* > */
/* >  P1, P2, and Q1 are represented as products of elementary reflectors. */
/* >  See CUNCSD2BY1 for details on generating P1, P2, and Q1 using CUNGQR */
/* >  and CUNGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cunbdb4_(integer *m, integer *p, integer *q, 
	doublecomplex *x11, integer *ldx11, doublecomplex *x21, integer *
	ldx21, doublereal *theta, doublereal *phi, doublecomplex *taup1, 
	doublecomplex *taup2, doublecomplex *tauq1, doublecomplex *phantom, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double atan2(doublereal, doublereal), cos(doublereal), sin(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static integer lworkmin, lworkopt;
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer childinfo;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static integer ilarf, llarf;
    extern /* Subroutine */ int csrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int cunbdb5_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *);
    static integer iorbdb5, lorbdb5;
    extern /* Subroutine */ int clarfgp_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
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

#line 257 "cunbdb4.f"
    /* Parameter adjustments */
#line 257 "cunbdb4.f"
    x11_dim1 = *ldx11;
#line 257 "cunbdb4.f"
    x11_offset = 1 + x11_dim1;
#line 257 "cunbdb4.f"
    x11 -= x11_offset;
#line 257 "cunbdb4.f"
    x21_dim1 = *ldx21;
#line 257 "cunbdb4.f"
    x21_offset = 1 + x21_dim1;
#line 257 "cunbdb4.f"
    x21 -= x21_offset;
#line 257 "cunbdb4.f"
    --theta;
#line 257 "cunbdb4.f"
    --phi;
#line 257 "cunbdb4.f"
    --taup1;
#line 257 "cunbdb4.f"
    --taup2;
#line 257 "cunbdb4.f"
    --tauq1;
#line 257 "cunbdb4.f"
    --phantom;
#line 257 "cunbdb4.f"
    --work;
#line 257 "cunbdb4.f"

#line 257 "cunbdb4.f"
    /* Function Body */
#line 257 "cunbdb4.f"
    *info = 0;
#line 258 "cunbdb4.f"
    lquery = *lwork == -1;

#line 260 "cunbdb4.f"
    if (*m < 0) {
#line 261 "cunbdb4.f"
	*info = -1;
#line 262 "cunbdb4.f"
    } else if (*p < *m - *q || *m - *p < *m - *q) {
#line 263 "cunbdb4.f"
	*info = -2;
#line 264 "cunbdb4.f"
    } else if (*q < *m - *q || *q > *m) {
#line 265 "cunbdb4.f"
	*info = -3;
#line 266 "cunbdb4.f"
    } else if (*ldx11 < max(1,*p)) {
#line 267 "cunbdb4.f"
	*info = -5;
#line 268 "cunbdb4.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 268 "cunbdb4.f"
	i__1 = 1, i__2 = *m - *p;
#line 268 "cunbdb4.f"
	if (*ldx21 < max(i__1,i__2)) {
#line 269 "cunbdb4.f"
	    *info = -7;
#line 270 "cunbdb4.f"
	}
#line 270 "cunbdb4.f"
    }

/*     Compute workspace */

#line 274 "cunbdb4.f"
    if (*info == 0) {
#line 275 "cunbdb4.f"
	ilarf = 2;
/* Computing MAX */
#line 276 "cunbdb4.f"
	i__1 = *q - 1, i__2 = *p - 1, i__1 = max(i__1,i__2), i__2 = *m - *p - 
		1;
#line 276 "cunbdb4.f"
	llarf = max(i__1,i__2);
#line 277 "cunbdb4.f"
	iorbdb5 = 2;
#line 278 "cunbdb4.f"
	lorbdb5 = *q;
#line 279 "cunbdb4.f"
	lworkopt = ilarf + llarf - 1;
/* Computing MAX */
#line 280 "cunbdb4.f"
	i__1 = lworkopt, i__2 = iorbdb5 + lorbdb5 - 1;
#line 280 "cunbdb4.f"
	lworkopt = max(i__1,i__2);
#line 281 "cunbdb4.f"
	lworkmin = lworkopt;
#line 282 "cunbdb4.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 283 "cunbdb4.f"
	if (*lwork < lworkmin && ! lquery) {
#line 284 "cunbdb4.f"
	    *info = -14;
#line 285 "cunbdb4.f"
	}
#line 286 "cunbdb4.f"
    }
#line 287 "cunbdb4.f"
    if (*info != 0) {
#line 288 "cunbdb4.f"
	i__1 = -(*info);
#line 288 "cunbdb4.f"
	xerbla_("CUNBDB4", &i__1, (ftnlen)7);
#line 289 "cunbdb4.f"
	return 0;
#line 290 "cunbdb4.f"
    } else if (lquery) {
#line 291 "cunbdb4.f"
	return 0;
#line 292 "cunbdb4.f"
    }

/*     Reduce columns 1, ..., M-Q of X11 and X21 */

#line 296 "cunbdb4.f"
    i__1 = *m - *q;
#line 296 "cunbdb4.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 298 "cunbdb4.f"
	if (i__ == 1) {
#line 299 "cunbdb4.f"
	    i__2 = *m;
#line 299 "cunbdb4.f"
	    for (j = 1; j <= i__2; ++j) {
#line 300 "cunbdb4.f"
		i__3 = j;
#line 300 "cunbdb4.f"
		phantom[i__3].r = 0., phantom[i__3].i = 0.;
#line 301 "cunbdb4.f"
	    }
#line 302 "cunbdb4.f"
	    i__2 = *m - *p;
#line 302 "cunbdb4.f"
	    cunbdb5_(p, &i__2, q, &phantom[1], &c__1, &phantom[*p + 1], &c__1,
		     &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &work[
		    iorbdb5], &lorbdb5, &childinfo);
#line 305 "cunbdb4.f"
	    cscal_(p, &c_b1, &phantom[1], &c__1);
#line 306 "cunbdb4.f"
	    clarfgp_(p, &phantom[1], &phantom[2], &c__1, &taup1[1]);
#line 307 "cunbdb4.f"
	    i__2 = *m - *p;
#line 307 "cunbdb4.f"
	    clarfgp_(&i__2, &phantom[*p + 1], &phantom[*p + 2], &c__1, &taup2[
		    1]);
#line 308 "cunbdb4.f"
	    theta[i__] = atan2((doublereal) phantom[1].r, (doublereal) 
		    phantom[*p + 1].r);
#line 309 "cunbdb4.f"
	    c__ = cos(theta[i__]);
#line 310 "cunbdb4.f"
	    s = sin(theta[i__]);
#line 311 "cunbdb4.f"
	    phantom[1].r = 1., phantom[1].i = 0.;
#line 312 "cunbdb4.f"
	    i__2 = *p + 1;
#line 312 "cunbdb4.f"
	    phantom[i__2].r = 1., phantom[i__2].i = 0.;
#line 313 "cunbdb4.f"
	    d_cnjg(&z__1, &taup1[1]);
#line 313 "cunbdb4.f"
	    clarf_("L", p, q, &phantom[1], &c__1, &z__1, &x11[x11_offset], 
		    ldx11, &work[ilarf], (ftnlen)1);
#line 315 "cunbdb4.f"
	    i__2 = *m - *p;
#line 315 "cunbdb4.f"
	    d_cnjg(&z__1, &taup2[1]);
#line 315 "cunbdb4.f"
	    clarf_("L", &i__2, q, &phantom[*p + 1], &c__1, &z__1, &x21[
		    x21_offset], ldx21, &work[ilarf], (ftnlen)1);
#line 317 "cunbdb4.f"
	} else {
#line 318 "cunbdb4.f"
	    i__2 = *p - i__ + 1;
#line 318 "cunbdb4.f"
	    i__3 = *m - *p - i__ + 1;
#line 318 "cunbdb4.f"
	    i__4 = *q - i__ + 1;
#line 318 "cunbdb4.f"
	    cunbdb5_(&i__2, &i__3, &i__4, &x11[i__ + (i__ - 1) * x11_dim1], &
		    c__1, &x21[i__ + (i__ - 1) * x21_dim1], &c__1, &x11[i__ + 
		    i__ * x11_dim1], ldx11, &x21[i__ + i__ * x21_dim1], ldx21,
		     &work[iorbdb5], &lorbdb5, &childinfo);
#line 321 "cunbdb4.f"
	    i__2 = *p - i__ + 1;
#line 321 "cunbdb4.f"
	    cscal_(&i__2, &c_b1, &x11[i__ + (i__ - 1) * x11_dim1], &c__1);
#line 322 "cunbdb4.f"
	    i__2 = *p - i__ + 1;
#line 322 "cunbdb4.f"
	    clarfgp_(&i__2, &x11[i__ + (i__ - 1) * x11_dim1], &x11[i__ + 1 + (
		    i__ - 1) * x11_dim1], &c__1, &taup1[i__]);
#line 323 "cunbdb4.f"
	    i__2 = *m - *p - i__ + 1;
#line 323 "cunbdb4.f"
	    clarfgp_(&i__2, &x21[i__ + (i__ - 1) * x21_dim1], &x21[i__ + 1 + (
		    i__ - 1) * x21_dim1], &c__1, &taup2[i__]);
#line 325 "cunbdb4.f"
	    theta[i__] = atan2((doublereal) x11[i__ + (i__ - 1) * x11_dim1].r,
		     (doublereal) x21[i__ + (i__ - 1) * x21_dim1].r);
#line 326 "cunbdb4.f"
	    c__ = cos(theta[i__]);
#line 327 "cunbdb4.f"
	    s = sin(theta[i__]);
#line 328 "cunbdb4.f"
	    i__2 = i__ + (i__ - 1) * x11_dim1;
#line 328 "cunbdb4.f"
	    x11[i__2].r = 1., x11[i__2].i = 0.;
#line 329 "cunbdb4.f"
	    i__2 = i__ + (i__ - 1) * x21_dim1;
#line 329 "cunbdb4.f"
	    x21[i__2].r = 1., x21[i__2].i = 0.;
#line 330 "cunbdb4.f"
	    i__2 = *p - i__ + 1;
#line 330 "cunbdb4.f"
	    i__3 = *q - i__ + 1;
#line 330 "cunbdb4.f"
	    d_cnjg(&z__1, &taup1[i__]);
#line 330 "cunbdb4.f"
	    clarf_("L", &i__2, &i__3, &x11[i__ + (i__ - 1) * x11_dim1], &c__1,
		     &z__1, &x11[i__ + i__ * x11_dim1], ldx11, &work[ilarf], (
		    ftnlen)1);
#line 332 "cunbdb4.f"
	    i__2 = *m - *p - i__ + 1;
#line 332 "cunbdb4.f"
	    i__3 = *q - i__ + 1;
#line 332 "cunbdb4.f"
	    d_cnjg(&z__1, &taup2[i__]);
#line 332 "cunbdb4.f"
	    clarf_("L", &i__2, &i__3, &x21[i__ + (i__ - 1) * x21_dim1], &c__1,
		     &z__1, &x21[i__ + i__ * x21_dim1], ldx21, &work[ilarf], (
		    ftnlen)1);
#line 334 "cunbdb4.f"
	}

#line 336 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 336 "cunbdb4.f"
	d__1 = -c__;
#line 336 "cunbdb4.f"
	csrot_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11, &x21[i__ + i__ * 
		x21_dim1], ldx21, &s, &d__1);
#line 337 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 337 "cunbdb4.f"
	clacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);
#line 338 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 338 "cunbdb4.f"
	clarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 1) * 
		x21_dim1], ldx21, &tauq1[i__]);
#line 339 "cunbdb4.f"
	i__2 = i__ + i__ * x21_dim1;
#line 339 "cunbdb4.f"
	c__ = x21[i__2].r;
#line 340 "cunbdb4.f"
	i__2 = i__ + i__ * x21_dim1;
#line 340 "cunbdb4.f"
	x21[i__2].r = 1., x21[i__2].i = 0.;
#line 341 "cunbdb4.f"
	i__2 = *p - i__;
#line 341 "cunbdb4.f"
	i__3 = *q - i__ + 1;
#line 341 "cunbdb4.f"
	clarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf], (
		ftnlen)1);
#line 343 "cunbdb4.f"
	i__2 = *m - *p - i__;
#line 343 "cunbdb4.f"
	i__3 = *q - i__ + 1;
#line 343 "cunbdb4.f"
	clarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[ilarf], (
		ftnlen)1);
#line 345 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 345 "cunbdb4.f"
	clacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);
#line 346 "cunbdb4.f"
	if (i__ < *m - *q) {
#line 347 "cunbdb4.f"
	    i__2 = *p - i__;
/* Computing 2nd power */
#line 347 "cunbdb4.f"
	    d__1 = scnrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 347 "cunbdb4.f"
	    i__3 = *m - *p - i__;
/* Computing 2nd power */
#line 347 "cunbdb4.f"
	    d__2 = scnrm2_(&i__3, &x21[i__ + 1 + i__ * x21_dim1], &c__1);
#line 347 "cunbdb4.f"
	    s = sqrt(d__1 * d__1 + d__2 * d__2);
#line 349 "cunbdb4.f"
	    phi[i__] = atan2(s, c__);
#line 350 "cunbdb4.f"
	}

#line 352 "cunbdb4.f"
    }

/*     Reduce the bottom-right portion of X11 to [ I 0 ] */

#line 356 "cunbdb4.f"
    i__1 = *p;
#line 356 "cunbdb4.f"
    for (i__ = *m - *q + 1; i__ <= i__1; ++i__) {
#line 357 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 357 "cunbdb4.f"
	clacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 358 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 358 "cunbdb4.f"
	clarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) * 
		x11_dim1], ldx11, &tauq1[i__]);
#line 359 "cunbdb4.f"
	i__2 = i__ + i__ * x11_dim1;
#line 359 "cunbdb4.f"
	x11[i__2].r = 1., x11[i__2].i = 0.;
#line 360 "cunbdb4.f"
	i__2 = *p - i__;
#line 360 "cunbdb4.f"
	i__3 = *q - i__ + 1;
#line 360 "cunbdb4.f"
	clarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf], (
		ftnlen)1);
#line 362 "cunbdb4.f"
	i__2 = *q - *p;
#line 362 "cunbdb4.f"
	i__3 = *q - i__ + 1;
#line 362 "cunbdb4.f"
	clarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x21[*m - *q + 1 + i__ * x21_dim1], ldx21, &work[ilarf],
		 (ftnlen)1);
#line 364 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 364 "cunbdb4.f"
	clacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 365 "cunbdb4.f"
    }

/*     Reduce the bottom-right portion of X21 to [ 0 I ] */

#line 369 "cunbdb4.f"
    i__1 = *q;
#line 369 "cunbdb4.f"
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
#line 370 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 370 "cunbdb4.f"
	clacgv_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], ldx21);
#line 371 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 371 "cunbdb4.f"
	clarfgp_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], &x21[*m - *
		q + i__ - *p + (i__ + 1) * x21_dim1], ldx21, &tauq1[i__]);
#line 373 "cunbdb4.f"
	i__2 = *m - *q + i__ - *p + i__ * x21_dim1;
#line 373 "cunbdb4.f"
	x21[i__2].r = 1., x21[i__2].i = 0.;
#line 374 "cunbdb4.f"
	i__2 = *q - i__;
#line 374 "cunbdb4.f"
	i__3 = *q - i__ + 1;
#line 374 "cunbdb4.f"
	clarf_("R", &i__2, &i__3, &x21[*m - *q + i__ - *p + i__ * x21_dim1], 
		ldx21, &tauq1[i__], &x21[*m - *q + i__ - *p + 1 + i__ * 
		x21_dim1], ldx21, &work[ilarf], (ftnlen)1);
#line 376 "cunbdb4.f"
	i__2 = *q - i__ + 1;
#line 376 "cunbdb4.f"
	clacgv_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], ldx21);
#line 377 "cunbdb4.f"
    }

#line 379 "cunbdb4.f"
    return 0;

/*     End of CUNBDB4 */

} /* cunbdb4_ */

