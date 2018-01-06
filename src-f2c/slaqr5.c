#line 1 "slaqr5.f"
/* slaqr5.f -- translated by f2c (version 20100827).
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

#line 1 "slaqr5.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLAQR5 performs a single small-bulge multi-shift QR sweep. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAQR5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqr5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqr5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqr5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, */
/*                          SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, */
/*                          LDU, NV, WV, LDWV, NH, WH, LDWH ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, */
/*      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               H( LDH, * ), SI( * ), SR( * ), U( LDU, * ), */
/*      $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SLAQR5, called by SLAQR0, performs a */
/* >    single small-bulge multi-shift QR sweep. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is logical scalar */
/* >             WANTT = .true. if the quasi-triangular Schur factor */
/* >             is being computed.  WANTT is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is logical scalar */
/* >             WANTZ = .true. if the orthogonal Schur factor is being */
/* >             computed.  WANTZ is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] KACC22 */
/* > \verbatim */
/* >          KACC22 is integer with value 0, 1, or 2. */
/* >             Specifies the computation mode of far-from-diagonal */
/* >             orthogonal updates. */
/* >        = 0: SLAQR5 does not accumulate reflections and does not */
/* >             use matrix-matrix multiply to update far-from-diagonal */
/* >             matrix entries. */
/* >        = 1: SLAQR5 accumulates reflections and uses matrix-matrix */
/* >             multiply to update the far-from-diagonal matrix entries. */
/* >        = 2: SLAQR5 accumulates reflections, uses matrix-matrix */
/* >             multiply to update the far-from-diagonal matrix entries, */
/* >             and takes advantage of 2-by-2 block structure during */
/* >             matrix multiplies. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is integer scalar */
/* >             N is the order of the Hessenberg matrix H upon which this */
/* >             subroutine operates. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* >          KTOP is integer scalar */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* >          KBOT is integer scalar */
/* >             These are the first and last rows and columns of an */
/* >             isolated diagonal block upon which the QR sweep is to be */
/* >             applied. It is assumed without a check that */
/* >                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0 */
/* >             and */
/* >                       either KBOT = N  or   H(KBOT+1,KBOT) = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NSHFTS */
/* > \verbatim */
/* >          NSHFTS is integer scalar */
/* >             NSHFTS gives the number of simultaneous shifts.  NSHFTS */
/* >             must be positive and even. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SR */
/* > \verbatim */
/* >          SR is REAL array of size (NSHFTS) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SI */
/* > \verbatim */
/* >          SI is REAL array of size (NSHFTS) */
/* >             SR contains the real parts and SI contains the imaginary */
/* >             parts of the NSHFTS shifts of origin that define the */
/* >             multi-shift QR sweep.  On output SR and SI may be */
/* >             reordered. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is REAL array of size (LDH,N) */
/* >             On input H contains a Hessenberg matrix.  On output a */
/* >             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied */
/* >             to the isolated diagonal block in rows and columns KTOP */
/* >             through KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is integer scalar */
/* >             LDH is the leading dimension of H just as declared in the */
/* >             calling procedure.  LDH.GE.MAX(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* >          ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* >          IHIZ is INTEGER */
/* >             Specify the rows of Z to which transformations must be */
/* >             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array of size (LDZ,IHI) */
/* >             If WANTZ = .TRUE., then the QR Sweep orthogonal */
/* >             similarity transformation is accumulated into */
/* >             Z(ILOZ:IHIZ,ILO:IHI) from the right. */
/* >             If WANTZ = .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is integer scalar */
/* >             LDA is the leading dimension of Z just as declared in */
/* >             the calling procedure. LDZ.GE.N. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is REAL array of size (LDV,NSHFTS/2) */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is integer scalar */
/* >             LDV is the leading dimension of V as declared in the */
/* >             calling procedure.  LDV.GE.3. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array of size */
/* >             (LDU,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is integer scalar */
/* >             LDU is the leading dimension of U just as declared in the */
/* >             in the calling subroutine.  LDU.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* >          NH is integer scalar */
/* >             NH is the number of columns in array WH available for */
/* >             workspace. NH.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WH */
/* > \verbatim */
/* >          WH is REAL array of size (LDWH,NH) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWH */
/* > \verbatim */
/* >          LDWH is integer scalar */
/* >             Leading dimension of WH just as declared in the */
/* >             calling procedure.  LDWH.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* >          NV is integer scalar */
/* >             NV is the number of rows in WV agailable for workspace. */
/* >             NV.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* >          WV is REAL array of size */
/* >             (LDWV,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* >          LDWV is integer scalar */
/* >             LDWV is the leading dimension of WV as declared in the */
/* >             in the calling subroutine.  LDWV.GE.NV. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */

/* > \par References: */
/*  ================ */
/* > */
/* >       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/* >       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3 */
/* >       Performance, SIAM Journal of Matrix Analysis, volume 23, pages */
/* >       929--947, 2002. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaqr5_(logical *wantt, logical *wantz, integer *kacc22, 
	integer *n, integer *ktop, integer *kbot, integer *nshfts, doublereal 
	*sr, doublereal *si, doublereal *h__, integer *ldh, integer *iloz, 
	integer *ihiz, doublereal *z__, integer *ldz, doublereal *v, integer *
	ldv, doublereal *u, integer *ldu, integer *nv, doublereal *wv, 
	integer *ldwv, integer *nh, doublereal *wh, integer *ldwh)
{
    /* System generated locals */
    integer h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, 
	    wh_offset, wv_dim1, wv_offset, z_dim1, z_offset, i__1, i__2, i__3,
	     i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static integer i__, j, k, m, i2, j2, i4, j4, k1;
    static doublereal h11, h12, h21, h22;
    static integer m22, ns, nu;
    static doublereal vt[3], scl;
    static integer kdu, kms;
    static doublereal ulp;
    static integer knz, kzs;
    static doublereal tst1, tst2, beta;
    static logical blk22, bmp22;
    static integer mend, jcol, jlen, jbot, mbot;
    static doublereal swap;
    static integer jtop, jrow, mtop;
    static doublereal alpha;
    static logical accum;
    static integer ndcol, incol;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer krcol, nbmps;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), slaqr1_(
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), slabad_(doublereal *, 
	    doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int slarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal safmax;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal refsum;
    static integer mstart;
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ================================================================ */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */

/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== If there are no shifts, then there is nothing to do. ==== */

#line 312 "slaqr5.f"
    /* Parameter adjustments */
#line 312 "slaqr5.f"
    --sr;
#line 312 "slaqr5.f"
    --si;
#line 312 "slaqr5.f"
    h_dim1 = *ldh;
#line 312 "slaqr5.f"
    h_offset = 1 + h_dim1;
#line 312 "slaqr5.f"
    h__ -= h_offset;
#line 312 "slaqr5.f"
    z_dim1 = *ldz;
#line 312 "slaqr5.f"
    z_offset = 1 + z_dim1;
#line 312 "slaqr5.f"
    z__ -= z_offset;
#line 312 "slaqr5.f"
    v_dim1 = *ldv;
#line 312 "slaqr5.f"
    v_offset = 1 + v_dim1;
#line 312 "slaqr5.f"
    v -= v_offset;
#line 312 "slaqr5.f"
    u_dim1 = *ldu;
#line 312 "slaqr5.f"
    u_offset = 1 + u_dim1;
#line 312 "slaqr5.f"
    u -= u_offset;
#line 312 "slaqr5.f"
    wv_dim1 = *ldwv;
#line 312 "slaqr5.f"
    wv_offset = 1 + wv_dim1;
#line 312 "slaqr5.f"
    wv -= wv_offset;
#line 312 "slaqr5.f"
    wh_dim1 = *ldwh;
#line 312 "slaqr5.f"
    wh_offset = 1 + wh_dim1;
#line 312 "slaqr5.f"
    wh -= wh_offset;
#line 312 "slaqr5.f"

#line 312 "slaqr5.f"
    /* Function Body */
#line 312 "slaqr5.f"
    if (*nshfts < 2) {
#line 312 "slaqr5.f"
	return 0;
#line 312 "slaqr5.f"
    }

/*     ==== If the active block is empty or 1-by-1, then there */
/*     .    is nothing to do. ==== */

#line 318 "slaqr5.f"
    if (*ktop >= *kbot) {
#line 318 "slaqr5.f"
	return 0;
#line 318 "slaqr5.f"
    }

/*     ==== Shuffle shifts into pairs of real shifts and pairs */
/*     .    of complex conjugate shifts assuming complex */
/*     .    conjugate shifts are already adjacent to one */
/*     .    another. ==== */

#line 326 "slaqr5.f"
    i__1 = *nshfts - 2;
#line 326 "slaqr5.f"
    for (i__ = 1; i__ <= i__1; i__ += 2) {
#line 327 "slaqr5.f"
	if (si[i__] != -si[i__ + 1]) {

#line 329 "slaqr5.f"
	    swap = sr[i__];
#line 330 "slaqr5.f"
	    sr[i__] = sr[i__ + 1];
#line 331 "slaqr5.f"
	    sr[i__ + 1] = sr[i__ + 2];
#line 332 "slaqr5.f"
	    sr[i__ + 2] = swap;

#line 334 "slaqr5.f"
	    swap = si[i__];
#line 335 "slaqr5.f"
	    si[i__] = si[i__ + 1];
#line 336 "slaqr5.f"
	    si[i__ + 1] = si[i__ + 2];
#line 337 "slaqr5.f"
	    si[i__ + 2] = swap;
#line 338 "slaqr5.f"
	}
#line 339 "slaqr5.f"
/* L10: */
#line 339 "slaqr5.f"
    }

/*     ==== NSHFTS is supposed to be even, but if it is odd, */
/*     .    then simply reduce it by one.  The shuffle above */
/*     .    ensures that the dropped shift is real and that */
/*     .    the remaining shifts are paired. ==== */

#line 346 "slaqr5.f"
    ns = *nshfts - *nshfts % 2;

/*     ==== Machine constants for deflation ==== */

#line 350 "slaqr5.f"
    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
#line 351 "slaqr5.f"
    safmax = 1. / safmin;
#line 352 "slaqr5.f"
    slabad_(&safmin, &safmax);
#line 353 "slaqr5.f"
    ulp = slamch_("PRECISION", (ftnlen)9);
#line 354 "slaqr5.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Use accumulated reflections to update far-from-diagonal */
/*     .    entries ? ==== */

#line 359 "slaqr5.f"
    accum = *kacc22 == 1 || *kacc22 == 2;

/*     ==== If so, exploit the 2-by-2 block structure? ==== */

#line 363 "slaqr5.f"
    blk22 = ns > 2 && *kacc22 == 2;

/*     ==== clear trash ==== */

#line 367 "slaqr5.f"
    if (*ktop + 2 <= *kbot) {
#line 367 "slaqr5.f"
	h__[*ktop + 2 + *ktop * h_dim1] = 0.;
#line 367 "slaqr5.f"
    }

/*     ==== NBMPS = number of 2-shift bulges in the chain ==== */

#line 372 "slaqr5.f"
    nbmps = ns / 2;

/*     ==== KDU = width of slab ==== */

#line 376 "slaqr5.f"
    kdu = nbmps * 6 - 3;

/*     ==== Create and chase chains of NBMPS bulges ==== */

#line 380 "slaqr5.f"
    i__1 = *kbot - 2;
#line 380 "slaqr5.f"
    i__2 = nbmps * 3 - 2;
#line 380 "slaqr5.f"
    for (incol = (1 - nbmps) * 3 + *ktop - 1; i__2 < 0 ? incol >= i__1 : 
	    incol <= i__1; incol += i__2) {
#line 381 "slaqr5.f"
	ndcol = incol + kdu;
#line 382 "slaqr5.f"
	if (accum) {
#line 382 "slaqr5.f"
	    slaset_("ALL", &kdu, &kdu, &c_b7, &c_b8, &u[u_offset], ldu, (
		    ftnlen)3);
#line 382 "slaqr5.f"
	}

/*        ==== Near-the-diagonal bulge chase.  The following loop */
/*        .    performs the near-the-diagonal part of a small bulge */
/*        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal */
/*        .    chunk extends from column INCOL to column NDCOL */
/*        .    (including both column INCOL and column NDCOL). The */
/*        .    following loop chases a 3*NBMPS column long chain of */
/*        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL */
/*        .    may be less than KTOP and and NDCOL may be greater than */
/*        .    KBOT indicating phantom columns from which to chase */
/*        .    bulges before they are actually introduced or to which */
/*        .    to chase bulges beyond column KBOT.)  ==== */

/* Computing MIN */
#line 397 "slaqr5.f"
	i__4 = incol + nbmps * 3 - 3, i__5 = *kbot - 2;
#line 397 "slaqr5.f"
	i__3 = min(i__4,i__5);
#line 397 "slaqr5.f"
	for (krcol = incol; krcol <= i__3; ++krcol) {

/*           ==== Bulges number MTOP to MBOT are active double implicit */
/*           .    shift bulges.  There may or may not also be small */
/*           .    2-by-2 bulge, if there is room.  The inactive bulges */
/*           .    (if any) must wait until the active bulges have moved */
/*           .    down the diagonal to make room.  The phantom matrix */
/*           .    paradigm described above helps keep track.  ==== */

/* Computing MAX */
#line 406 "slaqr5.f"
	    i__4 = 1, i__5 = (*ktop - 1 - krcol + 2) / 3 + 1;
#line 406 "slaqr5.f"
	    mtop = max(i__4,i__5);
/* Computing MIN */
#line 407 "slaqr5.f"
	    i__4 = nbmps, i__5 = (*kbot - krcol) / 3;
#line 407 "slaqr5.f"
	    mbot = min(i__4,i__5);
#line 408 "slaqr5.f"
	    m22 = mbot + 1;
#line 409 "slaqr5.f"
	    bmp22 = mbot < nbmps && krcol + (m22 - 1) * 3 == *kbot - 2;

/*           ==== Generate reflections to chase the chain right */
/*           .    one column.  (The minimum value of K is KTOP-1.) ==== */

#line 415 "slaqr5.f"
	    i__4 = mbot;
#line 415 "slaqr5.f"
	    for (m = mtop; m <= i__4; ++m) {
#line 416 "slaqr5.f"
		k = krcol + (m - 1) * 3;
#line 417 "slaqr5.f"
		if (k == *ktop - 1) {
#line 418 "slaqr5.f"
		    slaqr1_(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &sr[(m 
			    << 1) - 1], &si[(m << 1) - 1], &sr[m * 2], &si[m *
			     2], &v[m * v_dim1 + 1]);
#line 421 "slaqr5.f"
		    alpha = v[m * v_dim1 + 1];
#line 422 "slaqr5.f"
		    slarfg_(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1, &v[m * 
			    v_dim1 + 1]);
#line 423 "slaqr5.f"
		} else {
#line 424 "slaqr5.f"
		    beta = h__[k + 1 + k * h_dim1];
#line 425 "slaqr5.f"
		    v[m * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
#line 426 "slaqr5.f"
		    v[m * v_dim1 + 3] = h__[k + 3 + k * h_dim1];
#line 427 "slaqr5.f"
		    slarfg_(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * 
			    v_dim1 + 1]);

/*                 ==== A Bulge may collapse because of vigilant */
/*                 .    deflation or destructive underflow.  In the */
/*                 .    underflow case, try the two-small-subdiagonals */
/*                 .    trick to try to reinflate the bulge.  ==== */

#line 434 "slaqr5.f"
		    if (h__[k + 3 + k * h_dim1] != 0. || h__[k + 3 + (k + 1) *
			     h_dim1] != 0. || h__[k + 3 + (k + 2) * h_dim1] ==
			     0.) {

/*                    ==== Typical case: not collapsed (yet). ==== */

#line 439 "slaqr5.f"
			h__[k + 1 + k * h_dim1] = beta;
#line 440 "slaqr5.f"
			h__[k + 2 + k * h_dim1] = 0.;
#line 441 "slaqr5.f"
			h__[k + 3 + k * h_dim1] = 0.;
#line 442 "slaqr5.f"
		    } else {

/*                    ==== Atypical case: collapsed.  Attempt to */
/*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K). */
/*                    .    If the fill resulting from the new */
/*                    .    reflector is too large, then abandon it. */
/*                    .    Otherwise, use the new one. ==== */

#line 450 "slaqr5.f"
			slaqr1_(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh, &
				sr[(m << 1) - 1], &si[(m << 1) - 1], &sr[m * 
				2], &si[m * 2], vt);
#line 453 "slaqr5.f"
			alpha = vt[0];
#line 454 "slaqr5.f"
			slarfg_(&c__3, &alpha, &vt[1], &c__1, vt);
#line 455 "slaqr5.f"
			refsum = vt[0] * (h__[k + 1 + k * h_dim1] + vt[1] * 
				h__[k + 2 + k * h_dim1]);

#line 458 "slaqr5.f"
			if ((d__1 = h__[k + 2 + k * h_dim1] - refsum * vt[1], 
				abs(d__1)) + (d__2 = refsum * vt[2], abs(d__2)
				) > ulp * ((d__3 = h__[k + k * h_dim1], abs(
				d__3)) + (d__4 = h__[k + 1 + (k + 1) * h_dim1]
				, abs(d__4)) + (d__5 = h__[k + 2 + (k + 2) * 
				h_dim1], abs(d__5)))) {

/*                       ==== Starting a new bulge here would */
/*                       .    create non-negligible fill.  Use */
/*                       .    the old one with trepidation. ==== */

#line 467 "slaqr5.f"
			    h__[k + 1 + k * h_dim1] = beta;
#line 468 "slaqr5.f"
			    h__[k + 2 + k * h_dim1] = 0.;
#line 469 "slaqr5.f"
			    h__[k + 3 + k * h_dim1] = 0.;
#line 470 "slaqr5.f"
			} else {

/*                       ==== Stating a new bulge here would */
/*                       .    create only negligible fill. */
/*                       .    Replace the old reflector with */
/*                       .    the new one. ==== */

#line 477 "slaqr5.f"
			    h__[k + 1 + k * h_dim1] -= refsum;
#line 478 "slaqr5.f"
			    h__[k + 2 + k * h_dim1] = 0.;
#line 479 "slaqr5.f"
			    h__[k + 3 + k * h_dim1] = 0.;
#line 480 "slaqr5.f"
			    v[m * v_dim1 + 1] = vt[0];
#line 481 "slaqr5.f"
			    v[m * v_dim1 + 2] = vt[1];
#line 482 "slaqr5.f"
			    v[m * v_dim1 + 3] = vt[2];
#line 483 "slaqr5.f"
			}
#line 484 "slaqr5.f"
		    }
#line 485 "slaqr5.f"
		}
#line 486 "slaqr5.f"
/* L20: */
#line 486 "slaqr5.f"
	    }

/*           ==== Generate a 2-by-2 reflection, if needed. ==== */

#line 490 "slaqr5.f"
	    k = krcol + (m22 - 1) * 3;
#line 491 "slaqr5.f"
	    if (bmp22) {
#line 492 "slaqr5.f"
		if (k == *ktop - 1) {
#line 493 "slaqr5.f"
		    slaqr1_(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh, &sr[(
			    m22 << 1) - 1], &si[(m22 << 1) - 1], &sr[m22 * 2],
			     &si[m22 * 2], &v[m22 * v_dim1 + 1]);
#line 496 "slaqr5.f"
		    beta = v[m22 * v_dim1 + 1];
#line 497 "slaqr5.f"
		    slarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 
			    * v_dim1 + 1]);
#line 498 "slaqr5.f"
		} else {
#line 499 "slaqr5.f"
		    beta = h__[k + 1 + k * h_dim1];
#line 500 "slaqr5.f"
		    v[m22 * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
#line 501 "slaqr5.f"
		    slarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 
			    * v_dim1 + 1]);
#line 502 "slaqr5.f"
		    h__[k + 1 + k * h_dim1] = beta;
#line 503 "slaqr5.f"
		    h__[k + 2 + k * h_dim1] = 0.;
#line 504 "slaqr5.f"
		}
#line 505 "slaqr5.f"
	    }

/*           ==== Multiply H by reflections from the left ==== */

#line 509 "slaqr5.f"
	    if (accum) {
#line 510 "slaqr5.f"
		jbot = min(ndcol,*kbot);
#line 511 "slaqr5.f"
	    } else if (*wantt) {
#line 512 "slaqr5.f"
		jbot = *n;
#line 513 "slaqr5.f"
	    } else {
#line 514 "slaqr5.f"
		jbot = *kbot;
#line 515 "slaqr5.f"
	    }
#line 516 "slaqr5.f"
	    i__4 = jbot;
#line 516 "slaqr5.f"
	    for (j = max(*ktop,krcol); j <= i__4; ++j) {
/* Computing MIN */
#line 517 "slaqr5.f"
		i__5 = mbot, i__6 = (j - krcol + 2) / 3;
#line 517 "slaqr5.f"
		mend = min(i__5,i__6);
#line 518 "slaqr5.f"
		i__5 = mend;
#line 518 "slaqr5.f"
		for (m = mtop; m <= i__5; ++m) {
#line 519 "slaqr5.f"
		    k = krcol + (m - 1) * 3;
#line 520 "slaqr5.f"
		    refsum = v[m * v_dim1 + 1] * (h__[k + 1 + j * h_dim1] + v[
			    m * v_dim1 + 2] * h__[k + 2 + j * h_dim1] + v[m * 
			    v_dim1 + 3] * h__[k + 3 + j * h_dim1]);
#line 522 "slaqr5.f"
		    h__[k + 1 + j * h_dim1] -= refsum;
#line 523 "slaqr5.f"
		    h__[k + 2 + j * h_dim1] -= refsum * v[m * v_dim1 + 2];
#line 524 "slaqr5.f"
		    h__[k + 3 + j * h_dim1] -= refsum * v[m * v_dim1 + 3];
#line 525 "slaqr5.f"
/* L30: */
#line 525 "slaqr5.f"
		}
#line 526 "slaqr5.f"
/* L40: */
#line 526 "slaqr5.f"
	    }
#line 527 "slaqr5.f"
	    if (bmp22) {
#line 528 "slaqr5.f"
		k = krcol + (m22 - 1) * 3;
/* Computing MAX */
#line 529 "slaqr5.f"
		i__4 = k + 1;
#line 529 "slaqr5.f"
		i__5 = jbot;
#line 529 "slaqr5.f"
		for (j = max(i__4,*ktop); j <= i__5; ++j) {
#line 530 "slaqr5.f"
		    refsum = v[m22 * v_dim1 + 1] * (h__[k + 1 + j * h_dim1] + 
			    v[m22 * v_dim1 + 2] * h__[k + 2 + j * h_dim1]);
#line 532 "slaqr5.f"
		    h__[k + 1 + j * h_dim1] -= refsum;
#line 533 "slaqr5.f"
		    h__[k + 2 + j * h_dim1] -= refsum * v[m22 * v_dim1 + 2];
#line 534 "slaqr5.f"
/* L50: */
#line 534 "slaqr5.f"
		}
#line 535 "slaqr5.f"
	    }

/*           ==== Multiply H by reflections from the right. */
/*           .    Delay filling in the last row until the */
/*           .    vigilant deflation check is complete. ==== */

#line 541 "slaqr5.f"
	    if (accum) {
#line 542 "slaqr5.f"
		jtop = max(*ktop,incol);
#line 543 "slaqr5.f"
	    } else if (*wantt) {
#line 544 "slaqr5.f"
		jtop = 1;
#line 545 "slaqr5.f"
	    } else {
#line 546 "slaqr5.f"
		jtop = *ktop;
#line 547 "slaqr5.f"
	    }
#line 548 "slaqr5.f"
	    i__5 = mbot;
#line 548 "slaqr5.f"
	    for (m = mtop; m <= i__5; ++m) {
#line 549 "slaqr5.f"
		if (v[m * v_dim1 + 1] != 0.) {
#line 550 "slaqr5.f"
		    k = krcol + (m - 1) * 3;
/* Computing MIN */
#line 551 "slaqr5.f"
		    i__6 = *kbot, i__7 = k + 3;
#line 551 "slaqr5.f"
		    i__4 = min(i__6,i__7);
#line 551 "slaqr5.f"
		    for (j = jtop; j <= i__4; ++j) {
#line 552 "slaqr5.f"
			refsum = v[m * v_dim1 + 1] * (h__[j + (k + 1) * 
				h_dim1] + v[m * v_dim1 + 2] * h__[j + (k + 2) 
				* h_dim1] + v[m * v_dim1 + 3] * h__[j + (k + 
				3) * h_dim1]);
#line 554 "slaqr5.f"
			h__[j + (k + 1) * h_dim1] -= refsum;
#line 555 "slaqr5.f"
			h__[j + (k + 2) * h_dim1] -= refsum * v[m * v_dim1 + 
				2];
#line 556 "slaqr5.f"
			h__[j + (k + 3) * h_dim1] -= refsum * v[m * v_dim1 + 
				3];
#line 557 "slaqr5.f"
/* L60: */
#line 557 "slaqr5.f"
		    }

#line 559 "slaqr5.f"
		    if (accum) {

/*                    ==== Accumulate U. (If necessary, update Z later */
/*                    .    with with an efficient matrix-matrix */
/*                    .    multiply.) ==== */

#line 565 "slaqr5.f"
			kms = k - incol;
/* Computing MAX */
#line 566 "slaqr5.f"
			i__4 = 1, i__6 = *ktop - incol;
#line 566 "slaqr5.f"
			i__7 = kdu;
#line 566 "slaqr5.f"
			for (j = max(i__4,i__6); j <= i__7; ++j) {
#line 567 "slaqr5.f"
			    refsum = v[m * v_dim1 + 1] * (u[j + (kms + 1) * 
				    u_dim1] + v[m * v_dim1 + 2] * u[j + (kms 
				    + 2) * u_dim1] + v[m * v_dim1 + 3] * u[j 
				    + (kms + 3) * u_dim1]);
#line 569 "slaqr5.f"
			    u[j + (kms + 1) * u_dim1] -= refsum;
#line 570 "slaqr5.f"
			    u[j + (kms + 2) * u_dim1] -= refsum * v[m * 
				    v_dim1 + 2];
#line 571 "slaqr5.f"
			    u[j + (kms + 3) * u_dim1] -= refsum * v[m * 
				    v_dim1 + 3];
#line 572 "slaqr5.f"
/* L70: */
#line 572 "slaqr5.f"
			}
#line 573 "slaqr5.f"
		    } else if (*wantz) {

/*                    ==== U is not accumulated, so update Z */
/*                    .    now by multiplying by reflections */
/*                    .    from the right. ==== */

#line 579 "slaqr5.f"
			i__7 = *ihiz;
#line 579 "slaqr5.f"
			for (j = *iloz; j <= i__7; ++j) {
#line 580 "slaqr5.f"
			    refsum = v[m * v_dim1 + 1] * (z__[j + (k + 1) * 
				    z_dim1] + v[m * v_dim1 + 2] * z__[j + (k 
				    + 2) * z_dim1] + v[m * v_dim1 + 3] * z__[
				    j + (k + 3) * z_dim1]);
#line 582 "slaqr5.f"
			    z__[j + (k + 1) * z_dim1] -= refsum;
#line 583 "slaqr5.f"
			    z__[j + (k + 2) * z_dim1] -= refsum * v[m * 
				    v_dim1 + 2];
#line 584 "slaqr5.f"
			    z__[j + (k + 3) * z_dim1] -= refsum * v[m * 
				    v_dim1 + 3];
#line 585 "slaqr5.f"
/* L80: */
#line 585 "slaqr5.f"
			}
#line 586 "slaqr5.f"
		    }
#line 587 "slaqr5.f"
		}
#line 588 "slaqr5.f"
/* L90: */
#line 588 "slaqr5.f"
	    }

/*           ==== Special case: 2-by-2 reflection (if needed) ==== */

#line 592 "slaqr5.f"
	    k = krcol + (m22 - 1) * 3;
#line 593 "slaqr5.f"
	    if (bmp22) {
#line 594 "slaqr5.f"
		if (v[m22 * v_dim1 + 1] != 0.) {
/* Computing MIN */
#line 595 "slaqr5.f"
		    i__7 = *kbot, i__4 = k + 3;
#line 595 "slaqr5.f"
		    i__5 = min(i__7,i__4);
#line 595 "slaqr5.f"
		    for (j = jtop; j <= i__5; ++j) {
#line 596 "slaqr5.f"
			refsum = v[m22 * v_dim1 + 1] * (h__[j + (k + 1) * 
				h_dim1] + v[m22 * v_dim1 + 2] * h__[j + (k + 
				2) * h_dim1]);
#line 598 "slaqr5.f"
			h__[j + (k + 1) * h_dim1] -= refsum;
#line 599 "slaqr5.f"
			h__[j + (k + 2) * h_dim1] -= refsum * v[m22 * v_dim1 
				+ 2];
#line 600 "slaqr5.f"
/* L100: */
#line 600 "slaqr5.f"
		    }

#line 602 "slaqr5.f"
		    if (accum) {
#line 603 "slaqr5.f"
			kms = k - incol;
/* Computing MAX */
#line 604 "slaqr5.f"
			i__5 = 1, i__7 = *ktop - incol;
#line 604 "slaqr5.f"
			i__4 = kdu;
#line 604 "slaqr5.f"
			for (j = max(i__5,i__7); j <= i__4; ++j) {
#line 605 "slaqr5.f"
			    refsum = v[m22 * v_dim1 + 1] * (u[j + (kms + 1) * 
				    u_dim1] + v[m22 * v_dim1 + 2] * u[j + (
				    kms + 2) * u_dim1]);
#line 607 "slaqr5.f"
			    u[j + (kms + 1) * u_dim1] -= refsum;
#line 608 "slaqr5.f"
			    u[j + (kms + 2) * u_dim1] -= refsum * v[m22 * 
				    v_dim1 + 2];
#line 610 "slaqr5.f"
/* L110: */
#line 610 "slaqr5.f"
			}
#line 611 "slaqr5.f"
		    } else if (*wantz) {
#line 612 "slaqr5.f"
			i__4 = *ihiz;
#line 612 "slaqr5.f"
			for (j = *iloz; j <= i__4; ++j) {
#line 613 "slaqr5.f"
			    refsum = v[m22 * v_dim1 + 1] * (z__[j + (k + 1) * 
				    z_dim1] + v[m22 * v_dim1 + 2] * z__[j + (
				    k + 2) * z_dim1]);
#line 615 "slaqr5.f"
			    z__[j + (k + 1) * z_dim1] -= refsum;
#line 616 "slaqr5.f"
			    z__[j + (k + 2) * z_dim1] -= refsum * v[m22 * 
				    v_dim1 + 2];
#line 617 "slaqr5.f"
/* L120: */
#line 617 "slaqr5.f"
			}
#line 618 "slaqr5.f"
		    }
#line 619 "slaqr5.f"
		}
#line 620 "slaqr5.f"
	    }

/*           ==== Vigilant deflation check ==== */

#line 624 "slaqr5.f"
	    mstart = mtop;
#line 625 "slaqr5.f"
	    if (krcol + (mstart - 1) * 3 < *ktop) {
#line 625 "slaqr5.f"
		++mstart;
#line 625 "slaqr5.f"
	    }
#line 627 "slaqr5.f"
	    mend = mbot;
#line 628 "slaqr5.f"
	    if (bmp22) {
#line 628 "slaqr5.f"
		++mend;
#line 628 "slaqr5.f"
	    }
#line 630 "slaqr5.f"
	    if (krcol == *kbot - 2) {
#line 630 "slaqr5.f"
		++mend;
#line 630 "slaqr5.f"
	    }
#line 632 "slaqr5.f"
	    i__4 = mend;
#line 632 "slaqr5.f"
	    for (m = mstart; m <= i__4; ++m) {
/* Computing MIN */
#line 633 "slaqr5.f"
		i__5 = *kbot - 1, i__7 = krcol + (m - 1) * 3;
#line 633 "slaqr5.f"
		k = min(i__5,i__7);

/*              ==== The following convergence test requires that */
/*              .    the tradition small-compared-to-nearby-diagonals */
/*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997) */
/*              .    criteria both be satisfied.  The latter improves */
/*              .    accuracy in some examples. Falling back on an */
/*              .    alternate convergence criterion when TST1 or TST2 */
/*              .    is zero (as done here) is traditional but probably */
/*              .    unnecessary. ==== */

#line 644 "slaqr5.f"
		if (h__[k + 1 + k * h_dim1] != 0.) {
#line 645 "slaqr5.f"
		    tst1 = (d__1 = h__[k + k * h_dim1], abs(d__1)) + (d__2 = 
			    h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
#line 646 "slaqr5.f"
		    if (tst1 == 0.) {
#line 647 "slaqr5.f"
			if (k >= *ktop + 1) {
#line 647 "slaqr5.f"
			    tst1 += (d__1 = h__[k + (k - 1) * h_dim1], abs(
				    d__1));
#line 647 "slaqr5.f"
			}
#line 649 "slaqr5.f"
			if (k >= *ktop + 2) {
#line 649 "slaqr5.f"
			    tst1 += (d__1 = h__[k + (k - 2) * h_dim1], abs(
				    d__1));
#line 649 "slaqr5.f"
			}
#line 651 "slaqr5.f"
			if (k >= *ktop + 3) {
#line 651 "slaqr5.f"
			    tst1 += (d__1 = h__[k + (k - 3) * h_dim1], abs(
				    d__1));
#line 651 "slaqr5.f"
			}
#line 653 "slaqr5.f"
			if (k <= *kbot - 2) {
#line 653 "slaqr5.f"
			    tst1 += (d__1 = h__[k + 2 + (k + 1) * h_dim1], 
				    abs(d__1));
#line 653 "slaqr5.f"
			}
#line 655 "slaqr5.f"
			if (k <= *kbot - 3) {
#line 655 "slaqr5.f"
			    tst1 += (d__1 = h__[k + 3 + (k + 1) * h_dim1], 
				    abs(d__1));
#line 655 "slaqr5.f"
			}
#line 657 "slaqr5.f"
			if (k <= *kbot - 4) {
#line 657 "slaqr5.f"
			    tst1 += (d__1 = h__[k + 4 + (k + 1) * h_dim1], 
				    abs(d__1));
#line 657 "slaqr5.f"
			}
#line 659 "slaqr5.f"
		    }
/* Computing MAX */
#line 660 "slaqr5.f"
		    d__2 = smlnum, d__3 = ulp * tst1;
#line 660 "slaqr5.f"
		    if ((d__1 = h__[k + 1 + k * h_dim1], abs(d__1)) <= max(
			    d__2,d__3)) {
/* Computing MAX */
#line 662 "slaqr5.f"
			d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)), 
				d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(
				d__2));
#line 662 "slaqr5.f"
			h12 = max(d__3,d__4);
/* Computing MIN */
#line 663 "slaqr5.f"
			d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)), 
				d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(
				d__2));
#line 663 "slaqr5.f"
			h21 = min(d__3,d__4);
/* Computing MAX */
#line 664 "slaqr5.f"
			d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(
				d__1)), d__4 = (d__2 = h__[k + k * h_dim1] - 
				h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
#line 664 "slaqr5.f"
			h11 = max(d__3,d__4);
/* Computing MIN */
#line 666 "slaqr5.f"
			d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(
				d__1)), d__4 = (d__2 = h__[k + k * h_dim1] - 
				h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
#line 666 "slaqr5.f"
			h22 = min(d__3,d__4);
#line 668 "slaqr5.f"
			scl = h11 + h12;
#line 669 "slaqr5.f"
			tst2 = h22 * (h11 / scl);

/* Computing MAX */
#line 671 "slaqr5.f"
			d__1 = smlnum, d__2 = ulp * tst2;
#line 671 "slaqr5.f"
			if (tst2 == 0. || h21 * (h12 / scl) <= max(d__1,d__2))
				 {
#line 671 "slaqr5.f"
			    h__[k + 1 + k * h_dim1] = 0.;
#line 671 "slaqr5.f"
			}
#line 673 "slaqr5.f"
		    }
#line 674 "slaqr5.f"
		}
#line 675 "slaqr5.f"
/* L130: */
#line 675 "slaqr5.f"
	    }

/*           ==== Fill in the last row of each bulge. ==== */

/* Computing MIN */
#line 679 "slaqr5.f"
	    i__4 = nbmps, i__5 = (*kbot - krcol - 1) / 3;
#line 679 "slaqr5.f"
	    mend = min(i__4,i__5);
#line 680 "slaqr5.f"
	    i__4 = mend;
#line 680 "slaqr5.f"
	    for (m = mtop; m <= i__4; ++m) {
#line 681 "slaqr5.f"
		k = krcol + (m - 1) * 3;
#line 682 "slaqr5.f"
		refsum = v[m * v_dim1 + 1] * v[m * v_dim1 + 3] * h__[k + 4 + (
			k + 3) * h_dim1];
#line 683 "slaqr5.f"
		h__[k + 4 + (k + 1) * h_dim1] = -refsum;
#line 684 "slaqr5.f"
		h__[k + 4 + (k + 2) * h_dim1] = -refsum * v[m * v_dim1 + 2];
#line 685 "slaqr5.f"
		h__[k + 4 + (k + 3) * h_dim1] -= refsum * v[m * v_dim1 + 3];
#line 686 "slaqr5.f"
/* L140: */
#line 686 "slaqr5.f"
	    }

/*           ==== End of near-the-diagonal bulge chase. ==== */

#line 690 "slaqr5.f"
/* L150: */
#line 690 "slaqr5.f"
	}

/*        ==== Use U (if accumulated) to update far-from-diagonal */
/*        .    entries in H.  If required, use U to update Z as */
/*        .    well. ==== */

#line 696 "slaqr5.f"
	if (accum) {
#line 697 "slaqr5.f"
	    if (*wantt) {
#line 698 "slaqr5.f"
		jtop = 1;
#line 699 "slaqr5.f"
		jbot = *n;
#line 700 "slaqr5.f"
	    } else {
#line 701 "slaqr5.f"
		jtop = *ktop;
#line 702 "slaqr5.f"
		jbot = *kbot;
#line 703 "slaqr5.f"
	    }
#line 704 "slaqr5.f"
	    if (! blk22 || incol < *ktop || ndcol > *kbot || ns <= 2) {

/*              ==== Updates not exploiting the 2-by-2 block */
/*              .    structure of U.  K1 and NU keep track of */
/*              .    the location and size of U in the special */
/*              .    cases of introducing bulges and chasing */
/*              .    bulges off the bottom.  In these special */
/*              .    cases and in case the number of shifts */
/*              .    is NS = 2, there is no 2-by-2 block */
/*              .    structure to exploit.  ==== */

/* Computing MAX */
#line 716 "slaqr5.f"
		i__3 = 1, i__4 = *ktop - incol;
#line 716 "slaqr5.f"
		k1 = max(i__3,i__4);
/* Computing MAX */
#line 717 "slaqr5.f"
		i__3 = 0, i__4 = ndcol - *kbot;
#line 717 "slaqr5.f"
		nu = kdu - max(i__3,i__4) - k1 + 1;

/*              ==== Horizontal Multiply ==== */

#line 721 "slaqr5.f"
		i__3 = jbot;
#line 721 "slaqr5.f"
		i__4 = *nh;
#line 721 "slaqr5.f"
		for (jcol = min(ndcol,*kbot) + 1; i__4 < 0 ? jcol >= i__3 : 
			jcol <= i__3; jcol += i__4) {
/* Computing MIN */
#line 722 "slaqr5.f"
		    i__5 = *nh, i__7 = jbot - jcol + 1;
#line 722 "slaqr5.f"
		    jlen = min(i__5,i__7);
#line 723 "slaqr5.f"
		    sgemm_("C", "N", &nu, &jlen, &nu, &c_b8, &u[k1 + k1 * 
			    u_dim1], ldu, &h__[incol + k1 + jcol * h_dim1], 
			    ldh, &c_b7, &wh[wh_offset], ldwh, (ftnlen)1, (
			    ftnlen)1);
#line 726 "slaqr5.f"
		    slacpy_("ALL", &nu, &jlen, &wh[wh_offset], ldwh, &h__[
			    incol + k1 + jcol * h_dim1], ldh, (ftnlen)3);
#line 728 "slaqr5.f"
/* L160: */
#line 728 "slaqr5.f"
		}

/*              ==== Vertical multiply ==== */

#line 732 "slaqr5.f"
		i__4 = max(*ktop,incol) - 1;
#line 732 "slaqr5.f"
		i__3 = *nv;
#line 732 "slaqr5.f"
		for (jrow = jtop; i__3 < 0 ? jrow >= i__4 : jrow <= i__4; 
			jrow += i__3) {
/* Computing MIN */
#line 733 "slaqr5.f"
		    i__5 = *nv, i__7 = max(*ktop,incol) - jrow;
#line 733 "slaqr5.f"
		    jlen = min(i__5,i__7);
#line 734 "slaqr5.f"
		    sgemm_("N", "N", &jlen, &nu, &nu, &c_b8, &h__[jrow + (
			    incol + k1) * h_dim1], ldh, &u[k1 + k1 * u_dim1], 
			    ldu, &c_b7, &wv[wv_offset], ldwv, (ftnlen)1, (
			    ftnlen)1);
#line 737 "slaqr5.f"
		    slacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &h__[
			    jrow + (incol + k1) * h_dim1], ldh, (ftnlen)3);
#line 739 "slaqr5.f"
/* L170: */
#line 739 "slaqr5.f"
		}

/*              ==== Z multiply (also vertical) ==== */

#line 743 "slaqr5.f"
		if (*wantz) {
#line 744 "slaqr5.f"
		    i__3 = *ihiz;
#line 744 "slaqr5.f"
		    i__4 = *nv;
#line 744 "slaqr5.f"
		    for (jrow = *iloz; i__4 < 0 ? jrow >= i__3 : jrow <= i__3;
			     jrow += i__4) {
/* Computing MIN */
#line 745 "slaqr5.f"
			i__5 = *nv, i__7 = *ihiz - jrow + 1;
#line 745 "slaqr5.f"
			jlen = min(i__5,i__7);
#line 746 "slaqr5.f"
			sgemm_("N", "N", &jlen, &nu, &nu, &c_b8, &z__[jrow + (
				incol + k1) * z_dim1], ldz, &u[k1 + k1 * 
				u_dim1], ldu, &c_b7, &wv[wv_offset], ldwv, (
				ftnlen)1, (ftnlen)1);
#line 749 "slaqr5.f"
			slacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &z__[
				jrow + (incol + k1) * z_dim1], ldz, (ftnlen)3)
				;
#line 751 "slaqr5.f"
/* L180: */
#line 751 "slaqr5.f"
		    }
#line 752 "slaqr5.f"
		}
#line 753 "slaqr5.f"
	    } else {

/*              ==== Updates exploiting U's 2-by-2 block structure. */
/*              .    (I2, I4, J2, J4 are the last rows and columns */
/*              .    of the blocks.) ==== */

#line 759 "slaqr5.f"
		i2 = (kdu + 1) / 2;
#line 760 "slaqr5.f"
		i4 = kdu;
#line 761 "slaqr5.f"
		j2 = i4 - i2;
#line 762 "slaqr5.f"
		j4 = kdu;

/*              ==== KZS and KNZ deal with the band of zeros */
/*              .    along the diagonal of one of the triangular */
/*              .    blocks. ==== */

#line 768 "slaqr5.f"
		kzs = j4 - j2 - (ns + 1);
#line 769 "slaqr5.f"
		knz = ns + 1;

/*              ==== Horizontal multiply ==== */

#line 773 "slaqr5.f"
		i__4 = jbot;
#line 773 "slaqr5.f"
		i__3 = *nh;
#line 773 "slaqr5.f"
		for (jcol = min(ndcol,*kbot) + 1; i__3 < 0 ? jcol >= i__4 : 
			jcol <= i__4; jcol += i__3) {
/* Computing MIN */
#line 774 "slaqr5.f"
		    i__5 = *nh, i__7 = jbot - jcol + 1;
#line 774 "slaqr5.f"
		    jlen = min(i__5,i__7);

/*                 ==== Copy bottom of H to top+KZS of scratch ==== */
/*                  (The first KZS rows get multiplied by zero.) ==== */

#line 779 "slaqr5.f"
		    slacpy_("ALL", &knz, &jlen, &h__[incol + 1 + j2 + jcol * 
			    h_dim1], ldh, &wh[kzs + 1 + wh_dim1], ldwh, (
			    ftnlen)3);

/*                 ==== Multiply by U21**T ==== */

#line 784 "slaqr5.f"
		    slaset_("ALL", &kzs, &jlen, &c_b7, &c_b7, &wh[wh_offset], 
			    ldwh, (ftnlen)3);
#line 785 "slaqr5.f"
		    strmm_("L", "U", "C", "N", &knz, &jlen, &c_b8, &u[j2 + 1 
			    + (kzs + 1) * u_dim1], ldu, &wh[kzs + 1 + wh_dim1]
			    , ldwh, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1);

/*                 ==== Multiply top of H by U11**T ==== */

#line 791 "slaqr5.f"
		    sgemm_("C", "N", &i2, &jlen, &j2, &c_b8, &u[u_offset], 
			    ldu, &h__[incol + 1 + jcol * h_dim1], ldh, &c_b8, 
			    &wh[wh_offset], ldwh, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy top of H to bottom of WH ==== */

#line 796 "slaqr5.f"
		    slacpy_("ALL", &j2, &jlen, &h__[incol + 1 + jcol * h_dim1]
			    , ldh, &wh[i2 + 1 + wh_dim1], ldwh, (ftnlen)3);

/*                 ==== Multiply by U21**T ==== */

#line 801 "slaqr5.f"
		    strmm_("L", "L", "C", "N", &j2, &jlen, &c_b8, &u[(i2 + 1) 
			    * u_dim1 + 1], ldu, &wh[i2 + 1 + wh_dim1], ldwh, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 ==== Multiply by U22 ==== */

#line 806 "slaqr5.f"
		    i__5 = i4 - i2;
#line 806 "slaqr5.f"
		    i__7 = j4 - j2;
#line 806 "slaqr5.f"
		    sgemm_("C", "N", &i__5, &jlen, &i__7, &c_b8, &u[j2 + 1 + (
			    i2 + 1) * u_dim1], ldu, &h__[incol + 1 + j2 + 
			    jcol * h_dim1], ldh, &c_b8, &wh[i2 + 1 + wh_dim1],
			     ldwh, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy it back ==== */

#line 813 "slaqr5.f"
		    slacpy_("ALL", &kdu, &jlen, &wh[wh_offset], ldwh, &h__[
			    incol + 1 + jcol * h_dim1], ldh, (ftnlen)3);
#line 815 "slaqr5.f"
/* L190: */
#line 815 "slaqr5.f"
		}

/*              ==== Vertical multiply ==== */

#line 819 "slaqr5.f"
		i__3 = max(incol,*ktop) - 1;
#line 819 "slaqr5.f"
		i__4 = *nv;
#line 819 "slaqr5.f"
		for (jrow = jtop; i__4 < 0 ? jrow >= i__3 : jrow <= i__3; 
			jrow += i__4) {
/* Computing MIN */
#line 820 "slaqr5.f"
		    i__5 = *nv, i__7 = max(incol,*ktop) - jrow;
#line 820 "slaqr5.f"
		    jlen = min(i__5,i__7);

/*                 ==== Copy right of H to scratch (the first KZS */
/*                 .    columns get multiplied by zero) ==== */

#line 825 "slaqr5.f"
		    slacpy_("ALL", &jlen, &knz, &h__[jrow + (incol + 1 + j2) *
			     h_dim1], ldh, &wv[(kzs + 1) * wv_dim1 + 1], ldwv,
			     (ftnlen)3);

/*                 ==== Multiply by U21 ==== */

#line 830 "slaqr5.f"
		    slaset_("ALL", &jlen, &kzs, &c_b7, &c_b7, &wv[wv_offset], 
			    ldwv, (ftnlen)3);
#line 831 "slaqr5.f"
		    strmm_("R", "U", "N", "N", &jlen, &knz, &c_b8, &u[j2 + 1 
			    + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) * 
			    wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1);

/*                 ==== Multiply by U11 ==== */

#line 837 "slaqr5.f"
		    sgemm_("N", "N", &jlen, &i2, &j2, &c_b8, &h__[jrow + (
			    incol + 1) * h_dim1], ldh, &u[u_offset], ldu, &
			    c_b8, &wv[wv_offset], ldwv, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy left of H to right of scratch ==== */

#line 843 "slaqr5.f"
		    slacpy_("ALL", &jlen, &j2, &h__[jrow + (incol + 1) * 
			    h_dim1], ldh, &wv[(i2 + 1) * wv_dim1 + 1], ldwv, (
			    ftnlen)3);

/*                 ==== Multiply by U21 ==== */

#line 848 "slaqr5.f"
		    i__5 = i4 - i2;
#line 848 "slaqr5.f"
		    strmm_("R", "L", "N", "N", &jlen, &i__5, &c_b8, &u[(i2 + 
			    1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * wv_dim1 + 1]
			    , ldwv, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1);

/*                 ==== Multiply by U22 ==== */

#line 853 "slaqr5.f"
		    i__5 = i4 - i2;
#line 853 "slaqr5.f"
		    i__7 = j4 - j2;
#line 853 "slaqr5.f"
		    sgemm_("N", "N", &jlen, &i__5, &i__7, &c_b8, &h__[jrow + (
			    incol + 1 + j2) * h_dim1], ldh, &u[j2 + 1 + (i2 + 
			    1) * u_dim1], ldu, &c_b8, &wv[(i2 + 1) * wv_dim1 
			    + 1], ldwv, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy it back ==== */

#line 860 "slaqr5.f"
		    slacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, &h__[
			    jrow + (incol + 1) * h_dim1], ldh, (ftnlen)3);
#line 862 "slaqr5.f"
/* L200: */
#line 862 "slaqr5.f"
		}

/*              ==== Multiply Z (also vertical) ==== */

#line 866 "slaqr5.f"
		if (*wantz) {
#line 867 "slaqr5.f"
		    i__4 = *ihiz;
#line 867 "slaqr5.f"
		    i__3 = *nv;
#line 867 "slaqr5.f"
		    for (jrow = *iloz; i__3 < 0 ? jrow >= i__4 : jrow <= i__4;
			     jrow += i__3) {
/* Computing MIN */
#line 868 "slaqr5.f"
			i__5 = *nv, i__7 = *ihiz - jrow + 1;
#line 868 "slaqr5.f"
			jlen = min(i__5,i__7);

/*                    ==== Copy right of Z to left of scratch (first */
/*                    .     KZS columns get multiplied by zero) ==== */

#line 873 "slaqr5.f"
			slacpy_("ALL", &jlen, &knz, &z__[jrow + (incol + 1 + 
				j2) * z_dim1], ldz, &wv[(kzs + 1) * wv_dim1 + 
				1], ldwv, (ftnlen)3);

/*                    ==== Multiply by U12 ==== */

#line 879 "slaqr5.f"
			slaset_("ALL", &jlen, &kzs, &c_b7, &c_b7, &wv[
				wv_offset], ldwv, (ftnlen)3);
#line 881 "slaqr5.f"
			strmm_("R", "U", "N", "N", &jlen, &knz, &c_b8, &u[j2 
				+ 1 + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) 
				* wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

/*                    ==== Multiply by U11 ==== */

#line 887 "slaqr5.f"
			sgemm_("N", "N", &jlen, &i2, &j2, &c_b8, &z__[jrow + (
				incol + 1) * z_dim1], ldz, &u[u_offset], ldu, 
				&c_b8, &wv[wv_offset], ldwv, (ftnlen)1, (
				ftnlen)1);

/*                    ==== Copy left of Z to right of scratch ==== */

#line 893 "slaqr5.f"
			slacpy_("ALL", &jlen, &j2, &z__[jrow + (incol + 1) * 
				z_dim1], ldz, &wv[(i2 + 1) * wv_dim1 + 1], 
				ldwv, (ftnlen)3);

/*                    ==== Multiply by U21 ==== */

#line 898 "slaqr5.f"
			i__5 = i4 - i2;
#line 898 "slaqr5.f"
			strmm_("R", "L", "N", "N", &jlen, &i__5, &c_b8, &u[(
				i2 + 1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * 
				wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

/*                    ==== Multiply by U22 ==== */

#line 904 "slaqr5.f"
			i__5 = i4 - i2;
#line 904 "slaqr5.f"
			i__7 = j4 - j2;
#line 904 "slaqr5.f"
			sgemm_("N", "N", &jlen, &i__5, &i__7, &c_b8, &z__[
				jrow + (incol + 1 + j2) * z_dim1], ldz, &u[j2 
				+ 1 + (i2 + 1) * u_dim1], ldu, &c_b8, &wv[(i2 
				+ 1) * wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)
				1);

/*                    ==== Copy the result back to Z ==== */

#line 911 "slaqr5.f"
			slacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, &
				z__[jrow + (incol + 1) * z_dim1], ldz, (
				ftnlen)3);
#line 913 "slaqr5.f"
/* L210: */
#line 913 "slaqr5.f"
		    }
#line 914 "slaqr5.f"
		}
#line 915 "slaqr5.f"
	    }
#line 916 "slaqr5.f"
	}
#line 917 "slaqr5.f"
/* L220: */
#line 917 "slaqr5.f"
    }

/*     ==== End of SLAQR5 ==== */

#line 921 "slaqr5.f"
    return 0;
} /* slaqr5_ */

