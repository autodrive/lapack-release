#line 1 "claqr5.f"
/* claqr5.f -- translated by f2c (version 20100827).
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

#line 1 "claqr5.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b CLAQR5 performs a single small-bulge multi-shift QR sweep. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAQR5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, */
/*                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, */
/*                          WV, LDWV, NH, WH, LDWH ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, */
/*      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), */
/*      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLAQR5 called by CLAQR0 performs a */
/* >    single small-bulge multi-shift QR sweep. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >             WANTT = .true. if the triangular Schur factor */
/* >             is being computed.  WANTT is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >             WANTZ = .true. if the unitary Schur factor is being */
/* >             computed.  WANTZ is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] KACC22 */
/* > \verbatim */
/* >          KACC22 is INTEGER with value 0, 1, or 2. */
/* >             Specifies the computation mode of far-from-diagonal */
/* >             orthogonal updates. */
/* >        = 0: CLAQR5 does not accumulate reflections and does not */
/* >             use matrix-matrix multiply to update far-from-diagonal */
/* >             matrix entries. */
/* >        = 1: CLAQR5 accumulates reflections and uses matrix-matrix */
/* >             multiply to update the far-from-diagonal matrix entries. */
/* >        = 2: CLAQR5 accumulates reflections, uses matrix-matrix */
/* >             multiply to update the far-from-diagonal matrix entries, */
/* >             and takes advantage of 2-by-2 block structure during */
/* >             matrix multiplies. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >             N is the order of the Hessenberg matrix H upon which this */
/* >             subroutine operates. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* >          KTOP is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* >          KBOT is INTEGER */
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
/* >          NSHFTS is INTEGER */
/* >             NSHFTS gives the number of simultaneous shifts.  NSHFTS */
/* >             must be positive and even. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* >          S is COMPLEX array, dimension (NSHFTS) */
/* >             S contains the shifts of origin that define the multi- */
/* >             shift QR sweep.  On output S may be reordered. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is COMPLEX array, dimension (LDH,N) */
/* >             On input H contains a Hessenberg matrix.  On output a */
/* >             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied */
/* >             to the isolated diagonal block in rows and columns KTOP */
/* >             through KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
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
/* >          Z is COMPLEX array, dimension (LDZ,IHIZ) */
/* >             If WANTZ = .TRUE., then the QR Sweep unitary */
/* >             similarity transformation is accumulated into */
/* >             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. */
/* >             If WANTZ = .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >             LDA is the leading dimension of Z just as declared in */
/* >             the calling procedure. LDZ.GE.N. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,NSHFTS/2) */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >             LDV is the leading dimension of V as declared in the */
/* >             calling procedure.  LDV.GE.3. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX array, dimension (LDU,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >             LDU is the leading dimension of U just as declared in the */
/* >             in the calling subroutine.  LDU.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* >          NH is INTEGER */
/* >             NH is the number of columns in array WH available for */
/* >             workspace. NH.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WH */
/* > \verbatim */
/* >          WH is COMPLEX array, dimension (LDWH,NH) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWH */
/* > \verbatim */
/* >          LDWH is INTEGER */
/* >             Leading dimension of WH just as declared in the */
/* >             calling procedure.  LDWH.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* >          NV is INTEGER */
/* >             NV is the number of rows in WV agailable for workspace. */
/* >             NV.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* >          WV is COMPLEX array, dimension (LDWV,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* >          LDWV is INTEGER */
/* >             LDWV is the leading dimension of WV as declared in the */
/* >             in the calling subroutine.  LDWV.GE.NV. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup complexOTHERauxiliary */

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
/* Subroutine */ int claqr5_(logical *wantt, logical *wantz, integer *kacc22, 
	integer *n, integer *ktop, integer *kbot, integer *nshfts, 
	doublecomplex *s, doublecomplex *h__, integer *ldh, integer *iloz, 
	integer *ihiz, doublecomplex *z__, integer *ldz, doublecomplex *v, 
	integer *ldv, doublecomplex *u, integer *ldu, integer *nv, 
	doublecomplex *wv, integer *ldwv, integer *nh, doublecomplex *wh, 
	integer *ldwh)
{
    /* System generated locals */
    integer h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, 
	    wh_offset, wv_dim1, wv_offset, z_dim1, z_offset, i__1, i__2, i__3,
	     i__4, i__5, i__6, i__7, i__8, i__9, i__10, i__11;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer j, k, m, i2, j2, i4, j4, k1;
    static doublereal h11, h12, h21, h22;
    static integer m22, ns, nu;
    static doublecomplex vt[3];
    static doublereal scl;
    static integer kdu, kms;
    static doublereal ulp;
    static integer knz, kzs;
    static doublereal tst1, tst2;
    static doublecomplex beta;
    static logical blk22, bmp22;
    static integer mend, jcol, jlen, jbot, mbot, jtop, jrow, mtop;
    static doublecomplex alpha;
    static logical accum;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer ndcol, incol, krcol, nbmps;
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    claqr1_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *), slabad_(doublereal *, 
	    doublereal *), clarfg_(integer *, doublecomplex *, doublecomplex *
	    , integer *, doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    claset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal safmin, safmax;
    static doublecomplex refsum;
    static integer mstart;
    static doublereal smlnum;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== If there are no shifts, then there is nothing to do. ==== */

#line 310 "claqr5.f"
    /* Parameter adjustments */
#line 310 "claqr5.f"
    --s;
#line 310 "claqr5.f"
    h_dim1 = *ldh;
#line 310 "claqr5.f"
    h_offset = 1 + h_dim1;
#line 310 "claqr5.f"
    h__ -= h_offset;
#line 310 "claqr5.f"
    z_dim1 = *ldz;
#line 310 "claqr5.f"
    z_offset = 1 + z_dim1;
#line 310 "claqr5.f"
    z__ -= z_offset;
#line 310 "claqr5.f"
    v_dim1 = *ldv;
#line 310 "claqr5.f"
    v_offset = 1 + v_dim1;
#line 310 "claqr5.f"
    v -= v_offset;
#line 310 "claqr5.f"
    u_dim1 = *ldu;
#line 310 "claqr5.f"
    u_offset = 1 + u_dim1;
#line 310 "claqr5.f"
    u -= u_offset;
#line 310 "claqr5.f"
    wv_dim1 = *ldwv;
#line 310 "claqr5.f"
    wv_offset = 1 + wv_dim1;
#line 310 "claqr5.f"
    wv -= wv_offset;
#line 310 "claqr5.f"
    wh_dim1 = *ldwh;
#line 310 "claqr5.f"
    wh_offset = 1 + wh_dim1;
#line 310 "claqr5.f"
    wh -= wh_offset;
#line 310 "claqr5.f"

#line 310 "claqr5.f"
    /* Function Body */
#line 310 "claqr5.f"
    if (*nshfts < 2) {
#line 310 "claqr5.f"
	return 0;
#line 310 "claqr5.f"
    }

/*     ==== If the active block is empty or 1-by-1, then there */
/*     .    is nothing to do. ==== */

#line 316 "claqr5.f"
    if (*ktop >= *kbot) {
#line 316 "claqr5.f"
	return 0;
#line 316 "claqr5.f"
    }

/*     ==== NSHFTS is supposed to be even, but if it is odd, */
/*     .    then simply reduce it by one.  ==== */

#line 322 "claqr5.f"
    ns = *nshfts - *nshfts % 2;

/*     ==== Machine constants for deflation ==== */

#line 326 "claqr5.f"
    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
#line 327 "claqr5.f"
    safmax = 1. / safmin;
#line 328 "claqr5.f"
    slabad_(&safmin, &safmax);
#line 329 "claqr5.f"
    ulp = slamch_("PRECISION", (ftnlen)9);
#line 330 "claqr5.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Use accumulated reflections to update far-from-diagonal */
/*     .    entries ? ==== */

#line 335 "claqr5.f"
    accum = *kacc22 == 1 || *kacc22 == 2;

/*     ==== If so, exploit the 2-by-2 block structure? ==== */

#line 339 "claqr5.f"
    blk22 = ns > 2 && *kacc22 == 2;

/*     ==== clear trash ==== */

#line 343 "claqr5.f"
    if (*ktop + 2 <= *kbot) {
#line 343 "claqr5.f"
	i__1 = *ktop + 2 + *ktop * h_dim1;
#line 343 "claqr5.f"
	h__[i__1].r = 0., h__[i__1].i = 0.;
#line 343 "claqr5.f"
    }

/*     ==== NBMPS = number of 2-shift bulges in the chain ==== */

#line 348 "claqr5.f"
    nbmps = ns / 2;

/*     ==== KDU = width of slab ==== */

#line 352 "claqr5.f"
    kdu = nbmps * 6 - 3;

/*     ==== Create and chase chains of NBMPS bulges ==== */

#line 356 "claqr5.f"
    i__1 = *kbot - 2;
#line 356 "claqr5.f"
    i__2 = nbmps * 3 - 2;
#line 356 "claqr5.f"
    for (incol = (1 - nbmps) * 3 + *ktop - 1; i__2 < 0 ? incol >= i__1 : 
	    incol <= i__1; incol += i__2) {
#line 357 "claqr5.f"
	ndcol = incol + kdu;
#line 358 "claqr5.f"
	if (accum) {
#line 358 "claqr5.f"
	    claset_("ALL", &kdu, &kdu, &c_b1, &c_b2, &u[u_offset], ldu, (
		    ftnlen)3);
#line 358 "claqr5.f"
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
#line 373 "claqr5.f"
	i__4 = incol + nbmps * 3 - 3, i__5 = *kbot - 2;
#line 373 "claqr5.f"
	i__3 = min(i__4,i__5);
#line 373 "claqr5.f"
	for (krcol = incol; krcol <= i__3; ++krcol) {

/*           ==== Bulges number MTOP to MBOT are active double implicit */
/*           .    shift bulges.  There may or may not also be small */
/*           .    2-by-2 bulge, if there is room.  The inactive bulges */
/*           .    (if any) must wait until the active bulges have moved */
/*           .    down the diagonal to make room.  The phantom matrix */
/*           .    paradigm described above helps keep track.  ==== */

/* Computing MAX */
#line 382 "claqr5.f"
	    i__4 = 1, i__5 = (*ktop - 1 - krcol + 2) / 3 + 1;
#line 382 "claqr5.f"
	    mtop = max(i__4,i__5);
/* Computing MIN */
#line 383 "claqr5.f"
	    i__4 = nbmps, i__5 = (*kbot - krcol) / 3;
#line 383 "claqr5.f"
	    mbot = min(i__4,i__5);
#line 384 "claqr5.f"
	    m22 = mbot + 1;
#line 385 "claqr5.f"
	    bmp22 = mbot < nbmps && krcol + (m22 - 1) * 3 == *kbot - 2;

/*           ==== Generate reflections to chase the chain right */
/*           .    one column.  (The minimum value of K is KTOP-1.) ==== */

#line 391 "claqr5.f"
	    i__4 = mbot;
#line 391 "claqr5.f"
	    for (m = mtop; m <= i__4; ++m) {
#line 392 "claqr5.f"
		k = krcol + (m - 1) * 3;
#line 393 "claqr5.f"
		if (k == *ktop - 1) {
#line 394 "claqr5.f"
		    claqr1_(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &s[(m <<
			     1) - 1], &s[m * 2], &v[m * v_dim1 + 1]);
#line 396 "claqr5.f"
		    i__5 = m * v_dim1 + 1;
#line 396 "claqr5.f"
		    alpha.r = v[i__5].r, alpha.i = v[i__5].i;
#line 397 "claqr5.f"
		    clarfg_(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1, &v[m * 
			    v_dim1 + 1]);
#line 398 "claqr5.f"
		} else {
#line 399 "claqr5.f"
		    i__5 = k + 1 + k * h_dim1;
#line 399 "claqr5.f"
		    beta.r = h__[i__5].r, beta.i = h__[i__5].i;
#line 400 "claqr5.f"
		    i__5 = m * v_dim1 + 2;
#line 400 "claqr5.f"
		    i__6 = k + 2 + k * h_dim1;
#line 400 "claqr5.f"
		    v[i__5].r = h__[i__6].r, v[i__5].i = h__[i__6].i;
#line 401 "claqr5.f"
		    i__5 = m * v_dim1 + 3;
#line 401 "claqr5.f"
		    i__6 = k + 3 + k * h_dim1;
#line 401 "claqr5.f"
		    v[i__5].r = h__[i__6].r, v[i__5].i = h__[i__6].i;
#line 402 "claqr5.f"
		    clarfg_(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * 
			    v_dim1 + 1]);

/*                 ==== A Bulge may collapse because of vigilant */
/*                 .    deflation or destructive underflow.  In the */
/*                 .    underflow case, try the two-small-subdiagonals */
/*                 .    trick to try to reinflate the bulge.  ==== */

#line 409 "claqr5.f"
		    i__5 = k + 3 + k * h_dim1;
#line 409 "claqr5.f"
		    i__6 = k + 3 + (k + 1) * h_dim1;
#line 409 "claqr5.f"
		    i__7 = k + 3 + (k + 2) * h_dim1;
#line 409 "claqr5.f"
		    if (h__[i__5].r != 0. || h__[i__5].i != 0. || (h__[i__6]
			    .r != 0. || h__[i__6].i != 0.) || h__[i__7].r == 
			    0. && h__[i__7].i == 0.) {

/*                    ==== Typical case: not collapsed (yet). ==== */

#line 414 "claqr5.f"
			i__5 = k + 1 + k * h_dim1;
#line 414 "claqr5.f"
			h__[i__5].r = beta.r, h__[i__5].i = beta.i;
#line 415 "claqr5.f"
			i__5 = k + 2 + k * h_dim1;
#line 415 "claqr5.f"
			h__[i__5].r = 0., h__[i__5].i = 0.;
#line 416 "claqr5.f"
			i__5 = k + 3 + k * h_dim1;
#line 416 "claqr5.f"
			h__[i__5].r = 0., h__[i__5].i = 0.;
#line 417 "claqr5.f"
		    } else {

/*                    ==== Atypical case: collapsed.  Attempt to */
/*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K). */
/*                    .    If the fill resulting from the new */
/*                    .    reflector is too large, then abandon it. */
/*                    .    Otherwise, use the new one. ==== */

#line 425 "claqr5.f"
			claqr1_(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh, &
				s[(m << 1) - 1], &s[m * 2], vt);
#line 427 "claqr5.f"
			alpha.r = vt[0].r, alpha.i = vt[0].i;
#line 428 "claqr5.f"
			clarfg_(&c__3, &alpha, &vt[1], &c__1, vt);
#line 429 "claqr5.f"
			d_cnjg(&z__2, vt);
#line 429 "claqr5.f"
			i__5 = k + 1 + k * h_dim1;
#line 429 "claqr5.f"
			d_cnjg(&z__5, &vt[1]);
#line 429 "claqr5.f"
			i__6 = k + 2 + k * h_dim1;
#line 429 "claqr5.f"
			z__4.r = z__5.r * h__[i__6].r - z__5.i * h__[i__6].i, 
				z__4.i = z__5.r * h__[i__6].i + z__5.i * h__[
				i__6].r;
#line 429 "claqr5.f"
			z__3.r = h__[i__5].r + z__4.r, z__3.i = h__[i__5].i + 
				z__4.i;
#line 429 "claqr5.f"
			z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
				z__2.r * z__3.i + z__2.i * z__3.r;
#line 429 "claqr5.f"
			refsum.r = z__1.r, refsum.i = z__1.i;

#line 433 "claqr5.f"
			i__5 = k + 2 + k * h_dim1;
#line 433 "claqr5.f"
			z__3.r = refsum.r * vt[1].r - refsum.i * vt[1].i, 
				z__3.i = refsum.r * vt[1].i + refsum.i * vt[1]
				.r;
#line 433 "claqr5.f"
			z__2.r = h__[i__5].r - z__3.r, z__2.i = h__[i__5].i - 
				z__3.i;
#line 433 "claqr5.f"
			z__1.r = z__2.r, z__1.i = z__2.i;
#line 433 "claqr5.f"
			z__5.r = refsum.r * vt[2].r - refsum.i * vt[2].i, 
				z__5.i = refsum.r * vt[2].i + refsum.i * vt[2]
				.r;
#line 433 "claqr5.f"
			z__4.r = z__5.r, z__4.i = z__5.i;
#line 433 "claqr5.f"
			i__6 = k + k * h_dim1;
#line 433 "claqr5.f"
			i__7 = k + 1 + (k + 1) * h_dim1;
#line 433 "claqr5.f"
			i__8 = k + 2 + (k + 2) * h_dim1;
#line 433 "claqr5.f"
			if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1)
				, abs(d__2)) + ((d__3 = z__4.r, abs(d__3)) + (
				d__4 = d_imag(&z__4), abs(d__4))) > ulp * ((
				d__5 = h__[i__6].r, abs(d__5)) + (d__6 = 
				d_imag(&h__[k + k * h_dim1]), abs(d__6)) + ((
				d__7 = h__[i__7].r, abs(d__7)) + (d__8 = 
				d_imag(&h__[k + 1 + (k + 1) * h_dim1]), abs(
				d__8))) + ((d__9 = h__[i__8].r, abs(d__9)) + (
				d__10 = d_imag(&h__[k + 2 + (k + 2) * h_dim1])
				, abs(d__10))))) {

/*                       ==== Starting a new bulge here would */
/*                       .    create non-negligible fill.  Use */
/*                       .    the old one with trepidation. ==== */

#line 442 "claqr5.f"
			    i__5 = k + 1 + k * h_dim1;
#line 442 "claqr5.f"
			    h__[i__5].r = beta.r, h__[i__5].i = beta.i;
#line 443 "claqr5.f"
			    i__5 = k + 2 + k * h_dim1;
#line 443 "claqr5.f"
			    h__[i__5].r = 0., h__[i__5].i = 0.;
#line 444 "claqr5.f"
			    i__5 = k + 3 + k * h_dim1;
#line 444 "claqr5.f"
			    h__[i__5].r = 0., h__[i__5].i = 0.;
#line 445 "claqr5.f"
			} else {

/*                       ==== Stating a new bulge here would */
/*                       .    create only negligible fill. */
/*                       .    Replace the old reflector with */
/*                       .    the new one. ==== */

#line 452 "claqr5.f"
			    i__5 = k + 1 + k * h_dim1;
#line 452 "claqr5.f"
			    i__6 = k + 1 + k * h_dim1;
#line 452 "claqr5.f"
			    z__1.r = h__[i__6].r - refsum.r, z__1.i = h__[
				    i__6].i - refsum.i;
#line 452 "claqr5.f"
			    h__[i__5].r = z__1.r, h__[i__5].i = z__1.i;
#line 453 "claqr5.f"
			    i__5 = k + 2 + k * h_dim1;
#line 453 "claqr5.f"
			    h__[i__5].r = 0., h__[i__5].i = 0.;
#line 454 "claqr5.f"
			    i__5 = k + 3 + k * h_dim1;
#line 454 "claqr5.f"
			    h__[i__5].r = 0., h__[i__5].i = 0.;
#line 455 "claqr5.f"
			    i__5 = m * v_dim1 + 1;
#line 455 "claqr5.f"
			    v[i__5].r = vt[0].r, v[i__5].i = vt[0].i;
#line 456 "claqr5.f"
			    i__5 = m * v_dim1 + 2;
#line 456 "claqr5.f"
			    v[i__5].r = vt[1].r, v[i__5].i = vt[1].i;
#line 457 "claqr5.f"
			    i__5 = m * v_dim1 + 3;
#line 457 "claqr5.f"
			    v[i__5].r = vt[2].r, v[i__5].i = vt[2].i;
#line 458 "claqr5.f"
			}
#line 459 "claqr5.f"
		    }
#line 460 "claqr5.f"
		}
#line 461 "claqr5.f"
/* L10: */
#line 461 "claqr5.f"
	    }

/*           ==== Generate a 2-by-2 reflection, if needed. ==== */

#line 465 "claqr5.f"
	    k = krcol + (m22 - 1) * 3;
#line 466 "claqr5.f"
	    if (bmp22) {
#line 467 "claqr5.f"
		if (k == *ktop - 1) {
#line 468 "claqr5.f"
		    claqr1_(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh, &s[(
			    m22 << 1) - 1], &s[m22 * 2], &v[m22 * v_dim1 + 1])
			    ;
#line 470 "claqr5.f"
		    i__4 = m22 * v_dim1 + 1;
#line 470 "claqr5.f"
		    beta.r = v[i__4].r, beta.i = v[i__4].i;
#line 471 "claqr5.f"
		    clarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 
			    * v_dim1 + 1]);
#line 472 "claqr5.f"
		} else {
#line 473 "claqr5.f"
		    i__4 = k + 1 + k * h_dim1;
#line 473 "claqr5.f"
		    beta.r = h__[i__4].r, beta.i = h__[i__4].i;
#line 474 "claqr5.f"
		    i__4 = m22 * v_dim1 + 2;
#line 474 "claqr5.f"
		    i__5 = k + 2 + k * h_dim1;
#line 474 "claqr5.f"
		    v[i__4].r = h__[i__5].r, v[i__4].i = h__[i__5].i;
#line 475 "claqr5.f"
		    clarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 
			    * v_dim1 + 1]);
#line 476 "claqr5.f"
		    i__4 = k + 1 + k * h_dim1;
#line 476 "claqr5.f"
		    h__[i__4].r = beta.r, h__[i__4].i = beta.i;
#line 477 "claqr5.f"
		    i__4 = k + 2 + k * h_dim1;
#line 477 "claqr5.f"
		    h__[i__4].r = 0., h__[i__4].i = 0.;
#line 478 "claqr5.f"
		}
#line 479 "claqr5.f"
	    }

/*           ==== Multiply H by reflections from the left ==== */

#line 483 "claqr5.f"
	    if (accum) {
#line 484 "claqr5.f"
		jbot = min(ndcol,*kbot);
#line 485 "claqr5.f"
	    } else if (*wantt) {
#line 486 "claqr5.f"
		jbot = *n;
#line 487 "claqr5.f"
	    } else {
#line 488 "claqr5.f"
		jbot = *kbot;
#line 489 "claqr5.f"
	    }
#line 490 "claqr5.f"
	    i__4 = jbot;
#line 490 "claqr5.f"
	    for (j = max(*ktop,krcol); j <= i__4; ++j) {
/* Computing MIN */
#line 491 "claqr5.f"
		i__5 = mbot, i__6 = (j - krcol + 2) / 3;
#line 491 "claqr5.f"
		mend = min(i__5,i__6);
#line 492 "claqr5.f"
		i__5 = mend;
#line 492 "claqr5.f"
		for (m = mtop; m <= i__5; ++m) {
#line 493 "claqr5.f"
		    k = krcol + (m - 1) * 3;
#line 494 "claqr5.f"
		    d_cnjg(&z__2, &v[m * v_dim1 + 1]);
#line 494 "claqr5.f"
		    i__6 = k + 1 + j * h_dim1;
#line 494 "claqr5.f"
		    d_cnjg(&z__6, &v[m * v_dim1 + 2]);
#line 494 "claqr5.f"
		    i__7 = k + 2 + j * h_dim1;
#line 494 "claqr5.f"
		    z__5.r = z__6.r * h__[i__7].r - z__6.i * h__[i__7].i, 
			    z__5.i = z__6.r * h__[i__7].i + z__6.i * h__[i__7]
			    .r;
#line 494 "claqr5.f"
		    z__4.r = h__[i__6].r + z__5.r, z__4.i = h__[i__6].i + 
			    z__5.i;
#line 494 "claqr5.f"
		    d_cnjg(&z__8, &v[m * v_dim1 + 3]);
#line 494 "claqr5.f"
		    i__8 = k + 3 + j * h_dim1;
#line 494 "claqr5.f"
		    z__7.r = z__8.r * h__[i__8].r - z__8.i * h__[i__8].i, 
			    z__7.i = z__8.r * h__[i__8].i + z__8.i * h__[i__8]
			    .r;
#line 494 "claqr5.f"
		    z__3.r = z__4.r + z__7.r, z__3.i = z__4.i + z__7.i;
#line 494 "claqr5.f"
		    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
			    z__2.r * z__3.i + z__2.i * z__3.r;
#line 494 "claqr5.f"
		    refsum.r = z__1.r, refsum.i = z__1.i;
#line 497 "claqr5.f"
		    i__6 = k + 1 + j * h_dim1;
#line 497 "claqr5.f"
		    i__7 = k + 1 + j * h_dim1;
#line 497 "claqr5.f"
		    z__1.r = h__[i__7].r - refsum.r, z__1.i = h__[i__7].i - 
			    refsum.i;
#line 497 "claqr5.f"
		    h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 498 "claqr5.f"
		    i__6 = k + 2 + j * h_dim1;
#line 498 "claqr5.f"
		    i__7 = k + 2 + j * h_dim1;
#line 498 "claqr5.f"
		    i__8 = m * v_dim1 + 2;
#line 498 "claqr5.f"
		    z__2.r = refsum.r * v[i__8].r - refsum.i * v[i__8].i, 
			    z__2.i = refsum.r * v[i__8].i + refsum.i * v[i__8]
			    .r;
#line 498 "claqr5.f"
		    z__1.r = h__[i__7].r - z__2.r, z__1.i = h__[i__7].i - 
			    z__2.i;
#line 498 "claqr5.f"
		    h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 499 "claqr5.f"
		    i__6 = k + 3 + j * h_dim1;
#line 499 "claqr5.f"
		    i__7 = k + 3 + j * h_dim1;
#line 499 "claqr5.f"
		    i__8 = m * v_dim1 + 3;
#line 499 "claqr5.f"
		    z__2.r = refsum.r * v[i__8].r - refsum.i * v[i__8].i, 
			    z__2.i = refsum.r * v[i__8].i + refsum.i * v[i__8]
			    .r;
#line 499 "claqr5.f"
		    z__1.r = h__[i__7].r - z__2.r, z__1.i = h__[i__7].i - 
			    z__2.i;
#line 499 "claqr5.f"
		    h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 500 "claqr5.f"
/* L20: */
#line 500 "claqr5.f"
		}
#line 501 "claqr5.f"
/* L30: */
#line 501 "claqr5.f"
	    }
#line 502 "claqr5.f"
	    if (bmp22) {
#line 503 "claqr5.f"
		k = krcol + (m22 - 1) * 3;
/* Computing MAX */
#line 504 "claqr5.f"
		i__4 = k + 1;
#line 504 "claqr5.f"
		i__5 = jbot;
#line 504 "claqr5.f"
		for (j = max(i__4,*ktop); j <= i__5; ++j) {
#line 505 "claqr5.f"
		    d_cnjg(&z__2, &v[m22 * v_dim1 + 1]);
#line 505 "claqr5.f"
		    i__4 = k + 1 + j * h_dim1;
#line 505 "claqr5.f"
		    d_cnjg(&z__5, &v[m22 * v_dim1 + 2]);
#line 505 "claqr5.f"
		    i__6 = k + 2 + j * h_dim1;
#line 505 "claqr5.f"
		    z__4.r = z__5.r * h__[i__6].r - z__5.i * h__[i__6].i, 
			    z__4.i = z__5.r * h__[i__6].i + z__5.i * h__[i__6]
			    .r;
#line 505 "claqr5.f"
		    z__3.r = h__[i__4].r + z__4.r, z__3.i = h__[i__4].i + 
			    z__4.i;
#line 505 "claqr5.f"
		    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
			    z__2.r * z__3.i + z__2.i * z__3.r;
#line 505 "claqr5.f"
		    refsum.r = z__1.r, refsum.i = z__1.i;
#line 508 "claqr5.f"
		    i__4 = k + 1 + j * h_dim1;
#line 508 "claqr5.f"
		    i__6 = k + 1 + j * h_dim1;
#line 508 "claqr5.f"
		    z__1.r = h__[i__6].r - refsum.r, z__1.i = h__[i__6].i - 
			    refsum.i;
#line 508 "claqr5.f"
		    h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 509 "claqr5.f"
		    i__4 = k + 2 + j * h_dim1;
#line 509 "claqr5.f"
		    i__6 = k + 2 + j * h_dim1;
#line 509 "claqr5.f"
		    i__7 = m22 * v_dim1 + 2;
#line 509 "claqr5.f"
		    z__2.r = refsum.r * v[i__7].r - refsum.i * v[i__7].i, 
			    z__2.i = refsum.r * v[i__7].i + refsum.i * v[i__7]
			    .r;
#line 509 "claqr5.f"
		    z__1.r = h__[i__6].r - z__2.r, z__1.i = h__[i__6].i - 
			    z__2.i;
#line 509 "claqr5.f"
		    h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 510 "claqr5.f"
/* L40: */
#line 510 "claqr5.f"
		}
#line 511 "claqr5.f"
	    }

/*           ==== Multiply H by reflections from the right. */
/*           .    Delay filling in the last row until the */
/*           .    vigilant deflation check is complete. ==== */

#line 517 "claqr5.f"
	    if (accum) {
#line 518 "claqr5.f"
		jtop = max(*ktop,incol);
#line 519 "claqr5.f"
	    } else if (*wantt) {
#line 520 "claqr5.f"
		jtop = 1;
#line 521 "claqr5.f"
	    } else {
#line 522 "claqr5.f"
		jtop = *ktop;
#line 523 "claqr5.f"
	    }
#line 524 "claqr5.f"
	    i__5 = mbot;
#line 524 "claqr5.f"
	    for (m = mtop; m <= i__5; ++m) {
#line 525 "claqr5.f"
		i__4 = m * v_dim1 + 1;
#line 525 "claqr5.f"
		if (v[i__4].r != 0. || v[i__4].i != 0.) {
#line 526 "claqr5.f"
		    k = krcol + (m - 1) * 3;
/* Computing MIN */
#line 527 "claqr5.f"
		    i__6 = *kbot, i__7 = k + 3;
#line 527 "claqr5.f"
		    i__4 = min(i__6,i__7);
#line 527 "claqr5.f"
		    for (j = jtop; j <= i__4; ++j) {
#line 528 "claqr5.f"
			i__6 = m * v_dim1 + 1;
#line 528 "claqr5.f"
			i__7 = j + (k + 1) * h_dim1;
#line 528 "claqr5.f"
			i__8 = m * v_dim1 + 2;
#line 528 "claqr5.f"
			i__9 = j + (k + 2) * h_dim1;
#line 528 "claqr5.f"
			z__4.r = v[i__8].r * h__[i__9].r - v[i__8].i * h__[
				i__9].i, z__4.i = v[i__8].r * h__[i__9].i + v[
				i__8].i * h__[i__9].r;
#line 528 "claqr5.f"
			z__3.r = h__[i__7].r + z__4.r, z__3.i = h__[i__7].i + 
				z__4.i;
#line 528 "claqr5.f"
			i__10 = m * v_dim1 + 3;
#line 528 "claqr5.f"
			i__11 = j + (k + 3) * h_dim1;
#line 528 "claqr5.f"
			z__5.r = v[i__10].r * h__[i__11].r - v[i__10].i * h__[
				i__11].i, z__5.i = v[i__10].r * h__[i__11].i 
				+ v[i__10].i * h__[i__11].r;
#line 528 "claqr5.f"
			z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
#line 528 "claqr5.f"
			z__1.r = v[i__6].r * z__2.r - v[i__6].i * z__2.i, 
				z__1.i = v[i__6].r * z__2.i + v[i__6].i * 
				z__2.r;
#line 528 "claqr5.f"
			refsum.r = z__1.r, refsum.i = z__1.i;
#line 530 "claqr5.f"
			i__6 = j + (k + 1) * h_dim1;
#line 530 "claqr5.f"
			i__7 = j + (k + 1) * h_dim1;
#line 530 "claqr5.f"
			z__1.r = h__[i__7].r - refsum.r, z__1.i = h__[i__7].i 
				- refsum.i;
#line 530 "claqr5.f"
			h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 531 "claqr5.f"
			i__6 = j + (k + 2) * h_dim1;
#line 531 "claqr5.f"
			i__7 = j + (k + 2) * h_dim1;
#line 531 "claqr5.f"
			d_cnjg(&z__3, &v[m * v_dim1 + 2]);
#line 531 "claqr5.f"
			z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				z__2.i = refsum.r * z__3.i + refsum.i * 
				z__3.r;
#line 531 "claqr5.f"
			z__1.r = h__[i__7].r - z__2.r, z__1.i = h__[i__7].i - 
				z__2.i;
#line 531 "claqr5.f"
			h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 533 "claqr5.f"
			i__6 = j + (k + 3) * h_dim1;
#line 533 "claqr5.f"
			i__7 = j + (k + 3) * h_dim1;
#line 533 "claqr5.f"
			d_cnjg(&z__3, &v[m * v_dim1 + 3]);
#line 533 "claqr5.f"
			z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				z__2.i = refsum.r * z__3.i + refsum.i * 
				z__3.r;
#line 533 "claqr5.f"
			z__1.r = h__[i__7].r - z__2.r, z__1.i = h__[i__7].i - 
				z__2.i;
#line 533 "claqr5.f"
			h__[i__6].r = z__1.r, h__[i__6].i = z__1.i;
#line 535 "claqr5.f"
/* L50: */
#line 535 "claqr5.f"
		    }

#line 537 "claqr5.f"
		    if (accum) {

/*                    ==== Accumulate U. (If necessary, update Z later */
/*                    .    with with an efficient matrix-matrix */
/*                    .    multiply.) ==== */

#line 543 "claqr5.f"
			kms = k - incol;
/* Computing MAX */
#line 544 "claqr5.f"
			i__4 = 1, i__6 = *ktop - incol;
#line 544 "claqr5.f"
			i__7 = kdu;
#line 544 "claqr5.f"
			for (j = max(i__4,i__6); j <= i__7; ++j) {
#line 545 "claqr5.f"
			    i__4 = m * v_dim1 + 1;
#line 545 "claqr5.f"
			    i__6 = j + (kms + 1) * u_dim1;
#line 545 "claqr5.f"
			    i__8 = m * v_dim1 + 2;
#line 545 "claqr5.f"
			    i__9 = j + (kms + 2) * u_dim1;
#line 545 "claqr5.f"
			    z__4.r = v[i__8].r * u[i__9].r - v[i__8].i * u[
				    i__9].i, z__4.i = v[i__8].r * u[i__9].i + 
				    v[i__8].i * u[i__9].r;
#line 545 "claqr5.f"
			    z__3.r = u[i__6].r + z__4.r, z__3.i = u[i__6].i + 
				    z__4.i;
#line 545 "claqr5.f"
			    i__10 = m * v_dim1 + 3;
#line 545 "claqr5.f"
			    i__11 = j + (kms + 3) * u_dim1;
#line 545 "claqr5.f"
			    z__5.r = v[i__10].r * u[i__11].r - v[i__10].i * u[
				    i__11].i, z__5.i = v[i__10].r * u[i__11]
				    .i + v[i__10].i * u[i__11].r;
#line 545 "claqr5.f"
			    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + 
				    z__5.i;
#line 545 "claqr5.f"
			    z__1.r = v[i__4].r * z__2.r - v[i__4].i * z__2.i, 
				    z__1.i = v[i__4].r * z__2.i + v[i__4].i * 
				    z__2.r;
#line 545 "claqr5.f"
			    refsum.r = z__1.r, refsum.i = z__1.i;
#line 547 "claqr5.f"
			    i__4 = j + (kms + 1) * u_dim1;
#line 547 "claqr5.f"
			    i__6 = j + (kms + 1) * u_dim1;
#line 547 "claqr5.f"
			    z__1.r = u[i__6].r - refsum.r, z__1.i = u[i__6].i 
				    - refsum.i;
#line 547 "claqr5.f"
			    u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 548 "claqr5.f"
			    i__4 = j + (kms + 2) * u_dim1;
#line 548 "claqr5.f"
			    i__6 = j + (kms + 2) * u_dim1;
#line 548 "claqr5.f"
			    d_cnjg(&z__3, &v[m * v_dim1 + 2]);
#line 548 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 548 "claqr5.f"
			    z__1.r = u[i__6].r - z__2.r, z__1.i = u[i__6].i - 
				    z__2.i;
#line 548 "claqr5.f"
			    u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 550 "claqr5.f"
			    i__4 = j + (kms + 3) * u_dim1;
#line 550 "claqr5.f"
			    i__6 = j + (kms + 3) * u_dim1;
#line 550 "claqr5.f"
			    d_cnjg(&z__3, &v[m * v_dim1 + 3]);
#line 550 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 550 "claqr5.f"
			    z__1.r = u[i__6].r - z__2.r, z__1.i = u[i__6].i - 
				    z__2.i;
#line 550 "claqr5.f"
			    u[i__4].r = z__1.r, u[i__4].i = z__1.i;
#line 552 "claqr5.f"
/* L60: */
#line 552 "claqr5.f"
			}
#line 553 "claqr5.f"
		    } else if (*wantz) {

/*                    ==== U is not accumulated, so update Z */
/*                    .    now by multiplying by reflections */
/*                    .    from the right. ==== */

#line 559 "claqr5.f"
			i__7 = *ihiz;
#line 559 "claqr5.f"
			for (j = *iloz; j <= i__7; ++j) {
#line 560 "claqr5.f"
			    i__4 = m * v_dim1 + 1;
#line 560 "claqr5.f"
			    i__6 = j + (k + 1) * z_dim1;
#line 560 "claqr5.f"
			    i__8 = m * v_dim1 + 2;
#line 560 "claqr5.f"
			    i__9 = j + (k + 2) * z_dim1;
#line 560 "claqr5.f"
			    z__4.r = v[i__8].r * z__[i__9].r - v[i__8].i * 
				    z__[i__9].i, z__4.i = v[i__8].r * z__[
				    i__9].i + v[i__8].i * z__[i__9].r;
#line 560 "claqr5.f"
			    z__3.r = z__[i__6].r + z__4.r, z__3.i = z__[i__6]
				    .i + z__4.i;
#line 560 "claqr5.f"
			    i__10 = m * v_dim1 + 3;
#line 560 "claqr5.f"
			    i__11 = j + (k + 3) * z_dim1;
#line 560 "claqr5.f"
			    z__5.r = v[i__10].r * z__[i__11].r - v[i__10].i * 
				    z__[i__11].i, z__5.i = v[i__10].r * z__[
				    i__11].i + v[i__10].i * z__[i__11].r;
#line 560 "claqr5.f"
			    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + 
				    z__5.i;
#line 560 "claqr5.f"
			    z__1.r = v[i__4].r * z__2.r - v[i__4].i * z__2.i, 
				    z__1.i = v[i__4].r * z__2.i + v[i__4].i * 
				    z__2.r;
#line 560 "claqr5.f"
			    refsum.r = z__1.r, refsum.i = z__1.i;
#line 562 "claqr5.f"
			    i__4 = j + (k + 1) * z_dim1;
#line 562 "claqr5.f"
			    i__6 = j + (k + 1) * z_dim1;
#line 562 "claqr5.f"
			    z__1.r = z__[i__6].r - refsum.r, z__1.i = z__[
				    i__6].i - refsum.i;
#line 562 "claqr5.f"
			    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 563 "claqr5.f"
			    i__4 = j + (k + 2) * z_dim1;
#line 563 "claqr5.f"
			    i__6 = j + (k + 2) * z_dim1;
#line 563 "claqr5.f"
			    d_cnjg(&z__3, &v[m * v_dim1 + 2]);
#line 563 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 563 "claqr5.f"
			    z__1.r = z__[i__6].r - z__2.r, z__1.i = z__[i__6]
				    .i - z__2.i;
#line 563 "claqr5.f"
			    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 565 "claqr5.f"
			    i__4 = j + (k + 3) * z_dim1;
#line 565 "claqr5.f"
			    i__6 = j + (k + 3) * z_dim1;
#line 565 "claqr5.f"
			    d_cnjg(&z__3, &v[m * v_dim1 + 3]);
#line 565 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 565 "claqr5.f"
			    z__1.r = z__[i__6].r - z__2.r, z__1.i = z__[i__6]
				    .i - z__2.i;
#line 565 "claqr5.f"
			    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 567 "claqr5.f"
/* L70: */
#line 567 "claqr5.f"
			}
#line 568 "claqr5.f"
		    }
#line 569 "claqr5.f"
		}
#line 570 "claqr5.f"
/* L80: */
#line 570 "claqr5.f"
	    }

/*           ==== Special case: 2-by-2 reflection (if needed) ==== */

#line 574 "claqr5.f"
	    k = krcol + (m22 - 1) * 3;
#line 575 "claqr5.f"
	    if (bmp22) {
#line 576 "claqr5.f"
		i__5 = m22 * v_dim1 + 1;
#line 576 "claqr5.f"
		if (v[i__5].r != 0. || v[i__5].i != 0.) {
/* Computing MIN */
#line 577 "claqr5.f"
		    i__7 = *kbot, i__4 = k + 3;
#line 577 "claqr5.f"
		    i__5 = min(i__7,i__4);
#line 577 "claqr5.f"
		    for (j = jtop; j <= i__5; ++j) {
#line 578 "claqr5.f"
			i__7 = m22 * v_dim1 + 1;
#line 578 "claqr5.f"
			i__4 = j + (k + 1) * h_dim1;
#line 578 "claqr5.f"
			i__6 = m22 * v_dim1 + 2;
#line 578 "claqr5.f"
			i__8 = j + (k + 2) * h_dim1;
#line 578 "claqr5.f"
			z__3.r = v[i__6].r * h__[i__8].r - v[i__6].i * h__[
				i__8].i, z__3.i = v[i__6].r * h__[i__8].i + v[
				i__6].i * h__[i__8].r;
#line 578 "claqr5.f"
			z__2.r = h__[i__4].r + z__3.r, z__2.i = h__[i__4].i + 
				z__3.i;
#line 578 "claqr5.f"
			z__1.r = v[i__7].r * z__2.r - v[i__7].i * z__2.i, 
				z__1.i = v[i__7].r * z__2.i + v[i__7].i * 
				z__2.r;
#line 578 "claqr5.f"
			refsum.r = z__1.r, refsum.i = z__1.i;
#line 580 "claqr5.f"
			i__7 = j + (k + 1) * h_dim1;
#line 580 "claqr5.f"
			i__4 = j + (k + 1) * h_dim1;
#line 580 "claqr5.f"
			z__1.r = h__[i__4].r - refsum.r, z__1.i = h__[i__4].i 
				- refsum.i;
#line 580 "claqr5.f"
			h__[i__7].r = z__1.r, h__[i__7].i = z__1.i;
#line 581 "claqr5.f"
			i__7 = j + (k + 2) * h_dim1;
#line 581 "claqr5.f"
			i__4 = j + (k + 2) * h_dim1;
#line 581 "claqr5.f"
			d_cnjg(&z__3, &v[m22 * v_dim1 + 2]);
#line 581 "claqr5.f"
			z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				z__2.i = refsum.r * z__3.i + refsum.i * 
				z__3.r;
#line 581 "claqr5.f"
			z__1.r = h__[i__4].r - z__2.r, z__1.i = h__[i__4].i - 
				z__2.i;
#line 581 "claqr5.f"
			h__[i__7].r = z__1.r, h__[i__7].i = z__1.i;
#line 583 "claqr5.f"
/* L90: */
#line 583 "claqr5.f"
		    }

#line 585 "claqr5.f"
		    if (accum) {
#line 586 "claqr5.f"
			kms = k - incol;
/* Computing MAX */
#line 587 "claqr5.f"
			i__5 = 1, i__7 = *ktop - incol;
#line 587 "claqr5.f"
			i__4 = kdu;
#line 587 "claqr5.f"
			for (j = max(i__5,i__7); j <= i__4; ++j) {
#line 588 "claqr5.f"
			    i__5 = m22 * v_dim1 + 1;
#line 588 "claqr5.f"
			    i__7 = j + (kms + 1) * u_dim1;
#line 588 "claqr5.f"
			    i__6 = m22 * v_dim1 + 2;
#line 588 "claqr5.f"
			    i__8 = j + (kms + 2) * u_dim1;
#line 588 "claqr5.f"
			    z__3.r = v[i__6].r * u[i__8].r - v[i__6].i * u[
				    i__8].i, z__3.i = v[i__6].r * u[i__8].i + 
				    v[i__6].i * u[i__8].r;
#line 588 "claqr5.f"
			    z__2.r = u[i__7].r + z__3.r, z__2.i = u[i__7].i + 
				    z__3.i;
#line 588 "claqr5.f"
			    z__1.r = v[i__5].r * z__2.r - v[i__5].i * z__2.i, 
				    z__1.i = v[i__5].r * z__2.i + v[i__5].i * 
				    z__2.r;
#line 588 "claqr5.f"
			    refsum.r = z__1.r, refsum.i = z__1.i;
#line 590 "claqr5.f"
			    i__5 = j + (kms + 1) * u_dim1;
#line 590 "claqr5.f"
			    i__7 = j + (kms + 1) * u_dim1;
#line 590 "claqr5.f"
			    z__1.r = u[i__7].r - refsum.r, z__1.i = u[i__7].i 
				    - refsum.i;
#line 590 "claqr5.f"
			    u[i__5].r = z__1.r, u[i__5].i = z__1.i;
#line 591 "claqr5.f"
			    i__5 = j + (kms + 2) * u_dim1;
#line 591 "claqr5.f"
			    i__7 = j + (kms + 2) * u_dim1;
#line 591 "claqr5.f"
			    d_cnjg(&z__3, &v[m22 * v_dim1 + 2]);
#line 591 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 591 "claqr5.f"
			    z__1.r = u[i__7].r - z__2.r, z__1.i = u[i__7].i - 
				    z__2.i;
#line 591 "claqr5.f"
			    u[i__5].r = z__1.r, u[i__5].i = z__1.i;
#line 593 "claqr5.f"
/* L100: */
#line 593 "claqr5.f"
			}
#line 594 "claqr5.f"
		    } else if (*wantz) {
#line 595 "claqr5.f"
			i__4 = *ihiz;
#line 595 "claqr5.f"
			for (j = *iloz; j <= i__4; ++j) {
#line 596 "claqr5.f"
			    i__5 = m22 * v_dim1 + 1;
#line 596 "claqr5.f"
			    i__7 = j + (k + 1) * z_dim1;
#line 596 "claqr5.f"
			    i__6 = m22 * v_dim1 + 2;
#line 596 "claqr5.f"
			    i__8 = j + (k + 2) * z_dim1;
#line 596 "claqr5.f"
			    z__3.r = v[i__6].r * z__[i__8].r - v[i__6].i * 
				    z__[i__8].i, z__3.i = v[i__6].r * z__[
				    i__8].i + v[i__6].i * z__[i__8].r;
#line 596 "claqr5.f"
			    z__2.r = z__[i__7].r + z__3.r, z__2.i = z__[i__7]
				    .i + z__3.i;
#line 596 "claqr5.f"
			    z__1.r = v[i__5].r * z__2.r - v[i__5].i * z__2.i, 
				    z__1.i = v[i__5].r * z__2.i + v[i__5].i * 
				    z__2.r;
#line 596 "claqr5.f"
			    refsum.r = z__1.r, refsum.i = z__1.i;
#line 598 "claqr5.f"
			    i__5 = j + (k + 1) * z_dim1;
#line 598 "claqr5.f"
			    i__7 = j + (k + 1) * z_dim1;
#line 598 "claqr5.f"
			    z__1.r = z__[i__7].r - refsum.r, z__1.i = z__[
				    i__7].i - refsum.i;
#line 598 "claqr5.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 599 "claqr5.f"
			    i__5 = j + (k + 2) * z_dim1;
#line 599 "claqr5.f"
			    i__7 = j + (k + 2) * z_dim1;
#line 599 "claqr5.f"
			    d_cnjg(&z__3, &v[m22 * v_dim1 + 2]);
#line 599 "claqr5.f"
			    z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, 
				    z__2.i = refsum.r * z__3.i + refsum.i * 
				    z__3.r;
#line 599 "claqr5.f"
			    z__1.r = z__[i__7].r - z__2.r, z__1.i = z__[i__7]
				    .i - z__2.i;
#line 599 "claqr5.f"
			    z__[i__5].r = z__1.r, z__[i__5].i = z__1.i;
#line 601 "claqr5.f"
/* L110: */
#line 601 "claqr5.f"
			}
#line 602 "claqr5.f"
		    }
#line 603 "claqr5.f"
		}
#line 604 "claqr5.f"
	    }

/*           ==== Vigilant deflation check ==== */

#line 608 "claqr5.f"
	    mstart = mtop;
#line 609 "claqr5.f"
	    if (krcol + (mstart - 1) * 3 < *ktop) {
#line 609 "claqr5.f"
		++mstart;
#line 609 "claqr5.f"
	    }
#line 611 "claqr5.f"
	    mend = mbot;
#line 612 "claqr5.f"
	    if (bmp22) {
#line 612 "claqr5.f"
		++mend;
#line 612 "claqr5.f"
	    }
#line 614 "claqr5.f"
	    if (krcol == *kbot - 2) {
#line 614 "claqr5.f"
		++mend;
#line 614 "claqr5.f"
	    }
#line 616 "claqr5.f"
	    i__4 = mend;
#line 616 "claqr5.f"
	    for (m = mstart; m <= i__4; ++m) {
/* Computing MIN */
#line 617 "claqr5.f"
		i__5 = *kbot - 1, i__7 = krcol + (m - 1) * 3;
#line 617 "claqr5.f"
		k = min(i__5,i__7);

/*              ==== The following convergence test requires that */
/*              .    the tradition small-compared-to-nearby-diagonals */
/*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997) */
/*              .    criteria both be satisfied.  The latter improves */
/*              .    accuracy in some examples. Falling back on an */
/*              .    alternate convergence criterion when TST1 or TST2 */
/*              .    is zero (as done here) is traditional but probably */
/*              .    unnecessary. ==== */

#line 628 "claqr5.f"
		i__5 = k + 1 + k * h_dim1;
#line 628 "claqr5.f"
		if (h__[i__5].r != 0. || h__[i__5].i != 0.) {
#line 629 "claqr5.f"
		    i__5 = k + k * h_dim1;
#line 629 "claqr5.f"
		    i__7 = k + 1 + (k + 1) * h_dim1;
#line 629 "claqr5.f"
		    tst1 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&
			    h__[k + k * h_dim1]), abs(d__2)) + ((d__3 = h__[
			    i__7].r, abs(d__3)) + (d__4 = d_imag(&h__[k + 1 + 
			    (k + 1) * h_dim1]), abs(d__4)));
#line 630 "claqr5.f"
		    if (tst1 == 0.) {
#line 631 "claqr5.f"
			if (k >= *ktop + 1) {
#line 631 "claqr5.f"
			    i__5 = k + (k - 1) * h_dim1;
#line 631 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + (k - 1) * h_dim1]), abs(
				    d__2));
#line 631 "claqr5.f"
			}
#line 633 "claqr5.f"
			if (k >= *ktop + 2) {
#line 633 "claqr5.f"
			    i__5 = k + (k - 2) * h_dim1;
#line 633 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + (k - 2) * h_dim1]), abs(
				    d__2));
#line 633 "claqr5.f"
			}
#line 635 "claqr5.f"
			if (k >= *ktop + 3) {
#line 635 "claqr5.f"
			    i__5 = k + (k - 3) * h_dim1;
#line 635 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + (k - 3) * h_dim1]), abs(
				    d__2));
#line 635 "claqr5.f"
			}
#line 637 "claqr5.f"
			if (k <= *kbot - 2) {
#line 637 "claqr5.f"
			    i__5 = k + 2 + (k + 1) * h_dim1;
#line 637 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + 2 + (k + 1) * h_dim1]), 
				    abs(d__2));
#line 637 "claqr5.f"
			}
#line 639 "claqr5.f"
			if (k <= *kbot - 3) {
#line 639 "claqr5.f"
			    i__5 = k + 3 + (k + 1) * h_dim1;
#line 639 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + 3 + (k + 1) * h_dim1]), 
				    abs(d__2));
#line 639 "claqr5.f"
			}
#line 641 "claqr5.f"
			if (k <= *kbot - 4) {
#line 641 "claqr5.f"
			    i__5 = k + 4 + (k + 1) * h_dim1;
#line 641 "claqr5.f"
			    tst1 += (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				    d_imag(&h__[k + 4 + (k + 1) * h_dim1]), 
				    abs(d__2));
#line 641 "claqr5.f"
			}
#line 643 "claqr5.f"
		    }
#line 644 "claqr5.f"
		    i__5 = k + 1 + k * h_dim1;
/* Computing MAX */
#line 644 "claqr5.f"
		    d__3 = smlnum, d__4 = ulp * tst1;
#line 644 "claqr5.f"
		    if ((d__1 = h__[i__5].r, abs(d__1)) + (d__2 = d_imag(&h__[
			    k + 1 + k * h_dim1]), abs(d__2)) <= max(d__3,d__4)
			    ) {
/* Computing MAX */
#line 646 "claqr5.f"
			i__5 = k + 1 + k * h_dim1;
#line 646 "claqr5.f"
			i__7 = k + (k + 1) * h_dim1;
#line 646 "claqr5.f"
			d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				d_imag(&h__[k + 1 + k * h_dim1]), abs(d__2)), 
				d__6 = (d__3 = h__[i__7].r, abs(d__3)) + (
				d__4 = d_imag(&h__[k + (k + 1) * h_dim1]), 
				abs(d__4));
#line 646 "claqr5.f"
			h12 = max(d__5,d__6);
/* Computing MIN */
#line 648 "claqr5.f"
			i__5 = k + 1 + k * h_dim1;
#line 648 "claqr5.f"
			i__7 = k + (k + 1) * h_dim1;
#line 648 "claqr5.f"
			d__5 = (d__1 = h__[i__5].r, abs(d__1)) + (d__2 = 
				d_imag(&h__[k + 1 + k * h_dim1]), abs(d__2)), 
				d__6 = (d__3 = h__[i__7].r, abs(d__3)) + (
				d__4 = d_imag(&h__[k + (k + 1) * h_dim1]), 
				abs(d__4));
#line 648 "claqr5.f"
			h21 = min(d__5,d__6);
#line 650 "claqr5.f"
			i__5 = k + k * h_dim1;
#line 650 "claqr5.f"
			i__7 = k + 1 + (k + 1) * h_dim1;
#line 650 "claqr5.f"
			z__2.r = h__[i__5].r - h__[i__7].r, z__2.i = h__[i__5]
				.i - h__[i__7].i;
#line 650 "claqr5.f"
			z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MAX */
#line 650 "claqr5.f"
			i__6 = k + 1 + (k + 1) * h_dim1;
#line 650 "claqr5.f"
			d__5 = (d__1 = h__[i__6].r, abs(d__1)) + (d__2 = 
				d_imag(&h__[k + 1 + (k + 1) * h_dim1]), abs(
				d__2)), d__6 = (d__3 = z__1.r, abs(d__3)) + (
				d__4 = d_imag(&z__1), abs(d__4));
#line 650 "claqr5.f"
			h11 = max(d__5,d__6);
#line 652 "claqr5.f"
			i__5 = k + k * h_dim1;
#line 652 "claqr5.f"
			i__7 = k + 1 + (k + 1) * h_dim1;
#line 652 "claqr5.f"
			z__2.r = h__[i__5].r - h__[i__7].r, z__2.i = h__[i__5]
				.i - h__[i__7].i;
#line 652 "claqr5.f"
			z__1.r = z__2.r, z__1.i = z__2.i;
/* Computing MIN */
#line 652 "claqr5.f"
			i__6 = k + 1 + (k + 1) * h_dim1;
#line 652 "claqr5.f"
			d__5 = (d__1 = h__[i__6].r, abs(d__1)) + (d__2 = 
				d_imag(&h__[k + 1 + (k + 1) * h_dim1]), abs(
				d__2)), d__6 = (d__3 = z__1.r, abs(d__3)) + (
				d__4 = d_imag(&z__1), abs(d__4));
#line 652 "claqr5.f"
			h22 = min(d__5,d__6);
#line 654 "claqr5.f"
			scl = h11 + h12;
#line 655 "claqr5.f"
			tst2 = h22 * (h11 / scl);

/* Computing MAX */
#line 657 "claqr5.f"
			d__1 = smlnum, d__2 = ulp * tst2;
#line 657 "claqr5.f"
			if (tst2 == 0. || h21 * (h12 / scl) <= max(d__1,d__2))
				 {
#line 657 "claqr5.f"
			    i__5 = k + 1 + k * h_dim1;
#line 657 "claqr5.f"
			    h__[i__5].r = 0., h__[i__5].i = 0.;
#line 657 "claqr5.f"
			}
#line 659 "claqr5.f"
		    }
#line 660 "claqr5.f"
		}
#line 661 "claqr5.f"
/* L120: */
#line 661 "claqr5.f"
	    }

/*           ==== Fill in the last row of each bulge. ==== */

/* Computing MIN */
#line 665 "claqr5.f"
	    i__4 = nbmps, i__5 = (*kbot - krcol - 1) / 3;
#line 665 "claqr5.f"
	    mend = min(i__4,i__5);
#line 666 "claqr5.f"
	    i__4 = mend;
#line 666 "claqr5.f"
	    for (m = mtop; m <= i__4; ++m) {
#line 667 "claqr5.f"
		k = krcol + (m - 1) * 3;
#line 668 "claqr5.f"
		i__5 = m * v_dim1 + 1;
#line 668 "claqr5.f"
		i__7 = m * v_dim1 + 3;
#line 668 "claqr5.f"
		z__2.r = v[i__5].r * v[i__7].r - v[i__5].i * v[i__7].i, 
			z__2.i = v[i__5].r * v[i__7].i + v[i__5].i * v[i__7]
			.r;
#line 668 "claqr5.f"
		i__6 = k + 4 + (k + 3) * h_dim1;
#line 668 "claqr5.f"
		z__1.r = z__2.r * h__[i__6].r - z__2.i * h__[i__6].i, z__1.i =
			 z__2.r * h__[i__6].i + z__2.i * h__[i__6].r;
#line 668 "claqr5.f"
		refsum.r = z__1.r, refsum.i = z__1.i;
#line 669 "claqr5.f"
		i__5 = k + 4 + (k + 1) * h_dim1;
#line 669 "claqr5.f"
		z__1.r = -refsum.r, z__1.i = -refsum.i;
#line 669 "claqr5.f"
		h__[i__5].r = z__1.r, h__[i__5].i = z__1.i;
#line 670 "claqr5.f"
		i__5 = k + 4 + (k + 2) * h_dim1;
#line 670 "claqr5.f"
		z__2.r = -refsum.r, z__2.i = -refsum.i;
#line 670 "claqr5.f"
		d_cnjg(&z__3, &v[m * v_dim1 + 2]);
#line 670 "claqr5.f"
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
#line 670 "claqr5.f"
		h__[i__5].r = z__1.r, h__[i__5].i = z__1.i;
#line 671 "claqr5.f"
		i__5 = k + 4 + (k + 3) * h_dim1;
#line 671 "claqr5.f"
		i__7 = k + 4 + (k + 3) * h_dim1;
#line 671 "claqr5.f"
		d_cnjg(&z__3, &v[m * v_dim1 + 3]);
#line 671 "claqr5.f"
		z__2.r = refsum.r * z__3.r - refsum.i * z__3.i, z__2.i = 
			refsum.r * z__3.i + refsum.i * z__3.r;
#line 671 "claqr5.f"
		z__1.r = h__[i__7].r - z__2.r, z__1.i = h__[i__7].i - z__2.i;
#line 671 "claqr5.f"
		h__[i__5].r = z__1.r, h__[i__5].i = z__1.i;
#line 672 "claqr5.f"
/* L130: */
#line 672 "claqr5.f"
	    }

/*           ==== End of near-the-diagonal bulge chase. ==== */

#line 676 "claqr5.f"
/* L140: */
#line 676 "claqr5.f"
	}

/*        ==== Use U (if accumulated) to update far-from-diagonal */
/*        .    entries in H.  If required, use U to update Z as */
/*        .    well. ==== */

#line 682 "claqr5.f"
	if (accum) {
#line 683 "claqr5.f"
	    if (*wantt) {
#line 684 "claqr5.f"
		jtop = 1;
#line 685 "claqr5.f"
		jbot = *n;
#line 686 "claqr5.f"
	    } else {
#line 687 "claqr5.f"
		jtop = *ktop;
#line 688 "claqr5.f"
		jbot = *kbot;
#line 689 "claqr5.f"
	    }
#line 690 "claqr5.f"
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
#line 702 "claqr5.f"
		i__3 = 1, i__4 = *ktop - incol;
#line 702 "claqr5.f"
		k1 = max(i__3,i__4);
/* Computing MAX */
#line 703 "claqr5.f"
		i__3 = 0, i__4 = ndcol - *kbot;
#line 703 "claqr5.f"
		nu = kdu - max(i__3,i__4) - k1 + 1;

/*              ==== Horizontal Multiply ==== */

#line 707 "claqr5.f"
		i__3 = jbot;
#line 707 "claqr5.f"
		i__4 = *nh;
#line 707 "claqr5.f"
		for (jcol = min(ndcol,*kbot) + 1; i__4 < 0 ? jcol >= i__3 : 
			jcol <= i__3; jcol += i__4) {
/* Computing MIN */
#line 708 "claqr5.f"
		    i__5 = *nh, i__7 = jbot - jcol + 1;
#line 708 "claqr5.f"
		    jlen = min(i__5,i__7);
#line 709 "claqr5.f"
		    cgemm_("C", "N", &nu, &jlen, &nu, &c_b2, &u[k1 + k1 * 
			    u_dim1], ldu, &h__[incol + k1 + jcol * h_dim1], 
			    ldh, &c_b1, &wh[wh_offset], ldwh, (ftnlen)1, (
			    ftnlen)1);
#line 712 "claqr5.f"
		    clacpy_("ALL", &nu, &jlen, &wh[wh_offset], ldwh, &h__[
			    incol + k1 + jcol * h_dim1], ldh, (ftnlen)3);
#line 714 "claqr5.f"
/* L150: */
#line 714 "claqr5.f"
		}

/*              ==== Vertical multiply ==== */

#line 718 "claqr5.f"
		i__4 = max(*ktop,incol) - 1;
#line 718 "claqr5.f"
		i__3 = *nv;
#line 718 "claqr5.f"
		for (jrow = jtop; i__3 < 0 ? jrow >= i__4 : jrow <= i__4; 
			jrow += i__3) {
/* Computing MIN */
#line 719 "claqr5.f"
		    i__5 = *nv, i__7 = max(*ktop,incol) - jrow;
#line 719 "claqr5.f"
		    jlen = min(i__5,i__7);
#line 720 "claqr5.f"
		    cgemm_("N", "N", &jlen, &nu, &nu, &c_b2, &h__[jrow + (
			    incol + k1) * h_dim1], ldh, &u[k1 + k1 * u_dim1], 
			    ldu, &c_b1, &wv[wv_offset], ldwv, (ftnlen)1, (
			    ftnlen)1);
#line 723 "claqr5.f"
		    clacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &h__[
			    jrow + (incol + k1) * h_dim1], ldh, (ftnlen)3);
#line 725 "claqr5.f"
/* L160: */
#line 725 "claqr5.f"
		}

/*              ==== Z multiply (also vertical) ==== */

#line 729 "claqr5.f"
		if (*wantz) {
#line 730 "claqr5.f"
		    i__3 = *ihiz;
#line 730 "claqr5.f"
		    i__4 = *nv;
#line 730 "claqr5.f"
		    for (jrow = *iloz; i__4 < 0 ? jrow >= i__3 : jrow <= i__3;
			     jrow += i__4) {
/* Computing MIN */
#line 731 "claqr5.f"
			i__5 = *nv, i__7 = *ihiz - jrow + 1;
#line 731 "claqr5.f"
			jlen = min(i__5,i__7);
#line 732 "claqr5.f"
			cgemm_("N", "N", &jlen, &nu, &nu, &c_b2, &z__[jrow + (
				incol + k1) * z_dim1], ldz, &u[k1 + k1 * 
				u_dim1], ldu, &c_b1, &wv[wv_offset], ldwv, (
				ftnlen)1, (ftnlen)1);
#line 735 "claqr5.f"
			clacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &z__[
				jrow + (incol + k1) * z_dim1], ldz, (ftnlen)3)
				;
#line 737 "claqr5.f"
/* L170: */
#line 737 "claqr5.f"
		    }
#line 738 "claqr5.f"
		}
#line 739 "claqr5.f"
	    } else {

/*              ==== Updates exploiting U's 2-by-2 block structure. */
/*              .    (I2, I4, J2, J4 are the last rows and columns */
/*              .    of the blocks.) ==== */

#line 745 "claqr5.f"
		i2 = (kdu + 1) / 2;
#line 746 "claqr5.f"
		i4 = kdu;
#line 747 "claqr5.f"
		j2 = i4 - i2;
#line 748 "claqr5.f"
		j4 = kdu;

/*              ==== KZS and KNZ deal with the band of zeros */
/*              .    along the diagonal of one of the triangular */
/*              .    blocks. ==== */

#line 754 "claqr5.f"
		kzs = j4 - j2 - (ns + 1);
#line 755 "claqr5.f"
		knz = ns + 1;

/*              ==== Horizontal multiply ==== */

#line 759 "claqr5.f"
		i__4 = jbot;
#line 759 "claqr5.f"
		i__3 = *nh;
#line 759 "claqr5.f"
		for (jcol = min(ndcol,*kbot) + 1; i__3 < 0 ? jcol >= i__4 : 
			jcol <= i__4; jcol += i__3) {
/* Computing MIN */
#line 760 "claqr5.f"
		    i__5 = *nh, i__7 = jbot - jcol + 1;
#line 760 "claqr5.f"
		    jlen = min(i__5,i__7);

/*                 ==== Copy bottom of H to top+KZS of scratch ==== */
/*                  (The first KZS rows get multiplied by zero.) ==== */

#line 765 "claqr5.f"
		    clacpy_("ALL", &knz, &jlen, &h__[incol + 1 + j2 + jcol * 
			    h_dim1], ldh, &wh[kzs + 1 + wh_dim1], ldwh, (
			    ftnlen)3);

/*                 ==== Multiply by U21**H ==== */

#line 770 "claqr5.f"
		    claset_("ALL", &kzs, &jlen, &c_b1, &c_b1, &wh[wh_offset], 
			    ldwh, (ftnlen)3);
#line 771 "claqr5.f"
		    ctrmm_("L", "U", "C", "N", &knz, &jlen, &c_b2, &u[j2 + 1 
			    + (kzs + 1) * u_dim1], ldu, &wh[kzs + 1 + wh_dim1]
			    , ldwh, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1);

/*                 ==== Multiply top of H by U11**H ==== */

#line 777 "claqr5.f"
		    cgemm_("C", "N", &i2, &jlen, &j2, &c_b2, &u[u_offset], 
			    ldu, &h__[incol + 1 + jcol * h_dim1], ldh, &c_b2, 
			    &wh[wh_offset], ldwh, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy top of H to bottom of WH ==== */

#line 782 "claqr5.f"
		    clacpy_("ALL", &j2, &jlen, &h__[incol + 1 + jcol * h_dim1]
			    , ldh, &wh[i2 + 1 + wh_dim1], ldwh, (ftnlen)3);

/*                 ==== Multiply by U21**H ==== */

#line 787 "claqr5.f"
		    ctrmm_("L", "L", "C", "N", &j2, &jlen, &c_b2, &u[(i2 + 1) 
			    * u_dim1 + 1], ldu, &wh[i2 + 1 + wh_dim1], ldwh, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 ==== Multiply by U22 ==== */

#line 792 "claqr5.f"
		    i__5 = i4 - i2;
#line 792 "claqr5.f"
		    i__7 = j4 - j2;
#line 792 "claqr5.f"
		    cgemm_("C", "N", &i__5, &jlen, &i__7, &c_b2, &u[j2 + 1 + (
			    i2 + 1) * u_dim1], ldu, &h__[incol + 1 + j2 + 
			    jcol * h_dim1], ldh, &c_b2, &wh[i2 + 1 + wh_dim1],
			     ldwh, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy it back ==== */

#line 799 "claqr5.f"
		    clacpy_("ALL", &kdu, &jlen, &wh[wh_offset], ldwh, &h__[
			    incol + 1 + jcol * h_dim1], ldh, (ftnlen)3);
#line 801 "claqr5.f"
/* L180: */
#line 801 "claqr5.f"
		}

/*              ==== Vertical multiply ==== */

#line 805 "claqr5.f"
		i__3 = max(incol,*ktop) - 1;
#line 805 "claqr5.f"
		i__4 = *nv;
#line 805 "claqr5.f"
		for (jrow = jtop; i__4 < 0 ? jrow >= i__3 : jrow <= i__3; 
			jrow += i__4) {
/* Computing MIN */
#line 806 "claqr5.f"
		    i__5 = *nv, i__7 = max(incol,*ktop) - jrow;
#line 806 "claqr5.f"
		    jlen = min(i__5,i__7);

/*                 ==== Copy right of H to scratch (the first KZS */
/*                 .    columns get multiplied by zero) ==== */

#line 811 "claqr5.f"
		    clacpy_("ALL", &jlen, &knz, &h__[jrow + (incol + 1 + j2) *
			     h_dim1], ldh, &wv[(kzs + 1) * wv_dim1 + 1], ldwv,
			     (ftnlen)3);

/*                 ==== Multiply by U21 ==== */

#line 816 "claqr5.f"
		    claset_("ALL", &jlen, &kzs, &c_b1, &c_b1, &wv[wv_offset], 
			    ldwv, (ftnlen)3);
#line 817 "claqr5.f"
		    ctrmm_("R", "U", "N", "N", &jlen, &knz, &c_b2, &u[j2 + 1 
			    + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) * 
			    wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1);

/*                 ==== Multiply by U11 ==== */

#line 823 "claqr5.f"
		    cgemm_("N", "N", &jlen, &i2, &j2, &c_b2, &h__[jrow + (
			    incol + 1) * h_dim1], ldh, &u[u_offset], ldu, &
			    c_b2, &wv[wv_offset], ldwv, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy left of H to right of scratch ==== */

#line 829 "claqr5.f"
		    clacpy_("ALL", &jlen, &j2, &h__[jrow + (incol + 1) * 
			    h_dim1], ldh, &wv[(i2 + 1) * wv_dim1 + 1], ldwv, (
			    ftnlen)3);

/*                 ==== Multiply by U21 ==== */

#line 834 "claqr5.f"
		    i__5 = i4 - i2;
#line 834 "claqr5.f"
		    ctrmm_("R", "L", "N", "N", &jlen, &i__5, &c_b2, &u[(i2 + 
			    1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * wv_dim1 + 1]
			    , ldwv, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			    1);

/*                 ==== Multiply by U22 ==== */

#line 839 "claqr5.f"
		    i__5 = i4 - i2;
#line 839 "claqr5.f"
		    i__7 = j4 - j2;
#line 839 "claqr5.f"
		    cgemm_("N", "N", &jlen, &i__5, &i__7, &c_b2, &h__[jrow + (
			    incol + 1 + j2) * h_dim1], ldh, &u[j2 + 1 + (i2 + 
			    1) * u_dim1], ldu, &c_b2, &wv[(i2 + 1) * wv_dim1 
			    + 1], ldwv, (ftnlen)1, (ftnlen)1);

/*                 ==== Copy it back ==== */

#line 846 "claqr5.f"
		    clacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, &h__[
			    jrow + (incol + 1) * h_dim1], ldh, (ftnlen)3);
#line 848 "claqr5.f"
/* L190: */
#line 848 "claqr5.f"
		}

/*              ==== Multiply Z (also vertical) ==== */

#line 852 "claqr5.f"
		if (*wantz) {
#line 853 "claqr5.f"
		    i__4 = *ihiz;
#line 853 "claqr5.f"
		    i__3 = *nv;
#line 853 "claqr5.f"
		    for (jrow = *iloz; i__3 < 0 ? jrow >= i__4 : jrow <= i__4;
			     jrow += i__3) {
/* Computing MIN */
#line 854 "claqr5.f"
			i__5 = *nv, i__7 = *ihiz - jrow + 1;
#line 854 "claqr5.f"
			jlen = min(i__5,i__7);

/*                    ==== Copy right of Z to left of scratch (first */
/*                    .     KZS columns get multiplied by zero) ==== */

#line 859 "claqr5.f"
			clacpy_("ALL", &jlen, &knz, &z__[jrow + (incol + 1 + 
				j2) * z_dim1], ldz, &wv[(kzs + 1) * wv_dim1 + 
				1], ldwv, (ftnlen)3);

/*                    ==== Multiply by U12 ==== */

#line 865 "claqr5.f"
			claset_("ALL", &jlen, &kzs, &c_b1, &c_b1, &wv[
				wv_offset], ldwv, (ftnlen)3);
#line 867 "claqr5.f"
			ctrmm_("R", "U", "N", "N", &jlen, &knz, &c_b2, &u[j2 
				+ 1 + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) 
				* wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

/*                    ==== Multiply by U11 ==== */

#line 873 "claqr5.f"
			cgemm_("N", "N", &jlen, &i2, &j2, &c_b2, &z__[jrow + (
				incol + 1) * z_dim1], ldz, &u[u_offset], ldu, 
				&c_b2, &wv[wv_offset], ldwv, (ftnlen)1, (
				ftnlen)1);

/*                    ==== Copy left of Z to right of scratch ==== */

#line 879 "claqr5.f"
			clacpy_("ALL", &jlen, &j2, &z__[jrow + (incol + 1) * 
				z_dim1], ldz, &wv[(i2 + 1) * wv_dim1 + 1], 
				ldwv, (ftnlen)3);

/*                    ==== Multiply by U21 ==== */

#line 884 "claqr5.f"
			i__5 = i4 - i2;
#line 884 "claqr5.f"
			ctrmm_("R", "L", "N", "N", &jlen, &i__5, &c_b2, &u[(
				i2 + 1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * 
				wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)1, (
				ftnlen)1, (ftnlen)1);

/*                    ==== Multiply by U22 ==== */

#line 890 "claqr5.f"
			i__5 = i4 - i2;
#line 890 "claqr5.f"
			i__7 = j4 - j2;
#line 890 "claqr5.f"
			cgemm_("N", "N", &jlen, &i__5, &i__7, &c_b2, &z__[
				jrow + (incol + 1 + j2) * z_dim1], ldz, &u[j2 
				+ 1 + (i2 + 1) * u_dim1], ldu, &c_b2, &wv[(i2 
				+ 1) * wv_dim1 + 1], ldwv, (ftnlen)1, (ftnlen)
				1);

/*                    ==== Copy the result back to Z ==== */

#line 897 "claqr5.f"
			clacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, &
				z__[jrow + (incol + 1) * z_dim1], ldz, (
				ftnlen)3);
#line 899 "claqr5.f"
/* L200: */
#line 899 "claqr5.f"
		    }
#line 900 "claqr5.f"
		}
#line 901 "claqr5.f"
	    }
#line 902 "claqr5.f"
	}
#line 903 "claqr5.f"
/* L210: */
#line 903 "claqr5.f"
    }

/*     ==== End of CLAQR5 ==== */

#line 907 "claqr5.f"
    return 0;
} /* claqr5_ */

