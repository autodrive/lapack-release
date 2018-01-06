#line 1 "dlaqr2.f"
/* dlaqr2.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqr2.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static logical c_true = TRUE_;

/* > \brief \b DLAQR2 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and d
eflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation). 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
/*                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, */
/*                          LDT, NV, WV, LDWV, WORK, LWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, */
/*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), */
/*      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLAQR2 is identical to DLAQR3 except that it avoids */
/* >    recursion by calling DLAHQR instead of DLAQR4. */
/* > */
/* >    Aggressive early deflation: */
/* > */
/* >    This subroutine accepts as input an upper Hessenberg matrix */
/* >    H and performs an orthogonal similarity transformation */
/* >    designed to detect and deflate fully converged eigenvalues from */
/* >    a trailing principal submatrix.  On output H has been over- */
/* >    written by a new Hessenberg matrix that is a perturbation of */
/* >    an orthogonal similarity transformation of H.  It is to be */
/* >    hoped that the final version of H has many zero subdiagonal */
/* >    entries. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          If .TRUE., then the Hessenberg matrix H is fully updated */
/* >          so that the quasi-triangular Schur factor may be */
/* >          computed (in cooperation with the calling subroutine). */
/* >          If .FALSE., then only enough of H is updated to preserve */
/* >          the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          If .TRUE., then the orthogonal matrix Z is updated so */
/* >          so that the orthogonal Schur factor may be computed */
/* >          (in cooperation with the calling subroutine). */
/* >          If .FALSE., then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H and (if WANTZ is .TRUE.) the */
/* >          order of the orthogonal matrix Z. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* >          KTOP is INTEGER */
/* >          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0. */
/* >          KBOT and KTOP together determine an isolated block */
/* >          along the diagonal of the Hessenberg matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* >          KBOT is INTEGER */
/* >          It is assumed without a check that either */
/* >          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together */
/* >          determine an isolated block along the diagonal of the */
/* >          Hessenberg matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* >          NW is INTEGER */
/* >          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1). */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION array, dimension (LDH,N) */
/* >          On input the initial N-by-N section of H stores the */
/* >          Hessenberg matrix undergoing aggressive early deflation. */
/* >          On output H has been transformed by an orthogonal */
/* >          similarity transformation, perturbed, and the returned */
/* >          to Hessenberg form that (it is to be hoped) has some */
/* >          zero subdiagonal entries. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          Leading dimension of H just as declared in the calling */
/* >          subroutine.  N .LE. LDH */
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
/* >          Specify the rows of Z to which transformations must be */
/* >          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ,N) */
/* >          IF WANTZ is .TRUE., then on output, the orthogonal */
/* >          similarity transformation mentioned above has been */
/* >          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right. */
/* >          If WANTZ is .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of Z just as declared in the */
/* >          calling subroutine.  1 .LE. LDZ. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* >          NS is INTEGER */
/* >          The number of unconverged (ie approximate) eigenvalues */
/* >          returned in SR and SI that may be used as shifts by the */
/* >          calling subroutine. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* >          ND is INTEGER */
/* >          The number of converged eigenvalues uncovered by this */
/* >          subroutine. */
/* > \endverbatim */
/* > */
/* > \param[out] SR */
/* > \verbatim */
/* >          SR is DOUBLE PRECISION array, dimension (KBOT) */
/* > \endverbatim */
/* > */
/* > \param[out] SI */
/* > \verbatim */
/* >          SI is DOUBLE PRECISION array, dimension (KBOT) */
/* >          On output, the real and imaginary parts of approximate */
/* >          eigenvalues that may be used for shifts are stored in */
/* >          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and */
/* >          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively. */
/* >          The real and imaginary parts of converged eigenvalues */
/* >          are stored in SR(KBOT-ND+1) through SR(KBOT) and */
/* >          SI(KBOT-ND+1) through SI(KBOT), respectively. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension (LDV,NW) */
/* >          An NW-by-NW work array. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of V just as declared in the */
/* >          calling subroutine.  NW .LE. LDV */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* >          NH is INTEGER */
/* >          The number of columns of T.  NH.GE.NW. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of T just as declared in the */
/* >          calling subroutine.  NW .LE. LDT */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* >          NV is INTEGER */
/* >          The number of rows of work array WV available for */
/* >          workspace.  NV.GE.NW. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* >          WV is DOUBLE PRECISION array, dimension (LDWV,NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* >          LDWV is INTEGER */
/* >          The leading dimension of W just as declared in the */
/* >          calling subroutine.  NW .LE. LDV */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* >          On exit, WORK(1) is set to an estimate of the optimal value */
/* >          of LWORK for the given values of N, NW, KTOP and KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the work array WORK.  LWORK = 2*NW */
/* >          suffices, but greater efficiency may result from larger */
/* >          values of LWORK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; DLAQR2 */
/* >          only estimates the optimal workspace size for the given */
/* >          values of N, NW, KTOP and KBOT.  The estimate is returned */
/* >          in WORK(1).  No error message related to LWORK is issued */
/* >          by XERBLA.  Neither H nor Z are accessed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaqr2_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, doublereal *h__, integer *
	ldh, integer *iloz, integer *ihiz, doublereal *z__, integer *ldz, 
	integer *ns, integer *nd, doublereal *sr, doublereal *si, doublereal *
	v, integer *ldv, integer *nh, doublereal *t, integer *ldt, integer *
	nv, doublereal *wv, integer *ldwv, doublereal *work, integer *lwork)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, 
	    wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, aa, bb, cc, dd, cs, sn;
    static integer jw;
    static doublereal evi, evk, foo;
    static integer kln;
    static doublereal tau, ulp;
    static integer lwk1, lwk2;
    static doublereal beta;
    static integer kend, kcol, info, ifst, ilst, ltop, krow;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), dgemm_(char *, char *, integer *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical bulge;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer infqr, kwtop;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlarfg_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *), dlahqr_(logical *, logical *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int dtrexc_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), dormhr_(char *, char *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static logical sorted;
    static doublereal smlnum;
    static integer lwkopt;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Estimate optimal workspace. ==== */

#line 325 "dlaqr2.f"
    /* Parameter adjustments */
#line 325 "dlaqr2.f"
    h_dim1 = *ldh;
#line 325 "dlaqr2.f"
    h_offset = 1 + h_dim1;
#line 325 "dlaqr2.f"
    h__ -= h_offset;
#line 325 "dlaqr2.f"
    z_dim1 = *ldz;
#line 325 "dlaqr2.f"
    z_offset = 1 + z_dim1;
#line 325 "dlaqr2.f"
    z__ -= z_offset;
#line 325 "dlaqr2.f"
    --sr;
#line 325 "dlaqr2.f"
    --si;
#line 325 "dlaqr2.f"
    v_dim1 = *ldv;
#line 325 "dlaqr2.f"
    v_offset = 1 + v_dim1;
#line 325 "dlaqr2.f"
    v -= v_offset;
#line 325 "dlaqr2.f"
    t_dim1 = *ldt;
#line 325 "dlaqr2.f"
    t_offset = 1 + t_dim1;
#line 325 "dlaqr2.f"
    t -= t_offset;
#line 325 "dlaqr2.f"
    wv_dim1 = *ldwv;
#line 325 "dlaqr2.f"
    wv_offset = 1 + wv_dim1;
#line 325 "dlaqr2.f"
    wv -= wv_offset;
#line 325 "dlaqr2.f"
    --work;
#line 325 "dlaqr2.f"

#line 325 "dlaqr2.f"
    /* Function Body */
/* Computing MIN */
#line 325 "dlaqr2.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 325 "dlaqr2.f"
    jw = min(i__1,i__2);
#line 326 "dlaqr2.f"
    if (jw <= 2) {
#line 327 "dlaqr2.f"
	lwkopt = 1;
#line 328 "dlaqr2.f"
    } else {

/*        ==== Workspace query call to DGEHRD ==== */

#line 332 "dlaqr2.f"
	i__1 = jw - 1;
#line 332 "dlaqr2.f"
	dgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &
		c_n1, &info);
#line 333 "dlaqr2.f"
	lwk1 = (integer) work[1];

/*        ==== Workspace query call to DORMHR ==== */

#line 337 "dlaqr2.f"
	i__1 = jw - 1;
#line 337 "dlaqr2.f"
	dormhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
		 &v[v_offset], ldv, &work[1], &c_n1, &info, (ftnlen)1, (
		ftnlen)1);
#line 339 "dlaqr2.f"
	lwk2 = (integer) work[1];

/*        ==== Optimal workspace ==== */

#line 343 "dlaqr2.f"
	lwkopt = jw + max(lwk1,lwk2);
#line 344 "dlaqr2.f"
    }

/*     ==== Quick return in case of workspace query. ==== */

#line 348 "dlaqr2.f"
    if (*lwork == -1) {
#line 349 "dlaqr2.f"
	work[1] = (doublereal) lwkopt;
#line 350 "dlaqr2.f"
	return 0;
#line 351 "dlaqr2.f"
    }

/*     ==== Nothing to do ... */
/*     ... for an empty active block ... ==== */
#line 355 "dlaqr2.f"
    *ns = 0;
#line 356 "dlaqr2.f"
    *nd = 0;
#line 357 "dlaqr2.f"
    work[1] = 1.;
#line 358 "dlaqr2.f"
    if (*ktop > *kbot) {
#line 358 "dlaqr2.f"
	return 0;
#line 358 "dlaqr2.f"
    }
/*     ... nor for an empty deflation window. ==== */
#line 361 "dlaqr2.f"
    if (*nw < 1) {
#line 361 "dlaqr2.f"
	return 0;
#line 361 "dlaqr2.f"
    }

/*     ==== Machine constants ==== */

#line 366 "dlaqr2.f"
    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
#line 367 "dlaqr2.f"
    safmax = 1. / safmin;
#line 368 "dlaqr2.f"
    dlabad_(&safmin, &safmax);
#line 369 "dlaqr2.f"
    ulp = dlamch_("PRECISION", (ftnlen)9);
#line 370 "dlaqr2.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Setup deflation window ==== */

/* Computing MIN */
#line 374 "dlaqr2.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 374 "dlaqr2.f"
    jw = min(i__1,i__2);
#line 375 "dlaqr2.f"
    kwtop = *kbot - jw + 1;
#line 376 "dlaqr2.f"
    if (kwtop == *ktop) {
#line 377 "dlaqr2.f"
	s = 0.;
#line 378 "dlaqr2.f"
    } else {
#line 379 "dlaqr2.f"
	s = h__[kwtop + (kwtop - 1) * h_dim1];
#line 380 "dlaqr2.f"
    }

#line 382 "dlaqr2.f"
    if (*kbot == kwtop) {

/*        ==== 1-by-1 deflation window: not much to do ==== */

#line 386 "dlaqr2.f"
	sr[kwtop] = h__[kwtop + kwtop * h_dim1];
#line 387 "dlaqr2.f"
	si[kwtop] = 0.;
#line 388 "dlaqr2.f"
	*ns = 1;
#line 389 "dlaqr2.f"
	*nd = 0;
/* Computing MAX */
#line 390 "dlaqr2.f"
	d__2 = smlnum, d__3 = ulp * (d__1 = h__[kwtop + kwtop * h_dim1], abs(
		d__1));
#line 390 "dlaqr2.f"
	if (abs(s) <= max(d__2,d__3)) {
#line 392 "dlaqr2.f"
	    *ns = 0;
#line 393 "dlaqr2.f"
	    *nd = 1;
#line 394 "dlaqr2.f"
	    if (kwtop > *ktop) {
#line 394 "dlaqr2.f"
		h__[kwtop + (kwtop - 1) * h_dim1] = 0.;
#line 394 "dlaqr2.f"
	    }
#line 396 "dlaqr2.f"
	}
#line 397 "dlaqr2.f"
	work[1] = 1.;
#line 398 "dlaqr2.f"
	return 0;
#line 399 "dlaqr2.f"
    }

/*     ==== Convert to spike-triangular form.  (In case of a */
/*     .    rare QR failure, this routine continues to do */
/*     .    aggressive early deflation using that part of */
/*     .    the deflation window that converged using INFQR */
/*     .    here and there to keep track.) ==== */

#line 407 "dlaqr2.f"
    dlacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], 
	    ldt, (ftnlen)1);
#line 408 "dlaqr2.f"
    i__1 = jw - 1;
#line 408 "dlaqr2.f"
    i__2 = *ldh + 1;
#line 408 "dlaqr2.f"
    i__3 = *ldt + 1;
#line 408 "dlaqr2.f"
    dcopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &
	    i__3);

#line 410 "dlaqr2.f"
    dlaset_("A", &jw, &jw, &c_b12, &c_b13, &v[v_offset], ldv, (ftnlen)1);
#line 411 "dlaqr2.f"
    dlahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[kwtop], 
	    &si[kwtop], &c__1, &jw, &v[v_offset], ldv, &infqr);

/*     ==== DTREXC needs a clean margin near the diagonal ==== */

#line 416 "dlaqr2.f"
    i__1 = jw - 3;
#line 416 "dlaqr2.f"
    for (j = 1; j <= i__1; ++j) {
#line 417 "dlaqr2.f"
	t[j + 2 + j * t_dim1] = 0.;
#line 418 "dlaqr2.f"
	t[j + 3 + j * t_dim1] = 0.;
#line 419 "dlaqr2.f"
/* L10: */
#line 419 "dlaqr2.f"
    }
#line 420 "dlaqr2.f"
    if (jw > 2) {
#line 420 "dlaqr2.f"
	t[jw + (jw - 2) * t_dim1] = 0.;
#line 420 "dlaqr2.f"
    }

/*     ==== Deflation detection loop ==== */

#line 425 "dlaqr2.f"
    *ns = jw;
#line 426 "dlaqr2.f"
    ilst = infqr + 1;
#line 427 "dlaqr2.f"
L20:
#line 428 "dlaqr2.f"
    if (ilst <= *ns) {
#line 429 "dlaqr2.f"
	if (*ns == 1) {
#line 430 "dlaqr2.f"
	    bulge = FALSE_;
#line 431 "dlaqr2.f"
	} else {
#line 432 "dlaqr2.f"
	    bulge = t[*ns + (*ns - 1) * t_dim1] != 0.;
#line 433 "dlaqr2.f"
	}

/*        ==== Small spike tip test for deflation ==== */

#line 437 "dlaqr2.f"
	if (! bulge) {

/*           ==== Real eigenvalue ==== */

#line 441 "dlaqr2.f"
	    foo = (d__1 = t[*ns + *ns * t_dim1], abs(d__1));
#line 442 "dlaqr2.f"
	    if (foo == 0.) {
#line 442 "dlaqr2.f"
		foo = abs(s);
#line 442 "dlaqr2.f"
	    }
/* Computing MAX */
#line 444 "dlaqr2.f"
	    d__2 = smlnum, d__3 = ulp * foo;
#line 444 "dlaqr2.f"
	    if ((d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)) <= max(d__2,d__3))
		     {

/*              ==== Deflatable ==== */

#line 448 "dlaqr2.f"
		--(*ns);
#line 449 "dlaqr2.f"
	    } else {

/*              ==== Undeflatable.   Move it up out of the way. */
/*              .    (DTREXC can not fail in this case.) ==== */

#line 454 "dlaqr2.f"
		ifst = *ns;
#line 455 "dlaqr2.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 457 "dlaqr2.f"
		++ilst;
#line 458 "dlaqr2.f"
	    }
#line 459 "dlaqr2.f"
	} else {

/*           ==== Complex conjugate pair ==== */

#line 463 "dlaqr2.f"
	    foo = (d__3 = t[*ns + *ns * t_dim1], abs(d__3)) + sqrt((d__1 = t[*
		    ns + (*ns - 1) * t_dim1], abs(d__1))) * sqrt((d__2 = t[*
		    ns - 1 + *ns * t_dim1], abs(d__2)));
#line 465 "dlaqr2.f"
	    if (foo == 0.) {
#line 465 "dlaqr2.f"
		foo = abs(s);
#line 465 "dlaqr2.f"
	    }
/* Computing MAX */
#line 467 "dlaqr2.f"
	    d__3 = (d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)), d__4 = (d__2 =
		     s * v[(*ns - 1) * v_dim1 + 1], abs(d__2));
/* Computing MAX */
#line 467 "dlaqr2.f"
	    d__5 = smlnum, d__6 = ulp * foo;
#line 467 "dlaqr2.f"
	    if (max(d__3,d__4) <= max(d__5,d__6)) {

/*              ==== Deflatable ==== */

#line 472 "dlaqr2.f"
		*ns += -2;
#line 473 "dlaqr2.f"
	    } else {

/*              ==== Undeflatable. Move them up out of the way. */
/*              .    Fortunately, DTREXC does the right thing with */
/*              .    ILST in case of a rare exchange failure. ==== */

#line 479 "dlaqr2.f"
		ifst = *ns;
#line 480 "dlaqr2.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 482 "dlaqr2.f"
		ilst += 2;
#line 483 "dlaqr2.f"
	    }
#line 484 "dlaqr2.f"
	}

/*        ==== End deflation detection loop ==== */

#line 488 "dlaqr2.f"
	goto L20;
#line 489 "dlaqr2.f"
    }

/*        ==== Return to Hessenberg form ==== */

#line 493 "dlaqr2.f"
    if (*ns == 0) {
#line 493 "dlaqr2.f"
	s = 0.;
#line 493 "dlaqr2.f"
    }

#line 496 "dlaqr2.f"
    if (*ns < jw) {

/*        ==== sorting diagonal blocks of T improves accuracy for */
/*        .    graded matrices.  Bubble sort deals well with */
/*        .    exchange failures. ==== */

#line 502 "dlaqr2.f"
	sorted = FALSE_;
#line 503 "dlaqr2.f"
	i__ = *ns + 1;
#line 504 "dlaqr2.f"
L30:
#line 505 "dlaqr2.f"
	if (sorted) {
#line 505 "dlaqr2.f"
	    goto L50;
#line 505 "dlaqr2.f"
	}
#line 507 "dlaqr2.f"
	sorted = TRUE_;

#line 509 "dlaqr2.f"
	kend = i__ - 1;
#line 510 "dlaqr2.f"
	i__ = infqr + 1;
#line 511 "dlaqr2.f"
	if (i__ == *ns) {
#line 512 "dlaqr2.f"
	    k = i__ + 1;
#line 513 "dlaqr2.f"
	} else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
#line 514 "dlaqr2.f"
	    k = i__ + 1;
#line 515 "dlaqr2.f"
	} else {
#line 516 "dlaqr2.f"
	    k = i__ + 2;
#line 517 "dlaqr2.f"
	}
#line 518 "dlaqr2.f"
L40:
#line 519 "dlaqr2.f"
	if (k <= kend) {
#line 520 "dlaqr2.f"
	    if (k == i__ + 1) {
#line 521 "dlaqr2.f"
		evi = (d__1 = t[i__ + i__ * t_dim1], abs(d__1));
#line 522 "dlaqr2.f"
	    } else {
#line 523 "dlaqr2.f"
		evi = (d__3 = t[i__ + i__ * t_dim1], abs(d__3)) + sqrt((d__1 =
			 t[i__ + 1 + i__ * t_dim1], abs(d__1))) * sqrt((d__2 =
			 t[i__ + (i__ + 1) * t_dim1], abs(d__2)));
#line 525 "dlaqr2.f"
	    }

#line 527 "dlaqr2.f"
	    if (k == kend) {
#line 528 "dlaqr2.f"
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
#line 529 "dlaqr2.f"
	    } else if (t[k + 1 + k * t_dim1] == 0.) {
#line 530 "dlaqr2.f"
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
#line 531 "dlaqr2.f"
	    } else {
#line 532 "dlaqr2.f"
		evk = (d__3 = t[k + k * t_dim1], abs(d__3)) + sqrt((d__1 = t[
			k + 1 + k * t_dim1], abs(d__1))) * sqrt((d__2 = t[k + 
			(k + 1) * t_dim1], abs(d__2)));
#line 534 "dlaqr2.f"
	    }

#line 536 "dlaqr2.f"
	    if (evi >= evk) {
#line 537 "dlaqr2.f"
		i__ = k;
#line 538 "dlaqr2.f"
	    } else {
#line 539 "dlaqr2.f"
		sorted = FALSE_;
#line 540 "dlaqr2.f"
		ifst = i__;
#line 541 "dlaqr2.f"
		ilst = k;
#line 542 "dlaqr2.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 544 "dlaqr2.f"
		if (info == 0) {
#line 545 "dlaqr2.f"
		    i__ = ilst;
#line 546 "dlaqr2.f"
		} else {
#line 547 "dlaqr2.f"
		    i__ = k;
#line 548 "dlaqr2.f"
		}
#line 549 "dlaqr2.f"
	    }
#line 550 "dlaqr2.f"
	    if (i__ == kend) {
#line 551 "dlaqr2.f"
		k = i__ + 1;
#line 552 "dlaqr2.f"
	    } else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
#line 553 "dlaqr2.f"
		k = i__ + 1;
#line 554 "dlaqr2.f"
	    } else {
#line 555 "dlaqr2.f"
		k = i__ + 2;
#line 556 "dlaqr2.f"
	    }
#line 557 "dlaqr2.f"
	    goto L40;
#line 558 "dlaqr2.f"
	}
#line 559 "dlaqr2.f"
	goto L30;
#line 560 "dlaqr2.f"
L50:
#line 561 "dlaqr2.f"
	;
#line 561 "dlaqr2.f"
    }

/*     ==== Restore shift/eigenvalue array from T ==== */

#line 565 "dlaqr2.f"
    i__ = jw;
#line 566 "dlaqr2.f"
L60:
#line 567 "dlaqr2.f"
    if (i__ >= infqr + 1) {
#line 568 "dlaqr2.f"
	if (i__ == infqr + 1) {
#line 569 "dlaqr2.f"
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
#line 570 "dlaqr2.f"
	    si[kwtop + i__ - 1] = 0.;
#line 571 "dlaqr2.f"
	    --i__;
#line 572 "dlaqr2.f"
	} else if (t[i__ + (i__ - 1) * t_dim1] == 0.) {
#line 573 "dlaqr2.f"
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
#line 574 "dlaqr2.f"
	    si[kwtop + i__ - 1] = 0.;
#line 575 "dlaqr2.f"
	    --i__;
#line 576 "dlaqr2.f"
	} else {
#line 577 "dlaqr2.f"
	    aa = t[i__ - 1 + (i__ - 1) * t_dim1];
#line 578 "dlaqr2.f"
	    cc = t[i__ + (i__ - 1) * t_dim1];
#line 579 "dlaqr2.f"
	    bb = t[i__ - 1 + i__ * t_dim1];
#line 580 "dlaqr2.f"
	    dd = t[i__ + i__ * t_dim1];
#line 581 "dlaqr2.f"
	    dlanv2_(&aa, &bb, &cc, &dd, &sr[kwtop + i__ - 2], &si[kwtop + i__ 
		    - 2], &sr[kwtop + i__ - 1], &si[kwtop + i__ - 1], &cs, &
		    sn);
#line 584 "dlaqr2.f"
	    i__ += -2;
#line 585 "dlaqr2.f"
	}
#line 586 "dlaqr2.f"
	goto L60;
#line 587 "dlaqr2.f"
    }

#line 589 "dlaqr2.f"
    if (*ns < jw || s == 0.) {
#line 590 "dlaqr2.f"
	if (*ns > 1 && s != 0.) {

/*           ==== Reflect spike back into lower triangle ==== */

#line 594 "dlaqr2.f"
	    dcopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
#line 595 "dlaqr2.f"
	    beta = work[1];
#line 596 "dlaqr2.f"
	    dlarfg_(ns, &beta, &work[2], &c__1, &tau);
#line 597 "dlaqr2.f"
	    work[1] = 1.;

#line 599 "dlaqr2.f"
	    i__1 = jw - 2;
#line 599 "dlaqr2.f"
	    i__2 = jw - 2;
#line 599 "dlaqr2.f"
	    dlaset_("L", &i__1, &i__2, &c_b12, &c_b12, &t[t_dim1 + 3], ldt, (
		    ftnlen)1);

#line 601 "dlaqr2.f"
	    dlarf_("L", ns, &jw, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 603 "dlaqr2.f"
	    dlarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 605 "dlaqr2.f"
	    dlarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &
		    work[jw + 1], (ftnlen)1);

#line 608 "dlaqr2.f"
	    i__1 = *lwork - jw;
#line 608 "dlaqr2.f"
	    dgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1]
		    , &i__1, &info);
#line 610 "dlaqr2.f"
	}

/*        ==== Copy updated reduced window into place ==== */

#line 614 "dlaqr2.f"
	if (kwtop > 1) {
#line 614 "dlaqr2.f"
	    h__[kwtop + (kwtop - 1) * h_dim1] = s * v[v_dim1 + 1];
#line 614 "dlaqr2.f"
	}
#line 616 "dlaqr2.f"
	dlacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1]
		, ldh, (ftnlen)1);
#line 617 "dlaqr2.f"
	i__1 = jw - 1;
#line 617 "dlaqr2.f"
	i__2 = *ldt + 1;
#line 617 "dlaqr2.f"
	i__3 = *ldh + 1;
#line 617 "dlaqr2.f"
	dcopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1],
		 &i__3);

/*        ==== Accumulate orthogonal matrix in order update */
/*        .    H and Z, if requested.  ==== */

#line 623 "dlaqr2.f"
	if (*ns > 1 && s != 0.) {
#line 623 "dlaqr2.f"
	    i__1 = *lwork - jw;
#line 623 "dlaqr2.f"
	    dormhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
		     &v[v_offset], ldv, &work[jw + 1], &i__1, &info, (ftnlen)
		    1, (ftnlen)1);
#line 623 "dlaqr2.f"
	}

/*        ==== Update vertical slab in H ==== */

#line 629 "dlaqr2.f"
	if (*wantt) {
#line 630 "dlaqr2.f"
	    ltop = 1;
#line 631 "dlaqr2.f"
	} else {
#line 632 "dlaqr2.f"
	    ltop = *ktop;
#line 633 "dlaqr2.f"
	}
#line 634 "dlaqr2.f"
	i__1 = kwtop - 1;
#line 634 "dlaqr2.f"
	i__2 = *nv;
#line 634 "dlaqr2.f"
	for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += 
		i__2) {
/* Computing MIN */
#line 635 "dlaqr2.f"
	    i__3 = *nv, i__4 = kwtop - krow;
#line 635 "dlaqr2.f"
	    kln = min(i__3,i__4);
#line 636 "dlaqr2.f"
	    dgemm_("N", "N", &kln, &jw, &jw, &c_b13, &h__[krow + kwtop * 
		    h_dim1], ldh, &v[v_offset], ldv, &c_b12, &wv[wv_offset], 
		    ldwv, (ftnlen)1, (ftnlen)1);
#line 638 "dlaqr2.f"
	    dlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * 
		    h_dim1], ldh, (ftnlen)1);
#line 639 "dlaqr2.f"
/* L70: */
#line 639 "dlaqr2.f"
	}

/*        ==== Update horizontal slab in H ==== */

#line 643 "dlaqr2.f"
	if (*wantt) {
#line 644 "dlaqr2.f"
	    i__2 = *n;
#line 644 "dlaqr2.f"
	    i__1 = *nh;
#line 644 "dlaqr2.f"
	    for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; 
		    kcol += i__1) {
/* Computing MIN */
#line 645 "dlaqr2.f"
		i__3 = *nh, i__4 = *n - kcol + 1;
#line 645 "dlaqr2.f"
		kln = min(i__3,i__4);
#line 646 "dlaqr2.f"
		dgemm_("C", "N", &jw, &kln, &jw, &c_b13, &v[v_offset], ldv, &
			h__[kwtop + kcol * h_dim1], ldh, &c_b12, &t[t_offset],
			 ldt, (ftnlen)1, (ftnlen)1);
#line 648 "dlaqr2.f"
		dlacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol *
			 h_dim1], ldh, (ftnlen)1);
#line 650 "dlaqr2.f"
/* L80: */
#line 650 "dlaqr2.f"
	    }
#line 651 "dlaqr2.f"
	}

/*        ==== Update vertical slab in Z ==== */

#line 655 "dlaqr2.f"
	if (*wantz) {
#line 656 "dlaqr2.f"
	    i__1 = *ihiz;
#line 656 "dlaqr2.f"
	    i__2 = *nv;
#line 656 "dlaqr2.f"
	    for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow +=
		     i__2) {
/* Computing MIN */
#line 657 "dlaqr2.f"
		i__3 = *nv, i__4 = *ihiz - krow + 1;
#line 657 "dlaqr2.f"
		kln = min(i__3,i__4);
#line 658 "dlaqr2.f"
		dgemm_("N", "N", &kln, &jw, &jw, &c_b13, &z__[krow + kwtop * 
			z_dim1], ldz, &v[v_offset], ldv, &c_b12, &wv[
			wv_offset], ldwv, (ftnlen)1, (ftnlen)1);
#line 660 "dlaqr2.f"
		dlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + 
			kwtop * z_dim1], ldz, (ftnlen)1);
#line 662 "dlaqr2.f"
/* L90: */
#line 662 "dlaqr2.f"
	    }
#line 663 "dlaqr2.f"
	}
#line 664 "dlaqr2.f"
    }

/*     ==== Return the number of deflations ... ==== */

#line 668 "dlaqr2.f"
    *nd = jw - *ns;

/*     ==== ... and the number of shifts. (Subtracting */
/*     .    INFQR from the spike length takes care */
/*     .    of the case of a rare QR failure while */
/*     .    calculating eigenvalues of the deflation */
/*     .    window.)  ==== */

#line 676 "dlaqr2.f"
    *ns -= infqr;

/*      ==== Return optimal workspace. ==== */

#line 680 "dlaqr2.f"
    work[1] = (doublereal) lwkopt;

/*     ==== End of DLAQR2 ==== */

#line 684 "dlaqr2.f"
    return 0;
} /* dlaqr2_ */

