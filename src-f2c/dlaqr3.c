#line 1 "dlaqr3.f"
/* dlaqr3.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqr3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__12 = 12;

/* > \brief \b DLAQR3 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and d
eflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation). 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQR3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
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
/* >    Aggressive early deflation: */
/* > */
/* >    DLAQR3 accepts as input an upper Hessenberg matrix */
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
/* >          LDH is integer */
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
/* >          LDZ is integer */
/* >          The leading dimension of Z just as declared in the */
/* >          calling subroutine.  1 .LE. LDZ. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* >          NS is integer */
/* >          The number of unconverged (ie approximate) eigenvalues */
/* >          returned in SR and SI that may be used as shifts by the */
/* >          calling subroutine. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* >          ND is integer */
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
/* >          LDV is integer scalar */
/* >          The leading dimension of V just as declared in the */
/* >          calling subroutine.  NW .LE. LDV */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* >          NH is integer scalar */
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
/* >          LDT is integer */
/* >          The leading dimension of T just as declared in the */
/* >          calling subroutine.  NW .LE. LDT */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* >          NV is integer */
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
/* >          LDWV is integer */
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
/* >          LWORK is integer */
/* >          The dimension of the work array WORK.  LWORK = 2*NW */
/* >          suffices, but greater efficiency may result from larger */
/* >          values of LWORK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; DLAQR3 */
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

/* > \date June 2016 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaqr3_(logical *wantt, logical *wantz, integer *n, 
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
    static integer lwk1, lwk2, lwk3;
    static doublereal beta;
    static integer kend, kcol, info, nmin, ifst, ilst, ltop, krow;
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
	    doublereal *, doublereal *, doublereal *, doublereal *), dlaqr4_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    dlabad_(doublereal *, doublereal *);
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
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dtrexc_(char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen),
	     dormhr_(char *, char *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical sorted;
    static doublereal smlnum;
    static integer lwkopt;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Estimate optimal workspace. ==== */

#line 324 "dlaqr3.f"
    /* Parameter adjustments */
#line 324 "dlaqr3.f"
    h_dim1 = *ldh;
#line 324 "dlaqr3.f"
    h_offset = 1 + h_dim1;
#line 324 "dlaqr3.f"
    h__ -= h_offset;
#line 324 "dlaqr3.f"
    z_dim1 = *ldz;
#line 324 "dlaqr3.f"
    z_offset = 1 + z_dim1;
#line 324 "dlaqr3.f"
    z__ -= z_offset;
#line 324 "dlaqr3.f"
    --sr;
#line 324 "dlaqr3.f"
    --si;
#line 324 "dlaqr3.f"
    v_dim1 = *ldv;
#line 324 "dlaqr3.f"
    v_offset = 1 + v_dim1;
#line 324 "dlaqr3.f"
    v -= v_offset;
#line 324 "dlaqr3.f"
    t_dim1 = *ldt;
#line 324 "dlaqr3.f"
    t_offset = 1 + t_dim1;
#line 324 "dlaqr3.f"
    t -= t_offset;
#line 324 "dlaqr3.f"
    wv_dim1 = *ldwv;
#line 324 "dlaqr3.f"
    wv_offset = 1 + wv_dim1;
#line 324 "dlaqr3.f"
    wv -= wv_offset;
#line 324 "dlaqr3.f"
    --work;
#line 324 "dlaqr3.f"

#line 324 "dlaqr3.f"
    /* Function Body */
/* Computing MIN */
#line 324 "dlaqr3.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 324 "dlaqr3.f"
    jw = min(i__1,i__2);
#line 325 "dlaqr3.f"
    if (jw <= 2) {
#line 326 "dlaqr3.f"
	lwkopt = 1;
#line 327 "dlaqr3.f"
    } else {

/*        ==== Workspace query call to DGEHRD ==== */

#line 331 "dlaqr3.f"
	i__1 = jw - 1;
#line 331 "dlaqr3.f"
	dgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &
		c_n1, &info);
#line 332 "dlaqr3.f"
	lwk1 = (integer) work[1];

/*        ==== Workspace query call to DORMHR ==== */

#line 336 "dlaqr3.f"
	i__1 = jw - 1;
#line 336 "dlaqr3.f"
	dormhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
		 &v[v_offset], ldv, &work[1], &c_n1, &info, (ftnlen)1, (
		ftnlen)1);
#line 338 "dlaqr3.f"
	lwk2 = (integer) work[1];

/*        ==== Workspace query call to DLAQR4 ==== */

#line 342 "dlaqr3.f"
	dlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[1], 
		&si[1], &c__1, &jw, &v[v_offset], ldv, &work[1], &c_n1, &
		infqr);
#line 344 "dlaqr3.f"
	lwk3 = (integer) work[1];

/*        ==== Optimal workspace ==== */

/* Computing MAX */
#line 348 "dlaqr3.f"
	i__1 = jw + max(lwk1,lwk2);
#line 348 "dlaqr3.f"
	lwkopt = max(i__1,lwk3);
#line 349 "dlaqr3.f"
    }

/*     ==== Quick return in case of workspace query. ==== */

#line 353 "dlaqr3.f"
    if (*lwork == -1) {
#line 354 "dlaqr3.f"
	work[1] = (doublereal) lwkopt;
#line 355 "dlaqr3.f"
	return 0;
#line 356 "dlaqr3.f"
    }

/*     ==== Nothing to do ... */
/*     ... for an empty active block ... ==== */
#line 360 "dlaqr3.f"
    *ns = 0;
#line 361 "dlaqr3.f"
    *nd = 0;
#line 362 "dlaqr3.f"
    work[1] = 1.;
#line 363 "dlaqr3.f"
    if (*ktop > *kbot) {
#line 363 "dlaqr3.f"
	return 0;
#line 363 "dlaqr3.f"
    }
/*     ... nor for an empty deflation window. ==== */
#line 366 "dlaqr3.f"
    if (*nw < 1) {
#line 366 "dlaqr3.f"
	return 0;
#line 366 "dlaqr3.f"
    }

/*     ==== Machine constants ==== */

#line 371 "dlaqr3.f"
    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
#line 372 "dlaqr3.f"
    safmax = 1. / safmin;
#line 373 "dlaqr3.f"
    dlabad_(&safmin, &safmax);
#line 374 "dlaqr3.f"
    ulp = dlamch_("PRECISION", (ftnlen)9);
#line 375 "dlaqr3.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Setup deflation window ==== */

/* Computing MIN */
#line 379 "dlaqr3.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 379 "dlaqr3.f"
    jw = min(i__1,i__2);
#line 380 "dlaqr3.f"
    kwtop = *kbot - jw + 1;
#line 381 "dlaqr3.f"
    if (kwtop == *ktop) {
#line 382 "dlaqr3.f"
	s = 0.;
#line 383 "dlaqr3.f"
    } else {
#line 384 "dlaqr3.f"
	s = h__[kwtop + (kwtop - 1) * h_dim1];
#line 385 "dlaqr3.f"
    }

#line 387 "dlaqr3.f"
    if (*kbot == kwtop) {

/*        ==== 1-by-1 deflation window: not much to do ==== */

#line 391 "dlaqr3.f"
	sr[kwtop] = h__[kwtop + kwtop * h_dim1];
#line 392 "dlaqr3.f"
	si[kwtop] = 0.;
#line 393 "dlaqr3.f"
	*ns = 1;
#line 394 "dlaqr3.f"
	*nd = 0;
/* Computing MAX */
#line 395 "dlaqr3.f"
	d__2 = smlnum, d__3 = ulp * (d__1 = h__[kwtop + kwtop * h_dim1], abs(
		d__1));
#line 395 "dlaqr3.f"
	if (abs(s) <= max(d__2,d__3)) {
#line 397 "dlaqr3.f"
	    *ns = 0;
#line 398 "dlaqr3.f"
	    *nd = 1;
#line 399 "dlaqr3.f"
	    if (kwtop > *ktop) {
#line 399 "dlaqr3.f"
		h__[kwtop + (kwtop - 1) * h_dim1] = 0.;
#line 399 "dlaqr3.f"
	    }
#line 401 "dlaqr3.f"
	}
#line 402 "dlaqr3.f"
	work[1] = 1.;
#line 403 "dlaqr3.f"
	return 0;
#line 404 "dlaqr3.f"
    }

/*     ==== Convert to spike-triangular form.  (In case of a */
/*     .    rare QR failure, this routine continues to do */
/*     .    aggressive early deflation using that part of */
/*     .    the deflation window that converged using INFQR */
/*     .    here and there to keep track.) ==== */

#line 412 "dlaqr3.f"
    dlacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], 
	    ldt, (ftnlen)1);
#line 413 "dlaqr3.f"
    i__1 = jw - 1;
#line 413 "dlaqr3.f"
    i__2 = *ldh + 1;
#line 413 "dlaqr3.f"
    i__3 = *ldt + 1;
#line 413 "dlaqr3.f"
    dcopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &
	    i__3);

#line 415 "dlaqr3.f"
    dlaset_("A", &jw, &jw, &c_b17, &c_b18, &v[v_offset], ldv, (ftnlen)1);
#line 416 "dlaqr3.f"
    nmin = ilaenv_(&c__12, "DLAQR3", "SV", &jw, &c__1, &jw, lwork, (ftnlen)6, 
	    (ftnlen)2);
#line 417 "dlaqr3.f"
    if (jw > nmin) {
#line 418 "dlaqr3.f"
	dlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[
		kwtop], &si[kwtop], &c__1, &jw, &v[v_offset], ldv, &work[1], 
		lwork, &infqr);
#line 420 "dlaqr3.f"
    } else {
#line 421 "dlaqr3.f"
	dlahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[
		kwtop], &si[kwtop], &c__1, &jw, &v[v_offset], ldv, &infqr);
#line 423 "dlaqr3.f"
    }

/*     ==== DTREXC needs a clean margin near the diagonal ==== */

#line 427 "dlaqr3.f"
    i__1 = jw - 3;
#line 427 "dlaqr3.f"
    for (j = 1; j <= i__1; ++j) {
#line 428 "dlaqr3.f"
	t[j + 2 + j * t_dim1] = 0.;
#line 429 "dlaqr3.f"
	t[j + 3 + j * t_dim1] = 0.;
#line 430 "dlaqr3.f"
/* L10: */
#line 430 "dlaqr3.f"
    }
#line 431 "dlaqr3.f"
    if (jw > 2) {
#line 431 "dlaqr3.f"
	t[jw + (jw - 2) * t_dim1] = 0.;
#line 431 "dlaqr3.f"
    }

/*     ==== Deflation detection loop ==== */

#line 436 "dlaqr3.f"
    *ns = jw;
#line 437 "dlaqr3.f"
    ilst = infqr + 1;
#line 438 "dlaqr3.f"
L20:
#line 439 "dlaqr3.f"
    if (ilst <= *ns) {
#line 440 "dlaqr3.f"
	if (*ns == 1) {
#line 441 "dlaqr3.f"
	    bulge = FALSE_;
#line 442 "dlaqr3.f"
	} else {
#line 443 "dlaqr3.f"
	    bulge = t[*ns + (*ns - 1) * t_dim1] != 0.;
#line 444 "dlaqr3.f"
	}

/*        ==== Small spike tip test for deflation ==== */

#line 448 "dlaqr3.f"
	if (! bulge) {

/*           ==== Real eigenvalue ==== */

#line 452 "dlaqr3.f"
	    foo = (d__1 = t[*ns + *ns * t_dim1], abs(d__1));
#line 453 "dlaqr3.f"
	    if (foo == 0.) {
#line 453 "dlaqr3.f"
		foo = abs(s);
#line 453 "dlaqr3.f"
	    }
/* Computing MAX */
#line 455 "dlaqr3.f"
	    d__2 = smlnum, d__3 = ulp * foo;
#line 455 "dlaqr3.f"
	    if ((d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)) <= max(d__2,d__3))
		     {

/*              ==== Deflatable ==== */

#line 459 "dlaqr3.f"
		--(*ns);
#line 460 "dlaqr3.f"
	    } else {

/*              ==== Undeflatable.   Move it up out of the way. */
/*              .    (DTREXC can not fail in this case.) ==== */

#line 465 "dlaqr3.f"
		ifst = *ns;
#line 466 "dlaqr3.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 468 "dlaqr3.f"
		++ilst;
#line 469 "dlaqr3.f"
	    }
#line 470 "dlaqr3.f"
	} else {

/*           ==== Complex conjugate pair ==== */

#line 474 "dlaqr3.f"
	    foo = (d__3 = t[*ns + *ns * t_dim1], abs(d__3)) + sqrt((d__1 = t[*
		    ns + (*ns - 1) * t_dim1], abs(d__1))) * sqrt((d__2 = t[*
		    ns - 1 + *ns * t_dim1], abs(d__2)));
#line 476 "dlaqr3.f"
	    if (foo == 0.) {
#line 476 "dlaqr3.f"
		foo = abs(s);
#line 476 "dlaqr3.f"
	    }
/* Computing MAX */
#line 478 "dlaqr3.f"
	    d__3 = (d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)), d__4 = (d__2 =
		     s * v[(*ns - 1) * v_dim1 + 1], abs(d__2));
/* Computing MAX */
#line 478 "dlaqr3.f"
	    d__5 = smlnum, d__6 = ulp * foo;
#line 478 "dlaqr3.f"
	    if (max(d__3,d__4) <= max(d__5,d__6)) {

/*              ==== Deflatable ==== */

#line 483 "dlaqr3.f"
		*ns += -2;
#line 484 "dlaqr3.f"
	    } else {

/*              ==== Undeflatable. Move them up out of the way. */
/*              .    Fortunately, DTREXC does the right thing with */
/*              .    ILST in case of a rare exchange failure. ==== */

#line 490 "dlaqr3.f"
		ifst = *ns;
#line 491 "dlaqr3.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 493 "dlaqr3.f"
		ilst += 2;
#line 494 "dlaqr3.f"
	    }
#line 495 "dlaqr3.f"
	}

/*        ==== End deflation detection loop ==== */

#line 499 "dlaqr3.f"
	goto L20;
#line 500 "dlaqr3.f"
    }

/*        ==== Return to Hessenberg form ==== */

#line 504 "dlaqr3.f"
    if (*ns == 0) {
#line 504 "dlaqr3.f"
	s = 0.;
#line 504 "dlaqr3.f"
    }

#line 507 "dlaqr3.f"
    if (*ns < jw) {

/*        ==== sorting diagonal blocks of T improves accuracy for */
/*        .    graded matrices.  Bubble sort deals well with */
/*        .    exchange failures. ==== */

#line 513 "dlaqr3.f"
	sorted = FALSE_;
#line 514 "dlaqr3.f"
	i__ = *ns + 1;
#line 515 "dlaqr3.f"
L30:
#line 516 "dlaqr3.f"
	if (sorted) {
#line 516 "dlaqr3.f"
	    goto L50;
#line 516 "dlaqr3.f"
	}
#line 518 "dlaqr3.f"
	sorted = TRUE_;

#line 520 "dlaqr3.f"
	kend = i__ - 1;
#line 521 "dlaqr3.f"
	i__ = infqr + 1;
#line 522 "dlaqr3.f"
	if (i__ == *ns) {
#line 523 "dlaqr3.f"
	    k = i__ + 1;
#line 524 "dlaqr3.f"
	} else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
#line 525 "dlaqr3.f"
	    k = i__ + 1;
#line 526 "dlaqr3.f"
	} else {
#line 527 "dlaqr3.f"
	    k = i__ + 2;
#line 528 "dlaqr3.f"
	}
#line 529 "dlaqr3.f"
L40:
#line 530 "dlaqr3.f"
	if (k <= kend) {
#line 531 "dlaqr3.f"
	    if (k == i__ + 1) {
#line 532 "dlaqr3.f"
		evi = (d__1 = t[i__ + i__ * t_dim1], abs(d__1));
#line 533 "dlaqr3.f"
	    } else {
#line 534 "dlaqr3.f"
		evi = (d__3 = t[i__ + i__ * t_dim1], abs(d__3)) + sqrt((d__1 =
			 t[i__ + 1 + i__ * t_dim1], abs(d__1))) * sqrt((d__2 =
			 t[i__ + (i__ + 1) * t_dim1], abs(d__2)));
#line 536 "dlaqr3.f"
	    }

#line 538 "dlaqr3.f"
	    if (k == kend) {
#line 539 "dlaqr3.f"
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
#line 540 "dlaqr3.f"
	    } else if (t[k + 1 + k * t_dim1] == 0.) {
#line 541 "dlaqr3.f"
		evk = (d__1 = t[k + k * t_dim1], abs(d__1));
#line 542 "dlaqr3.f"
	    } else {
#line 543 "dlaqr3.f"
		evk = (d__3 = t[k + k * t_dim1], abs(d__3)) + sqrt((d__1 = t[
			k + 1 + k * t_dim1], abs(d__1))) * sqrt((d__2 = t[k + 
			(k + 1) * t_dim1], abs(d__2)));
#line 545 "dlaqr3.f"
	    }

#line 547 "dlaqr3.f"
	    if (evi >= evk) {
#line 548 "dlaqr3.f"
		i__ = k;
#line 549 "dlaqr3.f"
	    } else {
#line 550 "dlaqr3.f"
		sorted = FALSE_;
#line 551 "dlaqr3.f"
		ifst = i__;
#line 552 "dlaqr3.f"
		ilst = k;
#line 553 "dlaqr3.f"
		dtrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &work[1], &info, (ftnlen)1);
#line 555 "dlaqr3.f"
		if (info == 0) {
#line 556 "dlaqr3.f"
		    i__ = ilst;
#line 557 "dlaqr3.f"
		} else {
#line 558 "dlaqr3.f"
		    i__ = k;
#line 559 "dlaqr3.f"
		}
#line 560 "dlaqr3.f"
	    }
#line 561 "dlaqr3.f"
	    if (i__ == kend) {
#line 562 "dlaqr3.f"
		k = i__ + 1;
#line 563 "dlaqr3.f"
	    } else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
#line 564 "dlaqr3.f"
		k = i__ + 1;
#line 565 "dlaqr3.f"
	    } else {
#line 566 "dlaqr3.f"
		k = i__ + 2;
#line 567 "dlaqr3.f"
	    }
#line 568 "dlaqr3.f"
	    goto L40;
#line 569 "dlaqr3.f"
	}
#line 570 "dlaqr3.f"
	goto L30;
#line 571 "dlaqr3.f"
L50:
#line 572 "dlaqr3.f"
	;
#line 572 "dlaqr3.f"
    }

/*     ==== Restore shift/eigenvalue array from T ==== */

#line 576 "dlaqr3.f"
    i__ = jw;
#line 577 "dlaqr3.f"
L60:
#line 578 "dlaqr3.f"
    if (i__ >= infqr + 1) {
#line 579 "dlaqr3.f"
	if (i__ == infqr + 1) {
#line 580 "dlaqr3.f"
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
#line 581 "dlaqr3.f"
	    si[kwtop + i__ - 1] = 0.;
#line 582 "dlaqr3.f"
	    --i__;
#line 583 "dlaqr3.f"
	} else if (t[i__ + (i__ - 1) * t_dim1] == 0.) {
#line 584 "dlaqr3.f"
	    sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
#line 585 "dlaqr3.f"
	    si[kwtop + i__ - 1] = 0.;
#line 586 "dlaqr3.f"
	    --i__;
#line 587 "dlaqr3.f"
	} else {
#line 588 "dlaqr3.f"
	    aa = t[i__ - 1 + (i__ - 1) * t_dim1];
#line 589 "dlaqr3.f"
	    cc = t[i__ + (i__ - 1) * t_dim1];
#line 590 "dlaqr3.f"
	    bb = t[i__ - 1 + i__ * t_dim1];
#line 591 "dlaqr3.f"
	    dd = t[i__ + i__ * t_dim1];
#line 592 "dlaqr3.f"
	    dlanv2_(&aa, &bb, &cc, &dd, &sr[kwtop + i__ - 2], &si[kwtop + i__ 
		    - 2], &sr[kwtop + i__ - 1], &si[kwtop + i__ - 1], &cs, &
		    sn);
#line 595 "dlaqr3.f"
	    i__ += -2;
#line 596 "dlaqr3.f"
	}
#line 597 "dlaqr3.f"
	goto L60;
#line 598 "dlaqr3.f"
    }

#line 600 "dlaqr3.f"
    if (*ns < jw || s == 0.) {
#line 601 "dlaqr3.f"
	if (*ns > 1 && s != 0.) {

/*           ==== Reflect spike back into lower triangle ==== */

#line 605 "dlaqr3.f"
	    dcopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
#line 606 "dlaqr3.f"
	    beta = work[1];
#line 607 "dlaqr3.f"
	    dlarfg_(ns, &beta, &work[2], &c__1, &tau);
#line 608 "dlaqr3.f"
	    work[1] = 1.;

#line 610 "dlaqr3.f"
	    i__1 = jw - 2;
#line 610 "dlaqr3.f"
	    i__2 = jw - 2;
#line 610 "dlaqr3.f"
	    dlaset_("L", &i__1, &i__2, &c_b17, &c_b17, &t[t_dim1 + 3], ldt, (
		    ftnlen)1);

#line 612 "dlaqr3.f"
	    dlarf_("L", ns, &jw, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 614 "dlaqr3.f"
	    dlarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 616 "dlaqr3.f"
	    dlarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &
		    work[jw + 1], (ftnlen)1);

#line 619 "dlaqr3.f"
	    i__1 = *lwork - jw;
#line 619 "dlaqr3.f"
	    dgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1]
		    , &i__1, &info);
#line 621 "dlaqr3.f"
	}

/*        ==== Copy updated reduced window into place ==== */

#line 625 "dlaqr3.f"
	if (kwtop > 1) {
#line 625 "dlaqr3.f"
	    h__[kwtop + (kwtop - 1) * h_dim1] = s * v[v_dim1 + 1];
#line 625 "dlaqr3.f"
	}
#line 627 "dlaqr3.f"
	dlacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1]
		, ldh, (ftnlen)1);
#line 628 "dlaqr3.f"
	i__1 = jw - 1;
#line 628 "dlaqr3.f"
	i__2 = *ldt + 1;
#line 628 "dlaqr3.f"
	i__3 = *ldh + 1;
#line 628 "dlaqr3.f"
	dcopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1],
		 &i__3);

/*        ==== Accumulate orthogonal matrix in order update */
/*        .    H and Z, if requested.  ==== */

#line 634 "dlaqr3.f"
	if (*ns > 1 && s != 0.) {
#line 634 "dlaqr3.f"
	    i__1 = *lwork - jw;
#line 634 "dlaqr3.f"
	    dormhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
		     &v[v_offset], ldv, &work[jw + 1], &i__1, &info, (ftnlen)
		    1, (ftnlen)1);
#line 634 "dlaqr3.f"
	}

/*        ==== Update vertical slab in H ==== */

#line 640 "dlaqr3.f"
	if (*wantt) {
#line 641 "dlaqr3.f"
	    ltop = 1;
#line 642 "dlaqr3.f"
	} else {
#line 643 "dlaqr3.f"
	    ltop = *ktop;
#line 644 "dlaqr3.f"
	}
#line 645 "dlaqr3.f"
	i__1 = kwtop - 1;
#line 645 "dlaqr3.f"
	i__2 = *nv;
#line 645 "dlaqr3.f"
	for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += 
		i__2) {
/* Computing MIN */
#line 646 "dlaqr3.f"
	    i__3 = *nv, i__4 = kwtop - krow;
#line 646 "dlaqr3.f"
	    kln = min(i__3,i__4);
#line 647 "dlaqr3.f"
	    dgemm_("N", "N", &kln, &jw, &jw, &c_b18, &h__[krow + kwtop * 
		    h_dim1], ldh, &v[v_offset], ldv, &c_b17, &wv[wv_offset], 
		    ldwv, (ftnlen)1, (ftnlen)1);
#line 649 "dlaqr3.f"
	    dlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * 
		    h_dim1], ldh, (ftnlen)1);
#line 650 "dlaqr3.f"
/* L70: */
#line 650 "dlaqr3.f"
	}

/*        ==== Update horizontal slab in H ==== */

#line 654 "dlaqr3.f"
	if (*wantt) {
#line 655 "dlaqr3.f"
	    i__2 = *n;
#line 655 "dlaqr3.f"
	    i__1 = *nh;
#line 655 "dlaqr3.f"
	    for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; 
		    kcol += i__1) {
/* Computing MIN */
#line 656 "dlaqr3.f"
		i__3 = *nh, i__4 = *n - kcol + 1;
#line 656 "dlaqr3.f"
		kln = min(i__3,i__4);
#line 657 "dlaqr3.f"
		dgemm_("C", "N", &jw, &kln, &jw, &c_b18, &v[v_offset], ldv, &
			h__[kwtop + kcol * h_dim1], ldh, &c_b17, &t[t_offset],
			 ldt, (ftnlen)1, (ftnlen)1);
#line 659 "dlaqr3.f"
		dlacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol *
			 h_dim1], ldh, (ftnlen)1);
#line 661 "dlaqr3.f"
/* L80: */
#line 661 "dlaqr3.f"
	    }
#line 662 "dlaqr3.f"
	}

/*        ==== Update vertical slab in Z ==== */

#line 666 "dlaqr3.f"
	if (*wantz) {
#line 667 "dlaqr3.f"
	    i__1 = *ihiz;
#line 667 "dlaqr3.f"
	    i__2 = *nv;
#line 667 "dlaqr3.f"
	    for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow +=
		     i__2) {
/* Computing MIN */
#line 668 "dlaqr3.f"
		i__3 = *nv, i__4 = *ihiz - krow + 1;
#line 668 "dlaqr3.f"
		kln = min(i__3,i__4);
#line 669 "dlaqr3.f"
		dgemm_("N", "N", &kln, &jw, &jw, &c_b18, &z__[krow + kwtop * 
			z_dim1], ldz, &v[v_offset], ldv, &c_b17, &wv[
			wv_offset], ldwv, (ftnlen)1, (ftnlen)1);
#line 671 "dlaqr3.f"
		dlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + 
			kwtop * z_dim1], ldz, (ftnlen)1);
#line 673 "dlaqr3.f"
/* L90: */
#line 673 "dlaqr3.f"
	    }
#line 674 "dlaqr3.f"
	}
#line 675 "dlaqr3.f"
    }

/*     ==== Return the number of deflations ... ==== */

#line 679 "dlaqr3.f"
    *nd = jw - *ns;

/*     ==== ... and the number of shifts. (Subtracting */
/*     .    INFQR from the spike length takes care */
/*     .    of the case of a rare QR failure while */
/*     .    calculating eigenvalues of the deflation */
/*     .    window.)  ==== */

#line 687 "dlaqr3.f"
    *ns -= infqr;

/*      ==== Return optimal workspace. ==== */

#line 691 "dlaqr3.f"
    work[1] = (doublereal) lwkopt;

/*     ==== End of DLAQR3 ==== */

#line 695 "dlaqr3.f"
    return 0;
} /* dlaqr3_ */

