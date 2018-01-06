#line 1 "zlaqr3.f"
/* zlaqr3.f -- translated by f2c (version 20100827).
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

#line 1 "zlaqr3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
static integer c__12 = 12;

/* > \brief \b ZLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and defl
ate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAQR3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
/*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, */
/*                          NV, WV, LDWV, WORK, LWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, */
/*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    Aggressive early deflation: */
/* > */
/* >    ZLAQR3 accepts as input an upper Hessenberg matrix */
/* >    H and performs an unitary similarity transformation */
/* >    designed to detect and deflate fully converged eigenvalues from */
/* >    a trailing principal submatrix.  On output H has been over- */
/* >    written by a new Hessenberg matrix that is a perturbation of */
/* >    an unitary similarity transformation of H.  It is to be */
/* >    hoped that the final version of H has many zero subdiagonal */
/* >    entries. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTT */
/* > \verbatim */
/* >          WANTT is LOGICAL */
/* >          If .TRUE., then the Hessenberg matrix H is fully updated */
/* >          so that the triangular Schur factor may be */
/* >          computed (in cooperation with the calling subroutine). */
/* >          If .FALSE., then only enough of H is updated to preserve */
/* >          the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          If .TRUE., then the unitary matrix Z is updated so */
/* >          so that the unitary Schur factor may be computed */
/* >          (in cooperation with the calling subroutine). */
/* >          If .FALSE., then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H and (if WANTZ is .TRUE.) the */
/* >          order of the unitary matrix Z. */
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
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
/* >          On input the initial N-by-N section of H stores the */
/* >          Hessenberg matrix undergoing aggressive early deflation. */
/* >          On output H has been transformed by a unitary */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ,N) */
/* >          IF WANTZ is .TRUE., then on output, the unitary */
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
/* > \param[out] SH */
/* > \verbatim */
/* >          SH is COMPLEX*16 array, dimension (KBOT) */
/* >          On output, approximate eigenvalues that may */
/* >          be used for shifts are stored in SH(KBOT-ND-NS+1) */
/* >          through SR(KBOT-ND).  Converged eigenvalues are */
/* >          stored in SH(KBOT-ND+1) through SH(KBOT). */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,NW) */
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
/* >          T is COMPLEX*16 array, dimension (LDT,NW) */
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
/* >          WV is COMPLEX*16 array, dimension (LDWV,NW) */
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
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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
/* >          If LWORK = -1, then a workspace query is assumed; ZLAQR3 */
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

/* > \ingroup complex16OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zlaqr3_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, doublecomplex *h__, 
	integer *ldh, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *ns, integer *nd, doublecomplex *sh, 
	doublecomplex *v, integer *ldv, integer *nh, doublecomplex *t, 
	integer *ldt, integer *nv, doublecomplex *wv, integer *ldwv, 
	doublecomplex *work, integer *lwork)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, 
	    wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex s;
    static integer jw;
    static doublereal foo;
    static integer kln;
    static doublecomplex tau;
    static integer knt;
    static doublereal ulp;
    static integer lwk1, lwk2, lwk3;
    static doublecomplex beta;
    static integer kcol, info, nmin, ifst, ilst, ltop, krow;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static integer infqr;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer kwtop;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), dlabad_(doublereal *, doublereal *), 
	    zlaqr4_(logical *, logical *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    );
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal safmax;
    extern /* Subroutine */ int zgehrd_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zlahqr_(logical *, 
	    logical *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, doublecomplex *,
	     integer *, integer *), zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int ztrexc_(char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, integer *, 
	    integer *, ftnlen);
    static integer lwkopt;
    extern /* Subroutine */ int zunmhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Estimate optimal workspace. ==== */

#line 323 "zlaqr3.f"
    /* Parameter adjustments */
#line 323 "zlaqr3.f"
    h_dim1 = *ldh;
#line 323 "zlaqr3.f"
    h_offset = 1 + h_dim1;
#line 323 "zlaqr3.f"
    h__ -= h_offset;
#line 323 "zlaqr3.f"
    z_dim1 = *ldz;
#line 323 "zlaqr3.f"
    z_offset = 1 + z_dim1;
#line 323 "zlaqr3.f"
    z__ -= z_offset;
#line 323 "zlaqr3.f"
    --sh;
#line 323 "zlaqr3.f"
    v_dim1 = *ldv;
#line 323 "zlaqr3.f"
    v_offset = 1 + v_dim1;
#line 323 "zlaqr3.f"
    v -= v_offset;
#line 323 "zlaqr3.f"
    t_dim1 = *ldt;
#line 323 "zlaqr3.f"
    t_offset = 1 + t_dim1;
#line 323 "zlaqr3.f"
    t -= t_offset;
#line 323 "zlaqr3.f"
    wv_dim1 = *ldwv;
#line 323 "zlaqr3.f"
    wv_offset = 1 + wv_dim1;
#line 323 "zlaqr3.f"
    wv -= wv_offset;
#line 323 "zlaqr3.f"
    --work;
#line 323 "zlaqr3.f"

#line 323 "zlaqr3.f"
    /* Function Body */
/* Computing MIN */
#line 323 "zlaqr3.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 323 "zlaqr3.f"
    jw = min(i__1,i__2);
#line 324 "zlaqr3.f"
    if (jw <= 2) {
#line 325 "zlaqr3.f"
	lwkopt = 1;
#line 326 "zlaqr3.f"
    } else {

/*        ==== Workspace query call to ZGEHRD ==== */

#line 330 "zlaqr3.f"
	i__1 = jw - 1;
#line 330 "zlaqr3.f"
	zgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &
		c_n1, &info);
#line 331 "zlaqr3.f"
	lwk1 = (integer) work[1].r;

/*        ==== Workspace query call to ZUNMHR ==== */

#line 335 "zlaqr3.f"
	i__1 = jw - 1;
#line 335 "zlaqr3.f"
	zunmhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
		 &v[v_offset], ldv, &work[1], &c_n1, &info, (ftnlen)1, (
		ftnlen)1);
#line 337 "zlaqr3.f"
	lwk2 = (integer) work[1].r;

/*        ==== Workspace query call to ZLAQR4 ==== */

#line 341 "zlaqr3.f"
	zlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[1], 
		&c__1, &jw, &v[v_offset], ldv, &work[1], &c_n1, &infqr);
#line 343 "zlaqr3.f"
	lwk3 = (integer) work[1].r;

/*        ==== Optimal workspace ==== */

/* Computing MAX */
#line 347 "zlaqr3.f"
	i__1 = jw + max(lwk1,lwk2);
#line 347 "zlaqr3.f"
	lwkopt = max(i__1,lwk3);
#line 348 "zlaqr3.f"
    }

/*     ==== Quick return in case of workspace query. ==== */

#line 352 "zlaqr3.f"
    if (*lwork == -1) {
#line 353 "zlaqr3.f"
	d__1 = (doublereal) lwkopt;
#line 353 "zlaqr3.f"
	z__1.r = d__1, z__1.i = 0.;
#line 353 "zlaqr3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 354 "zlaqr3.f"
	return 0;
#line 355 "zlaqr3.f"
    }

/*     ==== Nothing to do ... */
/*     ... for an empty active block ... ==== */
#line 359 "zlaqr3.f"
    *ns = 0;
#line 360 "zlaqr3.f"
    *nd = 0;
#line 361 "zlaqr3.f"
    work[1].r = 1., work[1].i = 0.;
#line 362 "zlaqr3.f"
    if (*ktop > *kbot) {
#line 362 "zlaqr3.f"
	return 0;
#line 362 "zlaqr3.f"
    }
/*     ... nor for an empty deflation window. ==== */
#line 365 "zlaqr3.f"
    if (*nw < 1) {
#line 365 "zlaqr3.f"
	return 0;
#line 365 "zlaqr3.f"
    }

/*     ==== Machine constants ==== */

#line 370 "zlaqr3.f"
    safmin = dlamch_("SAFE MINIMUM", (ftnlen)12);
#line 371 "zlaqr3.f"
    safmax = 1. / safmin;
#line 372 "zlaqr3.f"
    dlabad_(&safmin, &safmax);
#line 373 "zlaqr3.f"
    ulp = dlamch_("PRECISION", (ftnlen)9);
#line 374 "zlaqr3.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Setup deflation window ==== */

/* Computing MIN */
#line 378 "zlaqr3.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 378 "zlaqr3.f"
    jw = min(i__1,i__2);
#line 379 "zlaqr3.f"
    kwtop = *kbot - jw + 1;
#line 380 "zlaqr3.f"
    if (kwtop == *ktop) {
#line 381 "zlaqr3.f"
	s.r = 0., s.i = 0.;
#line 382 "zlaqr3.f"
    } else {
#line 383 "zlaqr3.f"
	i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 383 "zlaqr3.f"
	s.r = h__[i__1].r, s.i = h__[i__1].i;
#line 384 "zlaqr3.f"
    }

#line 386 "zlaqr3.f"
    if (*kbot == kwtop) {

/*        ==== 1-by-1 deflation window: not much to do ==== */

#line 390 "zlaqr3.f"
	i__1 = kwtop;
#line 390 "zlaqr3.f"
	i__2 = kwtop + kwtop * h_dim1;
#line 390 "zlaqr3.f"
	sh[i__1].r = h__[i__2].r, sh[i__1].i = h__[i__2].i;
#line 391 "zlaqr3.f"
	*ns = 1;
#line 392 "zlaqr3.f"
	*nd = 0;
/* Computing MAX */
#line 393 "zlaqr3.f"
	i__1 = kwtop + kwtop * h_dim1;
#line 393 "zlaqr3.f"
	d__5 = smlnum, d__6 = ulp * ((d__1 = h__[i__1].r, abs(d__1)) + (d__2 =
		 d_imag(&h__[kwtop + kwtop * h_dim1]), abs(d__2)));
#line 393 "zlaqr3.f"
	if ((d__3 = s.r, abs(d__3)) + (d__4 = d_imag(&s), abs(d__4)) <= max(
		d__5,d__6)) {
#line 395 "zlaqr3.f"
	    *ns = 0;
#line 396 "zlaqr3.f"
	    *nd = 1;
#line 397 "zlaqr3.f"
	    if (kwtop > *ktop) {
#line 397 "zlaqr3.f"
		i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 397 "zlaqr3.f"
		h__[i__1].r = 0., h__[i__1].i = 0.;
#line 397 "zlaqr3.f"
	    }
#line 399 "zlaqr3.f"
	}
#line 400 "zlaqr3.f"
	work[1].r = 1., work[1].i = 0.;
#line 401 "zlaqr3.f"
	return 0;
#line 402 "zlaqr3.f"
    }

/*     ==== Convert to spike-triangular form.  (In case of a */
/*     .    rare QR failure, this routine continues to do */
/*     .    aggressive early deflation using that part of */
/*     .    the deflation window that converged using INFQR */
/*     .    here and there to keep track.) ==== */

#line 410 "zlaqr3.f"
    zlacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], 
	    ldt, (ftnlen)1);
#line 411 "zlaqr3.f"
    i__1 = jw - 1;
#line 411 "zlaqr3.f"
    i__2 = *ldh + 1;
#line 411 "zlaqr3.f"
    i__3 = *ldt + 1;
#line 411 "zlaqr3.f"
    zcopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &
	    i__3);

#line 413 "zlaqr3.f"
    zlaset_("A", &jw, &jw, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 414 "zlaqr3.f"
    nmin = ilaenv_(&c__12, "ZLAQR3", "SV", &jw, &c__1, &jw, lwork, (ftnlen)6, 
	    (ftnlen)2);
#line 415 "zlaqr3.f"
    if (jw > nmin) {
#line 416 "zlaqr3.f"
	zlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[
		kwtop], &c__1, &jw, &v[v_offset], ldv, &work[1], lwork, &
		infqr);
#line 418 "zlaqr3.f"
    } else {
#line 419 "zlaqr3.f"
	zlahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[
		kwtop], &c__1, &jw, &v[v_offset], ldv, &infqr);
#line 421 "zlaqr3.f"
    }

/*     ==== Deflation detection loop ==== */

#line 425 "zlaqr3.f"
    *ns = jw;
#line 426 "zlaqr3.f"
    ilst = infqr + 1;
#line 427 "zlaqr3.f"
    i__1 = jw;
#line 427 "zlaqr3.f"
    for (knt = infqr + 1; knt <= i__1; ++knt) {

/*        ==== Small spike tip deflation test ==== */

#line 431 "zlaqr3.f"
	i__2 = *ns + *ns * t_dim1;
#line 431 "zlaqr3.f"
	foo = (d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[*ns + *ns * 
		t_dim1]), abs(d__2));
#line 432 "zlaqr3.f"
	if (foo == 0.) {
#line 432 "zlaqr3.f"
	    foo = (d__1 = s.r, abs(d__1)) + (d__2 = d_imag(&s), abs(d__2));
#line 432 "zlaqr3.f"
	}
#line 434 "zlaqr3.f"
	i__2 = *ns * v_dim1 + 1;
/* Computing MAX */
#line 434 "zlaqr3.f"
	d__5 = smlnum, d__6 = ulp * foo;
#line 434 "zlaqr3.f"
	if (((d__1 = s.r, abs(d__1)) + (d__2 = d_imag(&s), abs(d__2))) * ((
		d__3 = v[i__2].r, abs(d__3)) + (d__4 = d_imag(&v[*ns * v_dim1 
		+ 1]), abs(d__4))) <= max(d__5,d__6)) {

/*           ==== One more converged eigenvalue ==== */

#line 439 "zlaqr3.f"
	    --(*ns);
#line 440 "zlaqr3.f"
	} else {

/*           ==== One undeflatable eigenvalue.  Move it up out of the */
/*           .    way.   (ZTREXC can not fail in this case.) ==== */

#line 445 "zlaqr3.f"
	    ifst = *ns;
#line 446 "zlaqr3.f"
	    ztrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &
		    ilst, &info, (ftnlen)1);
#line 447 "zlaqr3.f"
	    ++ilst;
#line 448 "zlaqr3.f"
	}
#line 449 "zlaqr3.f"
/* L10: */
#line 449 "zlaqr3.f"
    }

/*        ==== Return to Hessenberg form ==== */

#line 453 "zlaqr3.f"
    if (*ns == 0) {
#line 453 "zlaqr3.f"
	s.r = 0., s.i = 0.;
#line 453 "zlaqr3.f"
    }

#line 456 "zlaqr3.f"
    if (*ns < jw) {

/*        ==== sorting the diagonal of T improves accuracy for */
/*        .    graded matrices.  ==== */

#line 461 "zlaqr3.f"
	i__1 = *ns;
#line 461 "zlaqr3.f"
	for (i__ = infqr + 1; i__ <= i__1; ++i__) {
#line 462 "zlaqr3.f"
	    ifst = i__;
#line 463 "zlaqr3.f"
	    i__2 = *ns;
#line 463 "zlaqr3.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 464 "zlaqr3.f"
		i__3 = j + j * t_dim1;
#line 464 "zlaqr3.f"
		i__4 = ifst + ifst * t_dim1;
#line 464 "zlaqr3.f"
		if ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[j + j * 
			t_dim1]), abs(d__2)) > (d__3 = t[i__4].r, abs(d__3)) 
			+ (d__4 = d_imag(&t[ifst + ifst * t_dim1]), abs(d__4))
			) {
#line 464 "zlaqr3.f"
		    ifst = j;
#line 464 "zlaqr3.f"
		}
#line 466 "zlaqr3.f"
/* L20: */
#line 466 "zlaqr3.f"
	    }
#line 467 "zlaqr3.f"
	    ilst = i__;
#line 468 "zlaqr3.f"
	    if (ifst != ilst) {
#line 468 "zlaqr3.f"
		ztrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &info, (ftnlen)1);
#line 468 "zlaqr3.f"
	    }
#line 470 "zlaqr3.f"
/* L30: */
#line 470 "zlaqr3.f"
	}
#line 471 "zlaqr3.f"
    }

/*     ==== Restore shift/eigenvalue array from T ==== */

#line 475 "zlaqr3.f"
    i__1 = jw;
#line 475 "zlaqr3.f"
    for (i__ = infqr + 1; i__ <= i__1; ++i__) {
#line 476 "zlaqr3.f"
	i__2 = kwtop + i__ - 1;
#line 476 "zlaqr3.f"
	i__3 = i__ + i__ * t_dim1;
#line 476 "zlaqr3.f"
	sh[i__2].r = t[i__3].r, sh[i__2].i = t[i__3].i;
#line 477 "zlaqr3.f"
/* L40: */
#line 477 "zlaqr3.f"
    }


#line 480 "zlaqr3.f"
    if (*ns < jw || s.r == 0. && s.i == 0.) {
#line 481 "zlaqr3.f"
	if (*ns > 1 && (s.r != 0. || s.i != 0.)) {

/*           ==== Reflect spike back into lower triangle ==== */

#line 485 "zlaqr3.f"
	    zcopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
#line 486 "zlaqr3.f"
	    i__1 = *ns;
#line 486 "zlaqr3.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 487 "zlaqr3.f"
		i__2 = i__;
#line 487 "zlaqr3.f"
		d_cnjg(&z__1, &work[i__]);
#line 487 "zlaqr3.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 488 "zlaqr3.f"
/* L50: */
#line 488 "zlaqr3.f"
	    }
#line 489 "zlaqr3.f"
	    beta.r = work[1].r, beta.i = work[1].i;
#line 490 "zlaqr3.f"
	    zlarfg_(ns, &beta, &work[2], &c__1, &tau);
#line 491 "zlaqr3.f"
	    work[1].r = 1., work[1].i = 0.;

#line 493 "zlaqr3.f"
	    i__1 = jw - 2;
#line 493 "zlaqr3.f"
	    i__2 = jw - 2;
#line 493 "zlaqr3.f"
	    zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &t[t_dim1 + 3], ldt, (
		    ftnlen)1);

#line 495 "zlaqr3.f"
	    d_cnjg(&z__1, &tau);
#line 495 "zlaqr3.f"
	    zlarf_("L", ns, &jw, &work[1], &c__1, &z__1, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 497 "zlaqr3.f"
	    zlarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 499 "zlaqr3.f"
	    zlarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &
		    work[jw + 1], (ftnlen)1);

#line 502 "zlaqr3.f"
	    i__1 = *lwork - jw;
#line 502 "zlaqr3.f"
	    zgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1]
		    , &i__1, &info);
#line 504 "zlaqr3.f"
	}

/*        ==== Copy updated reduced window into place ==== */

#line 508 "zlaqr3.f"
	if (kwtop > 1) {
#line 508 "zlaqr3.f"
	    i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 508 "zlaqr3.f"
	    d_cnjg(&z__2, &v[v_dim1 + 1]);
#line 508 "zlaqr3.f"
	    z__1.r = s.r * z__2.r - s.i * z__2.i, z__1.i = s.r * z__2.i + s.i 
		    * z__2.r;
#line 508 "zlaqr3.f"
	    h__[i__1].r = z__1.r, h__[i__1].i = z__1.i;
#line 508 "zlaqr3.f"
	}
#line 510 "zlaqr3.f"
	zlacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1]
		, ldh, (ftnlen)1);
#line 511 "zlaqr3.f"
	i__1 = jw - 1;
#line 511 "zlaqr3.f"
	i__2 = *ldt + 1;
#line 511 "zlaqr3.f"
	i__3 = *ldh + 1;
#line 511 "zlaqr3.f"
	zcopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1],
		 &i__3);

/*        ==== Accumulate orthogonal matrix in order update */
/*        .    H and Z, if requested.  ==== */

#line 517 "zlaqr3.f"
	if (*ns > 1 && (s.r != 0. || s.i != 0.)) {
#line 517 "zlaqr3.f"
	    i__1 = *lwork - jw;
#line 517 "zlaqr3.f"
	    zunmhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
		     &v[v_offset], ldv, &work[jw + 1], &i__1, &info, (ftnlen)
		    1, (ftnlen)1);
#line 517 "zlaqr3.f"
	}

/*        ==== Update vertical slab in H ==== */

#line 523 "zlaqr3.f"
	if (*wantt) {
#line 524 "zlaqr3.f"
	    ltop = 1;
#line 525 "zlaqr3.f"
	} else {
#line 526 "zlaqr3.f"
	    ltop = *ktop;
#line 527 "zlaqr3.f"
	}
#line 528 "zlaqr3.f"
	i__1 = kwtop - 1;
#line 528 "zlaqr3.f"
	i__2 = *nv;
#line 528 "zlaqr3.f"
	for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += 
		i__2) {
/* Computing MIN */
#line 529 "zlaqr3.f"
	    i__3 = *nv, i__4 = kwtop - krow;
#line 529 "zlaqr3.f"
	    kln = min(i__3,i__4);
#line 530 "zlaqr3.f"
	    zgemm_("N", "N", &kln, &jw, &jw, &c_b2, &h__[krow + kwtop * 
		    h_dim1], ldh, &v[v_offset], ldv, &c_b1, &wv[wv_offset], 
		    ldwv, (ftnlen)1, (ftnlen)1);
#line 532 "zlaqr3.f"
	    zlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * 
		    h_dim1], ldh, (ftnlen)1);
#line 533 "zlaqr3.f"
/* L60: */
#line 533 "zlaqr3.f"
	}

/*        ==== Update horizontal slab in H ==== */

#line 537 "zlaqr3.f"
	if (*wantt) {
#line 538 "zlaqr3.f"
	    i__2 = *n;
#line 538 "zlaqr3.f"
	    i__1 = *nh;
#line 538 "zlaqr3.f"
	    for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; 
		    kcol += i__1) {
/* Computing MIN */
#line 539 "zlaqr3.f"
		i__3 = *nh, i__4 = *n - kcol + 1;
#line 539 "zlaqr3.f"
		kln = min(i__3,i__4);
#line 540 "zlaqr3.f"
		zgemm_("C", "N", &jw, &kln, &jw, &c_b2, &v[v_offset], ldv, &
			h__[kwtop + kcol * h_dim1], ldh, &c_b1, &t[t_offset], 
			ldt, (ftnlen)1, (ftnlen)1);
#line 542 "zlaqr3.f"
		zlacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol *
			 h_dim1], ldh, (ftnlen)1);
#line 544 "zlaqr3.f"
/* L70: */
#line 544 "zlaqr3.f"
	    }
#line 545 "zlaqr3.f"
	}

/*        ==== Update vertical slab in Z ==== */

#line 549 "zlaqr3.f"
	if (*wantz) {
#line 550 "zlaqr3.f"
	    i__1 = *ihiz;
#line 550 "zlaqr3.f"
	    i__2 = *nv;
#line 550 "zlaqr3.f"
	    for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow +=
		     i__2) {
/* Computing MIN */
#line 551 "zlaqr3.f"
		i__3 = *nv, i__4 = *ihiz - krow + 1;
#line 551 "zlaqr3.f"
		kln = min(i__3,i__4);
#line 552 "zlaqr3.f"
		zgemm_("N", "N", &kln, &jw, &jw, &c_b2, &z__[krow + kwtop * 
			z_dim1], ldz, &v[v_offset], ldv, &c_b1, &wv[wv_offset]
			, ldwv, (ftnlen)1, (ftnlen)1);
#line 554 "zlaqr3.f"
		zlacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + 
			kwtop * z_dim1], ldz, (ftnlen)1);
#line 556 "zlaqr3.f"
/* L80: */
#line 556 "zlaqr3.f"
	    }
#line 557 "zlaqr3.f"
	}
#line 558 "zlaqr3.f"
    }

/*     ==== Return the number of deflations ... ==== */

#line 562 "zlaqr3.f"
    *nd = jw - *ns;

/*     ==== ... and the number of shifts. (Subtracting */
/*     .    INFQR from the spike length takes care */
/*     .    of the case of a rare QR failure while */
/*     .    calculating eigenvalues of the deflation */
/*     .    window.)  ==== */

#line 570 "zlaqr3.f"
    *ns -= infqr;

/*      ==== Return optimal workspace. ==== */

#line 574 "zlaqr3.f"
    d__1 = (doublereal) lwkopt;
#line 574 "zlaqr3.f"
    z__1.r = d__1, z__1.i = 0.;
#line 574 "zlaqr3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

/*     ==== End of ZLAQR3 ==== */

#line 578 "zlaqr3.f"
    return 0;
} /* zlaqr3_ */

