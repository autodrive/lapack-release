#line 1 "claqr2.f"
/* claqr2.f -- translated by f2c (version 20100827).
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

#line 1 "claqr2.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;

/* > \brief \b CLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to detect and defl
ate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAQR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, */
/*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, */
/*                          NV, WV, LDWV, WORK, LWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, */
/*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW */
/*       LOGICAL            WANTT, WANTZ */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CLAQR2 is identical to CLAQR3 except that it avoids */
/* >    recursion by calling CLAHQR instead of CLAQR4. */
/* > */
/* >    Aggressive early deflation: */
/* > */
/* >    This subroutine accepts as input an upper Hessenberg matrix */
/* >    H and performs an unitary similarity transformation */
/* >    designed to detect and deflate fully converged eigenvalues from */
/* >    a trailing principal submatrix.  On output H has been over- */
/* >    written by a new Hessenberg matrix that is a perturbation of */
/* >    an unitary similarity transformation of H.  It is to be */
/* >    hoped that the final version of H has many zero subdiagonal */
/* >    entries. */
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
/* >          H is COMPLEX array, dimension (LDH,N) */
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
/* >          Z is COMPLEX array, dimension (LDZ,N) */
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
/* >          SH is COMPLEX array, dimension (KBOT) */
/* >          On output, approximate eigenvalues that may */
/* >          be used for shifts are stored in SH(KBOT-ND-NS+1) */
/* >          through SR(KBOT-ND).  Converged eigenvalues are */
/* >          stored in SH(KBOT-ND+1) through SH(KBOT). */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,NW) */
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
/* >          T is COMPLEX array, dimension (LDT,NW) */
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
/* >          WV is COMPLEX array, dimension (LDWV,NW) */
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
/* >          WORK is COMPLEX array, dimension (LWORK) */
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
/* >          If LWORK = -1, then a workspace query is assumed; CLAQR2 */
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

/* > \ingroup complexOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >       Karen Braman and Ralph Byers, Department of Mathematics, */
/* >       University of Kansas, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int claqr2_(logical *wantt, logical *wantz, integer *n, 
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
    static integer lwk1, lwk2;
    static doublecomplex beta;
    static integer kcol, info, ifst, ilst, ltop, krow;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), cgemm_(char *, char *, 
	    integer *, integer *, integer *, doublecomplex *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen, ftnlen), ccopy_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer infqr, kwtop;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *), cgehrd_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), clarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clahqr_(logical *, logical *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, integer *, integer *), 
	    clacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), claset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen);
    static doublereal safmin, safmax;
    extern /* Subroutine */ int ctrexc_(char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, integer *, integer *, 
	    integer *, ftnlen), cunmhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Estimate optimal workspace. ==== */

#line 323 "claqr2.f"
    /* Parameter adjustments */
#line 323 "claqr2.f"
    h_dim1 = *ldh;
#line 323 "claqr2.f"
    h_offset = 1 + h_dim1;
#line 323 "claqr2.f"
    h__ -= h_offset;
#line 323 "claqr2.f"
    z_dim1 = *ldz;
#line 323 "claqr2.f"
    z_offset = 1 + z_dim1;
#line 323 "claqr2.f"
    z__ -= z_offset;
#line 323 "claqr2.f"
    --sh;
#line 323 "claqr2.f"
    v_dim1 = *ldv;
#line 323 "claqr2.f"
    v_offset = 1 + v_dim1;
#line 323 "claqr2.f"
    v -= v_offset;
#line 323 "claqr2.f"
    t_dim1 = *ldt;
#line 323 "claqr2.f"
    t_offset = 1 + t_dim1;
#line 323 "claqr2.f"
    t -= t_offset;
#line 323 "claqr2.f"
    wv_dim1 = *ldwv;
#line 323 "claqr2.f"
    wv_offset = 1 + wv_dim1;
#line 323 "claqr2.f"
    wv -= wv_offset;
#line 323 "claqr2.f"
    --work;
#line 323 "claqr2.f"

#line 323 "claqr2.f"
    /* Function Body */
/* Computing MIN */
#line 323 "claqr2.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 323 "claqr2.f"
    jw = min(i__1,i__2);
#line 324 "claqr2.f"
    if (jw <= 2) {
#line 325 "claqr2.f"
	lwkopt = 1;
#line 326 "claqr2.f"
    } else {

/*        ==== Workspace query call to CGEHRD ==== */

#line 330 "claqr2.f"
	i__1 = jw - 1;
#line 330 "claqr2.f"
	cgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &
		c_n1, &info);
#line 331 "claqr2.f"
	lwk1 = (integer) work[1].r;

/*        ==== Workspace query call to CUNMHR ==== */

#line 335 "claqr2.f"
	i__1 = jw - 1;
#line 335 "claqr2.f"
	cunmhr_("R", "N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1],
		 &v[v_offset], ldv, &work[1], &c_n1, &info, (ftnlen)1, (
		ftnlen)1);
#line 337 "claqr2.f"
	lwk2 = (integer) work[1].r;

/*        ==== Optimal workspace ==== */

#line 341 "claqr2.f"
	lwkopt = jw + max(lwk1,lwk2);
#line 342 "claqr2.f"
    }

/*     ==== Quick return in case of workspace query. ==== */

#line 346 "claqr2.f"
    if (*lwork == -1) {
#line 347 "claqr2.f"
	d__1 = (doublereal) lwkopt;
#line 347 "claqr2.f"
	z__1.r = d__1, z__1.i = 0.;
#line 347 "claqr2.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 348 "claqr2.f"
	return 0;
#line 349 "claqr2.f"
    }

/*     ==== Nothing to do ... */
/*     ... for an empty active block ... ==== */
#line 353 "claqr2.f"
    *ns = 0;
#line 354 "claqr2.f"
    *nd = 0;
#line 355 "claqr2.f"
    work[1].r = 1., work[1].i = 0.;
#line 356 "claqr2.f"
    if (*ktop > *kbot) {
#line 356 "claqr2.f"
	return 0;
#line 356 "claqr2.f"
    }
/*     ... nor for an empty deflation window. ==== */
#line 359 "claqr2.f"
    if (*nw < 1) {
#line 359 "claqr2.f"
	return 0;
#line 359 "claqr2.f"
    }

/*     ==== Machine constants ==== */

#line 364 "claqr2.f"
    safmin = slamch_("SAFE MINIMUM", (ftnlen)12);
#line 365 "claqr2.f"
    safmax = 1. / safmin;
#line 366 "claqr2.f"
    slabad_(&safmin, &safmax);
#line 367 "claqr2.f"
    ulp = slamch_("PRECISION", (ftnlen)9);
#line 368 "claqr2.f"
    smlnum = safmin * ((doublereal) (*n) / ulp);

/*     ==== Setup deflation window ==== */

/* Computing MIN */
#line 372 "claqr2.f"
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
#line 372 "claqr2.f"
    jw = min(i__1,i__2);
#line 373 "claqr2.f"
    kwtop = *kbot - jw + 1;
#line 374 "claqr2.f"
    if (kwtop == *ktop) {
#line 375 "claqr2.f"
	s.r = 0., s.i = 0.;
#line 376 "claqr2.f"
    } else {
#line 377 "claqr2.f"
	i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 377 "claqr2.f"
	s.r = h__[i__1].r, s.i = h__[i__1].i;
#line 378 "claqr2.f"
    }

#line 380 "claqr2.f"
    if (*kbot == kwtop) {

/*        ==== 1-by-1 deflation window: not much to do ==== */

#line 384 "claqr2.f"
	i__1 = kwtop;
#line 384 "claqr2.f"
	i__2 = kwtop + kwtop * h_dim1;
#line 384 "claqr2.f"
	sh[i__1].r = h__[i__2].r, sh[i__1].i = h__[i__2].i;
#line 385 "claqr2.f"
	*ns = 1;
#line 386 "claqr2.f"
	*nd = 0;
/* Computing MAX */
#line 387 "claqr2.f"
	i__1 = kwtop + kwtop * h_dim1;
#line 387 "claqr2.f"
	d__5 = smlnum, d__6 = ulp * ((d__1 = h__[i__1].r, abs(d__1)) + (d__2 =
		 d_imag(&h__[kwtop + kwtop * h_dim1]), abs(d__2)));
#line 387 "claqr2.f"
	if ((d__3 = s.r, abs(d__3)) + (d__4 = d_imag(&s), abs(d__4)) <= max(
		d__5,d__6)) {
#line 389 "claqr2.f"
	    *ns = 0;
#line 390 "claqr2.f"
	    *nd = 1;
#line 391 "claqr2.f"
	    if (kwtop > *ktop) {
#line 391 "claqr2.f"
		i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 391 "claqr2.f"
		h__[i__1].r = 0., h__[i__1].i = 0.;
#line 391 "claqr2.f"
	    }
#line 393 "claqr2.f"
	}
#line 394 "claqr2.f"
	work[1].r = 1., work[1].i = 0.;
#line 395 "claqr2.f"
	return 0;
#line 396 "claqr2.f"
    }

/*     ==== Convert to spike-triangular form.  (In case of a */
/*     .    rare QR failure, this routine continues to do */
/*     .    aggressive early deflation using that part of */
/*     .    the deflation window that converged using INFQR */
/*     .    here and there to keep track.) ==== */

#line 404 "claqr2.f"
    clacpy_("U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], 
	    ldt, (ftnlen)1);
#line 405 "claqr2.f"
    i__1 = jw - 1;
#line 405 "claqr2.f"
    i__2 = *ldh + 1;
#line 405 "claqr2.f"
    i__3 = *ldt + 1;
#line 405 "claqr2.f"
    ccopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &
	    i__3);

#line 407 "claqr2.f"
    claset_("A", &jw, &jw, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 408 "claqr2.f"
    clahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sh[kwtop], 
	    &c__1, &jw, &v[v_offset], ldv, &infqr);

/*     ==== Deflation detection loop ==== */

#line 413 "claqr2.f"
    *ns = jw;
#line 414 "claqr2.f"
    ilst = infqr + 1;
#line 415 "claqr2.f"
    i__1 = jw;
#line 415 "claqr2.f"
    for (knt = infqr + 1; knt <= i__1; ++knt) {

/*        ==== Small spike tip deflation test ==== */

#line 419 "claqr2.f"
	i__2 = *ns + *ns * t_dim1;
#line 419 "claqr2.f"
	foo = (d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[*ns + *ns * 
		t_dim1]), abs(d__2));
#line 420 "claqr2.f"
	if (foo == 0.) {
#line 420 "claqr2.f"
	    foo = (d__1 = s.r, abs(d__1)) + (d__2 = d_imag(&s), abs(d__2));
#line 420 "claqr2.f"
	}
#line 422 "claqr2.f"
	i__2 = *ns * v_dim1 + 1;
/* Computing MAX */
#line 422 "claqr2.f"
	d__5 = smlnum, d__6 = ulp * foo;
#line 422 "claqr2.f"
	if (((d__1 = s.r, abs(d__1)) + (d__2 = d_imag(&s), abs(d__2))) * ((
		d__3 = v[i__2].r, abs(d__3)) + (d__4 = d_imag(&v[*ns * v_dim1 
		+ 1]), abs(d__4))) <= max(d__5,d__6)) {

/*           ==== One more converged eigenvalue ==== */

#line 427 "claqr2.f"
	    --(*ns);
#line 428 "claqr2.f"
	} else {

/*           ==== One undeflatable eigenvalue.  Move it up out of the */
/*           .    way.   (CTREXC can not fail in this case.) ==== */

#line 433 "claqr2.f"
	    ifst = *ns;
#line 434 "claqr2.f"
	    ctrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &
		    ilst, &info, (ftnlen)1);
#line 435 "claqr2.f"
	    ++ilst;
#line 436 "claqr2.f"
	}
#line 437 "claqr2.f"
/* L10: */
#line 437 "claqr2.f"
    }

/*        ==== Return to Hessenberg form ==== */

#line 441 "claqr2.f"
    if (*ns == 0) {
#line 441 "claqr2.f"
	s.r = 0., s.i = 0.;
#line 441 "claqr2.f"
    }

#line 444 "claqr2.f"
    if (*ns < jw) {

/*        ==== sorting the diagonal of T improves accuracy for */
/*        .    graded matrices.  ==== */

#line 449 "claqr2.f"
	i__1 = *ns;
#line 449 "claqr2.f"
	for (i__ = infqr + 1; i__ <= i__1; ++i__) {
#line 450 "claqr2.f"
	    ifst = i__;
#line 451 "claqr2.f"
	    i__2 = *ns;
#line 451 "claqr2.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 452 "claqr2.f"
		i__3 = j + j * t_dim1;
#line 452 "claqr2.f"
		i__4 = ifst + ifst * t_dim1;
#line 452 "claqr2.f"
		if ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[j + j * 
			t_dim1]), abs(d__2)) > (d__3 = t[i__4].r, abs(d__3)) 
			+ (d__4 = d_imag(&t[ifst + ifst * t_dim1]), abs(d__4))
			) {
#line 452 "claqr2.f"
		    ifst = j;
#line 452 "claqr2.f"
		}
#line 454 "claqr2.f"
/* L20: */
#line 454 "claqr2.f"
	    }
#line 455 "claqr2.f"
	    ilst = i__;
#line 456 "claqr2.f"
	    if (ifst != ilst) {
#line 456 "claqr2.f"
		ctrexc_("V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst,
			 &ilst, &info, (ftnlen)1);
#line 456 "claqr2.f"
	    }
#line 458 "claqr2.f"
/* L30: */
#line 458 "claqr2.f"
	}
#line 459 "claqr2.f"
    }

/*     ==== Restore shift/eigenvalue array from T ==== */

#line 463 "claqr2.f"
    i__1 = jw;
#line 463 "claqr2.f"
    for (i__ = infqr + 1; i__ <= i__1; ++i__) {
#line 464 "claqr2.f"
	i__2 = kwtop + i__ - 1;
#line 464 "claqr2.f"
	i__3 = i__ + i__ * t_dim1;
#line 464 "claqr2.f"
	sh[i__2].r = t[i__3].r, sh[i__2].i = t[i__3].i;
#line 465 "claqr2.f"
/* L40: */
#line 465 "claqr2.f"
    }


#line 468 "claqr2.f"
    if (*ns < jw || s.r == 0. && s.i == 0.) {
#line 469 "claqr2.f"
	if (*ns > 1 && (s.r != 0. || s.i != 0.)) {

/*           ==== Reflect spike back into lower triangle ==== */

#line 473 "claqr2.f"
	    ccopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
#line 474 "claqr2.f"
	    i__1 = *ns;
#line 474 "claqr2.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 475 "claqr2.f"
		i__2 = i__;
#line 475 "claqr2.f"
		d_cnjg(&z__1, &work[i__]);
#line 475 "claqr2.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 476 "claqr2.f"
/* L50: */
#line 476 "claqr2.f"
	    }
#line 477 "claqr2.f"
	    beta.r = work[1].r, beta.i = work[1].i;
#line 478 "claqr2.f"
	    clarfg_(ns, &beta, &work[2], &c__1, &tau);
#line 479 "claqr2.f"
	    work[1].r = 1., work[1].i = 0.;

#line 481 "claqr2.f"
	    i__1 = jw - 2;
#line 481 "claqr2.f"
	    i__2 = jw - 2;
#line 481 "claqr2.f"
	    claset_("L", &i__1, &i__2, &c_b1, &c_b1, &t[t_dim1 + 3], ldt, (
		    ftnlen)1);

#line 483 "claqr2.f"
	    d_cnjg(&z__1, &tau);
#line 483 "claqr2.f"
	    clarf_("L", ns, &jw, &work[1], &c__1, &z__1, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 485 "claqr2.f"
	    clarf_("R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &
		    work[jw + 1], (ftnlen)1);
#line 487 "claqr2.f"
	    clarf_("R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &
		    work[jw + 1], (ftnlen)1);

#line 490 "claqr2.f"
	    i__1 = *lwork - jw;
#line 490 "claqr2.f"
	    cgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1]
		    , &i__1, &info);
#line 492 "claqr2.f"
	}

/*        ==== Copy updated reduced window into place ==== */

#line 496 "claqr2.f"
	if (kwtop > 1) {
#line 496 "claqr2.f"
	    i__1 = kwtop + (kwtop - 1) * h_dim1;
#line 496 "claqr2.f"
	    d_cnjg(&z__2, &v[v_dim1 + 1]);
#line 496 "claqr2.f"
	    z__1.r = s.r * z__2.r - s.i * z__2.i, z__1.i = s.r * z__2.i + s.i 
		    * z__2.r;
#line 496 "claqr2.f"
	    h__[i__1].r = z__1.r, h__[i__1].i = z__1.i;
#line 496 "claqr2.f"
	}
#line 498 "claqr2.f"
	clacpy_("U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1]
		, ldh, (ftnlen)1);
#line 499 "claqr2.f"
	i__1 = jw - 1;
#line 499 "claqr2.f"
	i__2 = *ldt + 1;
#line 499 "claqr2.f"
	i__3 = *ldh + 1;
#line 499 "claqr2.f"
	ccopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1],
		 &i__3);

/*        ==== Accumulate orthogonal matrix in order update */
/*        .    H and Z, if requested.  ==== */

#line 505 "claqr2.f"
	if (*ns > 1 && (s.r != 0. || s.i != 0.)) {
#line 505 "claqr2.f"
	    i__1 = *lwork - jw;
#line 505 "claqr2.f"
	    cunmhr_("R", "N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1],
		     &v[v_offset], ldv, &work[jw + 1], &i__1, &info, (ftnlen)
		    1, (ftnlen)1);
#line 505 "claqr2.f"
	}

/*        ==== Update vertical slab in H ==== */

#line 511 "claqr2.f"
	if (*wantt) {
#line 512 "claqr2.f"
	    ltop = 1;
#line 513 "claqr2.f"
	} else {
#line 514 "claqr2.f"
	    ltop = *ktop;
#line 515 "claqr2.f"
	}
#line 516 "claqr2.f"
	i__1 = kwtop - 1;
#line 516 "claqr2.f"
	i__2 = *nv;
#line 516 "claqr2.f"
	for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += 
		i__2) {
/* Computing MIN */
#line 517 "claqr2.f"
	    i__3 = *nv, i__4 = kwtop - krow;
#line 517 "claqr2.f"
	    kln = min(i__3,i__4);
#line 518 "claqr2.f"
	    cgemm_("N", "N", &kln, &jw, &jw, &c_b2, &h__[krow + kwtop * 
		    h_dim1], ldh, &v[v_offset], ldv, &c_b1, &wv[wv_offset], 
		    ldwv, (ftnlen)1, (ftnlen)1);
#line 520 "claqr2.f"
	    clacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * 
		    h_dim1], ldh, (ftnlen)1);
#line 521 "claqr2.f"
/* L60: */
#line 521 "claqr2.f"
	}

/*        ==== Update horizontal slab in H ==== */

#line 525 "claqr2.f"
	if (*wantt) {
#line 526 "claqr2.f"
	    i__2 = *n;
#line 526 "claqr2.f"
	    i__1 = *nh;
#line 526 "claqr2.f"
	    for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; 
		    kcol += i__1) {
/* Computing MIN */
#line 527 "claqr2.f"
		i__3 = *nh, i__4 = *n - kcol + 1;
#line 527 "claqr2.f"
		kln = min(i__3,i__4);
#line 528 "claqr2.f"
		cgemm_("C", "N", &jw, &kln, &jw, &c_b2, &v[v_offset], ldv, &
			h__[kwtop + kcol * h_dim1], ldh, &c_b1, &t[t_offset], 
			ldt, (ftnlen)1, (ftnlen)1);
#line 530 "claqr2.f"
		clacpy_("A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol *
			 h_dim1], ldh, (ftnlen)1);
#line 532 "claqr2.f"
/* L70: */
#line 532 "claqr2.f"
	    }
#line 533 "claqr2.f"
	}

/*        ==== Update vertical slab in Z ==== */

#line 537 "claqr2.f"
	if (*wantz) {
#line 538 "claqr2.f"
	    i__1 = *ihiz;
#line 538 "claqr2.f"
	    i__2 = *nv;
#line 538 "claqr2.f"
	    for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow +=
		     i__2) {
/* Computing MIN */
#line 539 "claqr2.f"
		i__3 = *nv, i__4 = *ihiz - krow + 1;
#line 539 "claqr2.f"
		kln = min(i__3,i__4);
#line 540 "claqr2.f"
		cgemm_("N", "N", &kln, &jw, &jw, &c_b2, &z__[krow + kwtop * 
			z_dim1], ldz, &v[v_offset], ldv, &c_b1, &wv[wv_offset]
			, ldwv, (ftnlen)1, (ftnlen)1);
#line 542 "claqr2.f"
		clacpy_("A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + 
			kwtop * z_dim1], ldz, (ftnlen)1);
#line 544 "claqr2.f"
/* L80: */
#line 544 "claqr2.f"
	    }
#line 545 "claqr2.f"
	}
#line 546 "claqr2.f"
    }

/*     ==== Return the number of deflations ... ==== */

#line 550 "claqr2.f"
    *nd = jw - *ns;

/*     ==== ... and the number of shifts. (Subtracting */
/*     .    INFQR from the spike length takes care */
/*     .    of the case of a rare QR failure while */
/*     .    calculating eigenvalues of the deflation */
/*     .    window.)  ==== */

#line 558 "claqr2.f"
    *ns -= infqr;

/*      ==== Return optimal workspace. ==== */

#line 562 "claqr2.f"
    d__1 = (doublereal) lwkopt;
#line 562 "claqr2.f"
    z__1.r = d__1, z__1.i = 0.;
#line 562 "claqr2.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

/*     ==== End of CLAQR2 ==== */

#line 566 "claqr2.f"
    return 0;
} /* claqr2_ */

