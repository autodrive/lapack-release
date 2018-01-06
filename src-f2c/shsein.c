#line 1 "shsein.f"
/* shsein.f -- translated by f2c (version 20100827).
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

#line 1 "shsein.f"
/* Table of constant values */

static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* > \brief \b SHSEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SHSEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/shsein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/shsein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/shsein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI, */
/*                          VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL, */
/*                          IFAILR, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EIGSRC, INITV, SIDE */
/*       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IFAILL( * ), IFAILR( * ) */
/*       REAL               H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WI( * ), WORK( * ), WR( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SHSEIN uses inverse iteration to find specified right and/or left */
/* > eigenvectors of a real upper Hessenberg matrix H. */
/* > */
/* > The right eigenvector x and the left eigenvector y of the matrix H */
/* > corresponding to an eigenvalue w are defined by: */
/* > */
/* >              H * x = w * x,     y**h * H = w * y**h */
/* > */
/* > where y**h denotes the conjugate transpose of the vector y. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R': compute right eigenvectors only; */
/* >          = 'L': compute left eigenvectors only; */
/* >          = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] EIGSRC */
/* > \verbatim */
/* >          EIGSRC is CHARACTER*1 */
/* >          Specifies the source of eigenvalues supplied in (WR,WI): */
/* >          = 'Q': the eigenvalues were found using SHSEQR; thus, if */
/* >                 H has zero subdiagonal elements, and so is */
/* >                 block-triangular, then the j-th eigenvalue can be */
/* >                 assumed to be an eigenvalue of the block containing */
/* >                 the j-th row/column.  This property allows SHSEIN to */
/* >                 perform inverse iteration on just one diagonal block. */
/* >          = 'N': no assumptions are made on the correspondence */
/* >                 between eigenvalues and diagonal blocks.  In this */
/* >                 case, SHSEIN must always perform inverse iteration */
/* >                 using the whole matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] INITV */
/* > \verbatim */
/* >          INITV is CHARACTER*1 */
/* >          = 'N': no initial vectors are supplied; */
/* >          = 'U': user-supplied initial vectors are stored in the arrays */
/* >                 VL and/or VR. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          Specifies the eigenvectors to be computed. To select the */
/* >          real eigenvector corresponding to a real eigenvalue WR(j), */
/* >          SELECT(j) must be set to .TRUE.. To select the complex */
/* >          eigenvector corresponding to a complex eigenvalue */
/* >          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)), */
/* >          either SELECT(j) or SELECT(j+1) or both must be set to */
/* >          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is */
/* >          .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix H.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is REAL array, dimension (LDH,N) */
/* >          The upper Hessenberg matrix H. */
/* >          If a NaN is detected in H, the routine will return with INFO=-6. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H.  LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] WR */
/* > \verbatim */
/* >          WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] WI */
/* > \verbatim */
/* >          WI is REAL array, dimension (N) */
/* > */
/* >          On entry, the real and imaginary parts of the eigenvalues of */
/* >          H; a complex conjugate pair of eigenvalues must be stored in */
/* >          consecutive elements of WR and WI. */
/* >          On exit, WR may have been altered since close eigenvalues */
/* >          are perturbed slightly in searching for independent */
/* >          eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is REAL array, dimension (LDVL,MM) */
/* >          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must */
/* >          contain starting vectors for the inverse iteration for the */
/* >          left eigenvectors; the starting vector for each eigenvector */
/* >          must be in the same column(s) in which the eigenvector will */
/* >          be stored. */
/* >          On exit, if SIDE = 'L' or 'B', the left eigenvectors */
/* >          specified by SELECT will be stored consecutively in the */
/* >          columns of VL, in the same order as their eigenvalues. A */
/* >          complex eigenvector corresponding to a complex eigenvalue is */
/* >          stored in two consecutive columns, the first holding the real */
/* >          part and the second the imaginary part. */
/* >          If SIDE = 'R', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL. */
/* >          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is REAL array, dimension (LDVR,MM) */
/* >          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must */
/* >          contain starting vectors for the inverse iteration for the */
/* >          right eigenvectors; the starting vector for each eigenvector */
/* >          must be in the same column(s) in which the eigenvector will */
/* >          be stored. */
/* >          On exit, if SIDE = 'R' or 'B', the right eigenvectors */
/* >          specified by SELECT will be stored consecutively in the */
/* >          columns of VR, in the same order as their eigenvalues. A */
/* >          complex eigenvector corresponding to a complex eigenvalue is */
/* >          stored in two consecutive columns, the first holding the real */
/* >          part and the second the imaginary part. */
/* >          If SIDE = 'L', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR. */
/* >          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* >          MM is INTEGER */
/* >          The number of columns in the arrays VL and/or VR. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of columns in the arrays VL and/or VR required to */
/* >          store the eigenvectors; each selected real eigenvector */
/* >          occupies one column and each selected complex eigenvector */
/* >          occupies two columns. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension ((N+2)*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILL */
/* > \verbatim */
/* >          IFAILL is INTEGER array, dimension (MM) */
/* >          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left */
/* >          eigenvector in the i-th column of VL (corresponding to the */
/* >          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the */
/* >          eigenvector converged satisfactorily. If the i-th and (i+1)th */
/* >          columns of VL hold a complex eigenvector, then IFAILL(i) and */
/* >          IFAILL(i+1) are set to the same value. */
/* >          If SIDE = 'R', IFAILL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILR */
/* > \verbatim */
/* >          IFAILR is INTEGER array, dimension (MM) */
/* >          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right */
/* >          eigenvector in the i-th column of VR (corresponding to the */
/* >          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the */
/* >          eigenvector converged satisfactorily. If the i-th and (i+1)th */
/* >          columns of VR hold a complex eigenvector, then IFAILR(i) and */
/* >          IFAILR(i+1) are set to the same value. */
/* >          If SIDE = 'L', IFAILR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, i is the number of eigenvectors which */
/* >                failed to converge; see IFAILL and IFAILR for further */
/* >                details. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Each eigenvector is normalized so that the element of largest */
/* >  magnitude has magnitude 1; here the magnitude of a complex number */
/* >  (x,y) is taken to be |x|+|y|. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int shsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
	ifaill, integer *ifailr, integer *info, ftnlen side_len, ftnlen 
	eigsrc_len, ftnlen initv_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, kl, kr, kln, ksi;
    static doublereal wki;
    static integer ksr;
    static doublereal ulp, wkr, eps3;
    static logical pair;
    static doublereal unfl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical leftv, bothv;
    static doublereal hnorm;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slaein_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *), xerbla_(char *, integer 
	    *, ftnlen);
    static doublereal bignum;
    extern doublereal slanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical sisnan_(doublereal *);
    static logical noinit;
    static integer ldwork;
    static logical rightv, fromqr;
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

/*     Decode and test the input parameters. */

#line 309 "shsein.f"
    /* Parameter adjustments */
#line 309 "shsein.f"
    --select;
#line 309 "shsein.f"
    h_dim1 = *ldh;
#line 309 "shsein.f"
    h_offset = 1 + h_dim1;
#line 309 "shsein.f"
    h__ -= h_offset;
#line 309 "shsein.f"
    --wr;
#line 309 "shsein.f"
    --wi;
#line 309 "shsein.f"
    vl_dim1 = *ldvl;
#line 309 "shsein.f"
    vl_offset = 1 + vl_dim1;
#line 309 "shsein.f"
    vl -= vl_offset;
#line 309 "shsein.f"
    vr_dim1 = *ldvr;
#line 309 "shsein.f"
    vr_offset = 1 + vr_dim1;
#line 309 "shsein.f"
    vr -= vr_offset;
#line 309 "shsein.f"
    --work;
#line 309 "shsein.f"
    --ifaill;
#line 309 "shsein.f"
    --ifailr;
#line 309 "shsein.f"

#line 309 "shsein.f"
    /* Function Body */
#line 309 "shsein.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 310 "shsein.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 311 "shsein.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 313 "shsein.f"
    fromqr = lsame_(eigsrc, "Q", (ftnlen)1, (ftnlen)1);

#line 315 "shsein.f"
    noinit = lsame_(initv, "N", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors, and standardize the array SELECT. */

#line 320 "shsein.f"
    *m = 0;
#line 321 "shsein.f"
    pair = FALSE_;
#line 322 "shsein.f"
    i__1 = *n;
#line 322 "shsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 323 "shsein.f"
	if (pair) {
#line 324 "shsein.f"
	    pair = FALSE_;
#line 325 "shsein.f"
	    select[k] = FALSE_;
#line 326 "shsein.f"
	} else {
#line 327 "shsein.f"
	    if (wi[k] == 0.) {
#line 328 "shsein.f"
		if (select[k]) {
#line 328 "shsein.f"
		    ++(*m);
#line 328 "shsein.f"
		}
#line 330 "shsein.f"
	    } else {
#line 331 "shsein.f"
		pair = TRUE_;
#line 332 "shsein.f"
		if (select[k] || select[k + 1]) {
#line 333 "shsein.f"
		    select[k] = TRUE_;
#line 334 "shsein.f"
		    *m += 2;
#line 335 "shsein.f"
		}
#line 336 "shsein.f"
	    }
#line 337 "shsein.f"
	}
#line 338 "shsein.f"
/* L10: */
#line 338 "shsein.f"
    }

#line 340 "shsein.f"
    *info = 0;
#line 341 "shsein.f"
    if (! rightv && ! leftv) {
#line 342 "shsein.f"
	*info = -1;
#line 343 "shsein.f"
    } else if (! fromqr && ! lsame_(eigsrc, "N", (ftnlen)1, (ftnlen)1)) {
#line 344 "shsein.f"
	*info = -2;
#line 345 "shsein.f"
    } else if (! noinit && ! lsame_(initv, "U", (ftnlen)1, (ftnlen)1)) {
#line 346 "shsein.f"
	*info = -3;
#line 347 "shsein.f"
    } else if (*n < 0) {
#line 348 "shsein.f"
	*info = -5;
#line 349 "shsein.f"
    } else if (*ldh < max(1,*n)) {
#line 350 "shsein.f"
	*info = -7;
#line 351 "shsein.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 352 "shsein.f"
	*info = -11;
#line 353 "shsein.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 354 "shsein.f"
	*info = -13;
#line 355 "shsein.f"
    } else if (*mm < *m) {
#line 356 "shsein.f"
	*info = -14;
#line 357 "shsein.f"
    }
#line 358 "shsein.f"
    if (*info != 0) {
#line 359 "shsein.f"
	i__1 = -(*info);
#line 359 "shsein.f"
	xerbla_("SHSEIN", &i__1, (ftnlen)6);
#line 360 "shsein.f"
	return 0;
#line 361 "shsein.f"
    }

/*     Quick return if possible. */

#line 365 "shsein.f"
    if (*n == 0) {
#line 365 "shsein.f"
	return 0;
#line 365 "shsein.f"
    }

/*     Set machine-dependent constants. */

#line 370 "shsein.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 371 "shsein.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 372 "shsein.f"
    smlnum = unfl * (*n / ulp);
#line 373 "shsein.f"
    bignum = (1. - ulp) / smlnum;

#line 375 "shsein.f"
    ldwork = *n + 1;

#line 377 "shsein.f"
    kl = 1;
#line 378 "shsein.f"
    kln = 0;
#line 379 "shsein.f"
    if (fromqr) {
#line 380 "shsein.f"
	kr = 0;
#line 381 "shsein.f"
    } else {
#line 382 "shsein.f"
	kr = *n;
#line 383 "shsein.f"
    }
#line 384 "shsein.f"
    ksr = 1;

#line 386 "shsein.f"
    i__1 = *n;
#line 386 "shsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 387 "shsein.f"
	if (select[k]) {

/*           Compute eigenvector(s) corresponding to W(K). */

#line 391 "shsein.f"
	    if (fromqr) {

/*              If affiliation of eigenvalues is known, check whether */
/*              the matrix splits. */

/*              Determine KL and KR such that 1 <= KL <= K <= KR <= N */
/*              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or */
/*              KR = N). */

/*              Then inverse iteration can be performed with the */
/*              submatrix H(KL:N,KL:N) for a left eigenvector, and with */
/*              the submatrix H(1:KR,1:KR) for a right eigenvector. */

#line 404 "shsein.f"
		i__2 = kl + 1;
#line 404 "shsein.f"
		for (i__ = k; i__ >= i__2; --i__) {
#line 405 "shsein.f"
		    if (h__[i__ + (i__ - 1) * h_dim1] == 0.) {
#line 405 "shsein.f"
			goto L30;
#line 405 "shsein.f"
		    }
#line 407 "shsein.f"
/* L20: */
#line 407 "shsein.f"
		}
#line 408 "shsein.f"
L30:
#line 409 "shsein.f"
		kl = i__;
#line 410 "shsein.f"
		if (k > kr) {
#line 411 "shsein.f"
		    i__2 = *n - 1;
#line 411 "shsein.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 412 "shsein.f"
			if (h__[i__ + 1 + i__ * h_dim1] == 0.) {
#line 412 "shsein.f"
			    goto L50;
#line 412 "shsein.f"
			}
#line 414 "shsein.f"
/* L40: */
#line 414 "shsein.f"
		    }
#line 415 "shsein.f"
L50:
#line 416 "shsein.f"
		    kr = i__;
#line 417 "shsein.f"
		}
#line 418 "shsein.f"
	    }

#line 420 "shsein.f"
	    if (kl != kln) {
#line 421 "shsein.f"
		kln = kl;

/*              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it */
/*              has not ben computed before. */

#line 426 "shsein.f"
		i__2 = kr - kl + 1;
#line 426 "shsein.f"
		hnorm = slanhs_("I", &i__2, &h__[kl + kl * h_dim1], ldh, &
			work[1], (ftnlen)1);
#line 427 "shsein.f"
		if (sisnan_(&hnorm)) {
#line 428 "shsein.f"
		    *info = -6;
#line 429 "shsein.f"
		    return 0;
#line 430 "shsein.f"
		} else if (hnorm > 0.) {
#line 431 "shsein.f"
		    eps3 = hnorm * ulp;
#line 432 "shsein.f"
		} else {
#line 433 "shsein.f"
		    eps3 = smlnum;
#line 434 "shsein.f"
		}
#line 435 "shsein.f"
	    }

/*           Perturb eigenvalue if it is close to any previous */
/*           selected eigenvalues affiliated to the submatrix */
/*           H(KL:KR,KL:KR). Close roots are modified by EPS3. */

#line 441 "shsein.f"
	    wkr = wr[k];
#line 442 "shsein.f"
	    wki = wi[k];
#line 443 "shsein.f"
L60:
#line 444 "shsein.f"
	    i__2 = kl;
#line 444 "shsein.f"
	    for (i__ = k - 1; i__ >= i__2; --i__) {
#line 445 "shsein.f"
		if (select[i__] && (d__1 = wr[i__] - wkr, abs(d__1)) + (d__2 =
			 wi[i__] - wki, abs(d__2)) < eps3) {
#line 447 "shsein.f"
		    wkr += eps3;
#line 448 "shsein.f"
		    goto L60;
#line 449 "shsein.f"
		}
#line 450 "shsein.f"
/* L70: */
#line 450 "shsein.f"
	    }
#line 451 "shsein.f"
	    wr[k] = wkr;

#line 453 "shsein.f"
	    pair = wki != 0.;
#line 454 "shsein.f"
	    if (pair) {
#line 455 "shsein.f"
		ksi = ksr + 1;
#line 456 "shsein.f"
	    } else {
#line 457 "shsein.f"
		ksi = ksr;
#line 458 "shsein.f"
	    }
#line 459 "shsein.f"
	    if (leftv) {

/*              Compute left eigenvector. */

#line 463 "shsein.f"
		i__2 = *n - kl + 1;
#line 463 "shsein.f"
		slaein_(&c_false, &noinit, &i__2, &h__[kl + kl * h_dim1], ldh,
			 &wkr, &wki, &vl[kl + ksr * vl_dim1], &vl[kl + ksi * 
			vl_dim1], &work[1], &ldwork, &work[*n * *n + *n + 1], 
			&eps3, &smlnum, &bignum, &iinfo);
#line 467 "shsein.f"
		if (iinfo > 0) {
#line 468 "shsein.f"
		    if (pair) {
#line 469 "shsein.f"
			*info += 2;
#line 470 "shsein.f"
		    } else {
#line 471 "shsein.f"
			++(*info);
#line 472 "shsein.f"
		    }
#line 473 "shsein.f"
		    ifaill[ksr] = k;
#line 474 "shsein.f"
		    ifaill[ksi] = k;
#line 475 "shsein.f"
		} else {
#line 476 "shsein.f"
		    ifaill[ksr] = 0;
#line 477 "shsein.f"
		    ifaill[ksi] = 0;
#line 478 "shsein.f"
		}
#line 479 "shsein.f"
		i__2 = kl - 1;
#line 479 "shsein.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 480 "shsein.f"
		    vl[i__ + ksr * vl_dim1] = 0.;
#line 481 "shsein.f"
/* L80: */
#line 481 "shsein.f"
		}
#line 482 "shsein.f"
		if (pair) {
#line 483 "shsein.f"
		    i__2 = kl - 1;
#line 483 "shsein.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 484 "shsein.f"
			vl[i__ + ksi * vl_dim1] = 0.;
#line 485 "shsein.f"
/* L90: */
#line 485 "shsein.f"
		    }
#line 486 "shsein.f"
		}
#line 487 "shsein.f"
	    }
#line 488 "shsein.f"
	    if (rightv) {

/*              Compute right eigenvector. */

#line 492 "shsein.f"
		slaein_(&c_true, &noinit, &kr, &h__[h_offset], ldh, &wkr, &
			wki, &vr[ksr * vr_dim1 + 1], &vr[ksi * vr_dim1 + 1], &
			work[1], &ldwork, &work[*n * *n + *n + 1], &eps3, &
			smlnum, &bignum, &iinfo);
#line 496 "shsein.f"
		if (iinfo > 0) {
#line 497 "shsein.f"
		    if (pair) {
#line 498 "shsein.f"
			*info += 2;
#line 499 "shsein.f"
		    } else {
#line 500 "shsein.f"
			++(*info);
#line 501 "shsein.f"
		    }
#line 502 "shsein.f"
		    ifailr[ksr] = k;
#line 503 "shsein.f"
		    ifailr[ksi] = k;
#line 504 "shsein.f"
		} else {
#line 505 "shsein.f"
		    ifailr[ksr] = 0;
#line 506 "shsein.f"
		    ifailr[ksi] = 0;
#line 507 "shsein.f"
		}
#line 508 "shsein.f"
		i__2 = *n;
#line 508 "shsein.f"
		for (i__ = kr + 1; i__ <= i__2; ++i__) {
#line 509 "shsein.f"
		    vr[i__ + ksr * vr_dim1] = 0.;
#line 510 "shsein.f"
/* L100: */
#line 510 "shsein.f"
		}
#line 511 "shsein.f"
		if (pair) {
#line 512 "shsein.f"
		    i__2 = *n;
#line 512 "shsein.f"
		    for (i__ = kr + 1; i__ <= i__2; ++i__) {
#line 513 "shsein.f"
			vr[i__ + ksi * vr_dim1] = 0.;
#line 514 "shsein.f"
/* L110: */
#line 514 "shsein.f"
		    }
#line 515 "shsein.f"
		}
#line 516 "shsein.f"
	    }

#line 518 "shsein.f"
	    if (pair) {
#line 519 "shsein.f"
		ksr += 2;
#line 520 "shsein.f"
	    } else {
#line 521 "shsein.f"
		++ksr;
#line 522 "shsein.f"
	    }
#line 523 "shsein.f"
	}
#line 524 "shsein.f"
/* L120: */
#line 524 "shsein.f"
    }

#line 526 "shsein.f"
    return 0;

/*     End of SHSEIN */

} /* shsein_ */

