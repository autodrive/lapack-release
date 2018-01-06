#line 1 "chsein.f"
/* chsein.f -- translated by f2c (version 20100827).
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

#line 1 "chsein.f"
/* Table of constant values */

static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* > \brief \b CHSEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHSEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chsein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chsein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chsein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, */
/*                          IFAILR, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EIGSRC, INITV, SIDE */
/*       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IFAILL( * ), IFAILR( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHSEIN uses inverse iteration to find specified right and/or left */
/* > eigenvectors of a complex upper Hessenberg matrix H. */
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
/* >          Specifies the source of eigenvalues supplied in W: */
/* >          = 'Q': the eigenvalues were found using CHSEQR; thus, if */
/* >                 H has zero subdiagonal elements, and so is */
/* >                 block-triangular, then the j-th eigenvalue can be */
/* >                 assumed to be an eigenvalue of the block containing */
/* >                 the j-th row/column.  This property allows CHSEIN to */
/* >                 perform inverse iteration on just one diagonal block. */
/* >          = 'N': no assumptions are made on the correspondence */
/* >                 between eigenvalues and diagonal blocks.  In this */
/* >                 case, CHSEIN must always perform inverse iteration */
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
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          Specifies the eigenvectors to be computed. To select the */
/* >          eigenvector corresponding to the eigenvalue W(j), */
/* >          SELECT(j) must be set to .TRUE.. */
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
/* >          H is COMPLEX array, dimension (LDH,N) */
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
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is COMPLEX array, dimension (N) */
/* >          On entry, the eigenvalues of H. */
/* >          On exit, the real parts of W may have been altered since */
/* >          close eigenvalues are perturbed slightly in searching for */
/* >          independent eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,MM) */
/* >          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must */
/* >          contain starting vectors for the inverse iteration for the */
/* >          left eigenvectors; the starting vector for each eigenvector */
/* >          must be in the same column in which the eigenvector will be */
/* >          stored. */
/* >          On exit, if SIDE = 'L' or 'B', the left eigenvectors */
/* >          specified by SELECT will be stored consecutively in the */
/* >          columns of VL, in the same order as their eigenvalues. */
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
/* >          VR is COMPLEX array, dimension (LDVR,MM) */
/* >          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must */
/* >          contain starting vectors for the inverse iteration for the */
/* >          right eigenvectors; the starting vector for each eigenvector */
/* >          must be in the same column in which the eigenvector will be */
/* >          stored. */
/* >          On exit, if SIDE = 'R' or 'B', the right eigenvectors */
/* >          specified by SELECT will be stored consecutively in the */
/* >          columns of VR, in the same order as their eigenvalues. */
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
/* >          store the eigenvectors (= the number of .TRUE. elements in */
/* >          SELECT). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILL */
/* > \verbatim */
/* >          IFAILL is INTEGER array, dimension (MM) */
/* >          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left */
/* >          eigenvector in the i-th column of VL (corresponding to the */
/* >          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the */
/* >          eigenvector converged satisfactorily. */
/* >          If SIDE = 'R', IFAILL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILR */
/* > \verbatim */
/* >          IFAILR is INTEGER array, dimension (MM) */
/* >          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right */
/* >          eigenvector in the i-th column of VR (corresponding to the */
/* >          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the */
/* >          eigenvector converged satisfactorily. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int chsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *
	w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
	 integer *mm, integer *m, doublecomplex *work, doublereal *rwork, 
	integer *ifaill, integer *ifailr, integer *info, ftnlen side_len, 
	ftnlen eigsrc_len, ftnlen initv_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, k, kl, kr, ks;
    static doublecomplex wk;
    static integer kln;
    static doublereal ulp, eps3, unfl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical leftv, bothv;
    static doublereal hnorm;
    extern /* Subroutine */ int claein_(logical *, logical *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen), clanhs_(char *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern logical sisnan_(doublereal *);
    static logical noinit;
    static integer ldwork;
    static logical rightv, fromqr;
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters. */

#line 300 "chsein.f"
    /* Parameter adjustments */
#line 300 "chsein.f"
    --select;
#line 300 "chsein.f"
    h_dim1 = *ldh;
#line 300 "chsein.f"
    h_offset = 1 + h_dim1;
#line 300 "chsein.f"
    h__ -= h_offset;
#line 300 "chsein.f"
    --w;
#line 300 "chsein.f"
    vl_dim1 = *ldvl;
#line 300 "chsein.f"
    vl_offset = 1 + vl_dim1;
#line 300 "chsein.f"
    vl -= vl_offset;
#line 300 "chsein.f"
    vr_dim1 = *ldvr;
#line 300 "chsein.f"
    vr_offset = 1 + vr_dim1;
#line 300 "chsein.f"
    vr -= vr_offset;
#line 300 "chsein.f"
    --work;
#line 300 "chsein.f"
    --rwork;
#line 300 "chsein.f"
    --ifaill;
#line 300 "chsein.f"
    --ifailr;
#line 300 "chsein.f"

#line 300 "chsein.f"
    /* Function Body */
#line 300 "chsein.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 301 "chsein.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 302 "chsein.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 304 "chsein.f"
    fromqr = lsame_(eigsrc, "Q", (ftnlen)1, (ftnlen)1);

#line 306 "chsein.f"
    noinit = lsame_(initv, "N", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors. */

#line 311 "chsein.f"
    *m = 0;
#line 312 "chsein.f"
    i__1 = *n;
#line 312 "chsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 313 "chsein.f"
	if (select[k]) {
#line 313 "chsein.f"
	    ++(*m);
#line 313 "chsein.f"
	}
#line 315 "chsein.f"
/* L10: */
#line 315 "chsein.f"
    }

#line 317 "chsein.f"
    *info = 0;
#line 318 "chsein.f"
    if (! rightv && ! leftv) {
#line 319 "chsein.f"
	*info = -1;
#line 320 "chsein.f"
    } else if (! fromqr && ! lsame_(eigsrc, "N", (ftnlen)1, (ftnlen)1)) {
#line 321 "chsein.f"
	*info = -2;
#line 322 "chsein.f"
    } else if (! noinit && ! lsame_(initv, "U", (ftnlen)1, (ftnlen)1)) {
#line 323 "chsein.f"
	*info = -3;
#line 324 "chsein.f"
    } else if (*n < 0) {
#line 325 "chsein.f"
	*info = -5;
#line 326 "chsein.f"
    } else if (*ldh < max(1,*n)) {
#line 327 "chsein.f"
	*info = -7;
#line 328 "chsein.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 329 "chsein.f"
	*info = -10;
#line 330 "chsein.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 331 "chsein.f"
	*info = -12;
#line 332 "chsein.f"
    } else if (*mm < *m) {
#line 333 "chsein.f"
	*info = -13;
#line 334 "chsein.f"
    }
#line 335 "chsein.f"
    if (*info != 0) {
#line 336 "chsein.f"
	i__1 = -(*info);
#line 336 "chsein.f"
	xerbla_("CHSEIN", &i__1, (ftnlen)6);
#line 337 "chsein.f"
	return 0;
#line 338 "chsein.f"
    }

/*     Quick return if possible. */

#line 342 "chsein.f"
    if (*n == 0) {
#line 342 "chsein.f"
	return 0;
#line 342 "chsein.f"
    }

/*     Set machine-dependent constants. */

#line 347 "chsein.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 348 "chsein.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 349 "chsein.f"
    smlnum = unfl * (*n / ulp);

#line 351 "chsein.f"
    ldwork = *n;

#line 353 "chsein.f"
    kl = 1;
#line 354 "chsein.f"
    kln = 0;
#line 355 "chsein.f"
    if (fromqr) {
#line 356 "chsein.f"
	kr = 0;
#line 357 "chsein.f"
    } else {
#line 358 "chsein.f"
	kr = *n;
#line 359 "chsein.f"
    }
#line 360 "chsein.f"
    ks = 1;

#line 362 "chsein.f"
    i__1 = *n;
#line 362 "chsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 363 "chsein.f"
	if (select[k]) {

/*           Compute eigenvector(s) corresponding to W(K). */

#line 367 "chsein.f"
	    if (fromqr) {

/*              If affiliation of eigenvalues is known, check whether */
/*              the matrix splits. */

/*              Determine KL and KR such that 1 <= KL <= K <= KR <= N */
/*              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or */
/*              KR = N). */

/*              Then inverse iteration can be performed with the */
/*              submatrix H(KL:N,KL:N) for a left eigenvector, and with */
/*              the submatrix H(1:KR,1:KR) for a right eigenvector. */

#line 380 "chsein.f"
		i__2 = kl + 1;
#line 380 "chsein.f"
		for (i__ = k; i__ >= i__2; --i__) {
#line 381 "chsein.f"
		    i__3 = i__ + (i__ - 1) * h_dim1;
#line 381 "chsein.f"
		    if (h__[i__3].r == 0. && h__[i__3].i == 0.) {
#line 381 "chsein.f"
			goto L30;
#line 381 "chsein.f"
		    }
#line 383 "chsein.f"
/* L20: */
#line 383 "chsein.f"
		}
#line 384 "chsein.f"
L30:
#line 385 "chsein.f"
		kl = i__;
#line 386 "chsein.f"
		if (k > kr) {
#line 387 "chsein.f"
		    i__2 = *n - 1;
#line 387 "chsein.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 388 "chsein.f"
			i__3 = i__ + 1 + i__ * h_dim1;
#line 388 "chsein.f"
			if (h__[i__3].r == 0. && h__[i__3].i == 0.) {
#line 388 "chsein.f"
			    goto L50;
#line 388 "chsein.f"
			}
#line 390 "chsein.f"
/* L40: */
#line 390 "chsein.f"
		    }
#line 391 "chsein.f"
L50:
#line 392 "chsein.f"
		    kr = i__;
#line 393 "chsein.f"
		}
#line 394 "chsein.f"
	    }

#line 396 "chsein.f"
	    if (kl != kln) {
#line 397 "chsein.f"
		kln = kl;

/*              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it */
/*              has not ben computed before. */

#line 402 "chsein.f"
		i__2 = kr - kl + 1;
#line 402 "chsein.f"
		hnorm = clanhs_("I", &i__2, &h__[kl + kl * h_dim1], ldh, &
			rwork[1], (ftnlen)1);
#line 403 "chsein.f"
		if (sisnan_(&hnorm)) {
#line 404 "chsein.f"
		    *info = -6;
#line 405 "chsein.f"
		    return 0;
#line 406 "chsein.f"
		} else if (hnorm > 0.) {
#line 407 "chsein.f"
		    eps3 = hnorm * ulp;
#line 408 "chsein.f"
		} else {
#line 409 "chsein.f"
		    eps3 = smlnum;
#line 410 "chsein.f"
		}
#line 411 "chsein.f"
	    }

/*           Perturb eigenvalue if it is close to any previous */
/*           selected eigenvalues affiliated to the submatrix */
/*           H(KL:KR,KL:KR). Close roots are modified by EPS3. */

#line 417 "chsein.f"
	    i__2 = k;
#line 417 "chsein.f"
	    wk.r = w[i__2].r, wk.i = w[i__2].i;
#line 418 "chsein.f"
L60:
#line 419 "chsein.f"
	    i__2 = kl;
#line 419 "chsein.f"
	    for (i__ = k - 1; i__ >= i__2; --i__) {
#line 420 "chsein.f"
		i__3 = i__;
#line 420 "chsein.f"
		z__2.r = w[i__3].r - wk.r, z__2.i = w[i__3].i - wk.i;
#line 420 "chsein.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 420 "chsein.f"
		if (select[i__] && (d__1 = z__1.r, abs(d__1)) + (d__2 = 
			d_imag(&z__1), abs(d__2)) < eps3) {
#line 421 "chsein.f"
		    z__1.r = wk.r + eps3, z__1.i = wk.i;
#line 421 "chsein.f"
		    wk.r = z__1.r, wk.i = z__1.i;
#line 422 "chsein.f"
		    goto L60;
#line 423 "chsein.f"
		}
#line 424 "chsein.f"
/* L70: */
#line 424 "chsein.f"
	    }
#line 425 "chsein.f"
	    i__2 = k;
#line 425 "chsein.f"
	    w[i__2].r = wk.r, w[i__2].i = wk.i;

#line 427 "chsein.f"
	    if (leftv) {

/*              Compute left eigenvector. */

#line 431 "chsein.f"
		i__2 = *n - kl + 1;
#line 431 "chsein.f"
		claein_(&c_false, &noinit, &i__2, &h__[kl + kl * h_dim1], ldh,
			 &wk, &vl[kl + ks * vl_dim1], &work[1], &ldwork, &
			rwork[1], &eps3, &smlnum, &iinfo);
#line 434 "chsein.f"
		if (iinfo > 0) {
#line 435 "chsein.f"
		    ++(*info);
#line 436 "chsein.f"
		    ifaill[ks] = k;
#line 437 "chsein.f"
		} else {
#line 438 "chsein.f"
		    ifaill[ks] = 0;
#line 439 "chsein.f"
		}
#line 440 "chsein.f"
		i__2 = kl - 1;
#line 440 "chsein.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 441 "chsein.f"
		    i__3 = i__ + ks * vl_dim1;
#line 441 "chsein.f"
		    vl[i__3].r = 0., vl[i__3].i = 0.;
#line 442 "chsein.f"
/* L80: */
#line 442 "chsein.f"
		}
#line 443 "chsein.f"
	    }
#line 444 "chsein.f"
	    if (rightv) {

/*              Compute right eigenvector. */

#line 448 "chsein.f"
		claein_(&c_true, &noinit, &kr, &h__[h_offset], ldh, &wk, &vr[
			ks * vr_dim1 + 1], &work[1], &ldwork, &rwork[1], &
			eps3, &smlnum, &iinfo);
#line 450 "chsein.f"
		if (iinfo > 0) {
#line 451 "chsein.f"
		    ++(*info);
#line 452 "chsein.f"
		    ifailr[ks] = k;
#line 453 "chsein.f"
		} else {
#line 454 "chsein.f"
		    ifailr[ks] = 0;
#line 455 "chsein.f"
		}
#line 456 "chsein.f"
		i__2 = *n;
#line 456 "chsein.f"
		for (i__ = kr + 1; i__ <= i__2; ++i__) {
#line 457 "chsein.f"
		    i__3 = i__ + ks * vr_dim1;
#line 457 "chsein.f"
		    vr[i__3].r = 0., vr[i__3].i = 0.;
#line 458 "chsein.f"
/* L90: */
#line 458 "chsein.f"
		}
#line 459 "chsein.f"
	    }
#line 460 "chsein.f"
	    ++ks;
#line 461 "chsein.f"
	}
#line 462 "chsein.f"
/* L100: */
#line 462 "chsein.f"
    }

#line 464 "chsein.f"
    return 0;

/*     End of CHSEIN */

} /* chsein_ */

