#line 1 "zhsein.f"
/* zhsein.f -- translated by f2c (version 20100827).
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

#line 1 "zhsein.f"
/* Table of constant values */

static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* > \brief \b ZHSEIN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHSEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhsein.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhsein.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhsein.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, */
/*                          IFAILR, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EIGSRC, INITV, SIDE */
/*       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IFAILL( * ), IFAILR( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   W( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHSEIN uses inverse iteration to find specified right and/or left */
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
/* >          = 'Q': the eigenvalues were found using ZHSEQR; thus, if */
/* >                 H has zero subdiagonal elements, and so is */
/* >                 block-triangular, then the j-th eigenvalue can be */
/* >                 assumed to be an eigenvalue of the block containing */
/* >                 the j-th row/column.  This property allows ZHSEIN to */
/* >                 perform inverse iteration on just one diagonal block. */
/* >          = 'N': no assumptions are made on the correspondence */
/* >                 between eigenvalues and diagonal blocks.  In this */
/* >                 case, ZHSEIN must always perform inverse iteration */
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
/* >          H is COMPLEX*16 array, dimension (LDH,N) */
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
/* >          W is COMPLEX*16 array, dimension (N) */
/* >          On entry, the eigenvalues of H. */
/* >          On exit, the real parts of W may have been altered since */
/* >          close eigenvalues are perturbed slightly in searching for */
/* >          independent eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,MM) */
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
/* >          VR is COMPLEX*16 array, dimension (LDVR,MM) */
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
/* >          WORK is COMPLEX*16 array, dimension (N*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup complex16OTHERcomputational */

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
/* Subroutine */ int zhsein_(char *side, char *eigsrc, char *initv, logical *
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
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlaein_(
	    logical *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal zlanhs_(char *, integer *, doublecomplex *, integer *, 
	    doublereal *, ftnlen);
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

#line 300 "zhsein.f"
    /* Parameter adjustments */
#line 300 "zhsein.f"
    --select;
#line 300 "zhsein.f"
    h_dim1 = *ldh;
#line 300 "zhsein.f"
    h_offset = 1 + h_dim1;
#line 300 "zhsein.f"
    h__ -= h_offset;
#line 300 "zhsein.f"
    --w;
#line 300 "zhsein.f"
    vl_dim1 = *ldvl;
#line 300 "zhsein.f"
    vl_offset = 1 + vl_dim1;
#line 300 "zhsein.f"
    vl -= vl_offset;
#line 300 "zhsein.f"
    vr_dim1 = *ldvr;
#line 300 "zhsein.f"
    vr_offset = 1 + vr_dim1;
#line 300 "zhsein.f"
    vr -= vr_offset;
#line 300 "zhsein.f"
    --work;
#line 300 "zhsein.f"
    --rwork;
#line 300 "zhsein.f"
    --ifaill;
#line 300 "zhsein.f"
    --ifailr;
#line 300 "zhsein.f"

#line 300 "zhsein.f"
    /* Function Body */
#line 300 "zhsein.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 301 "zhsein.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 302 "zhsein.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 304 "zhsein.f"
    fromqr = lsame_(eigsrc, "Q", (ftnlen)1, (ftnlen)1);

#line 306 "zhsein.f"
    noinit = lsame_(initv, "N", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors. */

#line 311 "zhsein.f"
    *m = 0;
#line 312 "zhsein.f"
    i__1 = *n;
#line 312 "zhsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 313 "zhsein.f"
	if (select[k]) {
#line 313 "zhsein.f"
	    ++(*m);
#line 313 "zhsein.f"
	}
#line 315 "zhsein.f"
/* L10: */
#line 315 "zhsein.f"
    }

#line 317 "zhsein.f"
    *info = 0;
#line 318 "zhsein.f"
    if (! rightv && ! leftv) {
#line 319 "zhsein.f"
	*info = -1;
#line 320 "zhsein.f"
    } else if (! fromqr && ! lsame_(eigsrc, "N", (ftnlen)1, (ftnlen)1)) {
#line 321 "zhsein.f"
	*info = -2;
#line 322 "zhsein.f"
    } else if (! noinit && ! lsame_(initv, "U", (ftnlen)1, (ftnlen)1)) {
#line 323 "zhsein.f"
	*info = -3;
#line 324 "zhsein.f"
    } else if (*n < 0) {
#line 325 "zhsein.f"
	*info = -5;
#line 326 "zhsein.f"
    } else if (*ldh < max(1,*n)) {
#line 327 "zhsein.f"
	*info = -7;
#line 328 "zhsein.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 329 "zhsein.f"
	*info = -10;
#line 330 "zhsein.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 331 "zhsein.f"
	*info = -12;
#line 332 "zhsein.f"
    } else if (*mm < *m) {
#line 333 "zhsein.f"
	*info = -13;
#line 334 "zhsein.f"
    }
#line 335 "zhsein.f"
    if (*info != 0) {
#line 336 "zhsein.f"
	i__1 = -(*info);
#line 336 "zhsein.f"
	xerbla_("ZHSEIN", &i__1, (ftnlen)6);
#line 337 "zhsein.f"
	return 0;
#line 338 "zhsein.f"
    }

/*     Quick return if possible. */

#line 342 "zhsein.f"
    if (*n == 0) {
#line 342 "zhsein.f"
	return 0;
#line 342 "zhsein.f"
    }

/*     Set machine-dependent constants. */

#line 347 "zhsein.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 348 "zhsein.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 349 "zhsein.f"
    smlnum = unfl * (*n / ulp);

#line 351 "zhsein.f"
    ldwork = *n;

#line 353 "zhsein.f"
    kl = 1;
#line 354 "zhsein.f"
    kln = 0;
#line 355 "zhsein.f"
    if (fromqr) {
#line 356 "zhsein.f"
	kr = 0;
#line 357 "zhsein.f"
    } else {
#line 358 "zhsein.f"
	kr = *n;
#line 359 "zhsein.f"
    }
#line 360 "zhsein.f"
    ks = 1;

#line 362 "zhsein.f"
    i__1 = *n;
#line 362 "zhsein.f"
    for (k = 1; k <= i__1; ++k) {
#line 363 "zhsein.f"
	if (select[k]) {

/*           Compute eigenvector(s) corresponding to W(K). */

#line 367 "zhsein.f"
	    if (fromqr) {

/*              If affiliation of eigenvalues is known, check whether */
/*              the matrix splits. */

/*              Determine KL and KR such that 1 <= KL <= K <= KR <= N */
/*              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or */
/*              KR = N). */

/*              Then inverse iteration can be performed with the */
/*              submatrix H(KL:N,KL:N) for a left eigenvector, and with */
/*              the submatrix H(1:KR,1:KR) for a right eigenvector. */

#line 380 "zhsein.f"
		i__2 = kl + 1;
#line 380 "zhsein.f"
		for (i__ = k; i__ >= i__2; --i__) {
#line 381 "zhsein.f"
		    i__3 = i__ + (i__ - 1) * h_dim1;
#line 381 "zhsein.f"
		    if (h__[i__3].r == 0. && h__[i__3].i == 0.) {
#line 381 "zhsein.f"
			goto L30;
#line 381 "zhsein.f"
		    }
#line 383 "zhsein.f"
/* L20: */
#line 383 "zhsein.f"
		}
#line 384 "zhsein.f"
L30:
#line 385 "zhsein.f"
		kl = i__;
#line 386 "zhsein.f"
		if (k > kr) {
#line 387 "zhsein.f"
		    i__2 = *n - 1;
#line 387 "zhsein.f"
		    for (i__ = k; i__ <= i__2; ++i__) {
#line 388 "zhsein.f"
			i__3 = i__ + 1 + i__ * h_dim1;
#line 388 "zhsein.f"
			if (h__[i__3].r == 0. && h__[i__3].i == 0.) {
#line 388 "zhsein.f"
			    goto L50;
#line 388 "zhsein.f"
			}
#line 390 "zhsein.f"
/* L40: */
#line 390 "zhsein.f"
		    }
#line 391 "zhsein.f"
L50:
#line 392 "zhsein.f"
		    kr = i__;
#line 393 "zhsein.f"
		}
#line 394 "zhsein.f"
	    }

#line 396 "zhsein.f"
	    if (kl != kln) {
#line 397 "zhsein.f"
		kln = kl;

/*              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it */
/*              has not ben computed before. */

#line 402 "zhsein.f"
		i__2 = kr - kl + 1;
#line 402 "zhsein.f"
		hnorm = zlanhs_("I", &i__2, &h__[kl + kl * h_dim1], ldh, &
			rwork[1], (ftnlen)1);
#line 403 "zhsein.f"
		if (disnan_(&hnorm)) {
#line 404 "zhsein.f"
		    *info = -6;
#line 405 "zhsein.f"
		    return 0;
#line 406 "zhsein.f"
		} else if (hnorm > 0.) {
#line 407 "zhsein.f"
		    eps3 = hnorm * ulp;
#line 408 "zhsein.f"
		} else {
#line 409 "zhsein.f"
		    eps3 = smlnum;
#line 410 "zhsein.f"
		}
#line 411 "zhsein.f"
	    }

/*           Perturb eigenvalue if it is close to any previous */
/*           selected eigenvalues affiliated to the submatrix */
/*           H(KL:KR,KL:KR). Close roots are modified by EPS3. */

#line 417 "zhsein.f"
	    i__2 = k;
#line 417 "zhsein.f"
	    wk.r = w[i__2].r, wk.i = w[i__2].i;
#line 418 "zhsein.f"
L60:
#line 419 "zhsein.f"
	    i__2 = kl;
#line 419 "zhsein.f"
	    for (i__ = k - 1; i__ >= i__2; --i__) {
#line 420 "zhsein.f"
		i__3 = i__;
#line 420 "zhsein.f"
		z__2.r = w[i__3].r - wk.r, z__2.i = w[i__3].i - wk.i;
#line 420 "zhsein.f"
		z__1.r = z__2.r, z__1.i = z__2.i;
#line 420 "zhsein.f"
		if (select[i__] && (d__1 = z__1.r, abs(d__1)) + (d__2 = 
			d_imag(&z__1), abs(d__2)) < eps3) {
#line 421 "zhsein.f"
		    z__1.r = wk.r + eps3, z__1.i = wk.i;
#line 421 "zhsein.f"
		    wk.r = z__1.r, wk.i = z__1.i;
#line 422 "zhsein.f"
		    goto L60;
#line 423 "zhsein.f"
		}
#line 424 "zhsein.f"
/* L70: */
#line 424 "zhsein.f"
	    }
#line 425 "zhsein.f"
	    i__2 = k;
#line 425 "zhsein.f"
	    w[i__2].r = wk.r, w[i__2].i = wk.i;

#line 427 "zhsein.f"
	    if (leftv) {

/*              Compute left eigenvector. */

#line 431 "zhsein.f"
		i__2 = *n - kl + 1;
#line 431 "zhsein.f"
		zlaein_(&c_false, &noinit, &i__2, &h__[kl + kl * h_dim1], ldh,
			 &wk, &vl[kl + ks * vl_dim1], &work[1], &ldwork, &
			rwork[1], &eps3, &smlnum, &iinfo);
#line 434 "zhsein.f"
		if (iinfo > 0) {
#line 435 "zhsein.f"
		    ++(*info);
#line 436 "zhsein.f"
		    ifaill[ks] = k;
#line 437 "zhsein.f"
		} else {
#line 438 "zhsein.f"
		    ifaill[ks] = 0;
#line 439 "zhsein.f"
		}
#line 440 "zhsein.f"
		i__2 = kl - 1;
#line 440 "zhsein.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 441 "zhsein.f"
		    i__3 = i__ + ks * vl_dim1;
#line 441 "zhsein.f"
		    vl[i__3].r = 0., vl[i__3].i = 0.;
#line 442 "zhsein.f"
/* L80: */
#line 442 "zhsein.f"
		}
#line 443 "zhsein.f"
	    }
#line 444 "zhsein.f"
	    if (rightv) {

/*              Compute right eigenvector. */

#line 448 "zhsein.f"
		zlaein_(&c_true, &noinit, &kr, &h__[h_offset], ldh, &wk, &vr[
			ks * vr_dim1 + 1], &work[1], &ldwork, &rwork[1], &
			eps3, &smlnum, &iinfo);
#line 450 "zhsein.f"
		if (iinfo > 0) {
#line 451 "zhsein.f"
		    ++(*info);
#line 452 "zhsein.f"
		    ifailr[ks] = k;
#line 453 "zhsein.f"
		} else {
#line 454 "zhsein.f"
		    ifailr[ks] = 0;
#line 455 "zhsein.f"
		}
#line 456 "zhsein.f"
		i__2 = *n;
#line 456 "zhsein.f"
		for (i__ = kr + 1; i__ <= i__2; ++i__) {
#line 457 "zhsein.f"
		    i__3 = i__ + ks * vr_dim1;
#line 457 "zhsein.f"
		    vr[i__3].r = 0., vr[i__3].i = 0.;
#line 458 "zhsein.f"
/* L90: */
#line 458 "zhsein.f"
		}
#line 459 "zhsein.f"
	    }
#line 460 "zhsein.f"
	    ++ks;
#line 461 "zhsein.f"
	}
#line 462 "zhsein.f"
/* L100: */
#line 462 "zhsein.f"
    }

#line 464 "zhsein.f"
    return 0;

/*     End of ZHSEIN */

} /* zhsein_ */

