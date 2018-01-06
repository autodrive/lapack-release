#line 1 "ctrevc.f"
/* ctrevc.f -- translated by f2c (version 20100827).
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

#line 1 "ctrevc.f"
/* Table of constant values */

static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CTREVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTREVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, MM, M, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTREVC computes some or all of the right and/or left eigenvectors of */
/* > a complex upper triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* >              T*x = w*x,     (y**H)*T = w*(y**H) */
/* > */
/* > where y**H denotes the conjugate transpose of the vector y. */
/* > The eigenvalues are not input to this routine, but are read directly */
/* > from the diagonal of T. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/* > input matrix.  If Q is the unitary factor that reduces a matrix A to */
/* > Schur form T, then Q*X and Q*Y are the matrices of right and left */
/* > eigenvectors of A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R':  compute right eigenvectors only; */
/* >          = 'L':  compute left eigenvectors only; */
/* >          = 'B':  compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A':  compute all right and/or left eigenvectors; */
/* >          = 'B':  compute all right and/or left eigenvectors, */
/* >                  backtransformed using the matrices supplied in */
/* >                  VR and/or VL; */
/* >          = 'S':  compute selected right and/or left eigenvectors, */
/* >                  as indicated by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/* >          computed. */
/* >          The eigenvector corresponding to the j-th eigenvalue is */
/* >          computed if SELECT(j) = .TRUE.. */
/* >          Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
/* >          The upper triangular matrix T.  T is modified, but restored */
/* >          on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is COMPLEX array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* >          Schur vectors returned by CHSEQR). */
/* >          On exit, if SIDE = 'L' or 'B', VL contains: */
/* >          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*Y; */
/* >          if HOWMNY = 'S', the left eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VL, in the same order as their */
/* >                           eigenvalues. */
/* >          Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL.  LDVL >= 1, and if */
/* >          SIDE = 'L' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is COMPLEX array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* >          Schur vectors returned by CHSEQR). */
/* >          On exit, if SIDE = 'R' or 'B', VR contains: */
/* >          if HOWMNY = 'A', the matrix X of right eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*X; */
/* >          if HOWMNY = 'S', the right eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VR, in the same order as their */
/* >                           eigenvalues. */
/* >          Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1, and if */
/* >          SIDE = 'R' or 'B'; LDVR >= N. */
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
/* >          The number of columns in the arrays VL and/or VR actually */
/* >          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/* >          is set to N.  Each selected eigenvector occupies one */
/* >          column. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The algorithm used in this program is basically backward (forward) */
/* >  substitution, with scaling to make the the code robust against */
/* >  possible overflow. */
/* > */
/* >  Each eigenvector is normalized so that the element of largest */
/* >  magnitude has magnitude 1; here the magnitude of a complex number */
/* >  (x,y) is taken to be |x| + |y|. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctrevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer 
	*m, doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	side_len, ftnlen howmny_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, ii, ki, is;
    static doublereal ulp;
    static logical allv;
    static doublereal unfl, ovfl, smin;
    static logical over;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal remax;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical leftv, bothv, somev;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern integer icamax_(integer *, doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen), 
	    clatrs_(char *, char *, char *, char *, integer *, doublecomplex *
	    , integer *, doublecomplex *, doublereal *, doublereal *, integer 
	    *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);
    static logical rightv;
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

/*     Decode and test the input parameters */

#line 274 "ctrevc.f"
    /* Parameter adjustments */
#line 274 "ctrevc.f"
    --select;
#line 274 "ctrevc.f"
    t_dim1 = *ldt;
#line 274 "ctrevc.f"
    t_offset = 1 + t_dim1;
#line 274 "ctrevc.f"
    t -= t_offset;
#line 274 "ctrevc.f"
    vl_dim1 = *ldvl;
#line 274 "ctrevc.f"
    vl_offset = 1 + vl_dim1;
#line 274 "ctrevc.f"
    vl -= vl_offset;
#line 274 "ctrevc.f"
    vr_dim1 = *ldvr;
#line 274 "ctrevc.f"
    vr_offset = 1 + vr_dim1;
#line 274 "ctrevc.f"
    vr -= vr_offset;
#line 274 "ctrevc.f"
    --work;
#line 274 "ctrevc.f"
    --rwork;
#line 274 "ctrevc.f"

#line 274 "ctrevc.f"
    /* Function Body */
#line 274 "ctrevc.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 275 "ctrevc.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 276 "ctrevc.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 278 "ctrevc.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 279 "ctrevc.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 280 "ctrevc.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

/*     Set M to the number of columns required to store the selected */
/*     eigenvectors. */

#line 285 "ctrevc.f"
    if (somev) {
#line 286 "ctrevc.f"
	*m = 0;
#line 287 "ctrevc.f"
	i__1 = *n;
#line 287 "ctrevc.f"
	for (j = 1; j <= i__1; ++j) {
#line 288 "ctrevc.f"
	    if (select[j]) {
#line 288 "ctrevc.f"
		++(*m);
#line 288 "ctrevc.f"
	    }
#line 290 "ctrevc.f"
/* L10: */
#line 290 "ctrevc.f"
	}
#line 291 "ctrevc.f"
    } else {
#line 292 "ctrevc.f"
	*m = *n;
#line 293 "ctrevc.f"
    }

#line 295 "ctrevc.f"
    *info = 0;
#line 296 "ctrevc.f"
    if (! rightv && ! leftv) {
#line 297 "ctrevc.f"
	*info = -1;
#line 298 "ctrevc.f"
    } else if (! allv && ! over && ! somev) {
#line 299 "ctrevc.f"
	*info = -2;
#line 300 "ctrevc.f"
    } else if (*n < 0) {
#line 301 "ctrevc.f"
	*info = -4;
#line 302 "ctrevc.f"
    } else if (*ldt < max(1,*n)) {
#line 303 "ctrevc.f"
	*info = -6;
#line 304 "ctrevc.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 305 "ctrevc.f"
	*info = -8;
#line 306 "ctrevc.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 307 "ctrevc.f"
	*info = -10;
#line 308 "ctrevc.f"
    } else if (*mm < *m) {
#line 309 "ctrevc.f"
	*info = -11;
#line 310 "ctrevc.f"
    }
#line 311 "ctrevc.f"
    if (*info != 0) {
#line 312 "ctrevc.f"
	i__1 = -(*info);
#line 312 "ctrevc.f"
	xerbla_("CTREVC", &i__1, (ftnlen)6);
#line 313 "ctrevc.f"
	return 0;
#line 314 "ctrevc.f"
    }

/*     Quick return if possible. */

#line 318 "ctrevc.f"
    if (*n == 0) {
#line 318 "ctrevc.f"
	return 0;
#line 318 "ctrevc.f"
    }

/*     Set the constants to control overflow. */

#line 323 "ctrevc.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 324 "ctrevc.f"
    ovfl = 1. / unfl;
#line 325 "ctrevc.f"
    slabad_(&unfl, &ovfl);
#line 326 "ctrevc.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 327 "ctrevc.f"
    smlnum = unfl * (*n / ulp);

/*     Store the diagonal elements of T in working array WORK. */

#line 331 "ctrevc.f"
    i__1 = *n;
#line 331 "ctrevc.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 332 "ctrevc.f"
	i__2 = i__ + *n;
#line 332 "ctrevc.f"
	i__3 = i__ + i__ * t_dim1;
#line 332 "ctrevc.f"
	work[i__2].r = t[i__3].r, work[i__2].i = t[i__3].i;
#line 333 "ctrevc.f"
/* L20: */
#line 333 "ctrevc.f"
    }

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 338 "ctrevc.f"
    rwork[1] = 0.;
#line 339 "ctrevc.f"
    i__1 = *n;
#line 339 "ctrevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 340 "ctrevc.f"
	i__2 = j - 1;
#line 340 "ctrevc.f"
	rwork[j] = scasum_(&i__2, &t[j * t_dim1 + 1], &c__1);
#line 341 "ctrevc.f"
/* L30: */
#line 341 "ctrevc.f"
    }

#line 343 "ctrevc.f"
    if (rightv) {

/*        Compute right eigenvectors. */

#line 347 "ctrevc.f"
	is = *m;
#line 348 "ctrevc.f"
	for (ki = *n; ki >= 1; --ki) {

#line 350 "ctrevc.f"
	    if (somev) {
#line 351 "ctrevc.f"
		if (! select[ki]) {
#line 351 "ctrevc.f"
		    goto L80;
#line 351 "ctrevc.f"
		}
#line 353 "ctrevc.f"
	    }
/* Computing MAX */
#line 354 "ctrevc.f"
	    i__1 = ki + ki * t_dim1;
#line 354 "ctrevc.f"
	    d__3 = ulp * ((d__1 = t[i__1].r, abs(d__1)) + (d__2 = d_imag(&t[
		    ki + ki * t_dim1]), abs(d__2)));
#line 354 "ctrevc.f"
	    smin = max(d__3,smlnum);

#line 356 "ctrevc.f"
	    work[1].r = 1., work[1].i = 0.;

/*           Form right-hand side. */

#line 360 "ctrevc.f"
	    i__1 = ki - 1;
#line 360 "ctrevc.f"
	    for (k = 1; k <= i__1; ++k) {
#line 361 "ctrevc.f"
		i__2 = k;
#line 361 "ctrevc.f"
		i__3 = k + ki * t_dim1;
#line 361 "ctrevc.f"
		z__1.r = -t[i__3].r, z__1.i = -t[i__3].i;
#line 361 "ctrevc.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 362 "ctrevc.f"
/* L40: */
#line 362 "ctrevc.f"
	    }

/*           Solve the triangular system: */
/*              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK. */

#line 367 "ctrevc.f"
	    i__1 = ki - 1;
#line 367 "ctrevc.f"
	    for (k = 1; k <= i__1; ++k) {
#line 368 "ctrevc.f"
		i__2 = k + k * t_dim1;
#line 368 "ctrevc.f"
		i__3 = k + k * t_dim1;
#line 368 "ctrevc.f"
		i__4 = ki + ki * t_dim1;
#line 368 "ctrevc.f"
		z__1.r = t[i__3].r - t[i__4].r, z__1.i = t[i__3].i - t[i__4]
			.i;
#line 368 "ctrevc.f"
		t[i__2].r = z__1.r, t[i__2].i = z__1.i;
#line 369 "ctrevc.f"
		i__2 = k + k * t_dim1;
#line 369 "ctrevc.f"
		if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * 
			t_dim1]), abs(d__2)) < smin) {
#line 369 "ctrevc.f"
		    i__3 = k + k * t_dim1;
#line 369 "ctrevc.f"
		    t[i__3].r = smin, t[i__3].i = 0.;
#line 369 "ctrevc.f"
		}
#line 371 "ctrevc.f"
/* L50: */
#line 371 "ctrevc.f"
	    }

#line 373 "ctrevc.f"
	    if (ki > 1) {
#line 374 "ctrevc.f"
		i__1 = ki - 1;
#line 374 "ctrevc.f"
		clatrs_("Upper", "No transpose", "Non-unit", "Y", &i__1, &t[
			t_offset], ldt, &work[1], &scale, &rwork[1], info, (
			ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 377 "ctrevc.f"
		i__1 = ki;
#line 377 "ctrevc.f"
		work[i__1].r = scale, work[i__1].i = 0.;
#line 378 "ctrevc.f"
	    }

/*           Copy the vector x or Q*x to VR and normalize. */

#line 382 "ctrevc.f"
	    if (! over) {
#line 383 "ctrevc.f"
		ccopy_(&ki, &work[1], &c__1, &vr[is * vr_dim1 + 1], &c__1);

#line 385 "ctrevc.f"
		ii = icamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 386 "ctrevc.f"
		i__1 = ii + is * vr_dim1;
#line 386 "ctrevc.f"
		remax = 1. / ((d__1 = vr[i__1].r, abs(d__1)) + (d__2 = d_imag(
			&vr[ii + is * vr_dim1]), abs(d__2)));
#line 387 "ctrevc.f"
		csscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 389 "ctrevc.f"
		i__1 = *n;
#line 389 "ctrevc.f"
		for (k = ki + 1; k <= i__1; ++k) {
#line 390 "ctrevc.f"
		    i__2 = k + is * vr_dim1;
#line 390 "ctrevc.f"
		    vr[i__2].r = 0., vr[i__2].i = 0.;
#line 391 "ctrevc.f"
/* L60: */
#line 391 "ctrevc.f"
		}
#line 392 "ctrevc.f"
	    } else {
#line 393 "ctrevc.f"
		if (ki > 1) {
#line 393 "ctrevc.f"
		    i__1 = ki - 1;
#line 393 "ctrevc.f"
		    z__1.r = scale, z__1.i = 0.;
#line 393 "ctrevc.f"
		    cgemv_("N", n, &i__1, &c_b2, &vr[vr_offset], ldvr, &work[
			    1], &c__1, &z__1, &vr[ki * vr_dim1 + 1], &c__1, (
			    ftnlen)1);
#line 393 "ctrevc.f"
		}

#line 397 "ctrevc.f"
		ii = icamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 398 "ctrevc.f"
		i__1 = ii + ki * vr_dim1;
#line 398 "ctrevc.f"
		remax = 1. / ((d__1 = vr[i__1].r, abs(d__1)) + (d__2 = d_imag(
			&vr[ii + ki * vr_dim1]), abs(d__2)));
#line 399 "ctrevc.f"
		csscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
#line 400 "ctrevc.f"
	    }

/*           Set back the original diagonal elements of T. */

#line 404 "ctrevc.f"
	    i__1 = ki - 1;
#line 404 "ctrevc.f"
	    for (k = 1; k <= i__1; ++k) {
#line 405 "ctrevc.f"
		i__2 = k + k * t_dim1;
#line 405 "ctrevc.f"
		i__3 = k + *n;
#line 405 "ctrevc.f"
		t[i__2].r = work[i__3].r, t[i__2].i = work[i__3].i;
#line 406 "ctrevc.f"
/* L70: */
#line 406 "ctrevc.f"
	    }

#line 408 "ctrevc.f"
	    --is;
#line 409 "ctrevc.f"
L80:
#line 409 "ctrevc.f"
	    ;
#line 409 "ctrevc.f"
	}
#line 410 "ctrevc.f"
    }

#line 412 "ctrevc.f"
    if (leftv) {

/*        Compute left eigenvectors. */

#line 416 "ctrevc.f"
	is = 1;
#line 417 "ctrevc.f"
	i__1 = *n;
#line 417 "ctrevc.f"
	for (ki = 1; ki <= i__1; ++ki) {

#line 419 "ctrevc.f"
	    if (somev) {
#line 420 "ctrevc.f"
		if (! select[ki]) {
#line 420 "ctrevc.f"
		    goto L130;
#line 420 "ctrevc.f"
		}
#line 422 "ctrevc.f"
	    }
/* Computing MAX */
#line 423 "ctrevc.f"
	    i__2 = ki + ki * t_dim1;
#line 423 "ctrevc.f"
	    d__3 = ulp * ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[
		    ki + ki * t_dim1]), abs(d__2)));
#line 423 "ctrevc.f"
	    smin = max(d__3,smlnum);

#line 425 "ctrevc.f"
	    i__2 = *n;
#line 425 "ctrevc.f"
	    work[i__2].r = 1., work[i__2].i = 0.;

/*           Form right-hand side. */

#line 429 "ctrevc.f"
	    i__2 = *n;
#line 429 "ctrevc.f"
	    for (k = ki + 1; k <= i__2; ++k) {
#line 430 "ctrevc.f"
		i__3 = k;
#line 430 "ctrevc.f"
		d_cnjg(&z__2, &t[ki + k * t_dim1]);
#line 430 "ctrevc.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 430 "ctrevc.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 431 "ctrevc.f"
/* L90: */
#line 431 "ctrevc.f"
	    }

/*           Solve the triangular system: */
/*              (T(KI+1:N,KI+1:N) - T(KI,KI))**H*X = SCALE*WORK. */

#line 436 "ctrevc.f"
	    i__2 = *n;
#line 436 "ctrevc.f"
	    for (k = ki + 1; k <= i__2; ++k) {
#line 437 "ctrevc.f"
		i__3 = k + k * t_dim1;
#line 437 "ctrevc.f"
		i__4 = k + k * t_dim1;
#line 437 "ctrevc.f"
		i__5 = ki + ki * t_dim1;
#line 437 "ctrevc.f"
		z__1.r = t[i__4].r - t[i__5].r, z__1.i = t[i__4].i - t[i__5]
			.i;
#line 437 "ctrevc.f"
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
#line 438 "ctrevc.f"
		i__3 = k + k * t_dim1;
#line 438 "ctrevc.f"
		if ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * 
			t_dim1]), abs(d__2)) < smin) {
#line 438 "ctrevc.f"
		    i__4 = k + k * t_dim1;
#line 438 "ctrevc.f"
		    t[i__4].r = smin, t[i__4].i = 0.;
#line 438 "ctrevc.f"
		}
#line 440 "ctrevc.f"
/* L100: */
#line 440 "ctrevc.f"
	    }

#line 442 "ctrevc.f"
	    if (ki < *n) {
#line 443 "ctrevc.f"
		i__2 = *n - ki;
#line 443 "ctrevc.f"
		clatrs_("Upper", "Conjugate transpose", "Non-unit", "Y", &
			i__2, &t[ki + 1 + (ki + 1) * t_dim1], ldt, &work[ki + 
			1], &scale, &rwork[1], info, (ftnlen)5, (ftnlen)19, (
			ftnlen)8, (ftnlen)1);
#line 446 "ctrevc.f"
		i__2 = ki;
#line 446 "ctrevc.f"
		work[i__2].r = scale, work[i__2].i = 0.;
#line 447 "ctrevc.f"
	    }

/*           Copy the vector x or Q*x to VL and normalize. */

#line 451 "ctrevc.f"
	    if (! over) {
#line 452 "ctrevc.f"
		i__2 = *n - ki + 1;
#line 452 "ctrevc.f"
		ccopy_(&i__2, &work[ki], &c__1, &vl[ki + is * vl_dim1], &c__1)
			;

#line 454 "ctrevc.f"
		i__2 = *n - ki + 1;
#line 454 "ctrevc.f"
		ii = icamax_(&i__2, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
#line 455 "ctrevc.f"
		i__2 = ii + is * vl_dim1;
#line 455 "ctrevc.f"
		remax = 1. / ((d__1 = vl[i__2].r, abs(d__1)) + (d__2 = d_imag(
			&vl[ii + is * vl_dim1]), abs(d__2)));
#line 456 "ctrevc.f"
		i__2 = *n - ki + 1;
#line 456 "ctrevc.f"
		csscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 458 "ctrevc.f"
		i__2 = ki - 1;
#line 458 "ctrevc.f"
		for (k = 1; k <= i__2; ++k) {
#line 459 "ctrevc.f"
		    i__3 = k + is * vl_dim1;
#line 459 "ctrevc.f"
		    vl[i__3].r = 0., vl[i__3].i = 0.;
#line 460 "ctrevc.f"
/* L110: */
#line 460 "ctrevc.f"
		}
#line 461 "ctrevc.f"
	    } else {
#line 462 "ctrevc.f"
		if (ki < *n) {
#line 462 "ctrevc.f"
		    i__2 = *n - ki;
#line 462 "ctrevc.f"
		    z__1.r = scale, z__1.i = 0.;
#line 462 "ctrevc.f"
		    cgemv_("N", n, &i__2, &c_b2, &vl[(ki + 1) * vl_dim1 + 1], 
			    ldvl, &work[ki + 1], &c__1, &z__1, &vl[ki * 
			    vl_dim1 + 1], &c__1, (ftnlen)1);
#line 462 "ctrevc.f"
		}

#line 467 "ctrevc.f"
		ii = icamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 468 "ctrevc.f"
		i__2 = ii + ki * vl_dim1;
#line 468 "ctrevc.f"
		remax = 1. / ((d__1 = vl[i__2].r, abs(d__1)) + (d__2 = d_imag(
			&vl[ii + ki * vl_dim1]), abs(d__2)));
#line 469 "ctrevc.f"
		csscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
#line 470 "ctrevc.f"
	    }

/*           Set back the original diagonal elements of T. */

#line 474 "ctrevc.f"
	    i__2 = *n;
#line 474 "ctrevc.f"
	    for (k = ki + 1; k <= i__2; ++k) {
#line 475 "ctrevc.f"
		i__3 = k + k * t_dim1;
#line 475 "ctrevc.f"
		i__4 = k + *n;
#line 475 "ctrevc.f"
		t[i__3].r = work[i__4].r, t[i__3].i = work[i__4].i;
#line 476 "ctrevc.f"
/* L120: */
#line 476 "ctrevc.f"
	    }

#line 478 "ctrevc.f"
	    ++is;
#line 479 "ctrevc.f"
L130:
#line 479 "ctrevc.f"
	    ;
#line 479 "ctrevc.f"
	}
#line 480 "ctrevc.f"
    }

#line 482 "ctrevc.f"
    return 0;

/*     End of CTREVC */

} /* ctrevc_ */

