#line 1 "dbdsqr.f"
/* dbdsqr.f -- translated by f2c (version 20100827).
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

#line 1 "dbdsqr.f"
/* Table of constant values */

static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;

/* > \brief \b DBDSQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DBDSQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbdsqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbdsqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/*                          LDU, C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DBDSQR computes the singular values and, optionally, the right and/or */
/* > left singular vectors from the singular value decomposition (SVD) of */
/* > a real N-by-N (upper or lower) bidiagonal matrix B using the implicit */
/* > zero-shift QR algorithm.  The SVD of B has the form */
/* > */
/* >    B = Q * S * P**T */
/* > */
/* > where S is the diagonal matrix of singular values, Q is an orthogonal */
/* > matrix of left singular vectors, and P is an orthogonal matrix of */
/* > right singular vectors.  If left singular vectors are requested, this */
/* > subroutine actually returns U*Q instead of Q, and, if right singular */
/* > vectors are requested, this subroutine returns P**T*VT instead of */
/* > P**T, for given real input matrices U and VT.  When U and VT are the */
/* > orthogonal matrices that reduce a general matrix A to bidiagonal */
/* > form:  A = U*B*VT, as computed by DGEBRD, then */
/* > */
/* >    A = (U*Q) * S * (P**T*VT) */
/* > */
/* > is the SVD of A.  Optionally, the subroutine may also compute Q**T*C */
/* > for a given real input matrix C. */
/* > */
/* > See "Computing  Small Singular Values of Bidiagonal Matrices With */
/* > Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan, */
/* > LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11, */
/* > no. 5, pp. 873-912, Sept 1990) and */
/* > "Accurate singular values and differential qd algorithms," by */
/* > B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics */
/* > Department, University of California at Berkeley, July 1992 */
/* > for a detailed description of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  B is upper bidiagonal; */
/* >          = 'L':  B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCVT */
/* > \verbatim */
/* >          NCVT is INTEGER */
/* >          The number of columns of the matrix VT. NCVT >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRU */
/* > \verbatim */
/* >          NRU is INTEGER */
/* >          The number of rows of the matrix U. NRU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* >          NCC is INTEGER */
/* >          The number of columns of the matrix C. NCC >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the bidiagonal matrix B. */
/* >          On exit, if INFO=0, the singular values of B in decreasing */
/* >          order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the N-1 offdiagonal elements of the bidiagonal */
/* >          matrix B. */
/* >          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E */
/* >          will contain the diagonal and superdiagonal elements of a */
/* >          bidiagonal matrix orthogonally equivalent to the one given */
/* >          as input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT, NCVT) */
/* >          On entry, an N-by-NCVT matrix VT. */
/* >          On exit, VT is overwritten by P**T * VT. */
/* >          Not referenced if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT. */
/* >          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU, N) */
/* >          On entry, an NRU-by-N matrix U. */
/* >          On exit, U is overwritten by U * Q. */
/* >          Not referenced if NRU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= max(1,NRU). */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC, NCC) */
/* >          On entry, an N-by-NCC matrix C. */
/* >          On exit, C is overwritten by Q**T * C. */
/* >          Not referenced if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. */
/* >          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  If INFO = -i, the i-th argument had an illegal value */
/* >          > 0: */
/* >             if NCVT = NRU = NCC = 0, */
/* >                = 1, a split was marked by a positive value in E */
/* >                = 2, current block of Z not diagonalized after 30*N */
/* >                     iterations (in inner while loop) */
/* >                = 3, termination criterion of outer while loop not met */
/* >                     (program created more than N unreduced blocks) */
/* >             else NCVT = NRU = NCC = 0, */
/* >                   the algorithm did not converge; D and E contain the */
/* >                   elements of a bidiagonal matrix which is orthogonally */
/* >                   similar to the input matrix B;  if INFO = i, i */
/* >                   elements of E have not converged to zero. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8))) */
/* >          TOLMUL controls the convergence criterion of the QR loop. */
/* >          If it is positive, TOLMUL*EPS is the desired relative */
/* >             precision in the computed singular values. */
/* >          If it is negative, abs(TOLMUL*EPS*sigma_max) is the */
/* >             desired absolute accuracy in the computed singular */
/* >             values (corresponds to relative accuracy */
/* >             abs(TOLMUL*EPS) in the largest singular value. */
/* >          abs(TOLMUL) should be between 1 and 1/EPS, and preferably */
/* >             between 10 (for fast convergence) and .1/EPS */
/* >             (for there to be some accuracy in the results). */
/* >          Default is to lose at either one eighth or 2 of the */
/* >             available decimal digits in each computed singular value */
/* >             (whichever is smaller). */
/* > */
/* >  MAXITR  INTEGER, default = 6 */
/* >          MAXITR controls the maximum number of passes of the */
/* >          algorithm through its inner loop. The algorithms stops */
/* >          (and so fails to converge) if the number of passes */
/* >          through the inner loop exceeds MAXITR*N**2. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	ldc, doublereal *work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), d_sign(
	    doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__, j, m;
    static doublereal r__, cs;
    static integer ll;
    static doublereal sn, mu;
    static integer nm1, nm12, nm13, lll;
    static doublereal eps, sll, tol, abse;
    static integer idir;
    static doublereal abss;
    static integer oldm;
    static doublereal cosl;
    static integer isub, iter;
    static doublereal unfl, sinl, cosr, smin, smax, sinr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlas2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal oldcs;
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer oldll;
    static doublereal shift, sigmn, oldsn;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxit;
    static doublereal sminl, sigmx;
    static logical lower;
    extern /* Subroutine */ int dlasq1_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *), dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal sminoa, thresh;
    static logical rotate;
    static doublereal tolmul;


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 292 "dbdsqr.f"
    /* Parameter adjustments */
#line 292 "dbdsqr.f"
    --d__;
#line 292 "dbdsqr.f"
    --e;
#line 292 "dbdsqr.f"
    vt_dim1 = *ldvt;
#line 292 "dbdsqr.f"
    vt_offset = 1 + vt_dim1;
#line 292 "dbdsqr.f"
    vt -= vt_offset;
#line 292 "dbdsqr.f"
    u_dim1 = *ldu;
#line 292 "dbdsqr.f"
    u_offset = 1 + u_dim1;
#line 292 "dbdsqr.f"
    u -= u_offset;
#line 292 "dbdsqr.f"
    c_dim1 = *ldc;
#line 292 "dbdsqr.f"
    c_offset = 1 + c_dim1;
#line 292 "dbdsqr.f"
    c__ -= c_offset;
#line 292 "dbdsqr.f"
    --work;
#line 292 "dbdsqr.f"

#line 292 "dbdsqr.f"
    /* Function Body */
#line 292 "dbdsqr.f"
    *info = 0;
#line 293 "dbdsqr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 294 "dbdsqr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
#line 295 "dbdsqr.f"
	*info = -1;
#line 296 "dbdsqr.f"
    } else if (*n < 0) {
#line 297 "dbdsqr.f"
	*info = -2;
#line 298 "dbdsqr.f"
    } else if (*ncvt < 0) {
#line 299 "dbdsqr.f"
	*info = -3;
#line 300 "dbdsqr.f"
    } else if (*nru < 0) {
#line 301 "dbdsqr.f"
	*info = -4;
#line 302 "dbdsqr.f"
    } else if (*ncc < 0) {
#line 303 "dbdsqr.f"
	*info = -5;
#line 304 "dbdsqr.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 306 "dbdsqr.f"
	*info = -9;
#line 307 "dbdsqr.f"
    } else if (*ldu < max(1,*nru)) {
#line 308 "dbdsqr.f"
	*info = -11;
#line 309 "dbdsqr.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 311 "dbdsqr.f"
	*info = -13;
#line 312 "dbdsqr.f"
    }
#line 313 "dbdsqr.f"
    if (*info != 0) {
#line 314 "dbdsqr.f"
	i__1 = -(*info);
#line 314 "dbdsqr.f"
	xerbla_("DBDSQR", &i__1, (ftnlen)6);
#line 315 "dbdsqr.f"
	return 0;
#line 316 "dbdsqr.f"
    }
#line 317 "dbdsqr.f"
    if (*n == 0) {
#line 317 "dbdsqr.f"
	return 0;
#line 317 "dbdsqr.f"
    }
#line 319 "dbdsqr.f"
    if (*n == 1) {
#line 319 "dbdsqr.f"
	goto L160;
#line 319 "dbdsqr.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 324 "dbdsqr.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

#line 328 "dbdsqr.f"
    if (! rotate) {
#line 329 "dbdsqr.f"
	dlasq1_(n, &d__[1], &e[1], &work[1], info);

/*     If INFO equals 2, dqds didn't finish, try to finish */

#line 333 "dbdsqr.f"
	if (*info != 2) {
#line 333 "dbdsqr.f"
	    return 0;
#line 333 "dbdsqr.f"
	}
#line 334 "dbdsqr.f"
	*info = 0;
#line 335 "dbdsqr.f"
    }

#line 337 "dbdsqr.f"
    nm1 = *n - 1;
#line 338 "dbdsqr.f"
    nm12 = nm1 + nm1;
#line 339 "dbdsqr.f"
    nm13 = nm12 + nm1;
#line 340 "dbdsqr.f"
    idir = 0;

/*     Get machine constants */

#line 344 "dbdsqr.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 345 "dbdsqr.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 350 "dbdsqr.f"
    if (lower) {
#line 351 "dbdsqr.f"
	i__1 = *n - 1;
#line 351 "dbdsqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 352 "dbdsqr.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 353 "dbdsqr.f"
	    d__[i__] = r__;
#line 354 "dbdsqr.f"
	    e[i__] = sn * d__[i__ + 1];
#line 355 "dbdsqr.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 356 "dbdsqr.f"
	    work[i__] = cs;
#line 357 "dbdsqr.f"
	    work[nm1 + i__] = sn;
#line 358 "dbdsqr.f"
/* L10: */
#line 358 "dbdsqr.f"
	}

/*        Update singular vectors if desired */

#line 362 "dbdsqr.f"
	if (*nru > 0) {
#line 362 "dbdsqr.f"
	    dlasr_("R", "V", "F", nru, n, &work[1], &work[*n], &u[u_offset], 
		    ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 362 "dbdsqr.f"
	}
#line 365 "dbdsqr.f"
	if (*ncc > 0) {
#line 365 "dbdsqr.f"
	    dlasr_("L", "V", "F", n, ncc, &work[1], &work[*n], &c__[c_offset],
		     ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 365 "dbdsqr.f"
	}
#line 368 "dbdsqr.f"
    }

/*     Compute singular values to relative accuracy TOL */
/*     (By setting TOL to be negative, algorithm will compute */
/*     singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */

/* Computing MAX */
/* Computing MIN */
#line 374 "dbdsqr.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b15);
#line 374 "dbdsqr.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 374 "dbdsqr.f"
    tolmul = max(d__1,d__2);
#line 375 "dbdsqr.f"
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

#line 379 "dbdsqr.f"
    smax = 0.;
#line 380 "dbdsqr.f"
    i__1 = *n;
#line 380 "dbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 381 "dbdsqr.f"
	d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
#line 381 "dbdsqr.f"
	smax = max(d__2,d__3);
#line 382 "dbdsqr.f"
/* L20: */
#line 382 "dbdsqr.f"
    }
#line 383 "dbdsqr.f"
    i__1 = *n - 1;
#line 383 "dbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 384 "dbdsqr.f"
	d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
#line 384 "dbdsqr.f"
	smax = max(d__2,d__3);
#line 385 "dbdsqr.f"
/* L30: */
#line 385 "dbdsqr.f"
    }
#line 386 "dbdsqr.f"
    sminl = 0.;
#line 387 "dbdsqr.f"
    if (tol >= 0.) {

/*        Relative accuracy desired */

#line 391 "dbdsqr.f"
	sminoa = abs(d__[1]);
#line 392 "dbdsqr.f"
	if (sminoa == 0.) {
#line 392 "dbdsqr.f"
	    goto L50;
#line 392 "dbdsqr.f"
	}
#line 394 "dbdsqr.f"
	mu = sminoa;
#line 395 "dbdsqr.f"
	i__1 = *n;
#line 395 "dbdsqr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 396 "dbdsqr.f"
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
#line 397 "dbdsqr.f"
	    sminoa = min(sminoa,mu);
#line 398 "dbdsqr.f"
	    if (sminoa == 0.) {
#line 398 "dbdsqr.f"
		goto L50;
#line 398 "dbdsqr.f"
	    }
#line 400 "dbdsqr.f"
/* L40: */
#line 400 "dbdsqr.f"
	}
#line 401 "dbdsqr.f"
L50:
#line 402 "dbdsqr.f"
	sminoa /= sqrt((doublereal) (*n));
/* Computing MAX */
#line 403 "dbdsqr.f"
	d__1 = tol * sminoa, d__2 = *n * 6 * *n * unfl;
#line 403 "dbdsqr.f"
	thresh = max(d__1,d__2);
#line 404 "dbdsqr.f"
    } else {

/*        Absolute accuracy desired */

/* Computing MAX */
#line 408 "dbdsqr.f"
	d__1 = abs(tol) * smax, d__2 = *n * 6 * *n * unfl;
#line 408 "dbdsqr.f"
	thresh = max(d__1,d__2);
#line 409 "dbdsqr.f"
    }

/*     Prepare for main iteration loop for the singular values */
/*     (MAXIT is the maximum number of passes through the inner */
/*     loop permitted before nonconvergence signalled.) */

#line 415 "dbdsqr.f"
    maxit = *n * 6 * *n;
#line 416 "dbdsqr.f"
    iter = 0;
#line 417 "dbdsqr.f"
    oldll = -1;
#line 418 "dbdsqr.f"
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

#line 422 "dbdsqr.f"
    m = *n;

/*     Begin main iteration loop */

#line 426 "dbdsqr.f"
L60:

/*     Check for convergence or exceeding iteration count */

#line 430 "dbdsqr.f"
    if (m <= 1) {
#line 430 "dbdsqr.f"
	goto L160;
#line 430 "dbdsqr.f"
    }
#line 432 "dbdsqr.f"
    if (iter > maxit) {
#line 432 "dbdsqr.f"
	goto L200;
#line 432 "dbdsqr.f"
    }

/*     Find diagonal block of matrix to work on */

#line 437 "dbdsqr.f"
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
#line 437 "dbdsqr.f"
	d__[m] = 0.;
#line 437 "dbdsqr.f"
    }
#line 439 "dbdsqr.f"
    smax = (d__1 = d__[m], abs(d__1));
#line 440 "dbdsqr.f"
    smin = smax;
#line 441 "dbdsqr.f"
    i__1 = m - 1;
#line 441 "dbdsqr.f"
    for (lll = 1; lll <= i__1; ++lll) {
#line 442 "dbdsqr.f"
	ll = m - lll;
#line 443 "dbdsqr.f"
	abss = (d__1 = d__[ll], abs(d__1));
#line 444 "dbdsqr.f"
	abse = (d__1 = e[ll], abs(d__1));
#line 445 "dbdsqr.f"
	if (tol < 0. && abss <= thresh) {
#line 445 "dbdsqr.f"
	    d__[ll] = 0.;
#line 445 "dbdsqr.f"
	}
#line 447 "dbdsqr.f"
	if (abse <= thresh) {
#line 447 "dbdsqr.f"
	    goto L80;
#line 447 "dbdsqr.f"
	}
#line 449 "dbdsqr.f"
	smin = min(smin,abss);
/* Computing MAX */
#line 450 "dbdsqr.f"
	d__1 = max(smax,abss);
#line 450 "dbdsqr.f"
	smax = max(d__1,abse);
#line 451 "dbdsqr.f"
/* L70: */
#line 451 "dbdsqr.f"
    }
#line 452 "dbdsqr.f"
    ll = 0;
#line 453 "dbdsqr.f"
    goto L90;
#line 454 "dbdsqr.f"
L80:
#line 455 "dbdsqr.f"
    e[ll] = 0.;

/*     Matrix splits since E(LL) = 0 */

#line 459 "dbdsqr.f"
    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop */

#line 463 "dbdsqr.f"
	--m;
#line 464 "dbdsqr.f"
	goto L60;
#line 465 "dbdsqr.f"
    }
#line 466 "dbdsqr.f"
L90:
#line 467 "dbdsqr.f"
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

#line 471 "dbdsqr.f"
    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

#line 475 "dbdsqr.f"
	dlasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
#line 477 "dbdsqr.f"
	d__[m - 1] = sigmx;
#line 478 "dbdsqr.f"
	e[m - 1] = 0.;
#line 479 "dbdsqr.f"
	d__[m] = sigmn;

/*        Compute singular vectors, if desired */

#line 483 "dbdsqr.f"
	if (*ncvt > 0) {
#line 483 "dbdsqr.f"
	    drot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
#line 483 "dbdsqr.f"
	}
#line 486 "dbdsqr.f"
	if (*nru > 0) {
#line 486 "dbdsqr.f"
	    drot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
#line 486 "dbdsqr.f"
	}
#line 488 "dbdsqr.f"
	if (*ncc > 0) {
#line 488 "dbdsqr.f"
	    drot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
#line 488 "dbdsqr.f"
	}
#line 491 "dbdsqr.f"
	m += -2;
#line 492 "dbdsqr.f"
	goto L60;
#line 493 "dbdsqr.f"
    }

/*     If working on new submatrix, choose shift direction */
/*     (from larger end diagonal element towards smaller) */

#line 498 "dbdsqr.f"
    if (ll > oldm || m < oldll) {
#line 499 "dbdsqr.f"
	if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {

/*           Chase bulge from top (big end) to bottom (small end) */

#line 503 "dbdsqr.f"
	    idir = 1;
#line 504 "dbdsqr.f"
	} else {

/*           Chase bulge from bottom (big end) to top (small end) */

#line 508 "dbdsqr.f"
	    idir = 2;
#line 509 "dbdsqr.f"
	}
#line 510 "dbdsqr.f"
    }

/*     Apply convergence tests */

#line 514 "dbdsqr.f"
    if (idir == 1) {

/*        Run convergence test in forward direction */
/*        First apply standard test to bottom of matrix */

#line 519 "dbdsqr.f"
	if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(
		d__1)) || tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) 
		{
#line 521 "dbdsqr.f"
	    e[m - 1] = 0.;
#line 522 "dbdsqr.f"
	    goto L60;
#line 523 "dbdsqr.f"
	}

#line 525 "dbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion forward */

#line 530 "dbdsqr.f"
	    mu = (d__1 = d__[ll], abs(d__1));
#line 531 "dbdsqr.f"
	    sminl = mu;
#line 532 "dbdsqr.f"
	    i__1 = m - 1;
#line 532 "dbdsqr.f"
	    for (lll = ll; lll <= i__1; ++lll) {
#line 533 "dbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 534 "dbdsqr.f"
		    e[lll] = 0.;
#line 535 "dbdsqr.f"
		    goto L60;
#line 536 "dbdsqr.f"
		}
#line 537 "dbdsqr.f"
		mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[
			lll], abs(d__1))));
#line 538 "dbdsqr.f"
		sminl = min(sminl,mu);
#line 539 "dbdsqr.f"
/* L100: */
#line 539 "dbdsqr.f"
	    }
#line 540 "dbdsqr.f"
	}

#line 542 "dbdsqr.f"
    } else {

/*        Run convergence test in backward direction */
/*        First apply standard test to top of matrix */

#line 547 "dbdsqr.f"
	if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)
		) || tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
#line 549 "dbdsqr.f"
	    e[ll] = 0.;
#line 550 "dbdsqr.f"
	    goto L60;
#line 551 "dbdsqr.f"
	}

#line 553 "dbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion backward */

#line 558 "dbdsqr.f"
	    mu = (d__1 = d__[m], abs(d__1));
#line 559 "dbdsqr.f"
	    sminl = mu;
#line 560 "dbdsqr.f"
	    i__1 = ll;
#line 560 "dbdsqr.f"
	    for (lll = m - 1; lll >= i__1; --lll) {
#line 561 "dbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 562 "dbdsqr.f"
		    e[lll] = 0.;
#line 563 "dbdsqr.f"
		    goto L60;
#line 564 "dbdsqr.f"
		}
#line 565 "dbdsqr.f"
		mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll]
			, abs(d__1))));
#line 566 "dbdsqr.f"
		sminl = min(sminl,mu);
#line 567 "dbdsqr.f"
/* L110: */
#line 567 "dbdsqr.f"
	    }
#line 568 "dbdsqr.f"
	}
#line 569 "dbdsqr.f"
    }
#line 570 "dbdsqr.f"
    oldll = ll;
#line 571 "dbdsqr.f"
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative */
/*     accuracy, and if so set the shift to zero. */

/* Computing MAX */
#line 576 "dbdsqr.f"
    d__1 = eps, d__2 = tol * .01;
#line 576 "dbdsqr.f"
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

#line 581 "dbdsqr.f"
	shift = 0.;
#line 582 "dbdsqr.f"
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

#line 586 "dbdsqr.f"
	if (idir == 1) {
#line 587 "dbdsqr.f"
	    sll = (d__1 = d__[ll], abs(d__1));
#line 588 "dbdsqr.f"
	    dlas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
#line 589 "dbdsqr.f"
	} else {
#line 590 "dbdsqr.f"
	    sll = (d__1 = d__[m], abs(d__1));
#line 591 "dbdsqr.f"
	    dlas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
#line 592 "dbdsqr.f"
	}

/*        Test if shift negligible, and if so set to zero */

#line 596 "dbdsqr.f"
	if (sll > 0.) {
/* Computing 2nd power */
#line 597 "dbdsqr.f"
	    d__1 = shift / sll;
#line 597 "dbdsqr.f"
	    if (d__1 * d__1 < eps) {
#line 597 "dbdsqr.f"
		shift = 0.;
#line 597 "dbdsqr.f"
	    }
#line 599 "dbdsqr.f"
	}
#line 600 "dbdsqr.f"
    }

/*     Increment iteration count */

#line 604 "dbdsqr.f"
    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

#line 608 "dbdsqr.f"
    if (shift == 0.) {
#line 609 "dbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 614 "dbdsqr.f"
	    cs = 1.;
#line 615 "dbdsqr.f"
	    oldcs = 1.;
#line 616 "dbdsqr.f"
	    i__1 = m - 1;
#line 616 "dbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 617 "dbdsqr.f"
		d__1 = d__[i__] * cs;
#line 617 "dbdsqr.f"
		dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 618 "dbdsqr.f"
		if (i__ > ll) {
#line 618 "dbdsqr.f"
		    e[i__ - 1] = oldsn * r__;
#line 618 "dbdsqr.f"
		}
#line 620 "dbdsqr.f"
		d__1 = oldcs * r__;
#line 620 "dbdsqr.f"
		d__2 = d__[i__ + 1] * sn;
#line 620 "dbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 621 "dbdsqr.f"
		work[i__ - ll + 1] = cs;
#line 622 "dbdsqr.f"
		work[i__ - ll + 1 + nm1] = sn;
#line 623 "dbdsqr.f"
		work[i__ - ll + 1 + nm12] = oldcs;
#line 624 "dbdsqr.f"
		work[i__ - ll + 1 + nm13] = oldsn;
#line 625 "dbdsqr.f"
/* L120: */
#line 625 "dbdsqr.f"
	    }
#line 626 "dbdsqr.f"
	    h__ = d__[m] * cs;
#line 627 "dbdsqr.f"
	    d__[m] = h__ * oldcs;
#line 628 "dbdsqr.f"
	    e[m - 1] = h__ * oldsn;

/*           Update singular vectors */

#line 632 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 632 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 632 "dbdsqr.f"
		dlasr_("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 632 "dbdsqr.f"
	    }
#line 635 "dbdsqr.f"
	    if (*nru > 0) {
#line 635 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 635 "dbdsqr.f"
		dlasr_("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 635 "dbdsqr.f"
	    }
#line 638 "dbdsqr.f"
	    if (*ncc > 0) {
#line 638 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 638 "dbdsqr.f"
		dlasr_("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 638 "dbdsqr.f"
	    }

/*           Test convergence */

#line 644 "dbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 644 "dbdsqr.f"
		e[m - 1] = 0.;
#line 644 "dbdsqr.f"
	    }

#line 647 "dbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 652 "dbdsqr.f"
	    cs = 1.;
#line 653 "dbdsqr.f"
	    oldcs = 1.;
#line 654 "dbdsqr.f"
	    i__1 = ll + 1;
#line 654 "dbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 655 "dbdsqr.f"
		d__1 = d__[i__] * cs;
#line 655 "dbdsqr.f"
		dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 656 "dbdsqr.f"
		if (i__ < m) {
#line 656 "dbdsqr.f"
		    e[i__] = oldsn * r__;
#line 656 "dbdsqr.f"
		}
#line 658 "dbdsqr.f"
		d__1 = oldcs * r__;
#line 658 "dbdsqr.f"
		d__2 = d__[i__ - 1] * sn;
#line 658 "dbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 659 "dbdsqr.f"
		work[i__ - ll] = cs;
#line 660 "dbdsqr.f"
		work[i__ - ll + nm1] = -sn;
#line 661 "dbdsqr.f"
		work[i__ - ll + nm12] = oldcs;
#line 662 "dbdsqr.f"
		work[i__ - ll + nm13] = -oldsn;
#line 663 "dbdsqr.f"
/* L130: */
#line 663 "dbdsqr.f"
	    }
#line 664 "dbdsqr.f"
	    h__ = d__[ll] * cs;
#line 665 "dbdsqr.f"
	    d__[ll] = h__ * oldcs;
#line 666 "dbdsqr.f"
	    e[ll] = h__ * oldsn;

/*           Update singular vectors */

#line 670 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 670 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 670 "dbdsqr.f"
		dlasr_("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 670 "dbdsqr.f"
	    }
#line 673 "dbdsqr.f"
	    if (*nru > 0) {
#line 673 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 673 "dbdsqr.f"
		dlasr_("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 673 "dbdsqr.f"
	    }
#line 676 "dbdsqr.f"
	    if (*ncc > 0) {
#line 676 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 676 "dbdsqr.f"
		dlasr_("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 676 "dbdsqr.f"
	    }

/*           Test convergence */

#line 682 "dbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 682 "dbdsqr.f"
		e[ll] = 0.;
#line 682 "dbdsqr.f"
	    }
#line 684 "dbdsqr.f"
	}
#line 685 "dbdsqr.f"
    } else {

/*        Use nonzero shift */

#line 689 "dbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 694 "dbdsqr.f"
	    f = ((d__1 = d__[ll], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[
		    ll]) + shift / d__[ll]);
#line 696 "dbdsqr.f"
	    g = e[ll];
#line 697 "dbdsqr.f"
	    i__1 = m - 1;
#line 697 "dbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 698 "dbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 699 "dbdsqr.f"
		if (i__ > ll) {
#line 699 "dbdsqr.f"
		    e[i__ - 1] = r__;
#line 699 "dbdsqr.f"
		}
#line 701 "dbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 702 "dbdsqr.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 703 "dbdsqr.f"
		g = sinr * d__[i__ + 1];
#line 704 "dbdsqr.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 705 "dbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 706 "dbdsqr.f"
		d__[i__] = r__;
#line 707 "dbdsqr.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 708 "dbdsqr.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 709 "dbdsqr.f"
		if (i__ < m - 1) {
#line 710 "dbdsqr.f"
		    g = sinl * e[i__ + 1];
#line 711 "dbdsqr.f"
		    e[i__ + 1] = cosl * e[i__ + 1];
#line 712 "dbdsqr.f"
		}
#line 713 "dbdsqr.f"
		work[i__ - ll + 1] = cosr;
#line 714 "dbdsqr.f"
		work[i__ - ll + 1 + nm1] = sinr;
#line 715 "dbdsqr.f"
		work[i__ - ll + 1 + nm12] = cosl;
#line 716 "dbdsqr.f"
		work[i__ - ll + 1 + nm13] = sinl;
#line 717 "dbdsqr.f"
/* L140: */
#line 717 "dbdsqr.f"
	    }
#line 718 "dbdsqr.f"
	    e[m - 1] = f;

/*           Update singular vectors */

#line 722 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 722 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 722 "dbdsqr.f"
		dlasr_("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 722 "dbdsqr.f"
	    }
#line 725 "dbdsqr.f"
	    if (*nru > 0) {
#line 725 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 725 "dbdsqr.f"
		dlasr_("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 725 "dbdsqr.f"
	    }
#line 728 "dbdsqr.f"
	    if (*ncc > 0) {
#line 728 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 728 "dbdsqr.f"
		dlasr_("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 728 "dbdsqr.f"
	    }

/*           Test convergence */

#line 734 "dbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 734 "dbdsqr.f"
		e[m - 1] = 0.;
#line 734 "dbdsqr.f"
	    }

#line 737 "dbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 742 "dbdsqr.f"
	    f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[m]
		    ) + shift / d__[m]);
#line 744 "dbdsqr.f"
	    g = e[m - 1];
#line 745 "dbdsqr.f"
	    i__1 = ll + 1;
#line 745 "dbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 746 "dbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 747 "dbdsqr.f"
		if (i__ < m) {
#line 747 "dbdsqr.f"
		    e[i__] = r__;
#line 747 "dbdsqr.f"
		}
#line 749 "dbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 750 "dbdsqr.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 751 "dbdsqr.f"
		g = sinr * d__[i__ - 1];
#line 752 "dbdsqr.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 753 "dbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 754 "dbdsqr.f"
		d__[i__] = r__;
#line 755 "dbdsqr.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 756 "dbdsqr.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 757 "dbdsqr.f"
		if (i__ > ll + 1) {
#line 758 "dbdsqr.f"
		    g = sinl * e[i__ - 2];
#line 759 "dbdsqr.f"
		    e[i__ - 2] = cosl * e[i__ - 2];
#line 760 "dbdsqr.f"
		}
#line 761 "dbdsqr.f"
		work[i__ - ll] = cosr;
#line 762 "dbdsqr.f"
		work[i__ - ll + nm1] = -sinr;
#line 763 "dbdsqr.f"
		work[i__ - ll + nm12] = cosl;
#line 764 "dbdsqr.f"
		work[i__ - ll + nm13] = -sinl;
#line 765 "dbdsqr.f"
/* L150: */
#line 765 "dbdsqr.f"
	    }
#line 766 "dbdsqr.f"
	    e[ll] = f;

/*           Test convergence */

#line 770 "dbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 770 "dbdsqr.f"
		e[ll] = 0.;
#line 770 "dbdsqr.f"
	    }

/*           Update singular vectors if desired */

#line 775 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 775 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 775 "dbdsqr.f"
		dlasr_("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 775 "dbdsqr.f"
	    }
#line 778 "dbdsqr.f"
	    if (*nru > 0) {
#line 778 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 778 "dbdsqr.f"
		dlasr_("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 778 "dbdsqr.f"
	    }
#line 781 "dbdsqr.f"
	    if (*ncc > 0) {
#line 781 "dbdsqr.f"
		i__1 = m - ll + 1;
#line 781 "dbdsqr.f"
		dlasr_("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 781 "dbdsqr.f"
	    }
#line 784 "dbdsqr.f"
	}
#line 785 "dbdsqr.f"
    }

/*     QR iteration finished, go back and check convergence */

#line 789 "dbdsqr.f"
    goto L60;

/*     All singular values converged, so make them positive */

#line 793 "dbdsqr.f"
L160:
#line 794 "dbdsqr.f"
    i__1 = *n;
#line 794 "dbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 795 "dbdsqr.f"
	if (d__[i__] < 0.) {
#line 796 "dbdsqr.f"
	    d__[i__] = -d__[i__];

/*           Change sign of singular vectors, if desired */

#line 800 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 800 "dbdsqr.f"
		dscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
#line 800 "dbdsqr.f"
	    }
#line 802 "dbdsqr.f"
	}
#line 803 "dbdsqr.f"
/* L170: */
#line 803 "dbdsqr.f"
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 808 "dbdsqr.f"
    i__1 = *n - 1;
#line 808 "dbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I) */

#line 812 "dbdsqr.f"
	isub = 1;
#line 813 "dbdsqr.f"
	smin = d__[1];
#line 814 "dbdsqr.f"
	i__2 = *n + 1 - i__;
#line 814 "dbdsqr.f"
	for (j = 2; j <= i__2; ++j) {
#line 815 "dbdsqr.f"
	    if (d__[j] <= smin) {
#line 816 "dbdsqr.f"
		isub = j;
#line 817 "dbdsqr.f"
		smin = d__[j];
#line 818 "dbdsqr.f"
	    }
#line 819 "dbdsqr.f"
/* L180: */
#line 819 "dbdsqr.f"
	}
#line 820 "dbdsqr.f"
	if (isub != *n + 1 - i__) {

/*           Swap singular values and vectors */

#line 824 "dbdsqr.f"
	    d__[isub] = d__[*n + 1 - i__];
#line 825 "dbdsqr.f"
	    d__[*n + 1 - i__] = smin;
#line 826 "dbdsqr.f"
	    if (*ncvt > 0) {
#line 826 "dbdsqr.f"
		dswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
#line 826 "dbdsqr.f"
	    }
#line 829 "dbdsqr.f"
	    if (*nru > 0) {
#line 829 "dbdsqr.f"
		dswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
#line 829 "dbdsqr.f"
	    }
#line 831 "dbdsqr.f"
	    if (*ncc > 0) {
#line 831 "dbdsqr.f"
		dswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
#line 831 "dbdsqr.f"
	    }
#line 833 "dbdsqr.f"
	}
#line 834 "dbdsqr.f"
/* L190: */
#line 834 "dbdsqr.f"
    }
#line 835 "dbdsqr.f"
    goto L220;

/*     Maximum number of iterations exceeded, failure to converge */

#line 839 "dbdsqr.f"
L200:
#line 840 "dbdsqr.f"
    *info = 0;
#line 841 "dbdsqr.f"
    i__1 = *n - 1;
#line 841 "dbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 842 "dbdsqr.f"
	if (e[i__] != 0.) {
#line 842 "dbdsqr.f"
	    ++(*info);
#line 842 "dbdsqr.f"
	}
#line 844 "dbdsqr.f"
/* L210: */
#line 844 "dbdsqr.f"
    }
#line 845 "dbdsqr.f"
L220:
#line 846 "dbdsqr.f"
    return 0;

/*     End of DBDSQR */

} /* dbdsqr_ */

