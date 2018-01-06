#line 1 "zbdsqr.f"
/* zbdsqr.f -- translated by f2c (version 20100827).
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

#line 1 "zbdsqr.f"
/* Table of constant values */

static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;

/* > \brief \b ZBDSQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZBDSQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zbdsqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zbdsqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zbdsqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/*                          LDU, C, LDC, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * ) */
/*       COMPLEX*16         C( LDC, * ), U( LDU, * ), VT( LDVT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZBDSQR computes the singular values and, optionally, the right and/or */
/* > left singular vectors from the singular value decomposition (SVD) of */
/* > a real N-by-N (upper or lower) bidiagonal matrix B using the implicit */
/* > zero-shift QR algorithm.  The SVD of B has the form */
/* > */
/* >    B = Q * S * P**H */
/* > */
/* > where S is the diagonal matrix of singular values, Q is an orthogonal */
/* > matrix of left singular vectors, and P is an orthogonal matrix of */
/* > right singular vectors.  If left singular vectors are requested, this */
/* > subroutine actually returns U*Q instead of Q, and, if right singular */
/* > vectors are requested, this subroutine returns P**H*VT instead of */
/* > P**H, for given complex input matrices U and VT.  When U and VT are */
/* > the unitary matrices that reduce a general matrix A to bidiagonal */
/* > form: A = U*B*VT, as computed by ZGEBRD, then */
/* > */
/* >    A = (U*Q) * S * (P**H*VT) */
/* > */
/* > is the SVD of A.  Optionally, the subroutine may also compute Q**H*C */
/* > for a given complex input matrix C. */
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
/* >          VT is COMPLEX*16 array, dimension (LDVT, NCVT) */
/* >          On entry, an N-by-NCVT matrix VT. */
/* >          On exit, VT is overwritten by P**H * VT. */
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
/* >          U is COMPLEX*16 array, dimension (LDU, N) */
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
/* >          C is COMPLEX*16 array, dimension (LDC, NCC) */
/* >          On entry, an N-by-NCC matrix C. */
/* >          On exit, C is overwritten by Q**H * C. */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* >          if NCVT = NRU = NCC = 0, (max(1, 4*N-4)) otherwise */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  If INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  the algorithm did not converge; D and E contain the */
/* >                elements of a bidiagonal matrix which is orthogonally */
/* >                similar to the input matrix B;  if INFO = i, i */
/* >                elements of E have not converged to zero. */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublecomplex *vt, 
	integer *ldvt, doublecomplex *u, integer *ldu, doublecomplex *c__, 
	integer *ldc, doublereal *rwork, integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int dlas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal oldcs;
    static integer oldll;
    static doublereal shift, sigmn, oldsn;
    static integer maxit;
    static doublereal sminl, sigmx;
    static logical lower;
    extern /* Subroutine */ int zlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), zdrot_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *)
	    , zswap_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), dlasq1_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen), zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static doublereal sminoa, thresh;
    static logical rotate;
    static doublereal tolmul;


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 285 "zbdsqr.f"
    /* Parameter adjustments */
#line 285 "zbdsqr.f"
    --d__;
#line 285 "zbdsqr.f"
    --e;
#line 285 "zbdsqr.f"
    vt_dim1 = *ldvt;
#line 285 "zbdsqr.f"
    vt_offset = 1 + vt_dim1;
#line 285 "zbdsqr.f"
    vt -= vt_offset;
#line 285 "zbdsqr.f"
    u_dim1 = *ldu;
#line 285 "zbdsqr.f"
    u_offset = 1 + u_dim1;
#line 285 "zbdsqr.f"
    u -= u_offset;
#line 285 "zbdsqr.f"
    c_dim1 = *ldc;
#line 285 "zbdsqr.f"
    c_offset = 1 + c_dim1;
#line 285 "zbdsqr.f"
    c__ -= c_offset;
#line 285 "zbdsqr.f"
    --rwork;
#line 285 "zbdsqr.f"

#line 285 "zbdsqr.f"
    /* Function Body */
#line 285 "zbdsqr.f"
    *info = 0;
#line 286 "zbdsqr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 287 "zbdsqr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
#line 288 "zbdsqr.f"
	*info = -1;
#line 289 "zbdsqr.f"
    } else if (*n < 0) {
#line 290 "zbdsqr.f"
	*info = -2;
#line 291 "zbdsqr.f"
    } else if (*ncvt < 0) {
#line 292 "zbdsqr.f"
	*info = -3;
#line 293 "zbdsqr.f"
    } else if (*nru < 0) {
#line 294 "zbdsqr.f"
	*info = -4;
#line 295 "zbdsqr.f"
    } else if (*ncc < 0) {
#line 296 "zbdsqr.f"
	*info = -5;
#line 297 "zbdsqr.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 299 "zbdsqr.f"
	*info = -9;
#line 300 "zbdsqr.f"
    } else if (*ldu < max(1,*nru)) {
#line 301 "zbdsqr.f"
	*info = -11;
#line 302 "zbdsqr.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 304 "zbdsqr.f"
	*info = -13;
#line 305 "zbdsqr.f"
    }
#line 306 "zbdsqr.f"
    if (*info != 0) {
#line 307 "zbdsqr.f"
	i__1 = -(*info);
#line 307 "zbdsqr.f"
	xerbla_("ZBDSQR", &i__1, (ftnlen)6);
#line 308 "zbdsqr.f"
	return 0;
#line 309 "zbdsqr.f"
    }
#line 310 "zbdsqr.f"
    if (*n == 0) {
#line 310 "zbdsqr.f"
	return 0;
#line 310 "zbdsqr.f"
    }
#line 312 "zbdsqr.f"
    if (*n == 1) {
#line 312 "zbdsqr.f"
	goto L160;
#line 312 "zbdsqr.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 317 "zbdsqr.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

#line 321 "zbdsqr.f"
    if (! rotate) {
#line 322 "zbdsqr.f"
	dlasq1_(n, &d__[1], &e[1], &rwork[1], info);

/*     If INFO equals 2, dqds didn't finish, try to finish */

#line 326 "zbdsqr.f"
	if (*info != 2) {
#line 326 "zbdsqr.f"
	    return 0;
#line 326 "zbdsqr.f"
	}
#line 327 "zbdsqr.f"
	*info = 0;
#line 328 "zbdsqr.f"
    }

#line 330 "zbdsqr.f"
    nm1 = *n - 1;
#line 331 "zbdsqr.f"
    nm12 = nm1 + nm1;
#line 332 "zbdsqr.f"
    nm13 = nm12 + nm1;
#line 333 "zbdsqr.f"
    idir = 0;

/*     Get machine constants */

#line 337 "zbdsqr.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 338 "zbdsqr.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 343 "zbdsqr.f"
    if (lower) {
#line 344 "zbdsqr.f"
	i__1 = *n - 1;
#line 344 "zbdsqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "zbdsqr.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 346 "zbdsqr.f"
	    d__[i__] = r__;
#line 347 "zbdsqr.f"
	    e[i__] = sn * d__[i__ + 1];
#line 348 "zbdsqr.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 349 "zbdsqr.f"
	    rwork[i__] = cs;
#line 350 "zbdsqr.f"
	    rwork[nm1 + i__] = sn;
#line 351 "zbdsqr.f"
/* L10: */
#line 351 "zbdsqr.f"
	}

/*        Update singular vectors if desired */

#line 355 "zbdsqr.f"
	if (*nru > 0) {
#line 355 "zbdsqr.f"
	    zlasr_("R", "V", "F", nru, n, &rwork[1], &rwork[*n], &u[u_offset],
		     ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 355 "zbdsqr.f"
	}
#line 358 "zbdsqr.f"
	if (*ncc > 0) {
#line 358 "zbdsqr.f"
	    zlasr_("L", "V", "F", n, ncc, &rwork[1], &rwork[*n], &c__[
		    c_offset], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 358 "zbdsqr.f"
	}
#line 361 "zbdsqr.f"
    }

/*     Compute singular values to relative accuracy TOL */
/*     (By setting TOL to be negative, algorithm will compute */
/*     singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */

/* Computing MAX */
/* Computing MIN */
#line 367 "zbdsqr.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b15);
#line 367 "zbdsqr.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 367 "zbdsqr.f"
    tolmul = max(d__1,d__2);
#line 368 "zbdsqr.f"
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

#line 372 "zbdsqr.f"
    smax = 0.;
#line 373 "zbdsqr.f"
    i__1 = *n;
#line 373 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 374 "zbdsqr.f"
	d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
#line 374 "zbdsqr.f"
	smax = max(d__2,d__3);
#line 375 "zbdsqr.f"
/* L20: */
#line 375 "zbdsqr.f"
    }
#line 376 "zbdsqr.f"
    i__1 = *n - 1;
#line 376 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 377 "zbdsqr.f"
	d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
#line 377 "zbdsqr.f"
	smax = max(d__2,d__3);
#line 378 "zbdsqr.f"
/* L30: */
#line 378 "zbdsqr.f"
    }
#line 379 "zbdsqr.f"
    sminl = 0.;
#line 380 "zbdsqr.f"
    if (tol >= 0.) {

/*        Relative accuracy desired */

#line 384 "zbdsqr.f"
	sminoa = abs(d__[1]);
#line 385 "zbdsqr.f"
	if (sminoa == 0.) {
#line 385 "zbdsqr.f"
	    goto L50;
#line 385 "zbdsqr.f"
	}
#line 387 "zbdsqr.f"
	mu = sminoa;
#line 388 "zbdsqr.f"
	i__1 = *n;
#line 388 "zbdsqr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 389 "zbdsqr.f"
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
#line 390 "zbdsqr.f"
	    sminoa = min(sminoa,mu);
#line 391 "zbdsqr.f"
	    if (sminoa == 0.) {
#line 391 "zbdsqr.f"
		goto L50;
#line 391 "zbdsqr.f"
	    }
#line 393 "zbdsqr.f"
/* L40: */
#line 393 "zbdsqr.f"
	}
#line 394 "zbdsqr.f"
L50:
#line 395 "zbdsqr.f"
	sminoa /= sqrt((doublereal) (*n));
/* Computing MAX */
#line 396 "zbdsqr.f"
	d__1 = tol * sminoa, d__2 = *n * 6 * *n * unfl;
#line 396 "zbdsqr.f"
	thresh = max(d__1,d__2);
#line 397 "zbdsqr.f"
    } else {

/*        Absolute accuracy desired */

/* Computing MAX */
#line 401 "zbdsqr.f"
	d__1 = abs(tol) * smax, d__2 = *n * 6 * *n * unfl;
#line 401 "zbdsqr.f"
	thresh = max(d__1,d__2);
#line 402 "zbdsqr.f"
    }

/*     Prepare for main iteration loop for the singular values */
/*     (MAXIT is the maximum number of passes through the inner */
/*     loop permitted before nonconvergence signalled.) */

#line 408 "zbdsqr.f"
    maxit = *n * 6 * *n;
#line 409 "zbdsqr.f"
    iter = 0;
#line 410 "zbdsqr.f"
    oldll = -1;
#line 411 "zbdsqr.f"
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

#line 415 "zbdsqr.f"
    m = *n;

/*     Begin main iteration loop */

#line 419 "zbdsqr.f"
L60:

/*     Check for convergence or exceeding iteration count */

#line 423 "zbdsqr.f"
    if (m <= 1) {
#line 423 "zbdsqr.f"
	goto L160;
#line 423 "zbdsqr.f"
    }
#line 425 "zbdsqr.f"
    if (iter > maxit) {
#line 425 "zbdsqr.f"
	goto L200;
#line 425 "zbdsqr.f"
    }

/*     Find diagonal block of matrix to work on */

#line 430 "zbdsqr.f"
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
#line 430 "zbdsqr.f"
	d__[m] = 0.;
#line 430 "zbdsqr.f"
    }
#line 432 "zbdsqr.f"
    smax = (d__1 = d__[m], abs(d__1));
#line 433 "zbdsqr.f"
    smin = smax;
#line 434 "zbdsqr.f"
    i__1 = m - 1;
#line 434 "zbdsqr.f"
    for (lll = 1; lll <= i__1; ++lll) {
#line 435 "zbdsqr.f"
	ll = m - lll;
#line 436 "zbdsqr.f"
	abss = (d__1 = d__[ll], abs(d__1));
#line 437 "zbdsqr.f"
	abse = (d__1 = e[ll], abs(d__1));
#line 438 "zbdsqr.f"
	if (tol < 0. && abss <= thresh) {
#line 438 "zbdsqr.f"
	    d__[ll] = 0.;
#line 438 "zbdsqr.f"
	}
#line 440 "zbdsqr.f"
	if (abse <= thresh) {
#line 440 "zbdsqr.f"
	    goto L80;
#line 440 "zbdsqr.f"
	}
#line 442 "zbdsqr.f"
	smin = min(smin,abss);
/* Computing MAX */
#line 443 "zbdsqr.f"
	d__1 = max(smax,abss);
#line 443 "zbdsqr.f"
	smax = max(d__1,abse);
#line 444 "zbdsqr.f"
/* L70: */
#line 444 "zbdsqr.f"
    }
#line 445 "zbdsqr.f"
    ll = 0;
#line 446 "zbdsqr.f"
    goto L90;
#line 447 "zbdsqr.f"
L80:
#line 448 "zbdsqr.f"
    e[ll] = 0.;

/*     Matrix splits since E(LL) = 0 */

#line 452 "zbdsqr.f"
    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop */

#line 456 "zbdsqr.f"
	--m;
#line 457 "zbdsqr.f"
	goto L60;
#line 458 "zbdsqr.f"
    }
#line 459 "zbdsqr.f"
L90:
#line 460 "zbdsqr.f"
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

#line 464 "zbdsqr.f"
    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

#line 468 "zbdsqr.f"
	dlasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
#line 470 "zbdsqr.f"
	d__[m - 1] = sigmx;
#line 471 "zbdsqr.f"
	e[m - 1] = 0.;
#line 472 "zbdsqr.f"
	d__[m] = sigmn;

/*        Compute singular vectors, if desired */

#line 476 "zbdsqr.f"
	if (*ncvt > 0) {
#line 476 "zbdsqr.f"
	    zdrot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
#line 476 "zbdsqr.f"
	}
#line 479 "zbdsqr.f"
	if (*nru > 0) {
#line 479 "zbdsqr.f"
	    zdrot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
#line 479 "zbdsqr.f"
	}
#line 481 "zbdsqr.f"
	if (*ncc > 0) {
#line 481 "zbdsqr.f"
	    zdrot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
#line 481 "zbdsqr.f"
	}
#line 484 "zbdsqr.f"
	m += -2;
#line 485 "zbdsqr.f"
	goto L60;
#line 486 "zbdsqr.f"
    }

/*     If working on new submatrix, choose shift direction */
/*     (from larger end diagonal element towards smaller) */

#line 491 "zbdsqr.f"
    if (ll > oldm || m < oldll) {
#line 492 "zbdsqr.f"
	if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {

/*           Chase bulge from top (big end) to bottom (small end) */

#line 496 "zbdsqr.f"
	    idir = 1;
#line 497 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom (big end) to top (small end) */

#line 501 "zbdsqr.f"
	    idir = 2;
#line 502 "zbdsqr.f"
	}
#line 503 "zbdsqr.f"
    }

/*     Apply convergence tests */

#line 507 "zbdsqr.f"
    if (idir == 1) {

/*        Run convergence test in forward direction */
/*        First apply standard test to bottom of matrix */

#line 512 "zbdsqr.f"
	if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(
		d__1)) || tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) 
		{
#line 514 "zbdsqr.f"
	    e[m - 1] = 0.;
#line 515 "zbdsqr.f"
	    goto L60;
#line 516 "zbdsqr.f"
	}

#line 518 "zbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion forward */

#line 523 "zbdsqr.f"
	    mu = (d__1 = d__[ll], abs(d__1));
#line 524 "zbdsqr.f"
	    sminl = mu;
#line 525 "zbdsqr.f"
	    i__1 = m - 1;
#line 525 "zbdsqr.f"
	    for (lll = ll; lll <= i__1; ++lll) {
#line 526 "zbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 527 "zbdsqr.f"
		    e[lll] = 0.;
#line 528 "zbdsqr.f"
		    goto L60;
#line 529 "zbdsqr.f"
		}
#line 530 "zbdsqr.f"
		mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[
			lll], abs(d__1))));
#line 531 "zbdsqr.f"
		sminl = min(sminl,mu);
#line 532 "zbdsqr.f"
/* L100: */
#line 532 "zbdsqr.f"
	    }
#line 533 "zbdsqr.f"
	}

#line 535 "zbdsqr.f"
    } else {

/*        Run convergence test in backward direction */
/*        First apply standard test to top of matrix */

#line 540 "zbdsqr.f"
	if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)
		) || tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
#line 542 "zbdsqr.f"
	    e[ll] = 0.;
#line 543 "zbdsqr.f"
	    goto L60;
#line 544 "zbdsqr.f"
	}

#line 546 "zbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion backward */

#line 551 "zbdsqr.f"
	    mu = (d__1 = d__[m], abs(d__1));
#line 552 "zbdsqr.f"
	    sminl = mu;
#line 553 "zbdsqr.f"
	    i__1 = ll;
#line 553 "zbdsqr.f"
	    for (lll = m - 1; lll >= i__1; --lll) {
#line 554 "zbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 555 "zbdsqr.f"
		    e[lll] = 0.;
#line 556 "zbdsqr.f"
		    goto L60;
#line 557 "zbdsqr.f"
		}
#line 558 "zbdsqr.f"
		mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll]
			, abs(d__1))));
#line 559 "zbdsqr.f"
		sminl = min(sminl,mu);
#line 560 "zbdsqr.f"
/* L110: */
#line 560 "zbdsqr.f"
	    }
#line 561 "zbdsqr.f"
	}
#line 562 "zbdsqr.f"
    }
#line 563 "zbdsqr.f"
    oldll = ll;
#line 564 "zbdsqr.f"
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative */
/*     accuracy, and if so set the shift to zero. */

/* Computing MAX */
#line 569 "zbdsqr.f"
    d__1 = eps, d__2 = tol * .01;
#line 569 "zbdsqr.f"
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

#line 574 "zbdsqr.f"
	shift = 0.;
#line 575 "zbdsqr.f"
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

#line 579 "zbdsqr.f"
	if (idir == 1) {
#line 580 "zbdsqr.f"
	    sll = (d__1 = d__[ll], abs(d__1));
#line 581 "zbdsqr.f"
	    dlas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
#line 582 "zbdsqr.f"
	} else {
#line 583 "zbdsqr.f"
	    sll = (d__1 = d__[m], abs(d__1));
#line 584 "zbdsqr.f"
	    dlas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
#line 585 "zbdsqr.f"
	}

/*        Test if shift negligible, and if so set to zero */

#line 589 "zbdsqr.f"
	if (sll > 0.) {
/* Computing 2nd power */
#line 590 "zbdsqr.f"
	    d__1 = shift / sll;
#line 590 "zbdsqr.f"
	    if (d__1 * d__1 < eps) {
#line 590 "zbdsqr.f"
		shift = 0.;
#line 590 "zbdsqr.f"
	    }
#line 592 "zbdsqr.f"
	}
#line 593 "zbdsqr.f"
    }

/*     Increment iteration count */

#line 597 "zbdsqr.f"
    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

#line 601 "zbdsqr.f"
    if (shift == 0.) {
#line 602 "zbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 607 "zbdsqr.f"
	    cs = 1.;
#line 608 "zbdsqr.f"
	    oldcs = 1.;
#line 609 "zbdsqr.f"
	    i__1 = m - 1;
#line 609 "zbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 610 "zbdsqr.f"
		d__1 = d__[i__] * cs;
#line 610 "zbdsqr.f"
		dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 611 "zbdsqr.f"
		if (i__ > ll) {
#line 611 "zbdsqr.f"
		    e[i__ - 1] = oldsn * r__;
#line 611 "zbdsqr.f"
		}
#line 613 "zbdsqr.f"
		d__1 = oldcs * r__;
#line 613 "zbdsqr.f"
		d__2 = d__[i__ + 1] * sn;
#line 613 "zbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 614 "zbdsqr.f"
		rwork[i__ - ll + 1] = cs;
#line 615 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sn;
#line 616 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = oldcs;
#line 617 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = oldsn;
#line 618 "zbdsqr.f"
/* L120: */
#line 618 "zbdsqr.f"
	    }
#line 619 "zbdsqr.f"
	    h__ = d__[m] * cs;
#line 620 "zbdsqr.f"
	    d__[m] = h__ * oldcs;
#line 621 "zbdsqr.f"
	    e[m - 1] = h__ * oldsn;

/*           Update singular vectors */

#line 625 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 625 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 625 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 625 "zbdsqr.f"
	    }
#line 628 "zbdsqr.f"
	    if (*nru > 0) {
#line 628 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 628 "zbdsqr.f"
		zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 628 "zbdsqr.f"
	    }
#line 631 "zbdsqr.f"
	    if (*ncc > 0) {
#line 631 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 631 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 631 "zbdsqr.f"
	    }

/*           Test convergence */

#line 637 "zbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 637 "zbdsqr.f"
		e[m - 1] = 0.;
#line 637 "zbdsqr.f"
	    }

#line 640 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 645 "zbdsqr.f"
	    cs = 1.;
#line 646 "zbdsqr.f"
	    oldcs = 1.;
#line 647 "zbdsqr.f"
	    i__1 = ll + 1;
#line 647 "zbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 648 "zbdsqr.f"
		d__1 = d__[i__] * cs;
#line 648 "zbdsqr.f"
		dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 649 "zbdsqr.f"
		if (i__ < m) {
#line 649 "zbdsqr.f"
		    e[i__] = oldsn * r__;
#line 649 "zbdsqr.f"
		}
#line 651 "zbdsqr.f"
		d__1 = oldcs * r__;
#line 651 "zbdsqr.f"
		d__2 = d__[i__ - 1] * sn;
#line 651 "zbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 652 "zbdsqr.f"
		rwork[i__ - ll] = cs;
#line 653 "zbdsqr.f"
		rwork[i__ - ll + nm1] = -sn;
#line 654 "zbdsqr.f"
		rwork[i__ - ll + nm12] = oldcs;
#line 655 "zbdsqr.f"
		rwork[i__ - ll + nm13] = -oldsn;
#line 656 "zbdsqr.f"
/* L130: */
#line 656 "zbdsqr.f"
	    }
#line 657 "zbdsqr.f"
	    h__ = d__[ll] * cs;
#line 658 "zbdsqr.f"
	    d__[ll] = h__ * oldcs;
#line 659 "zbdsqr.f"
	    e[ll] = h__ * oldsn;

/*           Update singular vectors */

#line 663 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 663 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 663 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 663 "zbdsqr.f"
	    }
#line 666 "zbdsqr.f"
	    if (*nru > 0) {
#line 666 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 666 "zbdsqr.f"
		zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 666 "zbdsqr.f"
	    }
#line 669 "zbdsqr.f"
	    if (*ncc > 0) {
#line 669 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 669 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 669 "zbdsqr.f"
	    }

/*           Test convergence */

#line 675 "zbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 675 "zbdsqr.f"
		e[ll] = 0.;
#line 675 "zbdsqr.f"
	    }
#line 677 "zbdsqr.f"
	}
#line 678 "zbdsqr.f"
    } else {

/*        Use nonzero shift */

#line 682 "zbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 687 "zbdsqr.f"
	    f = ((d__1 = d__[ll], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[
		    ll]) + shift / d__[ll]);
#line 689 "zbdsqr.f"
	    g = e[ll];
#line 690 "zbdsqr.f"
	    i__1 = m - 1;
#line 690 "zbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 691 "zbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 692 "zbdsqr.f"
		if (i__ > ll) {
#line 692 "zbdsqr.f"
		    e[i__ - 1] = r__;
#line 692 "zbdsqr.f"
		}
#line 694 "zbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 695 "zbdsqr.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 696 "zbdsqr.f"
		g = sinr * d__[i__ + 1];
#line 697 "zbdsqr.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 698 "zbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 699 "zbdsqr.f"
		d__[i__] = r__;
#line 700 "zbdsqr.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 701 "zbdsqr.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 702 "zbdsqr.f"
		if (i__ < m - 1) {
#line 703 "zbdsqr.f"
		    g = sinl * e[i__ + 1];
#line 704 "zbdsqr.f"
		    e[i__ + 1] = cosl * e[i__ + 1];
#line 705 "zbdsqr.f"
		}
#line 706 "zbdsqr.f"
		rwork[i__ - ll + 1] = cosr;
#line 707 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sinr;
#line 708 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = cosl;
#line 709 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = sinl;
#line 710 "zbdsqr.f"
/* L140: */
#line 710 "zbdsqr.f"
	    }
#line 711 "zbdsqr.f"
	    e[m - 1] = f;

/*           Update singular vectors */

#line 715 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 715 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 715 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 715 "zbdsqr.f"
	    }
#line 718 "zbdsqr.f"
	    if (*nru > 0) {
#line 718 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 718 "zbdsqr.f"
		zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 718 "zbdsqr.f"
	    }
#line 721 "zbdsqr.f"
	    if (*ncc > 0) {
#line 721 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 721 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 721 "zbdsqr.f"
	    }

/*           Test convergence */

#line 727 "zbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 727 "zbdsqr.f"
		e[m - 1] = 0.;
#line 727 "zbdsqr.f"
	    }

#line 730 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 735 "zbdsqr.f"
	    f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[m]
		    ) + shift / d__[m]);
#line 737 "zbdsqr.f"
	    g = e[m - 1];
#line 738 "zbdsqr.f"
	    i__1 = ll + 1;
#line 738 "zbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 739 "zbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 740 "zbdsqr.f"
		if (i__ < m) {
#line 740 "zbdsqr.f"
		    e[i__] = r__;
#line 740 "zbdsqr.f"
		}
#line 742 "zbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 743 "zbdsqr.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 744 "zbdsqr.f"
		g = sinr * d__[i__ - 1];
#line 745 "zbdsqr.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 746 "zbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 747 "zbdsqr.f"
		d__[i__] = r__;
#line 748 "zbdsqr.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 749 "zbdsqr.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 750 "zbdsqr.f"
		if (i__ > ll + 1) {
#line 751 "zbdsqr.f"
		    g = sinl * e[i__ - 2];
#line 752 "zbdsqr.f"
		    e[i__ - 2] = cosl * e[i__ - 2];
#line 753 "zbdsqr.f"
		}
#line 754 "zbdsqr.f"
		rwork[i__ - ll] = cosr;
#line 755 "zbdsqr.f"
		rwork[i__ - ll + nm1] = -sinr;
#line 756 "zbdsqr.f"
		rwork[i__ - ll + nm12] = cosl;
#line 757 "zbdsqr.f"
		rwork[i__ - ll + nm13] = -sinl;
#line 758 "zbdsqr.f"
/* L150: */
#line 758 "zbdsqr.f"
	    }
#line 759 "zbdsqr.f"
	    e[ll] = f;

/*           Test convergence */

#line 763 "zbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 763 "zbdsqr.f"
		e[ll] = 0.;
#line 763 "zbdsqr.f"
	    }

/*           Update singular vectors if desired */

#line 768 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 768 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 768 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 768 "zbdsqr.f"
	    }
#line 771 "zbdsqr.f"
	    if (*nru > 0) {
#line 771 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 771 "zbdsqr.f"
		zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 771 "zbdsqr.f"
	    }
#line 774 "zbdsqr.f"
	    if (*ncc > 0) {
#line 774 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 774 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 774 "zbdsqr.f"
	    }
#line 777 "zbdsqr.f"
	}
#line 778 "zbdsqr.f"
    }

/*     QR iteration finished, go back and check convergence */

#line 782 "zbdsqr.f"
    goto L60;

/*     All singular values converged, so make them positive */

#line 786 "zbdsqr.f"
L160:
#line 787 "zbdsqr.f"
    i__1 = *n;
#line 787 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 788 "zbdsqr.f"
	if (d__[i__] < 0.) {
#line 789 "zbdsqr.f"
	    d__[i__] = -d__[i__];

/*           Change sign of singular vectors, if desired */

#line 793 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 793 "zbdsqr.f"
		zdscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
#line 793 "zbdsqr.f"
	    }
#line 795 "zbdsqr.f"
	}
#line 796 "zbdsqr.f"
/* L170: */
#line 796 "zbdsqr.f"
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 801 "zbdsqr.f"
    i__1 = *n - 1;
#line 801 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I) */

#line 805 "zbdsqr.f"
	isub = 1;
#line 806 "zbdsqr.f"
	smin = d__[1];
#line 807 "zbdsqr.f"
	i__2 = *n + 1 - i__;
#line 807 "zbdsqr.f"
	for (j = 2; j <= i__2; ++j) {
#line 808 "zbdsqr.f"
	    if (d__[j] <= smin) {
#line 809 "zbdsqr.f"
		isub = j;
#line 810 "zbdsqr.f"
		smin = d__[j];
#line 811 "zbdsqr.f"
	    }
#line 812 "zbdsqr.f"
/* L180: */
#line 812 "zbdsqr.f"
	}
#line 813 "zbdsqr.f"
	if (isub != *n + 1 - i__) {

/*           Swap singular values and vectors */

#line 817 "zbdsqr.f"
	    d__[isub] = d__[*n + 1 - i__];
#line 818 "zbdsqr.f"
	    d__[*n + 1 - i__] = smin;
#line 819 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 819 "zbdsqr.f"
		zswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
#line 819 "zbdsqr.f"
	    }
#line 822 "zbdsqr.f"
	    if (*nru > 0) {
#line 822 "zbdsqr.f"
		zswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
#line 822 "zbdsqr.f"
	    }
#line 824 "zbdsqr.f"
	    if (*ncc > 0) {
#line 824 "zbdsqr.f"
		zswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
#line 824 "zbdsqr.f"
	    }
#line 826 "zbdsqr.f"
	}
#line 827 "zbdsqr.f"
/* L190: */
#line 827 "zbdsqr.f"
    }
#line 828 "zbdsqr.f"
    goto L220;

/*     Maximum number of iterations exceeded, failure to converge */

#line 832 "zbdsqr.f"
L200:
#line 833 "zbdsqr.f"
    *info = 0;
#line 834 "zbdsqr.f"
    i__1 = *n - 1;
#line 834 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 835 "zbdsqr.f"
	if (e[i__] != 0.) {
#line 835 "zbdsqr.f"
	    ++(*info);
#line 835 "zbdsqr.f"
	}
#line 837 "zbdsqr.f"
/* L210: */
#line 837 "zbdsqr.f"
    }
#line 838 "zbdsqr.f"
L220:
#line 839 "zbdsqr.f"
    return 0;

/*     End of ZBDSQR */

} /* zbdsqr_ */

