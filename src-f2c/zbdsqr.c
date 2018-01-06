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
/* >          RWORK is DOUBLE PRECISION array, dimension (4*N) */
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

/* > \date November 2015 */

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


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 284 "zbdsqr.f"
    /* Parameter adjustments */
#line 284 "zbdsqr.f"
    --d__;
#line 284 "zbdsqr.f"
    --e;
#line 284 "zbdsqr.f"
    vt_dim1 = *ldvt;
#line 284 "zbdsqr.f"
    vt_offset = 1 + vt_dim1;
#line 284 "zbdsqr.f"
    vt -= vt_offset;
#line 284 "zbdsqr.f"
    u_dim1 = *ldu;
#line 284 "zbdsqr.f"
    u_offset = 1 + u_dim1;
#line 284 "zbdsqr.f"
    u -= u_offset;
#line 284 "zbdsqr.f"
    c_dim1 = *ldc;
#line 284 "zbdsqr.f"
    c_offset = 1 + c_dim1;
#line 284 "zbdsqr.f"
    c__ -= c_offset;
#line 284 "zbdsqr.f"
    --rwork;
#line 284 "zbdsqr.f"

#line 284 "zbdsqr.f"
    /* Function Body */
#line 284 "zbdsqr.f"
    *info = 0;
#line 285 "zbdsqr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 286 "zbdsqr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
#line 287 "zbdsqr.f"
	*info = -1;
#line 288 "zbdsqr.f"
    } else if (*n < 0) {
#line 289 "zbdsqr.f"
	*info = -2;
#line 290 "zbdsqr.f"
    } else if (*ncvt < 0) {
#line 291 "zbdsqr.f"
	*info = -3;
#line 292 "zbdsqr.f"
    } else if (*nru < 0) {
#line 293 "zbdsqr.f"
	*info = -4;
#line 294 "zbdsqr.f"
    } else if (*ncc < 0) {
#line 295 "zbdsqr.f"
	*info = -5;
#line 296 "zbdsqr.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 298 "zbdsqr.f"
	*info = -9;
#line 299 "zbdsqr.f"
    } else if (*ldu < max(1,*nru)) {
#line 300 "zbdsqr.f"
	*info = -11;
#line 301 "zbdsqr.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 303 "zbdsqr.f"
	*info = -13;
#line 304 "zbdsqr.f"
    }
#line 305 "zbdsqr.f"
    if (*info != 0) {
#line 306 "zbdsqr.f"
	i__1 = -(*info);
#line 306 "zbdsqr.f"
	xerbla_("ZBDSQR", &i__1, (ftnlen)6);
#line 307 "zbdsqr.f"
	return 0;
#line 308 "zbdsqr.f"
    }
#line 309 "zbdsqr.f"
    if (*n == 0) {
#line 309 "zbdsqr.f"
	return 0;
#line 309 "zbdsqr.f"
    }
#line 311 "zbdsqr.f"
    if (*n == 1) {
#line 311 "zbdsqr.f"
	goto L160;
#line 311 "zbdsqr.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 316 "zbdsqr.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

#line 320 "zbdsqr.f"
    if (! rotate) {
#line 321 "zbdsqr.f"
	dlasq1_(n, &d__[1], &e[1], &rwork[1], info);

/*     If INFO equals 2, dqds didn't finish, try to finish */

#line 325 "zbdsqr.f"
	if (*info != 2) {
#line 325 "zbdsqr.f"
	    return 0;
#line 325 "zbdsqr.f"
	}
#line 326 "zbdsqr.f"
	*info = 0;
#line 327 "zbdsqr.f"
    }

#line 329 "zbdsqr.f"
    nm1 = *n - 1;
#line 330 "zbdsqr.f"
    nm12 = nm1 + nm1;
#line 331 "zbdsqr.f"
    nm13 = nm12 + nm1;
#line 332 "zbdsqr.f"
    idir = 0;

/*     Get machine constants */

#line 336 "zbdsqr.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 337 "zbdsqr.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 342 "zbdsqr.f"
    if (lower) {
#line 343 "zbdsqr.f"
	i__1 = *n - 1;
#line 343 "zbdsqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 344 "zbdsqr.f"
	    dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 345 "zbdsqr.f"
	    d__[i__] = r__;
#line 346 "zbdsqr.f"
	    e[i__] = sn * d__[i__ + 1];
#line 347 "zbdsqr.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 348 "zbdsqr.f"
	    rwork[i__] = cs;
#line 349 "zbdsqr.f"
	    rwork[nm1 + i__] = sn;
#line 350 "zbdsqr.f"
/* L10: */
#line 350 "zbdsqr.f"
	}

/*        Update singular vectors if desired */

#line 354 "zbdsqr.f"
	if (*nru > 0) {
#line 354 "zbdsqr.f"
	    zlasr_("R", "V", "F", nru, n, &rwork[1], &rwork[*n], &u[u_offset],
		     ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 354 "zbdsqr.f"
	}
#line 357 "zbdsqr.f"
	if (*ncc > 0) {
#line 357 "zbdsqr.f"
	    zlasr_("L", "V", "F", n, ncc, &rwork[1], &rwork[*n], &c__[
		    c_offset], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 357 "zbdsqr.f"
	}
#line 360 "zbdsqr.f"
    }

/*     Compute singular values to relative accuracy TOL */
/*     (By setting TOL to be negative, algorithm will compute */
/*     singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */

/* Computing MAX */
/* Computing MIN */
#line 366 "zbdsqr.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b15);
#line 366 "zbdsqr.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 366 "zbdsqr.f"
    tolmul = max(d__1,d__2);
#line 367 "zbdsqr.f"
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

#line 371 "zbdsqr.f"
    smax = 0.;
#line 372 "zbdsqr.f"
    i__1 = *n;
#line 372 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 373 "zbdsqr.f"
	d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
#line 373 "zbdsqr.f"
	smax = max(d__2,d__3);
#line 374 "zbdsqr.f"
/* L20: */
#line 374 "zbdsqr.f"
    }
#line 375 "zbdsqr.f"
    i__1 = *n - 1;
#line 375 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 376 "zbdsqr.f"
	d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
#line 376 "zbdsqr.f"
	smax = max(d__2,d__3);
#line 377 "zbdsqr.f"
/* L30: */
#line 377 "zbdsqr.f"
    }
#line 378 "zbdsqr.f"
    sminl = 0.;
#line 379 "zbdsqr.f"
    if (tol >= 0.) {

/*        Relative accuracy desired */

#line 383 "zbdsqr.f"
	sminoa = abs(d__[1]);
#line 384 "zbdsqr.f"
	if (sminoa == 0.) {
#line 384 "zbdsqr.f"
	    goto L50;
#line 384 "zbdsqr.f"
	}
#line 386 "zbdsqr.f"
	mu = sminoa;
#line 387 "zbdsqr.f"
	i__1 = *n;
#line 387 "zbdsqr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 388 "zbdsqr.f"
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
#line 389 "zbdsqr.f"
	    sminoa = min(sminoa,mu);
#line 390 "zbdsqr.f"
	    if (sminoa == 0.) {
#line 390 "zbdsqr.f"
		goto L50;
#line 390 "zbdsqr.f"
	    }
#line 392 "zbdsqr.f"
/* L40: */
#line 392 "zbdsqr.f"
	}
#line 393 "zbdsqr.f"
L50:
#line 394 "zbdsqr.f"
	sminoa /= sqrt((doublereal) (*n));
/* Computing MAX */
#line 395 "zbdsqr.f"
	d__1 = tol * sminoa, d__2 = *n * 6 * *n * unfl;
#line 395 "zbdsqr.f"
	thresh = max(d__1,d__2);
#line 396 "zbdsqr.f"
    } else {

/*        Absolute accuracy desired */

/* Computing MAX */
#line 400 "zbdsqr.f"
	d__1 = abs(tol) * smax, d__2 = *n * 6 * *n * unfl;
#line 400 "zbdsqr.f"
	thresh = max(d__1,d__2);
#line 401 "zbdsqr.f"
    }

/*     Prepare for main iteration loop for the singular values */
/*     (MAXIT is the maximum number of passes through the inner */
/*     loop permitted before nonconvergence signalled.) */

#line 407 "zbdsqr.f"
    maxit = *n * 6 * *n;
#line 408 "zbdsqr.f"
    iter = 0;
#line 409 "zbdsqr.f"
    oldll = -1;
#line 410 "zbdsqr.f"
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

#line 414 "zbdsqr.f"
    m = *n;

/*     Begin main iteration loop */

#line 418 "zbdsqr.f"
L60:

/*     Check for convergence or exceeding iteration count */

#line 422 "zbdsqr.f"
    if (m <= 1) {
#line 422 "zbdsqr.f"
	goto L160;
#line 422 "zbdsqr.f"
    }
#line 424 "zbdsqr.f"
    if (iter > maxit) {
#line 424 "zbdsqr.f"
	goto L200;
#line 424 "zbdsqr.f"
    }

/*     Find diagonal block of matrix to work on */

#line 429 "zbdsqr.f"
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
#line 429 "zbdsqr.f"
	d__[m] = 0.;
#line 429 "zbdsqr.f"
    }
#line 431 "zbdsqr.f"
    smax = (d__1 = d__[m], abs(d__1));
#line 432 "zbdsqr.f"
    smin = smax;
#line 433 "zbdsqr.f"
    i__1 = m - 1;
#line 433 "zbdsqr.f"
    for (lll = 1; lll <= i__1; ++lll) {
#line 434 "zbdsqr.f"
	ll = m - lll;
#line 435 "zbdsqr.f"
	abss = (d__1 = d__[ll], abs(d__1));
#line 436 "zbdsqr.f"
	abse = (d__1 = e[ll], abs(d__1));
#line 437 "zbdsqr.f"
	if (tol < 0. && abss <= thresh) {
#line 437 "zbdsqr.f"
	    d__[ll] = 0.;
#line 437 "zbdsqr.f"
	}
#line 439 "zbdsqr.f"
	if (abse <= thresh) {
#line 439 "zbdsqr.f"
	    goto L80;
#line 439 "zbdsqr.f"
	}
#line 441 "zbdsqr.f"
	smin = min(smin,abss);
/* Computing MAX */
#line 442 "zbdsqr.f"
	d__1 = max(smax,abss);
#line 442 "zbdsqr.f"
	smax = max(d__1,abse);
#line 443 "zbdsqr.f"
/* L70: */
#line 443 "zbdsqr.f"
    }
#line 444 "zbdsqr.f"
    ll = 0;
#line 445 "zbdsqr.f"
    goto L90;
#line 446 "zbdsqr.f"
L80:
#line 447 "zbdsqr.f"
    e[ll] = 0.;

/*     Matrix splits since E(LL) = 0 */

#line 451 "zbdsqr.f"
    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop */

#line 455 "zbdsqr.f"
	--m;
#line 456 "zbdsqr.f"
	goto L60;
#line 457 "zbdsqr.f"
    }
#line 458 "zbdsqr.f"
L90:
#line 459 "zbdsqr.f"
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

#line 463 "zbdsqr.f"
    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

#line 467 "zbdsqr.f"
	dlasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
#line 469 "zbdsqr.f"
	d__[m - 1] = sigmx;
#line 470 "zbdsqr.f"
	e[m - 1] = 0.;
#line 471 "zbdsqr.f"
	d__[m] = sigmn;

/*        Compute singular vectors, if desired */

#line 475 "zbdsqr.f"
	if (*ncvt > 0) {
#line 475 "zbdsqr.f"
	    zdrot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
#line 475 "zbdsqr.f"
	}
#line 478 "zbdsqr.f"
	if (*nru > 0) {
#line 478 "zbdsqr.f"
	    zdrot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
#line 478 "zbdsqr.f"
	}
#line 480 "zbdsqr.f"
	if (*ncc > 0) {
#line 480 "zbdsqr.f"
	    zdrot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
#line 480 "zbdsqr.f"
	}
#line 483 "zbdsqr.f"
	m += -2;
#line 484 "zbdsqr.f"
	goto L60;
#line 485 "zbdsqr.f"
    }

/*     If working on new submatrix, choose shift direction */
/*     (from larger end diagonal element towards smaller) */

#line 490 "zbdsqr.f"
    if (ll > oldm || m < oldll) {
#line 491 "zbdsqr.f"
	if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {

/*           Chase bulge from top (big end) to bottom (small end) */

#line 495 "zbdsqr.f"
	    idir = 1;
#line 496 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom (big end) to top (small end) */

#line 500 "zbdsqr.f"
	    idir = 2;
#line 501 "zbdsqr.f"
	}
#line 502 "zbdsqr.f"
    }

/*     Apply convergence tests */

#line 506 "zbdsqr.f"
    if (idir == 1) {

/*        Run convergence test in forward direction */
/*        First apply standard test to bottom of matrix */

#line 511 "zbdsqr.f"
	if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(
		d__1)) || tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) 
		{
#line 513 "zbdsqr.f"
	    e[m - 1] = 0.;
#line 514 "zbdsqr.f"
	    goto L60;
#line 515 "zbdsqr.f"
	}

#line 517 "zbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion forward */

#line 522 "zbdsqr.f"
	    mu = (d__1 = d__[ll], abs(d__1));
#line 523 "zbdsqr.f"
	    sminl = mu;
#line 524 "zbdsqr.f"
	    i__1 = m - 1;
#line 524 "zbdsqr.f"
	    for (lll = ll; lll <= i__1; ++lll) {
#line 525 "zbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 526 "zbdsqr.f"
		    e[lll] = 0.;
#line 527 "zbdsqr.f"
		    goto L60;
#line 528 "zbdsqr.f"
		}
#line 529 "zbdsqr.f"
		mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[
			lll], abs(d__1))));
#line 530 "zbdsqr.f"
		sminl = min(sminl,mu);
#line 531 "zbdsqr.f"
/* L100: */
#line 531 "zbdsqr.f"
	    }
#line 532 "zbdsqr.f"
	}

#line 534 "zbdsqr.f"
    } else {

/*        Run convergence test in backward direction */
/*        First apply standard test to top of matrix */

#line 539 "zbdsqr.f"
	if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)
		) || tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
#line 541 "zbdsqr.f"
	    e[ll] = 0.;
#line 542 "zbdsqr.f"
	    goto L60;
#line 543 "zbdsqr.f"
	}

#line 545 "zbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion backward */

#line 550 "zbdsqr.f"
	    mu = (d__1 = d__[m], abs(d__1));
#line 551 "zbdsqr.f"
	    sminl = mu;
#line 552 "zbdsqr.f"
	    i__1 = ll;
#line 552 "zbdsqr.f"
	    for (lll = m - 1; lll >= i__1; --lll) {
#line 553 "zbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 554 "zbdsqr.f"
		    e[lll] = 0.;
#line 555 "zbdsqr.f"
		    goto L60;
#line 556 "zbdsqr.f"
		}
#line 557 "zbdsqr.f"
		mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll]
			, abs(d__1))));
#line 558 "zbdsqr.f"
		sminl = min(sminl,mu);
#line 559 "zbdsqr.f"
/* L110: */
#line 559 "zbdsqr.f"
	    }
#line 560 "zbdsqr.f"
	}
#line 561 "zbdsqr.f"
    }
#line 562 "zbdsqr.f"
    oldll = ll;
#line 563 "zbdsqr.f"
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative */
/*     accuracy, and if so set the shift to zero. */

/* Computing MAX */
#line 568 "zbdsqr.f"
    d__1 = eps, d__2 = tol * .01;
#line 568 "zbdsqr.f"
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

#line 573 "zbdsqr.f"
	shift = 0.;
#line 574 "zbdsqr.f"
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

#line 578 "zbdsqr.f"
	if (idir == 1) {
#line 579 "zbdsqr.f"
	    sll = (d__1 = d__[ll], abs(d__1));
#line 580 "zbdsqr.f"
	    dlas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
#line 581 "zbdsqr.f"
	} else {
#line 582 "zbdsqr.f"
	    sll = (d__1 = d__[m], abs(d__1));
#line 583 "zbdsqr.f"
	    dlas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
#line 584 "zbdsqr.f"
	}

/*        Test if shift negligible, and if so set to zero */

#line 588 "zbdsqr.f"
	if (sll > 0.) {
/* Computing 2nd power */
#line 589 "zbdsqr.f"
	    d__1 = shift / sll;
#line 589 "zbdsqr.f"
	    if (d__1 * d__1 < eps) {
#line 589 "zbdsqr.f"
		shift = 0.;
#line 589 "zbdsqr.f"
	    }
#line 591 "zbdsqr.f"
	}
#line 592 "zbdsqr.f"
    }

/*     Increment iteration count */

#line 596 "zbdsqr.f"
    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

#line 600 "zbdsqr.f"
    if (shift == 0.) {
#line 601 "zbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 606 "zbdsqr.f"
	    cs = 1.;
#line 607 "zbdsqr.f"
	    oldcs = 1.;
#line 608 "zbdsqr.f"
	    i__1 = m - 1;
#line 608 "zbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 609 "zbdsqr.f"
		d__1 = d__[i__] * cs;
#line 609 "zbdsqr.f"
		dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 610 "zbdsqr.f"
		if (i__ > ll) {
#line 610 "zbdsqr.f"
		    e[i__ - 1] = oldsn * r__;
#line 610 "zbdsqr.f"
		}
#line 612 "zbdsqr.f"
		d__1 = oldcs * r__;
#line 612 "zbdsqr.f"
		d__2 = d__[i__ + 1] * sn;
#line 612 "zbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 613 "zbdsqr.f"
		rwork[i__ - ll + 1] = cs;
#line 614 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sn;
#line 615 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = oldcs;
#line 616 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = oldsn;
#line 617 "zbdsqr.f"
/* L120: */
#line 617 "zbdsqr.f"
	    }
#line 618 "zbdsqr.f"
	    h__ = d__[m] * cs;
#line 619 "zbdsqr.f"
	    d__[m] = h__ * oldcs;
#line 620 "zbdsqr.f"
	    e[m - 1] = h__ * oldsn;

/*           Update singular vectors */

#line 624 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 624 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 624 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 624 "zbdsqr.f"
	    }
#line 627 "zbdsqr.f"
	    if (*nru > 0) {
#line 627 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 627 "zbdsqr.f"
		zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 627 "zbdsqr.f"
	    }
#line 630 "zbdsqr.f"
	    if (*ncc > 0) {
#line 630 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 630 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 630 "zbdsqr.f"
	    }

/*           Test convergence */

#line 636 "zbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 636 "zbdsqr.f"
		e[m - 1] = 0.;
#line 636 "zbdsqr.f"
	    }

#line 639 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 644 "zbdsqr.f"
	    cs = 1.;
#line 645 "zbdsqr.f"
	    oldcs = 1.;
#line 646 "zbdsqr.f"
	    i__1 = ll + 1;
#line 646 "zbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 647 "zbdsqr.f"
		d__1 = d__[i__] * cs;
#line 647 "zbdsqr.f"
		dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 648 "zbdsqr.f"
		if (i__ < m) {
#line 648 "zbdsqr.f"
		    e[i__] = oldsn * r__;
#line 648 "zbdsqr.f"
		}
#line 650 "zbdsqr.f"
		d__1 = oldcs * r__;
#line 650 "zbdsqr.f"
		d__2 = d__[i__ - 1] * sn;
#line 650 "zbdsqr.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 651 "zbdsqr.f"
		rwork[i__ - ll] = cs;
#line 652 "zbdsqr.f"
		rwork[i__ - ll + nm1] = -sn;
#line 653 "zbdsqr.f"
		rwork[i__ - ll + nm12] = oldcs;
#line 654 "zbdsqr.f"
		rwork[i__ - ll + nm13] = -oldsn;
#line 655 "zbdsqr.f"
/* L130: */
#line 655 "zbdsqr.f"
	    }
#line 656 "zbdsqr.f"
	    h__ = d__[ll] * cs;
#line 657 "zbdsqr.f"
	    d__[ll] = h__ * oldcs;
#line 658 "zbdsqr.f"
	    e[ll] = h__ * oldsn;

/*           Update singular vectors */

#line 662 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 662 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 662 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 662 "zbdsqr.f"
	    }
#line 665 "zbdsqr.f"
	    if (*nru > 0) {
#line 665 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 665 "zbdsqr.f"
		zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 665 "zbdsqr.f"
	    }
#line 668 "zbdsqr.f"
	    if (*ncc > 0) {
#line 668 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 668 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 668 "zbdsqr.f"
	    }

/*           Test convergence */

#line 674 "zbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 674 "zbdsqr.f"
		e[ll] = 0.;
#line 674 "zbdsqr.f"
	    }
#line 676 "zbdsqr.f"
	}
#line 677 "zbdsqr.f"
    } else {

/*        Use nonzero shift */

#line 681 "zbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 686 "zbdsqr.f"
	    f = ((d__1 = d__[ll], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[
		    ll]) + shift / d__[ll]);
#line 688 "zbdsqr.f"
	    g = e[ll];
#line 689 "zbdsqr.f"
	    i__1 = m - 1;
#line 689 "zbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 690 "zbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 691 "zbdsqr.f"
		if (i__ > ll) {
#line 691 "zbdsqr.f"
		    e[i__ - 1] = r__;
#line 691 "zbdsqr.f"
		}
#line 693 "zbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 694 "zbdsqr.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 695 "zbdsqr.f"
		g = sinr * d__[i__ + 1];
#line 696 "zbdsqr.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 697 "zbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 698 "zbdsqr.f"
		d__[i__] = r__;
#line 699 "zbdsqr.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 700 "zbdsqr.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 701 "zbdsqr.f"
		if (i__ < m - 1) {
#line 702 "zbdsqr.f"
		    g = sinl * e[i__ + 1];
#line 703 "zbdsqr.f"
		    e[i__ + 1] = cosl * e[i__ + 1];
#line 704 "zbdsqr.f"
		}
#line 705 "zbdsqr.f"
		rwork[i__ - ll + 1] = cosr;
#line 706 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sinr;
#line 707 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = cosl;
#line 708 "zbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = sinl;
#line 709 "zbdsqr.f"
/* L140: */
#line 709 "zbdsqr.f"
	    }
#line 710 "zbdsqr.f"
	    e[m - 1] = f;

/*           Update singular vectors */

#line 714 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 714 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 714 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 714 "zbdsqr.f"
	    }
#line 717 "zbdsqr.f"
	    if (*nru > 0) {
#line 717 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 717 "zbdsqr.f"
		zlasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 717 "zbdsqr.f"
	    }
#line 720 "zbdsqr.f"
	    if (*ncc > 0) {
#line 720 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 720 "zbdsqr.f"
		zlasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 720 "zbdsqr.f"
	    }

/*           Test convergence */

#line 726 "zbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 726 "zbdsqr.f"
		e[m - 1] = 0.;
#line 726 "zbdsqr.f"
	    }

#line 729 "zbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 734 "zbdsqr.f"
	    f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[m]
		    ) + shift / d__[m]);
#line 736 "zbdsqr.f"
	    g = e[m - 1];
#line 737 "zbdsqr.f"
	    i__1 = ll + 1;
#line 737 "zbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 738 "zbdsqr.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 739 "zbdsqr.f"
		if (i__ < m) {
#line 739 "zbdsqr.f"
		    e[i__] = r__;
#line 739 "zbdsqr.f"
		}
#line 741 "zbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 742 "zbdsqr.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 743 "zbdsqr.f"
		g = sinr * d__[i__ - 1];
#line 744 "zbdsqr.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 745 "zbdsqr.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 746 "zbdsqr.f"
		d__[i__] = r__;
#line 747 "zbdsqr.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 748 "zbdsqr.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 749 "zbdsqr.f"
		if (i__ > ll + 1) {
#line 750 "zbdsqr.f"
		    g = sinl * e[i__ - 2];
#line 751 "zbdsqr.f"
		    e[i__ - 2] = cosl * e[i__ - 2];
#line 752 "zbdsqr.f"
		}
#line 753 "zbdsqr.f"
		rwork[i__ - ll] = cosr;
#line 754 "zbdsqr.f"
		rwork[i__ - ll + nm1] = -sinr;
#line 755 "zbdsqr.f"
		rwork[i__ - ll + nm12] = cosl;
#line 756 "zbdsqr.f"
		rwork[i__ - ll + nm13] = -sinl;
#line 757 "zbdsqr.f"
/* L150: */
#line 757 "zbdsqr.f"
	    }
#line 758 "zbdsqr.f"
	    e[ll] = f;

/*           Test convergence */

#line 762 "zbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 762 "zbdsqr.f"
		e[ll] = 0.;
#line 762 "zbdsqr.f"
	    }

/*           Update singular vectors if desired */

#line 767 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 767 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 767 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 767 "zbdsqr.f"
	    }
#line 770 "zbdsqr.f"
	    if (*nru > 0) {
#line 770 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 770 "zbdsqr.f"
		zlasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 770 "zbdsqr.f"
	    }
#line 773 "zbdsqr.f"
	    if (*ncc > 0) {
#line 773 "zbdsqr.f"
		i__1 = m - ll + 1;
#line 773 "zbdsqr.f"
		zlasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 773 "zbdsqr.f"
	    }
#line 776 "zbdsqr.f"
	}
#line 777 "zbdsqr.f"
    }

/*     QR iteration finished, go back and check convergence */

#line 781 "zbdsqr.f"
    goto L60;

/*     All singular values converged, so make them positive */

#line 785 "zbdsqr.f"
L160:
#line 786 "zbdsqr.f"
    i__1 = *n;
#line 786 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 787 "zbdsqr.f"
	if (d__[i__] < 0.) {
#line 788 "zbdsqr.f"
	    d__[i__] = -d__[i__];

/*           Change sign of singular vectors, if desired */

#line 792 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 792 "zbdsqr.f"
		zdscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
#line 792 "zbdsqr.f"
	    }
#line 794 "zbdsqr.f"
	}
#line 795 "zbdsqr.f"
/* L170: */
#line 795 "zbdsqr.f"
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 800 "zbdsqr.f"
    i__1 = *n - 1;
#line 800 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I) */

#line 804 "zbdsqr.f"
	isub = 1;
#line 805 "zbdsqr.f"
	smin = d__[1];
#line 806 "zbdsqr.f"
	i__2 = *n + 1 - i__;
#line 806 "zbdsqr.f"
	for (j = 2; j <= i__2; ++j) {
#line 807 "zbdsqr.f"
	    if (d__[j] <= smin) {
#line 808 "zbdsqr.f"
		isub = j;
#line 809 "zbdsqr.f"
		smin = d__[j];
#line 810 "zbdsqr.f"
	    }
#line 811 "zbdsqr.f"
/* L180: */
#line 811 "zbdsqr.f"
	}
#line 812 "zbdsqr.f"
	if (isub != *n + 1 - i__) {

/*           Swap singular values and vectors */

#line 816 "zbdsqr.f"
	    d__[isub] = d__[*n + 1 - i__];
#line 817 "zbdsqr.f"
	    d__[*n + 1 - i__] = smin;
#line 818 "zbdsqr.f"
	    if (*ncvt > 0) {
#line 818 "zbdsqr.f"
		zswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
#line 818 "zbdsqr.f"
	    }
#line 821 "zbdsqr.f"
	    if (*nru > 0) {
#line 821 "zbdsqr.f"
		zswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
#line 821 "zbdsqr.f"
	    }
#line 823 "zbdsqr.f"
	    if (*ncc > 0) {
#line 823 "zbdsqr.f"
		zswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
#line 823 "zbdsqr.f"
	    }
#line 825 "zbdsqr.f"
	}
#line 826 "zbdsqr.f"
/* L190: */
#line 826 "zbdsqr.f"
    }
#line 827 "zbdsqr.f"
    goto L220;

/*     Maximum number of iterations exceeded, failure to converge */

#line 831 "zbdsqr.f"
L200:
#line 832 "zbdsqr.f"
    *info = 0;
#line 833 "zbdsqr.f"
    i__1 = *n - 1;
#line 833 "zbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 834 "zbdsqr.f"
	if (e[i__] != 0.) {
#line 834 "zbdsqr.f"
	    ++(*info);
#line 834 "zbdsqr.f"
	}
#line 836 "zbdsqr.f"
/* L210: */
#line 836 "zbdsqr.f"
    }
#line 837 "zbdsqr.f"
L220:
#line 838 "zbdsqr.f"
    return 0;

/*     End of ZBDSQR */

} /* zbdsqr_ */

