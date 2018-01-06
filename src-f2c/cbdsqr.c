#line 1 "cbdsqr.f"
/* cbdsqr.f -- translated by f2c (version 20100827).
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

#line 1 "cbdsqr.f"
/* Table of constant values */

static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;

/* > \brief \b CBDSQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CBDSQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cbdsqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cbdsqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cbdsqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/*                          LDU, C, LDC, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ), RWORK( * ) */
/*       COMPLEX            C( LDC, * ), U( LDU, * ), VT( LDVT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CBDSQR computes the singular values and, optionally, the right and/or */
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
/* > form: A = U*B*VT, as computed by CGEBRD, then */
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
/* >          D is REAL array, dimension (N) */
/* >          On entry, the n diagonal elements of the bidiagonal matrix B. */
/* >          On exit, if INFO=0, the singular values of B in decreasing */
/* >          order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
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
/* >          VT is COMPLEX array, dimension (LDVT, NCVT) */
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
/* >          U is COMPLEX array, dimension (LDU, N) */
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
/* >          C is COMPLEX array, dimension (LDC, NCC) */
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
/* >          RWORK is REAL array, dimension (4*N) */
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
/* >  TOLMUL  REAL, default = max(10,min(100,EPS**(-1/8))) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
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
    extern /* Subroutine */ int slas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal oldcs;
    extern /* Subroutine */ int clasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen);
    static integer oldll;
    static doublereal shift, sigmn, oldsn;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer maxit;
    static doublereal sminl, sigmx;
    static logical lower;
    extern /* Subroutine */ int csrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *), slasq1_(
	    integer *, doublereal *, doublereal *, doublereal *, integer *), 
	    slasv2_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen);
    static doublereal sminoa;
    extern /* Subroutine */ int slartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal thresh;
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

#line 284 "cbdsqr.f"
    /* Parameter adjustments */
#line 284 "cbdsqr.f"
    --d__;
#line 284 "cbdsqr.f"
    --e;
#line 284 "cbdsqr.f"
    vt_dim1 = *ldvt;
#line 284 "cbdsqr.f"
    vt_offset = 1 + vt_dim1;
#line 284 "cbdsqr.f"
    vt -= vt_offset;
#line 284 "cbdsqr.f"
    u_dim1 = *ldu;
#line 284 "cbdsqr.f"
    u_offset = 1 + u_dim1;
#line 284 "cbdsqr.f"
    u -= u_offset;
#line 284 "cbdsqr.f"
    c_dim1 = *ldc;
#line 284 "cbdsqr.f"
    c_offset = 1 + c_dim1;
#line 284 "cbdsqr.f"
    c__ -= c_offset;
#line 284 "cbdsqr.f"
    --rwork;
#line 284 "cbdsqr.f"

#line 284 "cbdsqr.f"
    /* Function Body */
#line 284 "cbdsqr.f"
    *info = 0;
#line 285 "cbdsqr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 286 "cbdsqr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
#line 287 "cbdsqr.f"
	*info = -1;
#line 288 "cbdsqr.f"
    } else if (*n < 0) {
#line 289 "cbdsqr.f"
	*info = -2;
#line 290 "cbdsqr.f"
    } else if (*ncvt < 0) {
#line 291 "cbdsqr.f"
	*info = -3;
#line 292 "cbdsqr.f"
    } else if (*nru < 0) {
#line 293 "cbdsqr.f"
	*info = -4;
#line 294 "cbdsqr.f"
    } else if (*ncc < 0) {
#line 295 "cbdsqr.f"
	*info = -5;
#line 296 "cbdsqr.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 298 "cbdsqr.f"
	*info = -9;
#line 299 "cbdsqr.f"
    } else if (*ldu < max(1,*nru)) {
#line 300 "cbdsqr.f"
	*info = -11;
#line 301 "cbdsqr.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 303 "cbdsqr.f"
	*info = -13;
#line 304 "cbdsqr.f"
    }
#line 305 "cbdsqr.f"
    if (*info != 0) {
#line 306 "cbdsqr.f"
	i__1 = -(*info);
#line 306 "cbdsqr.f"
	xerbla_("CBDSQR", &i__1, (ftnlen)6);
#line 307 "cbdsqr.f"
	return 0;
#line 308 "cbdsqr.f"
    }
#line 309 "cbdsqr.f"
    if (*n == 0) {
#line 309 "cbdsqr.f"
	return 0;
#line 309 "cbdsqr.f"
    }
#line 311 "cbdsqr.f"
    if (*n == 1) {
#line 311 "cbdsqr.f"
	goto L160;
#line 311 "cbdsqr.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 316 "cbdsqr.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

#line 320 "cbdsqr.f"
    if (! rotate) {
#line 321 "cbdsqr.f"
	slasq1_(n, &d__[1], &e[1], &rwork[1], info);

/*     If INFO equals 2, dqds didn't finish, try to finish */

#line 325 "cbdsqr.f"
	if (*info != 2) {
#line 325 "cbdsqr.f"
	    return 0;
#line 325 "cbdsqr.f"
	}
#line 326 "cbdsqr.f"
	*info = 0;
#line 327 "cbdsqr.f"
    }

#line 329 "cbdsqr.f"
    nm1 = *n - 1;
#line 330 "cbdsqr.f"
    nm12 = nm1 + nm1;
#line 331 "cbdsqr.f"
    nm13 = nm12 + nm1;
#line 332 "cbdsqr.f"
    idir = 0;

/*     Get machine constants */

#line 336 "cbdsqr.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 337 "cbdsqr.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 342 "cbdsqr.f"
    if (lower) {
#line 343 "cbdsqr.f"
	i__1 = *n - 1;
#line 343 "cbdsqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 344 "cbdsqr.f"
	    slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 345 "cbdsqr.f"
	    d__[i__] = r__;
#line 346 "cbdsqr.f"
	    e[i__] = sn * d__[i__ + 1];
#line 347 "cbdsqr.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 348 "cbdsqr.f"
	    rwork[i__] = cs;
#line 349 "cbdsqr.f"
	    rwork[nm1 + i__] = sn;
#line 350 "cbdsqr.f"
/* L10: */
#line 350 "cbdsqr.f"
	}

/*        Update singular vectors if desired */

#line 354 "cbdsqr.f"
	if (*nru > 0) {
#line 354 "cbdsqr.f"
	    clasr_("R", "V", "F", nru, n, &rwork[1], &rwork[*n], &u[u_offset],
		     ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 354 "cbdsqr.f"
	}
#line 357 "cbdsqr.f"
	if (*ncc > 0) {
#line 357 "cbdsqr.f"
	    clasr_("L", "V", "F", n, ncc, &rwork[1], &rwork[*n], &c__[
		    c_offset], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 357 "cbdsqr.f"
	}
#line 360 "cbdsqr.f"
    }

/*     Compute singular values to relative accuracy TOL */
/*     (By setting TOL to be negative, algorithm will compute */
/*     singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */

/* Computing MAX */
/* Computing MIN */
#line 366 "cbdsqr.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b15);
#line 366 "cbdsqr.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 366 "cbdsqr.f"
    tolmul = max(d__1,d__2);
#line 367 "cbdsqr.f"
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

#line 371 "cbdsqr.f"
    smax = 0.;
#line 372 "cbdsqr.f"
    i__1 = *n;
#line 372 "cbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 373 "cbdsqr.f"
	d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
#line 373 "cbdsqr.f"
	smax = max(d__2,d__3);
#line 374 "cbdsqr.f"
/* L20: */
#line 374 "cbdsqr.f"
    }
#line 375 "cbdsqr.f"
    i__1 = *n - 1;
#line 375 "cbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 376 "cbdsqr.f"
	d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
#line 376 "cbdsqr.f"
	smax = max(d__2,d__3);
#line 377 "cbdsqr.f"
/* L30: */
#line 377 "cbdsqr.f"
    }
#line 378 "cbdsqr.f"
    sminl = 0.;
#line 379 "cbdsqr.f"
    if (tol >= 0.) {

/*        Relative accuracy desired */

#line 383 "cbdsqr.f"
	sminoa = abs(d__[1]);
#line 384 "cbdsqr.f"
	if (sminoa == 0.) {
#line 384 "cbdsqr.f"
	    goto L50;
#line 384 "cbdsqr.f"
	}
#line 386 "cbdsqr.f"
	mu = sminoa;
#line 387 "cbdsqr.f"
	i__1 = *n;
#line 387 "cbdsqr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 388 "cbdsqr.f"
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
#line 389 "cbdsqr.f"
	    sminoa = min(sminoa,mu);
#line 390 "cbdsqr.f"
	    if (sminoa == 0.) {
#line 390 "cbdsqr.f"
		goto L50;
#line 390 "cbdsqr.f"
	    }
#line 392 "cbdsqr.f"
/* L40: */
#line 392 "cbdsqr.f"
	}
#line 393 "cbdsqr.f"
L50:
#line 394 "cbdsqr.f"
	sminoa /= sqrt((doublereal) (*n));
/* Computing MAX */
#line 395 "cbdsqr.f"
	d__1 = tol * sminoa, d__2 = *n * 6 * *n * unfl;
#line 395 "cbdsqr.f"
	thresh = max(d__1,d__2);
#line 396 "cbdsqr.f"
    } else {

/*        Absolute accuracy desired */

/* Computing MAX */
#line 400 "cbdsqr.f"
	d__1 = abs(tol) * smax, d__2 = *n * 6 * *n * unfl;
#line 400 "cbdsqr.f"
	thresh = max(d__1,d__2);
#line 401 "cbdsqr.f"
    }

/*     Prepare for main iteration loop for the singular values */
/*     (MAXIT is the maximum number of passes through the inner */
/*     loop permitted before nonconvergence signalled.) */

#line 407 "cbdsqr.f"
    maxit = *n * 6 * *n;
#line 408 "cbdsqr.f"
    iter = 0;
#line 409 "cbdsqr.f"
    oldll = -1;
#line 410 "cbdsqr.f"
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

#line 414 "cbdsqr.f"
    m = *n;

/*     Begin main iteration loop */

#line 418 "cbdsqr.f"
L60:

/*     Check for convergence or exceeding iteration count */

#line 422 "cbdsqr.f"
    if (m <= 1) {
#line 422 "cbdsqr.f"
	goto L160;
#line 422 "cbdsqr.f"
    }
#line 424 "cbdsqr.f"
    if (iter > maxit) {
#line 424 "cbdsqr.f"
	goto L200;
#line 424 "cbdsqr.f"
    }

/*     Find diagonal block of matrix to work on */

#line 429 "cbdsqr.f"
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
#line 429 "cbdsqr.f"
	d__[m] = 0.;
#line 429 "cbdsqr.f"
    }
#line 431 "cbdsqr.f"
    smax = (d__1 = d__[m], abs(d__1));
#line 432 "cbdsqr.f"
    smin = smax;
#line 433 "cbdsqr.f"
    i__1 = m - 1;
#line 433 "cbdsqr.f"
    for (lll = 1; lll <= i__1; ++lll) {
#line 434 "cbdsqr.f"
	ll = m - lll;
#line 435 "cbdsqr.f"
	abss = (d__1 = d__[ll], abs(d__1));
#line 436 "cbdsqr.f"
	abse = (d__1 = e[ll], abs(d__1));
#line 437 "cbdsqr.f"
	if (tol < 0. && abss <= thresh) {
#line 437 "cbdsqr.f"
	    d__[ll] = 0.;
#line 437 "cbdsqr.f"
	}
#line 439 "cbdsqr.f"
	if (abse <= thresh) {
#line 439 "cbdsqr.f"
	    goto L80;
#line 439 "cbdsqr.f"
	}
#line 441 "cbdsqr.f"
	smin = min(smin,abss);
/* Computing MAX */
#line 442 "cbdsqr.f"
	d__1 = max(smax,abss);
#line 442 "cbdsqr.f"
	smax = max(d__1,abse);
#line 443 "cbdsqr.f"
/* L70: */
#line 443 "cbdsqr.f"
    }
#line 444 "cbdsqr.f"
    ll = 0;
#line 445 "cbdsqr.f"
    goto L90;
#line 446 "cbdsqr.f"
L80:
#line 447 "cbdsqr.f"
    e[ll] = 0.;

/*     Matrix splits since E(LL) = 0 */

#line 451 "cbdsqr.f"
    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop */

#line 455 "cbdsqr.f"
	--m;
#line 456 "cbdsqr.f"
	goto L60;
#line 457 "cbdsqr.f"
    }
#line 458 "cbdsqr.f"
L90:
#line 459 "cbdsqr.f"
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

#line 463 "cbdsqr.f"
    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

#line 467 "cbdsqr.f"
	slasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
#line 469 "cbdsqr.f"
	d__[m - 1] = sigmx;
#line 470 "cbdsqr.f"
	e[m - 1] = 0.;
#line 471 "cbdsqr.f"
	d__[m] = sigmn;

/*        Compute singular vectors, if desired */

#line 475 "cbdsqr.f"
	if (*ncvt > 0) {
#line 475 "cbdsqr.f"
	    csrot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
#line 475 "cbdsqr.f"
	}
#line 478 "cbdsqr.f"
	if (*nru > 0) {
#line 478 "cbdsqr.f"
	    csrot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
#line 478 "cbdsqr.f"
	}
#line 480 "cbdsqr.f"
	if (*ncc > 0) {
#line 480 "cbdsqr.f"
	    csrot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
#line 480 "cbdsqr.f"
	}
#line 483 "cbdsqr.f"
	m += -2;
#line 484 "cbdsqr.f"
	goto L60;
#line 485 "cbdsqr.f"
    }

/*     If working on new submatrix, choose shift direction */
/*     (from larger end diagonal element towards smaller) */

#line 490 "cbdsqr.f"
    if (ll > oldm || m < oldll) {
#line 491 "cbdsqr.f"
	if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {

/*           Chase bulge from top (big end) to bottom (small end) */

#line 495 "cbdsqr.f"
	    idir = 1;
#line 496 "cbdsqr.f"
	} else {

/*           Chase bulge from bottom (big end) to top (small end) */

#line 500 "cbdsqr.f"
	    idir = 2;
#line 501 "cbdsqr.f"
	}
#line 502 "cbdsqr.f"
    }

/*     Apply convergence tests */

#line 506 "cbdsqr.f"
    if (idir == 1) {

/*        Run convergence test in forward direction */
/*        First apply standard test to bottom of matrix */

#line 511 "cbdsqr.f"
	if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(
		d__1)) || tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) 
		{
#line 513 "cbdsqr.f"
	    e[m - 1] = 0.;
#line 514 "cbdsqr.f"
	    goto L60;
#line 515 "cbdsqr.f"
	}

#line 517 "cbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion forward */

#line 522 "cbdsqr.f"
	    mu = (d__1 = d__[ll], abs(d__1));
#line 523 "cbdsqr.f"
	    sminl = mu;
#line 524 "cbdsqr.f"
	    i__1 = m - 1;
#line 524 "cbdsqr.f"
	    for (lll = ll; lll <= i__1; ++lll) {
#line 525 "cbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 526 "cbdsqr.f"
		    e[lll] = 0.;
#line 527 "cbdsqr.f"
		    goto L60;
#line 528 "cbdsqr.f"
		}
#line 529 "cbdsqr.f"
		mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[
			lll], abs(d__1))));
#line 530 "cbdsqr.f"
		sminl = min(sminl,mu);
#line 531 "cbdsqr.f"
/* L100: */
#line 531 "cbdsqr.f"
	    }
#line 532 "cbdsqr.f"
	}

#line 534 "cbdsqr.f"
    } else {

/*        Run convergence test in backward direction */
/*        First apply standard test to top of matrix */

#line 539 "cbdsqr.f"
	if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)
		) || tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
#line 541 "cbdsqr.f"
	    e[ll] = 0.;
#line 542 "cbdsqr.f"
	    goto L60;
#line 543 "cbdsqr.f"
	}

#line 545 "cbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion backward */

#line 550 "cbdsqr.f"
	    mu = (d__1 = d__[m], abs(d__1));
#line 551 "cbdsqr.f"
	    sminl = mu;
#line 552 "cbdsqr.f"
	    i__1 = ll;
#line 552 "cbdsqr.f"
	    for (lll = m - 1; lll >= i__1; --lll) {
#line 553 "cbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 554 "cbdsqr.f"
		    e[lll] = 0.;
#line 555 "cbdsqr.f"
		    goto L60;
#line 556 "cbdsqr.f"
		}
#line 557 "cbdsqr.f"
		mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll]
			, abs(d__1))));
#line 558 "cbdsqr.f"
		sminl = min(sminl,mu);
#line 559 "cbdsqr.f"
/* L110: */
#line 559 "cbdsqr.f"
	    }
#line 560 "cbdsqr.f"
	}
#line 561 "cbdsqr.f"
    }
#line 562 "cbdsqr.f"
    oldll = ll;
#line 563 "cbdsqr.f"
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative */
/*     accuracy, and if so set the shift to zero. */

/* Computing MAX */
#line 568 "cbdsqr.f"
    d__1 = eps, d__2 = tol * .01;
#line 568 "cbdsqr.f"
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

#line 573 "cbdsqr.f"
	shift = 0.;
#line 574 "cbdsqr.f"
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

#line 578 "cbdsqr.f"
	if (idir == 1) {
#line 579 "cbdsqr.f"
	    sll = (d__1 = d__[ll], abs(d__1));
#line 580 "cbdsqr.f"
	    slas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
#line 581 "cbdsqr.f"
	} else {
#line 582 "cbdsqr.f"
	    sll = (d__1 = d__[m], abs(d__1));
#line 583 "cbdsqr.f"
	    slas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
#line 584 "cbdsqr.f"
	}

/*        Test if shift negligible, and if so set to zero */

#line 588 "cbdsqr.f"
	if (sll > 0.) {
/* Computing 2nd power */
#line 589 "cbdsqr.f"
	    d__1 = shift / sll;
#line 589 "cbdsqr.f"
	    if (d__1 * d__1 < eps) {
#line 589 "cbdsqr.f"
		shift = 0.;
#line 589 "cbdsqr.f"
	    }
#line 591 "cbdsqr.f"
	}
#line 592 "cbdsqr.f"
    }

/*     Increment iteration count */

#line 596 "cbdsqr.f"
    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

#line 600 "cbdsqr.f"
    if (shift == 0.) {
#line 601 "cbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 606 "cbdsqr.f"
	    cs = 1.;
#line 607 "cbdsqr.f"
	    oldcs = 1.;
#line 608 "cbdsqr.f"
	    i__1 = m - 1;
#line 608 "cbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 609 "cbdsqr.f"
		d__1 = d__[i__] * cs;
#line 609 "cbdsqr.f"
		slartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 610 "cbdsqr.f"
		if (i__ > ll) {
#line 610 "cbdsqr.f"
		    e[i__ - 1] = oldsn * r__;
#line 610 "cbdsqr.f"
		}
#line 612 "cbdsqr.f"
		d__1 = oldcs * r__;
#line 612 "cbdsqr.f"
		d__2 = d__[i__ + 1] * sn;
#line 612 "cbdsqr.f"
		slartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 613 "cbdsqr.f"
		rwork[i__ - ll + 1] = cs;
#line 614 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sn;
#line 615 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = oldcs;
#line 616 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = oldsn;
#line 617 "cbdsqr.f"
/* L120: */
#line 617 "cbdsqr.f"
	    }
#line 618 "cbdsqr.f"
	    h__ = d__[m] * cs;
#line 619 "cbdsqr.f"
	    d__[m] = h__ * oldcs;
#line 620 "cbdsqr.f"
	    e[m - 1] = h__ * oldsn;

/*           Update singular vectors */

#line 624 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 624 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 624 "cbdsqr.f"
		clasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 624 "cbdsqr.f"
	    }
#line 627 "cbdsqr.f"
	    if (*nru > 0) {
#line 627 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 627 "cbdsqr.f"
		clasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 627 "cbdsqr.f"
	    }
#line 630 "cbdsqr.f"
	    if (*ncc > 0) {
#line 630 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 630 "cbdsqr.f"
		clasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 630 "cbdsqr.f"
	    }

/*           Test convergence */

#line 636 "cbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 636 "cbdsqr.f"
		e[m - 1] = 0.;
#line 636 "cbdsqr.f"
	    }

#line 639 "cbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 644 "cbdsqr.f"
	    cs = 1.;
#line 645 "cbdsqr.f"
	    oldcs = 1.;
#line 646 "cbdsqr.f"
	    i__1 = ll + 1;
#line 646 "cbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 647 "cbdsqr.f"
		d__1 = d__[i__] * cs;
#line 647 "cbdsqr.f"
		slartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 648 "cbdsqr.f"
		if (i__ < m) {
#line 648 "cbdsqr.f"
		    e[i__] = oldsn * r__;
#line 648 "cbdsqr.f"
		}
#line 650 "cbdsqr.f"
		d__1 = oldcs * r__;
#line 650 "cbdsqr.f"
		d__2 = d__[i__ - 1] * sn;
#line 650 "cbdsqr.f"
		slartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 651 "cbdsqr.f"
		rwork[i__ - ll] = cs;
#line 652 "cbdsqr.f"
		rwork[i__ - ll + nm1] = -sn;
#line 653 "cbdsqr.f"
		rwork[i__ - ll + nm12] = oldcs;
#line 654 "cbdsqr.f"
		rwork[i__ - ll + nm13] = -oldsn;
#line 655 "cbdsqr.f"
/* L130: */
#line 655 "cbdsqr.f"
	    }
#line 656 "cbdsqr.f"
	    h__ = d__[ll] * cs;
#line 657 "cbdsqr.f"
	    d__[ll] = h__ * oldcs;
#line 658 "cbdsqr.f"
	    e[ll] = h__ * oldsn;

/*           Update singular vectors */

#line 662 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 662 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 662 "cbdsqr.f"
		clasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 662 "cbdsqr.f"
	    }
#line 665 "cbdsqr.f"
	    if (*nru > 0) {
#line 665 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 665 "cbdsqr.f"
		clasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 665 "cbdsqr.f"
	    }
#line 668 "cbdsqr.f"
	    if (*ncc > 0) {
#line 668 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 668 "cbdsqr.f"
		clasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 668 "cbdsqr.f"
	    }

/*           Test convergence */

#line 674 "cbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 674 "cbdsqr.f"
		e[ll] = 0.;
#line 674 "cbdsqr.f"
	    }
#line 676 "cbdsqr.f"
	}
#line 677 "cbdsqr.f"
    } else {

/*        Use nonzero shift */

#line 681 "cbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 686 "cbdsqr.f"
	    f = ((d__1 = d__[ll], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[
		    ll]) + shift / d__[ll]);
#line 688 "cbdsqr.f"
	    g = e[ll];
#line 689 "cbdsqr.f"
	    i__1 = m - 1;
#line 689 "cbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 690 "cbdsqr.f"
		slartg_(&f, &g, &cosr, &sinr, &r__);
#line 691 "cbdsqr.f"
		if (i__ > ll) {
#line 691 "cbdsqr.f"
		    e[i__ - 1] = r__;
#line 691 "cbdsqr.f"
		}
#line 693 "cbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 694 "cbdsqr.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 695 "cbdsqr.f"
		g = sinr * d__[i__ + 1];
#line 696 "cbdsqr.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 697 "cbdsqr.f"
		slartg_(&f, &g, &cosl, &sinl, &r__);
#line 698 "cbdsqr.f"
		d__[i__] = r__;
#line 699 "cbdsqr.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 700 "cbdsqr.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 701 "cbdsqr.f"
		if (i__ < m - 1) {
#line 702 "cbdsqr.f"
		    g = sinl * e[i__ + 1];
#line 703 "cbdsqr.f"
		    e[i__ + 1] = cosl * e[i__ + 1];
#line 704 "cbdsqr.f"
		}
#line 705 "cbdsqr.f"
		rwork[i__ - ll + 1] = cosr;
#line 706 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm1] = sinr;
#line 707 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm12] = cosl;
#line 708 "cbdsqr.f"
		rwork[i__ - ll + 1 + nm13] = sinl;
#line 709 "cbdsqr.f"
/* L140: */
#line 709 "cbdsqr.f"
	    }
#line 710 "cbdsqr.f"
	    e[m - 1] = f;

/*           Update singular vectors */

#line 714 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 714 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 714 "cbdsqr.f"
		clasr_("L", "V", "F", &i__1, ncvt, &rwork[1], &rwork[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 714 "cbdsqr.f"
	    }
#line 717 "cbdsqr.f"
	    if (*nru > 0) {
#line 717 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 717 "cbdsqr.f"
		clasr_("R", "V", "F", nru, &i__1, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 717 "cbdsqr.f"
	    }
#line 720 "cbdsqr.f"
	    if (*ncc > 0) {
#line 720 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 720 "cbdsqr.f"
		clasr_("L", "V", "F", &i__1, ncc, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)
			1, (ftnlen)1);
#line 720 "cbdsqr.f"
	    }

/*           Test convergence */

#line 726 "cbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 726 "cbdsqr.f"
		e[m - 1] = 0.;
#line 726 "cbdsqr.f"
	    }

#line 729 "cbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 734 "cbdsqr.f"
	    f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[m]
		    ) + shift / d__[m]);
#line 736 "cbdsqr.f"
	    g = e[m - 1];
#line 737 "cbdsqr.f"
	    i__1 = ll + 1;
#line 737 "cbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 738 "cbdsqr.f"
		slartg_(&f, &g, &cosr, &sinr, &r__);
#line 739 "cbdsqr.f"
		if (i__ < m) {
#line 739 "cbdsqr.f"
		    e[i__] = r__;
#line 739 "cbdsqr.f"
		}
#line 741 "cbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 742 "cbdsqr.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 743 "cbdsqr.f"
		g = sinr * d__[i__ - 1];
#line 744 "cbdsqr.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 745 "cbdsqr.f"
		slartg_(&f, &g, &cosl, &sinl, &r__);
#line 746 "cbdsqr.f"
		d__[i__] = r__;
#line 747 "cbdsqr.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 748 "cbdsqr.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 749 "cbdsqr.f"
		if (i__ > ll + 1) {
#line 750 "cbdsqr.f"
		    g = sinl * e[i__ - 2];
#line 751 "cbdsqr.f"
		    e[i__ - 2] = cosl * e[i__ - 2];
#line 752 "cbdsqr.f"
		}
#line 753 "cbdsqr.f"
		rwork[i__ - ll] = cosr;
#line 754 "cbdsqr.f"
		rwork[i__ - ll + nm1] = -sinr;
#line 755 "cbdsqr.f"
		rwork[i__ - ll + nm12] = cosl;
#line 756 "cbdsqr.f"
		rwork[i__ - ll + nm13] = -sinl;
#line 757 "cbdsqr.f"
/* L150: */
#line 757 "cbdsqr.f"
	    }
#line 758 "cbdsqr.f"
	    e[ll] = f;

/*           Test convergence */

#line 762 "cbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 762 "cbdsqr.f"
		e[ll] = 0.;
#line 762 "cbdsqr.f"
	    }

/*           Update singular vectors if desired */

#line 767 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 767 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 767 "cbdsqr.f"
		clasr_("L", "V", "B", &i__1, ncvt, &rwork[nm12 + 1], &rwork[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 767 "cbdsqr.f"
	    }
#line 770 "cbdsqr.f"
	    if (*nru > 0) {
#line 770 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 770 "cbdsqr.f"
		clasr_("R", "V", "B", nru, &i__1, &rwork[1], &rwork[*n], &u[
			ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 770 "cbdsqr.f"
	    }
#line 773 "cbdsqr.f"
	    if (*ncc > 0) {
#line 773 "cbdsqr.f"
		i__1 = m - ll + 1;
#line 773 "cbdsqr.f"
		clasr_("L", "V", "B", &i__1, ncc, &rwork[1], &rwork[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 773 "cbdsqr.f"
	    }
#line 776 "cbdsqr.f"
	}
#line 777 "cbdsqr.f"
    }

/*     QR iteration finished, go back and check convergence */

#line 781 "cbdsqr.f"
    goto L60;

/*     All singular values converged, so make them positive */

#line 785 "cbdsqr.f"
L160:
#line 786 "cbdsqr.f"
    i__1 = *n;
#line 786 "cbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 787 "cbdsqr.f"
	if (d__[i__] < 0.) {
#line 788 "cbdsqr.f"
	    d__[i__] = -d__[i__];

/*           Change sign of singular vectors, if desired */

#line 792 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 792 "cbdsqr.f"
		csscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
#line 792 "cbdsqr.f"
	    }
#line 794 "cbdsqr.f"
	}
#line 795 "cbdsqr.f"
/* L170: */
#line 795 "cbdsqr.f"
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 800 "cbdsqr.f"
    i__1 = *n - 1;
#line 800 "cbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I) */

#line 804 "cbdsqr.f"
	isub = 1;
#line 805 "cbdsqr.f"
	smin = d__[1];
#line 806 "cbdsqr.f"
	i__2 = *n + 1 - i__;
#line 806 "cbdsqr.f"
	for (j = 2; j <= i__2; ++j) {
#line 807 "cbdsqr.f"
	    if (d__[j] <= smin) {
#line 808 "cbdsqr.f"
		isub = j;
#line 809 "cbdsqr.f"
		smin = d__[j];
#line 810 "cbdsqr.f"
	    }
#line 811 "cbdsqr.f"
/* L180: */
#line 811 "cbdsqr.f"
	}
#line 812 "cbdsqr.f"
	if (isub != *n + 1 - i__) {

/*           Swap singular values and vectors */

#line 816 "cbdsqr.f"
	    d__[isub] = d__[*n + 1 - i__];
#line 817 "cbdsqr.f"
	    d__[*n + 1 - i__] = smin;
#line 818 "cbdsqr.f"
	    if (*ncvt > 0) {
#line 818 "cbdsqr.f"
		cswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
#line 818 "cbdsqr.f"
	    }
#line 821 "cbdsqr.f"
	    if (*nru > 0) {
#line 821 "cbdsqr.f"
		cswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
#line 821 "cbdsqr.f"
	    }
#line 823 "cbdsqr.f"
	    if (*ncc > 0) {
#line 823 "cbdsqr.f"
		cswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
#line 823 "cbdsqr.f"
	    }
#line 825 "cbdsqr.f"
	}
#line 826 "cbdsqr.f"
/* L190: */
#line 826 "cbdsqr.f"
    }
#line 827 "cbdsqr.f"
    goto L220;

/*     Maximum number of iterations exceeded, failure to converge */

#line 831 "cbdsqr.f"
L200:
#line 832 "cbdsqr.f"
    *info = 0;
#line 833 "cbdsqr.f"
    i__1 = *n - 1;
#line 833 "cbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 834 "cbdsqr.f"
	if (e[i__] != 0.) {
#line 834 "cbdsqr.f"
	    ++(*info);
#line 834 "cbdsqr.f"
	}
#line 836 "cbdsqr.f"
/* L210: */
#line 836 "cbdsqr.f"
    }
#line 837 "cbdsqr.f"
L220:
#line 838 "cbdsqr.f"
    return 0;

/*     End of CBDSQR */

} /* cbdsqr_ */

