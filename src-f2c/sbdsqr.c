#line 1 "sbdsqr.f"
/* sbdsqr.f -- translated by f2c (version 20100827).
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

#line 1 "sbdsqr.f"
/* Table of constant values */

static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;

/* > \brief \b SBDSQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SBDSQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sbdsqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sbdsqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sbdsqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, */
/*                          LDU, C, LDC, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( LDC, * ), D( * ), E( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SBDSQR computes the singular values and, optionally, the right and/or */
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
/* > form:  A = U*B*VT, as computed by SGEBRD, then */
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
/* >          VT is REAL array, dimension (LDVT, NCVT) */
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
/* >          U is REAL array, dimension (LDU, N) */
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
/* >          C is REAL array, dimension (LDC, NCC) */
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
/* >          WORK is REAL array, dimension (4*N) */
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

/* > \par Note: */
/*  =========== */
/* > */
/* > \verbatim */
/* >  Bug report from Cezary Dendek. */
/* >  On March 23rd 2017, the INTEGER variable MAXIT = MAXITR*N**2 is */
/* >  removed since it can overflow pretty easily (for N larger or equal */
/* >  than 18,919). We instead use MAXITDIVN = MAXITR*N. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
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
    static integer iterdivn;
    static doublereal f, g, h__;
    static integer i__, j, m;
    static doublereal r__;
    static integer maxitdivn;
    static doublereal cs;
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
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), slas2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal oldcs;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer oldll;
    static doublereal shift, sigmn, oldsn, sminl;
    extern /* Subroutine */ int slasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal sigmx;
    static logical lower;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), slasq1_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *), slasv2_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal sminoa;
    extern /* Subroutine */ int slartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal thresh;
    static logical rotate;
    static doublereal tolmul;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 302 "sbdsqr.f"
    /* Parameter adjustments */
#line 302 "sbdsqr.f"
    --d__;
#line 302 "sbdsqr.f"
    --e;
#line 302 "sbdsqr.f"
    vt_dim1 = *ldvt;
#line 302 "sbdsqr.f"
    vt_offset = 1 + vt_dim1;
#line 302 "sbdsqr.f"
    vt -= vt_offset;
#line 302 "sbdsqr.f"
    u_dim1 = *ldu;
#line 302 "sbdsqr.f"
    u_offset = 1 + u_dim1;
#line 302 "sbdsqr.f"
    u -= u_offset;
#line 302 "sbdsqr.f"
    c_dim1 = *ldc;
#line 302 "sbdsqr.f"
    c_offset = 1 + c_dim1;
#line 302 "sbdsqr.f"
    c__ -= c_offset;
#line 302 "sbdsqr.f"
    --work;
#line 302 "sbdsqr.f"

#line 302 "sbdsqr.f"
    /* Function Body */
#line 302 "sbdsqr.f"
    *info = 0;
#line 303 "sbdsqr.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 304 "sbdsqr.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lower) {
#line 305 "sbdsqr.f"
	*info = -1;
#line 306 "sbdsqr.f"
    } else if (*n < 0) {
#line 307 "sbdsqr.f"
	*info = -2;
#line 308 "sbdsqr.f"
    } else if (*ncvt < 0) {
#line 309 "sbdsqr.f"
	*info = -3;
#line 310 "sbdsqr.f"
    } else if (*nru < 0) {
#line 311 "sbdsqr.f"
	*info = -4;
#line 312 "sbdsqr.f"
    } else if (*ncc < 0) {
#line 313 "sbdsqr.f"
	*info = -5;
#line 314 "sbdsqr.f"
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
#line 316 "sbdsqr.f"
	*info = -9;
#line 317 "sbdsqr.f"
    } else if (*ldu < max(1,*nru)) {
#line 318 "sbdsqr.f"
	*info = -11;
#line 319 "sbdsqr.f"
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
#line 321 "sbdsqr.f"
	*info = -13;
#line 322 "sbdsqr.f"
    }
#line 323 "sbdsqr.f"
    if (*info != 0) {
#line 324 "sbdsqr.f"
	i__1 = -(*info);
#line 324 "sbdsqr.f"
	xerbla_("SBDSQR", &i__1, (ftnlen)6);
#line 325 "sbdsqr.f"
	return 0;
#line 326 "sbdsqr.f"
    }
#line 327 "sbdsqr.f"
    if (*n == 0) {
#line 327 "sbdsqr.f"
	return 0;
#line 327 "sbdsqr.f"
    }
#line 329 "sbdsqr.f"
    if (*n == 1) {
#line 329 "sbdsqr.f"
	goto L160;
#line 329 "sbdsqr.f"
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

#line 334 "sbdsqr.f"
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

#line 338 "sbdsqr.f"
    if (! rotate) {
#line 339 "sbdsqr.f"
	slasq1_(n, &d__[1], &e[1], &work[1], info);

/*     If INFO equals 2, dqds didn't finish, try to finish */

#line 343 "sbdsqr.f"
	if (*info != 2) {
#line 343 "sbdsqr.f"
	    return 0;
#line 343 "sbdsqr.f"
	}
#line 344 "sbdsqr.f"
	*info = 0;
#line 345 "sbdsqr.f"
    }

#line 347 "sbdsqr.f"
    nm1 = *n - 1;
#line 348 "sbdsqr.f"
    nm12 = nm1 + nm1;
#line 349 "sbdsqr.f"
    nm13 = nm12 + nm1;
#line 350 "sbdsqr.f"
    idir = 0;

/*     Get machine constants */

#line 354 "sbdsqr.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 355 "sbdsqr.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */

#line 360 "sbdsqr.f"
    if (lower) {
#line 361 "sbdsqr.f"
	i__1 = *n - 1;
#line 361 "sbdsqr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 362 "sbdsqr.f"
	    slartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
#line 363 "sbdsqr.f"
	    d__[i__] = r__;
#line 364 "sbdsqr.f"
	    e[i__] = sn * d__[i__ + 1];
#line 365 "sbdsqr.f"
	    d__[i__ + 1] = cs * d__[i__ + 1];
#line 366 "sbdsqr.f"
	    work[i__] = cs;
#line 367 "sbdsqr.f"
	    work[nm1 + i__] = sn;
#line 368 "sbdsqr.f"
/* L10: */
#line 368 "sbdsqr.f"
	}

/*        Update singular vectors if desired */

#line 372 "sbdsqr.f"
	if (*nru > 0) {
#line 372 "sbdsqr.f"
	    slasr_("R", "V", "F", nru, n, &work[1], &work[*n], &u[u_offset], 
		    ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 372 "sbdsqr.f"
	}
#line 375 "sbdsqr.f"
	if (*ncc > 0) {
#line 375 "sbdsqr.f"
	    slasr_("L", "V", "F", n, ncc, &work[1], &work[*n], &c__[c_offset],
		     ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 375 "sbdsqr.f"
	}
#line 378 "sbdsqr.f"
    }

/*     Compute singular values to relative accuracy TOL */
/*     (By setting TOL to be negative, algorithm will compute */
/*     singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */

/* Computing MAX */
/* Computing MIN */
#line 384 "sbdsqr.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b15);
#line 384 "sbdsqr.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 384 "sbdsqr.f"
    tolmul = max(d__1,d__2);
#line 385 "sbdsqr.f"
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

#line 389 "sbdsqr.f"
    smax = 0.;
#line 390 "sbdsqr.f"
    i__1 = *n;
#line 390 "sbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 391 "sbdsqr.f"
	d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
#line 391 "sbdsqr.f"
	smax = max(d__2,d__3);
#line 392 "sbdsqr.f"
/* L20: */
#line 392 "sbdsqr.f"
    }
#line 393 "sbdsqr.f"
    i__1 = *n - 1;
#line 393 "sbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 394 "sbdsqr.f"
	d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
#line 394 "sbdsqr.f"
	smax = max(d__2,d__3);
#line 395 "sbdsqr.f"
/* L30: */
#line 395 "sbdsqr.f"
    }
#line 396 "sbdsqr.f"
    sminl = 0.;
#line 397 "sbdsqr.f"
    if (tol >= 0.) {

/*        Relative accuracy desired */

#line 401 "sbdsqr.f"
	sminoa = abs(d__[1]);
#line 402 "sbdsqr.f"
	if (sminoa == 0.) {
#line 402 "sbdsqr.f"
	    goto L50;
#line 402 "sbdsqr.f"
	}
#line 404 "sbdsqr.f"
	mu = sminoa;
#line 405 "sbdsqr.f"
	i__1 = *n;
#line 405 "sbdsqr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 406 "sbdsqr.f"
	    mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1]
		    , abs(d__1))));
#line 407 "sbdsqr.f"
	    sminoa = min(sminoa,mu);
#line 408 "sbdsqr.f"
	    if (sminoa == 0.) {
#line 408 "sbdsqr.f"
		goto L50;
#line 408 "sbdsqr.f"
	    }
#line 410 "sbdsqr.f"
/* L40: */
#line 410 "sbdsqr.f"
	}
#line 411 "sbdsqr.f"
L50:
#line 412 "sbdsqr.f"
	sminoa /= sqrt((doublereal) (*n));
/* Computing MAX */
#line 413 "sbdsqr.f"
	d__1 = tol * sminoa, d__2 = *n * (*n * unfl) * 6;
#line 413 "sbdsqr.f"
	thresh = max(d__1,d__2);
#line 414 "sbdsqr.f"
    } else {

/*        Absolute accuracy desired */

/* Computing MAX */
#line 418 "sbdsqr.f"
	d__1 = abs(tol) * smax, d__2 = *n * (*n * unfl) * 6;
#line 418 "sbdsqr.f"
	thresh = max(d__1,d__2);
#line 419 "sbdsqr.f"
    }

/*     Prepare for main iteration loop for the singular values */
/*     (MAXIT is the maximum number of passes through the inner */
/*     loop permitted before nonconvergence signalled.) */

#line 425 "sbdsqr.f"
    maxitdivn = *n * 6;
#line 426 "sbdsqr.f"
    iterdivn = 0;
#line 427 "sbdsqr.f"
    iter = -1;
#line 428 "sbdsqr.f"
    oldll = -1;
#line 429 "sbdsqr.f"
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

#line 433 "sbdsqr.f"
    m = *n;

/*     Begin main iteration loop */

#line 437 "sbdsqr.f"
L60:

/*     Check for convergence or exceeding iteration count */

#line 441 "sbdsqr.f"
    if (m <= 1) {
#line 441 "sbdsqr.f"
	goto L160;
#line 441 "sbdsqr.f"
    }

#line 444 "sbdsqr.f"
    if (iter >= *n) {
#line 445 "sbdsqr.f"
	iter -= *n;
#line 446 "sbdsqr.f"
	++iterdivn;
#line 447 "sbdsqr.f"
	if (iterdivn >= maxitdivn) {
#line 447 "sbdsqr.f"
	    goto L200;
#line 447 "sbdsqr.f"
	}
#line 449 "sbdsqr.f"
    }

/*     Find diagonal block of matrix to work on */

#line 453 "sbdsqr.f"
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
#line 453 "sbdsqr.f"
	d__[m] = 0.;
#line 453 "sbdsqr.f"
    }
#line 455 "sbdsqr.f"
    smax = (d__1 = d__[m], abs(d__1));
#line 456 "sbdsqr.f"
    smin = smax;
#line 457 "sbdsqr.f"
    i__1 = m - 1;
#line 457 "sbdsqr.f"
    for (lll = 1; lll <= i__1; ++lll) {
#line 458 "sbdsqr.f"
	ll = m - lll;
#line 459 "sbdsqr.f"
	abss = (d__1 = d__[ll], abs(d__1));
#line 460 "sbdsqr.f"
	abse = (d__1 = e[ll], abs(d__1));
#line 461 "sbdsqr.f"
	if (tol < 0. && abss <= thresh) {
#line 461 "sbdsqr.f"
	    d__[ll] = 0.;
#line 461 "sbdsqr.f"
	}
#line 463 "sbdsqr.f"
	if (abse <= thresh) {
#line 463 "sbdsqr.f"
	    goto L80;
#line 463 "sbdsqr.f"
	}
#line 465 "sbdsqr.f"
	smin = min(smin,abss);
/* Computing MAX */
#line 466 "sbdsqr.f"
	d__1 = max(smax,abss);
#line 466 "sbdsqr.f"
	smax = max(d__1,abse);
#line 467 "sbdsqr.f"
/* L70: */
#line 467 "sbdsqr.f"
    }
#line 468 "sbdsqr.f"
    ll = 0;
#line 469 "sbdsqr.f"
    goto L90;
#line 470 "sbdsqr.f"
L80:
#line 471 "sbdsqr.f"
    e[ll] = 0.;

/*     Matrix splits since E(LL) = 0 */

#line 475 "sbdsqr.f"
    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop */

#line 479 "sbdsqr.f"
	--m;
#line 480 "sbdsqr.f"
	goto L60;
#line 481 "sbdsqr.f"
    }
#line 482 "sbdsqr.f"
L90:
#line 483 "sbdsqr.f"
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

#line 487 "sbdsqr.f"
    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

#line 491 "sbdsqr.f"
	slasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
#line 493 "sbdsqr.f"
	d__[m - 1] = sigmx;
#line 494 "sbdsqr.f"
	e[m - 1] = 0.;
#line 495 "sbdsqr.f"
	d__[m] = sigmn;

/*        Compute singular vectors, if desired */

#line 499 "sbdsqr.f"
	if (*ncvt > 0) {
#line 499 "sbdsqr.f"
	    srot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
#line 499 "sbdsqr.f"
	}
#line 502 "sbdsqr.f"
	if (*nru > 0) {
#line 502 "sbdsqr.f"
	    srot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
#line 502 "sbdsqr.f"
	}
#line 504 "sbdsqr.f"
	if (*ncc > 0) {
#line 504 "sbdsqr.f"
	    srot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
#line 504 "sbdsqr.f"
	}
#line 507 "sbdsqr.f"
	m += -2;
#line 508 "sbdsqr.f"
	goto L60;
#line 509 "sbdsqr.f"
    }

/*     If working on new submatrix, choose shift direction */
/*     (from larger end diagonal element towards smaller) */

#line 514 "sbdsqr.f"
    if (ll > oldm || m < oldll) {
#line 515 "sbdsqr.f"
	if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {

/*           Chase bulge from top (big end) to bottom (small end) */

#line 519 "sbdsqr.f"
	    idir = 1;
#line 520 "sbdsqr.f"
	} else {

/*           Chase bulge from bottom (big end) to top (small end) */

#line 524 "sbdsqr.f"
	    idir = 2;
#line 525 "sbdsqr.f"
	}
#line 526 "sbdsqr.f"
    }

/*     Apply convergence tests */

#line 530 "sbdsqr.f"
    if (idir == 1) {

/*        Run convergence test in forward direction */
/*        First apply standard test to bottom of matrix */

#line 535 "sbdsqr.f"
	if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(
		d__1)) || tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) 
		{
#line 537 "sbdsqr.f"
	    e[m - 1] = 0.;
#line 538 "sbdsqr.f"
	    goto L60;
#line 539 "sbdsqr.f"
	}

#line 541 "sbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion forward */

#line 546 "sbdsqr.f"
	    mu = (d__1 = d__[ll], abs(d__1));
#line 547 "sbdsqr.f"
	    sminl = mu;
#line 548 "sbdsqr.f"
	    i__1 = m - 1;
#line 548 "sbdsqr.f"
	    for (lll = ll; lll <= i__1; ++lll) {
#line 549 "sbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 550 "sbdsqr.f"
		    e[lll] = 0.;
#line 551 "sbdsqr.f"
		    goto L60;
#line 552 "sbdsqr.f"
		}
#line 553 "sbdsqr.f"
		mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[
			lll], abs(d__1))));
#line 554 "sbdsqr.f"
		sminl = min(sminl,mu);
#line 555 "sbdsqr.f"
/* L100: */
#line 555 "sbdsqr.f"
	    }
#line 556 "sbdsqr.f"
	}

#line 558 "sbdsqr.f"
    } else {

/*        Run convergence test in backward direction */
/*        First apply standard test to top of matrix */

#line 563 "sbdsqr.f"
	if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)
		) || tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
#line 565 "sbdsqr.f"
	    e[ll] = 0.;
#line 566 "sbdsqr.f"
	    goto L60;
#line 567 "sbdsqr.f"
	}

#line 569 "sbdsqr.f"
	if (tol >= 0.) {

/*           If relative accuracy desired, */
/*           apply convergence criterion backward */

#line 574 "sbdsqr.f"
	    mu = (d__1 = d__[m], abs(d__1));
#line 575 "sbdsqr.f"
	    sminl = mu;
#line 576 "sbdsqr.f"
	    i__1 = ll;
#line 576 "sbdsqr.f"
	    for (lll = m - 1; lll >= i__1; --lll) {
#line 577 "sbdsqr.f"
		if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
#line 578 "sbdsqr.f"
		    e[lll] = 0.;
#line 579 "sbdsqr.f"
		    goto L60;
#line 580 "sbdsqr.f"
		}
#line 581 "sbdsqr.f"
		mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll]
			, abs(d__1))));
#line 582 "sbdsqr.f"
		sminl = min(sminl,mu);
#line 583 "sbdsqr.f"
/* L110: */
#line 583 "sbdsqr.f"
	    }
#line 584 "sbdsqr.f"
	}
#line 585 "sbdsqr.f"
    }
#line 586 "sbdsqr.f"
    oldll = ll;
#line 587 "sbdsqr.f"
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative */
/*     accuracy, and if so set the shift to zero. */

/* Computing MAX */
#line 592 "sbdsqr.f"
    d__1 = eps, d__2 = tol * .01;
#line 592 "sbdsqr.f"
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

#line 597 "sbdsqr.f"
	shift = 0.;
#line 598 "sbdsqr.f"
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

#line 602 "sbdsqr.f"
	if (idir == 1) {
#line 603 "sbdsqr.f"
	    sll = (d__1 = d__[ll], abs(d__1));
#line 604 "sbdsqr.f"
	    slas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
#line 605 "sbdsqr.f"
	} else {
#line 606 "sbdsqr.f"
	    sll = (d__1 = d__[m], abs(d__1));
#line 607 "sbdsqr.f"
	    slas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
#line 608 "sbdsqr.f"
	}

/*        Test if shift negligible, and if so set to zero */

#line 612 "sbdsqr.f"
	if (sll > 0.) {
/* Computing 2nd power */
#line 613 "sbdsqr.f"
	    d__1 = shift / sll;
#line 613 "sbdsqr.f"
	    if (d__1 * d__1 < eps) {
#line 613 "sbdsqr.f"
		shift = 0.;
#line 613 "sbdsqr.f"
	    }
#line 615 "sbdsqr.f"
	}
#line 616 "sbdsqr.f"
    }

/*     Increment iteration count */

#line 620 "sbdsqr.f"
    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

#line 624 "sbdsqr.f"
    if (shift == 0.) {
#line 625 "sbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 630 "sbdsqr.f"
	    cs = 1.;
#line 631 "sbdsqr.f"
	    oldcs = 1.;
#line 632 "sbdsqr.f"
	    i__1 = m - 1;
#line 632 "sbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 633 "sbdsqr.f"
		d__1 = d__[i__] * cs;
#line 633 "sbdsqr.f"
		slartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 634 "sbdsqr.f"
		if (i__ > ll) {
#line 634 "sbdsqr.f"
		    e[i__ - 1] = oldsn * r__;
#line 634 "sbdsqr.f"
		}
#line 636 "sbdsqr.f"
		d__1 = oldcs * r__;
#line 636 "sbdsqr.f"
		d__2 = d__[i__ + 1] * sn;
#line 636 "sbdsqr.f"
		slartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 637 "sbdsqr.f"
		work[i__ - ll + 1] = cs;
#line 638 "sbdsqr.f"
		work[i__ - ll + 1 + nm1] = sn;
#line 639 "sbdsqr.f"
		work[i__ - ll + 1 + nm12] = oldcs;
#line 640 "sbdsqr.f"
		work[i__ - ll + 1 + nm13] = oldsn;
#line 641 "sbdsqr.f"
/* L120: */
#line 641 "sbdsqr.f"
	    }
#line 642 "sbdsqr.f"
	    h__ = d__[m] * cs;
#line 643 "sbdsqr.f"
	    d__[m] = h__ * oldcs;
#line 644 "sbdsqr.f"
	    e[m - 1] = h__ * oldsn;

/*           Update singular vectors */

#line 648 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 648 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 648 "sbdsqr.f"
		slasr_("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 648 "sbdsqr.f"
	    }
#line 651 "sbdsqr.f"
	    if (*nru > 0) {
#line 651 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 651 "sbdsqr.f"
		slasr_("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 651 "sbdsqr.f"
	    }
#line 654 "sbdsqr.f"
	    if (*ncc > 0) {
#line 654 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 654 "sbdsqr.f"
		slasr_("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 654 "sbdsqr.f"
	    }

/*           Test convergence */

#line 660 "sbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 660 "sbdsqr.f"
		e[m - 1] = 0.;
#line 660 "sbdsqr.f"
	    }

#line 663 "sbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 668 "sbdsqr.f"
	    cs = 1.;
#line 669 "sbdsqr.f"
	    oldcs = 1.;
#line 670 "sbdsqr.f"
	    i__1 = ll + 1;
#line 670 "sbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 671 "sbdsqr.f"
		d__1 = d__[i__] * cs;
#line 671 "sbdsqr.f"
		slartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 672 "sbdsqr.f"
		if (i__ < m) {
#line 672 "sbdsqr.f"
		    e[i__] = oldsn * r__;
#line 672 "sbdsqr.f"
		}
#line 674 "sbdsqr.f"
		d__1 = oldcs * r__;
#line 674 "sbdsqr.f"
		d__2 = d__[i__ - 1] * sn;
#line 674 "sbdsqr.f"
		slartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 675 "sbdsqr.f"
		work[i__ - ll] = cs;
#line 676 "sbdsqr.f"
		work[i__ - ll + nm1] = -sn;
#line 677 "sbdsqr.f"
		work[i__ - ll + nm12] = oldcs;
#line 678 "sbdsqr.f"
		work[i__ - ll + nm13] = -oldsn;
#line 679 "sbdsqr.f"
/* L130: */
#line 679 "sbdsqr.f"
	    }
#line 680 "sbdsqr.f"
	    h__ = d__[ll] * cs;
#line 681 "sbdsqr.f"
	    d__[ll] = h__ * oldcs;
#line 682 "sbdsqr.f"
	    e[ll] = h__ * oldsn;

/*           Update singular vectors */

#line 686 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 686 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 686 "sbdsqr.f"
		slasr_("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 686 "sbdsqr.f"
	    }
#line 689 "sbdsqr.f"
	    if (*nru > 0) {
#line 689 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 689 "sbdsqr.f"
		slasr_("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 689 "sbdsqr.f"
	    }
#line 692 "sbdsqr.f"
	    if (*ncc > 0) {
#line 692 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 692 "sbdsqr.f"
		slasr_("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 692 "sbdsqr.f"
	    }

/*           Test convergence */

#line 698 "sbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 698 "sbdsqr.f"
		e[ll] = 0.;
#line 698 "sbdsqr.f"
	    }
#line 700 "sbdsqr.f"
	}
#line 701 "sbdsqr.f"
    } else {

/*        Use nonzero shift */

#line 705 "sbdsqr.f"
	if (idir == 1) {

/*           Chase bulge from top to bottom */
/*           Save cosines and sines for later singular vector updates */

#line 710 "sbdsqr.f"
	    f = ((d__1 = d__[ll], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[
		    ll]) + shift / d__[ll]);
#line 712 "sbdsqr.f"
	    g = e[ll];
#line 713 "sbdsqr.f"
	    i__1 = m - 1;
#line 713 "sbdsqr.f"
	    for (i__ = ll; i__ <= i__1; ++i__) {
#line 714 "sbdsqr.f"
		slartg_(&f, &g, &cosr, &sinr, &r__);
#line 715 "sbdsqr.f"
		if (i__ > ll) {
#line 715 "sbdsqr.f"
		    e[i__ - 1] = r__;
#line 715 "sbdsqr.f"
		}
#line 717 "sbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 718 "sbdsqr.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 719 "sbdsqr.f"
		g = sinr * d__[i__ + 1];
#line 720 "sbdsqr.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 721 "sbdsqr.f"
		slartg_(&f, &g, &cosl, &sinl, &r__);
#line 722 "sbdsqr.f"
		d__[i__] = r__;
#line 723 "sbdsqr.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 724 "sbdsqr.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 725 "sbdsqr.f"
		if (i__ < m - 1) {
#line 726 "sbdsqr.f"
		    g = sinl * e[i__ + 1];
#line 727 "sbdsqr.f"
		    e[i__ + 1] = cosl * e[i__ + 1];
#line 728 "sbdsqr.f"
		}
#line 729 "sbdsqr.f"
		work[i__ - ll + 1] = cosr;
#line 730 "sbdsqr.f"
		work[i__ - ll + 1 + nm1] = sinr;
#line 731 "sbdsqr.f"
		work[i__ - ll + 1 + nm12] = cosl;
#line 732 "sbdsqr.f"
		work[i__ - ll + 1 + nm13] = sinl;
#line 733 "sbdsqr.f"
/* L140: */
#line 733 "sbdsqr.f"
	    }
#line 734 "sbdsqr.f"
	    e[m - 1] = f;

/*           Update singular vectors */

#line 738 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 738 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 738 "sbdsqr.f"
		slasr_("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 738 "sbdsqr.f"
	    }
#line 741 "sbdsqr.f"
	    if (*nru > 0) {
#line 741 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 741 "sbdsqr.f"
		slasr_("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 741 "sbdsqr.f"
	    }
#line 744 "sbdsqr.f"
	    if (*ncc > 0) {
#line 744 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 744 "sbdsqr.f"
		slasr_("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 744 "sbdsqr.f"
	    }

/*           Test convergence */

#line 750 "sbdsqr.f"
	    if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
#line 750 "sbdsqr.f"
		e[m - 1] = 0.;
#line 750 "sbdsqr.f"
	    }

#line 753 "sbdsqr.f"
	} else {

/*           Chase bulge from bottom to top */
/*           Save cosines and sines for later singular vector updates */

#line 758 "sbdsqr.f"
	    f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_sign(&c_b49, &d__[m]
		    ) + shift / d__[m]);
#line 760 "sbdsqr.f"
	    g = e[m - 1];
#line 761 "sbdsqr.f"
	    i__1 = ll + 1;
#line 761 "sbdsqr.f"
	    for (i__ = m; i__ >= i__1; --i__) {
#line 762 "sbdsqr.f"
		slartg_(&f, &g, &cosr, &sinr, &r__);
#line 763 "sbdsqr.f"
		if (i__ < m) {
#line 763 "sbdsqr.f"
		    e[i__] = r__;
#line 763 "sbdsqr.f"
		}
#line 765 "sbdsqr.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 766 "sbdsqr.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 767 "sbdsqr.f"
		g = sinr * d__[i__ - 1];
#line 768 "sbdsqr.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 769 "sbdsqr.f"
		slartg_(&f, &g, &cosl, &sinl, &r__);
#line 770 "sbdsqr.f"
		d__[i__] = r__;
#line 771 "sbdsqr.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 772 "sbdsqr.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 773 "sbdsqr.f"
		if (i__ > ll + 1) {
#line 774 "sbdsqr.f"
		    g = sinl * e[i__ - 2];
#line 775 "sbdsqr.f"
		    e[i__ - 2] = cosl * e[i__ - 2];
#line 776 "sbdsqr.f"
		}
#line 777 "sbdsqr.f"
		work[i__ - ll] = cosr;
#line 778 "sbdsqr.f"
		work[i__ - ll + nm1] = -sinr;
#line 779 "sbdsqr.f"
		work[i__ - ll + nm12] = cosl;
#line 780 "sbdsqr.f"
		work[i__ - ll + nm13] = -sinl;
#line 781 "sbdsqr.f"
/* L150: */
#line 781 "sbdsqr.f"
	    }
#line 782 "sbdsqr.f"
	    e[ll] = f;

/*           Test convergence */

#line 786 "sbdsqr.f"
	    if ((d__1 = e[ll], abs(d__1)) <= thresh) {
#line 786 "sbdsqr.f"
		e[ll] = 0.;
#line 786 "sbdsqr.f"
	    }

/*           Update singular vectors if desired */

#line 791 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 791 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 791 "sbdsqr.f"
		slasr_("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 791 "sbdsqr.f"
	    }
#line 794 "sbdsqr.f"
	    if (*nru > 0) {
#line 794 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 794 "sbdsqr.f"
		slasr_("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 794 "sbdsqr.f"
	    }
#line 797 "sbdsqr.f"
	    if (*ncc > 0) {
#line 797 "sbdsqr.f"
		i__1 = m - ll + 1;
#line 797 "sbdsqr.f"
		slasr_("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 797 "sbdsqr.f"
	    }
#line 800 "sbdsqr.f"
	}
#line 801 "sbdsqr.f"
    }

/*     QR iteration finished, go back and check convergence */

#line 805 "sbdsqr.f"
    goto L60;

/*     All singular values converged, so make them positive */

#line 809 "sbdsqr.f"
L160:
#line 810 "sbdsqr.f"
    i__1 = *n;
#line 810 "sbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 811 "sbdsqr.f"
	if (d__[i__] < 0.) {
#line 812 "sbdsqr.f"
	    d__[i__] = -d__[i__];

/*           Change sign of singular vectors, if desired */

#line 816 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 816 "sbdsqr.f"
		sscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
#line 816 "sbdsqr.f"
	    }
#line 818 "sbdsqr.f"
	}
#line 819 "sbdsqr.f"
/* L170: */
#line 819 "sbdsqr.f"
    }

/*     Sort the singular values into decreasing order (insertion sort on */
/*     singular values, but only one transposition per singular vector) */

#line 824 "sbdsqr.f"
    i__1 = *n - 1;
#line 824 "sbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for smallest D(I) */

#line 828 "sbdsqr.f"
	isub = 1;
#line 829 "sbdsqr.f"
	smin = d__[1];
#line 830 "sbdsqr.f"
	i__2 = *n + 1 - i__;
#line 830 "sbdsqr.f"
	for (j = 2; j <= i__2; ++j) {
#line 831 "sbdsqr.f"
	    if (d__[j] <= smin) {
#line 832 "sbdsqr.f"
		isub = j;
#line 833 "sbdsqr.f"
		smin = d__[j];
#line 834 "sbdsqr.f"
	    }
#line 835 "sbdsqr.f"
/* L180: */
#line 835 "sbdsqr.f"
	}
#line 836 "sbdsqr.f"
	if (isub != *n + 1 - i__) {

/*           Swap singular values and vectors */

#line 840 "sbdsqr.f"
	    d__[isub] = d__[*n + 1 - i__];
#line 841 "sbdsqr.f"
	    d__[*n + 1 - i__] = smin;
#line 842 "sbdsqr.f"
	    if (*ncvt > 0) {
#line 842 "sbdsqr.f"
		sswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
#line 842 "sbdsqr.f"
	    }
#line 845 "sbdsqr.f"
	    if (*nru > 0) {
#line 845 "sbdsqr.f"
		sswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
#line 845 "sbdsqr.f"
	    }
#line 847 "sbdsqr.f"
	    if (*ncc > 0) {
#line 847 "sbdsqr.f"
		sswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
#line 847 "sbdsqr.f"
	    }
#line 849 "sbdsqr.f"
	}
#line 850 "sbdsqr.f"
/* L190: */
#line 850 "sbdsqr.f"
    }
#line 851 "sbdsqr.f"
    goto L220;

/*     Maximum number of iterations exceeded, failure to converge */

#line 855 "sbdsqr.f"
L200:
#line 856 "sbdsqr.f"
    *info = 0;
#line 857 "sbdsqr.f"
    i__1 = *n - 1;
#line 857 "sbdsqr.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 858 "sbdsqr.f"
	if (e[i__] != 0.) {
#line 858 "sbdsqr.f"
	    ++(*info);
#line 858 "sbdsqr.f"
	}
#line 860 "sbdsqr.f"
/* L210: */
#line 860 "sbdsqr.f"
    }
#line 861 "sbdsqr.f"
L220:
#line 862 "sbdsqr.f"
    return 0;

/*     End of SBDSQR */

} /* sbdsqr_ */

