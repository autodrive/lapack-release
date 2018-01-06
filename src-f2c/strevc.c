#line 1 "strevc.f"
/* strevc.f -- translated by f2c (version 20100827).
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

#line 1 "strevc.f"
/* Table of constant values */

static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b22 = 1.;
static doublereal c_b25 = 0.;
static integer c__2 = 2;
static logical c_true = TRUE_;

/* > \brief \b STREVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STREVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, MM, M, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL               T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STREVC computes some or all of the right and/or left eigenvectors of */
/* > a real upper quasi-triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* >    T*x = w*x,     (y**T)*T = w*(y**T) */
/* > */
/* > where y**T denotes the transpose of y. */
/* > The eigenvalues are not input to this routine, but are read directly */
/* > from the diagonal blocks of T. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/* > input matrix.  If Q is the orthogonal factor that reduces a matrix */
/* > A to Schur form T, then Q*X and Q*Y are the matrices of right and */
/* > left eigenvectors of A. */
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
/* >                  backtransformed by the matrices in VR and/or VL; */
/* >          = 'S':  compute selected right and/or left eigenvectors, */
/* >                  as indicated by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/* >          computed. */
/* >          If w(j) is a real eigenvalue, the corresponding real */
/* >          eigenvector is computed if SELECT(j) is .TRUE.. */
/* >          If w(j) and w(j+1) are the real and imaginary parts of a */
/* >          complex eigenvalue, the corresponding complex eigenvector is */
/* >          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and */
/* >          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to */
/* >          .FALSE.. */
/* >          Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,N) */
/* >          The upper quasi-triangular matrix T in Schur canonical form. */
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
/* >          VL is REAL array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of Schur vectors returned by SHSEQR). */
/* >          On exit, if SIDE = 'L' or 'B', VL contains: */
/* >          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*Y; */
/* >          if HOWMNY = 'S', the left eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VL, in the same order as their */
/* >                           eigenvalues. */
/* >          A complex eigenvector corresponding to a complex eigenvalue */
/* >          is stored in two consecutive columns, the first holding the */
/* >          real part, and the second the imaginary part. */
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
/* >          VR is REAL array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of Schur vectors returned by SHSEQR). */
/* >          On exit, if SIDE = 'R' or 'B', VR contains: */
/* >          if HOWMNY = 'A', the matrix X of right eigenvectors of T; */
/* >          if HOWMNY = 'B', the matrix Q*X; */
/* >          if HOWMNY = 'S', the right eigenvectors of T specified by */
/* >                           SELECT, stored consecutively in the columns */
/* >                           of VR, in the same order as their */
/* >                           eigenvalues. */
/* >          A complex eigenvector corresponding to a complex eigenvalue */
/* >          is stored in two consecutive columns, the first holding the */
/* >          real part and the second the imaginary part. */
/* >          Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR.  LDVR >= 1, and if */
/* >          SIDE = 'R' or 'B', LDVR >= N. */
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
/* >          used to store the eigenvectors. */
/* >          If HOWMNY = 'A' or 'B', M is set to N. */
/* >          Each selected real eigenvector occupies one column and each */
/* >          selected complex eigenvector occupies two columns. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (3*N) */
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

/* > \ingroup realOTHERcomputational */

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
/* Subroutine */ int strevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, 
	doublereal *work, integer *info, ftnlen side_len, ftnlen howmny_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal x[4]	/* was [2][2] */;
    static integer j1, j2, n2, ii, ki, ip, is;
    static doublereal wi, wr, rec, ulp, beta, emax;
    static logical pair, allv;
    static integer ierr;
    static doublereal unfl, ovfl, smin;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical over;
    static doublereal vmax;
    static integer jnxt;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal remax;
    static logical leftv;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical bothv;
    static doublereal vcrit;
    static logical somev;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal xnorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), slaln2_(logical *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , integer *), slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
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
/*     .. Local Arrays .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 273 "strevc.f"
    /* Parameter adjustments */
#line 273 "strevc.f"
    --select;
#line 273 "strevc.f"
    t_dim1 = *ldt;
#line 273 "strevc.f"
    t_offset = 1 + t_dim1;
#line 273 "strevc.f"
    t -= t_offset;
#line 273 "strevc.f"
    vl_dim1 = *ldvl;
#line 273 "strevc.f"
    vl_offset = 1 + vl_dim1;
#line 273 "strevc.f"
    vl -= vl_offset;
#line 273 "strevc.f"
    vr_dim1 = *ldvr;
#line 273 "strevc.f"
    vr_offset = 1 + vr_dim1;
#line 273 "strevc.f"
    vr -= vr_offset;
#line 273 "strevc.f"
    --work;
#line 273 "strevc.f"

#line 273 "strevc.f"
    /* Function Body */
#line 273 "strevc.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 274 "strevc.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 275 "strevc.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 277 "strevc.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 278 "strevc.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 279 "strevc.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 281 "strevc.f"
    *info = 0;
#line 282 "strevc.f"
    if (! rightv && ! leftv) {
#line 283 "strevc.f"
	*info = -1;
#line 284 "strevc.f"
    } else if (! allv && ! over && ! somev) {
#line 285 "strevc.f"
	*info = -2;
#line 286 "strevc.f"
    } else if (*n < 0) {
#line 287 "strevc.f"
	*info = -4;
#line 288 "strevc.f"
    } else if (*ldt < max(1,*n)) {
#line 289 "strevc.f"
	*info = -6;
#line 290 "strevc.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 291 "strevc.f"
	*info = -8;
#line 292 "strevc.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 293 "strevc.f"
	*info = -10;
#line 294 "strevc.f"
    } else {

/*        Set M to the number of columns required to store the selected */
/*        eigenvectors, standardize the array SELECT if necessary, and */
/*        test MM. */

#line 300 "strevc.f"
	if (somev) {
#line 301 "strevc.f"
	    *m = 0;
#line 302 "strevc.f"
	    pair = FALSE_;
#line 303 "strevc.f"
	    i__1 = *n;
#line 303 "strevc.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "strevc.f"
		if (pair) {
#line 305 "strevc.f"
		    pair = FALSE_;
#line 306 "strevc.f"
		    select[j] = FALSE_;
#line 307 "strevc.f"
		} else {
#line 308 "strevc.f"
		    if (j < *n) {
#line 309 "strevc.f"
			if (t[j + 1 + j * t_dim1] == 0.) {
#line 310 "strevc.f"
			    if (select[j]) {
#line 310 "strevc.f"
				++(*m);
#line 310 "strevc.f"
			    }
#line 312 "strevc.f"
			} else {
#line 313 "strevc.f"
			    pair = TRUE_;
#line 314 "strevc.f"
			    if (select[j] || select[j + 1]) {
#line 315 "strevc.f"
				select[j] = TRUE_;
#line 316 "strevc.f"
				*m += 2;
#line 317 "strevc.f"
			    }
#line 318 "strevc.f"
			}
#line 319 "strevc.f"
		    } else {
#line 320 "strevc.f"
			if (select[*n]) {
#line 320 "strevc.f"
			    ++(*m);
#line 320 "strevc.f"
			}
#line 322 "strevc.f"
		    }
#line 323 "strevc.f"
		}
#line 324 "strevc.f"
/* L10: */
#line 324 "strevc.f"
	    }
#line 325 "strevc.f"
	} else {
#line 326 "strevc.f"
	    *m = *n;
#line 327 "strevc.f"
	}

#line 329 "strevc.f"
	if (*mm < *m) {
#line 330 "strevc.f"
	    *info = -11;
#line 331 "strevc.f"
	}
#line 332 "strevc.f"
    }
#line 333 "strevc.f"
    if (*info != 0) {
#line 334 "strevc.f"
	i__1 = -(*info);
#line 334 "strevc.f"
	xerbla_("STREVC", &i__1, (ftnlen)6);
#line 335 "strevc.f"
	return 0;
#line 336 "strevc.f"
    }

/*     Quick return if possible. */

#line 340 "strevc.f"
    if (*n == 0) {
#line 340 "strevc.f"
	return 0;
#line 340 "strevc.f"
    }

/*     Set the constants to control overflow. */

#line 345 "strevc.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 346 "strevc.f"
    ovfl = 1. / unfl;
#line 347 "strevc.f"
    slabad_(&unfl, &ovfl);
#line 348 "strevc.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 349 "strevc.f"
    smlnum = unfl * (*n / ulp);
#line 350 "strevc.f"
    bignum = (1. - ulp) / smlnum;

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 355 "strevc.f"
    work[1] = 0.;
#line 356 "strevc.f"
    i__1 = *n;
#line 356 "strevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 357 "strevc.f"
	work[j] = 0.;
#line 358 "strevc.f"
	i__2 = j - 1;
#line 358 "strevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 359 "strevc.f"
	    work[j] += (d__1 = t[i__ + j * t_dim1], abs(d__1));
#line 360 "strevc.f"
/* L20: */
#line 360 "strevc.f"
	}
#line 361 "strevc.f"
/* L30: */
#line 361 "strevc.f"
    }

/*     Index IP is used to specify the real or complex eigenvalue: */
/*       IP = 0, real eigenvalue, */
/*            1, first of conjugate complex pair: (wr,wi) */
/*           -1, second of conjugate complex pair: (wr,wi) */

#line 368 "strevc.f"
    n2 = *n << 1;

#line 370 "strevc.f"
    if (rightv) {

/*        Compute right eigenvectors. */

#line 374 "strevc.f"
	ip = 0;
#line 375 "strevc.f"
	is = *m;
#line 376 "strevc.f"
	for (ki = *n; ki >= 1; --ki) {

#line 378 "strevc.f"
	    if (ip == 1) {
#line 378 "strevc.f"
		goto L130;
#line 378 "strevc.f"
	    }
#line 380 "strevc.f"
	    if (ki == 1) {
#line 380 "strevc.f"
		goto L40;
#line 380 "strevc.f"
	    }
#line 382 "strevc.f"
	    if (t[ki + (ki - 1) * t_dim1] == 0.) {
#line 382 "strevc.f"
		goto L40;
#line 382 "strevc.f"
	    }
#line 384 "strevc.f"
	    ip = -1;

#line 386 "strevc.f"
L40:
#line 387 "strevc.f"
	    if (somev) {
#line 388 "strevc.f"
		if (ip == 0) {
#line 389 "strevc.f"
		    if (! select[ki]) {
#line 389 "strevc.f"
			goto L130;
#line 389 "strevc.f"
		    }
#line 391 "strevc.f"
		} else {
#line 392 "strevc.f"
		    if (! select[ki - 1]) {
#line 392 "strevc.f"
			goto L130;
#line 392 "strevc.f"
		    }
#line 394 "strevc.f"
		}
#line 395 "strevc.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 399 "strevc.f"
	    wr = t[ki + ki * t_dim1];
#line 400 "strevc.f"
	    wi = 0.;
#line 401 "strevc.f"
	    if (ip != 0) {
#line 401 "strevc.f"
		wi = sqrt((d__1 = t[ki + (ki - 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki - 1 + ki * t_dim1], abs(d__2)));
#line 401 "strevc.f"
	    }
/* Computing MAX */
#line 404 "strevc.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 404 "strevc.f"
	    smin = max(d__1,smlnum);

#line 406 "strevc.f"
	    if (ip == 0) {

/*              Real right eigenvector */

#line 410 "strevc.f"
		work[ki + *n] = 1.;

/*              Form right-hand side */

#line 414 "strevc.f"
		i__1 = ki - 1;
#line 414 "strevc.f"
		for (k = 1; k <= i__1; ++k) {
#line 415 "strevc.f"
		    work[k + *n] = -t[k + ki * t_dim1];
#line 416 "strevc.f"
/* L50: */
#line 416 "strevc.f"
		}

/*              Solve the upper quasi-triangular system: */
/*                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */

#line 421 "strevc.f"
		jnxt = ki - 1;
#line 422 "strevc.f"
		for (j = ki - 1; j >= 1; --j) {
#line 423 "strevc.f"
		    if (j > jnxt) {
#line 423 "strevc.f"
			goto L60;
#line 423 "strevc.f"
		    }
#line 425 "strevc.f"
		    j1 = j;
#line 426 "strevc.f"
		    j2 = j;
#line 427 "strevc.f"
		    jnxt = j - 1;
#line 428 "strevc.f"
		    if (j > 1) {
#line 429 "strevc.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 430 "strevc.f"
			    j1 = j - 1;
#line 431 "strevc.f"
			    jnxt = j - 2;
#line 432 "strevc.f"
			}
#line 433 "strevc.f"
		    }

#line 435 "strevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 439 "strevc.f"
			slaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale X(1,1) to avoid overflow when updating */
/*                    the right-hand side. */

#line 446 "strevc.f"
			if (xnorm > 1.) {
#line 447 "strevc.f"
			    if (work[j] > bignum / xnorm) {
#line 448 "strevc.f"
				x[0] /= xnorm;
#line 449 "strevc.f"
				scale /= xnorm;
#line 450 "strevc.f"
			    }
#line 451 "strevc.f"
			}

/*                    Scale if necessary */

#line 455 "strevc.f"
			if (scale != 1.) {
#line 455 "strevc.f"
			    sscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 455 "strevc.f"
			}
#line 457 "strevc.f"
			work[j + *n] = x[0];

/*                    Update right-hand side */

#line 461 "strevc.f"
			i__1 = j - 1;
#line 461 "strevc.f"
			d__1 = -x[0];
#line 461 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);

#line 464 "strevc.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 468 "strevc.f"
			slaln2_(&c_false, &c__2, &c__1, &smin, &c_b22, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, &
				work[j - 1 + *n], n, &wr, &c_b25, x, &c__2, &
				scale, &xnorm, &ierr);

/*                    Scale X(1,1) and X(2,1) to avoid overflow when */
/*                    updating the right-hand side. */

#line 476 "strevc.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 477 "strevc.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 477 "strevc.f"
			    beta = max(d__1,d__2);
#line 478 "strevc.f"
			    if (beta > bignum / xnorm) {
#line 479 "strevc.f"
				x[0] /= xnorm;
#line 480 "strevc.f"
				x[1] /= xnorm;
#line 481 "strevc.f"
				scale /= xnorm;
#line 482 "strevc.f"
			    }
#line 483 "strevc.f"
			}

/*                    Scale if necessary */

#line 487 "strevc.f"
			if (scale != 1.) {
#line 487 "strevc.f"
			    sscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 487 "strevc.f"
			}
#line 489 "strevc.f"
			work[j - 1 + *n] = x[0];
#line 490 "strevc.f"
			work[j + *n] = x[1];

/*                    Update right-hand side */

#line 494 "strevc.f"
			i__1 = j - 2;
#line 494 "strevc.f"
			d__1 = -x[0];
#line 494 "strevc.f"
			saxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[*n + 1], &c__1);
#line 496 "strevc.f"
			i__1 = j - 2;
#line 496 "strevc.f"
			d__1 = -x[1];
#line 496 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 498 "strevc.f"
		    }
#line 499 "strevc.f"
L60:
#line 499 "strevc.f"
		    ;
#line 499 "strevc.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 503 "strevc.f"
		if (! over) {
#line 504 "strevc.f"
		    scopy_(&ki, &work[*n + 1], &c__1, &vr[is * vr_dim1 + 1], &
			    c__1);

#line 506 "strevc.f"
		    ii = isamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 507 "strevc.f"
		    remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs(d__1));
#line 508 "strevc.f"
		    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 510 "strevc.f"
		    i__1 = *n;
#line 510 "strevc.f"
		    for (k = ki + 1; k <= i__1; ++k) {
#line 511 "strevc.f"
			vr[k + is * vr_dim1] = 0.;
#line 512 "strevc.f"
/* L70: */
#line 512 "strevc.f"
		    }
#line 513 "strevc.f"
		} else {
#line 514 "strevc.f"
		    if (ki > 1) {
#line 514 "strevc.f"
			i__1 = ki - 1;
#line 514 "strevc.f"
			sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[*n + 1], &c__1, &work[ki + *n], &vr[ki * 
				vr_dim1 + 1], &c__1, (ftnlen)1);
#line 514 "strevc.f"
		    }

#line 519 "strevc.f"
		    ii = isamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 520 "strevc.f"
		    remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs(d__1));
#line 521 "strevc.f"
		    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
#line 522 "strevc.f"
		}

#line 524 "strevc.f"
	    } else {

/*              Complex right eigenvector. */

/*              Initial solve */
/*                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0. */
/*                [ (T(KI,KI-1)   T(KI,KI)   )               ] */

#line 532 "strevc.f"
		if ((d__1 = t[ki - 1 + ki * t_dim1], abs(d__1)) >= (d__2 = t[
			ki + (ki - 1) * t_dim1], abs(d__2))) {
#line 533 "strevc.f"
		    work[ki - 1 + *n] = 1.;
#line 534 "strevc.f"
		    work[ki + n2] = wi / t[ki - 1 + ki * t_dim1];
#line 535 "strevc.f"
		} else {
#line 536 "strevc.f"
		    work[ki - 1 + *n] = -wi / t[ki + (ki - 1) * t_dim1];
#line 537 "strevc.f"
		    work[ki + n2] = 1.;
#line 538 "strevc.f"
		}
#line 539 "strevc.f"
		work[ki + *n] = 0.;
#line 540 "strevc.f"
		work[ki - 1 + n2] = 0.;

/*              Form right-hand side */

#line 544 "strevc.f"
		i__1 = ki - 2;
#line 544 "strevc.f"
		for (k = 1; k <= i__1; ++k) {
#line 545 "strevc.f"
		    work[k + *n] = -work[ki - 1 + *n] * t[k + (ki - 1) * 
			    t_dim1];
#line 546 "strevc.f"
		    work[k + n2] = -work[ki + n2] * t[k + ki * t_dim1];
#line 547 "strevc.f"
/* L80: */
#line 547 "strevc.f"
		}

/*              Solve upper quasi-triangular system: */
/*              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */

#line 552 "strevc.f"
		jnxt = ki - 2;
#line 553 "strevc.f"
		for (j = ki - 2; j >= 1; --j) {
#line 554 "strevc.f"
		    if (j > jnxt) {
#line 554 "strevc.f"
			goto L90;
#line 554 "strevc.f"
		    }
#line 556 "strevc.f"
		    j1 = j;
#line 557 "strevc.f"
		    j2 = j;
#line 558 "strevc.f"
		    jnxt = j - 1;
#line 559 "strevc.f"
		    if (j > 1) {
#line 560 "strevc.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 561 "strevc.f"
			    j1 = j - 1;
#line 562 "strevc.f"
			    jnxt = j - 2;
#line 563 "strevc.f"
			}
#line 564 "strevc.f"
		    }

#line 566 "strevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 570 "strevc.f"
			slaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &wi, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale X(1,1) and X(1,2) to avoid overflow when */
/*                    updating the right-hand side. */

#line 577 "strevc.f"
			if (xnorm > 1.) {
#line 578 "strevc.f"
			    if (work[j] > bignum / xnorm) {
#line 579 "strevc.f"
				x[0] /= xnorm;
#line 580 "strevc.f"
				x[2] /= xnorm;
#line 581 "strevc.f"
				scale /= xnorm;
#line 582 "strevc.f"
			    }
#line 583 "strevc.f"
			}

/*                    Scale if necessary */

#line 587 "strevc.f"
			if (scale != 1.) {
#line 588 "strevc.f"
			    sscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 589 "strevc.f"
			    sscal_(&ki, &scale, &work[n2 + 1], &c__1);
#line 590 "strevc.f"
			}
#line 591 "strevc.f"
			work[j + *n] = x[0];
#line 592 "strevc.f"
			work[j + n2] = x[2];

/*                    Update the right-hand side */

#line 596 "strevc.f"
			i__1 = j - 1;
#line 596 "strevc.f"
			d__1 = -x[0];
#line 596 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 598 "strevc.f"
			i__1 = j - 1;
#line 598 "strevc.f"
			d__1 = -x[2];
#line 598 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				n2 + 1], &c__1);

#line 601 "strevc.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 605 "strevc.f"
			slaln2_(&c_false, &c__2, &c__2, &smin, &c_b22, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, &
				work[j - 1 + *n], n, &wr, &wi, x, &c__2, &
				scale, &xnorm, &ierr);

/*                    Scale X to avoid overflow when updating */
/*                    the right-hand side. */

#line 613 "strevc.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 614 "strevc.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 614 "strevc.f"
			    beta = max(d__1,d__2);
#line 615 "strevc.f"
			    if (beta > bignum / xnorm) {
#line 616 "strevc.f"
				rec = 1. / xnorm;
#line 617 "strevc.f"
				x[0] *= rec;
#line 618 "strevc.f"
				x[2] *= rec;
#line 619 "strevc.f"
				x[1] *= rec;
#line 620 "strevc.f"
				x[3] *= rec;
#line 621 "strevc.f"
				scale *= rec;
#line 622 "strevc.f"
			    }
#line 623 "strevc.f"
			}

/*                    Scale if necessary */

#line 627 "strevc.f"
			if (scale != 1.) {
#line 628 "strevc.f"
			    sscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 629 "strevc.f"
			    sscal_(&ki, &scale, &work[n2 + 1], &c__1);
#line 630 "strevc.f"
			}
#line 631 "strevc.f"
			work[j - 1 + *n] = x[0];
#line 632 "strevc.f"
			work[j + *n] = x[1];
#line 633 "strevc.f"
			work[j - 1 + n2] = x[2];
#line 634 "strevc.f"
			work[j + n2] = x[3];

/*                    Update the right-hand side */

#line 638 "strevc.f"
			i__1 = j - 2;
#line 638 "strevc.f"
			d__1 = -x[0];
#line 638 "strevc.f"
			saxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[*n + 1], &c__1);
#line 640 "strevc.f"
			i__1 = j - 2;
#line 640 "strevc.f"
			d__1 = -x[1];
#line 640 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 642 "strevc.f"
			i__1 = j - 2;
#line 642 "strevc.f"
			d__1 = -x[2];
#line 642 "strevc.f"
			saxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[n2 + 1], &c__1);
#line 644 "strevc.f"
			i__1 = j - 2;
#line 644 "strevc.f"
			d__1 = -x[3];
#line 644 "strevc.f"
			saxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				n2 + 1], &c__1);
#line 646 "strevc.f"
		    }
#line 647 "strevc.f"
L90:
#line 647 "strevc.f"
		    ;
#line 647 "strevc.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 651 "strevc.f"
		if (! over) {
#line 652 "strevc.f"
		    scopy_(&ki, &work[*n + 1], &c__1, &vr[(is - 1) * vr_dim1 
			    + 1], &c__1);
#line 653 "strevc.f"
		    scopy_(&ki, &work[n2 + 1], &c__1, &vr[is * vr_dim1 + 1], &
			    c__1);

#line 655 "strevc.f"
		    emax = 0.;
#line 656 "strevc.f"
		    i__1 = ki;
#line 656 "strevc.f"
		    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 657 "strevc.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + is * vr_dim1], 
				abs(d__2));
#line 657 "strevc.f"
			emax = max(d__3,d__4);
#line 659 "strevc.f"
/* L100: */
#line 659 "strevc.f"
		    }

#line 661 "strevc.f"
		    remax = 1. / emax;
#line 662 "strevc.f"
		    sscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
#line 663 "strevc.f"
		    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 665 "strevc.f"
		    i__1 = *n;
#line 665 "strevc.f"
		    for (k = ki + 1; k <= i__1; ++k) {
#line 666 "strevc.f"
			vr[k + (is - 1) * vr_dim1] = 0.;
#line 667 "strevc.f"
			vr[k + is * vr_dim1] = 0.;
#line 668 "strevc.f"
/* L110: */
#line 668 "strevc.f"
		    }

#line 670 "strevc.f"
		} else {

#line 672 "strevc.f"
		    if (ki > 2) {
#line 673 "strevc.f"
			i__1 = ki - 2;
#line 673 "strevc.f"
			sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[*n + 1], &c__1, &work[ki - 1 + *n], &vr[(
				ki - 1) * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 676 "strevc.f"
			i__1 = ki - 2;
#line 676 "strevc.f"
			sgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[n2 + 1], &c__1, &work[ki + n2], &vr[ki * 
				vr_dim1 + 1], &c__1, (ftnlen)1);
#line 679 "strevc.f"
		    } else {
#line 680 "strevc.f"
			sscal_(n, &work[ki - 1 + *n], &vr[(ki - 1) * vr_dim1 
				+ 1], &c__1);
#line 681 "strevc.f"
			sscal_(n, &work[ki + n2], &vr[ki * vr_dim1 + 1], &
				c__1);
#line 682 "strevc.f"
		    }

#line 684 "strevc.f"
		    emax = 0.;
#line 685 "strevc.f"
		    i__1 = *n;
#line 685 "strevc.f"
		    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 686 "strevc.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + ki * vr_dim1], 
				abs(d__2));
#line 686 "strevc.f"
			emax = max(d__3,d__4);
#line 688 "strevc.f"
/* L120: */
#line 688 "strevc.f"
		    }
#line 689 "strevc.f"
		    remax = 1. / emax;
#line 690 "strevc.f"
		    sscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
#line 691 "strevc.f"
		    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
#line 692 "strevc.f"
		}
#line 693 "strevc.f"
	    }

#line 695 "strevc.f"
	    --is;
#line 696 "strevc.f"
	    if (ip != 0) {
#line 696 "strevc.f"
		--is;
#line 696 "strevc.f"
	    }
#line 698 "strevc.f"
L130:
#line 699 "strevc.f"
	    if (ip == 1) {
#line 699 "strevc.f"
		ip = 0;
#line 699 "strevc.f"
	    }
#line 701 "strevc.f"
	    if (ip == -1) {
#line 701 "strevc.f"
		ip = 1;
#line 701 "strevc.f"
	    }
#line 703 "strevc.f"
/* L140: */
#line 703 "strevc.f"
	}
#line 704 "strevc.f"
    }

#line 706 "strevc.f"
    if (leftv) {

/*        Compute left eigenvectors. */

#line 710 "strevc.f"
	ip = 0;
#line 711 "strevc.f"
	is = 1;
#line 712 "strevc.f"
	i__1 = *n;
#line 712 "strevc.f"
	for (ki = 1; ki <= i__1; ++ki) {

#line 714 "strevc.f"
	    if (ip == -1) {
#line 714 "strevc.f"
		goto L250;
#line 714 "strevc.f"
	    }
#line 716 "strevc.f"
	    if (ki == *n) {
#line 716 "strevc.f"
		goto L150;
#line 716 "strevc.f"
	    }
#line 718 "strevc.f"
	    if (t[ki + 1 + ki * t_dim1] == 0.) {
#line 718 "strevc.f"
		goto L150;
#line 718 "strevc.f"
	    }
#line 720 "strevc.f"
	    ip = 1;

#line 722 "strevc.f"
L150:
#line 723 "strevc.f"
	    if (somev) {
#line 724 "strevc.f"
		if (! select[ki]) {
#line 724 "strevc.f"
		    goto L250;
#line 724 "strevc.f"
		}
#line 726 "strevc.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 730 "strevc.f"
	    wr = t[ki + ki * t_dim1];
#line 731 "strevc.f"
	    wi = 0.;
#line 732 "strevc.f"
	    if (ip != 0) {
#line 732 "strevc.f"
		wi = sqrt((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki + 1 + ki * t_dim1], abs(d__2)));
#line 732 "strevc.f"
	    }
/* Computing MAX */
#line 735 "strevc.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 735 "strevc.f"
	    smin = max(d__1,smlnum);

#line 737 "strevc.f"
	    if (ip == 0) {

/*              Real left eigenvector. */

#line 741 "strevc.f"
		work[ki + *n] = 1.;

/*              Form right-hand side */

#line 745 "strevc.f"
		i__2 = *n;
#line 745 "strevc.f"
		for (k = ki + 1; k <= i__2; ++k) {
#line 746 "strevc.f"
		    work[k + *n] = -t[ki + k * t_dim1];
#line 747 "strevc.f"
/* L160: */
#line 747 "strevc.f"
		}

/*              Solve the quasi-triangular system: */
/*                 (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK */

#line 752 "strevc.f"
		vmax = 1.;
#line 753 "strevc.f"
		vcrit = bignum;

#line 755 "strevc.f"
		jnxt = ki + 1;
#line 756 "strevc.f"
		i__2 = *n;
#line 756 "strevc.f"
		for (j = ki + 1; j <= i__2; ++j) {
#line 757 "strevc.f"
		    if (j < jnxt) {
#line 757 "strevc.f"
			goto L170;
#line 757 "strevc.f"
		    }
#line 759 "strevc.f"
		    j1 = j;
#line 760 "strevc.f"
		    j2 = j;
#line 761 "strevc.f"
		    jnxt = j + 1;
#line 762 "strevc.f"
		    if (j < *n) {
#line 763 "strevc.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 764 "strevc.f"
			    j2 = j + 1;
#line 765 "strevc.f"
			    jnxt = j + 2;
#line 766 "strevc.f"
			}
#line 767 "strevc.f"
		    }

#line 769 "strevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

#line 776 "strevc.f"
			if (work[j] > vcrit) {
#line 777 "strevc.f"
			    rec = 1. / vmax;
#line 778 "strevc.f"
			    i__3 = *n - ki + 1;
#line 778 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 779 "strevc.f"
			    vmax = 1.;
#line 780 "strevc.f"
			    vcrit = bignum;
#line 781 "strevc.f"
			}

#line 783 "strevc.f"
			i__3 = j - ki - 1;
#line 783 "strevc.f"
			work[j + *n] -= sdot_(&i__3, &t[ki + 1 + j * t_dim1], 
				&c__1, &work[ki + 1 + *n], &c__1);

/*                    Solve (T(J,J)-WR)**T*X = WORK */

#line 789 "strevc.f"
			slaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale if necessary */

#line 795 "strevc.f"
			if (scale != 1.) {
#line 795 "strevc.f"
			    i__3 = *n - ki + 1;
#line 795 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 795 "strevc.f"
			}
#line 797 "strevc.f"
			work[j + *n] = x[0];
/* Computing MAX */
#line 798 "strevc.f"
			d__2 = (d__1 = work[j + *n], abs(d__1));
#line 798 "strevc.f"
			vmax = max(d__2,vmax);
#line 799 "strevc.f"
			vcrit = bignum / vmax;

#line 801 "strevc.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

/* Computing MAX */
#line 808 "strevc.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 808 "strevc.f"
			beta = max(d__1,d__2);
#line 809 "strevc.f"
			if (beta > vcrit) {
#line 810 "strevc.f"
			    rec = 1. / vmax;
#line 811 "strevc.f"
			    i__3 = *n - ki + 1;
#line 811 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 812 "strevc.f"
			    vmax = 1.;
#line 813 "strevc.f"
			    vcrit = bignum;
#line 814 "strevc.f"
			}

#line 816 "strevc.f"
			i__3 = j - ki - 1;
#line 816 "strevc.f"
			work[j + *n] -= sdot_(&i__3, &t[ki + 1 + j * t_dim1], 
				&c__1, &work[ki + 1 + *n], &c__1);

#line 820 "strevc.f"
			i__3 = j - ki - 1;
#line 820 "strevc.f"
			work[j + 1 + *n] -= sdot_(&i__3, &t[ki + 1 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 1 + *n], &c__1);

/*                    Solve */
/*                      [T(J,J)-WR   T(J,J+1)     ]**T* X = SCALE*( WORK1 ) */
/*                      [T(J+1,J)    T(J+1,J+1)-WR]               ( WORK2 ) */

#line 828 "strevc.f"
			slaln2_(&c_true, &c__2, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale if necessary */

#line 834 "strevc.f"
			if (scale != 1.) {
#line 834 "strevc.f"
			    i__3 = *n - ki + 1;
#line 834 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 834 "strevc.f"
			}
#line 836 "strevc.f"
			work[j + *n] = x[0];
#line 837 "strevc.f"
			work[j + 1 + *n] = x[1];

/* Computing MAX */
#line 839 "strevc.f"
			d__3 = (d__1 = work[j + *n], abs(d__1)), d__4 = (d__2 
				= work[j + 1 + *n], abs(d__2)), d__3 = max(
				d__3,d__4);
#line 839 "strevc.f"
			vmax = max(d__3,vmax);
#line 841 "strevc.f"
			vcrit = bignum / vmax;

#line 843 "strevc.f"
		    }
#line 844 "strevc.f"
L170:
#line 844 "strevc.f"
		    ;
#line 844 "strevc.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 848 "strevc.f"
		if (! over) {
#line 849 "strevc.f"
		    i__2 = *n - ki + 1;
#line 849 "strevc.f"
		    scopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);

#line 851 "strevc.f"
		    i__2 = *n - ki + 1;
#line 851 "strevc.f"
		    ii = isamax_(&i__2, &vl[ki + is * vl_dim1], &c__1) + ki - 
			    1;
#line 852 "strevc.f"
		    remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs(d__1));
#line 853 "strevc.f"
		    i__2 = *n - ki + 1;
#line 853 "strevc.f"
		    sscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 855 "strevc.f"
		    i__2 = ki - 1;
#line 855 "strevc.f"
		    for (k = 1; k <= i__2; ++k) {
#line 856 "strevc.f"
			vl[k + is * vl_dim1] = 0.;
#line 857 "strevc.f"
/* L180: */
#line 857 "strevc.f"
		    }

#line 859 "strevc.f"
		} else {

#line 861 "strevc.f"
		    if (ki < *n) {
#line 861 "strevc.f"
			i__2 = *n - ki;
#line 861 "strevc.f"
			sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 1) * vl_dim1 
				+ 1], ldvl, &work[ki + 1 + *n], &c__1, &work[
				ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (
				ftnlen)1);
#line 861 "strevc.f"
		    }

#line 866 "strevc.f"
		    ii = isamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 867 "strevc.f"
		    remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs(d__1));
#line 868 "strevc.f"
		    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

#line 870 "strevc.f"
		}

#line 872 "strevc.f"
	    } else {

/*              Complex left eigenvector. */

/*               Initial solve: */
/*                 ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0. */
/*                 ((T(KI+1,KI) T(KI+1,KI+1))                ) */

#line 880 "strevc.f"
		if ((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1)) >= (d__2 = 
			t[ki + 1 + ki * t_dim1], abs(d__2))) {
#line 881 "strevc.f"
		    work[ki + *n] = wi / t[ki + (ki + 1) * t_dim1];
#line 882 "strevc.f"
		    work[ki + 1 + n2] = 1.;
#line 883 "strevc.f"
		} else {
#line 884 "strevc.f"
		    work[ki + *n] = 1.;
#line 885 "strevc.f"
		    work[ki + 1 + n2] = -wi / t[ki + 1 + ki * t_dim1];
#line 886 "strevc.f"
		}
#line 887 "strevc.f"
		work[ki + 1 + *n] = 0.;
#line 888 "strevc.f"
		work[ki + n2] = 0.;

/*              Form right-hand side */

#line 892 "strevc.f"
		i__2 = *n;
#line 892 "strevc.f"
		for (k = ki + 2; k <= i__2; ++k) {
#line 893 "strevc.f"
		    work[k + *n] = -work[ki + *n] * t[ki + k * t_dim1];
#line 894 "strevc.f"
		    work[k + n2] = -work[ki + 1 + n2] * t[ki + 1 + k * t_dim1]
			    ;
#line 895 "strevc.f"
/* L190: */
#line 895 "strevc.f"
		}

/*              Solve complex quasi-triangular system: */
/*              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */

#line 900 "strevc.f"
		vmax = 1.;
#line 901 "strevc.f"
		vcrit = bignum;

#line 903 "strevc.f"
		jnxt = ki + 2;
#line 904 "strevc.f"
		i__2 = *n;
#line 904 "strevc.f"
		for (j = ki + 2; j <= i__2; ++j) {
#line 905 "strevc.f"
		    if (j < jnxt) {
#line 905 "strevc.f"
			goto L200;
#line 905 "strevc.f"
		    }
#line 907 "strevc.f"
		    j1 = j;
#line 908 "strevc.f"
		    j2 = j;
#line 909 "strevc.f"
		    jnxt = j + 1;
#line 910 "strevc.f"
		    if (j < *n) {
#line 911 "strevc.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 912 "strevc.f"
			    j2 = j + 1;
#line 913 "strevc.f"
			    jnxt = j + 2;
#line 914 "strevc.f"
			}
#line 915 "strevc.f"
		    }

#line 917 "strevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when */
/*                    forming the right-hand side elements. */

#line 924 "strevc.f"
			if (work[j] > vcrit) {
#line 925 "strevc.f"
			    rec = 1. / vmax;
#line 926 "strevc.f"
			    i__3 = *n - ki + 1;
#line 926 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 927 "strevc.f"
			    i__3 = *n - ki + 1;
#line 927 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + n2], &c__1);
#line 928 "strevc.f"
			    vmax = 1.;
#line 929 "strevc.f"
			    vcrit = bignum;
#line 930 "strevc.f"
			}

#line 932 "strevc.f"
			i__3 = j - ki - 2;
#line 932 "strevc.f"
			work[j + *n] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + *n], &c__1);
#line 935 "strevc.f"
			i__3 = j - ki - 2;
#line 935 "strevc.f"
			work[j + n2] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + n2], &c__1);

/*                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */

#line 941 "strevc.f"
			d__1 = -wi;
#line 941 "strevc.f"
			slaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &d__1, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale if necessary */

#line 947 "strevc.f"
			if (scale != 1.) {
#line 948 "strevc.f"
			    i__3 = *n - ki + 1;
#line 948 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 949 "strevc.f"
			    i__3 = *n - ki + 1;
#line 949 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + n2], &c__1);
#line 950 "strevc.f"
			}
#line 951 "strevc.f"
			work[j + *n] = x[0];
#line 952 "strevc.f"
			work[j + n2] = x[2];
/* Computing MAX */
#line 953 "strevc.f"
			d__3 = (d__1 = work[j + *n], abs(d__1)), d__4 = (d__2 
				= work[j + n2], abs(d__2)), d__3 = max(d__3,
				d__4);
#line 953 "strevc.f"
			vmax = max(d__3,vmax);
#line 955 "strevc.f"
			vcrit = bignum / vmax;

#line 957 "strevc.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side elements. */

/* Computing MAX */
#line 964 "strevc.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 964 "strevc.f"
			beta = max(d__1,d__2);
#line 965 "strevc.f"
			if (beta > vcrit) {
#line 966 "strevc.f"
			    rec = 1. / vmax;
#line 967 "strevc.f"
			    i__3 = *n - ki + 1;
#line 967 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 968 "strevc.f"
			    i__3 = *n - ki + 1;
#line 968 "strevc.f"
			    sscal_(&i__3, &rec, &work[ki + n2], &c__1);
#line 969 "strevc.f"
			    vmax = 1.;
#line 970 "strevc.f"
			    vcrit = bignum;
#line 971 "strevc.f"
			}

#line 973 "strevc.f"
			i__3 = j - ki - 2;
#line 973 "strevc.f"
			work[j + *n] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + *n], &c__1);

#line 977 "strevc.f"
			i__3 = j - ki - 2;
#line 977 "strevc.f"
			work[j + n2] -= sdot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + n2], &c__1);

#line 981 "strevc.f"
			i__3 = j - ki - 2;
#line 981 "strevc.f"
			work[j + 1 + *n] -= sdot_(&i__3, &t[ki + 2 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 2 + *n], &c__1);

#line 985 "strevc.f"
			i__3 = j - ki - 2;
#line 985 "strevc.f"
			work[j + 1 + n2] -= sdot_(&i__3, &t[ki + 2 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 2 + n2], &c__1);

/*                    Solve 2-by-2 complex linear equation */
/*                      ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B */
/*                      ([T(j+1,j) T(j+1,j+1)]               ) */

#line 993 "strevc.f"
			d__1 = -wi;
#line 993 "strevc.f"
			slaln2_(&c_true, &c__2, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &d__1, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale if necessary */

#line 999 "strevc.f"
			if (scale != 1.) {
#line 1000 "strevc.f"
			    i__3 = *n - ki + 1;
#line 1000 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 1001 "strevc.f"
			    i__3 = *n - ki + 1;
#line 1001 "strevc.f"
			    sscal_(&i__3, &scale, &work[ki + n2], &c__1);
#line 1002 "strevc.f"
			}
#line 1003 "strevc.f"
			work[j + *n] = x[0];
#line 1004 "strevc.f"
			work[j + n2] = x[2];
#line 1005 "strevc.f"
			work[j + 1 + *n] = x[1];
#line 1006 "strevc.f"
			work[j + 1 + n2] = x[3];
/* Computing MAX */
#line 1007 "strevc.f"
			d__1 = abs(x[0]), d__2 = abs(x[2]), d__1 = max(d__1,
				d__2), d__2 = abs(x[1]), d__1 = max(d__1,d__2)
				, d__2 = abs(x[3]), d__1 = max(d__1,d__2);
#line 1007 "strevc.f"
			vmax = max(d__1,vmax);
#line 1009 "strevc.f"
			vcrit = bignum / vmax;

#line 1011 "strevc.f"
		    }
#line 1012 "strevc.f"
L200:
#line 1012 "strevc.f"
		    ;
#line 1012 "strevc.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 1016 "strevc.f"
		if (! over) {
#line 1017 "strevc.f"
		    i__2 = *n - ki + 1;
#line 1017 "strevc.f"
		    scopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);
#line 1018 "strevc.f"
		    i__2 = *n - ki + 1;
#line 1018 "strevc.f"
		    scopy_(&i__2, &work[ki + n2], &c__1, &vl[ki + (is + 1) * 
			    vl_dim1], &c__1);

#line 1021 "strevc.f"
		    emax = 0.;
#line 1022 "strevc.f"
		    i__2 = *n;
#line 1022 "strevc.f"
		    for (k = ki; k <= i__2; ++k) {
/* Computing MAX */
#line 1023 "strevc.f"
			d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (is + 1) * vl_dim1], 
				abs(d__2));
#line 1023 "strevc.f"
			emax = max(d__3,d__4);
#line 1025 "strevc.f"
/* L220: */
#line 1025 "strevc.f"
		    }
#line 1026 "strevc.f"
		    remax = 1. / emax;
#line 1027 "strevc.f"
		    i__2 = *n - ki + 1;
#line 1027 "strevc.f"
		    sscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);
#line 1028 "strevc.f"
		    i__2 = *n - ki + 1;
#line 1028 "strevc.f"
		    sscal_(&i__2, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1)
			    ;

#line 1030 "strevc.f"
		    i__2 = ki - 1;
#line 1030 "strevc.f"
		    for (k = 1; k <= i__2; ++k) {
#line 1031 "strevc.f"
			vl[k + is * vl_dim1] = 0.;
#line 1032 "strevc.f"
			vl[k + (is + 1) * vl_dim1] = 0.;
#line 1033 "strevc.f"
/* L230: */
#line 1033 "strevc.f"
		    }
#line 1034 "strevc.f"
		} else {
#line 1035 "strevc.f"
		    if (ki < *n - 1) {
#line 1036 "strevc.f"
			i__2 = *n - ki - 1;
#line 1036 "strevc.f"
			sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + *n], &c__1, &work[
				ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (
				ftnlen)1);
#line 1039 "strevc.f"
			i__2 = *n - ki - 1;
#line 1039 "strevc.f"
			sgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + n2], &c__1, &work[
				ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1042 "strevc.f"
		    } else {
#line 1043 "strevc.f"
			sscal_(n, &work[ki + *n], &vl[ki * vl_dim1 + 1], &
				c__1);
#line 1044 "strevc.f"
			sscal_(n, &work[ki + 1 + n2], &vl[(ki + 1) * vl_dim1 
				+ 1], &c__1);
#line 1045 "strevc.f"
		    }

#line 1047 "strevc.f"
		    emax = 0.;
#line 1048 "strevc.f"
		    i__2 = *n;
#line 1048 "strevc.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 1049 "strevc.f"
			d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (ki + 1) * vl_dim1], 
				abs(d__2));
#line 1049 "strevc.f"
			emax = max(d__3,d__4);
#line 1051 "strevc.f"
/* L240: */
#line 1051 "strevc.f"
		    }
#line 1052 "strevc.f"
		    remax = 1. / emax;
#line 1053 "strevc.f"
		    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
#line 1054 "strevc.f"
		    sscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);

#line 1056 "strevc.f"
		}

#line 1058 "strevc.f"
	    }

#line 1060 "strevc.f"
	    ++is;
#line 1061 "strevc.f"
	    if (ip != 0) {
#line 1061 "strevc.f"
		++is;
#line 1061 "strevc.f"
	    }
#line 1063 "strevc.f"
L250:
#line 1064 "strevc.f"
	    if (ip == -1) {
#line 1064 "strevc.f"
		ip = 0;
#line 1064 "strevc.f"
	    }
#line 1066 "strevc.f"
	    if (ip == 1) {
#line 1066 "strevc.f"
		ip = -1;
#line 1066 "strevc.f"
	    }

#line 1069 "strevc.f"
/* L260: */
#line 1069 "strevc.f"
	}

#line 1071 "strevc.f"
    }

#line 1073 "strevc.f"
    return 0;

/*     End of STREVC */

} /* strevc_ */

