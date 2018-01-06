#line 1 "dtrevc.f"
/* dtrevc.f -- translated by f2c (version 20100827).
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

#line 1 "dtrevc.f"
/* Table of constant values */

static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b22 = 1.;
static doublereal c_b25 = 0.;
static integer c__2 = 2;
static logical c_true = TRUE_;

/* > \brief \b DTREVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTREVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/*                          LDVR, MM, M, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTREVC computes some or all of the right and/or left eigenvectors of */
/* > a real upper quasi-triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* >    T*x = w*x,     (y**H)*T = w*(y**H) */
/* > */
/* > where y**H denotes the conjugate transpose of y. */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,N) */
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
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of Schur vectors returned by DHSEQR). */
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
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of Schur vectors returned by DHSEQR). */
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
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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

/* > \date November 2017 */

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dtrevc_(char *side, char *howmny, logical *select, 
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
    static logical pair;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical allv;
    static integer ierr;
    static doublereal unfl, ovfl, smin;
    static logical over;
    static doublereal vmax;
    static integer jnxt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal remax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical leftv, bothv;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal vcrit;
    static logical somev;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical rightv;
    static doublereal smlnum;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

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

#line 273 "dtrevc.f"
    /* Parameter adjustments */
#line 273 "dtrevc.f"
    --select;
#line 273 "dtrevc.f"
    t_dim1 = *ldt;
#line 273 "dtrevc.f"
    t_offset = 1 + t_dim1;
#line 273 "dtrevc.f"
    t -= t_offset;
#line 273 "dtrevc.f"
    vl_dim1 = *ldvl;
#line 273 "dtrevc.f"
    vl_offset = 1 + vl_dim1;
#line 273 "dtrevc.f"
    vl -= vl_offset;
#line 273 "dtrevc.f"
    vr_dim1 = *ldvr;
#line 273 "dtrevc.f"
    vr_offset = 1 + vr_dim1;
#line 273 "dtrevc.f"
    vr -= vr_offset;
#line 273 "dtrevc.f"
    --work;
#line 273 "dtrevc.f"

#line 273 "dtrevc.f"
    /* Function Body */
#line 273 "dtrevc.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 274 "dtrevc.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 275 "dtrevc.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 277 "dtrevc.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 278 "dtrevc.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 279 "dtrevc.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 281 "dtrevc.f"
    *info = 0;
#line 282 "dtrevc.f"
    if (! rightv && ! leftv) {
#line 283 "dtrevc.f"
	*info = -1;
#line 284 "dtrevc.f"
    } else if (! allv && ! over && ! somev) {
#line 285 "dtrevc.f"
	*info = -2;
#line 286 "dtrevc.f"
    } else if (*n < 0) {
#line 287 "dtrevc.f"
	*info = -4;
#line 288 "dtrevc.f"
    } else if (*ldt < max(1,*n)) {
#line 289 "dtrevc.f"
	*info = -6;
#line 290 "dtrevc.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 291 "dtrevc.f"
	*info = -8;
#line 292 "dtrevc.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 293 "dtrevc.f"
	*info = -10;
#line 294 "dtrevc.f"
    } else {

/*        Set M to the number of columns required to store the selected */
/*        eigenvectors, standardize the array SELECT if necessary, and */
/*        test MM. */

#line 300 "dtrevc.f"
	if (somev) {
#line 301 "dtrevc.f"
	    *m = 0;
#line 302 "dtrevc.f"
	    pair = FALSE_;
#line 303 "dtrevc.f"
	    i__1 = *n;
#line 303 "dtrevc.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "dtrevc.f"
		if (pair) {
#line 305 "dtrevc.f"
		    pair = FALSE_;
#line 306 "dtrevc.f"
		    select[j] = FALSE_;
#line 307 "dtrevc.f"
		} else {
#line 308 "dtrevc.f"
		    if (j < *n) {
#line 309 "dtrevc.f"
			if (t[j + 1 + j * t_dim1] == 0.) {
#line 310 "dtrevc.f"
			    if (select[j]) {
#line 310 "dtrevc.f"
				++(*m);
#line 310 "dtrevc.f"
			    }
#line 312 "dtrevc.f"
			} else {
#line 313 "dtrevc.f"
			    pair = TRUE_;
#line 314 "dtrevc.f"
			    if (select[j] || select[j + 1]) {
#line 315 "dtrevc.f"
				select[j] = TRUE_;
#line 316 "dtrevc.f"
				*m += 2;
#line 317 "dtrevc.f"
			    }
#line 318 "dtrevc.f"
			}
#line 319 "dtrevc.f"
		    } else {
#line 320 "dtrevc.f"
			if (select[*n]) {
#line 320 "dtrevc.f"
			    ++(*m);
#line 320 "dtrevc.f"
			}
#line 322 "dtrevc.f"
		    }
#line 323 "dtrevc.f"
		}
#line 324 "dtrevc.f"
/* L10: */
#line 324 "dtrevc.f"
	    }
#line 325 "dtrevc.f"
	} else {
#line 326 "dtrevc.f"
	    *m = *n;
#line 327 "dtrevc.f"
	}

#line 329 "dtrevc.f"
	if (*mm < *m) {
#line 330 "dtrevc.f"
	    *info = -11;
#line 331 "dtrevc.f"
	}
#line 332 "dtrevc.f"
    }
#line 333 "dtrevc.f"
    if (*info != 0) {
#line 334 "dtrevc.f"
	i__1 = -(*info);
#line 334 "dtrevc.f"
	xerbla_("DTREVC", &i__1, (ftnlen)6);
#line 335 "dtrevc.f"
	return 0;
#line 336 "dtrevc.f"
    }

/*     Quick return if possible. */

#line 340 "dtrevc.f"
    if (*n == 0) {
#line 340 "dtrevc.f"
	return 0;
#line 340 "dtrevc.f"
    }

/*     Set the constants to control overflow. */

#line 345 "dtrevc.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 346 "dtrevc.f"
    ovfl = 1. / unfl;
#line 347 "dtrevc.f"
    dlabad_(&unfl, &ovfl);
#line 348 "dtrevc.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 349 "dtrevc.f"
    smlnum = unfl * (*n / ulp);
#line 350 "dtrevc.f"
    bignum = (1. - ulp) / smlnum;

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 355 "dtrevc.f"
    work[1] = 0.;
#line 356 "dtrevc.f"
    i__1 = *n;
#line 356 "dtrevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 357 "dtrevc.f"
	work[j] = 0.;
#line 358 "dtrevc.f"
	i__2 = j - 1;
#line 358 "dtrevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 359 "dtrevc.f"
	    work[j] += (d__1 = t[i__ + j * t_dim1], abs(d__1));
#line 360 "dtrevc.f"
/* L20: */
#line 360 "dtrevc.f"
	}
#line 361 "dtrevc.f"
/* L30: */
#line 361 "dtrevc.f"
    }

/*     Index IP is used to specify the real or complex eigenvalue: */
/*       IP = 0, real eigenvalue, */
/*            1, first of conjugate complex pair: (wr,wi) */
/*           -1, second of conjugate complex pair: (wr,wi) */

#line 368 "dtrevc.f"
    n2 = *n << 1;

#line 370 "dtrevc.f"
    if (rightv) {

/*        Compute right eigenvectors. */

#line 374 "dtrevc.f"
	ip = 0;
#line 375 "dtrevc.f"
	is = *m;
#line 376 "dtrevc.f"
	for (ki = *n; ki >= 1; --ki) {

#line 378 "dtrevc.f"
	    if (ip == 1) {
#line 378 "dtrevc.f"
		goto L130;
#line 378 "dtrevc.f"
	    }
#line 380 "dtrevc.f"
	    if (ki == 1) {
#line 380 "dtrevc.f"
		goto L40;
#line 380 "dtrevc.f"
	    }
#line 382 "dtrevc.f"
	    if (t[ki + (ki - 1) * t_dim1] == 0.) {
#line 382 "dtrevc.f"
		goto L40;
#line 382 "dtrevc.f"
	    }
#line 384 "dtrevc.f"
	    ip = -1;

#line 386 "dtrevc.f"
L40:
#line 387 "dtrevc.f"
	    if (somev) {
#line 388 "dtrevc.f"
		if (ip == 0) {
#line 389 "dtrevc.f"
		    if (! select[ki]) {
#line 389 "dtrevc.f"
			goto L130;
#line 389 "dtrevc.f"
		    }
#line 391 "dtrevc.f"
		} else {
#line 392 "dtrevc.f"
		    if (! select[ki - 1]) {
#line 392 "dtrevc.f"
			goto L130;
#line 392 "dtrevc.f"
		    }
#line 394 "dtrevc.f"
		}
#line 395 "dtrevc.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 399 "dtrevc.f"
	    wr = t[ki + ki * t_dim1];
#line 400 "dtrevc.f"
	    wi = 0.;
#line 401 "dtrevc.f"
	    if (ip != 0) {
#line 401 "dtrevc.f"
		wi = sqrt((d__1 = t[ki + (ki - 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki - 1 + ki * t_dim1], abs(d__2)));
#line 401 "dtrevc.f"
	    }
/* Computing MAX */
#line 404 "dtrevc.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 404 "dtrevc.f"
	    smin = max(d__1,smlnum);

#line 406 "dtrevc.f"
	    if (ip == 0) {

/*              Real right eigenvector */

#line 410 "dtrevc.f"
		work[ki + *n] = 1.;

/*              Form right-hand side */

#line 414 "dtrevc.f"
		i__1 = ki - 1;
#line 414 "dtrevc.f"
		for (k = 1; k <= i__1; ++k) {
#line 415 "dtrevc.f"
		    work[k + *n] = -t[k + ki * t_dim1];
#line 416 "dtrevc.f"
/* L50: */
#line 416 "dtrevc.f"
		}

/*              Solve the upper quasi-triangular system: */
/*                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */

#line 421 "dtrevc.f"
		jnxt = ki - 1;
#line 422 "dtrevc.f"
		for (j = ki - 1; j >= 1; --j) {
#line 423 "dtrevc.f"
		    if (j > jnxt) {
#line 423 "dtrevc.f"
			goto L60;
#line 423 "dtrevc.f"
		    }
#line 425 "dtrevc.f"
		    j1 = j;
#line 426 "dtrevc.f"
		    j2 = j;
#line 427 "dtrevc.f"
		    jnxt = j - 1;
#line 428 "dtrevc.f"
		    if (j > 1) {
#line 429 "dtrevc.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 430 "dtrevc.f"
			    j1 = j - 1;
#line 431 "dtrevc.f"
			    jnxt = j - 2;
#line 432 "dtrevc.f"
			}
#line 433 "dtrevc.f"
		    }

#line 435 "dtrevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 439 "dtrevc.f"
			dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale X(1,1) to avoid overflow when updating */
/*                    the right-hand side. */

#line 446 "dtrevc.f"
			if (xnorm > 1.) {
#line 447 "dtrevc.f"
			    if (work[j] > bignum / xnorm) {
#line 448 "dtrevc.f"
				x[0] /= xnorm;
#line 449 "dtrevc.f"
				scale /= xnorm;
#line 450 "dtrevc.f"
			    }
#line 451 "dtrevc.f"
			}

/*                    Scale if necessary */

#line 455 "dtrevc.f"
			if (scale != 1.) {
#line 455 "dtrevc.f"
			    dscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 455 "dtrevc.f"
			}
#line 457 "dtrevc.f"
			work[j + *n] = x[0];

/*                    Update right-hand side */

#line 461 "dtrevc.f"
			i__1 = j - 1;
#line 461 "dtrevc.f"
			d__1 = -x[0];
#line 461 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);

#line 464 "dtrevc.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 468 "dtrevc.f"
			dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b22, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, &
				work[j - 1 + *n], n, &wr, &c_b25, x, &c__2, &
				scale, &xnorm, &ierr);

/*                    Scale X(1,1) and X(2,1) to avoid overflow when */
/*                    updating the right-hand side. */

#line 476 "dtrevc.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 477 "dtrevc.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 477 "dtrevc.f"
			    beta = max(d__1,d__2);
#line 478 "dtrevc.f"
			    if (beta > bignum / xnorm) {
#line 479 "dtrevc.f"
				x[0] /= xnorm;
#line 480 "dtrevc.f"
				x[1] /= xnorm;
#line 481 "dtrevc.f"
				scale /= xnorm;
#line 482 "dtrevc.f"
			    }
#line 483 "dtrevc.f"
			}

/*                    Scale if necessary */

#line 487 "dtrevc.f"
			if (scale != 1.) {
#line 487 "dtrevc.f"
			    dscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 487 "dtrevc.f"
			}
#line 489 "dtrevc.f"
			work[j - 1 + *n] = x[0];
#line 490 "dtrevc.f"
			work[j + *n] = x[1];

/*                    Update right-hand side */

#line 494 "dtrevc.f"
			i__1 = j - 2;
#line 494 "dtrevc.f"
			d__1 = -x[0];
#line 494 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[*n + 1], &c__1);
#line 496 "dtrevc.f"
			i__1 = j - 2;
#line 496 "dtrevc.f"
			d__1 = -x[1];
#line 496 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 498 "dtrevc.f"
		    }
#line 499 "dtrevc.f"
L60:
#line 499 "dtrevc.f"
		    ;
#line 499 "dtrevc.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 503 "dtrevc.f"
		if (! over) {
#line 504 "dtrevc.f"
		    dcopy_(&ki, &work[*n + 1], &c__1, &vr[is * vr_dim1 + 1], &
			    c__1);

#line 506 "dtrevc.f"
		    ii = idamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 507 "dtrevc.f"
		    remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs(d__1));
#line 508 "dtrevc.f"
		    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 510 "dtrevc.f"
		    i__1 = *n;
#line 510 "dtrevc.f"
		    for (k = ki + 1; k <= i__1; ++k) {
#line 511 "dtrevc.f"
			vr[k + is * vr_dim1] = 0.;
#line 512 "dtrevc.f"
/* L70: */
#line 512 "dtrevc.f"
		    }
#line 513 "dtrevc.f"
		} else {
#line 514 "dtrevc.f"
		    if (ki > 1) {
#line 514 "dtrevc.f"
			i__1 = ki - 1;
#line 514 "dtrevc.f"
			dgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[*n + 1], &c__1, &work[ki + *n], &vr[ki * 
				vr_dim1 + 1], &c__1, (ftnlen)1);
#line 514 "dtrevc.f"
		    }

#line 519 "dtrevc.f"
		    ii = idamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 520 "dtrevc.f"
		    remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs(d__1));
#line 521 "dtrevc.f"
		    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
#line 522 "dtrevc.f"
		}

#line 524 "dtrevc.f"
	    } else {

/*              Complex right eigenvector. */

/*              Initial solve */
/*                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0. */
/*                [ (T(KI,KI-1)   T(KI,KI)   )               ] */

#line 532 "dtrevc.f"
		if ((d__1 = t[ki - 1 + ki * t_dim1], abs(d__1)) >= (d__2 = t[
			ki + (ki - 1) * t_dim1], abs(d__2))) {
#line 533 "dtrevc.f"
		    work[ki - 1 + *n] = 1.;
#line 534 "dtrevc.f"
		    work[ki + n2] = wi / t[ki - 1 + ki * t_dim1];
#line 535 "dtrevc.f"
		} else {
#line 536 "dtrevc.f"
		    work[ki - 1 + *n] = -wi / t[ki + (ki - 1) * t_dim1];
#line 537 "dtrevc.f"
		    work[ki + n2] = 1.;
#line 538 "dtrevc.f"
		}
#line 539 "dtrevc.f"
		work[ki + *n] = 0.;
#line 540 "dtrevc.f"
		work[ki - 1 + n2] = 0.;

/*              Form right-hand side */

#line 544 "dtrevc.f"
		i__1 = ki - 2;
#line 544 "dtrevc.f"
		for (k = 1; k <= i__1; ++k) {
#line 545 "dtrevc.f"
		    work[k + *n] = -work[ki - 1 + *n] * t[k + (ki - 1) * 
			    t_dim1];
#line 546 "dtrevc.f"
		    work[k + n2] = -work[ki + n2] * t[k + ki * t_dim1];
#line 547 "dtrevc.f"
/* L80: */
#line 547 "dtrevc.f"
		}

/*              Solve upper quasi-triangular system: */
/*              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */

#line 552 "dtrevc.f"
		jnxt = ki - 2;
#line 553 "dtrevc.f"
		for (j = ki - 2; j >= 1; --j) {
#line 554 "dtrevc.f"
		    if (j > jnxt) {
#line 554 "dtrevc.f"
			goto L90;
#line 554 "dtrevc.f"
		    }
#line 556 "dtrevc.f"
		    j1 = j;
#line 557 "dtrevc.f"
		    j2 = j;
#line 558 "dtrevc.f"
		    jnxt = j - 1;
#line 559 "dtrevc.f"
		    if (j > 1) {
#line 560 "dtrevc.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 561 "dtrevc.f"
			    j1 = j - 1;
#line 562 "dtrevc.f"
			    jnxt = j - 2;
#line 563 "dtrevc.f"
			}
#line 564 "dtrevc.f"
		    }

#line 566 "dtrevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 570 "dtrevc.f"
			dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &wi, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale X(1,1) and X(1,2) to avoid overflow when */
/*                    updating the right-hand side. */

#line 577 "dtrevc.f"
			if (xnorm > 1.) {
#line 578 "dtrevc.f"
			    if (work[j] > bignum / xnorm) {
#line 579 "dtrevc.f"
				x[0] /= xnorm;
#line 580 "dtrevc.f"
				x[2] /= xnorm;
#line 581 "dtrevc.f"
				scale /= xnorm;
#line 582 "dtrevc.f"
			    }
#line 583 "dtrevc.f"
			}

/*                    Scale if necessary */

#line 587 "dtrevc.f"
			if (scale != 1.) {
#line 588 "dtrevc.f"
			    dscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 589 "dtrevc.f"
			    dscal_(&ki, &scale, &work[n2 + 1], &c__1);
#line 590 "dtrevc.f"
			}
#line 591 "dtrevc.f"
			work[j + *n] = x[0];
#line 592 "dtrevc.f"
			work[j + n2] = x[2];

/*                    Update the right-hand side */

#line 596 "dtrevc.f"
			i__1 = j - 1;
#line 596 "dtrevc.f"
			d__1 = -x[0];
#line 596 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 598 "dtrevc.f"
			i__1 = j - 1;
#line 598 "dtrevc.f"
			d__1 = -x[2];
#line 598 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				n2 + 1], &c__1);

#line 601 "dtrevc.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 605 "dtrevc.f"
			dlaln2_(&c_false, &c__2, &c__2, &smin, &c_b22, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b22, &c_b22, &
				work[j - 1 + *n], n, &wr, &wi, x, &c__2, &
				scale, &xnorm, &ierr);

/*                    Scale X to avoid overflow when updating */
/*                    the right-hand side. */

#line 613 "dtrevc.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 614 "dtrevc.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 614 "dtrevc.f"
			    beta = max(d__1,d__2);
#line 615 "dtrevc.f"
			    if (beta > bignum / xnorm) {
#line 616 "dtrevc.f"
				rec = 1. / xnorm;
#line 617 "dtrevc.f"
				x[0] *= rec;
#line 618 "dtrevc.f"
				x[2] *= rec;
#line 619 "dtrevc.f"
				x[1] *= rec;
#line 620 "dtrevc.f"
				x[3] *= rec;
#line 621 "dtrevc.f"
				scale *= rec;
#line 622 "dtrevc.f"
			    }
#line 623 "dtrevc.f"
			}

/*                    Scale if necessary */

#line 627 "dtrevc.f"
			if (scale != 1.) {
#line 628 "dtrevc.f"
			    dscal_(&ki, &scale, &work[*n + 1], &c__1);
#line 629 "dtrevc.f"
			    dscal_(&ki, &scale, &work[n2 + 1], &c__1);
#line 630 "dtrevc.f"
			}
#line 631 "dtrevc.f"
			work[j - 1 + *n] = x[0];
#line 632 "dtrevc.f"
			work[j + *n] = x[1];
#line 633 "dtrevc.f"
			work[j - 1 + n2] = x[2];
#line 634 "dtrevc.f"
			work[j + n2] = x[3];

/*                    Update the right-hand side */

#line 638 "dtrevc.f"
			i__1 = j - 2;
#line 638 "dtrevc.f"
			d__1 = -x[0];
#line 638 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[*n + 1], &c__1);
#line 640 "dtrevc.f"
			i__1 = j - 2;
#line 640 "dtrevc.f"
			d__1 = -x[1];
#line 640 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				*n + 1], &c__1);
#line 642 "dtrevc.f"
			i__1 = j - 2;
#line 642 "dtrevc.f"
			d__1 = -x[2];
#line 642 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[n2 + 1], &c__1);
#line 644 "dtrevc.f"
			i__1 = j - 2;
#line 644 "dtrevc.f"
			d__1 = -x[3];
#line 644 "dtrevc.f"
			daxpy_(&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				n2 + 1], &c__1);
#line 646 "dtrevc.f"
		    }
#line 647 "dtrevc.f"
L90:
#line 647 "dtrevc.f"
		    ;
#line 647 "dtrevc.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 651 "dtrevc.f"
		if (! over) {
#line 652 "dtrevc.f"
		    dcopy_(&ki, &work[*n + 1], &c__1, &vr[(is - 1) * vr_dim1 
			    + 1], &c__1);
#line 653 "dtrevc.f"
		    dcopy_(&ki, &work[n2 + 1], &c__1, &vr[is * vr_dim1 + 1], &
			    c__1);

#line 655 "dtrevc.f"
		    emax = 0.;
#line 656 "dtrevc.f"
		    i__1 = ki;
#line 656 "dtrevc.f"
		    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 657 "dtrevc.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + is * vr_dim1], 
				abs(d__2));
#line 657 "dtrevc.f"
			emax = max(d__3,d__4);
#line 659 "dtrevc.f"
/* L100: */
#line 659 "dtrevc.f"
		    }

#line 661 "dtrevc.f"
		    remax = 1. / emax;
#line 662 "dtrevc.f"
		    dscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
#line 663 "dtrevc.f"
		    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 665 "dtrevc.f"
		    i__1 = *n;
#line 665 "dtrevc.f"
		    for (k = ki + 1; k <= i__1; ++k) {
#line 666 "dtrevc.f"
			vr[k + (is - 1) * vr_dim1] = 0.;
#line 667 "dtrevc.f"
			vr[k + is * vr_dim1] = 0.;
#line 668 "dtrevc.f"
/* L110: */
#line 668 "dtrevc.f"
		    }

#line 670 "dtrevc.f"
		} else {

#line 672 "dtrevc.f"
		    if (ki > 2) {
#line 673 "dtrevc.f"
			i__1 = ki - 2;
#line 673 "dtrevc.f"
			dgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[*n + 1], &c__1, &work[ki - 1 + *n], &vr[(
				ki - 1) * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 676 "dtrevc.f"
			i__1 = ki - 2;
#line 676 "dtrevc.f"
			dgemv_("N", n, &i__1, &c_b22, &vr[vr_offset], ldvr, &
				work[n2 + 1], &c__1, &work[ki + n2], &vr[ki * 
				vr_dim1 + 1], &c__1, (ftnlen)1);
#line 679 "dtrevc.f"
		    } else {
#line 680 "dtrevc.f"
			dscal_(n, &work[ki - 1 + *n], &vr[(ki - 1) * vr_dim1 
				+ 1], &c__1);
#line 681 "dtrevc.f"
			dscal_(n, &work[ki + n2], &vr[ki * vr_dim1 + 1], &
				c__1);
#line 682 "dtrevc.f"
		    }

#line 684 "dtrevc.f"
		    emax = 0.;
#line 685 "dtrevc.f"
		    i__1 = *n;
#line 685 "dtrevc.f"
		    for (k = 1; k <= i__1; ++k) {
/* Computing MAX */
#line 686 "dtrevc.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + ki * vr_dim1], 
				abs(d__2));
#line 686 "dtrevc.f"
			emax = max(d__3,d__4);
#line 688 "dtrevc.f"
/* L120: */
#line 688 "dtrevc.f"
		    }
#line 689 "dtrevc.f"
		    remax = 1. / emax;
#line 690 "dtrevc.f"
		    dscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
#line 691 "dtrevc.f"
		    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
#line 692 "dtrevc.f"
		}
#line 693 "dtrevc.f"
	    }

#line 695 "dtrevc.f"
	    --is;
#line 696 "dtrevc.f"
	    if (ip != 0) {
#line 696 "dtrevc.f"
		--is;
#line 696 "dtrevc.f"
	    }
#line 698 "dtrevc.f"
L130:
#line 699 "dtrevc.f"
	    if (ip == 1) {
#line 699 "dtrevc.f"
		ip = 0;
#line 699 "dtrevc.f"
	    }
#line 701 "dtrevc.f"
	    if (ip == -1) {
#line 701 "dtrevc.f"
		ip = 1;
#line 701 "dtrevc.f"
	    }
#line 703 "dtrevc.f"
/* L140: */
#line 703 "dtrevc.f"
	}
#line 704 "dtrevc.f"
    }

#line 706 "dtrevc.f"
    if (leftv) {

/*        Compute left eigenvectors. */

#line 710 "dtrevc.f"
	ip = 0;
#line 711 "dtrevc.f"
	is = 1;
#line 712 "dtrevc.f"
	i__1 = *n;
#line 712 "dtrevc.f"
	for (ki = 1; ki <= i__1; ++ki) {

#line 714 "dtrevc.f"
	    if (ip == -1) {
#line 714 "dtrevc.f"
		goto L250;
#line 714 "dtrevc.f"
	    }
#line 716 "dtrevc.f"
	    if (ki == *n) {
#line 716 "dtrevc.f"
		goto L150;
#line 716 "dtrevc.f"
	    }
#line 718 "dtrevc.f"
	    if (t[ki + 1 + ki * t_dim1] == 0.) {
#line 718 "dtrevc.f"
		goto L150;
#line 718 "dtrevc.f"
	    }
#line 720 "dtrevc.f"
	    ip = 1;

#line 722 "dtrevc.f"
L150:
#line 723 "dtrevc.f"
	    if (somev) {
#line 724 "dtrevc.f"
		if (! select[ki]) {
#line 724 "dtrevc.f"
		    goto L250;
#line 724 "dtrevc.f"
		}
#line 726 "dtrevc.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 730 "dtrevc.f"
	    wr = t[ki + ki * t_dim1];
#line 731 "dtrevc.f"
	    wi = 0.;
#line 732 "dtrevc.f"
	    if (ip != 0) {
#line 732 "dtrevc.f"
		wi = sqrt((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki + 1 + ki * t_dim1], abs(d__2)));
#line 732 "dtrevc.f"
	    }
/* Computing MAX */
#line 735 "dtrevc.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 735 "dtrevc.f"
	    smin = max(d__1,smlnum);

#line 737 "dtrevc.f"
	    if (ip == 0) {

/*              Real left eigenvector. */

#line 741 "dtrevc.f"
		work[ki + *n] = 1.;

/*              Form right-hand side */

#line 745 "dtrevc.f"
		i__2 = *n;
#line 745 "dtrevc.f"
		for (k = ki + 1; k <= i__2; ++k) {
#line 746 "dtrevc.f"
		    work[k + *n] = -t[ki + k * t_dim1];
#line 747 "dtrevc.f"
/* L160: */
#line 747 "dtrevc.f"
		}

/*              Solve the quasi-triangular system: */
/*                 (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK */

#line 752 "dtrevc.f"
		vmax = 1.;
#line 753 "dtrevc.f"
		vcrit = bignum;

#line 755 "dtrevc.f"
		jnxt = ki + 1;
#line 756 "dtrevc.f"
		i__2 = *n;
#line 756 "dtrevc.f"
		for (j = ki + 1; j <= i__2; ++j) {
#line 757 "dtrevc.f"
		    if (j < jnxt) {
#line 757 "dtrevc.f"
			goto L170;
#line 757 "dtrevc.f"
		    }
#line 759 "dtrevc.f"
		    j1 = j;
#line 760 "dtrevc.f"
		    j2 = j;
#line 761 "dtrevc.f"
		    jnxt = j + 1;
#line 762 "dtrevc.f"
		    if (j < *n) {
#line 763 "dtrevc.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 764 "dtrevc.f"
			    j2 = j + 1;
#line 765 "dtrevc.f"
			    jnxt = j + 2;
#line 766 "dtrevc.f"
			}
#line 767 "dtrevc.f"
		    }

#line 769 "dtrevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

#line 776 "dtrevc.f"
			if (work[j] > vcrit) {
#line 777 "dtrevc.f"
			    rec = 1. / vmax;
#line 778 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 778 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 779 "dtrevc.f"
			    vmax = 1.;
#line 780 "dtrevc.f"
			    vcrit = bignum;
#line 781 "dtrevc.f"
			}

#line 783 "dtrevc.f"
			i__3 = j - ki - 1;
#line 783 "dtrevc.f"
			work[j + *n] -= ddot_(&i__3, &t[ki + 1 + j * t_dim1], 
				&c__1, &work[ki + 1 + *n], &c__1);

/*                    Solve (T(J,J)-WR)**T*X = WORK */

#line 789 "dtrevc.f"
			dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale if necessary */

#line 795 "dtrevc.f"
			if (scale != 1.) {
#line 795 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 795 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 795 "dtrevc.f"
			}
#line 797 "dtrevc.f"
			work[j + *n] = x[0];
/* Computing MAX */
#line 798 "dtrevc.f"
			d__2 = (d__1 = work[j + *n], abs(d__1));
#line 798 "dtrevc.f"
			vmax = max(d__2,vmax);
#line 799 "dtrevc.f"
			vcrit = bignum / vmax;

#line 801 "dtrevc.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

/* Computing MAX */
#line 808 "dtrevc.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 808 "dtrevc.f"
			beta = max(d__1,d__2);
#line 809 "dtrevc.f"
			if (beta > vcrit) {
#line 810 "dtrevc.f"
			    rec = 1. / vmax;
#line 811 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 811 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 812 "dtrevc.f"
			    vmax = 1.;
#line 813 "dtrevc.f"
			    vcrit = bignum;
#line 814 "dtrevc.f"
			}

#line 816 "dtrevc.f"
			i__3 = j - ki - 1;
#line 816 "dtrevc.f"
			work[j + *n] -= ddot_(&i__3, &t[ki + 1 + j * t_dim1], 
				&c__1, &work[ki + 1 + *n], &c__1);

#line 820 "dtrevc.f"
			i__3 = j - ki - 1;
#line 820 "dtrevc.f"
			work[j + 1 + *n] -= ddot_(&i__3, &t[ki + 1 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 1 + *n], &c__1);

/*                    Solve */
/*                      [T(J,J)-WR   T(J,J+1)     ]**T * X = SCALE*( WORK1 ) */
/*                      [T(J+1,J)    T(J+1,J+1)-WR]                ( WORK2 ) */

#line 828 "dtrevc.f"
			dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &c_b25, x, &c__2, &scale, &xnorm, 
				&ierr);

/*                    Scale if necessary */

#line 834 "dtrevc.f"
			if (scale != 1.) {
#line 834 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 834 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 834 "dtrevc.f"
			}
#line 836 "dtrevc.f"
			work[j + *n] = x[0];
#line 837 "dtrevc.f"
			work[j + 1 + *n] = x[1];

/* Computing MAX */
#line 839 "dtrevc.f"
			d__3 = (d__1 = work[j + *n], abs(d__1)), d__4 = (d__2 
				= work[j + 1 + *n], abs(d__2)), d__3 = max(
				d__3,d__4);
#line 839 "dtrevc.f"
			vmax = max(d__3,vmax);
#line 841 "dtrevc.f"
			vcrit = bignum / vmax;

#line 843 "dtrevc.f"
		    }
#line 844 "dtrevc.f"
L170:
#line 844 "dtrevc.f"
		    ;
#line 844 "dtrevc.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 848 "dtrevc.f"
		if (! over) {
#line 849 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 849 "dtrevc.f"
		    dcopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);

#line 851 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 851 "dtrevc.f"
		    ii = idamax_(&i__2, &vl[ki + is * vl_dim1], &c__1) + ki - 
			    1;
#line 852 "dtrevc.f"
		    remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs(d__1));
#line 853 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 853 "dtrevc.f"
		    dscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 855 "dtrevc.f"
		    i__2 = ki - 1;
#line 855 "dtrevc.f"
		    for (k = 1; k <= i__2; ++k) {
#line 856 "dtrevc.f"
			vl[k + is * vl_dim1] = 0.;
#line 857 "dtrevc.f"
/* L180: */
#line 857 "dtrevc.f"
		    }

#line 859 "dtrevc.f"
		} else {

#line 861 "dtrevc.f"
		    if (ki < *n) {
#line 861 "dtrevc.f"
			i__2 = *n - ki;
#line 861 "dtrevc.f"
			dgemv_("N", n, &i__2, &c_b22, &vl[(ki + 1) * vl_dim1 
				+ 1], ldvl, &work[ki + 1 + *n], &c__1, &work[
				ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (
				ftnlen)1);
#line 861 "dtrevc.f"
		    }

#line 866 "dtrevc.f"
		    ii = idamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 867 "dtrevc.f"
		    remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs(d__1));
#line 868 "dtrevc.f"
		    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

#line 870 "dtrevc.f"
		}

#line 872 "dtrevc.f"
	    } else {

/*              Complex left eigenvector. */

/*               Initial solve: */
/*                 ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0. */
/*                 ((T(KI+1,KI) T(KI+1,KI+1))                ) */

#line 880 "dtrevc.f"
		if ((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1)) >= (d__2 = 
			t[ki + 1 + ki * t_dim1], abs(d__2))) {
#line 881 "dtrevc.f"
		    work[ki + *n] = wi / t[ki + (ki + 1) * t_dim1];
#line 882 "dtrevc.f"
		    work[ki + 1 + n2] = 1.;
#line 883 "dtrevc.f"
		} else {
#line 884 "dtrevc.f"
		    work[ki + *n] = 1.;
#line 885 "dtrevc.f"
		    work[ki + 1 + n2] = -wi / t[ki + 1 + ki * t_dim1];
#line 886 "dtrevc.f"
		}
#line 887 "dtrevc.f"
		work[ki + 1 + *n] = 0.;
#line 888 "dtrevc.f"
		work[ki + n2] = 0.;

/*              Form right-hand side */

#line 892 "dtrevc.f"
		i__2 = *n;
#line 892 "dtrevc.f"
		for (k = ki + 2; k <= i__2; ++k) {
#line 893 "dtrevc.f"
		    work[k + *n] = -work[ki + *n] * t[ki + k * t_dim1];
#line 894 "dtrevc.f"
		    work[k + n2] = -work[ki + 1 + n2] * t[ki + 1 + k * t_dim1]
			    ;
#line 895 "dtrevc.f"
/* L190: */
#line 895 "dtrevc.f"
		}

/*              Solve complex quasi-triangular system: */
/*              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */

#line 900 "dtrevc.f"
		vmax = 1.;
#line 901 "dtrevc.f"
		vcrit = bignum;

#line 903 "dtrevc.f"
		jnxt = ki + 2;
#line 904 "dtrevc.f"
		i__2 = *n;
#line 904 "dtrevc.f"
		for (j = ki + 2; j <= i__2; ++j) {
#line 905 "dtrevc.f"
		    if (j < jnxt) {
#line 905 "dtrevc.f"
			goto L200;
#line 905 "dtrevc.f"
		    }
#line 907 "dtrevc.f"
		    j1 = j;
#line 908 "dtrevc.f"
		    j2 = j;
#line 909 "dtrevc.f"
		    jnxt = j + 1;
#line 910 "dtrevc.f"
		    if (j < *n) {
#line 911 "dtrevc.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 912 "dtrevc.f"
			    j2 = j + 1;
#line 913 "dtrevc.f"
			    jnxt = j + 2;
#line 914 "dtrevc.f"
			}
#line 915 "dtrevc.f"
		    }

#line 917 "dtrevc.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when */
/*                    forming the right-hand side elements. */

#line 924 "dtrevc.f"
			if (work[j] > vcrit) {
#line 925 "dtrevc.f"
			    rec = 1. / vmax;
#line 926 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 926 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 927 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 927 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + n2], &c__1);
#line 928 "dtrevc.f"
			    vmax = 1.;
#line 929 "dtrevc.f"
			    vcrit = bignum;
#line 930 "dtrevc.f"
			}

#line 932 "dtrevc.f"
			i__3 = j - ki - 2;
#line 932 "dtrevc.f"
			work[j + *n] -= ddot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + *n], &c__1);
#line 935 "dtrevc.f"
			i__3 = j - ki - 2;
#line 935 "dtrevc.f"
			work[j + n2] -= ddot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + n2], &c__1);

/*                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */

#line 941 "dtrevc.f"
			d__1 = -wi;
#line 941 "dtrevc.f"
			dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &d__1, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale if necessary */

#line 947 "dtrevc.f"
			if (scale != 1.) {
#line 948 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 948 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 949 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 949 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + n2], &c__1);
#line 950 "dtrevc.f"
			}
#line 951 "dtrevc.f"
			work[j + *n] = x[0];
#line 952 "dtrevc.f"
			work[j + n2] = x[2];
/* Computing MAX */
#line 953 "dtrevc.f"
			d__3 = (d__1 = work[j + *n], abs(d__1)), d__4 = (d__2 
				= work[j + n2], abs(d__2)), d__3 = max(d__3,
				d__4);
#line 953 "dtrevc.f"
			vmax = max(d__3,vmax);
#line 955 "dtrevc.f"
			vcrit = bignum / vmax;

#line 957 "dtrevc.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side elements. */

/* Computing MAX */
#line 964 "dtrevc.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 964 "dtrevc.f"
			beta = max(d__1,d__2);
#line 965 "dtrevc.f"
			if (beta > vcrit) {
#line 966 "dtrevc.f"
			    rec = 1. / vmax;
#line 967 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 967 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + *n], &c__1);
#line 968 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 968 "dtrevc.f"
			    dscal_(&i__3, &rec, &work[ki + n2], &c__1);
#line 969 "dtrevc.f"
			    vmax = 1.;
#line 970 "dtrevc.f"
			    vcrit = bignum;
#line 971 "dtrevc.f"
			}

#line 973 "dtrevc.f"
			i__3 = j - ki - 2;
#line 973 "dtrevc.f"
			work[j + *n] -= ddot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + *n], &c__1);

#line 977 "dtrevc.f"
			i__3 = j - ki - 2;
#line 977 "dtrevc.f"
			work[j + n2] -= ddot_(&i__3, &t[ki + 2 + j * t_dim1], 
				&c__1, &work[ki + 2 + n2], &c__1);

#line 981 "dtrevc.f"
			i__3 = j - ki - 2;
#line 981 "dtrevc.f"
			work[j + 1 + *n] -= ddot_(&i__3, &t[ki + 2 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 2 + *n], &c__1);

#line 985 "dtrevc.f"
			i__3 = j - ki - 2;
#line 985 "dtrevc.f"
			work[j + 1 + n2] -= ddot_(&i__3, &t[ki + 2 + (j + 1) *
				 t_dim1], &c__1, &work[ki + 2 + n2], &c__1);

/*                    Solve 2-by-2 complex linear equation */
/*                      ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B */
/*                      ([T(j+1,j) T(j+1,j+1)]               ) */

#line 993 "dtrevc.f"
			d__1 = -wi;
#line 993 "dtrevc.f"
			dlaln2_(&c_true, &c__2, &c__2, &smin, &c_b22, &t[j + 
				j * t_dim1], ldt, &c_b22, &c_b22, &work[j + *
				n], n, &wr, &d__1, x, &c__2, &scale, &xnorm, &
				ierr);

/*                    Scale if necessary */

#line 999 "dtrevc.f"
			if (scale != 1.) {
#line 1000 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 1000 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + *n], &c__1);
#line 1001 "dtrevc.f"
			    i__3 = *n - ki + 1;
#line 1001 "dtrevc.f"
			    dscal_(&i__3, &scale, &work[ki + n2], &c__1);
#line 1002 "dtrevc.f"
			}
#line 1003 "dtrevc.f"
			work[j + *n] = x[0];
#line 1004 "dtrevc.f"
			work[j + n2] = x[2];
#line 1005 "dtrevc.f"
			work[j + 1 + *n] = x[1];
#line 1006 "dtrevc.f"
			work[j + 1 + n2] = x[3];
/* Computing MAX */
#line 1007 "dtrevc.f"
			d__1 = abs(x[0]), d__2 = abs(x[2]), d__1 = max(d__1,
				d__2), d__2 = abs(x[1]), d__1 = max(d__1,d__2)
				, d__2 = abs(x[3]), d__1 = max(d__1,d__2);
#line 1007 "dtrevc.f"
			vmax = max(d__1,vmax);
#line 1009 "dtrevc.f"
			vcrit = bignum / vmax;

#line 1011 "dtrevc.f"
		    }
#line 1012 "dtrevc.f"
L200:
#line 1012 "dtrevc.f"
		    ;
#line 1012 "dtrevc.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 1016 "dtrevc.f"
		if (! over) {
#line 1017 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 1017 "dtrevc.f"
		    dcopy_(&i__2, &work[ki + *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);
#line 1018 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 1018 "dtrevc.f"
		    dcopy_(&i__2, &work[ki + n2], &c__1, &vl[ki + (is + 1) * 
			    vl_dim1], &c__1);

#line 1021 "dtrevc.f"
		    emax = 0.;
#line 1022 "dtrevc.f"
		    i__2 = *n;
#line 1022 "dtrevc.f"
		    for (k = ki; k <= i__2; ++k) {
/* Computing MAX */
#line 1023 "dtrevc.f"
			d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (is + 1) * vl_dim1], 
				abs(d__2));
#line 1023 "dtrevc.f"
			emax = max(d__3,d__4);
#line 1025 "dtrevc.f"
/* L220: */
#line 1025 "dtrevc.f"
		    }
#line 1026 "dtrevc.f"
		    remax = 1. / emax;
#line 1027 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 1027 "dtrevc.f"
		    dscal_(&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);
#line 1028 "dtrevc.f"
		    i__2 = *n - ki + 1;
#line 1028 "dtrevc.f"
		    dscal_(&i__2, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1)
			    ;

#line 1030 "dtrevc.f"
		    i__2 = ki - 1;
#line 1030 "dtrevc.f"
		    for (k = 1; k <= i__2; ++k) {
#line 1031 "dtrevc.f"
			vl[k + is * vl_dim1] = 0.;
#line 1032 "dtrevc.f"
			vl[k + (is + 1) * vl_dim1] = 0.;
#line 1033 "dtrevc.f"
/* L230: */
#line 1033 "dtrevc.f"
		    }
#line 1034 "dtrevc.f"
		} else {
#line 1035 "dtrevc.f"
		    if (ki < *n - 1) {
#line 1036 "dtrevc.f"
			i__2 = *n - ki - 1;
#line 1036 "dtrevc.f"
			dgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + *n], &c__1, &work[
				ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (
				ftnlen)1);
#line 1039 "dtrevc.f"
			i__2 = *n - ki - 1;
#line 1039 "dtrevc.f"
			dgemv_("N", n, &i__2, &c_b22, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + n2], &c__1, &work[
				ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1042 "dtrevc.f"
		    } else {
#line 1043 "dtrevc.f"
			dscal_(n, &work[ki + *n], &vl[ki * vl_dim1 + 1], &
				c__1);
#line 1044 "dtrevc.f"
			dscal_(n, &work[ki + 1 + n2], &vl[(ki + 1) * vl_dim1 
				+ 1], &c__1);
#line 1045 "dtrevc.f"
		    }

#line 1047 "dtrevc.f"
		    emax = 0.;
#line 1048 "dtrevc.f"
		    i__2 = *n;
#line 1048 "dtrevc.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 1049 "dtrevc.f"
			d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (ki + 1) * vl_dim1], 
				abs(d__2));
#line 1049 "dtrevc.f"
			emax = max(d__3,d__4);
#line 1051 "dtrevc.f"
/* L240: */
#line 1051 "dtrevc.f"
		    }
#line 1052 "dtrevc.f"
		    remax = 1. / emax;
#line 1053 "dtrevc.f"
		    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
#line 1054 "dtrevc.f"
		    dscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);

#line 1056 "dtrevc.f"
		}

#line 1058 "dtrevc.f"
	    }

#line 1060 "dtrevc.f"
	    ++is;
#line 1061 "dtrevc.f"
	    if (ip != 0) {
#line 1061 "dtrevc.f"
		++is;
#line 1061 "dtrevc.f"
	    }
#line 1063 "dtrevc.f"
L250:
#line 1064 "dtrevc.f"
	    if (ip == -1) {
#line 1064 "dtrevc.f"
		ip = 0;
#line 1064 "dtrevc.f"
	    }
#line 1066 "dtrevc.f"
	    if (ip == 1) {
#line 1066 "dtrevc.f"
		ip = -1;
#line 1066 "dtrevc.f"
	    }

#line 1069 "dtrevc.f"
/* L260: */
#line 1069 "dtrevc.f"
	}

#line 1071 "dtrevc.f"
    }

#line 1073 "dtrevc.f"
    return 0;

/*     End of DTREVC */

} /* dtrevc_ */

