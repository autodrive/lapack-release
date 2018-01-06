#line 1 "dtrevc3.f"
/* dtrevc3.f -- translated by f2c (version 20100827).
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

#line 1 "dtrevc3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static logical c_false = FALSE_;
static doublereal c_b29 = 1.;
static logical c_true = TRUE_;

/* > \brief \b DTREVC3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTREVC3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrevc3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrevc3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrevc3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, */
/*                           VR, LDVR, MM, M, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N */
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
/* > DTREVC3 computes some or all of the right and/or left eigenvectors of */
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
/* > input matrix. If Q is the orthogonal factor that reduces a matrix */
/* > A to Schur form T, then Q*X and Q*Y are the matrices of right and */
/* > left eigenvectors of A. */
/* > */
/* > This uses a Level 3 BLAS version of the back transformation. */
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
/* >          The leading dimension of the array VL. */
/* >          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N. */
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
/* >          The leading dimension of the array VR. */
/* >          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of array WORK. LWORK >= max(1,3*N). */
/* >          For optimum performance, LWORK >= N + 2*N*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/*  @precisions fortran d -> s */

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
/* Subroutine */ int dtrevc3_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, 
	doublereal *work, integer *lwork, integer *info, ftnlen side_len, 
	ftnlen howmny_len)
{
    /* System generated locals */
    address a__1[2];
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1[2],
	     i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal x[4]	/* was [2][2] */;
    static integer j1, j2, iscomplex[128], nb, ii, ki, ip, is, iv;
    static doublereal wi, wr;
    static integer ki2;
    static doublereal rec, ulp, beta, emax;
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
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
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
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal bignum;
    static logical rightv;
    static integer maxwrk;
    static doublereal smlnum;
    static logical lquery;


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

#line 296 "dtrevc3.f"
    /* Parameter adjustments */
#line 296 "dtrevc3.f"
    --select;
#line 296 "dtrevc3.f"
    t_dim1 = *ldt;
#line 296 "dtrevc3.f"
    t_offset = 1 + t_dim1;
#line 296 "dtrevc3.f"
    t -= t_offset;
#line 296 "dtrevc3.f"
    vl_dim1 = *ldvl;
#line 296 "dtrevc3.f"
    vl_offset = 1 + vl_dim1;
#line 296 "dtrevc3.f"
    vl -= vl_offset;
#line 296 "dtrevc3.f"
    vr_dim1 = *ldvr;
#line 296 "dtrevc3.f"
    vr_offset = 1 + vr_dim1;
#line 296 "dtrevc3.f"
    vr -= vr_offset;
#line 296 "dtrevc3.f"
    --work;
#line 296 "dtrevc3.f"

#line 296 "dtrevc3.f"
    /* Function Body */
#line 296 "dtrevc3.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 297 "dtrevc3.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 298 "dtrevc3.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 300 "dtrevc3.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 301 "dtrevc3.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 302 "dtrevc3.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 304 "dtrevc3.f"
    *info = 0;
/* Writing concatenation */
#line 305 "dtrevc3.f"
    i__1[0] = 1, a__1[0] = side;
#line 305 "dtrevc3.f"
    i__1[1] = 1, a__1[1] = howmny;
#line 305 "dtrevc3.f"
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 305 "dtrevc3.f"
    nb = ilaenv_(&c__1, "DTREVC", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
#line 306 "dtrevc3.f"
    maxwrk = *n + (*n << 1) * nb;
#line 307 "dtrevc3.f"
    work[1] = (doublereal) maxwrk;
#line 308 "dtrevc3.f"
    lquery = *lwork == -1;
#line 309 "dtrevc3.f"
    if (! rightv && ! leftv) {
#line 310 "dtrevc3.f"
	*info = -1;
#line 311 "dtrevc3.f"
    } else if (! allv && ! over && ! somev) {
#line 312 "dtrevc3.f"
	*info = -2;
#line 313 "dtrevc3.f"
    } else if (*n < 0) {
#line 314 "dtrevc3.f"
	*info = -4;
#line 315 "dtrevc3.f"
    } else if (*ldt < max(1,*n)) {
#line 316 "dtrevc3.f"
	*info = -6;
#line 317 "dtrevc3.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 318 "dtrevc3.f"
	*info = -8;
#line 319 "dtrevc3.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 320 "dtrevc3.f"
	*info = -10;
#line 321 "dtrevc3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 321 "dtrevc3.f"
	i__2 = 1, i__3 = *n * 3;
#line 321 "dtrevc3.f"
	if (*lwork < max(i__2,i__3) && ! lquery) {
#line 322 "dtrevc3.f"
	    *info = -14;
#line 323 "dtrevc3.f"
	} else {

/*        Set M to the number of columns required to store the selected */
/*        eigenvectors, standardize the array SELECT if necessary, and */
/*        test MM. */

#line 329 "dtrevc3.f"
	    if (somev) {
#line 330 "dtrevc3.f"
		*m = 0;
#line 331 "dtrevc3.f"
		pair = FALSE_;
#line 332 "dtrevc3.f"
		i__2 = *n;
#line 332 "dtrevc3.f"
		for (j = 1; j <= i__2; ++j) {
#line 333 "dtrevc3.f"
		    if (pair) {
#line 334 "dtrevc3.f"
			pair = FALSE_;
#line 335 "dtrevc3.f"
			select[j] = FALSE_;
#line 336 "dtrevc3.f"
		    } else {
#line 337 "dtrevc3.f"
			if (j < *n) {
#line 338 "dtrevc3.f"
			    if (t[j + 1 + j * t_dim1] == 0.) {
#line 339 "dtrevc3.f"
				if (select[j]) {
#line 339 "dtrevc3.f"
				    ++(*m);
#line 339 "dtrevc3.f"
				}
#line 341 "dtrevc3.f"
			    } else {
#line 342 "dtrevc3.f"
				pair = TRUE_;
#line 343 "dtrevc3.f"
				if (select[j] || select[j + 1]) {
#line 344 "dtrevc3.f"
				    select[j] = TRUE_;
#line 345 "dtrevc3.f"
				    *m += 2;
#line 346 "dtrevc3.f"
				}
#line 347 "dtrevc3.f"
			    }
#line 348 "dtrevc3.f"
			} else {
#line 349 "dtrevc3.f"
			    if (select[*n]) {
#line 349 "dtrevc3.f"
				++(*m);
#line 349 "dtrevc3.f"
			    }
#line 351 "dtrevc3.f"
			}
#line 352 "dtrevc3.f"
		    }
#line 353 "dtrevc3.f"
/* L10: */
#line 353 "dtrevc3.f"
		}
#line 354 "dtrevc3.f"
	    } else {
#line 355 "dtrevc3.f"
		*m = *n;
#line 356 "dtrevc3.f"
	    }

#line 358 "dtrevc3.f"
	    if (*mm < *m) {
#line 359 "dtrevc3.f"
		*info = -11;
#line 360 "dtrevc3.f"
	    }
#line 361 "dtrevc3.f"
	}
#line 361 "dtrevc3.f"
    }
#line 362 "dtrevc3.f"
    if (*info != 0) {
#line 363 "dtrevc3.f"
	i__2 = -(*info);
#line 363 "dtrevc3.f"
	xerbla_("DTREVC3", &i__2, (ftnlen)7);
#line 364 "dtrevc3.f"
	return 0;
#line 365 "dtrevc3.f"
    } else if (lquery) {
#line 366 "dtrevc3.f"
	return 0;
#line 367 "dtrevc3.f"
    }

/*     Quick return if possible. */

#line 371 "dtrevc3.f"
    if (*n == 0) {
#line 371 "dtrevc3.f"
	return 0;
#line 371 "dtrevc3.f"
    }

/*     Use blocked version of back-transformation if sufficient workspace. */
/*     Zero-out the workspace to avoid potential NaN propagation. */

#line 377 "dtrevc3.f"
    if (over && *lwork >= *n + (*n << 4)) {
#line 378 "dtrevc3.f"
	nb = (*lwork - *n) / (*n << 1);
#line 379 "dtrevc3.f"
	nb = min(nb,128);
#line 380 "dtrevc3.f"
	i__2 = (nb << 1) + 1;
#line 380 "dtrevc3.f"
	dlaset_("F", n, &i__2, &c_b17, &c_b17, &work[1], n, (ftnlen)1);
#line 381 "dtrevc3.f"
    } else {
#line 382 "dtrevc3.f"
	nb = 1;
#line 383 "dtrevc3.f"
    }

/*     Set the constants to control overflow. */

#line 387 "dtrevc3.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 388 "dtrevc3.f"
    ovfl = 1. / unfl;
#line 389 "dtrevc3.f"
    dlabad_(&unfl, &ovfl);
#line 390 "dtrevc3.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 391 "dtrevc3.f"
    smlnum = unfl * (*n / ulp);
#line 392 "dtrevc3.f"
    bignum = (1. - ulp) / smlnum;

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 397 "dtrevc3.f"
    work[1] = 0.;
#line 398 "dtrevc3.f"
    i__2 = *n;
#line 398 "dtrevc3.f"
    for (j = 2; j <= i__2; ++j) {
#line 399 "dtrevc3.f"
	work[j] = 0.;
#line 400 "dtrevc3.f"
	i__3 = j - 1;
#line 400 "dtrevc3.f"
	for (i__ = 1; i__ <= i__3; ++i__) {
#line 401 "dtrevc3.f"
	    work[j] += (d__1 = t[i__ + j * t_dim1], abs(d__1));
#line 402 "dtrevc3.f"
/* L20: */
#line 402 "dtrevc3.f"
	}
#line 403 "dtrevc3.f"
/* L30: */
#line 403 "dtrevc3.f"
    }

/*     Index IP is used to specify the real or complex eigenvalue: */
/*       IP = 0, real eigenvalue, */
/*            1, first  of conjugate complex pair: (wr,wi) */
/*           -1, second of conjugate complex pair: (wr,wi) */
/*       ISCOMPLEX array stores IP for each column in current block. */

#line 411 "dtrevc3.f"
    if (rightv) {

/*        ============================================================ */
/*        Compute right eigenvectors. */

/*        IV is index of column in current block. */
/*        For complex right vector, uses IV-1 for real part and IV for complex part. */
/*        Non-blocked version always uses IV=2; */
/*        blocked     version starts with IV=NB, goes down to 1 or 2. */
/*        (Note the "0-th" column is used for 1-norms computed above.) */
#line 421 "dtrevc3.f"
	iv = 2;
#line 422 "dtrevc3.f"
	if (nb > 2) {
#line 423 "dtrevc3.f"
	    iv = nb;
#line 424 "dtrevc3.f"
	}
#line 426 "dtrevc3.f"
	ip = 0;
#line 427 "dtrevc3.f"
	is = *m;
#line 428 "dtrevc3.f"
	for (ki = *n; ki >= 1; --ki) {
#line 429 "dtrevc3.f"
	    if (ip == -1) {
/*              previous iteration (ki+1) was second of conjugate pair, */
/*              so this ki is first of conjugate pair; skip to end of loop */
#line 432 "dtrevc3.f"
		ip = 1;
#line 433 "dtrevc3.f"
		goto L140;
#line 434 "dtrevc3.f"
	    } else if (ki == 1) {
/*              last column, so this ki must be real eigenvalue */
#line 436 "dtrevc3.f"
		ip = 0;
#line 437 "dtrevc3.f"
	    } else if (t[ki + (ki - 1) * t_dim1] == 0.) {
/*              zero on sub-diagonal, so this ki is real eigenvalue */
#line 439 "dtrevc3.f"
		ip = 0;
#line 440 "dtrevc3.f"
	    } else {
/*              non-zero on sub-diagonal, so this ki is second of conjugate pair */
#line 442 "dtrevc3.f"
		ip = -1;
#line 443 "dtrevc3.f"
	    }
#line 445 "dtrevc3.f"
	    if (somev) {
#line 446 "dtrevc3.f"
		if (ip == 0) {
#line 447 "dtrevc3.f"
		    if (! select[ki]) {
#line 447 "dtrevc3.f"
			goto L140;
#line 447 "dtrevc3.f"
		    }
#line 449 "dtrevc3.f"
		} else {
#line 450 "dtrevc3.f"
		    if (! select[ki - 1]) {
#line 450 "dtrevc3.f"
			goto L140;
#line 450 "dtrevc3.f"
		    }
#line 452 "dtrevc3.f"
		}
#line 453 "dtrevc3.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 457 "dtrevc3.f"
	    wr = t[ki + ki * t_dim1];
#line 458 "dtrevc3.f"
	    wi = 0.;
#line 459 "dtrevc3.f"
	    if (ip != 0) {
#line 459 "dtrevc3.f"
		wi = sqrt((d__1 = t[ki + (ki - 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki - 1 + ki * t_dim1], abs(d__2)));
#line 459 "dtrevc3.f"
	    }
/* Computing MAX */
#line 462 "dtrevc3.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 462 "dtrevc3.f"
	    smin = max(d__1,smlnum);

#line 464 "dtrevc3.f"
	    if (ip == 0) {

/*              -------------------------------------------------------- */
/*              Real right eigenvector */

#line 469 "dtrevc3.f"
		work[ki + iv * *n] = 1.;

/*              Form right-hand side. */

#line 473 "dtrevc3.f"
		i__2 = ki - 1;
#line 473 "dtrevc3.f"
		for (k = 1; k <= i__2; ++k) {
#line 474 "dtrevc3.f"
		    work[k + iv * *n] = -t[k + ki * t_dim1];
#line 475 "dtrevc3.f"
/* L50: */
#line 475 "dtrevc3.f"
		}

/*              Solve upper quasi-triangular system: */
/*              [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK. */

#line 480 "dtrevc3.f"
		jnxt = ki - 1;
#line 481 "dtrevc3.f"
		for (j = ki - 1; j >= 1; --j) {
#line 482 "dtrevc3.f"
		    if (j > jnxt) {
#line 482 "dtrevc3.f"
			goto L60;
#line 482 "dtrevc3.f"
		    }
#line 484 "dtrevc3.f"
		    j1 = j;
#line 485 "dtrevc3.f"
		    j2 = j;
#line 486 "dtrevc3.f"
		    jnxt = j - 1;
#line 487 "dtrevc3.f"
		    if (j > 1) {
#line 488 "dtrevc3.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 489 "dtrevc3.f"
			    j1 = j - 1;
#line 490 "dtrevc3.f"
			    jnxt = j - 2;
#line 491 "dtrevc3.f"
			}
#line 492 "dtrevc3.f"
		    }

#line 494 "dtrevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 498 "dtrevc3.f"
			dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale X(1,1) to avoid overflow when updating */
/*                    the right-hand side. */

#line 505 "dtrevc3.f"
			if (xnorm > 1.) {
#line 506 "dtrevc3.f"
			    if (work[j] > bignum / xnorm) {
#line 507 "dtrevc3.f"
				x[0] /= xnorm;
#line 508 "dtrevc3.f"
				scale /= xnorm;
#line 509 "dtrevc3.f"
			    }
#line 510 "dtrevc3.f"
			}

/*                    Scale if necessary */

#line 514 "dtrevc3.f"
			if (scale != 1.) {
#line 514 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 514 "dtrevc3.f"
			}
#line 516 "dtrevc3.f"
			work[j + iv * *n] = x[0];

/*                    Update right-hand side */

#line 520 "dtrevc3.f"
			i__2 = j - 1;
#line 520 "dtrevc3.f"
			d__1 = -x[0];
#line 520 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);

#line 523 "dtrevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 527 "dtrevc3.f"
			dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b29, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b29, &c_b29, &
				work[j - 1 + iv * *n], n, &wr, &c_b17, x, &
				c__2, &scale, &xnorm, &ierr);

/*                    Scale X(1,1) and X(2,1) to avoid overflow when */
/*                    updating the right-hand side. */

#line 535 "dtrevc3.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 536 "dtrevc3.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 536 "dtrevc3.f"
			    beta = max(d__1,d__2);
#line 537 "dtrevc3.f"
			    if (beta > bignum / xnorm) {
#line 538 "dtrevc3.f"
				x[0] /= xnorm;
#line 539 "dtrevc3.f"
				x[1] /= xnorm;
#line 540 "dtrevc3.f"
				scale /= xnorm;
#line 541 "dtrevc3.f"
			    }
#line 542 "dtrevc3.f"
			}

/*                    Scale if necessary */

#line 546 "dtrevc3.f"
			if (scale != 1.) {
#line 546 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 546 "dtrevc3.f"
			}
#line 548 "dtrevc3.f"
			work[j - 1 + iv * *n] = x[0];
#line 549 "dtrevc3.f"
			work[j + iv * *n] = x[1];

/*                    Update right-hand side */

#line 553 "dtrevc3.f"
			i__2 = j - 2;
#line 553 "dtrevc3.f"
			d__1 = -x[0];
#line 553 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[iv * *n + 1], &c__1);
#line 555 "dtrevc3.f"
			i__2 = j - 2;
#line 555 "dtrevc3.f"
			d__1 = -x[1];
#line 555 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);
#line 557 "dtrevc3.f"
		    }
#line 558 "dtrevc3.f"
L60:
#line 558 "dtrevc3.f"
		    ;
#line 558 "dtrevc3.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 562 "dtrevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VR and normalize. */
#line 565 "dtrevc3.f"
		    dcopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 
			    1], &c__1);

#line 567 "dtrevc3.f"
		    ii = idamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 568 "dtrevc3.f"
		    remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs(d__1));
#line 569 "dtrevc3.f"
		    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 571 "dtrevc3.f"
		    i__2 = *n;
#line 571 "dtrevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 572 "dtrevc3.f"
			vr[k + is * vr_dim1] = 0.;
#line 573 "dtrevc3.f"
/* L70: */
#line 573 "dtrevc3.f"
		    }

#line 575 "dtrevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 578 "dtrevc3.f"
		    if (ki > 1) {
#line 578 "dtrevc3.f"
			i__2 = ki - 1;
#line 578 "dtrevc3.f"
			dgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[iv * *n + 1], &c__1, &work[ki + iv * *n],
				 &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 578 "dtrevc3.f"
		    }

#line 583 "dtrevc3.f"
		    ii = idamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 584 "dtrevc3.f"
		    remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs(d__1));
#line 585 "dtrevc3.f"
		    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);

#line 587 "dtrevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out below vector */
#line 591 "dtrevc3.f"
		    i__2 = *n;
#line 591 "dtrevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 592 "dtrevc3.f"
			work[k + iv * *n] = 0.;
#line 593 "dtrevc3.f"
		    }
#line 594 "dtrevc3.f"
		    iscomplex[iv - 1] = ip;
/*                 back-transform and normalization is done below */
#line 596 "dtrevc3.f"
		}
#line 597 "dtrevc3.f"
	    } else {

/*              -------------------------------------------------------- */
/*              Complex right eigenvector. */

/*              Initial solve */
/*              [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0. */
/*              [ ( T(KI,  KI-1) T(KI,  KI) )               ] */

#line 606 "dtrevc3.f"
		if ((d__1 = t[ki - 1 + ki * t_dim1], abs(d__1)) >= (d__2 = t[
			ki + (ki - 1) * t_dim1], abs(d__2))) {
#line 607 "dtrevc3.f"
		    work[ki - 1 + (iv - 1) * *n] = 1.;
#line 608 "dtrevc3.f"
		    work[ki + iv * *n] = wi / t[ki - 1 + ki * t_dim1];
#line 609 "dtrevc3.f"
		} else {
#line 610 "dtrevc3.f"
		    work[ki - 1 + (iv - 1) * *n] = -wi / t[ki + (ki - 1) * 
			    t_dim1];
#line 611 "dtrevc3.f"
		    work[ki + iv * *n] = 1.;
#line 612 "dtrevc3.f"
		}
#line 613 "dtrevc3.f"
		work[ki + (iv - 1) * *n] = 0.;
#line 614 "dtrevc3.f"
		work[ki - 1 + iv * *n] = 0.;

/*              Form right-hand side. */

#line 618 "dtrevc3.f"
		i__2 = ki - 2;
#line 618 "dtrevc3.f"
		for (k = 1; k <= i__2; ++k) {
#line 619 "dtrevc3.f"
		    work[k + (iv - 1) * *n] = -work[ki - 1 + (iv - 1) * *n] * 
			    t[k + (ki - 1) * t_dim1];
#line 620 "dtrevc3.f"
		    work[k + iv * *n] = -work[ki + iv * *n] * t[k + ki * 
			    t_dim1];
#line 621 "dtrevc3.f"
/* L80: */
#line 621 "dtrevc3.f"
		}

/*              Solve upper quasi-triangular system: */
/*              [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2) */

#line 626 "dtrevc3.f"
		jnxt = ki - 2;
#line 627 "dtrevc3.f"
		for (j = ki - 2; j >= 1; --j) {
#line 628 "dtrevc3.f"
		    if (j > jnxt) {
#line 628 "dtrevc3.f"
			goto L90;
#line 628 "dtrevc3.f"
		    }
#line 630 "dtrevc3.f"
		    j1 = j;
#line 631 "dtrevc3.f"
		    j2 = j;
#line 632 "dtrevc3.f"
		    jnxt = j - 1;
#line 633 "dtrevc3.f"
		    if (j > 1) {
#line 634 "dtrevc3.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 635 "dtrevc3.f"
			    j1 = j - 1;
#line 636 "dtrevc3.f"
			    jnxt = j - 2;
#line 637 "dtrevc3.f"
			}
#line 638 "dtrevc3.f"
		    }

#line 640 "dtrevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 644 "dtrevc3.f"
			dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + (
				iv - 1) * *n], n, &wr, &wi, x, &c__2, &scale, 
				&xnorm, &ierr);

/*                    Scale X(1,1) and X(1,2) to avoid overflow when */
/*                    updating the right-hand side. */

#line 651 "dtrevc3.f"
			if (xnorm > 1.) {
#line 652 "dtrevc3.f"
			    if (work[j] > bignum / xnorm) {
#line 653 "dtrevc3.f"
				x[0] /= xnorm;
#line 654 "dtrevc3.f"
				x[2] /= xnorm;
#line 655 "dtrevc3.f"
				scale /= xnorm;
#line 656 "dtrevc3.f"
			    }
#line 657 "dtrevc3.f"
			}

/*                    Scale if necessary */

#line 661 "dtrevc3.f"
			if (scale != 1.) {
#line 662 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &
				    c__1);
#line 663 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 664 "dtrevc3.f"
			}
#line 665 "dtrevc3.f"
			work[j + (iv - 1) * *n] = x[0];
#line 666 "dtrevc3.f"
			work[j + iv * *n] = x[2];

/*                    Update the right-hand side */

#line 670 "dtrevc3.f"
			i__2 = j - 1;
#line 670 "dtrevc3.f"
			d__1 = -x[0];
#line 670 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				(iv - 1) * *n + 1], &c__1);
#line 672 "dtrevc3.f"
			i__2 = j - 1;
#line 672 "dtrevc3.f"
			d__1 = -x[2];
#line 672 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);

#line 675 "dtrevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 679 "dtrevc3.f"
			dlaln2_(&c_false, &c__2, &c__2, &smin, &c_b29, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b29, &c_b29, &
				work[j - 1 + (iv - 1) * *n], n, &wr, &wi, x, &
				c__2, &scale, &xnorm, &ierr);

/*                    Scale X to avoid overflow when updating */
/*                    the right-hand side. */

#line 687 "dtrevc3.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 688 "dtrevc3.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 688 "dtrevc3.f"
			    beta = max(d__1,d__2);
#line 689 "dtrevc3.f"
			    if (beta > bignum / xnorm) {
#line 690 "dtrevc3.f"
				rec = 1. / xnorm;
#line 691 "dtrevc3.f"
				x[0] *= rec;
#line 692 "dtrevc3.f"
				x[2] *= rec;
#line 693 "dtrevc3.f"
				x[1] *= rec;
#line 694 "dtrevc3.f"
				x[3] *= rec;
#line 695 "dtrevc3.f"
				scale *= rec;
#line 696 "dtrevc3.f"
			    }
#line 697 "dtrevc3.f"
			}

/*                    Scale if necessary */

#line 701 "dtrevc3.f"
			if (scale != 1.) {
#line 702 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &
				    c__1);
#line 703 "dtrevc3.f"
			    dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 704 "dtrevc3.f"
			}
#line 705 "dtrevc3.f"
			work[j - 1 + (iv - 1) * *n] = x[0];
#line 706 "dtrevc3.f"
			work[j + (iv - 1) * *n] = x[1];
#line 707 "dtrevc3.f"
			work[j - 1 + iv * *n] = x[2];
#line 708 "dtrevc3.f"
			work[j + iv * *n] = x[3];

/*                    Update the right-hand side */

#line 712 "dtrevc3.f"
			i__2 = j - 2;
#line 712 "dtrevc3.f"
			d__1 = -x[0];
#line 712 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[(iv - 1) * *n + 1], &c__1);
#line 714 "dtrevc3.f"
			i__2 = j - 2;
#line 714 "dtrevc3.f"
			d__1 = -x[1];
#line 714 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				(iv - 1) * *n + 1], &c__1);
#line 716 "dtrevc3.f"
			i__2 = j - 2;
#line 716 "dtrevc3.f"
			d__1 = -x[2];
#line 716 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[iv * *n + 1], &c__1);
#line 718 "dtrevc3.f"
			i__2 = j - 2;
#line 718 "dtrevc3.f"
			d__1 = -x[3];
#line 718 "dtrevc3.f"
			daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);
#line 720 "dtrevc3.f"
		    }
#line 721 "dtrevc3.f"
L90:
#line 721 "dtrevc3.f"
		    ;
#line 721 "dtrevc3.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 725 "dtrevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VR and normalize. */
#line 728 "dtrevc3.f"
		    dcopy_(&ki, &work[(iv - 1) * *n + 1], &c__1, &vr[(is - 1) 
			    * vr_dim1 + 1], &c__1);
#line 729 "dtrevc3.f"
		    dcopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 
			    1], &c__1);

#line 731 "dtrevc3.f"
		    emax = 0.;
#line 732 "dtrevc3.f"
		    i__2 = ki;
#line 732 "dtrevc3.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 733 "dtrevc3.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + is * vr_dim1], 
				abs(d__2));
#line 733 "dtrevc3.f"
			emax = max(d__3,d__4);
#line 735 "dtrevc3.f"
/* L100: */
#line 735 "dtrevc3.f"
		    }
#line 736 "dtrevc3.f"
		    remax = 1. / emax;
#line 737 "dtrevc3.f"
		    dscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
#line 738 "dtrevc3.f"
		    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 740 "dtrevc3.f"
		    i__2 = *n;
#line 740 "dtrevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 741 "dtrevc3.f"
			vr[k + (is - 1) * vr_dim1] = 0.;
#line 742 "dtrevc3.f"
			vr[k + is * vr_dim1] = 0.;
#line 743 "dtrevc3.f"
/* L110: */
#line 743 "dtrevc3.f"
		    }

#line 745 "dtrevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 748 "dtrevc3.f"
		    if (ki > 2) {
#line 749 "dtrevc3.f"
			i__2 = ki - 2;
#line 749 "dtrevc3.f"
			dgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[(iv - 1) * *n + 1], &c__1, &work[ki - 1 
				+ (iv - 1) * *n], &vr[(ki - 1) * vr_dim1 + 1],
				 &c__1, (ftnlen)1);
#line 752 "dtrevc3.f"
			i__2 = ki - 2;
#line 752 "dtrevc3.f"
			dgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[iv * *n + 1], &c__1, &work[ki + iv * *n],
				 &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 755 "dtrevc3.f"
		    } else {
#line 756 "dtrevc3.f"
			dscal_(n, &work[ki - 1 + (iv - 1) * *n], &vr[(ki - 1) 
				* vr_dim1 + 1], &c__1);
#line 757 "dtrevc3.f"
			dscal_(n, &work[ki + iv * *n], &vr[ki * vr_dim1 + 1], 
				&c__1);
#line 758 "dtrevc3.f"
		    }

#line 760 "dtrevc3.f"
		    emax = 0.;
#line 761 "dtrevc3.f"
		    i__2 = *n;
#line 761 "dtrevc3.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 762 "dtrevc3.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + ki * vr_dim1], 
				abs(d__2));
#line 762 "dtrevc3.f"
			emax = max(d__3,d__4);
#line 764 "dtrevc3.f"
/* L120: */
#line 764 "dtrevc3.f"
		    }
#line 765 "dtrevc3.f"
		    remax = 1. / emax;
#line 766 "dtrevc3.f"
		    dscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
#line 767 "dtrevc3.f"
		    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);

#line 769 "dtrevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out below vector */
#line 773 "dtrevc3.f"
		    i__2 = *n;
#line 773 "dtrevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 774 "dtrevc3.f"
			work[k + (iv - 1) * *n] = 0.;
#line 775 "dtrevc3.f"
			work[k + iv * *n] = 0.;
#line 776 "dtrevc3.f"
		    }
#line 777 "dtrevc3.f"
		    iscomplex[iv - 2] = -ip;
#line 778 "dtrevc3.f"
		    iscomplex[iv - 1] = ip;
#line 779 "dtrevc3.f"
		    --iv;
/*                 back-transform and normalization is done below */
#line 781 "dtrevc3.f"
		}
#line 782 "dtrevc3.f"
	    }
#line 784 "dtrevc3.f"
	    if (nb > 1) {
/*              -------------------------------------------------------- */
/*              Blocked version of back-transform */
/*              For complex case, KI2 includes both vectors (KI-1 and KI) */
#line 788 "dtrevc3.f"
		if (ip == 0) {
#line 789 "dtrevc3.f"
		    ki2 = ki;
#line 790 "dtrevc3.f"
		} else {
#line 791 "dtrevc3.f"
		    ki2 = ki - 1;
#line 792 "dtrevc3.f"
		}
/*              Columns IV:NB of work are valid vectors. */
/*              When the number of vectors stored reaches NB-1 or NB, */
/*              or if this was last vector, do the GEMM */
#line 797 "dtrevc3.f"
		if (iv <= 2 || ki2 == 1) {
#line 798 "dtrevc3.f"
		    i__2 = nb - iv + 1;
#line 798 "dtrevc3.f"
		    i__3 = ki2 + nb - iv;
#line 798 "dtrevc3.f"
		    dgemm_("N", "N", n, &i__2, &i__3, &c_b29, &vr[vr_offset], 
			    ldvr, &work[iv * *n + 1], n, &c_b17, &work[(nb + 
			    iv) * *n + 1], n, (ftnlen)1, (ftnlen)1);
/*                 normalize vectors */
#line 804 "dtrevc3.f"
		    i__2 = nb;
#line 804 "dtrevc3.f"
		    for (k = iv; k <= i__2; ++k) {
#line 805 "dtrevc3.f"
			if (iscomplex[k - 1] == 0) {
/*                       real eigenvector */
#line 807 "dtrevc3.f"
			    ii = idamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 808 "dtrevc3.f"
			    remax = 1. / (d__1 = work[ii + (nb + k) * *n], 
				    abs(d__1));
#line 809 "dtrevc3.f"
			} else if (iscomplex[k - 1] == 1) {
/*                       first eigenvector of conjugate pair */
#line 811 "dtrevc3.f"
			    emax = 0.;
#line 812 "dtrevc3.f"
			    i__3 = *n;
#line 812 "dtrevc3.f"
			    for (ii = 1; ii <= i__3; ++ii) {
/* Computing MAX */
#line 813 "dtrevc3.f"
				d__3 = emax, d__4 = (d__1 = work[ii + (nb + k)
					 * *n], abs(d__1)) + (d__2 = work[ii 
					+ (nb + k + 1) * *n], abs(d__2));
#line 813 "dtrevc3.f"
				emax = max(d__3,d__4);
#line 816 "dtrevc3.f"
			    }
#line 817 "dtrevc3.f"
			    remax = 1. / emax;
/*                    else if ISCOMPLEX(K).EQ.-1 */
/*                       second eigenvector of conjugate pair */
/*                       reuse same REMAX as previous K */
#line 821 "dtrevc3.f"
			}
#line 822 "dtrevc3.f"
			dscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 823 "dtrevc3.f"
		    }
#line 824 "dtrevc3.f"
		    i__2 = nb - iv + 1;
#line 824 "dtrevc3.f"
		    dlacpy_("F", n, &i__2, &work[(nb + iv) * *n + 1], n, &vr[
			    ki2 * vr_dim1 + 1], ldvr, (ftnlen)1);
#line 827 "dtrevc3.f"
		    iv = nb;
#line 828 "dtrevc3.f"
		} else {
#line 829 "dtrevc3.f"
		    --iv;
#line 830 "dtrevc3.f"
		}
#line 831 "dtrevc3.f"
	    }

/* blocked back-transform */
#line 833 "dtrevc3.f"
	    --is;
#line 834 "dtrevc3.f"
	    if (ip != 0) {
#line 834 "dtrevc3.f"
		--is;
#line 834 "dtrevc3.f"
	    }
#line 836 "dtrevc3.f"
L140:
#line 836 "dtrevc3.f"
	    ;
#line 836 "dtrevc3.f"
	}
#line 837 "dtrevc3.f"
    }
#line 839 "dtrevc3.f"
    if (leftv) {

/*        ============================================================ */
/*        Compute left eigenvectors. */

/*        IV is index of column in current block. */
/*        For complex left vector, uses IV for real part and IV+1 for complex part. */
/*        Non-blocked version always uses IV=1; */
/*        blocked     version starts with IV=1, goes up to NB-1 or NB. */
/*        (Note the "0-th" column is used for 1-norms computed above.) */
#line 849 "dtrevc3.f"
	iv = 1;
#line 850 "dtrevc3.f"
	ip = 0;
#line 851 "dtrevc3.f"
	is = 1;
#line 852 "dtrevc3.f"
	i__2 = *n;
#line 852 "dtrevc3.f"
	for (ki = 1; ki <= i__2; ++ki) {
#line 853 "dtrevc3.f"
	    if (ip == 1) {
/*              previous iteration (ki-1) was first of conjugate pair, */
/*              so this ki is second of conjugate pair; skip to end of loop */
#line 856 "dtrevc3.f"
		ip = -1;
#line 857 "dtrevc3.f"
		goto L260;
#line 858 "dtrevc3.f"
	    } else if (ki == *n) {
/*              last column, so this ki must be real eigenvalue */
#line 860 "dtrevc3.f"
		ip = 0;
#line 861 "dtrevc3.f"
	    } else if (t[ki + 1 + ki * t_dim1] == 0.) {
/*              zero on sub-diagonal, so this ki is real eigenvalue */
#line 863 "dtrevc3.f"
		ip = 0;
#line 864 "dtrevc3.f"
	    } else {
/*              non-zero on sub-diagonal, so this ki is first of conjugate pair */
#line 866 "dtrevc3.f"
		ip = 1;
#line 867 "dtrevc3.f"
	    }

#line 869 "dtrevc3.f"
	    if (somev) {
#line 870 "dtrevc3.f"
		if (! select[ki]) {
#line 870 "dtrevc3.f"
		    goto L260;
#line 870 "dtrevc3.f"
		}
#line 872 "dtrevc3.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 876 "dtrevc3.f"
	    wr = t[ki + ki * t_dim1];
#line 877 "dtrevc3.f"
	    wi = 0.;
#line 878 "dtrevc3.f"
	    if (ip != 0) {
#line 878 "dtrevc3.f"
		wi = sqrt((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki + 1 + ki * t_dim1], abs(d__2)));
#line 878 "dtrevc3.f"
	    }
/* Computing MAX */
#line 881 "dtrevc3.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 881 "dtrevc3.f"
	    smin = max(d__1,smlnum);

#line 883 "dtrevc3.f"
	    if (ip == 0) {

/*              -------------------------------------------------------- */
/*              Real left eigenvector */

#line 888 "dtrevc3.f"
		work[ki + iv * *n] = 1.;

/*              Form right-hand side. */

#line 892 "dtrevc3.f"
		i__3 = *n;
#line 892 "dtrevc3.f"
		for (k = ki + 1; k <= i__3; ++k) {
#line 893 "dtrevc3.f"
		    work[k + iv * *n] = -t[ki + k * t_dim1];
#line 894 "dtrevc3.f"
/* L160: */
#line 894 "dtrevc3.f"
		}

/*              Solve transposed quasi-triangular system: */
/*              [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK */

#line 899 "dtrevc3.f"
		vmax = 1.;
#line 900 "dtrevc3.f"
		vcrit = bignum;

#line 902 "dtrevc3.f"
		jnxt = ki + 1;
#line 903 "dtrevc3.f"
		i__3 = *n;
#line 903 "dtrevc3.f"
		for (j = ki + 1; j <= i__3; ++j) {
#line 904 "dtrevc3.f"
		    if (j < jnxt) {
#line 904 "dtrevc3.f"
			goto L170;
#line 904 "dtrevc3.f"
		    }
#line 906 "dtrevc3.f"
		    j1 = j;
#line 907 "dtrevc3.f"
		    j2 = j;
#line 908 "dtrevc3.f"
		    jnxt = j + 1;
#line 909 "dtrevc3.f"
		    if (j < *n) {
#line 910 "dtrevc3.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 911 "dtrevc3.f"
			    j2 = j + 1;
#line 912 "dtrevc3.f"
			    jnxt = j + 2;
#line 913 "dtrevc3.f"
			}
#line 914 "dtrevc3.f"
		    }

#line 916 "dtrevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

#line 923 "dtrevc3.f"
			if (work[j] > vcrit) {
#line 924 "dtrevc3.f"
			    rec = 1. / vmax;
#line 925 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 925 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 926 "dtrevc3.f"
			    vmax = 1.;
#line 927 "dtrevc3.f"
			    vcrit = bignum;
#line 928 "dtrevc3.f"
			}

#line 930 "dtrevc3.f"
			i__4 = j - ki - 1;
#line 930 "dtrevc3.f"
			work[j + iv * *n] -= ddot_(&i__4, &t[ki + 1 + j * 
				t_dim1], &c__1, &work[ki + 1 + iv * *n], &
				c__1);

/*                    Solve [ T(J,J) - WR ]**T * X = WORK */

#line 936 "dtrevc3.f"
			dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 942 "dtrevc3.f"
			if (scale != 1.) {
#line 942 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 942 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 942 "dtrevc3.f"
			}
#line 944 "dtrevc3.f"
			work[j + iv * *n] = x[0];
/* Computing MAX */
#line 945 "dtrevc3.f"
			d__2 = (d__1 = work[j + iv * *n], abs(d__1));
#line 945 "dtrevc3.f"
			vmax = max(d__2,vmax);
#line 946 "dtrevc3.f"
			vcrit = bignum / vmax;

#line 948 "dtrevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

/* Computing MAX */
#line 955 "dtrevc3.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 955 "dtrevc3.f"
			beta = max(d__1,d__2);
#line 956 "dtrevc3.f"
			if (beta > vcrit) {
#line 957 "dtrevc3.f"
			    rec = 1. / vmax;
#line 958 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 958 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 959 "dtrevc3.f"
			    vmax = 1.;
#line 960 "dtrevc3.f"
			    vcrit = bignum;
#line 961 "dtrevc3.f"
			}

#line 963 "dtrevc3.f"
			i__4 = j - ki - 1;
#line 963 "dtrevc3.f"
			work[j + iv * *n] -= ddot_(&i__4, &t[ki + 1 + j * 
				t_dim1], &c__1, &work[ki + 1 + iv * *n], &
				c__1);

#line 967 "dtrevc3.f"
			i__4 = j - ki - 1;
#line 967 "dtrevc3.f"
			work[j + 1 + iv * *n] -= ddot_(&i__4, &t[ki + 1 + (j 
				+ 1) * t_dim1], &c__1, &work[ki + 1 + iv * *n]
				, &c__1);

/*                    Solve */
/*                    [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 ) */
/*                    [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 ) */

#line 975 "dtrevc3.f"
			dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 981 "dtrevc3.f"
			if (scale != 1.) {
#line 981 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 981 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 981 "dtrevc3.f"
			}
#line 983 "dtrevc3.f"
			work[j + iv * *n] = x[0];
#line 984 "dtrevc3.f"
			work[j + 1 + iv * *n] = x[1];

/* Computing MAX */
#line 986 "dtrevc3.f"
			d__3 = (d__1 = work[j + iv * *n], abs(d__1)), d__4 = (
				d__2 = work[j + 1 + iv * *n], abs(d__2)), 
				d__3 = max(d__3,d__4);
#line 986 "dtrevc3.f"
			vmax = max(d__3,vmax);
#line 988 "dtrevc3.f"
			vcrit = bignum / vmax;

#line 990 "dtrevc3.f"
		    }
#line 991 "dtrevc3.f"
L170:
#line 991 "dtrevc3.f"
		    ;
#line 991 "dtrevc3.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 995 "dtrevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VL and normalize. */
#line 998 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 998 "dtrevc3.f"
		    dcopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);

#line 1001 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1001 "dtrevc3.f"
		    ii = idamax_(&i__3, &vl[ki + is * vl_dim1], &c__1) + ki - 
			    1;
#line 1002 "dtrevc3.f"
		    remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs(d__1));
#line 1003 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1003 "dtrevc3.f"
		    dscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 1005 "dtrevc3.f"
		    i__3 = ki - 1;
#line 1005 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1006 "dtrevc3.f"
			vl[k + is * vl_dim1] = 0.;
#line 1007 "dtrevc3.f"
/* L180: */
#line 1007 "dtrevc3.f"
		    }

#line 1009 "dtrevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 1012 "dtrevc3.f"
		    if (ki < *n) {
#line 1012 "dtrevc3.f"
			i__3 = *n - ki;
#line 1012 "dtrevc3.f"
			dgemv_("N", n, &i__3, &c_b29, &vl[(ki + 1) * vl_dim1 
				+ 1], ldvl, &work[ki + 1 + iv * *n], &c__1, &
				work[ki + iv * *n], &vl[ki * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1012 "dtrevc3.f"
		    }

#line 1018 "dtrevc3.f"
		    ii = idamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 1019 "dtrevc3.f"
		    remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs(d__1));
#line 1020 "dtrevc3.f"
		    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

#line 1022 "dtrevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out above vector */
/*                 could go from KI-NV+1 to KI-1 */
#line 1027 "dtrevc3.f"
		    i__3 = ki - 1;
#line 1027 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1028 "dtrevc3.f"
			work[k + iv * *n] = 0.;
#line 1029 "dtrevc3.f"
		    }
#line 1030 "dtrevc3.f"
		    iscomplex[iv - 1] = ip;
/*                 back-transform and normalization is done below */
#line 1032 "dtrevc3.f"
		}
#line 1033 "dtrevc3.f"
	    } else {

/*              -------------------------------------------------------- */
/*              Complex left eigenvector. */

/*              Initial solve: */
/*              [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0. */
/*              [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ] */

#line 1042 "dtrevc3.f"
		if ((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1)) >= (d__2 = 
			t[ki + 1 + ki * t_dim1], abs(d__2))) {
#line 1043 "dtrevc3.f"
		    work[ki + iv * *n] = wi / t[ki + (ki + 1) * t_dim1];
#line 1044 "dtrevc3.f"
		    work[ki + 1 + (iv + 1) * *n] = 1.;
#line 1045 "dtrevc3.f"
		} else {
#line 1046 "dtrevc3.f"
		    work[ki + iv * *n] = 1.;
#line 1047 "dtrevc3.f"
		    work[ki + 1 + (iv + 1) * *n] = -wi / t[ki + 1 + ki * 
			    t_dim1];
#line 1048 "dtrevc3.f"
		}
#line 1049 "dtrevc3.f"
		work[ki + 1 + iv * *n] = 0.;
#line 1050 "dtrevc3.f"
		work[ki + (iv + 1) * *n] = 0.;

/*              Form right-hand side. */

#line 1054 "dtrevc3.f"
		i__3 = *n;
#line 1054 "dtrevc3.f"
		for (k = ki + 2; k <= i__3; ++k) {
#line 1055 "dtrevc3.f"
		    work[k + iv * *n] = -work[ki + iv * *n] * t[ki + k * 
			    t_dim1];
#line 1056 "dtrevc3.f"
		    work[k + (iv + 1) * *n] = -work[ki + 1 + (iv + 1) * *n] * 
			    t[ki + 1 + k * t_dim1];
#line 1057 "dtrevc3.f"
/* L190: */
#line 1057 "dtrevc3.f"
		}

/*              Solve transposed quasi-triangular system: */
/*              [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2 */

#line 1062 "dtrevc3.f"
		vmax = 1.;
#line 1063 "dtrevc3.f"
		vcrit = bignum;

#line 1065 "dtrevc3.f"
		jnxt = ki + 2;
#line 1066 "dtrevc3.f"
		i__3 = *n;
#line 1066 "dtrevc3.f"
		for (j = ki + 2; j <= i__3; ++j) {
#line 1067 "dtrevc3.f"
		    if (j < jnxt) {
#line 1067 "dtrevc3.f"
			goto L200;
#line 1067 "dtrevc3.f"
		    }
#line 1069 "dtrevc3.f"
		    j1 = j;
#line 1070 "dtrevc3.f"
		    j2 = j;
#line 1071 "dtrevc3.f"
		    jnxt = j + 1;
#line 1072 "dtrevc3.f"
		    if (j < *n) {
#line 1073 "dtrevc3.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 1074 "dtrevc3.f"
			    j2 = j + 1;
#line 1075 "dtrevc3.f"
			    jnxt = j + 2;
#line 1076 "dtrevc3.f"
			}
#line 1077 "dtrevc3.f"
		    }

#line 1079 "dtrevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when */
/*                    forming the right-hand side elements. */

#line 1086 "dtrevc3.f"
			if (work[j] > vcrit) {
#line 1087 "dtrevc3.f"
			    rec = 1. / vmax;
#line 1088 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1088 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 1089 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1089 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1090 "dtrevc3.f"
			    vmax = 1.;
#line 1091 "dtrevc3.f"
			    vcrit = bignum;
#line 1092 "dtrevc3.f"
			}

#line 1094 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1094 "dtrevc3.f"
			work[j + iv * *n] -= ddot_(&i__4, &t[ki + 2 + j * 
				t_dim1], &c__1, &work[ki + 2 + iv * *n], &
				c__1);
#line 1097 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1097 "dtrevc3.f"
			work[j + (iv + 1) * *n] -= ddot_(&i__4, &t[ki + 2 + j 
				* t_dim1], &c__1, &work[ki + 2 + (iv + 1) * *
				n], &c__1);

/*                    Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2 */

#line 1103 "dtrevc3.f"
			d__1 = -wi;
#line 1103 "dtrevc3.f"
			dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &d__1, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 1109 "dtrevc3.f"
			if (scale != 1.) {
#line 1110 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1110 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 1111 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1111 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1112 "dtrevc3.f"
			}
#line 1113 "dtrevc3.f"
			work[j + iv * *n] = x[0];
#line 1114 "dtrevc3.f"
			work[j + (iv + 1) * *n] = x[2];
/* Computing MAX */
#line 1115 "dtrevc3.f"
			d__3 = (d__1 = work[j + iv * *n], abs(d__1)), d__4 = (
				d__2 = work[j + (iv + 1) * *n], abs(d__2)), 
				d__3 = max(d__3,d__4);
#line 1115 "dtrevc3.f"
			vmax = max(d__3,vmax);
#line 1117 "dtrevc3.f"
			vcrit = bignum / vmax;

#line 1119 "dtrevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side elements. */

/* Computing MAX */
#line 1126 "dtrevc3.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 1126 "dtrevc3.f"
			beta = max(d__1,d__2);
#line 1127 "dtrevc3.f"
			if (beta > vcrit) {
#line 1128 "dtrevc3.f"
			    rec = 1. / vmax;
#line 1129 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1129 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 1130 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1130 "dtrevc3.f"
			    dscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1131 "dtrevc3.f"
			    vmax = 1.;
#line 1132 "dtrevc3.f"
			    vcrit = bignum;
#line 1133 "dtrevc3.f"
			}

#line 1135 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1135 "dtrevc3.f"
			work[j + iv * *n] -= ddot_(&i__4, &t[ki + 2 + j * 
				t_dim1], &c__1, &work[ki + 2 + iv * *n], &
				c__1);

#line 1139 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1139 "dtrevc3.f"
			work[j + (iv + 1) * *n] -= ddot_(&i__4, &t[ki + 2 + j 
				* t_dim1], &c__1, &work[ki + 2 + (iv + 1) * *
				n], &c__1);

#line 1143 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1143 "dtrevc3.f"
			work[j + 1 + iv * *n] -= ddot_(&i__4, &t[ki + 2 + (j 
				+ 1) * t_dim1], &c__1, &work[ki + 2 + iv * *n]
				, &c__1);

#line 1147 "dtrevc3.f"
			i__4 = j - ki - 2;
#line 1147 "dtrevc3.f"
			work[j + 1 + (iv + 1) * *n] -= ddot_(&i__4, &t[ki + 2 
				+ (j + 1) * t_dim1], &c__1, &work[ki + 2 + (
				iv + 1) * *n], &c__1);

/*                    Solve 2-by-2 complex linear equation */
/*                    [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B */
/*                    [ (T(j+1,j) T(j+1,j+1))                  ] */

#line 1155 "dtrevc3.f"
			d__1 = -wi;
#line 1155 "dtrevc3.f"
			dlaln2_(&c_true, &c__2, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &d__1, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 1161 "dtrevc3.f"
			if (scale != 1.) {
#line 1162 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1162 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 1163 "dtrevc3.f"
			    i__4 = *n - ki + 1;
#line 1163 "dtrevc3.f"
			    dscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1164 "dtrevc3.f"
			}
#line 1165 "dtrevc3.f"
			work[j + iv * *n] = x[0];
#line 1166 "dtrevc3.f"
			work[j + (iv + 1) * *n] = x[2];
#line 1167 "dtrevc3.f"
			work[j + 1 + iv * *n] = x[1];
#line 1168 "dtrevc3.f"
			work[j + 1 + (iv + 1) * *n] = x[3];
/* Computing MAX */
#line 1169 "dtrevc3.f"
			d__1 = abs(x[0]), d__2 = abs(x[2]), d__1 = max(d__1,
				d__2), d__2 = abs(x[1]), d__1 = max(d__1,d__2)
				, d__2 = abs(x[3]), d__1 = max(d__1,d__2);
#line 1169 "dtrevc3.f"
			vmax = max(d__1,vmax);
#line 1172 "dtrevc3.f"
			vcrit = bignum / vmax;

#line 1174 "dtrevc3.f"
		    }
#line 1175 "dtrevc3.f"
L200:
#line 1175 "dtrevc3.f"
		    ;
#line 1175 "dtrevc3.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 1179 "dtrevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VL and normalize. */
#line 1182 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1182 "dtrevc3.f"
		    dcopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);
#line 1184 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1184 "dtrevc3.f"
		    dcopy_(&i__3, &work[ki + (iv + 1) * *n], &c__1, &vl[ki + (
			    is + 1) * vl_dim1], &c__1);

#line 1187 "dtrevc3.f"
		    emax = 0.;
#line 1188 "dtrevc3.f"
		    i__3 = *n;
#line 1188 "dtrevc3.f"
		    for (k = ki; k <= i__3; ++k) {
/* Computing MAX */
#line 1189 "dtrevc3.f"
			d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (is + 1) * vl_dim1], 
				abs(d__2));
#line 1189 "dtrevc3.f"
			emax = max(d__3,d__4);
#line 1191 "dtrevc3.f"
/* L220: */
#line 1191 "dtrevc3.f"
		    }
#line 1192 "dtrevc3.f"
		    remax = 1. / emax;
#line 1193 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1193 "dtrevc3.f"
		    dscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);
#line 1194 "dtrevc3.f"
		    i__3 = *n - ki + 1;
#line 1194 "dtrevc3.f"
		    dscal_(&i__3, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1)
			    ;

#line 1196 "dtrevc3.f"
		    i__3 = ki - 1;
#line 1196 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1197 "dtrevc3.f"
			vl[k + is * vl_dim1] = 0.;
#line 1198 "dtrevc3.f"
			vl[k + (is + 1) * vl_dim1] = 0.;
#line 1199 "dtrevc3.f"
/* L230: */
#line 1199 "dtrevc3.f"
		    }

#line 1201 "dtrevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 1204 "dtrevc3.f"
		    if (ki < *n - 1) {
#line 1205 "dtrevc3.f"
			i__3 = *n - ki - 1;
#line 1205 "dtrevc3.f"
			dgemv_("N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + iv * *n], &c__1, &
				work[ki + iv * *n], &vl[ki * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1210 "dtrevc3.f"
			i__3 = *n - ki - 1;
#line 1210 "dtrevc3.f"
			dgemv_("N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + (iv + 1) * *n], &
				c__1, &work[ki + 1 + (iv + 1) * *n], &vl[(ki 
				+ 1) * vl_dim1 + 1], &c__1, (ftnlen)1);
#line 1215 "dtrevc3.f"
		    } else {
#line 1216 "dtrevc3.f"
			dscal_(n, &work[ki + iv * *n], &vl[ki * vl_dim1 + 1], 
				&c__1);
#line 1217 "dtrevc3.f"
			dscal_(n, &work[ki + 1 + (iv + 1) * *n], &vl[(ki + 1) 
				* vl_dim1 + 1], &c__1);
#line 1218 "dtrevc3.f"
		    }

#line 1220 "dtrevc3.f"
		    emax = 0.;
#line 1221 "dtrevc3.f"
		    i__3 = *n;
#line 1221 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1222 "dtrevc3.f"
			d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (ki + 1) * vl_dim1], 
				abs(d__2));
#line 1222 "dtrevc3.f"
			emax = max(d__3,d__4);
#line 1224 "dtrevc3.f"
/* L240: */
#line 1224 "dtrevc3.f"
		    }
#line 1225 "dtrevc3.f"
		    remax = 1. / emax;
#line 1226 "dtrevc3.f"
		    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
#line 1227 "dtrevc3.f"
		    dscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);

#line 1229 "dtrevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out above vector */
/*                 could go from KI-NV+1 to KI-1 */
#line 1234 "dtrevc3.f"
		    i__3 = ki - 1;
#line 1234 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1235 "dtrevc3.f"
			work[k + iv * *n] = 0.;
#line 1236 "dtrevc3.f"
			work[k + (iv + 1) * *n] = 0.;
#line 1237 "dtrevc3.f"
		    }
#line 1238 "dtrevc3.f"
		    iscomplex[iv - 1] = ip;
#line 1239 "dtrevc3.f"
		    iscomplex[iv] = -ip;
#line 1240 "dtrevc3.f"
		    ++iv;
/*                 back-transform and normalization is done below */
#line 1242 "dtrevc3.f"
		}
#line 1243 "dtrevc3.f"
	    }
#line 1245 "dtrevc3.f"
	    if (nb > 1) {
/*              -------------------------------------------------------- */
/*              Blocked version of back-transform */
/*              For complex case, KI2 includes both vectors (KI and KI+1) */
#line 1249 "dtrevc3.f"
		if (ip == 0) {
#line 1250 "dtrevc3.f"
		    ki2 = ki;
#line 1251 "dtrevc3.f"
		} else {
#line 1252 "dtrevc3.f"
		    ki2 = ki + 1;
#line 1253 "dtrevc3.f"
		}
/*              Columns 1:IV of work are valid vectors. */
/*              When the number of vectors stored reaches NB-1 or NB, */
/*              or if this was last vector, do the GEMM */
#line 1258 "dtrevc3.f"
		if (iv >= nb - 1 || ki2 == *n) {
#line 1259 "dtrevc3.f"
		    i__3 = *n - ki2 + iv;
#line 1259 "dtrevc3.f"
		    dgemm_("N", "N", n, &iv, &i__3, &c_b29, &vl[(ki2 - iv + 1)
			     * vl_dim1 + 1], ldvl, &work[ki2 - iv + 1 + *n], 
			    n, &c_b17, &work[(nb + 1) * *n + 1], n, (ftnlen)1,
			     (ftnlen)1);
/*                 normalize vectors */
#line 1265 "dtrevc3.f"
		    i__3 = iv;
#line 1265 "dtrevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1266 "dtrevc3.f"
			if (iscomplex[k - 1] == 0) {
/*                       real eigenvector */
#line 1268 "dtrevc3.f"
			    ii = idamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 1269 "dtrevc3.f"
			    remax = 1. / (d__1 = work[ii + (nb + k) * *n], 
				    abs(d__1));
#line 1270 "dtrevc3.f"
			} else if (iscomplex[k - 1] == 1) {
/*                       first eigenvector of conjugate pair */
#line 1272 "dtrevc3.f"
			    emax = 0.;
#line 1273 "dtrevc3.f"
			    i__4 = *n;
#line 1273 "dtrevc3.f"
			    for (ii = 1; ii <= i__4; ++ii) {
/* Computing MAX */
#line 1274 "dtrevc3.f"
				d__3 = emax, d__4 = (d__1 = work[ii + (nb + k)
					 * *n], abs(d__1)) + (d__2 = work[ii 
					+ (nb + k + 1) * *n], abs(d__2));
#line 1274 "dtrevc3.f"
				emax = max(d__3,d__4);
#line 1277 "dtrevc3.f"
			    }
#line 1278 "dtrevc3.f"
			    remax = 1. / emax;
/*                    else if ISCOMPLEX(K).EQ.-1 */
/*                       second eigenvector of conjugate pair */
/*                       reuse same REMAX as previous K */
#line 1282 "dtrevc3.f"
			}
#line 1283 "dtrevc3.f"
			dscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 1284 "dtrevc3.f"
		    }
#line 1285 "dtrevc3.f"
		    dlacpy_("F", n, &iv, &work[(nb + 1) * *n + 1], n, &vl[(
			    ki2 - iv + 1) * vl_dim1 + 1], ldvl, (ftnlen)1);
#line 1288 "dtrevc3.f"
		    iv = 1;
#line 1289 "dtrevc3.f"
		} else {
#line 1290 "dtrevc3.f"
		    ++iv;
#line 1291 "dtrevc3.f"
		}
#line 1292 "dtrevc3.f"
	    }

/* blocked back-transform */
#line 1294 "dtrevc3.f"
	    ++is;
#line 1295 "dtrevc3.f"
	    if (ip != 0) {
#line 1295 "dtrevc3.f"
		++is;
#line 1295 "dtrevc3.f"
	    }
#line 1297 "dtrevc3.f"
L260:
#line 1297 "dtrevc3.f"
	    ;
#line 1297 "dtrevc3.f"
	}
#line 1298 "dtrevc3.f"
    }

#line 1300 "dtrevc3.f"
    return 0;

/*     End of DTREVC3 */

} /* dtrevc3_ */

