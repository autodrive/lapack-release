#line 1 "strevc3.f"
/* strevc3.f -- translated by f2c (version 20100827).
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

#line 1 "strevc3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static logical c_false = FALSE_;
static doublereal c_b29 = 1.;
static logical c_true = TRUE_;

/* > \brief \b STREVC3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STREVC3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, */
/*                           VR, LDVR, MM, M, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STREVC3 computes some or all of the right and/or left eigenvectors of */
/* > a real upper quasi-triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR. */
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
/* >          The leading dimension of the array VL. */
/* >          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N. */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/*  @generated from dtrevc3.f, fortran d -> s, Tue Apr 19 01:47:44 2016 */

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
/* Subroutine */ int strevc3_(char *side, char *howmny, logical *select, 
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
	    integer *), sgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
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
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal bignum;
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
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

#line 296 "strevc3.f"
    /* Parameter adjustments */
#line 296 "strevc3.f"
    --select;
#line 296 "strevc3.f"
    t_dim1 = *ldt;
#line 296 "strevc3.f"
    t_offset = 1 + t_dim1;
#line 296 "strevc3.f"
    t -= t_offset;
#line 296 "strevc3.f"
    vl_dim1 = *ldvl;
#line 296 "strevc3.f"
    vl_offset = 1 + vl_dim1;
#line 296 "strevc3.f"
    vl -= vl_offset;
#line 296 "strevc3.f"
    vr_dim1 = *ldvr;
#line 296 "strevc3.f"
    vr_offset = 1 + vr_dim1;
#line 296 "strevc3.f"
    vr -= vr_offset;
#line 296 "strevc3.f"
    --work;
#line 296 "strevc3.f"

#line 296 "strevc3.f"
    /* Function Body */
#line 296 "strevc3.f"
    bothv = lsame_(side, "B", (ftnlen)1, (ftnlen)1);
#line 297 "strevc3.f"
    rightv = lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bothv;
#line 298 "strevc3.f"
    leftv = lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bothv;

#line 300 "strevc3.f"
    allv = lsame_(howmny, "A", (ftnlen)1, (ftnlen)1);
#line 301 "strevc3.f"
    over = lsame_(howmny, "B", (ftnlen)1, (ftnlen)1);
#line 302 "strevc3.f"
    somev = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 304 "strevc3.f"
    *info = 0;
/* Writing concatenation */
#line 305 "strevc3.f"
    i__1[0] = 1, a__1[0] = side;
#line 305 "strevc3.f"
    i__1[1] = 1, a__1[1] = howmny;
#line 305 "strevc3.f"
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 305 "strevc3.f"
    nb = ilaenv_(&c__1, "STREVC", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
#line 306 "strevc3.f"
    maxwrk = *n + (*n << 1) * nb;
#line 307 "strevc3.f"
    work[1] = (doublereal) maxwrk;
#line 308 "strevc3.f"
    lquery = *lwork == -1;
#line 309 "strevc3.f"
    if (! rightv && ! leftv) {
#line 310 "strevc3.f"
	*info = -1;
#line 311 "strevc3.f"
    } else if (! allv && ! over && ! somev) {
#line 312 "strevc3.f"
	*info = -2;
#line 313 "strevc3.f"
    } else if (*n < 0) {
#line 314 "strevc3.f"
	*info = -4;
#line 315 "strevc3.f"
    } else if (*ldt < max(1,*n)) {
#line 316 "strevc3.f"
	*info = -6;
#line 317 "strevc3.f"
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
#line 318 "strevc3.f"
	*info = -8;
#line 319 "strevc3.f"
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
#line 320 "strevc3.f"
	*info = -10;
#line 321 "strevc3.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 321 "strevc3.f"
	i__2 = 1, i__3 = *n * 3;
#line 321 "strevc3.f"
	if (*lwork < max(i__2,i__3) && ! lquery) {
#line 322 "strevc3.f"
	    *info = -14;
#line 323 "strevc3.f"
	} else {

/*        Set M to the number of columns required to store the selected */
/*        eigenvectors, standardize the array SELECT if necessary, and */
/*        test MM. */

#line 329 "strevc3.f"
	    if (somev) {
#line 330 "strevc3.f"
		*m = 0;
#line 331 "strevc3.f"
		pair = FALSE_;
#line 332 "strevc3.f"
		i__2 = *n;
#line 332 "strevc3.f"
		for (j = 1; j <= i__2; ++j) {
#line 333 "strevc3.f"
		    if (pair) {
#line 334 "strevc3.f"
			pair = FALSE_;
#line 335 "strevc3.f"
			select[j] = FALSE_;
#line 336 "strevc3.f"
		    } else {
#line 337 "strevc3.f"
			if (j < *n) {
#line 338 "strevc3.f"
			    if (t[j + 1 + j * t_dim1] == 0.) {
#line 339 "strevc3.f"
				if (select[j]) {
#line 339 "strevc3.f"
				    ++(*m);
#line 339 "strevc3.f"
				}
#line 341 "strevc3.f"
			    } else {
#line 342 "strevc3.f"
				pair = TRUE_;
#line 343 "strevc3.f"
				if (select[j] || select[j + 1]) {
#line 344 "strevc3.f"
				    select[j] = TRUE_;
#line 345 "strevc3.f"
				    *m += 2;
#line 346 "strevc3.f"
				}
#line 347 "strevc3.f"
			    }
#line 348 "strevc3.f"
			} else {
#line 349 "strevc3.f"
			    if (select[*n]) {
#line 349 "strevc3.f"
				++(*m);
#line 349 "strevc3.f"
			    }
#line 351 "strevc3.f"
			}
#line 352 "strevc3.f"
		    }
#line 353 "strevc3.f"
/* L10: */
#line 353 "strevc3.f"
		}
#line 354 "strevc3.f"
	    } else {
#line 355 "strevc3.f"
		*m = *n;
#line 356 "strevc3.f"
	    }

#line 358 "strevc3.f"
	    if (*mm < *m) {
#line 359 "strevc3.f"
		*info = -11;
#line 360 "strevc3.f"
	    }
#line 361 "strevc3.f"
	}
#line 361 "strevc3.f"
    }
#line 362 "strevc3.f"
    if (*info != 0) {
#line 363 "strevc3.f"
	i__2 = -(*info);
#line 363 "strevc3.f"
	xerbla_("STREVC3", &i__2, (ftnlen)7);
#line 364 "strevc3.f"
	return 0;
#line 365 "strevc3.f"
    } else if (lquery) {
#line 366 "strevc3.f"
	return 0;
#line 367 "strevc3.f"
    }

/*     Quick return if possible. */

#line 371 "strevc3.f"
    if (*n == 0) {
#line 371 "strevc3.f"
	return 0;
#line 371 "strevc3.f"
    }

/*     Use blocked version of back-transformation if sufficient workspace. */
/*     Zero-out the workspace to avoid potential NaN propagation. */

#line 377 "strevc3.f"
    if (over && *lwork >= *n + (*n << 4)) {
#line 378 "strevc3.f"
	nb = (*lwork - *n) / (*n << 1);
#line 379 "strevc3.f"
	nb = min(nb,128);
#line 380 "strevc3.f"
	i__2 = (nb << 1) + 1;
#line 380 "strevc3.f"
	slaset_("F", n, &i__2, &c_b17, &c_b17, &work[1], n, (ftnlen)1);
#line 381 "strevc3.f"
    } else {
#line 382 "strevc3.f"
	nb = 1;
#line 383 "strevc3.f"
    }

/*     Set the constants to control overflow. */

#line 387 "strevc3.f"
    unfl = slamch_("Safe minimum", (ftnlen)12);
#line 388 "strevc3.f"
    ovfl = 1. / unfl;
#line 389 "strevc3.f"
    slabad_(&unfl, &ovfl);
#line 390 "strevc3.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 391 "strevc3.f"
    smlnum = unfl * (*n / ulp);
#line 392 "strevc3.f"
    bignum = (1. - ulp) / smlnum;

/*     Compute 1-norm of each column of strictly upper triangular */
/*     part of T to control overflow in triangular solver. */

#line 397 "strevc3.f"
    work[1] = 0.;
#line 398 "strevc3.f"
    i__2 = *n;
#line 398 "strevc3.f"
    for (j = 2; j <= i__2; ++j) {
#line 399 "strevc3.f"
	work[j] = 0.;
#line 400 "strevc3.f"
	i__3 = j - 1;
#line 400 "strevc3.f"
	for (i__ = 1; i__ <= i__3; ++i__) {
#line 401 "strevc3.f"
	    work[j] += (d__1 = t[i__ + j * t_dim1], abs(d__1));
#line 402 "strevc3.f"
/* L20: */
#line 402 "strevc3.f"
	}
#line 403 "strevc3.f"
/* L30: */
#line 403 "strevc3.f"
    }

/*     Index IP is used to specify the real or complex eigenvalue: */
/*       IP = 0, real eigenvalue, */
/*            1, first  of conjugate complex pair: (wr,wi) */
/*           -1, second of conjugate complex pair: (wr,wi) */
/*       ISCOMPLEX array stores IP for each column in current block. */

#line 411 "strevc3.f"
    if (rightv) {

/*        ============================================================ */
/*        Compute right eigenvectors. */

/*        IV is index of column in current block. */
/*        For complex right vector, uses IV-1 for real part and IV for complex part. */
/*        Non-blocked version always uses IV=2; */
/*        blocked     version starts with IV=NB, goes down to 1 or 2. */
/*        (Note the "0-th" column is used for 1-norms computed above.) */
#line 421 "strevc3.f"
	iv = 2;
#line 422 "strevc3.f"
	if (nb > 2) {
#line 423 "strevc3.f"
	    iv = nb;
#line 424 "strevc3.f"
	}
#line 426 "strevc3.f"
	ip = 0;
#line 427 "strevc3.f"
	is = *m;
#line 428 "strevc3.f"
	for (ki = *n; ki >= 1; --ki) {
#line 429 "strevc3.f"
	    if (ip == -1) {
/*              previous iteration (ki+1) was second of conjugate pair, */
/*              so this ki is first of conjugate pair; skip to end of loop */
#line 432 "strevc3.f"
		ip = 1;
#line 433 "strevc3.f"
		goto L140;
#line 434 "strevc3.f"
	    } else if (ki == 1) {
/*              last column, so this ki must be real eigenvalue */
#line 436 "strevc3.f"
		ip = 0;
#line 437 "strevc3.f"
	    } else if (t[ki + (ki - 1) * t_dim1] == 0.) {
/*              zero on sub-diagonal, so this ki is real eigenvalue */
#line 439 "strevc3.f"
		ip = 0;
#line 440 "strevc3.f"
	    } else {
/*              non-zero on sub-diagonal, so this ki is second of conjugate pair */
#line 442 "strevc3.f"
		ip = -1;
#line 443 "strevc3.f"
	    }
#line 445 "strevc3.f"
	    if (somev) {
#line 446 "strevc3.f"
		if (ip == 0) {
#line 447 "strevc3.f"
		    if (! select[ki]) {
#line 447 "strevc3.f"
			goto L140;
#line 447 "strevc3.f"
		    }
#line 449 "strevc3.f"
		} else {
#line 450 "strevc3.f"
		    if (! select[ki - 1]) {
#line 450 "strevc3.f"
			goto L140;
#line 450 "strevc3.f"
		    }
#line 452 "strevc3.f"
		}
#line 453 "strevc3.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 457 "strevc3.f"
	    wr = t[ki + ki * t_dim1];
#line 458 "strevc3.f"
	    wi = 0.;
#line 459 "strevc3.f"
	    if (ip != 0) {
#line 459 "strevc3.f"
		wi = sqrt((d__1 = t[ki + (ki - 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki - 1 + ki * t_dim1], abs(d__2)));
#line 459 "strevc3.f"
	    }
/* Computing MAX */
#line 462 "strevc3.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 462 "strevc3.f"
	    smin = max(d__1,smlnum);

#line 464 "strevc3.f"
	    if (ip == 0) {

/*              -------------------------------------------------------- */
/*              Real right eigenvector */

#line 469 "strevc3.f"
		work[ki + iv * *n] = 1.;

/*              Form right-hand side. */

#line 473 "strevc3.f"
		i__2 = ki - 1;
#line 473 "strevc3.f"
		for (k = 1; k <= i__2; ++k) {
#line 474 "strevc3.f"
		    work[k + iv * *n] = -t[k + ki * t_dim1];
#line 475 "strevc3.f"
/* L50: */
#line 475 "strevc3.f"
		}

/*              Solve upper quasi-triangular system: */
/*              [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK. */

#line 480 "strevc3.f"
		jnxt = ki - 1;
#line 481 "strevc3.f"
		for (j = ki - 1; j >= 1; --j) {
#line 482 "strevc3.f"
		    if (j > jnxt) {
#line 482 "strevc3.f"
			goto L60;
#line 482 "strevc3.f"
		    }
#line 484 "strevc3.f"
		    j1 = j;
#line 485 "strevc3.f"
		    j2 = j;
#line 486 "strevc3.f"
		    jnxt = j - 1;
#line 487 "strevc3.f"
		    if (j > 1) {
#line 488 "strevc3.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 489 "strevc3.f"
			    j1 = j - 1;
#line 490 "strevc3.f"
			    jnxt = j - 2;
#line 491 "strevc3.f"
			}
#line 492 "strevc3.f"
		    }

#line 494 "strevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 498 "strevc3.f"
			slaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale X(1,1) to avoid overflow when updating */
/*                    the right-hand side. */

#line 505 "strevc3.f"
			if (xnorm > 1.) {
#line 506 "strevc3.f"
			    if (work[j] > bignum / xnorm) {
#line 507 "strevc3.f"
				x[0] /= xnorm;
#line 508 "strevc3.f"
				scale /= xnorm;
#line 509 "strevc3.f"
			    }
#line 510 "strevc3.f"
			}

/*                    Scale if necessary */

#line 514 "strevc3.f"
			if (scale != 1.) {
#line 514 "strevc3.f"
			    sscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 514 "strevc3.f"
			}
#line 516 "strevc3.f"
			work[j + iv * *n] = x[0];

/*                    Update right-hand side */

#line 520 "strevc3.f"
			i__2 = j - 1;
#line 520 "strevc3.f"
			d__1 = -x[0];
#line 520 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);

#line 523 "strevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 527 "strevc3.f"
			slaln2_(&c_false, &c__2, &c__1, &smin, &c_b29, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b29, &c_b29, &
				work[j - 1 + iv * *n], n, &wr, &c_b17, x, &
				c__2, &scale, &xnorm, &ierr);

/*                    Scale X(1,1) and X(2,1) to avoid overflow when */
/*                    updating the right-hand side. */

#line 535 "strevc3.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 536 "strevc3.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 536 "strevc3.f"
			    beta = max(d__1,d__2);
#line 537 "strevc3.f"
			    if (beta > bignum / xnorm) {
#line 538 "strevc3.f"
				x[0] /= xnorm;
#line 539 "strevc3.f"
				x[1] /= xnorm;
#line 540 "strevc3.f"
				scale /= xnorm;
#line 541 "strevc3.f"
			    }
#line 542 "strevc3.f"
			}

/*                    Scale if necessary */

#line 546 "strevc3.f"
			if (scale != 1.) {
#line 546 "strevc3.f"
			    sscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 546 "strevc3.f"
			}
#line 548 "strevc3.f"
			work[j - 1 + iv * *n] = x[0];
#line 549 "strevc3.f"
			work[j + iv * *n] = x[1];

/*                    Update right-hand side */

#line 553 "strevc3.f"
			i__2 = j - 2;
#line 553 "strevc3.f"
			d__1 = -x[0];
#line 553 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[iv * *n + 1], &c__1);
#line 555 "strevc3.f"
			i__2 = j - 2;
#line 555 "strevc3.f"
			d__1 = -x[1];
#line 555 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);
#line 557 "strevc3.f"
		    }
#line 558 "strevc3.f"
L60:
#line 558 "strevc3.f"
		    ;
#line 558 "strevc3.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 562 "strevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VR and normalize. */
#line 565 "strevc3.f"
		    scopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 
			    1], &c__1);

#line 567 "strevc3.f"
		    ii = isamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
#line 568 "strevc3.f"
		    remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs(d__1));
#line 569 "strevc3.f"
		    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 571 "strevc3.f"
		    i__2 = *n;
#line 571 "strevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 572 "strevc3.f"
			vr[k + is * vr_dim1] = 0.;
#line 573 "strevc3.f"
/* L70: */
#line 573 "strevc3.f"
		    }

#line 575 "strevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 578 "strevc3.f"
		    if (ki > 1) {
#line 578 "strevc3.f"
			i__2 = ki - 1;
#line 578 "strevc3.f"
			sgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[iv * *n + 1], &c__1, &work[ki + iv * *n],
				 &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 578 "strevc3.f"
		    }

#line 583 "strevc3.f"
		    ii = isamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
#line 584 "strevc3.f"
		    remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs(d__1));
#line 585 "strevc3.f"
		    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);

#line 587 "strevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out below vector */
#line 591 "strevc3.f"
		    i__2 = *n;
#line 591 "strevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 592 "strevc3.f"
			work[k + iv * *n] = 0.;
#line 593 "strevc3.f"
		    }
#line 594 "strevc3.f"
		    iscomplex[iv - 1] = ip;
/*                 back-transform and normalization is done below */
#line 596 "strevc3.f"
		}
#line 597 "strevc3.f"
	    } else {

/*              -------------------------------------------------------- */
/*              Complex right eigenvector. */

/*              Initial solve */
/*              [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0. */
/*              [ ( T(KI,  KI-1) T(KI,  KI) )               ] */

#line 606 "strevc3.f"
		if ((d__1 = t[ki - 1 + ki * t_dim1], abs(d__1)) >= (d__2 = t[
			ki + (ki - 1) * t_dim1], abs(d__2))) {
#line 607 "strevc3.f"
		    work[ki - 1 + (iv - 1) * *n] = 1.;
#line 608 "strevc3.f"
		    work[ki + iv * *n] = wi / t[ki - 1 + ki * t_dim1];
#line 609 "strevc3.f"
		} else {
#line 610 "strevc3.f"
		    work[ki - 1 + (iv - 1) * *n] = -wi / t[ki + (ki - 1) * 
			    t_dim1];
#line 611 "strevc3.f"
		    work[ki + iv * *n] = 1.;
#line 612 "strevc3.f"
		}
#line 613 "strevc3.f"
		work[ki + (iv - 1) * *n] = 0.;
#line 614 "strevc3.f"
		work[ki - 1 + iv * *n] = 0.;

/*              Form right-hand side. */

#line 618 "strevc3.f"
		i__2 = ki - 2;
#line 618 "strevc3.f"
		for (k = 1; k <= i__2; ++k) {
#line 619 "strevc3.f"
		    work[k + (iv - 1) * *n] = -work[ki - 1 + (iv - 1) * *n] * 
			    t[k + (ki - 1) * t_dim1];
#line 620 "strevc3.f"
		    work[k + iv * *n] = -work[ki + iv * *n] * t[k + ki * 
			    t_dim1];
#line 621 "strevc3.f"
/* L80: */
#line 621 "strevc3.f"
		}

/*              Solve upper quasi-triangular system: */
/*              [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2) */

#line 626 "strevc3.f"
		jnxt = ki - 2;
#line 627 "strevc3.f"
		for (j = ki - 2; j >= 1; --j) {
#line 628 "strevc3.f"
		    if (j > jnxt) {
#line 628 "strevc3.f"
			goto L90;
#line 628 "strevc3.f"
		    }
#line 630 "strevc3.f"
		    j1 = j;
#line 631 "strevc3.f"
		    j2 = j;
#line 632 "strevc3.f"
		    jnxt = j - 1;
#line 633 "strevc3.f"
		    if (j > 1) {
#line 634 "strevc3.f"
			if (t[j + (j - 1) * t_dim1] != 0.) {
#line 635 "strevc3.f"
			    j1 = j - 1;
#line 636 "strevc3.f"
			    jnxt = j - 2;
#line 637 "strevc3.f"
			}
#line 638 "strevc3.f"
		    }

#line 640 "strevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

#line 644 "strevc3.f"
			slaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + (
				iv - 1) * *n], n, &wr, &wi, x, &c__2, &scale, 
				&xnorm, &ierr);

/*                    Scale X(1,1) and X(1,2) to avoid overflow when */
/*                    updating the right-hand side. */

#line 651 "strevc3.f"
			if (xnorm > 1.) {
#line 652 "strevc3.f"
			    if (work[j] > bignum / xnorm) {
#line 653 "strevc3.f"
				x[0] /= xnorm;
#line 654 "strevc3.f"
				x[2] /= xnorm;
#line 655 "strevc3.f"
				scale /= xnorm;
#line 656 "strevc3.f"
			    }
#line 657 "strevc3.f"
			}

/*                    Scale if necessary */

#line 661 "strevc3.f"
			if (scale != 1.) {
#line 662 "strevc3.f"
			    sscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &
				    c__1);
#line 663 "strevc3.f"
			    sscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 664 "strevc3.f"
			}
#line 665 "strevc3.f"
			work[j + (iv - 1) * *n] = x[0];
#line 666 "strevc3.f"
			work[j + iv * *n] = x[2];

/*                    Update the right-hand side */

#line 670 "strevc3.f"
			i__2 = j - 1;
#line 670 "strevc3.f"
			d__1 = -x[0];
#line 670 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				(iv - 1) * *n + 1], &c__1);
#line 672 "strevc3.f"
			i__2 = j - 1;
#line 672 "strevc3.f"
			d__1 = -x[2];
#line 672 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);

#line 675 "strevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

#line 679 "strevc3.f"
			slaln2_(&c_false, &c__2, &c__2, &smin, &c_b29, &t[j - 
				1 + (j - 1) * t_dim1], ldt, &c_b29, &c_b29, &
				work[j - 1 + (iv - 1) * *n], n, &wr, &wi, x, &
				c__2, &scale, &xnorm, &ierr);

/*                    Scale X to avoid overflow when updating */
/*                    the right-hand side. */

#line 687 "strevc3.f"
			if (xnorm > 1.) {
/* Computing MAX */
#line 688 "strevc3.f"
			    d__1 = work[j - 1], d__2 = work[j];
#line 688 "strevc3.f"
			    beta = max(d__1,d__2);
#line 689 "strevc3.f"
			    if (beta > bignum / xnorm) {
#line 690 "strevc3.f"
				rec = 1. / xnorm;
#line 691 "strevc3.f"
				x[0] *= rec;
#line 692 "strevc3.f"
				x[2] *= rec;
#line 693 "strevc3.f"
				x[1] *= rec;
#line 694 "strevc3.f"
				x[3] *= rec;
#line 695 "strevc3.f"
				scale *= rec;
#line 696 "strevc3.f"
			    }
#line 697 "strevc3.f"
			}

/*                    Scale if necessary */

#line 701 "strevc3.f"
			if (scale != 1.) {
#line 702 "strevc3.f"
			    sscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &
				    c__1);
#line 703 "strevc3.f"
			    sscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
#line 704 "strevc3.f"
			}
#line 705 "strevc3.f"
			work[j - 1 + (iv - 1) * *n] = x[0];
#line 706 "strevc3.f"
			work[j + (iv - 1) * *n] = x[1];
#line 707 "strevc3.f"
			work[j - 1 + iv * *n] = x[2];
#line 708 "strevc3.f"
			work[j + iv * *n] = x[3];

/*                    Update the right-hand side */

#line 712 "strevc3.f"
			i__2 = j - 2;
#line 712 "strevc3.f"
			d__1 = -x[0];
#line 712 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[(iv - 1) * *n + 1], &c__1);
#line 714 "strevc3.f"
			i__2 = j - 2;
#line 714 "strevc3.f"
			d__1 = -x[1];
#line 714 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				(iv - 1) * *n + 1], &c__1);
#line 716 "strevc3.f"
			i__2 = j - 2;
#line 716 "strevc3.f"
			d__1 = -x[2];
#line 716 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, 
				&work[iv * *n + 1], &c__1);
#line 718 "strevc3.f"
			i__2 = j - 2;
#line 718 "strevc3.f"
			d__1 = -x[3];
#line 718 "strevc3.f"
			saxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[
				iv * *n + 1], &c__1);
#line 720 "strevc3.f"
		    }
#line 721 "strevc3.f"
L90:
#line 721 "strevc3.f"
		    ;
#line 721 "strevc3.f"
		}

/*              Copy the vector x or Q*x to VR and normalize. */

#line 725 "strevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VR and normalize. */
#line 728 "strevc3.f"
		    scopy_(&ki, &work[(iv - 1) * *n + 1], &c__1, &vr[(is - 1) 
			    * vr_dim1 + 1], &c__1);
#line 729 "strevc3.f"
		    scopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 
			    1], &c__1);

#line 731 "strevc3.f"
		    emax = 0.;
#line 732 "strevc3.f"
		    i__2 = ki;
#line 732 "strevc3.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 733 "strevc3.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + is * vr_dim1], 
				abs(d__2));
#line 733 "strevc3.f"
			emax = max(d__3,d__4);
#line 735 "strevc3.f"
/* L100: */
#line 735 "strevc3.f"
		    }
#line 736 "strevc3.f"
		    remax = 1. / emax;
#line 737 "strevc3.f"
		    sscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
#line 738 "strevc3.f"
		    sscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

#line 740 "strevc3.f"
		    i__2 = *n;
#line 740 "strevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 741 "strevc3.f"
			vr[k + (is - 1) * vr_dim1] = 0.;
#line 742 "strevc3.f"
			vr[k + is * vr_dim1] = 0.;
#line 743 "strevc3.f"
/* L110: */
#line 743 "strevc3.f"
		    }

#line 745 "strevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 748 "strevc3.f"
		    if (ki > 2) {
#line 749 "strevc3.f"
			i__2 = ki - 2;
#line 749 "strevc3.f"
			sgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[(iv - 1) * *n + 1], &c__1, &work[ki - 1 
				+ (iv - 1) * *n], &vr[(ki - 1) * vr_dim1 + 1],
				 &c__1, (ftnlen)1);
#line 752 "strevc3.f"
			i__2 = ki - 2;
#line 752 "strevc3.f"
			sgemv_("N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &
				work[iv * *n + 1], &c__1, &work[ki + iv * *n],
				 &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
#line 755 "strevc3.f"
		    } else {
#line 756 "strevc3.f"
			sscal_(n, &work[ki - 1 + (iv - 1) * *n], &vr[(ki - 1) 
				* vr_dim1 + 1], &c__1);
#line 757 "strevc3.f"
			sscal_(n, &work[ki + iv * *n], &vr[ki * vr_dim1 + 1], 
				&c__1);
#line 758 "strevc3.f"
		    }

#line 760 "strevc3.f"
		    emax = 0.;
#line 761 "strevc3.f"
		    i__2 = *n;
#line 761 "strevc3.f"
		    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 762 "strevc3.f"
			d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1]
				, abs(d__1)) + (d__2 = vr[k + ki * vr_dim1], 
				abs(d__2));
#line 762 "strevc3.f"
			emax = max(d__3,d__4);
#line 764 "strevc3.f"
/* L120: */
#line 764 "strevc3.f"
		    }
#line 765 "strevc3.f"
		    remax = 1. / emax;
#line 766 "strevc3.f"
		    sscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
#line 767 "strevc3.f"
		    sscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);

#line 769 "strevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out below vector */
#line 773 "strevc3.f"
		    i__2 = *n;
#line 773 "strevc3.f"
		    for (k = ki + 1; k <= i__2; ++k) {
#line 774 "strevc3.f"
			work[k + (iv - 1) * *n] = 0.;
#line 775 "strevc3.f"
			work[k + iv * *n] = 0.;
#line 776 "strevc3.f"
		    }
#line 777 "strevc3.f"
		    iscomplex[iv - 2] = -ip;
#line 778 "strevc3.f"
		    iscomplex[iv - 1] = ip;
#line 779 "strevc3.f"
		    --iv;
/*                 back-transform and normalization is done below */
#line 781 "strevc3.f"
		}
#line 782 "strevc3.f"
	    }
#line 784 "strevc3.f"
	    if (nb > 1) {
/*              -------------------------------------------------------- */
/*              Blocked version of back-transform */
/*              For complex case, KI2 includes both vectors (KI-1 and KI) */
#line 788 "strevc3.f"
		if (ip == 0) {
#line 789 "strevc3.f"
		    ki2 = ki;
#line 790 "strevc3.f"
		} else {
#line 791 "strevc3.f"
		    ki2 = ki - 1;
#line 792 "strevc3.f"
		}
/*              Columns IV:NB of work are valid vectors. */
/*              When the number of vectors stored reaches NB-1 or NB, */
/*              or if this was last vector, do the GEMM */
#line 797 "strevc3.f"
		if (iv <= 2 || ki2 == 1) {
#line 798 "strevc3.f"
		    i__2 = nb - iv + 1;
#line 798 "strevc3.f"
		    i__3 = ki2 + nb - iv;
#line 798 "strevc3.f"
		    sgemm_("N", "N", n, &i__2, &i__3, &c_b29, &vr[vr_offset], 
			    ldvr, &work[iv * *n + 1], n, &c_b17, &work[(nb + 
			    iv) * *n + 1], n, (ftnlen)1, (ftnlen)1);
/*                 normalize vectors */
#line 804 "strevc3.f"
		    i__2 = nb;
#line 804 "strevc3.f"
		    for (k = iv; k <= i__2; ++k) {
#line 805 "strevc3.f"
			if (iscomplex[k - 1] == 0) {
/*                       real eigenvector */
#line 807 "strevc3.f"
			    ii = isamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 808 "strevc3.f"
			    remax = 1. / (d__1 = work[ii + (nb + k) * *n], 
				    abs(d__1));
#line 809 "strevc3.f"
			} else if (iscomplex[k - 1] == 1) {
/*                       first eigenvector of conjugate pair */
#line 811 "strevc3.f"
			    emax = 0.;
#line 812 "strevc3.f"
			    i__3 = *n;
#line 812 "strevc3.f"
			    for (ii = 1; ii <= i__3; ++ii) {
/* Computing MAX */
#line 813 "strevc3.f"
				d__3 = emax, d__4 = (d__1 = work[ii + (nb + k)
					 * *n], abs(d__1)) + (d__2 = work[ii 
					+ (nb + k + 1) * *n], abs(d__2));
#line 813 "strevc3.f"
				emax = max(d__3,d__4);
#line 816 "strevc3.f"
			    }
#line 817 "strevc3.f"
			    remax = 1. / emax;
/*                    else if ISCOMPLEX(K).EQ.-1 */
/*                       second eigenvector of conjugate pair */
/*                       reuse same REMAX as previous K */
#line 821 "strevc3.f"
			}
#line 822 "strevc3.f"
			sscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 823 "strevc3.f"
		    }
#line 824 "strevc3.f"
		    i__2 = nb - iv + 1;
#line 824 "strevc3.f"
		    slacpy_("F", n, &i__2, &work[(nb + iv) * *n + 1], n, &vr[
			    ki2 * vr_dim1 + 1], ldvr, (ftnlen)1);
#line 827 "strevc3.f"
		    iv = nb;
#line 828 "strevc3.f"
		} else {
#line 829 "strevc3.f"
		    --iv;
#line 830 "strevc3.f"
		}
#line 831 "strevc3.f"
	    }

/* blocked back-transform */
#line 833 "strevc3.f"
	    --is;
#line 834 "strevc3.f"
	    if (ip != 0) {
#line 834 "strevc3.f"
		--is;
#line 834 "strevc3.f"
	    }
#line 836 "strevc3.f"
L140:
#line 836 "strevc3.f"
	    ;
#line 836 "strevc3.f"
	}
#line 837 "strevc3.f"
    }
#line 839 "strevc3.f"
    if (leftv) {

/*        ============================================================ */
/*        Compute left eigenvectors. */

/*        IV is index of column in current block. */
/*        For complex left vector, uses IV for real part and IV+1 for complex part. */
/*        Non-blocked version always uses IV=1; */
/*        blocked     version starts with IV=1, goes up to NB-1 or NB. */
/*        (Note the "0-th" column is used for 1-norms computed above.) */
#line 849 "strevc3.f"
	iv = 1;
#line 850 "strevc3.f"
	ip = 0;
#line 851 "strevc3.f"
	is = 1;
#line 852 "strevc3.f"
	i__2 = *n;
#line 852 "strevc3.f"
	for (ki = 1; ki <= i__2; ++ki) {
#line 853 "strevc3.f"
	    if (ip == 1) {
/*              previous iteration (ki-1) was first of conjugate pair, */
/*              so this ki is second of conjugate pair; skip to end of loop */
#line 856 "strevc3.f"
		ip = -1;
#line 857 "strevc3.f"
		goto L260;
#line 858 "strevc3.f"
	    } else if (ki == *n) {
/*              last column, so this ki must be real eigenvalue */
#line 860 "strevc3.f"
		ip = 0;
#line 861 "strevc3.f"
	    } else if (t[ki + 1 + ki * t_dim1] == 0.) {
/*              zero on sub-diagonal, so this ki is real eigenvalue */
#line 863 "strevc3.f"
		ip = 0;
#line 864 "strevc3.f"
	    } else {
/*              non-zero on sub-diagonal, so this ki is first of conjugate pair */
#line 866 "strevc3.f"
		ip = 1;
#line 867 "strevc3.f"
	    }

#line 869 "strevc3.f"
	    if (somev) {
#line 870 "strevc3.f"
		if (! select[ki]) {
#line 870 "strevc3.f"
		    goto L260;
#line 870 "strevc3.f"
		}
#line 872 "strevc3.f"
	    }

/*           Compute the KI-th eigenvalue (WR,WI). */

#line 876 "strevc3.f"
	    wr = t[ki + ki * t_dim1];
#line 877 "strevc3.f"
	    wi = 0.;
#line 878 "strevc3.f"
	    if (ip != 0) {
#line 878 "strevc3.f"
		wi = sqrt((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1))) * 
			sqrt((d__2 = t[ki + 1 + ki * t_dim1], abs(d__2)));
#line 878 "strevc3.f"
	    }
/* Computing MAX */
#line 881 "strevc3.f"
	    d__1 = ulp * (abs(wr) + abs(wi));
#line 881 "strevc3.f"
	    smin = max(d__1,smlnum);

#line 883 "strevc3.f"
	    if (ip == 0) {

/*              -------------------------------------------------------- */
/*              Real left eigenvector */

#line 888 "strevc3.f"
		work[ki + iv * *n] = 1.;

/*              Form right-hand side. */

#line 892 "strevc3.f"
		i__3 = *n;
#line 892 "strevc3.f"
		for (k = ki + 1; k <= i__3; ++k) {
#line 893 "strevc3.f"
		    work[k + iv * *n] = -t[ki + k * t_dim1];
#line 894 "strevc3.f"
/* L160: */
#line 894 "strevc3.f"
		}

/*              Solve transposed quasi-triangular system: */
/*              [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK */

#line 899 "strevc3.f"
		vmax = 1.;
#line 900 "strevc3.f"
		vcrit = bignum;

#line 902 "strevc3.f"
		jnxt = ki + 1;
#line 903 "strevc3.f"
		i__3 = *n;
#line 903 "strevc3.f"
		for (j = ki + 1; j <= i__3; ++j) {
#line 904 "strevc3.f"
		    if (j < jnxt) {
#line 904 "strevc3.f"
			goto L170;
#line 904 "strevc3.f"
		    }
#line 906 "strevc3.f"
		    j1 = j;
#line 907 "strevc3.f"
		    j2 = j;
#line 908 "strevc3.f"
		    jnxt = j + 1;
#line 909 "strevc3.f"
		    if (j < *n) {
#line 910 "strevc3.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 911 "strevc3.f"
			    j2 = j + 1;
#line 912 "strevc3.f"
			    jnxt = j + 2;
#line 913 "strevc3.f"
			}
#line 914 "strevc3.f"
		    }

#line 916 "strevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

#line 923 "strevc3.f"
			if (work[j] > vcrit) {
#line 924 "strevc3.f"
			    rec = 1. / vmax;
#line 925 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 925 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 926 "strevc3.f"
			    vmax = 1.;
#line 927 "strevc3.f"
			    vcrit = bignum;
#line 928 "strevc3.f"
			}

#line 930 "strevc3.f"
			i__4 = j - ki - 1;
#line 930 "strevc3.f"
			work[j + iv * *n] -= sdot_(&i__4, &t[ki + 1 + j * 
				t_dim1], &c__1, &work[ki + 1 + iv * *n], &
				c__1);

/*                    Solve [ T(J,J) - WR ]**T * X = WORK */

#line 936 "strevc3.f"
			slaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 942 "strevc3.f"
			if (scale != 1.) {
#line 942 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 942 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 942 "strevc3.f"
			}
#line 944 "strevc3.f"
			work[j + iv * *n] = x[0];
/* Computing MAX */
#line 945 "strevc3.f"
			d__2 = (d__1 = work[j + iv * *n], abs(d__1));
#line 945 "strevc3.f"
			vmax = max(d__2,vmax);
#line 946 "strevc3.f"
			vcrit = bignum / vmax;

#line 948 "strevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side. */

/* Computing MAX */
#line 955 "strevc3.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 955 "strevc3.f"
			beta = max(d__1,d__2);
#line 956 "strevc3.f"
			if (beta > vcrit) {
#line 957 "strevc3.f"
			    rec = 1. / vmax;
#line 958 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 958 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 959 "strevc3.f"
			    vmax = 1.;
#line 960 "strevc3.f"
			    vcrit = bignum;
#line 961 "strevc3.f"
			}

#line 963 "strevc3.f"
			i__4 = j - ki - 1;
#line 963 "strevc3.f"
			work[j + iv * *n] -= sdot_(&i__4, &t[ki + 1 + j * 
				t_dim1], &c__1, &work[ki + 1 + iv * *n], &
				c__1);

#line 967 "strevc3.f"
			i__4 = j - ki - 1;
#line 967 "strevc3.f"
			work[j + 1 + iv * *n] -= sdot_(&i__4, &t[ki + 1 + (j 
				+ 1) * t_dim1], &c__1, &work[ki + 1 + iv * *n]
				, &c__1);

/*                    Solve */
/*                    [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 ) */
/*                    [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 ) */

#line 975 "strevc3.f"
			slaln2_(&c_true, &c__2, &c__1, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &c_b17, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 981 "strevc3.f"
			if (scale != 1.) {
#line 981 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 981 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 981 "strevc3.f"
			}
#line 983 "strevc3.f"
			work[j + iv * *n] = x[0];
#line 984 "strevc3.f"
			work[j + 1 + iv * *n] = x[1];

/* Computing MAX */
#line 986 "strevc3.f"
			d__3 = (d__1 = work[j + iv * *n], abs(d__1)), d__4 = (
				d__2 = work[j + 1 + iv * *n], abs(d__2)), 
				d__3 = max(d__3,d__4);
#line 986 "strevc3.f"
			vmax = max(d__3,vmax);
#line 988 "strevc3.f"
			vcrit = bignum / vmax;

#line 990 "strevc3.f"
		    }
#line 991 "strevc3.f"
L170:
#line 991 "strevc3.f"
		    ;
#line 991 "strevc3.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 995 "strevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VL and normalize. */
#line 998 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 998 "strevc3.f"
		    scopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);

#line 1001 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1001 "strevc3.f"
		    ii = isamax_(&i__3, &vl[ki + is * vl_dim1], &c__1) + ki - 
			    1;
#line 1002 "strevc3.f"
		    remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs(d__1));
#line 1003 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1003 "strevc3.f"
		    sscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);

#line 1005 "strevc3.f"
		    i__3 = ki - 1;
#line 1005 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1006 "strevc3.f"
			vl[k + is * vl_dim1] = 0.;
#line 1007 "strevc3.f"
/* L180: */
#line 1007 "strevc3.f"
		    }

#line 1009 "strevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 1012 "strevc3.f"
		    if (ki < *n) {
#line 1012 "strevc3.f"
			i__3 = *n - ki;
#line 1012 "strevc3.f"
			sgemv_("N", n, &i__3, &c_b29, &vl[(ki + 1) * vl_dim1 
				+ 1], ldvl, &work[ki + 1 + iv * *n], &c__1, &
				work[ki + iv * *n], &vl[ki * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1012 "strevc3.f"
		    }

#line 1018 "strevc3.f"
		    ii = isamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
#line 1019 "strevc3.f"
		    remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs(d__1));
#line 1020 "strevc3.f"
		    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

#line 1022 "strevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out above vector */
/*                 could go from KI-NV+1 to KI-1 */
#line 1027 "strevc3.f"
		    i__3 = ki - 1;
#line 1027 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1028 "strevc3.f"
			work[k + iv * *n] = 0.;
#line 1029 "strevc3.f"
		    }
#line 1030 "strevc3.f"
		    iscomplex[iv - 1] = ip;
/*                 back-transform and normalization is done below */
#line 1032 "strevc3.f"
		}
#line 1033 "strevc3.f"
	    } else {

/*              -------------------------------------------------------- */
/*              Complex left eigenvector. */

/*              Initial solve: */
/*              [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0. */
/*              [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ] */

#line 1042 "strevc3.f"
		if ((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1)) >= (d__2 = 
			t[ki + 1 + ki * t_dim1], abs(d__2))) {
#line 1043 "strevc3.f"
		    work[ki + iv * *n] = wi / t[ki + (ki + 1) * t_dim1];
#line 1044 "strevc3.f"
		    work[ki + 1 + (iv + 1) * *n] = 1.;
#line 1045 "strevc3.f"
		} else {
#line 1046 "strevc3.f"
		    work[ki + iv * *n] = 1.;
#line 1047 "strevc3.f"
		    work[ki + 1 + (iv + 1) * *n] = -wi / t[ki + 1 + ki * 
			    t_dim1];
#line 1048 "strevc3.f"
		}
#line 1049 "strevc3.f"
		work[ki + 1 + iv * *n] = 0.;
#line 1050 "strevc3.f"
		work[ki + (iv + 1) * *n] = 0.;

/*              Form right-hand side. */

#line 1054 "strevc3.f"
		i__3 = *n;
#line 1054 "strevc3.f"
		for (k = ki + 2; k <= i__3; ++k) {
#line 1055 "strevc3.f"
		    work[k + iv * *n] = -work[ki + iv * *n] * t[ki + k * 
			    t_dim1];
#line 1056 "strevc3.f"
		    work[k + (iv + 1) * *n] = -work[ki + 1 + (iv + 1) * *n] * 
			    t[ki + 1 + k * t_dim1];
#line 1057 "strevc3.f"
/* L190: */
#line 1057 "strevc3.f"
		}

/*              Solve transposed quasi-triangular system: */
/*              [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2 */

#line 1062 "strevc3.f"
		vmax = 1.;
#line 1063 "strevc3.f"
		vcrit = bignum;

#line 1065 "strevc3.f"
		jnxt = ki + 2;
#line 1066 "strevc3.f"
		i__3 = *n;
#line 1066 "strevc3.f"
		for (j = ki + 2; j <= i__3; ++j) {
#line 1067 "strevc3.f"
		    if (j < jnxt) {
#line 1067 "strevc3.f"
			goto L200;
#line 1067 "strevc3.f"
		    }
#line 1069 "strevc3.f"
		    j1 = j;
#line 1070 "strevc3.f"
		    j2 = j;
#line 1071 "strevc3.f"
		    jnxt = j + 1;
#line 1072 "strevc3.f"
		    if (j < *n) {
#line 1073 "strevc3.f"
			if (t[j + 1 + j * t_dim1] != 0.) {
#line 1074 "strevc3.f"
			    j2 = j + 1;
#line 1075 "strevc3.f"
			    jnxt = j + 2;
#line 1076 "strevc3.f"
			}
#line 1077 "strevc3.f"
		    }

#line 1079 "strevc3.f"
		    if (j1 == j2) {

/*                    1-by-1 diagonal block */

/*                    Scale if necessary to avoid overflow when */
/*                    forming the right-hand side elements. */

#line 1086 "strevc3.f"
			if (work[j] > vcrit) {
#line 1087 "strevc3.f"
			    rec = 1. / vmax;
#line 1088 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1088 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 1089 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1089 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1090 "strevc3.f"
			    vmax = 1.;
#line 1091 "strevc3.f"
			    vcrit = bignum;
#line 1092 "strevc3.f"
			}

#line 1094 "strevc3.f"
			i__4 = j - ki - 2;
#line 1094 "strevc3.f"
			work[j + iv * *n] -= sdot_(&i__4, &t[ki + 2 + j * 
				t_dim1], &c__1, &work[ki + 2 + iv * *n], &
				c__1);
#line 1097 "strevc3.f"
			i__4 = j - ki - 2;
#line 1097 "strevc3.f"
			work[j + (iv + 1) * *n] -= sdot_(&i__4, &t[ki + 2 + j 
				* t_dim1], &c__1, &work[ki + 2 + (iv + 1) * *
				n], &c__1);

/*                    Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2 */

#line 1103 "strevc3.f"
			d__1 = -wi;
#line 1103 "strevc3.f"
			slaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &d__1, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 1109 "strevc3.f"
			if (scale != 1.) {
#line 1110 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1110 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 1111 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1111 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1112 "strevc3.f"
			}
#line 1113 "strevc3.f"
			work[j + iv * *n] = x[0];
#line 1114 "strevc3.f"
			work[j + (iv + 1) * *n] = x[2];
/* Computing MAX */
#line 1115 "strevc3.f"
			d__3 = (d__1 = work[j + iv * *n], abs(d__1)), d__4 = (
				d__2 = work[j + (iv + 1) * *n], abs(d__2)), 
				d__3 = max(d__3,d__4);
#line 1115 "strevc3.f"
			vmax = max(d__3,vmax);
#line 1117 "strevc3.f"
			vcrit = bignum / vmax;

#line 1119 "strevc3.f"
		    } else {

/*                    2-by-2 diagonal block */

/*                    Scale if necessary to avoid overflow when forming */
/*                    the right-hand side elements. */

/* Computing MAX */
#line 1126 "strevc3.f"
			d__1 = work[j], d__2 = work[j + 1];
#line 1126 "strevc3.f"
			beta = max(d__1,d__2);
#line 1127 "strevc3.f"
			if (beta > vcrit) {
#line 1128 "strevc3.f"
			    rec = 1. / vmax;
#line 1129 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1129 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
#line 1130 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1130 "strevc3.f"
			    sscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1131 "strevc3.f"
			    vmax = 1.;
#line 1132 "strevc3.f"
			    vcrit = bignum;
#line 1133 "strevc3.f"
			}

#line 1135 "strevc3.f"
			i__4 = j - ki - 2;
#line 1135 "strevc3.f"
			work[j + iv * *n] -= sdot_(&i__4, &t[ki + 2 + j * 
				t_dim1], &c__1, &work[ki + 2 + iv * *n], &
				c__1);

#line 1139 "strevc3.f"
			i__4 = j - ki - 2;
#line 1139 "strevc3.f"
			work[j + (iv + 1) * *n] -= sdot_(&i__4, &t[ki + 2 + j 
				* t_dim1], &c__1, &work[ki + 2 + (iv + 1) * *
				n], &c__1);

#line 1143 "strevc3.f"
			i__4 = j - ki - 2;
#line 1143 "strevc3.f"
			work[j + 1 + iv * *n] -= sdot_(&i__4, &t[ki + 2 + (j 
				+ 1) * t_dim1], &c__1, &work[ki + 2 + iv * *n]
				, &c__1);

#line 1147 "strevc3.f"
			i__4 = j - ki - 2;
#line 1147 "strevc3.f"
			work[j + 1 + (iv + 1) * *n] -= sdot_(&i__4, &t[ki + 2 
				+ (j + 1) * t_dim1], &c__1, &work[ki + 2 + (
				iv + 1) * *n], &c__1);

/*                    Solve 2-by-2 complex linear equation */
/*                    [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B */
/*                    [ (T(j+1,j) T(j+1,j+1))                  ] */

#line 1155 "strevc3.f"
			d__1 = -wi;
#line 1155 "strevc3.f"
			slaln2_(&c_true, &c__2, &c__2, &smin, &c_b29, &t[j + 
				j * t_dim1], ldt, &c_b29, &c_b29, &work[j + 
				iv * *n], n, &wr, &d__1, x, &c__2, &scale, &
				xnorm, &ierr);

/*                    Scale if necessary */

#line 1161 "strevc3.f"
			if (scale != 1.) {
#line 1162 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1162 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
#line 1163 "strevc3.f"
			    i__4 = *n - ki + 1;
#line 1163 "strevc3.f"
			    sscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &
				    c__1);
#line 1164 "strevc3.f"
			}
#line 1165 "strevc3.f"
			work[j + iv * *n] = x[0];
#line 1166 "strevc3.f"
			work[j + (iv + 1) * *n] = x[2];
#line 1167 "strevc3.f"
			work[j + 1 + iv * *n] = x[1];
#line 1168 "strevc3.f"
			work[j + 1 + (iv + 1) * *n] = x[3];
/* Computing MAX */
#line 1169 "strevc3.f"
			d__1 = abs(x[0]), d__2 = abs(x[2]), d__1 = max(d__1,
				d__2), d__2 = abs(x[1]), d__1 = max(d__1,d__2)
				, d__2 = abs(x[3]), d__1 = max(d__1,d__2);
#line 1169 "strevc3.f"
			vmax = max(d__1,vmax);
#line 1172 "strevc3.f"
			vcrit = bignum / vmax;

#line 1174 "strevc3.f"
		    }
#line 1175 "strevc3.f"
L200:
#line 1175 "strevc3.f"
		    ;
#line 1175 "strevc3.f"
		}

/*              Copy the vector x or Q*x to VL and normalize. */

#line 1179 "strevc3.f"
		if (! over) {
/*                 ------------------------------ */
/*                 no back-transform: copy x to VL and normalize. */
#line 1182 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1182 "strevc3.f"
		    scopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * 
			    vl_dim1], &c__1);
#line 1184 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1184 "strevc3.f"
		    scopy_(&i__3, &work[ki + (iv + 1) * *n], &c__1, &vl[ki + (
			    is + 1) * vl_dim1], &c__1);

#line 1187 "strevc3.f"
		    emax = 0.;
#line 1188 "strevc3.f"
		    i__3 = *n;
#line 1188 "strevc3.f"
		    for (k = ki; k <= i__3; ++k) {
/* Computing MAX */
#line 1189 "strevc3.f"
			d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (is + 1) * vl_dim1], 
				abs(d__2));
#line 1189 "strevc3.f"
			emax = max(d__3,d__4);
#line 1191 "strevc3.f"
/* L220: */
#line 1191 "strevc3.f"
		    }
#line 1192 "strevc3.f"
		    remax = 1. / emax;
#line 1193 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1193 "strevc3.f"
		    sscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);
#line 1194 "strevc3.f"
		    i__3 = *n - ki + 1;
#line 1194 "strevc3.f"
		    sscal_(&i__3, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1)
			    ;

#line 1196 "strevc3.f"
		    i__3 = ki - 1;
#line 1196 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1197 "strevc3.f"
			vl[k + is * vl_dim1] = 0.;
#line 1198 "strevc3.f"
			vl[k + (is + 1) * vl_dim1] = 0.;
#line 1199 "strevc3.f"
/* L230: */
#line 1199 "strevc3.f"
		    }

#line 1201 "strevc3.f"
		} else if (nb == 1) {
/*                 ------------------------------ */
/*                 version 1: back-transform each vector with GEMV, Q*x. */
#line 1204 "strevc3.f"
		    if (ki < *n - 1) {
#line 1205 "strevc3.f"
			i__3 = *n - ki - 1;
#line 1205 "strevc3.f"
			sgemv_("N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + iv * *n], &c__1, &
				work[ki + iv * *n], &vl[ki * vl_dim1 + 1], &
				c__1, (ftnlen)1);
#line 1210 "strevc3.f"
			i__3 = *n - ki - 1;
#line 1210 "strevc3.f"
			sgemv_("N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 
				+ 1], ldvl, &work[ki + 2 + (iv + 1) * *n], &
				c__1, &work[ki + 1 + (iv + 1) * *n], &vl[(ki 
				+ 1) * vl_dim1 + 1], &c__1, (ftnlen)1);
#line 1215 "strevc3.f"
		    } else {
#line 1216 "strevc3.f"
			sscal_(n, &work[ki + iv * *n], &vl[ki * vl_dim1 + 1], 
				&c__1);
#line 1217 "strevc3.f"
			sscal_(n, &work[ki + 1 + (iv + 1) * *n], &vl[(ki + 1) 
				* vl_dim1 + 1], &c__1);
#line 1218 "strevc3.f"
		    }

#line 1220 "strevc3.f"
		    emax = 0.;
#line 1221 "strevc3.f"
		    i__3 = *n;
#line 1221 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1222 "strevc3.f"
			d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs(
				d__1)) + (d__2 = vl[k + (ki + 1) * vl_dim1], 
				abs(d__2));
#line 1222 "strevc3.f"
			emax = max(d__3,d__4);
#line 1224 "strevc3.f"
/* L240: */
#line 1224 "strevc3.f"
		    }
#line 1225 "strevc3.f"
		    remax = 1. / emax;
#line 1226 "strevc3.f"
		    sscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
#line 1227 "strevc3.f"
		    sscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);

#line 1229 "strevc3.f"
		} else {
/*                 ------------------------------ */
/*                 version 2: back-transform block of vectors with GEMM */
/*                 zero out above vector */
/*                 could go from KI-NV+1 to KI-1 */
#line 1234 "strevc3.f"
		    i__3 = ki - 1;
#line 1234 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1235 "strevc3.f"
			work[k + iv * *n] = 0.;
#line 1236 "strevc3.f"
			work[k + (iv + 1) * *n] = 0.;
#line 1237 "strevc3.f"
		    }
#line 1238 "strevc3.f"
		    iscomplex[iv - 1] = ip;
#line 1239 "strevc3.f"
		    iscomplex[iv] = -ip;
#line 1240 "strevc3.f"
		    ++iv;
/*                 back-transform and normalization is done below */
#line 1242 "strevc3.f"
		}
#line 1243 "strevc3.f"
	    }
#line 1245 "strevc3.f"
	    if (nb > 1) {
/*              -------------------------------------------------------- */
/*              Blocked version of back-transform */
/*              For complex case, KI2 includes both vectors (KI and KI+1) */
#line 1249 "strevc3.f"
		if (ip == 0) {
#line 1250 "strevc3.f"
		    ki2 = ki;
#line 1251 "strevc3.f"
		} else {
#line 1252 "strevc3.f"
		    ki2 = ki + 1;
#line 1253 "strevc3.f"
		}
/*              Columns 1:IV of work are valid vectors. */
/*              When the number of vectors stored reaches NB-1 or NB, */
/*              or if this was last vector, do the GEMM */
#line 1258 "strevc3.f"
		if (iv >= nb - 1 || ki2 == *n) {
#line 1259 "strevc3.f"
		    i__3 = *n - ki2 + iv;
#line 1259 "strevc3.f"
		    sgemm_("N", "N", n, &iv, &i__3, &c_b29, &vl[(ki2 - iv + 1)
			     * vl_dim1 + 1], ldvl, &work[ki2 - iv + 1 + *n], 
			    n, &c_b17, &work[(nb + 1) * *n + 1], n, (ftnlen)1,
			     (ftnlen)1);
/*                 normalize vectors */
#line 1265 "strevc3.f"
		    i__3 = iv;
#line 1265 "strevc3.f"
		    for (k = 1; k <= i__3; ++k) {
#line 1266 "strevc3.f"
			if (iscomplex[k - 1] == 0) {
/*                       real eigenvector */
#line 1268 "strevc3.f"
			    ii = isamax_(n, &work[(nb + k) * *n + 1], &c__1);
#line 1269 "strevc3.f"
			    remax = 1. / (d__1 = work[ii + (nb + k) * *n], 
				    abs(d__1));
#line 1270 "strevc3.f"
			} else if (iscomplex[k - 1] == 1) {
/*                       first eigenvector of conjugate pair */
#line 1272 "strevc3.f"
			    emax = 0.;
#line 1273 "strevc3.f"
			    i__4 = *n;
#line 1273 "strevc3.f"
			    for (ii = 1; ii <= i__4; ++ii) {
/* Computing MAX */
#line 1274 "strevc3.f"
				d__3 = emax, d__4 = (d__1 = work[ii + (nb + k)
					 * *n], abs(d__1)) + (d__2 = work[ii 
					+ (nb + k + 1) * *n], abs(d__2));
#line 1274 "strevc3.f"
				emax = max(d__3,d__4);
#line 1277 "strevc3.f"
			    }
#line 1278 "strevc3.f"
			    remax = 1. / emax;
/*                    else if ISCOMPLEX(K).EQ.-1 */
/*                       second eigenvector of conjugate pair */
/*                       reuse same REMAX as previous K */
#line 1282 "strevc3.f"
			}
#line 1283 "strevc3.f"
			sscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
#line 1284 "strevc3.f"
		    }
#line 1285 "strevc3.f"
		    slacpy_("F", n, &iv, &work[(nb + 1) * *n + 1], n, &vl[(
			    ki2 - iv + 1) * vl_dim1 + 1], ldvl, (ftnlen)1);
#line 1288 "strevc3.f"
		    iv = 1;
#line 1289 "strevc3.f"
		} else {
#line 1290 "strevc3.f"
		    ++iv;
#line 1291 "strevc3.f"
		}
#line 1292 "strevc3.f"
	    }

/* blocked back-transform */
#line 1294 "strevc3.f"
	    ++is;
#line 1295 "strevc3.f"
	    if (ip != 0) {
#line 1295 "strevc3.f"
		++is;
#line 1295 "strevc3.f"
	    }
#line 1297 "strevc3.f"
L260:
#line 1297 "strevc3.f"
	    ;
#line 1297 "strevc3.f"
	}
#line 1298 "strevc3.f"
    }

#line 1300 "strevc3.f"
    return 0;

/*     End of STREVC3 */

} /* strevc3_ */

