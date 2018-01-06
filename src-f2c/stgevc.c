#line 1 "stgevc.f"
/* stgevc.f -- translated by f2c (version 20100827).
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

#line 1 "stgevc.f"
/* Table of constant values */

static logical c_true = TRUE_;
static integer c__2 = 2;
static doublereal c_b34 = 1.;
static integer c__1 = 1;
static doublereal c_b36 = 0.;
static logical c_false = FALSE_;

/* > \brief \b STGEVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGEVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGEVC computes some or all of the right and/or left eigenvectors of */
/* > a pair of real matrices (S,P), where S is a quasi-triangular matrix */
/* > and P is upper triangular.  Matrix pairs of this type are produced by */
/* > the generalized Schur factorization of a matrix pair (A,B): */
/* > */
/* >    A = Q*S*Z**T,  B = Q*P*Z**T */
/* > */
/* > as computed by SGGHRD + SHGEQZ. */
/* > */
/* > The right eigenvector x and the left eigenvector y of (S,P) */
/* > corresponding to an eigenvalue w are defined by: */
/* > */
/* >    S*x = w*P*x,  (y**H)*S = w*(y**H)*P, */
/* > */
/* > where y**H denotes the conjugate tranpose of y. */
/* > The eigenvalues are not input to this routine, but are computed */
/* > directly from the diagonal blocks of S and P. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of (S,P), or the products Z*X and/or Q*Y, */
/* > where Z and Q are input matrices. */
/* > If Q and Z are the orthogonal factors from the generalized Schur */
/* > factorization of a matrix pair (A,B), then Z*X and Q*Y */
/* > are the matrices of right and left eigenvectors of (A,B). */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'R': compute right eigenvectors only; */
/* >          = 'L': compute left eigenvectors only; */
/* >          = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A': compute all right and/or left eigenvectors; */
/* >          = 'B': compute all right and/or left eigenvectors, */
/* >                 backtransformed by the matrices in VR and/or VL; */
/* >          = 'S': compute selected right and/or left eigenvectors, */
/* >                 specified by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY='S', SELECT specifies the eigenvectors to be */
/* >          computed.  If w(j) is a real eigenvalue, the corresponding */
/* >          real eigenvector is computed if SELECT(j) is .TRUE.. */
/* >          If w(j) and w(j+1) are the real and imaginary parts of a */
/* >          complex eigenvalue, the corresponding complex eigenvector */
/* >          is computed if either SELECT(j) or SELECT(j+1) is .TRUE., */
/* >          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is */
/* >          set to .FALSE.. */
/* >          Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices S and P.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL array, dimension (LDS,N) */
/* >          The upper quasi-triangular matrix S from a generalized Schur */
/* >          factorization, as computed by SHGEQZ. */
/* > \endverbatim */
/* > */
/* > \param[in] LDS */
/* > \verbatim */
/* >          LDS is INTEGER */
/* >          The leading dimension of array S.  LDS >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is REAL array, dimension (LDP,N) */
/* >          The upper triangular matrix P from a generalized Schur */
/* >          factorization, as computed by SHGEQZ. */
/* >          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks */
/* >          of S must be in positive diagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDP */
/* > \verbatim */
/* >          LDP is INTEGER */
/* >          The leading dimension of array P.  LDP >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* >          VL is REAL array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of left Schur vectors returned by SHGEQZ). */
/* >          On exit, if SIDE = 'L' or 'B', VL contains: */
/* >          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P); */
/* >          if HOWMNY = 'B', the matrix Q*Y; */
/* >          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by */
/* >                      SELECT, stored consecutively in the columns of */
/* >                      VL, in the same order as their eigenvalues. */
/* > */
/* >          A complex eigenvector corresponding to a complex eigenvalue */
/* >          is stored in two consecutive columns, the first holding the */
/* >          real part, and the second the imaginary part. */
/* > */
/* >          Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of array VL.  LDVL >= 1, and if */
/* >          SIDE = 'L' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* >          VR is REAL array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Z (usually the orthogonal matrix Z */
/* >          of right Schur vectors returned by SHGEQZ). */
/* > */
/* >          On exit, if SIDE = 'R' or 'B', VR contains: */
/* >          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P); */
/* >          if HOWMNY = 'B' or 'b', the matrix Z*X; */
/* >          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P) */
/* >                      specified by SELECT, stored consecutively in the */
/* >                      columns of VR, in the same order as their */
/* >                      eigenvalues. */
/* > */
/* >          A complex eigenvector corresponding to a complex eigenvalue */
/* >          is stored in two consecutive columns, the first holding the */
/* >          real part and the second the imaginary part. */
/* > */
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
/* >          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/* >          is set to N.  Each selected real eigenvector occupies one */
/* >          column and each selected complex eigenvector occupies two */
/* >          columns. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (6*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex */
/* >                eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Allocation of workspace: */
/* >  ---------- -- --------- */
/* > */
/* >     WORK( j ) = 1-norm of j-th column of A, above the diagonal */
/* >     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal */
/* >     WORK( 2*N+1:3*N ) = real part of eigenvector */
/* >     WORK( 3*N+1:4*N ) = imaginary part of eigenvector */
/* >     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector */
/* >     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector */
/* > */
/* >  Rowwise vs. columnwise solution methods: */
/* >  ------- --  ---------- -------- ------- */
/* > */
/* >  Finding a generalized eigenvector consists basically of solving the */
/* >  singular triangular system */
/* > */
/* >   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left) */
/* > */
/* >  Consider finding the i-th right eigenvector (assume all eigenvalues */
/* >  are real). The equation to be solved is: */
/* >       n                   i */
/* >  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1 */
/* >      k=j                 k=j */
/* > */
/* >  where  C = (A - w B)  (The components v(i+1:n) are 0.) */
/* > */
/* >  The "rowwise" method is: */
/* > */
/* >  (1)  v(i) := 1 */
/* >  for j = i-1,. . .,1: */
/* >                          i */
/* >      (2) compute  s = - sum C(j,k) v(k)   and */
/* >                        k=j+1 */
/* > */
/* >      (3) v(j) := s / C(j,j) */
/* > */
/* >  Step 2 is sometimes called the "dot product" step, since it is an */
/* >  inner product between the j-th row and the portion of the eigenvector */
/* >  that has been computed so far. */
/* > */
/* >  The "columnwise" method consists basically in doing the sums */
/* >  for all the rows in parallel.  As each v(j) is computed, the */
/* >  contribution of v(j) times the j-th column of C is added to the */
/* >  partial sums.  Since FORTRAN arrays are stored columnwise, this has */
/* >  the advantage that at each step, the elements of C that are accessed */
/* >  are adjacent to one another, whereas with the rowwise method, the */
/* >  elements accessed at a step are spaced LDS (and LDP) words apart. */
/* > */
/* >  When finding left eigenvectors, the matrix in question is the */
/* >  transpose of the one in storage, so the rowwise method then */
/* >  actually accesses columns of A and B at each step, and so is the */
/* >  preferred method. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *s, integer *lds, doublereal *p, integer *ldp, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	*mm, integer *m, doublereal *work, integer *info, ftnlen side_len, 
	ftnlen howmny_len)
{
    /* System generated locals */
    integer p_dim1, p_offset, s_dim1, s_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, j, ja, jc, je, na, im, jr, jw, nw;
    static doublereal big;
    static logical lsa, lsb;
    static doublereal ulp, sum[4]	/* was [2][2] */;
    static integer ibeg, ieig, iend;
    static doublereal dmin__, temp, xmax, sump[4]	/* was [2][2] */, 
	    sums[4]	/* was [2][2] */, cim2a, cim2b, cre2a, cre2b;
    extern /* Subroutine */ int slag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal temp2, bdiag[2], acoef, scale;
    static logical ilall;
    static integer iside;
    static doublereal sbeta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical il2by2;
    static integer iinfo;
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), slaln2_(logical *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *);
    static doublereal temp2i, temp2r;
    static logical ilabad, ilbbad;
    static doublereal acoefa, bcoefa, cimaga, cimagb;
    static logical ilback;
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    static doublereal bcoefi, ascale, bscale, creala, crealb, bcoefr;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal salfar, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal xscale, bignum;
    static logical ilcomp, ilcplx;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer ihwmny;


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test the input parameters */

#line 352 "stgevc.f"
    /* Parameter adjustments */
#line 352 "stgevc.f"
    --select;
#line 352 "stgevc.f"
    s_dim1 = *lds;
#line 352 "stgevc.f"
    s_offset = 1 + s_dim1;
#line 352 "stgevc.f"
    s -= s_offset;
#line 352 "stgevc.f"
    p_dim1 = *ldp;
#line 352 "stgevc.f"
    p_offset = 1 + p_dim1;
#line 352 "stgevc.f"
    p -= p_offset;
#line 352 "stgevc.f"
    vl_dim1 = *ldvl;
#line 352 "stgevc.f"
    vl_offset = 1 + vl_dim1;
#line 352 "stgevc.f"
    vl -= vl_offset;
#line 352 "stgevc.f"
    vr_dim1 = *ldvr;
#line 352 "stgevc.f"
    vr_offset = 1 + vr_dim1;
#line 352 "stgevc.f"
    vr -= vr_offset;
#line 352 "stgevc.f"
    --work;
#line 352 "stgevc.f"

#line 352 "stgevc.f"
    /* Function Body */
#line 352 "stgevc.f"
    if (lsame_(howmny, "A", (ftnlen)1, (ftnlen)1)) {
#line 353 "stgevc.f"
	ihwmny = 1;
#line 354 "stgevc.f"
	ilall = TRUE_;
#line 355 "stgevc.f"
	ilback = FALSE_;
#line 356 "stgevc.f"
    } else if (lsame_(howmny, "S", (ftnlen)1, (ftnlen)1)) {
#line 357 "stgevc.f"
	ihwmny = 2;
#line 358 "stgevc.f"
	ilall = FALSE_;
#line 359 "stgevc.f"
	ilback = FALSE_;
#line 360 "stgevc.f"
    } else if (lsame_(howmny, "B", (ftnlen)1, (ftnlen)1)) {
#line 361 "stgevc.f"
	ihwmny = 3;
#line 362 "stgevc.f"
	ilall = TRUE_;
#line 363 "stgevc.f"
	ilback = TRUE_;
#line 364 "stgevc.f"
    } else {
#line 365 "stgevc.f"
	ihwmny = -1;
#line 366 "stgevc.f"
	ilall = TRUE_;
#line 367 "stgevc.f"
    }

#line 369 "stgevc.f"
    if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 370 "stgevc.f"
	iside = 1;
#line 371 "stgevc.f"
	compl = FALSE_;
#line 372 "stgevc.f"
	compr = TRUE_;
#line 373 "stgevc.f"
    } else if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 374 "stgevc.f"
	iside = 2;
#line 375 "stgevc.f"
	compl = TRUE_;
#line 376 "stgevc.f"
	compr = FALSE_;
#line 377 "stgevc.f"
    } else if (lsame_(side, "B", (ftnlen)1, (ftnlen)1)) {
#line 378 "stgevc.f"
	iside = 3;
#line 379 "stgevc.f"
	compl = TRUE_;
#line 380 "stgevc.f"
	compr = TRUE_;
#line 381 "stgevc.f"
    } else {
#line 382 "stgevc.f"
	iside = -1;
#line 383 "stgevc.f"
    }

#line 385 "stgevc.f"
    *info = 0;
#line 386 "stgevc.f"
    if (iside < 0) {
#line 387 "stgevc.f"
	*info = -1;
#line 388 "stgevc.f"
    } else if (ihwmny < 0) {
#line 389 "stgevc.f"
	*info = -2;
#line 390 "stgevc.f"
    } else if (*n < 0) {
#line 391 "stgevc.f"
	*info = -4;
#line 392 "stgevc.f"
    } else if (*lds < max(1,*n)) {
#line 393 "stgevc.f"
	*info = -6;
#line 394 "stgevc.f"
    } else if (*ldp < max(1,*n)) {
#line 395 "stgevc.f"
	*info = -8;
#line 396 "stgevc.f"
    }
#line 397 "stgevc.f"
    if (*info != 0) {
#line 398 "stgevc.f"
	i__1 = -(*info);
#line 398 "stgevc.f"
	xerbla_("STGEVC", &i__1, (ftnlen)6);
#line 399 "stgevc.f"
	return 0;
#line 400 "stgevc.f"
    }

/*     Count the number of eigenvectors to be computed */

#line 404 "stgevc.f"
    if (! ilall) {
#line 405 "stgevc.f"
	im = 0;
#line 406 "stgevc.f"
	ilcplx = FALSE_;
#line 407 "stgevc.f"
	i__1 = *n;
#line 407 "stgevc.f"
	for (j = 1; j <= i__1; ++j) {
#line 408 "stgevc.f"
	    if (ilcplx) {
#line 409 "stgevc.f"
		ilcplx = FALSE_;
#line 410 "stgevc.f"
		goto L10;
#line 411 "stgevc.f"
	    }
#line 412 "stgevc.f"
	    if (j < *n) {
#line 413 "stgevc.f"
		if (s[j + 1 + j * s_dim1] != 0.) {
#line 413 "stgevc.f"
		    ilcplx = TRUE_;
#line 413 "stgevc.f"
		}
#line 415 "stgevc.f"
	    }
#line 416 "stgevc.f"
	    if (ilcplx) {
#line 417 "stgevc.f"
		if (select[j] || select[j + 1]) {
#line 417 "stgevc.f"
		    im += 2;
#line 417 "stgevc.f"
		}
#line 419 "stgevc.f"
	    } else {
#line 420 "stgevc.f"
		if (select[j]) {
#line 420 "stgevc.f"
		    ++im;
#line 420 "stgevc.f"
		}
#line 422 "stgevc.f"
	    }
#line 423 "stgevc.f"
L10:
#line 423 "stgevc.f"
	    ;
#line 423 "stgevc.f"
	}
#line 424 "stgevc.f"
    } else {
#line 425 "stgevc.f"
	im = *n;
#line 426 "stgevc.f"
    }

/*     Check 2-by-2 diagonal blocks of A, B */

#line 430 "stgevc.f"
    ilabad = FALSE_;
#line 431 "stgevc.f"
    ilbbad = FALSE_;
#line 432 "stgevc.f"
    i__1 = *n - 1;
#line 432 "stgevc.f"
    for (j = 1; j <= i__1; ++j) {
#line 433 "stgevc.f"
	if (s[j + 1 + j * s_dim1] != 0.) {
#line 434 "stgevc.f"
	    if (p[j + j * p_dim1] == 0. || p[j + 1 + (j + 1) * p_dim1] == 0. 
		    || p[j + (j + 1) * p_dim1] != 0.) {
#line 434 "stgevc.f"
		ilbbad = TRUE_;
#line 434 "stgevc.f"
	    }
#line 436 "stgevc.f"
	    if (j < *n - 1) {
#line 437 "stgevc.f"
		if (s[j + 2 + (j + 1) * s_dim1] != 0.) {
#line 437 "stgevc.f"
		    ilabad = TRUE_;
#line 437 "stgevc.f"
		}
#line 439 "stgevc.f"
	    }
#line 440 "stgevc.f"
	}
#line 441 "stgevc.f"
/* L20: */
#line 441 "stgevc.f"
    }

#line 443 "stgevc.f"
    if (ilabad) {
#line 444 "stgevc.f"
	*info = -5;
#line 445 "stgevc.f"
    } else if (ilbbad) {
#line 446 "stgevc.f"
	*info = -7;
#line 447 "stgevc.f"
    } else if (compl && *ldvl < *n || *ldvl < 1) {
#line 448 "stgevc.f"
	*info = -10;
#line 449 "stgevc.f"
    } else if (compr && *ldvr < *n || *ldvr < 1) {
#line 450 "stgevc.f"
	*info = -12;
#line 451 "stgevc.f"
    } else if (*mm < im) {
#line 452 "stgevc.f"
	*info = -13;
#line 453 "stgevc.f"
    }
#line 454 "stgevc.f"
    if (*info != 0) {
#line 455 "stgevc.f"
	i__1 = -(*info);
#line 455 "stgevc.f"
	xerbla_("STGEVC", &i__1, (ftnlen)6);
#line 456 "stgevc.f"
	return 0;
#line 457 "stgevc.f"
    }

/*     Quick return if possible */

#line 461 "stgevc.f"
    *m = im;
#line 462 "stgevc.f"
    if (*n == 0) {
#line 462 "stgevc.f"
	return 0;
#line 462 "stgevc.f"
    }

/*     Machine Constants */

#line 467 "stgevc.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 468 "stgevc.f"
    big = 1. / safmin;
#line 469 "stgevc.f"
    slabad_(&safmin, &big);
#line 470 "stgevc.f"
    ulp = slamch_("Epsilon", (ftnlen)7) * slamch_("Base", (ftnlen)4);
#line 471 "stgevc.f"
    small = safmin * *n / ulp;
#line 472 "stgevc.f"
    big = 1. / small;
#line 473 "stgevc.f"
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part (i.e., excluding all elements belonging to the diagonal */
/*     blocks) of A and B to check for possible overflow in the */
/*     triangular solver. */

#line 480 "stgevc.f"
    anorm = (d__1 = s[s_dim1 + 1], abs(d__1));
#line 481 "stgevc.f"
    if (*n > 1) {
#line 481 "stgevc.f"
	anorm += (d__1 = s[s_dim1 + 2], abs(d__1));
#line 481 "stgevc.f"
    }
#line 483 "stgevc.f"
    bnorm = (d__1 = p[p_dim1 + 1], abs(d__1));
#line 484 "stgevc.f"
    work[1] = 0.;
#line 485 "stgevc.f"
    work[*n + 1] = 0.;

#line 487 "stgevc.f"
    i__1 = *n;
#line 487 "stgevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 488 "stgevc.f"
	temp = 0.;
#line 489 "stgevc.f"
	temp2 = 0.;
#line 490 "stgevc.f"
	if (s[j + (j - 1) * s_dim1] == 0.) {
#line 491 "stgevc.f"
	    iend = j - 1;
#line 492 "stgevc.f"
	} else {
#line 493 "stgevc.f"
	    iend = j - 2;
#line 494 "stgevc.f"
	}
#line 495 "stgevc.f"
	i__2 = iend;
#line 495 "stgevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 496 "stgevc.f"
	    temp += (d__1 = s[i__ + j * s_dim1], abs(d__1));
#line 497 "stgevc.f"
	    temp2 += (d__1 = p[i__ + j * p_dim1], abs(d__1));
#line 498 "stgevc.f"
/* L30: */
#line 498 "stgevc.f"
	}
#line 499 "stgevc.f"
	work[j] = temp;
#line 500 "stgevc.f"
	work[*n + j] = temp2;
/* Computing MIN */
#line 501 "stgevc.f"
	i__3 = j + 1;
#line 501 "stgevc.f"
	i__2 = min(i__3,*n);
#line 501 "stgevc.f"
	for (i__ = iend + 1; i__ <= i__2; ++i__) {
#line 502 "stgevc.f"
	    temp += (d__1 = s[i__ + j * s_dim1], abs(d__1));
#line 503 "stgevc.f"
	    temp2 += (d__1 = p[i__ + j * p_dim1], abs(d__1));
#line 504 "stgevc.f"
/* L40: */
#line 504 "stgevc.f"
	}
#line 505 "stgevc.f"
	anorm = max(anorm,temp);
#line 506 "stgevc.f"
	bnorm = max(bnorm,temp2);
#line 507 "stgevc.f"
/* L50: */
#line 507 "stgevc.f"
    }

#line 509 "stgevc.f"
    ascale = 1. / max(anorm,safmin);
#line 510 "stgevc.f"
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

#line 514 "stgevc.f"
    if (compl) {
#line 515 "stgevc.f"
	ieig = 0;

/*        Main loop over eigenvalues */

#line 519 "stgevc.f"
	ilcplx = FALSE_;
#line 520 "stgevc.f"
	i__1 = *n;
#line 520 "stgevc.f"
	for (je = 1; je <= i__1; ++je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or */
/*           (b) this would be the second of a complex pair. */
/*           Check for complex eigenvalue, so as to be sure of which */
/*           entry(-ies) of SELECT to look at. */

#line 527 "stgevc.f"
	    if (ilcplx) {
#line 528 "stgevc.f"
		ilcplx = FALSE_;
#line 529 "stgevc.f"
		goto L220;
#line 530 "stgevc.f"
	    }
#line 531 "stgevc.f"
	    nw = 1;
#line 532 "stgevc.f"
	    if (je < *n) {
#line 533 "stgevc.f"
		if (s[je + 1 + je * s_dim1] != 0.) {
#line 534 "stgevc.f"
		    ilcplx = TRUE_;
#line 535 "stgevc.f"
		    nw = 2;
#line 536 "stgevc.f"
		}
#line 537 "stgevc.f"
	    }
#line 538 "stgevc.f"
	    if (ilall) {
#line 539 "stgevc.f"
		ilcomp = TRUE_;
#line 540 "stgevc.f"
	    } else if (ilcplx) {
#line 541 "stgevc.f"
		ilcomp = select[je] || select[je + 1];
#line 542 "stgevc.f"
	    } else {
#line 543 "stgevc.f"
		ilcomp = select[je];
#line 544 "stgevc.f"
	    }
#line 545 "stgevc.f"
	    if (! ilcomp) {
#line 545 "stgevc.f"
		goto L220;
#line 545 "stgevc.f"
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or */
/*           (c) complex eigenvalue. */

#line 551 "stgevc.f"
	    if (! ilcplx) {
#line 552 "stgevc.f"
		if ((d__1 = s[je + je * s_dim1], abs(d__1)) <= safmin && (
			d__2 = p[je + je * p_dim1], abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 557 "stgevc.f"
		    ++ieig;
#line 558 "stgevc.f"
		    i__2 = *n;
#line 558 "stgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 559 "stgevc.f"
			vl[jr + ieig * vl_dim1] = 0.;
#line 560 "stgevc.f"
/* L60: */
#line 560 "stgevc.f"
		    }
#line 561 "stgevc.f"
		    vl[ieig + ieig * vl_dim1] = 1.;
#line 562 "stgevc.f"
		    goto L220;
#line 563 "stgevc.f"
		}
#line 564 "stgevc.f"
	    }

/*           Clear vector */

#line 568 "stgevc.f"
	    i__2 = nw * *n;
#line 568 "stgevc.f"
	    for (jr = 1; jr <= i__2; ++jr) {
#line 569 "stgevc.f"
		work[(*n << 1) + jr] = 0.;
#line 570 "stgevc.f"
/* L70: */
#line 570 "stgevc.f"
	    }
/*                                                 T */
/*           Compute coefficients in  ( a A - b B )  y = 0 */
/*              a  is  ACOEF */
/*              b  is  BCOEFR + i*BCOEFI */

#line 576 "stgevc.f"
	    if (! ilcplx) {

/*              Real eigenvalue */

/* Computing MAX */
#line 580 "stgevc.f"
		d__3 = (d__1 = s[je + je * s_dim1], abs(d__1)) * ascale, d__4 
			= (d__2 = p[je + je * p_dim1], abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
#line 580 "stgevc.f"
		temp = 1. / max(d__3,safmin);
#line 582 "stgevc.f"
		salfar = temp * s[je + je * s_dim1] * ascale;
#line 583 "stgevc.f"
		sbeta = temp * p[je + je * p_dim1] * bscale;
#line 584 "stgevc.f"
		acoef = sbeta * ascale;
#line 585 "stgevc.f"
		bcoefr = salfar * bscale;
#line 586 "stgevc.f"
		bcoefi = 0.;

/*              Scale to avoid underflow */

#line 590 "stgevc.f"
		scale = 1.;
#line 591 "stgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
#line 592 "stgevc.f"
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
#line 594 "stgevc.f"
		if (lsa) {
#line 594 "stgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 594 "stgevc.f"
		}
#line 596 "stgevc.f"
		if (lsb) {
/* Computing MAX */
#line 596 "stgevc.f"
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
#line 596 "stgevc.f"
		    scale = max(d__1,d__2);
#line 596 "stgevc.f"
		}
#line 599 "stgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 600 "stgevc.f"
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
#line 600 "stgevc.f"
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
#line 600 "stgevc.f"
		    scale = min(d__1,d__2);
#line 603 "stgevc.f"
		    if (lsa) {
#line 604 "stgevc.f"
			acoef = ascale * (scale * sbeta);
#line 605 "stgevc.f"
		    } else {
#line 606 "stgevc.f"
			acoef = scale * acoef;
#line 607 "stgevc.f"
		    }
#line 608 "stgevc.f"
		    if (lsb) {
#line 609 "stgevc.f"
			bcoefr = bscale * (scale * salfar);
#line 610 "stgevc.f"
		    } else {
#line 611 "stgevc.f"
			bcoefr = scale * bcoefr;
#line 612 "stgevc.f"
		    }
#line 613 "stgevc.f"
		}
#line 614 "stgevc.f"
		acoefa = abs(acoef);
#line 615 "stgevc.f"
		bcoefa = abs(bcoefr);

/*              First component is 1 */

#line 619 "stgevc.f"
		work[(*n << 1) + je] = 1.;
#line 620 "stgevc.f"
		xmax = 1.;
#line 621 "stgevc.f"
	    } else {

/*              Complex eigenvalue */

#line 625 "stgevc.f"
		d__1 = safmin * 100.;
#line 625 "stgevc.f"
		slag2_(&s[je + je * s_dim1], lds, &p[je + je * p_dim1], ldp, &
			d__1, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
#line 628 "stgevc.f"
		bcoefi = -bcoefi;
#line 629 "stgevc.f"
		if (bcoefi == 0.) {
#line 630 "stgevc.f"
		    *info = je;
#line 631 "stgevc.f"
		    return 0;
#line 632 "stgevc.f"
		}

/*              Scale to avoid over/underflow */

#line 636 "stgevc.f"
		acoefa = abs(acoef);
#line 637 "stgevc.f"
		bcoefa = abs(bcoefr) + abs(bcoefi);
#line 638 "stgevc.f"
		scale = 1.;
#line 639 "stgevc.f"
		if (acoefa * ulp < safmin && acoefa >= safmin) {
#line 639 "stgevc.f"
		    scale = safmin / ulp / acoefa;
#line 639 "stgevc.f"
		}
#line 641 "stgevc.f"
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
#line 641 "stgevc.f"
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
#line 641 "stgevc.f"
		    scale = max(d__1,d__2);
#line 641 "stgevc.f"
		}
#line 643 "stgevc.f"
		if (safmin * acoefa > ascale) {
#line 643 "stgevc.f"
		    scale = ascale / (safmin * acoefa);
#line 643 "stgevc.f"
		}
#line 645 "stgevc.f"
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
#line 645 "stgevc.f"
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
#line 645 "stgevc.f"
		    scale = min(d__1,d__2);
#line 645 "stgevc.f"
		}
#line 647 "stgevc.f"
		if (scale != 1.) {
#line 648 "stgevc.f"
		    acoef = scale * acoef;
#line 649 "stgevc.f"
		    acoefa = abs(acoef);
#line 650 "stgevc.f"
		    bcoefr = scale * bcoefr;
#line 651 "stgevc.f"
		    bcoefi = scale * bcoefi;
#line 652 "stgevc.f"
		    bcoefa = abs(bcoefr) + abs(bcoefi);
#line 653 "stgevc.f"
		}

/*              Compute first two components of eigenvector */

#line 657 "stgevc.f"
		temp = acoef * s[je + 1 + je * s_dim1];
#line 658 "stgevc.f"
		temp2r = acoef * s[je + je * s_dim1] - bcoefr * p[je + je * 
			p_dim1];
#line 659 "stgevc.f"
		temp2i = -bcoefi * p[je + je * p_dim1];
#line 660 "stgevc.f"
		if (abs(temp) > abs(temp2r) + abs(temp2i)) {
#line 661 "stgevc.f"
		    work[(*n << 1) + je] = 1.;
#line 662 "stgevc.f"
		    work[*n * 3 + je] = 0.;
#line 663 "stgevc.f"
		    work[(*n << 1) + je + 1] = -temp2r / temp;
#line 664 "stgevc.f"
		    work[*n * 3 + je + 1] = -temp2i / temp;
#line 665 "stgevc.f"
		} else {
#line 666 "stgevc.f"
		    work[(*n << 1) + je + 1] = 1.;
#line 667 "stgevc.f"
		    work[*n * 3 + je + 1] = 0.;
#line 668 "stgevc.f"
		    temp = acoef * s[je + (je + 1) * s_dim1];
#line 669 "stgevc.f"
		    work[(*n << 1) + je] = (bcoefr * p[je + 1 + (je + 1) * 
			    p_dim1] - acoef * s[je + 1 + (je + 1) * s_dim1]) /
			     temp;
#line 671 "stgevc.f"
		    work[*n * 3 + je] = bcoefi * p[je + 1 + (je + 1) * p_dim1]
			     / temp;
#line 672 "stgevc.f"
		}
/* Computing MAX */
#line 673 "stgevc.f"
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je + 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je + 1], abs(d__4));
#line 673 "stgevc.f"
		xmax = max(d__5,d__6);
#line 675 "stgevc.f"
	    }

/* Computing MAX */
#line 677 "stgevc.f"
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
#line 677 "stgevc.f"
	    dmin__ = max(d__1,safmin);

/*                                           T */
/*           Triangular solve of  (a A - b B)  y = 0 */

/*                                   T */
/*           (rowwise in  (a A - b B) , or columnwise in (a A - b B) ) */

#line 685 "stgevc.f"
	    il2by2 = FALSE_;

#line 687 "stgevc.f"
	    i__2 = *n;
#line 687 "stgevc.f"
	    for (j = je + nw; j <= i__2; ++j) {
#line 688 "stgevc.f"
		if (il2by2) {
#line 689 "stgevc.f"
		    il2by2 = FALSE_;
#line 690 "stgevc.f"
		    goto L160;
#line 691 "stgevc.f"
		}

#line 693 "stgevc.f"
		na = 1;
#line 694 "stgevc.f"
		bdiag[0] = p[j + j * p_dim1];
#line 695 "stgevc.f"
		if (j < *n) {
#line 696 "stgevc.f"
		    if (s[j + 1 + j * s_dim1] != 0.) {
#line 697 "stgevc.f"
			il2by2 = TRUE_;
#line 698 "stgevc.f"
			bdiag[1] = p[j + 1 + (j + 1) * p_dim1];
#line 699 "stgevc.f"
			na = 2;
#line 700 "stgevc.f"
		    }
#line 701 "stgevc.f"
		}

/*              Check whether scaling is necessary for dot products */

#line 705 "stgevc.f"
		xscale = 1. / max(1.,xmax);
/* Computing MAX */
#line 706 "stgevc.f"
		d__1 = work[j], d__2 = work[*n + j], d__1 = max(d__1,d__2), 
			d__2 = acoefa * work[j] + bcoefa * work[*n + j];
#line 706 "stgevc.f"
		temp = max(d__1,d__2);
#line 708 "stgevc.f"
		if (il2by2) {
/* Computing MAX */
#line 708 "stgevc.f"
		    d__1 = temp, d__2 = work[j + 1], d__1 = max(d__1,d__2), 
			    d__2 = work[*n + j + 1], d__1 = max(d__1,d__2), 
			    d__2 = acoefa * work[j + 1] + bcoefa * work[*n + 
			    j + 1];
#line 708 "stgevc.f"
		    temp = max(d__1,d__2);
#line 708 "stgevc.f"
		}
#line 711 "stgevc.f"
		if (temp > bignum * xscale) {
#line 712 "stgevc.f"
		    i__3 = nw - 1;
#line 712 "stgevc.f"
		    for (jw = 0; jw <= i__3; ++jw) {
#line 713 "stgevc.f"
			i__4 = j - 1;
#line 713 "stgevc.f"
			for (jr = je; jr <= i__4; ++jr) {
#line 714 "stgevc.f"
			    work[(jw + 2) * *n + jr] = xscale * work[(jw + 2) 
				    * *n + jr];
#line 716 "stgevc.f"
/* L80: */
#line 716 "stgevc.f"
			}
#line 717 "stgevc.f"
/* L90: */
#line 717 "stgevc.f"
		    }
#line 718 "stgevc.f"
		    xmax *= xscale;
#line 719 "stgevc.f"
		}

/*              Compute dot products */

/*                    j-1 */
/*              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k) */
/*                    k=je */

/*              To reduce the op count, this is done as */

/*              _        j-1                  _        j-1 */
/*              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) ) */
/*                       k=je                          k=je */

/*              which may cause underflow problems if A or B are close */
/*              to underflow.  (E.g., less than SMALL.) */


#line 737 "stgevc.f"
		i__3 = nw;
#line 737 "stgevc.f"
		for (jw = 1; jw <= i__3; ++jw) {
#line 738 "stgevc.f"
		    i__4 = na;
#line 738 "stgevc.f"
		    for (ja = 1; ja <= i__4; ++ja) {
#line 739 "stgevc.f"
			sums[ja + (jw << 1) - 3] = 0.;
#line 740 "stgevc.f"
			sump[ja + (jw << 1) - 3] = 0.;

#line 742 "stgevc.f"
			i__5 = j - 1;
#line 742 "stgevc.f"
			for (jr = je; jr <= i__5; ++jr) {
#line 743 "stgevc.f"
			    sums[ja + (jw << 1) - 3] += s[jr + (j + ja - 1) * 
				    s_dim1] * work[(jw + 1) * *n + jr];
#line 746 "stgevc.f"
			    sump[ja + (jw << 1) - 3] += p[jr + (j + ja - 1) * 
				    p_dim1] * work[(jw + 1) * *n + jr];
#line 749 "stgevc.f"
/* L100: */
#line 749 "stgevc.f"
			}
#line 750 "stgevc.f"
/* L110: */
#line 750 "stgevc.f"
		    }
#line 751 "stgevc.f"
/* L120: */
#line 751 "stgevc.f"
		}

#line 753 "stgevc.f"
		i__3 = na;
#line 753 "stgevc.f"
		for (ja = 1; ja <= i__3; ++ja) {
#line 754 "stgevc.f"
		    if (ilcplx) {
#line 755 "stgevc.f"
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[
				ja - 1] - bcoefi * sump[ja + 1];
#line 758 "stgevc.f"
			sum[ja + 1] = -acoef * sums[ja + 1] + bcoefr * sump[
				ja + 1] + bcoefi * sump[ja - 1];
#line 761 "stgevc.f"
		    } else {
#line 762 "stgevc.f"
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[
				ja - 1];
#line 764 "stgevc.f"
		    }
#line 765 "stgevc.f"
/* L130: */
#line 765 "stgevc.f"
		}

/*                                  T */
/*              Solve  ( a A - b B )  y = SUM(,) */
/*              with scaling and perturbation of the denominator */

#line 771 "stgevc.f"
		slaln2_(&c_true, &na, &nw, &dmin__, &acoef, &s[j + j * s_dim1]
			, lds, bdiag, &bdiag[1], sum, &c__2, &bcoefr, &bcoefi,
			 &work[(*n << 1) + j], n, &scale, &temp, &iinfo);
#line 775 "stgevc.f"
		if (scale < 1.) {
#line 776 "stgevc.f"
		    i__3 = nw - 1;
#line 776 "stgevc.f"
		    for (jw = 0; jw <= i__3; ++jw) {
#line 777 "stgevc.f"
			i__4 = j - 1;
#line 777 "stgevc.f"
			for (jr = je; jr <= i__4; ++jr) {
#line 778 "stgevc.f"
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
#line 780 "stgevc.f"
/* L140: */
#line 780 "stgevc.f"
			}
#line 781 "stgevc.f"
/* L150: */
#line 781 "stgevc.f"
		    }
#line 782 "stgevc.f"
		    xmax = scale * xmax;
#line 783 "stgevc.f"
		}
#line 784 "stgevc.f"
		xmax = max(xmax,temp);
#line 785 "stgevc.f"
L160:
#line 785 "stgevc.f"
		;
#line 785 "stgevc.f"
	    }

/*           Copy eigenvector to VL, back transforming if */
/*           HOWMNY='B'. */

#line 790 "stgevc.f"
	    ++ieig;
#line 791 "stgevc.f"
	    if (ilback) {
#line 792 "stgevc.f"
		i__2 = nw - 1;
#line 792 "stgevc.f"
		for (jw = 0; jw <= i__2; ++jw) {
#line 793 "stgevc.f"
		    i__3 = *n + 1 - je;
#line 793 "stgevc.f"
		    sgemv_("N", n, &i__3, &c_b34, &vl[je * vl_dim1 + 1], ldvl,
			     &work[(jw + 2) * *n + je], &c__1, &c_b36, &work[(
			    jw + 4) * *n + 1], &c__1, (ftnlen)1);
#line 796 "stgevc.f"
/* L170: */
#line 796 "stgevc.f"
		}
#line 797 "stgevc.f"
		slacpy_(" ", n, &nw, &work[(*n << 2) + 1], n, &vl[je * 
			vl_dim1 + 1], ldvl, (ftnlen)1);
#line 799 "stgevc.f"
		ibeg = 1;
#line 800 "stgevc.f"
	    } else {
#line 801 "stgevc.f"
		slacpy_(" ", n, &nw, &work[(*n << 1) + 1], n, &vl[ieig * 
			vl_dim1 + 1], ldvl, (ftnlen)1);
#line 803 "stgevc.f"
		ibeg = je;
#line 804 "stgevc.f"
	    }

/*           Scale eigenvector */

#line 808 "stgevc.f"
	    xmax = 0.;
#line 809 "stgevc.f"
	    if (ilcplx) {
#line 810 "stgevc.f"
		i__2 = *n;
#line 810 "stgevc.f"
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
#line 811 "stgevc.f"
		    d__3 = xmax, d__4 = (d__1 = vl[j + ieig * vl_dim1], abs(
			    d__1)) + (d__2 = vl[j + (ieig + 1) * vl_dim1], 
			    abs(d__2));
#line 811 "stgevc.f"
		    xmax = max(d__3,d__4);
#line 813 "stgevc.f"
/* L180: */
#line 813 "stgevc.f"
		}
#line 814 "stgevc.f"
	    } else {
#line 815 "stgevc.f"
		i__2 = *n;
#line 815 "stgevc.f"
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
#line 816 "stgevc.f"
		    d__2 = xmax, d__3 = (d__1 = vl[j + ieig * vl_dim1], abs(
			    d__1));
#line 816 "stgevc.f"
		    xmax = max(d__2,d__3);
#line 817 "stgevc.f"
/* L190: */
#line 817 "stgevc.f"
		}
#line 818 "stgevc.f"
	    }

#line 820 "stgevc.f"
	    if (xmax > safmin) {
#line 821 "stgevc.f"
		xscale = 1. / xmax;

#line 823 "stgevc.f"
		i__2 = nw - 1;
#line 823 "stgevc.f"
		for (jw = 0; jw <= i__2; ++jw) {
#line 824 "stgevc.f"
		    i__3 = *n;
#line 824 "stgevc.f"
		    for (jr = ibeg; jr <= i__3; ++jr) {
#line 825 "stgevc.f"
			vl[jr + (ieig + jw) * vl_dim1] = xscale * vl[jr + (
				ieig + jw) * vl_dim1];
#line 826 "stgevc.f"
/* L200: */
#line 826 "stgevc.f"
		    }
#line 827 "stgevc.f"
/* L210: */
#line 827 "stgevc.f"
		}
#line 828 "stgevc.f"
	    }
#line 829 "stgevc.f"
	    ieig = ieig + nw - 1;

#line 831 "stgevc.f"
L220:
#line 831 "stgevc.f"
	    ;
#line 831 "stgevc.f"
	}
#line 832 "stgevc.f"
    }

/*     Right eigenvectors */

#line 836 "stgevc.f"
    if (compr) {
#line 837 "stgevc.f"
	ieig = im + 1;

/*        Main loop over eigenvalues */

#line 841 "stgevc.f"
	ilcplx = FALSE_;
#line 842 "stgevc.f"
	for (je = *n; je >= 1; --je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or */
/*           (b) this would be the second of a complex pair. */
/*           Check for complex eigenvalue, so as to be sure of which */
/*           entry(-ies) of SELECT to look at -- if complex, SELECT(JE) */
/*           or SELECT(JE-1). */
/*           If this is a complex pair, the 2-by-2 diagonal block */
/*           corresponding to the eigenvalue is in rows/columns JE-1:JE */

#line 852 "stgevc.f"
	    if (ilcplx) {
#line 853 "stgevc.f"
		ilcplx = FALSE_;
#line 854 "stgevc.f"
		goto L500;
#line 855 "stgevc.f"
	    }
#line 856 "stgevc.f"
	    nw = 1;
#line 857 "stgevc.f"
	    if (je > 1) {
#line 858 "stgevc.f"
		if (s[je + (je - 1) * s_dim1] != 0.) {
#line 859 "stgevc.f"
		    ilcplx = TRUE_;
#line 860 "stgevc.f"
		    nw = 2;
#line 861 "stgevc.f"
		}
#line 862 "stgevc.f"
	    }
#line 863 "stgevc.f"
	    if (ilall) {
#line 864 "stgevc.f"
		ilcomp = TRUE_;
#line 865 "stgevc.f"
	    } else if (ilcplx) {
#line 866 "stgevc.f"
		ilcomp = select[je] || select[je - 1];
#line 867 "stgevc.f"
	    } else {
#line 868 "stgevc.f"
		ilcomp = select[je];
#line 869 "stgevc.f"
	    }
#line 870 "stgevc.f"
	    if (! ilcomp) {
#line 870 "stgevc.f"
		goto L500;
#line 870 "stgevc.f"
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or */
/*           (c) complex eigenvalue. */

#line 876 "stgevc.f"
	    if (! ilcplx) {
#line 877 "stgevc.f"
		if ((d__1 = s[je + je * s_dim1], abs(d__1)) <= safmin && (
			d__2 = p[je + je * p_dim1], abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- unit eigenvector */

#line 882 "stgevc.f"
		    --ieig;
#line 883 "stgevc.f"
		    i__1 = *n;
#line 883 "stgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 884 "stgevc.f"
			vr[jr + ieig * vr_dim1] = 0.;
#line 885 "stgevc.f"
/* L230: */
#line 885 "stgevc.f"
		    }
#line 886 "stgevc.f"
		    vr[ieig + ieig * vr_dim1] = 1.;
#line 887 "stgevc.f"
		    goto L500;
#line 888 "stgevc.f"
		}
#line 889 "stgevc.f"
	    }

/*           Clear vector */

#line 893 "stgevc.f"
	    i__1 = nw - 1;
#line 893 "stgevc.f"
	    for (jw = 0; jw <= i__1; ++jw) {
#line 894 "stgevc.f"
		i__2 = *n;
#line 894 "stgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 895 "stgevc.f"
		    work[(jw + 2) * *n + jr] = 0.;
#line 896 "stgevc.f"
/* L240: */
#line 896 "stgevc.f"
		}
#line 897 "stgevc.f"
/* L250: */
#line 897 "stgevc.f"
	    }

/*           Compute coefficients in  ( a A - b B ) x = 0 */
/*              a  is  ACOEF */
/*              b  is  BCOEFR + i*BCOEFI */

#line 903 "stgevc.f"
	    if (! ilcplx) {

/*              Real eigenvalue */

/* Computing MAX */
#line 907 "stgevc.f"
		d__3 = (d__1 = s[je + je * s_dim1], abs(d__1)) * ascale, d__4 
			= (d__2 = p[je + je * p_dim1], abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
#line 907 "stgevc.f"
		temp = 1. / max(d__3,safmin);
#line 909 "stgevc.f"
		salfar = temp * s[je + je * s_dim1] * ascale;
#line 910 "stgevc.f"
		sbeta = temp * p[je + je * p_dim1] * bscale;
#line 911 "stgevc.f"
		acoef = sbeta * ascale;
#line 912 "stgevc.f"
		bcoefr = salfar * bscale;
#line 913 "stgevc.f"
		bcoefi = 0.;

/*              Scale to avoid underflow */

#line 917 "stgevc.f"
		scale = 1.;
#line 918 "stgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
#line 919 "stgevc.f"
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
#line 921 "stgevc.f"
		if (lsa) {
#line 921 "stgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 921 "stgevc.f"
		}
#line 923 "stgevc.f"
		if (lsb) {
/* Computing MAX */
#line 923 "stgevc.f"
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
#line 923 "stgevc.f"
		    scale = max(d__1,d__2);
#line 923 "stgevc.f"
		}
#line 926 "stgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 927 "stgevc.f"
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
#line 927 "stgevc.f"
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
#line 927 "stgevc.f"
		    scale = min(d__1,d__2);
#line 930 "stgevc.f"
		    if (lsa) {
#line 931 "stgevc.f"
			acoef = ascale * (scale * sbeta);
#line 932 "stgevc.f"
		    } else {
#line 933 "stgevc.f"
			acoef = scale * acoef;
#line 934 "stgevc.f"
		    }
#line 935 "stgevc.f"
		    if (lsb) {
#line 936 "stgevc.f"
			bcoefr = bscale * (scale * salfar);
#line 937 "stgevc.f"
		    } else {
#line 938 "stgevc.f"
			bcoefr = scale * bcoefr;
#line 939 "stgevc.f"
		    }
#line 940 "stgevc.f"
		}
#line 941 "stgevc.f"
		acoefa = abs(acoef);
#line 942 "stgevc.f"
		bcoefa = abs(bcoefr);

/*              First component is 1 */

#line 946 "stgevc.f"
		work[(*n << 1) + je] = 1.;
#line 947 "stgevc.f"
		xmax = 1.;

/*              Compute contribution from column JE of A and B to sum */
/*              (See "Further Details", above.) */

#line 952 "stgevc.f"
		i__1 = je - 1;
#line 952 "stgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 953 "stgevc.f"
		    work[(*n << 1) + jr] = bcoefr * p[jr + je * p_dim1] - 
			    acoef * s[jr + je * s_dim1];
#line 955 "stgevc.f"
/* L260: */
#line 955 "stgevc.f"
		}
#line 956 "stgevc.f"
	    } else {

/*              Complex eigenvalue */

#line 960 "stgevc.f"
		d__1 = safmin * 100.;
#line 960 "stgevc.f"
		slag2_(&s[je - 1 + (je - 1) * s_dim1], lds, &p[je - 1 + (je - 
			1) * p_dim1], ldp, &d__1, &acoef, &temp, &bcoefr, &
			temp2, &bcoefi);
#line 963 "stgevc.f"
		if (bcoefi == 0.) {
#line 964 "stgevc.f"
		    *info = je - 1;
#line 965 "stgevc.f"
		    return 0;
#line 966 "stgevc.f"
		}

/*              Scale to avoid over/underflow */

#line 970 "stgevc.f"
		acoefa = abs(acoef);
#line 971 "stgevc.f"
		bcoefa = abs(bcoefr) + abs(bcoefi);
#line 972 "stgevc.f"
		scale = 1.;
#line 973 "stgevc.f"
		if (acoefa * ulp < safmin && acoefa >= safmin) {
#line 973 "stgevc.f"
		    scale = safmin / ulp / acoefa;
#line 973 "stgevc.f"
		}
#line 975 "stgevc.f"
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
#line 975 "stgevc.f"
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
#line 975 "stgevc.f"
		    scale = max(d__1,d__2);
#line 975 "stgevc.f"
		}
#line 977 "stgevc.f"
		if (safmin * acoefa > ascale) {
#line 977 "stgevc.f"
		    scale = ascale / (safmin * acoefa);
#line 977 "stgevc.f"
		}
#line 979 "stgevc.f"
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
#line 979 "stgevc.f"
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
#line 979 "stgevc.f"
		    scale = min(d__1,d__2);
#line 979 "stgevc.f"
		}
#line 981 "stgevc.f"
		if (scale != 1.) {
#line 982 "stgevc.f"
		    acoef = scale * acoef;
#line 983 "stgevc.f"
		    acoefa = abs(acoef);
#line 984 "stgevc.f"
		    bcoefr = scale * bcoefr;
#line 985 "stgevc.f"
		    bcoefi = scale * bcoefi;
#line 986 "stgevc.f"
		    bcoefa = abs(bcoefr) + abs(bcoefi);
#line 987 "stgevc.f"
		}

/*              Compute first two components of eigenvector */
/*              and contribution to sums */

#line 992 "stgevc.f"
		temp = acoef * s[je + (je - 1) * s_dim1];
#line 993 "stgevc.f"
		temp2r = acoef * s[je + je * s_dim1] - bcoefr * p[je + je * 
			p_dim1];
#line 994 "stgevc.f"
		temp2i = -bcoefi * p[je + je * p_dim1];
#line 995 "stgevc.f"
		if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
#line 996 "stgevc.f"
		    work[(*n << 1) + je] = 1.;
#line 997 "stgevc.f"
		    work[*n * 3 + je] = 0.;
#line 998 "stgevc.f"
		    work[(*n << 1) + je - 1] = -temp2r / temp;
#line 999 "stgevc.f"
		    work[*n * 3 + je - 1] = -temp2i / temp;
#line 1000 "stgevc.f"
		} else {
#line 1001 "stgevc.f"
		    work[(*n << 1) + je - 1] = 1.;
#line 1002 "stgevc.f"
		    work[*n * 3 + je - 1] = 0.;
#line 1003 "stgevc.f"
		    temp = acoef * s[je - 1 + je * s_dim1];
#line 1004 "stgevc.f"
		    work[(*n << 1) + je] = (bcoefr * p[je - 1 + (je - 1) * 
			    p_dim1] - acoef * s[je - 1 + (je - 1) * s_dim1]) /
			     temp;
#line 1006 "stgevc.f"
		    work[*n * 3 + je] = bcoefi * p[je - 1 + (je - 1) * p_dim1]
			     / temp;
#line 1007 "stgevc.f"
		}

/* Computing MAX */
#line 1009 "stgevc.f"
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je - 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je - 1], abs(d__4));
#line 1009 "stgevc.f"
		xmax = max(d__5,d__6);

/*              Compute contribution from columns JE and JE-1 */
/*              of A and B to the sums. */

#line 1015 "stgevc.f"
		creala = acoef * work[(*n << 1) + je - 1];
#line 1016 "stgevc.f"
		cimaga = acoef * work[*n * 3 + je - 1];
#line 1017 "stgevc.f"
		crealb = bcoefr * work[(*n << 1) + je - 1] - bcoefi * work[*n 
			* 3 + je - 1];
#line 1019 "stgevc.f"
		cimagb = bcoefi * work[(*n << 1) + je - 1] + bcoefr * work[*n 
			* 3 + je - 1];
#line 1021 "stgevc.f"
		cre2a = acoef * work[(*n << 1) + je];
#line 1022 "stgevc.f"
		cim2a = acoef * work[*n * 3 + je];
#line 1023 "stgevc.f"
		cre2b = bcoefr * work[(*n << 1) + je] - bcoefi * work[*n * 3 
			+ je];
#line 1024 "stgevc.f"
		cim2b = bcoefi * work[(*n << 1) + je] + bcoefr * work[*n * 3 
			+ je];
#line 1025 "stgevc.f"
		i__1 = je - 2;
#line 1025 "stgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 1026 "stgevc.f"
		    work[(*n << 1) + jr] = -creala * s[jr + (je - 1) * s_dim1]
			     + crealb * p[jr + (je - 1) * p_dim1] - cre2a * s[
			    jr + je * s_dim1] + cre2b * p[jr + je * p_dim1];
#line 1029 "stgevc.f"
		    work[*n * 3 + jr] = -cimaga * s[jr + (je - 1) * s_dim1] + 
			    cimagb * p[jr + (je - 1) * p_dim1] - cim2a * s[jr 
			    + je * s_dim1] + cim2b * p[jr + je * p_dim1];
#line 1032 "stgevc.f"
/* L270: */
#line 1032 "stgevc.f"
		}
#line 1033 "stgevc.f"
	    }

/* Computing MAX */
#line 1035 "stgevc.f"
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
#line 1035 "stgevc.f"
	    dmin__ = max(d__1,safmin);

/*           Columnwise triangular solve of  (a A - b B)  x = 0 */

#line 1039 "stgevc.f"
	    il2by2 = FALSE_;
#line 1040 "stgevc.f"
	    for (j = je - nw; j >= 1; --j) {

/*              If a 2-by-2 block, is in position j-1:j, wait until */
/*              next iteration to process it (when it will be j:j+1) */

#line 1045 "stgevc.f"
		if (! il2by2 && j > 1) {
#line 1046 "stgevc.f"
		    if (s[j + (j - 1) * s_dim1] != 0.) {
#line 1047 "stgevc.f"
			il2by2 = TRUE_;
#line 1048 "stgevc.f"
			goto L370;
#line 1049 "stgevc.f"
		    }
#line 1050 "stgevc.f"
		}
#line 1051 "stgevc.f"
		bdiag[0] = p[j + j * p_dim1];
#line 1052 "stgevc.f"
		if (il2by2) {
#line 1053 "stgevc.f"
		    na = 2;
#line 1054 "stgevc.f"
		    bdiag[1] = p[j + 1 + (j + 1) * p_dim1];
#line 1055 "stgevc.f"
		} else {
#line 1056 "stgevc.f"
		    na = 1;
#line 1057 "stgevc.f"
		}

/*              Compute x(j) (and x(j+1), if 2-by-2 block) */

#line 1061 "stgevc.f"
		slaln2_(&c_false, &na, &nw, &dmin__, &acoef, &s[j + j * 
			s_dim1], lds, bdiag, &bdiag[1], &work[(*n << 1) + j], 
			n, &bcoefr, &bcoefi, sum, &c__2, &scale, &temp, &
			iinfo);
#line 1065 "stgevc.f"
		if (scale < 1.) {

#line 1067 "stgevc.f"
		    i__1 = nw - 1;
#line 1067 "stgevc.f"
		    for (jw = 0; jw <= i__1; ++jw) {
#line 1068 "stgevc.f"
			i__2 = je;
#line 1068 "stgevc.f"
			for (jr = 1; jr <= i__2; ++jr) {
#line 1069 "stgevc.f"
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
#line 1071 "stgevc.f"
/* L280: */
#line 1071 "stgevc.f"
			}
#line 1072 "stgevc.f"
/* L290: */
#line 1072 "stgevc.f"
		    }
#line 1073 "stgevc.f"
		}
/* Computing MAX */
#line 1074 "stgevc.f"
		d__1 = scale * xmax;
#line 1074 "stgevc.f"
		xmax = max(d__1,temp);

#line 1076 "stgevc.f"
		i__1 = nw;
#line 1076 "stgevc.f"
		for (jw = 1; jw <= i__1; ++jw) {
#line 1077 "stgevc.f"
		    i__2 = na;
#line 1077 "stgevc.f"
		    for (ja = 1; ja <= i__2; ++ja) {
#line 1078 "stgevc.f"
			work[(jw + 1) * *n + j + ja - 1] = sum[ja + (jw << 1) 
				- 3];
#line 1079 "stgevc.f"
/* L300: */
#line 1079 "stgevc.f"
		    }
#line 1080 "stgevc.f"
/* L310: */
#line 1080 "stgevc.f"
		}

/*              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

#line 1084 "stgevc.f"
		if (j > 1) {

/*                 Check whether scaling is necessary for sum. */

#line 1088 "stgevc.f"
		    xscale = 1. / max(1.,xmax);
#line 1089 "stgevc.f"
		    temp = acoefa * work[j] + bcoefa * work[*n + j];
#line 1090 "stgevc.f"
		    if (il2by2) {
/* Computing MAX */
#line 1090 "stgevc.f"
			d__1 = temp, d__2 = acoefa * work[j + 1] + bcoefa * 
				work[*n + j + 1];
#line 1090 "stgevc.f"
			temp = max(d__1,d__2);
#line 1090 "stgevc.f"
		    }
/* Computing MAX */
#line 1093 "stgevc.f"
		    d__1 = max(temp,acoefa);
#line 1093 "stgevc.f"
		    temp = max(d__1,bcoefa);
#line 1094 "stgevc.f"
		    if (temp > bignum * xscale) {

#line 1096 "stgevc.f"
			i__1 = nw - 1;
#line 1096 "stgevc.f"
			for (jw = 0; jw <= i__1; ++jw) {
#line 1097 "stgevc.f"
			    i__2 = je;
#line 1097 "stgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1098 "stgevc.f"
				work[(jw + 2) * *n + jr] = xscale * work[(jw 
					+ 2) * *n + jr];
#line 1100 "stgevc.f"
/* L320: */
#line 1100 "stgevc.f"
			    }
#line 1101 "stgevc.f"
/* L330: */
#line 1101 "stgevc.f"
			}
#line 1102 "stgevc.f"
			xmax *= xscale;
#line 1103 "stgevc.f"
		    }

/*                 Compute the contributions of the off-diagonals of */
/*                 column j (and j+1, if 2-by-2 block) of A and B to the */
/*                 sums. */


#line 1110 "stgevc.f"
		    i__1 = na;
#line 1110 "stgevc.f"
		    for (ja = 1; ja <= i__1; ++ja) {
#line 1111 "stgevc.f"
			if (ilcplx) {
#line 1112 "stgevc.f"
			    creala = acoef * work[(*n << 1) + j + ja - 1];
#line 1113 "stgevc.f"
			    cimaga = acoef * work[*n * 3 + j + ja - 1];
#line 1114 "stgevc.f"
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1] - 
				    bcoefi * work[*n * 3 + j + ja - 1];
#line 1116 "stgevc.f"
			    cimagb = bcoefi * work[(*n << 1) + j + ja - 1] + 
				    bcoefr * work[*n * 3 + j + ja - 1];
#line 1118 "stgevc.f"
			    i__2 = j - 1;
#line 1118 "stgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1119 "stgevc.f"
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * s[jr + (j + ja - 1) * s_dim1]
					 + crealb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1122 "stgevc.f"
				work[*n * 3 + jr] = work[*n * 3 + jr] - 
					cimaga * s[jr + (j + ja - 1) * s_dim1]
					 + cimagb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1125 "stgevc.f"
/* L340: */
#line 1125 "stgevc.f"
			    }
#line 1126 "stgevc.f"
			} else {
#line 1127 "stgevc.f"
			    creala = acoef * work[(*n << 1) + j + ja - 1];
#line 1128 "stgevc.f"
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1];
#line 1129 "stgevc.f"
			    i__2 = j - 1;
#line 1129 "stgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1130 "stgevc.f"
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * s[jr + (j + ja - 1) * s_dim1]
					 + crealb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1133 "stgevc.f"
/* L350: */
#line 1133 "stgevc.f"
			    }
#line 1134 "stgevc.f"
			}
#line 1135 "stgevc.f"
/* L360: */
#line 1135 "stgevc.f"
		    }
#line 1136 "stgevc.f"
		}

#line 1138 "stgevc.f"
		il2by2 = FALSE_;
#line 1139 "stgevc.f"
L370:
#line 1139 "stgevc.f"
		;
#line 1139 "stgevc.f"
	    }

/*           Copy eigenvector to VR, back transforming if */
/*           HOWMNY='B'. */

#line 1144 "stgevc.f"
	    ieig -= nw;
#line 1145 "stgevc.f"
	    if (ilback) {

#line 1147 "stgevc.f"
		i__1 = nw - 1;
#line 1147 "stgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1148 "stgevc.f"
		    i__2 = *n;
#line 1148 "stgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1149 "stgevc.f"
			work[(jw + 4) * *n + jr] = work[(jw + 2) * *n + 1] * 
				vr[jr + vr_dim1];
#line 1151 "stgevc.f"
/* L380: */
#line 1151 "stgevc.f"
		    }

/*                 A series of compiler directives to defeat */
/*                 vectorization for the next loop */


#line 1157 "stgevc.f"
		    i__2 = je;
#line 1157 "stgevc.f"
		    for (jc = 2; jc <= i__2; ++jc) {
#line 1158 "stgevc.f"
			i__3 = *n;
#line 1158 "stgevc.f"
			for (jr = 1; jr <= i__3; ++jr) {
#line 1159 "stgevc.f"
			    work[(jw + 4) * *n + jr] += work[(jw + 2) * *n + 
				    jc] * vr[jr + jc * vr_dim1];
#line 1161 "stgevc.f"
/* L390: */
#line 1161 "stgevc.f"
			}
#line 1162 "stgevc.f"
/* L400: */
#line 1162 "stgevc.f"
		    }
#line 1163 "stgevc.f"
/* L410: */
#line 1163 "stgevc.f"
		}

#line 1165 "stgevc.f"
		i__1 = nw - 1;
#line 1165 "stgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1166 "stgevc.f"
		    i__2 = *n;
#line 1166 "stgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1167 "stgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = work[(jw + 4) * *n + 
				jr];
#line 1168 "stgevc.f"
/* L420: */
#line 1168 "stgevc.f"
		    }
#line 1169 "stgevc.f"
/* L430: */
#line 1169 "stgevc.f"
		}

#line 1171 "stgevc.f"
		iend = *n;
#line 1172 "stgevc.f"
	    } else {
#line 1173 "stgevc.f"
		i__1 = nw - 1;
#line 1173 "stgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1174 "stgevc.f"
		    i__2 = *n;
#line 1174 "stgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1175 "stgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = work[(jw + 2) * *n + 
				jr];
#line 1176 "stgevc.f"
/* L440: */
#line 1176 "stgevc.f"
		    }
#line 1177 "stgevc.f"
/* L450: */
#line 1177 "stgevc.f"
		}

#line 1179 "stgevc.f"
		iend = je;
#line 1180 "stgevc.f"
	    }

/*           Scale eigenvector */

#line 1184 "stgevc.f"
	    xmax = 0.;
#line 1185 "stgevc.f"
	    if (ilcplx) {
#line 1186 "stgevc.f"
		i__1 = iend;
#line 1186 "stgevc.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 1187 "stgevc.f"
		    d__3 = xmax, d__4 = (d__1 = vr[j + ieig * vr_dim1], abs(
			    d__1)) + (d__2 = vr[j + (ieig + 1) * vr_dim1], 
			    abs(d__2));
#line 1187 "stgevc.f"
		    xmax = max(d__3,d__4);
#line 1189 "stgevc.f"
/* L460: */
#line 1189 "stgevc.f"
		}
#line 1190 "stgevc.f"
	    } else {
#line 1191 "stgevc.f"
		i__1 = iend;
#line 1191 "stgevc.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 1192 "stgevc.f"
		    d__2 = xmax, d__3 = (d__1 = vr[j + ieig * vr_dim1], abs(
			    d__1));
#line 1192 "stgevc.f"
		    xmax = max(d__2,d__3);
#line 1193 "stgevc.f"
/* L470: */
#line 1193 "stgevc.f"
		}
#line 1194 "stgevc.f"
	    }

#line 1196 "stgevc.f"
	    if (xmax > safmin) {
#line 1197 "stgevc.f"
		xscale = 1. / xmax;
#line 1198 "stgevc.f"
		i__1 = nw - 1;
#line 1198 "stgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1199 "stgevc.f"
		    i__2 = iend;
#line 1199 "stgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1200 "stgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = xscale * vr[jr + (
				ieig + jw) * vr_dim1];
#line 1201 "stgevc.f"
/* L480: */
#line 1201 "stgevc.f"
		    }
#line 1202 "stgevc.f"
/* L490: */
#line 1202 "stgevc.f"
		}
#line 1203 "stgevc.f"
	    }
#line 1204 "stgevc.f"
L500:
#line 1204 "stgevc.f"
	    ;
#line 1204 "stgevc.f"
	}
#line 1205 "stgevc.f"
    }

#line 1207 "stgevc.f"
    return 0;

/*     End of STGEVC */

} /* stgevc_ */

