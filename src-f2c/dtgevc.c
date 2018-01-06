#line 1 "dtgevc.f"
/* dtgevc.f -- translated by f2c (version 20100827).
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

#line 1 "dtgevc.f"
/* Table of constant values */

static logical c_true = TRUE_;
static integer c__2 = 2;
static doublereal c_b34 = 1.;
static integer c__1 = 1;
static doublereal c_b36 = 0.;
static logical c_false = FALSE_;

/* > \brief \b DTGEVC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGEVC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgevc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgevc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgevc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, */
/*                          LDVL, VR, LDVR, MM, M, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, SIDE */
/*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       DOUBLE PRECISION   P( LDP, * ), S( LDS, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */



/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGEVC computes some or all of the right and/or left eigenvectors of */
/* > a pair of real matrices (S,P), where S is a quasi-triangular matrix */
/* > and P is upper triangular.  Matrix pairs of this type are produced by */
/* > the generalized Schur factorization of a matrix pair (A,B): */
/* > */
/* >    A = Q*S*Z**T,  B = Q*P*Z**T */
/* > */
/* > as computed by DGGHRD + DHGEQZ. */
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
/* >          S is DOUBLE PRECISION array, dimension (LDS,N) */
/* >          The upper quasi-triangular matrix S from a generalized Schur */
/* >          factorization, as computed by DHGEQZ. */
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
/* >          P is DOUBLE PRECISION array, dimension (LDP,N) */
/* >          The upper triangular matrix P from a generalized Schur */
/* >          factorization, as computed by DHGEQZ. */
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
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,MM) */
/* >          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* >          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
/* >          of left Schur vectors returned by DHGEQZ). */
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
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,MM) */
/* >          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* >          contain an N-by-N matrix Z (usually the orthogonal matrix Z */
/* >          of right Schur vectors returned by DHGEQZ). */
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
/* >          WORK is DOUBLE PRECISION array, dimension (6*N) */
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

/* > \date December 2016 */

/* > \ingroup doubleGEcomputational */

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
/* Subroutine */ int dtgevc_(char *side, char *howmny, logical *select, 
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
	    sums[4]	/* was [2][2] */;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal cim2a, cim2b, cre2a, cre2b, temp2, bdiag[2], acoef, 
	    scale;
    static logical ilall;
    static integer iside;
    static doublereal sbeta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical il2by2;
    static integer iinfo;
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    static doublereal temp2i;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    static doublereal temp2r;
    static logical ilabad, ilbbad;
    static doublereal acoefa, bcoefa, cimaga, cimagb;
    static logical ilback;
    static doublereal bcoefi, ascale, bscale, creala, crealb;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bcoefr, salfar, safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal xscale, bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ilcomp, ilcplx;
    static integer ihwmny;


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

#line 352 "dtgevc.f"
    /* Parameter adjustments */
#line 352 "dtgevc.f"
    --select;
#line 352 "dtgevc.f"
    s_dim1 = *lds;
#line 352 "dtgevc.f"
    s_offset = 1 + s_dim1;
#line 352 "dtgevc.f"
    s -= s_offset;
#line 352 "dtgevc.f"
    p_dim1 = *ldp;
#line 352 "dtgevc.f"
    p_offset = 1 + p_dim1;
#line 352 "dtgevc.f"
    p -= p_offset;
#line 352 "dtgevc.f"
    vl_dim1 = *ldvl;
#line 352 "dtgevc.f"
    vl_offset = 1 + vl_dim1;
#line 352 "dtgevc.f"
    vl -= vl_offset;
#line 352 "dtgevc.f"
    vr_dim1 = *ldvr;
#line 352 "dtgevc.f"
    vr_offset = 1 + vr_dim1;
#line 352 "dtgevc.f"
    vr -= vr_offset;
#line 352 "dtgevc.f"
    --work;
#line 352 "dtgevc.f"

#line 352 "dtgevc.f"
    /* Function Body */
#line 352 "dtgevc.f"
    if (lsame_(howmny, "A", (ftnlen)1, (ftnlen)1)) {
#line 353 "dtgevc.f"
	ihwmny = 1;
#line 354 "dtgevc.f"
	ilall = TRUE_;
#line 355 "dtgevc.f"
	ilback = FALSE_;
#line 356 "dtgevc.f"
    } else if (lsame_(howmny, "S", (ftnlen)1, (ftnlen)1)) {
#line 357 "dtgevc.f"
	ihwmny = 2;
#line 358 "dtgevc.f"
	ilall = FALSE_;
#line 359 "dtgevc.f"
	ilback = FALSE_;
#line 360 "dtgevc.f"
    } else if (lsame_(howmny, "B", (ftnlen)1, (ftnlen)1)) {
#line 361 "dtgevc.f"
	ihwmny = 3;
#line 362 "dtgevc.f"
	ilall = TRUE_;
#line 363 "dtgevc.f"
	ilback = TRUE_;
#line 364 "dtgevc.f"
    } else {
#line 365 "dtgevc.f"
	ihwmny = -1;
#line 366 "dtgevc.f"
	ilall = TRUE_;
#line 367 "dtgevc.f"
    }

#line 369 "dtgevc.f"
    if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 370 "dtgevc.f"
	iside = 1;
#line 371 "dtgevc.f"
	compl = FALSE_;
#line 372 "dtgevc.f"
	compr = TRUE_;
#line 373 "dtgevc.f"
    } else if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 374 "dtgevc.f"
	iside = 2;
#line 375 "dtgevc.f"
	compl = TRUE_;
#line 376 "dtgevc.f"
	compr = FALSE_;
#line 377 "dtgevc.f"
    } else if (lsame_(side, "B", (ftnlen)1, (ftnlen)1)) {
#line 378 "dtgevc.f"
	iside = 3;
#line 379 "dtgevc.f"
	compl = TRUE_;
#line 380 "dtgevc.f"
	compr = TRUE_;
#line 381 "dtgevc.f"
    } else {
#line 382 "dtgevc.f"
	iside = -1;
#line 383 "dtgevc.f"
    }

#line 385 "dtgevc.f"
    *info = 0;
#line 386 "dtgevc.f"
    if (iside < 0) {
#line 387 "dtgevc.f"
	*info = -1;
#line 388 "dtgevc.f"
    } else if (ihwmny < 0) {
#line 389 "dtgevc.f"
	*info = -2;
#line 390 "dtgevc.f"
    } else if (*n < 0) {
#line 391 "dtgevc.f"
	*info = -4;
#line 392 "dtgevc.f"
    } else if (*lds < max(1,*n)) {
#line 393 "dtgevc.f"
	*info = -6;
#line 394 "dtgevc.f"
    } else if (*ldp < max(1,*n)) {
#line 395 "dtgevc.f"
	*info = -8;
#line 396 "dtgevc.f"
    }
#line 397 "dtgevc.f"
    if (*info != 0) {
#line 398 "dtgevc.f"
	i__1 = -(*info);
#line 398 "dtgevc.f"
	xerbla_("DTGEVC", &i__1, (ftnlen)6);
#line 399 "dtgevc.f"
	return 0;
#line 400 "dtgevc.f"
    }

/*     Count the number of eigenvectors to be computed */

#line 404 "dtgevc.f"
    if (! ilall) {
#line 405 "dtgevc.f"
	im = 0;
#line 406 "dtgevc.f"
	ilcplx = FALSE_;
#line 407 "dtgevc.f"
	i__1 = *n;
#line 407 "dtgevc.f"
	for (j = 1; j <= i__1; ++j) {
#line 408 "dtgevc.f"
	    if (ilcplx) {
#line 409 "dtgevc.f"
		ilcplx = FALSE_;
#line 410 "dtgevc.f"
		goto L10;
#line 411 "dtgevc.f"
	    }
#line 412 "dtgevc.f"
	    if (j < *n) {
#line 413 "dtgevc.f"
		if (s[j + 1 + j * s_dim1] != 0.) {
#line 413 "dtgevc.f"
		    ilcplx = TRUE_;
#line 413 "dtgevc.f"
		}
#line 415 "dtgevc.f"
	    }
#line 416 "dtgevc.f"
	    if (ilcplx) {
#line 417 "dtgevc.f"
		if (select[j] || select[j + 1]) {
#line 417 "dtgevc.f"
		    im += 2;
#line 417 "dtgevc.f"
		}
#line 419 "dtgevc.f"
	    } else {
#line 420 "dtgevc.f"
		if (select[j]) {
#line 420 "dtgevc.f"
		    ++im;
#line 420 "dtgevc.f"
		}
#line 422 "dtgevc.f"
	    }
#line 423 "dtgevc.f"
L10:
#line 423 "dtgevc.f"
	    ;
#line 423 "dtgevc.f"
	}
#line 424 "dtgevc.f"
    } else {
#line 425 "dtgevc.f"
	im = *n;
#line 426 "dtgevc.f"
    }

/*     Check 2-by-2 diagonal blocks of A, B */

#line 430 "dtgevc.f"
    ilabad = FALSE_;
#line 431 "dtgevc.f"
    ilbbad = FALSE_;
#line 432 "dtgevc.f"
    i__1 = *n - 1;
#line 432 "dtgevc.f"
    for (j = 1; j <= i__1; ++j) {
#line 433 "dtgevc.f"
	if (s[j + 1 + j * s_dim1] != 0.) {
#line 434 "dtgevc.f"
	    if (p[j + j * p_dim1] == 0. || p[j + 1 + (j + 1) * p_dim1] == 0. 
		    || p[j + (j + 1) * p_dim1] != 0.) {
#line 434 "dtgevc.f"
		ilbbad = TRUE_;
#line 434 "dtgevc.f"
	    }
#line 436 "dtgevc.f"
	    if (j < *n - 1) {
#line 437 "dtgevc.f"
		if (s[j + 2 + (j + 1) * s_dim1] != 0.) {
#line 437 "dtgevc.f"
		    ilabad = TRUE_;
#line 437 "dtgevc.f"
		}
#line 439 "dtgevc.f"
	    }
#line 440 "dtgevc.f"
	}
#line 441 "dtgevc.f"
/* L20: */
#line 441 "dtgevc.f"
    }

#line 443 "dtgevc.f"
    if (ilabad) {
#line 444 "dtgevc.f"
	*info = -5;
#line 445 "dtgevc.f"
    } else if (ilbbad) {
#line 446 "dtgevc.f"
	*info = -7;
#line 447 "dtgevc.f"
    } else if (compl && *ldvl < *n || *ldvl < 1) {
#line 448 "dtgevc.f"
	*info = -10;
#line 449 "dtgevc.f"
    } else if (compr && *ldvr < *n || *ldvr < 1) {
#line 450 "dtgevc.f"
	*info = -12;
#line 451 "dtgevc.f"
    } else if (*mm < im) {
#line 452 "dtgevc.f"
	*info = -13;
#line 453 "dtgevc.f"
    }
#line 454 "dtgevc.f"
    if (*info != 0) {
#line 455 "dtgevc.f"
	i__1 = -(*info);
#line 455 "dtgevc.f"
	xerbla_("DTGEVC", &i__1, (ftnlen)6);
#line 456 "dtgevc.f"
	return 0;
#line 457 "dtgevc.f"
    }

/*     Quick return if possible */

#line 461 "dtgevc.f"
    *m = im;
#line 462 "dtgevc.f"
    if (*n == 0) {
#line 462 "dtgevc.f"
	return 0;
#line 462 "dtgevc.f"
    }

/*     Machine Constants */

#line 467 "dtgevc.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 468 "dtgevc.f"
    big = 1. / safmin;
#line 469 "dtgevc.f"
    dlabad_(&safmin, &big);
#line 470 "dtgevc.f"
    ulp = dlamch_("Epsilon", (ftnlen)7) * dlamch_("Base", (ftnlen)4);
#line 471 "dtgevc.f"
    small = safmin * *n / ulp;
#line 472 "dtgevc.f"
    big = 1. / small;
#line 473 "dtgevc.f"
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part (i.e., excluding all elements belonging to the diagonal */
/*     blocks) of A and B to check for possible overflow in the */
/*     triangular solver. */

#line 480 "dtgevc.f"
    anorm = (d__1 = s[s_dim1 + 1], abs(d__1));
#line 481 "dtgevc.f"
    if (*n > 1) {
#line 481 "dtgevc.f"
	anorm += (d__1 = s[s_dim1 + 2], abs(d__1));
#line 481 "dtgevc.f"
    }
#line 483 "dtgevc.f"
    bnorm = (d__1 = p[p_dim1 + 1], abs(d__1));
#line 484 "dtgevc.f"
    work[1] = 0.;
#line 485 "dtgevc.f"
    work[*n + 1] = 0.;

#line 487 "dtgevc.f"
    i__1 = *n;
#line 487 "dtgevc.f"
    for (j = 2; j <= i__1; ++j) {
#line 488 "dtgevc.f"
	temp = 0.;
#line 489 "dtgevc.f"
	temp2 = 0.;
#line 490 "dtgevc.f"
	if (s[j + (j - 1) * s_dim1] == 0.) {
#line 491 "dtgevc.f"
	    iend = j - 1;
#line 492 "dtgevc.f"
	} else {
#line 493 "dtgevc.f"
	    iend = j - 2;
#line 494 "dtgevc.f"
	}
#line 495 "dtgevc.f"
	i__2 = iend;
#line 495 "dtgevc.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 496 "dtgevc.f"
	    temp += (d__1 = s[i__ + j * s_dim1], abs(d__1));
#line 497 "dtgevc.f"
	    temp2 += (d__1 = p[i__ + j * p_dim1], abs(d__1));
#line 498 "dtgevc.f"
/* L30: */
#line 498 "dtgevc.f"
	}
#line 499 "dtgevc.f"
	work[j] = temp;
#line 500 "dtgevc.f"
	work[*n + j] = temp2;
/* Computing MIN */
#line 501 "dtgevc.f"
	i__3 = j + 1;
#line 501 "dtgevc.f"
	i__2 = min(i__3,*n);
#line 501 "dtgevc.f"
	for (i__ = iend + 1; i__ <= i__2; ++i__) {
#line 502 "dtgevc.f"
	    temp += (d__1 = s[i__ + j * s_dim1], abs(d__1));
#line 503 "dtgevc.f"
	    temp2 += (d__1 = p[i__ + j * p_dim1], abs(d__1));
#line 504 "dtgevc.f"
/* L40: */
#line 504 "dtgevc.f"
	}
#line 505 "dtgevc.f"
	anorm = max(anorm,temp);
#line 506 "dtgevc.f"
	bnorm = max(bnorm,temp2);
#line 507 "dtgevc.f"
/* L50: */
#line 507 "dtgevc.f"
    }

#line 509 "dtgevc.f"
    ascale = 1. / max(anorm,safmin);
#line 510 "dtgevc.f"
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

#line 514 "dtgevc.f"
    if (compl) {
#line 515 "dtgevc.f"
	ieig = 0;

/*        Main loop over eigenvalues */

#line 519 "dtgevc.f"
	ilcplx = FALSE_;
#line 520 "dtgevc.f"
	i__1 = *n;
#line 520 "dtgevc.f"
	for (je = 1; je <= i__1; ++je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or */
/*           (b) this would be the second of a complex pair. */
/*           Check for complex eigenvalue, so as to be sure of which */
/*           entry(-ies) of SELECT to look at. */

#line 527 "dtgevc.f"
	    if (ilcplx) {
#line 528 "dtgevc.f"
		ilcplx = FALSE_;
#line 529 "dtgevc.f"
		goto L220;
#line 530 "dtgevc.f"
	    }
#line 531 "dtgevc.f"
	    nw = 1;
#line 532 "dtgevc.f"
	    if (je < *n) {
#line 533 "dtgevc.f"
		if (s[je + 1 + je * s_dim1] != 0.) {
#line 534 "dtgevc.f"
		    ilcplx = TRUE_;
#line 535 "dtgevc.f"
		    nw = 2;
#line 536 "dtgevc.f"
		}
#line 537 "dtgevc.f"
	    }
#line 538 "dtgevc.f"
	    if (ilall) {
#line 539 "dtgevc.f"
		ilcomp = TRUE_;
#line 540 "dtgevc.f"
	    } else if (ilcplx) {
#line 541 "dtgevc.f"
		ilcomp = select[je] || select[je + 1];
#line 542 "dtgevc.f"
	    } else {
#line 543 "dtgevc.f"
		ilcomp = select[je];
#line 544 "dtgevc.f"
	    }
#line 545 "dtgevc.f"
	    if (! ilcomp) {
#line 545 "dtgevc.f"
		goto L220;
#line 545 "dtgevc.f"
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or */
/*           (c) complex eigenvalue. */

#line 551 "dtgevc.f"
	    if (! ilcplx) {
#line 552 "dtgevc.f"
		if ((d__1 = s[je + je * s_dim1], abs(d__1)) <= safmin && (
			d__2 = p[je + je * p_dim1], abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

#line 557 "dtgevc.f"
		    ++ieig;
#line 558 "dtgevc.f"
		    i__2 = *n;
#line 558 "dtgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 559 "dtgevc.f"
			vl[jr + ieig * vl_dim1] = 0.;
#line 560 "dtgevc.f"
/* L60: */
#line 560 "dtgevc.f"
		    }
#line 561 "dtgevc.f"
		    vl[ieig + ieig * vl_dim1] = 1.;
#line 562 "dtgevc.f"
		    goto L220;
#line 563 "dtgevc.f"
		}
#line 564 "dtgevc.f"
	    }

/*           Clear vector */

#line 568 "dtgevc.f"
	    i__2 = nw * *n;
#line 568 "dtgevc.f"
	    for (jr = 1; jr <= i__2; ++jr) {
#line 569 "dtgevc.f"
		work[(*n << 1) + jr] = 0.;
#line 570 "dtgevc.f"
/* L70: */
#line 570 "dtgevc.f"
	    }
/*                                                 T */
/*           Compute coefficients in  ( a A - b B )  y = 0 */
/*              a  is  ACOEF */
/*              b  is  BCOEFR + i*BCOEFI */

#line 576 "dtgevc.f"
	    if (! ilcplx) {

/*              Real eigenvalue */

/* Computing MAX */
#line 580 "dtgevc.f"
		d__3 = (d__1 = s[je + je * s_dim1], abs(d__1)) * ascale, d__4 
			= (d__2 = p[je + je * p_dim1], abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
#line 580 "dtgevc.f"
		temp = 1. / max(d__3,safmin);
#line 582 "dtgevc.f"
		salfar = temp * s[je + je * s_dim1] * ascale;
#line 583 "dtgevc.f"
		sbeta = temp * p[je + je * p_dim1] * bscale;
#line 584 "dtgevc.f"
		acoef = sbeta * ascale;
#line 585 "dtgevc.f"
		bcoefr = salfar * bscale;
#line 586 "dtgevc.f"
		bcoefi = 0.;

/*              Scale to avoid underflow */

#line 590 "dtgevc.f"
		scale = 1.;
#line 591 "dtgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
#line 592 "dtgevc.f"
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
#line 594 "dtgevc.f"
		if (lsa) {
#line 594 "dtgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 594 "dtgevc.f"
		}
#line 596 "dtgevc.f"
		if (lsb) {
/* Computing MAX */
#line 596 "dtgevc.f"
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
#line 596 "dtgevc.f"
		    scale = max(d__1,d__2);
#line 596 "dtgevc.f"
		}
#line 599 "dtgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 600 "dtgevc.f"
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
#line 600 "dtgevc.f"
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
#line 600 "dtgevc.f"
		    scale = min(d__1,d__2);
#line 603 "dtgevc.f"
		    if (lsa) {
#line 604 "dtgevc.f"
			acoef = ascale * (scale * sbeta);
#line 605 "dtgevc.f"
		    } else {
#line 606 "dtgevc.f"
			acoef = scale * acoef;
#line 607 "dtgevc.f"
		    }
#line 608 "dtgevc.f"
		    if (lsb) {
#line 609 "dtgevc.f"
			bcoefr = bscale * (scale * salfar);
#line 610 "dtgevc.f"
		    } else {
#line 611 "dtgevc.f"
			bcoefr = scale * bcoefr;
#line 612 "dtgevc.f"
		    }
#line 613 "dtgevc.f"
		}
#line 614 "dtgevc.f"
		acoefa = abs(acoef);
#line 615 "dtgevc.f"
		bcoefa = abs(bcoefr);

/*              First component is 1 */

#line 619 "dtgevc.f"
		work[(*n << 1) + je] = 1.;
#line 620 "dtgevc.f"
		xmax = 1.;
#line 621 "dtgevc.f"
	    } else {

/*              Complex eigenvalue */

#line 625 "dtgevc.f"
		d__1 = safmin * 100.;
#line 625 "dtgevc.f"
		dlag2_(&s[je + je * s_dim1], lds, &p[je + je * p_dim1], ldp, &
			d__1, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
#line 628 "dtgevc.f"
		bcoefi = -bcoefi;
#line 629 "dtgevc.f"
		if (bcoefi == 0.) {
#line 630 "dtgevc.f"
		    *info = je;
#line 631 "dtgevc.f"
		    return 0;
#line 632 "dtgevc.f"
		}

/*              Scale to avoid over/underflow */

#line 636 "dtgevc.f"
		acoefa = abs(acoef);
#line 637 "dtgevc.f"
		bcoefa = abs(bcoefr) + abs(bcoefi);
#line 638 "dtgevc.f"
		scale = 1.;
#line 639 "dtgevc.f"
		if (acoefa * ulp < safmin && acoefa >= safmin) {
#line 639 "dtgevc.f"
		    scale = safmin / ulp / acoefa;
#line 639 "dtgevc.f"
		}
#line 641 "dtgevc.f"
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
#line 641 "dtgevc.f"
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
#line 641 "dtgevc.f"
		    scale = max(d__1,d__2);
#line 641 "dtgevc.f"
		}
#line 643 "dtgevc.f"
		if (safmin * acoefa > ascale) {
#line 643 "dtgevc.f"
		    scale = ascale / (safmin * acoefa);
#line 643 "dtgevc.f"
		}
#line 645 "dtgevc.f"
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
#line 645 "dtgevc.f"
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
#line 645 "dtgevc.f"
		    scale = min(d__1,d__2);
#line 645 "dtgevc.f"
		}
#line 647 "dtgevc.f"
		if (scale != 1.) {
#line 648 "dtgevc.f"
		    acoef = scale * acoef;
#line 649 "dtgevc.f"
		    acoefa = abs(acoef);
#line 650 "dtgevc.f"
		    bcoefr = scale * bcoefr;
#line 651 "dtgevc.f"
		    bcoefi = scale * bcoefi;
#line 652 "dtgevc.f"
		    bcoefa = abs(bcoefr) + abs(bcoefi);
#line 653 "dtgevc.f"
		}

/*              Compute first two components of eigenvector */

#line 657 "dtgevc.f"
		temp = acoef * s[je + 1 + je * s_dim1];
#line 658 "dtgevc.f"
		temp2r = acoef * s[je + je * s_dim1] - bcoefr * p[je + je * 
			p_dim1];
#line 659 "dtgevc.f"
		temp2i = -bcoefi * p[je + je * p_dim1];
#line 660 "dtgevc.f"
		if (abs(temp) > abs(temp2r) + abs(temp2i)) {
#line 661 "dtgevc.f"
		    work[(*n << 1) + je] = 1.;
#line 662 "dtgevc.f"
		    work[*n * 3 + je] = 0.;
#line 663 "dtgevc.f"
		    work[(*n << 1) + je + 1] = -temp2r / temp;
#line 664 "dtgevc.f"
		    work[*n * 3 + je + 1] = -temp2i / temp;
#line 665 "dtgevc.f"
		} else {
#line 666 "dtgevc.f"
		    work[(*n << 1) + je + 1] = 1.;
#line 667 "dtgevc.f"
		    work[*n * 3 + je + 1] = 0.;
#line 668 "dtgevc.f"
		    temp = acoef * s[je + (je + 1) * s_dim1];
#line 669 "dtgevc.f"
		    work[(*n << 1) + je] = (bcoefr * p[je + 1 + (je + 1) * 
			    p_dim1] - acoef * s[je + 1 + (je + 1) * s_dim1]) /
			     temp;
#line 671 "dtgevc.f"
		    work[*n * 3 + je] = bcoefi * p[je + 1 + (je + 1) * p_dim1]
			     / temp;
#line 672 "dtgevc.f"
		}
/* Computing MAX */
#line 673 "dtgevc.f"
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je + 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je + 1], abs(d__4));
#line 673 "dtgevc.f"
		xmax = max(d__5,d__6);
#line 675 "dtgevc.f"
	    }

/* Computing MAX */
#line 677 "dtgevc.f"
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
#line 677 "dtgevc.f"
	    dmin__ = max(d__1,safmin);

/*                                           T */
/*           Triangular solve of  (a A - b B)  y = 0 */

/*                                   T */
/*           (rowwise in  (a A - b B) , or columnwise in (a A - b B) ) */

#line 685 "dtgevc.f"
	    il2by2 = FALSE_;

#line 687 "dtgevc.f"
	    i__2 = *n;
#line 687 "dtgevc.f"
	    for (j = je + nw; j <= i__2; ++j) {
#line 688 "dtgevc.f"
		if (il2by2) {
#line 689 "dtgevc.f"
		    il2by2 = FALSE_;
#line 690 "dtgevc.f"
		    goto L160;
#line 691 "dtgevc.f"
		}

#line 693 "dtgevc.f"
		na = 1;
#line 694 "dtgevc.f"
		bdiag[0] = p[j + j * p_dim1];
#line 695 "dtgevc.f"
		if (j < *n) {
#line 696 "dtgevc.f"
		    if (s[j + 1 + j * s_dim1] != 0.) {
#line 697 "dtgevc.f"
			il2by2 = TRUE_;
#line 698 "dtgevc.f"
			bdiag[1] = p[j + 1 + (j + 1) * p_dim1];
#line 699 "dtgevc.f"
			na = 2;
#line 700 "dtgevc.f"
		    }
#line 701 "dtgevc.f"
		}

/*              Check whether scaling is necessary for dot products */

#line 705 "dtgevc.f"
		xscale = 1. / max(1.,xmax);
/* Computing MAX */
#line 706 "dtgevc.f"
		d__1 = work[j], d__2 = work[*n + j], d__1 = max(d__1,d__2), 
			d__2 = acoefa * work[j] + bcoefa * work[*n + j];
#line 706 "dtgevc.f"
		temp = max(d__1,d__2);
#line 708 "dtgevc.f"
		if (il2by2) {
/* Computing MAX */
#line 708 "dtgevc.f"
		    d__1 = temp, d__2 = work[j + 1], d__1 = max(d__1,d__2), 
			    d__2 = work[*n + j + 1], d__1 = max(d__1,d__2), 
			    d__2 = acoefa * work[j + 1] + bcoefa * work[*n + 
			    j + 1];
#line 708 "dtgevc.f"
		    temp = max(d__1,d__2);
#line 708 "dtgevc.f"
		}
#line 711 "dtgevc.f"
		if (temp > bignum * xscale) {
#line 712 "dtgevc.f"
		    i__3 = nw - 1;
#line 712 "dtgevc.f"
		    for (jw = 0; jw <= i__3; ++jw) {
#line 713 "dtgevc.f"
			i__4 = j - 1;
#line 713 "dtgevc.f"
			for (jr = je; jr <= i__4; ++jr) {
#line 714 "dtgevc.f"
			    work[(jw + 2) * *n + jr] = xscale * work[(jw + 2) 
				    * *n + jr];
#line 716 "dtgevc.f"
/* L80: */
#line 716 "dtgevc.f"
			}
#line 717 "dtgevc.f"
/* L90: */
#line 717 "dtgevc.f"
		    }
#line 718 "dtgevc.f"
		    xmax *= xscale;
#line 719 "dtgevc.f"
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


#line 737 "dtgevc.f"
		i__3 = nw;
#line 737 "dtgevc.f"
		for (jw = 1; jw <= i__3; ++jw) {
#line 738 "dtgevc.f"
		    i__4 = na;
#line 738 "dtgevc.f"
		    for (ja = 1; ja <= i__4; ++ja) {
#line 739 "dtgevc.f"
			sums[ja + (jw << 1) - 3] = 0.;
#line 740 "dtgevc.f"
			sump[ja + (jw << 1) - 3] = 0.;

#line 742 "dtgevc.f"
			i__5 = j - 1;
#line 742 "dtgevc.f"
			for (jr = je; jr <= i__5; ++jr) {
#line 743 "dtgevc.f"
			    sums[ja + (jw << 1) - 3] += s[jr + (j + ja - 1) * 
				    s_dim1] * work[(jw + 1) * *n + jr];
#line 746 "dtgevc.f"
			    sump[ja + (jw << 1) - 3] += p[jr + (j + ja - 1) * 
				    p_dim1] * work[(jw + 1) * *n + jr];
#line 749 "dtgevc.f"
/* L100: */
#line 749 "dtgevc.f"
			}
#line 750 "dtgevc.f"
/* L110: */
#line 750 "dtgevc.f"
		    }
#line 751 "dtgevc.f"
/* L120: */
#line 751 "dtgevc.f"
		}

#line 753 "dtgevc.f"
		i__3 = na;
#line 753 "dtgevc.f"
		for (ja = 1; ja <= i__3; ++ja) {
#line 754 "dtgevc.f"
		    if (ilcplx) {
#line 755 "dtgevc.f"
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[
				ja - 1] - bcoefi * sump[ja + 1];
#line 758 "dtgevc.f"
			sum[ja + 1] = -acoef * sums[ja + 1] + bcoefr * sump[
				ja + 1] + bcoefi * sump[ja - 1];
#line 761 "dtgevc.f"
		    } else {
#line 762 "dtgevc.f"
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[
				ja - 1];
#line 764 "dtgevc.f"
		    }
#line 765 "dtgevc.f"
/* L130: */
#line 765 "dtgevc.f"
		}

/*                                  T */
/*              Solve  ( a A - b B )  y = SUM(,) */
/*              with scaling and perturbation of the denominator */

#line 771 "dtgevc.f"
		dlaln2_(&c_true, &na, &nw, &dmin__, &acoef, &s[j + j * s_dim1]
			, lds, bdiag, &bdiag[1], sum, &c__2, &bcoefr, &bcoefi,
			 &work[(*n << 1) + j], n, &scale, &temp, &iinfo);
#line 775 "dtgevc.f"
		if (scale < 1.) {
#line 776 "dtgevc.f"
		    i__3 = nw - 1;
#line 776 "dtgevc.f"
		    for (jw = 0; jw <= i__3; ++jw) {
#line 777 "dtgevc.f"
			i__4 = j - 1;
#line 777 "dtgevc.f"
			for (jr = je; jr <= i__4; ++jr) {
#line 778 "dtgevc.f"
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
#line 780 "dtgevc.f"
/* L140: */
#line 780 "dtgevc.f"
			}
#line 781 "dtgevc.f"
/* L150: */
#line 781 "dtgevc.f"
		    }
#line 782 "dtgevc.f"
		    xmax = scale * xmax;
#line 783 "dtgevc.f"
		}
#line 784 "dtgevc.f"
		xmax = max(xmax,temp);
#line 785 "dtgevc.f"
L160:
#line 785 "dtgevc.f"
		;
#line 785 "dtgevc.f"
	    }

/*           Copy eigenvector to VL, back transforming if */
/*           HOWMNY='B'. */

#line 790 "dtgevc.f"
	    ++ieig;
#line 791 "dtgevc.f"
	    if (ilback) {
#line 792 "dtgevc.f"
		i__2 = nw - 1;
#line 792 "dtgevc.f"
		for (jw = 0; jw <= i__2; ++jw) {
#line 793 "dtgevc.f"
		    i__3 = *n + 1 - je;
#line 793 "dtgevc.f"
		    dgemv_("N", n, &i__3, &c_b34, &vl[je * vl_dim1 + 1], ldvl,
			     &work[(jw + 2) * *n + je], &c__1, &c_b36, &work[(
			    jw + 4) * *n + 1], &c__1, (ftnlen)1);
#line 796 "dtgevc.f"
/* L170: */
#line 796 "dtgevc.f"
		}
#line 797 "dtgevc.f"
		dlacpy_(" ", n, &nw, &work[(*n << 2) + 1], n, &vl[je * 
			vl_dim1 + 1], ldvl, (ftnlen)1);
#line 799 "dtgevc.f"
		ibeg = 1;
#line 800 "dtgevc.f"
	    } else {
#line 801 "dtgevc.f"
		dlacpy_(" ", n, &nw, &work[(*n << 1) + 1], n, &vl[ieig * 
			vl_dim1 + 1], ldvl, (ftnlen)1);
#line 803 "dtgevc.f"
		ibeg = je;
#line 804 "dtgevc.f"
	    }

/*           Scale eigenvector */

#line 808 "dtgevc.f"
	    xmax = 0.;
#line 809 "dtgevc.f"
	    if (ilcplx) {
#line 810 "dtgevc.f"
		i__2 = *n;
#line 810 "dtgevc.f"
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
#line 811 "dtgevc.f"
		    d__3 = xmax, d__4 = (d__1 = vl[j + ieig * vl_dim1], abs(
			    d__1)) + (d__2 = vl[j + (ieig + 1) * vl_dim1], 
			    abs(d__2));
#line 811 "dtgevc.f"
		    xmax = max(d__3,d__4);
#line 813 "dtgevc.f"
/* L180: */
#line 813 "dtgevc.f"
		}
#line 814 "dtgevc.f"
	    } else {
#line 815 "dtgevc.f"
		i__2 = *n;
#line 815 "dtgevc.f"
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
#line 816 "dtgevc.f"
		    d__2 = xmax, d__3 = (d__1 = vl[j + ieig * vl_dim1], abs(
			    d__1));
#line 816 "dtgevc.f"
		    xmax = max(d__2,d__3);
#line 817 "dtgevc.f"
/* L190: */
#line 817 "dtgevc.f"
		}
#line 818 "dtgevc.f"
	    }

#line 820 "dtgevc.f"
	    if (xmax > safmin) {
#line 821 "dtgevc.f"
		xscale = 1. / xmax;

#line 823 "dtgevc.f"
		i__2 = nw - 1;
#line 823 "dtgevc.f"
		for (jw = 0; jw <= i__2; ++jw) {
#line 824 "dtgevc.f"
		    i__3 = *n;
#line 824 "dtgevc.f"
		    for (jr = ibeg; jr <= i__3; ++jr) {
#line 825 "dtgevc.f"
			vl[jr + (ieig + jw) * vl_dim1] = xscale * vl[jr + (
				ieig + jw) * vl_dim1];
#line 826 "dtgevc.f"
/* L200: */
#line 826 "dtgevc.f"
		    }
#line 827 "dtgevc.f"
/* L210: */
#line 827 "dtgevc.f"
		}
#line 828 "dtgevc.f"
	    }
#line 829 "dtgevc.f"
	    ieig = ieig + nw - 1;

#line 831 "dtgevc.f"
L220:
#line 831 "dtgevc.f"
	    ;
#line 831 "dtgevc.f"
	}
#line 832 "dtgevc.f"
    }

/*     Right eigenvectors */

#line 836 "dtgevc.f"
    if (compr) {
#line 837 "dtgevc.f"
	ieig = im + 1;

/*        Main loop over eigenvalues */

#line 841 "dtgevc.f"
	ilcplx = FALSE_;
#line 842 "dtgevc.f"
	for (je = *n; je >= 1; --je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or */
/*           (b) this would be the second of a complex pair. */
/*           Check for complex eigenvalue, so as to be sure of which */
/*           entry(-ies) of SELECT to look at -- if complex, SELECT(JE) */
/*           or SELECT(JE-1). */
/*           If this is a complex pair, the 2-by-2 diagonal block */
/*           corresponding to the eigenvalue is in rows/columns JE-1:JE */

#line 852 "dtgevc.f"
	    if (ilcplx) {
#line 853 "dtgevc.f"
		ilcplx = FALSE_;
#line 854 "dtgevc.f"
		goto L500;
#line 855 "dtgevc.f"
	    }
#line 856 "dtgevc.f"
	    nw = 1;
#line 857 "dtgevc.f"
	    if (je > 1) {
#line 858 "dtgevc.f"
		if (s[je + (je - 1) * s_dim1] != 0.) {
#line 859 "dtgevc.f"
		    ilcplx = TRUE_;
#line 860 "dtgevc.f"
		    nw = 2;
#line 861 "dtgevc.f"
		}
#line 862 "dtgevc.f"
	    }
#line 863 "dtgevc.f"
	    if (ilall) {
#line 864 "dtgevc.f"
		ilcomp = TRUE_;
#line 865 "dtgevc.f"
	    } else if (ilcplx) {
#line 866 "dtgevc.f"
		ilcomp = select[je] || select[je - 1];
#line 867 "dtgevc.f"
	    } else {
#line 868 "dtgevc.f"
		ilcomp = select[je];
#line 869 "dtgevc.f"
	    }
#line 870 "dtgevc.f"
	    if (! ilcomp) {
#line 870 "dtgevc.f"
		goto L500;
#line 870 "dtgevc.f"
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or */
/*           (c) complex eigenvalue. */

#line 876 "dtgevc.f"
	    if (! ilcplx) {
#line 877 "dtgevc.f"
		if ((d__1 = s[je + je * s_dim1], abs(d__1)) <= safmin && (
			d__2 = p[je + je * p_dim1], abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- unit eigenvector */

#line 882 "dtgevc.f"
		    --ieig;
#line 883 "dtgevc.f"
		    i__1 = *n;
#line 883 "dtgevc.f"
		    for (jr = 1; jr <= i__1; ++jr) {
#line 884 "dtgevc.f"
			vr[jr + ieig * vr_dim1] = 0.;
#line 885 "dtgevc.f"
/* L230: */
#line 885 "dtgevc.f"
		    }
#line 886 "dtgevc.f"
		    vr[ieig + ieig * vr_dim1] = 1.;
#line 887 "dtgevc.f"
		    goto L500;
#line 888 "dtgevc.f"
		}
#line 889 "dtgevc.f"
	    }

/*           Clear vector */

#line 893 "dtgevc.f"
	    i__1 = nw - 1;
#line 893 "dtgevc.f"
	    for (jw = 0; jw <= i__1; ++jw) {
#line 894 "dtgevc.f"
		i__2 = *n;
#line 894 "dtgevc.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 895 "dtgevc.f"
		    work[(jw + 2) * *n + jr] = 0.;
#line 896 "dtgevc.f"
/* L240: */
#line 896 "dtgevc.f"
		}
#line 897 "dtgevc.f"
/* L250: */
#line 897 "dtgevc.f"
	    }

/*           Compute coefficients in  ( a A - b B ) x = 0 */
/*              a  is  ACOEF */
/*              b  is  BCOEFR + i*BCOEFI */

#line 903 "dtgevc.f"
	    if (! ilcplx) {

/*              Real eigenvalue */

/* Computing MAX */
#line 907 "dtgevc.f"
		d__3 = (d__1 = s[je + je * s_dim1], abs(d__1)) * ascale, d__4 
			= (d__2 = p[je + je * p_dim1], abs(d__2)) * bscale, 
			d__3 = max(d__3,d__4);
#line 907 "dtgevc.f"
		temp = 1. / max(d__3,safmin);
#line 909 "dtgevc.f"
		salfar = temp * s[je + je * s_dim1] * ascale;
#line 910 "dtgevc.f"
		sbeta = temp * p[je + je * p_dim1] * bscale;
#line 911 "dtgevc.f"
		acoef = sbeta * ascale;
#line 912 "dtgevc.f"
		bcoefr = salfar * bscale;
#line 913 "dtgevc.f"
		bcoefi = 0.;

/*              Scale to avoid underflow */

#line 917 "dtgevc.f"
		scale = 1.;
#line 918 "dtgevc.f"
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
#line 919 "dtgevc.f"
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
#line 921 "dtgevc.f"
		if (lsa) {
#line 921 "dtgevc.f"
		    scale = small / abs(sbeta) * min(anorm,big);
#line 921 "dtgevc.f"
		}
#line 923 "dtgevc.f"
		if (lsb) {
/* Computing MAX */
#line 923 "dtgevc.f"
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
#line 923 "dtgevc.f"
		    scale = max(d__1,d__2);
#line 923 "dtgevc.f"
		}
#line 926 "dtgevc.f"
		if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
#line 927 "dtgevc.f"
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
#line 927 "dtgevc.f"
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
#line 927 "dtgevc.f"
		    scale = min(d__1,d__2);
#line 930 "dtgevc.f"
		    if (lsa) {
#line 931 "dtgevc.f"
			acoef = ascale * (scale * sbeta);
#line 932 "dtgevc.f"
		    } else {
#line 933 "dtgevc.f"
			acoef = scale * acoef;
#line 934 "dtgevc.f"
		    }
#line 935 "dtgevc.f"
		    if (lsb) {
#line 936 "dtgevc.f"
			bcoefr = bscale * (scale * salfar);
#line 937 "dtgevc.f"
		    } else {
#line 938 "dtgevc.f"
			bcoefr = scale * bcoefr;
#line 939 "dtgevc.f"
		    }
#line 940 "dtgevc.f"
		}
#line 941 "dtgevc.f"
		acoefa = abs(acoef);
#line 942 "dtgevc.f"
		bcoefa = abs(bcoefr);

/*              First component is 1 */

#line 946 "dtgevc.f"
		work[(*n << 1) + je] = 1.;
#line 947 "dtgevc.f"
		xmax = 1.;

/*              Compute contribution from column JE of A and B to sum */
/*              (See "Further Details", above.) */

#line 952 "dtgevc.f"
		i__1 = je - 1;
#line 952 "dtgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 953 "dtgevc.f"
		    work[(*n << 1) + jr] = bcoefr * p[jr + je * p_dim1] - 
			    acoef * s[jr + je * s_dim1];
#line 955 "dtgevc.f"
/* L260: */
#line 955 "dtgevc.f"
		}
#line 956 "dtgevc.f"
	    } else {

/*              Complex eigenvalue */

#line 960 "dtgevc.f"
		d__1 = safmin * 100.;
#line 960 "dtgevc.f"
		dlag2_(&s[je - 1 + (je - 1) * s_dim1], lds, &p[je - 1 + (je - 
			1) * p_dim1], ldp, &d__1, &acoef, &temp, &bcoefr, &
			temp2, &bcoefi);
#line 963 "dtgevc.f"
		if (bcoefi == 0.) {
#line 964 "dtgevc.f"
		    *info = je - 1;
#line 965 "dtgevc.f"
		    return 0;
#line 966 "dtgevc.f"
		}

/*              Scale to avoid over/underflow */

#line 970 "dtgevc.f"
		acoefa = abs(acoef);
#line 971 "dtgevc.f"
		bcoefa = abs(bcoefr) + abs(bcoefi);
#line 972 "dtgevc.f"
		scale = 1.;
#line 973 "dtgevc.f"
		if (acoefa * ulp < safmin && acoefa >= safmin) {
#line 973 "dtgevc.f"
		    scale = safmin / ulp / acoefa;
#line 973 "dtgevc.f"
		}
#line 975 "dtgevc.f"
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
#line 975 "dtgevc.f"
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
#line 975 "dtgevc.f"
		    scale = max(d__1,d__2);
#line 975 "dtgevc.f"
		}
#line 977 "dtgevc.f"
		if (safmin * acoefa > ascale) {
#line 977 "dtgevc.f"
		    scale = ascale / (safmin * acoefa);
#line 977 "dtgevc.f"
		}
#line 979 "dtgevc.f"
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
#line 979 "dtgevc.f"
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
#line 979 "dtgevc.f"
		    scale = min(d__1,d__2);
#line 979 "dtgevc.f"
		}
#line 981 "dtgevc.f"
		if (scale != 1.) {
#line 982 "dtgevc.f"
		    acoef = scale * acoef;
#line 983 "dtgevc.f"
		    acoefa = abs(acoef);
#line 984 "dtgevc.f"
		    bcoefr = scale * bcoefr;
#line 985 "dtgevc.f"
		    bcoefi = scale * bcoefi;
#line 986 "dtgevc.f"
		    bcoefa = abs(bcoefr) + abs(bcoefi);
#line 987 "dtgevc.f"
		}

/*              Compute first two components of eigenvector */
/*              and contribution to sums */

#line 992 "dtgevc.f"
		temp = acoef * s[je + (je - 1) * s_dim1];
#line 993 "dtgevc.f"
		temp2r = acoef * s[je + je * s_dim1] - bcoefr * p[je + je * 
			p_dim1];
#line 994 "dtgevc.f"
		temp2i = -bcoefi * p[je + je * p_dim1];
#line 995 "dtgevc.f"
		if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
#line 996 "dtgevc.f"
		    work[(*n << 1) + je] = 1.;
#line 997 "dtgevc.f"
		    work[*n * 3 + je] = 0.;
#line 998 "dtgevc.f"
		    work[(*n << 1) + je - 1] = -temp2r / temp;
#line 999 "dtgevc.f"
		    work[*n * 3 + je - 1] = -temp2i / temp;
#line 1000 "dtgevc.f"
		} else {
#line 1001 "dtgevc.f"
		    work[(*n << 1) + je - 1] = 1.;
#line 1002 "dtgevc.f"
		    work[*n * 3 + je - 1] = 0.;
#line 1003 "dtgevc.f"
		    temp = acoef * s[je - 1 + je * s_dim1];
#line 1004 "dtgevc.f"
		    work[(*n << 1) + je] = (bcoefr * p[je - 1 + (je - 1) * 
			    p_dim1] - acoef * s[je - 1 + (je - 1) * s_dim1]) /
			     temp;
#line 1006 "dtgevc.f"
		    work[*n * 3 + je] = bcoefi * p[je - 1 + (je - 1) * p_dim1]
			     / temp;
#line 1007 "dtgevc.f"
		}

/* Computing MAX */
#line 1009 "dtgevc.f"
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je - 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je - 1], abs(d__4));
#line 1009 "dtgevc.f"
		xmax = max(d__5,d__6);

/*              Compute contribution from columns JE and JE-1 */
/*              of A and B to the sums. */

#line 1015 "dtgevc.f"
		creala = acoef * work[(*n << 1) + je - 1];
#line 1016 "dtgevc.f"
		cimaga = acoef * work[*n * 3 + je - 1];
#line 1017 "dtgevc.f"
		crealb = bcoefr * work[(*n << 1) + je - 1] - bcoefi * work[*n 
			* 3 + je - 1];
#line 1019 "dtgevc.f"
		cimagb = bcoefi * work[(*n << 1) + je - 1] + bcoefr * work[*n 
			* 3 + je - 1];
#line 1021 "dtgevc.f"
		cre2a = acoef * work[(*n << 1) + je];
#line 1022 "dtgevc.f"
		cim2a = acoef * work[*n * 3 + je];
#line 1023 "dtgevc.f"
		cre2b = bcoefr * work[(*n << 1) + je] - bcoefi * work[*n * 3 
			+ je];
#line 1024 "dtgevc.f"
		cim2b = bcoefi * work[(*n << 1) + je] + bcoefr * work[*n * 3 
			+ je];
#line 1025 "dtgevc.f"
		i__1 = je - 2;
#line 1025 "dtgevc.f"
		for (jr = 1; jr <= i__1; ++jr) {
#line 1026 "dtgevc.f"
		    work[(*n << 1) + jr] = -creala * s[jr + (je - 1) * s_dim1]
			     + crealb * p[jr + (je - 1) * p_dim1] - cre2a * s[
			    jr + je * s_dim1] + cre2b * p[jr + je * p_dim1];
#line 1029 "dtgevc.f"
		    work[*n * 3 + jr] = -cimaga * s[jr + (je - 1) * s_dim1] + 
			    cimagb * p[jr + (je - 1) * p_dim1] - cim2a * s[jr 
			    + je * s_dim1] + cim2b * p[jr + je * p_dim1];
#line 1032 "dtgevc.f"
/* L270: */
#line 1032 "dtgevc.f"
		}
#line 1033 "dtgevc.f"
	    }

/* Computing MAX */
#line 1035 "dtgevc.f"
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
#line 1035 "dtgevc.f"
	    dmin__ = max(d__1,safmin);

/*           Columnwise triangular solve of  (a A - b B)  x = 0 */

#line 1039 "dtgevc.f"
	    il2by2 = FALSE_;
#line 1040 "dtgevc.f"
	    for (j = je - nw; j >= 1; --j) {

/*              If a 2-by-2 block, is in position j-1:j, wait until */
/*              next iteration to process it (when it will be j:j+1) */

#line 1045 "dtgevc.f"
		if (! il2by2 && j > 1) {
#line 1046 "dtgevc.f"
		    if (s[j + (j - 1) * s_dim1] != 0.) {
#line 1047 "dtgevc.f"
			il2by2 = TRUE_;
#line 1048 "dtgevc.f"
			goto L370;
#line 1049 "dtgevc.f"
		    }
#line 1050 "dtgevc.f"
		}
#line 1051 "dtgevc.f"
		bdiag[0] = p[j + j * p_dim1];
#line 1052 "dtgevc.f"
		if (il2by2) {
#line 1053 "dtgevc.f"
		    na = 2;
#line 1054 "dtgevc.f"
		    bdiag[1] = p[j + 1 + (j + 1) * p_dim1];
#line 1055 "dtgevc.f"
		} else {
#line 1056 "dtgevc.f"
		    na = 1;
#line 1057 "dtgevc.f"
		}

/*              Compute x(j) (and x(j+1), if 2-by-2 block) */

#line 1061 "dtgevc.f"
		dlaln2_(&c_false, &na, &nw, &dmin__, &acoef, &s[j + j * 
			s_dim1], lds, bdiag, &bdiag[1], &work[(*n << 1) + j], 
			n, &bcoefr, &bcoefi, sum, &c__2, &scale, &temp, &
			iinfo);
#line 1065 "dtgevc.f"
		if (scale < 1.) {

#line 1067 "dtgevc.f"
		    i__1 = nw - 1;
#line 1067 "dtgevc.f"
		    for (jw = 0; jw <= i__1; ++jw) {
#line 1068 "dtgevc.f"
			i__2 = je;
#line 1068 "dtgevc.f"
			for (jr = 1; jr <= i__2; ++jr) {
#line 1069 "dtgevc.f"
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
#line 1071 "dtgevc.f"
/* L280: */
#line 1071 "dtgevc.f"
			}
#line 1072 "dtgevc.f"
/* L290: */
#line 1072 "dtgevc.f"
		    }
#line 1073 "dtgevc.f"
		}
/* Computing MAX */
#line 1074 "dtgevc.f"
		d__1 = scale * xmax;
#line 1074 "dtgevc.f"
		xmax = max(d__1,temp);

#line 1076 "dtgevc.f"
		i__1 = nw;
#line 1076 "dtgevc.f"
		for (jw = 1; jw <= i__1; ++jw) {
#line 1077 "dtgevc.f"
		    i__2 = na;
#line 1077 "dtgevc.f"
		    for (ja = 1; ja <= i__2; ++ja) {
#line 1078 "dtgevc.f"
			work[(jw + 1) * *n + j + ja - 1] = sum[ja + (jw << 1) 
				- 3];
#line 1079 "dtgevc.f"
/* L300: */
#line 1079 "dtgevc.f"
		    }
#line 1080 "dtgevc.f"
/* L310: */
#line 1080 "dtgevc.f"
		}

/*              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

#line 1084 "dtgevc.f"
		if (j > 1) {

/*                 Check whether scaling is necessary for sum. */

#line 1088 "dtgevc.f"
		    xscale = 1. / max(1.,xmax);
#line 1089 "dtgevc.f"
		    temp = acoefa * work[j] + bcoefa * work[*n + j];
#line 1090 "dtgevc.f"
		    if (il2by2) {
/* Computing MAX */
#line 1090 "dtgevc.f"
			d__1 = temp, d__2 = acoefa * work[j + 1] + bcoefa * 
				work[*n + j + 1];
#line 1090 "dtgevc.f"
			temp = max(d__1,d__2);
#line 1090 "dtgevc.f"
		    }
/* Computing MAX */
#line 1093 "dtgevc.f"
		    d__1 = max(temp,acoefa);
#line 1093 "dtgevc.f"
		    temp = max(d__1,bcoefa);
#line 1094 "dtgevc.f"
		    if (temp > bignum * xscale) {

#line 1096 "dtgevc.f"
			i__1 = nw - 1;
#line 1096 "dtgevc.f"
			for (jw = 0; jw <= i__1; ++jw) {
#line 1097 "dtgevc.f"
			    i__2 = je;
#line 1097 "dtgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1098 "dtgevc.f"
				work[(jw + 2) * *n + jr] = xscale * work[(jw 
					+ 2) * *n + jr];
#line 1100 "dtgevc.f"
/* L320: */
#line 1100 "dtgevc.f"
			    }
#line 1101 "dtgevc.f"
/* L330: */
#line 1101 "dtgevc.f"
			}
#line 1102 "dtgevc.f"
			xmax *= xscale;
#line 1103 "dtgevc.f"
		    }

/*                 Compute the contributions of the off-diagonals of */
/*                 column j (and j+1, if 2-by-2 block) of A and B to the */
/*                 sums. */


#line 1110 "dtgevc.f"
		    i__1 = na;
#line 1110 "dtgevc.f"
		    for (ja = 1; ja <= i__1; ++ja) {
#line 1111 "dtgevc.f"
			if (ilcplx) {
#line 1112 "dtgevc.f"
			    creala = acoef * work[(*n << 1) + j + ja - 1];
#line 1113 "dtgevc.f"
			    cimaga = acoef * work[*n * 3 + j + ja - 1];
#line 1114 "dtgevc.f"
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1] - 
				    bcoefi * work[*n * 3 + j + ja - 1];
#line 1116 "dtgevc.f"
			    cimagb = bcoefi * work[(*n << 1) + j + ja - 1] + 
				    bcoefr * work[*n * 3 + j + ja - 1];
#line 1118 "dtgevc.f"
			    i__2 = j - 1;
#line 1118 "dtgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1119 "dtgevc.f"
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * s[jr + (j + ja - 1) * s_dim1]
					 + crealb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1122 "dtgevc.f"
				work[*n * 3 + jr] = work[*n * 3 + jr] - 
					cimaga * s[jr + (j + ja - 1) * s_dim1]
					 + cimagb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1125 "dtgevc.f"
/* L340: */
#line 1125 "dtgevc.f"
			    }
#line 1126 "dtgevc.f"
			} else {
#line 1127 "dtgevc.f"
			    creala = acoef * work[(*n << 1) + j + ja - 1];
#line 1128 "dtgevc.f"
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1];
#line 1129 "dtgevc.f"
			    i__2 = j - 1;
#line 1129 "dtgevc.f"
			    for (jr = 1; jr <= i__2; ++jr) {
#line 1130 "dtgevc.f"
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * s[jr + (j + ja - 1) * s_dim1]
					 + crealb * p[jr + (j + ja - 1) * 
					p_dim1];
#line 1133 "dtgevc.f"
/* L350: */
#line 1133 "dtgevc.f"
			    }
#line 1134 "dtgevc.f"
			}
#line 1135 "dtgevc.f"
/* L360: */
#line 1135 "dtgevc.f"
		    }
#line 1136 "dtgevc.f"
		}

#line 1138 "dtgevc.f"
		il2by2 = FALSE_;
#line 1139 "dtgevc.f"
L370:
#line 1139 "dtgevc.f"
		;
#line 1139 "dtgevc.f"
	    }

/*           Copy eigenvector to VR, back transforming if */
/*           HOWMNY='B'. */

#line 1144 "dtgevc.f"
	    ieig -= nw;
#line 1145 "dtgevc.f"
	    if (ilback) {

#line 1147 "dtgevc.f"
		i__1 = nw - 1;
#line 1147 "dtgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1148 "dtgevc.f"
		    i__2 = *n;
#line 1148 "dtgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1149 "dtgevc.f"
			work[(jw + 4) * *n + jr] = work[(jw + 2) * *n + 1] * 
				vr[jr + vr_dim1];
#line 1151 "dtgevc.f"
/* L380: */
#line 1151 "dtgevc.f"
		    }

/*                 A series of compiler directives to defeat */
/*                 vectorization for the next loop */


#line 1157 "dtgevc.f"
		    i__2 = je;
#line 1157 "dtgevc.f"
		    for (jc = 2; jc <= i__2; ++jc) {
#line 1158 "dtgevc.f"
			i__3 = *n;
#line 1158 "dtgevc.f"
			for (jr = 1; jr <= i__3; ++jr) {
#line 1159 "dtgevc.f"
			    work[(jw + 4) * *n + jr] += work[(jw + 2) * *n + 
				    jc] * vr[jr + jc * vr_dim1];
#line 1161 "dtgevc.f"
/* L390: */
#line 1161 "dtgevc.f"
			}
#line 1162 "dtgevc.f"
/* L400: */
#line 1162 "dtgevc.f"
		    }
#line 1163 "dtgevc.f"
/* L410: */
#line 1163 "dtgevc.f"
		}

#line 1165 "dtgevc.f"
		i__1 = nw - 1;
#line 1165 "dtgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1166 "dtgevc.f"
		    i__2 = *n;
#line 1166 "dtgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1167 "dtgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = work[(jw + 4) * *n + 
				jr];
#line 1168 "dtgevc.f"
/* L420: */
#line 1168 "dtgevc.f"
		    }
#line 1169 "dtgevc.f"
/* L430: */
#line 1169 "dtgevc.f"
		}

#line 1171 "dtgevc.f"
		iend = *n;
#line 1172 "dtgevc.f"
	    } else {
#line 1173 "dtgevc.f"
		i__1 = nw - 1;
#line 1173 "dtgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1174 "dtgevc.f"
		    i__2 = *n;
#line 1174 "dtgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1175 "dtgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = work[(jw + 2) * *n + 
				jr];
#line 1176 "dtgevc.f"
/* L440: */
#line 1176 "dtgevc.f"
		    }
#line 1177 "dtgevc.f"
/* L450: */
#line 1177 "dtgevc.f"
		}

#line 1179 "dtgevc.f"
		iend = je;
#line 1180 "dtgevc.f"
	    }

/*           Scale eigenvector */

#line 1184 "dtgevc.f"
	    xmax = 0.;
#line 1185 "dtgevc.f"
	    if (ilcplx) {
#line 1186 "dtgevc.f"
		i__1 = iend;
#line 1186 "dtgevc.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 1187 "dtgevc.f"
		    d__3 = xmax, d__4 = (d__1 = vr[j + ieig * vr_dim1], abs(
			    d__1)) + (d__2 = vr[j + (ieig + 1) * vr_dim1], 
			    abs(d__2));
#line 1187 "dtgevc.f"
		    xmax = max(d__3,d__4);
#line 1189 "dtgevc.f"
/* L460: */
#line 1189 "dtgevc.f"
		}
#line 1190 "dtgevc.f"
	    } else {
#line 1191 "dtgevc.f"
		i__1 = iend;
#line 1191 "dtgevc.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 1192 "dtgevc.f"
		    d__2 = xmax, d__3 = (d__1 = vr[j + ieig * vr_dim1], abs(
			    d__1));
#line 1192 "dtgevc.f"
		    xmax = max(d__2,d__3);
#line 1193 "dtgevc.f"
/* L470: */
#line 1193 "dtgevc.f"
		}
#line 1194 "dtgevc.f"
	    }

#line 1196 "dtgevc.f"
	    if (xmax > safmin) {
#line 1197 "dtgevc.f"
		xscale = 1. / xmax;
#line 1198 "dtgevc.f"
		i__1 = nw - 1;
#line 1198 "dtgevc.f"
		for (jw = 0; jw <= i__1; ++jw) {
#line 1199 "dtgevc.f"
		    i__2 = iend;
#line 1199 "dtgevc.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 1200 "dtgevc.f"
			vr[jr + (ieig + jw) * vr_dim1] = xscale * vr[jr + (
				ieig + jw) * vr_dim1];
#line 1201 "dtgevc.f"
/* L480: */
#line 1201 "dtgevc.f"
		    }
#line 1202 "dtgevc.f"
/* L490: */
#line 1202 "dtgevc.f"
		}
#line 1203 "dtgevc.f"
	    }
#line 1204 "dtgevc.f"
L500:
#line 1204 "dtgevc.f"
	    ;
#line 1204 "dtgevc.f"
	}
#line 1205 "dtgevc.f"
    }

#line 1207 "dtgevc.f"
    return 0;

/*     End of DTGEVC */

} /* dtgevc_ */

