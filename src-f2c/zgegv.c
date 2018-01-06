#line 1 "zgegv.f"
/* zgegv.f -- translated by f2c (version 20100827).
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

#line 1 "zgegv.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b29 = 1.;

/* > \brief <b> ZGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgegv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgegv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgegv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, */
/*                         VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine ZGGEV. */
/* > */
/* > ZGEGV computes the eigenvalues and, optionally, the left and/or right */
/* > eigenvectors of a complex matrix pair (A,B). */
/* > Given two square matrices A and B, */
/* > the generalized nonsymmetric eigenvalue problem (GNEP) is to find the */
/* > eigenvalues lambda and corresponding (non-zero) eigenvectors x such */
/* > that */
/* >    A*x = lambda*B*x. */
/* > */
/* > An alternate form is to find the eigenvalues mu and corresponding */
/* > eigenvectors y such that */
/* >    mu*A*y = B*y. */
/* > */
/* > These two forms are equivalent with mu = 1/lambda and x = y if */
/* > neither lambda nor mu is zero.  In order to deal with the case that */
/* > lambda or mu is zero or small, two values alpha and beta are returned */
/* > for each eigenvalue, such that lambda = alpha/beta and */
/* > mu = beta/alpha. */
/* > */
/* > The vectors x and y in the above equations are right eigenvectors of */
/* > the matrix pair (A,B).  Vectors u and v satisfying */
/* >    u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B */
/* > are left eigenvectors of (A,B). */
/* > */
/* > Note: this routine performs "full balancing" on A and B */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBVL */
/* > \verbatim */
/* >          JOBVL is CHARACTER*1 */
/* >          = 'N':  do not compute the left generalized eigenvectors; */
/* >          = 'V':  compute the left generalized eigenvectors (returned */
/* >                  in VL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* >          JOBVR is CHARACTER*1 */
/* >          = 'N':  do not compute the right generalized eigenvectors; */
/* >          = 'V':  compute the right generalized eigenvectors (returned */
/* >                  in VR). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A, B, VL, and VR.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the matrix A. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit A */
/* >          contains the Schur form of A from the generalized Schur */
/* >          factorization of the pair (A,B) after balancing.  If no */
/* >          eigenvectors were computed, then only the diagonal elements */
/* >          of the Schur form will be correct.  See ZGGHRD and ZHGEQZ */
/* >          for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
/* >          On entry, the matrix B. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the */
/* >          upper triangular matrix obtained from B in the generalized */
/* >          Schur factorization of the pair (A,B) after balancing. */
/* >          If no eigenvectors were computed, then only the diagonal */
/* >          elements of B will be correct.  See ZGGHRD and ZHGEQZ for */
/* >          details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 array, dimension (N) */
/* >          The complex scalars alpha that define the eigenvalues of */
/* >          GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
/* >          The complex scalars beta that define the eigenvalues of GNEP. */
/* > */
/* >          Together, the quantities alpha = ALPHA(j) and beta = BETA(j) */
/* >          represent the j-th eigenvalue of the matrix pair (A,B), in */
/* >          one of the forms lambda = alpha/beta or mu = beta/alpha. */
/* >          Since either lambda or mu may overflow, they should not, */
/* >          in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored */
/* >          in the columns of VL, in the same order as their eigenvalues. */
/* >          Each eigenvector is scaled so that its largest component has */
/* >          abs(real part) + abs(imag. part) = 1, except for eigenvectors */
/* >          corresponding to an eigenvalue with alpha = beta = 0, which */
/* >          are set to zero. */
/* >          Not referenced if JOBVL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the matrix VL. LDVL >= 1, and */
/* >          if JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* >          VR is COMPLEX*16 array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors x(j) are stored */
/* >          in the columns of VR, in the same order as their eigenvalues. */
/* >          Each eigenvector is scaled so that its largest component has */
/* >          abs(real part) + abs(imag. part) = 1, except for eigenvectors */
/* >          corresponding to an eigenvalue with alpha = beta = 0, which */
/* >          are set to zero. */
/* >          Not referenced if JOBVR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the matrix VR. LDVR >= 1, and */
/* >          if JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,2*N). */
/* >          For good performance, LWORK must generally be larger. */
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for ZGEQRF, ZUNMQR, and ZUNGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for ZGEQRF, ZUNMQR, and ZUNGQR; */
/* >          The optimal LWORK is  MAX( 2*N, N*(NB+1) ). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          =1,...,N: */
/* >                The QZ iteration failed.  No eigenvectors have been */
/* >                calculated, but ALPHA(j) and BETA(j) should be */
/* >                correct for j=INFO+1,...,N. */
/* >          > N:  errors that usually indicate LAPACK problems: */
/* >                =N+1: error return from ZGGBAL */
/* >                =N+2: error return from ZGEQRF */
/* >                =N+3: error return from ZUNMQR */
/* >                =N+4: error return from ZUNGQR */
/* >                =N+5: error return from ZGGHRD */
/* >                =N+6: error return from ZHGEQZ (other than failed */
/* >                                               iteration) */
/* >                =N+7: error return from ZTGEVC */
/* >                =N+8: error return from ZGGBAK (computing VL) */
/* >                =N+9: error return from ZGGBAK (computing VR) */
/* >                =N+10: error return from ZLASCL (various calls) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16GEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Balancing */
/* >  --------- */
/* > */
/* >  This driver calls ZGGBAL to both permute and scale rows and columns */
/* >  of A and B.  The permutations PL and PR are chosen so that PL*A*PR */
/* >  and PL*B*R will be upper triangular except for the diagonal blocks */
/* >  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as */
/* >  possible.  The diagonal scaling matrices DL and DR are chosen so */
/* >  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to */
/* >  one (except for the elements that start out zero.) */
/* > */
/* >  After the eigenvalues and eigenvectors of the balanced matrices */
/* >  have been computed, ZGGBAK transforms the eigenvectors back to what */
/* >  they would have been (in perfect arithmetic) if they had not been */
/* >  balanced. */
/* > */
/* >  Contents of A and B on Exit */
/* >  -------- -- - --- - -- ---- */
/* > */
/* >  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or */
/* >  both), then on exit the arrays A and B will contain the complex Schur */
/* >  form[*] of the "balanced" versions of A and B.  If no eigenvectors */
/* >  are computed, then only the diagonal blocks will be correct. */
/* > */
/* >  [*] In other words, upper triangular form. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info, ftnlen jobvl_len, ftnlen 
	jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer jc, nb, in, jr, nb1, nb2, nb3, ihi, ilo;
    static doublereal eps;
    static logical ilv;
    static doublereal absb, anrm, bnrm;
    static integer itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static integer lopt;
    static doublereal anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ileft, iinfo, icols, iwork, irows;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal salfai;
    extern /* Subroutine */ int zggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen, ftnlen), zggbal_(char *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static doublereal salfar, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer ijobvl, iright;
    static logical ilimit;
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static integer lwkmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), ztgevc_(
	    char *, char *, logical *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublecomplex *,
	     doublereal *, integer *, ftnlen, ftnlen), zhgeqz_(char *, char *,
	     char *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer irwork, lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmqr_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 347 "zgegv.f"
    /* Parameter adjustments */
#line 347 "zgegv.f"
    a_dim1 = *lda;
#line 347 "zgegv.f"
    a_offset = 1 + a_dim1;
#line 347 "zgegv.f"
    a -= a_offset;
#line 347 "zgegv.f"
    b_dim1 = *ldb;
#line 347 "zgegv.f"
    b_offset = 1 + b_dim1;
#line 347 "zgegv.f"
    b -= b_offset;
#line 347 "zgegv.f"
    --alpha;
#line 347 "zgegv.f"
    --beta;
#line 347 "zgegv.f"
    vl_dim1 = *ldvl;
#line 347 "zgegv.f"
    vl_offset = 1 + vl_dim1;
#line 347 "zgegv.f"
    vl -= vl_offset;
#line 347 "zgegv.f"
    vr_dim1 = *ldvr;
#line 347 "zgegv.f"
    vr_offset = 1 + vr_dim1;
#line 347 "zgegv.f"
    vr -= vr_offset;
#line 347 "zgegv.f"
    --work;
#line 347 "zgegv.f"
    --rwork;
#line 347 "zgegv.f"

#line 347 "zgegv.f"
    /* Function Body */
#line 347 "zgegv.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 348 "zgegv.f"
	ijobvl = 1;
#line 349 "zgegv.f"
	ilvl = FALSE_;
#line 350 "zgegv.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 351 "zgegv.f"
	ijobvl = 2;
#line 352 "zgegv.f"
	ilvl = TRUE_;
#line 353 "zgegv.f"
    } else {
#line 354 "zgegv.f"
	ijobvl = -1;
#line 355 "zgegv.f"
	ilvl = FALSE_;
#line 356 "zgegv.f"
    }

#line 358 "zgegv.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 359 "zgegv.f"
	ijobvr = 1;
#line 360 "zgegv.f"
	ilvr = FALSE_;
#line 361 "zgegv.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 362 "zgegv.f"
	ijobvr = 2;
#line 363 "zgegv.f"
	ilvr = TRUE_;
#line 364 "zgegv.f"
    } else {
#line 365 "zgegv.f"
	ijobvr = -1;
#line 366 "zgegv.f"
	ilvr = FALSE_;
#line 367 "zgegv.f"
    }
#line 368 "zgegv.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

/* Computing MAX */
#line 372 "zgegv.f"
    i__1 = *n << 1;
#line 372 "zgegv.f"
    lwkmin = max(i__1,1);
#line 373 "zgegv.f"
    lwkopt = lwkmin;
#line 374 "zgegv.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 375 "zgegv.f"
    lquery = *lwork == -1;
#line 376 "zgegv.f"
    *info = 0;
#line 377 "zgegv.f"
    if (ijobvl <= 0) {
#line 378 "zgegv.f"
	*info = -1;
#line 379 "zgegv.f"
    } else if (ijobvr <= 0) {
#line 380 "zgegv.f"
	*info = -2;
#line 381 "zgegv.f"
    } else if (*n < 0) {
#line 382 "zgegv.f"
	*info = -3;
#line 383 "zgegv.f"
    } else if (*lda < max(1,*n)) {
#line 384 "zgegv.f"
	*info = -5;
#line 385 "zgegv.f"
    } else if (*ldb < max(1,*n)) {
#line 386 "zgegv.f"
	*info = -7;
#line 387 "zgegv.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 388 "zgegv.f"
	*info = -11;
#line 389 "zgegv.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 390 "zgegv.f"
	*info = -13;
#line 391 "zgegv.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 392 "zgegv.f"
	*info = -15;
#line 393 "zgegv.f"
    }

#line 395 "zgegv.f"
    if (*info == 0) {
#line 396 "zgegv.f"
	nb1 = ilaenv_(&c__1, "ZGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 397 "zgegv.f"
	nb2 = ilaenv_(&c__1, "ZUNMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 398 "zgegv.f"
	nb3 = ilaenv_(&c__1, "ZUNGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 399 "zgegv.f"
	i__1 = max(nb1,nb2);
#line 399 "zgegv.f"
	nb = max(i__1,nb3);
/* Computing MAX */
#line 400 "zgegv.f"
	i__1 = *n << 1, i__2 = *n * (nb + 1);
#line 400 "zgegv.f"
	lopt = max(i__1,i__2);
#line 401 "zgegv.f"
	work[1].r = (doublereal) lopt, work[1].i = 0.;
#line 402 "zgegv.f"
    }

#line 404 "zgegv.f"
    if (*info != 0) {
#line 405 "zgegv.f"
	i__1 = -(*info);
#line 405 "zgegv.f"
	xerbla_("ZGEGV ", &i__1, (ftnlen)6);
#line 406 "zgegv.f"
	return 0;
#line 407 "zgegv.f"
    } else if (lquery) {
#line 408 "zgegv.f"
	return 0;
#line 409 "zgegv.f"
    }

/*     Quick return if possible */

#line 413 "zgegv.f"
    if (*n == 0) {
#line 413 "zgegv.f"
	return 0;
#line 413 "zgegv.f"
    }

/*     Get machine constants */

#line 418 "zgegv.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 419 "zgegv.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 420 "zgegv.f"
    safmin += safmin;
#line 421 "zgegv.f"
    safmax = 1. / safmin;

/*     Scale A */

#line 425 "zgegv.f"
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 426 "zgegv.f"
    anrm1 = anrm;
#line 427 "zgegv.f"
    anrm2 = 1.;
#line 428 "zgegv.f"
    if (anrm < 1.) {
#line 429 "zgegv.f"
	if (safmax * anrm < 1.) {
#line 430 "zgegv.f"
	    anrm1 = safmin;
#line 431 "zgegv.f"
	    anrm2 = safmax * anrm;
#line 432 "zgegv.f"
	}
#line 433 "zgegv.f"
    }

#line 435 "zgegv.f"
    if (anrm > 0.) {
#line 436 "zgegv.f"
	zlascl_("G", &c_n1, &c_n1, &anrm, &c_b29, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 437 "zgegv.f"
	if (iinfo != 0) {
#line 438 "zgegv.f"
	    *info = *n + 10;
#line 439 "zgegv.f"
	    return 0;
#line 440 "zgegv.f"
	}
#line 441 "zgegv.f"
    }

/*     Scale B */

#line 445 "zgegv.f"
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);
#line 446 "zgegv.f"
    bnrm1 = bnrm;
#line 447 "zgegv.f"
    bnrm2 = 1.;
#line 448 "zgegv.f"
    if (bnrm < 1.) {
#line 449 "zgegv.f"
	if (safmax * bnrm < 1.) {
#line 450 "zgegv.f"
	    bnrm1 = safmin;
#line 451 "zgegv.f"
	    bnrm2 = safmax * bnrm;
#line 452 "zgegv.f"
	}
#line 453 "zgegv.f"
    }

#line 455 "zgegv.f"
    if (bnrm > 0.) {
#line 456 "zgegv.f"
	zlascl_("G", &c_n1, &c_n1, &bnrm, &c_b29, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 457 "zgegv.f"
	if (iinfo != 0) {
#line 458 "zgegv.f"
	    *info = *n + 10;
#line 459 "zgegv.f"
	    return 0;
#line 460 "zgegv.f"
	}
#line 461 "zgegv.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     Also "balance" the matrix. */

#line 466 "zgegv.f"
    ileft = 1;
#line 467 "zgegv.f"
    iright = *n + 1;
#line 468 "zgegv.f"
    irwork = iright + *n;
#line 469 "zgegv.f"
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[
	    ileft], &rwork[iright], &rwork[irwork], &iinfo, (ftnlen)1);
#line 471 "zgegv.f"
    if (iinfo != 0) {
#line 472 "zgegv.f"
	*info = *n + 1;
#line 473 "zgegv.f"
	goto L80;
#line 474 "zgegv.f"
    }

/*     Reduce B to triangular form, and initialize VL and/or VR */

#line 478 "zgegv.f"
    irows = ihi + 1 - ilo;
#line 479 "zgegv.f"
    if (ilv) {
#line 480 "zgegv.f"
	icols = *n + 1 - ilo;
#line 481 "zgegv.f"
    } else {
#line 482 "zgegv.f"
	icols = irows;
#line 483 "zgegv.f"
    }
#line 484 "zgegv.f"
    itau = 1;
#line 485 "zgegv.f"
    iwork = itau + irows;
#line 486 "zgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 486 "zgegv.f"
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 488 "zgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 488 "zgegv.f"
	i__3 = iwork;
#line 488 "zgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 488 "zgegv.f"
	lwkopt = max(i__1,i__2);
#line 488 "zgegv.f"
    }
#line 490 "zgegv.f"
    if (iinfo != 0) {
#line 491 "zgegv.f"
	*info = *n + 2;
#line 492 "zgegv.f"
	goto L80;
#line 493 "zgegv.f"
    }

#line 495 "zgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 495 "zgegv.f"
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 498 "zgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 498 "zgegv.f"
	i__3 = iwork;
#line 498 "zgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 498 "zgegv.f"
	lwkopt = max(i__1,i__2);
#line 498 "zgegv.f"
    }
#line 500 "zgegv.f"
    if (iinfo != 0) {
#line 501 "zgegv.f"
	*info = *n + 3;
#line 502 "zgegv.f"
	goto L80;
#line 503 "zgegv.f"
    }

#line 505 "zgegv.f"
    if (ilvl) {
#line 506 "zgegv.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl, (ftnlen)4);
#line 507 "zgegv.f"
	i__1 = irows - 1;
#line 507 "zgegv.f"
	i__2 = irows - 1;
#line 507 "zgegv.f"
	zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[ilo + 
		1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 509 "zgegv.f"
	i__1 = *lwork + 1 - iwork;
#line 509 "zgegv.f"
	zungqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwork], &i__1, &iinfo);
#line 512 "zgegv.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 512 "zgegv.f"
	    i__3 = iwork;
#line 512 "zgegv.f"
	    i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 512 "zgegv.f"
	    lwkopt = max(i__1,i__2);
#line 512 "zgegv.f"
	}
#line 514 "zgegv.f"
	if (iinfo != 0) {
#line 515 "zgegv.f"
	    *info = *n + 4;
#line 516 "zgegv.f"
	    goto L80;
#line 517 "zgegv.f"
	}
#line 518 "zgegv.f"
    }

#line 520 "zgegv.f"
    if (ilvr) {
#line 520 "zgegv.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr, (ftnlen)4);
#line 520 "zgegv.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 525 "zgegv.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 529 "zgegv.f"
	zgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 531 "zgegv.f"
    } else {
#line 532 "zgegv.f"
	zgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 534 "zgegv.f"
    }
#line 535 "zgegv.f"
    if (iinfo != 0) {
#line 536 "zgegv.f"
	*info = *n + 5;
#line 537 "zgegv.f"
	goto L80;
#line 538 "zgegv.f"
    }

/*     Perform QZ algorithm */

#line 542 "zgegv.f"
    iwork = itau;
#line 543 "zgegv.f"
    if (ilv) {
#line 544 "zgegv.f"
	*(unsigned char *)chtemp = 'S';
#line 545 "zgegv.f"
    } else {
#line 546 "zgegv.f"
	*(unsigned char *)chtemp = 'E';
#line 547 "zgegv.f"
    }
#line 548 "zgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 548 "zgegv.f"
    zhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[
	    vr_offset], ldvr, &work[iwork], &i__1, &rwork[irwork], &iinfo, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 551 "zgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 551 "zgegv.f"
	i__3 = iwork;
#line 551 "zgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[i__3].r + iwork - 1;
#line 551 "zgegv.f"
	lwkopt = max(i__1,i__2);
#line 551 "zgegv.f"
    }
#line 553 "zgegv.f"
    if (iinfo != 0) {
#line 554 "zgegv.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 555 "zgegv.f"
	    *info = iinfo;
#line 556 "zgegv.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 557 "zgegv.f"
	    *info = iinfo - *n;
#line 558 "zgegv.f"
	} else {
#line 559 "zgegv.f"
	    *info = *n + 6;
#line 560 "zgegv.f"
	}
#line 561 "zgegv.f"
	goto L80;
#line 562 "zgegv.f"
    }

#line 564 "zgegv.f"
    if (ilv) {

/*        Compute Eigenvectors */

#line 568 "zgegv.f"
	if (ilvl) {
#line 569 "zgegv.f"
	    if (ilvr) {
#line 570 "zgegv.f"
		*(unsigned char *)chtemp = 'B';
#line 571 "zgegv.f"
	    } else {
#line 572 "zgegv.f"
		*(unsigned char *)chtemp = 'L';
#line 573 "zgegv.f"
	    }
#line 574 "zgegv.f"
	} else {
#line 575 "zgegv.f"
	    *(unsigned char *)chtemp = 'R';
#line 576 "zgegv.f"
	}

#line 578 "zgegv.f"
	ztgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwork], &rwork[irwork], &iinfo, (ftnlen)1, (ftnlen)1);
#line 581 "zgegv.f"
	if (iinfo != 0) {
#line 582 "zgegv.f"
	    *info = *n + 7;
#line 583 "zgegv.f"
	    goto L80;
#line 584 "zgegv.f"
	}

/*        Undo balancing on VL and VR, rescale */

#line 588 "zgegv.f"
	if (ilvl) {
#line 589 "zgegv.f"
	    zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vl[vl_offset], ldvl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 591 "zgegv.f"
	    if (iinfo != 0) {
#line 592 "zgegv.f"
		*info = *n + 8;
#line 593 "zgegv.f"
		goto L80;
#line 594 "zgegv.f"
	    }
#line 595 "zgegv.f"
	    i__1 = *n;
#line 595 "zgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 596 "zgegv.f"
		temp = 0.;
#line 597 "zgegv.f"
		i__2 = *n;
#line 597 "zgegv.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 598 "zgegv.f"
		    i__3 = jr + jc * vl_dim1;
#line 598 "zgegv.f"
		    d__3 = temp, d__4 = (d__1 = vl[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vl[jr + jc * vl_dim1]), abs(d__2));
#line 598 "zgegv.f"
		    temp = max(d__3,d__4);
#line 599 "zgegv.f"
/* L10: */
#line 599 "zgegv.f"
		}
#line 600 "zgegv.f"
		if (temp < safmin) {
#line 600 "zgegv.f"
		    goto L30;
#line 600 "zgegv.f"
		}
#line 602 "zgegv.f"
		temp = 1. / temp;
#line 603 "zgegv.f"
		i__2 = *n;
#line 603 "zgegv.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 604 "zgegv.f"
		    i__3 = jr + jc * vl_dim1;
#line 604 "zgegv.f"
		    i__4 = jr + jc * vl_dim1;
#line 604 "zgegv.f"
		    z__1.r = temp * vl[i__4].r, z__1.i = temp * vl[i__4].i;
#line 604 "zgegv.f"
		    vl[i__3].r = z__1.r, vl[i__3].i = z__1.i;
#line 605 "zgegv.f"
/* L20: */
#line 605 "zgegv.f"
		}
#line 606 "zgegv.f"
L30:
#line 606 "zgegv.f"
		;
#line 606 "zgegv.f"
	    }
#line 607 "zgegv.f"
	}
#line 608 "zgegv.f"
	if (ilvr) {
#line 609 "zgegv.f"
	    zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
		     &vr[vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 611 "zgegv.f"
	    if (iinfo != 0) {
#line 612 "zgegv.f"
		*info = *n + 9;
#line 613 "zgegv.f"
		goto L80;
#line 614 "zgegv.f"
	    }
#line 615 "zgegv.f"
	    i__1 = *n;
#line 615 "zgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 616 "zgegv.f"
		temp = 0.;
#line 617 "zgegv.f"
		i__2 = *n;
#line 617 "zgegv.f"
		for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 618 "zgegv.f"
		    i__3 = jr + jc * vr_dim1;
#line 618 "zgegv.f"
		    d__3 = temp, d__4 = (d__1 = vr[i__3].r, abs(d__1)) + (
			    d__2 = d_imag(&vr[jr + jc * vr_dim1]), abs(d__2));
#line 618 "zgegv.f"
		    temp = max(d__3,d__4);
#line 619 "zgegv.f"
/* L40: */
#line 619 "zgegv.f"
		}
#line 620 "zgegv.f"
		if (temp < safmin) {
#line 620 "zgegv.f"
		    goto L60;
#line 620 "zgegv.f"
		}
#line 622 "zgegv.f"
		temp = 1. / temp;
#line 623 "zgegv.f"
		i__2 = *n;
#line 623 "zgegv.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 624 "zgegv.f"
		    i__3 = jr + jc * vr_dim1;
#line 624 "zgegv.f"
		    i__4 = jr + jc * vr_dim1;
#line 624 "zgegv.f"
		    z__1.r = temp * vr[i__4].r, z__1.i = temp * vr[i__4].i;
#line 624 "zgegv.f"
		    vr[i__3].r = z__1.r, vr[i__3].i = z__1.i;
#line 625 "zgegv.f"
/* L50: */
#line 625 "zgegv.f"
		}
#line 626 "zgegv.f"
L60:
#line 626 "zgegv.f"
		;
#line 626 "zgegv.f"
	    }
#line 627 "zgegv.f"
	}

/*        End of eigenvector calculation */

#line 631 "zgegv.f"
    }

/*     Undo scaling in alpha, beta */

/*     Note: this does not give the alpha and beta for the unscaled */
/*     problem. */

/*     Un-scaling is limited to avoid underflow in alpha and beta */
/*     if they are significant. */

#line 641 "zgegv.f"
    i__1 = *n;
#line 641 "zgegv.f"
    for (jc = 1; jc <= i__1; ++jc) {
#line 642 "zgegv.f"
	i__2 = jc;
#line 642 "zgegv.f"
	absar = (d__1 = alpha[i__2].r, abs(d__1));
#line 643 "zgegv.f"
	absai = (d__1 = d_imag(&alpha[jc]), abs(d__1));
#line 644 "zgegv.f"
	i__2 = jc;
#line 644 "zgegv.f"
	absb = (d__1 = beta[i__2].r, abs(d__1));
#line 645 "zgegv.f"
	i__2 = jc;
#line 645 "zgegv.f"
	salfar = anrm * alpha[i__2].r;
#line 646 "zgegv.f"
	salfai = anrm * d_imag(&alpha[jc]);
#line 647 "zgegv.f"
	i__2 = jc;
#line 647 "zgegv.f"
	sbeta = bnrm * beta[i__2].r;
#line 648 "zgegv.f"
	ilimit = FALSE_;
#line 649 "zgegv.f"
	scale = 1.;

/*        Check for significant underflow in imaginary part of ALPHA */

/* Computing MAX */
#line 653 "zgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 653 "zgegv.f"
	if (abs(salfai) < safmin && absai >= max(d__1,d__2)) {
#line 655 "zgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
#line 656 "zgegv.f"
	    d__1 = safmin, d__2 = anrm2 * absai;
#line 656 "zgegv.f"
	    scale = safmin / anrm1 / max(d__1,d__2);
#line 657 "zgegv.f"
	}

/*        Check for significant underflow in real part of ALPHA */

/* Computing MAX */
#line 661 "zgegv.f"
	d__1 = safmin, d__2 = eps * absai, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 661 "zgegv.f"
	if (abs(salfar) < safmin && absar >= max(d__1,d__2)) {
#line 663 "zgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 664 "zgegv.f"
	    d__3 = safmin, d__4 = anrm2 * absar;
#line 664 "zgegv.f"
	    d__1 = scale, d__2 = safmin / anrm1 / max(d__3,d__4);
#line 664 "zgegv.f"
	    scale = max(d__1,d__2);
#line 666 "zgegv.f"
	}

/*        Check for significant underflow in BETA */

/* Computing MAX */
#line 670 "zgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absai;
#line 670 "zgegv.f"
	if (abs(sbeta) < safmin && absb >= max(d__1,d__2)) {
#line 672 "zgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 673 "zgegv.f"
	    d__3 = safmin, d__4 = bnrm2 * absb;
#line 673 "zgegv.f"
	    d__1 = scale, d__2 = safmin / bnrm1 / max(d__3,d__4);
#line 673 "zgegv.f"
	    scale = max(d__1,d__2);
#line 675 "zgegv.f"
	}

/*        Check for possible overflow when limiting scaling */

#line 679 "zgegv.f"
	if (ilimit) {
/* Computing MAX */
#line 680 "zgegv.f"
	    d__1 = abs(salfar), d__2 = abs(salfai), d__1 = max(d__1,d__2), 
		    d__2 = abs(sbeta);
#line 680 "zgegv.f"
	    temp = scale * safmin * max(d__1,d__2);
#line 682 "zgegv.f"
	    if (temp > 1.) {
#line 682 "zgegv.f"
		scale /= temp;
#line 682 "zgegv.f"
	    }
#line 684 "zgegv.f"
	    if (scale < 1.) {
#line 684 "zgegv.f"
		ilimit = FALSE_;
#line 684 "zgegv.f"
	    }
#line 686 "zgegv.f"
	}

/*        Recompute un-scaled ALPHA, BETA if necessary. */

#line 690 "zgegv.f"
	if (ilimit) {
#line 691 "zgegv.f"
	    i__2 = jc;
#line 691 "zgegv.f"
	    salfar = scale * alpha[i__2].r * anrm;
#line 692 "zgegv.f"
	    salfai = scale * d_imag(&alpha[jc]) * anrm;
#line 693 "zgegv.f"
	    i__2 = jc;
#line 693 "zgegv.f"
	    z__2.r = scale * beta[i__2].r, z__2.i = scale * beta[i__2].i;
#line 693 "zgegv.f"
	    z__1.r = bnrm * z__2.r, z__1.i = bnrm * z__2.i;
#line 693 "zgegv.f"
	    sbeta = z__1.r;
#line 694 "zgegv.f"
	}
#line 695 "zgegv.f"
	i__2 = jc;
#line 695 "zgegv.f"
	z__1.r = salfar, z__1.i = salfai;
#line 695 "zgegv.f"
	alpha[i__2].r = z__1.r, alpha[i__2].i = z__1.i;
#line 696 "zgegv.f"
	i__2 = jc;
#line 696 "zgegv.f"
	beta[i__2].r = sbeta, beta[i__2].i = 0.;
#line 697 "zgegv.f"
/* L70: */
#line 697 "zgegv.f"
    }

#line 699 "zgegv.f"
L80:
#line 700 "zgegv.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

#line 702 "zgegv.f"
    return 0;

/*     End of ZGEGV */

} /* zgegv_ */

