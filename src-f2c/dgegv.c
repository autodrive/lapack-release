#line 1 "dgegv.f"
/* dgegv.f -- translated by f2c (version 20100827).
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

#line 1 "dgegv.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b27 = 1.;
static doublereal c_b38 = 0.;

/* > \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgegv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgegv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgegv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, */
/*                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGGEV. */
/* > */
/* > DGEGV computes the eigenvalues and, optionally, the left and/or right */
/* > eigenvectors of a real matrix pair (A,B). */
/* > Given two square matrices A and B, */
/* > the generalized nonsymmetric eigenvalue problem (GNEP) is to find the */
/* > eigenvalues lambda and corresponding (non-zero) eigenvectors x such */
/* > that */
/* > */
/* >    A*x = lambda*B*x. */
/* > */
/* > An alternate form is to find the eigenvalues mu and corresponding */
/* > eigenvectors y such that */
/* > */
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
/* > */
/* >    u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B */
/* > */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
/* >          On entry, the matrix A. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit A */
/* >          contains the real Schur form of A from the generalized Schur */
/* >          factorization of the pair (A,B) after balancing. */
/* >          If no eigenvectors were computed, then only the diagonal */
/* >          blocks from the Schur form will be correct.  See DGGHRD and */
/* >          DHGEQZ for details. */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, the matrix B. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the */
/* >          upper triangular matrix obtained from B in the generalized */
/* >          Schur factorization of the pair (A,B) after balancing. */
/* >          If no eigenvectors were computed, then only those elements of */
/* >          B corresponding to the diagonal blocks from the Schur form of */
/* >          A will be correct.  See DGGHRD and DHGEQZ for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* >          ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* >          The real parts of each scalar alpha defining an eigenvalue of */
/* >          GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* >          The imaginary parts of each scalar alpha defining an */
/* >          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th */
/* >          eigenvalue is real; if positive, then the j-th and */
/* >          (j+1)-st eigenvalues are a complex conjugate pair, with */
/* >          ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
/* >          The scalars beta that define the eigenvalues of GNEP. */
/* > */
/* >          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* >          beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* >          pair (A,B), in one of the forms lambda = alpha/beta or */
/* >          mu = beta/alpha.  Since either lambda or mu may overflow, */
/* >          they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,N) */
/* >          If JOBVL = 'V', the left eigenvectors u(j) are stored */
/* >          in the columns of VL, in the same order as their eigenvalues. */
/* >          If the j-th eigenvalue is real, then u(j) = VL(:,j). */
/* >          If the j-th and (j+1)-st eigenvalues form a complex conjugate */
/* >          pair, then */
/* >             u(j) = VL(:,j) + i*VL(:,j+1) */
/* >          and */
/* >            u(j+1) = VL(:,j) - i*VL(:,j+1). */
/* > */
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
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,N) */
/* >          If JOBVR = 'V', the right eigenvectors x(j) are stored */
/* >          in the columns of VR, in the same order as their eigenvalues. */
/* >          If the j-th eigenvalue is real, then x(j) = VR(:,j). */
/* >          If the j-th and (j+1)-st eigenvalues form a complex conjugate */
/* >          pair, then */
/* >            x(j) = VR(:,j) + i*VR(:,j+1) */
/* >          and */
/* >            x(j+1) = VR(:,j) - i*VR(:,j+1). */
/* > */
/* >          Each eigenvector is scaled so that its largest component has */
/* >          abs(real part) + abs(imag. part) = 1, except for eigenvalues */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,8*N). */
/* >          For good performance, LWORK must generally be larger. */
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR; */
/* >          The optimal LWORK is: */
/* >              2*N + MAX( 6*N, N*(NB+1) ). */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          = 1,...,N: */
/* >                The QZ iteration failed.  No eigenvectors have been */
/* >                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) */
/* >                should be correct for j=INFO+1,...,N. */
/* >          > N:  errors that usually indicate LAPACK problems: */
/* >                =N+1: error return from DGGBAL */
/* >                =N+2: error return from DGEQRF */
/* >                =N+3: error return from DORMQR */
/* >                =N+4: error return from DORGQR */
/* >                =N+5: error return from DGGHRD */
/* >                =N+6: error return from DHGEQZ (other than failed */
/* >                                                iteration) */
/* >                =N+7: error return from DTGEVC */
/* >                =N+8: error return from DGGBAK (computing VL) */
/* >                =N+9: error return from DGGBAK (computing VR) */
/* >                =N+10: error return from DLASCL (various calls) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Balancing */
/* >  --------- */
/* > */
/* >  This driver calls DGGBAL to both permute and scale rows and columns */
/* >  of A and B.  The permutations PL and PR are chosen so that PL*A*PR */
/* >  and PL*B*R will be upper triangular except for the diagonal blocks */
/* >  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as */
/* >  possible.  The diagonal scaling matrices DL and DR are chosen so */
/* >  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to */
/* >  one (except for the elements that start out zero.) */
/* > */
/* >  After the eigenvalues and eigenvectors of the balanced matrices */
/* >  have been computed, DGGBAK transforms the eigenvectors back to what */
/* >  they would have been (in perfect arithmetic) if they had not been */
/* >  balanced. */
/* > */
/* >  Contents of A and B on Exit */
/* >  -------- -- - --- - -- ---- */
/* > */
/* >  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or */
/* >  both), then on exit the arrays A and B will contain the real Schur */
/* >  form[*] of the "balanced" versions of A and B.  If no eigenvectors */
/* >  are computed, then only the diagonal blocks will be correct. */
/* > */
/* >  [*] See DHGEQZ, DGEGS, or read the book "Matrix Computations", */
/* >      by Golub & van Loan, pub. by Johns Hopkins U. Press. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info, ftnlen jobvl_len, ftnlen jobvr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

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
    extern /* Subroutine */ int dggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal salfai;
    extern /* Subroutine */ int dgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlascl_(char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal salfar;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dtgevc_(char *, char *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer ijobvl, iright;
    static logical ilimit;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal onepls;
    static integer lwkmin;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer lwkopt;
    static logical lquery;


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
/*     .. Executable Statements .. */

/*     Decode the input arguments */

#line 360 "dgegv.f"
    /* Parameter adjustments */
#line 360 "dgegv.f"
    a_dim1 = *lda;
#line 360 "dgegv.f"
    a_offset = 1 + a_dim1;
#line 360 "dgegv.f"
    a -= a_offset;
#line 360 "dgegv.f"
    b_dim1 = *ldb;
#line 360 "dgegv.f"
    b_offset = 1 + b_dim1;
#line 360 "dgegv.f"
    b -= b_offset;
#line 360 "dgegv.f"
    --alphar;
#line 360 "dgegv.f"
    --alphai;
#line 360 "dgegv.f"
    --beta;
#line 360 "dgegv.f"
    vl_dim1 = *ldvl;
#line 360 "dgegv.f"
    vl_offset = 1 + vl_dim1;
#line 360 "dgegv.f"
    vl -= vl_offset;
#line 360 "dgegv.f"
    vr_dim1 = *ldvr;
#line 360 "dgegv.f"
    vr_offset = 1 + vr_dim1;
#line 360 "dgegv.f"
    vr -= vr_offset;
#line 360 "dgegv.f"
    --work;
#line 360 "dgegv.f"

#line 360 "dgegv.f"
    /* Function Body */
#line 360 "dgegv.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 361 "dgegv.f"
	ijobvl = 1;
#line 362 "dgegv.f"
	ilvl = FALSE_;
#line 363 "dgegv.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 364 "dgegv.f"
	ijobvl = 2;
#line 365 "dgegv.f"
	ilvl = TRUE_;
#line 366 "dgegv.f"
    } else {
#line 367 "dgegv.f"
	ijobvl = -1;
#line 368 "dgegv.f"
	ilvl = FALSE_;
#line 369 "dgegv.f"
    }

#line 371 "dgegv.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 372 "dgegv.f"
	ijobvr = 1;
#line 373 "dgegv.f"
	ilvr = FALSE_;
#line 374 "dgegv.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 375 "dgegv.f"
	ijobvr = 2;
#line 376 "dgegv.f"
	ilvr = TRUE_;
#line 377 "dgegv.f"
    } else {
#line 378 "dgegv.f"
	ijobvr = -1;
#line 379 "dgegv.f"
	ilvr = FALSE_;
#line 380 "dgegv.f"
    }
#line 381 "dgegv.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

/* Computing MAX */
#line 385 "dgegv.f"
    i__1 = *n << 3;
#line 385 "dgegv.f"
    lwkmin = max(i__1,1);
#line 386 "dgegv.f"
    lwkopt = lwkmin;
#line 387 "dgegv.f"
    work[1] = (doublereal) lwkopt;
#line 388 "dgegv.f"
    lquery = *lwork == -1;
#line 389 "dgegv.f"
    *info = 0;
#line 390 "dgegv.f"
    if (ijobvl <= 0) {
#line 391 "dgegv.f"
	*info = -1;
#line 392 "dgegv.f"
    } else if (ijobvr <= 0) {
#line 393 "dgegv.f"
	*info = -2;
#line 394 "dgegv.f"
    } else if (*n < 0) {
#line 395 "dgegv.f"
	*info = -3;
#line 396 "dgegv.f"
    } else if (*lda < max(1,*n)) {
#line 397 "dgegv.f"
	*info = -5;
#line 398 "dgegv.f"
    } else if (*ldb < max(1,*n)) {
#line 399 "dgegv.f"
	*info = -7;
#line 400 "dgegv.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 401 "dgegv.f"
	*info = -12;
#line 402 "dgegv.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 403 "dgegv.f"
	*info = -14;
#line 404 "dgegv.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 405 "dgegv.f"
	*info = -16;
#line 406 "dgegv.f"
    }

#line 408 "dgegv.f"
    if (*info == 0) {
#line 409 "dgegv.f"
	nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 410 "dgegv.f"
	nb2 = ilaenv_(&c__1, "DORMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 411 "dgegv.f"
	nb3 = ilaenv_(&c__1, "DORGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 412 "dgegv.f"
	i__1 = max(nb1,nb2);
#line 412 "dgegv.f"
	nb = max(i__1,nb3);
/* Computing MAX */
#line 413 "dgegv.f"
	i__1 = *n * 6, i__2 = *n * (nb + 1);
#line 413 "dgegv.f"
	lopt = (*n << 1) + max(i__1,i__2);
#line 414 "dgegv.f"
	work[1] = (doublereal) lopt;
#line 415 "dgegv.f"
    }

#line 417 "dgegv.f"
    if (*info != 0) {
#line 418 "dgegv.f"
	i__1 = -(*info);
#line 418 "dgegv.f"
	xerbla_("DGEGV ", &i__1, (ftnlen)6);
#line 419 "dgegv.f"
	return 0;
#line 420 "dgegv.f"
    } else if (lquery) {
#line 421 "dgegv.f"
	return 0;
#line 422 "dgegv.f"
    }

/*     Quick return if possible */

#line 426 "dgegv.f"
    if (*n == 0) {
#line 426 "dgegv.f"
	return 0;
#line 426 "dgegv.f"
    }

/*     Get machine constants */

#line 431 "dgegv.f"
    eps = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 432 "dgegv.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 433 "dgegv.f"
    safmin += safmin;
#line 434 "dgegv.f"
    safmax = 1. / safmin;
#line 435 "dgegv.f"
    onepls = eps * 4 + 1.;

/*     Scale A */

#line 439 "dgegv.f"
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 440 "dgegv.f"
    anrm1 = anrm;
#line 441 "dgegv.f"
    anrm2 = 1.;
#line 442 "dgegv.f"
    if (anrm < 1.) {
#line 443 "dgegv.f"
	if (safmax * anrm < 1.) {
#line 444 "dgegv.f"
	    anrm1 = safmin;
#line 445 "dgegv.f"
	    anrm2 = safmax * anrm;
#line 446 "dgegv.f"
	}
#line 447 "dgegv.f"
    }

#line 449 "dgegv.f"
    if (anrm > 0.) {
#line 450 "dgegv.f"
	dlascl_("G", &c_n1, &c_n1, &anrm, &c_b27, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 451 "dgegv.f"
	if (iinfo != 0) {
#line 452 "dgegv.f"
	    *info = *n + 10;
#line 453 "dgegv.f"
	    return 0;
#line 454 "dgegv.f"
	}
#line 455 "dgegv.f"
    }

/*     Scale B */

#line 459 "dgegv.f"
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 460 "dgegv.f"
    bnrm1 = bnrm;
#line 461 "dgegv.f"
    bnrm2 = 1.;
#line 462 "dgegv.f"
    if (bnrm < 1.) {
#line 463 "dgegv.f"
	if (safmax * bnrm < 1.) {
#line 464 "dgegv.f"
	    bnrm1 = safmin;
#line 465 "dgegv.f"
	    bnrm2 = safmax * bnrm;
#line 466 "dgegv.f"
	}
#line 467 "dgegv.f"
    }

#line 469 "dgegv.f"
    if (bnrm > 0.) {
#line 470 "dgegv.f"
	dlascl_("G", &c_n1, &c_n1, &bnrm, &c_b27, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 471 "dgegv.f"
	if (iinfo != 0) {
#line 472 "dgegv.f"
	    *info = *n + 10;
#line 473 "dgegv.f"
	    return 0;
#line 474 "dgegv.f"
	}
#line 475 "dgegv.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     Workspace layout:  (8*N words -- "work" requires 6*N words) */
/*        left_permutation, right_permutation, work... */

#line 481 "dgegv.f"
    ileft = 1;
#line 482 "dgegv.f"
    iright = *n + 1;
#line 483 "dgegv.f"
    iwork = iright + *n;
#line 484 "dgegv.f"
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwork], &iinfo, (ftnlen)1);
#line 486 "dgegv.f"
    if (iinfo != 0) {
#line 487 "dgegv.f"
	*info = *n + 1;
#line 488 "dgegv.f"
	goto L120;
#line 489 "dgegv.f"
    }

/*     Reduce B to triangular form, and initialize VL and/or VR */
/*     Workspace layout:  ("work..." must have at least N words) */
/*        left_permutation, right_permutation, tau, work... */

#line 495 "dgegv.f"
    irows = ihi + 1 - ilo;
#line 496 "dgegv.f"
    if (ilv) {
#line 497 "dgegv.f"
	icols = *n + 1 - ilo;
#line 498 "dgegv.f"
    } else {
#line 499 "dgegv.f"
	icols = irows;
#line 500 "dgegv.f"
    }
#line 501 "dgegv.f"
    itau = iwork;
#line 502 "dgegv.f"
    iwork = itau + irows;
#line 503 "dgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 503 "dgegv.f"
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 505 "dgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 505 "dgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 505 "dgegv.f"
	lwkopt = max(i__1,i__2);
#line 505 "dgegv.f"
    }
#line 507 "dgegv.f"
    if (iinfo != 0) {
#line 508 "dgegv.f"
	*info = *n + 2;
#line 509 "dgegv.f"
	goto L120;
#line 510 "dgegv.f"
    }

#line 512 "dgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 512 "dgegv.f"
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 515 "dgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 515 "dgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 515 "dgegv.f"
	lwkopt = max(i__1,i__2);
#line 515 "dgegv.f"
    }
#line 517 "dgegv.f"
    if (iinfo != 0) {
#line 518 "dgegv.f"
	*info = *n + 3;
#line 519 "dgegv.f"
	goto L120;
#line 520 "dgegv.f"
    }

#line 522 "dgegv.f"
    if (ilvl) {
#line 523 "dgegv.f"
	dlaset_("Full", n, n, &c_b38, &c_b27, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 524 "dgegv.f"
	i__1 = irows - 1;
#line 524 "dgegv.f"
	i__2 = irows - 1;
#line 524 "dgegv.f"
	dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[ilo + 
		1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 526 "dgegv.f"
	i__1 = *lwork + 1 - iwork;
#line 526 "dgegv.f"
	dorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwork], &i__1, &iinfo);
#line 529 "dgegv.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 529 "dgegv.f"
	    i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 529 "dgegv.f"
	    lwkopt = max(i__1,i__2);
#line 529 "dgegv.f"
	}
#line 531 "dgegv.f"
	if (iinfo != 0) {
#line 532 "dgegv.f"
	    *info = *n + 4;
#line 533 "dgegv.f"
	    goto L120;
#line 534 "dgegv.f"
	}
#line 535 "dgegv.f"
    }

#line 537 "dgegv.f"
    if (ilvr) {
#line 537 "dgegv.f"
	dlaset_("Full", n, n, &c_b38, &c_b27, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 537 "dgegv.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 542 "dgegv.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 546 "dgegv.f"
	dgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 548 "dgegv.f"
    } else {
#line 549 "dgegv.f"
	dgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 551 "dgegv.f"
    }
#line 552 "dgegv.f"
    if (iinfo != 0) {
#line 553 "dgegv.f"
	*info = *n + 5;
#line 554 "dgegv.f"
	goto L120;
#line 555 "dgegv.f"
    }

/*     Perform QZ algorithm */
/*     Workspace layout:  ("work..." must have at least 1 word) */
/*        left_permutation, right_permutation, work... */

#line 561 "dgegv.f"
    iwork = itau;
#line 562 "dgegv.f"
    if (ilv) {
#line 563 "dgegv.f"
	*(unsigned char *)chtemp = 'S';
#line 564 "dgegv.f"
    } else {
#line 565 "dgegv.f"
	*(unsigned char *)chtemp = 'E';
#line 566 "dgegv.f"
    }
#line 567 "dgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 567 "dgegv.f"
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwork], &i__1, &iinfo, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1);
#line 570 "dgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 570 "dgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 570 "dgegv.f"
	lwkopt = max(i__1,i__2);
#line 570 "dgegv.f"
    }
#line 572 "dgegv.f"
    if (iinfo != 0) {
#line 573 "dgegv.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 574 "dgegv.f"
	    *info = iinfo;
#line 575 "dgegv.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 576 "dgegv.f"
	    *info = iinfo - *n;
#line 577 "dgegv.f"
	} else {
#line 578 "dgegv.f"
	    *info = *n + 6;
#line 579 "dgegv.f"
	}
#line 580 "dgegv.f"
	goto L120;
#line 581 "dgegv.f"
    }

#line 583 "dgegv.f"
    if (ilv) {

/*        Compute Eigenvectors  (DTGEVC requires 6*N words of workspace) */

#line 587 "dgegv.f"
	if (ilvl) {
#line 588 "dgegv.f"
	    if (ilvr) {
#line 589 "dgegv.f"
		*(unsigned char *)chtemp = 'B';
#line 590 "dgegv.f"
	    } else {
#line 591 "dgegv.f"
		*(unsigned char *)chtemp = 'L';
#line 592 "dgegv.f"
	    }
#line 593 "dgegv.f"
	} else {
#line 594 "dgegv.f"
	    *(unsigned char *)chtemp = 'R';
#line 595 "dgegv.f"
	}

#line 597 "dgegv.f"
	dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwork], &iinfo, (ftnlen)1, (ftnlen)1);
#line 599 "dgegv.f"
	if (iinfo != 0) {
#line 600 "dgegv.f"
	    *info = *n + 7;
#line 601 "dgegv.f"
	    goto L120;
#line 602 "dgegv.f"
	}

/*        Undo balancing on VL and VR, rescale */

#line 606 "dgegv.f"
	if (ilvl) {
#line 607 "dgegv.f"
	    dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 609 "dgegv.f"
	    if (iinfo != 0) {
#line 610 "dgegv.f"
		*info = *n + 8;
#line 611 "dgegv.f"
		goto L120;
#line 612 "dgegv.f"
	    }
#line 613 "dgegv.f"
	    i__1 = *n;
#line 613 "dgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 614 "dgegv.f"
		if (alphai[jc] < 0.) {
#line 614 "dgegv.f"
		    goto L50;
#line 614 "dgegv.f"
		}
#line 616 "dgegv.f"
		temp = 0.;
#line 617 "dgegv.f"
		if (alphai[jc] == 0.) {
#line 618 "dgegv.f"
		    i__2 = *n;
#line 618 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 619 "dgegv.f"
			d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1));
#line 619 "dgegv.f"
			temp = max(d__2,d__3);
#line 620 "dgegv.f"
/* L10: */
#line 620 "dgegv.f"
		    }
#line 621 "dgegv.f"
		} else {
#line 622 "dgegv.f"
		    i__2 = *n;
#line 622 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 623 "dgegv.f"
			d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1)) + (d__2 = vl[jr + (jc + 1) * 
				vl_dim1], abs(d__2));
#line 623 "dgegv.f"
			temp = max(d__3,d__4);
#line 625 "dgegv.f"
/* L20: */
#line 625 "dgegv.f"
		    }
#line 626 "dgegv.f"
		}
#line 627 "dgegv.f"
		if (temp < safmin) {
#line 627 "dgegv.f"
		    goto L50;
#line 627 "dgegv.f"
		}
#line 629 "dgegv.f"
		temp = 1. / temp;
#line 630 "dgegv.f"
		if (alphai[jc] == 0.) {
#line 631 "dgegv.f"
		    i__2 = *n;
#line 631 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 632 "dgegv.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 633 "dgegv.f"
/* L30: */
#line 633 "dgegv.f"
		    }
#line 634 "dgegv.f"
		} else {
#line 635 "dgegv.f"
		    i__2 = *n;
#line 635 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 636 "dgegv.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 637 "dgegv.f"
			vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 638 "dgegv.f"
/* L40: */
#line 638 "dgegv.f"
		    }
#line 639 "dgegv.f"
		}
#line 640 "dgegv.f"
L50:
#line 640 "dgegv.f"
		;
#line 640 "dgegv.f"
	    }
#line 641 "dgegv.f"
	}
#line 642 "dgegv.f"
	if (ilvr) {
#line 643 "dgegv.f"
	    dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 645 "dgegv.f"
	    if (iinfo != 0) {
#line 646 "dgegv.f"
		*info = *n + 9;
#line 647 "dgegv.f"
		goto L120;
#line 648 "dgegv.f"
	    }
#line 649 "dgegv.f"
	    i__1 = *n;
#line 649 "dgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 650 "dgegv.f"
		if (alphai[jc] < 0.) {
#line 650 "dgegv.f"
		    goto L100;
#line 650 "dgegv.f"
		}
#line 652 "dgegv.f"
		temp = 0.;
#line 653 "dgegv.f"
		if (alphai[jc] == 0.) {
#line 654 "dgegv.f"
		    i__2 = *n;
#line 654 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 655 "dgegv.f"
			d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1));
#line 655 "dgegv.f"
			temp = max(d__2,d__3);
#line 656 "dgegv.f"
/* L60: */
#line 656 "dgegv.f"
		    }
#line 657 "dgegv.f"
		} else {
#line 658 "dgegv.f"
		    i__2 = *n;
#line 658 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 659 "dgegv.f"
			d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1)) + (d__2 = vr[jr + (jc + 1) * 
				vr_dim1], abs(d__2));
#line 659 "dgegv.f"
			temp = max(d__3,d__4);
#line 661 "dgegv.f"
/* L70: */
#line 661 "dgegv.f"
		    }
#line 662 "dgegv.f"
		}
#line 663 "dgegv.f"
		if (temp < safmin) {
#line 663 "dgegv.f"
		    goto L100;
#line 663 "dgegv.f"
		}
#line 665 "dgegv.f"
		temp = 1. / temp;
#line 666 "dgegv.f"
		if (alphai[jc] == 0.) {
#line 667 "dgegv.f"
		    i__2 = *n;
#line 667 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 668 "dgegv.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 669 "dgegv.f"
/* L80: */
#line 669 "dgegv.f"
		    }
#line 670 "dgegv.f"
		} else {
#line 671 "dgegv.f"
		    i__2 = *n;
#line 671 "dgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 672 "dgegv.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 673 "dgegv.f"
			vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 674 "dgegv.f"
/* L90: */
#line 674 "dgegv.f"
		    }
#line 675 "dgegv.f"
		}
#line 676 "dgegv.f"
L100:
#line 676 "dgegv.f"
		;
#line 676 "dgegv.f"
	    }
#line 677 "dgegv.f"
	}

/*        End of eigenvector calculation */

#line 681 "dgegv.f"
    }

/*     Undo scaling in alpha, beta */

/*     Note: this does not give the alpha and beta for the unscaled */
/*     problem. */

/*     Un-scaling is limited to avoid underflow in alpha and beta */
/*     if they are significant. */

#line 691 "dgegv.f"
    i__1 = *n;
#line 691 "dgegv.f"
    for (jc = 1; jc <= i__1; ++jc) {
#line 692 "dgegv.f"
	absar = (d__1 = alphar[jc], abs(d__1));
#line 693 "dgegv.f"
	absai = (d__1 = alphai[jc], abs(d__1));
#line 694 "dgegv.f"
	absb = (d__1 = beta[jc], abs(d__1));
#line 695 "dgegv.f"
	salfar = anrm * alphar[jc];
#line 696 "dgegv.f"
	salfai = anrm * alphai[jc];
#line 697 "dgegv.f"
	sbeta = bnrm * beta[jc];
#line 698 "dgegv.f"
	ilimit = FALSE_;
#line 699 "dgegv.f"
	scale = 1.;

/*        Check for significant underflow in ALPHAI */

/* Computing MAX */
#line 703 "dgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 703 "dgegv.f"
	if (abs(salfai) < safmin && absai >= max(d__1,d__2)) {
#line 705 "dgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
#line 706 "dgegv.f"
	    d__1 = onepls * safmin, d__2 = anrm2 * absai;
#line 706 "dgegv.f"
	    scale = onepls * safmin / anrm1 / max(d__1,d__2);

#line 709 "dgegv.f"
	} else if (salfai == 0.) {

/*           If insignificant underflow in ALPHAI, then make the */
/*           conjugate eigenvalue real. */

#line 714 "dgegv.f"
	    if (alphai[jc] < 0. && jc > 1) {
#line 715 "dgegv.f"
		alphai[jc - 1] = 0.;
#line 716 "dgegv.f"
	    } else if (alphai[jc] > 0. && jc < *n) {
#line 717 "dgegv.f"
		alphai[jc + 1] = 0.;
#line 718 "dgegv.f"
	    }
#line 719 "dgegv.f"
	}

/*        Check for significant underflow in ALPHAR */

/* Computing MAX */
#line 723 "dgegv.f"
	d__1 = safmin, d__2 = eps * absai, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 723 "dgegv.f"
	if (abs(salfar) < safmin && absar >= max(d__1,d__2)) {
#line 725 "dgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 726 "dgegv.f"
	    d__3 = onepls * safmin, d__4 = anrm2 * absar;
#line 726 "dgegv.f"
	    d__1 = scale, d__2 = onepls * safmin / anrm1 / max(d__3,d__4);
#line 726 "dgegv.f"
	    scale = max(d__1,d__2);
#line 728 "dgegv.f"
	}

/*        Check for significant underflow in BETA */

/* Computing MAX */
#line 732 "dgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absai;
#line 732 "dgegv.f"
	if (abs(sbeta) < safmin && absb >= max(d__1,d__2)) {
#line 734 "dgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 735 "dgegv.f"
	    d__3 = onepls * safmin, d__4 = bnrm2 * absb;
#line 735 "dgegv.f"
	    d__1 = scale, d__2 = onepls * safmin / bnrm1 / max(d__3,d__4);
#line 735 "dgegv.f"
	    scale = max(d__1,d__2);
#line 737 "dgegv.f"
	}

/*        Check for possible overflow when limiting scaling */

#line 741 "dgegv.f"
	if (ilimit) {
/* Computing MAX */
#line 742 "dgegv.f"
	    d__1 = abs(salfar), d__2 = abs(salfai), d__1 = max(d__1,d__2), 
		    d__2 = abs(sbeta);
#line 742 "dgegv.f"
	    temp = scale * safmin * max(d__1,d__2);
#line 744 "dgegv.f"
	    if (temp > 1.) {
#line 744 "dgegv.f"
		scale /= temp;
#line 744 "dgegv.f"
	    }
#line 746 "dgegv.f"
	    if (scale < 1.) {
#line 746 "dgegv.f"
		ilimit = FALSE_;
#line 746 "dgegv.f"
	    }
#line 748 "dgegv.f"
	}

/*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */

#line 752 "dgegv.f"
	if (ilimit) {
#line 753 "dgegv.f"
	    salfar = scale * alphar[jc] * anrm;
#line 754 "dgegv.f"
	    salfai = scale * alphai[jc] * anrm;
#line 755 "dgegv.f"
	    sbeta = scale * beta[jc] * bnrm;
#line 756 "dgegv.f"
	}
#line 757 "dgegv.f"
	alphar[jc] = salfar;
#line 758 "dgegv.f"
	alphai[jc] = salfai;
#line 759 "dgegv.f"
	beta[jc] = sbeta;
#line 760 "dgegv.f"
/* L110: */
#line 760 "dgegv.f"
    }

#line 762 "dgegv.f"
L120:
#line 763 "dgegv.f"
    work[1] = (doublereal) lwkopt;

#line 765 "dgegv.f"
    return 0;

/*     End of DGEGV */

} /* dgegv_ */

