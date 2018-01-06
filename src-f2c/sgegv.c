#line 1 "sgegv.f"
/* sgegv.f -- translated by f2c (version 20100827).
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

#line 1 "sgegv.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b27 = 1.;
static doublereal c_b38 = 0.;

/* > \brief <b> SGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgegv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgegv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgegv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, */
/*                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBVL, JOBVR */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine SGGEV. */
/* > */
/* > SGEGV computes the eigenvalues and, optionally, the left and/or right */
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
/* >          A is REAL array, dimension (LDA, N) */
/* >          On entry, the matrix A. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit A */
/* >          contains the real Schur form of A from the generalized Schur */
/* >          factorization of the pair (A,B) after balancing. */
/* >          If no eigenvectors were computed, then only the diagonal */
/* >          blocks from the Schur form will be correct.  See SGGHRD and */
/* >          SHGEQZ for details. */
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
/* >          B is REAL array, dimension (LDB, N) */
/* >          On entry, the matrix B. */
/* >          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the */
/* >          upper triangular matrix obtained from B in the generalized */
/* >          Schur factorization of the pair (A,B) after balancing. */
/* >          If no eigenvectors were computed, then only those elements of */
/* >          B corresponding to the diagonal blocks from the Schur form of */
/* >          A will be correct.  See SGGHRD and SHGEQZ for details. */
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
/* >          ALPHAR is REAL array, dimension (N) */
/* >          The real parts of each scalar alpha defining an eigenvalue of */
/* >          GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is REAL array, dimension (N) */
/* >          The imaginary parts of each scalar alpha defining an */
/* >          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th */
/* >          eigenvalue is real; if positive, then the j-th and */
/* >          (j+1)-st eigenvalues are a complex conjugate pair, with */
/* >          ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is REAL array, dimension (N) */
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
/* >          VL is REAL array, dimension (LDVL,N) */
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
/* >          VR is REAL array, dimension (LDVR,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,8*N). */
/* >          For good performance, LWORK must generally be larger. */
/* >          To compute the optimal value of LWORK, call ILAENV to get */
/* >          blocksizes (for SGEQRF, SORMQR, and SORGQR.)  Then compute: */
/* >          NB  -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR; */
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
/* >                =N+1: error return from SGGBAL */
/* >                =N+2: error return from SGEQRF */
/* >                =N+3: error return from SORMQR */
/* >                =N+4: error return from SORGQR */
/* >                =N+5: error return from SGGHRD */
/* >                =N+6: error return from SHGEQZ (other than failed */
/* >                                                iteration) */
/* >                =N+7: error return from STGEVC */
/* >                =N+8: error return from SGGBAK (computing VL) */
/* >                =N+9: error return from SGGBAK (computing VR) */
/* >                =N+10: error return from SLASCL (various calls) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realGEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Balancing */
/* >  --------- */
/* > */
/* >  This driver calls SGGBAL to both permute and scale rows and columns */
/* >  of A and B.  The permutations PL and PR are chosen so that PL*A*PR */
/* >  and PL*B*R will be upper triangular except for the diagonal blocks */
/* >  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as */
/* >  possible.  The diagonal scaling matrices DL and DR are chosen so */
/* >  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to */
/* >  one (except for the elements that start out zero.) */
/* > */
/* >  After the eigenvalues and eigenvectors of the balanced matrices */
/* >  have been computed, SGGBAK transforms the eigenvectors back to what */
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
/* >  [*] See SHGEQZ, SGEGS, or read the book "Matrix Computations", */
/* >      by Golub & van Loan, pub. by Johns Hopkins U. Press. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
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
    static doublereal salfai;
    extern /* Subroutine */ int sggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), sggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static doublereal salfar;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int sgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvl, iright;
    static logical ilimit;
    extern /* Subroutine */ int sgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer ijobvr;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), stgevc_(char *, char *, logical 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal onepls;
    static integer lwkmin;
    extern /* Subroutine */ int shgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), sorgqr_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int sormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);


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

#line 360 "sgegv.f"
    /* Parameter adjustments */
#line 360 "sgegv.f"
    a_dim1 = *lda;
#line 360 "sgegv.f"
    a_offset = 1 + a_dim1;
#line 360 "sgegv.f"
    a -= a_offset;
#line 360 "sgegv.f"
    b_dim1 = *ldb;
#line 360 "sgegv.f"
    b_offset = 1 + b_dim1;
#line 360 "sgegv.f"
    b -= b_offset;
#line 360 "sgegv.f"
    --alphar;
#line 360 "sgegv.f"
    --alphai;
#line 360 "sgegv.f"
    --beta;
#line 360 "sgegv.f"
    vl_dim1 = *ldvl;
#line 360 "sgegv.f"
    vl_offset = 1 + vl_dim1;
#line 360 "sgegv.f"
    vl -= vl_offset;
#line 360 "sgegv.f"
    vr_dim1 = *ldvr;
#line 360 "sgegv.f"
    vr_offset = 1 + vr_dim1;
#line 360 "sgegv.f"
    vr -= vr_offset;
#line 360 "sgegv.f"
    --work;
#line 360 "sgegv.f"

#line 360 "sgegv.f"
    /* Function Body */
#line 360 "sgegv.f"
    if (lsame_(jobvl, "N", (ftnlen)1, (ftnlen)1)) {
#line 361 "sgegv.f"
	ijobvl = 1;
#line 362 "sgegv.f"
	ilvl = FALSE_;
#line 363 "sgegv.f"
    } else if (lsame_(jobvl, "V", (ftnlen)1, (ftnlen)1)) {
#line 364 "sgegv.f"
	ijobvl = 2;
#line 365 "sgegv.f"
	ilvl = TRUE_;
#line 366 "sgegv.f"
    } else {
#line 367 "sgegv.f"
	ijobvl = -1;
#line 368 "sgegv.f"
	ilvl = FALSE_;
#line 369 "sgegv.f"
    }

#line 371 "sgegv.f"
    if (lsame_(jobvr, "N", (ftnlen)1, (ftnlen)1)) {
#line 372 "sgegv.f"
	ijobvr = 1;
#line 373 "sgegv.f"
	ilvr = FALSE_;
#line 374 "sgegv.f"
    } else if (lsame_(jobvr, "V", (ftnlen)1, (ftnlen)1)) {
#line 375 "sgegv.f"
	ijobvr = 2;
#line 376 "sgegv.f"
	ilvr = TRUE_;
#line 377 "sgegv.f"
    } else {
#line 378 "sgegv.f"
	ijobvr = -1;
#line 379 "sgegv.f"
	ilvr = FALSE_;
#line 380 "sgegv.f"
    }
#line 381 "sgegv.f"
    ilv = ilvl || ilvr;

/*     Test the input arguments */

/* Computing MAX */
#line 385 "sgegv.f"
    i__1 = *n << 3;
#line 385 "sgegv.f"
    lwkmin = max(i__1,1);
#line 386 "sgegv.f"
    lwkopt = lwkmin;
#line 387 "sgegv.f"
    work[1] = (doublereal) lwkopt;
#line 388 "sgegv.f"
    lquery = *lwork == -1;
#line 389 "sgegv.f"
    *info = 0;
#line 390 "sgegv.f"
    if (ijobvl <= 0) {
#line 391 "sgegv.f"
	*info = -1;
#line 392 "sgegv.f"
    } else if (ijobvr <= 0) {
#line 393 "sgegv.f"
	*info = -2;
#line 394 "sgegv.f"
    } else if (*n < 0) {
#line 395 "sgegv.f"
	*info = -3;
#line 396 "sgegv.f"
    } else if (*lda < max(1,*n)) {
#line 397 "sgegv.f"
	*info = -5;
#line 398 "sgegv.f"
    } else if (*ldb < max(1,*n)) {
#line 399 "sgegv.f"
	*info = -7;
#line 400 "sgegv.f"
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
#line 401 "sgegv.f"
	*info = -12;
#line 402 "sgegv.f"
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
#line 403 "sgegv.f"
	*info = -14;
#line 404 "sgegv.f"
    } else if (*lwork < lwkmin && ! lquery) {
#line 405 "sgegv.f"
	*info = -16;
#line 406 "sgegv.f"
    }

#line 408 "sgegv.f"
    if (*info == 0) {
#line 409 "sgegv.f"
	nb1 = ilaenv_(&c__1, "SGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 410 "sgegv.f"
	nb2 = ilaenv_(&c__1, "SORMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 411 "sgegv.f"
	nb3 = ilaenv_(&c__1, "SORGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
#line 412 "sgegv.f"
	i__1 = max(nb1,nb2);
#line 412 "sgegv.f"
	nb = max(i__1,nb3);
/* Computing MAX */
#line 413 "sgegv.f"
	i__1 = *n * 6, i__2 = *n * (nb + 1);
#line 413 "sgegv.f"
	lopt = (*n << 1) + max(i__1,i__2);
#line 414 "sgegv.f"
	work[1] = (doublereal) lopt;
#line 415 "sgegv.f"
    }

#line 417 "sgegv.f"
    if (*info != 0) {
#line 418 "sgegv.f"
	i__1 = -(*info);
#line 418 "sgegv.f"
	xerbla_("SGEGV ", &i__1, (ftnlen)6);
#line 419 "sgegv.f"
	return 0;
#line 420 "sgegv.f"
    } else if (lquery) {
#line 421 "sgegv.f"
	return 0;
#line 422 "sgegv.f"
    }

/*     Quick return if possible */

#line 426 "sgegv.f"
    if (*n == 0) {
#line 426 "sgegv.f"
	return 0;
#line 426 "sgegv.f"
    }

/*     Get machine constants */

#line 431 "sgegv.f"
    eps = slamch_("E", (ftnlen)1) * slamch_("B", (ftnlen)1);
#line 432 "sgegv.f"
    safmin = slamch_("S", (ftnlen)1);
#line 433 "sgegv.f"
    safmin += safmin;
#line 434 "sgegv.f"
    safmax = 1. / safmin;
#line 435 "sgegv.f"
    onepls = eps * 4 + 1.;

/*     Scale A */

#line 439 "sgegv.f"
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 440 "sgegv.f"
    anrm1 = anrm;
#line 441 "sgegv.f"
    anrm2 = 1.;
#line 442 "sgegv.f"
    if (anrm < 1.) {
#line 443 "sgegv.f"
	if (safmax * anrm < 1.) {
#line 444 "sgegv.f"
	    anrm1 = safmin;
#line 445 "sgegv.f"
	    anrm2 = safmax * anrm;
#line 446 "sgegv.f"
	}
#line 447 "sgegv.f"
    }

#line 449 "sgegv.f"
    if (anrm > 0.) {
#line 450 "sgegv.f"
	slascl_("G", &c_n1, &c_n1, &anrm, &c_b27, n, n, &a[a_offset], lda, &
		iinfo, (ftnlen)1);
#line 451 "sgegv.f"
	if (iinfo != 0) {
#line 452 "sgegv.f"
	    *info = *n + 10;
#line 453 "sgegv.f"
	    return 0;
#line 454 "sgegv.f"
	}
#line 455 "sgegv.f"
    }

/*     Scale B */

#line 459 "sgegv.f"
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1], (ftnlen)1);
#line 460 "sgegv.f"
    bnrm1 = bnrm;
#line 461 "sgegv.f"
    bnrm2 = 1.;
#line 462 "sgegv.f"
    if (bnrm < 1.) {
#line 463 "sgegv.f"
	if (safmax * bnrm < 1.) {
#line 464 "sgegv.f"
	    bnrm1 = safmin;
#line 465 "sgegv.f"
	    bnrm2 = safmax * bnrm;
#line 466 "sgegv.f"
	}
#line 467 "sgegv.f"
    }

#line 469 "sgegv.f"
    if (bnrm > 0.) {
#line 470 "sgegv.f"
	slascl_("G", &c_n1, &c_n1, &bnrm, &c_b27, n, n, &b[b_offset], ldb, &
		iinfo, (ftnlen)1);
#line 471 "sgegv.f"
	if (iinfo != 0) {
#line 472 "sgegv.f"
	    *info = *n + 10;
#line 473 "sgegv.f"
	    return 0;
#line 474 "sgegv.f"
	}
#line 475 "sgegv.f"
    }

/*     Permute the matrix to make it more nearly triangular */
/*     Workspace layout:  (8*N words -- "work" requires 6*N words) */
/*        left_permutation, right_permutation, work... */

#line 481 "sgegv.f"
    ileft = 1;
#line 482 "sgegv.f"
    iright = *n + 1;
#line 483 "sgegv.f"
    iwork = iright + *n;
#line 484 "sgegv.f"
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwork], &iinfo, (ftnlen)1);
#line 486 "sgegv.f"
    if (iinfo != 0) {
#line 487 "sgegv.f"
	*info = *n + 1;
#line 488 "sgegv.f"
	goto L120;
#line 489 "sgegv.f"
    }

/*     Reduce B to triangular form, and initialize VL and/or VR */
/*     Workspace layout:  ("work..." must have at least N words) */
/*        left_permutation, right_permutation, tau, work... */

#line 495 "sgegv.f"
    irows = ihi + 1 - ilo;
#line 496 "sgegv.f"
    if (ilv) {
#line 497 "sgegv.f"
	icols = *n + 1 - ilo;
#line 498 "sgegv.f"
    } else {
#line 499 "sgegv.f"
	icols = irows;
#line 500 "sgegv.f"
    }
#line 501 "sgegv.f"
    itau = iwork;
#line 502 "sgegv.f"
    iwork = itau + irows;
#line 503 "sgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 503 "sgegv.f"
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[
	    iwork], &i__1, &iinfo);
#line 505 "sgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 505 "sgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 505 "sgegv.f"
	lwkopt = max(i__1,i__2);
#line 505 "sgegv.f"
    }
#line 507 "sgegv.f"
    if (iinfo != 0) {
#line 508 "sgegv.f"
	*info = *n + 2;
#line 509 "sgegv.f"
	goto L120;
#line 510 "sgegv.f"
    }

#line 512 "sgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 512 "sgegv.f"
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &
	    work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &
	    iinfo, (ftnlen)1, (ftnlen)1);
#line 515 "sgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 515 "sgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 515 "sgegv.f"
	lwkopt = max(i__1,i__2);
#line 515 "sgegv.f"
    }
#line 517 "sgegv.f"
    if (iinfo != 0) {
#line 518 "sgegv.f"
	*info = *n + 3;
#line 519 "sgegv.f"
	goto L120;
#line 520 "sgegv.f"
    }

#line 522 "sgegv.f"
    if (ilvl) {
#line 523 "sgegv.f"
	slaset_("Full", n, n, &c_b38, &c_b27, &vl[vl_offset], ldvl, (ftnlen)4)
		;
#line 524 "sgegv.f"
	i__1 = irows - 1;
#line 524 "sgegv.f"
	i__2 = irows - 1;
#line 524 "sgegv.f"
	slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[ilo + 
		1 + ilo * vl_dim1], ldvl, (ftnlen)1);
#line 526 "sgegv.f"
	i__1 = *lwork + 1 - iwork;
#line 526 "sgegv.f"
	sorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[
		itau], &work[iwork], &i__1, &iinfo);
#line 529 "sgegv.f"
	if (iinfo >= 0) {
/* Computing MAX */
#line 529 "sgegv.f"
	    i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 529 "sgegv.f"
	    lwkopt = max(i__1,i__2);
#line 529 "sgegv.f"
	}
#line 531 "sgegv.f"
	if (iinfo != 0) {
#line 532 "sgegv.f"
	    *info = *n + 4;
#line 533 "sgegv.f"
	    goto L120;
#line 534 "sgegv.f"
	}
#line 535 "sgegv.f"
    }

#line 537 "sgegv.f"
    if (ilvr) {
#line 537 "sgegv.f"
	slaset_("Full", n, n, &c_b38, &c_b27, &vr[vr_offset], ldvr, (ftnlen)4)
		;
#line 537 "sgegv.f"
    }

/*     Reduce to generalized Hessenberg form */

#line 542 "sgegv.f"
    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

#line 546 "sgegv.f"
	sgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo, (
		ftnlen)1, (ftnlen)1);
#line 548 "sgegv.f"
    } else {
#line 549 "sgegv.f"
	sgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, 
		&b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[
		vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 551 "sgegv.f"
    }
#line 552 "sgegv.f"
    if (iinfo != 0) {
#line 553 "sgegv.f"
	*info = *n + 5;
#line 554 "sgegv.f"
	goto L120;
#line 555 "sgegv.f"
    }

/*     Perform QZ algorithm */
/*     Workspace layout:  ("work..." must have at least 1 word) */
/*        left_permutation, right_permutation, work... */

#line 561 "sgegv.f"
    iwork = itau;
#line 562 "sgegv.f"
    if (ilv) {
#line 563 "sgegv.f"
	*(unsigned char *)chtemp = 'S';
#line 564 "sgegv.f"
    } else {
#line 565 "sgegv.f"
	*(unsigned char *)chtemp = 'E';
#line 566 "sgegv.f"
    }
#line 567 "sgegv.f"
    i__1 = *lwork + 1 - iwork;
#line 567 "sgegv.f"
    shgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwork], &i__1, &iinfo, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1);
#line 570 "sgegv.f"
    if (iinfo >= 0) {
/* Computing MAX */
#line 570 "sgegv.f"
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
#line 570 "sgegv.f"
	lwkopt = max(i__1,i__2);
#line 570 "sgegv.f"
    }
#line 572 "sgegv.f"
    if (iinfo != 0) {
#line 573 "sgegv.f"
	if (iinfo > 0 && iinfo <= *n) {
#line 574 "sgegv.f"
	    *info = iinfo;
#line 575 "sgegv.f"
	} else if (iinfo > *n && iinfo <= *n << 1) {
#line 576 "sgegv.f"
	    *info = iinfo - *n;
#line 577 "sgegv.f"
	} else {
#line 578 "sgegv.f"
	    *info = *n + 6;
#line 579 "sgegv.f"
	}
#line 580 "sgegv.f"
	goto L120;
#line 581 "sgegv.f"
    }

#line 583 "sgegv.f"
    if (ilv) {

/*        Compute Eigenvectors  (STGEVC requires 6*N words of workspace) */

#line 587 "sgegv.f"
	if (ilvl) {
#line 588 "sgegv.f"
	    if (ilvr) {
#line 589 "sgegv.f"
		*(unsigned char *)chtemp = 'B';
#line 590 "sgegv.f"
	    } else {
#line 591 "sgegv.f"
		*(unsigned char *)chtemp = 'L';
#line 592 "sgegv.f"
	    }
#line 593 "sgegv.f"
	} else {
#line 594 "sgegv.f"
	    *(unsigned char *)chtemp = 'R';
#line 595 "sgegv.f"
	}

#line 597 "sgegv.f"
	stgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwork], &iinfo, (ftnlen)1, (ftnlen)1);
#line 599 "sgegv.f"
	if (iinfo != 0) {
#line 600 "sgegv.f"
	    *info = *n + 7;
#line 601 "sgegv.f"
	    goto L120;
#line 602 "sgegv.f"
	}

/*        Undo balancing on VL and VR, rescale */

#line 606 "sgegv.f"
	if (ilvl) {
#line 607 "sgegv.f"
	    sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &iinfo, (ftnlen)1, (ftnlen)1);
#line 609 "sgegv.f"
	    if (iinfo != 0) {
#line 610 "sgegv.f"
		*info = *n + 8;
#line 611 "sgegv.f"
		goto L120;
#line 612 "sgegv.f"
	    }
#line 613 "sgegv.f"
	    i__1 = *n;
#line 613 "sgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 614 "sgegv.f"
		if (alphai[jc] < 0.) {
#line 614 "sgegv.f"
		    goto L50;
#line 614 "sgegv.f"
		}
#line 616 "sgegv.f"
		temp = 0.;
#line 617 "sgegv.f"
		if (alphai[jc] == 0.) {
#line 618 "sgegv.f"
		    i__2 = *n;
#line 618 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 619 "sgegv.f"
			d__2 = temp, d__3 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1));
#line 619 "sgegv.f"
			temp = max(d__2,d__3);
#line 620 "sgegv.f"
/* L10: */
#line 620 "sgegv.f"
		    }
#line 621 "sgegv.f"
		} else {
#line 622 "sgegv.f"
		    i__2 = *n;
#line 622 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 623 "sgegv.f"
			d__3 = temp, d__4 = (d__1 = vl[jr + jc * vl_dim1], 
				abs(d__1)) + (d__2 = vl[jr + (jc + 1) * 
				vl_dim1], abs(d__2));
#line 623 "sgegv.f"
			temp = max(d__3,d__4);
#line 625 "sgegv.f"
/* L20: */
#line 625 "sgegv.f"
		    }
#line 626 "sgegv.f"
		}
#line 627 "sgegv.f"
		if (temp < safmin) {
#line 627 "sgegv.f"
		    goto L50;
#line 627 "sgegv.f"
		}
#line 629 "sgegv.f"
		temp = 1. / temp;
#line 630 "sgegv.f"
		if (alphai[jc] == 0.) {
#line 631 "sgegv.f"
		    i__2 = *n;
#line 631 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 632 "sgegv.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 633 "sgegv.f"
/* L30: */
#line 633 "sgegv.f"
		    }
#line 634 "sgegv.f"
		} else {
#line 635 "sgegv.f"
		    i__2 = *n;
#line 635 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 636 "sgegv.f"
			vl[jr + jc * vl_dim1] *= temp;
#line 637 "sgegv.f"
			vl[jr + (jc + 1) * vl_dim1] *= temp;
#line 638 "sgegv.f"
/* L40: */
#line 638 "sgegv.f"
		    }
#line 639 "sgegv.f"
		}
#line 640 "sgegv.f"
L50:
#line 640 "sgegv.f"
		;
#line 640 "sgegv.f"
	    }
#line 641 "sgegv.f"
	}
#line 642 "sgegv.f"
	if (ilvr) {
#line 643 "sgegv.f"
	    sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &iinfo, (ftnlen)1, (ftnlen)1);
#line 645 "sgegv.f"
	    if (iinfo != 0) {
#line 646 "sgegv.f"
		*info = *n + 9;
#line 647 "sgegv.f"
		goto L120;
#line 648 "sgegv.f"
	    }
#line 649 "sgegv.f"
	    i__1 = *n;
#line 649 "sgegv.f"
	    for (jc = 1; jc <= i__1; ++jc) {
#line 650 "sgegv.f"
		if (alphai[jc] < 0.) {
#line 650 "sgegv.f"
		    goto L100;
#line 650 "sgegv.f"
		}
#line 652 "sgegv.f"
		temp = 0.;
#line 653 "sgegv.f"
		if (alphai[jc] == 0.) {
#line 654 "sgegv.f"
		    i__2 = *n;
#line 654 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 655 "sgegv.f"
			d__2 = temp, d__3 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1));
#line 655 "sgegv.f"
			temp = max(d__2,d__3);
#line 656 "sgegv.f"
/* L60: */
#line 656 "sgegv.f"
		    }
#line 657 "sgegv.f"
		} else {
#line 658 "sgegv.f"
		    i__2 = *n;
#line 658 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
#line 659 "sgegv.f"
			d__3 = temp, d__4 = (d__1 = vr[jr + jc * vr_dim1], 
				abs(d__1)) + (d__2 = vr[jr + (jc + 1) * 
				vr_dim1], abs(d__2));
#line 659 "sgegv.f"
			temp = max(d__3,d__4);
#line 661 "sgegv.f"
/* L70: */
#line 661 "sgegv.f"
		    }
#line 662 "sgegv.f"
		}
#line 663 "sgegv.f"
		if (temp < safmin) {
#line 663 "sgegv.f"
		    goto L100;
#line 663 "sgegv.f"
		}
#line 665 "sgegv.f"
		temp = 1. / temp;
#line 666 "sgegv.f"
		if (alphai[jc] == 0.) {
#line 667 "sgegv.f"
		    i__2 = *n;
#line 667 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 668 "sgegv.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 669 "sgegv.f"
/* L80: */
#line 669 "sgegv.f"
		    }
#line 670 "sgegv.f"
		} else {
#line 671 "sgegv.f"
		    i__2 = *n;
#line 671 "sgegv.f"
		    for (jr = 1; jr <= i__2; ++jr) {
#line 672 "sgegv.f"
			vr[jr + jc * vr_dim1] *= temp;
#line 673 "sgegv.f"
			vr[jr + (jc + 1) * vr_dim1] *= temp;
#line 674 "sgegv.f"
/* L90: */
#line 674 "sgegv.f"
		    }
#line 675 "sgegv.f"
		}
#line 676 "sgegv.f"
L100:
#line 676 "sgegv.f"
		;
#line 676 "sgegv.f"
	    }
#line 677 "sgegv.f"
	}

/*        End of eigenvector calculation */

#line 681 "sgegv.f"
    }

/*     Undo scaling in alpha, beta */

/*     Note: this does not give the alpha and beta for the unscaled */
/*     problem. */

/*     Un-scaling is limited to avoid underflow in alpha and beta */
/*     if they are significant. */

#line 691 "sgegv.f"
    i__1 = *n;
#line 691 "sgegv.f"
    for (jc = 1; jc <= i__1; ++jc) {
#line 692 "sgegv.f"
	absar = (d__1 = alphar[jc], abs(d__1));
#line 693 "sgegv.f"
	absai = (d__1 = alphai[jc], abs(d__1));
#line 694 "sgegv.f"
	absb = (d__1 = beta[jc], abs(d__1));
#line 695 "sgegv.f"
	salfar = anrm * alphar[jc];
#line 696 "sgegv.f"
	salfai = anrm * alphai[jc];
#line 697 "sgegv.f"
	sbeta = bnrm * beta[jc];
#line 698 "sgegv.f"
	ilimit = FALSE_;
#line 699 "sgegv.f"
	scale = 1.;

/*        Check for significant underflow in ALPHAI */

/* Computing MAX */
#line 703 "sgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 703 "sgegv.f"
	if (abs(salfai) < safmin && absai >= max(d__1,d__2)) {
#line 705 "sgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
#line 706 "sgegv.f"
	    d__1 = onepls * safmin, d__2 = anrm2 * absai;
#line 706 "sgegv.f"
	    scale = onepls * safmin / anrm1 / max(d__1,d__2);

#line 709 "sgegv.f"
	} else if (salfai == 0.) {

/*           If insignificant underflow in ALPHAI, then make the */
/*           conjugate eigenvalue real. */

#line 714 "sgegv.f"
	    if (alphai[jc] < 0. && jc > 1) {
#line 715 "sgegv.f"
		alphai[jc - 1] = 0.;
#line 716 "sgegv.f"
	    } else if (alphai[jc] > 0. && jc < *n) {
#line 717 "sgegv.f"
		alphai[jc + 1] = 0.;
#line 718 "sgegv.f"
	    }
#line 719 "sgegv.f"
	}

/*        Check for significant underflow in ALPHAR */

/* Computing MAX */
#line 723 "sgegv.f"
	d__1 = safmin, d__2 = eps * absai, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
#line 723 "sgegv.f"
	if (abs(salfar) < safmin && absar >= max(d__1,d__2)) {
#line 725 "sgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 726 "sgegv.f"
	    d__3 = onepls * safmin, d__4 = anrm2 * absar;
#line 726 "sgegv.f"
	    d__1 = scale, d__2 = onepls * safmin / anrm1 / max(d__3,d__4);
#line 726 "sgegv.f"
	    scale = max(d__1,d__2);
#line 728 "sgegv.f"
	}

/*        Check for significant underflow in BETA */

/* Computing MAX */
#line 732 "sgegv.f"
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absai;
#line 732 "sgegv.f"
	if (abs(sbeta) < safmin && absb >= max(d__1,d__2)) {
#line 734 "sgegv.f"
	    ilimit = TRUE_;
/* Computing MAX */
/* Computing MAX */
#line 735 "sgegv.f"
	    d__3 = onepls * safmin, d__4 = bnrm2 * absb;
#line 735 "sgegv.f"
	    d__1 = scale, d__2 = onepls * safmin / bnrm1 / max(d__3,d__4);
#line 735 "sgegv.f"
	    scale = max(d__1,d__2);
#line 737 "sgegv.f"
	}

/*        Check for possible overflow when limiting scaling */

#line 741 "sgegv.f"
	if (ilimit) {
/* Computing MAX */
#line 742 "sgegv.f"
	    d__1 = abs(salfar), d__2 = abs(salfai), d__1 = max(d__1,d__2), 
		    d__2 = abs(sbeta);
#line 742 "sgegv.f"
	    temp = scale * safmin * max(d__1,d__2);
#line 744 "sgegv.f"
	    if (temp > 1.) {
#line 744 "sgegv.f"
		scale /= temp;
#line 744 "sgegv.f"
	    }
#line 746 "sgegv.f"
	    if (scale < 1.) {
#line 746 "sgegv.f"
		ilimit = FALSE_;
#line 746 "sgegv.f"
	    }
#line 748 "sgegv.f"
	}

/*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */

#line 752 "sgegv.f"
	if (ilimit) {
#line 753 "sgegv.f"
	    salfar = scale * alphar[jc] * anrm;
#line 754 "sgegv.f"
	    salfai = scale * alphai[jc] * anrm;
#line 755 "sgegv.f"
	    sbeta = scale * beta[jc] * bnrm;
#line 756 "sgegv.f"
	}
#line 757 "sgegv.f"
	alphar[jc] = salfar;
#line 758 "sgegv.f"
	alphai[jc] = salfai;
#line 759 "sgegv.f"
	beta[jc] = sbeta;
#line 760 "sgegv.f"
/* L110: */
#line 760 "sgegv.f"
    }

#line 762 "sgegv.f"
L120:
#line 763 "sgegv.f"
    work[1] = (doublereal) lwkopt;

#line 765 "sgegv.f"
    return 0;

/*     End of SGEGV */

} /* sgegv_ */

