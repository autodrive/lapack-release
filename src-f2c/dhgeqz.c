#line 1 "dhgeqz.f"
/* dhgeqz.f -- translated by f2c (version 20100827).
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

#line 1 "dhgeqz.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static integer c__1 = 1;
static integer c__3 = 3;

/* > \brief \b DHGEQZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DHGEQZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhgeqz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhgeqz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhgeqz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, */
/*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, */
/*                          LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ, JOB */
/*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ), */
/*      $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DHGEQZ computes the eigenvalues of a real matrix pair (H,T), */
/* > where H is an upper Hessenberg matrix and T is upper triangular, */
/* > using the double-shift QZ method. */
/* > Matrix pairs of this type are produced by the reduction to */
/* > generalized upper Hessenberg form of a real matrix pair (A,B): */
/* > */
/* >    A = Q1*H*Z1**T,  B = Q1*T*Z1**T, */
/* > */
/* > as computed by DGGHRD. */
/* > */
/* > If JOB='S', then the Hessenberg-triangular pair (H,T) is */
/* > also reduced to generalized Schur form, */
/* > */
/* >    H = Q*S*Z**T,  T = Q*P*Z**T, */
/* > */
/* > where Q and Z are orthogonal matrices, P is an upper triangular */
/* > matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2 */
/* > diagonal blocks. */
/* > */
/* > The 1-by-1 blocks correspond to real eigenvalues of the matrix pair */
/* > (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of */
/* > eigenvalues. */
/* > */
/* > Additionally, the 2-by-2 upper triangular diagonal blocks of P */
/* > corresponding to 2-by-2 blocks of S are reduced to positive diagonal */
/* > form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0, */
/* > P(j,j) > 0, and P(j+1,j+1) > 0. */
/* > */
/* > Optionally, the orthogonal matrix Q from the generalized Schur */
/* > factorization may be postmultiplied into an input matrix Q1, and the */
/* > orthogonal matrix Z may be postmultiplied into an input matrix Z1. */
/* > If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced */
/* > the matrix pair (A,B) to generalized upper Hessenberg form, then the */
/* > output matrices Q1*Q and Z1*Z are the orthogonal factors from the */
/* > generalized Schur factorization of (A,B): */
/* > */
/* >    A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T. */
/* > */
/* > To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently, */
/* > of (A,B)) are computed as a pair of values (alpha,beta), where alpha is */
/* > complex and beta real. */
/* > If beta is nonzero, lambda = alpha / beta is an eigenvalue of the */
/* > generalized nonsymmetric eigenvalue problem (GNEP) */
/* >    A*x = lambda*B*x */
/* > and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the */
/* > alternate form of the GNEP */
/* >    mu*A*y = B*y. */
/* > Real eigenvalues can be read directly from the generalized Schur */
/* > form: */
/* >   alpha = S(i,i), beta = P(i,i). */
/* > */
/* > Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix */
/* >      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973), */
/* >      pp. 241--256. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          = 'E': Compute eigenvalues only; */
/* >          = 'S': Compute eigenvalues and the Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'N': Left Schur vectors (Q) are not computed; */
/* >          = 'I': Q is initialized to the unit matrix and the matrix Q */
/* >                 of left Schur vectors of (H,T) is returned; */
/* >          = 'V': Q must contain an orthogonal matrix Q1 on entry and */
/* >                 the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N': Right Schur vectors (Z) are not computed; */
/* >          = 'I': Z is initialized to the unit matrix and the matrix Z */
/* >                 of right Schur vectors of (H,T) is returned; */
/* >          = 'V': Z must contain an orthogonal matrix Z1 on entry and */
/* >                 the product Z1*Z is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices H, T, Q, and Z.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >          ILO and IHI mark the rows and columns of H which are in */
/* >          Hessenberg form.  It is assumed that A is already upper */
/* >          triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/* >          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION array, dimension (LDH, N) */
/* >          On entry, the N-by-N upper Hessenberg matrix H. */
/* >          On exit, if JOB = 'S', H contains the upper quasi-triangular */
/* >          matrix S from the generalized Schur factorization. */
/* >          If JOB = 'E', the diagonal blocks of H match those of S, but */
/* >          the rest of H is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* >          LDH is INTEGER */
/* >          The leading dimension of the array H.  LDH >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT, N) */
/* >          On entry, the N-by-N upper triangular matrix T. */
/* >          On exit, if JOB = 'S', T contains the upper triangular */
/* >          matrix P from the generalized Schur factorization; */
/* >          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S */
/* >          are reduced to positive diagonal form, i.e., if H(j+1,j) is */
/* >          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and */
/* >          T(j+1,j+1) > 0. */
/* >          If JOB = 'E', the diagonal blocks of T match those of P, but */
/* >          the rest of T is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* >          ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* >          The real parts of each scalar alpha defining an eigenvalue */
/* >          of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* >          The imaginary parts of each scalar alpha defining an */
/* >          eigenvalue of GNEP. */
/* >          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if */
/* >          positive, then the j-th and (j+1)-st eigenvalues are a */
/* >          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
/* >          The scalars beta that define the eigenvalues of GNEP. */
/* >          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* >          beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* >          pair (A,B), in one of the forms lambda = alpha/beta or */
/* >          mu = beta/alpha.  Since either lambda or mu may overflow, */
/* >          they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* >          On entry, if COMPZ = 'V', the orthogonal matrix Q1 used in */
/* >          the reduction of (A,B) to generalized Hessenberg form. */
/* >          On exit, if COMPZ = 'I', the orthogonal matrix of left Schur */
/* >          vectors of (H,T), and if COMPZ = 'V', the orthogonal matrix */
/* >          of left Schur vectors of (A,B). */
/* >          Not referenced if COMPZ = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= 1. */
/* >          If COMPQ='V' or 'I', then LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in */
/* >          the reduction of (A,B) to generalized Hessenberg form. */
/* >          On exit, if COMPZ = 'I', the orthogonal matrix of */
/* >          right Schur vectors of (H,T), and if COMPZ = 'V', the */
/* >          orthogonal matrix of right Schur vectors of (A,B). */
/* >          Not referenced if COMPZ = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1. */
/* >          If COMPZ='V' or 'I', then LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK.  LWORK >= max(1,N). */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          = 1,...,N: the QZ iteration did not converge.  (H,T) is not */
/* >                     in Schur form, but ALPHAR(i), ALPHAI(i), and */
/* >                     BETA(i), i=INFO+1,...,N should be correct. */
/* >          = N+1,...,2*N: the shift calculation failed.  (H,T) is not */
/* >                     in Schur form, but ALPHAR(i), ALPHAI(i), and */
/* >                     BETA(i), i=INFO-N+1,...,N should be correct. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2013 */

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Iteration counters: */
/* > */
/* >  JITER  -- counts iterations. */
/* >  IITER  -- counts iterations run since ILAST was last */
/* >            changed.  This is therefore reset only when a 1-by-1 or */
/* >            2-by-2 block deflates off the bottom. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*t, integer *ldt, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	doublereal *work, integer *lwork, integer *info, ftnlen job_len, 
	ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, q_dim1, q_offset, t_dim1, t_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublereal s, v[3], s1, s2, t1, u1, u2, a11, a12, a21, a22, b11, 
	    b22, c12, c21;
    static integer jc;
    static doublereal an, bn, cl, cq, cr;
    static integer in;
    static doublereal u12, w11, w12, w21;
    static integer jr;
    static doublereal cz, w22, sl, wi, sr, vs, wr, b1a, b2a, a1i, a2i, b1i, 
	    b2i, a1r, a2r, b1r, b2r, wr2, ad11, ad12, ad21, ad22, c11i, c22i;
    static integer jch;
    static doublereal c11r, c22r;
    static logical ilq;
    static doublereal u12l, tau, sqi;
    static logical ilz;
    static doublereal ulp, sqr, szi, szr, ad11l, ad12l, ad21l, ad22l, ad32l, 
	    wabs, atol, btol, temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlag2_(
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal temp2, s1inv, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iiter, ilast, jiter;
    static doublereal anorm, bnorm;
    static integer maxit;
    static doublereal tempi, tempr;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static logical ilazr2;
    static doublereal ascale, bscale;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal eshift;
    static logical ilschr;
    static integer icompq, ilastm, ischur;
    static logical ilazro;
    static integer icompz, ifirst, ifrstm, istart;
    static logical ilpivt, lquery;


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*    $                     SAFETY = 1.0E+0 ) */
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

/*     Decode JOB, COMPQ, COMPZ */

#line 366 "dhgeqz.f"
    /* Parameter adjustments */
#line 366 "dhgeqz.f"
    h_dim1 = *ldh;
#line 366 "dhgeqz.f"
    h_offset = 1 + h_dim1;
#line 366 "dhgeqz.f"
    h__ -= h_offset;
#line 366 "dhgeqz.f"
    t_dim1 = *ldt;
#line 366 "dhgeqz.f"
    t_offset = 1 + t_dim1;
#line 366 "dhgeqz.f"
    t -= t_offset;
#line 366 "dhgeqz.f"
    --alphar;
#line 366 "dhgeqz.f"
    --alphai;
#line 366 "dhgeqz.f"
    --beta;
#line 366 "dhgeqz.f"
    q_dim1 = *ldq;
#line 366 "dhgeqz.f"
    q_offset = 1 + q_dim1;
#line 366 "dhgeqz.f"
    q -= q_offset;
#line 366 "dhgeqz.f"
    z_dim1 = *ldz;
#line 366 "dhgeqz.f"
    z_offset = 1 + z_dim1;
#line 366 "dhgeqz.f"
    z__ -= z_offset;
#line 366 "dhgeqz.f"
    --work;
#line 366 "dhgeqz.f"

#line 366 "dhgeqz.f"
    /* Function Body */
#line 366 "dhgeqz.f"
    if (lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
#line 367 "dhgeqz.f"
	ilschr = FALSE_;
#line 368 "dhgeqz.f"
	ischur = 1;
#line 369 "dhgeqz.f"
    } else if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 370 "dhgeqz.f"
	ilschr = TRUE_;
#line 371 "dhgeqz.f"
	ischur = 2;
#line 372 "dhgeqz.f"
    } else {
#line 373 "dhgeqz.f"
	ischur = 0;
#line 374 "dhgeqz.f"
    }

#line 376 "dhgeqz.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 377 "dhgeqz.f"
	ilq = FALSE_;
#line 378 "dhgeqz.f"
	icompq = 1;
#line 379 "dhgeqz.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 380 "dhgeqz.f"
	ilq = TRUE_;
#line 381 "dhgeqz.f"
	icompq = 2;
#line 382 "dhgeqz.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 383 "dhgeqz.f"
	ilq = TRUE_;
#line 384 "dhgeqz.f"
	icompq = 3;
#line 385 "dhgeqz.f"
    } else {
#line 386 "dhgeqz.f"
	icompq = 0;
#line 387 "dhgeqz.f"
    }

#line 389 "dhgeqz.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 390 "dhgeqz.f"
	ilz = FALSE_;
#line 391 "dhgeqz.f"
	icompz = 1;
#line 392 "dhgeqz.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 393 "dhgeqz.f"
	ilz = TRUE_;
#line 394 "dhgeqz.f"
	icompz = 2;
#line 395 "dhgeqz.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 396 "dhgeqz.f"
	ilz = TRUE_;
#line 397 "dhgeqz.f"
	icompz = 3;
#line 398 "dhgeqz.f"
    } else {
#line 399 "dhgeqz.f"
	icompz = 0;
#line 400 "dhgeqz.f"
    }

/*     Check Argument Values */

#line 404 "dhgeqz.f"
    *info = 0;
#line 405 "dhgeqz.f"
    work[1] = (doublereal) max(1,*n);
#line 406 "dhgeqz.f"
    lquery = *lwork == -1;
#line 407 "dhgeqz.f"
    if (ischur == 0) {
#line 408 "dhgeqz.f"
	*info = -1;
#line 409 "dhgeqz.f"
    } else if (icompq == 0) {
#line 410 "dhgeqz.f"
	*info = -2;
#line 411 "dhgeqz.f"
    } else if (icompz == 0) {
#line 412 "dhgeqz.f"
	*info = -3;
#line 413 "dhgeqz.f"
    } else if (*n < 0) {
#line 414 "dhgeqz.f"
	*info = -4;
#line 415 "dhgeqz.f"
    } else if (*ilo < 1) {
#line 416 "dhgeqz.f"
	*info = -5;
#line 417 "dhgeqz.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 418 "dhgeqz.f"
	*info = -6;
#line 419 "dhgeqz.f"
    } else if (*ldh < *n) {
#line 420 "dhgeqz.f"
	*info = -8;
#line 421 "dhgeqz.f"
    } else if (*ldt < *n) {
#line 422 "dhgeqz.f"
	*info = -10;
#line 423 "dhgeqz.f"
    } else if (*ldq < 1 || ilq && *ldq < *n) {
#line 424 "dhgeqz.f"
	*info = -15;
#line 425 "dhgeqz.f"
    } else if (*ldz < 1 || ilz && *ldz < *n) {
#line 426 "dhgeqz.f"
	*info = -17;
#line 427 "dhgeqz.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 428 "dhgeqz.f"
	*info = -19;
#line 429 "dhgeqz.f"
    }
#line 430 "dhgeqz.f"
    if (*info != 0) {
#line 431 "dhgeqz.f"
	i__1 = -(*info);
#line 431 "dhgeqz.f"
	xerbla_("DHGEQZ", &i__1, (ftnlen)6);
#line 432 "dhgeqz.f"
	return 0;
#line 433 "dhgeqz.f"
    } else if (lquery) {
#line 434 "dhgeqz.f"
	return 0;
#line 435 "dhgeqz.f"
    }

/*     Quick return if possible */

#line 439 "dhgeqz.f"
    if (*n <= 0) {
#line 440 "dhgeqz.f"
	work[1] = 1.;
#line 441 "dhgeqz.f"
	return 0;
#line 442 "dhgeqz.f"
    }

/*     Initialize Q and Z */

#line 446 "dhgeqz.f"
    if (icompq == 3) {
#line 446 "dhgeqz.f"
	dlaset_("Full", n, n, &c_b12, &c_b13, &q[q_offset], ldq, (ftnlen)4);
#line 446 "dhgeqz.f"
    }
#line 448 "dhgeqz.f"
    if (icompz == 3) {
#line 448 "dhgeqz.f"
	dlaset_("Full", n, n, &c_b12, &c_b13, &z__[z_offset], ldz, (ftnlen)4);
#line 448 "dhgeqz.f"
    }

/*     Machine Constants */

#line 453 "dhgeqz.f"
    in = *ihi + 1 - *ilo;
#line 454 "dhgeqz.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 455 "dhgeqz.f"
    safmax = 1. / safmin;
#line 456 "dhgeqz.f"
    ulp = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 457 "dhgeqz.f"
    anorm = dlanhs_("F", &in, &h__[*ilo + *ilo * h_dim1], ldh, &work[1], (
	    ftnlen)1);
#line 458 "dhgeqz.f"
    bnorm = dlanhs_("F", &in, &t[*ilo + *ilo * t_dim1], ldt, &work[1], (
	    ftnlen)1);
/* Computing MAX */
#line 459 "dhgeqz.f"
    d__1 = safmin, d__2 = ulp * anorm;
#line 459 "dhgeqz.f"
    atol = max(d__1,d__2);
/* Computing MAX */
#line 460 "dhgeqz.f"
    d__1 = safmin, d__2 = ulp * bnorm;
#line 460 "dhgeqz.f"
    btol = max(d__1,d__2);
#line 461 "dhgeqz.f"
    ascale = 1. / max(safmin,anorm);
#line 462 "dhgeqz.f"
    bscale = 1. / max(safmin,bnorm);

/*     Set Eigenvalues IHI+1:N */

#line 466 "dhgeqz.f"
    i__1 = *n;
#line 466 "dhgeqz.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 467 "dhgeqz.f"
	if (t[j + j * t_dim1] < 0.) {
#line 468 "dhgeqz.f"
	    if (ilschr) {
#line 469 "dhgeqz.f"
		i__2 = j;
#line 469 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 470 "dhgeqz.f"
		    h__[jr + j * h_dim1] = -h__[jr + j * h_dim1];
#line 471 "dhgeqz.f"
		    t[jr + j * t_dim1] = -t[jr + j * t_dim1];
#line 472 "dhgeqz.f"
/* L10: */
#line 472 "dhgeqz.f"
		}
#line 473 "dhgeqz.f"
	    } else {
#line 474 "dhgeqz.f"
		h__[j + j * h_dim1] = -h__[j + j * h_dim1];
#line 475 "dhgeqz.f"
		t[j + j * t_dim1] = -t[j + j * t_dim1];
#line 476 "dhgeqz.f"
	    }
#line 477 "dhgeqz.f"
	    if (ilz) {
#line 478 "dhgeqz.f"
		i__2 = *n;
#line 478 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 479 "dhgeqz.f"
		    z__[jr + j * z_dim1] = -z__[jr + j * z_dim1];
#line 480 "dhgeqz.f"
/* L20: */
#line 480 "dhgeqz.f"
		}
#line 481 "dhgeqz.f"
	    }
#line 482 "dhgeqz.f"
	}
#line 483 "dhgeqz.f"
	alphar[j] = h__[j + j * h_dim1];
#line 484 "dhgeqz.f"
	alphai[j] = 0.;
#line 485 "dhgeqz.f"
	beta[j] = t[j + j * t_dim1];
#line 486 "dhgeqz.f"
/* L30: */
#line 486 "dhgeqz.f"
    }

/*     If IHI < ILO, skip QZ steps */

#line 490 "dhgeqz.f"
    if (*ihi < *ilo) {
#line 490 "dhgeqz.f"
	goto L380;
#line 490 "dhgeqz.f"
    }

/*     MAIN QZ ITERATION LOOP */

/*     Initialize dynamic indices */

/*     Eigenvalues ILAST+1:N have been found. */
/*        Column operations modify rows IFRSTM:whatever. */
/*        Row operations modify columns whatever:ILASTM. */

/*     If only eigenvalues are being computed, then */
/*        IFRSTM is the row of the last splitting row above row ILAST; */
/*        this is always at least ILO. */
/*     IITER counts iterations since the last eigenvalue was found, */
/*        to tell when to use an extraordinary shift. */
/*     MAXIT is the maximum number of QZ sweeps allowed. */

#line 508 "dhgeqz.f"
    ilast = *ihi;
#line 509 "dhgeqz.f"
    if (ilschr) {
#line 510 "dhgeqz.f"
	ifrstm = 1;
#line 511 "dhgeqz.f"
	ilastm = *n;
#line 512 "dhgeqz.f"
    } else {
#line 513 "dhgeqz.f"
	ifrstm = *ilo;
#line 514 "dhgeqz.f"
	ilastm = *ihi;
#line 515 "dhgeqz.f"
    }
#line 516 "dhgeqz.f"
    iiter = 0;
#line 517 "dhgeqz.f"
    eshift = 0.;
#line 518 "dhgeqz.f"
    maxit = (*ihi - *ilo + 1) * 30;

#line 520 "dhgeqz.f"
    i__1 = maxit;
#line 520 "dhgeqz.f"
    for (jiter = 1; jiter <= i__1; ++jiter) {

/*        Split the matrix if possible. */

/*        Two tests: */
/*           1: H(j,j-1)=0  or  j=ILO */
/*           2: T(j,j)=0 */

#line 528 "dhgeqz.f"
	if (ilast == *ilo) {

/*           Special case: j=ILAST */

#line 532 "dhgeqz.f"
	    goto L80;
#line 533 "dhgeqz.f"
	} else {
#line 534 "dhgeqz.f"
	    if ((d__1 = h__[ilast + (ilast - 1) * h_dim1], abs(d__1)) <= atol)
		     {
#line 535 "dhgeqz.f"
		h__[ilast + (ilast - 1) * h_dim1] = 0.;
#line 536 "dhgeqz.f"
		goto L80;
#line 537 "dhgeqz.f"
	    }
#line 538 "dhgeqz.f"
	}

#line 540 "dhgeqz.f"
	if ((d__1 = t[ilast + ilast * t_dim1], abs(d__1)) <= btol) {
#line 541 "dhgeqz.f"
	    t[ilast + ilast * t_dim1] = 0.;
#line 542 "dhgeqz.f"
	    goto L70;
#line 543 "dhgeqz.f"
	}

/*        General case: j<ILAST */

#line 547 "dhgeqz.f"
	i__2 = *ilo;
#line 547 "dhgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {

/*           Test 1: for H(j,j-1)=0 or j=ILO */

#line 551 "dhgeqz.f"
	    if (j == *ilo) {
#line 552 "dhgeqz.f"
		ilazro = TRUE_;
#line 553 "dhgeqz.f"
	    } else {
#line 554 "dhgeqz.f"
		if ((d__1 = h__[j + (j - 1) * h_dim1], abs(d__1)) <= atol) {
#line 555 "dhgeqz.f"
		    h__[j + (j - 1) * h_dim1] = 0.;
#line 556 "dhgeqz.f"
		    ilazro = TRUE_;
#line 557 "dhgeqz.f"
		} else {
#line 558 "dhgeqz.f"
		    ilazro = FALSE_;
#line 559 "dhgeqz.f"
		}
#line 560 "dhgeqz.f"
	    }

/*           Test 2: for T(j,j)=0 */

#line 564 "dhgeqz.f"
	    if ((d__1 = t[j + j * t_dim1], abs(d__1)) < btol) {
#line 565 "dhgeqz.f"
		t[j + j * t_dim1] = 0.;

/*              Test 1a: Check for 2 consecutive small subdiagonals in A */

#line 569 "dhgeqz.f"
		ilazr2 = FALSE_;
#line 570 "dhgeqz.f"
		if (! ilazro) {
#line 571 "dhgeqz.f"
		    temp = (d__1 = h__[j + (j - 1) * h_dim1], abs(d__1));
#line 572 "dhgeqz.f"
		    temp2 = (d__1 = h__[j + j * h_dim1], abs(d__1));
#line 573 "dhgeqz.f"
		    tempr = max(temp,temp2);
#line 574 "dhgeqz.f"
		    if (tempr < 1. && tempr != 0.) {
#line 575 "dhgeqz.f"
			temp /= tempr;
#line 576 "dhgeqz.f"
			temp2 /= tempr;
#line 577 "dhgeqz.f"
		    }
#line 578 "dhgeqz.f"
		    if (temp * (ascale * (d__1 = h__[j + 1 + j * h_dim1], abs(
			    d__1))) <= temp2 * (ascale * atol)) {
#line 578 "dhgeqz.f"
			ilazr2 = TRUE_;
#line 578 "dhgeqz.f"
		    }
#line 580 "dhgeqz.f"
		}

/*              If both tests pass (1 & 2), i.e., the leading diagonal */
/*              element of B in the block is zero, split a 1x1 block off */
/*              at the top. (I.e., at the J-th row/column) The leading */
/*              diagonal element of the remainder can also be zero, so */
/*              this may have to be done repeatedly. */

#line 588 "dhgeqz.f"
		if (ilazro || ilazr2) {
#line 589 "dhgeqz.f"
		    i__3 = ilast - 1;
#line 589 "dhgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 590 "dhgeqz.f"
			temp = h__[jch + jch * h_dim1];
#line 591 "dhgeqz.f"
			dlartg_(&temp, &h__[jch + 1 + jch * h_dim1], &c__, &s,
				 &h__[jch + jch * h_dim1]);
#line 593 "dhgeqz.f"
			h__[jch + 1 + jch * h_dim1] = 0.;
#line 594 "dhgeqz.f"
			i__4 = ilastm - jch;
#line 594 "dhgeqz.f"
			drot_(&i__4, &h__[jch + (jch + 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch + 1) * h_dim1], ldh, &c__, 
				&s);
#line 596 "dhgeqz.f"
			i__4 = ilastm - jch;
#line 596 "dhgeqz.f"
			drot_(&i__4, &t[jch + (jch + 1) * t_dim1], ldt, &t[
				jch + 1 + (jch + 1) * t_dim1], ldt, &c__, &s);
#line 598 "dhgeqz.f"
			if (ilq) {
#line 598 "dhgeqz.f"
			    drot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &s);
#line 598 "dhgeqz.f"
			}
#line 601 "dhgeqz.f"
			if (ilazr2) {
#line 601 "dhgeqz.f"
			    h__[jch + (jch - 1) * h_dim1] *= c__;
#line 601 "dhgeqz.f"
			}
#line 603 "dhgeqz.f"
			ilazr2 = FALSE_;
#line 604 "dhgeqz.f"
			if ((d__1 = t[jch + 1 + (jch + 1) * t_dim1], abs(d__1)
				) >= btol) {
#line 605 "dhgeqz.f"
			    if (jch + 1 >= ilast) {
#line 606 "dhgeqz.f"
				goto L80;
#line 607 "dhgeqz.f"
			    } else {
#line 608 "dhgeqz.f"
				ifirst = jch + 1;
#line 609 "dhgeqz.f"
				goto L110;
#line 610 "dhgeqz.f"
			    }
#line 611 "dhgeqz.f"
			}
#line 612 "dhgeqz.f"
			t[jch + 1 + (jch + 1) * t_dim1] = 0.;
#line 613 "dhgeqz.f"
/* L40: */
#line 613 "dhgeqz.f"
		    }
#line 614 "dhgeqz.f"
		    goto L70;
#line 615 "dhgeqz.f"
		} else {

/*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST) */
/*                 Then process as in the case T(ILAST,ILAST)=0 */

#line 620 "dhgeqz.f"
		    i__3 = ilast - 1;
#line 620 "dhgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 621 "dhgeqz.f"
			temp = t[jch + (jch + 1) * t_dim1];
#line 622 "dhgeqz.f"
			dlartg_(&temp, &t[jch + 1 + (jch + 1) * t_dim1], &c__,
				 &s, &t[jch + (jch + 1) * t_dim1]);
#line 624 "dhgeqz.f"
			t[jch + 1 + (jch + 1) * t_dim1] = 0.;
#line 625 "dhgeqz.f"
			if (jch < ilastm - 1) {
#line 625 "dhgeqz.f"
			    i__4 = ilastm - jch - 1;
#line 625 "dhgeqz.f"
			    drot_(&i__4, &t[jch + (jch + 2) * t_dim1], ldt, &
				    t[jch + 1 + (jch + 2) * t_dim1], ldt, &
				    c__, &s);
#line 625 "dhgeqz.f"
			}
#line 628 "dhgeqz.f"
			i__4 = ilastm - jch + 2;
#line 628 "dhgeqz.f"
			drot_(&i__4, &h__[jch + (jch - 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch - 1) * h_dim1], ldh, &c__, 
				&s);
#line 630 "dhgeqz.f"
			if (ilq) {
#line 630 "dhgeqz.f"
			    drot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &s);
#line 630 "dhgeqz.f"
			}
#line 633 "dhgeqz.f"
			temp = h__[jch + 1 + jch * h_dim1];
#line 634 "dhgeqz.f"
			dlartg_(&temp, &h__[jch + 1 + (jch - 1) * h_dim1], &
				c__, &s, &h__[jch + 1 + jch * h_dim1]);
#line 636 "dhgeqz.f"
			h__[jch + 1 + (jch - 1) * h_dim1] = 0.;
#line 637 "dhgeqz.f"
			i__4 = jch + 1 - ifrstm;
#line 637 "dhgeqz.f"
			drot_(&i__4, &h__[ifrstm + jch * h_dim1], &c__1, &h__[
				ifrstm + (jch - 1) * h_dim1], &c__1, &c__, &s)
				;
#line 639 "dhgeqz.f"
			i__4 = jch - ifrstm;
#line 639 "dhgeqz.f"
			drot_(&i__4, &t[ifrstm + jch * t_dim1], &c__1, &t[
				ifrstm + (jch - 1) * t_dim1], &c__1, &c__, &s)
				;
#line 641 "dhgeqz.f"
			if (ilz) {
#line 641 "dhgeqz.f"
			    drot_(n, &z__[jch * z_dim1 + 1], &c__1, &z__[(jch 
				    - 1) * z_dim1 + 1], &c__1, &c__, &s);
#line 641 "dhgeqz.f"
			}
#line 644 "dhgeqz.f"
/* L50: */
#line 644 "dhgeqz.f"
		    }
#line 645 "dhgeqz.f"
		    goto L70;
#line 646 "dhgeqz.f"
		}
#line 647 "dhgeqz.f"
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

#line 651 "dhgeqz.f"
		ifirst = j;
#line 652 "dhgeqz.f"
		goto L110;
#line 653 "dhgeqz.f"
	    }

/*           Neither test passed -- try next J */

#line 657 "dhgeqz.f"
/* L60: */
#line 657 "dhgeqz.f"
	}

/*        (Drop-through is "impossible") */

#line 661 "dhgeqz.f"
	*info = *n + 1;
#line 662 "dhgeqz.f"
	goto L420;

/*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a */
/*        1x1 block. */

#line 667 "dhgeqz.f"
L70:
#line 668 "dhgeqz.f"
	temp = h__[ilast + ilast * h_dim1];
#line 669 "dhgeqz.f"
	dlartg_(&temp, &h__[ilast + (ilast - 1) * h_dim1], &c__, &s, &h__[
		ilast + ilast * h_dim1]);
#line 671 "dhgeqz.f"
	h__[ilast + (ilast - 1) * h_dim1] = 0.;
#line 672 "dhgeqz.f"
	i__2 = ilast - ifrstm;
#line 672 "dhgeqz.f"
	drot_(&i__2, &h__[ifrstm + ilast * h_dim1], &c__1, &h__[ifrstm + (
		ilast - 1) * h_dim1], &c__1, &c__, &s);
#line 674 "dhgeqz.f"
	i__2 = ilast - ifrstm;
#line 674 "dhgeqz.f"
	drot_(&i__2, &t[ifrstm + ilast * t_dim1], &c__1, &t[ifrstm + (ilast - 
		1) * t_dim1], &c__1, &c__, &s);
#line 676 "dhgeqz.f"
	if (ilz) {
#line 676 "dhgeqz.f"
	    drot_(n, &z__[ilast * z_dim1 + 1], &c__1, &z__[(ilast - 1) * 
		    z_dim1 + 1], &c__1, &c__, &s);
#line 676 "dhgeqz.f"
	}

/*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI, */
/*                              and BETA */

#line 682 "dhgeqz.f"
L80:
#line 683 "dhgeqz.f"
	if (t[ilast + ilast * t_dim1] < 0.) {
#line 684 "dhgeqz.f"
	    if (ilschr) {
#line 685 "dhgeqz.f"
		i__2 = ilast;
#line 685 "dhgeqz.f"
		for (j = ifrstm; j <= i__2; ++j) {
#line 686 "dhgeqz.f"
		    h__[j + ilast * h_dim1] = -h__[j + ilast * h_dim1];
#line 687 "dhgeqz.f"
		    t[j + ilast * t_dim1] = -t[j + ilast * t_dim1];
#line 688 "dhgeqz.f"
/* L90: */
#line 688 "dhgeqz.f"
		}
#line 689 "dhgeqz.f"
	    } else {
#line 690 "dhgeqz.f"
		h__[ilast + ilast * h_dim1] = -h__[ilast + ilast * h_dim1];
#line 691 "dhgeqz.f"
		t[ilast + ilast * t_dim1] = -t[ilast + ilast * t_dim1];
#line 692 "dhgeqz.f"
	    }
#line 693 "dhgeqz.f"
	    if (ilz) {
#line 694 "dhgeqz.f"
		i__2 = *n;
#line 694 "dhgeqz.f"
		for (j = 1; j <= i__2; ++j) {
#line 695 "dhgeqz.f"
		    z__[j + ilast * z_dim1] = -z__[j + ilast * z_dim1];
#line 696 "dhgeqz.f"
/* L100: */
#line 696 "dhgeqz.f"
		}
#line 697 "dhgeqz.f"
	    }
#line 698 "dhgeqz.f"
	}
#line 699 "dhgeqz.f"
	alphar[ilast] = h__[ilast + ilast * h_dim1];
#line 700 "dhgeqz.f"
	alphai[ilast] = 0.;
#line 701 "dhgeqz.f"
	beta[ilast] = t[ilast + ilast * t_dim1];

/*        Go to next block -- exit if finished. */

#line 705 "dhgeqz.f"
	--ilast;
#line 706 "dhgeqz.f"
	if (ilast < *ilo) {
#line 706 "dhgeqz.f"
	    goto L380;
#line 706 "dhgeqz.f"
	}

/*        Reset counters */

#line 711 "dhgeqz.f"
	iiter = 0;
#line 712 "dhgeqz.f"
	eshift = 0.;
#line 713 "dhgeqz.f"
	if (! ilschr) {
#line 714 "dhgeqz.f"
	    ilastm = ilast;
#line 715 "dhgeqz.f"
	    if (ifrstm > ilast) {
#line 715 "dhgeqz.f"
		ifrstm = *ilo;
#line 715 "dhgeqz.f"
	    }
#line 717 "dhgeqz.f"
	}
#line 718 "dhgeqz.f"
	goto L350;

/*        QZ step */

/*        This iteration only involves rows/columns IFIRST:ILAST. We */
/*        assume IFIRST < ILAST, and that the diagonal of B is non-zero. */

#line 725 "dhgeqz.f"
L110:
#line 726 "dhgeqz.f"
	++iiter;
#line 727 "dhgeqz.f"
	if (! ilschr) {
#line 728 "dhgeqz.f"
	    ifrstm = ifirst;
#line 729 "dhgeqz.f"
	}

/*        Compute single shifts. */

/*        At this point, IFIRST < ILAST, and the diagonal elements of */
/*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in */
/*        magnitude) */

#line 737 "dhgeqz.f"
	if (iiter / 10 * 10 == iiter) {

/*           Exceptional shift.  Chosen for no particularly good reason. */
/*           (Single shift only.) */

#line 742 "dhgeqz.f"
	    if ((doublereal) maxit * safmin * (d__1 = h__[ilast + (ilast - 1) 
		    * h_dim1], abs(d__1)) < (d__2 = t[ilast - 1 + (ilast - 1) 
		    * t_dim1], abs(d__2))) {
#line 744 "dhgeqz.f"
		eshift = h__[ilast + (ilast - 1) * h_dim1] / t[ilast - 1 + (
			ilast - 1) * t_dim1];
#line 746 "dhgeqz.f"
	    } else {
#line 747 "dhgeqz.f"
		eshift += 1. / (safmin * (doublereal) maxit);
#line 748 "dhgeqz.f"
	    }
#line 749 "dhgeqz.f"
	    s1 = 1.;
#line 750 "dhgeqz.f"
	    wr = eshift;

#line 752 "dhgeqz.f"
	} else {

/*           Shifts based on the generalized eigenvalues of the */
/*           bottom-right 2x2 block of A and B. The first eigenvalue */
/*           returned by DLAG2 is the Wilkinson shift (AEP p.512), */

#line 758 "dhgeqz.f"
	    d__1 = safmin * 100.;
#line 758 "dhgeqz.f"
	    dlag2_(&h__[ilast - 1 + (ilast - 1) * h_dim1], ldh, &t[ilast - 1 
		    + (ilast - 1) * t_dim1], ldt, &d__1, &s1, &s2, &wr, &wr2, 
		    &wi);

#line 762 "dhgeqz.f"
	    if ((d__1 = wr / s1 * t[ilast + ilast * t_dim1] - h__[ilast + 
		    ilast * h_dim1], abs(d__1)) > (d__2 = wr2 / s2 * t[ilast 
		    + ilast * t_dim1] - h__[ilast + ilast * h_dim1], abs(d__2)
		    )) {
#line 765 "dhgeqz.f"
		temp = wr;
#line 766 "dhgeqz.f"
		wr = wr2;
#line 767 "dhgeqz.f"
		wr2 = temp;
#line 768 "dhgeqz.f"
		temp = s1;
#line 769 "dhgeqz.f"
		s1 = s2;
#line 770 "dhgeqz.f"
		s2 = temp;
#line 771 "dhgeqz.f"
	    }
/* Computing MAX */
/* Computing MAX */
#line 772 "dhgeqz.f"
	    d__3 = 1., d__4 = abs(wr), d__3 = max(d__3,d__4), d__4 = abs(wi);
#line 772 "dhgeqz.f"
	    d__1 = s1, d__2 = safmin * max(d__3,d__4);
#line 772 "dhgeqz.f"
	    temp = max(d__1,d__2);
#line 773 "dhgeqz.f"
	    if (wi != 0.) {
#line 773 "dhgeqz.f"
		goto L200;
#line 773 "dhgeqz.f"
	    }
#line 775 "dhgeqz.f"
	}

/*        Fiddle with shift to avoid overflow */

#line 779 "dhgeqz.f"
	temp = min(ascale,1.) * (safmax * .5);
#line 780 "dhgeqz.f"
	if (s1 > temp) {
#line 781 "dhgeqz.f"
	    scale = temp / s1;
#line 782 "dhgeqz.f"
	} else {
#line 783 "dhgeqz.f"
	    scale = 1.;
#line 784 "dhgeqz.f"
	}

#line 786 "dhgeqz.f"
	temp = min(bscale,1.) * (safmax * .5);
#line 787 "dhgeqz.f"
	if (abs(wr) > temp) {
/* Computing MIN */
#line 787 "dhgeqz.f"
	    d__1 = scale, d__2 = temp / abs(wr);
#line 787 "dhgeqz.f"
	    scale = min(d__1,d__2);
#line 787 "dhgeqz.f"
	}
#line 789 "dhgeqz.f"
	s1 = scale * s1;
#line 790 "dhgeqz.f"
	wr = scale * wr;

/*        Now check for two consecutive small subdiagonals. */

#line 794 "dhgeqz.f"
	i__2 = ifirst + 1;
#line 794 "dhgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {
#line 795 "dhgeqz.f"
	    istart = j;
#line 796 "dhgeqz.f"
	    temp = (d__1 = s1 * h__[j + (j - 1) * h_dim1], abs(d__1));
#line 797 "dhgeqz.f"
	    temp2 = (d__1 = s1 * h__[j + j * h_dim1] - wr * t[j + j * t_dim1],
		     abs(d__1));
#line 798 "dhgeqz.f"
	    tempr = max(temp,temp2);
#line 799 "dhgeqz.f"
	    if (tempr < 1. && tempr != 0.) {
#line 800 "dhgeqz.f"
		temp /= tempr;
#line 801 "dhgeqz.f"
		temp2 /= tempr;
#line 802 "dhgeqz.f"
	    }
#line 803 "dhgeqz.f"
	    if ((d__1 = ascale * h__[j + 1 + j * h_dim1] * temp, abs(d__1)) <=
		     ascale * atol * temp2) {
#line 803 "dhgeqz.f"
		goto L130;
#line 803 "dhgeqz.f"
	    }
#line 805 "dhgeqz.f"
/* L120: */
#line 805 "dhgeqz.f"
	}

#line 807 "dhgeqz.f"
	istart = ifirst;
#line 808 "dhgeqz.f"
L130:

/*        Do an implicit single-shift QZ sweep. */

/*        Initial Q */

#line 814 "dhgeqz.f"
	temp = s1 * h__[istart + istart * h_dim1] - wr * t[istart + istart * 
		t_dim1];
#line 815 "dhgeqz.f"
	temp2 = s1 * h__[istart + 1 + istart * h_dim1];
#line 816 "dhgeqz.f"
	dlartg_(&temp, &temp2, &c__, &s, &tempr);

/*        Sweep */

#line 820 "dhgeqz.f"
	i__2 = ilast - 1;
#line 820 "dhgeqz.f"
	for (j = istart; j <= i__2; ++j) {
#line 821 "dhgeqz.f"
	    if (j > istart) {
#line 822 "dhgeqz.f"
		temp = h__[j + (j - 1) * h_dim1];
#line 823 "dhgeqz.f"
		dlartg_(&temp, &h__[j + 1 + (j - 1) * h_dim1], &c__, &s, &h__[
			j + (j - 1) * h_dim1]);
#line 824 "dhgeqz.f"
		h__[j + 1 + (j - 1) * h_dim1] = 0.;
#line 825 "dhgeqz.f"
	    }

#line 827 "dhgeqz.f"
	    i__3 = ilastm;
#line 827 "dhgeqz.f"
	    for (jc = j; jc <= i__3; ++jc) {
#line 828 "dhgeqz.f"
		temp = c__ * h__[j + jc * h_dim1] + s * h__[j + 1 + jc * 
			h_dim1];
#line 829 "dhgeqz.f"
		h__[j + 1 + jc * h_dim1] = -s * h__[j + jc * h_dim1] + c__ * 
			h__[j + 1 + jc * h_dim1];
#line 830 "dhgeqz.f"
		h__[j + jc * h_dim1] = temp;
#line 831 "dhgeqz.f"
		temp2 = c__ * t[j + jc * t_dim1] + s * t[j + 1 + jc * t_dim1];
#line 832 "dhgeqz.f"
		t[j + 1 + jc * t_dim1] = -s * t[j + jc * t_dim1] + c__ * t[j 
			+ 1 + jc * t_dim1];
#line 833 "dhgeqz.f"
		t[j + jc * t_dim1] = temp2;
#line 834 "dhgeqz.f"
/* L140: */
#line 834 "dhgeqz.f"
	    }
#line 835 "dhgeqz.f"
	    if (ilq) {
#line 836 "dhgeqz.f"
		i__3 = *n;
#line 836 "dhgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 837 "dhgeqz.f"
		    temp = c__ * q[jr + j * q_dim1] + s * q[jr + (j + 1) * 
			    q_dim1];
#line 838 "dhgeqz.f"
		    q[jr + (j + 1) * q_dim1] = -s * q[jr + j * q_dim1] + c__ *
			     q[jr + (j + 1) * q_dim1];
#line 839 "dhgeqz.f"
		    q[jr + j * q_dim1] = temp;
#line 840 "dhgeqz.f"
/* L150: */
#line 840 "dhgeqz.f"
		}
#line 841 "dhgeqz.f"
	    }

#line 843 "dhgeqz.f"
	    temp = t[j + 1 + (j + 1) * t_dim1];
#line 844 "dhgeqz.f"
	    dlartg_(&temp, &t[j + 1 + j * t_dim1], &c__, &s, &t[j + 1 + (j + 
		    1) * t_dim1]);
#line 845 "dhgeqz.f"
	    t[j + 1 + j * t_dim1] = 0.;

/* Computing MIN */
#line 847 "dhgeqz.f"
	    i__4 = j + 2;
#line 847 "dhgeqz.f"
	    i__3 = min(i__4,ilast);
#line 847 "dhgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 848 "dhgeqz.f"
		temp = c__ * h__[jr + (j + 1) * h_dim1] + s * h__[jr + j * 
			h_dim1];
#line 849 "dhgeqz.f"
		h__[jr + j * h_dim1] = -s * h__[jr + (j + 1) * h_dim1] + c__ *
			 h__[jr + j * h_dim1];
#line 850 "dhgeqz.f"
		h__[jr + (j + 1) * h_dim1] = temp;
#line 851 "dhgeqz.f"
/* L160: */
#line 851 "dhgeqz.f"
	    }
#line 852 "dhgeqz.f"
	    i__3 = j;
#line 852 "dhgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 853 "dhgeqz.f"
		temp = c__ * t[jr + (j + 1) * t_dim1] + s * t[jr + j * t_dim1]
			;
#line 854 "dhgeqz.f"
		t[jr + j * t_dim1] = -s * t[jr + (j + 1) * t_dim1] + c__ * t[
			jr + j * t_dim1];
#line 855 "dhgeqz.f"
		t[jr + (j + 1) * t_dim1] = temp;
#line 856 "dhgeqz.f"
/* L170: */
#line 856 "dhgeqz.f"
	    }
#line 857 "dhgeqz.f"
	    if (ilz) {
#line 858 "dhgeqz.f"
		i__3 = *n;
#line 858 "dhgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 859 "dhgeqz.f"
		    temp = c__ * z__[jr + (j + 1) * z_dim1] + s * z__[jr + j *
			     z_dim1];
#line 860 "dhgeqz.f"
		    z__[jr + j * z_dim1] = -s * z__[jr + (j + 1) * z_dim1] + 
			    c__ * z__[jr + j * z_dim1];
#line 861 "dhgeqz.f"
		    z__[jr + (j + 1) * z_dim1] = temp;
#line 862 "dhgeqz.f"
/* L180: */
#line 862 "dhgeqz.f"
		}
#line 863 "dhgeqz.f"
	    }
#line 864 "dhgeqz.f"
/* L190: */
#line 864 "dhgeqz.f"
	}

#line 866 "dhgeqz.f"
	goto L350;

/*        Use Francis double-shift */

/*        Note: the Francis double-shift should work with real shifts, */
/*              but only if the block is at least 3x3. */
/*              This code may break if this point is reached with */
/*              a 2x2 block with real eigenvalues. */

#line 875 "dhgeqz.f"
L200:
#line 876 "dhgeqz.f"
	if (ifirst + 1 == ilast) {

/*           Special case -- 2x2 block with complex eigenvectors */

/*           Step 1: Standardize, that is, rotate so that */

/*                       ( B11  0  ) */
/*                   B = (         )  with B11 non-negative. */
/*                       (  0  B22 ) */

#line 886 "dhgeqz.f"
	    dlasv2_(&t[ilast - 1 + (ilast - 1) * t_dim1], &t[ilast - 1 + 
		    ilast * t_dim1], &t[ilast + ilast * t_dim1], &b22, &b11, &
		    sr, &cr, &sl, &cl);

#line 889 "dhgeqz.f"
	    if (b11 < 0.) {
#line 890 "dhgeqz.f"
		cr = -cr;
#line 891 "dhgeqz.f"
		sr = -sr;
#line 892 "dhgeqz.f"
		b11 = -b11;
#line 893 "dhgeqz.f"
		b22 = -b22;
#line 894 "dhgeqz.f"
	    }

#line 896 "dhgeqz.f"
	    i__2 = ilastm + 1 - ifirst;
#line 896 "dhgeqz.f"
	    drot_(&i__2, &h__[ilast - 1 + (ilast - 1) * h_dim1], ldh, &h__[
		    ilast + (ilast - 1) * h_dim1], ldh, &cl, &sl);
#line 898 "dhgeqz.f"
	    i__2 = ilast + 1 - ifrstm;
#line 898 "dhgeqz.f"
	    drot_(&i__2, &h__[ifrstm + (ilast - 1) * h_dim1], &c__1, &h__[
		    ifrstm + ilast * h_dim1], &c__1, &cr, &sr);

#line 901 "dhgeqz.f"
	    if (ilast < ilastm) {
#line 901 "dhgeqz.f"
		i__2 = ilastm - ilast;
#line 901 "dhgeqz.f"
		drot_(&i__2, &t[ilast - 1 + (ilast + 1) * t_dim1], ldt, &t[
			ilast + (ilast + 1) * t_dim1], ldt, &cl, &sl);
#line 901 "dhgeqz.f"
	    }
#line 904 "dhgeqz.f"
	    if (ifrstm < ilast - 1) {
#line 904 "dhgeqz.f"
		i__2 = ifirst - ifrstm;
#line 904 "dhgeqz.f"
		drot_(&i__2, &t[ifrstm + (ilast - 1) * t_dim1], &c__1, &t[
			ifrstm + ilast * t_dim1], &c__1, &cr, &sr);
#line 904 "dhgeqz.f"
	    }

#line 908 "dhgeqz.f"
	    if (ilq) {
#line 908 "dhgeqz.f"
		drot_(n, &q[(ilast - 1) * q_dim1 + 1], &c__1, &q[ilast * 
			q_dim1 + 1], &c__1, &cl, &sl);
#line 908 "dhgeqz.f"
	    }
#line 911 "dhgeqz.f"
	    if (ilz) {
#line 911 "dhgeqz.f"
		drot_(n, &z__[(ilast - 1) * z_dim1 + 1], &c__1, &z__[ilast * 
			z_dim1 + 1], &c__1, &cr, &sr);
#line 911 "dhgeqz.f"
	    }

#line 915 "dhgeqz.f"
	    t[ilast - 1 + (ilast - 1) * t_dim1] = b11;
#line 916 "dhgeqz.f"
	    t[ilast - 1 + ilast * t_dim1] = 0.;
#line 917 "dhgeqz.f"
	    t[ilast + (ilast - 1) * t_dim1] = 0.;
#line 918 "dhgeqz.f"
	    t[ilast + ilast * t_dim1] = b22;

/*           If B22 is negative, negate column ILAST */

#line 922 "dhgeqz.f"
	    if (b22 < 0.) {
#line 923 "dhgeqz.f"
		i__2 = ilast;
#line 923 "dhgeqz.f"
		for (j = ifrstm; j <= i__2; ++j) {
#line 924 "dhgeqz.f"
		    h__[j + ilast * h_dim1] = -h__[j + ilast * h_dim1];
#line 925 "dhgeqz.f"
		    t[j + ilast * t_dim1] = -t[j + ilast * t_dim1];
#line 926 "dhgeqz.f"
/* L210: */
#line 926 "dhgeqz.f"
		}

#line 928 "dhgeqz.f"
		if (ilz) {
#line 929 "dhgeqz.f"
		    i__2 = *n;
#line 929 "dhgeqz.f"
		    for (j = 1; j <= i__2; ++j) {
#line 930 "dhgeqz.f"
			z__[j + ilast * z_dim1] = -z__[j + ilast * z_dim1];
#line 931 "dhgeqz.f"
/* L220: */
#line 931 "dhgeqz.f"
		    }
#line 932 "dhgeqz.f"
		}
#line 933 "dhgeqz.f"
		b22 = -b22;
#line 934 "dhgeqz.f"
	    }

/*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.) */

/*           Recompute shift */

#line 940 "dhgeqz.f"
	    d__1 = safmin * 100.;
#line 940 "dhgeqz.f"
	    dlag2_(&h__[ilast - 1 + (ilast - 1) * h_dim1], ldh, &t[ilast - 1 
		    + (ilast - 1) * t_dim1], ldt, &d__1, &s1, &temp, &wr, &
		    temp2, &wi);

/*           If standardization has perturbed the shift onto real line, */
/*           do another (real single-shift) QR step. */

#line 947 "dhgeqz.f"
	    if (wi == 0.) {
#line 947 "dhgeqz.f"
		goto L350;
#line 947 "dhgeqz.f"
	    }
#line 949 "dhgeqz.f"
	    s1inv = 1. / s1;

/*           Do EISPACK (QZVAL) computation of alpha and beta */

#line 953 "dhgeqz.f"
	    a11 = h__[ilast - 1 + (ilast - 1) * h_dim1];
#line 954 "dhgeqz.f"
	    a21 = h__[ilast + (ilast - 1) * h_dim1];
#line 955 "dhgeqz.f"
	    a12 = h__[ilast - 1 + ilast * h_dim1];
#line 956 "dhgeqz.f"
	    a22 = h__[ilast + ilast * h_dim1];

/*           Compute complex Givens rotation on right */
/*           (Assume some element of C = (sA - wB) > unfl ) */
/*                            __ */
/*           (sA - wB) ( CZ   -SZ ) */
/*                     ( SZ    CZ ) */

#line 964 "dhgeqz.f"
	    c11r = s1 * a11 - wr * b11;
#line 965 "dhgeqz.f"
	    c11i = -wi * b11;
#line 966 "dhgeqz.f"
	    c12 = s1 * a12;
#line 967 "dhgeqz.f"
	    c21 = s1 * a21;
#line 968 "dhgeqz.f"
	    c22r = s1 * a22 - wr * b22;
#line 969 "dhgeqz.f"
	    c22i = -wi * b22;

#line 971 "dhgeqz.f"
	    if (abs(c11r) + abs(c11i) + abs(c12) > abs(c21) + abs(c22r) + abs(
		    c22i)) {
#line 973 "dhgeqz.f"
		t1 = dlapy3_(&c12, &c11r, &c11i);
#line 974 "dhgeqz.f"
		cz = c12 / t1;
#line 975 "dhgeqz.f"
		szr = -c11r / t1;
#line 976 "dhgeqz.f"
		szi = -c11i / t1;
#line 977 "dhgeqz.f"
	    } else {
#line 978 "dhgeqz.f"
		cz = dlapy2_(&c22r, &c22i);
#line 979 "dhgeqz.f"
		if (cz <= safmin) {
#line 980 "dhgeqz.f"
		    cz = 0.;
#line 981 "dhgeqz.f"
		    szr = 1.;
#line 982 "dhgeqz.f"
		    szi = 0.;
#line 983 "dhgeqz.f"
		} else {
#line 984 "dhgeqz.f"
		    tempr = c22r / cz;
#line 985 "dhgeqz.f"
		    tempi = c22i / cz;
#line 986 "dhgeqz.f"
		    t1 = dlapy2_(&cz, &c21);
#line 987 "dhgeqz.f"
		    cz /= t1;
#line 988 "dhgeqz.f"
		    szr = -c21 * tempr / t1;
#line 989 "dhgeqz.f"
		    szi = c21 * tempi / t1;
#line 990 "dhgeqz.f"
		}
#line 991 "dhgeqz.f"
	    }

/*           Compute Givens rotation on left */

/*           (  CQ   SQ ) */
/*           (  __      )  A or B */
/*           ( -SQ   CQ ) */

#line 999 "dhgeqz.f"
	    an = abs(a11) + abs(a12) + abs(a21) + abs(a22);
#line 1000 "dhgeqz.f"
	    bn = abs(b11) + abs(b22);
#line 1001 "dhgeqz.f"
	    wabs = abs(wr) + abs(wi);
#line 1002 "dhgeqz.f"
	    if (s1 * an > wabs * bn) {
#line 1003 "dhgeqz.f"
		cq = cz * b11;
#line 1004 "dhgeqz.f"
		sqr = szr * b22;
#line 1005 "dhgeqz.f"
		sqi = -szi * b22;
#line 1006 "dhgeqz.f"
	    } else {
#line 1007 "dhgeqz.f"
		a1r = cz * a11 + szr * a12;
#line 1008 "dhgeqz.f"
		a1i = szi * a12;
#line 1009 "dhgeqz.f"
		a2r = cz * a21 + szr * a22;
#line 1010 "dhgeqz.f"
		a2i = szi * a22;
#line 1011 "dhgeqz.f"
		cq = dlapy2_(&a1r, &a1i);
#line 1012 "dhgeqz.f"
		if (cq <= safmin) {
#line 1013 "dhgeqz.f"
		    cq = 0.;
#line 1014 "dhgeqz.f"
		    sqr = 1.;
#line 1015 "dhgeqz.f"
		    sqi = 0.;
#line 1016 "dhgeqz.f"
		} else {
#line 1017 "dhgeqz.f"
		    tempr = a1r / cq;
#line 1018 "dhgeqz.f"
		    tempi = a1i / cq;
#line 1019 "dhgeqz.f"
		    sqr = tempr * a2r + tempi * a2i;
#line 1020 "dhgeqz.f"
		    sqi = tempi * a2r - tempr * a2i;
#line 1021 "dhgeqz.f"
		}
#line 1022 "dhgeqz.f"
	    }
#line 1023 "dhgeqz.f"
	    t1 = dlapy3_(&cq, &sqr, &sqi);
#line 1024 "dhgeqz.f"
	    cq /= t1;
#line 1025 "dhgeqz.f"
	    sqr /= t1;
#line 1026 "dhgeqz.f"
	    sqi /= t1;

/*           Compute diagonal elements of QBZ */

#line 1030 "dhgeqz.f"
	    tempr = sqr * szr - sqi * szi;
#line 1031 "dhgeqz.f"
	    tempi = sqr * szi + sqi * szr;
#line 1032 "dhgeqz.f"
	    b1r = cq * cz * b11 + tempr * b22;
#line 1033 "dhgeqz.f"
	    b1i = tempi * b22;
#line 1034 "dhgeqz.f"
	    b1a = dlapy2_(&b1r, &b1i);
#line 1035 "dhgeqz.f"
	    b2r = cq * cz * b22 + tempr * b11;
#line 1036 "dhgeqz.f"
	    b2i = -tempi * b11;
#line 1037 "dhgeqz.f"
	    b2a = dlapy2_(&b2r, &b2i);

/*           Normalize so beta > 0, and Im( alpha1 ) > 0 */

#line 1041 "dhgeqz.f"
	    beta[ilast - 1] = b1a;
#line 1042 "dhgeqz.f"
	    beta[ilast] = b2a;
#line 1043 "dhgeqz.f"
	    alphar[ilast - 1] = wr * b1a * s1inv;
#line 1044 "dhgeqz.f"
	    alphai[ilast - 1] = wi * b1a * s1inv;
#line 1045 "dhgeqz.f"
	    alphar[ilast] = wr * b2a * s1inv;
#line 1046 "dhgeqz.f"
	    alphai[ilast] = -(wi * b2a) * s1inv;

/*           Step 3: Go to next block -- exit if finished. */

#line 1050 "dhgeqz.f"
	    ilast = ifirst - 1;
#line 1051 "dhgeqz.f"
	    if (ilast < *ilo) {
#line 1051 "dhgeqz.f"
		goto L380;
#line 1051 "dhgeqz.f"
	    }

/*           Reset counters */

#line 1056 "dhgeqz.f"
	    iiter = 0;
#line 1057 "dhgeqz.f"
	    eshift = 0.;
#line 1058 "dhgeqz.f"
	    if (! ilschr) {
#line 1059 "dhgeqz.f"
		ilastm = ilast;
#line 1060 "dhgeqz.f"
		if (ifrstm > ilast) {
#line 1060 "dhgeqz.f"
		    ifrstm = *ilo;
#line 1060 "dhgeqz.f"
		}
#line 1062 "dhgeqz.f"
	    }
#line 1063 "dhgeqz.f"
	    goto L350;
#line 1064 "dhgeqz.f"
	} else {

/*           Usual case: 3x3 or larger block, using Francis implicit */
/*                       double-shift */

/*                                    2 */
/*           Eigenvalue equation is  w  - c w + d = 0, */

/*                                         -1 2        -1 */
/*           so compute 1st column of  (A B  )  - c A B   + d */
/*           using the formula in QZIT (from EISPACK) */

/*           We assume that the block is at least 3x3 */

#line 1078 "dhgeqz.f"
	    ad11 = ascale * h__[ilast - 1 + (ilast - 1) * h_dim1] / (bscale * 
		    t[ilast - 1 + (ilast - 1) * t_dim1]);
#line 1080 "dhgeqz.f"
	    ad21 = ascale * h__[ilast + (ilast - 1) * h_dim1] / (bscale * t[
		    ilast - 1 + (ilast - 1) * t_dim1]);
#line 1082 "dhgeqz.f"
	    ad12 = ascale * h__[ilast - 1 + ilast * h_dim1] / (bscale * t[
		    ilast + ilast * t_dim1]);
#line 1084 "dhgeqz.f"
	    ad22 = ascale * h__[ilast + ilast * h_dim1] / (bscale * t[ilast + 
		    ilast * t_dim1]);
#line 1086 "dhgeqz.f"
	    u12 = t[ilast - 1 + ilast * t_dim1] / t[ilast + ilast * t_dim1];
#line 1087 "dhgeqz.f"
	    ad11l = ascale * h__[ifirst + ifirst * h_dim1] / (bscale * t[
		    ifirst + ifirst * t_dim1]);
#line 1089 "dhgeqz.f"
	    ad21l = ascale * h__[ifirst + 1 + ifirst * h_dim1] / (bscale * t[
		    ifirst + ifirst * t_dim1]);
#line 1091 "dhgeqz.f"
	    ad12l = ascale * h__[ifirst + (ifirst + 1) * h_dim1] / (bscale * 
		    t[ifirst + 1 + (ifirst + 1) * t_dim1]);
#line 1093 "dhgeqz.f"
	    ad22l = ascale * h__[ifirst + 1 + (ifirst + 1) * h_dim1] / (
		    bscale * t[ifirst + 1 + (ifirst + 1) * t_dim1]);
#line 1095 "dhgeqz.f"
	    ad32l = ascale * h__[ifirst + 2 + (ifirst + 1) * h_dim1] / (
		    bscale * t[ifirst + 1 + (ifirst + 1) * t_dim1]);
#line 1097 "dhgeqz.f"
	    u12l = t[ifirst + (ifirst + 1) * t_dim1] / t[ifirst + 1 + (ifirst 
		    + 1) * t_dim1];

#line 1099 "dhgeqz.f"
	    v[0] = (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * u12 
		    * ad11l + (ad12l - ad11l * u12l) * ad21l;
#line 1101 "dhgeqz.f"
	    v[1] = (ad22l - ad11l - ad21l * u12l - (ad11 - ad11l) - (ad22 - 
		    ad11l) + ad21 * u12) * ad21l;
#line 1103 "dhgeqz.f"
	    v[2] = ad32l * ad21l;

#line 1105 "dhgeqz.f"
	    istart = ifirst;

#line 1107 "dhgeqz.f"
	    dlarfg_(&c__3, v, &v[1], &c__1, &tau);
#line 1108 "dhgeqz.f"
	    v[0] = 1.;

/*           Sweep */

#line 1112 "dhgeqz.f"
	    i__2 = ilast - 2;
#line 1112 "dhgeqz.f"
	    for (j = istart; j <= i__2; ++j) {

/*              All but last elements: use 3x3 Householder transforms. */

/*              Zero (j-1)st column of A */

#line 1118 "dhgeqz.f"
		if (j > istart) {
#line 1119 "dhgeqz.f"
		    v[0] = h__[j + (j - 1) * h_dim1];
#line 1120 "dhgeqz.f"
		    v[1] = h__[j + 1 + (j - 1) * h_dim1];
#line 1121 "dhgeqz.f"
		    v[2] = h__[j + 2 + (j - 1) * h_dim1];

#line 1123 "dhgeqz.f"
		    dlarfg_(&c__3, &h__[j + (j - 1) * h_dim1], &v[1], &c__1, &
			    tau);
#line 1124 "dhgeqz.f"
		    v[0] = 1.;
#line 1125 "dhgeqz.f"
		    h__[j + 1 + (j - 1) * h_dim1] = 0.;
#line 1126 "dhgeqz.f"
		    h__[j + 2 + (j - 1) * h_dim1] = 0.;
#line 1127 "dhgeqz.f"
		}

#line 1129 "dhgeqz.f"
		i__3 = ilastm;
#line 1129 "dhgeqz.f"
		for (jc = j; jc <= i__3; ++jc) {
#line 1130 "dhgeqz.f"
		    temp = tau * (h__[j + jc * h_dim1] + v[1] * h__[j + 1 + 
			    jc * h_dim1] + v[2] * h__[j + 2 + jc * h_dim1]);
#line 1132 "dhgeqz.f"
		    h__[j + jc * h_dim1] -= temp;
#line 1133 "dhgeqz.f"
		    h__[j + 1 + jc * h_dim1] -= temp * v[1];
#line 1134 "dhgeqz.f"
		    h__[j + 2 + jc * h_dim1] -= temp * v[2];
#line 1135 "dhgeqz.f"
		    temp2 = tau * (t[j + jc * t_dim1] + v[1] * t[j + 1 + jc * 
			    t_dim1] + v[2] * t[j + 2 + jc * t_dim1]);
#line 1137 "dhgeqz.f"
		    t[j + jc * t_dim1] -= temp2;
#line 1138 "dhgeqz.f"
		    t[j + 1 + jc * t_dim1] -= temp2 * v[1];
#line 1139 "dhgeqz.f"
		    t[j + 2 + jc * t_dim1] -= temp2 * v[2];
#line 1140 "dhgeqz.f"
/* L230: */
#line 1140 "dhgeqz.f"
		}
#line 1141 "dhgeqz.f"
		if (ilq) {
#line 1142 "dhgeqz.f"
		    i__3 = *n;
#line 1142 "dhgeqz.f"
		    for (jr = 1; jr <= i__3; ++jr) {
#line 1143 "dhgeqz.f"
			temp = tau * (q[jr + j * q_dim1] + v[1] * q[jr + (j + 
				1) * q_dim1] + v[2] * q[jr + (j + 2) * q_dim1]
				);
#line 1145 "dhgeqz.f"
			q[jr + j * q_dim1] -= temp;
#line 1146 "dhgeqz.f"
			q[jr + (j + 1) * q_dim1] -= temp * v[1];
#line 1147 "dhgeqz.f"
			q[jr + (j + 2) * q_dim1] -= temp * v[2];
#line 1148 "dhgeqz.f"
/* L240: */
#line 1148 "dhgeqz.f"
		    }
#line 1149 "dhgeqz.f"
		}

/*              Zero j-th column of B (see DLAGBC for details) */

/*              Swap rows to pivot */

#line 1155 "dhgeqz.f"
		ilpivt = FALSE_;
/* Computing MAX */
#line 1156 "dhgeqz.f"
		d__3 = (d__1 = t[j + 1 + (j + 1) * t_dim1], abs(d__1)), d__4 =
			 (d__2 = t[j + 1 + (j + 2) * t_dim1], abs(d__2));
#line 1156 "dhgeqz.f"
		temp = max(d__3,d__4);
/* Computing MAX */
#line 1157 "dhgeqz.f"
		d__3 = (d__1 = t[j + 2 + (j + 1) * t_dim1], abs(d__1)), d__4 =
			 (d__2 = t[j + 2 + (j + 2) * t_dim1], abs(d__2));
#line 1157 "dhgeqz.f"
		temp2 = max(d__3,d__4);
#line 1158 "dhgeqz.f"
		if (max(temp,temp2) < safmin) {
#line 1159 "dhgeqz.f"
		    scale = 0.;
#line 1160 "dhgeqz.f"
		    u1 = 1.;
#line 1161 "dhgeqz.f"
		    u2 = 0.;
#line 1162 "dhgeqz.f"
		    goto L250;
#line 1163 "dhgeqz.f"
		} else if (temp >= temp2) {
#line 1164 "dhgeqz.f"
		    w11 = t[j + 1 + (j + 1) * t_dim1];
#line 1165 "dhgeqz.f"
		    w21 = t[j + 2 + (j + 1) * t_dim1];
#line 1166 "dhgeqz.f"
		    w12 = t[j + 1 + (j + 2) * t_dim1];
#line 1167 "dhgeqz.f"
		    w22 = t[j + 2 + (j + 2) * t_dim1];
#line 1168 "dhgeqz.f"
		    u1 = t[j + 1 + j * t_dim1];
#line 1169 "dhgeqz.f"
		    u2 = t[j + 2 + j * t_dim1];
#line 1170 "dhgeqz.f"
		} else {
#line 1171 "dhgeqz.f"
		    w21 = t[j + 1 + (j + 1) * t_dim1];
#line 1172 "dhgeqz.f"
		    w11 = t[j + 2 + (j + 1) * t_dim1];
#line 1173 "dhgeqz.f"
		    w22 = t[j + 1 + (j + 2) * t_dim1];
#line 1174 "dhgeqz.f"
		    w12 = t[j + 2 + (j + 2) * t_dim1];
#line 1175 "dhgeqz.f"
		    u2 = t[j + 1 + j * t_dim1];
#line 1176 "dhgeqz.f"
		    u1 = t[j + 2 + j * t_dim1];
#line 1177 "dhgeqz.f"
		}

/*              Swap columns if nec. */

#line 1181 "dhgeqz.f"
		if (abs(w12) > abs(w11)) {
#line 1182 "dhgeqz.f"
		    ilpivt = TRUE_;
#line 1183 "dhgeqz.f"
		    temp = w12;
#line 1184 "dhgeqz.f"
		    temp2 = w22;
#line 1185 "dhgeqz.f"
		    w12 = w11;
#line 1186 "dhgeqz.f"
		    w22 = w21;
#line 1187 "dhgeqz.f"
		    w11 = temp;
#line 1188 "dhgeqz.f"
		    w21 = temp2;
#line 1189 "dhgeqz.f"
		}

/*              LU-factor */

#line 1193 "dhgeqz.f"
		temp = w21 / w11;
#line 1194 "dhgeqz.f"
		u2 -= temp * u1;
#line 1195 "dhgeqz.f"
		w22 -= temp * w12;
#line 1196 "dhgeqz.f"
		w21 = 0.;

/*              Compute SCALE */

#line 1200 "dhgeqz.f"
		scale = 1.;
#line 1201 "dhgeqz.f"
		if (abs(w22) < safmin) {
#line 1202 "dhgeqz.f"
		    scale = 0.;
#line 1203 "dhgeqz.f"
		    u2 = 1.;
#line 1204 "dhgeqz.f"
		    u1 = -w12 / w11;
#line 1205 "dhgeqz.f"
		    goto L250;
#line 1206 "dhgeqz.f"
		}
#line 1207 "dhgeqz.f"
		if (abs(w22) < abs(u2)) {
#line 1207 "dhgeqz.f"
		    scale = (d__1 = w22 / u2, abs(d__1));
#line 1207 "dhgeqz.f"
		}
#line 1209 "dhgeqz.f"
		if (abs(w11) < abs(u1)) {
/* Computing MIN */
#line 1209 "dhgeqz.f"
		    d__2 = scale, d__3 = (d__1 = w11 / u1, abs(d__1));
#line 1209 "dhgeqz.f"
		    scale = min(d__2,d__3);
#line 1209 "dhgeqz.f"
		}

/*              Solve */

#line 1214 "dhgeqz.f"
		u2 = scale * u2 / w22;
#line 1215 "dhgeqz.f"
		u1 = (scale * u1 - w12 * u2) / w11;

#line 1217 "dhgeqz.f"
L250:
#line 1218 "dhgeqz.f"
		if (ilpivt) {
#line 1219 "dhgeqz.f"
		    temp = u2;
#line 1220 "dhgeqz.f"
		    u2 = u1;
#line 1221 "dhgeqz.f"
		    u1 = temp;
#line 1222 "dhgeqz.f"
		}

/*              Compute Householder Vector */

/* Computing 2nd power */
#line 1226 "dhgeqz.f"
		d__1 = scale;
/* Computing 2nd power */
#line 1226 "dhgeqz.f"
		d__2 = u1;
/* Computing 2nd power */
#line 1226 "dhgeqz.f"
		d__3 = u2;
#line 1226 "dhgeqz.f"
		t1 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
#line 1227 "dhgeqz.f"
		tau = scale / t1 + 1.;
#line 1228 "dhgeqz.f"
		vs = -1. / (scale + t1);
#line 1229 "dhgeqz.f"
		v[0] = 1.;
#line 1230 "dhgeqz.f"
		v[1] = vs * u1;
#line 1231 "dhgeqz.f"
		v[2] = vs * u2;

/*              Apply transformations from the right. */

/* Computing MIN */
#line 1235 "dhgeqz.f"
		i__4 = j + 3;
#line 1235 "dhgeqz.f"
		i__3 = min(i__4,ilast);
#line 1235 "dhgeqz.f"
		for (jr = ifrstm; jr <= i__3; ++jr) {
#line 1236 "dhgeqz.f"
		    temp = tau * (h__[jr + j * h_dim1] + v[1] * h__[jr + (j + 
			    1) * h_dim1] + v[2] * h__[jr + (j + 2) * h_dim1]);
#line 1238 "dhgeqz.f"
		    h__[jr + j * h_dim1] -= temp;
#line 1239 "dhgeqz.f"
		    h__[jr + (j + 1) * h_dim1] -= temp * v[1];
#line 1240 "dhgeqz.f"
		    h__[jr + (j + 2) * h_dim1] -= temp * v[2];
#line 1241 "dhgeqz.f"
/* L260: */
#line 1241 "dhgeqz.f"
		}
#line 1242 "dhgeqz.f"
		i__3 = j + 2;
#line 1242 "dhgeqz.f"
		for (jr = ifrstm; jr <= i__3; ++jr) {
#line 1243 "dhgeqz.f"
		    temp = tau * (t[jr + j * t_dim1] + v[1] * t[jr + (j + 1) *
			     t_dim1] + v[2] * t[jr + (j + 2) * t_dim1]);
#line 1245 "dhgeqz.f"
		    t[jr + j * t_dim1] -= temp;
#line 1246 "dhgeqz.f"
		    t[jr + (j + 1) * t_dim1] -= temp * v[1];
#line 1247 "dhgeqz.f"
		    t[jr + (j + 2) * t_dim1] -= temp * v[2];
#line 1248 "dhgeqz.f"
/* L270: */
#line 1248 "dhgeqz.f"
		}
#line 1249 "dhgeqz.f"
		if (ilz) {
#line 1250 "dhgeqz.f"
		    i__3 = *n;
#line 1250 "dhgeqz.f"
		    for (jr = 1; jr <= i__3; ++jr) {
#line 1251 "dhgeqz.f"
			temp = tau * (z__[jr + j * z_dim1] + v[1] * z__[jr + (
				j + 1) * z_dim1] + v[2] * z__[jr + (j + 2) * 
				z_dim1]);
#line 1253 "dhgeqz.f"
			z__[jr + j * z_dim1] -= temp;
#line 1254 "dhgeqz.f"
			z__[jr + (j + 1) * z_dim1] -= temp * v[1];
#line 1255 "dhgeqz.f"
			z__[jr + (j + 2) * z_dim1] -= temp * v[2];
#line 1256 "dhgeqz.f"
/* L280: */
#line 1256 "dhgeqz.f"
		    }
#line 1257 "dhgeqz.f"
		}
#line 1258 "dhgeqz.f"
		t[j + 1 + j * t_dim1] = 0.;
#line 1259 "dhgeqz.f"
		t[j + 2 + j * t_dim1] = 0.;
#line 1260 "dhgeqz.f"
/* L290: */
#line 1260 "dhgeqz.f"
	    }

/*           Last elements: Use Givens rotations */

/*           Rotations from the left */

#line 1266 "dhgeqz.f"
	    j = ilast - 1;
#line 1267 "dhgeqz.f"
	    temp = h__[j + (j - 1) * h_dim1];
#line 1268 "dhgeqz.f"
	    dlartg_(&temp, &h__[j + 1 + (j - 1) * h_dim1], &c__, &s, &h__[j + 
		    (j - 1) * h_dim1]);
#line 1269 "dhgeqz.f"
	    h__[j + 1 + (j - 1) * h_dim1] = 0.;

#line 1271 "dhgeqz.f"
	    i__2 = ilastm;
#line 1271 "dhgeqz.f"
	    for (jc = j; jc <= i__2; ++jc) {
#line 1272 "dhgeqz.f"
		temp = c__ * h__[j + jc * h_dim1] + s * h__[j + 1 + jc * 
			h_dim1];
#line 1273 "dhgeqz.f"
		h__[j + 1 + jc * h_dim1] = -s * h__[j + jc * h_dim1] + c__ * 
			h__[j + 1 + jc * h_dim1];
#line 1274 "dhgeqz.f"
		h__[j + jc * h_dim1] = temp;
#line 1275 "dhgeqz.f"
		temp2 = c__ * t[j + jc * t_dim1] + s * t[j + 1 + jc * t_dim1];
#line 1276 "dhgeqz.f"
		t[j + 1 + jc * t_dim1] = -s * t[j + jc * t_dim1] + c__ * t[j 
			+ 1 + jc * t_dim1];
#line 1277 "dhgeqz.f"
		t[j + jc * t_dim1] = temp2;
#line 1278 "dhgeqz.f"
/* L300: */
#line 1278 "dhgeqz.f"
	    }
#line 1279 "dhgeqz.f"
	    if (ilq) {
#line 1280 "dhgeqz.f"
		i__2 = *n;
#line 1280 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 1281 "dhgeqz.f"
		    temp = c__ * q[jr + j * q_dim1] + s * q[jr + (j + 1) * 
			    q_dim1];
#line 1282 "dhgeqz.f"
		    q[jr + (j + 1) * q_dim1] = -s * q[jr + j * q_dim1] + c__ *
			     q[jr + (j + 1) * q_dim1];
#line 1283 "dhgeqz.f"
		    q[jr + j * q_dim1] = temp;
#line 1284 "dhgeqz.f"
/* L310: */
#line 1284 "dhgeqz.f"
		}
#line 1285 "dhgeqz.f"
	    }

/*           Rotations from the right. */

#line 1289 "dhgeqz.f"
	    temp = t[j + 1 + (j + 1) * t_dim1];
#line 1290 "dhgeqz.f"
	    dlartg_(&temp, &t[j + 1 + j * t_dim1], &c__, &s, &t[j + 1 + (j + 
		    1) * t_dim1]);
#line 1291 "dhgeqz.f"
	    t[j + 1 + j * t_dim1] = 0.;

#line 1293 "dhgeqz.f"
	    i__2 = ilast;
#line 1293 "dhgeqz.f"
	    for (jr = ifrstm; jr <= i__2; ++jr) {
#line 1294 "dhgeqz.f"
		temp = c__ * h__[jr + (j + 1) * h_dim1] + s * h__[jr + j * 
			h_dim1];
#line 1295 "dhgeqz.f"
		h__[jr + j * h_dim1] = -s * h__[jr + (j + 1) * h_dim1] + c__ *
			 h__[jr + j * h_dim1];
#line 1296 "dhgeqz.f"
		h__[jr + (j + 1) * h_dim1] = temp;
#line 1297 "dhgeqz.f"
/* L320: */
#line 1297 "dhgeqz.f"
	    }
#line 1298 "dhgeqz.f"
	    i__2 = ilast - 1;
#line 1298 "dhgeqz.f"
	    for (jr = ifrstm; jr <= i__2; ++jr) {
#line 1299 "dhgeqz.f"
		temp = c__ * t[jr + (j + 1) * t_dim1] + s * t[jr + j * t_dim1]
			;
#line 1300 "dhgeqz.f"
		t[jr + j * t_dim1] = -s * t[jr + (j + 1) * t_dim1] + c__ * t[
			jr + j * t_dim1];
#line 1301 "dhgeqz.f"
		t[jr + (j + 1) * t_dim1] = temp;
#line 1302 "dhgeqz.f"
/* L330: */
#line 1302 "dhgeqz.f"
	    }
#line 1303 "dhgeqz.f"
	    if (ilz) {
#line 1304 "dhgeqz.f"
		i__2 = *n;
#line 1304 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 1305 "dhgeqz.f"
		    temp = c__ * z__[jr + (j + 1) * z_dim1] + s * z__[jr + j *
			     z_dim1];
#line 1306 "dhgeqz.f"
		    z__[jr + j * z_dim1] = -s * z__[jr + (j + 1) * z_dim1] + 
			    c__ * z__[jr + j * z_dim1];
#line 1307 "dhgeqz.f"
		    z__[jr + (j + 1) * z_dim1] = temp;
#line 1308 "dhgeqz.f"
/* L340: */
#line 1308 "dhgeqz.f"
		}
#line 1309 "dhgeqz.f"
	    }

/*           End of Double-Shift code */

#line 1313 "dhgeqz.f"
	}

#line 1315 "dhgeqz.f"
	goto L350;

/*        End of iteration loop */

#line 1319 "dhgeqz.f"
L350:
#line 1320 "dhgeqz.f"
/* L360: */
#line 1320 "dhgeqz.f"
	;
#line 1320 "dhgeqz.f"
    }

/*     Drop-through = non-convergence */

#line 1324 "dhgeqz.f"
    *info = ilast;
#line 1325 "dhgeqz.f"
    goto L420;

/*     Successful completion of all QZ steps */

#line 1329 "dhgeqz.f"
L380:

/*     Set Eigenvalues 1:ILO-1 */

#line 1333 "dhgeqz.f"
    i__1 = *ilo - 1;
#line 1333 "dhgeqz.f"
    for (j = 1; j <= i__1; ++j) {
#line 1334 "dhgeqz.f"
	if (t[j + j * t_dim1] < 0.) {
#line 1335 "dhgeqz.f"
	    if (ilschr) {
#line 1336 "dhgeqz.f"
		i__2 = j;
#line 1336 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 1337 "dhgeqz.f"
		    h__[jr + j * h_dim1] = -h__[jr + j * h_dim1];
#line 1338 "dhgeqz.f"
		    t[jr + j * t_dim1] = -t[jr + j * t_dim1];
#line 1339 "dhgeqz.f"
/* L390: */
#line 1339 "dhgeqz.f"
		}
#line 1340 "dhgeqz.f"
	    } else {
#line 1341 "dhgeqz.f"
		h__[j + j * h_dim1] = -h__[j + j * h_dim1];
#line 1342 "dhgeqz.f"
		t[j + j * t_dim1] = -t[j + j * t_dim1];
#line 1343 "dhgeqz.f"
	    }
#line 1344 "dhgeqz.f"
	    if (ilz) {
#line 1345 "dhgeqz.f"
		i__2 = *n;
#line 1345 "dhgeqz.f"
		for (jr = 1; jr <= i__2; ++jr) {
#line 1346 "dhgeqz.f"
		    z__[jr + j * z_dim1] = -z__[jr + j * z_dim1];
#line 1347 "dhgeqz.f"
/* L400: */
#line 1347 "dhgeqz.f"
		}
#line 1348 "dhgeqz.f"
	    }
#line 1349 "dhgeqz.f"
	}
#line 1350 "dhgeqz.f"
	alphar[j] = h__[j + j * h_dim1];
#line 1351 "dhgeqz.f"
	alphai[j] = 0.;
#line 1352 "dhgeqz.f"
	beta[j] = t[j + j * t_dim1];
#line 1353 "dhgeqz.f"
/* L410: */
#line 1353 "dhgeqz.f"
    }

/*     Normal Termination */

#line 1357 "dhgeqz.f"
    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size */

#line 1361 "dhgeqz.f"
L420:
#line 1362 "dhgeqz.f"
    work[1] = (doublereal) (*n);
#line 1363 "dhgeqz.f"
    return 0;

/*     End of DHGEQZ */

} /* dhgeqz_ */

