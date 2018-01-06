#line 1 "chgeqz.f"
/* chgeqz.f -- translated by f2c (version 20100827).
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

#line 1 "chgeqz.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b CHGEQZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHGEQZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chgeqz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chgeqz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chgeqz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, */
/*                          ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ, JOB */
/*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            ALPHA( * ), BETA( * ), H( LDH, * ), */
/*      $                   Q( LDQ, * ), T( LDT, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHGEQZ computes the eigenvalues of a complex matrix pair (H,T), */
/* > where H is an upper Hessenberg matrix and T is upper triangular, */
/* > using the single-shift QZ method. */
/* > Matrix pairs of this type are produced by the reduction to */
/* > generalized upper Hessenberg form of a complex matrix pair (A,B): */
/* > */
/* >    A = Q1*H*Z1**H,  B = Q1*T*Z1**H, */
/* > */
/* > as computed by CGGHRD. */
/* > */
/* > If JOB='S', then the Hessenberg-triangular pair (H,T) is */
/* > also reduced to generalized Schur form, */
/* > */
/* >    H = Q*S*Z**H,  T = Q*P*Z**H, */
/* > */
/* > where Q and Z are unitary matrices and S and P are upper triangular. */
/* > */
/* > Optionally, the unitary matrix Q from the generalized Schur */
/* > factorization may be postmultiplied into an input matrix Q1, and the */
/* > unitary matrix Z may be postmultiplied into an input matrix Z1. */
/* > If Q1 and Z1 are the unitary matrices from CGGHRD that reduced */
/* > the matrix pair (A,B) to generalized Hessenberg form, then the output */
/* > matrices Q1*Q and Z1*Z are the unitary factors from the generalized */
/* > Schur factorization of (A,B): */
/* > */
/* >    A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H. */
/* > */
/* > To avoid overflow, eigenvalues of the matrix pair (H,T) */
/* > (equivalently, of (A,B)) are computed as a pair of complex values */
/* > (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an */
/* > eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP) */
/* >    A*x = lambda*B*x */
/* > and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the */
/* > alternate form of the GNEP */
/* >    mu*A*y = B*y. */
/* > The values of alpha and beta for the i-th eigenvalue can be read */
/* > directly from the generalized Schur form:  alpha = S(i,i), */
/* > beta = P(i,i). */
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
/* >          = 'S': Computer eigenvalues and the Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'N': Left Schur vectors (Q) are not computed; */
/* >          = 'I': Q is initialized to the unit matrix and the matrix Q */
/* >                 of left Schur vectors of (H,T) is returned; */
/* >          = 'V': Q must contain a unitary matrix Q1 on entry and */
/* >                 the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N': Right Schur vectors (Z) are not computed; */
/* >          = 'I': Q is initialized to the unit matrix and the matrix Z */
/* >                 of right Schur vectors of (H,T) is returned; */
/* >          = 'V': Z must contain a unitary matrix Z1 on entry and */
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
/* >          H is COMPLEX array, dimension (LDH, N) */
/* >          On entry, the N-by-N upper Hessenberg matrix H. */
/* >          On exit, if JOB = 'S', H contains the upper triangular */
/* >          matrix S from the generalized Schur factorization. */
/* >          If JOB = 'E', the diagonal of H matches that of S, but */
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
/* >          T is COMPLEX array, dimension (LDT, N) */
/* >          On entry, the N-by-N upper triangular matrix T. */
/* >          On exit, if JOB = 'S', T contains the upper triangular */
/* >          matrix P from the generalized Schur factorization. */
/* >          If JOB = 'E', the diagonal of T matches that of P, but */
/* >          the rest of T is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX array, dimension (N) */
/* >          The complex scalars alpha that define the eigenvalues of */
/* >          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur */
/* >          factorization. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX array, dimension (N) */
/* >          The real non-negative scalars beta that define the */
/* >          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized */
/* >          Schur factorization. */
/* > */
/* >          Together, the quantities alpha = ALPHA(j) and beta = BETA(j) */
/* >          represent the j-th eigenvalue of the matrix pair (A,B), in */
/* >          one of the forms lambda = alpha/beta or mu = beta/alpha. */
/* >          Since either lambda or mu may overflow, they should not, */
/* >          in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ, N) */
/* >          On entry, if COMPZ = 'V', the unitary matrix Q1 used in the */
/* >          reduction of (A,B) to generalized Hessenberg form. */
/* >          On exit, if COMPZ = 'I', the unitary matrix of left Schur */
/* >          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of */
/* >          left Schur vectors of (A,B). */
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
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the */
/* >          reduction of (A,B) to generalized Hessenberg form. */
/* >          On exit, if COMPZ = 'I', the unitary matrix of right Schur */
/* >          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of */
/* >          right Schur vectors of (A,B). */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          = 1,...,N: the QZ iteration did not converge.  (H,T) is not */
/* >                     in Schur form, but ALPHA(i) and BETA(i), */
/* >                     i=INFO+1,...,N should be correct. */
/* >          = N+1,...,2*N: the shift calculation failed.  (H,T) is not */
/* >                     in Schur form, but ALPHA(i) and BETA(i), */
/* >                     i=INFO-N+1,...,N should be correct. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup complexGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  We assume that complex ABS works as long as its value is less than */
/* >  overflow. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *t, integer *ldt, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	info, ftnlen job_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, q_dim1, q_offset, t_dim1, t_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), pow_zi(
	    doublecomplex *, doublecomplex *, integer *), z_sqrt(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublecomplex s, t1;
    static integer jc, in;
    static doublecomplex u12;
    static integer jr;
    static doublecomplex ad11, ad12, ad21, ad22;
    static integer jch;
    static logical ilq, ilz;
    static doublereal ulp;
    static doublecomplex abi22;
    static doublereal absb, atol, btol, temp;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal temp2;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex ctemp;
    static integer iiter, ilast, jiter;
    static doublereal anorm, bnorm;
    static integer maxit;
    static doublecomplex shift;
    static doublereal tempr;
    static doublecomplex ctemp2, ctemp3;
    static logical ilazr2;
    static doublereal ascale, bscale;
    static doublecomplex signbc;
    extern doublereal slamch_(char *, ftnlen), clanhs_(char *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublecomplex eshift;
    static logical ilschr;
    static integer icompq, ilastm;
    static doublecomplex rtdisc;
    static integer ischur;
    static logical ilazro;
    static integer icompz, ifirst, ifrstm, istart;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode JOB, COMPQ, COMPZ */

#line 346 "chgeqz.f"
    /* Parameter adjustments */
#line 346 "chgeqz.f"
    h_dim1 = *ldh;
#line 346 "chgeqz.f"
    h_offset = 1 + h_dim1;
#line 346 "chgeqz.f"
    h__ -= h_offset;
#line 346 "chgeqz.f"
    t_dim1 = *ldt;
#line 346 "chgeqz.f"
    t_offset = 1 + t_dim1;
#line 346 "chgeqz.f"
    t -= t_offset;
#line 346 "chgeqz.f"
    --alpha;
#line 346 "chgeqz.f"
    --beta;
#line 346 "chgeqz.f"
    q_dim1 = *ldq;
#line 346 "chgeqz.f"
    q_offset = 1 + q_dim1;
#line 346 "chgeqz.f"
    q -= q_offset;
#line 346 "chgeqz.f"
    z_dim1 = *ldz;
#line 346 "chgeqz.f"
    z_offset = 1 + z_dim1;
#line 346 "chgeqz.f"
    z__ -= z_offset;
#line 346 "chgeqz.f"
    --work;
#line 346 "chgeqz.f"
    --rwork;
#line 346 "chgeqz.f"

#line 346 "chgeqz.f"
    /* Function Body */
#line 346 "chgeqz.f"
    if (lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
#line 347 "chgeqz.f"
	ilschr = FALSE_;
#line 348 "chgeqz.f"
	ischur = 1;
#line 349 "chgeqz.f"
    } else if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 350 "chgeqz.f"
	ilschr = TRUE_;
#line 351 "chgeqz.f"
	ischur = 2;
#line 352 "chgeqz.f"
    } else {
#line 353 "chgeqz.f"
	ischur = 0;
#line 354 "chgeqz.f"
    }

#line 356 "chgeqz.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 357 "chgeqz.f"
	ilq = FALSE_;
#line 358 "chgeqz.f"
	icompq = 1;
#line 359 "chgeqz.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 360 "chgeqz.f"
	ilq = TRUE_;
#line 361 "chgeqz.f"
	icompq = 2;
#line 362 "chgeqz.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 363 "chgeqz.f"
	ilq = TRUE_;
#line 364 "chgeqz.f"
	icompq = 3;
#line 365 "chgeqz.f"
    } else {
#line 366 "chgeqz.f"
	icompq = 0;
#line 367 "chgeqz.f"
    }

#line 369 "chgeqz.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 370 "chgeqz.f"
	ilz = FALSE_;
#line 371 "chgeqz.f"
	icompz = 1;
#line 372 "chgeqz.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 373 "chgeqz.f"
	ilz = TRUE_;
#line 374 "chgeqz.f"
	icompz = 2;
#line 375 "chgeqz.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 376 "chgeqz.f"
	ilz = TRUE_;
#line 377 "chgeqz.f"
	icompz = 3;
#line 378 "chgeqz.f"
    } else {
#line 379 "chgeqz.f"
	icompz = 0;
#line 380 "chgeqz.f"
    }

/*     Check Argument Values */

#line 384 "chgeqz.f"
    *info = 0;
#line 385 "chgeqz.f"
    i__1 = max(1,*n);
#line 385 "chgeqz.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 386 "chgeqz.f"
    lquery = *lwork == -1;
#line 387 "chgeqz.f"
    if (ischur == 0) {
#line 388 "chgeqz.f"
	*info = -1;
#line 389 "chgeqz.f"
    } else if (icompq == 0) {
#line 390 "chgeqz.f"
	*info = -2;
#line 391 "chgeqz.f"
    } else if (icompz == 0) {
#line 392 "chgeqz.f"
	*info = -3;
#line 393 "chgeqz.f"
    } else if (*n < 0) {
#line 394 "chgeqz.f"
	*info = -4;
#line 395 "chgeqz.f"
    } else if (*ilo < 1) {
#line 396 "chgeqz.f"
	*info = -5;
#line 397 "chgeqz.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 398 "chgeqz.f"
	*info = -6;
#line 399 "chgeqz.f"
    } else if (*ldh < *n) {
#line 400 "chgeqz.f"
	*info = -8;
#line 401 "chgeqz.f"
    } else if (*ldt < *n) {
#line 402 "chgeqz.f"
	*info = -10;
#line 403 "chgeqz.f"
    } else if (*ldq < 1 || ilq && *ldq < *n) {
#line 404 "chgeqz.f"
	*info = -14;
#line 405 "chgeqz.f"
    } else if (*ldz < 1 || ilz && *ldz < *n) {
#line 406 "chgeqz.f"
	*info = -16;
#line 407 "chgeqz.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 408 "chgeqz.f"
	*info = -18;
#line 409 "chgeqz.f"
    }
#line 410 "chgeqz.f"
    if (*info != 0) {
#line 411 "chgeqz.f"
	i__1 = -(*info);
#line 411 "chgeqz.f"
	xerbla_("CHGEQZ", &i__1, (ftnlen)6);
#line 412 "chgeqz.f"
	return 0;
#line 413 "chgeqz.f"
    } else if (lquery) {
#line 414 "chgeqz.f"
	return 0;
#line 415 "chgeqz.f"
    }

/*     Quick return if possible */

/*     WORK( 1 ) = CMPLX( 1 ) */
#line 420 "chgeqz.f"
    if (*n <= 0) {
#line 421 "chgeqz.f"
	work[1].r = 1., work[1].i = 0.;
#line 422 "chgeqz.f"
	return 0;
#line 423 "chgeqz.f"
    }

/*     Initialize Q and Z */

#line 427 "chgeqz.f"
    if (icompq == 3) {
#line 427 "chgeqz.f"
	claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 427 "chgeqz.f"
    }
#line 429 "chgeqz.f"
    if (icompz == 3) {
#line 429 "chgeqz.f"
	claset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 429 "chgeqz.f"
    }

/*     Machine Constants */

#line 434 "chgeqz.f"
    in = *ihi + 1 - *ilo;
#line 435 "chgeqz.f"
    safmin = slamch_("S", (ftnlen)1);
#line 436 "chgeqz.f"
    ulp = slamch_("E", (ftnlen)1) * slamch_("B", (ftnlen)1);
#line 437 "chgeqz.f"
    anorm = clanhs_("F", &in, &h__[*ilo + *ilo * h_dim1], ldh, &rwork[1], (
	    ftnlen)1);
#line 438 "chgeqz.f"
    bnorm = clanhs_("F", &in, &t[*ilo + *ilo * t_dim1], ldt, &rwork[1], (
	    ftnlen)1);
/* Computing MAX */
#line 439 "chgeqz.f"
    d__1 = safmin, d__2 = ulp * anorm;
#line 439 "chgeqz.f"
    atol = max(d__1,d__2);
/* Computing MAX */
#line 440 "chgeqz.f"
    d__1 = safmin, d__2 = ulp * bnorm;
#line 440 "chgeqz.f"
    btol = max(d__1,d__2);
#line 441 "chgeqz.f"
    ascale = 1. / max(safmin,anorm);
#line 442 "chgeqz.f"
    bscale = 1. / max(safmin,bnorm);


/*     Set Eigenvalues IHI+1:N */

#line 447 "chgeqz.f"
    i__1 = *n;
#line 447 "chgeqz.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 448 "chgeqz.f"
	absb = z_abs(&t[j + j * t_dim1]);
#line 449 "chgeqz.f"
	if (absb > safmin) {
#line 450 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 450 "chgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 450 "chgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 450 "chgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 451 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 451 "chgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 452 "chgeqz.f"
	    if (ilschr) {
#line 453 "chgeqz.f"
		i__2 = j - 1;
#line 453 "chgeqz.f"
		cscal_(&i__2, &signbc, &t[j * t_dim1 + 1], &c__1);
#line 454 "chgeqz.f"
		cscal_(&j, &signbc, &h__[j * h_dim1 + 1], &c__1);
#line 455 "chgeqz.f"
	    } else {
#line 456 "chgeqz.f"
		i__2 = j + j * h_dim1;
#line 456 "chgeqz.f"
		i__3 = j + j * h_dim1;
#line 456 "chgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 456 "chgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 457 "chgeqz.f"
	    }
#line 458 "chgeqz.f"
	    if (ilz) {
#line 458 "chgeqz.f"
		cscal_(n, &signbc, &z__[j * z_dim1 + 1], &c__1);
#line 458 "chgeqz.f"
	    }
#line 460 "chgeqz.f"
	} else {
#line 461 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 461 "chgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 462 "chgeqz.f"
	}
#line 463 "chgeqz.f"
	i__2 = j;
#line 463 "chgeqz.f"
	i__3 = j + j * h_dim1;
#line 463 "chgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 464 "chgeqz.f"
	i__2 = j;
#line 464 "chgeqz.f"
	i__3 = j + j * t_dim1;
#line 464 "chgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;
#line 465 "chgeqz.f"
/* L10: */
#line 465 "chgeqz.f"
    }

/*     If IHI < ILO, skip QZ steps */

#line 469 "chgeqz.f"
    if (*ihi < *ilo) {
#line 469 "chgeqz.f"
	goto L190;
#line 469 "chgeqz.f"
    }

/*     MAIN QZ ITERATION LOOP */

/*     Initialize dynamic indices */

/*     Eigenvalues ILAST+1:N have been found. */
/*        Column operations modify rows IFRSTM:whatever */
/*        Row operations modify columns whatever:ILASTM */

/*     If only eigenvalues are being computed, then */
/*        IFRSTM is the row of the last splitting row above row ILAST; */
/*        this is always at least ILO. */
/*     IITER counts iterations since the last eigenvalue was found, */
/*        to tell when to use an extraordinary shift. */
/*     MAXIT is the maximum number of QZ sweeps allowed. */

#line 487 "chgeqz.f"
    ilast = *ihi;
#line 488 "chgeqz.f"
    if (ilschr) {
#line 489 "chgeqz.f"
	ifrstm = 1;
#line 490 "chgeqz.f"
	ilastm = *n;
#line 491 "chgeqz.f"
    } else {
#line 492 "chgeqz.f"
	ifrstm = *ilo;
#line 493 "chgeqz.f"
	ilastm = *ihi;
#line 494 "chgeqz.f"
    }
#line 495 "chgeqz.f"
    iiter = 0;
#line 496 "chgeqz.f"
    eshift.r = 0., eshift.i = 0.;
#line 497 "chgeqz.f"
    maxit = (*ihi - *ilo + 1) * 30;

#line 499 "chgeqz.f"
    i__1 = maxit;
#line 499 "chgeqz.f"
    for (jiter = 1; jiter <= i__1; ++jiter) {

/*        Check for too many iterations. */

#line 503 "chgeqz.f"
	if (jiter > maxit) {
#line 503 "chgeqz.f"
	    goto L180;
#line 503 "chgeqz.f"
	}

/*        Split the matrix if possible. */

/*        Two tests: */
/*           1: H(j,j-1)=0  or  j=ILO */
/*           2: T(j,j)=0 */

/*        Special case: j=ILAST */

#line 514 "chgeqz.f"
	if (ilast == *ilo) {
#line 515 "chgeqz.f"
	    goto L60;
#line 516 "chgeqz.f"
	} else {
#line 517 "chgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 517 "chgeqz.f"
	    if ((d__1 = h__[i__2].r, abs(d__1)) + (d__2 = d_imag(&h__[ilast + 
		    (ilast - 1) * h_dim1]), abs(d__2)) <= atol) {
#line 518 "chgeqz.f"
		i__2 = ilast + (ilast - 1) * h_dim1;
#line 518 "chgeqz.f"
		h__[i__2].r = 0., h__[i__2].i = 0.;
#line 519 "chgeqz.f"
		goto L60;
#line 520 "chgeqz.f"
	    }
#line 521 "chgeqz.f"
	}

#line 523 "chgeqz.f"
	if (z_abs(&t[ilast + ilast * t_dim1]) <= btol) {
#line 524 "chgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 524 "chgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 525 "chgeqz.f"
	    goto L50;
#line 526 "chgeqz.f"
	}

/*        General case: j<ILAST */

#line 530 "chgeqz.f"
	i__2 = *ilo;
#line 530 "chgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {

/*           Test 1: for H(j,j-1)=0 or j=ILO */

#line 534 "chgeqz.f"
	    if (j == *ilo) {
#line 535 "chgeqz.f"
		ilazro = TRUE_;
#line 536 "chgeqz.f"
	    } else {
#line 537 "chgeqz.f"
		i__3 = j + (j - 1) * h_dim1;
#line 537 "chgeqz.f"
		if ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[j + 
			(j - 1) * h_dim1]), abs(d__2)) <= atol) {
#line 538 "chgeqz.f"
		    i__3 = j + (j - 1) * h_dim1;
#line 538 "chgeqz.f"
		    h__[i__3].r = 0., h__[i__3].i = 0.;
#line 539 "chgeqz.f"
		    ilazro = TRUE_;
#line 540 "chgeqz.f"
		} else {
#line 541 "chgeqz.f"
		    ilazro = FALSE_;
#line 542 "chgeqz.f"
		}
#line 543 "chgeqz.f"
	    }

/*           Test 2: for T(j,j)=0 */

#line 547 "chgeqz.f"
	    if (z_abs(&t[j + j * t_dim1]) < btol) {
#line 548 "chgeqz.f"
		i__3 = j + j * t_dim1;
#line 548 "chgeqz.f"
		t[i__3].r = 0., t[i__3].i = 0.;

/*              Test 1a: Check for 2 consecutive small subdiagonals in A */

#line 552 "chgeqz.f"
		ilazr2 = FALSE_;
#line 553 "chgeqz.f"
		if (! ilazro) {
#line 554 "chgeqz.f"
		    i__3 = j + (j - 1) * h_dim1;
#line 554 "chgeqz.f"
		    i__4 = j + 1 + j * h_dim1;
#line 554 "chgeqz.f"
		    i__5 = j + j * h_dim1;
#line 554 "chgeqz.f"
		    if (((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
			    h__[j + (j - 1) * h_dim1]), abs(d__2))) * (ascale 
			    * ((d__3 = h__[i__4].r, abs(d__3)) + (d__4 = 
			    d_imag(&h__[j + 1 + j * h_dim1]), abs(d__4)))) <= 
			    ((d__5 = h__[i__5].r, abs(d__5)) + (d__6 = d_imag(
			    &h__[j + j * h_dim1]), abs(d__6))) * (ascale * 
			    atol)) {
#line 554 "chgeqz.f"
			ilazr2 = TRUE_;
#line 554 "chgeqz.f"
		    }
#line 557 "chgeqz.f"
		}

/*              If both tests pass (1 & 2), i.e., the leading diagonal */
/*              element of B in the block is zero, split a 1x1 block off */
/*              at the top. (I.e., at the J-th row/column) The leading */
/*              diagonal element of the remainder can also be zero, so */
/*              this may have to be done repeatedly. */

#line 565 "chgeqz.f"
		if (ilazro || ilazr2) {
#line 566 "chgeqz.f"
		    i__3 = ilast - 1;
#line 566 "chgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 567 "chgeqz.f"
			i__4 = jch + jch * h_dim1;
#line 567 "chgeqz.f"
			ctemp.r = h__[i__4].r, ctemp.i = h__[i__4].i;
#line 568 "chgeqz.f"
			clartg_(&ctemp, &h__[jch + 1 + jch * h_dim1], &c__, &
				s, &h__[jch + jch * h_dim1]);
#line 570 "chgeqz.f"
			i__4 = jch + 1 + jch * h_dim1;
#line 570 "chgeqz.f"
			h__[i__4].r = 0., h__[i__4].i = 0.;
#line 571 "chgeqz.f"
			i__4 = ilastm - jch;
#line 571 "chgeqz.f"
			crot_(&i__4, &h__[jch + (jch + 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch + 1) * h_dim1], ldh, &c__, 
				&s);
#line 573 "chgeqz.f"
			i__4 = ilastm - jch;
#line 573 "chgeqz.f"
			crot_(&i__4, &t[jch + (jch + 1) * t_dim1], ldt, &t[
				jch + 1 + (jch + 1) * t_dim1], ldt, &c__, &s);
#line 575 "chgeqz.f"
			if (ilq) {
#line 575 "chgeqz.f"
			    d_cnjg(&z__1, &s);
#line 575 "chgeqz.f"
			    crot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &z__1);
#line 575 "chgeqz.f"
			}
#line 578 "chgeqz.f"
			if (ilazr2) {
#line 578 "chgeqz.f"
			    i__4 = jch + (jch - 1) * h_dim1;
#line 578 "chgeqz.f"
			    i__5 = jch + (jch - 1) * h_dim1;
#line 578 "chgeqz.f"
			    z__1.r = c__ * h__[i__5].r, z__1.i = c__ * h__[
				    i__5].i;
#line 578 "chgeqz.f"
			    h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 578 "chgeqz.f"
			}
#line 580 "chgeqz.f"
			ilazr2 = FALSE_;
#line 581 "chgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 581 "chgeqz.f"
			if ((d__1 = t[i__4].r, abs(d__1)) + (d__2 = d_imag(&t[
				jch + 1 + (jch + 1) * t_dim1]), abs(d__2)) >= 
				btol) {
#line 582 "chgeqz.f"
			    if (jch + 1 >= ilast) {
#line 583 "chgeqz.f"
				goto L60;
#line 584 "chgeqz.f"
			    } else {
#line 585 "chgeqz.f"
				ifirst = jch + 1;
#line 586 "chgeqz.f"
				goto L70;
#line 587 "chgeqz.f"
			    }
#line 588 "chgeqz.f"
			}
#line 589 "chgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 589 "chgeqz.f"
			t[i__4].r = 0., t[i__4].i = 0.;
#line 590 "chgeqz.f"
/* L20: */
#line 590 "chgeqz.f"
		    }
#line 591 "chgeqz.f"
		    goto L50;
#line 592 "chgeqz.f"
		} else {

/*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST) */
/*                 Then process as in the case T(ILAST,ILAST)=0 */

#line 597 "chgeqz.f"
		    i__3 = ilast - 1;
#line 597 "chgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 598 "chgeqz.f"
			i__4 = jch + (jch + 1) * t_dim1;
#line 598 "chgeqz.f"
			ctemp.r = t[i__4].r, ctemp.i = t[i__4].i;
#line 599 "chgeqz.f"
			clartg_(&ctemp, &t[jch + 1 + (jch + 1) * t_dim1], &
				c__, &s, &t[jch + (jch + 1) * t_dim1]);
#line 601 "chgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 601 "chgeqz.f"
			t[i__4].r = 0., t[i__4].i = 0.;
#line 602 "chgeqz.f"
			if (jch < ilastm - 1) {
#line 602 "chgeqz.f"
			    i__4 = ilastm - jch - 1;
#line 602 "chgeqz.f"
			    crot_(&i__4, &t[jch + (jch + 2) * t_dim1], ldt, &
				    t[jch + 1 + (jch + 2) * t_dim1], ldt, &
				    c__, &s);
#line 602 "chgeqz.f"
			}
#line 605 "chgeqz.f"
			i__4 = ilastm - jch + 2;
#line 605 "chgeqz.f"
			crot_(&i__4, &h__[jch + (jch - 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch - 1) * h_dim1], ldh, &c__, 
				&s);
#line 607 "chgeqz.f"
			if (ilq) {
#line 607 "chgeqz.f"
			    d_cnjg(&z__1, &s);
#line 607 "chgeqz.f"
			    crot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &z__1);
#line 607 "chgeqz.f"
			}
#line 610 "chgeqz.f"
			i__4 = jch + 1 + jch * h_dim1;
#line 610 "chgeqz.f"
			ctemp.r = h__[i__4].r, ctemp.i = h__[i__4].i;
#line 611 "chgeqz.f"
			clartg_(&ctemp, &h__[jch + 1 + (jch - 1) * h_dim1], &
				c__, &s, &h__[jch + 1 + jch * h_dim1]);
#line 613 "chgeqz.f"
			i__4 = jch + 1 + (jch - 1) * h_dim1;
#line 613 "chgeqz.f"
			h__[i__4].r = 0., h__[i__4].i = 0.;
#line 614 "chgeqz.f"
			i__4 = jch + 1 - ifrstm;
#line 614 "chgeqz.f"
			crot_(&i__4, &h__[ifrstm + jch * h_dim1], &c__1, &h__[
				ifrstm + (jch - 1) * h_dim1], &c__1, &c__, &s)
				;
#line 616 "chgeqz.f"
			i__4 = jch - ifrstm;
#line 616 "chgeqz.f"
			crot_(&i__4, &t[ifrstm + jch * t_dim1], &c__1, &t[
				ifrstm + (jch - 1) * t_dim1], &c__1, &c__, &s)
				;
#line 618 "chgeqz.f"
			if (ilz) {
#line 618 "chgeqz.f"
			    crot_(n, &z__[jch * z_dim1 + 1], &c__1, &z__[(jch 
				    - 1) * z_dim1 + 1], &c__1, &c__, &s);
#line 618 "chgeqz.f"
			}
#line 621 "chgeqz.f"
/* L30: */
#line 621 "chgeqz.f"
		    }
#line 622 "chgeqz.f"
		    goto L50;
#line 623 "chgeqz.f"
		}
#line 624 "chgeqz.f"
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

#line 628 "chgeqz.f"
		ifirst = j;
#line 629 "chgeqz.f"
		goto L70;
#line 630 "chgeqz.f"
	    }

/*           Neither test passed -- try next J */

#line 634 "chgeqz.f"
/* L40: */
#line 634 "chgeqz.f"
	}

/*        (Drop-through is "impossible") */

#line 638 "chgeqz.f"
	*info = (*n << 1) + 1;
#line 639 "chgeqz.f"
	goto L210;

/*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a */
/*        1x1 block. */

#line 644 "chgeqz.f"
L50:
#line 645 "chgeqz.f"
	i__2 = ilast + ilast * h_dim1;
#line 645 "chgeqz.f"
	ctemp.r = h__[i__2].r, ctemp.i = h__[i__2].i;
#line 646 "chgeqz.f"
	clartg_(&ctemp, &h__[ilast + (ilast - 1) * h_dim1], &c__, &s, &h__[
		ilast + ilast * h_dim1]);
#line 648 "chgeqz.f"
	i__2 = ilast + (ilast - 1) * h_dim1;
#line 648 "chgeqz.f"
	h__[i__2].r = 0., h__[i__2].i = 0.;
#line 649 "chgeqz.f"
	i__2 = ilast - ifrstm;
#line 649 "chgeqz.f"
	crot_(&i__2, &h__[ifrstm + ilast * h_dim1], &c__1, &h__[ifrstm + (
		ilast - 1) * h_dim1], &c__1, &c__, &s);
#line 651 "chgeqz.f"
	i__2 = ilast - ifrstm;
#line 651 "chgeqz.f"
	crot_(&i__2, &t[ifrstm + ilast * t_dim1], &c__1, &t[ifrstm + (ilast - 
		1) * t_dim1], &c__1, &c__, &s);
#line 653 "chgeqz.f"
	if (ilz) {
#line 653 "chgeqz.f"
	    crot_(n, &z__[ilast * z_dim1 + 1], &c__1, &z__[(ilast - 1) * 
		    z_dim1 + 1], &c__1, &c__, &s);
#line 653 "chgeqz.f"
	}

/*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA */

#line 658 "chgeqz.f"
L60:
#line 659 "chgeqz.f"
	absb = z_abs(&t[ilast + ilast * t_dim1]);
#line 660 "chgeqz.f"
	if (absb > safmin) {
#line 661 "chgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 661 "chgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 661 "chgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 661 "chgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 662 "chgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 662 "chgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 663 "chgeqz.f"
	    if (ilschr) {
#line 664 "chgeqz.f"
		i__2 = ilast - ifrstm;
#line 664 "chgeqz.f"
		cscal_(&i__2, &signbc, &t[ifrstm + ilast * t_dim1], &c__1);
#line 665 "chgeqz.f"
		i__2 = ilast + 1 - ifrstm;
#line 665 "chgeqz.f"
		cscal_(&i__2, &signbc, &h__[ifrstm + ilast * h_dim1], &c__1);
#line 667 "chgeqz.f"
	    } else {
#line 668 "chgeqz.f"
		i__2 = ilast + ilast * h_dim1;
#line 668 "chgeqz.f"
		i__3 = ilast + ilast * h_dim1;
#line 668 "chgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 668 "chgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 669 "chgeqz.f"
	    }
#line 670 "chgeqz.f"
	    if (ilz) {
#line 670 "chgeqz.f"
		cscal_(n, &signbc, &z__[ilast * z_dim1 + 1], &c__1);
#line 670 "chgeqz.f"
	    }
#line 672 "chgeqz.f"
	} else {
#line 673 "chgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 673 "chgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 674 "chgeqz.f"
	}
#line 675 "chgeqz.f"
	i__2 = ilast;
#line 675 "chgeqz.f"
	i__3 = ilast + ilast * h_dim1;
#line 675 "chgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 676 "chgeqz.f"
	i__2 = ilast;
#line 676 "chgeqz.f"
	i__3 = ilast + ilast * t_dim1;
#line 676 "chgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;

/*        Go to next block -- exit if finished. */

#line 680 "chgeqz.f"
	--ilast;
#line 681 "chgeqz.f"
	if (ilast < *ilo) {
#line 681 "chgeqz.f"
	    goto L190;
#line 681 "chgeqz.f"
	}

/*        Reset counters */

#line 686 "chgeqz.f"
	iiter = 0;
#line 687 "chgeqz.f"
	eshift.r = 0., eshift.i = 0.;
#line 688 "chgeqz.f"
	if (! ilschr) {
#line 689 "chgeqz.f"
	    ilastm = ilast;
#line 690 "chgeqz.f"
	    if (ifrstm > ilast) {
#line 690 "chgeqz.f"
		ifrstm = *ilo;
#line 690 "chgeqz.f"
	    }
#line 692 "chgeqz.f"
	}
#line 693 "chgeqz.f"
	goto L160;

/*        QZ step */

/*        This iteration only involves rows/columns IFIRST:ILAST.  We */
/*        assume IFIRST < ILAST, and that the diagonal of B is non-zero. */

#line 700 "chgeqz.f"
L70:
#line 701 "chgeqz.f"
	++iiter;
#line 702 "chgeqz.f"
	if (! ilschr) {
#line 703 "chgeqz.f"
	    ifrstm = ifirst;
#line 704 "chgeqz.f"
	}

/*        Compute the Shift. */

/*        At this point, IFIRST < ILAST, and the diagonal elements of */
/*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in */
/*        magnitude) */

#line 712 "chgeqz.f"
	if (iiter / 10 * 10 != iiter) {

/*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of */
/*           the bottom-right 2x2 block of A inv(B) which is nearest to */
/*           the bottom-right element. */

/*           We factor B as U*D, where U has unit diagonals, and */
/*           compute (A*inv(D))*inv(U). */

#line 721 "chgeqz.f"
	    i__2 = ilast - 1 + ilast * t_dim1;
#line 721 "chgeqz.f"
	    z__2.r = bscale * t[i__2].r, z__2.i = bscale * t[i__2].i;
#line 721 "chgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 721 "chgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 721 "chgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 721 "chgeqz.f"
	    u12.r = z__1.r, u12.i = z__1.i;
#line 723 "chgeqz.f"
	    i__2 = ilast - 1 + (ilast - 1) * h_dim1;
#line 723 "chgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 723 "chgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 723 "chgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 723 "chgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 723 "chgeqz.f"
	    ad11.r = z__1.r, ad11.i = z__1.i;
#line 725 "chgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 725 "chgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 725 "chgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 725 "chgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 725 "chgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 725 "chgeqz.f"
	    ad21.r = z__1.r, ad21.i = z__1.i;
#line 727 "chgeqz.f"
	    i__2 = ilast - 1 + ilast * h_dim1;
#line 727 "chgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 727 "chgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 727 "chgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 727 "chgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 727 "chgeqz.f"
	    ad12.r = z__1.r, ad12.i = z__1.i;
#line 729 "chgeqz.f"
	    i__2 = ilast + ilast * h_dim1;
#line 729 "chgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 729 "chgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 729 "chgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 729 "chgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 729 "chgeqz.f"
	    ad22.r = z__1.r, ad22.i = z__1.i;
#line 731 "chgeqz.f"
	    z__2.r = u12.r * ad21.r - u12.i * ad21.i, z__2.i = u12.r * ad21.i 
		    + u12.i * ad21.r;
#line 731 "chgeqz.f"
	    z__1.r = ad22.r - z__2.r, z__1.i = ad22.i - z__2.i;
#line 731 "chgeqz.f"
	    abi22.r = z__1.r, abi22.i = z__1.i;

#line 733 "chgeqz.f"
	    z__2.r = ad11.r + abi22.r, z__2.i = ad11.i + abi22.i;
#line 733 "chgeqz.f"
	    z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
#line 733 "chgeqz.f"
	    t1.r = z__1.r, t1.i = z__1.i;
#line 734 "chgeqz.f"
	    pow_zi(&z__4, &t1, &c__2);
#line 734 "chgeqz.f"
	    z__5.r = ad12.r * ad21.r - ad12.i * ad21.i, z__5.i = ad12.r * 
		    ad21.i + ad12.i * ad21.r;
#line 734 "chgeqz.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 734 "chgeqz.f"
	    z__6.r = ad11.r * ad22.r - ad11.i * ad22.i, z__6.i = ad11.r * 
		    ad22.i + ad11.i * ad22.r;
#line 734 "chgeqz.f"
	    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 734 "chgeqz.f"
	    z_sqrt(&z__1, &z__2);
#line 734 "chgeqz.f"
	    rtdisc.r = z__1.r, rtdisc.i = z__1.i;
#line 735 "chgeqz.f"
	    z__1.r = t1.r - abi22.r, z__1.i = t1.i - abi22.i;
#line 735 "chgeqz.f"
	    z__2.r = t1.r - abi22.r, z__2.i = t1.i - abi22.i;
#line 735 "chgeqz.f"
	    temp = z__1.r * rtdisc.r + d_imag(&z__2) * d_imag(&rtdisc);
#line 737 "chgeqz.f"
	    if (temp <= 0.) {
#line 738 "chgeqz.f"
		z__1.r = t1.r + rtdisc.r, z__1.i = t1.i + rtdisc.i;
#line 738 "chgeqz.f"
		shift.r = z__1.r, shift.i = z__1.i;
#line 739 "chgeqz.f"
	    } else {
#line 740 "chgeqz.f"
		z__1.r = t1.r - rtdisc.r, z__1.i = t1.i - rtdisc.i;
#line 740 "chgeqz.f"
		shift.r = z__1.r, shift.i = z__1.i;
#line 741 "chgeqz.f"
	    }
#line 742 "chgeqz.f"
	} else {

/*           Exceptional shift.  Chosen for no particularly good reason. */

#line 746 "chgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 746 "chgeqz.f"
	    z__3.r = ascale * h__[i__2].r, z__3.i = ascale * h__[i__2].i;
#line 746 "chgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 746 "chgeqz.f"
	    z__4.r = bscale * t[i__3].r, z__4.i = bscale * t[i__3].i;
#line 746 "chgeqz.f"
	    z_div(&z__2, &z__3, &z__4);
#line 746 "chgeqz.f"
	    z__1.r = eshift.r + z__2.r, z__1.i = eshift.i + z__2.i;
#line 746 "chgeqz.f"
	    eshift.r = z__1.r, eshift.i = z__1.i;
#line 748 "chgeqz.f"
	    shift.r = eshift.r, shift.i = eshift.i;
#line 749 "chgeqz.f"
	}

/*        Now check for two consecutive small subdiagonals. */

#line 753 "chgeqz.f"
	i__2 = ifirst + 1;
#line 753 "chgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {
#line 754 "chgeqz.f"
	    istart = j;
#line 755 "chgeqz.f"
	    i__3 = j + j * h_dim1;
#line 755 "chgeqz.f"
	    z__2.r = ascale * h__[i__3].r, z__2.i = ascale * h__[i__3].i;
#line 755 "chgeqz.f"
	    i__4 = j + j * t_dim1;
#line 755 "chgeqz.f"
	    z__4.r = bscale * t[i__4].r, z__4.i = bscale * t[i__4].i;
#line 755 "chgeqz.f"
	    z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		    z__4.i + shift.i * z__4.r;
#line 755 "chgeqz.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 755 "chgeqz.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 756 "chgeqz.f"
	    temp = (d__1 = ctemp.r, abs(d__1)) + (d__2 = d_imag(&ctemp), abs(
		    d__2));
#line 757 "chgeqz.f"
	    i__3 = j + 1 + j * h_dim1;
#line 757 "chgeqz.f"
	    temp2 = ascale * ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&h__[j + 1 + j * h_dim1]), abs(d__2)));
#line 758 "chgeqz.f"
	    tempr = max(temp,temp2);
#line 759 "chgeqz.f"
	    if (tempr < 1. && tempr != 0.) {
#line 760 "chgeqz.f"
		temp /= tempr;
#line 761 "chgeqz.f"
		temp2 /= tempr;
#line 762 "chgeqz.f"
	    }
#line 763 "chgeqz.f"
	    i__3 = j + (j - 1) * h_dim1;
#line 763 "chgeqz.f"
	    if (((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[j + (j 
		    - 1) * h_dim1]), abs(d__2))) * temp2 <= temp * atol) {
#line 763 "chgeqz.f"
		goto L90;
#line 763 "chgeqz.f"
	    }
#line 765 "chgeqz.f"
/* L80: */
#line 765 "chgeqz.f"
	}

#line 767 "chgeqz.f"
	istart = ifirst;
#line 768 "chgeqz.f"
	i__2 = ifirst + ifirst * h_dim1;
#line 768 "chgeqz.f"
	z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 768 "chgeqz.f"
	i__3 = ifirst + ifirst * t_dim1;
#line 768 "chgeqz.f"
	z__4.r = bscale * t[i__3].r, z__4.i = bscale * t[i__3].i;
#line 768 "chgeqz.f"
	z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		z__4.i + shift.i * z__4.r;
#line 768 "chgeqz.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 768 "chgeqz.f"
	ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 770 "chgeqz.f"
L90:

/*        Do an implicit-shift QZ sweep. */

/*        Initial Q */

#line 776 "chgeqz.f"
	i__2 = istart + 1 + istart * h_dim1;
#line 776 "chgeqz.f"
	z__1.r = ascale * h__[i__2].r, z__1.i = ascale * h__[i__2].i;
#line 776 "chgeqz.f"
	ctemp2.r = z__1.r, ctemp2.i = z__1.i;
#line 777 "chgeqz.f"
	clartg_(&ctemp, &ctemp2, &c__, &s, &ctemp3);

/*        Sweep */

#line 781 "chgeqz.f"
	i__2 = ilast - 1;
#line 781 "chgeqz.f"
	for (j = istart; j <= i__2; ++j) {
#line 782 "chgeqz.f"
	    if (j > istart) {
#line 783 "chgeqz.f"
		i__3 = j + (j - 1) * h_dim1;
#line 783 "chgeqz.f"
		ctemp.r = h__[i__3].r, ctemp.i = h__[i__3].i;
#line 784 "chgeqz.f"
		clartg_(&ctemp, &h__[j + 1 + (j - 1) * h_dim1], &c__, &s, &
			h__[j + (j - 1) * h_dim1]);
#line 785 "chgeqz.f"
		i__3 = j + 1 + (j - 1) * h_dim1;
#line 785 "chgeqz.f"
		h__[i__3].r = 0., h__[i__3].i = 0.;
#line 786 "chgeqz.f"
	    }

#line 788 "chgeqz.f"
	    i__3 = ilastm;
#line 788 "chgeqz.f"
	    for (jc = j; jc <= i__3; ++jc) {
#line 789 "chgeqz.f"
		i__4 = j + jc * h_dim1;
#line 789 "chgeqz.f"
		z__2.r = c__ * h__[i__4].r, z__2.i = c__ * h__[i__4].i;
#line 789 "chgeqz.f"
		i__5 = j + 1 + jc * h_dim1;
#line 789 "chgeqz.f"
		z__3.r = s.r * h__[i__5].r - s.i * h__[i__5].i, z__3.i = s.r *
			 h__[i__5].i + s.i * h__[i__5].r;
#line 789 "chgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 789 "chgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 790 "chgeqz.f"
		i__4 = j + 1 + jc * h_dim1;
#line 790 "chgeqz.f"
		d_cnjg(&z__4, &s);
#line 790 "chgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 790 "chgeqz.f"
		i__5 = j + jc * h_dim1;
#line 790 "chgeqz.f"
		z__2.r = z__3.r * h__[i__5].r - z__3.i * h__[i__5].i, z__2.i =
			 z__3.r * h__[i__5].i + z__3.i * h__[i__5].r;
#line 790 "chgeqz.f"
		i__6 = j + 1 + jc * h_dim1;
#line 790 "chgeqz.f"
		z__5.r = c__ * h__[i__6].r, z__5.i = c__ * h__[i__6].i;
#line 790 "chgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 790 "chgeqz.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 791 "chgeqz.f"
		i__4 = j + jc * h_dim1;
#line 791 "chgeqz.f"
		h__[i__4].r = ctemp.r, h__[i__4].i = ctemp.i;
#line 792 "chgeqz.f"
		i__4 = j + jc * t_dim1;
#line 792 "chgeqz.f"
		z__2.r = c__ * t[i__4].r, z__2.i = c__ * t[i__4].i;
#line 792 "chgeqz.f"
		i__5 = j + 1 + jc * t_dim1;
#line 792 "chgeqz.f"
		z__3.r = s.r * t[i__5].r - s.i * t[i__5].i, z__3.i = s.r * t[
			i__5].i + s.i * t[i__5].r;
#line 792 "chgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 792 "chgeqz.f"
		ctemp2.r = z__1.r, ctemp2.i = z__1.i;
#line 793 "chgeqz.f"
		i__4 = j + 1 + jc * t_dim1;
#line 793 "chgeqz.f"
		d_cnjg(&z__4, &s);
#line 793 "chgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 793 "chgeqz.f"
		i__5 = j + jc * t_dim1;
#line 793 "chgeqz.f"
		z__2.r = z__3.r * t[i__5].r - z__3.i * t[i__5].i, z__2.i = 
			z__3.r * t[i__5].i + z__3.i * t[i__5].r;
#line 793 "chgeqz.f"
		i__6 = j + 1 + jc * t_dim1;
#line 793 "chgeqz.f"
		z__5.r = c__ * t[i__6].r, z__5.i = c__ * t[i__6].i;
#line 793 "chgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 793 "chgeqz.f"
		t[i__4].r = z__1.r, t[i__4].i = z__1.i;
#line 794 "chgeqz.f"
		i__4 = j + jc * t_dim1;
#line 794 "chgeqz.f"
		t[i__4].r = ctemp2.r, t[i__4].i = ctemp2.i;
#line 795 "chgeqz.f"
/* L100: */
#line 795 "chgeqz.f"
	    }
#line 796 "chgeqz.f"
	    if (ilq) {
#line 797 "chgeqz.f"
		i__3 = *n;
#line 797 "chgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 798 "chgeqz.f"
		    i__4 = jr + j * q_dim1;
#line 798 "chgeqz.f"
		    z__2.r = c__ * q[i__4].r, z__2.i = c__ * q[i__4].i;
#line 798 "chgeqz.f"
		    d_cnjg(&z__4, &s);
#line 798 "chgeqz.f"
		    i__5 = jr + (j + 1) * q_dim1;
#line 798 "chgeqz.f"
		    z__3.r = z__4.r * q[i__5].r - z__4.i * q[i__5].i, z__3.i =
			     z__4.r * q[i__5].i + z__4.i * q[i__5].r;
#line 798 "chgeqz.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 798 "chgeqz.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 799 "chgeqz.f"
		    i__4 = jr + (j + 1) * q_dim1;
#line 799 "chgeqz.f"
		    z__3.r = -s.r, z__3.i = -s.i;
#line 799 "chgeqz.f"
		    i__5 = jr + j * q_dim1;
#line 799 "chgeqz.f"
		    z__2.r = z__3.r * q[i__5].r - z__3.i * q[i__5].i, z__2.i =
			     z__3.r * q[i__5].i + z__3.i * q[i__5].r;
#line 799 "chgeqz.f"
		    i__6 = jr + (j + 1) * q_dim1;
#line 799 "chgeqz.f"
		    z__4.r = c__ * q[i__6].r, z__4.i = c__ * q[i__6].i;
#line 799 "chgeqz.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 799 "chgeqz.f"
		    q[i__4].r = z__1.r, q[i__4].i = z__1.i;
#line 800 "chgeqz.f"
		    i__4 = jr + j * q_dim1;
#line 800 "chgeqz.f"
		    q[i__4].r = ctemp.r, q[i__4].i = ctemp.i;
#line 801 "chgeqz.f"
/* L110: */
#line 801 "chgeqz.f"
		}
#line 802 "chgeqz.f"
	    }

#line 804 "chgeqz.f"
	    i__3 = j + 1 + (j + 1) * t_dim1;
#line 804 "chgeqz.f"
	    ctemp.r = t[i__3].r, ctemp.i = t[i__3].i;
#line 805 "chgeqz.f"
	    clartg_(&ctemp, &t[j + 1 + j * t_dim1], &c__, &s, &t[j + 1 + (j + 
		    1) * t_dim1]);
#line 806 "chgeqz.f"
	    i__3 = j + 1 + j * t_dim1;
#line 806 "chgeqz.f"
	    t[i__3].r = 0., t[i__3].i = 0.;

/* Computing MIN */
#line 808 "chgeqz.f"
	    i__4 = j + 2;
#line 808 "chgeqz.f"
	    i__3 = min(i__4,ilast);
#line 808 "chgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 809 "chgeqz.f"
		i__4 = jr + (j + 1) * h_dim1;
#line 809 "chgeqz.f"
		z__2.r = c__ * h__[i__4].r, z__2.i = c__ * h__[i__4].i;
#line 809 "chgeqz.f"
		i__5 = jr + j * h_dim1;
#line 809 "chgeqz.f"
		z__3.r = s.r * h__[i__5].r - s.i * h__[i__5].i, z__3.i = s.r *
			 h__[i__5].i + s.i * h__[i__5].r;
#line 809 "chgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 809 "chgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 810 "chgeqz.f"
		i__4 = jr + j * h_dim1;
#line 810 "chgeqz.f"
		d_cnjg(&z__4, &s);
#line 810 "chgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 810 "chgeqz.f"
		i__5 = jr + (j + 1) * h_dim1;
#line 810 "chgeqz.f"
		z__2.r = z__3.r * h__[i__5].r - z__3.i * h__[i__5].i, z__2.i =
			 z__3.r * h__[i__5].i + z__3.i * h__[i__5].r;
#line 810 "chgeqz.f"
		i__6 = jr + j * h_dim1;
#line 810 "chgeqz.f"
		z__5.r = c__ * h__[i__6].r, z__5.i = c__ * h__[i__6].i;
#line 810 "chgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 810 "chgeqz.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 811 "chgeqz.f"
		i__4 = jr + (j + 1) * h_dim1;
#line 811 "chgeqz.f"
		h__[i__4].r = ctemp.r, h__[i__4].i = ctemp.i;
#line 812 "chgeqz.f"
/* L120: */
#line 812 "chgeqz.f"
	    }
#line 813 "chgeqz.f"
	    i__3 = j;
#line 813 "chgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 814 "chgeqz.f"
		i__4 = jr + (j + 1) * t_dim1;
#line 814 "chgeqz.f"
		z__2.r = c__ * t[i__4].r, z__2.i = c__ * t[i__4].i;
#line 814 "chgeqz.f"
		i__5 = jr + j * t_dim1;
#line 814 "chgeqz.f"
		z__3.r = s.r * t[i__5].r - s.i * t[i__5].i, z__3.i = s.r * t[
			i__5].i + s.i * t[i__5].r;
#line 814 "chgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 814 "chgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 815 "chgeqz.f"
		i__4 = jr + j * t_dim1;
#line 815 "chgeqz.f"
		d_cnjg(&z__4, &s);
#line 815 "chgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 815 "chgeqz.f"
		i__5 = jr + (j + 1) * t_dim1;
#line 815 "chgeqz.f"
		z__2.r = z__3.r * t[i__5].r - z__3.i * t[i__5].i, z__2.i = 
			z__3.r * t[i__5].i + z__3.i * t[i__5].r;
#line 815 "chgeqz.f"
		i__6 = jr + j * t_dim1;
#line 815 "chgeqz.f"
		z__5.r = c__ * t[i__6].r, z__5.i = c__ * t[i__6].i;
#line 815 "chgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 815 "chgeqz.f"
		t[i__4].r = z__1.r, t[i__4].i = z__1.i;
#line 816 "chgeqz.f"
		i__4 = jr + (j + 1) * t_dim1;
#line 816 "chgeqz.f"
		t[i__4].r = ctemp.r, t[i__4].i = ctemp.i;
#line 817 "chgeqz.f"
/* L130: */
#line 817 "chgeqz.f"
	    }
#line 818 "chgeqz.f"
	    if (ilz) {
#line 819 "chgeqz.f"
		i__3 = *n;
#line 819 "chgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 820 "chgeqz.f"
		    i__4 = jr + (j + 1) * z_dim1;
#line 820 "chgeqz.f"
		    z__2.r = c__ * z__[i__4].r, z__2.i = c__ * z__[i__4].i;
#line 820 "chgeqz.f"
		    i__5 = jr + j * z_dim1;
#line 820 "chgeqz.f"
		    z__3.r = s.r * z__[i__5].r - s.i * z__[i__5].i, z__3.i = 
			    s.r * z__[i__5].i + s.i * z__[i__5].r;
#line 820 "chgeqz.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 820 "chgeqz.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 821 "chgeqz.f"
		    i__4 = jr + j * z_dim1;
#line 821 "chgeqz.f"
		    d_cnjg(&z__4, &s);
#line 821 "chgeqz.f"
		    z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 821 "chgeqz.f"
		    i__5 = jr + (j + 1) * z_dim1;
#line 821 "chgeqz.f"
		    z__2.r = z__3.r * z__[i__5].r - z__3.i * z__[i__5].i, 
			    z__2.i = z__3.r * z__[i__5].i + z__3.i * z__[i__5]
			    .r;
#line 821 "chgeqz.f"
		    i__6 = jr + j * z_dim1;
#line 821 "chgeqz.f"
		    z__5.r = c__ * z__[i__6].r, z__5.i = c__ * z__[i__6].i;
#line 821 "chgeqz.f"
		    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 821 "chgeqz.f"
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 822 "chgeqz.f"
		    i__4 = jr + (j + 1) * z_dim1;
#line 822 "chgeqz.f"
		    z__[i__4].r = ctemp.r, z__[i__4].i = ctemp.i;
#line 823 "chgeqz.f"
/* L140: */
#line 823 "chgeqz.f"
		}
#line 824 "chgeqz.f"
	    }
#line 825 "chgeqz.f"
/* L150: */
#line 825 "chgeqz.f"
	}

#line 827 "chgeqz.f"
L160:

#line 829 "chgeqz.f"
/* L170: */
#line 829 "chgeqz.f"
	;
#line 829 "chgeqz.f"
    }

/*     Drop-through = non-convergence */

#line 833 "chgeqz.f"
L180:
#line 834 "chgeqz.f"
    *info = ilast;
#line 835 "chgeqz.f"
    goto L210;

/*     Successful completion of all QZ steps */

#line 839 "chgeqz.f"
L190:

/*     Set Eigenvalues 1:ILO-1 */

#line 843 "chgeqz.f"
    i__1 = *ilo - 1;
#line 843 "chgeqz.f"
    for (j = 1; j <= i__1; ++j) {
#line 844 "chgeqz.f"
	absb = z_abs(&t[j + j * t_dim1]);
#line 845 "chgeqz.f"
	if (absb > safmin) {
#line 846 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 846 "chgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 846 "chgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 846 "chgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 847 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 847 "chgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 848 "chgeqz.f"
	    if (ilschr) {
#line 849 "chgeqz.f"
		i__2 = j - 1;
#line 849 "chgeqz.f"
		cscal_(&i__2, &signbc, &t[j * t_dim1 + 1], &c__1);
#line 850 "chgeqz.f"
		cscal_(&j, &signbc, &h__[j * h_dim1 + 1], &c__1);
#line 851 "chgeqz.f"
	    } else {
#line 852 "chgeqz.f"
		i__2 = j + j * h_dim1;
#line 852 "chgeqz.f"
		i__3 = j + j * h_dim1;
#line 852 "chgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 852 "chgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 853 "chgeqz.f"
	    }
#line 854 "chgeqz.f"
	    if (ilz) {
#line 854 "chgeqz.f"
		cscal_(n, &signbc, &z__[j * z_dim1 + 1], &c__1);
#line 854 "chgeqz.f"
	    }
#line 856 "chgeqz.f"
	} else {
#line 857 "chgeqz.f"
	    i__2 = j + j * t_dim1;
#line 857 "chgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 858 "chgeqz.f"
	}
#line 859 "chgeqz.f"
	i__2 = j;
#line 859 "chgeqz.f"
	i__3 = j + j * h_dim1;
#line 859 "chgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 860 "chgeqz.f"
	i__2 = j;
#line 860 "chgeqz.f"
	i__3 = j + j * t_dim1;
#line 860 "chgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;
#line 861 "chgeqz.f"
/* L200: */
#line 861 "chgeqz.f"
    }

/*     Normal Termination */

#line 865 "chgeqz.f"
    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size */

#line 869 "chgeqz.f"
L210:
#line 870 "chgeqz.f"
    z__1.r = (doublereal) (*n), z__1.i = 0.;
#line 870 "chgeqz.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 871 "chgeqz.f"
    return 0;

/*     End of CHGEQZ */

} /* chgeqz_ */

