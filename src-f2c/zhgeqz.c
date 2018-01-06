#line 1 "zhgeqz.f"
/* zhgeqz.f -- translated by f2c (version 20100827).
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

#line 1 "zhgeqz.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b ZHGEQZ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHGEQZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhgeqz.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhgeqz.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhgeqz.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT, */
/*                          ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ, JOB */
/*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ), */
/*      $                   Q( LDQ, * ), T( LDT, * ), WORK( * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T), */
/* > where H is an upper Hessenberg matrix and T is upper triangular, */
/* > using the single-shift QZ method. */
/* > Matrix pairs of this type are produced by the reduction to */
/* > generalized upper Hessenberg form of a complex matrix pair (A,B): */
/* > */
/* >    A = Q1*H*Z1**H,  B = Q1*T*Z1**H, */
/* > */
/* > as computed by ZGGHRD. */
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
/* > If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced */
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
/* >          H is COMPLEX*16 array, dimension (LDH, N) */
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
/* >          T is COMPLEX*16 array, dimension (LDT, N) */
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
/* >          ALPHA is COMPLEX*16 array, dimension (N) */
/* >          The complex scalars alpha that define the eigenvalues of */
/* >          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur */
/* >          factorization. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
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
/* >          Q is COMPLEX*16 array, dimension (LDQ, N) */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup complex16GEcomputational */

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
/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, integer *n, 
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
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex ctemp;
    static integer iiter, ilast, jiter;
    static doublereal anorm, bnorm;
    static integer maxit;
    static doublecomplex shift;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal tempr;
    static doublecomplex ctemp2, ctemp3;
    static logical ilazr2;
    static doublereal ascale, bscale;
    extern doublereal dlamch_(char *, ftnlen);
    static doublecomplex signbc;
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublecomplex eshift;
    static logical ilschr;
    static integer icompq, ilastm;
    static doublecomplex rtdisc;
    static integer ischur;
    extern doublereal zlanhs_(char *, integer *, doublecomplex *, integer *, 
	    doublereal *, ftnlen);
    static logical ilazro;
    static integer icompz, ifirst;
    extern /* Subroutine */ int zlartg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *);
    static integer ifrstm;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static integer istart;
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

#line 347 "zhgeqz.f"
    /* Parameter adjustments */
#line 347 "zhgeqz.f"
    h_dim1 = *ldh;
#line 347 "zhgeqz.f"
    h_offset = 1 + h_dim1;
#line 347 "zhgeqz.f"
    h__ -= h_offset;
#line 347 "zhgeqz.f"
    t_dim1 = *ldt;
#line 347 "zhgeqz.f"
    t_offset = 1 + t_dim1;
#line 347 "zhgeqz.f"
    t -= t_offset;
#line 347 "zhgeqz.f"
    --alpha;
#line 347 "zhgeqz.f"
    --beta;
#line 347 "zhgeqz.f"
    q_dim1 = *ldq;
#line 347 "zhgeqz.f"
    q_offset = 1 + q_dim1;
#line 347 "zhgeqz.f"
    q -= q_offset;
#line 347 "zhgeqz.f"
    z_dim1 = *ldz;
#line 347 "zhgeqz.f"
    z_offset = 1 + z_dim1;
#line 347 "zhgeqz.f"
    z__ -= z_offset;
#line 347 "zhgeqz.f"
    --work;
#line 347 "zhgeqz.f"
    --rwork;
#line 347 "zhgeqz.f"

#line 347 "zhgeqz.f"
    /* Function Body */
#line 347 "zhgeqz.f"
    if (lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
#line 348 "zhgeqz.f"
	ilschr = FALSE_;
#line 349 "zhgeqz.f"
	ischur = 1;
#line 350 "zhgeqz.f"
    } else if (lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {
#line 351 "zhgeqz.f"
	ilschr = TRUE_;
#line 352 "zhgeqz.f"
	ischur = 2;
#line 353 "zhgeqz.f"
    } else {
#line 354 "zhgeqz.f"
	ischur = 0;
#line 355 "zhgeqz.f"
    }

#line 357 "zhgeqz.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 358 "zhgeqz.f"
	ilq = FALSE_;
#line 359 "zhgeqz.f"
	icompq = 1;
#line 360 "zhgeqz.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 361 "zhgeqz.f"
	ilq = TRUE_;
#line 362 "zhgeqz.f"
	icompq = 2;
#line 363 "zhgeqz.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 364 "zhgeqz.f"
	ilq = TRUE_;
#line 365 "zhgeqz.f"
	icompq = 3;
#line 366 "zhgeqz.f"
    } else {
#line 367 "zhgeqz.f"
	icompq = 0;
#line 368 "zhgeqz.f"
    }

#line 370 "zhgeqz.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 371 "zhgeqz.f"
	ilz = FALSE_;
#line 372 "zhgeqz.f"
	icompz = 1;
#line 373 "zhgeqz.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 374 "zhgeqz.f"
	ilz = TRUE_;
#line 375 "zhgeqz.f"
	icompz = 2;
#line 376 "zhgeqz.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 377 "zhgeqz.f"
	ilz = TRUE_;
#line 378 "zhgeqz.f"
	icompz = 3;
#line 379 "zhgeqz.f"
    } else {
#line 380 "zhgeqz.f"
	icompz = 0;
#line 381 "zhgeqz.f"
    }

/*     Check Argument Values */

#line 385 "zhgeqz.f"
    *info = 0;
#line 386 "zhgeqz.f"
    i__1 = max(1,*n);
#line 386 "zhgeqz.f"
    work[1].r = (doublereal) i__1, work[1].i = 0.;
#line 387 "zhgeqz.f"
    lquery = *lwork == -1;
#line 388 "zhgeqz.f"
    if (ischur == 0) {
#line 389 "zhgeqz.f"
	*info = -1;
#line 390 "zhgeqz.f"
    } else if (icompq == 0) {
#line 391 "zhgeqz.f"
	*info = -2;
#line 392 "zhgeqz.f"
    } else if (icompz == 0) {
#line 393 "zhgeqz.f"
	*info = -3;
#line 394 "zhgeqz.f"
    } else if (*n < 0) {
#line 395 "zhgeqz.f"
	*info = -4;
#line 396 "zhgeqz.f"
    } else if (*ilo < 1) {
#line 397 "zhgeqz.f"
	*info = -5;
#line 398 "zhgeqz.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 399 "zhgeqz.f"
	*info = -6;
#line 400 "zhgeqz.f"
    } else if (*ldh < *n) {
#line 401 "zhgeqz.f"
	*info = -8;
#line 402 "zhgeqz.f"
    } else if (*ldt < *n) {
#line 403 "zhgeqz.f"
	*info = -10;
#line 404 "zhgeqz.f"
    } else if (*ldq < 1 || ilq && *ldq < *n) {
#line 405 "zhgeqz.f"
	*info = -14;
#line 406 "zhgeqz.f"
    } else if (*ldz < 1 || ilz && *ldz < *n) {
#line 407 "zhgeqz.f"
	*info = -16;
#line 408 "zhgeqz.f"
    } else if (*lwork < max(1,*n) && ! lquery) {
#line 409 "zhgeqz.f"
	*info = -18;
#line 410 "zhgeqz.f"
    }
#line 411 "zhgeqz.f"
    if (*info != 0) {
#line 412 "zhgeqz.f"
	i__1 = -(*info);
#line 412 "zhgeqz.f"
	xerbla_("ZHGEQZ", &i__1, (ftnlen)6);
#line 413 "zhgeqz.f"
	return 0;
#line 414 "zhgeqz.f"
    } else if (lquery) {
#line 415 "zhgeqz.f"
	return 0;
#line 416 "zhgeqz.f"
    }

/*     Quick return if possible */

/*     WORK( 1 ) = CMPLX( 1 ) */
#line 421 "zhgeqz.f"
    if (*n <= 0) {
#line 422 "zhgeqz.f"
	work[1].r = 1., work[1].i = 0.;
#line 423 "zhgeqz.f"
	return 0;
#line 424 "zhgeqz.f"
    }

/*     Initialize Q and Z */

#line 428 "zhgeqz.f"
    if (icompq == 3) {
#line 428 "zhgeqz.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 428 "zhgeqz.f"
    }
#line 430 "zhgeqz.f"
    if (icompz == 3) {
#line 430 "zhgeqz.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &z__[z_offset], ldz, (ftnlen)4);
#line 430 "zhgeqz.f"
    }

/*     Machine Constants */

#line 435 "zhgeqz.f"
    in = *ihi + 1 - *ilo;
#line 436 "zhgeqz.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 437 "zhgeqz.f"
    ulp = dlamch_("E", (ftnlen)1) * dlamch_("B", (ftnlen)1);
#line 438 "zhgeqz.f"
    anorm = zlanhs_("F", &in, &h__[*ilo + *ilo * h_dim1], ldh, &rwork[1], (
	    ftnlen)1);
#line 439 "zhgeqz.f"
    bnorm = zlanhs_("F", &in, &t[*ilo + *ilo * t_dim1], ldt, &rwork[1], (
	    ftnlen)1);
/* Computing MAX */
#line 440 "zhgeqz.f"
    d__1 = safmin, d__2 = ulp * anorm;
#line 440 "zhgeqz.f"
    atol = max(d__1,d__2);
/* Computing MAX */
#line 441 "zhgeqz.f"
    d__1 = safmin, d__2 = ulp * bnorm;
#line 441 "zhgeqz.f"
    btol = max(d__1,d__2);
#line 442 "zhgeqz.f"
    ascale = 1. / max(safmin,anorm);
#line 443 "zhgeqz.f"
    bscale = 1. / max(safmin,bnorm);


/*     Set Eigenvalues IHI+1:N */

#line 448 "zhgeqz.f"
    i__1 = *n;
#line 448 "zhgeqz.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 449 "zhgeqz.f"
	absb = z_abs(&t[j + j * t_dim1]);
#line 450 "zhgeqz.f"
	if (absb > safmin) {
#line 451 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 451 "zhgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 451 "zhgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 451 "zhgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 452 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 452 "zhgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 453 "zhgeqz.f"
	    if (ilschr) {
#line 454 "zhgeqz.f"
		i__2 = j - 1;
#line 454 "zhgeqz.f"
		zscal_(&i__2, &signbc, &t[j * t_dim1 + 1], &c__1);
#line 455 "zhgeqz.f"
		zscal_(&j, &signbc, &h__[j * h_dim1 + 1], &c__1);
#line 456 "zhgeqz.f"
	    } else {
#line 457 "zhgeqz.f"
		i__2 = j + j * h_dim1;
#line 457 "zhgeqz.f"
		i__3 = j + j * h_dim1;
#line 457 "zhgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 457 "zhgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 458 "zhgeqz.f"
	    }
#line 459 "zhgeqz.f"
	    if (ilz) {
#line 459 "zhgeqz.f"
		zscal_(n, &signbc, &z__[j * z_dim1 + 1], &c__1);
#line 459 "zhgeqz.f"
	    }
#line 461 "zhgeqz.f"
	} else {
#line 462 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 462 "zhgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 463 "zhgeqz.f"
	}
#line 464 "zhgeqz.f"
	i__2 = j;
#line 464 "zhgeqz.f"
	i__3 = j + j * h_dim1;
#line 464 "zhgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 465 "zhgeqz.f"
	i__2 = j;
#line 465 "zhgeqz.f"
	i__3 = j + j * t_dim1;
#line 465 "zhgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;
#line 466 "zhgeqz.f"
/* L10: */
#line 466 "zhgeqz.f"
    }

/*     If IHI < ILO, skip QZ steps */

#line 470 "zhgeqz.f"
    if (*ihi < *ilo) {
#line 470 "zhgeqz.f"
	goto L190;
#line 470 "zhgeqz.f"
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

#line 488 "zhgeqz.f"
    ilast = *ihi;
#line 489 "zhgeqz.f"
    if (ilschr) {
#line 490 "zhgeqz.f"
	ifrstm = 1;
#line 491 "zhgeqz.f"
	ilastm = *n;
#line 492 "zhgeqz.f"
    } else {
#line 493 "zhgeqz.f"
	ifrstm = *ilo;
#line 494 "zhgeqz.f"
	ilastm = *ihi;
#line 495 "zhgeqz.f"
    }
#line 496 "zhgeqz.f"
    iiter = 0;
#line 497 "zhgeqz.f"
    eshift.r = 0., eshift.i = 0.;
#line 498 "zhgeqz.f"
    maxit = (*ihi - *ilo + 1) * 30;

#line 500 "zhgeqz.f"
    i__1 = maxit;
#line 500 "zhgeqz.f"
    for (jiter = 1; jiter <= i__1; ++jiter) {

/*        Check for too many iterations. */

#line 504 "zhgeqz.f"
	if (jiter > maxit) {
#line 504 "zhgeqz.f"
	    goto L180;
#line 504 "zhgeqz.f"
	}

/*        Split the matrix if possible. */

/*        Two tests: */
/*           1: H(j,j-1)=0  or  j=ILO */
/*           2: T(j,j)=0 */

/*        Special case: j=ILAST */

#line 515 "zhgeqz.f"
	if (ilast == *ilo) {
#line 516 "zhgeqz.f"
	    goto L60;
#line 517 "zhgeqz.f"
	} else {
#line 518 "zhgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 518 "zhgeqz.f"
	    if ((d__1 = h__[i__2].r, abs(d__1)) + (d__2 = d_imag(&h__[ilast + 
		    (ilast - 1) * h_dim1]), abs(d__2)) <= atol) {
#line 519 "zhgeqz.f"
		i__2 = ilast + (ilast - 1) * h_dim1;
#line 519 "zhgeqz.f"
		h__[i__2].r = 0., h__[i__2].i = 0.;
#line 520 "zhgeqz.f"
		goto L60;
#line 521 "zhgeqz.f"
	    }
#line 522 "zhgeqz.f"
	}

#line 524 "zhgeqz.f"
	if (z_abs(&t[ilast + ilast * t_dim1]) <= btol) {
#line 525 "zhgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 525 "zhgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 526 "zhgeqz.f"
	    goto L50;
#line 527 "zhgeqz.f"
	}

/*        General case: j<ILAST */

#line 531 "zhgeqz.f"
	i__2 = *ilo;
#line 531 "zhgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {

/*           Test 1: for H(j,j-1)=0 or j=ILO */

#line 535 "zhgeqz.f"
	    if (j == *ilo) {
#line 536 "zhgeqz.f"
		ilazro = TRUE_;
#line 537 "zhgeqz.f"
	    } else {
#line 538 "zhgeqz.f"
		i__3 = j + (j - 1) * h_dim1;
#line 538 "zhgeqz.f"
		if ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[j + 
			(j - 1) * h_dim1]), abs(d__2)) <= atol) {
#line 539 "zhgeqz.f"
		    i__3 = j + (j - 1) * h_dim1;
#line 539 "zhgeqz.f"
		    h__[i__3].r = 0., h__[i__3].i = 0.;
#line 540 "zhgeqz.f"
		    ilazro = TRUE_;
#line 541 "zhgeqz.f"
		} else {
#line 542 "zhgeqz.f"
		    ilazro = FALSE_;
#line 543 "zhgeqz.f"
		}
#line 544 "zhgeqz.f"
	    }

/*           Test 2: for T(j,j)=0 */

#line 548 "zhgeqz.f"
	    if (z_abs(&t[j + j * t_dim1]) < btol) {
#line 549 "zhgeqz.f"
		i__3 = j + j * t_dim1;
#line 549 "zhgeqz.f"
		t[i__3].r = 0., t[i__3].i = 0.;

/*              Test 1a: Check for 2 consecutive small subdiagonals in A */

#line 553 "zhgeqz.f"
		ilazr2 = FALSE_;
#line 554 "zhgeqz.f"
		if (! ilazro) {
#line 555 "zhgeqz.f"
		    i__3 = j + (j - 1) * h_dim1;
#line 555 "zhgeqz.f"
		    i__4 = j + 1 + j * h_dim1;
#line 555 "zhgeqz.f"
		    i__5 = j + j * h_dim1;
#line 555 "zhgeqz.f"
		    if (((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
			    h__[j + (j - 1) * h_dim1]), abs(d__2))) * (ascale 
			    * ((d__3 = h__[i__4].r, abs(d__3)) + (d__4 = 
			    d_imag(&h__[j + 1 + j * h_dim1]), abs(d__4)))) <= 
			    ((d__5 = h__[i__5].r, abs(d__5)) + (d__6 = d_imag(
			    &h__[j + j * h_dim1]), abs(d__6))) * (ascale * 
			    atol)) {
#line 555 "zhgeqz.f"
			ilazr2 = TRUE_;
#line 555 "zhgeqz.f"
		    }
#line 558 "zhgeqz.f"
		}

/*              If both tests pass (1 & 2), i.e., the leading diagonal */
/*              element of B in the block is zero, split a 1x1 block off */
/*              at the top. (I.e., at the J-th row/column) The leading */
/*              diagonal element of the remainder can also be zero, so */
/*              this may have to be done repeatedly. */

#line 566 "zhgeqz.f"
		if (ilazro || ilazr2) {
#line 567 "zhgeqz.f"
		    i__3 = ilast - 1;
#line 567 "zhgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 568 "zhgeqz.f"
			i__4 = jch + jch * h_dim1;
#line 568 "zhgeqz.f"
			ctemp.r = h__[i__4].r, ctemp.i = h__[i__4].i;
#line 569 "zhgeqz.f"
			zlartg_(&ctemp, &h__[jch + 1 + jch * h_dim1], &c__, &
				s, &h__[jch + jch * h_dim1]);
#line 571 "zhgeqz.f"
			i__4 = jch + 1 + jch * h_dim1;
#line 571 "zhgeqz.f"
			h__[i__4].r = 0., h__[i__4].i = 0.;
#line 572 "zhgeqz.f"
			i__4 = ilastm - jch;
#line 572 "zhgeqz.f"
			zrot_(&i__4, &h__[jch + (jch + 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch + 1) * h_dim1], ldh, &c__, 
				&s);
#line 574 "zhgeqz.f"
			i__4 = ilastm - jch;
#line 574 "zhgeqz.f"
			zrot_(&i__4, &t[jch + (jch + 1) * t_dim1], ldt, &t[
				jch + 1 + (jch + 1) * t_dim1], ldt, &c__, &s);
#line 576 "zhgeqz.f"
			if (ilq) {
#line 576 "zhgeqz.f"
			    d_cnjg(&z__1, &s);
#line 576 "zhgeqz.f"
			    zrot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &z__1);
#line 576 "zhgeqz.f"
			}
#line 579 "zhgeqz.f"
			if (ilazr2) {
#line 579 "zhgeqz.f"
			    i__4 = jch + (jch - 1) * h_dim1;
#line 579 "zhgeqz.f"
			    i__5 = jch + (jch - 1) * h_dim1;
#line 579 "zhgeqz.f"
			    z__1.r = c__ * h__[i__5].r, z__1.i = c__ * h__[
				    i__5].i;
#line 579 "zhgeqz.f"
			    h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 579 "zhgeqz.f"
			}
#line 581 "zhgeqz.f"
			ilazr2 = FALSE_;
#line 582 "zhgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 582 "zhgeqz.f"
			if ((d__1 = t[i__4].r, abs(d__1)) + (d__2 = d_imag(&t[
				jch + 1 + (jch + 1) * t_dim1]), abs(d__2)) >= 
				btol) {
#line 583 "zhgeqz.f"
			    if (jch + 1 >= ilast) {
#line 584 "zhgeqz.f"
				goto L60;
#line 585 "zhgeqz.f"
			    } else {
#line 586 "zhgeqz.f"
				ifirst = jch + 1;
#line 587 "zhgeqz.f"
				goto L70;
#line 588 "zhgeqz.f"
			    }
#line 589 "zhgeqz.f"
			}
#line 590 "zhgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 590 "zhgeqz.f"
			t[i__4].r = 0., t[i__4].i = 0.;
#line 591 "zhgeqz.f"
/* L20: */
#line 591 "zhgeqz.f"
		    }
#line 592 "zhgeqz.f"
		    goto L50;
#line 593 "zhgeqz.f"
		} else {

/*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST) */
/*                 Then process as in the case T(ILAST,ILAST)=0 */

#line 598 "zhgeqz.f"
		    i__3 = ilast - 1;
#line 598 "zhgeqz.f"
		    for (jch = j; jch <= i__3; ++jch) {
#line 599 "zhgeqz.f"
			i__4 = jch + (jch + 1) * t_dim1;
#line 599 "zhgeqz.f"
			ctemp.r = t[i__4].r, ctemp.i = t[i__4].i;
#line 600 "zhgeqz.f"
			zlartg_(&ctemp, &t[jch + 1 + (jch + 1) * t_dim1], &
				c__, &s, &t[jch + (jch + 1) * t_dim1]);
#line 602 "zhgeqz.f"
			i__4 = jch + 1 + (jch + 1) * t_dim1;
#line 602 "zhgeqz.f"
			t[i__4].r = 0., t[i__4].i = 0.;
#line 603 "zhgeqz.f"
			if (jch < ilastm - 1) {
#line 603 "zhgeqz.f"
			    i__4 = ilastm - jch - 1;
#line 603 "zhgeqz.f"
			    zrot_(&i__4, &t[jch + (jch + 2) * t_dim1], ldt, &
				    t[jch + 1 + (jch + 2) * t_dim1], ldt, &
				    c__, &s);
#line 603 "zhgeqz.f"
			}
#line 606 "zhgeqz.f"
			i__4 = ilastm - jch + 2;
#line 606 "zhgeqz.f"
			zrot_(&i__4, &h__[jch + (jch - 1) * h_dim1], ldh, &
				h__[jch + 1 + (jch - 1) * h_dim1], ldh, &c__, 
				&s);
#line 608 "zhgeqz.f"
			if (ilq) {
#line 608 "zhgeqz.f"
			    d_cnjg(&z__1, &s);
#line 608 "zhgeqz.f"
			    zrot_(n, &q[jch * q_dim1 + 1], &c__1, &q[(jch + 1)
				     * q_dim1 + 1], &c__1, &c__, &z__1);
#line 608 "zhgeqz.f"
			}
#line 611 "zhgeqz.f"
			i__4 = jch + 1 + jch * h_dim1;
#line 611 "zhgeqz.f"
			ctemp.r = h__[i__4].r, ctemp.i = h__[i__4].i;
#line 612 "zhgeqz.f"
			zlartg_(&ctemp, &h__[jch + 1 + (jch - 1) * h_dim1], &
				c__, &s, &h__[jch + 1 + jch * h_dim1]);
#line 614 "zhgeqz.f"
			i__4 = jch + 1 + (jch - 1) * h_dim1;
#line 614 "zhgeqz.f"
			h__[i__4].r = 0., h__[i__4].i = 0.;
#line 615 "zhgeqz.f"
			i__4 = jch + 1 - ifrstm;
#line 615 "zhgeqz.f"
			zrot_(&i__4, &h__[ifrstm + jch * h_dim1], &c__1, &h__[
				ifrstm + (jch - 1) * h_dim1], &c__1, &c__, &s)
				;
#line 617 "zhgeqz.f"
			i__4 = jch - ifrstm;
#line 617 "zhgeqz.f"
			zrot_(&i__4, &t[ifrstm + jch * t_dim1], &c__1, &t[
				ifrstm + (jch - 1) * t_dim1], &c__1, &c__, &s)
				;
#line 619 "zhgeqz.f"
			if (ilz) {
#line 619 "zhgeqz.f"
			    zrot_(n, &z__[jch * z_dim1 + 1], &c__1, &z__[(jch 
				    - 1) * z_dim1 + 1], &c__1, &c__, &s);
#line 619 "zhgeqz.f"
			}
#line 622 "zhgeqz.f"
/* L30: */
#line 622 "zhgeqz.f"
		    }
#line 623 "zhgeqz.f"
		    goto L50;
#line 624 "zhgeqz.f"
		}
#line 625 "zhgeqz.f"
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

#line 629 "zhgeqz.f"
		ifirst = j;
#line 630 "zhgeqz.f"
		goto L70;
#line 631 "zhgeqz.f"
	    }

/*           Neither test passed -- try next J */

#line 635 "zhgeqz.f"
/* L40: */
#line 635 "zhgeqz.f"
	}

/*        (Drop-through is "impossible") */

#line 639 "zhgeqz.f"
	*info = (*n << 1) + 1;
#line 640 "zhgeqz.f"
	goto L210;

/*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a */
/*        1x1 block. */

#line 645 "zhgeqz.f"
L50:
#line 646 "zhgeqz.f"
	i__2 = ilast + ilast * h_dim1;
#line 646 "zhgeqz.f"
	ctemp.r = h__[i__2].r, ctemp.i = h__[i__2].i;
#line 647 "zhgeqz.f"
	zlartg_(&ctemp, &h__[ilast + (ilast - 1) * h_dim1], &c__, &s, &h__[
		ilast + ilast * h_dim1]);
#line 649 "zhgeqz.f"
	i__2 = ilast + (ilast - 1) * h_dim1;
#line 649 "zhgeqz.f"
	h__[i__2].r = 0., h__[i__2].i = 0.;
#line 650 "zhgeqz.f"
	i__2 = ilast - ifrstm;
#line 650 "zhgeqz.f"
	zrot_(&i__2, &h__[ifrstm + ilast * h_dim1], &c__1, &h__[ifrstm + (
		ilast - 1) * h_dim1], &c__1, &c__, &s);
#line 652 "zhgeqz.f"
	i__2 = ilast - ifrstm;
#line 652 "zhgeqz.f"
	zrot_(&i__2, &t[ifrstm + ilast * t_dim1], &c__1, &t[ifrstm + (ilast - 
		1) * t_dim1], &c__1, &c__, &s);
#line 654 "zhgeqz.f"
	if (ilz) {
#line 654 "zhgeqz.f"
	    zrot_(n, &z__[ilast * z_dim1 + 1], &c__1, &z__[(ilast - 1) * 
		    z_dim1 + 1], &c__1, &c__, &s);
#line 654 "zhgeqz.f"
	}

/*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA */

#line 659 "zhgeqz.f"
L60:
#line 660 "zhgeqz.f"
	absb = z_abs(&t[ilast + ilast * t_dim1]);
#line 661 "zhgeqz.f"
	if (absb > safmin) {
#line 662 "zhgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 662 "zhgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 662 "zhgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 662 "zhgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 663 "zhgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 663 "zhgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 664 "zhgeqz.f"
	    if (ilschr) {
#line 665 "zhgeqz.f"
		i__2 = ilast - ifrstm;
#line 665 "zhgeqz.f"
		zscal_(&i__2, &signbc, &t[ifrstm + ilast * t_dim1], &c__1);
#line 666 "zhgeqz.f"
		i__2 = ilast + 1 - ifrstm;
#line 666 "zhgeqz.f"
		zscal_(&i__2, &signbc, &h__[ifrstm + ilast * h_dim1], &c__1);
#line 668 "zhgeqz.f"
	    } else {
#line 669 "zhgeqz.f"
		i__2 = ilast + ilast * h_dim1;
#line 669 "zhgeqz.f"
		i__3 = ilast + ilast * h_dim1;
#line 669 "zhgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 669 "zhgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 670 "zhgeqz.f"
	    }
#line 671 "zhgeqz.f"
	    if (ilz) {
#line 671 "zhgeqz.f"
		zscal_(n, &signbc, &z__[ilast * z_dim1 + 1], &c__1);
#line 671 "zhgeqz.f"
	    }
#line 673 "zhgeqz.f"
	} else {
#line 674 "zhgeqz.f"
	    i__2 = ilast + ilast * t_dim1;
#line 674 "zhgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 675 "zhgeqz.f"
	}
#line 676 "zhgeqz.f"
	i__2 = ilast;
#line 676 "zhgeqz.f"
	i__3 = ilast + ilast * h_dim1;
#line 676 "zhgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 677 "zhgeqz.f"
	i__2 = ilast;
#line 677 "zhgeqz.f"
	i__3 = ilast + ilast * t_dim1;
#line 677 "zhgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;

/*        Go to next block -- exit if finished. */

#line 681 "zhgeqz.f"
	--ilast;
#line 682 "zhgeqz.f"
	if (ilast < *ilo) {
#line 682 "zhgeqz.f"
	    goto L190;
#line 682 "zhgeqz.f"
	}

/*        Reset counters */

#line 687 "zhgeqz.f"
	iiter = 0;
#line 688 "zhgeqz.f"
	eshift.r = 0., eshift.i = 0.;
#line 689 "zhgeqz.f"
	if (! ilschr) {
#line 690 "zhgeqz.f"
	    ilastm = ilast;
#line 691 "zhgeqz.f"
	    if (ifrstm > ilast) {
#line 691 "zhgeqz.f"
		ifrstm = *ilo;
#line 691 "zhgeqz.f"
	    }
#line 693 "zhgeqz.f"
	}
#line 694 "zhgeqz.f"
	goto L160;

/*        QZ step */

/*        This iteration only involves rows/columns IFIRST:ILAST.  We */
/*        assume IFIRST < ILAST, and that the diagonal of B is non-zero. */

#line 701 "zhgeqz.f"
L70:
#line 702 "zhgeqz.f"
	++iiter;
#line 703 "zhgeqz.f"
	if (! ilschr) {
#line 704 "zhgeqz.f"
	    ifrstm = ifirst;
#line 705 "zhgeqz.f"
	}

/*        Compute the Shift. */

/*        At this point, IFIRST < ILAST, and the diagonal elements of */
/*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in */
/*        magnitude) */

#line 713 "zhgeqz.f"
	if (iiter / 10 * 10 != iiter) {

/*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of */
/*           the bottom-right 2x2 block of A inv(B) which is nearest to */
/*           the bottom-right element. */

/*           We factor B as U*D, where U has unit diagonals, and */
/*           compute (A*inv(D))*inv(U). */

#line 722 "zhgeqz.f"
	    i__2 = ilast - 1 + ilast * t_dim1;
#line 722 "zhgeqz.f"
	    z__2.r = bscale * t[i__2].r, z__2.i = bscale * t[i__2].i;
#line 722 "zhgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 722 "zhgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 722 "zhgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 722 "zhgeqz.f"
	    u12.r = z__1.r, u12.i = z__1.i;
#line 724 "zhgeqz.f"
	    i__2 = ilast - 1 + (ilast - 1) * h_dim1;
#line 724 "zhgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 724 "zhgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 724 "zhgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 724 "zhgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 724 "zhgeqz.f"
	    ad11.r = z__1.r, ad11.i = z__1.i;
#line 726 "zhgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 726 "zhgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 726 "zhgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 726 "zhgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 726 "zhgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 726 "zhgeqz.f"
	    ad21.r = z__1.r, ad21.i = z__1.i;
#line 728 "zhgeqz.f"
	    i__2 = ilast - 1 + ilast * h_dim1;
#line 728 "zhgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 728 "zhgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 728 "zhgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 728 "zhgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 728 "zhgeqz.f"
	    ad12.r = z__1.r, ad12.i = z__1.i;
#line 730 "zhgeqz.f"
	    i__2 = ilast + ilast * h_dim1;
#line 730 "zhgeqz.f"
	    z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 730 "zhgeqz.f"
	    i__3 = ilast + ilast * t_dim1;
#line 730 "zhgeqz.f"
	    z__3.r = bscale * t[i__3].r, z__3.i = bscale * t[i__3].i;
#line 730 "zhgeqz.f"
	    z_div(&z__1, &z__2, &z__3);
#line 730 "zhgeqz.f"
	    ad22.r = z__1.r, ad22.i = z__1.i;
#line 732 "zhgeqz.f"
	    z__2.r = u12.r * ad21.r - u12.i * ad21.i, z__2.i = u12.r * ad21.i 
		    + u12.i * ad21.r;
#line 732 "zhgeqz.f"
	    z__1.r = ad22.r - z__2.r, z__1.i = ad22.i - z__2.i;
#line 732 "zhgeqz.f"
	    abi22.r = z__1.r, abi22.i = z__1.i;

#line 734 "zhgeqz.f"
	    z__2.r = ad11.r + abi22.r, z__2.i = ad11.i + abi22.i;
#line 734 "zhgeqz.f"
	    z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
#line 734 "zhgeqz.f"
	    t1.r = z__1.r, t1.i = z__1.i;
#line 735 "zhgeqz.f"
	    pow_zi(&z__4, &t1, &c__2);
#line 735 "zhgeqz.f"
	    z__5.r = ad12.r * ad21.r - ad12.i * ad21.i, z__5.i = ad12.r * 
		    ad21.i + ad12.i * ad21.r;
#line 735 "zhgeqz.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 735 "zhgeqz.f"
	    z__6.r = ad11.r * ad22.r - ad11.i * ad22.i, z__6.i = ad11.r * 
		    ad22.i + ad11.i * ad22.r;
#line 735 "zhgeqz.f"
	    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 735 "zhgeqz.f"
	    z_sqrt(&z__1, &z__2);
#line 735 "zhgeqz.f"
	    rtdisc.r = z__1.r, rtdisc.i = z__1.i;
#line 736 "zhgeqz.f"
	    z__1.r = t1.r - abi22.r, z__1.i = t1.i - abi22.i;
#line 736 "zhgeqz.f"
	    z__2.r = t1.r - abi22.r, z__2.i = t1.i - abi22.i;
#line 736 "zhgeqz.f"
	    temp = z__1.r * rtdisc.r + d_imag(&z__2) * d_imag(&rtdisc);
#line 738 "zhgeqz.f"
	    if (temp <= 0.) {
#line 739 "zhgeqz.f"
		z__1.r = t1.r + rtdisc.r, z__1.i = t1.i + rtdisc.i;
#line 739 "zhgeqz.f"
		shift.r = z__1.r, shift.i = z__1.i;
#line 740 "zhgeqz.f"
	    } else {
#line 741 "zhgeqz.f"
		z__1.r = t1.r - rtdisc.r, z__1.i = t1.i - rtdisc.i;
#line 741 "zhgeqz.f"
		shift.r = z__1.r, shift.i = z__1.i;
#line 742 "zhgeqz.f"
	    }
#line 743 "zhgeqz.f"
	} else {

/*           Exceptional shift.  Chosen for no particularly good reason. */

#line 747 "zhgeqz.f"
	    i__2 = ilast + (ilast - 1) * h_dim1;
#line 747 "zhgeqz.f"
	    z__3.r = ascale * h__[i__2].r, z__3.i = ascale * h__[i__2].i;
#line 747 "zhgeqz.f"
	    i__3 = ilast - 1 + (ilast - 1) * t_dim1;
#line 747 "zhgeqz.f"
	    z__4.r = bscale * t[i__3].r, z__4.i = bscale * t[i__3].i;
#line 747 "zhgeqz.f"
	    z_div(&z__2, &z__3, &z__4);
#line 747 "zhgeqz.f"
	    z__1.r = eshift.r + z__2.r, z__1.i = eshift.i + z__2.i;
#line 747 "zhgeqz.f"
	    eshift.r = z__1.r, eshift.i = z__1.i;
#line 749 "zhgeqz.f"
	    shift.r = eshift.r, shift.i = eshift.i;
#line 750 "zhgeqz.f"
	}

/*        Now check for two consecutive small subdiagonals. */

#line 754 "zhgeqz.f"
	i__2 = ifirst + 1;
#line 754 "zhgeqz.f"
	for (j = ilast - 1; j >= i__2; --j) {
#line 755 "zhgeqz.f"
	    istart = j;
#line 756 "zhgeqz.f"
	    i__3 = j + j * h_dim1;
#line 756 "zhgeqz.f"
	    z__2.r = ascale * h__[i__3].r, z__2.i = ascale * h__[i__3].i;
#line 756 "zhgeqz.f"
	    i__4 = j + j * t_dim1;
#line 756 "zhgeqz.f"
	    z__4.r = bscale * t[i__4].r, z__4.i = bscale * t[i__4].i;
#line 756 "zhgeqz.f"
	    z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		    z__4.i + shift.i * z__4.r;
#line 756 "zhgeqz.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 756 "zhgeqz.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 757 "zhgeqz.f"
	    temp = (d__1 = ctemp.r, abs(d__1)) + (d__2 = d_imag(&ctemp), abs(
		    d__2));
#line 758 "zhgeqz.f"
	    i__3 = j + 1 + j * h_dim1;
#line 758 "zhgeqz.f"
	    temp2 = ascale * ((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&h__[j + 1 + j * h_dim1]), abs(d__2)));
#line 759 "zhgeqz.f"
	    tempr = max(temp,temp2);
#line 760 "zhgeqz.f"
	    if (tempr < 1. && tempr != 0.) {
#line 761 "zhgeqz.f"
		temp /= tempr;
#line 762 "zhgeqz.f"
		temp2 /= tempr;
#line 763 "zhgeqz.f"
	    }
#line 764 "zhgeqz.f"
	    i__3 = j + (j - 1) * h_dim1;
#line 764 "zhgeqz.f"
	    if (((d__1 = h__[i__3].r, abs(d__1)) + (d__2 = d_imag(&h__[j + (j 
		    - 1) * h_dim1]), abs(d__2))) * temp2 <= temp * atol) {
#line 764 "zhgeqz.f"
		goto L90;
#line 764 "zhgeqz.f"
	    }
#line 766 "zhgeqz.f"
/* L80: */
#line 766 "zhgeqz.f"
	}

#line 768 "zhgeqz.f"
	istart = ifirst;
#line 769 "zhgeqz.f"
	i__2 = ifirst + ifirst * h_dim1;
#line 769 "zhgeqz.f"
	z__2.r = ascale * h__[i__2].r, z__2.i = ascale * h__[i__2].i;
#line 769 "zhgeqz.f"
	i__3 = ifirst + ifirst * t_dim1;
#line 769 "zhgeqz.f"
	z__4.r = bscale * t[i__3].r, z__4.i = bscale * t[i__3].i;
#line 769 "zhgeqz.f"
	z__3.r = shift.r * z__4.r - shift.i * z__4.i, z__3.i = shift.r * 
		z__4.i + shift.i * z__4.r;
#line 769 "zhgeqz.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 769 "zhgeqz.f"
	ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 771 "zhgeqz.f"
L90:

/*        Do an implicit-shift QZ sweep. */

/*        Initial Q */

#line 777 "zhgeqz.f"
	i__2 = istart + 1 + istart * h_dim1;
#line 777 "zhgeqz.f"
	z__1.r = ascale * h__[i__2].r, z__1.i = ascale * h__[i__2].i;
#line 777 "zhgeqz.f"
	ctemp2.r = z__1.r, ctemp2.i = z__1.i;
#line 778 "zhgeqz.f"
	zlartg_(&ctemp, &ctemp2, &c__, &s, &ctemp3);

/*        Sweep */

#line 782 "zhgeqz.f"
	i__2 = ilast - 1;
#line 782 "zhgeqz.f"
	for (j = istart; j <= i__2; ++j) {
#line 783 "zhgeqz.f"
	    if (j > istart) {
#line 784 "zhgeqz.f"
		i__3 = j + (j - 1) * h_dim1;
#line 784 "zhgeqz.f"
		ctemp.r = h__[i__3].r, ctemp.i = h__[i__3].i;
#line 785 "zhgeqz.f"
		zlartg_(&ctemp, &h__[j + 1 + (j - 1) * h_dim1], &c__, &s, &
			h__[j + (j - 1) * h_dim1]);
#line 786 "zhgeqz.f"
		i__3 = j + 1 + (j - 1) * h_dim1;
#line 786 "zhgeqz.f"
		h__[i__3].r = 0., h__[i__3].i = 0.;
#line 787 "zhgeqz.f"
	    }

#line 789 "zhgeqz.f"
	    i__3 = ilastm;
#line 789 "zhgeqz.f"
	    for (jc = j; jc <= i__3; ++jc) {
#line 790 "zhgeqz.f"
		i__4 = j + jc * h_dim1;
#line 790 "zhgeqz.f"
		z__2.r = c__ * h__[i__4].r, z__2.i = c__ * h__[i__4].i;
#line 790 "zhgeqz.f"
		i__5 = j + 1 + jc * h_dim1;
#line 790 "zhgeqz.f"
		z__3.r = s.r * h__[i__5].r - s.i * h__[i__5].i, z__3.i = s.r *
			 h__[i__5].i + s.i * h__[i__5].r;
#line 790 "zhgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 790 "zhgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 791 "zhgeqz.f"
		i__4 = j + 1 + jc * h_dim1;
#line 791 "zhgeqz.f"
		d_cnjg(&z__4, &s);
#line 791 "zhgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 791 "zhgeqz.f"
		i__5 = j + jc * h_dim1;
#line 791 "zhgeqz.f"
		z__2.r = z__3.r * h__[i__5].r - z__3.i * h__[i__5].i, z__2.i =
			 z__3.r * h__[i__5].i + z__3.i * h__[i__5].r;
#line 791 "zhgeqz.f"
		i__6 = j + 1 + jc * h_dim1;
#line 791 "zhgeqz.f"
		z__5.r = c__ * h__[i__6].r, z__5.i = c__ * h__[i__6].i;
#line 791 "zhgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 791 "zhgeqz.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 792 "zhgeqz.f"
		i__4 = j + jc * h_dim1;
#line 792 "zhgeqz.f"
		h__[i__4].r = ctemp.r, h__[i__4].i = ctemp.i;
#line 793 "zhgeqz.f"
		i__4 = j + jc * t_dim1;
#line 793 "zhgeqz.f"
		z__2.r = c__ * t[i__4].r, z__2.i = c__ * t[i__4].i;
#line 793 "zhgeqz.f"
		i__5 = j + 1 + jc * t_dim1;
#line 793 "zhgeqz.f"
		z__3.r = s.r * t[i__5].r - s.i * t[i__5].i, z__3.i = s.r * t[
			i__5].i + s.i * t[i__5].r;
#line 793 "zhgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 793 "zhgeqz.f"
		ctemp2.r = z__1.r, ctemp2.i = z__1.i;
#line 794 "zhgeqz.f"
		i__4 = j + 1 + jc * t_dim1;
#line 794 "zhgeqz.f"
		d_cnjg(&z__4, &s);
#line 794 "zhgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 794 "zhgeqz.f"
		i__5 = j + jc * t_dim1;
#line 794 "zhgeqz.f"
		z__2.r = z__3.r * t[i__5].r - z__3.i * t[i__5].i, z__2.i = 
			z__3.r * t[i__5].i + z__3.i * t[i__5].r;
#line 794 "zhgeqz.f"
		i__6 = j + 1 + jc * t_dim1;
#line 794 "zhgeqz.f"
		z__5.r = c__ * t[i__6].r, z__5.i = c__ * t[i__6].i;
#line 794 "zhgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 794 "zhgeqz.f"
		t[i__4].r = z__1.r, t[i__4].i = z__1.i;
#line 795 "zhgeqz.f"
		i__4 = j + jc * t_dim1;
#line 795 "zhgeqz.f"
		t[i__4].r = ctemp2.r, t[i__4].i = ctemp2.i;
#line 796 "zhgeqz.f"
/* L100: */
#line 796 "zhgeqz.f"
	    }
#line 797 "zhgeqz.f"
	    if (ilq) {
#line 798 "zhgeqz.f"
		i__3 = *n;
#line 798 "zhgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 799 "zhgeqz.f"
		    i__4 = jr + j * q_dim1;
#line 799 "zhgeqz.f"
		    z__2.r = c__ * q[i__4].r, z__2.i = c__ * q[i__4].i;
#line 799 "zhgeqz.f"
		    d_cnjg(&z__4, &s);
#line 799 "zhgeqz.f"
		    i__5 = jr + (j + 1) * q_dim1;
#line 799 "zhgeqz.f"
		    z__3.r = z__4.r * q[i__5].r - z__4.i * q[i__5].i, z__3.i =
			     z__4.r * q[i__5].i + z__4.i * q[i__5].r;
#line 799 "zhgeqz.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 799 "zhgeqz.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 800 "zhgeqz.f"
		    i__4 = jr + (j + 1) * q_dim1;
#line 800 "zhgeqz.f"
		    z__3.r = -s.r, z__3.i = -s.i;
#line 800 "zhgeqz.f"
		    i__5 = jr + j * q_dim1;
#line 800 "zhgeqz.f"
		    z__2.r = z__3.r * q[i__5].r - z__3.i * q[i__5].i, z__2.i =
			     z__3.r * q[i__5].i + z__3.i * q[i__5].r;
#line 800 "zhgeqz.f"
		    i__6 = jr + (j + 1) * q_dim1;
#line 800 "zhgeqz.f"
		    z__4.r = c__ * q[i__6].r, z__4.i = c__ * q[i__6].i;
#line 800 "zhgeqz.f"
		    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 800 "zhgeqz.f"
		    q[i__4].r = z__1.r, q[i__4].i = z__1.i;
#line 801 "zhgeqz.f"
		    i__4 = jr + j * q_dim1;
#line 801 "zhgeqz.f"
		    q[i__4].r = ctemp.r, q[i__4].i = ctemp.i;
#line 802 "zhgeqz.f"
/* L110: */
#line 802 "zhgeqz.f"
		}
#line 803 "zhgeqz.f"
	    }

#line 805 "zhgeqz.f"
	    i__3 = j + 1 + (j + 1) * t_dim1;
#line 805 "zhgeqz.f"
	    ctemp.r = t[i__3].r, ctemp.i = t[i__3].i;
#line 806 "zhgeqz.f"
	    zlartg_(&ctemp, &t[j + 1 + j * t_dim1], &c__, &s, &t[j + 1 + (j + 
		    1) * t_dim1]);
#line 807 "zhgeqz.f"
	    i__3 = j + 1 + j * t_dim1;
#line 807 "zhgeqz.f"
	    t[i__3].r = 0., t[i__3].i = 0.;

/* Computing MIN */
#line 809 "zhgeqz.f"
	    i__4 = j + 2;
#line 809 "zhgeqz.f"
	    i__3 = min(i__4,ilast);
#line 809 "zhgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 810 "zhgeqz.f"
		i__4 = jr + (j + 1) * h_dim1;
#line 810 "zhgeqz.f"
		z__2.r = c__ * h__[i__4].r, z__2.i = c__ * h__[i__4].i;
#line 810 "zhgeqz.f"
		i__5 = jr + j * h_dim1;
#line 810 "zhgeqz.f"
		z__3.r = s.r * h__[i__5].r - s.i * h__[i__5].i, z__3.i = s.r *
			 h__[i__5].i + s.i * h__[i__5].r;
#line 810 "zhgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 810 "zhgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 811 "zhgeqz.f"
		i__4 = jr + j * h_dim1;
#line 811 "zhgeqz.f"
		d_cnjg(&z__4, &s);
#line 811 "zhgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 811 "zhgeqz.f"
		i__5 = jr + (j + 1) * h_dim1;
#line 811 "zhgeqz.f"
		z__2.r = z__3.r * h__[i__5].r - z__3.i * h__[i__5].i, z__2.i =
			 z__3.r * h__[i__5].i + z__3.i * h__[i__5].r;
#line 811 "zhgeqz.f"
		i__6 = jr + j * h_dim1;
#line 811 "zhgeqz.f"
		z__5.r = c__ * h__[i__6].r, z__5.i = c__ * h__[i__6].i;
#line 811 "zhgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 811 "zhgeqz.f"
		h__[i__4].r = z__1.r, h__[i__4].i = z__1.i;
#line 812 "zhgeqz.f"
		i__4 = jr + (j + 1) * h_dim1;
#line 812 "zhgeqz.f"
		h__[i__4].r = ctemp.r, h__[i__4].i = ctemp.i;
#line 813 "zhgeqz.f"
/* L120: */
#line 813 "zhgeqz.f"
	    }
#line 814 "zhgeqz.f"
	    i__3 = j;
#line 814 "zhgeqz.f"
	    for (jr = ifrstm; jr <= i__3; ++jr) {
#line 815 "zhgeqz.f"
		i__4 = jr + (j + 1) * t_dim1;
#line 815 "zhgeqz.f"
		z__2.r = c__ * t[i__4].r, z__2.i = c__ * t[i__4].i;
#line 815 "zhgeqz.f"
		i__5 = jr + j * t_dim1;
#line 815 "zhgeqz.f"
		z__3.r = s.r * t[i__5].r - s.i * t[i__5].i, z__3.i = s.r * t[
			i__5].i + s.i * t[i__5].r;
#line 815 "zhgeqz.f"
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 815 "zhgeqz.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 816 "zhgeqz.f"
		i__4 = jr + j * t_dim1;
#line 816 "zhgeqz.f"
		d_cnjg(&z__4, &s);
#line 816 "zhgeqz.f"
		z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 816 "zhgeqz.f"
		i__5 = jr + (j + 1) * t_dim1;
#line 816 "zhgeqz.f"
		z__2.r = z__3.r * t[i__5].r - z__3.i * t[i__5].i, z__2.i = 
			z__3.r * t[i__5].i + z__3.i * t[i__5].r;
#line 816 "zhgeqz.f"
		i__6 = jr + j * t_dim1;
#line 816 "zhgeqz.f"
		z__5.r = c__ * t[i__6].r, z__5.i = c__ * t[i__6].i;
#line 816 "zhgeqz.f"
		z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 816 "zhgeqz.f"
		t[i__4].r = z__1.r, t[i__4].i = z__1.i;
#line 817 "zhgeqz.f"
		i__4 = jr + (j + 1) * t_dim1;
#line 817 "zhgeqz.f"
		t[i__4].r = ctemp.r, t[i__4].i = ctemp.i;
#line 818 "zhgeqz.f"
/* L130: */
#line 818 "zhgeqz.f"
	    }
#line 819 "zhgeqz.f"
	    if (ilz) {
#line 820 "zhgeqz.f"
		i__3 = *n;
#line 820 "zhgeqz.f"
		for (jr = 1; jr <= i__3; ++jr) {
#line 821 "zhgeqz.f"
		    i__4 = jr + (j + 1) * z_dim1;
#line 821 "zhgeqz.f"
		    z__2.r = c__ * z__[i__4].r, z__2.i = c__ * z__[i__4].i;
#line 821 "zhgeqz.f"
		    i__5 = jr + j * z_dim1;
#line 821 "zhgeqz.f"
		    z__3.r = s.r * z__[i__5].r - s.i * z__[i__5].i, z__3.i = 
			    s.r * z__[i__5].i + s.i * z__[i__5].r;
#line 821 "zhgeqz.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 821 "zhgeqz.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 822 "zhgeqz.f"
		    i__4 = jr + j * z_dim1;
#line 822 "zhgeqz.f"
		    d_cnjg(&z__4, &s);
#line 822 "zhgeqz.f"
		    z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 822 "zhgeqz.f"
		    i__5 = jr + (j + 1) * z_dim1;
#line 822 "zhgeqz.f"
		    z__2.r = z__3.r * z__[i__5].r - z__3.i * z__[i__5].i, 
			    z__2.i = z__3.r * z__[i__5].i + z__3.i * z__[i__5]
			    .r;
#line 822 "zhgeqz.f"
		    i__6 = jr + j * z_dim1;
#line 822 "zhgeqz.f"
		    z__5.r = c__ * z__[i__6].r, z__5.i = c__ * z__[i__6].i;
#line 822 "zhgeqz.f"
		    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 822 "zhgeqz.f"
		    z__[i__4].r = z__1.r, z__[i__4].i = z__1.i;
#line 823 "zhgeqz.f"
		    i__4 = jr + (j + 1) * z_dim1;
#line 823 "zhgeqz.f"
		    z__[i__4].r = ctemp.r, z__[i__4].i = ctemp.i;
#line 824 "zhgeqz.f"
/* L140: */
#line 824 "zhgeqz.f"
		}
#line 825 "zhgeqz.f"
	    }
#line 826 "zhgeqz.f"
/* L150: */
#line 826 "zhgeqz.f"
	}

#line 828 "zhgeqz.f"
L160:

#line 830 "zhgeqz.f"
/* L170: */
#line 830 "zhgeqz.f"
	;
#line 830 "zhgeqz.f"
    }

/*     Drop-through = non-convergence */

#line 834 "zhgeqz.f"
L180:
#line 835 "zhgeqz.f"
    *info = ilast;
#line 836 "zhgeqz.f"
    goto L210;

/*     Successful completion of all QZ steps */

#line 840 "zhgeqz.f"
L190:

/*     Set Eigenvalues 1:ILO-1 */

#line 844 "zhgeqz.f"
    i__1 = *ilo - 1;
#line 844 "zhgeqz.f"
    for (j = 1; j <= i__1; ++j) {
#line 845 "zhgeqz.f"
	absb = z_abs(&t[j + j * t_dim1]);
#line 846 "zhgeqz.f"
	if (absb > safmin) {
#line 847 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 847 "zhgeqz.f"
	    z__2.r = t[i__2].r / absb, z__2.i = t[i__2].i / absb;
#line 847 "zhgeqz.f"
	    d_cnjg(&z__1, &z__2);
#line 847 "zhgeqz.f"
	    signbc.r = z__1.r, signbc.i = z__1.i;
#line 848 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 848 "zhgeqz.f"
	    t[i__2].r = absb, t[i__2].i = 0.;
#line 849 "zhgeqz.f"
	    if (ilschr) {
#line 850 "zhgeqz.f"
		i__2 = j - 1;
#line 850 "zhgeqz.f"
		zscal_(&i__2, &signbc, &t[j * t_dim1 + 1], &c__1);
#line 851 "zhgeqz.f"
		zscal_(&j, &signbc, &h__[j * h_dim1 + 1], &c__1);
#line 852 "zhgeqz.f"
	    } else {
#line 853 "zhgeqz.f"
		i__2 = j + j * h_dim1;
#line 853 "zhgeqz.f"
		i__3 = j + j * h_dim1;
#line 853 "zhgeqz.f"
		z__1.r = h__[i__3].r * signbc.r - h__[i__3].i * signbc.i, 
			z__1.i = h__[i__3].r * signbc.i + h__[i__3].i * 
			signbc.r;
#line 853 "zhgeqz.f"
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
#line 854 "zhgeqz.f"
	    }
#line 855 "zhgeqz.f"
	    if (ilz) {
#line 855 "zhgeqz.f"
		zscal_(n, &signbc, &z__[j * z_dim1 + 1], &c__1);
#line 855 "zhgeqz.f"
	    }
#line 857 "zhgeqz.f"
	} else {
#line 858 "zhgeqz.f"
	    i__2 = j + j * t_dim1;
#line 858 "zhgeqz.f"
	    t[i__2].r = 0., t[i__2].i = 0.;
#line 859 "zhgeqz.f"
	}
#line 860 "zhgeqz.f"
	i__2 = j;
#line 860 "zhgeqz.f"
	i__3 = j + j * h_dim1;
#line 860 "zhgeqz.f"
	alpha[i__2].r = h__[i__3].r, alpha[i__2].i = h__[i__3].i;
#line 861 "zhgeqz.f"
	i__2 = j;
#line 861 "zhgeqz.f"
	i__3 = j + j * t_dim1;
#line 861 "zhgeqz.f"
	beta[i__2].r = t[i__3].r, beta[i__2].i = t[i__3].i;
#line 862 "zhgeqz.f"
/* L200: */
#line 862 "zhgeqz.f"
    }

/*     Normal Termination */

#line 866 "zhgeqz.f"
    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size */

#line 870 "zhgeqz.f"
L210:
#line 871 "zhgeqz.f"
    z__1.r = (doublereal) (*n), z__1.i = 0.;
#line 871 "zhgeqz.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 872 "zhgeqz.f"
    return 0;

/*     End of ZHGEQZ */

} /* zhgeqz_ */

