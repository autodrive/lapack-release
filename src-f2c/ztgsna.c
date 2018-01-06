#line 1 "ztgsna.f"
/* ztgsna.f -- translated by f2c (version 20100827).
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

#line 1 "ztgsna.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b19 = {1.,0.};
static doublecomplex c_b20 = {0.,0.};
static logical c_false = FALSE_;
static integer c__3 = 3;

/* > \brief \b ZTGSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTGSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, */
/*                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   DIF( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or eigenvectors of a matrix pair (A, B). */
/* > */
/* > (A, B) must be in generalized Schur canonical form, that is, A and */
/* > B are both upper triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies whether condition numbers are required for */
/* >          eigenvalues (S) or eigenvectors (DIF): */
/* >          = 'E': for eigenvalues only (S); */
/* >          = 'V': for eigenvectors only (DIF); */
/* >          = 'B': for both eigenvalues and eigenvectors (S and DIF). */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* >          HOWMNY is CHARACTER*1 */
/* >          = 'A': compute condition numbers for all eigenpairs; */
/* >          = 'S': compute condition numbers for selected eigenpairs */
/* >                 specified by the array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          If HOWMNY = 'S', SELECT specifies the eigenpairs for which */
/* >          condition numbers are required. To select condition numbers */
/* >          for the corresponding j-th eigenvalue and/or eigenvector, */
/* >          SELECT(j) must be set to .TRUE.. */
/* >          If HOWMNY = 'A', SELECT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the square matrix pair (A, B). N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The upper triangular matrix A in the pair (A,B). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          The upper triangular matrix B in the pair (A, B). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >          VL is COMPLEX*16 array, dimension (LDVL,M) */
/* >          IF JOB = 'E' or 'B', VL must contain left eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT.  The eigenvectors must be stored in consecutive */
/* >          columns of VL, as returned by ZTGEVC. */
/* >          If JOB = 'V', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL. LDVL >= 1; and */
/* >          If JOB = 'E' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] VR */
/* > \verbatim */
/* >          VR is COMPLEX*16 array, dimension (LDVR,M) */
/* >          IF JOB = 'E' or 'B', VR must contain right eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT.  The eigenvectors must be stored in consecutive */
/* >          columns of VR, as returned by ZTGEVC. */
/* >          If JOB = 'V', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR. LDVR >= 1; */
/* >          If JOB = 'E' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* >          selected eigenvalues, stored in consecutive elements of the */
/* >          array. */
/* >          If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the selected eigenvectors, stored in consecutive */
/* >          elements of the array. */
/* >          If the eigenvalues cannot be reordered to compute DIF(j), */
/* >          DIF(j) is set to 0; this can only occur when the true value */
/* >          would be very small anyway. */
/* >          For each eigenvalue/vector specified by SELECT, DIF stores */
/* >          a Frobenius norm-based estimate of Difl. */
/* >          If JOB = 'E', DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* >          MM is INTEGER */
/* >          The number of elements in the arrays S and DIF. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of elements of the arrays S and DIF used to store */
/* >          the specified condition numbers; for each selected eigenvalue */
/* >          one element is used. If HOWMNY = 'A', M is set to N. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N). */
/* >          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N+2) */
/* >          If JOB = 'E', IWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: Successful exit */
/* >          < 0: If INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The reciprocal of the condition number of the i-th generalized */
/* >  eigenvalue w = (a, b) is defined as */
/* > */
/* >          S(I) = (|v**HAu|**2 + |v**HBu|**2)**(1/2) / (norm(u)*norm(v)) */
/* > */
/* >  where u and v are the right and left eigenvectors of (A, B) */
/* >  corresponding to w; |z| denotes the absolute value of the complex */
/* >  number, and norm(u) denotes the 2-norm of the vector u. The pair */
/* >  (a, b) corresponds to an eigenvalue w = a/b (= v**HAu/v**HBu) of the */
/* >  matrix pair (A, B). If both a and b equal zero, then (A,B) is */
/* >  singular and S(I) = -1 is returned. */
/* > */
/* >  An approximate error bound on the chordal distance between the i-th */
/* >  computed generalized eigenvalue w and the corresponding exact */
/* >  eigenvalue lambda is */
/* > */
/* >          chord(w, lambda) <=   EPS * norm(A, B) / S(I), */
/* > */
/* >  where EPS is the machine precision. */
/* > */
/* >  The reciprocal of the condition number of the right eigenvector u */
/* >  and left eigenvector v corresponding to the generalized eigenvalue w */
/* >  is defined as follows. Suppose */
/* > */
/* >                   (A, B) = ( a   *  ) ( b  *  )  1 */
/* >                            ( 0  A22 ),( 0 B22 )  n-1 */
/* >                              1  n-1     1 n-1 */
/* > */
/* >  Then the reciprocal condition number DIF(I) is */
/* > */
/* >          Difl[(a, b), (A22, B22)]  = sigma-min( Zl ) */
/* > */
/* >  where sigma-min(Zl) denotes the smallest singular value of */
/* > */
/* >         Zl = [ kron(a, In-1) -kron(1, A22) ] */
/* >              [ kron(b, In-1) -kron(1, B22) ]. */
/* > */
/* >  Here In-1 is the identity matrix of size n-1 and X**H is the conjugate */
/* >  transpose of X. kron(X, Y) is the Kronecker product between the */
/* >  matrices X and Y. */
/* > */
/* >  We approximate the smallest singular value of Zl with an upper */
/* >  bound. This is done by ZLATDF. */
/* > */
/* >  An approximate error bound for a computed eigenvector VL(i) or */
/* >  VR(i) is given by */
/* > */
/* >                      EPS * norm(A, B) / DIF(i). */
/* > */
/* >  See ref. [2-3] for more details and further references. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > */
/* >  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified */
/* >      Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* >      Estimation: Theory, Algorithms and Software, Report */
/* >      UMINF - 94.04, Department of Computing Science, Umea University, */
/* >      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87. */
/* >      To appear in Numerical Algorithms, 1996. */
/* > */
/* >  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* >      Note 75. */
/* >      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztgsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *info, 
	ftnlen job_len, ftnlen howmny_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, k, n1, n2, ks;
    static doublereal eps, cond;
    static integer ierr, ifst;
    static doublereal lnrm;
    static doublecomplex yhax, yhbx;
    static integer ilst;
    static doublereal rnrm, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer lwmin;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical wants;
    static doublecomplex dummy[1];
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    static doublecomplex dummy1[1];
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical wantbh, wantdf, somcon;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    ztgexc_(logical *, logical *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *);
    static doublereal smlnum;
    static logical lquery;
    extern /* Subroutine */ int ztgsyl_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, integer *,
	     integer *, ftnlen);


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

/*     Decode and test the input parameters */

#line 363 "ztgsna.f"
    /* Parameter adjustments */
#line 363 "ztgsna.f"
    --select;
#line 363 "ztgsna.f"
    a_dim1 = *lda;
#line 363 "ztgsna.f"
    a_offset = 1 + a_dim1;
#line 363 "ztgsna.f"
    a -= a_offset;
#line 363 "ztgsna.f"
    b_dim1 = *ldb;
#line 363 "ztgsna.f"
    b_offset = 1 + b_dim1;
#line 363 "ztgsna.f"
    b -= b_offset;
#line 363 "ztgsna.f"
    vl_dim1 = *ldvl;
#line 363 "ztgsna.f"
    vl_offset = 1 + vl_dim1;
#line 363 "ztgsna.f"
    vl -= vl_offset;
#line 363 "ztgsna.f"
    vr_dim1 = *ldvr;
#line 363 "ztgsna.f"
    vr_offset = 1 + vr_dim1;
#line 363 "ztgsna.f"
    vr -= vr_offset;
#line 363 "ztgsna.f"
    --s;
#line 363 "ztgsna.f"
    --dif;
#line 363 "ztgsna.f"
    --work;
#line 363 "ztgsna.f"
    --iwork;
#line 363 "ztgsna.f"

#line 363 "ztgsna.f"
    /* Function Body */
#line 363 "ztgsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 364 "ztgsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 365 "ztgsna.f"
    wantdf = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 367 "ztgsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 369 "ztgsna.f"
    *info = 0;
#line 370 "ztgsna.f"
    lquery = *lwork == -1;

#line 372 "ztgsna.f"
    if (! wants && ! wantdf) {
#line 373 "ztgsna.f"
	*info = -1;
#line 374 "ztgsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 375 "ztgsna.f"
	*info = -2;
#line 376 "ztgsna.f"
    } else if (*n < 0) {
#line 377 "ztgsna.f"
	*info = -4;
#line 378 "ztgsna.f"
    } else if (*lda < max(1,*n)) {
#line 379 "ztgsna.f"
	*info = -6;
#line 380 "ztgsna.f"
    } else if (*ldb < max(1,*n)) {
#line 381 "ztgsna.f"
	*info = -8;
#line 382 "ztgsna.f"
    } else if (wants && *ldvl < *n) {
#line 383 "ztgsna.f"
	*info = -10;
#line 384 "ztgsna.f"
    } else if (wants && *ldvr < *n) {
#line 385 "ztgsna.f"
	*info = -12;
#line 386 "ztgsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 391 "ztgsna.f"
	if (somcon) {
#line 392 "ztgsna.f"
	    *m = 0;
#line 393 "ztgsna.f"
	    i__1 = *n;
#line 393 "ztgsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 394 "ztgsna.f"
		if (select[k]) {
#line 394 "ztgsna.f"
		    ++(*m);
#line 394 "ztgsna.f"
		}
#line 396 "ztgsna.f"
/* L10: */
#line 396 "ztgsna.f"
	    }
#line 397 "ztgsna.f"
	} else {
#line 398 "ztgsna.f"
	    *m = *n;
#line 399 "ztgsna.f"
	}

#line 401 "ztgsna.f"
	if (*n == 0) {
#line 402 "ztgsna.f"
	    lwmin = 1;
#line 403 "ztgsna.f"
	} else if (lsame_(job, "V", (ftnlen)1, (ftnlen)1) || lsame_(job, 
		"B", (ftnlen)1, (ftnlen)1)) {
#line 404 "ztgsna.f"
	    lwmin = (*n << 1) * *n;
#line 405 "ztgsna.f"
	} else {
#line 406 "ztgsna.f"
	    lwmin = *n;
#line 407 "ztgsna.f"
	}
#line 408 "ztgsna.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 410 "ztgsna.f"
	if (*mm < *m) {
#line 411 "ztgsna.f"
	    *info = -15;
#line 412 "ztgsna.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 413 "ztgsna.f"
	    *info = -18;
#line 414 "ztgsna.f"
	}
#line 415 "ztgsna.f"
    }

#line 417 "ztgsna.f"
    if (*info != 0) {
#line 418 "ztgsna.f"
	i__1 = -(*info);
#line 418 "ztgsna.f"
	xerbla_("ZTGSNA", &i__1, (ftnlen)6);
#line 419 "ztgsna.f"
	return 0;
#line 420 "ztgsna.f"
    } else if (lquery) {
#line 421 "ztgsna.f"
	return 0;
#line 422 "ztgsna.f"
    }

/*     Quick return if possible */

#line 426 "ztgsna.f"
    if (*n == 0) {
#line 426 "ztgsna.f"
	return 0;
#line 426 "ztgsna.f"
    }

/*     Get machine constants */

#line 431 "ztgsna.f"
    eps = dlamch_("P", (ftnlen)1);
#line 432 "ztgsna.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 433 "ztgsna.f"
    bignum = 1. / smlnum;
#line 434 "ztgsna.f"
    dlabad_(&smlnum, &bignum);
#line 435 "ztgsna.f"
    ks = 0;
#line 436 "ztgsna.f"
    i__1 = *n;
#line 436 "ztgsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 441 "ztgsna.f"
	if (somcon) {
#line 442 "ztgsna.f"
	    if (! select[k]) {
#line 442 "ztgsna.f"
		goto L20;
#line 442 "ztgsna.f"
	    }
#line 444 "ztgsna.f"
	}

#line 446 "ztgsna.f"
	++ks;

#line 448 "ztgsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 453 "ztgsna.f"
	    rnrm = dznrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 454 "ztgsna.f"
	    lnrm = dznrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 455 "ztgsna.f"
	    zgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 + 1]
		    , &c__1, &c_b20, &work[1], &c__1, (ftnlen)1);
#line 457 "ztgsna.f"
	    zdotc_(&z__1, n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
#line 457 "ztgsna.f"
	    yhax.r = z__1.r, yhax.i = z__1.i;
#line 458 "ztgsna.f"
	    zgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 + 1]
		    , &c__1, &c_b20, &work[1], &c__1, (ftnlen)1);
#line 460 "ztgsna.f"
	    zdotc_(&z__1, n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
#line 460 "ztgsna.f"
	    yhbx.r = z__1.r, yhbx.i = z__1.i;
#line 461 "ztgsna.f"
	    d__1 = z_abs(&yhax);
#line 461 "ztgsna.f"
	    d__2 = z_abs(&yhbx);
#line 461 "ztgsna.f"
	    cond = dlapy2_(&d__1, &d__2);
#line 462 "ztgsna.f"
	    if (cond == 0.) {
#line 463 "ztgsna.f"
		s[ks] = -1.;
#line 464 "ztgsna.f"
	    } else {
#line 465 "ztgsna.f"
		s[ks] = cond / (rnrm * lnrm);
#line 466 "ztgsna.f"
	    }
#line 467 "ztgsna.f"
	}

#line 469 "ztgsna.f"
	if (wantdf) {
#line 470 "ztgsna.f"
	    if (*n == 1) {
#line 471 "ztgsna.f"
		d__1 = z_abs(&a[a_dim1 + 1]);
#line 471 "ztgsna.f"
		d__2 = z_abs(&b[b_dim1 + 1]);
#line 471 "ztgsna.f"
		dif[ks] = dlapy2_(&d__1, &d__2);
#line 472 "ztgsna.f"
	    } else {

/*              Estimate the reciprocal condition number of the k-th */
/*              eigenvectors. */

/*              Copy the matrix (A, B) to the array WORK and move the */
/*              (k,k)th pair to the (1,1) position. */

#line 480 "ztgsna.f"
		zlacpy_("Full", n, n, &a[a_offset], lda, &work[1], n, (ftnlen)
			4);
#line 481 "ztgsna.f"
		zlacpy_("Full", n, n, &b[b_offset], ldb, &work[*n * *n + 1], 
			n, (ftnlen)4);
#line 482 "ztgsna.f"
		ifst = k;
#line 483 "ztgsna.f"
		ilst = 1;

#line 485 "ztgsna.f"
		ztgexc_(&c_false, &c_false, n, &work[1], n, &work[*n * *n + 1]
			, n, dummy, &c__1, dummy1, &c__1, &ifst, &ilst, &ierr)
			;

#line 488 "ztgsna.f"
		if (ierr > 0) {

/*                 Ill-conditioned problem - swap rejected. */

#line 492 "ztgsna.f"
		    dif[ks] = 0.;
#line 493 "ztgsna.f"
		} else {

/*                 Reordering successful, solve generalized Sylvester */
/*                 equation for R and L, */
/*                            A22 * R - L * A11 = A12 */
/*                            B22 * R - L * B11 = B12, */
/*                 and compute estimate of Difl[(A11,B11), (A22, B22)]. */

#line 501 "ztgsna.f"
		    n1 = 1;
#line 502 "ztgsna.f"
		    n2 = *n - n1;
#line 503 "ztgsna.f"
		    i__ = *n * *n + 1;
#line 504 "ztgsna.f"
		    ztgsyl_("N", &c__3, &n2, &n1, &work[*n * n1 + n1 + 1], n, 
			    &work[1], n, &work[n1 + 1], n, &work[*n * n1 + n1 
			    + i__], n, &work[i__], n, &work[n1 + i__], n, &
			    scale, &dif[ks], dummy, &c__1, &iwork[1], &ierr, (
			    ftnlen)1);
#line 509 "ztgsna.f"
		}
#line 510 "ztgsna.f"
	    }
#line 511 "ztgsna.f"
	}

#line 513 "ztgsna.f"
L20:
#line 513 "ztgsna.f"
	;
#line 513 "ztgsna.f"
    }
#line 514 "ztgsna.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 515 "ztgsna.f"
    return 0;

/*     End of ZTGSNA */

} /* ztgsna_ */

