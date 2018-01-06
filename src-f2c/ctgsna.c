#line 1 "ctgsna.f"
/* ctgsna.f -- translated by f2c (version 20100827).
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

#line 1 "ctgsna.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b19 = {1.,0.};
static doublecomplex c_b20 = {0.,0.};
static logical c_false = FALSE_;
static integer c__3 = 3;

/* > \brief \b CTGSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, */
/*                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       REAL               DIF( * ), S( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), VL( LDVL, * ), */
/*      $                   VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGSNA estimates reciprocal condition numbers for specified */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          VL is COMPLEX array, dimension (LDVL,M) */
/* >          IF JOB = 'E' or 'B', VL must contain left eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT.  The eigenvectors must be stored in consecutive */
/* >          columns of VL, as returned by CTGEVC. */
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
/* >          VR is COMPLEX array, dimension (LDVR,M) */
/* >          IF JOB = 'E' or 'B', VR must contain right eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT.  The eigenvectors must be stored in consecutive */
/* >          columns of VR, as returned by CTGEVC. */
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
/* >          S is REAL array, dimension (MM) */
/* >          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* >          selected eigenvalues, stored in consecutive elements of the */
/* >          array. */
/* >          If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is REAL array, dimension (MM) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complexOTHERcomputational */

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
/* >  bound. This is done by CLATDF. */
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
/* Subroutine */ int ctgsna_(char *job, char *howmny, logical *select, 
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
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer lwmin;
    static logical wants;
    static doublecomplex dummy[1];
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *), slapy2_(
	    doublereal *, doublereal *);
    static doublecomplex dummy1[1];
    extern /* Subroutine */ int slabad_(doublereal *, doublereal *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    ctgexc_(logical *, logical *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical wantbh, wantdf, somcon;
    extern /* Subroutine */ int ctgsyl_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, integer *,
	     integer *, ftnlen);
    static doublereal smlnum;
    static logical lquery;


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

#line 363 "ctgsna.f"
    /* Parameter adjustments */
#line 363 "ctgsna.f"
    --select;
#line 363 "ctgsna.f"
    a_dim1 = *lda;
#line 363 "ctgsna.f"
    a_offset = 1 + a_dim1;
#line 363 "ctgsna.f"
    a -= a_offset;
#line 363 "ctgsna.f"
    b_dim1 = *ldb;
#line 363 "ctgsna.f"
    b_offset = 1 + b_dim1;
#line 363 "ctgsna.f"
    b -= b_offset;
#line 363 "ctgsna.f"
    vl_dim1 = *ldvl;
#line 363 "ctgsna.f"
    vl_offset = 1 + vl_dim1;
#line 363 "ctgsna.f"
    vl -= vl_offset;
#line 363 "ctgsna.f"
    vr_dim1 = *ldvr;
#line 363 "ctgsna.f"
    vr_offset = 1 + vr_dim1;
#line 363 "ctgsna.f"
    vr -= vr_offset;
#line 363 "ctgsna.f"
    --s;
#line 363 "ctgsna.f"
    --dif;
#line 363 "ctgsna.f"
    --work;
#line 363 "ctgsna.f"
    --iwork;
#line 363 "ctgsna.f"

#line 363 "ctgsna.f"
    /* Function Body */
#line 363 "ctgsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 364 "ctgsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 365 "ctgsna.f"
    wantdf = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 367 "ctgsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 369 "ctgsna.f"
    *info = 0;
#line 370 "ctgsna.f"
    lquery = *lwork == -1;

#line 372 "ctgsna.f"
    if (! wants && ! wantdf) {
#line 373 "ctgsna.f"
	*info = -1;
#line 374 "ctgsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 375 "ctgsna.f"
	*info = -2;
#line 376 "ctgsna.f"
    } else if (*n < 0) {
#line 377 "ctgsna.f"
	*info = -4;
#line 378 "ctgsna.f"
    } else if (*lda < max(1,*n)) {
#line 379 "ctgsna.f"
	*info = -6;
#line 380 "ctgsna.f"
    } else if (*ldb < max(1,*n)) {
#line 381 "ctgsna.f"
	*info = -8;
#line 382 "ctgsna.f"
    } else if (wants && *ldvl < *n) {
#line 383 "ctgsna.f"
	*info = -10;
#line 384 "ctgsna.f"
    } else if (wants && *ldvr < *n) {
#line 385 "ctgsna.f"
	*info = -12;
#line 386 "ctgsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 391 "ctgsna.f"
	if (somcon) {
#line 392 "ctgsna.f"
	    *m = 0;
#line 393 "ctgsna.f"
	    i__1 = *n;
#line 393 "ctgsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 394 "ctgsna.f"
		if (select[k]) {
#line 394 "ctgsna.f"
		    ++(*m);
#line 394 "ctgsna.f"
		}
#line 396 "ctgsna.f"
/* L10: */
#line 396 "ctgsna.f"
	    }
#line 397 "ctgsna.f"
	} else {
#line 398 "ctgsna.f"
	    *m = *n;
#line 399 "ctgsna.f"
	}

#line 401 "ctgsna.f"
	if (*n == 0) {
#line 402 "ctgsna.f"
	    lwmin = 1;
#line 403 "ctgsna.f"
	} else if (lsame_(job, "V", (ftnlen)1, (ftnlen)1) || lsame_(job, 
		"B", (ftnlen)1, (ftnlen)1)) {
#line 404 "ctgsna.f"
	    lwmin = (*n << 1) * *n;
#line 405 "ctgsna.f"
	} else {
#line 406 "ctgsna.f"
	    lwmin = *n;
#line 407 "ctgsna.f"
	}
#line 408 "ctgsna.f"
	work[1].r = (doublereal) lwmin, work[1].i = 0.;

#line 410 "ctgsna.f"
	if (*mm < *m) {
#line 411 "ctgsna.f"
	    *info = -15;
#line 412 "ctgsna.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 413 "ctgsna.f"
	    *info = -18;
#line 414 "ctgsna.f"
	}
#line 415 "ctgsna.f"
    }

#line 417 "ctgsna.f"
    if (*info != 0) {
#line 418 "ctgsna.f"
	i__1 = -(*info);
#line 418 "ctgsna.f"
	xerbla_("CTGSNA", &i__1, (ftnlen)6);
#line 419 "ctgsna.f"
	return 0;
#line 420 "ctgsna.f"
    } else if (lquery) {
#line 421 "ctgsna.f"
	return 0;
#line 422 "ctgsna.f"
    }

/*     Quick return if possible */

#line 426 "ctgsna.f"
    if (*n == 0) {
#line 426 "ctgsna.f"
	return 0;
#line 426 "ctgsna.f"
    }

/*     Get machine constants */

#line 431 "ctgsna.f"
    eps = slamch_("P", (ftnlen)1);
#line 432 "ctgsna.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 433 "ctgsna.f"
    bignum = 1. / smlnum;
#line 434 "ctgsna.f"
    slabad_(&smlnum, &bignum);
#line 435 "ctgsna.f"
    ks = 0;
#line 436 "ctgsna.f"
    i__1 = *n;
#line 436 "ctgsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 441 "ctgsna.f"
	if (somcon) {
#line 442 "ctgsna.f"
	    if (! select[k]) {
#line 442 "ctgsna.f"
		goto L20;
#line 442 "ctgsna.f"
	    }
#line 444 "ctgsna.f"
	}

#line 446 "ctgsna.f"
	++ks;

#line 448 "ctgsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 453 "ctgsna.f"
	    rnrm = scnrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 454 "ctgsna.f"
	    lnrm = scnrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 455 "ctgsna.f"
	    cgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 + 1]
		    , &c__1, &c_b20, &work[1], &c__1, (ftnlen)1);
#line 457 "ctgsna.f"
	    cdotc_(&z__1, n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
#line 457 "ctgsna.f"
	    yhax.r = z__1.r, yhax.i = z__1.i;
#line 458 "ctgsna.f"
	    cgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 + 1]
		    , &c__1, &c_b20, &work[1], &c__1, (ftnlen)1);
#line 460 "ctgsna.f"
	    cdotc_(&z__1, n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
#line 460 "ctgsna.f"
	    yhbx.r = z__1.r, yhbx.i = z__1.i;
#line 461 "ctgsna.f"
	    d__1 = z_abs(&yhax);
#line 461 "ctgsna.f"
	    d__2 = z_abs(&yhbx);
#line 461 "ctgsna.f"
	    cond = slapy2_(&d__1, &d__2);
#line 462 "ctgsna.f"
	    if (cond == 0.) {
#line 463 "ctgsna.f"
		s[ks] = -1.;
#line 464 "ctgsna.f"
	    } else {
#line 465 "ctgsna.f"
		s[ks] = cond / (rnrm * lnrm);
#line 466 "ctgsna.f"
	    }
#line 467 "ctgsna.f"
	}

#line 469 "ctgsna.f"
	if (wantdf) {
#line 470 "ctgsna.f"
	    if (*n == 1) {
#line 471 "ctgsna.f"
		d__1 = z_abs(&a[a_dim1 + 1]);
#line 471 "ctgsna.f"
		d__2 = z_abs(&b[b_dim1 + 1]);
#line 471 "ctgsna.f"
		dif[ks] = slapy2_(&d__1, &d__2);
#line 472 "ctgsna.f"
	    } else {

/*              Estimate the reciprocal condition number of the k-th */
/*              eigenvectors. */

/*              Copy the matrix (A, B) to the array WORK and move the */
/*              (k,k)th pair to the (1,1) position. */

#line 480 "ctgsna.f"
		clacpy_("Full", n, n, &a[a_offset], lda, &work[1], n, (ftnlen)
			4);
#line 481 "ctgsna.f"
		clacpy_("Full", n, n, &b[b_offset], ldb, &work[*n * *n + 1], 
			n, (ftnlen)4);
#line 482 "ctgsna.f"
		ifst = k;
#line 483 "ctgsna.f"
		ilst = 1;

#line 485 "ctgsna.f"
		ctgexc_(&c_false, &c_false, n, &work[1], n, &work[*n * *n + 1]
			, n, dummy, &c__1, dummy1, &c__1, &ifst, &ilst, &ierr)
			;

#line 488 "ctgsna.f"
		if (ierr > 0) {

/*                 Ill-conditioned problem - swap rejected. */

#line 492 "ctgsna.f"
		    dif[ks] = 0.;
#line 493 "ctgsna.f"
		} else {

/*                 Reordering successful, solve generalized Sylvester */
/*                 equation for R and L, */
/*                            A22 * R - L * A11 = A12 */
/*                            B22 * R - L * B11 = B12, */
/*                 and compute estimate of Difl[(A11,B11), (A22, B22)]. */

#line 501 "ctgsna.f"
		    n1 = 1;
#line 502 "ctgsna.f"
		    n2 = *n - n1;
#line 503 "ctgsna.f"
		    i__ = *n * *n + 1;
#line 504 "ctgsna.f"
		    ctgsyl_("N", &c__3, &n2, &n1, &work[*n * n1 + n1 + 1], n, 
			    &work[1], n, &work[n1 + 1], n, &work[*n * n1 + n1 
			    + i__], n, &work[i__], n, &work[n1 + i__], n, &
			    scale, &dif[ks], dummy, &c__1, &iwork[1], &ierr, (
			    ftnlen)1);
#line 509 "ctgsna.f"
		}
#line 510 "ctgsna.f"
	    }
#line 511 "ctgsna.f"
	}

#line 513 "ctgsna.f"
L20:
#line 513 "ctgsna.f"
	;
#line 513 "ctgsna.f"
    }
#line 514 "ctgsna.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 515 "ctgsna.f"
    return 0;

/*     End of CTGSNA */

} /* ctgsna_ */

