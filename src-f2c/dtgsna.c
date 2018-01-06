#line 1 "dtgsna.f"
/* dtgsna.f -- translated by f2c (version 20100827).
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

#line 1 "dtgsna.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = 1.;
static doublereal c_b21 = 0.;
static integer c__2 = 2;
static logical c_false = FALSE_;
static integer c__3 = 3;

/* > \brief \b DTGSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, */
/*                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), DIF( * ), S( * ), */
/*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or eigenvectors of a matrix pair (A, B) in */
/* > generalized real Schur canonical form (or of any matrix pair */
/* > (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where */
/* > Z**T denotes the transpose of Z. */
/* > */
/* > (A, B) must be in generalized real Schur form (as returned by DGGES), */
/* > i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal */
/* > blocks. B is upper triangular. */
/* > */
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
/* >          for the eigenpair corresponding to a real eigenvalue w(j), */
/* >          SELECT(j) must be set to .TRUE.. To select condition numbers */
/* >          corresponding to a complex conjugate pair of eigenvalues w(j) */
/* >          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be */
/* >          set to .TRUE.. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The upper quasi-triangular matrix A in the pair (A,B). */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          The upper triangular matrix B in the pair (A,B). */
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
/* >          VL is DOUBLE PRECISION array, dimension (LDVL,M) */
/* >          If JOB = 'E' or 'B', VL must contain left eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT. The eigenvectors must be stored in consecutive */
/* >          columns of VL, as returned by DTGEVC. */
/* >          If JOB = 'V', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* >          LDVL is INTEGER */
/* >          The leading dimension of the array VL. LDVL >= 1. */
/* >          If JOB = 'E' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] VR */
/* > \verbatim */
/* >          VR is DOUBLE PRECISION array, dimension (LDVR,M) */
/* >          If JOB = 'E' or 'B', VR must contain right eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT. The eigenvectors must be stored in consecutive */
/* >          columns ov VR, as returned by DTGEVC. */
/* >          If JOB = 'V', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* >          LDVR is INTEGER */
/* >          The leading dimension of the array VR. LDVR >= 1. */
/* >          If JOB = 'E' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
/* >          selected eigenvalues, stored in consecutive elements of the */
/* >          array. For a complex conjugate pair of eigenvalues two */
/* >          consecutive elements of S are set to the same value. Thus */
/* >          S(j), DIF(j), and the j-th columns of VL and VR all */
/* >          correspond to the same eigenpair (but not in general the */
/* >          j-th eigenpair, unless all eigenpairs are selected). */
/* >          If JOB = 'V', S is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is DOUBLE PRECISION array, dimension (MM) */
/* >          If JOB = 'V' or 'B', the estimated reciprocal condition */
/* >          numbers of the selected eigenvectors, stored in consecutive */
/* >          elements of the array. For a complex eigenvector two */
/* >          consecutive elements of DIF are set to the same value. If */
/* >          the eigenvalues cannot be reordered to compute DIF(j), DIF(j) */
/* >          is set to 0; this can only occur when the true value would be */
/* >          very small anyway. */
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
/* >          the specified condition numbers; for each selected real */
/* >          eigenvalue one element is used, and for each selected complex */
/* >          conjugate pair of eigenvalues, two elements are used. */
/* >          If HOWMNY = 'A', M is set to N. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,N). */
/* >          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N + 6) */
/* >          If JOB = 'E', IWORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          =0: Successful exit */
/* >          <0: If INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The reciprocal of the condition number of a generalized eigenvalue */
/* >  w = (a, b) is defined as */
/* > */
/* >       S(w) = (|u**TAv|**2 + |u**TBv|**2)**(1/2) / (norm(u)*norm(v)) */
/* > */
/* >  where u and v are the left and right eigenvectors of (A, B) */
/* >  corresponding to w; |z| denotes the absolute value of the complex */
/* >  number, and norm(u) denotes the 2-norm of the vector u. */
/* >  The pair (a, b) corresponds to an eigenvalue w = a/b (= u**TAv/u**TBv) */
/* >  of the matrix pair (A, B). If both a and b equal zero, then (A B) is */
/* >  singular and S(I) = -1 is returned. */
/* > */
/* >  An approximate error bound on the chordal distance between the i-th */
/* >  computed generalized eigenvalue w and the corresponding exact */
/* >  eigenvalue lambda is */
/* > */
/* >       chord(w, lambda) <= EPS * norm(A, B) / S(I) */
/* > */
/* >  where EPS is the machine precision. */
/* > */
/* >  The reciprocal of the condition number DIF(i) of right eigenvector u */
/* >  and left eigenvector v corresponding to the generalized eigenvalue w */
/* >  is defined as follows: */
/* > */
/* >  a) If the i-th eigenvalue w = (a,b) is real */
/* > */
/* >     Suppose U and V are orthogonal transformations such that */
/* > */
/* >              U**T*(A, B)*V  = (S, T) = ( a   *  ) ( b  *  )  1 */
/* >                                        ( 0  S22 ),( 0 T22 )  n-1 */
/* >                                          1  n-1     1 n-1 */
/* > */
/* >     Then the reciprocal condition number DIF(i) is */
/* > */
/* >                Difl((a, b), (S22, T22)) = sigma-min( Zl ), */
/* > */
/* >     where sigma-min(Zl) denotes the smallest singular value of the */
/* >     2(n-1)-by-2(n-1) matrix */
/* > */
/* >         Zl = [ kron(a, In-1)  -kron(1, S22) ] */
/* >              [ kron(b, In-1)  -kron(1, T22) ] . */
/* > */
/* >     Here In-1 is the identity matrix of size n-1. kron(X, Y) is the */
/* >     Kronecker product between the matrices X and Y. */
/* > */
/* >     Note that if the default method for computing DIF(i) is wanted */
/* >     (see DLATDF), then the parameter DIFDRI (see below) should be */
/* >     changed from 3 to 4 (routine DLATDF(IJOB = 2 will be used)). */
/* >     See DTGSYL for more details. */
/* > */
/* >  b) If the i-th and (i+1)-th eigenvalues are complex conjugate pair, */
/* > */
/* >     Suppose U and V are orthogonal transformations such that */
/* > */
/* >              U**T*(A, B)*V = (S, T) = ( S11  *   ) ( T11  *  )  2 */
/* >                                       ( 0    S22 ),( 0    T22) n-2 */
/* >                                         2    n-2     2    n-2 */
/* > */
/* >     and (S11, T11) corresponds to the complex conjugate eigenvalue */
/* >     pair (w, conjg(w)). There exist unitary matrices U1 and V1 such */
/* >     that */
/* > */
/* >       U1**T*S11*V1 = ( s11 s12 ) and U1**T*T11*V1 = ( t11 t12 ) */
/* >                      (  0  s22 )                    (  0  t22 ) */
/* > */
/* >     where the generalized eigenvalues w = s11/t11 and */
/* >     conjg(w) = s22/t22. */
/* > */
/* >     Then the reciprocal condition number DIF(i) is bounded by */
/* > */
/* >         min( d1, max( 1, |real(s11)/real(s22)| )*d2 ) */
/* > */
/* >     where, d1 = Difl((s11, t11), (s22, t22)) = sigma-min(Z1), where */
/* >     Z1 is the complex 2-by-2 matrix */
/* > */
/* >              Z1 =  [ s11  -s22 ] */
/* >                    [ t11  -t22 ], */
/* > */
/* >     This is done by computing (using real arithmetic) the */
/* >     roots of the characteristical polynomial det(Z1**T * Z1 - lambda I), */
/* >     where Z1**T denotes the transpose of Z1 and det(X) denotes */
/* >     the determinant of X. */
/* > */
/* >     and d2 is an upper bound on Difl((S11, T11), (S22, T22)), i.e. an */
/* >     upper bound on sigma-min(Z2), where Z2 is (2n-2)-by-(2n-2) */
/* > */
/* >              Z2 = [ kron(S11**T, In-2)  -kron(I2, S22) ] */
/* >                   [ kron(T11**T, In-2)  -kron(I2, T22) ] */
/* > */
/* >     Note that if the default method for computing DIF is wanted (see */
/* >     DLATDF), then the parameter DIFDRI (see below) should be changed */
/* >     from 3 to 4 (routine DLATDF(IJOB = 2 will be used)). See DTGSYL */
/* >     for more details. */
/* > */
/* >  For each eigenvalue/vector specified by SELECT, DIF stores a */
/* >  Frobenius norm-based estimate of Difl. */
/* > */
/* >  An approximate error bound for the i-th computed eigenvector VL(i) or */
/* >  VR(i) is given by */
/* > */
/* >             EPS * norm(A, B) / DIF(i). */
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
/* >      Estimation: Theory, Algorithms and Software, */
/* >      Report UMINF - 94.04, Department of Computing Science, Umea */
/* >      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working */
/* >      Note 87. To appear in Numerical Algorithms, 1996. */
/* > */
/* >  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* >      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22, */
/* >      No 1, 1996. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtgsna_(char *job, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
	work, integer *lwork, integer *iwork, integer *info, ftnlen job_len, 
	ftnlen howmny_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal c1, c2;
    static integer n1, n2, ks, iz;
    static doublereal eps, beta, cond;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical pair;
    static integer ierr;
    static doublereal uhav, uhbv;
    static integer ifst;
    static doublereal lnrm;
    static integer ilst;
    static doublereal rnrm;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal root1, root2, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal uhavi, uhbvi, tmpii;
    static integer lwmin;
    static logical wants;
    static doublereal tmpir, tmpri, dummy[1], tmprr;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal dummy1[1];
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal alphai, alphar;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtgexc_(logical *, logical *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *);
    static logical wantbh, wantdf, somcon;
    static doublereal alprqt;
    extern /* Subroutine */ int dtgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
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

#line 435 "dtgsna.f"
    /* Parameter adjustments */
#line 435 "dtgsna.f"
    --select;
#line 435 "dtgsna.f"
    a_dim1 = *lda;
#line 435 "dtgsna.f"
    a_offset = 1 + a_dim1;
#line 435 "dtgsna.f"
    a -= a_offset;
#line 435 "dtgsna.f"
    b_dim1 = *ldb;
#line 435 "dtgsna.f"
    b_offset = 1 + b_dim1;
#line 435 "dtgsna.f"
    b -= b_offset;
#line 435 "dtgsna.f"
    vl_dim1 = *ldvl;
#line 435 "dtgsna.f"
    vl_offset = 1 + vl_dim1;
#line 435 "dtgsna.f"
    vl -= vl_offset;
#line 435 "dtgsna.f"
    vr_dim1 = *ldvr;
#line 435 "dtgsna.f"
    vr_offset = 1 + vr_dim1;
#line 435 "dtgsna.f"
    vr -= vr_offset;
#line 435 "dtgsna.f"
    --s;
#line 435 "dtgsna.f"
    --dif;
#line 435 "dtgsna.f"
    --work;
#line 435 "dtgsna.f"
    --iwork;
#line 435 "dtgsna.f"

#line 435 "dtgsna.f"
    /* Function Body */
#line 435 "dtgsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 436 "dtgsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 437 "dtgsna.f"
    wantdf = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 439 "dtgsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 441 "dtgsna.f"
    *info = 0;
#line 442 "dtgsna.f"
    lquery = *lwork == -1;

#line 444 "dtgsna.f"
    if (! wants && ! wantdf) {
#line 445 "dtgsna.f"
	*info = -1;
#line 446 "dtgsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 447 "dtgsna.f"
	*info = -2;
#line 448 "dtgsna.f"
    } else if (*n < 0) {
#line 449 "dtgsna.f"
	*info = -4;
#line 450 "dtgsna.f"
    } else if (*lda < max(1,*n)) {
#line 451 "dtgsna.f"
	*info = -6;
#line 452 "dtgsna.f"
    } else if (*ldb < max(1,*n)) {
#line 453 "dtgsna.f"
	*info = -8;
#line 454 "dtgsna.f"
    } else if (wants && *ldvl < *n) {
#line 455 "dtgsna.f"
	*info = -10;
#line 456 "dtgsna.f"
    } else if (wants && *ldvr < *n) {
#line 457 "dtgsna.f"
	*info = -12;
#line 458 "dtgsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 463 "dtgsna.f"
	if (somcon) {
#line 464 "dtgsna.f"
	    *m = 0;
#line 465 "dtgsna.f"
	    pair = FALSE_;
#line 466 "dtgsna.f"
	    i__1 = *n;
#line 466 "dtgsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 467 "dtgsna.f"
		if (pair) {
#line 468 "dtgsna.f"
		    pair = FALSE_;
#line 469 "dtgsna.f"
		} else {
#line 470 "dtgsna.f"
		    if (k < *n) {
#line 471 "dtgsna.f"
			if (a[k + 1 + k * a_dim1] == 0.) {
#line 472 "dtgsna.f"
			    if (select[k]) {
#line 472 "dtgsna.f"
				++(*m);
#line 472 "dtgsna.f"
			    }
#line 474 "dtgsna.f"
			} else {
#line 475 "dtgsna.f"
			    pair = TRUE_;
#line 476 "dtgsna.f"
			    if (select[k] || select[k + 1]) {
#line 476 "dtgsna.f"
				*m += 2;
#line 476 "dtgsna.f"
			    }
#line 478 "dtgsna.f"
			}
#line 479 "dtgsna.f"
		    } else {
#line 480 "dtgsna.f"
			if (select[*n]) {
#line 480 "dtgsna.f"
			    ++(*m);
#line 480 "dtgsna.f"
			}
#line 482 "dtgsna.f"
		    }
#line 483 "dtgsna.f"
		}
#line 484 "dtgsna.f"
/* L10: */
#line 484 "dtgsna.f"
	    }
#line 485 "dtgsna.f"
	} else {
#line 486 "dtgsna.f"
	    *m = *n;
#line 487 "dtgsna.f"
	}

#line 489 "dtgsna.f"
	if (*n == 0) {
#line 490 "dtgsna.f"
	    lwmin = 1;
#line 491 "dtgsna.f"
	} else if (lsame_(job, "V", (ftnlen)1, (ftnlen)1) || lsame_(job, 
		"B", (ftnlen)1, (ftnlen)1)) {
#line 492 "dtgsna.f"
	    lwmin = (*n << 1) * (*n + 2) + 16;
#line 493 "dtgsna.f"
	} else {
#line 494 "dtgsna.f"
	    lwmin = *n;
#line 495 "dtgsna.f"
	}
#line 496 "dtgsna.f"
	work[1] = (doublereal) lwmin;

#line 498 "dtgsna.f"
	if (*mm < *m) {
#line 499 "dtgsna.f"
	    *info = -15;
#line 500 "dtgsna.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 501 "dtgsna.f"
	    *info = -18;
#line 502 "dtgsna.f"
	}
#line 503 "dtgsna.f"
    }

#line 505 "dtgsna.f"
    if (*info != 0) {
#line 506 "dtgsna.f"
	i__1 = -(*info);
#line 506 "dtgsna.f"
	xerbla_("DTGSNA", &i__1, (ftnlen)6);
#line 507 "dtgsna.f"
	return 0;
#line 508 "dtgsna.f"
    } else if (lquery) {
#line 509 "dtgsna.f"
	return 0;
#line 510 "dtgsna.f"
    }

/*     Quick return if possible */

#line 514 "dtgsna.f"
    if (*n == 0) {
#line 514 "dtgsna.f"
	return 0;
#line 514 "dtgsna.f"
    }

/*     Get machine constants */

#line 519 "dtgsna.f"
    eps = dlamch_("P", (ftnlen)1);
#line 520 "dtgsna.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 521 "dtgsna.f"
    ks = 0;
#line 522 "dtgsna.f"
    pair = FALSE_;

#line 524 "dtgsna.f"
    i__1 = *n;
#line 524 "dtgsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block. */

#line 528 "dtgsna.f"
	if (pair) {
#line 529 "dtgsna.f"
	    pair = FALSE_;
#line 530 "dtgsna.f"
	    goto L20;
#line 531 "dtgsna.f"
	} else {
#line 532 "dtgsna.f"
	    if (k < *n) {
#line 532 "dtgsna.f"
		pair = a[k + 1 + k * a_dim1] != 0.;
#line 532 "dtgsna.f"
	    }
#line 534 "dtgsna.f"
	}

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 539 "dtgsna.f"
	if (somcon) {
#line 540 "dtgsna.f"
	    if (pair) {
#line 541 "dtgsna.f"
		if (! select[k] && ! select[k + 1]) {
#line 541 "dtgsna.f"
		    goto L20;
#line 541 "dtgsna.f"
		}
#line 543 "dtgsna.f"
	    } else {
#line 544 "dtgsna.f"
		if (! select[k]) {
#line 544 "dtgsna.f"
		    goto L20;
#line 544 "dtgsna.f"
		}
#line 546 "dtgsna.f"
	    }
#line 547 "dtgsna.f"
	}

#line 549 "dtgsna.f"
	++ks;

#line 551 "dtgsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 556 "dtgsna.f"
	    if (pair) {

/*              Complex eigenvalue pair. */

#line 560 "dtgsna.f"
		d__1 = dnrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 560 "dtgsna.f"
		d__2 = dnrm2_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
#line 560 "dtgsna.f"
		rnrm = dlapy2_(&d__1, &d__2);
#line 562 "dtgsna.f"
		d__1 = dnrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 562 "dtgsna.f"
		d__2 = dnrm2_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
#line 562 "dtgsna.f"
		lnrm = dlapy2_(&d__1, &d__2);
#line 564 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 566 "dtgsna.f"
		tmprr = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 567 "dtgsna.f"
		tmpri = ddot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 568 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)
			1);
#line 570 "dtgsna.f"
		tmpii = ddot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 571 "dtgsna.f"
		tmpir = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 572 "dtgsna.f"
		uhav = tmprr + tmpii;
#line 573 "dtgsna.f"
		uhavi = tmpir - tmpri;
#line 574 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 576 "dtgsna.f"
		tmprr = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 577 "dtgsna.f"
		tmpri = ddot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 578 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)
			1);
#line 580 "dtgsna.f"
		tmpii = ddot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 581 "dtgsna.f"
		tmpir = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 582 "dtgsna.f"
		uhbv = tmprr + tmpii;
#line 583 "dtgsna.f"
		uhbvi = tmpir - tmpri;
#line 584 "dtgsna.f"
		uhav = dlapy2_(&uhav, &uhavi);
#line 585 "dtgsna.f"
		uhbv = dlapy2_(&uhbv, &uhbvi);
#line 586 "dtgsna.f"
		cond = dlapy2_(&uhav, &uhbv);
#line 587 "dtgsna.f"
		s[ks] = cond / (rnrm * lnrm);
#line 588 "dtgsna.f"
		s[ks + 1] = s[ks];

#line 590 "dtgsna.f"
	    } else {

/*              Real eigenvalue. */

#line 594 "dtgsna.f"
		rnrm = dnrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 595 "dtgsna.f"
		lnrm = dnrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 596 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 598 "dtgsna.f"
		uhav = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1)
			;
#line 599 "dtgsna.f"
		dgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 601 "dtgsna.f"
		uhbv = ddot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1)
			;
#line 602 "dtgsna.f"
		cond = dlapy2_(&uhav, &uhbv);
#line 603 "dtgsna.f"
		if (cond == 0.) {
#line 604 "dtgsna.f"
		    s[ks] = -1.;
#line 605 "dtgsna.f"
		} else {
#line 606 "dtgsna.f"
		    s[ks] = cond / (rnrm * lnrm);
#line 607 "dtgsna.f"
		}
#line 608 "dtgsna.f"
	    }
#line 609 "dtgsna.f"
	}

#line 611 "dtgsna.f"
	if (wantdf) {
#line 612 "dtgsna.f"
	    if (*n == 1) {
#line 613 "dtgsna.f"
		dif[ks] = dlapy2_(&a[a_dim1 + 1], &b[b_dim1 + 1]);
#line 614 "dtgsna.f"
		goto L20;
#line 615 "dtgsna.f"
	    }

/*           Estimate the reciprocal condition number of the k-th */
/*           eigenvectors. */
#line 619 "dtgsna.f"
	    if (pair) {

/*              Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)). */
/*              Compute the eigenvalue(s) at position K. */

#line 624 "dtgsna.f"
		work[1] = a[k + k * a_dim1];
#line 625 "dtgsna.f"
		work[2] = a[k + 1 + k * a_dim1];
#line 626 "dtgsna.f"
		work[3] = a[k + (k + 1) * a_dim1];
#line 627 "dtgsna.f"
		work[4] = a[k + 1 + (k + 1) * a_dim1];
#line 628 "dtgsna.f"
		work[5] = b[k + k * b_dim1];
#line 629 "dtgsna.f"
		work[6] = b[k + 1 + k * b_dim1];
#line 630 "dtgsna.f"
		work[7] = b[k + (k + 1) * b_dim1];
#line 631 "dtgsna.f"
		work[8] = b[k + 1 + (k + 1) * b_dim1];
#line 632 "dtgsna.f"
		d__1 = smlnum * eps;
#line 632 "dtgsna.f"
		dlag2_(&work[1], &c__2, &work[5], &c__2, &d__1, &beta, dummy1,
			 &alphar, dummy, &alphai);
#line 634 "dtgsna.f"
		alprqt = 1.;
#line 635 "dtgsna.f"
		c1 = (alphar * alphar + alphai * alphai + beta * beta) * 2.;
#line 636 "dtgsna.f"
		c2 = beta * 4. * beta * alphai * alphai;
#line 637 "dtgsna.f"
		root1 = c1 + sqrt(c1 * c1 - c2 * 4.);
#line 638 "dtgsna.f"
		root2 = c2 / root1;
#line 639 "dtgsna.f"
		root1 /= 2.;
/* Computing MIN */
#line 640 "dtgsna.f"
		d__1 = sqrt(root1), d__2 = sqrt(root2);
#line 640 "dtgsna.f"
		cond = min(d__1,d__2);
#line 641 "dtgsna.f"
	    }

/*           Copy the matrix (A, B) to the array WORK and swap the */
/*           diagonal block beginning at A(k,k) to the (1,1) position. */

#line 646 "dtgsna.f"
	    dlacpy_("Full", n, n, &a[a_offset], lda, &work[1], n, (ftnlen)4);
#line 647 "dtgsna.f"
	    dlacpy_("Full", n, n, &b[b_offset], ldb, &work[*n * *n + 1], n, (
		    ftnlen)4);
#line 648 "dtgsna.f"
	    ifst = k;
#line 649 "dtgsna.f"
	    ilst = 1;

#line 651 "dtgsna.f"
	    i__2 = *lwork - (*n << 1) * *n;
#line 651 "dtgsna.f"
	    dtgexc_(&c_false, &c_false, n, &work[1], n, &work[*n * *n + 1], n,
		     dummy, &c__1, dummy1, &c__1, &ifst, &ilst, &work[(*n * *
		    n << 1) + 1], &i__2, &ierr);

#line 655 "dtgsna.f"
	    if (ierr > 0) {

/*              Ill-conditioned problem - swap rejected. */

#line 659 "dtgsna.f"
		dif[ks] = 0.;
#line 660 "dtgsna.f"
	    } else {

/*              Reordering successful, solve generalized Sylvester */
/*              equation for R and L, */
/*                         A22 * R - L * A11 = A12 */
/*                         B22 * R - L * B11 = B12, */
/*              and compute estimate of Difl((A11,B11), (A22, B22)). */

#line 668 "dtgsna.f"
		n1 = 1;
#line 669 "dtgsna.f"
		if (work[2] != 0.) {
#line 669 "dtgsna.f"
		    n1 = 2;
#line 669 "dtgsna.f"
		}
#line 671 "dtgsna.f"
		n2 = *n - n1;
#line 672 "dtgsna.f"
		if (n2 == 0) {
#line 673 "dtgsna.f"
		    dif[ks] = cond;
#line 674 "dtgsna.f"
		} else {
#line 675 "dtgsna.f"
		    i__ = *n * *n + 1;
#line 676 "dtgsna.f"
		    iz = (*n << 1) * *n + 1;
#line 677 "dtgsna.f"
		    i__2 = *lwork - (*n << 1) * *n;
#line 677 "dtgsna.f"
		    dtgsyl_("N", &c__3, &n2, &n1, &work[*n * n1 + n1 + 1], n, 
			    &work[1], n, &work[n1 + 1], n, &work[*n * n1 + n1 
			    + i__], n, &work[i__], n, &work[n1 + i__], n, &
			    scale, &dif[ks], &work[iz + 1], &i__2, &iwork[1], 
			    &ierr, (ftnlen)1);

#line 683 "dtgsna.f"
		    if (pair) {
/* Computing MIN */
#line 683 "dtgsna.f"
			d__1 = max(1.,alprqt) * dif[ks];
#line 683 "dtgsna.f"
			dif[ks] = min(d__1,cond);
#line 683 "dtgsna.f"
		    }
#line 686 "dtgsna.f"
		}
#line 687 "dtgsna.f"
	    }
#line 688 "dtgsna.f"
	    if (pair) {
#line 688 "dtgsna.f"
		dif[ks + 1] = dif[ks];
#line 688 "dtgsna.f"
	    }
#line 690 "dtgsna.f"
	}
#line 691 "dtgsna.f"
	if (pair) {
#line 691 "dtgsna.f"
	    ++ks;
#line 691 "dtgsna.f"
	}

#line 694 "dtgsna.f"
L20:
#line 694 "dtgsna.f"
	;
#line 694 "dtgsna.f"
    }
#line 695 "dtgsna.f"
    work[1] = (doublereal) lwmin;
#line 696 "dtgsna.f"
    return 0;

/*     End of DTGSNA */

} /* dtgsna_ */

