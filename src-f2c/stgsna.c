#line 1 "stgsna.f"
/* stgsna.f -- translated by f2c (version 20100827).
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

#line 1 "stgsna.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = 1.;
static doublereal c_b21 = 0.;
static integer c__2 = 2;
static logical c_false = FALSE_;
static integer c__3 = 3;

/* > \brief \b STGSNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGSNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, */
/*                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          HOWMNY, JOB */
/*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), DIF( * ), S( * ), */
/*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSNA estimates reciprocal condition numbers for specified */
/* > eigenvalues and/or eigenvectors of a matrix pair (A, B) in */
/* > generalized real Schur canonical form (or of any matrix pair */
/* > (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where */
/* > Z**T denotes the transpose of Z. */
/* > */
/* > (A, B) must be in generalized real Schur form (as returned by SGGES), */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          B is REAL array, dimension (LDB,N) */
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
/* >          VL is REAL array, dimension (LDVL,M) */
/* >          If JOB = 'E' or 'B', VL must contain left eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT. The eigenvectors must be stored in consecutive */
/* >          columns of VL, as returned by STGEVC. */
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
/* >          VR is REAL array, dimension (LDVR,M) */
/* >          If JOB = 'E' or 'B', VR must contain right eigenvectors of */
/* >          (A, B), corresponding to the eigenpairs specified by HOWMNY */
/* >          and SELECT. The eigenvectors must be stored in consecutive */
/* >          columns ov VR, as returned by STGEVC. */
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
/* >          S is REAL array, dimension (MM) */
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
/* >          DIF is REAL array, dimension (MM) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup realOTHERcomputational */

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
/* >     (see SLATDF), then the parameter DIFDRI (see below) should be */
/* >     changed from 3 to 4 (routine SLATDF(IJOB = 2 will be used)). */
/* >     See STGSYL for more details. */
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
/* >     SLATDF), then the parameter DIFDRI (see below) should be changed */
/* >     from 3 to 4 (routine SLATDF(IJOB = 2 will be used)). See STGSYL */
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
/* Subroutine */ int stgsna_(char *job, char *howmny, logical *select, 
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
    static logical pair;
    static integer ierr;
    static doublereal uhav, uhbv;
    static integer ifst;
    static doublereal lnrm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ilst;
    static doublereal rnrm;
    extern /* Subroutine */ int slag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal root1, root2, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal uhavi, uhbvi;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal tmpii;
    static integer lwmin;
    static logical wants;
    static doublereal tmpir, tmpri, dummy[1], tmprr;
    extern doublereal slapy2_(doublereal *, doublereal *);
    static doublereal dummy1[1], alphai, alphar;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical wantbh, wantdf;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    stgexc_(logical *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *);
    static logical somcon;
    static doublereal alprqt, smlnum;
    static logical lquery;
    extern /* Subroutine */ int stgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);


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

#line 435 "stgsna.f"
    /* Parameter adjustments */
#line 435 "stgsna.f"
    --select;
#line 435 "stgsna.f"
    a_dim1 = *lda;
#line 435 "stgsna.f"
    a_offset = 1 + a_dim1;
#line 435 "stgsna.f"
    a -= a_offset;
#line 435 "stgsna.f"
    b_dim1 = *ldb;
#line 435 "stgsna.f"
    b_offset = 1 + b_dim1;
#line 435 "stgsna.f"
    b -= b_offset;
#line 435 "stgsna.f"
    vl_dim1 = *ldvl;
#line 435 "stgsna.f"
    vl_offset = 1 + vl_dim1;
#line 435 "stgsna.f"
    vl -= vl_offset;
#line 435 "stgsna.f"
    vr_dim1 = *ldvr;
#line 435 "stgsna.f"
    vr_offset = 1 + vr_dim1;
#line 435 "stgsna.f"
    vr -= vr_offset;
#line 435 "stgsna.f"
    --s;
#line 435 "stgsna.f"
    --dif;
#line 435 "stgsna.f"
    --work;
#line 435 "stgsna.f"
    --iwork;
#line 435 "stgsna.f"

#line 435 "stgsna.f"
    /* Function Body */
#line 435 "stgsna.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 436 "stgsna.f"
    wants = lsame_(job, "E", (ftnlen)1, (ftnlen)1) || wantbh;
#line 437 "stgsna.f"
    wantdf = lsame_(job, "V", (ftnlen)1, (ftnlen)1) || wantbh;

#line 439 "stgsna.f"
    somcon = lsame_(howmny, "S", (ftnlen)1, (ftnlen)1);

#line 441 "stgsna.f"
    *info = 0;
#line 442 "stgsna.f"
    lquery = *lwork == -1;

#line 444 "stgsna.f"
    if (! wants && ! wantdf) {
#line 445 "stgsna.f"
	*info = -1;
#line 446 "stgsna.f"
    } else if (! lsame_(howmny, "A", (ftnlen)1, (ftnlen)1) && ! somcon) {
#line 447 "stgsna.f"
	*info = -2;
#line 448 "stgsna.f"
    } else if (*n < 0) {
#line 449 "stgsna.f"
	*info = -4;
#line 450 "stgsna.f"
    } else if (*lda < max(1,*n)) {
#line 451 "stgsna.f"
	*info = -6;
#line 452 "stgsna.f"
    } else if (*ldb < max(1,*n)) {
#line 453 "stgsna.f"
	*info = -8;
#line 454 "stgsna.f"
    } else if (wants && *ldvl < *n) {
#line 455 "stgsna.f"
	*info = -10;
#line 456 "stgsna.f"
    } else if (wants && *ldvr < *n) {
#line 457 "stgsna.f"
	*info = -12;
#line 458 "stgsna.f"
    } else {

/*        Set M to the number of eigenpairs for which condition numbers */
/*        are required, and test MM. */

#line 463 "stgsna.f"
	if (somcon) {
#line 464 "stgsna.f"
	    *m = 0;
#line 465 "stgsna.f"
	    pair = FALSE_;
#line 466 "stgsna.f"
	    i__1 = *n;
#line 466 "stgsna.f"
	    for (k = 1; k <= i__1; ++k) {
#line 467 "stgsna.f"
		if (pair) {
#line 468 "stgsna.f"
		    pair = FALSE_;
#line 469 "stgsna.f"
		} else {
#line 470 "stgsna.f"
		    if (k < *n) {
#line 471 "stgsna.f"
			if (a[k + 1 + k * a_dim1] == 0.) {
#line 472 "stgsna.f"
			    if (select[k]) {
#line 472 "stgsna.f"
				++(*m);
#line 472 "stgsna.f"
			    }
#line 474 "stgsna.f"
			} else {
#line 475 "stgsna.f"
			    pair = TRUE_;
#line 476 "stgsna.f"
			    if (select[k] || select[k + 1]) {
#line 476 "stgsna.f"
				*m += 2;
#line 476 "stgsna.f"
			    }
#line 478 "stgsna.f"
			}
#line 479 "stgsna.f"
		    } else {
#line 480 "stgsna.f"
			if (select[*n]) {
#line 480 "stgsna.f"
			    ++(*m);
#line 480 "stgsna.f"
			}
#line 482 "stgsna.f"
		    }
#line 483 "stgsna.f"
		}
#line 484 "stgsna.f"
/* L10: */
#line 484 "stgsna.f"
	    }
#line 485 "stgsna.f"
	} else {
#line 486 "stgsna.f"
	    *m = *n;
#line 487 "stgsna.f"
	}

#line 489 "stgsna.f"
	if (*n == 0) {
#line 490 "stgsna.f"
	    lwmin = 1;
#line 491 "stgsna.f"
	} else if (lsame_(job, "V", (ftnlen)1, (ftnlen)1) || lsame_(job, 
		"B", (ftnlen)1, (ftnlen)1)) {
#line 492 "stgsna.f"
	    lwmin = (*n << 1) * (*n + 2) + 16;
#line 493 "stgsna.f"
	} else {
#line 494 "stgsna.f"
	    lwmin = *n;
#line 495 "stgsna.f"
	}
#line 496 "stgsna.f"
	work[1] = (doublereal) lwmin;

#line 498 "stgsna.f"
	if (*mm < *m) {
#line 499 "stgsna.f"
	    *info = -15;
#line 500 "stgsna.f"
	} else if (*lwork < lwmin && ! lquery) {
#line 501 "stgsna.f"
	    *info = -18;
#line 502 "stgsna.f"
	}
#line 503 "stgsna.f"
    }

#line 505 "stgsna.f"
    if (*info != 0) {
#line 506 "stgsna.f"
	i__1 = -(*info);
#line 506 "stgsna.f"
	xerbla_("STGSNA", &i__1, (ftnlen)6);
#line 507 "stgsna.f"
	return 0;
#line 508 "stgsna.f"
    } else if (lquery) {
#line 509 "stgsna.f"
	return 0;
#line 510 "stgsna.f"
    }

/*     Quick return if possible */

#line 514 "stgsna.f"
    if (*n == 0) {
#line 514 "stgsna.f"
	return 0;
#line 514 "stgsna.f"
    }

/*     Get machine constants */

#line 519 "stgsna.f"
    eps = slamch_("P", (ftnlen)1);
#line 520 "stgsna.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 521 "stgsna.f"
    ks = 0;
#line 522 "stgsna.f"
    pair = FALSE_;

#line 524 "stgsna.f"
    i__1 = *n;
#line 524 "stgsna.f"
    for (k = 1; k <= i__1; ++k) {

/*        Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block. */

#line 528 "stgsna.f"
	if (pair) {
#line 529 "stgsna.f"
	    pair = FALSE_;
#line 530 "stgsna.f"
	    goto L20;
#line 531 "stgsna.f"
	} else {
#line 532 "stgsna.f"
	    if (k < *n) {
#line 532 "stgsna.f"
		pair = a[k + 1 + k * a_dim1] != 0.;
#line 532 "stgsna.f"
	    }
#line 534 "stgsna.f"
	}

/*        Determine whether condition numbers are required for the k-th */
/*        eigenpair. */

#line 539 "stgsna.f"
	if (somcon) {
#line 540 "stgsna.f"
	    if (pair) {
#line 541 "stgsna.f"
		if (! select[k] && ! select[k + 1]) {
#line 541 "stgsna.f"
		    goto L20;
#line 541 "stgsna.f"
		}
#line 543 "stgsna.f"
	    } else {
#line 544 "stgsna.f"
		if (! select[k]) {
#line 544 "stgsna.f"
		    goto L20;
#line 544 "stgsna.f"
		}
#line 546 "stgsna.f"
	    }
#line 547 "stgsna.f"
	}

#line 549 "stgsna.f"
	++ks;

#line 551 "stgsna.f"
	if (wants) {

/*           Compute the reciprocal condition number of the k-th */
/*           eigenvalue. */

#line 556 "stgsna.f"
	    if (pair) {

/*              Complex eigenvalue pair. */

#line 560 "stgsna.f"
		d__1 = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 560 "stgsna.f"
		d__2 = snrm2_(n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
#line 560 "stgsna.f"
		rnrm = slapy2_(&d__1, &d__2);
#line 562 "stgsna.f"
		d__1 = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 562 "stgsna.f"
		d__2 = snrm2_(n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
#line 562 "stgsna.f"
		lnrm = slapy2_(&d__1, &d__2);
#line 564 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 566 "stgsna.f"
		tmprr = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 567 "stgsna.f"
		tmpri = sdot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 568 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)
			1);
#line 570 "stgsna.f"
		tmpii = sdot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 571 "stgsna.f"
		tmpir = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 572 "stgsna.f"
		uhav = tmprr + tmpii;
#line 573 "stgsna.f"
		uhavi = tmpir - tmpri;
#line 574 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 576 "stgsna.f"
		tmprr = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 577 "stgsna.f"
		tmpri = sdot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 578 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[(ks + 1) * 
			vr_dim1 + 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)
			1);
#line 580 "stgsna.f"
		tmpii = sdot_(n, &work[1], &c__1, &vl[(ks + 1) * vl_dim1 + 1],
			 &c__1);
#line 581 "stgsna.f"
		tmpir = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &
			c__1);
#line 582 "stgsna.f"
		uhbv = tmprr + tmpii;
#line 583 "stgsna.f"
		uhbvi = tmpir - tmpri;
#line 584 "stgsna.f"
		uhav = slapy2_(&uhav, &uhavi);
#line 585 "stgsna.f"
		uhbv = slapy2_(&uhbv, &uhbvi);
#line 586 "stgsna.f"
		cond = slapy2_(&uhav, &uhbv);
#line 587 "stgsna.f"
		s[ks] = cond / (rnrm * lnrm);
#line 588 "stgsna.f"
		s[ks + 1] = s[ks];

#line 590 "stgsna.f"
	    } else {

/*              Real eigenvalue. */

#line 594 "stgsna.f"
		rnrm = snrm2_(n, &vr[ks * vr_dim1 + 1], &c__1);
#line 595 "stgsna.f"
		lnrm = snrm2_(n, &vl[ks * vl_dim1 + 1], &c__1);
#line 596 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &a[a_offset], lda, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 598 "stgsna.f"
		uhav = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1)
			;
#line 599 "stgsna.f"
		sgemv_("N", n, n, &c_b19, &b[b_offset], ldb, &vr[ks * vr_dim1 
			+ 1], &c__1, &c_b21, &work[1], &c__1, (ftnlen)1);
#line 601 "stgsna.f"
		uhbv = sdot_(n, &work[1], &c__1, &vl[ks * vl_dim1 + 1], &c__1)
			;
#line 602 "stgsna.f"
		cond = slapy2_(&uhav, &uhbv);
#line 603 "stgsna.f"
		if (cond == 0.) {
#line 604 "stgsna.f"
		    s[ks] = -1.;
#line 605 "stgsna.f"
		} else {
#line 606 "stgsna.f"
		    s[ks] = cond / (rnrm * lnrm);
#line 607 "stgsna.f"
		}
#line 608 "stgsna.f"
	    }
#line 609 "stgsna.f"
	}

#line 611 "stgsna.f"
	if (wantdf) {
#line 612 "stgsna.f"
	    if (*n == 1) {
#line 613 "stgsna.f"
		dif[ks] = slapy2_(&a[a_dim1 + 1], &b[b_dim1 + 1]);
#line 614 "stgsna.f"
		goto L20;
#line 615 "stgsna.f"
	    }

/*           Estimate the reciprocal condition number of the k-th */
/*           eigenvectors. */
#line 619 "stgsna.f"
	    if (pair) {

/*              Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)). */
/*              Compute the eigenvalue(s) at position K. */

#line 624 "stgsna.f"
		work[1] = a[k + k * a_dim1];
#line 625 "stgsna.f"
		work[2] = a[k + 1 + k * a_dim1];
#line 626 "stgsna.f"
		work[3] = a[k + (k + 1) * a_dim1];
#line 627 "stgsna.f"
		work[4] = a[k + 1 + (k + 1) * a_dim1];
#line 628 "stgsna.f"
		work[5] = b[k + k * b_dim1];
#line 629 "stgsna.f"
		work[6] = b[k + 1 + k * b_dim1];
#line 630 "stgsna.f"
		work[7] = b[k + (k + 1) * b_dim1];
#line 631 "stgsna.f"
		work[8] = b[k + 1 + (k + 1) * b_dim1];
#line 632 "stgsna.f"
		d__1 = smlnum * eps;
#line 632 "stgsna.f"
		slag2_(&work[1], &c__2, &work[5], &c__2, &d__1, &beta, dummy1,
			 &alphar, dummy, &alphai);
#line 634 "stgsna.f"
		alprqt = 1.;
#line 635 "stgsna.f"
		c1 = (alphar * alphar + alphai * alphai + beta * beta) * 2.;
#line 636 "stgsna.f"
		c2 = beta * 4. * beta * alphai * alphai;
#line 637 "stgsna.f"
		root1 = c1 + sqrt(c1 * c1 - c2 * 4.);
#line 638 "stgsna.f"
		root2 = c2 / root1;
#line 639 "stgsna.f"
		root1 /= 2.;
/* Computing MIN */
#line 640 "stgsna.f"
		d__1 = sqrt(root1), d__2 = sqrt(root2);
#line 640 "stgsna.f"
		cond = min(d__1,d__2);
#line 641 "stgsna.f"
	    }

/*           Copy the matrix (A, B) to the array WORK and swap the */
/*           diagonal block beginning at A(k,k) to the (1,1) position. */

#line 646 "stgsna.f"
	    slacpy_("Full", n, n, &a[a_offset], lda, &work[1], n, (ftnlen)4);
#line 647 "stgsna.f"
	    slacpy_("Full", n, n, &b[b_offset], ldb, &work[*n * *n + 1], n, (
		    ftnlen)4);
#line 648 "stgsna.f"
	    ifst = k;
#line 649 "stgsna.f"
	    ilst = 1;

#line 651 "stgsna.f"
	    i__2 = *lwork - (*n << 1) * *n;
#line 651 "stgsna.f"
	    stgexc_(&c_false, &c_false, n, &work[1], n, &work[*n * *n + 1], n,
		     dummy, &c__1, dummy1, &c__1, &ifst, &ilst, &work[(*n * *
		    n << 1) + 1], &i__2, &ierr);

#line 655 "stgsna.f"
	    if (ierr > 0) {

/*              Ill-conditioned problem - swap rejected. */

#line 659 "stgsna.f"
		dif[ks] = 0.;
#line 660 "stgsna.f"
	    } else {

/*              Reordering successful, solve generalized Sylvester */
/*              equation for R and L, */
/*                         A22 * R - L * A11 = A12 */
/*                         B22 * R - L * B11 = B12, */
/*              and compute estimate of Difl((A11,B11), (A22, B22)). */

#line 668 "stgsna.f"
		n1 = 1;
#line 669 "stgsna.f"
		if (work[2] != 0.) {
#line 669 "stgsna.f"
		    n1 = 2;
#line 669 "stgsna.f"
		}
#line 671 "stgsna.f"
		n2 = *n - n1;
#line 672 "stgsna.f"
		if (n2 == 0) {
#line 673 "stgsna.f"
		    dif[ks] = cond;
#line 674 "stgsna.f"
		} else {
#line 675 "stgsna.f"
		    i__ = *n * *n + 1;
#line 676 "stgsna.f"
		    iz = (*n << 1) * *n + 1;
#line 677 "stgsna.f"
		    i__2 = *lwork - (*n << 1) * *n;
#line 677 "stgsna.f"
		    stgsyl_("N", &c__3, &n2, &n1, &work[*n * n1 + n1 + 1], n, 
			    &work[1], n, &work[n1 + 1], n, &work[*n * n1 + n1 
			    + i__], n, &work[i__], n, &work[n1 + i__], n, &
			    scale, &dif[ks], &work[iz + 1], &i__2, &iwork[1], 
			    &ierr, (ftnlen)1);

#line 683 "stgsna.f"
		    if (pair) {
/* Computing MIN */
#line 683 "stgsna.f"
			d__1 = max(1.,alprqt) * dif[ks];
#line 683 "stgsna.f"
			dif[ks] = min(d__1,cond);
#line 683 "stgsna.f"
		    }
#line 686 "stgsna.f"
		}
#line 687 "stgsna.f"
	    }
#line 688 "stgsna.f"
	    if (pair) {
#line 688 "stgsna.f"
		dif[ks + 1] = dif[ks];
#line 688 "stgsna.f"
	    }
#line 690 "stgsna.f"
	}
#line 691 "stgsna.f"
	if (pair) {
#line 691 "stgsna.f"
	    ++ks;
#line 691 "stgsna.f"
	}

#line 694 "stgsna.f"
L20:
#line 694 "stgsna.f"
	;
#line 694 "stgsna.f"
    }
#line 695 "stgsna.f"
    work[1] = (doublereal) lwmin;
#line 696 "stgsna.f"
    return 0;

/*     End of STGSNA */

} /* stgsna_ */

