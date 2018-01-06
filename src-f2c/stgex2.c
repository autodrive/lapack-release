#line 1 "stgex2.f"
/* stgex2.f -- translated by f2c (version 20100827).
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

#line 1 "stgex2.f"
/* Table of constant values */

static integer c__4 = 4;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b42 = 1.;
static doublereal c_b48 = -1.;
static integer c__0 = 0;

/* > \brief \b STGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogon
al equivalence transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STGEX2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgex2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgex2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgex2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/*                          LDZ, J1, N1, N2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22) */
/* > of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair */
/* > (A, B) by an orthogonal equivalence transformation. */
/* > */
/* > (A, B) must be in generalized real Schur canonical form (as returned */
/* > by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2 */
/* > diagonal blocks. B is upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* >        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T */
/* >        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] WANTQ */
/* > \verbatim */
/* >          WANTQ is LOGICAL */
/* >          .TRUE. : update the left transformation matrix Q; */
/* >          .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          .TRUE. : update the right transformation matrix Z; */
/* >          .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the matrix A in the pair (A, B). */
/* >          On exit, the updated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the matrix B in the pair (A, B). */
/* >          On exit, the updated matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ,N) */
/* >          On entry, if WANTQ = .TRUE., the orthogonal matrix Q. */
/* >          On exit, the updated matrix Q. */
/* >          Not referenced if WANTQ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= 1. */
/* >          If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ,N) */
/* >          On entry, if WANTZ =.TRUE., the orthogonal matrix Z. */
/* >          On exit, the updated matrix Z. */
/* >          Not referenced if WANTZ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1. */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] J1 */
/* > \verbatim */
/* >          J1 is INTEGER */
/* >          The index to the first block (A11, B11). 1 <= J1 <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          The order of the first block (A11, B11). N1 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* >          The order of the second block (A22, B22). N2 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)). */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >=  MAX( N*(N2+N1), (N2+N1)*(N2+N1)*2 ) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >            =0: Successful exit */
/* >            >0: If INFO = 1, the transformed matrix (A, B) would be */
/* >                too far from generalized Schur form; the blocks are */
/* >                not swapped and (A, B) and (Q, Z) are unchanged. */
/* >                The problem of swapping is too ill-conditioned. */
/* >            <0: If INFO = -16: LWORK is too small. Appropriate value */
/* >                for LWORK is returned in WORK(1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup realGEauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  In the current code both weak and strong stability tests are */
/* >  performed. The user can omit the strong stability test by changing */
/* >  the internal logical parameter WANDS to .FALSE.. See ref. [2] for */
/* >  details. */

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
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stgex2_(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
	n1, integer *n2, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer i__, m;
    static doublereal s[16]	/* was [4][4] */, t[16]	/* was [4][4] */, be[
	    2], ai[2], ar[2], sa, sb, li[16]	/* was [4][4] */, ir[16]	
	    /* was [4][4] */, ss, ws, eps;
    static logical weak;
    static doublereal ddum;
    static integer idum;
    static doublereal taul[4], dsum, taur[4], scpy[16]	/* was [4][4] */, 
	    tcpy[16]	/* was [4][4] */;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal scale, bqra21, brqa21;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal licop[16]	/* was [4][4] */;
    static integer linfo;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal ircop[16]	/* was [4][4] */, dnorm;
    static integer iwork[4];
    extern /* Subroutine */ int slagv2_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *), sgeqr2_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *), sgerq2_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *), sorg2r_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), sorgr2_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    sorm2r_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), sormr2_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal dscale;
    extern /* Subroutine */ int stgsy2_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal thresh;
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    slassq_(integer *, doublereal *, integer *, doublereal *, 
	    doublereal *);
    static doublereal smlnum;
    static logical strong;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */
/*  Replaced various illegal calls to SCOPY by calls to SLASET, or by DO */
/*  loops. Sven Hammarling, 1/5/02. */

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

#line 280 "stgex2.f"
    /* Parameter adjustments */
#line 280 "stgex2.f"
    a_dim1 = *lda;
#line 280 "stgex2.f"
    a_offset = 1 + a_dim1;
#line 280 "stgex2.f"
    a -= a_offset;
#line 280 "stgex2.f"
    b_dim1 = *ldb;
#line 280 "stgex2.f"
    b_offset = 1 + b_dim1;
#line 280 "stgex2.f"
    b -= b_offset;
#line 280 "stgex2.f"
    q_dim1 = *ldq;
#line 280 "stgex2.f"
    q_offset = 1 + q_dim1;
#line 280 "stgex2.f"
    q -= q_offset;
#line 280 "stgex2.f"
    z_dim1 = *ldz;
#line 280 "stgex2.f"
    z_offset = 1 + z_dim1;
#line 280 "stgex2.f"
    z__ -= z_offset;
#line 280 "stgex2.f"
    --work;
#line 280 "stgex2.f"

#line 280 "stgex2.f"
    /* Function Body */
#line 280 "stgex2.f"
    *info = 0;

/*     Quick return if possible */

#line 284 "stgex2.f"
    if (*n <= 1 || *n1 <= 0 || *n2 <= 0) {
#line 284 "stgex2.f"
	return 0;
#line 284 "stgex2.f"
    }
#line 286 "stgex2.f"
    if (*n1 > *n || *j1 + *n1 > *n) {
#line 286 "stgex2.f"
	return 0;
#line 286 "stgex2.f"
    }
#line 288 "stgex2.f"
    m = *n1 + *n2;
/* Computing MAX */
#line 289 "stgex2.f"
    i__1 = *n * m, i__2 = m * m << 1;
#line 289 "stgex2.f"
    if (*lwork < max(i__1,i__2)) {
#line 290 "stgex2.f"
	*info = -16;
/* Computing MAX */
#line 291 "stgex2.f"
	i__1 = *n * m, i__2 = m * m << 1;
#line 291 "stgex2.f"
	work[1] = (doublereal) max(i__1,i__2);
#line 292 "stgex2.f"
	return 0;
#line 293 "stgex2.f"
    }

#line 295 "stgex2.f"
    weak = FALSE_;
#line 296 "stgex2.f"
    strong = FALSE_;

/*     Make a local copy of selected block */

#line 300 "stgex2.f"
    slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, li, &c__4, (ftnlen)4);
#line 301 "stgex2.f"
    slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, ir, &c__4, (ftnlen)4);
#line 302 "stgex2.f"
    slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, s, &c__4, (ftnlen)4);
#line 303 "stgex2.f"
    slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, t, &c__4, (ftnlen)4);

/*     Compute threshold for testing acceptance of swapping. */

#line 307 "stgex2.f"
    eps = slamch_("P", (ftnlen)1);
#line 308 "stgex2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 309 "stgex2.f"
    dscale = 0.;
#line 310 "stgex2.f"
    dsum = 1.;
#line 311 "stgex2.f"
    slacpy_("Full", &m, &m, s, &c__4, &work[1], &m, (ftnlen)4);
#line 312 "stgex2.f"
    i__1 = m * m;
#line 312 "stgex2.f"
    slassq_(&i__1, &work[1], &c__1, &dscale, &dsum);
#line 313 "stgex2.f"
    slacpy_("Full", &m, &m, t, &c__4, &work[1], &m, (ftnlen)4);
#line 314 "stgex2.f"
    i__1 = m * m;
#line 314 "stgex2.f"
    slassq_(&i__1, &work[1], &c__1, &dscale, &dsum);
#line 315 "stgex2.f"
    dnorm = dscale * sqrt(dsum);

/*     THRES has been changed from */
/*        THRESH = MAX( TEN*EPS*SA, SMLNUM ) */
/*     to */
/*        THRESH = MAX( TWENTY*EPS*SA, SMLNUM ) */
/*     on 04/01/10. */
/*     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by */
/*     Jim Demmel and Guillaume Revy. See forum post 1783. */

/* Computing MAX */
#line 325 "stgex2.f"
    d__1 = eps * 20. * dnorm;
#line 325 "stgex2.f"
    thresh = max(d__1,smlnum);

#line 327 "stgex2.f"
    if (m == 2) {

/*        CASE 1: Swap 1-by-1 and 1-by-1 blocks. */

/*        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks */
/*        using Givens rotations and perform the swap tentatively. */

#line 334 "stgex2.f"
	f = s[5] * t[0] - t[5] * s[0];
#line 335 "stgex2.f"
	g = s[5] * t[4] - t[5] * s[4];
#line 336 "stgex2.f"
	sb = abs(t[5]);
#line 337 "stgex2.f"
	sa = abs(s[5]);
#line 338 "stgex2.f"
	slartg_(&f, &g, &ir[4], ir, &ddum);
#line 339 "stgex2.f"
	ir[1] = -ir[4];
#line 340 "stgex2.f"
	ir[5] = ir[0];
#line 341 "stgex2.f"
	srot_(&c__2, s, &c__1, &s[4], &c__1, ir, &ir[1]);
#line 343 "stgex2.f"
	srot_(&c__2, t, &c__1, &t[4], &c__1, ir, &ir[1]);
#line 345 "stgex2.f"
	if (sa >= sb) {
#line 346 "stgex2.f"
	    slartg_(s, &s[1], li, &li[1], &ddum);
#line 348 "stgex2.f"
	} else {
#line 349 "stgex2.f"
	    slartg_(t, &t[1], li, &li[1], &ddum);
#line 351 "stgex2.f"
	}
#line 352 "stgex2.f"
	srot_(&c__2, s, &c__4, &s[1], &c__4, li, &li[1]);
#line 354 "stgex2.f"
	srot_(&c__2, t, &c__4, &t[1], &c__4, li, &li[1]);
#line 356 "stgex2.f"
	li[5] = li[0];
#line 357 "stgex2.f"
	li[4] = -li[1];

/*        Weak stability test: */
/*           |S21| + |T21| <= O(EPS * F-norm((S, T))) */

#line 362 "stgex2.f"
	ws = abs(s[1]) + abs(t[1]);
#line 363 "stgex2.f"
	weak = ws <= thresh;
#line 364 "stgex2.f"
	if (! weak) {
#line 364 "stgex2.f"
	    goto L70;
#line 364 "stgex2.f"
	}

#line 367 "stgex2.f"
	if (TRUE_) {

/*           Strong stability test: */
/*           F-norm((A-QL**T*S*QR, B-QL**T*T*QR)) <= O(EPS*F-norm((A, B))) */

#line 372 "stgex2.f"
	    slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, &work[m * m 
		    + 1], &m, (ftnlen)4);
#line 374 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, &
		    work[1], &m, (ftnlen)1, (ftnlen)1);
#line 376 "stgex2.f"
	    sgemm_("N", "T", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, &
		    c_b42, &work[m * m + 1], &m, (ftnlen)1, (ftnlen)1);
#line 378 "stgex2.f"
	    dscale = 0.;
#line 379 "stgex2.f"
	    dsum = 1.;
#line 380 "stgex2.f"
	    i__1 = m * m;
#line 380 "stgex2.f"
	    slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);

#line 382 "stgex2.f"
	    slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, &work[m * m 
		    + 1], &m, (ftnlen)4);
#line 384 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, &
		    work[1], &m, (ftnlen)1, (ftnlen)1);
#line 386 "stgex2.f"
	    sgemm_("N", "T", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, &
		    c_b42, &work[m * m + 1], &m, (ftnlen)1, (ftnlen)1);
#line 388 "stgex2.f"
	    i__1 = m * m;
#line 388 "stgex2.f"
	    slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
#line 389 "stgex2.f"
	    ss = dscale * sqrt(dsum);
#line 390 "stgex2.f"
	    strong = ss <= thresh;
#line 391 "stgex2.f"
	    if (! strong) {
#line 391 "stgex2.f"
		goto L70;
#line 391 "stgex2.f"
	    }
#line 393 "stgex2.f"
	}

/*        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and */
/*               (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)). */

#line 398 "stgex2.f"
	i__1 = *j1 + 1;
#line 398 "stgex2.f"
	srot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[(*j1 + 1) * a_dim1 + 1], 
		&c__1, ir, &ir[1]);
#line 400 "stgex2.f"
	i__1 = *j1 + 1;
#line 400 "stgex2.f"
	srot_(&i__1, &b[*j1 * b_dim1 + 1], &c__1, &b[(*j1 + 1) * b_dim1 + 1], 
		&c__1, ir, &ir[1]);
#line 402 "stgex2.f"
	i__1 = *n - *j1 + 1;
#line 402 "stgex2.f"
	srot_(&i__1, &a[*j1 + *j1 * a_dim1], lda, &a[*j1 + 1 + *j1 * a_dim1], 
		lda, li, &li[1]);
#line 404 "stgex2.f"
	i__1 = *n - *j1 + 1;
#line 404 "stgex2.f"
	srot_(&i__1, &b[*j1 + *j1 * b_dim1], ldb, &b[*j1 + 1 + *j1 * b_dim1], 
		ldb, li, &li[1]);

/*        Set  N1-by-N2 (2,1) - blocks to ZERO. */

#line 409 "stgex2.f"
	a[*j1 + 1 + *j1 * a_dim1] = 0.;
#line 410 "stgex2.f"
	b[*j1 + 1 + *j1 * b_dim1] = 0.;

/*        Accumulate transformations into Q and Z if requested. */

#line 414 "stgex2.f"
	if (*wantz) {
#line 414 "stgex2.f"
	    srot_(n, &z__[*j1 * z_dim1 + 1], &c__1, &z__[(*j1 + 1) * z_dim1 + 
		    1], &c__1, ir, &ir[1]);
#line 414 "stgex2.f"
	}
#line 417 "stgex2.f"
	if (*wantq) {
#line 417 "stgex2.f"
	    srot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[(*j1 + 1) * q_dim1 + 1], 
		    &c__1, li, &li[1]);
#line 417 "stgex2.f"
	}

/*        Exit with INFO = 0 if swap was successfully performed. */

#line 423 "stgex2.f"
	return 0;

#line 425 "stgex2.f"
    } else {

/*        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2 */
/*                and 2-by-2 blocks. */

/*        Solve the generalized Sylvester equation */
/*                 S11 * R - L * S22 = SCALE * S12 */
/*                 T11 * R - L * T22 = SCALE * T12 */
/*        for R and L. Solutions in LI and IR. */

#line 435 "stgex2.f"
	slacpy_("Full", n1, n2, &t[(*n1 + 1 << 2) - 4], &c__4, li, &c__4, (
		ftnlen)4);
#line 436 "stgex2.f"
	slacpy_("Full", n1, n2, &s[(*n1 + 1 << 2) - 4], &c__4, &ir[*n2 + 1 + (
		*n1 + 1 << 2) - 5], &c__4, (ftnlen)4);
#line 438 "stgex2.f"
	stgsy2_("N", &c__0, n1, n2, s, &c__4, &s[*n1 + 1 + (*n1 + 1 << 2) - 5]
		, &c__4, &ir[*n2 + 1 + (*n1 + 1 << 2) - 5], &c__4, t, &c__4, &
		t[*n1 + 1 + (*n1 + 1 << 2) - 5], &c__4, li, &c__4, &scale, &
		dsum, &dscale, iwork, &idum, &linfo, (ftnlen)1);

/*        Compute orthogonal matrix QL: */

/*                    QL**T * LI = [ TL ] */
/*                                 [ 0  ] */
/*        where */
/*                    LI =  [      -L              ] */
/*                          [ SCALE * identity(N2) ] */

#line 451 "stgex2.f"
	i__1 = *n2;
#line 451 "stgex2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 452 "stgex2.f"
	    sscal_(n1, &c_b48, &li[(i__ << 2) - 4], &c__1);
#line 453 "stgex2.f"
	    li[*n1 + i__ + (i__ << 2) - 5] = scale;
#line 454 "stgex2.f"
/* L10: */
#line 454 "stgex2.f"
	}
#line 455 "stgex2.f"
	sgeqr2_(&m, n2, li, &c__4, taul, &work[1], &linfo);
#line 456 "stgex2.f"
	if (linfo != 0) {
#line 456 "stgex2.f"
	    goto L70;
#line 456 "stgex2.f"
	}
#line 458 "stgex2.f"
	sorg2r_(&m, &m, n2, li, &c__4, taul, &work[1], &linfo);
#line 459 "stgex2.f"
	if (linfo != 0) {
#line 459 "stgex2.f"
	    goto L70;
#line 459 "stgex2.f"
	}

/*        Compute orthogonal matrix RQ: */

/*                    IR * RQ**T =   [ 0  TR], */

/*         where IR = [ SCALE * identity(N1), R ] */

#line 468 "stgex2.f"
	i__1 = *n1;
#line 468 "stgex2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 469 "stgex2.f"
	    ir[*n2 + i__ + (i__ << 2) - 5] = scale;
#line 470 "stgex2.f"
/* L20: */
#line 470 "stgex2.f"
	}
#line 471 "stgex2.f"
	sgerq2_(n1, &m, &ir[*n2], &c__4, taur, &work[1], &linfo);
#line 472 "stgex2.f"
	if (linfo != 0) {
#line 472 "stgex2.f"
	    goto L70;
#line 472 "stgex2.f"
	}
#line 474 "stgex2.f"
	sorgr2_(&m, &m, n1, ir, &c__4, taur, &work[1], &linfo);
#line 475 "stgex2.f"
	if (linfo != 0) {
#line 475 "stgex2.f"
	    goto L70;
#line 475 "stgex2.f"
	}

/*        Perform the swapping tentatively: */

#line 480 "stgex2.f"
	sgemm_("T", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, &
		work[1], &m, (ftnlen)1, (ftnlen)1);
#line 482 "stgex2.f"
	sgemm_("N", "T", &m, &m, &m, &c_b42, &work[1], &m, ir, &c__4, &c_b5, 
		s, &c__4, (ftnlen)1, (ftnlen)1);
#line 484 "stgex2.f"
	sgemm_("T", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, &
		work[1], &m, (ftnlen)1, (ftnlen)1);
#line 486 "stgex2.f"
	sgemm_("N", "T", &m, &m, &m, &c_b42, &work[1], &m, ir, &c__4, &c_b5, 
		t, &c__4, (ftnlen)1, (ftnlen)1);
#line 488 "stgex2.f"
	slacpy_("F", &m, &m, s, &c__4, scpy, &c__4, (ftnlen)1);
#line 489 "stgex2.f"
	slacpy_("F", &m, &m, t, &c__4, tcpy, &c__4, (ftnlen)1);
#line 490 "stgex2.f"
	slacpy_("F", &m, &m, ir, &c__4, ircop, &c__4, (ftnlen)1);
#line 491 "stgex2.f"
	slacpy_("F", &m, &m, li, &c__4, licop, &c__4, (ftnlen)1);

/*        Triangularize the B-part by an RQ factorization. */
/*        Apply transformation (from left) to A-part, giving S. */

#line 496 "stgex2.f"
	sgerq2_(&m, &m, t, &c__4, taur, &work[1], &linfo);
#line 497 "stgex2.f"
	if (linfo != 0) {
#line 497 "stgex2.f"
	    goto L70;
#line 497 "stgex2.f"
	}
#line 499 "stgex2.f"
	sormr2_("R", "T", &m, &m, &m, t, &c__4, taur, s, &c__4, &work[1], &
		linfo, (ftnlen)1, (ftnlen)1);
#line 501 "stgex2.f"
	if (linfo != 0) {
#line 501 "stgex2.f"
	    goto L70;
#line 501 "stgex2.f"
	}
#line 503 "stgex2.f"
	sormr2_("L", "N", &m, &m, &m, t, &c__4, taur, ir, &c__4, &work[1], &
		linfo, (ftnlen)1, (ftnlen)1);
#line 505 "stgex2.f"
	if (linfo != 0) {
#line 505 "stgex2.f"
	    goto L70;
#line 505 "stgex2.f"
	}

/*        Compute F-norm(S21) in BRQA21. (T21 is 0.) */

#line 510 "stgex2.f"
	dscale = 0.;
#line 511 "stgex2.f"
	dsum = 1.;
#line 512 "stgex2.f"
	i__1 = *n2;
#line 512 "stgex2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 513 "stgex2.f"
	    slassq_(n1, &s[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &dsum);
#line 514 "stgex2.f"
/* L30: */
#line 514 "stgex2.f"
	}
#line 515 "stgex2.f"
	brqa21 = dscale * sqrt(dsum);

/*        Triangularize the B-part by a QR factorization. */
/*        Apply transformation (from right) to A-part, giving S. */

#line 520 "stgex2.f"
	sgeqr2_(&m, &m, tcpy, &c__4, taul, &work[1], &linfo);
#line 521 "stgex2.f"
	if (linfo != 0) {
#line 521 "stgex2.f"
	    goto L70;
#line 521 "stgex2.f"
	}
#line 523 "stgex2.f"
	sorm2r_("L", "T", &m, &m, &m, tcpy, &c__4, taul, scpy, &c__4, &work[1]
		, info, (ftnlen)1, (ftnlen)1);
#line 525 "stgex2.f"
	sorm2r_("R", "N", &m, &m, &m, tcpy, &c__4, taul, licop, &c__4, &work[
		1], info, (ftnlen)1, (ftnlen)1);
#line 527 "stgex2.f"
	if (linfo != 0) {
#line 527 "stgex2.f"
	    goto L70;
#line 527 "stgex2.f"
	}

/*        Compute F-norm(S21) in BQRA21. (T21 is 0.) */

#line 532 "stgex2.f"
	dscale = 0.;
#line 533 "stgex2.f"
	dsum = 1.;
#line 534 "stgex2.f"
	i__1 = *n2;
#line 534 "stgex2.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 535 "stgex2.f"
	    slassq_(n1, &scpy[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &
		    dsum);
#line 536 "stgex2.f"
/* L40: */
#line 536 "stgex2.f"
	}
#line 537 "stgex2.f"
	bqra21 = dscale * sqrt(dsum);

/*        Decide which method to use. */
/*          Weak stability test: */
/*             F-norm(S21) <= O(EPS * F-norm((S, T))) */

#line 543 "stgex2.f"
	if (bqra21 <= brqa21 && bqra21 <= thresh) {
#line 544 "stgex2.f"
	    slacpy_("F", &m, &m, scpy, &c__4, s, &c__4, (ftnlen)1);
#line 545 "stgex2.f"
	    slacpy_("F", &m, &m, tcpy, &c__4, t, &c__4, (ftnlen)1);
#line 546 "stgex2.f"
	    slacpy_("F", &m, &m, ircop, &c__4, ir, &c__4, (ftnlen)1);
#line 547 "stgex2.f"
	    slacpy_("F", &m, &m, licop, &c__4, li, &c__4, (ftnlen)1);
#line 548 "stgex2.f"
	} else if (brqa21 >= thresh) {
#line 549 "stgex2.f"
	    goto L70;
#line 550 "stgex2.f"
	}

/*        Set lower triangle of B-part to zero */

#line 554 "stgex2.f"
	i__1 = m - 1;
#line 554 "stgex2.f"
	i__2 = m - 1;
#line 554 "stgex2.f"
	slaset_("Lower", &i__1, &i__2, &c_b5, &c_b5, &t[1], &c__4, (ftnlen)5);

#line 556 "stgex2.f"
	if (TRUE_) {

/*           Strong stability test: */
/*              F-norm((A-QL*S*QR**T, B-QL*T*QR**T)) <= O(EPS*F-norm((A,B))) */

#line 561 "stgex2.f"
	    slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, &work[m * m 
		    + 1], &m, (ftnlen)4);
#line 563 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, &
		    work[1], &m, (ftnlen)1, (ftnlen)1);
#line 565 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, &
		    c_b42, &work[m * m + 1], &m, (ftnlen)1, (ftnlen)1);
#line 567 "stgex2.f"
	    dscale = 0.;
#line 568 "stgex2.f"
	    dsum = 1.;
#line 569 "stgex2.f"
	    i__1 = m * m;
#line 569 "stgex2.f"
	    slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);

#line 571 "stgex2.f"
	    slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, &work[m * m 
		    + 1], &m, (ftnlen)4);
#line 573 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, &
		    work[1], &m, (ftnlen)1, (ftnlen)1);
#line 575 "stgex2.f"
	    sgemm_("N", "N", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, &
		    c_b42, &work[m * m + 1], &m, (ftnlen)1, (ftnlen)1);
#line 577 "stgex2.f"
	    i__1 = m * m;
#line 577 "stgex2.f"
	    slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
#line 578 "stgex2.f"
	    ss = dscale * sqrt(dsum);
#line 579 "stgex2.f"
	    strong = ss <= thresh;
#line 580 "stgex2.f"
	    if (! strong) {
#line 580 "stgex2.f"
		goto L70;
#line 580 "stgex2.f"
	    }

#line 583 "stgex2.f"
	}

/*        If the swap is accepted ("weakly" and "strongly"), apply the */
/*        transformations and set N1-by-N2 (2,1)-block to zero. */

#line 588 "stgex2.f"
	slaset_("Full", n1, n2, &c_b5, &c_b5, &s[*n2], &c__4, (ftnlen)4);

/*        copy back M-by-M diagonal block starting at index J1 of (A, B) */

#line 592 "stgex2.f"
	slacpy_("F", &m, &m, s, &c__4, &a[*j1 + *j1 * a_dim1], lda, (ftnlen)1)
		;
#line 593 "stgex2.f"
	slacpy_("F", &m, &m, t, &c__4, &b[*j1 + *j1 * b_dim1], ldb, (ftnlen)1)
		;
#line 594 "stgex2.f"
	slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, t, &c__4, (ftnlen)4);

/*        Standardize existing 2-by-2 blocks. */

#line 598 "stgex2.f"
	slaset_("Full", &m, &m, &c_b5, &c_b5, &work[1], &m, (ftnlen)4);
#line 599 "stgex2.f"
	work[1] = 1.;
#line 600 "stgex2.f"
	t[0] = 1.;
#line 601 "stgex2.f"
	idum = *lwork - m * m - 2;
#line 602 "stgex2.f"
	if (*n2 > 1) {
#line 603 "stgex2.f"
	    slagv2_(&a[*j1 + *j1 * a_dim1], lda, &b[*j1 + *j1 * b_dim1], ldb, 
		    ar, ai, be, &work[1], &work[2], t, &t[1]);
#line 605 "stgex2.f"
	    work[m + 1] = -work[2];
#line 606 "stgex2.f"
	    work[m + 2] = work[1];
#line 607 "stgex2.f"
	    t[*n2 + (*n2 << 2) - 5] = t[0];
#line 608 "stgex2.f"
	    t[4] = -t[1];
#line 609 "stgex2.f"
	}
#line 610 "stgex2.f"
	work[m * m] = 1.;
#line 611 "stgex2.f"
	t[m + (m << 2) - 5] = 1.;

#line 613 "stgex2.f"
	if (*n1 > 1) {
#line 614 "stgex2.f"
	    slagv2_(&a[*j1 + *n2 + (*j1 + *n2) * a_dim1], lda, &b[*j1 + *n2 + 
		    (*j1 + *n2) * b_dim1], ldb, taur, taul, &work[m * m + 1], 
		    &work[*n2 * m + *n2 + 1], &work[*n2 * m + *n2 + 2], &t[*
		    n2 + 1 + (*n2 + 1 << 2) - 5], &t[m + (m - 1 << 2) - 5]);
#line 618 "stgex2.f"
	    work[m * m] = work[*n2 * m + *n2 + 1];
#line 619 "stgex2.f"
	    work[m * m - 1] = -work[*n2 * m + *n2 + 2];
#line 620 "stgex2.f"
	    t[m + (m << 2) - 5] = t[*n2 + 1 + (*n2 + 1 << 2) - 5];
#line 621 "stgex2.f"
	    t[m - 1 + (m << 2) - 5] = -t[m + (m - 1 << 2) - 5];
#line 622 "stgex2.f"
	}
#line 623 "stgex2.f"
	sgemm_("T", "N", n2, n1, n2, &c_b42, &work[1], &m, &a[*j1 + (*j1 + *
		n2) * a_dim1], lda, &c_b5, &work[m * m + 1], n2, (ftnlen)1, (
		ftnlen)1);
#line 625 "stgex2.f"
	slacpy_("Full", n2, n1, &work[m * m + 1], n2, &a[*j1 + (*j1 + *n2) * 
		a_dim1], lda, (ftnlen)4);
#line 627 "stgex2.f"
	sgemm_("T", "N", n2, n1, n2, &c_b42, &work[1], &m, &b[*j1 + (*j1 + *
		n2) * b_dim1], ldb, &c_b5, &work[m * m + 1], n2, (ftnlen)1, (
		ftnlen)1);
#line 629 "stgex2.f"
	slacpy_("Full", n2, n1, &work[m * m + 1], n2, &b[*j1 + (*j1 + *n2) * 
		b_dim1], ldb, (ftnlen)4);
#line 631 "stgex2.f"
	sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, &work[1], &m, &c_b5, &
		work[m * m + 1], &m, (ftnlen)1, (ftnlen)1);
#line 633 "stgex2.f"
	slacpy_("Full", &m, &m, &work[m * m + 1], &m, li, &c__4, (ftnlen)4);
#line 634 "stgex2.f"
	sgemm_("N", "N", n2, n1, n1, &c_b42, &a[*j1 + (*j1 + *n2) * a_dim1], 
		lda, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &c__4, &c_b5, &work[1],
		 n2, (ftnlen)1, (ftnlen)1);
#line 636 "stgex2.f"
	slacpy_("Full", n2, n1, &work[1], n2, &a[*j1 + (*j1 + *n2) * a_dim1], 
		lda, (ftnlen)4);
#line 637 "stgex2.f"
	sgemm_("N", "N", n2, n1, n1, &c_b42, &b[*j1 + (*j1 + *n2) * b_dim1], 
		ldb, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &c__4, &c_b5, &work[1],
		 n2, (ftnlen)1, (ftnlen)1);
#line 639 "stgex2.f"
	slacpy_("Full", n2, n1, &work[1], n2, &b[*j1 + (*j1 + *n2) * b_dim1], 
		ldb, (ftnlen)4);
#line 640 "stgex2.f"
	sgemm_("T", "N", &m, &m, &m, &c_b42, ir, &c__4, t, &c__4, &c_b5, &
		work[1], &m, (ftnlen)1, (ftnlen)1);
#line 642 "stgex2.f"
	slacpy_("Full", &m, &m, &work[1], &m, ir, &c__4, (ftnlen)4);

/*        Accumulate transformations into Q and Z if requested. */

#line 646 "stgex2.f"
	if (*wantq) {
#line 647 "stgex2.f"
	    sgemm_("N", "N", n, &m, &m, &c_b42, &q[*j1 * q_dim1 + 1], ldq, li,
		     &c__4, &c_b5, &work[1], n, (ftnlen)1, (ftnlen)1);
#line 649 "stgex2.f"
	    slacpy_("Full", n, &m, &work[1], n, &q[*j1 * q_dim1 + 1], ldq, (
		    ftnlen)4);

#line 651 "stgex2.f"
	}

#line 653 "stgex2.f"
	if (*wantz) {
#line 654 "stgex2.f"
	    sgemm_("N", "N", n, &m, &m, &c_b42, &z__[*j1 * z_dim1 + 1], ldz, 
		    ir, &c__4, &c_b5, &work[1], n, (ftnlen)1, (ftnlen)1);
#line 656 "stgex2.f"
	    slacpy_("Full", n, &m, &work[1], n, &z__[*j1 * z_dim1 + 1], ldz, (
		    ftnlen)4);

#line 658 "stgex2.f"
	}

/*        Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and */
/*                (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)). */

#line 663 "stgex2.f"
	i__ = *j1 + m;
#line 664 "stgex2.f"
	if (i__ <= *n) {
#line 665 "stgex2.f"
	    i__1 = *n - i__ + 1;
#line 665 "stgex2.f"
	    sgemm_("T", "N", &m, &i__1, &m, &c_b42, li, &c__4, &a[*j1 + i__ * 
		    a_dim1], lda, &c_b5, &work[1], &m, (ftnlen)1, (ftnlen)1);
#line 667 "stgex2.f"
	    i__1 = *n - i__ + 1;
#line 667 "stgex2.f"
	    slacpy_("Full", &m, &i__1, &work[1], &m, &a[*j1 + i__ * a_dim1], 
		    lda, (ftnlen)4);
#line 668 "stgex2.f"
	    i__1 = *n - i__ + 1;
#line 668 "stgex2.f"
	    sgemm_("T", "N", &m, &i__1, &m, &c_b42, li, &c__4, &b[*j1 + i__ * 
		    b_dim1], ldb, &c_b5, &work[1], &m, (ftnlen)1, (ftnlen)1);
#line 670 "stgex2.f"
	    i__1 = *n - i__ + 1;
#line 670 "stgex2.f"
	    slacpy_("Full", &m, &i__1, &work[1], &m, &b[*j1 + i__ * b_dim1], 
		    ldb, (ftnlen)4);
#line 671 "stgex2.f"
	}
#line 672 "stgex2.f"
	i__ = *j1 - 1;
#line 673 "stgex2.f"
	if (i__ > 0) {
#line 674 "stgex2.f"
	    sgemm_("N", "N", &i__, &m, &m, &c_b42, &a[*j1 * a_dim1 + 1], lda, 
		    ir, &c__4, &c_b5, &work[1], &i__, (ftnlen)1, (ftnlen)1);
#line 676 "stgex2.f"
	    slacpy_("Full", &i__, &m, &work[1], &i__, &a[*j1 * a_dim1 + 1], 
		    lda, (ftnlen)4);
#line 677 "stgex2.f"
	    sgemm_("N", "N", &i__, &m, &m, &c_b42, &b[*j1 * b_dim1 + 1], ldb, 
		    ir, &c__4, &c_b5, &work[1], &i__, (ftnlen)1, (ftnlen)1);
#line 679 "stgex2.f"
	    slacpy_("Full", &i__, &m, &work[1], &i__, &b[*j1 * b_dim1 + 1], 
		    ldb, (ftnlen)4);
#line 680 "stgex2.f"
	}

/*        Exit with INFO = 0 if swap was successfully performed. */

#line 684 "stgex2.f"
	return 0;

#line 686 "stgex2.f"
    }

/*     Exit with INFO = 1 if swap was rejected. */

#line 690 "stgex2.f"
L70:

#line 692 "stgex2.f"
    *info = 1;
#line 693 "stgex2.f"
    return 0;

/*     End of STGEX2 */

} /* stgex2_ */

