#line 1 "ctgex2.f"
/* ctgex2.f -- translated by f2c (version 20100827).
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

#line 1 "ctgex2.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* > \brief \b CTGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an unitary 
equivalence transformation. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGEX2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgex2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgex2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgex2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/*                          LDZ, J1, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22) */
/* > in an upper triangular matrix pair (A, B) by an unitary equivalence */
/* > transformation. */
/* > */
/* > (A, B) must be in generalized Schur canonical form, that is, A and */
/* > B are both upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* >        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H */
/* >        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          Q is COMPLEX array, dimension (LDQ,N) */
/* >          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit, */
/* >          the updated matrix Q. */
/* >          Not referenced if WANTQ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= 1; */
/* >          If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ,N) */
/* >          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit, */
/* >          the updated matrix Z. */
/* >          Not referenced if WANTZ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1; */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] J1 */
/* > \verbatim */
/* >          J1 is INTEGER */
/* >          The index to the first block (A11, B11). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           =0:  Successful exit. */
/* >           =1:  The transformed matrix pair (A, B) would be too far */
/* >                from generalized Schur form; the problem is ill- */
/* >                conditioned. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complexGEauxiliary */

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
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \n */
/* >  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified */
/* >      Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* >      Estimation: Theory, Algorithms and Software, Report UMINF-94.04, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, 1994. Also as LAPACK Working Note 87. To appear in */
/* >      Numerical Algorithms, 1996. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctgex2_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *j1, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex f, g;
    static integer i__, m;
    static doublecomplex s[4]	/* was [2][2] */, t[4]	/* was [2][2] */;
    static doublereal cq, sa, sb, cz;
    static doublecomplex sq;
    static doublereal ss, ws;
    static doublecomplex sz;
    static doublereal eps, sum;
    static logical weak;
    static doublecomplex cdum;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublecomplex work[8];
    static doublereal scale;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), classq_(integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *);
    static doublereal thresh, smlnum;
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

#line 242 "ctgex2.f"
    /* Parameter adjustments */
#line 242 "ctgex2.f"
    a_dim1 = *lda;
#line 242 "ctgex2.f"
    a_offset = 1 + a_dim1;
#line 242 "ctgex2.f"
    a -= a_offset;
#line 242 "ctgex2.f"
    b_dim1 = *ldb;
#line 242 "ctgex2.f"
    b_offset = 1 + b_dim1;
#line 242 "ctgex2.f"
    b -= b_offset;
#line 242 "ctgex2.f"
    q_dim1 = *ldq;
#line 242 "ctgex2.f"
    q_offset = 1 + q_dim1;
#line 242 "ctgex2.f"
    q -= q_offset;
#line 242 "ctgex2.f"
    z_dim1 = *ldz;
#line 242 "ctgex2.f"
    z_offset = 1 + z_dim1;
#line 242 "ctgex2.f"
    z__ -= z_offset;
#line 242 "ctgex2.f"

#line 242 "ctgex2.f"
    /* Function Body */
#line 242 "ctgex2.f"
    *info = 0;

/*     Quick return if possible */

#line 246 "ctgex2.f"
    if (*n <= 1) {
#line 246 "ctgex2.f"
	return 0;
#line 246 "ctgex2.f"
    }

#line 249 "ctgex2.f"
    m = 2;
#line 250 "ctgex2.f"
    weak = FALSE_;
#line 251 "ctgex2.f"
    strong = FALSE_;

/*     Make a local copy of selected block in (A, B) */

#line 255 "ctgex2.f"
    clacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, s, &c__2, (ftnlen)4);
#line 256 "ctgex2.f"
    clacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, t, &c__2, (ftnlen)4);

/*     Compute the threshold for testing the acceptance of swapping. */

#line 260 "ctgex2.f"
    eps = slamch_("P", (ftnlen)1);
#line 261 "ctgex2.f"
    smlnum = slamch_("S", (ftnlen)1) / eps;
#line 262 "ctgex2.f"
    scale = 0.;
#line 263 "ctgex2.f"
    sum = 1.;
#line 264 "ctgex2.f"
    clacpy_("Full", &m, &m, s, &c__2, work, &m, (ftnlen)4);
#line 265 "ctgex2.f"
    clacpy_("Full", &m, &m, t, &c__2, &work[m * m], &m, (ftnlen)4);
#line 266 "ctgex2.f"
    i__1 = (m << 1) * m;
#line 266 "ctgex2.f"
    classq_(&i__1, work, &c__1, &scale, &sum);
#line 267 "ctgex2.f"
    sa = scale * sqrt(sum);

/*     THRES has been changed from */
/*        THRESH = MAX( TEN*EPS*SA, SMLNUM ) */
/*     to */
/*        THRESH = MAX( TWENTY*EPS*SA, SMLNUM ) */
/*     on 04/01/10. */
/*     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by */
/*     Jim Demmel and Guillaume Revy. See forum post 1783. */

/* Computing MAX */
#line 277 "ctgex2.f"
    d__1 = eps * 20. * sa;
#line 277 "ctgex2.f"
    thresh = max(d__1,smlnum);

/*     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks */
/*     using Givens rotations and perform the swap tentatively. */

#line 282 "ctgex2.f"
    z__2.r = s[3].r * t[0].r - s[3].i * t[0].i, z__2.i = s[3].r * t[0].i + s[
	    3].i * t[0].r;
#line 282 "ctgex2.f"
    z__3.r = t[3].r * s[0].r - t[3].i * s[0].i, z__3.i = t[3].r * s[0].i + t[
	    3].i * s[0].r;
#line 282 "ctgex2.f"
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 282 "ctgex2.f"
    f.r = z__1.r, f.i = z__1.i;
#line 283 "ctgex2.f"
    z__2.r = s[3].r * t[2].r - s[3].i * t[2].i, z__2.i = s[3].r * t[2].i + s[
	    3].i * t[2].r;
#line 283 "ctgex2.f"
    z__3.r = t[3].r * s[2].r - t[3].i * s[2].i, z__3.i = t[3].r * s[2].i + t[
	    3].i * s[2].r;
#line 283 "ctgex2.f"
    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 283 "ctgex2.f"
    g.r = z__1.r, g.i = z__1.i;
#line 284 "ctgex2.f"
    sa = z_abs(&s[3]);
#line 285 "ctgex2.f"
    sb = z_abs(&t[3]);
#line 286 "ctgex2.f"
    clartg_(&g, &f, &cz, &sz, &cdum);
#line 287 "ctgex2.f"
    z__1.r = -sz.r, z__1.i = -sz.i;
#line 287 "ctgex2.f"
    sz.r = z__1.r, sz.i = z__1.i;
#line 288 "ctgex2.f"
    d_cnjg(&z__1, &sz);
#line 288 "ctgex2.f"
    crot_(&c__2, s, &c__1, &s[2], &c__1, &cz, &z__1);
#line 289 "ctgex2.f"
    d_cnjg(&z__1, &sz);
#line 289 "ctgex2.f"
    crot_(&c__2, t, &c__1, &t[2], &c__1, &cz, &z__1);
#line 290 "ctgex2.f"
    if (sa >= sb) {
#line 291 "ctgex2.f"
	clartg_(s, &s[1], &cq, &sq, &cdum);
#line 292 "ctgex2.f"
    } else {
#line 293 "ctgex2.f"
	clartg_(t, &t[1], &cq, &sq, &cdum);
#line 294 "ctgex2.f"
    }
#line 295 "ctgex2.f"
    crot_(&c__2, s, &c__2, &s[1], &c__2, &cq, &sq);
#line 296 "ctgex2.f"
    crot_(&c__2, t, &c__2, &t[1], &c__2, &cq, &sq);

/*     Weak stability test: |S21| + |T21| <= O(EPS F-norm((S, T))) */

#line 300 "ctgex2.f"
    ws = z_abs(&s[1]) + z_abs(&t[1]);
#line 301 "ctgex2.f"
    weak = ws <= thresh;
#line 302 "ctgex2.f"
    if (! weak) {
#line 302 "ctgex2.f"
	goto L20;
#line 302 "ctgex2.f"
    }

#line 305 "ctgex2.f"
    if (TRUE_) {

/*        Strong stability test: */
/*           F-norm((A-QL**H*S*QR, B-QL**H*T*QR)) <= O(EPS*F-norm((A, B))) */

#line 310 "ctgex2.f"
	clacpy_("Full", &m, &m, s, &c__2, work, &m, (ftnlen)4);
#line 311 "ctgex2.f"
	clacpy_("Full", &m, &m, t, &c__2, &work[m * m], &m, (ftnlen)4);
#line 312 "ctgex2.f"
	d_cnjg(&z__2, &sz);
#line 312 "ctgex2.f"
	z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 312 "ctgex2.f"
	crot_(&c__2, work, &c__1, &work[2], &c__1, &cz, &z__1);
#line 313 "ctgex2.f"
	d_cnjg(&z__2, &sz);
#line 313 "ctgex2.f"
	z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 313 "ctgex2.f"
	crot_(&c__2, &work[4], &c__1, &work[6], &c__1, &cz, &z__1);
#line 314 "ctgex2.f"
	z__1.r = -sq.r, z__1.i = -sq.i;
#line 314 "ctgex2.f"
	crot_(&c__2, work, &c__2, &work[1], &c__2, &cq, &z__1);
#line 315 "ctgex2.f"
	z__1.r = -sq.r, z__1.i = -sq.i;
#line 315 "ctgex2.f"
	crot_(&c__2, &work[4], &c__2, &work[5], &c__2, &cq, &z__1);
#line 316 "ctgex2.f"
	for (i__ = 1; i__ <= 2; ++i__) {
#line 317 "ctgex2.f"
	    i__1 = i__ - 1;
#line 317 "ctgex2.f"
	    i__2 = i__ - 1;
#line 317 "ctgex2.f"
	    i__3 = *j1 + i__ - 1 + *j1 * a_dim1;
#line 317 "ctgex2.f"
	    z__1.r = work[i__2].r - a[i__3].r, z__1.i = work[i__2].i - a[i__3]
		    .i;
#line 317 "ctgex2.f"
	    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 318 "ctgex2.f"
	    i__1 = i__ + 1;
#line 318 "ctgex2.f"
	    i__2 = i__ + 1;
#line 318 "ctgex2.f"
	    i__3 = *j1 + i__ - 1 + (*j1 + 1) * a_dim1;
#line 318 "ctgex2.f"
	    z__1.r = work[i__2].r - a[i__3].r, z__1.i = work[i__2].i - a[i__3]
		    .i;
#line 318 "ctgex2.f"
	    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 319 "ctgex2.f"
	    i__1 = i__ + 3;
#line 319 "ctgex2.f"
	    i__2 = i__ + 3;
#line 319 "ctgex2.f"
	    i__3 = *j1 + i__ - 1 + *j1 * b_dim1;
#line 319 "ctgex2.f"
	    z__1.r = work[i__2].r - b[i__3].r, z__1.i = work[i__2].i - b[i__3]
		    .i;
#line 319 "ctgex2.f"
	    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 320 "ctgex2.f"
	    i__1 = i__ + 5;
#line 320 "ctgex2.f"
	    i__2 = i__ + 5;
#line 320 "ctgex2.f"
	    i__3 = *j1 + i__ - 1 + (*j1 + 1) * b_dim1;
#line 320 "ctgex2.f"
	    z__1.r = work[i__2].r - b[i__3].r, z__1.i = work[i__2].i - b[i__3]
		    .i;
#line 320 "ctgex2.f"
	    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 321 "ctgex2.f"
/* L10: */
#line 321 "ctgex2.f"
	}
#line 322 "ctgex2.f"
	scale = 0.;
#line 323 "ctgex2.f"
	sum = 1.;
#line 324 "ctgex2.f"
	i__1 = (m << 1) * m;
#line 324 "ctgex2.f"
	classq_(&i__1, work, &c__1, &scale, &sum);
#line 325 "ctgex2.f"
	ss = scale * sqrt(sum);
#line 326 "ctgex2.f"
	strong = ss <= thresh;
#line 327 "ctgex2.f"
	if (! strong) {
#line 327 "ctgex2.f"
	    goto L20;
#line 327 "ctgex2.f"
	}
#line 329 "ctgex2.f"
    }

/*     If the swap is accepted ("weakly" and "strongly"), apply the */
/*     equivalence transformations to the original matrix pair (A,B) */

#line 334 "ctgex2.f"
    i__1 = *j1 + 1;
#line 334 "ctgex2.f"
    d_cnjg(&z__1, &sz);
#line 334 "ctgex2.f"
    crot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[(*j1 + 1) * a_dim1 + 1], &
	    c__1, &cz, &z__1);
#line 335 "ctgex2.f"
    i__1 = *j1 + 1;
#line 335 "ctgex2.f"
    d_cnjg(&z__1, &sz);
#line 335 "ctgex2.f"
    crot_(&i__1, &b[*j1 * b_dim1 + 1], &c__1, &b[(*j1 + 1) * b_dim1 + 1], &
	    c__1, &cz, &z__1);
#line 336 "ctgex2.f"
    i__1 = *n - *j1 + 1;
#line 336 "ctgex2.f"
    crot_(&i__1, &a[*j1 + *j1 * a_dim1], lda, &a[*j1 + 1 + *j1 * a_dim1], lda,
	     &cq, &sq);
#line 337 "ctgex2.f"
    i__1 = *n - *j1 + 1;
#line 337 "ctgex2.f"
    crot_(&i__1, &b[*j1 + *j1 * b_dim1], ldb, &b[*j1 + 1 + *j1 * b_dim1], ldb,
	     &cq, &sq);

/*     Set  N1 by N2 (2,1) blocks to 0 */

#line 341 "ctgex2.f"
    i__1 = *j1 + 1 + *j1 * a_dim1;
#line 341 "ctgex2.f"
    a[i__1].r = 0., a[i__1].i = 0.;
#line 342 "ctgex2.f"
    i__1 = *j1 + 1 + *j1 * b_dim1;
#line 342 "ctgex2.f"
    b[i__1].r = 0., b[i__1].i = 0.;

/*     Accumulate transformations into Q and Z if requested. */

#line 346 "ctgex2.f"
    if (*wantz) {
#line 346 "ctgex2.f"
	d_cnjg(&z__1, &sz);
#line 346 "ctgex2.f"
	crot_(n, &z__[*j1 * z_dim1 + 1], &c__1, &z__[(*j1 + 1) * z_dim1 + 1], 
		&c__1, &cz, &z__1);
#line 346 "ctgex2.f"
    }
#line 348 "ctgex2.f"
    if (*wantq) {
#line 348 "ctgex2.f"
	d_cnjg(&z__1, &sq);
#line 348 "ctgex2.f"
	crot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[(*j1 + 1) * q_dim1 + 1], &
		c__1, &cq, &z__1);
#line 348 "ctgex2.f"
    }

/*     Exit with INFO = 0 if swap was successfully performed. */

#line 353 "ctgex2.f"
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

#line 357 "ctgex2.f"
L20:
#line 358 "ctgex2.f"
    *info = 1;
#line 359 "ctgex2.f"
    return 0;

/*     End of CTGEX2 */

} /* ctgex2_ */

