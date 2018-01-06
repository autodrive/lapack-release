#line 1 "cgghd3.f"
/* cgghd3.f -- translated by f2c (version 20100827).
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

#line 1 "cgghd3.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__16 = 16;

/* > \brief \b CGGHD3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGHD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgghd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgghd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgghd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*        SUBROUTINE CGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*       $                   LDQ, Z, LDZ, WORK, LWORK, INFO ) */

/*        .. Scalar Arguments .. */
/*        CHARACTER          COMPQ, COMPZ */
/*        INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK */
/*        .. */
/*        .. Array Arguments .. */
/*        COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*       $                   Z( LDZ, * ), WORK( * ) */
/*        .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CGGHD3 reduces a pair of complex matrices (A,B) to generalized upper */
/* > Hessenberg form using unitary transformations, where A is a */
/* > general matrix and B is upper triangular.  The form of the */
/* > generalized eigenvalue problem is */
/* >    A*x = lambda*B*x, */
/* > and B is typically made upper triangular by computing its QR */
/* > factorization and moving the unitary matrix Q to the left side */
/* > of the equation. */
/* > */
/* > This subroutine simultaneously reduces A to a Hessenberg matrix H: */
/* >    Q**H*A*Z = H */
/* > and transforms B to another upper triangular matrix T: */
/* >    Q**H*B*Z = T */
/* > in order to reduce the problem to its standard form */
/* >    H*y = lambda*T*y */
/* > where y = Z**H*x. */
/* > */
/* > The unitary matrices Q and Z are determined as products of Givens */
/* > rotations.  They may either be formed explicitly, or they may be */
/* > postmultiplied into input matrices Q1 and Z1, so that */
/* > */
/* >      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H */
/* > */
/* >      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H */
/* > */
/* > If Q1 is the unitary matrix from the QR factorization of B in the */
/* > original equation A*x = lambda*B*x, then CGGHD3 reduces the original */
/* > problem to generalized Hessenberg form. */
/* > */
/* > This is a blocked variant of CGGHRD, using matrix-matrix */
/* > multiplications for parts of the computation to enhance performance. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'N': do not compute Q; */
/* >          = 'I': Q is initialized to the unit matrix, and the */
/* >                 unitary matrix Q is returned; */
/* >          = 'V': Q must contain a unitary matrix Q1 on entry, */
/* >                 and the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N': do not compute Z; */
/* >          = 'I': Z is initialized to the unit matrix, and the */
/* >                 unitary matrix Z is returned; */
/* >          = 'V': Z must contain a unitary matrix Z1 on entry, */
/* >                 and the product Z1*Z is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
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
/* > */
/* >          ILO and IHI mark the rows and columns of A which are to be */
/* >          reduced.  It is assumed that A is already upper triangular */
/* >          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are */
/* >          normally set by a previous call to CGGBAL; otherwise they */
/* >          should be set to 1 and N respectively. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          On entry, the N-by-N general matrix to be reduced. */
/* >          On exit, the upper triangle and the first subdiagonal of A */
/* >          are overwritten with the upper Hessenberg matrix H, and the */
/* >          rest is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB, N) */
/* >          On entry, the N-by-N upper triangular matrix B. */
/* >          On exit, the upper triangular matrix T = Q**H B Z.  The */
/* >          elements below the diagonal are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ, N) */
/* >          On entry, if COMPQ = 'V', the unitary matrix Q1, typically */
/* >          from the QR factorization of B. */
/* >          On exit, if COMPQ='I', the unitary matrix Q, and if */
/* >          COMPQ = 'V', the product Q1*Q. */
/* >          Not referenced if COMPQ='N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the unitary matrix Z1. */
/* >          On exit, if COMPZ='I', the unitary matrix Z, and if */
/* >          COMPZ = 'V', the product Z1*Z. */
/* >          Not referenced if COMPZ='N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. */
/* >          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (LWORK) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in]  LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK.  LWORK >= 1. */
/* >          For optimum performance LWORK >= 6*N*NB, where NB is the */
/* >          optimal blocksize. */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date January 2015 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  This routine reduces A to Hessenberg form and maintains B in */
/* >  using a blocked variant of Moler and Stewart's original algorithm, */
/* >  as described by Kagstrom, Kressner, Quintana-Orti, and Quintana-Orti */
/* >  (BIT 2008). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgghd3_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, integer *info, 
	ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublecomplex s, c1, c2;
    static integer j0;
    static doublecomplex s1, s2;
    static integer nb, jj, nh, nx, pw, nnb, len, top, ppw, n2nb;
    static logical blk22;
    static integer cola, jcol, ierr;
    static doublecomplex temp;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer jrow, topq, ppwo;
    static doublecomplex temp1, temp2, temp3;
    static integer kacc22;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int cunm22_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen)
	    ;
    static doublecomplex ctemp;
    static integer nblst;
    static logical initq, wantq;
    extern /* Subroutine */ int ctrmv_(char *, char *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static logical initz, wantz;
    static char compq2[1], compz2[1];
    extern /* Subroutine */ int cgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), clacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     January 2015 */



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
/*     .. Executable Statements .. */

/*     Decode and test the input parameters. */

#line 284 "cgghd3.f"
    /* Parameter adjustments */
#line 284 "cgghd3.f"
    a_dim1 = *lda;
#line 284 "cgghd3.f"
    a_offset = 1 + a_dim1;
#line 284 "cgghd3.f"
    a -= a_offset;
#line 284 "cgghd3.f"
    b_dim1 = *ldb;
#line 284 "cgghd3.f"
    b_offset = 1 + b_dim1;
#line 284 "cgghd3.f"
    b -= b_offset;
#line 284 "cgghd3.f"
    q_dim1 = *ldq;
#line 284 "cgghd3.f"
    q_offset = 1 + q_dim1;
#line 284 "cgghd3.f"
    q -= q_offset;
#line 284 "cgghd3.f"
    z_dim1 = *ldz;
#line 284 "cgghd3.f"
    z_offset = 1 + z_dim1;
#line 284 "cgghd3.f"
    z__ -= z_offset;
#line 284 "cgghd3.f"
    --work;
#line 284 "cgghd3.f"

#line 284 "cgghd3.f"
    /* Function Body */
#line 284 "cgghd3.f"
    *info = 0;
#line 285 "cgghd3.f"
    nb = ilaenv_(&c__1, "CGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)
	    1);
/* Computing MAX */
#line 286 "cgghd3.f"
    i__1 = *n * 6 * nb;
#line 286 "cgghd3.f"
    lwkopt = max(i__1,1);
#line 287 "cgghd3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 287 "cgghd3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 288 "cgghd3.f"
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 289 "cgghd3.f"
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 290 "cgghd3.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 291 "cgghd3.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);
#line 292 "cgghd3.f"
    lquery = *lwork == -1;

#line 294 "cgghd3.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 295 "cgghd3.f"
	*info = -1;
#line 296 "cgghd3.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 297 "cgghd3.f"
	*info = -2;
#line 298 "cgghd3.f"
    } else if (*n < 0) {
#line 299 "cgghd3.f"
	*info = -3;
#line 300 "cgghd3.f"
    } else if (*ilo < 1) {
#line 301 "cgghd3.f"
	*info = -4;
#line 302 "cgghd3.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 303 "cgghd3.f"
	*info = -5;
#line 304 "cgghd3.f"
    } else if (*lda < max(1,*n)) {
#line 305 "cgghd3.f"
	*info = -7;
#line 306 "cgghd3.f"
    } else if (*ldb < max(1,*n)) {
#line 307 "cgghd3.f"
	*info = -9;
#line 308 "cgghd3.f"
    } else if (wantq && *ldq < *n || *ldq < 1) {
#line 309 "cgghd3.f"
	*info = -11;
#line 310 "cgghd3.f"
    } else if (wantz && *ldz < *n || *ldz < 1) {
#line 311 "cgghd3.f"
	*info = -13;
#line 312 "cgghd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 313 "cgghd3.f"
	*info = -15;
#line 314 "cgghd3.f"
    }
#line 315 "cgghd3.f"
    if (*info != 0) {
#line 316 "cgghd3.f"
	i__1 = -(*info);
#line 316 "cgghd3.f"
	xerbla_("CGGHD3", &i__1, (ftnlen)6);
#line 317 "cgghd3.f"
	return 0;
#line 318 "cgghd3.f"
    } else if (lquery) {
#line 319 "cgghd3.f"
	return 0;
#line 320 "cgghd3.f"
    }

/*     Initialize Q and Z if desired. */

#line 324 "cgghd3.f"
    if (initq) {
#line 324 "cgghd3.f"
	claset_("All", n, n, &c_b2, &c_b1, &q[q_offset], ldq, (ftnlen)3);
#line 324 "cgghd3.f"
    }
#line 326 "cgghd3.f"
    if (initz) {
#line 326 "cgghd3.f"
	claset_("All", n, n, &c_b2, &c_b1, &z__[z_offset], ldz, (ftnlen)3);
#line 326 "cgghd3.f"
    }

/*     Zero out lower triangle of B. */

#line 331 "cgghd3.f"
    if (*n > 1) {
#line 331 "cgghd3.f"
	i__1 = *n - 1;
#line 331 "cgghd3.f"
	i__2 = *n - 1;
#line 331 "cgghd3.f"
	claset_("Lower", &i__1, &i__2, &c_b2, &c_b2, &b[b_dim1 + 2], ldb, (
		ftnlen)5);
#line 331 "cgghd3.f"
    }

/*     Quick return if possible */

#line 336 "cgghd3.f"
    nh = *ihi - *ilo + 1;
#line 337 "cgghd3.f"
    if (nh <= 1) {
#line 338 "cgghd3.f"
	work[1].r = 1., work[1].i = 0.;
#line 339 "cgghd3.f"
	return 0;
#line 340 "cgghd3.f"
    }

/*     Determine the blocksize. */

#line 344 "cgghd3.f"
    nbmin = ilaenv_(&c__2, "CGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 345 "cgghd3.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to use unblocked instead of blocked code. */

/* Computing MAX */
#line 349 "cgghd3.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "CGGHD3", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 349 "cgghd3.f"
	nx = max(i__1,i__2);
#line 350 "cgghd3.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code. */

#line 354 "cgghd3.f"
	    if (*lwork < lwkopt) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code. */

/* Computing MAX */
#line 360 "cgghd3.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "CGGHD3", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 360 "cgghd3.f"
		nbmin = max(i__1,i__2);
#line 362 "cgghd3.f"
		if (*lwork >= *n * 6 * nbmin) {
#line 363 "cgghd3.f"
		    nb = *lwork / (*n * 6);
#line 364 "cgghd3.f"
		} else {
#line 365 "cgghd3.f"
		    nb = 1;
#line 366 "cgghd3.f"
		}
#line 367 "cgghd3.f"
	    }
#line 368 "cgghd3.f"
	}
#line 369 "cgghd3.f"
    }

#line 371 "cgghd3.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 375 "cgghd3.f"
	jcol = *ilo;

#line 377 "cgghd3.f"
    } else {

/*        Use blocked code */

#line 381 "cgghd3.f"
	kacc22 = ilaenv_(&c__16, "CGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 382 "cgghd3.f"
	blk22 = kacc22 == 2;
#line 383 "cgghd3.f"
	i__1 = *ihi - 2;
#line 383 "cgghd3.f"
	i__2 = nb;
#line 383 "cgghd3.f"
	for (jcol = *ilo; i__2 < 0 ? jcol >= i__1 : jcol <= i__1; jcol += 
		i__2) {
/* Computing MIN */
#line 384 "cgghd3.f"
	    i__3 = nb, i__4 = *ihi - jcol - 1;
#line 384 "cgghd3.f"
	    nnb = min(i__3,i__4);

/*           Initialize small unitary factors that will hold the */
/*           accumulated Givens rotations in workspace. */
/*           N2NB   denotes the number of 2*NNB-by-2*NNB factors */
/*           NBLST  denotes the (possibly smaller) order of the last */
/*                  factor. */

#line 392 "cgghd3.f"
	    n2nb = (*ihi - jcol - 1) / nnb - 1;
#line 393 "cgghd3.f"
	    nblst = *ihi - jcol - n2nb * nnb;
#line 394 "cgghd3.f"
	    claset_("All", &nblst, &nblst, &c_b2, &c_b1, &work[1], &nblst, (
		    ftnlen)3);
#line 395 "cgghd3.f"
	    pw = nblst * nblst + 1;
#line 396 "cgghd3.f"
	    i__3 = n2nb;
#line 396 "cgghd3.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 397 "cgghd3.f"
		i__4 = nnb << 1;
#line 397 "cgghd3.f"
		i__5 = nnb << 1;
#line 397 "cgghd3.f"
		i__6 = nnb << 1;
#line 397 "cgghd3.f"
		claset_("All", &i__4, &i__5, &c_b2, &c_b1, &work[pw], &i__6, (
			ftnlen)3);
#line 399 "cgghd3.f"
		pw += (nnb << 2) * nnb;
#line 400 "cgghd3.f"
	    }

/*           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form. */

#line 404 "cgghd3.f"
	    i__3 = jcol + nnb - 1;
#line 404 "cgghd3.f"
	    for (j = jcol; j <= i__3; ++j) {

/*              Reduce Jth column of A. Store cosines and sines in Jth */
/*              column of A and B, respectively. */

#line 409 "cgghd3.f"
		i__4 = j + 2;
#line 409 "cgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 410 "cgghd3.f"
		    i__5 = i__ - 1 + j * a_dim1;
#line 410 "cgghd3.f"
		    temp.r = a[i__5].r, temp.i = a[i__5].i;
#line 411 "cgghd3.f"
		    clartg_(&temp, &a[i__ + j * a_dim1], &c__, &s, &a[i__ - 1 
			    + j * a_dim1]);
#line 412 "cgghd3.f"
		    i__5 = i__ + j * a_dim1;
#line 412 "cgghd3.f"
		    z__1.r = c__, z__1.i = 0.;
#line 412 "cgghd3.f"
		    a[i__5].r = z__1.r, a[i__5].i = z__1.i;
#line 413 "cgghd3.f"
		    i__5 = i__ + j * b_dim1;
#line 413 "cgghd3.f"
		    b[i__5].r = s.r, b[i__5].i = s.i;
#line 414 "cgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 418 "cgghd3.f"
		ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 419 "cgghd3.f"
		len = j + 2 - jcol;
#line 420 "cgghd3.f"
		jrow = j + n2nb * nnb + 2;
#line 421 "cgghd3.f"
		i__4 = jrow;
#line 421 "cgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 422 "cgghd3.f"
		    i__5 = i__ + j * a_dim1;
#line 422 "cgghd3.f"
		    ctemp.r = a[i__5].r, ctemp.i = a[i__5].i;
#line 423 "cgghd3.f"
		    i__5 = i__ + j * b_dim1;
#line 423 "cgghd3.f"
		    s.r = b[i__5].r, s.i = b[i__5].i;
#line 424 "cgghd3.f"
		    i__5 = ppw + len - 1;
#line 424 "cgghd3.f"
		    for (jj = ppw; jj <= i__5; ++jj) {
#line 425 "cgghd3.f"
			i__6 = jj + nblst;
#line 425 "cgghd3.f"
			temp.r = work[i__6].r, temp.i = work[i__6].i;
#line 426 "cgghd3.f"
			i__6 = jj + nblst;
#line 426 "cgghd3.f"
			z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, z__2.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 426 "cgghd3.f"
			i__7 = jj;
#line 426 "cgghd3.f"
			z__3.r = s.r * work[i__7].r - s.i * work[i__7].i, 
				z__3.i = s.r * work[i__7].i + s.i * work[i__7]
				.r;
#line 426 "cgghd3.f"
			z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 426 "cgghd3.f"
			work[i__6].r = z__1.r, work[i__6].i = z__1.i;
#line 427 "cgghd3.f"
			i__6 = jj;
#line 427 "cgghd3.f"
			d_cnjg(&z__3, &s);
#line 427 "cgghd3.f"
			z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
				z__3.r * temp.i + z__3.i * temp.r;
#line 427 "cgghd3.f"
			i__7 = jj;
#line 427 "cgghd3.f"
			z__4.r = ctemp.r * work[i__7].r - ctemp.i * work[i__7]
				.i, z__4.i = ctemp.r * work[i__7].i + ctemp.i 
				* work[i__7].r;
#line 427 "cgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 427 "cgghd3.f"
			work[i__6].r = z__1.r, work[i__6].i = z__1.i;
#line 428 "cgghd3.f"
		    }
#line 429 "cgghd3.f"
		    ++len;
#line 430 "cgghd3.f"
		    ppw = ppw - nblst - 1;
#line 431 "cgghd3.f"
		}

#line 433 "cgghd3.f"
		ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + nnb;
#line 434 "cgghd3.f"
		j0 = jrow - nnb;
#line 435 "cgghd3.f"
		i__4 = j + 2;
#line 435 "cgghd3.f"
		i__5 = -nnb;
#line 435 "cgghd3.f"
		for (jrow = j0; i__5 < 0 ? jrow >= i__4 : jrow <= i__4; jrow 
			+= i__5) {
#line 436 "cgghd3.f"
		    ppw = ppwo;
#line 437 "cgghd3.f"
		    len = j + 2 - jcol;
#line 438 "cgghd3.f"
		    i__6 = jrow;
#line 438 "cgghd3.f"
		    for (i__ = jrow + nnb - 1; i__ >= i__6; --i__) {
#line 439 "cgghd3.f"
			i__7 = i__ + j * a_dim1;
#line 439 "cgghd3.f"
			ctemp.r = a[i__7].r, ctemp.i = a[i__7].i;
#line 440 "cgghd3.f"
			i__7 = i__ + j * b_dim1;
#line 440 "cgghd3.f"
			s.r = b[i__7].r, s.i = b[i__7].i;
#line 441 "cgghd3.f"
			i__7 = ppw + len - 1;
#line 441 "cgghd3.f"
			for (jj = ppw; jj <= i__7; ++jj) {
#line 442 "cgghd3.f"
			    i__8 = jj + (nnb << 1);
#line 442 "cgghd3.f"
			    temp.r = work[i__8].r, temp.i = work[i__8].i;
#line 443 "cgghd3.f"
			    i__8 = jj + (nnb << 1);
#line 443 "cgghd3.f"
			    z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
				    z__2.i = ctemp.r * temp.i + ctemp.i * 
				    temp.r;
#line 443 "cgghd3.f"
			    i__9 = jj;
#line 443 "cgghd3.f"
			    z__3.r = s.r * work[i__9].r - s.i * work[i__9].i, 
				    z__3.i = s.r * work[i__9].i + s.i * work[
				    i__9].r;
#line 443 "cgghd3.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 443 "cgghd3.f"
			    work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 444 "cgghd3.f"
			    i__8 = jj;
#line 444 "cgghd3.f"
			    d_cnjg(&z__3, &s);
#line 444 "cgghd3.f"
			    z__2.r = z__3.r * temp.r - z__3.i * temp.i, 
				    z__2.i = z__3.r * temp.i + z__3.i * 
				    temp.r;
#line 444 "cgghd3.f"
			    i__9 = jj;
#line 444 "cgghd3.f"
			    z__4.r = ctemp.r * work[i__9].r - ctemp.i * work[
				    i__9].i, z__4.i = ctemp.r * work[i__9].i 
				    + ctemp.i * work[i__9].r;
#line 444 "cgghd3.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 444 "cgghd3.f"
			    work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 445 "cgghd3.f"
			}
#line 446 "cgghd3.f"
			++len;
#line 447 "cgghd3.f"
			ppw = ppw - (nnb << 1) - 1;
#line 448 "cgghd3.f"
		    }
#line 449 "cgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 450 "cgghd3.f"
		}

/*              TOP denotes the number of top rows in A and B that will */
/*              not be updated during the next steps. */

#line 455 "cgghd3.f"
		if (jcol <= 2) {
#line 456 "cgghd3.f"
		    top = 0;
#line 457 "cgghd3.f"
		} else {
#line 458 "cgghd3.f"
		    top = jcol;
#line 459 "cgghd3.f"
		}

/*              Propagate transformations through B and replace stored */
/*              left sines/cosines by right sines/cosines. */

#line 464 "cgghd3.f"
		i__5 = j + 1;
#line 464 "cgghd3.f"
		for (jj = *n; jj >= i__5; --jj) {

/*                 Update JJth column of B. */

/* Computing MIN */
#line 468 "cgghd3.f"
		    i__4 = jj + 1;
#line 468 "cgghd3.f"
		    i__6 = j + 2;
#line 468 "cgghd3.f"
		    for (i__ = min(i__4,*ihi); i__ >= i__6; --i__) {
#line 469 "cgghd3.f"
			i__4 = i__ + j * a_dim1;
#line 469 "cgghd3.f"
			ctemp.r = a[i__4].r, ctemp.i = a[i__4].i;
#line 470 "cgghd3.f"
			i__4 = i__ + j * b_dim1;
#line 470 "cgghd3.f"
			s.r = b[i__4].r, s.i = b[i__4].i;
#line 471 "cgghd3.f"
			i__4 = i__ + jj * b_dim1;
#line 471 "cgghd3.f"
			temp.r = b[i__4].r, temp.i = b[i__4].i;
#line 472 "cgghd3.f"
			i__4 = i__ + jj * b_dim1;
#line 472 "cgghd3.f"
			z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, z__2.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 472 "cgghd3.f"
			d_cnjg(&z__4, &s);
#line 472 "cgghd3.f"
			i__7 = i__ - 1 + jj * b_dim1;
#line 472 "cgghd3.f"
			z__3.r = z__4.r * b[i__7].r - z__4.i * b[i__7].i, 
				z__3.i = z__4.r * b[i__7].i + z__4.i * b[i__7]
				.r;
#line 472 "cgghd3.f"
			z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 472 "cgghd3.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 473 "cgghd3.f"
			i__4 = i__ - 1 + jj * b_dim1;
#line 473 "cgghd3.f"
			z__2.r = s.r * temp.r - s.i * temp.i, z__2.i = s.r * 
				temp.i + s.i * temp.r;
#line 473 "cgghd3.f"
			i__7 = i__ - 1 + jj * b_dim1;
#line 473 "cgghd3.f"
			z__3.r = ctemp.r * b[i__7].r - ctemp.i * b[i__7].i, 
				z__3.i = ctemp.r * b[i__7].i + ctemp.i * b[
				i__7].r;
#line 473 "cgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 473 "cgghd3.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 474 "cgghd3.f"
		    }

/*                 Annihilate B( JJ+1, JJ ). */

#line 478 "cgghd3.f"
		    if (jj < *ihi) {
#line 479 "cgghd3.f"
			i__6 = jj + 1 + (jj + 1) * b_dim1;
#line 479 "cgghd3.f"
			temp.r = b[i__6].r, temp.i = b[i__6].i;
#line 480 "cgghd3.f"
			clartg_(&temp, &b[jj + 1 + jj * b_dim1], &c__, &s, &b[
				jj + 1 + (jj + 1) * b_dim1]);
#line 482 "cgghd3.f"
			i__6 = jj + 1 + jj * b_dim1;
#line 482 "cgghd3.f"
			b[i__6].r = 0., b[i__6].i = 0.;
#line 483 "cgghd3.f"
			i__6 = jj - top;
#line 483 "cgghd3.f"
			crot_(&i__6, &b[top + 1 + (jj + 1) * b_dim1], &c__1, &
				b[top + 1 + jj * b_dim1], &c__1, &c__, &s);
#line 485 "cgghd3.f"
			i__6 = jj + 1 + j * a_dim1;
#line 485 "cgghd3.f"
			z__1.r = c__, z__1.i = 0.;
#line 485 "cgghd3.f"
			a[i__6].r = z__1.r, a[i__6].i = z__1.i;
#line 486 "cgghd3.f"
			i__6 = jj + 1 + j * b_dim1;
#line 486 "cgghd3.f"
			d_cnjg(&z__2, &s);
#line 486 "cgghd3.f"
			z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 486 "cgghd3.f"
			b[i__6].r = z__1.r, b[i__6].i = z__1.i;
#line 487 "cgghd3.f"
		    }
#line 488 "cgghd3.f"
		}

/*              Update A by transformations from right. */

#line 492 "cgghd3.f"
		jj = (*ihi - j - 1) % 3;
#line 493 "cgghd3.f"
		i__5 = jj + 1;
#line 493 "cgghd3.f"
		for (i__ = *ihi - j - 3; i__ >= i__5; i__ += -3) {
#line 494 "cgghd3.f"
		    i__6 = j + 1 + i__ + j * a_dim1;
#line 494 "cgghd3.f"
		    ctemp.r = a[i__6].r, ctemp.i = a[i__6].i;
#line 495 "cgghd3.f"
		    i__6 = j + 1 + i__ + j * b_dim1;
#line 495 "cgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 495 "cgghd3.f"
		    s.r = z__1.r, s.i = z__1.i;
#line 496 "cgghd3.f"
		    i__6 = j + 2 + i__ + j * a_dim1;
#line 496 "cgghd3.f"
		    c1.r = a[i__6].r, c1.i = a[i__6].i;
#line 497 "cgghd3.f"
		    i__6 = j + 2 + i__ + j * b_dim1;
#line 497 "cgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 497 "cgghd3.f"
		    s1.r = z__1.r, s1.i = z__1.i;
#line 498 "cgghd3.f"
		    i__6 = j + 3 + i__ + j * a_dim1;
#line 498 "cgghd3.f"
		    c2.r = a[i__6].r, c2.i = a[i__6].i;
#line 499 "cgghd3.f"
		    i__6 = j + 3 + i__ + j * b_dim1;
#line 499 "cgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 499 "cgghd3.f"
		    s2.r = z__1.r, s2.i = z__1.i;

#line 501 "cgghd3.f"
		    i__6 = *ihi;
#line 501 "cgghd3.f"
		    for (k = top + 1; k <= i__6; ++k) {
#line 502 "cgghd3.f"
			i__4 = k + (j + i__) * a_dim1;
#line 502 "cgghd3.f"
			temp.r = a[i__4].r, temp.i = a[i__4].i;
#line 503 "cgghd3.f"
			i__4 = k + (j + i__ + 1) * a_dim1;
#line 503 "cgghd3.f"
			temp1.r = a[i__4].r, temp1.i = a[i__4].i;
#line 504 "cgghd3.f"
			i__4 = k + (j + i__ + 2) * a_dim1;
#line 504 "cgghd3.f"
			temp2.r = a[i__4].r, temp2.i = a[i__4].i;
#line 505 "cgghd3.f"
			i__4 = k + (j + i__ + 3) * a_dim1;
#line 505 "cgghd3.f"
			temp3.r = a[i__4].r, temp3.i = a[i__4].i;
#line 506 "cgghd3.f"
			i__4 = k + (j + i__ + 3) * a_dim1;
#line 506 "cgghd3.f"
			z__2.r = c2.r * temp3.r - c2.i * temp3.i, z__2.i = 
				c2.r * temp3.i + c2.i * temp3.r;
#line 506 "cgghd3.f"
			d_cnjg(&z__4, &s2);
#line 506 "cgghd3.f"
			z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, z__3.i =
				 z__4.r * temp2.i + z__4.i * temp2.r;
#line 506 "cgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 506 "cgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 507 "cgghd3.f"
			z__3.r = -s2.r, z__3.i = -s2.i;
#line 507 "cgghd3.f"
			z__2.r = z__3.r * temp3.r - z__3.i * temp3.i, z__2.i =
				 z__3.r * temp3.i + z__3.i * temp3.r;
#line 507 "cgghd3.f"
			z__4.r = c2.r * temp2.r - c2.i * temp2.i, z__4.i = 
				c2.r * temp2.i + c2.i * temp2.r;
#line 507 "cgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 507 "cgghd3.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 508 "cgghd3.f"
			i__4 = k + (j + i__ + 2) * a_dim1;
#line 508 "cgghd3.f"
			z__2.r = c1.r * temp2.r - c1.i * temp2.i, z__2.i = 
				c1.r * temp2.i + c1.i * temp2.r;
#line 508 "cgghd3.f"
			d_cnjg(&z__4, &s1);
#line 508 "cgghd3.f"
			z__3.r = z__4.r * temp1.r - z__4.i * temp1.i, z__3.i =
				 z__4.r * temp1.i + z__4.i * temp1.r;
#line 508 "cgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 508 "cgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 509 "cgghd3.f"
			z__3.r = -s1.r, z__3.i = -s1.i;
#line 509 "cgghd3.f"
			z__2.r = z__3.r * temp2.r - z__3.i * temp2.i, z__2.i =
				 z__3.r * temp2.i + z__3.i * temp2.r;
#line 509 "cgghd3.f"
			z__4.r = c1.r * temp1.r - c1.i * temp1.i, z__4.i = 
				c1.r * temp1.i + c1.i * temp1.r;
#line 509 "cgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 509 "cgghd3.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 510 "cgghd3.f"
			i__4 = k + (j + i__ + 1) * a_dim1;
#line 510 "cgghd3.f"
			z__2.r = ctemp.r * temp1.r - ctemp.i * temp1.i, 
				z__2.i = ctemp.r * temp1.i + ctemp.i * 
				temp1.r;
#line 510 "cgghd3.f"
			d_cnjg(&z__4, &s);
#line 510 "cgghd3.f"
			z__3.r = z__4.r * temp.r - z__4.i * temp.i, z__3.i = 
				z__4.r * temp.i + z__4.i * temp.r;
#line 510 "cgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 510 "cgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 511 "cgghd3.f"
			i__4 = k + (j + i__) * a_dim1;
#line 511 "cgghd3.f"
			z__3.r = -s.r, z__3.i = -s.i;
#line 511 "cgghd3.f"
			z__2.r = z__3.r * temp1.r - z__3.i * temp1.i, z__2.i =
				 z__3.r * temp1.i + z__3.i * temp1.r;
#line 511 "cgghd3.f"
			z__4.r = ctemp.r * temp.r - ctemp.i * temp.i, z__4.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 511 "cgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 511 "cgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 512 "cgghd3.f"
		    }
#line 513 "cgghd3.f"
		}

#line 515 "cgghd3.f"
		if (jj > 0) {
#line 516 "cgghd3.f"
		    for (i__ = jj; i__ >= 1; --i__) {
#line 517 "cgghd3.f"
			i__5 = j + 1 + i__ + j * a_dim1;
#line 517 "cgghd3.f"
			c__ = a[i__5].r;
#line 518 "cgghd3.f"
			i__5 = *ihi - top;
#line 518 "cgghd3.f"
			d_cnjg(&z__2, &b[j + 1 + i__ + j * b_dim1]);
#line 518 "cgghd3.f"
			z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 518 "cgghd3.f"
			crot_(&i__5, &a[top + 1 + (j + i__ + 1) * a_dim1], &
				c__1, &a[top + 1 + (j + i__) * a_dim1], &c__1,
				 &c__, &z__1);
#line 521 "cgghd3.f"
		    }
#line 522 "cgghd3.f"
		}

/*              Update (J+1)th column of A by transformations from left. */

#line 526 "cgghd3.f"
		if (j < jcol + nnb - 1) {
#line 527 "cgghd3.f"
		    len = j + 1 - jcol;

/*                 Multiply with the trailing accumulated unitary */
/*                 matrix, which takes the form */

/*                        [  U11  U12  ] */
/*                    U = [            ], */
/*                        [  U21  U22  ] */

/*                 where U21 is a LEN-by-LEN matrix and U12 is lower */
/*                 triangular. */

#line 539 "cgghd3.f"
		    jrow = *ihi - nblst + 1;
#line 540 "cgghd3.f"
		    cgemv_("Conjugate", &nblst, &len, &c_b1, &work[1], &nblst,
			     &a[jrow + (j + 1) * a_dim1], &c__1, &c_b2, &work[
			    pw], &c__1, (ftnlen)9);
#line 543 "cgghd3.f"
		    ppw = pw + len;
#line 544 "cgghd3.f"
		    i__5 = jrow + nblst - len - 1;
#line 544 "cgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 545 "cgghd3.f"
			i__6 = ppw;
#line 545 "cgghd3.f"
			i__4 = i__ + (j + 1) * a_dim1;
#line 545 "cgghd3.f"
			work[i__6].r = a[i__4].r, work[i__6].i = a[i__4].i;
#line 546 "cgghd3.f"
			++ppw;
#line 547 "cgghd3.f"
		    }
#line 548 "cgghd3.f"
		    i__5 = nblst - len;
#line 548 "cgghd3.f"
		    ctrmv_("Lower", "Conjugate", "Non-unit", &i__5, &work[len 
			    * nblst + 1], &nblst, &work[pw + len], &c__1, (
			    ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 551 "cgghd3.f"
		    i__5 = nblst - len;
#line 551 "cgghd3.f"
		    cgemv_("Conjugate", &len, &i__5, &c_b1, &work[(len + 1) * 
			    nblst - len + 1], &nblst, &a[jrow + nblst - len + 
			    (j + 1) * a_dim1], &c__1, &c_b1, &work[pw + len], 
			    &c__1, (ftnlen)9);
#line 555 "cgghd3.f"
		    ppw = pw;
#line 556 "cgghd3.f"
		    i__5 = jrow + nblst - 1;
#line 556 "cgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 557 "cgghd3.f"
			i__6 = i__ + (j + 1) * a_dim1;
#line 557 "cgghd3.f"
			i__4 = ppw;
#line 557 "cgghd3.f"
			a[i__6].r = work[i__4].r, a[i__6].i = work[i__4].i;
#line 558 "cgghd3.f"
			++ppw;
#line 559 "cgghd3.f"
		    }

/*                 Multiply with the other accumulated unitary */
/*                 matrices, which take the form */

/*                        [  U11  U12   0  ] */
/*                        [                ] */
/*                    U = [  U21  U22   0  ], */
/*                        [                ] */
/*                        [   0    0    I  ] */

/*                 where I denotes the (NNB-LEN)-by-(NNB-LEN) identity */
/*                 matrix, U21 is a LEN-by-LEN upper triangular matrix */
/*                 and U12 is an NNB-by-NNB lower triangular matrix. */

#line 574 "cgghd3.f"
		    ppwo = nblst * nblst + 1;
#line 575 "cgghd3.f"
		    j0 = jrow - nnb;
#line 576 "cgghd3.f"
		    i__5 = jcol + 1;
#line 576 "cgghd3.f"
		    i__6 = -nnb;
#line 576 "cgghd3.f"
		    for (jrow = j0; i__6 < 0 ? jrow >= i__5 : jrow <= i__5; 
			    jrow += i__6) {
#line 577 "cgghd3.f"
			ppw = pw + len;
#line 578 "cgghd3.f"
			i__4 = jrow + nnb - 1;
#line 578 "cgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 579 "cgghd3.f"
			    i__7 = ppw;
#line 579 "cgghd3.f"
			    i__8 = i__ + (j + 1) * a_dim1;
#line 579 "cgghd3.f"
			    work[i__7].r = a[i__8].r, work[i__7].i = a[i__8]
				    .i;
#line 580 "cgghd3.f"
			    ++ppw;
#line 581 "cgghd3.f"
			}
#line 582 "cgghd3.f"
			ppw = pw;
#line 583 "cgghd3.f"
			i__4 = jrow + nnb + len - 1;
#line 583 "cgghd3.f"
			for (i__ = jrow + nnb; i__ <= i__4; ++i__) {
#line 584 "cgghd3.f"
			    i__7 = ppw;
#line 584 "cgghd3.f"
			    i__8 = i__ + (j + 1) * a_dim1;
#line 584 "cgghd3.f"
			    work[i__7].r = a[i__8].r, work[i__7].i = a[i__8]
				    .i;
#line 585 "cgghd3.f"
			    ++ppw;
#line 586 "cgghd3.f"
			}
#line 587 "cgghd3.f"
			i__4 = nnb << 1;
#line 587 "cgghd3.f"
			ctrmv_("Upper", "Conjugate", "Non-unit", &len, &work[
				ppwo + nnb], &i__4, &work[pw], &c__1, (ftnlen)
				5, (ftnlen)9, (ftnlen)8);
#line 590 "cgghd3.f"
			i__4 = nnb << 1;
#line 590 "cgghd3.f"
			ctrmv_("Lower", "Conjugate", "Non-unit", &nnb, &work[
				ppwo + (len << 1) * nnb], &i__4, &work[pw + 
				len], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 593 "cgghd3.f"
			i__4 = nnb << 1;
#line 593 "cgghd3.f"
			cgemv_("Conjugate", &nnb, &len, &c_b1, &work[ppwo], &
				i__4, &a[jrow + (j + 1) * a_dim1], &c__1, &
				c_b1, &work[pw], &c__1, (ftnlen)9);
#line 596 "cgghd3.f"
			i__4 = nnb << 1;
#line 596 "cgghd3.f"
			cgemv_("Conjugate", &len, &nnb, &c_b1, &work[ppwo + (
				len << 1) * nnb + nnb], &i__4, &a[jrow + nnb 
				+ (j + 1) * a_dim1], &c__1, &c_b1, &work[pw + 
				len], &c__1, (ftnlen)9);
#line 600 "cgghd3.f"
			ppw = pw;
#line 601 "cgghd3.f"
			i__4 = jrow + len + nnb - 1;
#line 601 "cgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 602 "cgghd3.f"
			    i__7 = i__ + (j + 1) * a_dim1;
#line 602 "cgghd3.f"
			    i__8 = ppw;
#line 602 "cgghd3.f"
			    a[i__7].r = work[i__8].r, a[i__7].i = work[i__8]
				    .i;
#line 603 "cgghd3.f"
			    ++ppw;
#line 604 "cgghd3.f"
			}
#line 605 "cgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 606 "cgghd3.f"
		    }
#line 607 "cgghd3.f"
		}
#line 608 "cgghd3.f"
	    }

/*           Apply accumulated unitary matrices to A. */

#line 612 "cgghd3.f"
	    cola = *n - jcol - nnb + 1;
#line 613 "cgghd3.f"
	    j = *ihi - nblst + 1;
#line 614 "cgghd3.f"
	    cgemm_("Conjugate", "No Transpose", &nblst, &cola, &nblst, &c_b1, 
		    &work[1], &nblst, &a[j + (jcol + nnb) * a_dim1], lda, &
		    c_b2, &work[pw], &nblst, (ftnlen)9, (ftnlen)12);
#line 618 "cgghd3.f"
	    clacpy_("All", &nblst, &cola, &work[pw], &nblst, &a[j + (jcol + 
		    nnb) * a_dim1], lda, (ftnlen)3);
#line 620 "cgghd3.f"
	    ppwo = nblst * nblst + 1;
#line 621 "cgghd3.f"
	    j0 = j - nnb;
#line 622 "cgghd3.f"
	    i__3 = jcol + 1;
#line 622 "cgghd3.f"
	    i__6 = -nnb;
#line 622 "cgghd3.f"
	    for (j = j0; i__6 < 0 ? j >= i__3 : j <= i__3; j += i__6) {
#line 623 "cgghd3.f"
		if (blk22) {

/*                 Exploit the structure of */

/*                        [  U11  U12  ] */
/*                    U = [            ] */
/*                        [  U21  U22  ], */

/*                 where all blocks are NNB-by-NNB, U21 is upper */
/*                 triangular and U12 is lower triangular. */

#line 634 "cgghd3.f"
		    i__5 = nnb << 1;
#line 634 "cgghd3.f"
		    i__4 = nnb << 1;
#line 634 "cgghd3.f"
		    i__7 = *lwork - pw + 1;
#line 634 "cgghd3.f"
		    cunm22_("Left", "Conjugate", &i__5, &cola, &nnb, &nnb, &
			    work[ppwo], &i__4, &a[j + (jcol + nnb) * a_dim1], 
			    lda, &work[pw], &i__7, &ierr, (ftnlen)4, (ftnlen)
			    9);
#line 638 "cgghd3.f"
		} else {

/*                 Ignore the structure of U. */

#line 642 "cgghd3.f"
		    i__5 = nnb << 1;
#line 642 "cgghd3.f"
		    i__4 = nnb << 1;
#line 642 "cgghd3.f"
		    i__7 = nnb << 1;
#line 642 "cgghd3.f"
		    i__8 = nnb << 1;
#line 642 "cgghd3.f"
		    cgemm_("Conjugate", "No Transpose", &i__5, &cola, &i__4, &
			    c_b1, &work[ppwo], &i__7, &a[j + (jcol + nnb) * 
			    a_dim1], lda, &c_b2, &work[pw], &i__8, (ftnlen)9, 
			    (ftnlen)12);
#line 646 "cgghd3.f"
		    i__5 = nnb << 1;
#line 646 "cgghd3.f"
		    i__4 = nnb << 1;
#line 646 "cgghd3.f"
		    clacpy_("All", &i__5, &cola, &work[pw], &i__4, &a[j + (
			    jcol + nnb) * a_dim1], lda, (ftnlen)3);
#line 648 "cgghd3.f"
		}
#line 649 "cgghd3.f"
		ppwo += (nnb << 2) * nnb;
#line 650 "cgghd3.f"
	    }

/*           Apply accumulated unitary matrices to Q. */

#line 654 "cgghd3.f"
	    if (wantq) {
#line 655 "cgghd3.f"
		j = *ihi - nblst + 1;
#line 656 "cgghd3.f"
		if (initq) {
/* Computing MAX */
#line 657 "cgghd3.f"
		    i__6 = 2, i__3 = j - jcol + 1;
#line 657 "cgghd3.f"
		    topq = max(i__6,i__3);
#line 658 "cgghd3.f"
		    nh = *ihi - topq + 1;
#line 659 "cgghd3.f"
		} else {
#line 660 "cgghd3.f"
		    topq = 1;
#line 661 "cgghd3.f"
		    nh = *n;
#line 662 "cgghd3.f"
		}
#line 663 "cgghd3.f"
		cgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b1, &q[topq + j * q_dim1], ldq, &work[1], &nblst, &
			c_b2, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 666 "cgghd3.f"
		clacpy_("All", &nh, &nblst, &work[pw], &nh, &q[topq + j * 
			q_dim1], ldq, (ftnlen)3);
#line 668 "cgghd3.f"
		ppwo = nblst * nblst + 1;
#line 669 "cgghd3.f"
		j0 = j - nnb;
#line 670 "cgghd3.f"
		i__6 = jcol + 1;
#line 670 "cgghd3.f"
		i__3 = -nnb;
#line 670 "cgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__6 : j <= i__6; j += i__3) {
#line 671 "cgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 672 "cgghd3.f"
			i__5 = 2, i__4 = j - jcol + 1;
#line 672 "cgghd3.f"
			topq = max(i__5,i__4);
#line 673 "cgghd3.f"
			nh = *ihi - topq + 1;
#line 674 "cgghd3.f"
		    }
#line 675 "cgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 679 "cgghd3.f"
			i__5 = nnb << 1;
#line 679 "cgghd3.f"
			i__4 = nnb << 1;
#line 679 "cgghd3.f"
			i__7 = *lwork - pw + 1;
#line 679 "cgghd3.f"
			cunm22_("Right", "No Transpose", &nh, &i__5, &nnb, &
				nnb, &work[ppwo], &i__4, &q[topq + j * q_dim1]
				, ldq, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 683 "cgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 687 "cgghd3.f"
			i__5 = nnb << 1;
#line 687 "cgghd3.f"
			i__4 = nnb << 1;
#line 687 "cgghd3.f"
			i__7 = nnb << 1;
#line 687 "cgghd3.f"
			cgemm_("No Transpose", "No Transpose", &nh, &i__5, &
				i__4, &c_b1, &q[topq + j * q_dim1], ldq, &
				work[ppwo], &i__7, &c_b2, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 691 "cgghd3.f"
			i__5 = nnb << 1;
#line 691 "cgghd3.f"
			clacpy_("All", &nh, &i__5, &work[pw], &nh, &q[topq + 
				j * q_dim1], ldq, (ftnlen)3);
#line 693 "cgghd3.f"
		    }
#line 694 "cgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 695 "cgghd3.f"
		}
#line 696 "cgghd3.f"
	    }

/*           Accumulate right Givens rotations if required. */

#line 700 "cgghd3.f"
	    if (wantz || top > 0) {

/*              Initialize small unitary factors that will hold the */
/*              accumulated Givens rotations in workspace. */

#line 705 "cgghd3.f"
		claset_("All", &nblst, &nblst, &c_b2, &c_b1, &work[1], &nblst,
			 (ftnlen)3);
#line 707 "cgghd3.f"
		pw = nblst * nblst + 1;
#line 708 "cgghd3.f"
		i__3 = n2nb;
#line 708 "cgghd3.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 709 "cgghd3.f"
		    i__6 = nnb << 1;
#line 709 "cgghd3.f"
		    i__5 = nnb << 1;
#line 709 "cgghd3.f"
		    i__4 = nnb << 1;
#line 709 "cgghd3.f"
		    claset_("All", &i__6, &i__5, &c_b2, &c_b1, &work[pw], &
			    i__4, (ftnlen)3);
#line 711 "cgghd3.f"
		    pw += (nnb << 2) * nnb;
#line 712 "cgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 716 "cgghd3.f"
		i__3 = jcol + nnb - 1;
#line 716 "cgghd3.f"
		for (j = jcol; j <= i__3; ++j) {
#line 717 "cgghd3.f"
		    ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 718 "cgghd3.f"
		    len = j + 2 - jcol;
#line 719 "cgghd3.f"
		    jrow = j + n2nb * nnb + 2;
#line 720 "cgghd3.f"
		    i__6 = jrow;
#line 720 "cgghd3.f"
		    for (i__ = *ihi; i__ >= i__6; --i__) {
#line 721 "cgghd3.f"
			i__5 = i__ + j * a_dim1;
#line 721 "cgghd3.f"
			ctemp.r = a[i__5].r, ctemp.i = a[i__5].i;
#line 722 "cgghd3.f"
			i__5 = i__ + j * a_dim1;
#line 722 "cgghd3.f"
			a[i__5].r = 0., a[i__5].i = 0.;
#line 723 "cgghd3.f"
			i__5 = i__ + j * b_dim1;
#line 723 "cgghd3.f"
			s.r = b[i__5].r, s.i = b[i__5].i;
#line 724 "cgghd3.f"
			i__5 = i__ + j * b_dim1;
#line 724 "cgghd3.f"
			b[i__5].r = 0., b[i__5].i = 0.;
#line 725 "cgghd3.f"
			i__5 = ppw + len - 1;
#line 725 "cgghd3.f"
			for (jj = ppw; jj <= i__5; ++jj) {
#line 726 "cgghd3.f"
			    i__4 = jj + nblst;
#line 726 "cgghd3.f"
			    temp.r = work[i__4].r, temp.i = work[i__4].i;
#line 727 "cgghd3.f"
			    i__4 = jj + nblst;
#line 727 "cgghd3.f"
			    z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
				    z__2.i = ctemp.r * temp.i + ctemp.i * 
				    temp.r;
#line 727 "cgghd3.f"
			    d_cnjg(&z__4, &s);
#line 727 "cgghd3.f"
			    i__7 = jj;
#line 727 "cgghd3.f"
			    z__3.r = z__4.r * work[i__7].r - z__4.i * work[
				    i__7].i, z__3.i = z__4.r * work[i__7].i + 
				    z__4.i * work[i__7].r;
#line 727 "cgghd3.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 727 "cgghd3.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 729 "cgghd3.f"
			    i__4 = jj;
#line 729 "cgghd3.f"
			    z__2.r = s.r * temp.r - s.i * temp.i, z__2.i = 
				    s.r * temp.i + s.i * temp.r;
#line 729 "cgghd3.f"
			    i__7 = jj;
#line 729 "cgghd3.f"
			    z__3.r = ctemp.r * work[i__7].r - ctemp.i * work[
				    i__7].i, z__3.i = ctemp.r * work[i__7].i 
				    + ctemp.i * work[i__7].r;
#line 729 "cgghd3.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 729 "cgghd3.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 730 "cgghd3.f"
			}
#line 731 "cgghd3.f"
			++len;
#line 732 "cgghd3.f"
			ppw = ppw - nblst - 1;
#line 733 "cgghd3.f"
		    }

#line 735 "cgghd3.f"
		    ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + 
			    nnb;
#line 736 "cgghd3.f"
		    j0 = jrow - nnb;
#line 737 "cgghd3.f"
		    i__6 = j + 2;
#line 737 "cgghd3.f"
		    i__5 = -nnb;
#line 737 "cgghd3.f"
		    for (jrow = j0; i__5 < 0 ? jrow >= i__6 : jrow <= i__6; 
			    jrow += i__5) {
#line 738 "cgghd3.f"
			ppw = ppwo;
#line 739 "cgghd3.f"
			len = j + 2 - jcol;
#line 740 "cgghd3.f"
			i__4 = jrow;
#line 740 "cgghd3.f"
			for (i__ = jrow + nnb - 1; i__ >= i__4; --i__) {
#line 741 "cgghd3.f"
			    i__7 = i__ + j * a_dim1;
#line 741 "cgghd3.f"
			    ctemp.r = a[i__7].r, ctemp.i = a[i__7].i;
#line 742 "cgghd3.f"
			    i__7 = i__ + j * a_dim1;
#line 742 "cgghd3.f"
			    a[i__7].r = 0., a[i__7].i = 0.;
#line 743 "cgghd3.f"
			    i__7 = i__ + j * b_dim1;
#line 743 "cgghd3.f"
			    s.r = b[i__7].r, s.i = b[i__7].i;
#line 744 "cgghd3.f"
			    i__7 = i__ + j * b_dim1;
#line 744 "cgghd3.f"
			    b[i__7].r = 0., b[i__7].i = 0.;
#line 745 "cgghd3.f"
			    i__7 = ppw + len - 1;
#line 745 "cgghd3.f"
			    for (jj = ppw; jj <= i__7; ++jj) {
#line 746 "cgghd3.f"
				i__8 = jj + (nnb << 1);
#line 746 "cgghd3.f"
				temp.r = work[i__8].r, temp.i = work[i__8].i;
#line 747 "cgghd3.f"
				i__8 = jj + (nnb << 1);
#line 747 "cgghd3.f"
				z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
					z__2.i = ctemp.r * temp.i + ctemp.i * 
					temp.r;
#line 747 "cgghd3.f"
				d_cnjg(&z__4, &s);
#line 747 "cgghd3.f"
				i__9 = jj;
#line 747 "cgghd3.f"
				z__3.r = z__4.r * work[i__9].r - z__4.i * 
					work[i__9].i, z__3.i = z__4.r * work[
					i__9].i + z__4.i * work[i__9].r;
#line 747 "cgghd3.f"
				z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
					z__3.i;
#line 747 "cgghd3.f"
				work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 749 "cgghd3.f"
				i__8 = jj;
#line 749 "cgghd3.f"
				z__2.r = s.r * temp.r - s.i * temp.i, z__2.i =
					 s.r * temp.i + s.i * temp.r;
#line 749 "cgghd3.f"
				i__9 = jj;
#line 749 "cgghd3.f"
				z__3.r = ctemp.r * work[i__9].r - ctemp.i * 
					work[i__9].i, z__3.i = ctemp.r * work[
					i__9].i + ctemp.i * work[i__9].r;
#line 749 "cgghd3.f"
				z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
					z__3.i;
#line 749 "cgghd3.f"
				work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 750 "cgghd3.f"
			    }
#line 751 "cgghd3.f"
			    ++len;
#line 752 "cgghd3.f"
			    ppw = ppw - (nnb << 1) - 1;
#line 753 "cgghd3.f"
			}
#line 754 "cgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 755 "cgghd3.f"
		    }
#line 756 "cgghd3.f"
		}
#line 757 "cgghd3.f"
	    } else {

#line 759 "cgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 759 "cgghd3.f"
		claset_("Lower", &i__3, &nnb, &c_b2, &c_b2, &a[jcol + 2 + 
			jcol * a_dim1], lda, (ftnlen)5);
#line 761 "cgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 761 "cgghd3.f"
		claset_("Lower", &i__3, &nnb, &c_b2, &c_b2, &b[jcol + 2 + 
			jcol * b_dim1], ldb, (ftnlen)5);
#line 763 "cgghd3.f"
	    }

/*           Apply accumulated unitary matrices to A and B. */

#line 767 "cgghd3.f"
	    if (top > 0) {
#line 768 "cgghd3.f"
		j = *ihi - nblst + 1;
#line 769 "cgghd3.f"
		cgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b1, &a[j * a_dim1 + 1], lda, &work[1], &nblst, &
			c_b2, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 772 "cgghd3.f"
		clacpy_("All", &top, &nblst, &work[pw], &top, &a[j * a_dim1 + 
			1], lda, (ftnlen)3);
#line 774 "cgghd3.f"
		ppwo = nblst * nblst + 1;
#line 775 "cgghd3.f"
		j0 = j - nnb;
#line 776 "cgghd3.f"
		i__3 = jcol + 1;
#line 776 "cgghd3.f"
		i__5 = -nnb;
#line 776 "cgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 777 "cgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 781 "cgghd3.f"
			i__6 = nnb << 1;
#line 781 "cgghd3.f"
			i__4 = nnb << 1;
#line 781 "cgghd3.f"
			i__7 = *lwork - pw + 1;
#line 781 "cgghd3.f"
			cunm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &a[j * a_dim1 + 1], 
				lda, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 785 "cgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 789 "cgghd3.f"
			i__6 = nnb << 1;
#line 789 "cgghd3.f"
			i__4 = nnb << 1;
#line 789 "cgghd3.f"
			i__7 = nnb << 1;
#line 789 "cgghd3.f"
			cgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b1, &a[j * a_dim1 + 1], lda, &work[
				ppwo], &i__7, &c_b2, &work[pw], &top, (ftnlen)
				12, (ftnlen)12);
#line 793 "cgghd3.f"
			i__6 = nnb << 1;
#line 793 "cgghd3.f"
			clacpy_("All", &top, &i__6, &work[pw], &top, &a[j * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 795 "cgghd3.f"
		    }
#line 796 "cgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 797 "cgghd3.f"
		}

#line 799 "cgghd3.f"
		j = *ihi - nblst + 1;
#line 800 "cgghd3.f"
		cgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b1, &b[j * b_dim1 + 1], ldb, &work[1], &nblst, &
			c_b2, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 803 "cgghd3.f"
		clacpy_("All", &top, &nblst, &work[pw], &top, &b[j * b_dim1 + 
			1], ldb, (ftnlen)3);
#line 805 "cgghd3.f"
		ppwo = nblst * nblst + 1;
#line 806 "cgghd3.f"
		j0 = j - nnb;
#line 807 "cgghd3.f"
		i__5 = jcol + 1;
#line 807 "cgghd3.f"
		i__3 = -nnb;
#line 807 "cgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__5 : j <= i__5; j += i__3) {
#line 808 "cgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 812 "cgghd3.f"
			i__6 = nnb << 1;
#line 812 "cgghd3.f"
			i__4 = nnb << 1;
#line 812 "cgghd3.f"
			i__7 = *lwork - pw + 1;
#line 812 "cgghd3.f"
			cunm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &b[j * b_dim1 + 1], 
				ldb, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 816 "cgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 820 "cgghd3.f"
			i__6 = nnb << 1;
#line 820 "cgghd3.f"
			i__4 = nnb << 1;
#line 820 "cgghd3.f"
			i__7 = nnb << 1;
#line 820 "cgghd3.f"
			cgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b1, &b[j * b_dim1 + 1], ldb, &work[
				ppwo], &i__7, &c_b2, &work[pw], &top, (ftnlen)
				12, (ftnlen)12);
#line 824 "cgghd3.f"
			i__6 = nnb << 1;
#line 824 "cgghd3.f"
			clacpy_("All", &top, &i__6, &work[pw], &top, &b[j * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 826 "cgghd3.f"
		    }
#line 827 "cgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 828 "cgghd3.f"
		}
#line 829 "cgghd3.f"
	    }

/*           Apply accumulated unitary matrices to Z. */

#line 833 "cgghd3.f"
	    if (wantz) {
#line 834 "cgghd3.f"
		j = *ihi - nblst + 1;
#line 835 "cgghd3.f"
		if (initq) {
/* Computing MAX */
#line 836 "cgghd3.f"
		    i__3 = 2, i__5 = j - jcol + 1;
#line 836 "cgghd3.f"
		    topq = max(i__3,i__5);
#line 837 "cgghd3.f"
		    nh = *ihi - topq + 1;
#line 838 "cgghd3.f"
		} else {
#line 839 "cgghd3.f"
		    topq = 1;
#line 840 "cgghd3.f"
		    nh = *n;
#line 841 "cgghd3.f"
		}
#line 842 "cgghd3.f"
		cgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b1, &z__[topq + j * z_dim1], ldz, &work[1], &nblst, 
			&c_b2, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 845 "cgghd3.f"
		clacpy_("All", &nh, &nblst, &work[pw], &nh, &z__[topq + j * 
			z_dim1], ldz, (ftnlen)3);
#line 847 "cgghd3.f"
		ppwo = nblst * nblst + 1;
#line 848 "cgghd3.f"
		j0 = j - nnb;
#line 849 "cgghd3.f"
		i__3 = jcol + 1;
#line 849 "cgghd3.f"
		i__5 = -nnb;
#line 849 "cgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 850 "cgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 851 "cgghd3.f"
			i__6 = 2, i__4 = j - jcol + 1;
#line 851 "cgghd3.f"
			topq = max(i__6,i__4);
#line 852 "cgghd3.f"
			nh = *ihi - topq + 1;
#line 853 "cgghd3.f"
		    }
#line 854 "cgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 858 "cgghd3.f"
			i__6 = nnb << 1;
#line 858 "cgghd3.f"
			i__4 = nnb << 1;
#line 858 "cgghd3.f"
			i__7 = *lwork - pw + 1;
#line 858 "cgghd3.f"
			cunm22_("Right", "No Transpose", &nh, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &z__[topq + j * 
				z_dim1], ldz, &work[pw], &i__7, &ierr, (
				ftnlen)5, (ftnlen)12);
#line 862 "cgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 866 "cgghd3.f"
			i__6 = nnb << 1;
#line 866 "cgghd3.f"
			i__4 = nnb << 1;
#line 866 "cgghd3.f"
			i__7 = nnb << 1;
#line 866 "cgghd3.f"
			cgemm_("No Transpose", "No Transpose", &nh, &i__6, &
				i__4, &c_b1, &z__[topq + j * z_dim1], ldz, &
				work[ppwo], &i__7, &c_b2, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 870 "cgghd3.f"
			i__6 = nnb << 1;
#line 870 "cgghd3.f"
			clacpy_("All", &nh, &i__6, &work[pw], &nh, &z__[topq 
				+ j * z_dim1], ldz, (ftnlen)3);
#line 872 "cgghd3.f"
		    }
#line 873 "cgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 874 "cgghd3.f"
		}
#line 875 "cgghd3.f"
	    }
#line 876 "cgghd3.f"
	}
#line 877 "cgghd3.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */
/*     Avoid re-initialization of modified Q and Z. */

#line 882 "cgghd3.f"
    *(unsigned char *)compq2 = *(unsigned char *)compq;
#line 883 "cgghd3.f"
    *(unsigned char *)compz2 = *(unsigned char *)compz;
#line 884 "cgghd3.f"
    if (jcol != *ilo) {
#line 885 "cgghd3.f"
	if (wantq) {
#line 885 "cgghd3.f"
	    *(unsigned char *)compq2 = 'V';
#line 885 "cgghd3.f"
	}
#line 887 "cgghd3.f"
	if (wantz) {
#line 887 "cgghd3.f"
	    *(unsigned char *)compz2 = 'V';
#line 887 "cgghd3.f"
	}
#line 889 "cgghd3.f"
    }

#line 891 "cgghd3.f"
    if (jcol < *ihi) {
#line 891 "cgghd3.f"
	cgghrd_(compq2, compz2, n, &jcol, ihi, &a[a_offset], lda, &b[b_offset]
		, ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &ierr, (ftnlen)
		1, (ftnlen)1);
#line 891 "cgghd3.f"
    }
#line 894 "cgghd3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 894 "cgghd3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 896 "cgghd3.f"
    return 0;

/*     End of CGGHD3 */

} /* cgghd3_ */

