#line 1 "zgghd3.f"
/* zgghd3.f -- translated by f2c (version 20100827).
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

#line 1 "zgghd3.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__16 = 16;

/* > \brief \b ZGGHD3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGHD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgghd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgghd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgghd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGHD3 reduces a pair of complex matrices (A,B) to generalized upper */
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
/* >      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H */
/* >      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H */
/* > If Q1 is the unitary matrix from the QR factorization of B in the */
/* > original equation A*x = lambda*B*x, then ZGGHD3 reduces the original */
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
/* >          normally set by a previous call to ZGGBAL; otherwise they */
/* >          should be set to 1 and N respectively. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB, N) */
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
/* >          Q is COMPLEX*16 array, dimension (LDQ, N) */
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
/* >          Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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

/* > \ingroup complex16OTHERcomputational */

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
/* Subroutine */ int zgghd3_(char *compq, char *compz, integer *n, integer *
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
    static integer jrow, topq, ppwo;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublecomplex temp1, temp2, temp3;
    static integer kacc22;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    static doublecomplex ctemp;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static integer nblst;
    static logical initq;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical wantq, initz;
    extern /* Subroutine */ int zunm22_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen)
	    ;
    static logical wantz;
    extern /* Subroutine */ int ztrmv_(char *, char *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static char compq2[1], compz2[1];
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), zlartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), zlacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen);
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.6.1) -- */
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

#line 278 "zgghd3.f"
    /* Parameter adjustments */
#line 278 "zgghd3.f"
    a_dim1 = *lda;
#line 278 "zgghd3.f"
    a_offset = 1 + a_dim1;
#line 278 "zgghd3.f"
    a -= a_offset;
#line 278 "zgghd3.f"
    b_dim1 = *ldb;
#line 278 "zgghd3.f"
    b_offset = 1 + b_dim1;
#line 278 "zgghd3.f"
    b -= b_offset;
#line 278 "zgghd3.f"
    q_dim1 = *ldq;
#line 278 "zgghd3.f"
    q_offset = 1 + q_dim1;
#line 278 "zgghd3.f"
    q -= q_offset;
#line 278 "zgghd3.f"
    z_dim1 = *ldz;
#line 278 "zgghd3.f"
    z_offset = 1 + z_dim1;
#line 278 "zgghd3.f"
    z__ -= z_offset;
#line 278 "zgghd3.f"
    --work;
#line 278 "zgghd3.f"

#line 278 "zgghd3.f"
    /* Function Body */
#line 278 "zgghd3.f"
    *info = 0;
#line 279 "zgghd3.f"
    nb = ilaenv_(&c__1, "ZGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)
	    1);
/* Computing MAX */
#line 280 "zgghd3.f"
    i__1 = *n * 6 * nb;
#line 280 "zgghd3.f"
    lwkopt = max(i__1,1);
#line 281 "zgghd3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 281 "zgghd3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 282 "zgghd3.f"
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 283 "zgghd3.f"
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 284 "zgghd3.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 285 "zgghd3.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);
#line 286 "zgghd3.f"
    lquery = *lwork == -1;

#line 288 "zgghd3.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 289 "zgghd3.f"
	*info = -1;
#line 290 "zgghd3.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 291 "zgghd3.f"
	*info = -2;
#line 292 "zgghd3.f"
    } else if (*n < 0) {
#line 293 "zgghd3.f"
	*info = -3;
#line 294 "zgghd3.f"
    } else if (*ilo < 1) {
#line 295 "zgghd3.f"
	*info = -4;
#line 296 "zgghd3.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 297 "zgghd3.f"
	*info = -5;
#line 298 "zgghd3.f"
    } else if (*lda < max(1,*n)) {
#line 299 "zgghd3.f"
	*info = -7;
#line 300 "zgghd3.f"
    } else if (*ldb < max(1,*n)) {
#line 301 "zgghd3.f"
	*info = -9;
#line 302 "zgghd3.f"
    } else if (wantq && *ldq < *n || *ldq < 1) {
#line 303 "zgghd3.f"
	*info = -11;
#line 304 "zgghd3.f"
    } else if (wantz && *ldz < *n || *ldz < 1) {
#line 305 "zgghd3.f"
	*info = -13;
#line 306 "zgghd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 307 "zgghd3.f"
	*info = -15;
#line 308 "zgghd3.f"
    }
#line 309 "zgghd3.f"
    if (*info != 0) {
#line 310 "zgghd3.f"
	i__1 = -(*info);
#line 310 "zgghd3.f"
	xerbla_("ZGGHD3", &i__1, (ftnlen)6);
#line 311 "zgghd3.f"
	return 0;
#line 312 "zgghd3.f"
    } else if (lquery) {
#line 313 "zgghd3.f"
	return 0;
#line 314 "zgghd3.f"
    }

/*     Initialize Q and Z if desired. */

#line 318 "zgghd3.f"
    if (initq) {
#line 318 "zgghd3.f"
	zlaset_("All", n, n, &c_b2, &c_b1, &q[q_offset], ldq, (ftnlen)3);
#line 318 "zgghd3.f"
    }
#line 320 "zgghd3.f"
    if (initz) {
#line 320 "zgghd3.f"
	zlaset_("All", n, n, &c_b2, &c_b1, &z__[z_offset], ldz, (ftnlen)3);
#line 320 "zgghd3.f"
    }

/*     Zero out lower triangle of B. */

#line 325 "zgghd3.f"
    if (*n > 1) {
#line 325 "zgghd3.f"
	i__1 = *n - 1;
#line 325 "zgghd3.f"
	i__2 = *n - 1;
#line 325 "zgghd3.f"
	zlaset_("Lower", &i__1, &i__2, &c_b2, &c_b2, &b[b_dim1 + 2], ldb, (
		ftnlen)5);
#line 325 "zgghd3.f"
    }

/*     Quick return if possible */

#line 330 "zgghd3.f"
    nh = *ihi - *ilo + 1;
#line 331 "zgghd3.f"
    if (nh <= 1) {
#line 332 "zgghd3.f"
	work[1].r = 1., work[1].i = 0.;
#line 333 "zgghd3.f"
	return 0;
#line 334 "zgghd3.f"
    }

/*     Determine the blocksize. */

#line 338 "zgghd3.f"
    nbmin = ilaenv_(&c__2, "ZGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 339 "zgghd3.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to use unblocked instead of blocked code. */

/* Computing MAX */
#line 343 "zgghd3.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "ZGGHD3", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 343 "zgghd3.f"
	nx = max(i__1,i__2);
#line 344 "zgghd3.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code. */

#line 348 "zgghd3.f"
	    if (*lwork < lwkopt) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code. */

/* Computing MAX */
#line 354 "zgghd3.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "ZGGHD3", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 354 "zgghd3.f"
		nbmin = max(i__1,i__2);
#line 356 "zgghd3.f"
		if (*lwork >= *n * 6 * nbmin) {
#line 357 "zgghd3.f"
		    nb = *lwork / (*n * 6);
#line 358 "zgghd3.f"
		} else {
#line 359 "zgghd3.f"
		    nb = 1;
#line 360 "zgghd3.f"
		}
#line 361 "zgghd3.f"
	    }
#line 362 "zgghd3.f"
	}
#line 363 "zgghd3.f"
    }

#line 365 "zgghd3.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 369 "zgghd3.f"
	jcol = *ilo;

#line 371 "zgghd3.f"
    } else {

/*        Use blocked code */

#line 375 "zgghd3.f"
	kacc22 = ilaenv_(&c__16, "ZGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 376 "zgghd3.f"
	blk22 = kacc22 == 2;
#line 377 "zgghd3.f"
	i__1 = *ihi - 2;
#line 377 "zgghd3.f"
	i__2 = nb;
#line 377 "zgghd3.f"
	for (jcol = *ilo; i__2 < 0 ? jcol >= i__1 : jcol <= i__1; jcol += 
		i__2) {
/* Computing MIN */
#line 378 "zgghd3.f"
	    i__3 = nb, i__4 = *ihi - jcol - 1;
#line 378 "zgghd3.f"
	    nnb = min(i__3,i__4);

/*           Initialize small unitary factors that will hold the */
/*           accumulated Givens rotations in workspace. */
/*           N2NB   denotes the number of 2*NNB-by-2*NNB factors */
/*           NBLST  denotes the (possibly smaller) order of the last */
/*                  factor. */

#line 386 "zgghd3.f"
	    n2nb = (*ihi - jcol - 1) / nnb - 1;
#line 387 "zgghd3.f"
	    nblst = *ihi - jcol - n2nb * nnb;
#line 388 "zgghd3.f"
	    zlaset_("All", &nblst, &nblst, &c_b2, &c_b1, &work[1], &nblst, (
		    ftnlen)3);
#line 389 "zgghd3.f"
	    pw = nblst * nblst + 1;
#line 390 "zgghd3.f"
	    i__3 = n2nb;
#line 390 "zgghd3.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 391 "zgghd3.f"
		i__4 = nnb << 1;
#line 391 "zgghd3.f"
		i__5 = nnb << 1;
#line 391 "zgghd3.f"
		i__6 = nnb << 1;
#line 391 "zgghd3.f"
		zlaset_("All", &i__4, &i__5, &c_b2, &c_b1, &work[pw], &i__6, (
			ftnlen)3);
#line 393 "zgghd3.f"
		pw += (nnb << 2) * nnb;
#line 394 "zgghd3.f"
	    }

/*           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form. */

#line 398 "zgghd3.f"
	    i__3 = jcol + nnb - 1;
#line 398 "zgghd3.f"
	    for (j = jcol; j <= i__3; ++j) {

/*              Reduce Jth column of A. Store cosines and sines in Jth */
/*              column of A and B, respectively. */

#line 403 "zgghd3.f"
		i__4 = j + 2;
#line 403 "zgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 404 "zgghd3.f"
		    i__5 = i__ - 1 + j * a_dim1;
#line 404 "zgghd3.f"
		    temp.r = a[i__5].r, temp.i = a[i__5].i;
#line 405 "zgghd3.f"
		    zlartg_(&temp, &a[i__ + j * a_dim1], &c__, &s, &a[i__ - 1 
			    + j * a_dim1]);
#line 406 "zgghd3.f"
		    i__5 = i__ + j * a_dim1;
#line 406 "zgghd3.f"
		    z__1.r = c__, z__1.i = 0.;
#line 406 "zgghd3.f"
		    a[i__5].r = z__1.r, a[i__5].i = z__1.i;
#line 407 "zgghd3.f"
		    i__5 = i__ + j * b_dim1;
#line 407 "zgghd3.f"
		    b[i__5].r = s.r, b[i__5].i = s.i;
#line 408 "zgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 412 "zgghd3.f"
		ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 413 "zgghd3.f"
		len = j + 2 - jcol;
#line 414 "zgghd3.f"
		jrow = j + n2nb * nnb + 2;
#line 415 "zgghd3.f"
		i__4 = jrow;
#line 415 "zgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 416 "zgghd3.f"
		    i__5 = i__ + j * a_dim1;
#line 416 "zgghd3.f"
		    ctemp.r = a[i__5].r, ctemp.i = a[i__5].i;
#line 417 "zgghd3.f"
		    i__5 = i__ + j * b_dim1;
#line 417 "zgghd3.f"
		    s.r = b[i__5].r, s.i = b[i__5].i;
#line 418 "zgghd3.f"
		    i__5 = ppw + len - 1;
#line 418 "zgghd3.f"
		    for (jj = ppw; jj <= i__5; ++jj) {
#line 419 "zgghd3.f"
			i__6 = jj + nblst;
#line 419 "zgghd3.f"
			temp.r = work[i__6].r, temp.i = work[i__6].i;
#line 420 "zgghd3.f"
			i__6 = jj + nblst;
#line 420 "zgghd3.f"
			z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, z__2.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 420 "zgghd3.f"
			i__7 = jj;
#line 420 "zgghd3.f"
			z__3.r = s.r * work[i__7].r - s.i * work[i__7].i, 
				z__3.i = s.r * work[i__7].i + s.i * work[i__7]
				.r;
#line 420 "zgghd3.f"
			z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 420 "zgghd3.f"
			work[i__6].r = z__1.r, work[i__6].i = z__1.i;
#line 421 "zgghd3.f"
			i__6 = jj;
#line 421 "zgghd3.f"
			d_cnjg(&z__3, &s);
#line 421 "zgghd3.f"
			z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
				z__3.r * temp.i + z__3.i * temp.r;
#line 421 "zgghd3.f"
			i__7 = jj;
#line 421 "zgghd3.f"
			z__4.r = ctemp.r * work[i__7].r - ctemp.i * work[i__7]
				.i, z__4.i = ctemp.r * work[i__7].i + ctemp.i 
				* work[i__7].r;
#line 421 "zgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 421 "zgghd3.f"
			work[i__6].r = z__1.r, work[i__6].i = z__1.i;
#line 422 "zgghd3.f"
		    }
#line 423 "zgghd3.f"
		    ++len;
#line 424 "zgghd3.f"
		    ppw = ppw - nblst - 1;
#line 425 "zgghd3.f"
		}

#line 427 "zgghd3.f"
		ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + nnb;
#line 428 "zgghd3.f"
		j0 = jrow - nnb;
#line 429 "zgghd3.f"
		i__4 = j + 2;
#line 429 "zgghd3.f"
		i__5 = -nnb;
#line 429 "zgghd3.f"
		for (jrow = j0; i__5 < 0 ? jrow >= i__4 : jrow <= i__4; jrow 
			+= i__5) {
#line 430 "zgghd3.f"
		    ppw = ppwo;
#line 431 "zgghd3.f"
		    len = j + 2 - jcol;
#line 432 "zgghd3.f"
		    i__6 = jrow;
#line 432 "zgghd3.f"
		    for (i__ = jrow + nnb - 1; i__ >= i__6; --i__) {
#line 433 "zgghd3.f"
			i__7 = i__ + j * a_dim1;
#line 433 "zgghd3.f"
			ctemp.r = a[i__7].r, ctemp.i = a[i__7].i;
#line 434 "zgghd3.f"
			i__7 = i__ + j * b_dim1;
#line 434 "zgghd3.f"
			s.r = b[i__7].r, s.i = b[i__7].i;
#line 435 "zgghd3.f"
			i__7 = ppw + len - 1;
#line 435 "zgghd3.f"
			for (jj = ppw; jj <= i__7; ++jj) {
#line 436 "zgghd3.f"
			    i__8 = jj + (nnb << 1);
#line 436 "zgghd3.f"
			    temp.r = work[i__8].r, temp.i = work[i__8].i;
#line 437 "zgghd3.f"
			    i__8 = jj + (nnb << 1);
#line 437 "zgghd3.f"
			    z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
				    z__2.i = ctemp.r * temp.i + ctemp.i * 
				    temp.r;
#line 437 "zgghd3.f"
			    i__9 = jj;
#line 437 "zgghd3.f"
			    z__3.r = s.r * work[i__9].r - s.i * work[i__9].i, 
				    z__3.i = s.r * work[i__9].i + s.i * work[
				    i__9].r;
#line 437 "zgghd3.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 437 "zgghd3.f"
			    work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 438 "zgghd3.f"
			    i__8 = jj;
#line 438 "zgghd3.f"
			    d_cnjg(&z__3, &s);
#line 438 "zgghd3.f"
			    z__2.r = z__3.r * temp.r - z__3.i * temp.i, 
				    z__2.i = z__3.r * temp.i + z__3.i * 
				    temp.r;
#line 438 "zgghd3.f"
			    i__9 = jj;
#line 438 "zgghd3.f"
			    z__4.r = ctemp.r * work[i__9].r - ctemp.i * work[
				    i__9].i, z__4.i = ctemp.r * work[i__9].i 
				    + ctemp.i * work[i__9].r;
#line 438 "zgghd3.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 438 "zgghd3.f"
			    work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 439 "zgghd3.f"
			}
#line 440 "zgghd3.f"
			++len;
#line 441 "zgghd3.f"
			ppw = ppw - (nnb << 1) - 1;
#line 442 "zgghd3.f"
		    }
#line 443 "zgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 444 "zgghd3.f"
		}

/*              TOP denotes the number of top rows in A and B that will */
/*              not be updated during the next steps. */

#line 449 "zgghd3.f"
		if (jcol <= 2) {
#line 450 "zgghd3.f"
		    top = 0;
#line 451 "zgghd3.f"
		} else {
#line 452 "zgghd3.f"
		    top = jcol;
#line 453 "zgghd3.f"
		}

/*              Propagate transformations through B and replace stored */
/*              left sines/cosines by right sines/cosines. */

#line 458 "zgghd3.f"
		i__5 = j + 1;
#line 458 "zgghd3.f"
		for (jj = *n; jj >= i__5; --jj) {

/*                 Update JJth column of B. */

/* Computing MIN */
#line 462 "zgghd3.f"
		    i__4 = jj + 1;
#line 462 "zgghd3.f"
		    i__6 = j + 2;
#line 462 "zgghd3.f"
		    for (i__ = min(i__4,*ihi); i__ >= i__6; --i__) {
#line 463 "zgghd3.f"
			i__4 = i__ + j * a_dim1;
#line 463 "zgghd3.f"
			ctemp.r = a[i__4].r, ctemp.i = a[i__4].i;
#line 464 "zgghd3.f"
			i__4 = i__ + j * b_dim1;
#line 464 "zgghd3.f"
			s.r = b[i__4].r, s.i = b[i__4].i;
#line 465 "zgghd3.f"
			i__4 = i__ + jj * b_dim1;
#line 465 "zgghd3.f"
			temp.r = b[i__4].r, temp.i = b[i__4].i;
#line 466 "zgghd3.f"
			i__4 = i__ + jj * b_dim1;
#line 466 "zgghd3.f"
			z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, z__2.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 466 "zgghd3.f"
			d_cnjg(&z__4, &s);
#line 466 "zgghd3.f"
			i__7 = i__ - 1 + jj * b_dim1;
#line 466 "zgghd3.f"
			z__3.r = z__4.r * b[i__7].r - z__4.i * b[i__7].i, 
				z__3.i = z__4.r * b[i__7].i + z__4.i * b[i__7]
				.r;
#line 466 "zgghd3.f"
			z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 466 "zgghd3.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 467 "zgghd3.f"
			i__4 = i__ - 1 + jj * b_dim1;
#line 467 "zgghd3.f"
			z__2.r = s.r * temp.r - s.i * temp.i, z__2.i = s.r * 
				temp.i + s.i * temp.r;
#line 467 "zgghd3.f"
			i__7 = i__ - 1 + jj * b_dim1;
#line 467 "zgghd3.f"
			z__3.r = ctemp.r * b[i__7].r - ctemp.i * b[i__7].i, 
				z__3.i = ctemp.r * b[i__7].i + ctemp.i * b[
				i__7].r;
#line 467 "zgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 467 "zgghd3.f"
			b[i__4].r = z__1.r, b[i__4].i = z__1.i;
#line 468 "zgghd3.f"
		    }

/*                 Annihilate B( JJ+1, JJ ). */

#line 472 "zgghd3.f"
		    if (jj < *ihi) {
#line 473 "zgghd3.f"
			i__6 = jj + 1 + (jj + 1) * b_dim1;
#line 473 "zgghd3.f"
			temp.r = b[i__6].r, temp.i = b[i__6].i;
#line 474 "zgghd3.f"
			zlartg_(&temp, &b[jj + 1 + jj * b_dim1], &c__, &s, &b[
				jj + 1 + (jj + 1) * b_dim1]);
#line 476 "zgghd3.f"
			i__6 = jj + 1 + jj * b_dim1;
#line 476 "zgghd3.f"
			b[i__6].r = 0., b[i__6].i = 0.;
#line 477 "zgghd3.f"
			i__6 = jj - top;
#line 477 "zgghd3.f"
			zrot_(&i__6, &b[top + 1 + (jj + 1) * b_dim1], &c__1, &
				b[top + 1 + jj * b_dim1], &c__1, &c__, &s);
#line 479 "zgghd3.f"
			i__6 = jj + 1 + j * a_dim1;
#line 479 "zgghd3.f"
			z__1.r = c__, z__1.i = 0.;
#line 479 "zgghd3.f"
			a[i__6].r = z__1.r, a[i__6].i = z__1.i;
#line 480 "zgghd3.f"
			i__6 = jj + 1 + j * b_dim1;
#line 480 "zgghd3.f"
			d_cnjg(&z__2, &s);
#line 480 "zgghd3.f"
			z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 480 "zgghd3.f"
			b[i__6].r = z__1.r, b[i__6].i = z__1.i;
#line 481 "zgghd3.f"
		    }
#line 482 "zgghd3.f"
		}

/*              Update A by transformations from right. */

#line 486 "zgghd3.f"
		jj = (*ihi - j - 1) % 3;
#line 487 "zgghd3.f"
		i__5 = jj + 1;
#line 487 "zgghd3.f"
		for (i__ = *ihi - j - 3; i__ >= i__5; i__ += -3) {
#line 488 "zgghd3.f"
		    i__6 = j + 1 + i__ + j * a_dim1;
#line 488 "zgghd3.f"
		    ctemp.r = a[i__6].r, ctemp.i = a[i__6].i;
#line 489 "zgghd3.f"
		    i__6 = j + 1 + i__ + j * b_dim1;
#line 489 "zgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 489 "zgghd3.f"
		    s.r = z__1.r, s.i = z__1.i;
#line 490 "zgghd3.f"
		    i__6 = j + 2 + i__ + j * a_dim1;
#line 490 "zgghd3.f"
		    c1.r = a[i__6].r, c1.i = a[i__6].i;
#line 491 "zgghd3.f"
		    i__6 = j + 2 + i__ + j * b_dim1;
#line 491 "zgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 491 "zgghd3.f"
		    s1.r = z__1.r, s1.i = z__1.i;
#line 492 "zgghd3.f"
		    i__6 = j + 3 + i__ + j * a_dim1;
#line 492 "zgghd3.f"
		    c2.r = a[i__6].r, c2.i = a[i__6].i;
#line 493 "zgghd3.f"
		    i__6 = j + 3 + i__ + j * b_dim1;
#line 493 "zgghd3.f"
		    z__1.r = -b[i__6].r, z__1.i = -b[i__6].i;
#line 493 "zgghd3.f"
		    s2.r = z__1.r, s2.i = z__1.i;

#line 495 "zgghd3.f"
		    i__6 = *ihi;
#line 495 "zgghd3.f"
		    for (k = top + 1; k <= i__6; ++k) {
#line 496 "zgghd3.f"
			i__4 = k + (j + i__) * a_dim1;
#line 496 "zgghd3.f"
			temp.r = a[i__4].r, temp.i = a[i__4].i;
#line 497 "zgghd3.f"
			i__4 = k + (j + i__ + 1) * a_dim1;
#line 497 "zgghd3.f"
			temp1.r = a[i__4].r, temp1.i = a[i__4].i;
#line 498 "zgghd3.f"
			i__4 = k + (j + i__ + 2) * a_dim1;
#line 498 "zgghd3.f"
			temp2.r = a[i__4].r, temp2.i = a[i__4].i;
#line 499 "zgghd3.f"
			i__4 = k + (j + i__ + 3) * a_dim1;
#line 499 "zgghd3.f"
			temp3.r = a[i__4].r, temp3.i = a[i__4].i;
#line 500 "zgghd3.f"
			i__4 = k + (j + i__ + 3) * a_dim1;
#line 500 "zgghd3.f"
			z__2.r = c2.r * temp3.r - c2.i * temp3.i, z__2.i = 
				c2.r * temp3.i + c2.i * temp3.r;
#line 500 "zgghd3.f"
			d_cnjg(&z__4, &s2);
#line 500 "zgghd3.f"
			z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, z__3.i =
				 z__4.r * temp2.i + z__4.i * temp2.r;
#line 500 "zgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 500 "zgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 501 "zgghd3.f"
			z__3.r = -s2.r, z__3.i = -s2.i;
#line 501 "zgghd3.f"
			z__2.r = z__3.r * temp3.r - z__3.i * temp3.i, z__2.i =
				 z__3.r * temp3.i + z__3.i * temp3.r;
#line 501 "zgghd3.f"
			z__4.r = c2.r * temp2.r - c2.i * temp2.i, z__4.i = 
				c2.r * temp2.i + c2.i * temp2.r;
#line 501 "zgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 501 "zgghd3.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 502 "zgghd3.f"
			i__4 = k + (j + i__ + 2) * a_dim1;
#line 502 "zgghd3.f"
			z__2.r = c1.r * temp2.r - c1.i * temp2.i, z__2.i = 
				c1.r * temp2.i + c1.i * temp2.r;
#line 502 "zgghd3.f"
			d_cnjg(&z__4, &s1);
#line 502 "zgghd3.f"
			z__3.r = z__4.r * temp1.r - z__4.i * temp1.i, z__3.i =
				 z__4.r * temp1.i + z__4.i * temp1.r;
#line 502 "zgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 502 "zgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 503 "zgghd3.f"
			z__3.r = -s1.r, z__3.i = -s1.i;
#line 503 "zgghd3.f"
			z__2.r = z__3.r * temp2.r - z__3.i * temp2.i, z__2.i =
				 z__3.r * temp2.i + z__3.i * temp2.r;
#line 503 "zgghd3.f"
			z__4.r = c1.r * temp1.r - c1.i * temp1.i, z__4.i = 
				c1.r * temp1.i + c1.i * temp1.r;
#line 503 "zgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 503 "zgghd3.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 504 "zgghd3.f"
			i__4 = k + (j + i__ + 1) * a_dim1;
#line 504 "zgghd3.f"
			z__2.r = ctemp.r * temp1.r - ctemp.i * temp1.i, 
				z__2.i = ctemp.r * temp1.i + ctemp.i * 
				temp1.r;
#line 504 "zgghd3.f"
			d_cnjg(&z__4, &s);
#line 504 "zgghd3.f"
			z__3.r = z__4.r * temp.r - z__4.i * temp.i, z__3.i = 
				z__4.r * temp.i + z__4.i * temp.r;
#line 504 "zgghd3.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 504 "zgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 505 "zgghd3.f"
			i__4 = k + (j + i__) * a_dim1;
#line 505 "zgghd3.f"
			z__3.r = -s.r, z__3.i = -s.i;
#line 505 "zgghd3.f"
			z__2.r = z__3.r * temp1.r - z__3.i * temp1.i, z__2.i =
				 z__3.r * temp1.i + z__3.i * temp1.r;
#line 505 "zgghd3.f"
			z__4.r = ctemp.r * temp.r - ctemp.i * temp.i, z__4.i =
				 ctemp.r * temp.i + ctemp.i * temp.r;
#line 505 "zgghd3.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 505 "zgghd3.f"
			a[i__4].r = z__1.r, a[i__4].i = z__1.i;
#line 506 "zgghd3.f"
		    }
#line 507 "zgghd3.f"
		}

#line 509 "zgghd3.f"
		if (jj > 0) {
#line 510 "zgghd3.f"
		    for (i__ = jj; i__ >= 1; --i__) {
#line 511 "zgghd3.f"
			i__5 = j + 1 + i__ + j * a_dim1;
#line 511 "zgghd3.f"
			c__ = a[i__5].r;
#line 512 "zgghd3.f"
			i__5 = *ihi - top;
#line 512 "zgghd3.f"
			d_cnjg(&z__2, &b[j + 1 + i__ + j * b_dim1]);
#line 512 "zgghd3.f"
			z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 512 "zgghd3.f"
			zrot_(&i__5, &a[top + 1 + (j + i__ + 1) * a_dim1], &
				c__1, &a[top + 1 + (j + i__) * a_dim1], &c__1,
				 &c__, &z__1);
#line 515 "zgghd3.f"
		    }
#line 516 "zgghd3.f"
		}

/*              Update (J+1)th column of A by transformations from left. */

#line 520 "zgghd3.f"
		if (j < jcol + nnb - 1) {
#line 521 "zgghd3.f"
		    len = j + 1 - jcol;

/*                 Multiply with the trailing accumulated unitary */
/*                 matrix, which takes the form */

/*                        [  U11  U12  ] */
/*                    U = [            ], */
/*                        [  U21  U22  ] */

/*                 where U21 is a LEN-by-LEN matrix and U12 is lower */
/*                 triangular. */

#line 533 "zgghd3.f"
		    jrow = *ihi - nblst + 1;
#line 534 "zgghd3.f"
		    zgemv_("Conjugate", &nblst, &len, &c_b1, &work[1], &nblst,
			     &a[jrow + (j + 1) * a_dim1], &c__1, &c_b2, &work[
			    pw], &c__1, (ftnlen)9);
#line 537 "zgghd3.f"
		    ppw = pw + len;
#line 538 "zgghd3.f"
		    i__5 = jrow + nblst - len - 1;
#line 538 "zgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 539 "zgghd3.f"
			i__6 = ppw;
#line 539 "zgghd3.f"
			i__4 = i__ + (j + 1) * a_dim1;
#line 539 "zgghd3.f"
			work[i__6].r = a[i__4].r, work[i__6].i = a[i__4].i;
#line 540 "zgghd3.f"
			++ppw;
#line 541 "zgghd3.f"
		    }
#line 542 "zgghd3.f"
		    i__5 = nblst - len;
#line 542 "zgghd3.f"
		    ztrmv_("Lower", "Conjugate", "Non-unit", &i__5, &work[len 
			    * nblst + 1], &nblst, &work[pw + len], &c__1, (
			    ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 545 "zgghd3.f"
		    i__5 = nblst - len;
#line 545 "zgghd3.f"
		    zgemv_("Conjugate", &len, &i__5, &c_b1, &work[(len + 1) * 
			    nblst - len + 1], &nblst, &a[jrow + nblst - len + 
			    (j + 1) * a_dim1], &c__1, &c_b1, &work[pw + len], 
			    &c__1, (ftnlen)9);
#line 549 "zgghd3.f"
		    ppw = pw;
#line 550 "zgghd3.f"
		    i__5 = jrow + nblst - 1;
#line 550 "zgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 551 "zgghd3.f"
			i__6 = i__ + (j + 1) * a_dim1;
#line 551 "zgghd3.f"
			i__4 = ppw;
#line 551 "zgghd3.f"
			a[i__6].r = work[i__4].r, a[i__6].i = work[i__4].i;
#line 552 "zgghd3.f"
			++ppw;
#line 553 "zgghd3.f"
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

#line 568 "zgghd3.f"
		    ppwo = nblst * nblst + 1;
#line 569 "zgghd3.f"
		    j0 = jrow - nnb;
#line 570 "zgghd3.f"
		    i__5 = jcol + 1;
#line 570 "zgghd3.f"
		    i__6 = -nnb;
#line 570 "zgghd3.f"
		    for (jrow = j0; i__6 < 0 ? jrow >= i__5 : jrow <= i__5; 
			    jrow += i__6) {
#line 571 "zgghd3.f"
			ppw = pw + len;
#line 572 "zgghd3.f"
			i__4 = jrow + nnb - 1;
#line 572 "zgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 573 "zgghd3.f"
			    i__7 = ppw;
#line 573 "zgghd3.f"
			    i__8 = i__ + (j + 1) * a_dim1;
#line 573 "zgghd3.f"
			    work[i__7].r = a[i__8].r, work[i__7].i = a[i__8]
				    .i;
#line 574 "zgghd3.f"
			    ++ppw;
#line 575 "zgghd3.f"
			}
#line 576 "zgghd3.f"
			ppw = pw;
#line 577 "zgghd3.f"
			i__4 = jrow + nnb + len - 1;
#line 577 "zgghd3.f"
			for (i__ = jrow + nnb; i__ <= i__4; ++i__) {
#line 578 "zgghd3.f"
			    i__7 = ppw;
#line 578 "zgghd3.f"
			    i__8 = i__ + (j + 1) * a_dim1;
#line 578 "zgghd3.f"
			    work[i__7].r = a[i__8].r, work[i__7].i = a[i__8]
				    .i;
#line 579 "zgghd3.f"
			    ++ppw;
#line 580 "zgghd3.f"
			}
#line 581 "zgghd3.f"
			i__4 = nnb << 1;
#line 581 "zgghd3.f"
			ztrmv_("Upper", "Conjugate", "Non-unit", &len, &work[
				ppwo + nnb], &i__4, &work[pw], &c__1, (ftnlen)
				5, (ftnlen)9, (ftnlen)8);
#line 584 "zgghd3.f"
			i__4 = nnb << 1;
#line 584 "zgghd3.f"
			ztrmv_("Lower", "Conjugate", "Non-unit", &nnb, &work[
				ppwo + (len << 1) * nnb], &i__4, &work[pw + 
				len], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 587 "zgghd3.f"
			i__4 = nnb << 1;
#line 587 "zgghd3.f"
			zgemv_("Conjugate", &nnb, &len, &c_b1, &work[ppwo], &
				i__4, &a[jrow + (j + 1) * a_dim1], &c__1, &
				c_b1, &work[pw], &c__1, (ftnlen)9);
#line 590 "zgghd3.f"
			i__4 = nnb << 1;
#line 590 "zgghd3.f"
			zgemv_("Conjugate", &len, &nnb, &c_b1, &work[ppwo + (
				len << 1) * nnb + nnb], &i__4, &a[jrow + nnb 
				+ (j + 1) * a_dim1], &c__1, &c_b1, &work[pw + 
				len], &c__1, (ftnlen)9);
#line 594 "zgghd3.f"
			ppw = pw;
#line 595 "zgghd3.f"
			i__4 = jrow + len + nnb - 1;
#line 595 "zgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 596 "zgghd3.f"
			    i__7 = i__ + (j + 1) * a_dim1;
#line 596 "zgghd3.f"
			    i__8 = ppw;
#line 596 "zgghd3.f"
			    a[i__7].r = work[i__8].r, a[i__7].i = work[i__8]
				    .i;
#line 597 "zgghd3.f"
			    ++ppw;
#line 598 "zgghd3.f"
			}
#line 599 "zgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 600 "zgghd3.f"
		    }
#line 601 "zgghd3.f"
		}
#line 602 "zgghd3.f"
	    }

/*           Apply accumulated unitary matrices to A. */

#line 606 "zgghd3.f"
	    cola = *n - jcol - nnb + 1;
#line 607 "zgghd3.f"
	    j = *ihi - nblst + 1;
#line 608 "zgghd3.f"
	    zgemm_("Conjugate", "No Transpose", &nblst, &cola, &nblst, &c_b1, 
		    &work[1], &nblst, &a[j + (jcol + nnb) * a_dim1], lda, &
		    c_b2, &work[pw], &nblst, (ftnlen)9, (ftnlen)12);
#line 612 "zgghd3.f"
	    zlacpy_("All", &nblst, &cola, &work[pw], &nblst, &a[j + (jcol + 
		    nnb) * a_dim1], lda, (ftnlen)3);
#line 614 "zgghd3.f"
	    ppwo = nblst * nblst + 1;
#line 615 "zgghd3.f"
	    j0 = j - nnb;
#line 616 "zgghd3.f"
	    i__3 = jcol + 1;
#line 616 "zgghd3.f"
	    i__6 = -nnb;
#line 616 "zgghd3.f"
	    for (j = j0; i__6 < 0 ? j >= i__3 : j <= i__3; j += i__6) {
#line 617 "zgghd3.f"
		if (blk22) {

/*                 Exploit the structure of */

/*                        [  U11  U12  ] */
/*                    U = [            ] */
/*                        [  U21  U22  ], */

/*                 where all blocks are NNB-by-NNB, U21 is upper */
/*                 triangular and U12 is lower triangular. */

#line 628 "zgghd3.f"
		    i__5 = nnb << 1;
#line 628 "zgghd3.f"
		    i__4 = nnb << 1;
#line 628 "zgghd3.f"
		    i__7 = *lwork - pw + 1;
#line 628 "zgghd3.f"
		    zunm22_("Left", "Conjugate", &i__5, &cola, &nnb, &nnb, &
			    work[ppwo], &i__4, &a[j + (jcol + nnb) * a_dim1], 
			    lda, &work[pw], &i__7, &ierr, (ftnlen)4, (ftnlen)
			    9);
#line 632 "zgghd3.f"
		} else {

/*                 Ignore the structure of U. */

#line 636 "zgghd3.f"
		    i__5 = nnb << 1;
#line 636 "zgghd3.f"
		    i__4 = nnb << 1;
#line 636 "zgghd3.f"
		    i__7 = nnb << 1;
#line 636 "zgghd3.f"
		    i__8 = nnb << 1;
#line 636 "zgghd3.f"
		    zgemm_("Conjugate", "No Transpose", &i__5, &cola, &i__4, &
			    c_b1, &work[ppwo], &i__7, &a[j + (jcol + nnb) * 
			    a_dim1], lda, &c_b2, &work[pw], &i__8, (ftnlen)9, 
			    (ftnlen)12);
#line 640 "zgghd3.f"
		    i__5 = nnb << 1;
#line 640 "zgghd3.f"
		    i__4 = nnb << 1;
#line 640 "zgghd3.f"
		    zlacpy_("All", &i__5, &cola, &work[pw], &i__4, &a[j + (
			    jcol + nnb) * a_dim1], lda, (ftnlen)3);
#line 642 "zgghd3.f"
		}
#line 643 "zgghd3.f"
		ppwo += (nnb << 2) * nnb;
#line 644 "zgghd3.f"
	    }

/*           Apply accumulated unitary matrices to Q. */

#line 648 "zgghd3.f"
	    if (wantq) {
#line 649 "zgghd3.f"
		j = *ihi - nblst + 1;
#line 650 "zgghd3.f"
		if (initq) {
/* Computing MAX */
#line 651 "zgghd3.f"
		    i__6 = 2, i__3 = j - jcol + 1;
#line 651 "zgghd3.f"
		    topq = max(i__6,i__3);
#line 652 "zgghd3.f"
		    nh = *ihi - topq + 1;
#line 653 "zgghd3.f"
		} else {
#line 654 "zgghd3.f"
		    topq = 1;
#line 655 "zgghd3.f"
		    nh = *n;
#line 656 "zgghd3.f"
		}
#line 657 "zgghd3.f"
		zgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b1, &q[topq + j * q_dim1], ldq, &work[1], &nblst, &
			c_b2, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 660 "zgghd3.f"
		zlacpy_("All", &nh, &nblst, &work[pw], &nh, &q[topq + j * 
			q_dim1], ldq, (ftnlen)3);
#line 662 "zgghd3.f"
		ppwo = nblst * nblst + 1;
#line 663 "zgghd3.f"
		j0 = j - nnb;
#line 664 "zgghd3.f"
		i__6 = jcol + 1;
#line 664 "zgghd3.f"
		i__3 = -nnb;
#line 664 "zgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__6 : j <= i__6; j += i__3) {
#line 665 "zgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 666 "zgghd3.f"
			i__5 = 2, i__4 = j - jcol + 1;
#line 666 "zgghd3.f"
			topq = max(i__5,i__4);
#line 667 "zgghd3.f"
			nh = *ihi - topq + 1;
#line 668 "zgghd3.f"
		    }
#line 669 "zgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 673 "zgghd3.f"
			i__5 = nnb << 1;
#line 673 "zgghd3.f"
			i__4 = nnb << 1;
#line 673 "zgghd3.f"
			i__7 = *lwork - pw + 1;
#line 673 "zgghd3.f"
			zunm22_("Right", "No Transpose", &nh, &i__5, &nnb, &
				nnb, &work[ppwo], &i__4, &q[topq + j * q_dim1]
				, ldq, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 677 "zgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 681 "zgghd3.f"
			i__5 = nnb << 1;
#line 681 "zgghd3.f"
			i__4 = nnb << 1;
#line 681 "zgghd3.f"
			i__7 = nnb << 1;
#line 681 "zgghd3.f"
			zgemm_("No Transpose", "No Transpose", &nh, &i__5, &
				i__4, &c_b1, &q[topq + j * q_dim1], ldq, &
				work[ppwo], &i__7, &c_b2, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 685 "zgghd3.f"
			i__5 = nnb << 1;
#line 685 "zgghd3.f"
			zlacpy_("All", &nh, &i__5, &work[pw], &nh, &q[topq + 
				j * q_dim1], ldq, (ftnlen)3);
#line 687 "zgghd3.f"
		    }
#line 688 "zgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 689 "zgghd3.f"
		}
#line 690 "zgghd3.f"
	    }

/*           Accumulate right Givens rotations if required. */

#line 694 "zgghd3.f"
	    if (wantz || top > 0) {

/*              Initialize small unitary factors that will hold the */
/*              accumulated Givens rotations in workspace. */

#line 699 "zgghd3.f"
		zlaset_("All", &nblst, &nblst, &c_b2, &c_b1, &work[1], &nblst,
			 (ftnlen)3);
#line 701 "zgghd3.f"
		pw = nblst * nblst + 1;
#line 702 "zgghd3.f"
		i__3 = n2nb;
#line 702 "zgghd3.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 703 "zgghd3.f"
		    i__6 = nnb << 1;
#line 703 "zgghd3.f"
		    i__5 = nnb << 1;
#line 703 "zgghd3.f"
		    i__4 = nnb << 1;
#line 703 "zgghd3.f"
		    zlaset_("All", &i__6, &i__5, &c_b2, &c_b1, &work[pw], &
			    i__4, (ftnlen)3);
#line 705 "zgghd3.f"
		    pw += (nnb << 2) * nnb;
#line 706 "zgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 710 "zgghd3.f"
		i__3 = jcol + nnb - 1;
#line 710 "zgghd3.f"
		for (j = jcol; j <= i__3; ++j) {
#line 711 "zgghd3.f"
		    ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 712 "zgghd3.f"
		    len = j + 2 - jcol;
#line 713 "zgghd3.f"
		    jrow = j + n2nb * nnb + 2;
#line 714 "zgghd3.f"
		    i__6 = jrow;
#line 714 "zgghd3.f"
		    for (i__ = *ihi; i__ >= i__6; --i__) {
#line 715 "zgghd3.f"
			i__5 = i__ + j * a_dim1;
#line 715 "zgghd3.f"
			ctemp.r = a[i__5].r, ctemp.i = a[i__5].i;
#line 716 "zgghd3.f"
			i__5 = i__ + j * a_dim1;
#line 716 "zgghd3.f"
			a[i__5].r = 0., a[i__5].i = 0.;
#line 717 "zgghd3.f"
			i__5 = i__ + j * b_dim1;
#line 717 "zgghd3.f"
			s.r = b[i__5].r, s.i = b[i__5].i;
#line 718 "zgghd3.f"
			i__5 = i__ + j * b_dim1;
#line 718 "zgghd3.f"
			b[i__5].r = 0., b[i__5].i = 0.;
#line 719 "zgghd3.f"
			i__5 = ppw + len - 1;
#line 719 "zgghd3.f"
			for (jj = ppw; jj <= i__5; ++jj) {
#line 720 "zgghd3.f"
			    i__4 = jj + nblst;
#line 720 "zgghd3.f"
			    temp.r = work[i__4].r, temp.i = work[i__4].i;
#line 721 "zgghd3.f"
			    i__4 = jj + nblst;
#line 721 "zgghd3.f"
			    z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
				    z__2.i = ctemp.r * temp.i + ctemp.i * 
				    temp.r;
#line 721 "zgghd3.f"
			    d_cnjg(&z__4, &s);
#line 721 "zgghd3.f"
			    i__7 = jj;
#line 721 "zgghd3.f"
			    z__3.r = z__4.r * work[i__7].r - z__4.i * work[
				    i__7].i, z__3.i = z__4.r * work[i__7].i + 
				    z__4.i * work[i__7].r;
#line 721 "zgghd3.f"
			    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
				    z__3.i;
#line 721 "zgghd3.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 723 "zgghd3.f"
			    i__4 = jj;
#line 723 "zgghd3.f"
			    z__2.r = s.r * temp.r - s.i * temp.i, z__2.i = 
				    s.r * temp.i + s.i * temp.r;
#line 723 "zgghd3.f"
			    i__7 = jj;
#line 723 "zgghd3.f"
			    z__3.r = ctemp.r * work[i__7].r - ctemp.i * work[
				    i__7].i, z__3.i = ctemp.r * work[i__7].i 
				    + ctemp.i * work[i__7].r;
#line 723 "zgghd3.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 723 "zgghd3.f"
			    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 724 "zgghd3.f"
			}
#line 725 "zgghd3.f"
			++len;
#line 726 "zgghd3.f"
			ppw = ppw - nblst - 1;
#line 727 "zgghd3.f"
		    }

#line 729 "zgghd3.f"
		    ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + 
			    nnb;
#line 730 "zgghd3.f"
		    j0 = jrow - nnb;
#line 731 "zgghd3.f"
		    i__6 = j + 2;
#line 731 "zgghd3.f"
		    i__5 = -nnb;
#line 731 "zgghd3.f"
		    for (jrow = j0; i__5 < 0 ? jrow >= i__6 : jrow <= i__6; 
			    jrow += i__5) {
#line 732 "zgghd3.f"
			ppw = ppwo;
#line 733 "zgghd3.f"
			len = j + 2 - jcol;
#line 734 "zgghd3.f"
			i__4 = jrow;
#line 734 "zgghd3.f"
			for (i__ = jrow + nnb - 1; i__ >= i__4; --i__) {
#line 735 "zgghd3.f"
			    i__7 = i__ + j * a_dim1;
#line 735 "zgghd3.f"
			    ctemp.r = a[i__7].r, ctemp.i = a[i__7].i;
#line 736 "zgghd3.f"
			    i__7 = i__ + j * a_dim1;
#line 736 "zgghd3.f"
			    a[i__7].r = 0., a[i__7].i = 0.;
#line 737 "zgghd3.f"
			    i__7 = i__ + j * b_dim1;
#line 737 "zgghd3.f"
			    s.r = b[i__7].r, s.i = b[i__7].i;
#line 738 "zgghd3.f"
			    i__7 = i__ + j * b_dim1;
#line 738 "zgghd3.f"
			    b[i__7].r = 0., b[i__7].i = 0.;
#line 739 "zgghd3.f"
			    i__7 = ppw + len - 1;
#line 739 "zgghd3.f"
			    for (jj = ppw; jj <= i__7; ++jj) {
#line 740 "zgghd3.f"
				i__8 = jj + (nnb << 1);
#line 740 "zgghd3.f"
				temp.r = work[i__8].r, temp.i = work[i__8].i;
#line 741 "zgghd3.f"
				i__8 = jj + (nnb << 1);
#line 741 "zgghd3.f"
				z__2.r = ctemp.r * temp.r - ctemp.i * temp.i, 
					z__2.i = ctemp.r * temp.i + ctemp.i * 
					temp.r;
#line 741 "zgghd3.f"
				d_cnjg(&z__4, &s);
#line 741 "zgghd3.f"
				i__9 = jj;
#line 741 "zgghd3.f"
				z__3.r = z__4.r * work[i__9].r - z__4.i * 
					work[i__9].i, z__3.i = z__4.r * work[
					i__9].i + z__4.i * work[i__9].r;
#line 741 "zgghd3.f"
				z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - 
					z__3.i;
#line 741 "zgghd3.f"
				work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 743 "zgghd3.f"
				i__8 = jj;
#line 743 "zgghd3.f"
				z__2.r = s.r * temp.r - s.i * temp.i, z__2.i =
					 s.r * temp.i + s.i * temp.r;
#line 743 "zgghd3.f"
				i__9 = jj;
#line 743 "zgghd3.f"
				z__3.r = ctemp.r * work[i__9].r - ctemp.i * 
					work[i__9].i, z__3.i = ctemp.r * work[
					i__9].i + ctemp.i * work[i__9].r;
#line 743 "zgghd3.f"
				z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
					z__3.i;
#line 743 "zgghd3.f"
				work[i__8].r = z__1.r, work[i__8].i = z__1.i;
#line 744 "zgghd3.f"
			    }
#line 745 "zgghd3.f"
			    ++len;
#line 746 "zgghd3.f"
			    ppw = ppw - (nnb << 1) - 1;
#line 747 "zgghd3.f"
			}
#line 748 "zgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 749 "zgghd3.f"
		    }
#line 750 "zgghd3.f"
		}
#line 751 "zgghd3.f"
	    } else {

#line 753 "zgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 753 "zgghd3.f"
		zlaset_("Lower", &i__3, &nnb, &c_b2, &c_b2, &a[jcol + 2 + 
			jcol * a_dim1], lda, (ftnlen)5);
#line 755 "zgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 755 "zgghd3.f"
		zlaset_("Lower", &i__3, &nnb, &c_b2, &c_b2, &b[jcol + 2 + 
			jcol * b_dim1], ldb, (ftnlen)5);
#line 757 "zgghd3.f"
	    }

/*           Apply accumulated unitary matrices to A and B. */

#line 761 "zgghd3.f"
	    if (top > 0) {
#line 762 "zgghd3.f"
		j = *ihi - nblst + 1;
#line 763 "zgghd3.f"
		zgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b1, &a[j * a_dim1 + 1], lda, &work[1], &nblst, &
			c_b2, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 766 "zgghd3.f"
		zlacpy_("All", &top, &nblst, &work[pw], &top, &a[j * a_dim1 + 
			1], lda, (ftnlen)3);
#line 768 "zgghd3.f"
		ppwo = nblst * nblst + 1;
#line 769 "zgghd3.f"
		j0 = j - nnb;
#line 770 "zgghd3.f"
		i__3 = jcol + 1;
#line 770 "zgghd3.f"
		i__5 = -nnb;
#line 770 "zgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 771 "zgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 775 "zgghd3.f"
			i__6 = nnb << 1;
#line 775 "zgghd3.f"
			i__4 = nnb << 1;
#line 775 "zgghd3.f"
			i__7 = *lwork - pw + 1;
#line 775 "zgghd3.f"
			zunm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &a[j * a_dim1 + 1], 
				lda, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 779 "zgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 783 "zgghd3.f"
			i__6 = nnb << 1;
#line 783 "zgghd3.f"
			i__4 = nnb << 1;
#line 783 "zgghd3.f"
			i__7 = nnb << 1;
#line 783 "zgghd3.f"
			zgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b1, &a[j * a_dim1 + 1], lda, &work[
				ppwo], &i__7, &c_b2, &work[pw], &top, (ftnlen)
				12, (ftnlen)12);
#line 787 "zgghd3.f"
			i__6 = nnb << 1;
#line 787 "zgghd3.f"
			zlacpy_("All", &top, &i__6, &work[pw], &top, &a[j * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 789 "zgghd3.f"
		    }
#line 790 "zgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 791 "zgghd3.f"
		}

#line 793 "zgghd3.f"
		j = *ihi - nblst + 1;
#line 794 "zgghd3.f"
		zgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b1, &b[j * b_dim1 + 1], ldb, &work[1], &nblst, &
			c_b2, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 797 "zgghd3.f"
		zlacpy_("All", &top, &nblst, &work[pw], &top, &b[j * b_dim1 + 
			1], ldb, (ftnlen)3);
#line 799 "zgghd3.f"
		ppwo = nblst * nblst + 1;
#line 800 "zgghd3.f"
		j0 = j - nnb;
#line 801 "zgghd3.f"
		i__5 = jcol + 1;
#line 801 "zgghd3.f"
		i__3 = -nnb;
#line 801 "zgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__5 : j <= i__5; j += i__3) {
#line 802 "zgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 806 "zgghd3.f"
			i__6 = nnb << 1;
#line 806 "zgghd3.f"
			i__4 = nnb << 1;
#line 806 "zgghd3.f"
			i__7 = *lwork - pw + 1;
#line 806 "zgghd3.f"
			zunm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &b[j * b_dim1 + 1], 
				ldb, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 810 "zgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 814 "zgghd3.f"
			i__6 = nnb << 1;
#line 814 "zgghd3.f"
			i__4 = nnb << 1;
#line 814 "zgghd3.f"
			i__7 = nnb << 1;
#line 814 "zgghd3.f"
			zgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b1, &b[j * b_dim1 + 1], ldb, &work[
				ppwo], &i__7, &c_b2, &work[pw], &top, (ftnlen)
				12, (ftnlen)12);
#line 818 "zgghd3.f"
			i__6 = nnb << 1;
#line 818 "zgghd3.f"
			zlacpy_("All", &top, &i__6, &work[pw], &top, &b[j * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 820 "zgghd3.f"
		    }
#line 821 "zgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 822 "zgghd3.f"
		}
#line 823 "zgghd3.f"
	    }

/*           Apply accumulated unitary matrices to Z. */

#line 827 "zgghd3.f"
	    if (wantz) {
#line 828 "zgghd3.f"
		j = *ihi - nblst + 1;
#line 829 "zgghd3.f"
		if (initq) {
/* Computing MAX */
#line 830 "zgghd3.f"
		    i__3 = 2, i__5 = j - jcol + 1;
#line 830 "zgghd3.f"
		    topq = max(i__3,i__5);
#line 831 "zgghd3.f"
		    nh = *ihi - topq + 1;
#line 832 "zgghd3.f"
		} else {
#line 833 "zgghd3.f"
		    topq = 1;
#line 834 "zgghd3.f"
		    nh = *n;
#line 835 "zgghd3.f"
		}
#line 836 "zgghd3.f"
		zgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b1, &z__[topq + j * z_dim1], ldz, &work[1], &nblst, 
			&c_b2, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 839 "zgghd3.f"
		zlacpy_("All", &nh, &nblst, &work[pw], &nh, &z__[topq + j * 
			z_dim1], ldz, (ftnlen)3);
#line 841 "zgghd3.f"
		ppwo = nblst * nblst + 1;
#line 842 "zgghd3.f"
		j0 = j - nnb;
#line 843 "zgghd3.f"
		i__3 = jcol + 1;
#line 843 "zgghd3.f"
		i__5 = -nnb;
#line 843 "zgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 844 "zgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 845 "zgghd3.f"
			i__6 = 2, i__4 = j - jcol + 1;
#line 845 "zgghd3.f"
			topq = max(i__6,i__4);
#line 846 "zgghd3.f"
			nh = *ihi - topq + 1;
#line 847 "zgghd3.f"
		    }
#line 848 "zgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 852 "zgghd3.f"
			i__6 = nnb << 1;
#line 852 "zgghd3.f"
			i__4 = nnb << 1;
#line 852 "zgghd3.f"
			i__7 = *lwork - pw + 1;
#line 852 "zgghd3.f"
			zunm22_("Right", "No Transpose", &nh, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &z__[topq + j * 
				z_dim1], ldz, &work[pw], &i__7, &ierr, (
				ftnlen)5, (ftnlen)12);
#line 856 "zgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 860 "zgghd3.f"
			i__6 = nnb << 1;
#line 860 "zgghd3.f"
			i__4 = nnb << 1;
#line 860 "zgghd3.f"
			i__7 = nnb << 1;
#line 860 "zgghd3.f"
			zgemm_("No Transpose", "No Transpose", &nh, &i__6, &
				i__4, &c_b1, &z__[topq + j * z_dim1], ldz, &
				work[ppwo], &i__7, &c_b2, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 864 "zgghd3.f"
			i__6 = nnb << 1;
#line 864 "zgghd3.f"
			zlacpy_("All", &nh, &i__6, &work[pw], &nh, &z__[topq 
				+ j * z_dim1], ldz, (ftnlen)3);
#line 866 "zgghd3.f"
		    }
#line 867 "zgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 868 "zgghd3.f"
		}
#line 869 "zgghd3.f"
	    }
#line 870 "zgghd3.f"
	}
#line 871 "zgghd3.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */
/*     Avoid re-initialization of modified Q and Z. */

#line 876 "zgghd3.f"
    *(unsigned char *)compq2 = *(unsigned char *)compq;
#line 877 "zgghd3.f"
    *(unsigned char *)compz2 = *(unsigned char *)compz;
#line 878 "zgghd3.f"
    if (jcol != *ilo) {
#line 879 "zgghd3.f"
	if (wantq) {
#line 879 "zgghd3.f"
	    *(unsigned char *)compq2 = 'V';
#line 879 "zgghd3.f"
	}
#line 881 "zgghd3.f"
	if (wantz) {
#line 881 "zgghd3.f"
	    *(unsigned char *)compz2 = 'V';
#line 881 "zgghd3.f"
	}
#line 883 "zgghd3.f"
    }

#line 885 "zgghd3.f"
    if (jcol < *ihi) {
#line 885 "zgghd3.f"
	zgghrd_(compq2, compz2, n, &jcol, ihi, &a[a_offset], lda, &b[b_offset]
		, ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &ierr, (ftnlen)
		1, (ftnlen)1);
#line 885 "zgghd3.f"
    }
#line 888 "zgghd3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 888 "zgghd3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;

#line 890 "zgghd3.f"
    return 0;

/*     End of ZGGHD3 */

} /* zgghd3_ */

