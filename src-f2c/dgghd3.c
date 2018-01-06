#line 1 "dgghd3.f"
/* dgghd3.f -- translated by f2c (version 20100827).
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

#line 1 "dgghd3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 0.;
static doublereal c_b15 = 1.;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__16 = 16;

/* > \brief \b DGGHD3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGHD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGHD3 reduces a pair of real matrices (A,B) to generalized upper */
/* > Hessenberg form using orthogonal transformations, where A is a */
/* > general matrix and B is upper triangular.  The form of the */
/* > generalized eigenvalue problem is */
/* >    A*x = lambda*B*x, */
/* > and B is typically made upper triangular by computing its QR */
/* > factorization and moving the orthogonal matrix Q to the left side */
/* > of the equation. */
/* > */
/* > This subroutine simultaneously reduces A to a Hessenberg matrix H: */
/* >    Q**T*A*Z = H */
/* > and transforms B to another upper triangular matrix T: */
/* >    Q**T*B*Z = T */
/* > in order to reduce the problem to its standard form */
/* >    H*y = lambda*T*y */
/* > where y = Z**T*x. */
/* > */
/* > The orthogonal matrices Q and Z are determined as products of Givens */
/* > rotations.  They may either be formed explicitly, or they may be */
/* > postmultiplied into input matrices Q1 and Z1, so that */
/* > */
/* >      Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T */
/* > */
/* >      Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T */
/* > */
/* > If Q1 is the orthogonal matrix from the QR factorization of B in the */
/* > original equation A*x = lambda*B*x, then DGGHD3 reduces the original */
/* > problem to generalized Hessenberg form. */
/* > */
/* > This is a blocked variant of DGGHRD, using matrix-matrix */
/* > multiplications for parts of the computation to enhance performance. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'N': do not compute Q; */
/* >          = 'I': Q is initialized to the unit matrix, and the */
/* >                 orthogonal matrix Q is returned; */
/* >          = 'V': Q must contain an orthogonal matrix Q1 on entry, */
/* >                 and the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N': do not compute Z; */
/* >          = 'I': Z is initialized to the unit matrix, and the */
/* >                 orthogonal matrix Z is returned; */
/* >          = 'V': Z must contain an orthogonal matrix Z1 on entry, */
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
/* >          normally set by a previous call to DGGBAL; otherwise they */
/* >          should be set to 1 and N respectively. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA, N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB, N) */
/* >          On entry, the N-by-N upper triangular matrix B. */
/* >          On exit, the upper triangular matrix T = Q**T B Z.  The */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* >          On entry, if COMPQ = 'V', the orthogonal matrix Q1, */
/* >          typically from the QR factorization of B. */
/* >          On exit, if COMPQ='I', the orthogonal matrix Q, and if */
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
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          On entry, if COMPZ = 'V', the orthogonal matrix Z1. */
/* >          On exit, if COMPZ='I', the orthogonal matrix Z, and if */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
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

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dgghd3_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, doublereal *work, integer *lwork, integer *info, ftnlen 
	compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublereal s, c1, c2;
    static integer j0;
    static doublereal s1, s2;
    static integer nb, jj, nh, nx, pw, nnb, len, top, ppw, n2nb;
    static logical blk22;
    static integer cola, jcol, ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jrow, topq, ppwo;
    static doublereal temp1, temp2, temp3;
    static integer kacc22;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int dorm22_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer nblst;
    static logical initq, wantq;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical initz, wantz;
    static char compq2[1], compz2[1];
    extern /* Subroutine */ int dgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlaset_(char *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, ftnlen), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
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

/* ===================================================================== */

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

#line 278 "dgghd3.f"
    /* Parameter adjustments */
#line 278 "dgghd3.f"
    a_dim1 = *lda;
#line 278 "dgghd3.f"
    a_offset = 1 + a_dim1;
#line 278 "dgghd3.f"
    a -= a_offset;
#line 278 "dgghd3.f"
    b_dim1 = *ldb;
#line 278 "dgghd3.f"
    b_offset = 1 + b_dim1;
#line 278 "dgghd3.f"
    b -= b_offset;
#line 278 "dgghd3.f"
    q_dim1 = *ldq;
#line 278 "dgghd3.f"
    q_offset = 1 + q_dim1;
#line 278 "dgghd3.f"
    q -= q_offset;
#line 278 "dgghd3.f"
    z_dim1 = *ldz;
#line 278 "dgghd3.f"
    z_offset = 1 + z_dim1;
#line 278 "dgghd3.f"
    z__ -= z_offset;
#line 278 "dgghd3.f"
    --work;
#line 278 "dgghd3.f"

#line 278 "dgghd3.f"
    /* Function Body */
#line 278 "dgghd3.f"
    *info = 0;
#line 279 "dgghd3.f"
    nb = ilaenv_(&c__1, "DGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)
	    1);
/* Computing MAX */
#line 280 "dgghd3.f"
    i__1 = *n * 6 * nb;
#line 280 "dgghd3.f"
    lwkopt = max(i__1,1);
#line 281 "dgghd3.f"
    work[1] = (doublereal) lwkopt;
#line 282 "dgghd3.f"
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 283 "dgghd3.f"
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 284 "dgghd3.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 285 "dgghd3.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);
#line 286 "dgghd3.f"
    lquery = *lwork == -1;

#line 288 "dgghd3.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 289 "dgghd3.f"
	*info = -1;
#line 290 "dgghd3.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 291 "dgghd3.f"
	*info = -2;
#line 292 "dgghd3.f"
    } else if (*n < 0) {
#line 293 "dgghd3.f"
	*info = -3;
#line 294 "dgghd3.f"
    } else if (*ilo < 1) {
#line 295 "dgghd3.f"
	*info = -4;
#line 296 "dgghd3.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 297 "dgghd3.f"
	*info = -5;
#line 298 "dgghd3.f"
    } else if (*lda < max(1,*n)) {
#line 299 "dgghd3.f"
	*info = -7;
#line 300 "dgghd3.f"
    } else if (*ldb < max(1,*n)) {
#line 301 "dgghd3.f"
	*info = -9;
#line 302 "dgghd3.f"
    } else if (wantq && *ldq < *n || *ldq < 1) {
#line 303 "dgghd3.f"
	*info = -11;
#line 304 "dgghd3.f"
    } else if (wantz && *ldz < *n || *ldz < 1) {
#line 305 "dgghd3.f"
	*info = -13;
#line 306 "dgghd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 307 "dgghd3.f"
	*info = -15;
#line 308 "dgghd3.f"
    }
#line 309 "dgghd3.f"
    if (*info != 0) {
#line 310 "dgghd3.f"
	i__1 = -(*info);
#line 310 "dgghd3.f"
	xerbla_("DGGHD3", &i__1, (ftnlen)6);
#line 311 "dgghd3.f"
	return 0;
#line 312 "dgghd3.f"
    } else if (lquery) {
#line 313 "dgghd3.f"
	return 0;
#line 314 "dgghd3.f"
    }

/*     Initialize Q and Z if desired. */

#line 318 "dgghd3.f"
    if (initq) {
#line 318 "dgghd3.f"
	dlaset_("All", n, n, &c_b14, &c_b15, &q[q_offset], ldq, (ftnlen)3);
#line 318 "dgghd3.f"
    }
#line 320 "dgghd3.f"
    if (initz) {
#line 320 "dgghd3.f"
	dlaset_("All", n, n, &c_b14, &c_b15, &z__[z_offset], ldz, (ftnlen)3);
#line 320 "dgghd3.f"
    }

/*     Zero out lower triangle of B. */

#line 325 "dgghd3.f"
    if (*n > 1) {
#line 325 "dgghd3.f"
	i__1 = *n - 1;
#line 325 "dgghd3.f"
	i__2 = *n - 1;
#line 325 "dgghd3.f"
	dlaset_("Lower", &i__1, &i__2, &c_b14, &c_b14, &b[b_dim1 + 2], ldb, (
		ftnlen)5);
#line 325 "dgghd3.f"
    }

/*     Quick return if possible */

#line 330 "dgghd3.f"
    nh = *ihi - *ilo + 1;
#line 331 "dgghd3.f"
    if (nh <= 1) {
#line 332 "dgghd3.f"
	work[1] = 1.;
#line 333 "dgghd3.f"
	return 0;
#line 334 "dgghd3.f"
    }

/*     Determine the blocksize. */

#line 338 "dgghd3.f"
    nbmin = ilaenv_(&c__2, "DGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 339 "dgghd3.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to use unblocked instead of blocked code. */

/* Computing MAX */
#line 343 "dgghd3.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "DGGHD3", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 343 "dgghd3.f"
	nx = max(i__1,i__2);
#line 344 "dgghd3.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code. */

#line 348 "dgghd3.f"
	    if (*lwork < lwkopt) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code. */

/* Computing MAX */
#line 354 "dgghd3.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGGHD3", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 354 "dgghd3.f"
		nbmin = max(i__1,i__2);
#line 356 "dgghd3.f"
		if (*lwork >= *n * 6 * nbmin) {
#line 357 "dgghd3.f"
		    nb = *lwork / (*n * 6);
#line 358 "dgghd3.f"
		} else {
#line 359 "dgghd3.f"
		    nb = 1;
#line 360 "dgghd3.f"
		}
#line 361 "dgghd3.f"
	    }
#line 362 "dgghd3.f"
	}
#line 363 "dgghd3.f"
    }

#line 365 "dgghd3.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 369 "dgghd3.f"
	jcol = *ilo;

#line 371 "dgghd3.f"
    } else {

/*        Use blocked code */

#line 375 "dgghd3.f"
	kacc22 = ilaenv_(&c__16, "DGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 376 "dgghd3.f"
	blk22 = kacc22 == 2;
#line 377 "dgghd3.f"
	i__1 = *ihi - 2;
#line 377 "dgghd3.f"
	i__2 = nb;
#line 377 "dgghd3.f"
	for (jcol = *ilo; i__2 < 0 ? jcol >= i__1 : jcol <= i__1; jcol += 
		i__2) {
/* Computing MIN */
#line 378 "dgghd3.f"
	    i__3 = nb, i__4 = *ihi - jcol - 1;
#line 378 "dgghd3.f"
	    nnb = min(i__3,i__4);

/*           Initialize small orthogonal factors that will hold the */
/*           accumulated Givens rotations in workspace. */
/*           N2NB   denotes the number of 2*NNB-by-2*NNB factors */
/*           NBLST  denotes the (possibly smaller) order of the last */
/*                  factor. */

#line 386 "dgghd3.f"
	    n2nb = (*ihi - jcol - 1) / nnb - 1;
#line 387 "dgghd3.f"
	    nblst = *ihi - jcol - n2nb * nnb;
#line 388 "dgghd3.f"
	    dlaset_("All", &nblst, &nblst, &c_b14, &c_b15, &work[1], &nblst, (
		    ftnlen)3);
#line 389 "dgghd3.f"
	    pw = nblst * nblst + 1;
#line 390 "dgghd3.f"
	    i__3 = n2nb;
#line 390 "dgghd3.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 391 "dgghd3.f"
		i__4 = nnb << 1;
#line 391 "dgghd3.f"
		i__5 = nnb << 1;
#line 391 "dgghd3.f"
		i__6 = nnb << 1;
#line 391 "dgghd3.f"
		dlaset_("All", &i__4, &i__5, &c_b14, &c_b15, &work[pw], &i__6,
			 (ftnlen)3);
#line 393 "dgghd3.f"
		pw += (nnb << 2) * nnb;
#line 394 "dgghd3.f"
	    }

/*           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form. */

#line 398 "dgghd3.f"
	    i__3 = jcol + nnb - 1;
#line 398 "dgghd3.f"
	    for (j = jcol; j <= i__3; ++j) {

/*              Reduce Jth column of A. Store cosines and sines in Jth */
/*              column of A and B, respectively. */

#line 403 "dgghd3.f"
		i__4 = j + 2;
#line 403 "dgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 404 "dgghd3.f"
		    temp = a[i__ - 1 + j * a_dim1];
#line 405 "dgghd3.f"
		    dlartg_(&temp, &a[i__ + j * a_dim1], &c__, &s, &a[i__ - 1 
			    + j * a_dim1]);
#line 406 "dgghd3.f"
		    a[i__ + j * a_dim1] = c__;
#line 407 "dgghd3.f"
		    b[i__ + j * b_dim1] = s;
#line 408 "dgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 412 "dgghd3.f"
		ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 413 "dgghd3.f"
		len = j + 2 - jcol;
#line 414 "dgghd3.f"
		jrow = j + n2nb * nnb + 2;
#line 415 "dgghd3.f"
		i__4 = jrow;
#line 415 "dgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 416 "dgghd3.f"
		    c__ = a[i__ + j * a_dim1];
#line 417 "dgghd3.f"
		    s = b[i__ + j * b_dim1];
#line 418 "dgghd3.f"
		    i__5 = ppw + len - 1;
#line 418 "dgghd3.f"
		    for (jj = ppw; jj <= i__5; ++jj) {
#line 419 "dgghd3.f"
			temp = work[jj + nblst];
#line 420 "dgghd3.f"
			work[jj + nblst] = c__ * temp - s * work[jj];
#line 421 "dgghd3.f"
			work[jj] = s * temp + c__ * work[jj];
#line 422 "dgghd3.f"
		    }
#line 423 "dgghd3.f"
		    ++len;
#line 424 "dgghd3.f"
		    ppw = ppw - nblst - 1;
#line 425 "dgghd3.f"
		}

#line 427 "dgghd3.f"
		ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + nnb;
#line 428 "dgghd3.f"
		j0 = jrow - nnb;
#line 429 "dgghd3.f"
		i__4 = j + 2;
#line 429 "dgghd3.f"
		i__5 = -nnb;
#line 429 "dgghd3.f"
		for (jrow = j0; i__5 < 0 ? jrow >= i__4 : jrow <= i__4; jrow 
			+= i__5) {
#line 430 "dgghd3.f"
		    ppw = ppwo;
#line 431 "dgghd3.f"
		    len = j + 2 - jcol;
#line 432 "dgghd3.f"
		    i__6 = jrow;
#line 432 "dgghd3.f"
		    for (i__ = jrow + nnb - 1; i__ >= i__6; --i__) {
#line 433 "dgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 434 "dgghd3.f"
			s = b[i__ + j * b_dim1];
#line 435 "dgghd3.f"
			i__7 = ppw + len - 1;
#line 435 "dgghd3.f"
			for (jj = ppw; jj <= i__7; ++jj) {
#line 436 "dgghd3.f"
			    temp = work[jj + (nnb << 1)];
#line 437 "dgghd3.f"
			    work[jj + (nnb << 1)] = c__ * temp - s * work[jj];
#line 438 "dgghd3.f"
			    work[jj] = s * temp + c__ * work[jj];
#line 439 "dgghd3.f"
			}
#line 440 "dgghd3.f"
			++len;
#line 441 "dgghd3.f"
			ppw = ppw - (nnb << 1) - 1;
#line 442 "dgghd3.f"
		    }
#line 443 "dgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 444 "dgghd3.f"
		}

/*              TOP denotes the number of top rows in A and B that will */
/*              not be updated during the next steps. */

#line 449 "dgghd3.f"
		if (jcol <= 2) {
#line 450 "dgghd3.f"
		    top = 0;
#line 451 "dgghd3.f"
		} else {
#line 452 "dgghd3.f"
		    top = jcol;
#line 453 "dgghd3.f"
		}

/*              Propagate transformations through B and replace stored */
/*              left sines/cosines by right sines/cosines. */

#line 458 "dgghd3.f"
		i__5 = j + 1;
#line 458 "dgghd3.f"
		for (jj = *n; jj >= i__5; --jj) {

/*                 Update JJth column of B. */

/* Computing MIN */
#line 462 "dgghd3.f"
		    i__4 = jj + 1;
#line 462 "dgghd3.f"
		    i__6 = j + 2;
#line 462 "dgghd3.f"
		    for (i__ = min(i__4,*ihi); i__ >= i__6; --i__) {
#line 463 "dgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 464 "dgghd3.f"
			s = b[i__ + j * b_dim1];
#line 465 "dgghd3.f"
			temp = b[i__ + jj * b_dim1];
#line 466 "dgghd3.f"
			b[i__ + jj * b_dim1] = c__ * temp - s * b[i__ - 1 + 
				jj * b_dim1];
#line 467 "dgghd3.f"
			b[i__ - 1 + jj * b_dim1] = s * temp + c__ * b[i__ - 1 
				+ jj * b_dim1];
#line 468 "dgghd3.f"
		    }

/*                 Annihilate B( JJ+1, JJ ). */

#line 472 "dgghd3.f"
		    if (jj < *ihi) {
#line 473 "dgghd3.f"
			temp = b[jj + 1 + (jj + 1) * b_dim1];
#line 474 "dgghd3.f"
			dlartg_(&temp, &b[jj + 1 + jj * b_dim1], &c__, &s, &b[
				jj + 1 + (jj + 1) * b_dim1]);
#line 476 "dgghd3.f"
			b[jj + 1 + jj * b_dim1] = 0.;
#line 477 "dgghd3.f"
			i__6 = jj - top;
#line 477 "dgghd3.f"
			drot_(&i__6, &b[top + 1 + (jj + 1) * b_dim1], &c__1, &
				b[top + 1 + jj * b_dim1], &c__1, &c__, &s);
#line 479 "dgghd3.f"
			a[jj + 1 + j * a_dim1] = c__;
#line 480 "dgghd3.f"
			b[jj + 1 + j * b_dim1] = -s;
#line 481 "dgghd3.f"
		    }
#line 482 "dgghd3.f"
		}

/*              Update A by transformations from right. */
/*              Explicit loop unrolling provides better performance */
/*              compared to DLASR. */
/*               CALL DLASR( 'Right', 'Variable', 'Backward', IHI-TOP, */
/*     $                     IHI-J, A( J+2, J ), B( J+2, J ), */
/*     $                     A( TOP+1, J+1 ), LDA ) */

#line 491 "dgghd3.f"
		jj = (*ihi - j - 1) % 3;
#line 492 "dgghd3.f"
		i__5 = jj + 1;
#line 492 "dgghd3.f"
		for (i__ = *ihi - j - 3; i__ >= i__5; i__ += -3) {
#line 493 "dgghd3.f"
		    c__ = a[j + 1 + i__ + j * a_dim1];
#line 494 "dgghd3.f"
		    s = -b[j + 1 + i__ + j * b_dim1];
#line 495 "dgghd3.f"
		    c1 = a[j + 2 + i__ + j * a_dim1];
#line 496 "dgghd3.f"
		    s1 = -b[j + 2 + i__ + j * b_dim1];
#line 497 "dgghd3.f"
		    c2 = a[j + 3 + i__ + j * a_dim1];
#line 498 "dgghd3.f"
		    s2 = -b[j + 3 + i__ + j * b_dim1];

#line 500 "dgghd3.f"
		    i__6 = *ihi;
#line 500 "dgghd3.f"
		    for (k = top + 1; k <= i__6; ++k) {
#line 501 "dgghd3.f"
			temp = a[k + (j + i__) * a_dim1];
#line 502 "dgghd3.f"
			temp1 = a[k + (j + i__ + 1) * a_dim1];
#line 503 "dgghd3.f"
			temp2 = a[k + (j + i__ + 2) * a_dim1];
#line 504 "dgghd3.f"
			temp3 = a[k + (j + i__ + 3) * a_dim1];
#line 505 "dgghd3.f"
			a[k + (j + i__ + 3) * a_dim1] = c2 * temp3 + s2 * 
				temp2;
#line 506 "dgghd3.f"
			temp2 = -s2 * temp3 + c2 * temp2;
#line 507 "dgghd3.f"
			a[k + (j + i__ + 2) * a_dim1] = c1 * temp2 + s1 * 
				temp1;
#line 508 "dgghd3.f"
			temp1 = -s1 * temp2 + c1 * temp1;
#line 509 "dgghd3.f"
			a[k + (j + i__ + 1) * a_dim1] = c__ * temp1 + s * 
				temp;
#line 510 "dgghd3.f"
			a[k + (j + i__) * a_dim1] = -s * temp1 + c__ * temp;
#line 511 "dgghd3.f"
		    }
#line 512 "dgghd3.f"
		}

#line 514 "dgghd3.f"
		if (jj > 0) {
#line 515 "dgghd3.f"
		    for (i__ = jj; i__ >= 1; --i__) {
#line 516 "dgghd3.f"
			i__5 = *ihi - top;
#line 516 "dgghd3.f"
			d__1 = -b[j + 1 + i__ + j * b_dim1];
#line 516 "dgghd3.f"
			drot_(&i__5, &a[top + 1 + (j + i__ + 1) * a_dim1], &
				c__1, &a[top + 1 + (j + i__) * a_dim1], &c__1,
				 &a[j + 1 + i__ + j * a_dim1], &d__1);
#line 519 "dgghd3.f"
		    }
#line 520 "dgghd3.f"
		}

/*              Update (J+1)th column of A by transformations from left. */

#line 524 "dgghd3.f"
		if (j < jcol + nnb - 1) {
#line 525 "dgghd3.f"
		    len = j + 1 - jcol;

/*                 Multiply with the trailing accumulated orthogonal */
/*                 matrix, which takes the form */

/*                        [  U11  U12  ] */
/*                    U = [            ], */
/*                        [  U21  U22  ] */

/*                 where U21 is a LEN-by-LEN matrix and U12 is lower */
/*                 triangular. */

#line 537 "dgghd3.f"
		    jrow = *ihi - nblst + 1;
#line 538 "dgghd3.f"
		    dgemv_("Transpose", &nblst, &len, &c_b15, &work[1], &
			    nblst, &a[jrow + (j + 1) * a_dim1], &c__1, &c_b14,
			     &work[pw], &c__1, (ftnlen)9);
#line 541 "dgghd3.f"
		    ppw = pw + len;
#line 542 "dgghd3.f"
		    i__5 = jrow + nblst - len - 1;
#line 542 "dgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 543 "dgghd3.f"
			work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 544 "dgghd3.f"
			++ppw;
#line 545 "dgghd3.f"
		    }
#line 546 "dgghd3.f"
		    i__5 = nblst - len;
#line 546 "dgghd3.f"
		    dtrmv_("Lower", "Transpose", "Non-unit", &i__5, &work[len 
			    * nblst + 1], &nblst, &work[pw + len], &c__1, (
			    ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 549 "dgghd3.f"
		    i__5 = nblst - len;
#line 549 "dgghd3.f"
		    dgemv_("Transpose", &len, &i__5, &c_b15, &work[(len + 1) *
			     nblst - len + 1], &nblst, &a[jrow + nblst - len 
			    + (j + 1) * a_dim1], &c__1, &c_b15, &work[pw + 
			    len], &c__1, (ftnlen)9);
#line 553 "dgghd3.f"
		    ppw = pw;
#line 554 "dgghd3.f"
		    i__5 = jrow + nblst - 1;
#line 554 "dgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 555 "dgghd3.f"
			a[i__ + (j + 1) * a_dim1] = work[ppw];
#line 556 "dgghd3.f"
			++ppw;
#line 557 "dgghd3.f"
		    }

/*                 Multiply with the other accumulated orthogonal */
/*                 matrices, which take the form */

/*                        [  U11  U12   0  ] */
/*                        [                ] */
/*                    U = [  U21  U22   0  ], */
/*                        [                ] */
/*                        [   0    0    I  ] */

/*                 where I denotes the (NNB-LEN)-by-(NNB-LEN) identity */
/*                 matrix, U21 is a LEN-by-LEN upper triangular matrix */
/*                 and U12 is an NNB-by-NNB lower triangular matrix. */

#line 572 "dgghd3.f"
		    ppwo = nblst * nblst + 1;
#line 573 "dgghd3.f"
		    j0 = jrow - nnb;
#line 574 "dgghd3.f"
		    i__5 = jcol + 1;
#line 574 "dgghd3.f"
		    i__6 = -nnb;
#line 574 "dgghd3.f"
		    for (jrow = j0; i__6 < 0 ? jrow >= i__5 : jrow <= i__5; 
			    jrow += i__6) {
#line 575 "dgghd3.f"
			ppw = pw + len;
#line 576 "dgghd3.f"
			i__4 = jrow + nnb - 1;
#line 576 "dgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 577 "dgghd3.f"
			    work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 578 "dgghd3.f"
			    ++ppw;
#line 579 "dgghd3.f"
			}
#line 580 "dgghd3.f"
			ppw = pw;
#line 581 "dgghd3.f"
			i__4 = jrow + nnb + len - 1;
#line 581 "dgghd3.f"
			for (i__ = jrow + nnb; i__ <= i__4; ++i__) {
#line 582 "dgghd3.f"
			    work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 583 "dgghd3.f"
			    ++ppw;
#line 584 "dgghd3.f"
			}
#line 585 "dgghd3.f"
			i__4 = nnb << 1;
#line 585 "dgghd3.f"
			dtrmv_("Upper", "Transpose", "Non-unit", &len, &work[
				ppwo + nnb], &i__4, &work[pw], &c__1, (ftnlen)
				5, (ftnlen)9, (ftnlen)8);
#line 588 "dgghd3.f"
			i__4 = nnb << 1;
#line 588 "dgghd3.f"
			dtrmv_("Lower", "Transpose", "Non-unit", &nnb, &work[
				ppwo + (len << 1) * nnb], &i__4, &work[pw + 
				len], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 591 "dgghd3.f"
			i__4 = nnb << 1;
#line 591 "dgghd3.f"
			dgemv_("Transpose", &nnb, &len, &c_b15, &work[ppwo], &
				i__4, &a[jrow + (j + 1) * a_dim1], &c__1, &
				c_b15, &work[pw], &c__1, (ftnlen)9);
#line 594 "dgghd3.f"
			i__4 = nnb << 1;
#line 594 "dgghd3.f"
			dgemv_("Transpose", &len, &nnb, &c_b15, &work[ppwo + (
				len << 1) * nnb + nnb], &i__4, &a[jrow + nnb 
				+ (j + 1) * a_dim1], &c__1, &c_b15, &work[pw 
				+ len], &c__1, (ftnlen)9);
#line 598 "dgghd3.f"
			ppw = pw;
#line 599 "dgghd3.f"
			i__4 = jrow + len + nnb - 1;
#line 599 "dgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 600 "dgghd3.f"
			    a[i__ + (j + 1) * a_dim1] = work[ppw];
#line 601 "dgghd3.f"
			    ++ppw;
#line 602 "dgghd3.f"
			}
#line 603 "dgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 604 "dgghd3.f"
		    }
#line 605 "dgghd3.f"
		}
#line 606 "dgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to A. */

#line 610 "dgghd3.f"
	    cola = *n - jcol - nnb + 1;
#line 611 "dgghd3.f"
	    j = *ihi - nblst + 1;
#line 612 "dgghd3.f"
	    dgemm_("Transpose", "No Transpose", &nblst, &cola, &nblst, &c_b15,
		     &work[1], &nblst, &a[j + (jcol + nnb) * a_dim1], lda, &
		    c_b14, &work[pw], &nblst, (ftnlen)9, (ftnlen)12);
#line 616 "dgghd3.f"
	    dlacpy_("All", &nblst, &cola, &work[pw], &nblst, &a[j + (jcol + 
		    nnb) * a_dim1], lda, (ftnlen)3);
#line 618 "dgghd3.f"
	    ppwo = nblst * nblst + 1;
#line 619 "dgghd3.f"
	    j0 = j - nnb;
#line 620 "dgghd3.f"
	    i__3 = jcol + 1;
#line 620 "dgghd3.f"
	    i__6 = -nnb;
#line 620 "dgghd3.f"
	    for (j = j0; i__6 < 0 ? j >= i__3 : j <= i__3; j += i__6) {
#line 621 "dgghd3.f"
		if (blk22) {

/*                 Exploit the structure of */

/*                        [  U11  U12  ] */
/*                    U = [            ] */
/*                        [  U21  U22  ], */

/*                 where all blocks are NNB-by-NNB, U21 is upper */
/*                 triangular and U12 is lower triangular. */

#line 632 "dgghd3.f"
		    i__5 = nnb << 1;
#line 632 "dgghd3.f"
		    i__4 = nnb << 1;
#line 632 "dgghd3.f"
		    i__7 = *lwork - pw + 1;
#line 632 "dgghd3.f"
		    dorm22_("Left", "Transpose", &i__5, &cola, &nnb, &nnb, &
			    work[ppwo], &i__4, &a[j + (jcol + nnb) * a_dim1], 
			    lda, &work[pw], &i__7, &ierr, (ftnlen)4, (ftnlen)
			    9);
#line 636 "dgghd3.f"
		} else {

/*                 Ignore the structure of U. */

#line 640 "dgghd3.f"
		    i__5 = nnb << 1;
#line 640 "dgghd3.f"
		    i__4 = nnb << 1;
#line 640 "dgghd3.f"
		    i__7 = nnb << 1;
#line 640 "dgghd3.f"
		    i__8 = nnb << 1;
#line 640 "dgghd3.f"
		    dgemm_("Transpose", "No Transpose", &i__5, &cola, &i__4, &
			    c_b15, &work[ppwo], &i__7, &a[j + (jcol + nnb) * 
			    a_dim1], lda, &c_b14, &work[pw], &i__8, (ftnlen)9,
			     (ftnlen)12);
#line 644 "dgghd3.f"
		    i__5 = nnb << 1;
#line 644 "dgghd3.f"
		    i__4 = nnb << 1;
#line 644 "dgghd3.f"
		    dlacpy_("All", &i__5, &cola, &work[pw], &i__4, &a[j + (
			    jcol + nnb) * a_dim1], lda, (ftnlen)3);
#line 646 "dgghd3.f"
		}
#line 647 "dgghd3.f"
		ppwo += (nnb << 2) * nnb;
#line 648 "dgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to Q. */

#line 652 "dgghd3.f"
	    if (wantq) {
#line 653 "dgghd3.f"
		j = *ihi - nblst + 1;
#line 654 "dgghd3.f"
		if (initq) {
/* Computing MAX */
#line 655 "dgghd3.f"
		    i__6 = 2, i__3 = j - jcol + 1;
#line 655 "dgghd3.f"
		    topq = max(i__6,i__3);
#line 656 "dgghd3.f"
		    nh = *ihi - topq + 1;
#line 657 "dgghd3.f"
		} else {
#line 658 "dgghd3.f"
		    topq = 1;
#line 659 "dgghd3.f"
		    nh = *n;
#line 660 "dgghd3.f"
		}
#line 661 "dgghd3.f"
		dgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b15, &q[topq + j * q_dim1], ldq, &work[1], &nblst, &
			c_b14, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 664 "dgghd3.f"
		dlacpy_("All", &nh, &nblst, &work[pw], &nh, &q[topq + j * 
			q_dim1], ldq, (ftnlen)3);
#line 666 "dgghd3.f"
		ppwo = nblst * nblst + 1;
#line 667 "dgghd3.f"
		j0 = j - nnb;
#line 668 "dgghd3.f"
		i__6 = jcol + 1;
#line 668 "dgghd3.f"
		i__3 = -nnb;
#line 668 "dgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__6 : j <= i__6; j += i__3) {
#line 669 "dgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 670 "dgghd3.f"
			i__5 = 2, i__4 = j - jcol + 1;
#line 670 "dgghd3.f"
			topq = max(i__5,i__4);
#line 671 "dgghd3.f"
			nh = *ihi - topq + 1;
#line 672 "dgghd3.f"
		    }
#line 673 "dgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 677 "dgghd3.f"
			i__5 = nnb << 1;
#line 677 "dgghd3.f"
			i__4 = nnb << 1;
#line 677 "dgghd3.f"
			i__7 = *lwork - pw + 1;
#line 677 "dgghd3.f"
			dorm22_("Right", "No Transpose", &nh, &i__5, &nnb, &
				nnb, &work[ppwo], &i__4, &q[topq + j * q_dim1]
				, ldq, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 681 "dgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 685 "dgghd3.f"
			i__5 = nnb << 1;
#line 685 "dgghd3.f"
			i__4 = nnb << 1;
#line 685 "dgghd3.f"
			i__7 = nnb << 1;
#line 685 "dgghd3.f"
			dgemm_("No Transpose", "No Transpose", &nh, &i__5, &
				i__4, &c_b15, &q[topq + j * q_dim1], ldq, &
				work[ppwo], &i__7, &c_b14, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 689 "dgghd3.f"
			i__5 = nnb << 1;
#line 689 "dgghd3.f"
			dlacpy_("All", &nh, &i__5, &work[pw], &nh, &q[topq + 
				j * q_dim1], ldq, (ftnlen)3);
#line 691 "dgghd3.f"
		    }
#line 692 "dgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 693 "dgghd3.f"
		}
#line 694 "dgghd3.f"
	    }

/*           Accumulate right Givens rotations if required. */

#line 698 "dgghd3.f"
	    if (wantz || top > 0) {

/*              Initialize small orthogonal factors that will hold the */
/*              accumulated Givens rotations in workspace. */

#line 703 "dgghd3.f"
		dlaset_("All", &nblst, &nblst, &c_b14, &c_b15, &work[1], &
			nblst, (ftnlen)3);
#line 705 "dgghd3.f"
		pw = nblst * nblst + 1;
#line 706 "dgghd3.f"
		i__3 = n2nb;
#line 706 "dgghd3.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 707 "dgghd3.f"
		    i__6 = nnb << 1;
#line 707 "dgghd3.f"
		    i__5 = nnb << 1;
#line 707 "dgghd3.f"
		    i__4 = nnb << 1;
#line 707 "dgghd3.f"
		    dlaset_("All", &i__6, &i__5, &c_b14, &c_b15, &work[pw], &
			    i__4, (ftnlen)3);
#line 709 "dgghd3.f"
		    pw += (nnb << 2) * nnb;
#line 710 "dgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 714 "dgghd3.f"
		i__3 = jcol + nnb - 1;
#line 714 "dgghd3.f"
		for (j = jcol; j <= i__3; ++j) {
#line 715 "dgghd3.f"
		    ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 716 "dgghd3.f"
		    len = j + 2 - jcol;
#line 717 "dgghd3.f"
		    jrow = j + n2nb * nnb + 2;
#line 718 "dgghd3.f"
		    i__6 = jrow;
#line 718 "dgghd3.f"
		    for (i__ = *ihi; i__ >= i__6; --i__) {
#line 719 "dgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 720 "dgghd3.f"
			a[i__ + j * a_dim1] = 0.;
#line 721 "dgghd3.f"
			s = b[i__ + j * b_dim1];
#line 722 "dgghd3.f"
			b[i__ + j * b_dim1] = 0.;
#line 723 "dgghd3.f"
			i__5 = ppw + len - 1;
#line 723 "dgghd3.f"
			for (jj = ppw; jj <= i__5; ++jj) {
#line 724 "dgghd3.f"
			    temp = work[jj + nblst];
#line 725 "dgghd3.f"
			    work[jj + nblst] = c__ * temp - s * work[jj];
#line 726 "dgghd3.f"
			    work[jj] = s * temp + c__ * work[jj];
#line 727 "dgghd3.f"
			}
#line 728 "dgghd3.f"
			++len;
#line 729 "dgghd3.f"
			ppw = ppw - nblst - 1;
#line 730 "dgghd3.f"
		    }

#line 732 "dgghd3.f"
		    ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + 
			    nnb;
#line 733 "dgghd3.f"
		    j0 = jrow - nnb;
#line 734 "dgghd3.f"
		    i__6 = j + 2;
#line 734 "dgghd3.f"
		    i__5 = -nnb;
#line 734 "dgghd3.f"
		    for (jrow = j0; i__5 < 0 ? jrow >= i__6 : jrow <= i__6; 
			    jrow += i__5) {
#line 735 "dgghd3.f"
			ppw = ppwo;
#line 736 "dgghd3.f"
			len = j + 2 - jcol;
#line 737 "dgghd3.f"
			i__4 = jrow;
#line 737 "dgghd3.f"
			for (i__ = jrow + nnb - 1; i__ >= i__4; --i__) {
#line 738 "dgghd3.f"
			    c__ = a[i__ + j * a_dim1];
#line 739 "dgghd3.f"
			    a[i__ + j * a_dim1] = 0.;
#line 740 "dgghd3.f"
			    s = b[i__ + j * b_dim1];
#line 741 "dgghd3.f"
			    b[i__ + j * b_dim1] = 0.;
#line 742 "dgghd3.f"
			    i__7 = ppw + len - 1;
#line 742 "dgghd3.f"
			    for (jj = ppw; jj <= i__7; ++jj) {
#line 743 "dgghd3.f"
				temp = work[jj + (nnb << 1)];
#line 744 "dgghd3.f"
				work[jj + (nnb << 1)] = c__ * temp - s * work[
					jj];
#line 745 "dgghd3.f"
				work[jj] = s * temp + c__ * work[jj];
#line 746 "dgghd3.f"
			    }
#line 747 "dgghd3.f"
			    ++len;
#line 748 "dgghd3.f"
			    ppw = ppw - (nnb << 1) - 1;
#line 749 "dgghd3.f"
			}
#line 750 "dgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 751 "dgghd3.f"
		    }
#line 752 "dgghd3.f"
		}
#line 753 "dgghd3.f"
	    } else {

#line 755 "dgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 755 "dgghd3.f"
		dlaset_("Lower", &i__3, &nnb, &c_b14, &c_b14, &a[jcol + 2 + 
			jcol * a_dim1], lda, (ftnlen)5);
#line 757 "dgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 757 "dgghd3.f"
		dlaset_("Lower", &i__3, &nnb, &c_b14, &c_b14, &b[jcol + 2 + 
			jcol * b_dim1], ldb, (ftnlen)5);
#line 759 "dgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to A and B. */

#line 763 "dgghd3.f"
	    if (top > 0) {
#line 764 "dgghd3.f"
		j = *ihi - nblst + 1;
#line 765 "dgghd3.f"
		dgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b15, &a[j * a_dim1 + 1], lda, &work[1], &nblst, &
			c_b14, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 768 "dgghd3.f"
		dlacpy_("All", &top, &nblst, &work[pw], &top, &a[j * a_dim1 + 
			1], lda, (ftnlen)3);
#line 770 "dgghd3.f"
		ppwo = nblst * nblst + 1;
#line 771 "dgghd3.f"
		j0 = j - nnb;
#line 772 "dgghd3.f"
		i__3 = jcol + 1;
#line 772 "dgghd3.f"
		i__5 = -nnb;
#line 772 "dgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 773 "dgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 777 "dgghd3.f"
			i__6 = nnb << 1;
#line 777 "dgghd3.f"
			i__4 = nnb << 1;
#line 777 "dgghd3.f"
			i__7 = *lwork - pw + 1;
#line 777 "dgghd3.f"
			dorm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &a[j * a_dim1 + 1], 
				lda, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 781 "dgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 785 "dgghd3.f"
			i__6 = nnb << 1;
#line 785 "dgghd3.f"
			i__4 = nnb << 1;
#line 785 "dgghd3.f"
			i__7 = nnb << 1;
#line 785 "dgghd3.f"
			dgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b15, &a[j * a_dim1 + 1], lda, &work[
				ppwo], &i__7, &c_b14, &work[pw], &top, (
				ftnlen)12, (ftnlen)12);
#line 789 "dgghd3.f"
			i__6 = nnb << 1;
#line 789 "dgghd3.f"
			dlacpy_("All", &top, &i__6, &work[pw], &top, &a[j * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 791 "dgghd3.f"
		    }
#line 792 "dgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 793 "dgghd3.f"
		}

#line 795 "dgghd3.f"
		j = *ihi - nblst + 1;
#line 796 "dgghd3.f"
		dgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b15, &b[j * b_dim1 + 1], ldb, &work[1], &nblst, &
			c_b14, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 799 "dgghd3.f"
		dlacpy_("All", &top, &nblst, &work[pw], &top, &b[j * b_dim1 + 
			1], ldb, (ftnlen)3);
#line 801 "dgghd3.f"
		ppwo = nblst * nblst + 1;
#line 802 "dgghd3.f"
		j0 = j - nnb;
#line 803 "dgghd3.f"
		i__5 = jcol + 1;
#line 803 "dgghd3.f"
		i__3 = -nnb;
#line 803 "dgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__5 : j <= i__5; j += i__3) {
#line 804 "dgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 808 "dgghd3.f"
			i__6 = nnb << 1;
#line 808 "dgghd3.f"
			i__4 = nnb << 1;
#line 808 "dgghd3.f"
			i__7 = *lwork - pw + 1;
#line 808 "dgghd3.f"
			dorm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &b[j * b_dim1 + 1], 
				ldb, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 812 "dgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 816 "dgghd3.f"
			i__6 = nnb << 1;
#line 816 "dgghd3.f"
			i__4 = nnb << 1;
#line 816 "dgghd3.f"
			i__7 = nnb << 1;
#line 816 "dgghd3.f"
			dgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b15, &b[j * b_dim1 + 1], ldb, &work[
				ppwo], &i__7, &c_b14, &work[pw], &top, (
				ftnlen)12, (ftnlen)12);
#line 820 "dgghd3.f"
			i__6 = nnb << 1;
#line 820 "dgghd3.f"
			dlacpy_("All", &top, &i__6, &work[pw], &top, &b[j * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 822 "dgghd3.f"
		    }
#line 823 "dgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 824 "dgghd3.f"
		}
#line 825 "dgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to Z. */

#line 829 "dgghd3.f"
	    if (wantz) {
#line 830 "dgghd3.f"
		j = *ihi - nblst + 1;
#line 831 "dgghd3.f"
		if (initq) {
/* Computing MAX */
#line 832 "dgghd3.f"
		    i__3 = 2, i__5 = j - jcol + 1;
#line 832 "dgghd3.f"
		    topq = max(i__3,i__5);
#line 833 "dgghd3.f"
		    nh = *ihi - topq + 1;
#line 834 "dgghd3.f"
		} else {
#line 835 "dgghd3.f"
		    topq = 1;
#line 836 "dgghd3.f"
		    nh = *n;
#line 837 "dgghd3.f"
		}
#line 838 "dgghd3.f"
		dgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b15, &z__[topq + j * z_dim1], ldz, &work[1], &nblst,
			 &c_b14, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 841 "dgghd3.f"
		dlacpy_("All", &nh, &nblst, &work[pw], &nh, &z__[topq + j * 
			z_dim1], ldz, (ftnlen)3);
#line 843 "dgghd3.f"
		ppwo = nblst * nblst + 1;
#line 844 "dgghd3.f"
		j0 = j - nnb;
#line 845 "dgghd3.f"
		i__3 = jcol + 1;
#line 845 "dgghd3.f"
		i__5 = -nnb;
#line 845 "dgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 846 "dgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 847 "dgghd3.f"
			i__6 = 2, i__4 = j - jcol + 1;
#line 847 "dgghd3.f"
			topq = max(i__6,i__4);
#line 848 "dgghd3.f"
			nh = *ihi - topq + 1;
#line 849 "dgghd3.f"
		    }
#line 850 "dgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 854 "dgghd3.f"
			i__6 = nnb << 1;
#line 854 "dgghd3.f"
			i__4 = nnb << 1;
#line 854 "dgghd3.f"
			i__7 = *lwork - pw + 1;
#line 854 "dgghd3.f"
			dorm22_("Right", "No Transpose", &nh, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &z__[topq + j * 
				z_dim1], ldz, &work[pw], &i__7, &ierr, (
				ftnlen)5, (ftnlen)12);
#line 858 "dgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 862 "dgghd3.f"
			i__6 = nnb << 1;
#line 862 "dgghd3.f"
			i__4 = nnb << 1;
#line 862 "dgghd3.f"
			i__7 = nnb << 1;
#line 862 "dgghd3.f"
			dgemm_("No Transpose", "No Transpose", &nh, &i__6, &
				i__4, &c_b15, &z__[topq + j * z_dim1], ldz, &
				work[ppwo], &i__7, &c_b14, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 866 "dgghd3.f"
			i__6 = nnb << 1;
#line 866 "dgghd3.f"
			dlacpy_("All", &nh, &i__6, &work[pw], &nh, &z__[topq 
				+ j * z_dim1], ldz, (ftnlen)3);
#line 868 "dgghd3.f"
		    }
#line 869 "dgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 870 "dgghd3.f"
		}
#line 871 "dgghd3.f"
	    }
#line 872 "dgghd3.f"
	}
#line 873 "dgghd3.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */
/*     Avoid re-initialization of modified Q and Z. */

#line 878 "dgghd3.f"
    *(unsigned char *)compq2 = *(unsigned char *)compq;
#line 879 "dgghd3.f"
    *(unsigned char *)compz2 = *(unsigned char *)compz;
#line 880 "dgghd3.f"
    if (jcol != *ilo) {
#line 881 "dgghd3.f"
	if (wantq) {
#line 881 "dgghd3.f"
	    *(unsigned char *)compq2 = 'V';
#line 881 "dgghd3.f"
	}
#line 883 "dgghd3.f"
	if (wantz) {
#line 883 "dgghd3.f"
	    *(unsigned char *)compz2 = 'V';
#line 883 "dgghd3.f"
	}
#line 885 "dgghd3.f"
    }

#line 887 "dgghd3.f"
    if (jcol < *ihi) {
#line 887 "dgghd3.f"
	dgghrd_(compq2, compz2, n, &jcol, ihi, &a[a_offset], lda, &b[b_offset]
		, ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &ierr, (ftnlen)
		1, (ftnlen)1);
#line 887 "dgghd3.f"
    }
#line 890 "dgghd3.f"
    work[1] = (doublereal) lwkopt;

#line 892 "dgghd3.f"
    return 0;

/*     End of DGGHD3 */

} /* dgghd3_ */

