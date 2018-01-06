#line 1 "sgghd3.f"
/* sgghd3.f -- translated by f2c (version 20100827).
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

#line 1 "sgghd3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 0.;
static doublereal c_b15 = 1.;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__16 = 16;

/* > \brief \b SGGHD3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgghd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgghd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgghd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGHD3( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGHD3 reduces a pair of real matrices (A,B) to generalized upper */
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
/* > original equation A*x = lambda*B*x, then SGGHD3 reduces the original */
/* > problem to generalized Hessenberg form. */
/* > */
/* > This is a blocked variant of SGGHRD, using matrix-matrix */
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
/* >          normally set by a previous call to SGGBAL; otherwise they */
/* >          should be set to 1 and N respectively. */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA, N) */
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
/* >          B is REAL array, dimension (LDB, N) */
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
/* >          Q is REAL array, dimension (LDQ, N) */
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
/* >          Z is REAL array, dimension (LDZ, N) */
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
/* >          WORK is REAL array, dimension (LWORK) */
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

/* > \ingroup realOTHERcomputational */

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
/* Subroutine */ int sgghd3_(char *compq, char *compz, integer *n, integer *
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
    static integer jrow, topq, ppwo;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal temp1, temp2, temp3;
    static integer kacc22;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static integer nblst;
    static logical initq;
    extern /* Subroutine */ int sorm22_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical wantq, initz, wantz;
    extern /* Subroutine */ int strmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static char compq2[1], compz2[1];
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), slaset_(char *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, ftnlen), slartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), slacpy_(char *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen);
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

#line 279 "sgghd3.f"
    /* Parameter adjustments */
#line 279 "sgghd3.f"
    a_dim1 = *lda;
#line 279 "sgghd3.f"
    a_offset = 1 + a_dim1;
#line 279 "sgghd3.f"
    a -= a_offset;
#line 279 "sgghd3.f"
    b_dim1 = *ldb;
#line 279 "sgghd3.f"
    b_offset = 1 + b_dim1;
#line 279 "sgghd3.f"
    b -= b_offset;
#line 279 "sgghd3.f"
    q_dim1 = *ldq;
#line 279 "sgghd3.f"
    q_offset = 1 + q_dim1;
#line 279 "sgghd3.f"
    q -= q_offset;
#line 279 "sgghd3.f"
    z_dim1 = *ldz;
#line 279 "sgghd3.f"
    z_offset = 1 + z_dim1;
#line 279 "sgghd3.f"
    z__ -= z_offset;
#line 279 "sgghd3.f"
    --work;
#line 279 "sgghd3.f"

#line 279 "sgghd3.f"
    /* Function Body */
#line 279 "sgghd3.f"
    *info = 0;
#line 280 "sgghd3.f"
    nb = ilaenv_(&c__1, "SGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)
	    1);
/* Computing MAX */
#line 281 "sgghd3.f"
    i__1 = *n * 6 * nb;
#line 281 "sgghd3.f"
    lwkopt = max(i__1,1);
#line 282 "sgghd3.f"
    work[1] = (doublereal) lwkopt;
#line 283 "sgghd3.f"
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 284 "sgghd3.f"
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 285 "sgghd3.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 286 "sgghd3.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);
#line 287 "sgghd3.f"
    lquery = *lwork == -1;

#line 289 "sgghd3.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 290 "sgghd3.f"
	*info = -1;
#line 291 "sgghd3.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 292 "sgghd3.f"
	*info = -2;
#line 293 "sgghd3.f"
    } else if (*n < 0) {
#line 294 "sgghd3.f"
	*info = -3;
#line 295 "sgghd3.f"
    } else if (*ilo < 1) {
#line 296 "sgghd3.f"
	*info = -4;
#line 297 "sgghd3.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 298 "sgghd3.f"
	*info = -5;
#line 299 "sgghd3.f"
    } else if (*lda < max(1,*n)) {
#line 300 "sgghd3.f"
	*info = -7;
#line 301 "sgghd3.f"
    } else if (*ldb < max(1,*n)) {
#line 302 "sgghd3.f"
	*info = -9;
#line 303 "sgghd3.f"
    } else if (wantq && *ldq < *n || *ldq < 1) {
#line 304 "sgghd3.f"
	*info = -11;
#line 305 "sgghd3.f"
    } else if (wantz && *ldz < *n || *ldz < 1) {
#line 306 "sgghd3.f"
	*info = -13;
#line 307 "sgghd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 308 "sgghd3.f"
	*info = -15;
#line 309 "sgghd3.f"
    }
#line 310 "sgghd3.f"
    if (*info != 0) {
#line 311 "sgghd3.f"
	i__1 = -(*info);
#line 311 "sgghd3.f"
	xerbla_("SGGHD3", &i__1, (ftnlen)6);
#line 312 "sgghd3.f"
	return 0;
#line 313 "sgghd3.f"
    } else if (lquery) {
#line 314 "sgghd3.f"
	return 0;
#line 315 "sgghd3.f"
    }

/*     Initialize Q and Z if desired. */

#line 319 "sgghd3.f"
    if (initq) {
#line 319 "sgghd3.f"
	slaset_("All", n, n, &c_b14, &c_b15, &q[q_offset], ldq, (ftnlen)3);
#line 319 "sgghd3.f"
    }
#line 321 "sgghd3.f"
    if (initz) {
#line 321 "sgghd3.f"
	slaset_("All", n, n, &c_b14, &c_b15, &z__[z_offset], ldz, (ftnlen)3);
#line 321 "sgghd3.f"
    }

/*     Zero out lower triangle of B. */

#line 326 "sgghd3.f"
    if (*n > 1) {
#line 326 "sgghd3.f"
	i__1 = *n - 1;
#line 326 "sgghd3.f"
	i__2 = *n - 1;
#line 326 "sgghd3.f"
	slaset_("Lower", &i__1, &i__2, &c_b14, &c_b14, &b[b_dim1 + 2], ldb, (
		ftnlen)5);
#line 326 "sgghd3.f"
    }

/*     Quick return if possible */

#line 331 "sgghd3.f"
    nh = *ihi - *ilo + 1;
#line 332 "sgghd3.f"
    if (nh <= 1) {
#line 333 "sgghd3.f"
	work[1] = 1.;
#line 334 "sgghd3.f"
	return 0;
#line 335 "sgghd3.f"
    }

/*     Determine the blocksize. */

#line 339 "sgghd3.f"
    nbmin = ilaenv_(&c__2, "SGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 340 "sgghd3.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to use unblocked instead of blocked code. */

/* Computing MAX */
#line 344 "sgghd3.f"
	i__1 = nb, i__2 = ilaenv_(&c__3, "SGGHD3", " ", n, ilo, ihi, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 344 "sgghd3.f"
	nx = max(i__1,i__2);
#line 345 "sgghd3.f"
	if (nx < nh) {

/*           Determine if workspace is large enough for blocked code. */

#line 349 "sgghd3.f"
	    if (*lwork < lwkopt) {

/*              Not enough workspace to use optimal NB:  determine the */
/*              minimum value of NB, and reduce NB or force use of */
/*              unblocked code. */

/* Computing MAX */
#line 355 "sgghd3.f"
		i__1 = 2, i__2 = ilaenv_(&c__2, "SGGHD3", " ", n, ilo, ihi, &
			c_n1, (ftnlen)6, (ftnlen)1);
#line 355 "sgghd3.f"
		nbmin = max(i__1,i__2);
#line 357 "sgghd3.f"
		if (*lwork >= *n * 6 * nbmin) {
#line 358 "sgghd3.f"
		    nb = *lwork / (*n * 6);
#line 359 "sgghd3.f"
		} else {
#line 360 "sgghd3.f"
		    nb = 1;
#line 361 "sgghd3.f"
		}
#line 362 "sgghd3.f"
	    }
#line 363 "sgghd3.f"
	}
#line 364 "sgghd3.f"
    }

#line 366 "sgghd3.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code below */

#line 370 "sgghd3.f"
	jcol = *ilo;

#line 372 "sgghd3.f"
    } else {

/*        Use blocked code */

#line 376 "sgghd3.f"
	kacc22 = ilaenv_(&c__16, "SGGHD3", " ", n, ilo, ihi, &c_n1, (ftnlen)6,
		 (ftnlen)1);
#line 377 "sgghd3.f"
	blk22 = kacc22 == 2;
#line 378 "sgghd3.f"
	i__1 = *ihi - 2;
#line 378 "sgghd3.f"
	i__2 = nb;
#line 378 "sgghd3.f"
	for (jcol = *ilo; i__2 < 0 ? jcol >= i__1 : jcol <= i__1; jcol += 
		i__2) {
/* Computing MIN */
#line 379 "sgghd3.f"
	    i__3 = nb, i__4 = *ihi - jcol - 1;
#line 379 "sgghd3.f"
	    nnb = min(i__3,i__4);

/*           Initialize small orthogonal factors that will hold the */
/*           accumulated Givens rotations in workspace. */
/*           N2NB   denotes the number of 2*NNB-by-2*NNB factors */
/*           NBLST  denotes the (possibly smaller) order of the last */
/*                  factor. */

#line 387 "sgghd3.f"
	    n2nb = (*ihi - jcol - 1) / nnb - 1;
#line 388 "sgghd3.f"
	    nblst = *ihi - jcol - n2nb * nnb;
#line 389 "sgghd3.f"
	    slaset_("All", &nblst, &nblst, &c_b14, &c_b15, &work[1], &nblst, (
		    ftnlen)3);
#line 390 "sgghd3.f"
	    pw = nblst * nblst + 1;
#line 391 "sgghd3.f"
	    i__3 = n2nb;
#line 391 "sgghd3.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 392 "sgghd3.f"
		i__4 = nnb << 1;
#line 392 "sgghd3.f"
		i__5 = nnb << 1;
#line 392 "sgghd3.f"
		i__6 = nnb << 1;
#line 392 "sgghd3.f"
		slaset_("All", &i__4, &i__5, &c_b14, &c_b15, &work[pw], &i__6,
			 (ftnlen)3);
#line 394 "sgghd3.f"
		pw += (nnb << 2) * nnb;
#line 395 "sgghd3.f"
	    }

/*           Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form. */

#line 399 "sgghd3.f"
	    i__3 = jcol + nnb - 1;
#line 399 "sgghd3.f"
	    for (j = jcol; j <= i__3; ++j) {

/*              Reduce Jth column of A. Store cosines and sines in Jth */
/*              column of A and B, respectively. */

#line 404 "sgghd3.f"
		i__4 = j + 2;
#line 404 "sgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 405 "sgghd3.f"
		    temp = a[i__ - 1 + j * a_dim1];
#line 406 "sgghd3.f"
		    slartg_(&temp, &a[i__ + j * a_dim1], &c__, &s, &a[i__ - 1 
			    + j * a_dim1]);
#line 407 "sgghd3.f"
		    a[i__ + j * a_dim1] = c__;
#line 408 "sgghd3.f"
		    b[i__ + j * b_dim1] = s;
#line 409 "sgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 413 "sgghd3.f"
		ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 414 "sgghd3.f"
		len = j + 2 - jcol;
#line 415 "sgghd3.f"
		jrow = j + n2nb * nnb + 2;
#line 416 "sgghd3.f"
		i__4 = jrow;
#line 416 "sgghd3.f"
		for (i__ = *ihi; i__ >= i__4; --i__) {
#line 417 "sgghd3.f"
		    c__ = a[i__ + j * a_dim1];
#line 418 "sgghd3.f"
		    s = b[i__ + j * b_dim1];
#line 419 "sgghd3.f"
		    i__5 = ppw + len - 1;
#line 419 "sgghd3.f"
		    for (jj = ppw; jj <= i__5; ++jj) {
#line 420 "sgghd3.f"
			temp = work[jj + nblst];
#line 421 "sgghd3.f"
			work[jj + nblst] = c__ * temp - s * work[jj];
#line 422 "sgghd3.f"
			work[jj] = s * temp + c__ * work[jj];
#line 423 "sgghd3.f"
		    }
#line 424 "sgghd3.f"
		    ++len;
#line 425 "sgghd3.f"
		    ppw = ppw - nblst - 1;
#line 426 "sgghd3.f"
		}

#line 428 "sgghd3.f"
		ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + nnb;
#line 429 "sgghd3.f"
		j0 = jrow - nnb;
#line 430 "sgghd3.f"
		i__4 = j + 2;
#line 430 "sgghd3.f"
		i__5 = -nnb;
#line 430 "sgghd3.f"
		for (jrow = j0; i__5 < 0 ? jrow >= i__4 : jrow <= i__4; jrow 
			+= i__5) {
#line 431 "sgghd3.f"
		    ppw = ppwo;
#line 432 "sgghd3.f"
		    len = j + 2 - jcol;
#line 433 "sgghd3.f"
		    i__6 = jrow;
#line 433 "sgghd3.f"
		    for (i__ = jrow + nnb - 1; i__ >= i__6; --i__) {
#line 434 "sgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 435 "sgghd3.f"
			s = b[i__ + j * b_dim1];
#line 436 "sgghd3.f"
			i__7 = ppw + len - 1;
#line 436 "sgghd3.f"
			for (jj = ppw; jj <= i__7; ++jj) {
#line 437 "sgghd3.f"
			    temp = work[jj + (nnb << 1)];
#line 438 "sgghd3.f"
			    work[jj + (nnb << 1)] = c__ * temp - s * work[jj];
#line 439 "sgghd3.f"
			    work[jj] = s * temp + c__ * work[jj];
#line 440 "sgghd3.f"
			}
#line 441 "sgghd3.f"
			++len;
#line 442 "sgghd3.f"
			ppw = ppw - (nnb << 1) - 1;
#line 443 "sgghd3.f"
		    }
#line 444 "sgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 445 "sgghd3.f"
		}

/*              TOP denotes the number of top rows in A and B that will */
/*              not be updated during the next steps. */

#line 450 "sgghd3.f"
		if (jcol <= 2) {
#line 451 "sgghd3.f"
		    top = 0;
#line 452 "sgghd3.f"
		} else {
#line 453 "sgghd3.f"
		    top = jcol;
#line 454 "sgghd3.f"
		}

/*              Propagate transformations through B and replace stored */
/*              left sines/cosines by right sines/cosines. */

#line 459 "sgghd3.f"
		i__5 = j + 1;
#line 459 "sgghd3.f"
		for (jj = *n; jj >= i__5; --jj) {

/*                 Update JJth column of B. */

/* Computing MIN */
#line 463 "sgghd3.f"
		    i__4 = jj + 1;
#line 463 "sgghd3.f"
		    i__6 = j + 2;
#line 463 "sgghd3.f"
		    for (i__ = min(i__4,*ihi); i__ >= i__6; --i__) {
#line 464 "sgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 465 "sgghd3.f"
			s = b[i__ + j * b_dim1];
#line 466 "sgghd3.f"
			temp = b[i__ + jj * b_dim1];
#line 467 "sgghd3.f"
			b[i__ + jj * b_dim1] = c__ * temp - s * b[i__ - 1 + 
				jj * b_dim1];
#line 468 "sgghd3.f"
			b[i__ - 1 + jj * b_dim1] = s * temp + c__ * b[i__ - 1 
				+ jj * b_dim1];
#line 469 "sgghd3.f"
		    }

/*                 Annihilate B( JJ+1, JJ ). */

#line 473 "sgghd3.f"
		    if (jj < *ihi) {
#line 474 "sgghd3.f"
			temp = b[jj + 1 + (jj + 1) * b_dim1];
#line 475 "sgghd3.f"
			slartg_(&temp, &b[jj + 1 + jj * b_dim1], &c__, &s, &b[
				jj + 1 + (jj + 1) * b_dim1]);
#line 477 "sgghd3.f"
			b[jj + 1 + jj * b_dim1] = 0.;
#line 478 "sgghd3.f"
			i__6 = jj - top;
#line 478 "sgghd3.f"
			srot_(&i__6, &b[top + 1 + (jj + 1) * b_dim1], &c__1, &
				b[top + 1 + jj * b_dim1], &c__1, &c__, &s);
#line 480 "sgghd3.f"
			a[jj + 1 + j * a_dim1] = c__;
#line 481 "sgghd3.f"
			b[jj + 1 + j * b_dim1] = -s;
#line 482 "sgghd3.f"
		    }
#line 483 "sgghd3.f"
		}

/*              Update A by transformations from right. */
/*              Explicit loop unrolling provides better performance */
/*              compared to SLASR. */
/*               CALL SLASR( 'Right', 'Variable', 'Backward', IHI-TOP, */
/*     $                     IHI-J, A( J+2, J ), B( J+2, J ), */
/*     $                     A( TOP+1, J+1 ), LDA ) */

#line 492 "sgghd3.f"
		jj = (*ihi - j - 1) % 3;
#line 493 "sgghd3.f"
		i__5 = jj + 1;
#line 493 "sgghd3.f"
		for (i__ = *ihi - j - 3; i__ >= i__5; i__ += -3) {
#line 494 "sgghd3.f"
		    c__ = a[j + 1 + i__ + j * a_dim1];
#line 495 "sgghd3.f"
		    s = -b[j + 1 + i__ + j * b_dim1];
#line 496 "sgghd3.f"
		    c1 = a[j + 2 + i__ + j * a_dim1];
#line 497 "sgghd3.f"
		    s1 = -b[j + 2 + i__ + j * b_dim1];
#line 498 "sgghd3.f"
		    c2 = a[j + 3 + i__ + j * a_dim1];
#line 499 "sgghd3.f"
		    s2 = -b[j + 3 + i__ + j * b_dim1];

#line 501 "sgghd3.f"
		    i__6 = *ihi;
#line 501 "sgghd3.f"
		    for (k = top + 1; k <= i__6; ++k) {
#line 502 "sgghd3.f"
			temp = a[k + (j + i__) * a_dim1];
#line 503 "sgghd3.f"
			temp1 = a[k + (j + i__ + 1) * a_dim1];
#line 504 "sgghd3.f"
			temp2 = a[k + (j + i__ + 2) * a_dim1];
#line 505 "sgghd3.f"
			temp3 = a[k + (j + i__ + 3) * a_dim1];
#line 506 "sgghd3.f"
			a[k + (j + i__ + 3) * a_dim1] = c2 * temp3 + s2 * 
				temp2;
#line 507 "sgghd3.f"
			temp2 = -s2 * temp3 + c2 * temp2;
#line 508 "sgghd3.f"
			a[k + (j + i__ + 2) * a_dim1] = c1 * temp2 + s1 * 
				temp1;
#line 509 "sgghd3.f"
			temp1 = -s1 * temp2 + c1 * temp1;
#line 510 "sgghd3.f"
			a[k + (j + i__ + 1) * a_dim1] = c__ * temp1 + s * 
				temp;
#line 511 "sgghd3.f"
			a[k + (j + i__) * a_dim1] = -s * temp1 + c__ * temp;
#line 512 "sgghd3.f"
		    }
#line 513 "sgghd3.f"
		}

#line 515 "sgghd3.f"
		if (jj > 0) {
#line 516 "sgghd3.f"
		    for (i__ = jj; i__ >= 1; --i__) {
#line 517 "sgghd3.f"
			i__5 = *ihi - top;
#line 517 "sgghd3.f"
			d__1 = -b[j + 1 + i__ + j * b_dim1];
#line 517 "sgghd3.f"
			srot_(&i__5, &a[top + 1 + (j + i__ + 1) * a_dim1], &
				c__1, &a[top + 1 + (j + i__) * a_dim1], &c__1,
				 &a[j + 1 + i__ + j * a_dim1], &d__1);
#line 520 "sgghd3.f"
		    }
#line 521 "sgghd3.f"
		}

/*              Update (J+1)th column of A by transformations from left. */

#line 525 "sgghd3.f"
		if (j < jcol + nnb - 1) {
#line 526 "sgghd3.f"
		    len = j + 1 - jcol;

/*                 Multiply with the trailing accumulated orthogonal */
/*                 matrix, which takes the form */

/*                        [  U11  U12  ] */
/*                    U = [            ], */
/*                        [  U21  U22  ] */

/*                 where U21 is a LEN-by-LEN matrix and U12 is lower */
/*                 triangular. */

#line 538 "sgghd3.f"
		    jrow = *ihi - nblst + 1;
#line 539 "sgghd3.f"
		    sgemv_("Transpose", &nblst, &len, &c_b15, &work[1], &
			    nblst, &a[jrow + (j + 1) * a_dim1], &c__1, &c_b14,
			     &work[pw], &c__1, (ftnlen)9);
#line 542 "sgghd3.f"
		    ppw = pw + len;
#line 543 "sgghd3.f"
		    i__5 = jrow + nblst - len - 1;
#line 543 "sgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 544 "sgghd3.f"
			work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 545 "sgghd3.f"
			++ppw;
#line 546 "sgghd3.f"
		    }
#line 547 "sgghd3.f"
		    i__5 = nblst - len;
#line 547 "sgghd3.f"
		    strmv_("Lower", "Transpose", "Non-unit", &i__5, &work[len 
			    * nblst + 1], &nblst, &work[pw + len], &c__1, (
			    ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 550 "sgghd3.f"
		    i__5 = nblst - len;
#line 550 "sgghd3.f"
		    sgemv_("Transpose", &len, &i__5, &c_b15, &work[(len + 1) *
			     nblst - len + 1], &nblst, &a[jrow + nblst - len 
			    + (j + 1) * a_dim1], &c__1, &c_b15, &work[pw + 
			    len], &c__1, (ftnlen)9);
#line 554 "sgghd3.f"
		    ppw = pw;
#line 555 "sgghd3.f"
		    i__5 = jrow + nblst - 1;
#line 555 "sgghd3.f"
		    for (i__ = jrow; i__ <= i__5; ++i__) {
#line 556 "sgghd3.f"
			a[i__ + (j + 1) * a_dim1] = work[ppw];
#line 557 "sgghd3.f"
			++ppw;
#line 558 "sgghd3.f"
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

#line 573 "sgghd3.f"
		    ppwo = nblst * nblst + 1;
#line 574 "sgghd3.f"
		    j0 = jrow - nnb;
#line 575 "sgghd3.f"
		    i__5 = jcol + 1;
#line 575 "sgghd3.f"
		    i__6 = -nnb;
#line 575 "sgghd3.f"
		    for (jrow = j0; i__6 < 0 ? jrow >= i__5 : jrow <= i__5; 
			    jrow += i__6) {
#line 576 "sgghd3.f"
			ppw = pw + len;
#line 577 "sgghd3.f"
			i__4 = jrow + nnb - 1;
#line 577 "sgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 578 "sgghd3.f"
			    work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 579 "sgghd3.f"
			    ++ppw;
#line 580 "sgghd3.f"
			}
#line 581 "sgghd3.f"
			ppw = pw;
#line 582 "sgghd3.f"
			i__4 = jrow + nnb + len - 1;
#line 582 "sgghd3.f"
			for (i__ = jrow + nnb; i__ <= i__4; ++i__) {
#line 583 "sgghd3.f"
			    work[ppw] = a[i__ + (j + 1) * a_dim1];
#line 584 "sgghd3.f"
			    ++ppw;
#line 585 "sgghd3.f"
			}
#line 586 "sgghd3.f"
			i__4 = nnb << 1;
#line 586 "sgghd3.f"
			strmv_("Upper", "Transpose", "Non-unit", &len, &work[
				ppwo + nnb], &i__4, &work[pw], &c__1, (ftnlen)
				5, (ftnlen)9, (ftnlen)8);
#line 589 "sgghd3.f"
			i__4 = nnb << 1;
#line 589 "sgghd3.f"
			strmv_("Lower", "Transpose", "Non-unit", &nnb, &work[
				ppwo + (len << 1) * nnb], &i__4, &work[pw + 
				len], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 592 "sgghd3.f"
			i__4 = nnb << 1;
#line 592 "sgghd3.f"
			sgemv_("Transpose", &nnb, &len, &c_b15, &work[ppwo], &
				i__4, &a[jrow + (j + 1) * a_dim1], &c__1, &
				c_b15, &work[pw], &c__1, (ftnlen)9);
#line 595 "sgghd3.f"
			i__4 = nnb << 1;
#line 595 "sgghd3.f"
			sgemv_("Transpose", &len, &nnb, &c_b15, &work[ppwo + (
				len << 1) * nnb + nnb], &i__4, &a[jrow + nnb 
				+ (j + 1) * a_dim1], &c__1, &c_b15, &work[pw 
				+ len], &c__1, (ftnlen)9);
#line 599 "sgghd3.f"
			ppw = pw;
#line 600 "sgghd3.f"
			i__4 = jrow + len + nnb - 1;
#line 600 "sgghd3.f"
			for (i__ = jrow; i__ <= i__4; ++i__) {
#line 601 "sgghd3.f"
			    a[i__ + (j + 1) * a_dim1] = work[ppw];
#line 602 "sgghd3.f"
			    ++ppw;
#line 603 "sgghd3.f"
			}
#line 604 "sgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 605 "sgghd3.f"
		    }
#line 606 "sgghd3.f"
		}
#line 607 "sgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to A. */

#line 611 "sgghd3.f"
	    cola = *n - jcol - nnb + 1;
#line 612 "sgghd3.f"
	    j = *ihi - nblst + 1;
#line 613 "sgghd3.f"
	    sgemm_("Transpose", "No Transpose", &nblst, &cola, &nblst, &c_b15,
		     &work[1], &nblst, &a[j + (jcol + nnb) * a_dim1], lda, &
		    c_b14, &work[pw], &nblst, (ftnlen)9, (ftnlen)12);
#line 617 "sgghd3.f"
	    slacpy_("All", &nblst, &cola, &work[pw], &nblst, &a[j + (jcol + 
		    nnb) * a_dim1], lda, (ftnlen)3);
#line 619 "sgghd3.f"
	    ppwo = nblst * nblst + 1;
#line 620 "sgghd3.f"
	    j0 = j - nnb;
#line 621 "sgghd3.f"
	    i__3 = jcol + 1;
#line 621 "sgghd3.f"
	    i__6 = -nnb;
#line 621 "sgghd3.f"
	    for (j = j0; i__6 < 0 ? j >= i__3 : j <= i__3; j += i__6) {
#line 622 "sgghd3.f"
		if (blk22) {

/*                 Exploit the structure of */

/*                        [  U11  U12  ] */
/*                    U = [            ] */
/*                        [  U21  U22  ], */

/*                 where all blocks are NNB-by-NNB, U21 is upper */
/*                 triangular and U12 is lower triangular. */

#line 633 "sgghd3.f"
		    i__5 = nnb << 1;
#line 633 "sgghd3.f"
		    i__4 = nnb << 1;
#line 633 "sgghd3.f"
		    i__7 = *lwork - pw + 1;
#line 633 "sgghd3.f"
		    sorm22_("Left", "Transpose", &i__5, &cola, &nnb, &nnb, &
			    work[ppwo], &i__4, &a[j + (jcol + nnb) * a_dim1], 
			    lda, &work[pw], &i__7, &ierr, (ftnlen)4, (ftnlen)
			    9);
#line 637 "sgghd3.f"
		} else {

/*                 Ignore the structure of U. */

#line 641 "sgghd3.f"
		    i__5 = nnb << 1;
#line 641 "sgghd3.f"
		    i__4 = nnb << 1;
#line 641 "sgghd3.f"
		    i__7 = nnb << 1;
#line 641 "sgghd3.f"
		    i__8 = nnb << 1;
#line 641 "sgghd3.f"
		    sgemm_("Transpose", "No Transpose", &i__5, &cola, &i__4, &
			    c_b15, &work[ppwo], &i__7, &a[j + (jcol + nnb) * 
			    a_dim1], lda, &c_b14, &work[pw], &i__8, (ftnlen)9,
			     (ftnlen)12);
#line 645 "sgghd3.f"
		    i__5 = nnb << 1;
#line 645 "sgghd3.f"
		    i__4 = nnb << 1;
#line 645 "sgghd3.f"
		    slacpy_("All", &i__5, &cola, &work[pw], &i__4, &a[j + (
			    jcol + nnb) * a_dim1], lda, (ftnlen)3);
#line 647 "sgghd3.f"
		}
#line 648 "sgghd3.f"
		ppwo += (nnb << 2) * nnb;
#line 649 "sgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to Q. */

#line 653 "sgghd3.f"
	    if (wantq) {
#line 654 "sgghd3.f"
		j = *ihi - nblst + 1;
#line 655 "sgghd3.f"
		if (initq) {
/* Computing MAX */
#line 656 "sgghd3.f"
		    i__6 = 2, i__3 = j - jcol + 1;
#line 656 "sgghd3.f"
		    topq = max(i__6,i__3);
#line 657 "sgghd3.f"
		    nh = *ihi - topq + 1;
#line 658 "sgghd3.f"
		} else {
#line 659 "sgghd3.f"
		    topq = 1;
#line 660 "sgghd3.f"
		    nh = *n;
#line 661 "sgghd3.f"
		}
#line 662 "sgghd3.f"
		sgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b15, &q[topq + j * q_dim1], ldq, &work[1], &nblst, &
			c_b14, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 665 "sgghd3.f"
		slacpy_("All", &nh, &nblst, &work[pw], &nh, &q[topq + j * 
			q_dim1], ldq, (ftnlen)3);
#line 667 "sgghd3.f"
		ppwo = nblst * nblst + 1;
#line 668 "sgghd3.f"
		j0 = j - nnb;
#line 669 "sgghd3.f"
		i__6 = jcol + 1;
#line 669 "sgghd3.f"
		i__3 = -nnb;
#line 669 "sgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__6 : j <= i__6; j += i__3) {
#line 670 "sgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 671 "sgghd3.f"
			i__5 = 2, i__4 = j - jcol + 1;
#line 671 "sgghd3.f"
			topq = max(i__5,i__4);
#line 672 "sgghd3.f"
			nh = *ihi - topq + 1;
#line 673 "sgghd3.f"
		    }
#line 674 "sgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 678 "sgghd3.f"
			i__5 = nnb << 1;
#line 678 "sgghd3.f"
			i__4 = nnb << 1;
#line 678 "sgghd3.f"
			i__7 = *lwork - pw + 1;
#line 678 "sgghd3.f"
			sorm22_("Right", "No Transpose", &nh, &i__5, &nnb, &
				nnb, &work[ppwo], &i__4, &q[topq + j * q_dim1]
				, ldq, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 682 "sgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 686 "sgghd3.f"
			i__5 = nnb << 1;
#line 686 "sgghd3.f"
			i__4 = nnb << 1;
#line 686 "sgghd3.f"
			i__7 = nnb << 1;
#line 686 "sgghd3.f"
			sgemm_("No Transpose", "No Transpose", &nh, &i__5, &
				i__4, &c_b15, &q[topq + j * q_dim1], ldq, &
				work[ppwo], &i__7, &c_b14, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 690 "sgghd3.f"
			i__5 = nnb << 1;
#line 690 "sgghd3.f"
			slacpy_("All", &nh, &i__5, &work[pw], &nh, &q[topq + 
				j * q_dim1], ldq, (ftnlen)3);
#line 692 "sgghd3.f"
		    }
#line 693 "sgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 694 "sgghd3.f"
		}
#line 695 "sgghd3.f"
	    }

/*           Accumulate right Givens rotations if required. */

#line 699 "sgghd3.f"
	    if (wantz || top > 0) {

/*              Initialize small orthogonal factors that will hold the */
/*              accumulated Givens rotations in workspace. */

#line 704 "sgghd3.f"
		slaset_("All", &nblst, &nblst, &c_b14, &c_b15, &work[1], &
			nblst, (ftnlen)3);
#line 706 "sgghd3.f"
		pw = nblst * nblst + 1;
#line 707 "sgghd3.f"
		i__3 = n2nb;
#line 707 "sgghd3.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 708 "sgghd3.f"
		    i__6 = nnb << 1;
#line 708 "sgghd3.f"
		    i__5 = nnb << 1;
#line 708 "sgghd3.f"
		    i__4 = nnb << 1;
#line 708 "sgghd3.f"
		    slaset_("All", &i__6, &i__5, &c_b14, &c_b15, &work[pw], &
			    i__4, (ftnlen)3);
#line 710 "sgghd3.f"
		    pw += (nnb << 2) * nnb;
#line 711 "sgghd3.f"
		}

/*              Accumulate Givens rotations into workspace array. */

#line 715 "sgghd3.f"
		i__3 = jcol + nnb - 1;
#line 715 "sgghd3.f"
		for (j = jcol; j <= i__3; ++j) {
#line 716 "sgghd3.f"
		    ppw = (nblst + 1) * (nblst - 2) - j + jcol + 1;
#line 717 "sgghd3.f"
		    len = j + 2 - jcol;
#line 718 "sgghd3.f"
		    jrow = j + n2nb * nnb + 2;
#line 719 "sgghd3.f"
		    i__6 = jrow;
#line 719 "sgghd3.f"
		    for (i__ = *ihi; i__ >= i__6; --i__) {
#line 720 "sgghd3.f"
			c__ = a[i__ + j * a_dim1];
#line 721 "sgghd3.f"
			a[i__ + j * a_dim1] = 0.;
#line 722 "sgghd3.f"
			s = b[i__ + j * b_dim1];
#line 723 "sgghd3.f"
			b[i__ + j * b_dim1] = 0.;
#line 724 "sgghd3.f"
			i__5 = ppw + len - 1;
#line 724 "sgghd3.f"
			for (jj = ppw; jj <= i__5; ++jj) {
#line 725 "sgghd3.f"
			    temp = work[jj + nblst];
#line 726 "sgghd3.f"
			    work[jj + nblst] = c__ * temp - s * work[jj];
#line 727 "sgghd3.f"
			    work[jj] = s * temp + c__ * work[jj];
#line 728 "sgghd3.f"
			}
#line 729 "sgghd3.f"
			++len;
#line 730 "sgghd3.f"
			ppw = ppw - nblst - 1;
#line 731 "sgghd3.f"
		    }

#line 733 "sgghd3.f"
		    ppwo = nblst * nblst + (nnb + j - jcol - 1 << 1) * nnb + 
			    nnb;
#line 734 "sgghd3.f"
		    j0 = jrow - nnb;
#line 735 "sgghd3.f"
		    i__6 = j + 2;
#line 735 "sgghd3.f"
		    i__5 = -nnb;
#line 735 "sgghd3.f"
		    for (jrow = j0; i__5 < 0 ? jrow >= i__6 : jrow <= i__6; 
			    jrow += i__5) {
#line 736 "sgghd3.f"
			ppw = ppwo;
#line 737 "sgghd3.f"
			len = j + 2 - jcol;
#line 738 "sgghd3.f"
			i__4 = jrow;
#line 738 "sgghd3.f"
			for (i__ = jrow + nnb - 1; i__ >= i__4; --i__) {
#line 739 "sgghd3.f"
			    c__ = a[i__ + j * a_dim1];
#line 740 "sgghd3.f"
			    a[i__ + j * a_dim1] = 0.;
#line 741 "sgghd3.f"
			    s = b[i__ + j * b_dim1];
#line 742 "sgghd3.f"
			    b[i__ + j * b_dim1] = 0.;
#line 743 "sgghd3.f"
			    i__7 = ppw + len - 1;
#line 743 "sgghd3.f"
			    for (jj = ppw; jj <= i__7; ++jj) {
#line 744 "sgghd3.f"
				temp = work[jj + (nnb << 1)];
#line 745 "sgghd3.f"
				work[jj + (nnb << 1)] = c__ * temp - s * work[
					jj];
#line 746 "sgghd3.f"
				work[jj] = s * temp + c__ * work[jj];
#line 747 "sgghd3.f"
			    }
#line 748 "sgghd3.f"
			    ++len;
#line 749 "sgghd3.f"
			    ppw = ppw - (nnb << 1) - 1;
#line 750 "sgghd3.f"
			}
#line 751 "sgghd3.f"
			ppwo += (nnb << 2) * nnb;
#line 752 "sgghd3.f"
		    }
#line 753 "sgghd3.f"
		}
#line 754 "sgghd3.f"
	    } else {

#line 756 "sgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 756 "sgghd3.f"
		slaset_("Lower", &i__3, &nnb, &c_b14, &c_b14, &a[jcol + 2 + 
			jcol * a_dim1], lda, (ftnlen)5);
#line 758 "sgghd3.f"
		i__3 = *ihi - jcol - 1;
#line 758 "sgghd3.f"
		slaset_("Lower", &i__3, &nnb, &c_b14, &c_b14, &b[jcol + 2 + 
			jcol * b_dim1], ldb, (ftnlen)5);
#line 760 "sgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to A and B. */

#line 764 "sgghd3.f"
	    if (top > 0) {
#line 765 "sgghd3.f"
		j = *ihi - nblst + 1;
#line 766 "sgghd3.f"
		sgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b15, &a[j * a_dim1 + 1], lda, &work[1], &nblst, &
			c_b14, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 769 "sgghd3.f"
		slacpy_("All", &top, &nblst, &work[pw], &top, &a[j * a_dim1 + 
			1], lda, (ftnlen)3);
#line 771 "sgghd3.f"
		ppwo = nblst * nblst + 1;
#line 772 "sgghd3.f"
		j0 = j - nnb;
#line 773 "sgghd3.f"
		i__3 = jcol + 1;
#line 773 "sgghd3.f"
		i__5 = -nnb;
#line 773 "sgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 774 "sgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 778 "sgghd3.f"
			i__6 = nnb << 1;
#line 778 "sgghd3.f"
			i__4 = nnb << 1;
#line 778 "sgghd3.f"
			i__7 = *lwork - pw + 1;
#line 778 "sgghd3.f"
			sorm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &a[j * a_dim1 + 1], 
				lda, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 782 "sgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 786 "sgghd3.f"
			i__6 = nnb << 1;
#line 786 "sgghd3.f"
			i__4 = nnb << 1;
#line 786 "sgghd3.f"
			i__7 = nnb << 1;
#line 786 "sgghd3.f"
			sgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b15, &a[j * a_dim1 + 1], lda, &work[
				ppwo], &i__7, &c_b14, &work[pw], &top, (
				ftnlen)12, (ftnlen)12);
#line 790 "sgghd3.f"
			i__6 = nnb << 1;
#line 790 "sgghd3.f"
			slacpy_("All", &top, &i__6, &work[pw], &top, &a[j * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 792 "sgghd3.f"
		    }
#line 793 "sgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 794 "sgghd3.f"
		}

#line 796 "sgghd3.f"
		j = *ihi - nblst + 1;
#line 797 "sgghd3.f"
		sgemm_("No Transpose", "No Transpose", &top, &nblst, &nblst, &
			c_b15, &b[j * b_dim1 + 1], ldb, &work[1], &nblst, &
			c_b14, &work[pw], &top, (ftnlen)12, (ftnlen)12);
#line 800 "sgghd3.f"
		slacpy_("All", &top, &nblst, &work[pw], &top, &b[j * b_dim1 + 
			1], ldb, (ftnlen)3);
#line 802 "sgghd3.f"
		ppwo = nblst * nblst + 1;
#line 803 "sgghd3.f"
		j0 = j - nnb;
#line 804 "sgghd3.f"
		i__5 = jcol + 1;
#line 804 "sgghd3.f"
		i__3 = -nnb;
#line 804 "sgghd3.f"
		for (j = j0; i__3 < 0 ? j >= i__5 : j <= i__5; j += i__3) {
#line 805 "sgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 809 "sgghd3.f"
			i__6 = nnb << 1;
#line 809 "sgghd3.f"
			i__4 = nnb << 1;
#line 809 "sgghd3.f"
			i__7 = *lwork - pw + 1;
#line 809 "sgghd3.f"
			sorm22_("Right", "No Transpose", &top, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &b[j * b_dim1 + 1], 
				ldb, &work[pw], &i__7, &ierr, (ftnlen)5, (
				ftnlen)12);
#line 813 "sgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 817 "sgghd3.f"
			i__6 = nnb << 1;
#line 817 "sgghd3.f"
			i__4 = nnb << 1;
#line 817 "sgghd3.f"
			i__7 = nnb << 1;
#line 817 "sgghd3.f"
			sgemm_("No Transpose", "No Transpose", &top, &i__6, &
				i__4, &c_b15, &b[j * b_dim1 + 1], ldb, &work[
				ppwo], &i__7, &c_b14, &work[pw], &top, (
				ftnlen)12, (ftnlen)12);
#line 821 "sgghd3.f"
			i__6 = nnb << 1;
#line 821 "sgghd3.f"
			slacpy_("All", &top, &i__6, &work[pw], &top, &b[j * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 823 "sgghd3.f"
		    }
#line 824 "sgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 825 "sgghd3.f"
		}
#line 826 "sgghd3.f"
	    }

/*           Apply accumulated orthogonal matrices to Z. */

#line 830 "sgghd3.f"
	    if (wantz) {
#line 831 "sgghd3.f"
		j = *ihi - nblst + 1;
#line 832 "sgghd3.f"
		if (initq) {
/* Computing MAX */
#line 833 "sgghd3.f"
		    i__3 = 2, i__5 = j - jcol + 1;
#line 833 "sgghd3.f"
		    topq = max(i__3,i__5);
#line 834 "sgghd3.f"
		    nh = *ihi - topq + 1;
#line 835 "sgghd3.f"
		} else {
#line 836 "sgghd3.f"
		    topq = 1;
#line 837 "sgghd3.f"
		    nh = *n;
#line 838 "sgghd3.f"
		}
#line 839 "sgghd3.f"
		sgemm_("No Transpose", "No Transpose", &nh, &nblst, &nblst, &
			c_b15, &z__[topq + j * z_dim1], ldz, &work[1], &nblst,
			 &c_b14, &work[pw], &nh, (ftnlen)12, (ftnlen)12);
#line 842 "sgghd3.f"
		slacpy_("All", &nh, &nblst, &work[pw], &nh, &z__[topq + j * 
			z_dim1], ldz, (ftnlen)3);
#line 844 "sgghd3.f"
		ppwo = nblst * nblst + 1;
#line 845 "sgghd3.f"
		j0 = j - nnb;
#line 846 "sgghd3.f"
		i__3 = jcol + 1;
#line 846 "sgghd3.f"
		i__5 = -nnb;
#line 846 "sgghd3.f"
		for (j = j0; i__5 < 0 ? j >= i__3 : j <= i__3; j += i__5) {
#line 847 "sgghd3.f"
		    if (initq) {
/* Computing MAX */
#line 848 "sgghd3.f"
			i__6 = 2, i__4 = j - jcol + 1;
#line 848 "sgghd3.f"
			topq = max(i__6,i__4);
#line 849 "sgghd3.f"
			nh = *ihi - topq + 1;
#line 850 "sgghd3.f"
		    }
#line 851 "sgghd3.f"
		    if (blk22) {

/*                    Exploit the structure of U. */

#line 855 "sgghd3.f"
			i__6 = nnb << 1;
#line 855 "sgghd3.f"
			i__4 = nnb << 1;
#line 855 "sgghd3.f"
			i__7 = *lwork - pw + 1;
#line 855 "sgghd3.f"
			sorm22_("Right", "No Transpose", &nh, &i__6, &nnb, &
				nnb, &work[ppwo], &i__4, &z__[topq + j * 
				z_dim1], ldz, &work[pw], &i__7, &ierr, (
				ftnlen)5, (ftnlen)12);
#line 859 "sgghd3.f"
		    } else {

/*                    Ignore the structure of U. */

#line 863 "sgghd3.f"
			i__6 = nnb << 1;
#line 863 "sgghd3.f"
			i__4 = nnb << 1;
#line 863 "sgghd3.f"
			i__7 = nnb << 1;
#line 863 "sgghd3.f"
			sgemm_("No Transpose", "No Transpose", &nh, &i__6, &
				i__4, &c_b15, &z__[topq + j * z_dim1], ldz, &
				work[ppwo], &i__7, &c_b14, &work[pw], &nh, (
				ftnlen)12, (ftnlen)12);
#line 867 "sgghd3.f"
			i__6 = nnb << 1;
#line 867 "sgghd3.f"
			slacpy_("All", &nh, &i__6, &work[pw], &nh, &z__[topq 
				+ j * z_dim1], ldz, (ftnlen)3);
#line 869 "sgghd3.f"
		    }
#line 870 "sgghd3.f"
		    ppwo += (nnb << 2) * nnb;
#line 871 "sgghd3.f"
		}
#line 872 "sgghd3.f"
	    }
#line 873 "sgghd3.f"
	}
#line 874 "sgghd3.f"
    }

/*     Use unblocked code to reduce the rest of the matrix */
/*     Avoid re-initialization of modified Q and Z. */

#line 879 "sgghd3.f"
    *(unsigned char *)compq2 = *(unsigned char *)compq;
#line 880 "sgghd3.f"
    *(unsigned char *)compz2 = *(unsigned char *)compz;
#line 881 "sgghd3.f"
    if (jcol != *ilo) {
#line 882 "sgghd3.f"
	if (wantq) {
#line 882 "sgghd3.f"
	    *(unsigned char *)compq2 = 'V';
#line 882 "sgghd3.f"
	}
#line 884 "sgghd3.f"
	if (wantz) {
#line 884 "sgghd3.f"
	    *(unsigned char *)compz2 = 'V';
#line 884 "sgghd3.f"
	}
#line 886 "sgghd3.f"
    }

#line 888 "sgghd3.f"
    if (jcol < *ihi) {
#line 888 "sgghd3.f"
	sgghrd_(compq2, compz2, n, &jcol, ihi, &a[a_offset], lda, &b[b_offset]
		, ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &ierr, (ftnlen)
		1, (ftnlen)1);
#line 888 "sgghd3.f"
    }
#line 891 "sgghd3.f"
    work[1] = (doublereal) lwkopt;

#line 893 "sgghd3.f"
    return 0;

/*     End of SGGHD3 */

} /* sgghd3_ */

