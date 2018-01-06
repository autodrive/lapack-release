#line 1 "dgghrd.f"
/* dgghrd.f -- translated by f2c (version 20100827).
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

#line 1 "dgghrd.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief \b DGGHRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGHRD reduces a pair of real matrices (A,B) to generalized upper */
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
/* > original equation A*x = lambda*B*x, then DGGHRD reduces the original */
/* > problem to generalized Hessenberg form. */
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

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  This routine reduces A to Hessenberg and B to triangular form by */
/* >  an unblocked reduction, as described in _Matrix_Computations_, */
/* >  by Golub and Van Loan (Johns Hopkins Press.) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, integer *info, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__, s;
    static logical ilq, ilz;
    static integer jcol;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jrow;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen);
    static integer icompq, icompz;


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode COMPQ */

#line 249 "dgghrd.f"
    /* Parameter adjustments */
#line 249 "dgghrd.f"
    a_dim1 = *lda;
#line 249 "dgghrd.f"
    a_offset = 1 + a_dim1;
#line 249 "dgghrd.f"
    a -= a_offset;
#line 249 "dgghrd.f"
    b_dim1 = *ldb;
#line 249 "dgghrd.f"
    b_offset = 1 + b_dim1;
#line 249 "dgghrd.f"
    b -= b_offset;
#line 249 "dgghrd.f"
    q_dim1 = *ldq;
#line 249 "dgghrd.f"
    q_offset = 1 + q_dim1;
#line 249 "dgghrd.f"
    q -= q_offset;
#line 249 "dgghrd.f"
    z_dim1 = *ldz;
#line 249 "dgghrd.f"
    z_offset = 1 + z_dim1;
#line 249 "dgghrd.f"
    z__ -= z_offset;
#line 249 "dgghrd.f"

#line 249 "dgghrd.f"
    /* Function Body */
#line 249 "dgghrd.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 250 "dgghrd.f"
	ilq = FALSE_;
#line 251 "dgghrd.f"
	icompq = 1;
#line 252 "dgghrd.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 253 "dgghrd.f"
	ilq = TRUE_;
#line 254 "dgghrd.f"
	icompq = 2;
#line 255 "dgghrd.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 256 "dgghrd.f"
	ilq = TRUE_;
#line 257 "dgghrd.f"
	icompq = 3;
#line 258 "dgghrd.f"
    } else {
#line 259 "dgghrd.f"
	icompq = 0;
#line 260 "dgghrd.f"
    }

/*     Decode COMPZ */

#line 264 "dgghrd.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 265 "dgghrd.f"
	ilz = FALSE_;
#line 266 "dgghrd.f"
	icompz = 1;
#line 267 "dgghrd.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 268 "dgghrd.f"
	ilz = TRUE_;
#line 269 "dgghrd.f"
	icompz = 2;
#line 270 "dgghrd.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 271 "dgghrd.f"
	ilz = TRUE_;
#line 272 "dgghrd.f"
	icompz = 3;
#line 273 "dgghrd.f"
    } else {
#line 274 "dgghrd.f"
	icompz = 0;
#line 275 "dgghrd.f"
    }

/*     Test the input parameters. */

#line 279 "dgghrd.f"
    *info = 0;
#line 280 "dgghrd.f"
    if (icompq <= 0) {
#line 281 "dgghrd.f"
	*info = -1;
#line 282 "dgghrd.f"
    } else if (icompz <= 0) {
#line 283 "dgghrd.f"
	*info = -2;
#line 284 "dgghrd.f"
    } else if (*n < 0) {
#line 285 "dgghrd.f"
	*info = -3;
#line 286 "dgghrd.f"
    } else if (*ilo < 1) {
#line 287 "dgghrd.f"
	*info = -4;
#line 288 "dgghrd.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 289 "dgghrd.f"
	*info = -5;
#line 290 "dgghrd.f"
    } else if (*lda < max(1,*n)) {
#line 291 "dgghrd.f"
	*info = -7;
#line 292 "dgghrd.f"
    } else if (*ldb < max(1,*n)) {
#line 293 "dgghrd.f"
	*info = -9;
#line 294 "dgghrd.f"
    } else if (ilq && *ldq < *n || *ldq < 1) {
#line 295 "dgghrd.f"
	*info = -11;
#line 296 "dgghrd.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 297 "dgghrd.f"
	*info = -13;
#line 298 "dgghrd.f"
    }
#line 299 "dgghrd.f"
    if (*info != 0) {
#line 300 "dgghrd.f"
	i__1 = -(*info);
#line 300 "dgghrd.f"
	xerbla_("DGGHRD", &i__1, (ftnlen)6);
#line 301 "dgghrd.f"
	return 0;
#line 302 "dgghrd.f"
    }

/*     Initialize Q and Z if desired. */

#line 306 "dgghrd.f"
    if (icompq == 3) {
#line 306 "dgghrd.f"
	dlaset_("Full", n, n, &c_b10, &c_b11, &q[q_offset], ldq, (ftnlen)4);
#line 306 "dgghrd.f"
    }
#line 308 "dgghrd.f"
    if (icompz == 3) {
#line 308 "dgghrd.f"
	dlaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz, (ftnlen)4);
#line 308 "dgghrd.f"
    }

/*     Quick return if possible */

#line 313 "dgghrd.f"
    if (*n <= 1) {
#line 313 "dgghrd.f"
	return 0;
#line 313 "dgghrd.f"
    }

/*     Zero out lower triangle of B */

#line 318 "dgghrd.f"
    i__1 = *n - 1;
#line 318 "dgghrd.f"
    for (jcol = 1; jcol <= i__1; ++jcol) {
#line 319 "dgghrd.f"
	i__2 = *n;
#line 319 "dgghrd.f"
	for (jrow = jcol + 1; jrow <= i__2; ++jrow) {
#line 320 "dgghrd.f"
	    b[jrow + jcol * b_dim1] = 0.;
#line 321 "dgghrd.f"
/* L10: */
#line 321 "dgghrd.f"
	}
#line 322 "dgghrd.f"
/* L20: */
#line 322 "dgghrd.f"
    }

/*     Reduce A and B */

#line 326 "dgghrd.f"
    i__1 = *ihi - 2;
#line 326 "dgghrd.f"
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

#line 328 "dgghrd.f"
	i__2 = jcol + 2;
#line 328 "dgghrd.f"
	for (jrow = *ihi; jrow >= i__2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL) */

#line 332 "dgghrd.f"
	    temp = a[jrow - 1 + jcol * a_dim1];
#line 333 "dgghrd.f"
	    dlartg_(&temp, &a[jrow + jcol * a_dim1], &c__, &s, &a[jrow - 1 + 
		    jcol * a_dim1]);
#line 335 "dgghrd.f"
	    a[jrow + jcol * a_dim1] = 0.;
#line 336 "dgghrd.f"
	    i__3 = *n - jcol;
#line 336 "dgghrd.f"
	    drot_(&i__3, &a[jrow - 1 + (jcol + 1) * a_dim1], lda, &a[jrow + (
		    jcol + 1) * a_dim1], lda, &c__, &s);
#line 338 "dgghrd.f"
	    i__3 = *n + 2 - jrow;
#line 338 "dgghrd.f"
	    drot_(&i__3, &b[jrow - 1 + (jrow - 1) * b_dim1], ldb, &b[jrow + (
		    jrow - 1) * b_dim1], ldb, &c__, &s);
#line 340 "dgghrd.f"
	    if (ilq) {
#line 340 "dgghrd.f"
		drot_(n, &q[(jrow - 1) * q_dim1 + 1], &c__1, &q[jrow * q_dim1 
			+ 1], &c__1, &c__, &s);
#line 340 "dgghrd.f"
	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1) */

#line 345 "dgghrd.f"
	    temp = b[jrow + jrow * b_dim1];
#line 346 "dgghrd.f"
	    dlartg_(&temp, &b[jrow + (jrow - 1) * b_dim1], &c__, &s, &b[jrow 
		    + jrow * b_dim1]);
#line 348 "dgghrd.f"
	    b[jrow + (jrow - 1) * b_dim1] = 0.;
#line 349 "dgghrd.f"
	    drot_(ihi, &a[jrow * a_dim1 + 1], &c__1, &a[(jrow - 1) * a_dim1 + 
		    1], &c__1, &c__, &s);
#line 350 "dgghrd.f"
	    i__3 = jrow - 1;
#line 350 "dgghrd.f"
	    drot_(&i__3, &b[jrow * b_dim1 + 1], &c__1, &b[(jrow - 1) * b_dim1 
		    + 1], &c__1, &c__, &s);
#line 352 "dgghrd.f"
	    if (ilz) {
#line 352 "dgghrd.f"
		drot_(n, &z__[jrow * z_dim1 + 1], &c__1, &z__[(jrow - 1) * 
			z_dim1 + 1], &c__1, &c__, &s);
#line 352 "dgghrd.f"
	    }
#line 354 "dgghrd.f"
/* L30: */
#line 354 "dgghrd.f"
	}
#line 355 "dgghrd.f"
/* L40: */
#line 355 "dgghrd.f"
    }

#line 357 "dgghrd.f"
    return 0;

/*     End of DGGHRD */

} /* dgghrd_ */

