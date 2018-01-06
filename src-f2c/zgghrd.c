#line 1 "zgghrd.f"
/* zgghrd.f -- translated by f2c (version 20100827).
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

#line 1 "zgghrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZGGHRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgghrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgghrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgghrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGHRD reduces a pair of complex matrices (A,B) to generalized upper */
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
/* > original equation A*x = lambda*B*x, then ZGGHRD reduces the original */
/* > problem to generalized Hessenberg form. */
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

/* > \date November 2015 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  This routine reduces A to Hessenberg and B to triangular form by */
/* >  an unblocked reduction, as described in _Matrix_Computations_, */
/* >  by Golub and van Loan (Johns Hopkins Press). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__, 
	integer *ldz, integer *info, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal c__;
    static doublecomplex s;
    static logical ilq, ilz;
    static integer jcol, jrow;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex ctemp;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icompq, icompz;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), zlartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 248 "zgghrd.f"
    /* Parameter adjustments */
#line 248 "zgghrd.f"
    a_dim1 = *lda;
#line 248 "zgghrd.f"
    a_offset = 1 + a_dim1;
#line 248 "zgghrd.f"
    a -= a_offset;
#line 248 "zgghrd.f"
    b_dim1 = *ldb;
#line 248 "zgghrd.f"
    b_offset = 1 + b_dim1;
#line 248 "zgghrd.f"
    b -= b_offset;
#line 248 "zgghrd.f"
    q_dim1 = *ldq;
#line 248 "zgghrd.f"
    q_offset = 1 + q_dim1;
#line 248 "zgghrd.f"
    q -= q_offset;
#line 248 "zgghrd.f"
    z_dim1 = *ldz;
#line 248 "zgghrd.f"
    z_offset = 1 + z_dim1;
#line 248 "zgghrd.f"
    z__ -= z_offset;
#line 248 "zgghrd.f"

#line 248 "zgghrd.f"
    /* Function Body */
#line 248 "zgghrd.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 249 "zgghrd.f"
	ilq = FALSE_;
#line 250 "zgghrd.f"
	icompq = 1;
#line 251 "zgghrd.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 252 "zgghrd.f"
	ilq = TRUE_;
#line 253 "zgghrd.f"
	icompq = 2;
#line 254 "zgghrd.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 255 "zgghrd.f"
	ilq = TRUE_;
#line 256 "zgghrd.f"
	icompq = 3;
#line 257 "zgghrd.f"
    } else {
#line 258 "zgghrd.f"
	icompq = 0;
#line 259 "zgghrd.f"
    }

/*     Decode COMPZ */

#line 263 "zgghrd.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 264 "zgghrd.f"
	ilz = FALSE_;
#line 265 "zgghrd.f"
	icompz = 1;
#line 266 "zgghrd.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 267 "zgghrd.f"
	ilz = TRUE_;
#line 268 "zgghrd.f"
	icompz = 2;
#line 269 "zgghrd.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 270 "zgghrd.f"
	ilz = TRUE_;
#line 271 "zgghrd.f"
	icompz = 3;
#line 272 "zgghrd.f"
    } else {
#line 273 "zgghrd.f"
	icompz = 0;
#line 274 "zgghrd.f"
    }

/*     Test the input parameters. */

#line 278 "zgghrd.f"
    *info = 0;
#line 279 "zgghrd.f"
    if (icompq <= 0) {
#line 280 "zgghrd.f"
	*info = -1;
#line 281 "zgghrd.f"
    } else if (icompz <= 0) {
#line 282 "zgghrd.f"
	*info = -2;
#line 283 "zgghrd.f"
    } else if (*n < 0) {
#line 284 "zgghrd.f"
	*info = -3;
#line 285 "zgghrd.f"
    } else if (*ilo < 1) {
#line 286 "zgghrd.f"
	*info = -4;
#line 287 "zgghrd.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 288 "zgghrd.f"
	*info = -5;
#line 289 "zgghrd.f"
    } else if (*lda < max(1,*n)) {
#line 290 "zgghrd.f"
	*info = -7;
#line 291 "zgghrd.f"
    } else if (*ldb < max(1,*n)) {
#line 292 "zgghrd.f"
	*info = -9;
#line 293 "zgghrd.f"
    } else if (ilq && *ldq < *n || *ldq < 1) {
#line 294 "zgghrd.f"
	*info = -11;
#line 295 "zgghrd.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 296 "zgghrd.f"
	*info = -13;
#line 297 "zgghrd.f"
    }
#line 298 "zgghrd.f"
    if (*info != 0) {
#line 299 "zgghrd.f"
	i__1 = -(*info);
#line 299 "zgghrd.f"
	xerbla_("ZGGHRD", &i__1, (ftnlen)6);
#line 300 "zgghrd.f"
	return 0;
#line 301 "zgghrd.f"
    }

/*     Initialize Q and Z if desired. */

#line 305 "zgghrd.f"
    if (icompq == 3) {
#line 305 "zgghrd.f"
	zlaset_("Full", n, n, &c_b2, &c_b1, &q[q_offset], ldq, (ftnlen)4);
#line 305 "zgghrd.f"
    }
#line 307 "zgghrd.f"
    if (icompz == 3) {
#line 307 "zgghrd.f"
	zlaset_("Full", n, n, &c_b2, &c_b1, &z__[z_offset], ldz, (ftnlen)4);
#line 307 "zgghrd.f"
    }

/*     Quick return if possible */

#line 312 "zgghrd.f"
    if (*n <= 1) {
#line 312 "zgghrd.f"
	return 0;
#line 312 "zgghrd.f"
    }

/*     Zero out lower triangle of B */

#line 317 "zgghrd.f"
    i__1 = *n - 1;
#line 317 "zgghrd.f"
    for (jcol = 1; jcol <= i__1; ++jcol) {
#line 318 "zgghrd.f"
	i__2 = *n;
#line 318 "zgghrd.f"
	for (jrow = jcol + 1; jrow <= i__2; ++jrow) {
#line 319 "zgghrd.f"
	    i__3 = jrow + jcol * b_dim1;
#line 319 "zgghrd.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 320 "zgghrd.f"
/* L10: */
#line 320 "zgghrd.f"
	}
#line 321 "zgghrd.f"
/* L20: */
#line 321 "zgghrd.f"
    }

/*     Reduce A and B */

#line 325 "zgghrd.f"
    i__1 = *ihi - 2;
#line 325 "zgghrd.f"
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

#line 327 "zgghrd.f"
	i__2 = jcol + 2;
#line 327 "zgghrd.f"
	for (jrow = *ihi; jrow >= i__2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL) */

#line 331 "zgghrd.f"
	    i__3 = jrow - 1 + jcol * a_dim1;
#line 331 "zgghrd.f"
	    ctemp.r = a[i__3].r, ctemp.i = a[i__3].i;
#line 332 "zgghrd.f"
	    zlartg_(&ctemp, &a[jrow + jcol * a_dim1], &c__, &s, &a[jrow - 1 + 
		    jcol * a_dim1]);
#line 334 "zgghrd.f"
	    i__3 = jrow + jcol * a_dim1;
#line 334 "zgghrd.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 335 "zgghrd.f"
	    i__3 = *n - jcol;
#line 335 "zgghrd.f"
	    zrot_(&i__3, &a[jrow - 1 + (jcol + 1) * a_dim1], lda, &a[jrow + (
		    jcol + 1) * a_dim1], lda, &c__, &s);
#line 337 "zgghrd.f"
	    i__3 = *n + 2 - jrow;
#line 337 "zgghrd.f"
	    zrot_(&i__3, &b[jrow - 1 + (jrow - 1) * b_dim1], ldb, &b[jrow + (
		    jrow - 1) * b_dim1], ldb, &c__, &s);
#line 339 "zgghrd.f"
	    if (ilq) {
#line 339 "zgghrd.f"
		d_cnjg(&z__1, &s);
#line 339 "zgghrd.f"
		zrot_(n, &q[(jrow - 1) * q_dim1 + 1], &c__1, &q[jrow * q_dim1 
			+ 1], &c__1, &c__, &z__1);
#line 339 "zgghrd.f"
	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1) */

#line 345 "zgghrd.f"
	    i__3 = jrow + jrow * b_dim1;
#line 345 "zgghrd.f"
	    ctemp.r = b[i__3].r, ctemp.i = b[i__3].i;
#line 346 "zgghrd.f"
	    zlartg_(&ctemp, &b[jrow + (jrow - 1) * b_dim1], &c__, &s, &b[jrow 
		    + jrow * b_dim1]);
#line 348 "zgghrd.f"
	    i__3 = jrow + (jrow - 1) * b_dim1;
#line 348 "zgghrd.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 349 "zgghrd.f"
	    zrot_(ihi, &a[jrow * a_dim1 + 1], &c__1, &a[(jrow - 1) * a_dim1 + 
		    1], &c__1, &c__, &s);
#line 350 "zgghrd.f"
	    i__3 = jrow - 1;
#line 350 "zgghrd.f"
	    zrot_(&i__3, &b[jrow * b_dim1 + 1], &c__1, &b[(jrow - 1) * b_dim1 
		    + 1], &c__1, &c__, &s);
#line 352 "zgghrd.f"
	    if (ilz) {
#line 352 "zgghrd.f"
		zrot_(n, &z__[jrow * z_dim1 + 1], &c__1, &z__[(jrow - 1) * 
			z_dim1 + 1], &c__1, &c__, &s);
#line 352 "zgghrd.f"
	    }
#line 354 "zgghrd.f"
/* L30: */
#line 354 "zgghrd.f"
	}
#line 355 "zgghrd.f"
/* L40: */
#line 355 "zgghrd.f"
    }

#line 357 "zgghrd.f"
    return 0;

/*     End of ZGGHRD */

} /* zgghrd_ */

