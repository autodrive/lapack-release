#line 1 "cgghrd.f"
/* cgghrd.f -- translated by f2c (version 20100827).
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

#line 1 "cgghrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CGGHRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgghrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgghrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgghrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/*                          LDQ, Z, LDZ, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ, COMPZ */
/*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N */
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
/* > CGGHRD reduces a pair of complex matrices (A,B) to generalized upper */
/* > Hessenberg form using unitary transformations, where A is a */
/* > general matrix and B is upper triangular.  The form of the generalized */
/* > eigenvalue problem is */
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
/* > original equation A*x = lambda*B*x, then CGGHRD reduces the original */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int cgghrd_(char *compq, char *compz, integer *n, integer *
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
    static integer jcol;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer jrow;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex ctemp;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), xerbla_(char *, integer *, 
	    ftnlen);
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

#line 248 "cgghrd.f"
    /* Parameter adjustments */
#line 248 "cgghrd.f"
    a_dim1 = *lda;
#line 248 "cgghrd.f"
    a_offset = 1 + a_dim1;
#line 248 "cgghrd.f"
    a -= a_offset;
#line 248 "cgghrd.f"
    b_dim1 = *ldb;
#line 248 "cgghrd.f"
    b_offset = 1 + b_dim1;
#line 248 "cgghrd.f"
    b -= b_offset;
#line 248 "cgghrd.f"
    q_dim1 = *ldq;
#line 248 "cgghrd.f"
    q_offset = 1 + q_dim1;
#line 248 "cgghrd.f"
    q -= q_offset;
#line 248 "cgghrd.f"
    z_dim1 = *ldz;
#line 248 "cgghrd.f"
    z_offset = 1 + z_dim1;
#line 248 "cgghrd.f"
    z__ -= z_offset;
#line 248 "cgghrd.f"

#line 248 "cgghrd.f"
    /* Function Body */
#line 248 "cgghrd.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 249 "cgghrd.f"
	ilq = FALSE_;
#line 250 "cgghrd.f"
	icompq = 1;
#line 251 "cgghrd.f"
    } else if (lsame_(compq, "V", (ftnlen)1, (ftnlen)1)) {
#line 252 "cgghrd.f"
	ilq = TRUE_;
#line 253 "cgghrd.f"
	icompq = 2;
#line 254 "cgghrd.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 255 "cgghrd.f"
	ilq = TRUE_;
#line 256 "cgghrd.f"
	icompq = 3;
#line 257 "cgghrd.f"
    } else {
#line 258 "cgghrd.f"
	icompq = 0;
#line 259 "cgghrd.f"
    }

/*     Decode COMPZ */

#line 263 "cgghrd.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 264 "cgghrd.f"
	ilz = FALSE_;
#line 265 "cgghrd.f"
	icompz = 1;
#line 266 "cgghrd.f"
    } else if (lsame_(compz, "V", (ftnlen)1, (ftnlen)1)) {
#line 267 "cgghrd.f"
	ilz = TRUE_;
#line 268 "cgghrd.f"
	icompz = 2;
#line 269 "cgghrd.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 270 "cgghrd.f"
	ilz = TRUE_;
#line 271 "cgghrd.f"
	icompz = 3;
#line 272 "cgghrd.f"
    } else {
#line 273 "cgghrd.f"
	icompz = 0;
#line 274 "cgghrd.f"
    }

/*     Test the input parameters. */

#line 278 "cgghrd.f"
    *info = 0;
#line 279 "cgghrd.f"
    if (icompq <= 0) {
#line 280 "cgghrd.f"
	*info = -1;
#line 281 "cgghrd.f"
    } else if (icompz <= 0) {
#line 282 "cgghrd.f"
	*info = -2;
#line 283 "cgghrd.f"
    } else if (*n < 0) {
#line 284 "cgghrd.f"
	*info = -3;
#line 285 "cgghrd.f"
    } else if (*ilo < 1) {
#line 286 "cgghrd.f"
	*info = -4;
#line 287 "cgghrd.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 288 "cgghrd.f"
	*info = -5;
#line 289 "cgghrd.f"
    } else if (*lda < max(1,*n)) {
#line 290 "cgghrd.f"
	*info = -7;
#line 291 "cgghrd.f"
    } else if (*ldb < max(1,*n)) {
#line 292 "cgghrd.f"
	*info = -9;
#line 293 "cgghrd.f"
    } else if (ilq && *ldq < *n || *ldq < 1) {
#line 294 "cgghrd.f"
	*info = -11;
#line 295 "cgghrd.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 296 "cgghrd.f"
	*info = -13;
#line 297 "cgghrd.f"
    }
#line 298 "cgghrd.f"
    if (*info != 0) {
#line 299 "cgghrd.f"
	i__1 = -(*info);
#line 299 "cgghrd.f"
	xerbla_("CGGHRD", &i__1, (ftnlen)6);
#line 300 "cgghrd.f"
	return 0;
#line 301 "cgghrd.f"
    }

/*     Initialize Q and Z if desired. */

#line 305 "cgghrd.f"
    if (icompq == 3) {
#line 305 "cgghrd.f"
	claset_("Full", n, n, &c_b2, &c_b1, &q[q_offset], ldq, (ftnlen)4);
#line 305 "cgghrd.f"
    }
#line 307 "cgghrd.f"
    if (icompz == 3) {
#line 307 "cgghrd.f"
	claset_("Full", n, n, &c_b2, &c_b1, &z__[z_offset], ldz, (ftnlen)4);
#line 307 "cgghrd.f"
    }

/*     Quick return if possible */

#line 312 "cgghrd.f"
    if (*n <= 1) {
#line 312 "cgghrd.f"
	return 0;
#line 312 "cgghrd.f"
    }

/*     Zero out lower triangle of B */

#line 317 "cgghrd.f"
    i__1 = *n - 1;
#line 317 "cgghrd.f"
    for (jcol = 1; jcol <= i__1; ++jcol) {
#line 318 "cgghrd.f"
	i__2 = *n;
#line 318 "cgghrd.f"
	for (jrow = jcol + 1; jrow <= i__2; ++jrow) {
#line 319 "cgghrd.f"
	    i__3 = jrow + jcol * b_dim1;
#line 319 "cgghrd.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 320 "cgghrd.f"
/* L10: */
#line 320 "cgghrd.f"
	}
#line 321 "cgghrd.f"
/* L20: */
#line 321 "cgghrd.f"
    }

/*     Reduce A and B */

#line 325 "cgghrd.f"
    i__1 = *ihi - 2;
#line 325 "cgghrd.f"
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

#line 327 "cgghrd.f"
	i__2 = jcol + 2;
#line 327 "cgghrd.f"
	for (jrow = *ihi; jrow >= i__2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL) */

#line 331 "cgghrd.f"
	    i__3 = jrow - 1 + jcol * a_dim1;
#line 331 "cgghrd.f"
	    ctemp.r = a[i__3].r, ctemp.i = a[i__3].i;
#line 332 "cgghrd.f"
	    clartg_(&ctemp, &a[jrow + jcol * a_dim1], &c__, &s, &a[jrow - 1 + 
		    jcol * a_dim1]);
#line 334 "cgghrd.f"
	    i__3 = jrow + jcol * a_dim1;
#line 334 "cgghrd.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 335 "cgghrd.f"
	    i__3 = *n - jcol;
#line 335 "cgghrd.f"
	    crot_(&i__3, &a[jrow - 1 + (jcol + 1) * a_dim1], lda, &a[jrow + (
		    jcol + 1) * a_dim1], lda, &c__, &s);
#line 337 "cgghrd.f"
	    i__3 = *n + 2 - jrow;
#line 337 "cgghrd.f"
	    crot_(&i__3, &b[jrow - 1 + (jrow - 1) * b_dim1], ldb, &b[jrow + (
		    jrow - 1) * b_dim1], ldb, &c__, &s);
#line 339 "cgghrd.f"
	    if (ilq) {
#line 339 "cgghrd.f"
		d_cnjg(&z__1, &s);
#line 339 "cgghrd.f"
		crot_(n, &q[(jrow - 1) * q_dim1 + 1], &c__1, &q[jrow * q_dim1 
			+ 1], &c__1, &c__, &z__1);
#line 339 "cgghrd.f"
	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1) */

#line 345 "cgghrd.f"
	    i__3 = jrow + jrow * b_dim1;
#line 345 "cgghrd.f"
	    ctemp.r = b[i__3].r, ctemp.i = b[i__3].i;
#line 346 "cgghrd.f"
	    clartg_(&ctemp, &b[jrow + (jrow - 1) * b_dim1], &c__, &s, &b[jrow 
		    + jrow * b_dim1]);
#line 348 "cgghrd.f"
	    i__3 = jrow + (jrow - 1) * b_dim1;
#line 348 "cgghrd.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 349 "cgghrd.f"
	    crot_(ihi, &a[jrow * a_dim1 + 1], &c__1, &a[(jrow - 1) * a_dim1 + 
		    1], &c__1, &c__, &s);
#line 350 "cgghrd.f"
	    i__3 = jrow - 1;
#line 350 "cgghrd.f"
	    crot_(&i__3, &b[jrow * b_dim1 + 1], &c__1, &b[(jrow - 1) * b_dim1 
		    + 1], &c__1, &c__, &s);
#line 352 "cgghrd.f"
	    if (ilz) {
#line 352 "cgghrd.f"
		crot_(n, &z__[jrow * z_dim1 + 1], &c__1, &z__[(jrow - 1) * 
			z_dim1 + 1], &c__1, &c__, &s);
#line 352 "cgghrd.f"
	    }
#line 354 "cgghrd.f"
/* L30: */
#line 354 "cgghrd.f"
	}
#line 355 "cgghrd.f"
/* L40: */
#line 355 "cgghrd.f"
    }

#line 357 "cgghrd.f"
    return 0;

/*     End of CGGHRD */

} /* cgghrd_ */

