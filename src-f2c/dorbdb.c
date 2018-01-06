#line 1 "dorbdb.f"
/* dorbdb.f -- translated by f2c (version 20100827).
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

#line 1 "dorbdb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DORBDB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORBDB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, */
/*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, */
/*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIGNS, TRANS */
/*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, */
/*      $                   Q */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   PHI( * ), THETA( * ) */
/*       DOUBLE PRECISION   TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), */
/*      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ), */
/*      $                   X21( LDX21, * ), X22( LDX22, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORBDB simultaneously bidiagonalizes the blocks of an M-by-M */
/* > partitioned orthogonal matrix X: */
/* > */
/* >                                 [ B11 | B12 0  0 ] */
/* >     [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**T */
/* > X = [-----------] = [---------] [----------------] [---------]   . */
/* >     [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ] */
/* >                                 [  0  |  0  0  I ] */
/* > */
/* > X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is */
/* > not the case, then X must be transposed and/or permuted. This can be */
/* > done in constant time using the TRANS and SIGNS options. See DORCSD */
/* > for details.) */
/* > */
/* > The orthogonal matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by- */
/* > (M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. They are */
/* > represented implicitly by Householder vectors. */
/* > */
/* > B11, B12, B21, and B22 are Q-by-Q bidiagonal matrices represented */
/* > implicitly by angles THETA, PHI. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER */
/* >          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major */
/* >                      order; */
/* >          otherwise:  X, U1, U2, V1T, and V2T are stored in column- */
/* >                      major order. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGNS */
/* > \verbatim */
/* >          SIGNS is CHARACTER */
/* >          = 'O':      The lower-left block is made nonpositive (the */
/* >                      "other" convention); */
/* >          otherwise:  The upper-right block is made nonpositive (the */
/* >                      "default" convention). */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows and columns in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of rows in X11 and X12. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >          The number of columns in X11 and X21. 0 <= Q <= */
/* >          MIN(P,M-P,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* >          X11 is DOUBLE PRECISION array, dimension (LDX11,Q) */
/* >          On entry, the top-left block of the orthogonal matrix to be */
/* >          reduced. On exit, the form depends on TRANS: */
/* >          If TRANS = 'N', then */
/* >             the columns of tril(X11) specify reflectors for P1, */
/* >             the rows of triu(X11,1) specify reflectors for Q1; */
/* >          else TRANS = 'T', and */
/* >             the rows of triu(X11) specify reflectors for P1, */
/* >             the columns of tril(X11,-1) specify reflectors for Q1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* >          LDX11 is INTEGER */
/* >          The leading dimension of X11. If TRANS = 'N', then LDX11 >= */
/* >          P; else LDX11 >= Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X12 */
/* > \verbatim */
/* >          X12 is DOUBLE PRECISION array, dimension (LDX12,M-Q) */
/* >          On entry, the top-right block of the orthogonal matrix to */
/* >          be reduced. On exit, the form depends on TRANS: */
/* >          If TRANS = 'N', then */
/* >             the rows of triu(X12) specify the first P reflectors for */
/* >             Q2; */
/* >          else TRANS = 'T', and */
/* >             the columns of tril(X12) specify the first P reflectors */
/* >             for Q2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX12 */
/* > \verbatim */
/* >          LDX12 is INTEGER */
/* >          The leading dimension of X12. If TRANS = 'N', then LDX12 >= */
/* >          P; else LDX11 >= M-Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* >          X21 is DOUBLE PRECISION array, dimension (LDX21,Q) */
/* >          On entry, the bottom-left block of the orthogonal matrix to */
/* >          be reduced. On exit, the form depends on TRANS: */
/* >          If TRANS = 'N', then */
/* >             the columns of tril(X21) specify reflectors for P2; */
/* >          else TRANS = 'T', and */
/* >             the rows of triu(X21) specify reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* >          LDX21 is INTEGER */
/* >          The leading dimension of X21. If TRANS = 'N', then LDX21 >= */
/* >          M-P; else LDX21 >= Q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X22 */
/* > \verbatim */
/* >          X22 is DOUBLE PRECISION array, dimension (LDX22,M-Q) */
/* >          On entry, the bottom-right block of the orthogonal matrix to */
/* >          be reduced. On exit, the form depends on TRANS: */
/* >          If TRANS = 'N', then */
/* >             the rows of triu(X22(Q+1:M-P,P+1:M-Q)) specify the last */
/* >             M-P-Q reflectors for Q2, */
/* >          else TRANS = 'T', and */
/* >             the columns of tril(X22(P+1:M-Q,Q+1:M-P)) specify the last */
/* >             M-P-Q reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX22 */
/* > \verbatim */
/* >          LDX22 is INTEGER */
/* >          The leading dimension of X22. If TRANS = 'N', then LDX22 >= */
/* >          M-P; else LDX22 >= M-Q. */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* >          THETA is DOUBLE PRECISION array, dimension (Q) */
/* >          The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* >          be computed from the angles THETA and PHI. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* >          PHI is DOUBLE PRECISION array, dimension (Q-1) */
/* >          The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* >          be computed from the angles THETA and PHI. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* >          TAUP1 is DOUBLE PRECISION array, dimension (P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is DOUBLE PRECISION array, dimension (M-P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is DOUBLE PRECISION array, dimension (Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ2 */
/* > \verbatim */
/* >          TAUQ2 is DOUBLE PRECISION array, dimension (M-Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= M-Q. */
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

/* > \date November 2013 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The bidiagonal blocks B11, B12, B21, and B22 are represented */
/* >  implicitly by angles THETA(1), ..., THETA(Q) and PHI(1), ..., */
/* >  PHI(Q-1). B11 and B21 are upper bidiagonal, while B21 and B22 are */
/* >  lower bidiagonal. Every entry in each bidiagonal band is a product */
/* >  of a sine or cosine of a THETA with a sine or cosine of a PHI. See */
/* >  [1] or DORCSD for details. */
/* > */
/* >  P1, P2, Q1, and Q2 are represented as products of elementary */
/* >  reflectors. See DORCSD for details on generating P1, P2, Q1, and Q2 */
/* >  using DORGQR and DORGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dorbdb_(char *trans, char *signs, integer *m, integer *p,
	 integer *q, doublereal *x11, integer *ldx11, doublereal *x12, 
	integer *ldx12, doublereal *x21, integer *ldx21, doublereal *x22, 
	integer *ldx22, doublereal *theta, doublereal *phi, doublereal *taup1,
	 doublereal *taup2, doublereal *tauq1, doublereal *tauq2, doublereal *
	work, integer *lwork, integer *info, ftnlen trans_len, ftnlen 
	signs_len)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x12_dim1, x12_offset, x21_dim1, x21_offset, 
	    x22_dim1, x22_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    static logical colmajor;
    static integer lworkmin, lworkopt, i__;
    static doublereal z1, z2, z3, z4;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    static logical lquery;
    extern /* Subroutine */ int dlarfgp_(integer *, doublereal *, doublereal *
	    , integer *, doublereal *);


/*  -- LAPACK computational routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2013 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ==================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions */
/*     .. */
/*     .. Executable Statements .. */

/*     Test input arguments */

#line 335 "dorbdb.f"
    /* Parameter adjustments */
#line 335 "dorbdb.f"
    x11_dim1 = *ldx11;
#line 335 "dorbdb.f"
    x11_offset = 1 + x11_dim1;
#line 335 "dorbdb.f"
    x11 -= x11_offset;
#line 335 "dorbdb.f"
    x12_dim1 = *ldx12;
#line 335 "dorbdb.f"
    x12_offset = 1 + x12_dim1;
#line 335 "dorbdb.f"
    x12 -= x12_offset;
#line 335 "dorbdb.f"
    x21_dim1 = *ldx21;
#line 335 "dorbdb.f"
    x21_offset = 1 + x21_dim1;
#line 335 "dorbdb.f"
    x21 -= x21_offset;
#line 335 "dorbdb.f"
    x22_dim1 = *ldx22;
#line 335 "dorbdb.f"
    x22_offset = 1 + x22_dim1;
#line 335 "dorbdb.f"
    x22 -= x22_offset;
#line 335 "dorbdb.f"
    --theta;
#line 335 "dorbdb.f"
    --phi;
#line 335 "dorbdb.f"
    --taup1;
#line 335 "dorbdb.f"
    --taup2;
#line 335 "dorbdb.f"
    --tauq1;
#line 335 "dorbdb.f"
    --tauq2;
#line 335 "dorbdb.f"
    --work;
#line 335 "dorbdb.f"

#line 335 "dorbdb.f"
    /* Function Body */
#line 335 "dorbdb.f"
    *info = 0;
#line 336 "dorbdb.f"
    colmajor = ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 337 "dorbdb.f"
    if (! lsame_(signs, "O", (ftnlen)1, (ftnlen)1)) {
#line 338 "dorbdb.f"
	z1 = 1.;
#line 339 "dorbdb.f"
	z2 = 1.;
#line 340 "dorbdb.f"
	z3 = 1.;
#line 341 "dorbdb.f"
	z4 = 1.;
#line 342 "dorbdb.f"
    } else {
#line 343 "dorbdb.f"
	z1 = 1.;
#line 344 "dorbdb.f"
	z2 = -1.;
#line 345 "dorbdb.f"
	z3 = 1.;
#line 346 "dorbdb.f"
	z4 = -1.;
#line 347 "dorbdb.f"
    }
#line 348 "dorbdb.f"
    lquery = *lwork == -1;

#line 350 "dorbdb.f"
    if (*m < 0) {
#line 351 "dorbdb.f"
	*info = -3;
#line 352 "dorbdb.f"
    } else if (*p < 0 || *p > *m) {
#line 353 "dorbdb.f"
	*info = -4;
#line 354 "dorbdb.f"
    } else if (*q < 0 || *q > *p || *q > *m - *p || *q > *m - *q) {
#line 356 "dorbdb.f"
	*info = -5;
#line 357 "dorbdb.f"
    } else if (colmajor && *ldx11 < max(1,*p)) {
#line 358 "dorbdb.f"
	*info = -7;
#line 359 "dorbdb.f"
    } else if (! colmajor && *ldx11 < max(1,*q)) {
#line 360 "dorbdb.f"
	*info = -7;
#line 361 "dorbdb.f"
    } else if (colmajor && *ldx12 < max(1,*p)) {
#line 362 "dorbdb.f"
	*info = -9;
#line 363 "dorbdb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 363 "dorbdb.f"
	i__1 = 1, i__2 = *m - *q;
#line 363 "dorbdb.f"
	if (! colmajor && *ldx12 < max(i__1,i__2)) {
#line 364 "dorbdb.f"
	    *info = -9;
#line 365 "dorbdb.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 365 "dorbdb.f"
	    i__1 = 1, i__2 = *m - *p;
#line 365 "dorbdb.f"
	    if (colmajor && *ldx21 < max(i__1,i__2)) {
#line 366 "dorbdb.f"
		*info = -11;
#line 367 "dorbdb.f"
	    } else if (! colmajor && *ldx21 < max(1,*q)) {
#line 368 "dorbdb.f"
		*info = -11;
#line 369 "dorbdb.f"
	    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 369 "dorbdb.f"
		i__1 = 1, i__2 = *m - *p;
#line 369 "dorbdb.f"
		if (colmajor && *ldx22 < max(i__1,i__2)) {
#line 370 "dorbdb.f"
		    *info = -13;
#line 371 "dorbdb.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 371 "dorbdb.f"
		    i__1 = 1, i__2 = *m - *q;
#line 371 "dorbdb.f"
		    if (! colmajor && *ldx22 < max(i__1,i__2)) {
#line 372 "dorbdb.f"
			*info = -13;
#line 373 "dorbdb.f"
		    }
#line 373 "dorbdb.f"
		}
#line 373 "dorbdb.f"
	    }
#line 373 "dorbdb.f"
	}
#line 373 "dorbdb.f"
    }

/*     Compute workspace */

#line 377 "dorbdb.f"
    if (*info == 0) {
#line 378 "dorbdb.f"
	lworkopt = *m - *q;
#line 379 "dorbdb.f"
	lworkmin = *m - *q;
#line 380 "dorbdb.f"
	work[1] = (doublereal) lworkopt;
#line 381 "dorbdb.f"
	if (*lwork < lworkmin && ! lquery) {
#line 382 "dorbdb.f"
	    *info = -21;
#line 383 "dorbdb.f"
	}
#line 384 "dorbdb.f"
    }
#line 385 "dorbdb.f"
    if (*info != 0) {
#line 386 "dorbdb.f"
	i__1 = -(*info);
#line 386 "dorbdb.f"
	xerbla_("xORBDB", &i__1, (ftnlen)6);
#line 387 "dorbdb.f"
	return 0;
#line 388 "dorbdb.f"
    } else if (lquery) {
#line 389 "dorbdb.f"
	return 0;
#line 390 "dorbdb.f"
    }

/*     Handle column-major and row-major separately */

#line 394 "dorbdb.f"
    if (colmajor) {

/*        Reduce columns 1, ..., Q of X11, X12, X21, and X22 */

#line 398 "dorbdb.f"
	i__1 = *q;
#line 398 "dorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 400 "dorbdb.f"
	    if (i__ == 1) {
#line 401 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 401 "dorbdb.f"
		dscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 402 "dorbdb.f"
	    } else {
#line 403 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 403 "dorbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 403 "dorbdb.f"
		dscal_(&i__2, &d__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 404 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 404 "dorbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 404 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x12[i__ + (i__ - 1) * x12_dim1], &c__1, 
			&x11[i__ + i__ * x11_dim1], &c__1);
#line 406 "dorbdb.f"
	    }
#line 407 "dorbdb.f"
	    if (i__ == 1) {
#line 408 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 408 "dorbdb.f"
		dscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], &c__1);
#line 409 "dorbdb.f"
	    } else {
#line 410 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 410 "dorbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 410 "dorbdb.f"
		dscal_(&i__2, &d__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 411 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 411 "dorbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 411 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x22[i__ + (i__ - 1) * x22_dim1], &c__1, 
			&x21[i__ + i__ * x21_dim1], &c__1);
#line 413 "dorbdb.f"
	    }

#line 415 "dorbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 415 "dorbdb.f"
	    i__3 = *p - i__ + 1;
#line 415 "dorbdb.f"
	    theta[i__] = atan2(dnrm2_(&i__2, &x21[i__ + i__ * x21_dim1], &
		    c__1), dnrm2_(&i__3, &x11[i__ + i__ * x11_dim1], &c__1));

#line 418 "dorbdb.f"
	    if (*p > i__) {
#line 419 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 419 "dorbdb.f"
		dlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + 
			i__ * x11_dim1], &c__1, &taup1[i__]);
#line 420 "dorbdb.f"
	    } else if (*p == i__) {
#line 421 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 421 "dorbdb.f"
		dlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + i__ * 
			x11_dim1], &c__1, &taup1[i__]);
#line 422 "dorbdb.f"
	    }
#line 423 "dorbdb.f"
	    x11[i__ + i__ * x11_dim1] = 1.;
#line 424 "dorbdb.f"
	    if (*m - *p > i__) {
#line 425 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 425 "dorbdb.f"
		dlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + 
			i__ * x21_dim1], &c__1, &taup2[i__]);
#line 427 "dorbdb.f"
	    } else if (*m - *p == i__) {
#line 428 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 428 "dorbdb.f"
		dlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], &c__1, &taup2[i__]);
#line 429 "dorbdb.f"
	    }
#line 430 "dorbdb.f"
	    x21[i__ + i__ * x21_dim1] = 1.;

#line 432 "dorbdb.f"
	    if (*q > i__) {
#line 433 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 433 "dorbdb.f"
		i__3 = *q - i__;
#line 433 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			taup1[i__], &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &
			work[1], (ftnlen)1);
#line 435 "dorbdb.f"
	    }
#line 436 "dorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 437 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 437 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 437 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 439 "dorbdb.f"
	    }
#line 440 "dorbdb.f"
	    if (*q > i__) {
#line 441 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 441 "dorbdb.f"
		i__3 = *q - i__;
#line 441 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			taup2[i__], &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &
			work[1], (ftnlen)1);
#line 443 "dorbdb.f"
	    }
#line 444 "dorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 445 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 445 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 445 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 447 "dorbdb.f"
	    }

#line 449 "dorbdb.f"
	    if (i__ < *q) {
#line 450 "dorbdb.f"
		i__2 = *q - i__;
#line 450 "dorbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 450 "dorbdb.f"
		dscal_(&i__2, &d__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 452 "dorbdb.f"
		i__2 = *q - i__;
#line 452 "dorbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 452 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, 
			&x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 454 "dorbdb.f"
	    }
#line 455 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 455 "dorbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 455 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 456 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 456 "dorbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 456 "dorbdb.f"
	    daxpy_(&i__2, &d__1, &x22[i__ + i__ * x22_dim1], ldx22, &x12[i__ 
		    + i__ * x12_dim1], ldx12);

#line 459 "dorbdb.f"
	    if (i__ < *q) {
#line 459 "dorbdb.f"
		i__2 = *q - i__;
#line 459 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 459 "dorbdb.f"
		phi[i__] = atan2(dnrm2_(&i__2, &x11[i__ + (i__ + 1) * 
			x11_dim1], ldx11), dnrm2_(&i__3, &x12[i__ + i__ * 
			x12_dim1], ldx12));
#line 459 "dorbdb.f"
	    }

#line 463 "dorbdb.f"
	    if (i__ < *q) {
#line 464 "dorbdb.f"
		if (*q - i__ == 1) {
#line 465 "dorbdb.f"
		    i__2 = *q - i__;
#line 465 "dorbdb.f"
		    dlarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__]);
#line 467 "dorbdb.f"
		} else {
#line 468 "dorbdb.f"
		    i__2 = *q - i__;
#line 468 "dorbdb.f"
		    dlarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 2) * x11_dim1], ldx11, &tauq1[i__]);
#line 470 "dorbdb.f"
		}
#line 471 "dorbdb.f"
		x11[i__ + (i__ + 1) * x11_dim1] = 1.;
#line 472 "dorbdb.f"
	    }
#line 473 "dorbdb.f"
	    if (*q + i__ - 1 < *m) {
#line 474 "dorbdb.f"
		if (*m - *q == i__) {
#line 475 "dorbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 475 "dorbdb.f"
		    dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 
			    i__ * x12_dim1], ldx12, &tauq2[i__]);
#line 477 "dorbdb.f"
		} else {
#line 478 "dorbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 478 "dorbdb.f"
		    dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (
			    i__ + 1) * x12_dim1], ldx12, &tauq2[i__]);
#line 480 "dorbdb.f"
		}
#line 481 "dorbdb.f"
	    }
#line 482 "dorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 484 "dorbdb.f"
	    if (i__ < *q) {
#line 485 "dorbdb.f"
		i__2 = *p - i__;
#line 485 "dorbdb.f"
		i__3 = *q - i__;
#line 485 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 487 "dorbdb.f"
		i__2 = *m - *p - i__;
#line 487 "dorbdb.f"
		i__3 = *q - i__;
#line 487 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 489 "dorbdb.f"
	    }
#line 490 "dorbdb.f"
	    if (*p > i__) {
#line 491 "dorbdb.f"
		i__2 = *p - i__;
#line 491 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 491 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 493 "dorbdb.f"
	    }
#line 494 "dorbdb.f"
	    if (*m - *p > i__) {
#line 495 "dorbdb.f"
		i__2 = *m - *p - i__;
#line 495 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 495 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[i__ + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 497 "dorbdb.f"
	    }

#line 499 "dorbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 503 "dorbdb.f"
	i__1 = *p;
#line 503 "dorbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 505 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 505 "dorbdb.f"
	    d__1 = -z1 * z4;
#line 505 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 506 "dorbdb.f"
	    if (i__ >= *m - *q) {
#line 507 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 507 "dorbdb.f"
		dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], ldx12, &tauq2[i__]);
#line 509 "dorbdb.f"
	    } else {
#line 510 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 510 "dorbdb.f"
		dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 
			1) * x12_dim1], ldx12, &tauq2[i__]);
#line 512 "dorbdb.f"
	    }
#line 513 "dorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 515 "dorbdb.f"
	    if (*p > i__) {
#line 516 "dorbdb.f"
		i__2 = *p - i__;
#line 516 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 516 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 518 "dorbdb.f"
	    }
#line 519 "dorbdb.f"
	    if (*m - *p - *q >= 1) {
#line 519 "dorbdb.f"
		i__2 = *m - *p - *q;
#line 519 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 519 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[*q + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 519 "dorbdb.f"
	    }

#line 523 "dorbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 527 "dorbdb.f"
	i__1 = *m - *p - *q;
#line 527 "dorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 529 "dorbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 529 "dorbdb.f"
	    d__1 = z2 * z4;
#line 529 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22);
#line 530 "dorbdb.f"
	    if (i__ == *m - *p - *q) {
#line 531 "dorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 531 "dorbdb.f"
		dlarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*
			q + i__ + (*p + i__) * x22_dim1], ldx22, &tauq2[*p + 
			i__]);
#line 533 "dorbdb.f"
	    } else {
#line 534 "dorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 534 "dorbdb.f"
		dlarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*
			q + i__ + (*p + i__ + 1) * x22_dim1], ldx22, &tauq2[*
			p + i__]);
#line 536 "dorbdb.f"
	    }
#line 537 "dorbdb.f"
	    x22[*q + i__ + (*p + i__) * x22_dim1] = 1.;
#line 538 "dorbdb.f"
	    if (i__ < *m - *p - *q) {
#line 539 "dorbdb.f"
		i__2 = *m - *p - *q - i__;
#line 539 "dorbdb.f"
		i__3 = *m - *p - *q - i__ + 1;
#line 539 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x22[*q + i__ + (*p + i__) * 
			x22_dim1], ldx22, &tauq2[*p + i__], &x22[*q + i__ + 1 
			+ (*p + i__) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 541 "dorbdb.f"
	    }

#line 543 "dorbdb.f"
	}

#line 545 "dorbdb.f"
    } else {

/*        Reduce columns 1, ..., Q of X11, X12, X21, X22 */

#line 549 "dorbdb.f"
	i__1 = *q;
#line 549 "dorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 551 "dorbdb.f"
	    if (i__ == 1) {
#line 552 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 552 "dorbdb.f"
		dscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 553 "dorbdb.f"
	    } else {
#line 554 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 554 "dorbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 554 "dorbdb.f"
		dscal_(&i__2, &d__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 555 "dorbdb.f"
		i__2 = *p - i__ + 1;
#line 555 "dorbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 555 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x12[i__ - 1 + i__ * x12_dim1], ldx12, &
			x11[i__ + i__ * x11_dim1], ldx11);
#line 557 "dorbdb.f"
	    }
#line 558 "dorbdb.f"
	    if (i__ == 1) {
#line 559 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 559 "dorbdb.f"
		dscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], ldx21);
#line 560 "dorbdb.f"
	    } else {
#line 561 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 561 "dorbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 561 "dorbdb.f"
		dscal_(&i__2, &d__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 562 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 562 "dorbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 562 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x22[i__ - 1 + i__ * x22_dim1], ldx22, &
			x21[i__ + i__ * x21_dim1], ldx21);
#line 564 "dorbdb.f"
	    }

#line 566 "dorbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 566 "dorbdb.f"
	    i__3 = *p - i__ + 1;
#line 566 "dorbdb.f"
	    theta[i__] = atan2(dnrm2_(&i__2, &x21[i__ + i__ * x21_dim1], 
		    ldx21), dnrm2_(&i__3, &x11[i__ + i__ * x11_dim1], ldx11));

#line 569 "dorbdb.f"
	    i__2 = *p - i__ + 1;
#line 569 "dorbdb.f"
	    dlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) *
		     x11_dim1], ldx11, &taup1[i__]);
#line 570 "dorbdb.f"
	    x11[i__ + i__ * x11_dim1] = 1.;
#line 571 "dorbdb.f"
	    if (i__ == *m - *p) {
#line 572 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 572 "dorbdb.f"
		dlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], ldx21, &taup2[i__]);
#line 574 "dorbdb.f"
	    } else {
#line 575 "dorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 575 "dorbdb.f"
		dlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 
			1) * x21_dim1], ldx21, &taup2[i__]);
#line 577 "dorbdb.f"
	    }
#line 578 "dorbdb.f"
	    x21[i__ + i__ * x21_dim1] = 1.;

#line 580 "dorbdb.f"
	    if (*q > i__) {
#line 581 "dorbdb.f"
		i__2 = *q - i__;
#line 581 "dorbdb.f"
		i__3 = *p - i__ + 1;
#line 581 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
			taup1[i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &
			work[1], (ftnlen)1);
#line 583 "dorbdb.f"
	    }
#line 584 "dorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 585 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 585 "dorbdb.f"
		i__3 = *p - i__ + 1;
#line 585 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
			taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 587 "dorbdb.f"
	    }
#line 588 "dorbdb.f"
	    if (*q > i__) {
#line 589 "dorbdb.f"
		i__2 = *q - i__;
#line 589 "dorbdb.f"
		i__3 = *m - *p - i__ + 1;
#line 589 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
			taup2[i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &
			work[1], (ftnlen)1);
#line 591 "dorbdb.f"
	    }
#line 592 "dorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 593 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 593 "dorbdb.f"
		i__3 = *m - *p - i__ + 1;
#line 593 "dorbdb.f"
		dlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
			taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 595 "dorbdb.f"
	    }

#line 597 "dorbdb.f"
	    if (i__ < *q) {
#line 598 "dorbdb.f"
		i__2 = *q - i__;
#line 598 "dorbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 598 "dorbdb.f"
		dscal_(&i__2, &d__1, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 599 "dorbdb.f"
		i__2 = *q - i__;
#line 599 "dorbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 599 "dorbdb.f"
		daxpy_(&i__2, &d__1, &x21[i__ + 1 + i__ * x21_dim1], &c__1, &
			x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 601 "dorbdb.f"
	    }
#line 602 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 602 "dorbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 602 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 603 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 603 "dorbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 603 "dorbdb.f"
	    daxpy_(&i__2, &d__1, &x22[i__ + i__ * x22_dim1], &c__1, &x12[i__ 
		    + i__ * x12_dim1], &c__1);

#line 606 "dorbdb.f"
	    if (i__ < *q) {
#line 606 "dorbdb.f"
		i__2 = *q - i__;
#line 606 "dorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 606 "dorbdb.f"
		phi[i__] = atan2(dnrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1],
			 &c__1), dnrm2_(&i__3, &x12[i__ + i__ * x12_dim1], &
			c__1));
#line 606 "dorbdb.f"
	    }

#line 610 "dorbdb.f"
	    if (i__ < *q) {
#line 611 "dorbdb.f"
		if (*q - i__ == 1) {
#line 612 "dorbdb.f"
		    i__2 = *q - i__;
#line 612 "dorbdb.f"
		    dlarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ 
			    + 1 + i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 614 "dorbdb.f"
		} else {
#line 615 "dorbdb.f"
		    i__2 = *q - i__;
#line 615 "dorbdb.f"
		    dlarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ 
			    + 2 + i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 617 "dorbdb.f"
		}
#line 618 "dorbdb.f"
		x11[i__ + 1 + i__ * x11_dim1] = 1.;
#line 619 "dorbdb.f"
	    }
#line 620 "dorbdb.f"
	    if (*m - *q > i__) {
#line 621 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 621 "dorbdb.f"
		dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + 
			i__ * x12_dim1], &c__1, &tauq2[i__]);
#line 623 "dorbdb.f"
	    } else {
#line 624 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 624 "dorbdb.f"
		dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], &c__1, &tauq2[i__]);
#line 626 "dorbdb.f"
	    }
#line 627 "dorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 629 "dorbdb.f"
	    if (i__ < *q) {
#line 630 "dorbdb.f"
		i__2 = *q - i__;
#line 630 "dorbdb.f"
		i__3 = *p - i__;
#line 630 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 632 "dorbdb.f"
		i__2 = *q - i__;
#line 632 "dorbdb.f"
		i__3 = *m - *p - i__;
#line 632 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 634 "dorbdb.f"
	    }
#line 635 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 635 "dorbdb.f"
	    i__3 = *p - i__;
#line 635 "dorbdb.f"
	    dlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
		    tauq2[i__], &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &
		    work[1], (ftnlen)1);
#line 637 "dorbdb.f"
	    if (*m - *p - i__ > 0) {
#line 638 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 638 "dorbdb.f"
		i__3 = *m - *p - i__;
#line 638 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x22[i__ + (i__ + 1) * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 640 "dorbdb.f"
	    }

#line 642 "dorbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 646 "dorbdb.f"
	i__1 = *p;
#line 646 "dorbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 648 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 648 "dorbdb.f"
	    d__1 = -z1 * z4;
#line 648 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 649 "dorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 649 "dorbdb.f"
	    dlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 650 "dorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 652 "dorbdb.f"
	    if (*p > i__) {
#line 653 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 653 "dorbdb.f"
		i__3 = *p - i__;
#line 653 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 655 "dorbdb.f"
	    }
#line 656 "dorbdb.f"
	    if (*m - *p - *q >= 1) {
#line 656 "dorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 656 "dorbdb.f"
		i__3 = *m - *p - *q;
#line 656 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x22[i__ + (*q + 1) * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 656 "dorbdb.f"
	    }

#line 660 "dorbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 664 "dorbdb.f"
	i__1 = *m - *p - *q;
#line 664 "dorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 666 "dorbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 666 "dorbdb.f"
	    d__1 = z2 * z4;
#line 666 "dorbdb.f"
	    dscal_(&i__2, &d__1, &x22[*p + i__ + (*q + i__) * x22_dim1], &
		    c__1);
#line 667 "dorbdb.f"
	    if (*m - *p - *q == i__) {
#line 668 "dorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 668 "dorbdb.f"
		dlarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*
			p + i__ + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + 
			i__]);
#line 670 "dorbdb.f"
	    } else {
#line 671 "dorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 671 "dorbdb.f"
		dlarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*
			p + i__ + 1 + (*q + i__) * x22_dim1], &c__1, &tauq2[*
			p + i__]);
#line 673 "dorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 673 "dorbdb.f"
		i__3 = *m - *p - *q - i__;
#line 673 "dorbdb.f"
		dlarf_("L", &i__2, &i__3, &x22[*p + i__ + (*q + i__) * 
			x22_dim1], &c__1, &tauq2[*p + i__], &x22[*p + i__ + (*
			q + i__ + 1) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 675 "dorbdb.f"
	    }
#line 676 "dorbdb.f"
	    x22[*p + i__ + (*q + i__) * x22_dim1] = 1.;

#line 678 "dorbdb.f"
	}

#line 680 "dorbdb.f"
    }

#line 682 "dorbdb.f"
    return 0;

/*     End of DORBDB */

} /* dorbdb_ */

