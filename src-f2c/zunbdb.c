#line 1 "zunbdb.f"
/* zunbdb.f -- translated by f2c (version 20100827).
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

#line 1 "zunbdb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZUNBDB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNBDB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunbdb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunbdb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunbdb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, */
/*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, */
/*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIGNS, TRANS */
/*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, */
/*      $                   Q */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   PHI( * ), THETA( * ) */
/*       COMPLEX*16         TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), */
/*      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ), */
/*      $                   X21( LDX21, * ), X22( LDX22, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNBDB simultaneously bidiagonalizes the blocks of an M-by-M */
/* > partitioned unitary matrix X: */
/* > */
/* >                                 [ B11 | B12 0  0 ] */
/* >     [ X11 | X12 ]   [ P1 |    ] [  0  |  0 -I  0 ] [ Q1 |    ]**H */
/* > X = [-----------] = [---------] [----------------] [---------]   . */
/* >     [ X21 | X22 ]   [    | P2 ] [ B21 | B22 0  0 ] [    | Q2 ] */
/* >                                 [  0  |  0  0  I ] */
/* > */
/* > X11 is P-by-Q. Q must be no larger than P, M-P, or M-Q. (If this is */
/* > not the case, then X must be transposed and/or permuted. This can be */
/* > done in constant time using the TRANS and SIGNS options. See ZUNCSD */
/* > for details.) */
/* > */
/* > The unitary matrices P1, P2, Q1, and Q2 are P-by-P, (M-P)-by- */
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
/* >          X11 is COMPLEX*16 array, dimension (LDX11,Q) */
/* >          On entry, the top-left block of the unitary matrix to be */
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
/* >          X12 is COMPLEX*16 array, dimension (LDX12,M-Q) */
/* >          On entry, the top-right block of the unitary matrix to */
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
/* >          X21 is COMPLEX*16 array, dimension (LDX21,Q) */
/* >          On entry, the bottom-left block of the unitary matrix to */
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
/* >          X22 is COMPLEX*16 array, dimension (LDX22,M-Q) */
/* >          On entry, the bottom-right block of the unitary matrix to */
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
/* >          TAUP1 is COMPLEX*16 array, dimension (P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is COMPLEX*16 array, dimension (M-P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is COMPLEX*16 array, dimension (Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ2 */
/* > \verbatim */
/* >          TAUQ2 is COMPLEX*16 array, dimension (M-Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
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

/* > \ingroup complex16OTHERcomputational */

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
/* >  [1] or ZUNCSD for details. */
/* > */
/* >  P1, P2, Q1, and Q2 are represented as products of elementary */
/* >  reflectors. See ZUNCSD for details on generating P1, P2, Q1, and Q2 */
/* >  using ZUNGQR and ZUNGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zunbdb_(char *trans, char *signs, integer *m, integer *p,
	 integer *q, doublecomplex *x11, integer *ldx11, doublecomplex *x12, 
	integer *ldx12, doublecomplex *x21, integer *ldx21, doublecomplex *
	x22, integer *ldx22, doublereal *theta, doublereal *phi, 
	doublecomplex *taup1, doublecomplex *taup2, doublecomplex *tauq1, 
	doublecomplex *tauq2, doublecomplex *work, integer *lwork, integer *
	info, ftnlen trans_len, ftnlen signs_len)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x12_dim1, x12_offset, x21_dim1, x21_offset, 
	    x22_dim1, x22_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), atan2(doublereal, doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static logical colmajor;
    static integer lworkmin, lworkopt, i__;
    static doublereal z1, z2, z3, z4;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlacgv_(
	    integer *, doublecomplex *, integer *);
    static logical lquery;
    extern /* Subroutine */ int zlarfgp_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);


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

#line 338 "zunbdb.f"
    /* Parameter adjustments */
#line 338 "zunbdb.f"
    x11_dim1 = *ldx11;
#line 338 "zunbdb.f"
    x11_offset = 1 + x11_dim1;
#line 338 "zunbdb.f"
    x11 -= x11_offset;
#line 338 "zunbdb.f"
    x12_dim1 = *ldx12;
#line 338 "zunbdb.f"
    x12_offset = 1 + x12_dim1;
#line 338 "zunbdb.f"
    x12 -= x12_offset;
#line 338 "zunbdb.f"
    x21_dim1 = *ldx21;
#line 338 "zunbdb.f"
    x21_offset = 1 + x21_dim1;
#line 338 "zunbdb.f"
    x21 -= x21_offset;
#line 338 "zunbdb.f"
    x22_dim1 = *ldx22;
#line 338 "zunbdb.f"
    x22_offset = 1 + x22_dim1;
#line 338 "zunbdb.f"
    x22 -= x22_offset;
#line 338 "zunbdb.f"
    --theta;
#line 338 "zunbdb.f"
    --phi;
#line 338 "zunbdb.f"
    --taup1;
#line 338 "zunbdb.f"
    --taup2;
#line 338 "zunbdb.f"
    --tauq1;
#line 338 "zunbdb.f"
    --tauq2;
#line 338 "zunbdb.f"
    --work;
#line 338 "zunbdb.f"

#line 338 "zunbdb.f"
    /* Function Body */
#line 338 "zunbdb.f"
    *info = 0;
#line 339 "zunbdb.f"
    colmajor = ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 340 "zunbdb.f"
    if (! lsame_(signs, "O", (ftnlen)1, (ftnlen)1)) {
#line 341 "zunbdb.f"
	z1 = 1.;
#line 342 "zunbdb.f"
	z2 = 1.;
#line 343 "zunbdb.f"
	z3 = 1.;
#line 344 "zunbdb.f"
	z4 = 1.;
#line 345 "zunbdb.f"
    } else {
#line 346 "zunbdb.f"
	z1 = 1.;
#line 347 "zunbdb.f"
	z2 = -1.;
#line 348 "zunbdb.f"
	z3 = 1.;
#line 349 "zunbdb.f"
	z4 = -1.;
#line 350 "zunbdb.f"
    }
#line 351 "zunbdb.f"
    lquery = *lwork == -1;

#line 353 "zunbdb.f"
    if (*m < 0) {
#line 354 "zunbdb.f"
	*info = -3;
#line 355 "zunbdb.f"
    } else if (*p < 0 || *p > *m) {
#line 356 "zunbdb.f"
	*info = -4;
#line 357 "zunbdb.f"
    } else if (*q < 0 || *q > *p || *q > *m - *p || *q > *m - *q) {
#line 359 "zunbdb.f"
	*info = -5;
#line 360 "zunbdb.f"
    } else if (colmajor && *ldx11 < max(1,*p)) {
#line 361 "zunbdb.f"
	*info = -7;
#line 362 "zunbdb.f"
    } else if (! colmajor && *ldx11 < max(1,*q)) {
#line 363 "zunbdb.f"
	*info = -7;
#line 364 "zunbdb.f"
    } else if (colmajor && *ldx12 < max(1,*p)) {
#line 365 "zunbdb.f"
	*info = -9;
#line 366 "zunbdb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 366 "zunbdb.f"
	i__1 = 1, i__2 = *m - *q;
#line 366 "zunbdb.f"
	if (! colmajor && *ldx12 < max(i__1,i__2)) {
#line 367 "zunbdb.f"
	    *info = -9;
#line 368 "zunbdb.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 368 "zunbdb.f"
	    i__1 = 1, i__2 = *m - *p;
#line 368 "zunbdb.f"
	    if (colmajor && *ldx21 < max(i__1,i__2)) {
#line 369 "zunbdb.f"
		*info = -11;
#line 370 "zunbdb.f"
	    } else if (! colmajor && *ldx21 < max(1,*q)) {
#line 371 "zunbdb.f"
		*info = -11;
#line 372 "zunbdb.f"
	    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 372 "zunbdb.f"
		i__1 = 1, i__2 = *m - *p;
#line 372 "zunbdb.f"
		if (colmajor && *ldx22 < max(i__1,i__2)) {
#line 373 "zunbdb.f"
		    *info = -13;
#line 374 "zunbdb.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 374 "zunbdb.f"
		    i__1 = 1, i__2 = *m - *q;
#line 374 "zunbdb.f"
		    if (! colmajor && *ldx22 < max(i__1,i__2)) {
#line 375 "zunbdb.f"
			*info = -13;
#line 376 "zunbdb.f"
		    }
#line 376 "zunbdb.f"
		}
#line 376 "zunbdb.f"
	    }
#line 376 "zunbdb.f"
	}
#line 376 "zunbdb.f"
    }

/*     Compute workspace */

#line 380 "zunbdb.f"
    if (*info == 0) {
#line 381 "zunbdb.f"
	lworkopt = *m - *q;
#line 382 "zunbdb.f"
	lworkmin = *m - *q;
#line 383 "zunbdb.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 384 "zunbdb.f"
	if (*lwork < lworkmin && ! lquery) {
#line 385 "zunbdb.f"
	    *info = -21;
#line 386 "zunbdb.f"
	}
#line 387 "zunbdb.f"
    }
#line 388 "zunbdb.f"
    if (*info != 0) {
#line 389 "zunbdb.f"
	i__1 = -(*info);
#line 389 "zunbdb.f"
	xerbla_("xORBDB", &i__1, (ftnlen)6);
#line 390 "zunbdb.f"
	return 0;
#line 391 "zunbdb.f"
    } else if (lquery) {
#line 392 "zunbdb.f"
	return 0;
#line 393 "zunbdb.f"
    }

/*     Handle column-major and row-major separately */

#line 397 "zunbdb.f"
    if (colmajor) {

/*        Reduce columns 1, ..., Q of X11, X12, X21, and X22 */

#line 401 "zunbdb.f"
	i__1 = *q;
#line 401 "zunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 403 "zunbdb.f"
	    if (i__ == 1) {
#line 404 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 404 "zunbdb.f"
		z__1.r = z1, z__1.i = 0.;
#line 404 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 405 "zunbdb.f"
	    } else {
#line 406 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 406 "zunbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 406 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 406 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 408 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 408 "zunbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 408 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 408 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x12[i__ + (i__ - 1) * x12_dim1], &c__1, 
			&x11[i__ + i__ * x11_dim1], &c__1);
#line 410 "zunbdb.f"
	    }
#line 411 "zunbdb.f"
	    if (i__ == 1) {
#line 412 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 412 "zunbdb.f"
		z__1.r = z2, z__1.i = 0.;
#line 412 "zunbdb.f"
		zscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 413 "zunbdb.f"
	    } else {
#line 414 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 414 "zunbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 414 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 414 "zunbdb.f"
		zscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 416 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 416 "zunbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 416 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 416 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x22[i__ + (i__ - 1) * x22_dim1], &c__1, 
			&x21[i__ + i__ * x21_dim1], &c__1);
#line 418 "zunbdb.f"
	    }

#line 420 "zunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 420 "zunbdb.f"
	    i__3 = *p - i__ + 1;
#line 420 "zunbdb.f"
	    theta[i__] = atan2(dznrm2_(&i__2, &x21[i__ + i__ * x21_dim1], &
		    c__1), dznrm2_(&i__3, &x11[i__ + i__ * x11_dim1], &c__1));

#line 423 "zunbdb.f"
	    if (*p > i__) {
#line 424 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 424 "zunbdb.f"
		zlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + 
			i__ * x11_dim1], &c__1, &taup1[i__]);
#line 425 "zunbdb.f"
	    } else if (*p == i__) {
#line 426 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 426 "zunbdb.f"
		zlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + i__ * 
			x11_dim1], &c__1, &taup1[i__]);
#line 427 "zunbdb.f"
	    }
#line 428 "zunbdb.f"
	    i__2 = i__ + i__ * x11_dim1;
#line 428 "zunbdb.f"
	    x11[i__2].r = 1., x11[i__2].i = 0.;
#line 429 "zunbdb.f"
	    if (*m - *p > i__) {
#line 430 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 430 "zunbdb.f"
		zlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + 
			i__ * x21_dim1], &c__1, &taup2[i__]);
#line 432 "zunbdb.f"
	    } else if (*m - *p == i__) {
#line 433 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 433 "zunbdb.f"
		zlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], &c__1, &taup2[i__]);
#line 435 "zunbdb.f"
	    }
#line 436 "zunbdb.f"
	    i__2 = i__ + i__ * x21_dim1;
#line 436 "zunbdb.f"
	    x21[i__2].r = 1., x21[i__2].i = 0.;

#line 438 "zunbdb.f"
	    if (*q > i__) {
#line 439 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 439 "zunbdb.f"
		i__3 = *q - i__;
#line 439 "zunbdb.f"
		d_cnjg(&z__1, &taup1[i__]);
#line 439 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			z__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &work[
			1], (ftnlen)1);
#line 441 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 441 "zunbdb.f"
		i__3 = *q - i__;
#line 441 "zunbdb.f"
		d_cnjg(&z__1, &taup2[i__]);
#line 441 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			z__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[
			1], (ftnlen)1);
#line 443 "zunbdb.f"
	    }
#line 444 "zunbdb.f"
	    if (*m - *q + 1 > i__) {
#line 445 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 445 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 445 "zunbdb.f"
		d_cnjg(&z__1, &taup1[i__]);
#line 445 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			z__1, &x12[i__ + i__ * x12_dim1], ldx12, &work[1], (
			ftnlen)1);
#line 447 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 447 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 447 "zunbdb.f"
		d_cnjg(&z__1, &taup2[i__]);
#line 447 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			z__1, &x22[i__ + i__ * x22_dim1], ldx22, &work[1], (
			ftnlen)1);
#line 449 "zunbdb.f"
	    }

#line 451 "zunbdb.f"
	    if (i__ < *q) {
#line 452 "zunbdb.f"
		i__2 = *q - i__;
#line 452 "zunbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 452 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 452 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 454 "zunbdb.f"
		i__2 = *q - i__;
#line 454 "zunbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 454 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 454 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, 
			&x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 456 "zunbdb.f"
	    }
#line 457 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 457 "zunbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 457 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 457 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 459 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 459 "zunbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 459 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 459 "zunbdb.f"
	    zaxpy_(&i__2, &z__1, &x22[i__ + i__ * x22_dim1], ldx22, &x12[i__ 
		    + i__ * x12_dim1], ldx12);

#line 462 "zunbdb.f"
	    if (i__ < *q) {
#line 462 "zunbdb.f"
		i__2 = *q - i__;
#line 462 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 462 "zunbdb.f"
		phi[i__] = atan2(dznrm2_(&i__2, &x11[i__ + (i__ + 1) * 
			x11_dim1], ldx11), dznrm2_(&i__3, &x12[i__ + i__ * 
			x12_dim1], ldx12));
#line 462 "zunbdb.f"
	    }

#line 466 "zunbdb.f"
	    if (i__ < *q) {
#line 467 "zunbdb.f"
		i__2 = *q - i__;
#line 467 "zunbdb.f"
		zlacgv_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 468 "zunbdb.f"
		if (i__ == *q - 1) {
#line 469 "zunbdb.f"
		    i__2 = *q - i__;
#line 469 "zunbdb.f"
		    zlarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__]);
#line 471 "zunbdb.f"
		} else {
#line 472 "zunbdb.f"
		    i__2 = *q - i__;
#line 472 "zunbdb.f"
		    zlarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 2) * x11_dim1], ldx11, &tauq1[i__]);
#line 474 "zunbdb.f"
		}
#line 475 "zunbdb.f"
		i__2 = i__ + (i__ + 1) * x11_dim1;
#line 475 "zunbdb.f"
		x11[i__2].r = 1., x11[i__2].i = 0.;
#line 476 "zunbdb.f"
	    }
#line 477 "zunbdb.f"
	    if (*m - *q + 1 > i__) {
#line 478 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 478 "zunbdb.f"
		zlacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);
#line 479 "zunbdb.f"
		if (*m - *q == i__) {
#line 480 "zunbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 480 "zunbdb.f"
		    zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 
			    i__ * x12_dim1], ldx12, &tauq2[i__]);
#line 482 "zunbdb.f"
		} else {
#line 483 "zunbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 483 "zunbdb.f"
		    zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (
			    i__ + 1) * x12_dim1], ldx12, &tauq2[i__]);
#line 485 "zunbdb.f"
		}
#line 486 "zunbdb.f"
	    }
#line 487 "zunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 487 "zunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 489 "zunbdb.f"
	    if (i__ < *q) {
#line 490 "zunbdb.f"
		i__2 = *p - i__;
#line 490 "zunbdb.f"
		i__3 = *q - i__;
#line 490 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 492 "zunbdb.f"
		i__2 = *m - *p - i__;
#line 492 "zunbdb.f"
		i__3 = *q - i__;
#line 492 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 494 "zunbdb.f"
	    }
#line 495 "zunbdb.f"
	    if (*p > i__) {
#line 496 "zunbdb.f"
		i__2 = *p - i__;
#line 496 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 496 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 498 "zunbdb.f"
	    }
#line 499 "zunbdb.f"
	    if (*m - *p > i__) {
#line 500 "zunbdb.f"
		i__2 = *m - *p - i__;
#line 500 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 500 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[i__ + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 502 "zunbdb.f"
	    }

#line 504 "zunbdb.f"
	    if (i__ < *q) {
#line 504 "zunbdb.f"
		i__2 = *q - i__;
#line 504 "zunbdb.f"
		zlacgv_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 504 "zunbdb.f"
	    }
#line 506 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 506 "zunbdb.f"
	    zlacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);

#line 508 "zunbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 512 "zunbdb.f"
	i__1 = *p;
#line 512 "zunbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 514 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 514 "zunbdb.f"
	    d__1 = -z1 * z4;
#line 514 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 514 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 516 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 516 "zunbdb.f"
	    zlacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);
#line 517 "zunbdb.f"
	    if (i__ >= *m - *q) {
#line 518 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 518 "zunbdb.f"
		zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], ldx12, &tauq2[i__]);
#line 520 "zunbdb.f"
	    } else {
#line 521 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 521 "zunbdb.f"
		zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 
			1) * x12_dim1], ldx12, &tauq2[i__]);
#line 523 "zunbdb.f"
	    }
#line 524 "zunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 524 "zunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 526 "zunbdb.f"
	    if (*p > i__) {
#line 527 "zunbdb.f"
		i__2 = *p - i__;
#line 527 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 527 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 529 "zunbdb.f"
	    }
#line 530 "zunbdb.f"
	    if (*m - *p - *q >= 1) {
#line 530 "zunbdb.f"
		i__2 = *m - *p - *q;
#line 530 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 530 "zunbdb.f"
		zlarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[*q + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 530 "zunbdb.f"
	    }

#line 534 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 534 "zunbdb.f"
	    zlacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);

#line 536 "zunbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 540 "zunbdb.f"
	i__1 = *m - *p - *q;
#line 540 "zunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 542 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 542 "zunbdb.f"
	    d__1 = z2 * z4;
#line 542 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 542 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22);
#line 544 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 544 "zunbdb.f"
	    zlacgv_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22);
#line 545 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 545 "zunbdb.f"
	    zlarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*q + 
		    i__ + (*p + i__ + 1) * x22_dim1], ldx22, &tauq2[*p + i__])
		    ;
#line 547 "zunbdb.f"
	    i__2 = *q + i__ + (*p + i__) * x22_dim1;
#line 547 "zunbdb.f"
	    x22[i__2].r = 1., x22[i__2].i = 0.;
#line 548 "zunbdb.f"
	    i__2 = *m - *p - *q - i__;
#line 548 "zunbdb.f"
	    i__3 = *m - *p - *q - i__ + 1;
#line 548 "zunbdb.f"
	    zlarf_("R", &i__2, &i__3, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22, &tauq2[*p + i__], &x22[*q + i__ + 1 + (*p + i__) * 
		    x22_dim1], ldx22, &work[1], (ftnlen)1);

#line 551 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 551 "zunbdb.f"
	    zlacgv_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22);

#line 553 "zunbdb.f"
	}

#line 555 "zunbdb.f"
    } else {

/*        Reduce columns 1, ..., Q of X11, X12, X21, X22 */

#line 559 "zunbdb.f"
	i__1 = *q;
#line 559 "zunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 561 "zunbdb.f"
	    if (i__ == 1) {
#line 562 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 562 "zunbdb.f"
		z__1.r = z1, z__1.i = 0.;
#line 562 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 564 "zunbdb.f"
	    } else {
#line 565 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 565 "zunbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 565 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 565 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 567 "zunbdb.f"
		i__2 = *p - i__ + 1;
#line 567 "zunbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 567 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 567 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x12[i__ - 1 + i__ * x12_dim1], ldx12, &
			x11[i__ + i__ * x11_dim1], ldx11);
#line 569 "zunbdb.f"
	    }
#line 570 "zunbdb.f"
	    if (i__ == 1) {
#line 571 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 571 "zunbdb.f"
		z__1.r = z2, z__1.i = 0.;
#line 571 "zunbdb.f"
		zscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 573 "zunbdb.f"
	    } else {
#line 574 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 574 "zunbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 574 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 574 "zunbdb.f"
		zscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 576 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 576 "zunbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 576 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 576 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x22[i__ - 1 + i__ * x22_dim1], ldx22, &
			x21[i__ + i__ * x21_dim1], ldx21);
#line 578 "zunbdb.f"
	    }

#line 580 "zunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 580 "zunbdb.f"
	    i__3 = *p - i__ + 1;
#line 580 "zunbdb.f"
	    theta[i__] = atan2(dznrm2_(&i__2, &x21[i__ + i__ * x21_dim1], 
		    ldx21), dznrm2_(&i__3, &x11[i__ + i__ * x11_dim1], ldx11))
		    ;

#line 583 "zunbdb.f"
	    i__2 = *p - i__ + 1;
#line 583 "zunbdb.f"
	    zlacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 584 "zunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 584 "zunbdb.f"
	    zlacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);

#line 586 "zunbdb.f"
	    i__2 = *p - i__ + 1;
#line 586 "zunbdb.f"
	    zlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) *
		     x11_dim1], ldx11, &taup1[i__]);
#line 587 "zunbdb.f"
	    i__2 = i__ + i__ * x11_dim1;
#line 587 "zunbdb.f"
	    x11[i__2].r = 1., x11[i__2].i = 0.;
#line 588 "zunbdb.f"
	    if (i__ == *m - *p) {
#line 589 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 589 "zunbdb.f"
		zlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], ldx21, &taup2[i__]);
#line 591 "zunbdb.f"
	    } else {
#line 592 "zunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 592 "zunbdb.f"
		zlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 
			1) * x21_dim1], ldx21, &taup2[i__]);
#line 594 "zunbdb.f"
	    }
#line 595 "zunbdb.f"
	    i__2 = i__ + i__ * x21_dim1;
#line 595 "zunbdb.f"
	    x21[i__2].r = 1., x21[i__2].i = 0.;

#line 597 "zunbdb.f"
	    i__2 = *q - i__;
#line 597 "zunbdb.f"
	    i__3 = *p - i__ + 1;
#line 597 "zunbdb.f"
	    zlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
		    taup1[i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[
		    1], (ftnlen)1);
#line 599 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 599 "zunbdb.f"
	    i__3 = *p - i__ + 1;
#line 599 "zunbdb.f"
	    zlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
		    taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[1], (
		    ftnlen)1);
#line 601 "zunbdb.f"
	    i__2 = *q - i__;
#line 601 "zunbdb.f"
	    i__3 = *m - *p - i__ + 1;
#line 601 "zunbdb.f"
	    zlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
		    taup2[i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[
		    1], (ftnlen)1);
#line 603 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 603 "zunbdb.f"
	    i__3 = *m - *p - i__ + 1;
#line 603 "zunbdb.f"
	    zlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
		    taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[1], (
		    ftnlen)1);

#line 606 "zunbdb.f"
	    i__2 = *p - i__ + 1;
#line 606 "zunbdb.f"
	    zlacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 607 "zunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 607 "zunbdb.f"
	    zlacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);

#line 609 "zunbdb.f"
	    if (i__ < *q) {
#line 610 "zunbdb.f"
		i__2 = *q - i__;
#line 610 "zunbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 610 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 610 "zunbdb.f"
		zscal_(&i__2, &z__1, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 612 "zunbdb.f"
		i__2 = *q - i__;
#line 612 "zunbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 612 "zunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 612 "zunbdb.f"
		zaxpy_(&i__2, &z__1, &x21[i__ + 1 + i__ * x21_dim1], &c__1, &
			x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 614 "zunbdb.f"
	    }
#line 615 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 615 "zunbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 615 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 615 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 617 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 617 "zunbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 617 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 617 "zunbdb.f"
	    zaxpy_(&i__2, &z__1, &x22[i__ + i__ * x22_dim1], &c__1, &x12[i__ 
		    + i__ * x12_dim1], &c__1);

#line 620 "zunbdb.f"
	    if (i__ < *q) {
#line 620 "zunbdb.f"
		i__2 = *q - i__;
#line 620 "zunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 620 "zunbdb.f"
		phi[i__] = atan2(dznrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1]
			, &c__1), dznrm2_(&i__3, &x12[i__ + i__ * x12_dim1], &
			c__1));
#line 620 "zunbdb.f"
	    }

#line 624 "zunbdb.f"
	    if (i__ < *q) {
#line 625 "zunbdb.f"
		i__2 = *q - i__;
#line 625 "zunbdb.f"
		zlarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ + 2 
			+ i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 626 "zunbdb.f"
		i__2 = i__ + 1 + i__ * x11_dim1;
#line 626 "zunbdb.f"
		x11[i__2].r = 1., x11[i__2].i = 0.;
#line 627 "zunbdb.f"
	    }
#line 628 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 628 "zunbdb.f"
	    zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 629 "zunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 629 "zunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 631 "zunbdb.f"
	    if (i__ < *q) {
#line 632 "zunbdb.f"
		i__2 = *q - i__;
#line 632 "zunbdb.f"
		i__3 = *p - i__;
#line 632 "zunbdb.f"
		d_cnjg(&z__1, &tauq1[i__]);
#line 632 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &z__1, &x11[i__ + 1 + (i__ + 1) * x11_dim1], 
			ldx11, &work[1], (ftnlen)1);
#line 634 "zunbdb.f"
		i__2 = *q - i__;
#line 634 "zunbdb.f"
		i__3 = *m - *p - i__;
#line 634 "zunbdb.f"
		d_cnjg(&z__1, &tauq1[i__]);
#line 634 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &z__1, &x21[i__ + 1 + (i__ + 1) * x21_dim1], 
			ldx21, &work[1], (ftnlen)1);
#line 636 "zunbdb.f"
	    }
#line 637 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 637 "zunbdb.f"
	    i__3 = *p - i__;
#line 637 "zunbdb.f"
	    d_cnjg(&z__1, &tauq2[i__]);
#line 637 "zunbdb.f"
	    zlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
		    z__1, &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[1], (
		    ftnlen)1);
#line 639 "zunbdb.f"
	    if (*m - *p > i__) {
#line 640 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 640 "zunbdb.f"
		i__3 = *m - *p - i__;
#line 640 "zunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 640 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x22[i__ + (i__ + 1) * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 642 "zunbdb.f"
	    }

#line 644 "zunbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 648 "zunbdb.f"
	i__1 = *p;
#line 648 "zunbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 650 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 650 "zunbdb.f"
	    d__1 = -z1 * z4;
#line 650 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 650 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 651 "zunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 651 "zunbdb.f"
	    zlarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 652 "zunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 652 "zunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 654 "zunbdb.f"
	    if (*p > i__) {
#line 655 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 655 "zunbdb.f"
		i__3 = *p - i__;
#line 655 "zunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 655 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 657 "zunbdb.f"
	    }
#line 658 "zunbdb.f"
	    if (*m - *p - *q >= 1) {
#line 658 "zunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 658 "zunbdb.f"
		i__3 = *m - *p - *q;
#line 658 "zunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 658 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x22[i__ + (*q + 1) * x22_dim1], ldx22, &work[1]
			, (ftnlen)1);
#line 658 "zunbdb.f"
	    }

#line 662 "zunbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 666 "zunbdb.f"
	i__1 = *m - *p - *q;
#line 666 "zunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 668 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 668 "zunbdb.f"
	    d__1 = z2 * z4;
#line 668 "zunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 668 "zunbdb.f"
	    zscal_(&i__2, &z__1, &x22[*p + i__ + (*q + i__) * x22_dim1], &
		    c__1);
#line 670 "zunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 670 "zunbdb.f"
	    zlarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*p + 
		    i__ + 1 + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + i__])
		    ;
#line 672 "zunbdb.f"
	    i__2 = *p + i__ + (*q + i__) * x22_dim1;
#line 672 "zunbdb.f"
	    x22[i__2].r = 1., x22[i__2].i = 0.;

#line 674 "zunbdb.f"
	    if (*m - *p - *q != i__) {
#line 675 "zunbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 675 "zunbdb.f"
		i__3 = *m - *p - *q - i__;
#line 675 "zunbdb.f"
		d_cnjg(&z__1, &tauq2[*p + i__]);
#line 675 "zunbdb.f"
		zlarf_("L", &i__2, &i__3, &x22[*p + i__ + (*q + i__) * 
			x22_dim1], &c__1, &z__1, &x22[*p + i__ + (*q + i__ + 
			1) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 678 "zunbdb.f"
	    }

#line 680 "zunbdb.f"
	}

#line 682 "zunbdb.f"
    }

#line 684 "zunbdb.f"
    return 0;

/*     End of ZUNBDB */

} /* zunbdb_ */

