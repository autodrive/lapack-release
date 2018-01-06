#line 1 "cunbdb.f"
/* cunbdb.f -- translated by f2c (version 20100827).
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

#line 1 "cunbdb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CUNBDB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNBDB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, */
/*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, */
/*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIGNS, TRANS */
/*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, */
/*      $                   Q */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               PHI( * ), THETA( * ) */
/*       COMPLEX            TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), */
/*      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ), */
/*      $                   X21( LDX21, * ), X22( LDX22, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNBDB simultaneously bidiagonalizes the blocks of an M-by-M */
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
/* > done in constant time using the TRANS and SIGNS options. See CUNCSD */
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
/* >          X11 is COMPLEX array, dimension (LDX11,Q) */
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
/* >          X12 is COMPLEX array, dimension (LDX12,M-Q) */
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
/* >          X21 is COMPLEX array, dimension (LDX21,Q) */
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
/* >          X22 is COMPLEX array, dimension (LDX22,M-Q) */
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
/* >          THETA is REAL array, dimension (Q) */
/* >          The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* >          be computed from the angles THETA and PHI. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* >          PHI is REAL array, dimension (Q-1) */
/* >          The entries of the bidiagonal blocks B11, B12, B21, B22 can */
/* >          be computed from the angles THETA and PHI. See Further */
/* >          Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* >          TAUP1 is COMPLEX array, dimension (P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is COMPLEX array, dimension (M-P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is COMPLEX array, dimension (Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ2 */
/* > \verbatim */
/* >          TAUQ2 is COMPLEX array, dimension (M-Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (LWORK) */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

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
/* >  [1] or CUNCSD for details. */
/* > */
/* >  P1, P2, Q1, and Q2 are represented as products of elementary */
/* >  reflectors. See CUNCSD for details on generating P1, P2, Q1, and Q2 */
/* >  using CUNGQR and CUNGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cunbdb_(char *trans, char *signs, integer *m, integer *p,
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
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int clarfgp_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 338 "cunbdb.f"
    /* Parameter adjustments */
#line 338 "cunbdb.f"
    x11_dim1 = *ldx11;
#line 338 "cunbdb.f"
    x11_offset = 1 + x11_dim1;
#line 338 "cunbdb.f"
    x11 -= x11_offset;
#line 338 "cunbdb.f"
    x12_dim1 = *ldx12;
#line 338 "cunbdb.f"
    x12_offset = 1 + x12_dim1;
#line 338 "cunbdb.f"
    x12 -= x12_offset;
#line 338 "cunbdb.f"
    x21_dim1 = *ldx21;
#line 338 "cunbdb.f"
    x21_offset = 1 + x21_dim1;
#line 338 "cunbdb.f"
    x21 -= x21_offset;
#line 338 "cunbdb.f"
    x22_dim1 = *ldx22;
#line 338 "cunbdb.f"
    x22_offset = 1 + x22_dim1;
#line 338 "cunbdb.f"
    x22 -= x22_offset;
#line 338 "cunbdb.f"
    --theta;
#line 338 "cunbdb.f"
    --phi;
#line 338 "cunbdb.f"
    --taup1;
#line 338 "cunbdb.f"
    --taup2;
#line 338 "cunbdb.f"
    --tauq1;
#line 338 "cunbdb.f"
    --tauq2;
#line 338 "cunbdb.f"
    --work;
#line 338 "cunbdb.f"

#line 338 "cunbdb.f"
    /* Function Body */
#line 338 "cunbdb.f"
    *info = 0;
#line 339 "cunbdb.f"
    colmajor = ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 340 "cunbdb.f"
    if (! lsame_(signs, "O", (ftnlen)1, (ftnlen)1)) {
#line 341 "cunbdb.f"
	z1 = 1.;
#line 342 "cunbdb.f"
	z2 = 1.;
#line 343 "cunbdb.f"
	z3 = 1.;
#line 344 "cunbdb.f"
	z4 = 1.;
#line 345 "cunbdb.f"
    } else {
#line 346 "cunbdb.f"
	z1 = 1.;
#line 347 "cunbdb.f"
	z2 = -1.;
#line 348 "cunbdb.f"
	z3 = 1.;
#line 349 "cunbdb.f"
	z4 = -1.;
#line 350 "cunbdb.f"
    }
#line 351 "cunbdb.f"
    lquery = *lwork == -1;

#line 353 "cunbdb.f"
    if (*m < 0) {
#line 354 "cunbdb.f"
	*info = -3;
#line 355 "cunbdb.f"
    } else if (*p < 0 || *p > *m) {
#line 356 "cunbdb.f"
	*info = -4;
#line 357 "cunbdb.f"
    } else if (*q < 0 || *q > *p || *q > *m - *p || *q > *m - *q) {
#line 359 "cunbdb.f"
	*info = -5;
#line 360 "cunbdb.f"
    } else if (colmajor && *ldx11 < max(1,*p)) {
#line 361 "cunbdb.f"
	*info = -7;
#line 362 "cunbdb.f"
    } else if (! colmajor && *ldx11 < max(1,*q)) {
#line 363 "cunbdb.f"
	*info = -7;
#line 364 "cunbdb.f"
    } else if (colmajor && *ldx12 < max(1,*p)) {
#line 365 "cunbdb.f"
	*info = -9;
#line 366 "cunbdb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 366 "cunbdb.f"
	i__1 = 1, i__2 = *m - *q;
#line 366 "cunbdb.f"
	if (! colmajor && *ldx12 < max(i__1,i__2)) {
#line 367 "cunbdb.f"
	    *info = -9;
#line 368 "cunbdb.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 368 "cunbdb.f"
	    i__1 = 1, i__2 = *m - *p;
#line 368 "cunbdb.f"
	    if (colmajor && *ldx21 < max(i__1,i__2)) {
#line 369 "cunbdb.f"
		*info = -11;
#line 370 "cunbdb.f"
	    } else if (! colmajor && *ldx21 < max(1,*q)) {
#line 371 "cunbdb.f"
		*info = -11;
#line 372 "cunbdb.f"
	    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 372 "cunbdb.f"
		i__1 = 1, i__2 = *m - *p;
#line 372 "cunbdb.f"
		if (colmajor && *ldx22 < max(i__1,i__2)) {
#line 373 "cunbdb.f"
		    *info = -13;
#line 374 "cunbdb.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 374 "cunbdb.f"
		    i__1 = 1, i__2 = *m - *q;
#line 374 "cunbdb.f"
		    if (! colmajor && *ldx22 < max(i__1,i__2)) {
#line 375 "cunbdb.f"
			*info = -13;
#line 376 "cunbdb.f"
		    }
#line 376 "cunbdb.f"
		}
#line 376 "cunbdb.f"
	    }
#line 376 "cunbdb.f"
	}
#line 376 "cunbdb.f"
    }

/*     Compute workspace */

#line 380 "cunbdb.f"
    if (*info == 0) {
#line 381 "cunbdb.f"
	lworkopt = *m - *q;
#line 382 "cunbdb.f"
	lworkmin = *m - *q;
#line 383 "cunbdb.f"
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
#line 384 "cunbdb.f"
	if (*lwork < lworkmin && ! lquery) {
#line 385 "cunbdb.f"
	    *info = -21;
#line 386 "cunbdb.f"
	}
#line 387 "cunbdb.f"
    }
#line 388 "cunbdb.f"
    if (*info != 0) {
#line 389 "cunbdb.f"
	i__1 = -(*info);
#line 389 "cunbdb.f"
	xerbla_("xORBDB", &i__1, (ftnlen)6);
#line 390 "cunbdb.f"
	return 0;
#line 391 "cunbdb.f"
    } else if (lquery) {
#line 392 "cunbdb.f"
	return 0;
#line 393 "cunbdb.f"
    }

/*     Handle column-major and row-major separately */

#line 397 "cunbdb.f"
    if (colmajor) {

/*        Reduce columns 1, ..., Q of X11, X12, X21, and X22 */

#line 401 "cunbdb.f"
	i__1 = *q;
#line 401 "cunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 403 "cunbdb.f"
	    if (i__ == 1) {
#line 404 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 404 "cunbdb.f"
		z__1.r = z1, z__1.i = 0.;
#line 404 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 405 "cunbdb.f"
	    } else {
#line 406 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 406 "cunbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 406 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 406 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 408 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 408 "cunbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 408 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 408 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x12[i__ + (i__ - 1) * x12_dim1], &c__1, 
			&x11[i__ + i__ * x11_dim1], &c__1);
#line 410 "cunbdb.f"
	    }
#line 411 "cunbdb.f"
	    if (i__ == 1) {
#line 412 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 412 "cunbdb.f"
		z__1.r = z2, z__1.i = 0.;
#line 412 "cunbdb.f"
		cscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 413 "cunbdb.f"
	    } else {
#line 414 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 414 "cunbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 414 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 414 "cunbdb.f"
		cscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 416 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 416 "cunbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 416 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 416 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x22[i__ + (i__ - 1) * x22_dim1], &c__1, 
			&x21[i__ + i__ * x21_dim1], &c__1);
#line 418 "cunbdb.f"
	    }

#line 420 "cunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 420 "cunbdb.f"
	    i__3 = *p - i__ + 1;
#line 420 "cunbdb.f"
	    theta[i__] = atan2(scnrm2_(&i__2, &x21[i__ + i__ * x21_dim1], &
		    c__1), scnrm2_(&i__3, &x11[i__ + i__ * x11_dim1], &c__1));

#line 423 "cunbdb.f"
	    if (*p > i__) {
#line 424 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 424 "cunbdb.f"
		clarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + 
			i__ * x11_dim1], &c__1, &taup1[i__]);
#line 425 "cunbdb.f"
	    } else if (*p == i__) {
#line 426 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 426 "cunbdb.f"
		clarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + i__ * 
			x11_dim1], &c__1, &taup1[i__]);
#line 427 "cunbdb.f"
	    }
#line 428 "cunbdb.f"
	    i__2 = i__ + i__ * x11_dim1;
#line 428 "cunbdb.f"
	    x11[i__2].r = 1., x11[i__2].i = 0.;
#line 429 "cunbdb.f"
	    if (*m - *p > i__) {
#line 430 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 430 "cunbdb.f"
		clarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + 
			i__ * x21_dim1], &c__1, &taup2[i__]);
#line 432 "cunbdb.f"
	    } else if (*m - *p == i__) {
#line 433 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 433 "cunbdb.f"
		clarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], &c__1, &taup2[i__]);
#line 435 "cunbdb.f"
	    }
#line 436 "cunbdb.f"
	    i__2 = i__ + i__ * x21_dim1;
#line 436 "cunbdb.f"
	    x21[i__2].r = 1., x21[i__2].i = 0.;

#line 438 "cunbdb.f"
	    if (*q > i__) {
#line 439 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 439 "cunbdb.f"
		i__3 = *q - i__;
#line 439 "cunbdb.f"
		d_cnjg(&z__1, &taup1[i__]);
#line 439 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			z__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &work[
			1], (ftnlen)1);
#line 441 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 441 "cunbdb.f"
		i__3 = *q - i__;
#line 441 "cunbdb.f"
		d_cnjg(&z__1, &taup2[i__]);
#line 441 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			z__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[
			1], (ftnlen)1);
#line 443 "cunbdb.f"
	    }
#line 444 "cunbdb.f"
	    if (*m - *q + 1 > i__) {
#line 445 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 445 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 445 "cunbdb.f"
		d_cnjg(&z__1, &taup1[i__]);
#line 445 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			z__1, &x12[i__ + i__ * x12_dim1], ldx12, &work[1], (
			ftnlen)1);
#line 447 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 447 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 447 "cunbdb.f"
		d_cnjg(&z__1, &taup2[i__]);
#line 447 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			z__1, &x22[i__ + i__ * x22_dim1], ldx22, &work[1], (
			ftnlen)1);
#line 449 "cunbdb.f"
	    }

#line 451 "cunbdb.f"
	    if (i__ < *q) {
#line 452 "cunbdb.f"
		i__2 = *q - i__;
#line 452 "cunbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 452 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 452 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 454 "cunbdb.f"
		i__2 = *q - i__;
#line 454 "cunbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 454 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 454 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, 
			&x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 456 "cunbdb.f"
	    }
#line 457 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 457 "cunbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 457 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 457 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 459 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 459 "cunbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 459 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 459 "cunbdb.f"
	    caxpy_(&i__2, &z__1, &x22[i__ + i__ * x22_dim1], ldx22, &x12[i__ 
		    + i__ * x12_dim1], ldx12);

#line 462 "cunbdb.f"
	    if (i__ < *q) {
#line 462 "cunbdb.f"
		i__2 = *q - i__;
#line 462 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 462 "cunbdb.f"
		phi[i__] = atan2(scnrm2_(&i__2, &x11[i__ + (i__ + 1) * 
			x11_dim1], ldx11), scnrm2_(&i__3, &x12[i__ + i__ * 
			x12_dim1], ldx12));
#line 462 "cunbdb.f"
	    }

#line 466 "cunbdb.f"
	    if (i__ < *q) {
#line 467 "cunbdb.f"
		i__2 = *q - i__;
#line 467 "cunbdb.f"
		clacgv_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 468 "cunbdb.f"
		if (i__ == *q - 1) {
#line 469 "cunbdb.f"
		    i__2 = *q - i__;
#line 469 "cunbdb.f"
		    clarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__]);
#line 471 "cunbdb.f"
		} else {
#line 472 "cunbdb.f"
		    i__2 = *q - i__;
#line 472 "cunbdb.f"
		    clarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 2) * x11_dim1], ldx11, &tauq1[i__]);
#line 474 "cunbdb.f"
		}
#line 475 "cunbdb.f"
		i__2 = i__ + (i__ + 1) * x11_dim1;
#line 475 "cunbdb.f"
		x11[i__2].r = 1., x11[i__2].i = 0.;
#line 476 "cunbdb.f"
	    }
#line 477 "cunbdb.f"
	    if (*m - *q + 1 > i__) {
#line 478 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 478 "cunbdb.f"
		clacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);
#line 479 "cunbdb.f"
		if (*m - *q == i__) {
#line 480 "cunbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 480 "cunbdb.f"
		    clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 
			    i__ * x12_dim1], ldx12, &tauq2[i__]);
#line 482 "cunbdb.f"
		} else {
#line 483 "cunbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 483 "cunbdb.f"
		    clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (
			    i__ + 1) * x12_dim1], ldx12, &tauq2[i__]);
#line 485 "cunbdb.f"
		}
#line 486 "cunbdb.f"
	    }
#line 487 "cunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 487 "cunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 489 "cunbdb.f"
	    if (i__ < *q) {
#line 490 "cunbdb.f"
		i__2 = *p - i__;
#line 490 "cunbdb.f"
		i__3 = *q - i__;
#line 490 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 492 "cunbdb.f"
		i__2 = *m - *p - i__;
#line 492 "cunbdb.f"
		i__3 = *q - i__;
#line 492 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 494 "cunbdb.f"
	    }
#line 495 "cunbdb.f"
	    if (*p > i__) {
#line 496 "cunbdb.f"
		i__2 = *p - i__;
#line 496 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 496 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 498 "cunbdb.f"
	    }
#line 499 "cunbdb.f"
	    if (*m - *p > i__) {
#line 500 "cunbdb.f"
		i__2 = *m - *p - i__;
#line 500 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 500 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[i__ + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 502 "cunbdb.f"
	    }

#line 504 "cunbdb.f"
	    if (i__ < *q) {
#line 504 "cunbdb.f"
		i__2 = *q - i__;
#line 504 "cunbdb.f"
		clacgv_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 504 "cunbdb.f"
	    }
#line 506 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 506 "cunbdb.f"
	    clacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);

#line 508 "cunbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 512 "cunbdb.f"
	i__1 = *p;
#line 512 "cunbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 514 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 514 "cunbdb.f"
	    d__1 = -z1 * z4;
#line 514 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 514 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 516 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 516 "cunbdb.f"
	    clacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);
#line 517 "cunbdb.f"
	    if (i__ >= *m - *q) {
#line 518 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 518 "cunbdb.f"
		clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], ldx12, &tauq2[i__]);
#line 520 "cunbdb.f"
	    } else {
#line 521 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 521 "cunbdb.f"
		clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 
			1) * x12_dim1], ldx12, &tauq2[i__]);
#line 523 "cunbdb.f"
	    }
#line 524 "cunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 524 "cunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 526 "cunbdb.f"
	    if (*p > i__) {
#line 527 "cunbdb.f"
		i__2 = *p - i__;
#line 527 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 527 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 529 "cunbdb.f"
	    }
#line 530 "cunbdb.f"
	    if (*m - *p - *q >= 1) {
#line 530 "cunbdb.f"
		i__2 = *m - *p - *q;
#line 530 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 530 "cunbdb.f"
		clarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[*q + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 530 "cunbdb.f"
	    }

#line 534 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 534 "cunbdb.f"
	    clacgv_(&i__2, &x12[i__ + i__ * x12_dim1], ldx12);

#line 536 "cunbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 540 "cunbdb.f"
	i__1 = *m - *p - *q;
#line 540 "cunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 542 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 542 "cunbdb.f"
	    d__1 = z2 * z4;
#line 542 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 542 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22);
#line 544 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 544 "cunbdb.f"
	    clacgv_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22);
#line 545 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 545 "cunbdb.f"
	    clarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*q + 
		    i__ + (*p + i__ + 1) * x22_dim1], ldx22, &tauq2[*p + i__])
		    ;
#line 547 "cunbdb.f"
	    i__2 = *q + i__ + (*p + i__) * x22_dim1;
#line 547 "cunbdb.f"
	    x22[i__2].r = 1., x22[i__2].i = 0.;
#line 548 "cunbdb.f"
	    i__2 = *m - *p - *q - i__;
#line 548 "cunbdb.f"
	    i__3 = *m - *p - *q - i__ + 1;
#line 548 "cunbdb.f"
	    clarf_("R", &i__2, &i__3, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22, &tauq2[*p + i__], &x22[*q + i__ + 1 + (*p + i__) * 
		    x22_dim1], ldx22, &work[1], (ftnlen)1);

#line 551 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 551 "cunbdb.f"
	    clacgv_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], ldx22);

#line 553 "cunbdb.f"
	}

#line 555 "cunbdb.f"
    } else {

/*        Reduce columns 1, ..., Q of X11, X12, X21, X22 */

#line 559 "cunbdb.f"
	i__1 = *q;
#line 559 "cunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 561 "cunbdb.f"
	    if (i__ == 1) {
#line 562 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 562 "cunbdb.f"
		z__1.r = z1, z__1.i = 0.;
#line 562 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 564 "cunbdb.f"
	    } else {
#line 565 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 565 "cunbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 565 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 565 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 567 "cunbdb.f"
		i__2 = *p - i__ + 1;
#line 567 "cunbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 567 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 567 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x12[i__ - 1 + i__ * x12_dim1], ldx12, &
			x11[i__ + i__ * x11_dim1], ldx11);
#line 569 "cunbdb.f"
	    }
#line 570 "cunbdb.f"
	    if (i__ == 1) {
#line 571 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 571 "cunbdb.f"
		z__1.r = z2, z__1.i = 0.;
#line 571 "cunbdb.f"
		cscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 573 "cunbdb.f"
	    } else {
#line 574 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 574 "cunbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 574 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 574 "cunbdb.f"
		cscal_(&i__2, &z__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 576 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 576 "cunbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 576 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 576 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x22[i__ - 1 + i__ * x22_dim1], ldx22, &
			x21[i__ + i__ * x21_dim1], ldx21);
#line 578 "cunbdb.f"
	    }

#line 580 "cunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 580 "cunbdb.f"
	    i__3 = *p - i__ + 1;
#line 580 "cunbdb.f"
	    theta[i__] = atan2(scnrm2_(&i__2, &x21[i__ + i__ * x21_dim1], 
		    ldx21), scnrm2_(&i__3, &x11[i__ + i__ * x11_dim1], ldx11))
		    ;

#line 583 "cunbdb.f"
	    i__2 = *p - i__ + 1;
#line 583 "cunbdb.f"
	    clacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 584 "cunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 584 "cunbdb.f"
	    clacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);

#line 586 "cunbdb.f"
	    i__2 = *p - i__ + 1;
#line 586 "cunbdb.f"
	    clarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) *
		     x11_dim1], ldx11, &taup1[i__]);
#line 587 "cunbdb.f"
	    i__2 = i__ + i__ * x11_dim1;
#line 587 "cunbdb.f"
	    x11[i__2].r = 1., x11[i__2].i = 0.;
#line 588 "cunbdb.f"
	    if (i__ == *m - *p) {
#line 589 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 589 "cunbdb.f"
		clarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], ldx21, &taup2[i__]);
#line 591 "cunbdb.f"
	    } else {
#line 592 "cunbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 592 "cunbdb.f"
		clarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 
			1) * x21_dim1], ldx21, &taup2[i__]);
#line 594 "cunbdb.f"
	    }
#line 595 "cunbdb.f"
	    i__2 = i__ + i__ * x21_dim1;
#line 595 "cunbdb.f"
	    x21[i__2].r = 1., x21[i__2].i = 0.;

#line 597 "cunbdb.f"
	    i__2 = *q - i__;
#line 597 "cunbdb.f"
	    i__3 = *p - i__ + 1;
#line 597 "cunbdb.f"
	    clarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
		    taup1[i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[
		    1], (ftnlen)1);
#line 599 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 599 "cunbdb.f"
	    i__3 = *p - i__ + 1;
#line 599 "cunbdb.f"
	    clarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
		    taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[1], (
		    ftnlen)1);
#line 601 "cunbdb.f"
	    i__2 = *q - i__;
#line 601 "cunbdb.f"
	    i__3 = *m - *p - i__ + 1;
#line 601 "cunbdb.f"
	    clarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
		    taup2[i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[
		    1], (ftnlen)1);
#line 603 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 603 "cunbdb.f"
	    i__3 = *m - *p - i__ + 1;
#line 603 "cunbdb.f"
	    clarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
		    taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[1], (
		    ftnlen)1);

#line 606 "cunbdb.f"
	    i__2 = *p - i__ + 1;
#line 606 "cunbdb.f"
	    clacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
#line 607 "cunbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 607 "cunbdb.f"
	    clacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);

#line 609 "cunbdb.f"
	    if (i__ < *q) {
#line 610 "cunbdb.f"
		i__2 = *q - i__;
#line 610 "cunbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 610 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 610 "cunbdb.f"
		cscal_(&i__2, &z__1, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 612 "cunbdb.f"
		i__2 = *q - i__;
#line 612 "cunbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 612 "cunbdb.f"
		z__1.r = d__1, z__1.i = 0.;
#line 612 "cunbdb.f"
		caxpy_(&i__2, &z__1, &x21[i__ + 1 + i__ * x21_dim1], &c__1, &
			x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 614 "cunbdb.f"
	    }
#line 615 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 615 "cunbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 615 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 615 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 617 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 617 "cunbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 617 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 617 "cunbdb.f"
	    caxpy_(&i__2, &z__1, &x22[i__ + i__ * x22_dim1], &c__1, &x12[i__ 
		    + i__ * x12_dim1], &c__1);

#line 620 "cunbdb.f"
	    if (i__ < *q) {
#line 620 "cunbdb.f"
		i__2 = *q - i__;
#line 620 "cunbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 620 "cunbdb.f"
		phi[i__] = atan2(scnrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1]
			, &c__1), scnrm2_(&i__3, &x12[i__ + i__ * x12_dim1], &
			c__1));
#line 620 "cunbdb.f"
	    }

#line 624 "cunbdb.f"
	    if (i__ < *q) {
#line 625 "cunbdb.f"
		i__2 = *q - i__;
#line 625 "cunbdb.f"
		clarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ + 2 
			+ i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 626 "cunbdb.f"
		i__2 = i__ + 1 + i__ * x11_dim1;
#line 626 "cunbdb.f"
		x11[i__2].r = 1., x11[i__2].i = 0.;
#line 627 "cunbdb.f"
	    }
#line 628 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 628 "cunbdb.f"
	    clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 629 "cunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 629 "cunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 631 "cunbdb.f"
	    if (i__ < *q) {
#line 632 "cunbdb.f"
		i__2 = *q - i__;
#line 632 "cunbdb.f"
		i__3 = *p - i__;
#line 632 "cunbdb.f"
		d_cnjg(&z__1, &tauq1[i__]);
#line 632 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &z__1, &x11[i__ + 1 + (i__ + 1) * x11_dim1], 
			ldx11, &work[1], (ftnlen)1);
#line 634 "cunbdb.f"
		i__2 = *q - i__;
#line 634 "cunbdb.f"
		i__3 = *m - *p - i__;
#line 634 "cunbdb.f"
		d_cnjg(&z__1, &tauq1[i__]);
#line 634 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &z__1, &x21[i__ + 1 + (i__ + 1) * x21_dim1], 
			ldx21, &work[1], (ftnlen)1);
#line 636 "cunbdb.f"
	    }
#line 637 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 637 "cunbdb.f"
	    i__3 = *p - i__;
#line 637 "cunbdb.f"
	    d_cnjg(&z__1, &tauq2[i__]);
#line 637 "cunbdb.f"
	    clarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
		    z__1, &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[1], (
		    ftnlen)1);
#line 640 "cunbdb.f"
	    if (*m - *p > i__) {
#line 641 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 641 "cunbdb.f"
		i__3 = *m - *p - i__;
#line 641 "cunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 641 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x22[i__ + (i__ + 1) * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 643 "cunbdb.f"
	    }
#line 644 "cunbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 648 "cunbdb.f"
	i__1 = *p;
#line 648 "cunbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 650 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 650 "cunbdb.f"
	    d__1 = -z1 * z4;
#line 650 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 650 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 651 "cunbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 651 "cunbdb.f"
	    clarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 652 "cunbdb.f"
	    i__2 = i__ + i__ * x12_dim1;
#line 652 "cunbdb.f"
	    x12[i__2].r = 1., x12[i__2].i = 0.;

#line 654 "cunbdb.f"
	    if (*p > i__) {
#line 655 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 655 "cunbdb.f"
		i__3 = *p - i__;
#line 655 "cunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 655 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 657 "cunbdb.f"
	    }
#line 658 "cunbdb.f"
	    if (*m - *p - *q >= 1) {
#line 658 "cunbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 658 "cunbdb.f"
		i__3 = *m - *p - *q;
#line 658 "cunbdb.f"
		d_cnjg(&z__1, &tauq2[i__]);
#line 658 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			z__1, &x22[i__ + (*q + 1) * x22_dim1], ldx22, &work[1]
			, (ftnlen)1);
#line 658 "cunbdb.f"
	    }

#line 662 "cunbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 666 "cunbdb.f"
	i__1 = *m - *p - *q;
#line 666 "cunbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 668 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 668 "cunbdb.f"
	    d__1 = z2 * z4;
#line 668 "cunbdb.f"
	    z__1.r = d__1, z__1.i = 0.;
#line 668 "cunbdb.f"
	    cscal_(&i__2, &z__1, &x22[*p + i__ + (*q + i__) * x22_dim1], &
		    c__1);
#line 670 "cunbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 670 "cunbdb.f"
	    clarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*p + 
		    i__ + 1 + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + i__])
		    ;
#line 672 "cunbdb.f"
	    i__2 = *p + i__ + (*q + i__) * x22_dim1;
#line 672 "cunbdb.f"
	    x22[i__2].r = 1., x22[i__2].i = 0.;
#line 673 "cunbdb.f"
	    if (*m - *p - *q != i__) {
#line 674 "cunbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 674 "cunbdb.f"
		i__3 = *m - *p - *q - i__;
#line 674 "cunbdb.f"
		d_cnjg(&z__1, &tauq2[*p + i__]);
#line 674 "cunbdb.f"
		clarf_("L", &i__2, &i__3, &x22[*p + i__ + (*q + i__) * 
			x22_dim1], &c__1, &z__1, &x22[*p + i__ + (*q + i__ + 
			1) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 677 "cunbdb.f"
	    }
#line 678 "cunbdb.f"
	}

#line 680 "cunbdb.f"
    }

#line 682 "cunbdb.f"
    return 0;

/*     End of CUNBDB */

} /* cunbdb_ */

