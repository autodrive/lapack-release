#line 1 "sorbdb.f"
/* sorbdb.f -- translated by f2c (version 20100827).
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

#line 1 "sorbdb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SORBDB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORBDB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, */
/*                          X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, */
/*                          TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIGNS, TRANS */
/*       INTEGER            INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, */
/*      $                   Q */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               PHI( * ), THETA( * ) */
/*       REAL               TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), */
/*      $                   WORK( * ), X11( LDX11, * ), X12( LDX12, * ), */
/*      $                   X21( LDX21, * ), X22( LDX22, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORBDB simultaneously bidiagonalizes the blocks of an M-by-M */
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
/* > done in constant time using the TRANS and SIGNS options. See SORCSD */
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
/* >          X11 is REAL array, dimension (LDX11,Q) */
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
/* >          X12 is REAL array, dimension (LDX12,M-Q) */
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
/* >          X21 is REAL array, dimension (LDX21,Q) */
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
/* >          X22 is REAL array, dimension (LDX22,M-Q) */
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
/* >          TAUP1 is REAL array, dimension (P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is REAL array, dimension (M-P) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is REAL array, dimension (Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ2 */
/* > \verbatim */
/* >          TAUQ2 is REAL array, dimension (M-Q) */
/* >          The scalar factors of the elementary reflectors that define */
/* >          Q2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (LWORK) */
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

/* > \date November 2015 */

/* > \ingroup realOTHERcomputational */

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
/* >  [1] or SORCSD for details. */
/* > */
/* >  P1, P2, Q1, and Q2 are represented as products of elementary */
/* >  reflectors. See SORCSD for details on generating P1, P2, Q1, and Q2 */
/* >  using SORGQR and SORGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sorbdb_(char *trans, char *signs, integer *m, integer *p,
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), slarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen), saxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int slarfgp_(integer *, doublereal *, doublereal *
	    , integer *, doublereal *);


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 335 "sorbdb.f"
    /* Parameter adjustments */
#line 335 "sorbdb.f"
    x11_dim1 = *ldx11;
#line 335 "sorbdb.f"
    x11_offset = 1 + x11_dim1;
#line 335 "sorbdb.f"
    x11 -= x11_offset;
#line 335 "sorbdb.f"
    x12_dim1 = *ldx12;
#line 335 "sorbdb.f"
    x12_offset = 1 + x12_dim1;
#line 335 "sorbdb.f"
    x12 -= x12_offset;
#line 335 "sorbdb.f"
    x21_dim1 = *ldx21;
#line 335 "sorbdb.f"
    x21_offset = 1 + x21_dim1;
#line 335 "sorbdb.f"
    x21 -= x21_offset;
#line 335 "sorbdb.f"
    x22_dim1 = *ldx22;
#line 335 "sorbdb.f"
    x22_offset = 1 + x22_dim1;
#line 335 "sorbdb.f"
    x22 -= x22_offset;
#line 335 "sorbdb.f"
    --theta;
#line 335 "sorbdb.f"
    --phi;
#line 335 "sorbdb.f"
    --taup1;
#line 335 "sorbdb.f"
    --taup2;
#line 335 "sorbdb.f"
    --tauq1;
#line 335 "sorbdb.f"
    --tauq2;
#line 335 "sorbdb.f"
    --work;
#line 335 "sorbdb.f"

#line 335 "sorbdb.f"
    /* Function Body */
#line 335 "sorbdb.f"
    *info = 0;
#line 336 "sorbdb.f"
    colmajor = ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 337 "sorbdb.f"
    if (! lsame_(signs, "O", (ftnlen)1, (ftnlen)1)) {
#line 338 "sorbdb.f"
	z1 = 1.;
#line 339 "sorbdb.f"
	z2 = 1.;
#line 340 "sorbdb.f"
	z3 = 1.;
#line 341 "sorbdb.f"
	z4 = 1.;
#line 342 "sorbdb.f"
    } else {
#line 343 "sorbdb.f"
	z1 = 1.;
#line 344 "sorbdb.f"
	z2 = -1.;
#line 345 "sorbdb.f"
	z3 = 1.;
#line 346 "sorbdb.f"
	z4 = -1.;
#line 347 "sorbdb.f"
    }
#line 348 "sorbdb.f"
    lquery = *lwork == -1;

#line 350 "sorbdb.f"
    if (*m < 0) {
#line 351 "sorbdb.f"
	*info = -3;
#line 352 "sorbdb.f"
    } else if (*p < 0 || *p > *m) {
#line 353 "sorbdb.f"
	*info = -4;
#line 354 "sorbdb.f"
    } else if (*q < 0 || *q > *p || *q > *m - *p || *q > *m - *q) {
#line 356 "sorbdb.f"
	*info = -5;
#line 357 "sorbdb.f"
    } else if (colmajor && *ldx11 < max(1,*p)) {
#line 358 "sorbdb.f"
	*info = -7;
#line 359 "sorbdb.f"
    } else if (! colmajor && *ldx11 < max(1,*q)) {
#line 360 "sorbdb.f"
	*info = -7;
#line 361 "sorbdb.f"
    } else if (colmajor && *ldx12 < max(1,*p)) {
#line 362 "sorbdb.f"
	*info = -9;
#line 363 "sorbdb.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 363 "sorbdb.f"
	i__1 = 1, i__2 = *m - *q;
#line 363 "sorbdb.f"
	if (! colmajor && *ldx12 < max(i__1,i__2)) {
#line 364 "sorbdb.f"
	    *info = -9;
#line 365 "sorbdb.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 365 "sorbdb.f"
	    i__1 = 1, i__2 = *m - *p;
#line 365 "sorbdb.f"
	    if (colmajor && *ldx21 < max(i__1,i__2)) {
#line 366 "sorbdb.f"
		*info = -11;
#line 367 "sorbdb.f"
	    } else if (! colmajor && *ldx21 < max(1,*q)) {
#line 368 "sorbdb.f"
		*info = -11;
#line 369 "sorbdb.f"
	    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 369 "sorbdb.f"
		i__1 = 1, i__2 = *m - *p;
#line 369 "sorbdb.f"
		if (colmajor && *ldx22 < max(i__1,i__2)) {
#line 370 "sorbdb.f"
		    *info = -13;
#line 371 "sorbdb.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 371 "sorbdb.f"
		    i__1 = 1, i__2 = *m - *q;
#line 371 "sorbdb.f"
		    if (! colmajor && *ldx22 < max(i__1,i__2)) {
#line 372 "sorbdb.f"
			*info = -13;
#line 373 "sorbdb.f"
		    }
#line 373 "sorbdb.f"
		}
#line 373 "sorbdb.f"
	    }
#line 373 "sorbdb.f"
	}
#line 373 "sorbdb.f"
    }

/*     Compute workspace */

#line 377 "sorbdb.f"
    if (*info == 0) {
#line 378 "sorbdb.f"
	lworkopt = *m - *q;
#line 379 "sorbdb.f"
	lworkmin = *m - *q;
#line 380 "sorbdb.f"
	work[1] = (doublereal) lworkopt;
#line 381 "sorbdb.f"
	if (*lwork < lworkmin && ! lquery) {
#line 382 "sorbdb.f"
	    *info = -21;
#line 383 "sorbdb.f"
	}
#line 384 "sorbdb.f"
    }
#line 385 "sorbdb.f"
    if (*info != 0) {
#line 386 "sorbdb.f"
	i__1 = -(*info);
#line 386 "sorbdb.f"
	xerbla_("xORBDB", &i__1, (ftnlen)6);
#line 387 "sorbdb.f"
	return 0;
#line 388 "sorbdb.f"
    } else if (lquery) {
#line 389 "sorbdb.f"
	return 0;
#line 390 "sorbdb.f"
    }

/*     Handle column-major and row-major separately */

#line 394 "sorbdb.f"
    if (colmajor) {

/*        Reduce columns 1, ..., Q of X11, X12, X21, and X22 */

#line 398 "sorbdb.f"
	i__1 = *q;
#line 398 "sorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 400 "sorbdb.f"
	    if (i__ == 1) {
#line 401 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 401 "sorbdb.f"
		sscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 402 "sorbdb.f"
	    } else {
#line 403 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 403 "sorbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 403 "sorbdb.f"
		sscal_(&i__2, &d__1, &x11[i__ + i__ * x11_dim1], &c__1);
#line 404 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 404 "sorbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 404 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x12[i__ + (i__ - 1) * x12_dim1], &c__1, 
			&x11[i__ + i__ * x11_dim1], &c__1);
#line 406 "sorbdb.f"
	    }
#line 407 "sorbdb.f"
	    if (i__ == 1) {
#line 408 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 408 "sorbdb.f"
		sscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], &c__1);
#line 409 "sorbdb.f"
	    } else {
#line 410 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 410 "sorbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 410 "sorbdb.f"
		sscal_(&i__2, &d__1, &x21[i__ + i__ * x21_dim1], &c__1);
#line 411 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 411 "sorbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 411 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x22[i__ + (i__ - 1) * x22_dim1], &c__1, 
			&x21[i__ + i__ * x21_dim1], &c__1);
#line 413 "sorbdb.f"
	    }

#line 415 "sorbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 415 "sorbdb.f"
	    i__3 = *p - i__ + 1;
#line 415 "sorbdb.f"
	    theta[i__] = atan2(snrm2_(&i__2, &x21[i__ + i__ * x21_dim1], &
		    c__1), snrm2_(&i__3, &x11[i__ + i__ * x11_dim1], &c__1));

#line 418 "sorbdb.f"
	    if (*p > i__) {
#line 419 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 419 "sorbdb.f"
		slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + 
			i__ * x11_dim1], &c__1, &taup1[i__]);
#line 420 "sorbdb.f"
	    } else if (*p == i__) {
#line 421 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 421 "sorbdb.f"
		slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + i__ * 
			x11_dim1], &c__1, &taup1[i__]);
#line 422 "sorbdb.f"
	    }
#line 423 "sorbdb.f"
	    x11[i__ + i__ * x11_dim1] = 1.;
#line 424 "sorbdb.f"
	    if (*m - *p > i__) {
#line 425 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 425 "sorbdb.f"
		slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + 
			i__ * x21_dim1], &c__1, &taup2[i__]);
#line 427 "sorbdb.f"
	    } else if (*m - *p == i__) {
#line 428 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 428 "sorbdb.f"
		slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], &c__1, &taup2[i__]);
#line 429 "sorbdb.f"
	    }
#line 430 "sorbdb.f"
	    x21[i__ + i__ * x21_dim1] = 1.;

#line 432 "sorbdb.f"
	    if (*q > i__) {
#line 433 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 433 "sorbdb.f"
		i__3 = *q - i__;
#line 433 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			taup1[i__], &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &
			work[1], (ftnlen)1);
#line 435 "sorbdb.f"
	    }
#line 436 "sorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 437 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 437 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 437 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &
			taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 439 "sorbdb.f"
	    }
#line 440 "sorbdb.f"
	    if (*q > i__) {
#line 441 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 441 "sorbdb.f"
		i__3 = *q - i__;
#line 441 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			taup2[i__], &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &
			work[1], (ftnlen)1);
#line 443 "sorbdb.f"
	    }
#line 444 "sorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 445 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 445 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 445 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &
			taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 447 "sorbdb.f"
	    }

#line 449 "sorbdb.f"
	    if (i__ < *q) {
#line 450 "sorbdb.f"
		i__2 = *q - i__;
#line 450 "sorbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 450 "sorbdb.f"
		sscal_(&i__2, &d__1, &x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 452 "sorbdb.f"
		i__2 = *q - i__;
#line 452 "sorbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 452 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, 
			&x11[i__ + (i__ + 1) * x11_dim1], ldx11);
#line 454 "sorbdb.f"
	    }
#line 455 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 455 "sorbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 455 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 456 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 456 "sorbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 456 "sorbdb.f"
	    saxpy_(&i__2, &d__1, &x22[i__ + i__ * x22_dim1], ldx22, &x12[i__ 
		    + i__ * x12_dim1], ldx12);

#line 459 "sorbdb.f"
	    if (i__ < *q) {
#line 459 "sorbdb.f"
		i__2 = *q - i__;
#line 459 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 459 "sorbdb.f"
		phi[i__] = atan2(snrm2_(&i__2, &x11[i__ + (i__ + 1) * 
			x11_dim1], ldx11), snrm2_(&i__3, &x12[i__ + i__ * 
			x12_dim1], ldx12));
#line 459 "sorbdb.f"
	    }

#line 463 "sorbdb.f"
	    if (i__ < *q) {
#line 464 "sorbdb.f"
		if (*q - i__ == 1) {
#line 465 "sorbdb.f"
		    i__2 = *q - i__;
#line 465 "sorbdb.f"
		    slarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 1) * x11_dim1], ldx11, &tauq1[i__]);
#line 467 "sorbdb.f"
		} else {
#line 468 "sorbdb.f"
		    i__2 = *q - i__;
#line 468 "sorbdb.f"
		    slarfgp_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], &x11[
			    i__ + (i__ + 2) * x11_dim1], ldx11, &tauq1[i__]);
#line 470 "sorbdb.f"
		}
#line 471 "sorbdb.f"
		x11[i__ + (i__ + 1) * x11_dim1] = 1.;
#line 472 "sorbdb.f"
	    }
#line 473 "sorbdb.f"
	    if (*q + i__ - 1 < *m) {
#line 474 "sorbdb.f"
		if (*m - *q == i__) {
#line 475 "sorbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 475 "sorbdb.f"
		    slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 
			    i__ * x12_dim1], ldx12, &tauq2[i__]);
#line 477 "sorbdb.f"
		} else {
#line 478 "sorbdb.f"
		    i__2 = *m - *q - i__ + 1;
#line 478 "sorbdb.f"
		    slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (
			    i__ + 1) * x12_dim1], ldx12, &tauq2[i__]);
#line 480 "sorbdb.f"
		}
#line 481 "sorbdb.f"
	    }
#line 482 "sorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 484 "sorbdb.f"
	    if (i__ < *q) {
#line 485 "sorbdb.f"
		i__2 = *p - i__;
#line 485 "sorbdb.f"
		i__3 = *q - i__;
#line 485 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 487 "sorbdb.f"
		i__2 = *m - *p - i__;
#line 487 "sorbdb.f"
		i__3 = *q - i__;
#line 487 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x11[i__ + (i__ + 1) * x11_dim1], 
			ldx11, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 489 "sorbdb.f"
	    }
#line 490 "sorbdb.f"
	    if (*p > i__) {
#line 491 "sorbdb.f"
		i__2 = *p - i__;
#line 491 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 491 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 493 "sorbdb.f"
	    }
#line 494 "sorbdb.f"
	    if (*m - *p > i__) {
#line 495 "sorbdb.f"
		i__2 = *m - *p - i__;
#line 495 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 495 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[i__ + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 497 "sorbdb.f"
	    }

#line 499 "sorbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 503 "sorbdb.f"
	i__1 = *p;
#line 503 "sorbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 505 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 505 "sorbdb.f"
	    d__1 = -z1 * z4;
#line 505 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], ldx12);
#line 506 "sorbdb.f"
	    if (i__ >= *m - *q) {
#line 507 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 507 "sorbdb.f"
		slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], ldx12, &tauq2[i__]);
#line 509 "sorbdb.f"
	    } else {
#line 510 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 510 "sorbdb.f"
		slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + (i__ + 
			1) * x12_dim1], ldx12, &tauq2[i__]);
#line 512 "sorbdb.f"
	    }
#line 513 "sorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 515 "sorbdb.f"
	    if (*p > i__) {
#line 516 "sorbdb.f"
		i__2 = *p - i__;
#line 516 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 516 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x12[i__ + 1 + i__ * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 518 "sorbdb.f"
	    }
#line 519 "sorbdb.f"
	    if (*m - *p - *q >= 1) {
#line 519 "sorbdb.f"
		i__2 = *m - *p - *q;
#line 519 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 519 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], ldx12, &
			tauq2[i__], &x22[*q + 1 + i__ * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 519 "sorbdb.f"
	    }

#line 523 "sorbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 527 "sorbdb.f"
	i__1 = *m - *p - *q;
#line 527 "sorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 529 "sorbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 529 "sorbdb.f"
	    d__1 = z2 * z4;
#line 529 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x22[*q + i__ + (*p + i__) * x22_dim1], 
		    ldx22);
#line 530 "sorbdb.f"
	    if (i__ == *m - *p - *q) {
#line 531 "sorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 531 "sorbdb.f"
		slarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*
			q + i__ + (*p + i__) * x22_dim1], ldx22, &tauq2[*p + 
			i__]);
#line 533 "sorbdb.f"
	    } else {
#line 534 "sorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 534 "sorbdb.f"
		slarfgp_(&i__2, &x22[*q + i__ + (*p + i__) * x22_dim1], &x22[*
			q + i__ + (*p + i__ + 1) * x22_dim1], ldx22, &tauq2[*
			p + i__]);
#line 536 "sorbdb.f"
	    }
#line 537 "sorbdb.f"
	    x22[*q + i__ + (*p + i__) * x22_dim1] = 1.;
#line 538 "sorbdb.f"
	    if (i__ < *m - *p - *q) {
#line 539 "sorbdb.f"
		i__2 = *m - *p - *q - i__;
#line 539 "sorbdb.f"
		i__3 = *m - *p - *q - i__ + 1;
#line 539 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x22[*q + i__ + (*p + i__) * 
			x22_dim1], ldx22, &tauq2[*p + i__], &x22[*q + i__ + 1 
			+ (*p + i__) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 541 "sorbdb.f"
	    }

#line 543 "sorbdb.f"
	}

#line 545 "sorbdb.f"
    } else {

/*        Reduce columns 1, ..., Q of X11, X12, X21, X22 */

#line 549 "sorbdb.f"
	i__1 = *q;
#line 549 "sorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 551 "sorbdb.f"
	    if (i__ == 1) {
#line 552 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 552 "sorbdb.f"
		sscal_(&i__2, &z1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 553 "sorbdb.f"
	    } else {
#line 554 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 554 "sorbdb.f"
		d__1 = z1 * cos(phi[i__ - 1]);
#line 554 "sorbdb.f"
		sscal_(&i__2, &d__1, &x11[i__ + i__ * x11_dim1], ldx11);
#line 555 "sorbdb.f"
		i__2 = *p - i__ + 1;
#line 555 "sorbdb.f"
		d__1 = -z1 * z3 * z4 * sin(phi[i__ - 1]);
#line 555 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x12[i__ - 1 + i__ * x12_dim1], ldx12, &
			x11[i__ + i__ * x11_dim1], ldx11);
#line 557 "sorbdb.f"
	    }
#line 558 "sorbdb.f"
	    if (i__ == 1) {
#line 559 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 559 "sorbdb.f"
		sscal_(&i__2, &z2, &x21[i__ + i__ * x21_dim1], ldx21);
#line 560 "sorbdb.f"
	    } else {
#line 561 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 561 "sorbdb.f"
		d__1 = z2 * cos(phi[i__ - 1]);
#line 561 "sorbdb.f"
		sscal_(&i__2, &d__1, &x21[i__ + i__ * x21_dim1], ldx21);
#line 562 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 562 "sorbdb.f"
		d__1 = -z2 * z3 * z4 * sin(phi[i__ - 1]);
#line 562 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x22[i__ - 1 + i__ * x22_dim1], ldx22, &
			x21[i__ + i__ * x21_dim1], ldx21);
#line 564 "sorbdb.f"
	    }

#line 566 "sorbdb.f"
	    i__2 = *m - *p - i__ + 1;
#line 566 "sorbdb.f"
	    i__3 = *p - i__ + 1;
#line 566 "sorbdb.f"
	    theta[i__] = atan2(snrm2_(&i__2, &x21[i__ + i__ * x21_dim1], 
		    ldx21), snrm2_(&i__3, &x11[i__ + i__ * x11_dim1], ldx11));

#line 569 "sorbdb.f"
	    i__2 = *p - i__ + 1;
#line 569 "sorbdb.f"
	    slarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) *
		     x11_dim1], ldx11, &taup1[i__]);
#line 570 "sorbdb.f"
	    x11[i__ + i__ * x11_dim1] = 1.;
#line 571 "sorbdb.f"
	    if (i__ == *m - *p) {
#line 572 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 572 "sorbdb.f"
		slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + i__ * 
			x21_dim1], ldx21, &taup2[i__]);
#line 574 "sorbdb.f"
	    } else {
#line 575 "sorbdb.f"
		i__2 = *m - *p - i__ + 1;
#line 575 "sorbdb.f"
		slarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 
			1) * x21_dim1], ldx21, &taup2[i__]);
#line 577 "sorbdb.f"
	    }
#line 578 "sorbdb.f"
	    x21[i__ + i__ * x21_dim1] = 1.;

#line 580 "sorbdb.f"
	    if (*q > i__) {
#line 581 "sorbdb.f"
		i__2 = *q - i__;
#line 581 "sorbdb.f"
		i__3 = *p - i__ + 1;
#line 581 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
			taup1[i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &
			work[1], (ftnlen)1);
#line 583 "sorbdb.f"
	    }
#line 584 "sorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 585 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 585 "sorbdb.f"
		i__3 = *p - i__ + 1;
#line 585 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &
			taup1[i__], &x12[i__ + i__ * x12_dim1], ldx12, &work[
			1], (ftnlen)1);
#line 587 "sorbdb.f"
	    }
#line 588 "sorbdb.f"
	    if (*q > i__) {
#line 589 "sorbdb.f"
		i__2 = *q - i__;
#line 589 "sorbdb.f"
		i__3 = *m - *p - i__ + 1;
#line 589 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
			taup2[i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &
			work[1], (ftnlen)1);
#line 591 "sorbdb.f"
	    }
#line 592 "sorbdb.f"
	    if (*m - *q + 1 > i__) {
#line 593 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 593 "sorbdb.f"
		i__3 = *m - *p - i__ + 1;
#line 593 "sorbdb.f"
		slarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &
			taup2[i__], &x22[i__ + i__ * x22_dim1], ldx22, &work[
			1], (ftnlen)1);
#line 595 "sorbdb.f"
	    }

#line 597 "sorbdb.f"
	    if (i__ < *q) {
#line 598 "sorbdb.f"
		i__2 = *q - i__;
#line 598 "sorbdb.f"
		d__1 = -z1 * z3 * sin(theta[i__]);
#line 598 "sorbdb.f"
		sscal_(&i__2, &d__1, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 599 "sorbdb.f"
		i__2 = *q - i__;
#line 599 "sorbdb.f"
		d__1 = z2 * z3 * cos(theta[i__]);
#line 599 "sorbdb.f"
		saxpy_(&i__2, &d__1, &x21[i__ + 1 + i__ * x21_dim1], &c__1, &
			x11[i__ + 1 + i__ * x11_dim1], &c__1);
#line 601 "sorbdb.f"
	    }
#line 602 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 602 "sorbdb.f"
	    d__1 = -z1 * z4 * sin(theta[i__]);
#line 602 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 603 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 603 "sorbdb.f"
	    d__1 = z2 * z4 * cos(theta[i__]);
#line 603 "sorbdb.f"
	    saxpy_(&i__2, &d__1, &x22[i__ + i__ * x22_dim1], &c__1, &x12[i__ 
		    + i__ * x12_dim1], &c__1);

#line 606 "sorbdb.f"
	    if (i__ < *q) {
#line 606 "sorbdb.f"
		i__2 = *q - i__;
#line 606 "sorbdb.f"
		i__3 = *m - *q - i__ + 1;
#line 606 "sorbdb.f"
		phi[i__] = atan2(snrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1],
			 &c__1), snrm2_(&i__3, &x12[i__ + i__ * x12_dim1], &
			c__1));
#line 606 "sorbdb.f"
	    }

#line 610 "sorbdb.f"
	    if (i__ < *q) {
#line 611 "sorbdb.f"
		if (*q - i__ == 1) {
#line 612 "sorbdb.f"
		    i__2 = *q - i__;
#line 612 "sorbdb.f"
		    slarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ 
			    + 1 + i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 614 "sorbdb.f"
		} else {
#line 615 "sorbdb.f"
		    i__2 = *q - i__;
#line 615 "sorbdb.f"
		    slarfgp_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &x11[i__ 
			    + 2 + i__ * x11_dim1], &c__1, &tauq1[i__]);
#line 617 "sorbdb.f"
		}
#line 618 "sorbdb.f"
		x11[i__ + 1 + i__ * x11_dim1] = 1.;
#line 619 "sorbdb.f"
	    }
#line 620 "sorbdb.f"
	    if (*m - *q > i__) {
#line 621 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 621 "sorbdb.f"
		slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + 
			i__ * x12_dim1], &c__1, &tauq2[i__]);
#line 623 "sorbdb.f"
	    } else {
#line 624 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 624 "sorbdb.f"
		slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + i__ * 
			x12_dim1], &c__1, &tauq2[i__]);
#line 626 "sorbdb.f"
	    }
#line 627 "sorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 629 "sorbdb.f"
	    if (i__ < *q) {
#line 630 "sorbdb.f"
		i__2 = *q - i__;
#line 630 "sorbdb.f"
		i__3 = *p - i__;
#line 630 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * 
			x11_dim1], ldx11, &work[1], (ftnlen)1);
#line 632 "sorbdb.f"
		i__2 = *q - i__;
#line 632 "sorbdb.f"
		i__3 = *m - *p - i__;
#line 632 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x11[i__ + 1 + i__ * x11_dim1], &
			c__1, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * 
			x21_dim1], ldx21, &work[1], (ftnlen)1);
#line 634 "sorbdb.f"
	    }
#line 635 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 635 "sorbdb.f"
	    i__3 = *p - i__;
#line 635 "sorbdb.f"
	    slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
		    tauq2[i__], &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &
		    work[1], (ftnlen)1);
#line 637 "sorbdb.f"
	    if (*m - *p - i__ > 0) {
#line 638 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 638 "sorbdb.f"
		i__3 = *m - *p - i__;
#line 638 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x22[i__ + (i__ + 1) * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 640 "sorbdb.f"
	    }

#line 642 "sorbdb.f"
	}

/*        Reduce columns Q + 1, ..., P of X12, X22 */

#line 646 "sorbdb.f"
	i__1 = *p;
#line 646 "sorbdb.f"
	for (i__ = *q + 1; i__ <= i__1; ++i__) {

#line 648 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 648 "sorbdb.f"
	    d__1 = -z1 * z4;
#line 648 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x12[i__ + i__ * x12_dim1], &c__1);
#line 649 "sorbdb.f"
	    i__2 = *m - *q - i__ + 1;
#line 649 "sorbdb.f"
	    slarfgp_(&i__2, &x12[i__ + i__ * x12_dim1], &x12[i__ + 1 + i__ * 
		    x12_dim1], &c__1, &tauq2[i__]);
#line 650 "sorbdb.f"
	    x12[i__ + i__ * x12_dim1] = 1.;

#line 652 "sorbdb.f"
	    if (*p > i__) {
#line 653 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 653 "sorbdb.f"
		i__3 = *p - i__;
#line 653 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x12[i__ + (i__ + 1) * x12_dim1], ldx12, &
			work[1], (ftnlen)1);
#line 655 "sorbdb.f"
	    }
#line 656 "sorbdb.f"
	    if (*m - *p - *q >= 1) {
#line 656 "sorbdb.f"
		i__2 = *m - *q - i__ + 1;
#line 656 "sorbdb.f"
		i__3 = *m - *p - *q;
#line 656 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x12[i__ + i__ * x12_dim1], &c__1, &
			tauq2[i__], &x22[i__ + (*q + 1) * x22_dim1], ldx22, &
			work[1], (ftnlen)1);
#line 656 "sorbdb.f"
	    }

#line 660 "sorbdb.f"
	}

/*        Reduce columns P + 1, ..., M - Q of X12, X22 */

#line 664 "sorbdb.f"
	i__1 = *m - *p - *q;
#line 664 "sorbdb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 666 "sorbdb.f"
	    i__2 = *m - *p - *q - i__ + 1;
#line 666 "sorbdb.f"
	    d__1 = z2 * z4;
#line 666 "sorbdb.f"
	    sscal_(&i__2, &d__1, &x22[*p + i__ + (*q + i__) * x22_dim1], &
		    c__1);
#line 667 "sorbdb.f"
	    if (*m - *p - *q == i__) {
#line 668 "sorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 668 "sorbdb.f"
		slarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*
			p + i__ + (*q + i__) * x22_dim1], &c__1, &tauq2[*p + 
			i__]);
#line 670 "sorbdb.f"
		x22[*p + i__ + (*q + i__) * x22_dim1] = 1.;
#line 671 "sorbdb.f"
	    } else {
#line 672 "sorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 672 "sorbdb.f"
		slarfgp_(&i__2, &x22[*p + i__ + (*q + i__) * x22_dim1], &x22[*
			p + i__ + 1 + (*q + i__) * x22_dim1], &c__1, &tauq2[*
			p + i__]);
#line 674 "sorbdb.f"
		x22[*p + i__ + (*q + i__) * x22_dim1] = 1.;
#line 675 "sorbdb.f"
		i__2 = *m - *p - *q - i__ + 1;
#line 675 "sorbdb.f"
		i__3 = *m - *p - *q - i__;
#line 675 "sorbdb.f"
		slarf_("L", &i__2, &i__3, &x22[*p + i__ + (*q + i__) * 
			x22_dim1], &c__1, &tauq2[*p + i__], &x22[*p + i__ + (*
			q + i__ + 1) * x22_dim1], ldx22, &work[1], (ftnlen)1);
#line 677 "sorbdb.f"
	    }


#line 680 "sorbdb.f"
	}

#line 682 "sorbdb.f"
    }

#line 684 "sorbdb.f"
    return 0;

/*     End of SORBDB */

} /* sorbdb_ */

