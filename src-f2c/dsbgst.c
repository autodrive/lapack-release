#line 1 "dsbgst.f"
/* dsbgst.f -- translated by f2c (version 20100827).
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

#line 1 "dsbgst.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b20 = -1.;

/* > \brief \b DSBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, */
/*                          LDX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBGST reduces a real symmetric-definite banded generalized */
/* > eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y, */
/* > such that C has the same bandwidth as A. */
/* > */
/* > B must have been previously factorized as S**T*S by DPBSTF, using a */
/* > split Cholesky factorization. A is overwritten by C = X**T*A*X, where */
/* > X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the */
/* > bandwidth of A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          = 'N':  do not form the transformation matrix X; */
/* >          = 'V':  form X. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KA */
/* > \verbatim */
/* >          KA is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KA >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of superdiagonals of the matrix B if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the transformed matrix X**T*A*X, stored in the same */
/* >          format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KA+1. */
/* > \endverbatim */
/* > */
/* > \param[in] BB */
/* > \verbatim */
/* >          BB is DOUBLE PRECISION array, dimension (LDBB,N) */
/* >          The banded factor S from the split Cholesky factorization of */
/* >          B, as returned by DPBSTF, stored in the first KB+1 rows of */
/* >          the array. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBB */
/* > \verbatim */
/* >          LDBB is INTEGER */
/* >          The leading dimension of the array BB.  LDBB >= KB+1. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,N) */
/* >          If VECT = 'V', the n-by-n matrix X. */
/* >          If VECT = 'N', the array X is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X. */
/* >          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
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

/*  ===================================================================== */
/* Subroutine */ int dsbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info, 
	ftnlen vect_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal t;
    static integer i0, i1, i2, j1, j2;
    static doublereal ra;
    static integer nr, nx, ka1, kb1;
    static doublereal ra1;
    static integer j1t, j2t;
    static doublereal bii;
    static integer kbt, nrt, inca;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), drot_(integer *, doublereal *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper, wantx;
    extern /* Subroutine */ int dlar2v_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen), dlargv_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static logical update;
    extern /* Subroutine */ int dlartv_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);


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

/*     Test the input parameters */

#line 203 "dsbgst.f"
    /* Parameter adjustments */
#line 203 "dsbgst.f"
    ab_dim1 = *ldab;
#line 203 "dsbgst.f"
    ab_offset = 1 + ab_dim1;
#line 203 "dsbgst.f"
    ab -= ab_offset;
#line 203 "dsbgst.f"
    bb_dim1 = *ldbb;
#line 203 "dsbgst.f"
    bb_offset = 1 + bb_dim1;
#line 203 "dsbgst.f"
    bb -= bb_offset;
#line 203 "dsbgst.f"
    x_dim1 = *ldx;
#line 203 "dsbgst.f"
    x_offset = 1 + x_dim1;
#line 203 "dsbgst.f"
    x -= x_offset;
#line 203 "dsbgst.f"
    --work;
#line 203 "dsbgst.f"

#line 203 "dsbgst.f"
    /* Function Body */
#line 203 "dsbgst.f"
    wantx = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 204 "dsbgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 205 "dsbgst.f"
    ka1 = *ka + 1;
#line 206 "dsbgst.f"
    kb1 = *kb + 1;
#line 207 "dsbgst.f"
    *info = 0;
#line 208 "dsbgst.f"
    if (! wantx && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 209 "dsbgst.f"
	*info = -1;
#line 210 "dsbgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 211 "dsbgst.f"
	*info = -2;
#line 212 "dsbgst.f"
    } else if (*n < 0) {
#line 213 "dsbgst.f"
	*info = -3;
#line 214 "dsbgst.f"
    } else if (*ka < 0) {
#line 215 "dsbgst.f"
	*info = -4;
#line 216 "dsbgst.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 217 "dsbgst.f"
	*info = -5;
#line 218 "dsbgst.f"
    } else if (*ldab < *ka + 1) {
#line 219 "dsbgst.f"
	*info = -7;
#line 220 "dsbgst.f"
    } else if (*ldbb < *kb + 1) {
#line 221 "dsbgst.f"
	*info = -9;
#line 222 "dsbgst.f"
    } else if (*ldx < 1 || wantx && *ldx < max(1,*n)) {
#line 223 "dsbgst.f"
	*info = -11;
#line 224 "dsbgst.f"
    }
#line 225 "dsbgst.f"
    if (*info != 0) {
#line 226 "dsbgst.f"
	i__1 = -(*info);
#line 226 "dsbgst.f"
	xerbla_("DSBGST", &i__1, (ftnlen)6);
#line 227 "dsbgst.f"
	return 0;
#line 228 "dsbgst.f"
    }

/*     Quick return if possible */

#line 232 "dsbgst.f"
    if (*n == 0) {
#line 232 "dsbgst.f"
	return 0;
#line 232 "dsbgst.f"
    }

#line 235 "dsbgst.f"
    inca = *ldab * ka1;

/*     Initialize X to the unit matrix, if needed */

#line 239 "dsbgst.f"
    if (wantx) {
#line 239 "dsbgst.f"
	dlaset_("Full", n, n, &c_b8, &c_b9, &x[x_offset], ldx, (ftnlen)4);
#line 239 "dsbgst.f"
    }

/*     Set M to the splitting point m. It must be the same value as is */
/*     used in DPBSTF. The chosen value allows the arrays WORK and RWORK */
/*     to be of dimension (N). */

#line 246 "dsbgst.f"
    m = (*n + *kb) / 2;

/*     The routine works in two phases, corresponding to the two halves */
/*     of the split Cholesky factorization of B as S**T*S where */

/*     S = ( U    ) */
/*         ( M  L ) */

/*     with U upper triangular of order m, and L lower triangular of */
/*     order n-m. S has the same bandwidth as B. */

/*     S is treated as a product of elementary matrices: */

/*     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n) */

/*     where S(i) is determined by the i-th row of S. */

/*     In phase 1, the index i takes the values n, n-1, ... , m+1; */
/*     in phase 2, it takes the values 1, 2, ... , m. */

/*     For each value of i, the current matrix A is updated by forming */
/*     inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside */
/*     the band of A. The bulge is then pushed down toward the bottom of */
/*     A in phase 1, and up toward the top of A in phase 2, by applying */
/*     plane rotations. */

/*     There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1 */
/*     of them are linearly independent, so annihilating a bulge requires */
/*     only 2*kb-1 plane rotations. The rotations are divided into a 1st */
/*     set of kb-1 rotations, and a 2nd set of kb rotations. */

/*     Wherever possible, rotations are generated and applied in vector */
/*     operations of length NR between the indices J1 and J2 (sometimes */
/*     replaced by modified values NRT, J1T or J2T). */

/*     The cosines and sines of the rotations are stored in the array */
/*     WORK. The cosines of the 1st set of rotations are stored in */
/*     elements n+2:n+m-kb-1 and the sines of the 1st set in elements */
/*     2:m-kb-1; the cosines of the 2nd set are stored in elements */
/*     n+m-kb+1:2*n and the sines of the second set in elements m-kb+1:n. */

/*     The bulges are not formed explicitly; nonzero elements outside the */
/*     band are created only when they are required for generating new */
/*     rotations; they are stored in the array WORK, in positions where */
/*     they are later overwritten by the sines of the rotations which */
/*     annihilate them. */

/*     **************************** Phase 1 ***************************** */

/*     The logical structure of this phase is: */

/*     UPDATE = .TRUE. */
/*     DO I = N, M + 1, -1 */
/*        use S(i) to update A and create a new bulge */
/*        apply rotations to push all bulges KA positions downward */
/*     END DO */
/*     UPDATE = .FALSE. */
/*     DO I = M + KA + 1, N - 1 */
/*        apply rotations to push all bulges KA positions downward */
/*     END DO */

/*     To avoid duplicating code, the two loops are merged. */

#line 309 "dsbgst.f"
    update = TRUE_;
#line 310 "dsbgst.f"
    i__ = *n + 1;
#line 311 "dsbgst.f"
L10:
#line 312 "dsbgst.f"
    if (update) {
#line 313 "dsbgst.f"
	--i__;
/* Computing MIN */
#line 314 "dsbgst.f"
	i__1 = *kb, i__2 = i__ - 1;
#line 314 "dsbgst.f"
	kbt = min(i__1,i__2);
#line 315 "dsbgst.f"
	i0 = i__ - 1;
/* Computing MIN */
#line 316 "dsbgst.f"
	i__1 = *n, i__2 = i__ + *ka;
#line 316 "dsbgst.f"
	i1 = min(i__1,i__2);
#line 317 "dsbgst.f"
	i2 = i__ - kbt + ka1;
#line 318 "dsbgst.f"
	if (i__ < m + 1) {
#line 319 "dsbgst.f"
	    update = FALSE_;
#line 320 "dsbgst.f"
	    ++i__;
#line 321 "dsbgst.f"
	    i0 = m;
#line 322 "dsbgst.f"
	    if (*ka == 0) {
#line 322 "dsbgst.f"
		goto L480;
#line 322 "dsbgst.f"
	    }
#line 324 "dsbgst.f"
	    goto L10;
#line 325 "dsbgst.f"
	}
#line 326 "dsbgst.f"
    } else {
#line 327 "dsbgst.f"
	i__ += *ka;
#line 328 "dsbgst.f"
	if (i__ > *n - 1) {
#line 328 "dsbgst.f"
	    goto L480;
#line 328 "dsbgst.f"
	}
#line 330 "dsbgst.f"
    }

#line 332 "dsbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 336 "dsbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 340 "dsbgst.f"
	    bii = bb[kb1 + i__ * bb_dim1];
#line 341 "dsbgst.f"
	    i__1 = i1;
#line 341 "dsbgst.f"
	    for (j = i__; j <= i__1; ++j) {
#line 342 "dsbgst.f"
		ab[i__ - j + ka1 + j * ab_dim1] /= bii;
#line 343 "dsbgst.f"
/* L20: */
#line 343 "dsbgst.f"
	    }
/* Computing MAX */
#line 344 "dsbgst.f"
	    i__1 = 1, i__2 = i__ - *ka;
#line 344 "dsbgst.f"
	    i__3 = i__;
#line 344 "dsbgst.f"
	    for (j = max(i__1,i__2); j <= i__3; ++j) {
#line 345 "dsbgst.f"
		ab[j - i__ + ka1 + i__ * ab_dim1] /= bii;
#line 346 "dsbgst.f"
/* L30: */
#line 346 "dsbgst.f"
	    }
#line 347 "dsbgst.f"
	    i__3 = i__ - 1;
#line 347 "dsbgst.f"
	    for (k = i__ - kbt; k <= i__3; ++k) {
#line 348 "dsbgst.f"
		i__1 = k;
#line 348 "dsbgst.f"
		for (j = i__ - kbt; j <= i__1; ++j) {
#line 349 "dsbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] = ab[j - k + ka1 + k * 
			    ab_dim1] - bb[j - i__ + kb1 + i__ * bb_dim1] * ab[
			    k - i__ + ka1 + i__ * ab_dim1] - bb[k - i__ + kb1 
			    + i__ * bb_dim1] * ab[j - i__ + ka1 + i__ * 
			    ab_dim1] + ab[ka1 + i__ * ab_dim1] * bb[j - i__ + 
			    kb1 + i__ * bb_dim1] * bb[k - i__ + kb1 + i__ * 
			    bb_dim1];
#line 354 "dsbgst.f"
/* L40: */
#line 354 "dsbgst.f"
		}
/* Computing MAX */
#line 355 "dsbgst.f"
		i__1 = 1, i__2 = i__ - *ka;
#line 355 "dsbgst.f"
		i__4 = i__ - kbt - 1;
#line 355 "dsbgst.f"
		for (j = max(i__1,i__2); j <= i__4; ++j) {
#line 356 "dsbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] -= bb[k - i__ + kb1 + i__ * 
			    bb_dim1] * ab[j - i__ + ka1 + i__ * ab_dim1];
#line 358 "dsbgst.f"
/* L50: */
#line 358 "dsbgst.f"
		}
#line 359 "dsbgst.f"
/* L60: */
#line 359 "dsbgst.f"
	    }
#line 360 "dsbgst.f"
	    i__3 = i1;
#line 360 "dsbgst.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MAX */
#line 361 "dsbgst.f"
		i__4 = j - *ka, i__1 = i__ - kbt;
#line 361 "dsbgst.f"
		i__2 = i__ - 1;
#line 361 "dsbgst.f"
		for (k = max(i__4,i__1); k <= i__2; ++k) {
#line 362 "dsbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] -= bb[k - i__ + kb1 + i__ * 
			    bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1];
#line 364 "dsbgst.f"
/* L70: */
#line 364 "dsbgst.f"
		}
#line 365 "dsbgst.f"
/* L80: */
#line 365 "dsbgst.f"
	    }

#line 367 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 371 "dsbgst.f"
		i__3 = *n - m;
#line 371 "dsbgst.f"
		d__1 = 1. / bii;
#line 371 "dsbgst.f"
		dscal_(&i__3, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 372 "dsbgst.f"
		if (kbt > 0) {
#line 372 "dsbgst.f"
		    i__3 = *n - m;
#line 372 "dsbgst.f"
		    dger_(&i__3, &kbt, &c_b20, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kb1 - kbt + i__ * bb_dim1], &c__1, &x[m 
			    + 1 + (i__ - kbt) * x_dim1], ldx);
#line 372 "dsbgst.f"
		}
#line 375 "dsbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 379 "dsbgst.f"
	    ra1 = ab[i__ - i1 + ka1 + i1 * ab_dim1];
#line 380 "dsbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 386 "dsbgst.f"
	i__3 = *kb - 1;
#line 386 "dsbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 387 "dsbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 392 "dsbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i,i-k+ka+1) */

#line 396 "dsbgst.f"
		    dlartg_(&ab[k + 1 + (i__ - k + *ka) * ab_dim1], &ra1, &
			    work[*n + i__ - k + *ka - m], &work[i__ - k + *ka 
			    - m], &ra);

/*                 create nonzero element a(i-k,i-k+ka+1) outside the */
/*                 band and store it in WORK(i-k) */

#line 403 "dsbgst.f"
		    t = -bb[kb1 - k + i__ * bb_dim1] * ra1;
#line 404 "dsbgst.f"
		    work[i__ - k] = work[*n + i__ - k + *ka - m] * t - work[
			    i__ - k + *ka - m] * ab[(i__ - k + *ka) * ab_dim1 
			    + 1];
#line 406 "dsbgst.f"
		    ab[(i__ - k + *ka) * ab_dim1 + 1] = work[i__ - k + *ka - 
			    m] * t + work[*n + i__ - k + *ka - m] * ab[(i__ - 
			    k + *ka) * ab_dim1 + 1];
#line 408 "dsbgst.f"
		    ra1 = ra;
#line 409 "dsbgst.f"
		}
#line 410 "dsbgst.f"
	    }
/* Computing MAX */
#line 411 "dsbgst.f"
	    i__2 = 1, i__4 = k - i0 + 2;
#line 411 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__2,i__4) * ka1;
#line 412 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 413 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 414 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 415 "dsbgst.f"
		i__2 = j2, i__4 = i__ + (*ka << 1) - k + 1;
#line 415 "dsbgst.f"
		j2t = max(i__2,i__4);
#line 416 "dsbgst.f"
	    } else {
#line 417 "dsbgst.f"
		j2t = j2;
#line 418 "dsbgst.f"
	    }
#line 419 "dsbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 420 "dsbgst.f"
	    i__2 = j1;
#line 420 "dsbgst.f"
	    i__4 = ka1;
#line 420 "dsbgst.f"
	    for (j = j2t; i__4 < 0 ? j >= i__2 : j <= i__2; j += i__4) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j-m) */

#line 425 "dsbgst.f"
		work[j - m] *= ab[(j + 1) * ab_dim1 + 1];
#line 426 "dsbgst.f"
		ab[(j + 1) * ab_dim1 + 1] = work[*n + j - m] * ab[(j + 1) * 
			ab_dim1 + 1];
#line 427 "dsbgst.f"
/* L90: */
#line 427 "dsbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 432 "dsbgst.f"
	    if (nrt > 0) {
#line 432 "dsbgst.f"
		dlargv_(&nrt, &ab[j2t * ab_dim1 + 1], &inca, &work[j2t - m], &
			ka1, &work[*n + j2t - m], &ka1);
#line 432 "dsbgst.f"
	    }
#line 435 "dsbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 439 "dsbgst.f"
		i__4 = *ka - 1;
#line 439 "dsbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 440 "dsbgst.f"
		    dlartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &work[*n + j2 - 
			    m], &work[j2 - m], &ka1);
#line 443 "dsbgst.f"
/* L100: */
#line 443 "dsbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 448 "dsbgst.f"
		dlar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &work[
			*n + j2 - m], &work[j2 - m], &ka1);

#line 452 "dsbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 456 "dsbgst.f"
	    i__4 = *kb - k + 1;
#line 456 "dsbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 457 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 458 "dsbgst.f"
		if (nrt > 0) {
#line 458 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2 - m], &work[j2 - m], &ka1);
#line 458 "dsbgst.f"
		}
#line 462 "dsbgst.f"
/* L110: */
#line 462 "dsbgst.f"
	    }

#line 464 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 468 "dsbgst.f"
		i__4 = j1;
#line 468 "dsbgst.f"
		i__2 = ka1;
#line 468 "dsbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__4 : j <= i__4; j += i__2) {
#line 469 "dsbgst.f"
		    i__1 = *n - m;
#line 469 "dsbgst.f"
		    drot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j - m], &work[j 
			    - m]);
#line 471 "dsbgst.f"
/* L120: */
#line 471 "dsbgst.f"
		}
#line 472 "dsbgst.f"
	    }
#line 473 "dsbgst.f"
/* L130: */
#line 473 "dsbgst.f"
	}

#line 475 "dsbgst.f"
	if (update) {
#line 476 "dsbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt,i-kbt+ka+1) outside the */
/*              band and store it in WORK(i-kbt) */

#line 481 "dsbgst.f"
		work[i__ - kbt] = -bb[kb1 - kbt + i__ * bb_dim1] * ra1;
#line 482 "dsbgst.f"
	    }
#line 483 "dsbgst.f"
	}

#line 485 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 486 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 487 "dsbgst.f"
		i__3 = 2, i__2 = k - i0 + 1;
#line 487 "dsbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 488 "dsbgst.f"
	    } else {
/* Computing MAX */
#line 489 "dsbgst.f"
		i__3 = 1, i__2 = k - i0 + 1;
#line 489 "dsbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 490 "dsbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 494 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 495 "dsbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 496 "dsbgst.f"
		if (nrt > 0) {
#line 496 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + (j2 - l + 1) * ab_dim1], &inca, &ab[
			    l + 1 + (j2 - l + 1) * ab_dim1], &inca, &work[*n 
			    + j2 - *ka], &work[j2 - *ka], &ka1);
#line 496 "dsbgst.f"
		}
#line 500 "dsbgst.f"
/* L140: */
#line 500 "dsbgst.f"
	    }
#line 501 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 502 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 503 "dsbgst.f"
	    i__3 = j2;
#line 503 "dsbgst.f"
	    i__2 = -ka1;
#line 503 "dsbgst.f"
	    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 504 "dsbgst.f"
		work[j] = work[j - *ka];
#line 505 "dsbgst.f"
		work[*n + j] = work[*n + j - *ka];
#line 506 "dsbgst.f"
/* L150: */
#line 506 "dsbgst.f"
	    }
#line 507 "dsbgst.f"
	    i__2 = j1;
#line 507 "dsbgst.f"
	    i__3 = ka1;
#line 507 "dsbgst.f"
	    for (j = j2; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j) */

#line 512 "dsbgst.f"
		work[j] *= ab[(j + 1) * ab_dim1 + 1];
#line 513 "dsbgst.f"
		ab[(j + 1) * ab_dim1 + 1] = work[*n + j] * ab[(j + 1) * 
			ab_dim1 + 1];
#line 514 "dsbgst.f"
/* L160: */
#line 514 "dsbgst.f"
	    }
#line 515 "dsbgst.f"
	    if (update) {
#line 516 "dsbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 516 "dsbgst.f"
		    work[i__ - k + *ka] = work[i__ - k];
#line 516 "dsbgst.f"
		}
#line 518 "dsbgst.f"
	    }
#line 519 "dsbgst.f"
/* L170: */
#line 519 "dsbgst.f"
	}

#line 521 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 522 "dsbgst.f"
	    i__3 = 1, i__2 = k - i0 + 1;
#line 522 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 523 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 524 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 525 "dsbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 530 "dsbgst.f"
		dlargv_(&nr, &ab[j2 * ab_dim1 + 1], &inca, &work[j2], &ka1, &
			work[*n + j2], &ka1);

/*              apply rotations in 2nd set from the right */

#line 535 "dsbgst.f"
		i__3 = *ka - 1;
#line 535 "dsbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 536 "dsbgst.f"
		    dlartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &work[*n + j2], 
			    &work[j2], &ka1);
#line 539 "dsbgst.f"
/* L180: */
#line 539 "dsbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 544 "dsbgst.f"
		dlar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &work[
			*n + j2], &work[j2], &ka1);

#line 548 "dsbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 552 "dsbgst.f"
	    i__3 = *kb - k + 1;
#line 552 "dsbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 553 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 554 "dsbgst.f"
		if (nrt > 0) {
#line 554 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2], &work[j2], &ka1);
#line 554 "dsbgst.f"
		}
#line 558 "dsbgst.f"
/* L190: */
#line 558 "dsbgst.f"
	    }

#line 560 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 564 "dsbgst.f"
		i__3 = j1;
#line 564 "dsbgst.f"
		i__2 = ka1;
#line 564 "dsbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 565 "dsbgst.f"
		    i__4 = *n - m;
#line 565 "dsbgst.f"
		    drot_(&i__4, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j], &work[j]);
#line 567 "dsbgst.f"
/* L200: */
#line 567 "dsbgst.f"
		}
#line 568 "dsbgst.f"
	    }
#line 569 "dsbgst.f"
/* L210: */
#line 569 "dsbgst.f"
	}

#line 571 "dsbgst.f"
	i__2 = *kb - 1;
#line 571 "dsbgst.f"
	for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 572 "dsbgst.f"
	    i__3 = 1, i__4 = k - i0 + 2;
#line 572 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__4) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 576 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 577 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 578 "dsbgst.f"
		if (nrt > 0) {
#line 578 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2 - m], &work[j2 - m], &ka1);
#line 578 "dsbgst.f"
		}
#line 582 "dsbgst.f"
/* L220: */
#line 582 "dsbgst.f"
	    }
#line 583 "dsbgst.f"
/* L230: */
#line 583 "dsbgst.f"
	}

#line 585 "dsbgst.f"
	if (*kb > 1) {
#line 586 "dsbgst.f"
	    i__2 = i__ - *kb + (*ka << 1) + 1;
#line 586 "dsbgst.f"
	    for (j = *n - 1; j >= i__2; --j) {
#line 587 "dsbgst.f"
		work[*n + j - m] = work[*n + j - *ka - m];
#line 588 "dsbgst.f"
		work[j - m] = work[j - *ka - m];
#line 589 "dsbgst.f"
/* L240: */
#line 589 "dsbgst.f"
	    }
#line 590 "dsbgst.f"
	}

#line 592 "dsbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 596 "dsbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 600 "dsbgst.f"
	    bii = bb[i__ * bb_dim1 + 1];
#line 601 "dsbgst.f"
	    i__2 = i1;
#line 601 "dsbgst.f"
	    for (j = i__; j <= i__2; ++j) {
#line 602 "dsbgst.f"
		ab[j - i__ + 1 + i__ * ab_dim1] /= bii;
#line 603 "dsbgst.f"
/* L250: */
#line 603 "dsbgst.f"
	    }
/* Computing MAX */
#line 604 "dsbgst.f"
	    i__2 = 1, i__3 = i__ - *ka;
#line 604 "dsbgst.f"
	    i__4 = i__;
#line 604 "dsbgst.f"
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 605 "dsbgst.f"
		ab[i__ - j + 1 + j * ab_dim1] /= bii;
#line 606 "dsbgst.f"
/* L260: */
#line 606 "dsbgst.f"
	    }
#line 607 "dsbgst.f"
	    i__4 = i__ - 1;
#line 607 "dsbgst.f"
	    for (k = i__ - kbt; k <= i__4; ++k) {
#line 608 "dsbgst.f"
		i__2 = k;
#line 608 "dsbgst.f"
		for (j = i__ - kbt; j <= i__2; ++j) {
#line 609 "dsbgst.f"
		    ab[k - j + 1 + j * ab_dim1] = ab[k - j + 1 + j * ab_dim1] 
			    - bb[i__ - j + 1 + j * bb_dim1] * ab[i__ - k + 1 
			    + k * ab_dim1] - bb[i__ - k + 1 + k * bb_dim1] * 
			    ab[i__ - j + 1 + j * ab_dim1] + ab[i__ * ab_dim1 
			    + 1] * bb[i__ - j + 1 + j * bb_dim1] * bb[i__ - k 
			    + 1 + k * bb_dim1];
#line 614 "dsbgst.f"
/* L270: */
#line 614 "dsbgst.f"
		}
/* Computing MAX */
#line 615 "dsbgst.f"
		i__2 = 1, i__3 = i__ - *ka;
#line 615 "dsbgst.f"
		i__1 = i__ - kbt - 1;
#line 615 "dsbgst.f"
		for (j = max(i__2,i__3); j <= i__1; ++j) {
#line 616 "dsbgst.f"
		    ab[k - j + 1 + j * ab_dim1] -= bb[i__ - k + 1 + k * 
			    bb_dim1] * ab[i__ - j + 1 + j * ab_dim1];
#line 618 "dsbgst.f"
/* L280: */
#line 618 "dsbgst.f"
		}
#line 619 "dsbgst.f"
/* L290: */
#line 619 "dsbgst.f"
	    }
#line 620 "dsbgst.f"
	    i__4 = i1;
#line 620 "dsbgst.f"
	    for (j = i__; j <= i__4; ++j) {
/* Computing MAX */
#line 621 "dsbgst.f"
		i__1 = j - *ka, i__2 = i__ - kbt;
#line 621 "dsbgst.f"
		i__3 = i__ - 1;
#line 621 "dsbgst.f"
		for (k = max(i__1,i__2); k <= i__3; ++k) {
#line 622 "dsbgst.f"
		    ab[j - k + 1 + k * ab_dim1] -= bb[i__ - k + 1 + k * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1];
#line 624 "dsbgst.f"
/* L300: */
#line 624 "dsbgst.f"
		}
#line 625 "dsbgst.f"
/* L310: */
#line 625 "dsbgst.f"
	    }

#line 627 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 631 "dsbgst.f"
		i__4 = *n - m;
#line 631 "dsbgst.f"
		d__1 = 1. / bii;
#line 631 "dsbgst.f"
		dscal_(&i__4, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 632 "dsbgst.f"
		if (kbt > 0) {
#line 632 "dsbgst.f"
		    i__4 = *n - m;
#line 632 "dsbgst.f"
		    i__3 = *ldbb - 1;
#line 632 "dsbgst.f"
		    dger_(&i__4, &kbt, &c_b20, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kbt + 1 + (i__ - kbt) * bb_dim1], &i__3,
			     &x[m + 1 + (i__ - kbt) * x_dim1], ldx);
#line 632 "dsbgst.f"
		}
#line 636 "dsbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 640 "dsbgst.f"
	    ra1 = ab[i1 - i__ + 1 + i__ * ab_dim1];
#line 641 "dsbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 647 "dsbgst.f"
	i__4 = *kb - 1;
#line 647 "dsbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 648 "dsbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 653 "dsbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i-k+ka+1,i) */

#line 657 "dsbgst.f"
		    dlartg_(&ab[ka1 - k + i__ * ab_dim1], &ra1, &work[*n + 
			    i__ - k + *ka - m], &work[i__ - k + *ka - m], &ra)
			    ;

/*                 create nonzero element a(i-k+ka+1,i-k) outside the */
/*                 band and store it in WORK(i-k) */

#line 663 "dsbgst.f"
		    t = -bb[k + 1 + (i__ - k) * bb_dim1] * ra1;
#line 664 "dsbgst.f"
		    work[i__ - k] = work[*n + i__ - k + *ka - m] * t - work[
			    i__ - k + *ka - m] * ab[ka1 + (i__ - k) * ab_dim1]
			    ;
#line 666 "dsbgst.f"
		    ab[ka1 + (i__ - k) * ab_dim1] = work[i__ - k + *ka - m] * 
			    t + work[*n + i__ - k + *ka - m] * ab[ka1 + (i__ 
			    - k) * ab_dim1];
#line 668 "dsbgst.f"
		    ra1 = ra;
#line 669 "dsbgst.f"
		}
#line 670 "dsbgst.f"
	    }
/* Computing MAX */
#line 671 "dsbgst.f"
	    i__3 = 1, i__1 = k - i0 + 2;
#line 671 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__1) * ka1;
#line 672 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 673 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 674 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 675 "dsbgst.f"
		i__3 = j2, i__1 = i__ + (*ka << 1) - k + 1;
#line 675 "dsbgst.f"
		j2t = max(i__3,i__1);
#line 676 "dsbgst.f"
	    } else {
#line 677 "dsbgst.f"
		j2t = j2;
#line 678 "dsbgst.f"
	    }
#line 679 "dsbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 680 "dsbgst.f"
	    i__3 = j1;
#line 680 "dsbgst.f"
	    i__1 = ka1;
#line 680 "dsbgst.f"
	    for (j = j2t; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j-m) */

#line 685 "dsbgst.f"
		work[j - m] *= ab[ka1 + (j - *ka + 1) * ab_dim1];
#line 686 "dsbgst.f"
		ab[ka1 + (j - *ka + 1) * ab_dim1] = work[*n + j - m] * ab[ka1 
			+ (j - *ka + 1) * ab_dim1];
#line 687 "dsbgst.f"
/* L320: */
#line 687 "dsbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 692 "dsbgst.f"
	    if (nrt > 0) {
#line 692 "dsbgst.f"
		dlargv_(&nrt, &ab[ka1 + (j2t - *ka) * ab_dim1], &inca, &work[
			j2t - m], &ka1, &work[*n + j2t - m], &ka1);
#line 692 "dsbgst.f"
	    }
#line 695 "dsbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 699 "dsbgst.f"
		i__1 = *ka - 1;
#line 699 "dsbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 700 "dsbgst.f"
		    dlartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &work[*n + j2 
			    - m], &work[j2 - m], &ka1);
#line 703 "dsbgst.f"
/* L330: */
#line 703 "dsbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 708 "dsbgst.f"
		dlar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &work[*n + j2 - m], 
			&work[j2 - m], &ka1);

#line 711 "dsbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 715 "dsbgst.f"
	    i__1 = *kb - k + 1;
#line 715 "dsbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 716 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 717 "dsbgst.f"
		if (nrt > 0) {
#line 717 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2 - m], &work[j2 - m], &ka1);
#line 717 "dsbgst.f"
		}
#line 721 "dsbgst.f"
/* L340: */
#line 721 "dsbgst.f"
	    }

#line 723 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 727 "dsbgst.f"
		i__1 = j1;
#line 727 "dsbgst.f"
		i__3 = ka1;
#line 727 "dsbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 728 "dsbgst.f"
		    i__2 = *n - m;
#line 728 "dsbgst.f"
		    drot_(&i__2, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j - m], &work[j 
			    - m]);
#line 730 "dsbgst.f"
/* L350: */
#line 730 "dsbgst.f"
		}
#line 731 "dsbgst.f"
	    }
#line 732 "dsbgst.f"
/* L360: */
#line 732 "dsbgst.f"
	}

#line 734 "dsbgst.f"
	if (update) {
#line 735 "dsbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt+ka+1,i-kbt) outside the */
/*              band and store it in WORK(i-kbt) */

#line 740 "dsbgst.f"
		work[i__ - kbt] = -bb[kbt + 1 + (i__ - kbt) * bb_dim1] * ra1;
#line 741 "dsbgst.f"
	    }
#line 742 "dsbgst.f"
	}

#line 744 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 745 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 746 "dsbgst.f"
		i__4 = 2, i__3 = k - i0 + 1;
#line 746 "dsbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 747 "dsbgst.f"
	    } else {
/* Computing MAX */
#line 748 "dsbgst.f"
		i__4 = 1, i__3 = k - i0 + 1;
#line 748 "dsbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 749 "dsbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 753 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 754 "dsbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 755 "dsbgst.f"
		if (nrt > 0) {
#line 755 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + (j2 - *ka) * ab_dim1], &
			    inca, &ab[ka1 - l + (j2 - *ka + 1) * ab_dim1], &
			    inca, &work[*n + j2 - *ka], &work[j2 - *ka], &ka1)
			    ;
#line 755 "dsbgst.f"
		}
#line 759 "dsbgst.f"
/* L370: */
#line 759 "dsbgst.f"
	    }
#line 760 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 761 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 762 "dsbgst.f"
	    i__4 = j2;
#line 762 "dsbgst.f"
	    i__3 = -ka1;
#line 762 "dsbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 763 "dsbgst.f"
		work[j] = work[j - *ka];
#line 764 "dsbgst.f"
		work[*n + j] = work[*n + j - *ka];
#line 765 "dsbgst.f"
/* L380: */
#line 765 "dsbgst.f"
	    }
#line 766 "dsbgst.f"
	    i__3 = j1;
#line 766 "dsbgst.f"
	    i__4 = ka1;
#line 766 "dsbgst.f"
	    for (j = j2; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j) */

#line 771 "dsbgst.f"
		work[j] *= ab[ka1 + (j - *ka + 1) * ab_dim1];
#line 772 "dsbgst.f"
		ab[ka1 + (j - *ka + 1) * ab_dim1] = work[*n + j] * ab[ka1 + (
			j - *ka + 1) * ab_dim1];
#line 773 "dsbgst.f"
/* L390: */
#line 773 "dsbgst.f"
	    }
#line 774 "dsbgst.f"
	    if (update) {
#line 775 "dsbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 775 "dsbgst.f"
		    work[i__ - k + *ka] = work[i__ - k];
#line 775 "dsbgst.f"
		}
#line 777 "dsbgst.f"
	    }
#line 778 "dsbgst.f"
/* L400: */
#line 778 "dsbgst.f"
	}

#line 780 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 781 "dsbgst.f"
	    i__4 = 1, i__3 = k - i0 + 1;
#line 781 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 782 "dsbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 783 "dsbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 784 "dsbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 789 "dsbgst.f"
		dlargv_(&nr, &ab[ka1 + (j2 - *ka) * ab_dim1], &inca, &work[j2]
			, &ka1, &work[*n + j2], &ka1);

/*              apply rotations in 2nd set from the left */

#line 794 "dsbgst.f"
		i__4 = *ka - 1;
#line 794 "dsbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 795 "dsbgst.f"
		    dlartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &work[*n + j2]
			    , &work[j2], &ka1);
#line 798 "dsbgst.f"
/* L410: */
#line 798 "dsbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 803 "dsbgst.f"
		dlar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &work[*n + j2], &
			work[j2], &ka1);

#line 806 "dsbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 810 "dsbgst.f"
	    i__4 = *kb - k + 1;
#line 810 "dsbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 811 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 812 "dsbgst.f"
		if (nrt > 0) {
#line 812 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2], &work[j2], &ka1);
#line 812 "dsbgst.f"
		}
#line 816 "dsbgst.f"
/* L420: */
#line 816 "dsbgst.f"
	    }

#line 818 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 822 "dsbgst.f"
		i__4 = j1;
#line 822 "dsbgst.f"
		i__3 = ka1;
#line 822 "dsbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 823 "dsbgst.f"
		    i__1 = *n - m;
#line 823 "dsbgst.f"
		    drot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j], &work[j]);
#line 825 "dsbgst.f"
/* L430: */
#line 825 "dsbgst.f"
		}
#line 826 "dsbgst.f"
	    }
#line 827 "dsbgst.f"
/* L440: */
#line 827 "dsbgst.f"
	}

#line 829 "dsbgst.f"
	i__3 = *kb - 1;
#line 829 "dsbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 830 "dsbgst.f"
	    i__4 = 1, i__1 = k - i0 + 2;
#line 830 "dsbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 834 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 835 "dsbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 836 "dsbgst.f"
		if (nrt > 0) {
#line 836 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2 - m], &work[j2 - m], &ka1);
#line 836 "dsbgst.f"
		}
#line 840 "dsbgst.f"
/* L450: */
#line 840 "dsbgst.f"
	    }
#line 841 "dsbgst.f"
/* L460: */
#line 841 "dsbgst.f"
	}

#line 843 "dsbgst.f"
	if (*kb > 1) {
#line 844 "dsbgst.f"
	    i__3 = i__ - *kb + (*ka << 1) + 1;
#line 844 "dsbgst.f"
	    for (j = *n - 1; j >= i__3; --j) {
#line 845 "dsbgst.f"
		work[*n + j - m] = work[*n + j - *ka - m];
#line 846 "dsbgst.f"
		work[j - m] = work[j - *ka - m];
#line 847 "dsbgst.f"
/* L470: */
#line 847 "dsbgst.f"
	    }
#line 848 "dsbgst.f"
	}

#line 850 "dsbgst.f"
    }

#line 852 "dsbgst.f"
    goto L10;

#line 854 "dsbgst.f"
L480:

/*     **************************** Phase 2 ***************************** */

/*     The logical structure of this phase is: */

/*     UPDATE = .TRUE. */
/*     DO I = 1, M */
/*        use S(i) to update A and create a new bulge */
/*        apply rotations to push all bulges KA positions upward */
/*     END DO */
/*     UPDATE = .FALSE. */
/*     DO I = M - KA - 1, 2, -1 */
/*        apply rotations to push all bulges KA positions upward */
/*     END DO */

/*     To avoid duplicating code, the two loops are merged. */

#line 872 "dsbgst.f"
    update = TRUE_;
#line 873 "dsbgst.f"
    i__ = 0;
#line 874 "dsbgst.f"
L490:
#line 875 "dsbgst.f"
    if (update) {
#line 876 "dsbgst.f"
	++i__;
/* Computing MIN */
#line 877 "dsbgst.f"
	i__3 = *kb, i__4 = m - i__;
#line 877 "dsbgst.f"
	kbt = min(i__3,i__4);
#line 878 "dsbgst.f"
	i0 = i__ + 1;
/* Computing MAX */
#line 879 "dsbgst.f"
	i__3 = 1, i__4 = i__ - *ka;
#line 879 "dsbgst.f"
	i1 = max(i__3,i__4);
#line 880 "dsbgst.f"
	i2 = i__ + kbt - ka1;
#line 881 "dsbgst.f"
	if (i__ > m) {
#line 882 "dsbgst.f"
	    update = FALSE_;
#line 883 "dsbgst.f"
	    --i__;
#line 884 "dsbgst.f"
	    i0 = m + 1;
#line 885 "dsbgst.f"
	    if (*ka == 0) {
#line 885 "dsbgst.f"
		return 0;
#line 885 "dsbgst.f"
	    }
#line 887 "dsbgst.f"
	    goto L490;
#line 888 "dsbgst.f"
	}
#line 889 "dsbgst.f"
    } else {
#line 890 "dsbgst.f"
	i__ -= *ka;
#line 891 "dsbgst.f"
	if (i__ < 2) {
#line 891 "dsbgst.f"
	    return 0;
#line 891 "dsbgst.f"
	}
#line 893 "dsbgst.f"
    }

#line 895 "dsbgst.f"
    if (i__ < m - kbt) {
#line 896 "dsbgst.f"
	nx = m;
#line 897 "dsbgst.f"
    } else {
#line 898 "dsbgst.f"
	nx = *n;
#line 899 "dsbgst.f"
    }

#line 901 "dsbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 905 "dsbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 909 "dsbgst.f"
	    bii = bb[kb1 + i__ * bb_dim1];
#line 910 "dsbgst.f"
	    i__3 = i__;
#line 910 "dsbgst.f"
	    for (j = i1; j <= i__3; ++j) {
#line 911 "dsbgst.f"
		ab[j - i__ + ka1 + i__ * ab_dim1] /= bii;
#line 912 "dsbgst.f"
/* L500: */
#line 912 "dsbgst.f"
	    }
/* Computing MIN */
#line 913 "dsbgst.f"
	    i__4 = *n, i__1 = i__ + *ka;
#line 913 "dsbgst.f"
	    i__3 = min(i__4,i__1);
#line 913 "dsbgst.f"
	    for (j = i__; j <= i__3; ++j) {
#line 914 "dsbgst.f"
		ab[i__ - j + ka1 + j * ab_dim1] /= bii;
#line 915 "dsbgst.f"
/* L510: */
#line 915 "dsbgst.f"
	    }
#line 916 "dsbgst.f"
	    i__3 = i__ + kbt;
#line 916 "dsbgst.f"
	    for (k = i__ + 1; k <= i__3; ++k) {
#line 917 "dsbgst.f"
		i__4 = i__ + kbt;
#line 917 "dsbgst.f"
		for (j = k; j <= i__4; ++j) {
#line 918 "dsbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] = ab[k - j + ka1 + j * 
			    ab_dim1] - bb[i__ - j + kb1 + j * bb_dim1] * ab[
			    i__ - k + ka1 + k * ab_dim1] - bb[i__ - k + kb1 + 
			    k * bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1] + 
			    ab[ka1 + i__ * ab_dim1] * bb[i__ - j + kb1 + j * 
			    bb_dim1] * bb[i__ - k + kb1 + k * bb_dim1];
#line 923 "dsbgst.f"
/* L520: */
#line 923 "dsbgst.f"
		}
/* Computing MIN */
#line 924 "dsbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 924 "dsbgst.f"
		i__4 = min(i__1,i__2);
#line 924 "dsbgst.f"
		for (j = i__ + kbt + 1; j <= i__4; ++j) {
#line 925 "dsbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] -= bb[i__ - k + kb1 + k * 
			    bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1];
#line 927 "dsbgst.f"
/* L530: */
#line 927 "dsbgst.f"
		}
#line 928 "dsbgst.f"
/* L540: */
#line 928 "dsbgst.f"
	    }
#line 929 "dsbgst.f"
	    i__3 = i__;
#line 929 "dsbgst.f"
	    for (j = i1; j <= i__3; ++j) {
/* Computing MIN */
#line 930 "dsbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 930 "dsbgst.f"
		i__4 = min(i__1,i__2);
#line 930 "dsbgst.f"
		for (k = i__ + 1; k <= i__4; ++k) {
#line 931 "dsbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] -= bb[i__ - k + kb1 + k * 
			    bb_dim1] * ab[j - i__ + ka1 + i__ * ab_dim1];
#line 933 "dsbgst.f"
/* L550: */
#line 933 "dsbgst.f"
		}
#line 934 "dsbgst.f"
/* L560: */
#line 934 "dsbgst.f"
	    }

#line 936 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 940 "dsbgst.f"
		d__1 = 1. / bii;
#line 940 "dsbgst.f"
		dscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 941 "dsbgst.f"
		if (kbt > 0) {
#line 941 "dsbgst.f"
		    i__3 = *ldbb - 1;
#line 941 "dsbgst.f"
		    dger_(&nx, &kbt, &c_b20, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    *kb + (i__ + 1) * bb_dim1], &i__3, &x[(i__ + 1) * 
			    x_dim1 + 1], ldx);
#line 941 "dsbgst.f"
		}
#line 944 "dsbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 948 "dsbgst.f"
	    ra1 = ab[i1 - i__ + ka1 + i__ * ab_dim1];
#line 949 "dsbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 954 "dsbgst.f"
	i__3 = *kb - 1;
#line 954 "dsbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 955 "dsbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 960 "dsbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i+k-ka-1,i) */

#line 964 "dsbgst.f"
		    dlartg_(&ab[k + 1 + i__ * ab_dim1], &ra1, &work[*n + i__ 
			    + k - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k-ka-1,i+k) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 970 "dsbgst.f"
		    t = -bb[kb1 - k + (i__ + k) * bb_dim1] * ra1;
#line 971 "dsbgst.f"
		    work[m - *kb + i__ + k] = work[*n + i__ + k - *ka] * t - 
			    work[i__ + k - *ka] * ab[(i__ + k) * ab_dim1 + 1];
#line 973 "dsbgst.f"
		    ab[(i__ + k) * ab_dim1 + 1] = work[i__ + k - *ka] * t + 
			    work[*n + i__ + k - *ka] * ab[(i__ + k) * ab_dim1 
			    + 1];
#line 975 "dsbgst.f"
		    ra1 = ra;
#line 976 "dsbgst.f"
		}
#line 977 "dsbgst.f"
	    }
/* Computing MAX */
#line 978 "dsbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 978 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;
#line 979 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 980 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 981 "dsbgst.f"
	    if (update) {
/* Computing MIN */
#line 982 "dsbgst.f"
		i__4 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 982 "dsbgst.f"
		j2t = min(i__4,i__1);
#line 983 "dsbgst.f"
	    } else {
#line 984 "dsbgst.f"
		j2t = j2;
#line 985 "dsbgst.f"
	    }
#line 986 "dsbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 987 "dsbgst.f"
	    i__4 = j2t;
#line 987 "dsbgst.f"
	    i__1 = ka1;
#line 987 "dsbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__4 : j <= i__4; j += i__1) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(j) */

#line 992 "dsbgst.f"
		work[j] *= ab[(j + *ka - 1) * ab_dim1 + 1];
#line 993 "dsbgst.f"
		ab[(j + *ka - 1) * ab_dim1 + 1] = work[*n + j] * ab[(j + *ka 
			- 1) * ab_dim1 + 1];
#line 994 "dsbgst.f"
/* L570: */
#line 994 "dsbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 999 "dsbgst.f"
	    if (nrt > 0) {
#line 999 "dsbgst.f"
		dlargv_(&nrt, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[j1],
			 &ka1, &work[*n + j1], &ka1);
#line 999 "dsbgst.f"
	    }
#line 1002 "dsbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 1006 "dsbgst.f"
		i__1 = *ka - 1;
#line 1006 "dsbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1007 "dsbgst.f"
		    dlartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &work[*n 
			    + j1], &work[j1], &ka1);
#line 1010 "dsbgst.f"
/* L580: */
#line 1010 "dsbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1015 "dsbgst.f"
		dlar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &work[*n + 
			j1], &work[j1], &ka1);

#line 1019 "dsbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 1023 "dsbgst.f"
	    i__1 = *kb - k + 1;
#line 1023 "dsbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1024 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1025 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1026 "dsbgst.f"
		if (nrt > 0) {
#line 1026 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + j1t], &
			    work[j1t], &ka1);
#line 1026 "dsbgst.f"
		}
#line 1030 "dsbgst.f"
/* L590: */
#line 1030 "dsbgst.f"
	    }

#line 1032 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1036 "dsbgst.f"
		i__1 = j2;
#line 1036 "dsbgst.f"
		i__4 = ka1;
#line 1036 "dsbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__1 : j <= i__1; j += i__4) {
#line 1037 "dsbgst.f"
		    drot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + j], &work[j]);
#line 1039 "dsbgst.f"
/* L600: */
#line 1039 "dsbgst.f"
		}
#line 1040 "dsbgst.f"
	    }
#line 1041 "dsbgst.f"
/* L610: */
#line 1041 "dsbgst.f"
	}

#line 1043 "dsbgst.f"
	if (update) {
#line 1044 "dsbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt-ka-1,i+kbt) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1049 "dsbgst.f"
		work[m - *kb + i__ + kbt] = -bb[kb1 - kbt + (i__ + kbt) * 
			bb_dim1] * ra1;
#line 1050 "dsbgst.f"
	    }
#line 1051 "dsbgst.f"
	}

#line 1053 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1054 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 1055 "dsbgst.f"
		i__3 = 2, i__4 = k + i0 - m;
#line 1055 "dsbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1056 "dsbgst.f"
	    } else {
/* Computing MAX */
#line 1057 "dsbgst.f"
		i__3 = 1, i__4 = k + i0 - m;
#line 1057 "dsbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1058 "dsbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 1062 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1063 "dsbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1064 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1065 "dsbgst.f"
		if (nrt > 0) {
#line 1065 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + (j1t + *ka) * ab_dim1], &inca, &ab[
			    l + 1 + (j1t + *ka - 1) * ab_dim1], &inca, &work[*
			    n + m - *kb + j1t + *ka], &work[m - *kb + j1t + *
			    ka], &ka1);
#line 1065 "dsbgst.f"
		}
#line 1070 "dsbgst.f"
/* L620: */
#line 1070 "dsbgst.f"
	    }
#line 1071 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1072 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1073 "dsbgst.f"
	    i__3 = j2;
#line 1073 "dsbgst.f"
	    i__4 = ka1;
#line 1073 "dsbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1074 "dsbgst.f"
		work[m - *kb + j] = work[m - *kb + j + *ka];
#line 1075 "dsbgst.f"
		work[*n + m - *kb + j] = work[*n + m - *kb + j + *ka];
#line 1076 "dsbgst.f"
/* L630: */
#line 1076 "dsbgst.f"
	    }
#line 1077 "dsbgst.f"
	    i__4 = j2;
#line 1077 "dsbgst.f"
	    i__3 = ka1;
#line 1077 "dsbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1082 "dsbgst.f"
		work[m - *kb + j] *= ab[(j + *ka - 1) * ab_dim1 + 1];
#line 1083 "dsbgst.f"
		ab[(j + *ka - 1) * ab_dim1 + 1] = work[*n + m - *kb + j] * ab[
			(j + *ka - 1) * ab_dim1 + 1];
#line 1084 "dsbgst.f"
/* L640: */
#line 1084 "dsbgst.f"
	    }
#line 1085 "dsbgst.f"
	    if (update) {
#line 1086 "dsbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1086 "dsbgst.f"
		    work[m - *kb + i__ + k - *ka] = work[m - *kb + i__ + k];
#line 1086 "dsbgst.f"
		}
#line 1088 "dsbgst.f"
	    }
#line 1089 "dsbgst.f"
/* L650: */
#line 1089 "dsbgst.f"
	}

#line 1091 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1092 "dsbgst.f"
	    i__3 = 1, i__4 = k + i0 - m;
#line 1092 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1093 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1094 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1095 "dsbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1100 "dsbgst.f"
		dlargv_(&nr, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[m - *
			kb + j1], &ka1, &work[*n + m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the left */

#line 1105 "dsbgst.f"
		i__3 = *ka - 1;
#line 1105 "dsbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 1106 "dsbgst.f"
		    dlartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &work[*n 
			    + m - *kb + j1], &work[m - *kb + j1], &ka1);
#line 1109 "dsbgst.f"
/* L660: */
#line 1109 "dsbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1114 "dsbgst.f"
		dlar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &work[*n + 
			m - *kb + j1], &work[m - *kb + j1], &ka1);

#line 1118 "dsbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 1122 "dsbgst.f"
	    i__3 = *kb - k + 1;
#line 1122 "dsbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 1123 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1124 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1125 "dsbgst.f"
		if (nrt > 0) {
#line 1125 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + m - *kb + 
			    j1t], &work[m - *kb + j1t], &ka1);
#line 1125 "dsbgst.f"
		}
#line 1130 "dsbgst.f"
/* L670: */
#line 1130 "dsbgst.f"
	    }

#line 1132 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1136 "dsbgst.f"
		i__3 = j2;
#line 1136 "dsbgst.f"
		i__4 = ka1;
#line 1136 "dsbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1137 "dsbgst.f"
		    drot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + m - *kb + j], &work[m - *
			    kb + j]);
#line 1139 "dsbgst.f"
/* L680: */
#line 1139 "dsbgst.f"
		}
#line 1140 "dsbgst.f"
	    }
#line 1141 "dsbgst.f"
/* L690: */
#line 1141 "dsbgst.f"
	}

#line 1143 "dsbgst.f"
	i__4 = *kb - 1;
#line 1143 "dsbgst.f"
	for (k = 1; k <= i__4; ++k) {
/* Computing MAX */
#line 1144 "dsbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1144 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 1148 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1149 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1150 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1151 "dsbgst.f"
		if (nrt > 0) {
#line 1151 "dsbgst.f"
		    dlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + j1t], &
			    work[j1t], &ka1);
#line 1151 "dsbgst.f"
		}
#line 1155 "dsbgst.f"
/* L700: */
#line 1155 "dsbgst.f"
	    }
#line 1156 "dsbgst.f"
/* L710: */
#line 1156 "dsbgst.f"
	}

#line 1158 "dsbgst.f"
	if (*kb > 1) {
/* Computing MIN */
#line 1159 "dsbgst.f"
	    i__3 = i__ + *kb;
#line 1159 "dsbgst.f"
	    i__4 = min(i__3,m) - (*ka << 1) - 1;
#line 1159 "dsbgst.f"
	    for (j = 2; j <= i__4; ++j) {
#line 1160 "dsbgst.f"
		work[*n + j] = work[*n + j + *ka];
#line 1161 "dsbgst.f"
		work[j] = work[j + *ka];
#line 1162 "dsbgst.f"
/* L720: */
#line 1162 "dsbgst.f"
	    }
#line 1163 "dsbgst.f"
	}

#line 1165 "dsbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 1169 "dsbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 1173 "dsbgst.f"
	    bii = bb[i__ * bb_dim1 + 1];
#line 1174 "dsbgst.f"
	    i__4 = i__;
#line 1174 "dsbgst.f"
	    for (j = i1; j <= i__4; ++j) {
#line 1175 "dsbgst.f"
		ab[i__ - j + 1 + j * ab_dim1] /= bii;
#line 1176 "dsbgst.f"
/* L730: */
#line 1176 "dsbgst.f"
	    }
/* Computing MIN */
#line 1177 "dsbgst.f"
	    i__3 = *n, i__1 = i__ + *ka;
#line 1177 "dsbgst.f"
	    i__4 = min(i__3,i__1);
#line 1177 "dsbgst.f"
	    for (j = i__; j <= i__4; ++j) {
#line 1178 "dsbgst.f"
		ab[j - i__ + 1 + i__ * ab_dim1] /= bii;
#line 1179 "dsbgst.f"
/* L740: */
#line 1179 "dsbgst.f"
	    }
#line 1180 "dsbgst.f"
	    i__4 = i__ + kbt;
#line 1180 "dsbgst.f"
	    for (k = i__ + 1; k <= i__4; ++k) {
#line 1181 "dsbgst.f"
		i__3 = i__ + kbt;
#line 1181 "dsbgst.f"
		for (j = k; j <= i__3; ++j) {
#line 1182 "dsbgst.f"
		    ab[j - k + 1 + k * ab_dim1] = ab[j - k + 1 + k * ab_dim1] 
			    - bb[j - i__ + 1 + i__ * bb_dim1] * ab[k - i__ + 
			    1 + i__ * ab_dim1] - bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1] + ab[
			    i__ * ab_dim1 + 1] * bb[j - i__ + 1 + i__ * 
			    bb_dim1] * bb[k - i__ + 1 + i__ * bb_dim1];
#line 1187 "dsbgst.f"
/* L750: */
#line 1187 "dsbgst.f"
		}
/* Computing MIN */
#line 1188 "dsbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 1188 "dsbgst.f"
		i__3 = min(i__1,i__2);
#line 1188 "dsbgst.f"
		for (j = i__ + kbt + 1; j <= i__3; ++j) {
#line 1189 "dsbgst.f"
		    ab[j - k + 1 + k * ab_dim1] -= bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1];
#line 1191 "dsbgst.f"
/* L760: */
#line 1191 "dsbgst.f"
		}
#line 1192 "dsbgst.f"
/* L770: */
#line 1192 "dsbgst.f"
	    }
#line 1193 "dsbgst.f"
	    i__4 = i__;
#line 1193 "dsbgst.f"
	    for (j = i1; j <= i__4; ++j) {
/* Computing MIN */
#line 1194 "dsbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 1194 "dsbgst.f"
		i__3 = min(i__1,i__2);
#line 1194 "dsbgst.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 1195 "dsbgst.f"
		    ab[k - j + 1 + j * ab_dim1] -= bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[i__ - j + 1 + j * ab_dim1];
#line 1197 "dsbgst.f"
/* L780: */
#line 1197 "dsbgst.f"
		}
#line 1198 "dsbgst.f"
/* L790: */
#line 1198 "dsbgst.f"
	    }

#line 1200 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 1204 "dsbgst.f"
		d__1 = 1. / bii;
#line 1204 "dsbgst.f"
		dscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 1205 "dsbgst.f"
		if (kbt > 0) {
#line 1205 "dsbgst.f"
		    dger_(&nx, &kbt, &c_b20, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    i__ * bb_dim1 + 2], &c__1, &x[(i__ + 1) * x_dim1 
			    + 1], ldx);
#line 1205 "dsbgst.f"
		}
#line 1208 "dsbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 1212 "dsbgst.f"
	    ra1 = ab[i__ - i1 + 1 + i1 * ab_dim1];
#line 1213 "dsbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 1218 "dsbgst.f"
	i__4 = *kb - 1;
#line 1218 "dsbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 1219 "dsbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 1224 "dsbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i,i+k-ka-1) */

#line 1228 "dsbgst.f"
		    dlartg_(&ab[ka1 - k + (i__ + k - *ka) * ab_dim1], &ra1, &
			    work[*n + i__ + k - *ka], &work[i__ + k - *ka], &
			    ra);

/*                 create nonzero element a(i+k,i+k-ka-1) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 1234 "dsbgst.f"
		    t = -bb[k + 1 + i__ * bb_dim1] * ra1;
#line 1235 "dsbgst.f"
		    work[m - *kb + i__ + k] = work[*n + i__ + k - *ka] * t - 
			    work[i__ + k - *ka] * ab[ka1 + (i__ + k - *ka) * 
			    ab_dim1];
#line 1237 "dsbgst.f"
		    ab[ka1 + (i__ + k - *ka) * ab_dim1] = work[i__ + k - *ka] 
			    * t + work[*n + i__ + k - *ka] * ab[ka1 + (i__ + 
			    k - *ka) * ab_dim1];
#line 1239 "dsbgst.f"
		    ra1 = ra;
#line 1240 "dsbgst.f"
		}
#line 1241 "dsbgst.f"
	    }
/* Computing MAX */
#line 1242 "dsbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1242 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;
#line 1243 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1244 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1245 "dsbgst.f"
	    if (update) {
/* Computing MIN */
#line 1246 "dsbgst.f"
		i__3 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1246 "dsbgst.f"
		j2t = min(i__3,i__1);
#line 1247 "dsbgst.f"
	    } else {
#line 1248 "dsbgst.f"
		j2t = j2;
#line 1249 "dsbgst.f"
	    }
#line 1250 "dsbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1251 "dsbgst.f"
	    i__3 = j2t;
#line 1251 "dsbgst.f"
	    i__1 = ka1;
#line 1251 "dsbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(j) */

#line 1256 "dsbgst.f"
		work[j] *= ab[ka1 + (j - 1) * ab_dim1];
#line 1257 "dsbgst.f"
		ab[ka1 + (j - 1) * ab_dim1] = work[*n + j] * ab[ka1 + (j - 1) 
			* ab_dim1];
#line 1258 "dsbgst.f"
/* L800: */
#line 1258 "dsbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1263 "dsbgst.f"
	    if (nrt > 0) {
#line 1263 "dsbgst.f"
		dlargv_(&nrt, &ab[ka1 + j1 * ab_dim1], &inca, &work[j1], &ka1,
			 &work[*n + j1], &ka1);
#line 1263 "dsbgst.f"
	    }
#line 1266 "dsbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 1270 "dsbgst.f"
		i__1 = *ka - 1;
#line 1270 "dsbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1271 "dsbgst.f"
		    dlartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &work[*n + j1], &
			    work[j1], &ka1);
#line 1273 "dsbgst.f"
/* L810: */
#line 1273 "dsbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1278 "dsbgst.f"
		dlar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &work[*n + j1]
			, &work[j1], &ka1);

#line 1282 "dsbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 1286 "dsbgst.f"
	    i__1 = *kb - k + 1;
#line 1286 "dsbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1287 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1288 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1289 "dsbgst.f"
		if (nrt > 0) {
#line 1289 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + j1t], &work[j1t], &ka1);
#line 1289 "dsbgst.f"
		}
#line 1293 "dsbgst.f"
/* L820: */
#line 1293 "dsbgst.f"
	    }

#line 1295 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1299 "dsbgst.f"
		i__1 = j2;
#line 1299 "dsbgst.f"
		i__3 = ka1;
#line 1299 "dsbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 1300 "dsbgst.f"
		    drot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + j], &work[j]);
#line 1302 "dsbgst.f"
/* L830: */
#line 1302 "dsbgst.f"
		}
#line 1303 "dsbgst.f"
	    }
#line 1304 "dsbgst.f"
/* L840: */
#line 1304 "dsbgst.f"
	}

#line 1306 "dsbgst.f"
	if (update) {
#line 1307 "dsbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt,i+kbt-ka-1) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1312 "dsbgst.f"
		work[m - *kb + i__ + kbt] = -bb[kbt + 1 + i__ * bb_dim1] * 
			ra1;
#line 1313 "dsbgst.f"
	    }
#line 1314 "dsbgst.f"
	}

#line 1316 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1317 "dsbgst.f"
	    if (update) {
/* Computing MAX */
#line 1318 "dsbgst.f"
		i__4 = 2, i__3 = k + i0 - m;
#line 1318 "dsbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1319 "dsbgst.f"
	    } else {
/* Computing MAX */
#line 1320 "dsbgst.f"
		i__4 = 1, i__3 = k + i0 - m;
#line 1320 "dsbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1321 "dsbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 1325 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1326 "dsbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1327 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1328 "dsbgst.f"
		if (nrt > 0) {
#line 1328 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + (j1t + l - 1) * ab_dim1], 
			    &inca, &ab[ka1 - l + (j1t + l - 1) * ab_dim1], &
			    inca, &work[*n + m - *kb + j1t + *ka], &work[m - *
			    kb + j1t + *ka], &ka1);
#line 1328 "dsbgst.f"
		}
#line 1333 "dsbgst.f"
/* L850: */
#line 1333 "dsbgst.f"
	    }
#line 1334 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1335 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1336 "dsbgst.f"
	    i__4 = j2;
#line 1336 "dsbgst.f"
	    i__3 = ka1;
#line 1336 "dsbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1337 "dsbgst.f"
		work[m - *kb + j] = work[m - *kb + j + *ka];
#line 1338 "dsbgst.f"
		work[*n + m - *kb + j] = work[*n + m - *kb + j + *ka];
#line 1339 "dsbgst.f"
/* L860: */
#line 1339 "dsbgst.f"
	    }
#line 1340 "dsbgst.f"
	    i__3 = j2;
#line 1340 "dsbgst.f"
	    i__4 = ka1;
#line 1340 "dsbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1345 "dsbgst.f"
		work[m - *kb + j] *= ab[ka1 + (j - 1) * ab_dim1];
#line 1346 "dsbgst.f"
		ab[ka1 + (j - 1) * ab_dim1] = work[*n + m - *kb + j] * ab[ka1 
			+ (j - 1) * ab_dim1];
#line 1347 "dsbgst.f"
/* L870: */
#line 1347 "dsbgst.f"
	    }
#line 1348 "dsbgst.f"
	    if (update) {
#line 1349 "dsbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1349 "dsbgst.f"
		    work[m - *kb + i__ + k - *ka] = work[m - *kb + i__ + k];
#line 1349 "dsbgst.f"
		}
#line 1351 "dsbgst.f"
	    }
#line 1352 "dsbgst.f"
/* L880: */
#line 1352 "dsbgst.f"
	}

#line 1354 "dsbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1355 "dsbgst.f"
	    i__4 = 1, i__3 = k + i0 - m;
#line 1355 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1356 "dsbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1357 "dsbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1358 "dsbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1363 "dsbgst.f"
		dlargv_(&nr, &ab[ka1 + j1 * ab_dim1], &inca, &work[m - *kb + 
			j1], &ka1, &work[*n + m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the right */

#line 1368 "dsbgst.f"
		i__4 = *ka - 1;
#line 1368 "dsbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 1369 "dsbgst.f"
		    dlartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &work[*n + m - *kb 
			    + j1], &work[m - *kb + j1], &ka1);
#line 1372 "dsbgst.f"
/* L890: */
#line 1372 "dsbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1377 "dsbgst.f"
		dlar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &work[*n + m 
			- *kb + j1], &work[m - *kb + j1], &ka1);

#line 1381 "dsbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 1385 "dsbgst.f"
	    i__4 = *kb - k + 1;
#line 1385 "dsbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 1386 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1387 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1388 "dsbgst.f"
		if (nrt > 0) {
#line 1388 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + m - *kb + j1t], &work[m - *kb 
			    + j1t], &ka1);
#line 1388 "dsbgst.f"
		}
#line 1393 "dsbgst.f"
/* L900: */
#line 1393 "dsbgst.f"
	    }

#line 1395 "dsbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1399 "dsbgst.f"
		i__4 = j2;
#line 1399 "dsbgst.f"
		i__3 = ka1;
#line 1399 "dsbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1400 "dsbgst.f"
		    drot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + m - *kb + j], &work[m - *
			    kb + j]);
#line 1402 "dsbgst.f"
/* L910: */
#line 1402 "dsbgst.f"
		}
#line 1403 "dsbgst.f"
	    }
#line 1404 "dsbgst.f"
/* L920: */
#line 1404 "dsbgst.f"
	}

#line 1406 "dsbgst.f"
	i__3 = *kb - 1;
#line 1406 "dsbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1407 "dsbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1407 "dsbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 1411 "dsbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1412 "dsbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1413 "dsbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1414 "dsbgst.f"
		if (nrt > 0) {
#line 1414 "dsbgst.f"
		    dlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + j1t], &work[j1t], &ka1);
#line 1414 "dsbgst.f"
		}
#line 1418 "dsbgst.f"
/* L930: */
#line 1418 "dsbgst.f"
	    }
#line 1419 "dsbgst.f"
/* L940: */
#line 1419 "dsbgst.f"
	}

#line 1421 "dsbgst.f"
	if (*kb > 1) {
/* Computing MIN */
#line 1422 "dsbgst.f"
	    i__4 = i__ + *kb;
#line 1422 "dsbgst.f"
	    i__3 = min(i__4,m) - (*ka << 1) - 1;
#line 1422 "dsbgst.f"
	    for (j = 2; j <= i__3; ++j) {
#line 1423 "dsbgst.f"
		work[*n + j] = work[*n + j + *ka];
#line 1424 "dsbgst.f"
		work[j] = work[j + *ka];
#line 1425 "dsbgst.f"
/* L950: */
#line 1425 "dsbgst.f"
	    }
#line 1426 "dsbgst.f"
	}

#line 1428 "dsbgst.f"
    }

#line 1430 "dsbgst.f"
    goto L490;

/*     End of DSBGST */

} /* dsbgst_ */

