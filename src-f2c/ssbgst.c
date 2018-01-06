#line 1 "ssbgst.f"
/* ssbgst.f -- translated by f2c (version 20100827).
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

#line 1 "ssbgst.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b20 = -1.;

/* > \brief \b SSBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSBGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, */
/*                          LDX, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBGST reduces a real symmetric-definite banded generalized */
/* > eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y, */
/* > such that C has the same bandwidth as A. */
/* > */
/* > B must have been previously factorized as S**T*S by SPBSTF, using a */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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
/* >          BB is REAL array, dimension (LDBB,N) */
/* >          The banded factor S from the split Cholesky factorization of */
/* >          B, as returned by SPBSTF, stored in the first KB+1 rows of */
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
/* >          X is REAL array, dimension (LDX,N) */
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
/* >          WORK is REAL array, dimension (2*N) */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssbgst_(char *vect, char *uplo, integer *n, integer *ka, 
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
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), srot_(integer *, doublereal *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper, wantx;
    extern /* Subroutine */ int slar2v_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    static logical update;
    extern /* Subroutine */ int slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    slartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), slargv_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), slartv_(
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);


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

#line 203 "ssbgst.f"
    /* Parameter adjustments */
#line 203 "ssbgst.f"
    ab_dim1 = *ldab;
#line 203 "ssbgst.f"
    ab_offset = 1 + ab_dim1;
#line 203 "ssbgst.f"
    ab -= ab_offset;
#line 203 "ssbgst.f"
    bb_dim1 = *ldbb;
#line 203 "ssbgst.f"
    bb_offset = 1 + bb_dim1;
#line 203 "ssbgst.f"
    bb -= bb_offset;
#line 203 "ssbgst.f"
    x_dim1 = *ldx;
#line 203 "ssbgst.f"
    x_offset = 1 + x_dim1;
#line 203 "ssbgst.f"
    x -= x_offset;
#line 203 "ssbgst.f"
    --work;
#line 203 "ssbgst.f"

#line 203 "ssbgst.f"
    /* Function Body */
#line 203 "ssbgst.f"
    wantx = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 204 "ssbgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 205 "ssbgst.f"
    ka1 = *ka + 1;
#line 206 "ssbgst.f"
    kb1 = *kb + 1;
#line 207 "ssbgst.f"
    *info = 0;
#line 208 "ssbgst.f"
    if (! wantx && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 209 "ssbgst.f"
	*info = -1;
#line 210 "ssbgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 211 "ssbgst.f"
	*info = -2;
#line 212 "ssbgst.f"
    } else if (*n < 0) {
#line 213 "ssbgst.f"
	*info = -3;
#line 214 "ssbgst.f"
    } else if (*ka < 0) {
#line 215 "ssbgst.f"
	*info = -4;
#line 216 "ssbgst.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 217 "ssbgst.f"
	*info = -5;
#line 218 "ssbgst.f"
    } else if (*ldab < *ka + 1) {
#line 219 "ssbgst.f"
	*info = -7;
#line 220 "ssbgst.f"
    } else if (*ldbb < *kb + 1) {
#line 221 "ssbgst.f"
	*info = -9;
#line 222 "ssbgst.f"
    } else if (*ldx < 1 || wantx && *ldx < max(1,*n)) {
#line 223 "ssbgst.f"
	*info = -11;
#line 224 "ssbgst.f"
    }
#line 225 "ssbgst.f"
    if (*info != 0) {
#line 226 "ssbgst.f"
	i__1 = -(*info);
#line 226 "ssbgst.f"
	xerbla_("SSBGST", &i__1, (ftnlen)6);
#line 227 "ssbgst.f"
	return 0;
#line 228 "ssbgst.f"
    }

/*     Quick return if possible */

#line 232 "ssbgst.f"
    if (*n == 0) {
#line 232 "ssbgst.f"
	return 0;
#line 232 "ssbgst.f"
    }

#line 235 "ssbgst.f"
    inca = *ldab * ka1;

/*     Initialize X to the unit matrix, if needed */

#line 239 "ssbgst.f"
    if (wantx) {
#line 239 "ssbgst.f"
	slaset_("Full", n, n, &c_b8, &c_b9, &x[x_offset], ldx, (ftnlen)4);
#line 239 "ssbgst.f"
    }

/*     Set M to the splitting point m. It must be the same value as is */
/*     used in SPBSTF. The chosen value allows the arrays WORK and RWORK */
/*     to be of dimension (N). */

#line 246 "ssbgst.f"
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

#line 309 "ssbgst.f"
    update = TRUE_;
#line 310 "ssbgst.f"
    i__ = *n + 1;
#line 311 "ssbgst.f"
L10:
#line 312 "ssbgst.f"
    if (update) {
#line 313 "ssbgst.f"
	--i__;
/* Computing MIN */
#line 314 "ssbgst.f"
	i__1 = *kb, i__2 = i__ - 1;
#line 314 "ssbgst.f"
	kbt = min(i__1,i__2);
#line 315 "ssbgst.f"
	i0 = i__ - 1;
/* Computing MIN */
#line 316 "ssbgst.f"
	i__1 = *n, i__2 = i__ + *ka;
#line 316 "ssbgst.f"
	i1 = min(i__1,i__2);
#line 317 "ssbgst.f"
	i2 = i__ - kbt + ka1;
#line 318 "ssbgst.f"
	if (i__ < m + 1) {
#line 319 "ssbgst.f"
	    update = FALSE_;
#line 320 "ssbgst.f"
	    ++i__;
#line 321 "ssbgst.f"
	    i0 = m;
#line 322 "ssbgst.f"
	    if (*ka == 0) {
#line 322 "ssbgst.f"
		goto L480;
#line 322 "ssbgst.f"
	    }
#line 324 "ssbgst.f"
	    goto L10;
#line 325 "ssbgst.f"
	}
#line 326 "ssbgst.f"
    } else {
#line 327 "ssbgst.f"
	i__ += *ka;
#line 328 "ssbgst.f"
	if (i__ > *n - 1) {
#line 328 "ssbgst.f"
	    goto L480;
#line 328 "ssbgst.f"
	}
#line 330 "ssbgst.f"
    }

#line 332 "ssbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 336 "ssbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 340 "ssbgst.f"
	    bii = bb[kb1 + i__ * bb_dim1];
#line 341 "ssbgst.f"
	    i__1 = i1;
#line 341 "ssbgst.f"
	    for (j = i__; j <= i__1; ++j) {
#line 342 "ssbgst.f"
		ab[i__ - j + ka1 + j * ab_dim1] /= bii;
#line 343 "ssbgst.f"
/* L20: */
#line 343 "ssbgst.f"
	    }
/* Computing MAX */
#line 344 "ssbgst.f"
	    i__1 = 1, i__2 = i__ - *ka;
#line 344 "ssbgst.f"
	    i__3 = i__;
#line 344 "ssbgst.f"
	    for (j = max(i__1,i__2); j <= i__3; ++j) {
#line 345 "ssbgst.f"
		ab[j - i__ + ka1 + i__ * ab_dim1] /= bii;
#line 346 "ssbgst.f"
/* L30: */
#line 346 "ssbgst.f"
	    }
#line 347 "ssbgst.f"
	    i__3 = i__ - 1;
#line 347 "ssbgst.f"
	    for (k = i__ - kbt; k <= i__3; ++k) {
#line 348 "ssbgst.f"
		i__1 = k;
#line 348 "ssbgst.f"
		for (j = i__ - kbt; j <= i__1; ++j) {
#line 349 "ssbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] = ab[j - k + ka1 + k * 
			    ab_dim1] - bb[j - i__ + kb1 + i__ * bb_dim1] * ab[
			    k - i__ + ka1 + i__ * ab_dim1] - bb[k - i__ + kb1 
			    + i__ * bb_dim1] * ab[j - i__ + ka1 + i__ * 
			    ab_dim1] + ab[ka1 + i__ * ab_dim1] * bb[j - i__ + 
			    kb1 + i__ * bb_dim1] * bb[k - i__ + kb1 + i__ * 
			    bb_dim1];
#line 354 "ssbgst.f"
/* L40: */
#line 354 "ssbgst.f"
		}
/* Computing MAX */
#line 355 "ssbgst.f"
		i__1 = 1, i__2 = i__ - *ka;
#line 355 "ssbgst.f"
		i__4 = i__ - kbt - 1;
#line 355 "ssbgst.f"
		for (j = max(i__1,i__2); j <= i__4; ++j) {
#line 356 "ssbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] -= bb[k - i__ + kb1 + i__ * 
			    bb_dim1] * ab[j - i__ + ka1 + i__ * ab_dim1];
#line 358 "ssbgst.f"
/* L50: */
#line 358 "ssbgst.f"
		}
#line 359 "ssbgst.f"
/* L60: */
#line 359 "ssbgst.f"
	    }
#line 360 "ssbgst.f"
	    i__3 = i1;
#line 360 "ssbgst.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MAX */
#line 361 "ssbgst.f"
		i__4 = j - *ka, i__1 = i__ - kbt;
#line 361 "ssbgst.f"
		i__2 = i__ - 1;
#line 361 "ssbgst.f"
		for (k = max(i__4,i__1); k <= i__2; ++k) {
#line 362 "ssbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] -= bb[k - i__ + kb1 + i__ * 
			    bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1];
#line 364 "ssbgst.f"
/* L70: */
#line 364 "ssbgst.f"
		}
#line 365 "ssbgst.f"
/* L80: */
#line 365 "ssbgst.f"
	    }

#line 367 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 371 "ssbgst.f"
		i__3 = *n - m;
#line 371 "ssbgst.f"
		d__1 = 1. / bii;
#line 371 "ssbgst.f"
		sscal_(&i__3, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 372 "ssbgst.f"
		if (kbt > 0) {
#line 372 "ssbgst.f"
		    i__3 = *n - m;
#line 372 "ssbgst.f"
		    sger_(&i__3, &kbt, &c_b20, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kb1 - kbt + i__ * bb_dim1], &c__1, &x[m 
			    + 1 + (i__ - kbt) * x_dim1], ldx);
#line 372 "ssbgst.f"
		}
#line 375 "ssbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 379 "ssbgst.f"
	    ra1 = ab[i__ - i1 + ka1 + i1 * ab_dim1];
#line 380 "ssbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 386 "ssbgst.f"
	i__3 = *kb - 1;
#line 386 "ssbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 387 "ssbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 392 "ssbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i,i-k+ka+1) */

#line 396 "ssbgst.f"
		    slartg_(&ab[k + 1 + (i__ - k + *ka) * ab_dim1], &ra1, &
			    work[*n + i__ - k + *ka - m], &work[i__ - k + *ka 
			    - m], &ra);

/*                 create nonzero element a(i-k,i-k+ka+1) outside the */
/*                 band and store it in WORK(i-k) */

#line 403 "ssbgst.f"
		    t = -bb[kb1 - k + i__ * bb_dim1] * ra1;
#line 404 "ssbgst.f"
		    work[i__ - k] = work[*n + i__ - k + *ka - m] * t - work[
			    i__ - k + *ka - m] * ab[(i__ - k + *ka) * ab_dim1 
			    + 1];
#line 406 "ssbgst.f"
		    ab[(i__ - k + *ka) * ab_dim1 + 1] = work[i__ - k + *ka - 
			    m] * t + work[*n + i__ - k + *ka - m] * ab[(i__ - 
			    k + *ka) * ab_dim1 + 1];
#line 408 "ssbgst.f"
		    ra1 = ra;
#line 409 "ssbgst.f"
		}
#line 410 "ssbgst.f"
	    }
/* Computing MAX */
#line 411 "ssbgst.f"
	    i__2 = 1, i__4 = k - i0 + 2;
#line 411 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__2,i__4) * ka1;
#line 412 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 413 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 414 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 415 "ssbgst.f"
		i__2 = j2, i__4 = i__ + (*ka << 1) - k + 1;
#line 415 "ssbgst.f"
		j2t = max(i__2,i__4);
#line 416 "ssbgst.f"
	    } else {
#line 417 "ssbgst.f"
		j2t = j2;
#line 418 "ssbgst.f"
	    }
#line 419 "ssbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 420 "ssbgst.f"
	    i__2 = j1;
#line 420 "ssbgst.f"
	    i__4 = ka1;
#line 420 "ssbgst.f"
	    for (j = j2t; i__4 < 0 ? j >= i__2 : j <= i__2; j += i__4) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j-m) */

#line 425 "ssbgst.f"
		work[j - m] *= ab[(j + 1) * ab_dim1 + 1];
#line 426 "ssbgst.f"
		ab[(j + 1) * ab_dim1 + 1] = work[*n + j - m] * ab[(j + 1) * 
			ab_dim1 + 1];
#line 427 "ssbgst.f"
/* L90: */
#line 427 "ssbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 432 "ssbgst.f"
	    if (nrt > 0) {
#line 432 "ssbgst.f"
		slargv_(&nrt, &ab[j2t * ab_dim1 + 1], &inca, &work[j2t - m], &
			ka1, &work[*n + j2t - m], &ka1);
#line 432 "ssbgst.f"
	    }
#line 435 "ssbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 439 "ssbgst.f"
		i__4 = *ka - 1;
#line 439 "ssbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 440 "ssbgst.f"
		    slartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &work[*n + j2 - 
			    m], &work[j2 - m], &ka1);
#line 443 "ssbgst.f"
/* L100: */
#line 443 "ssbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 448 "ssbgst.f"
		slar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &work[
			*n + j2 - m], &work[j2 - m], &ka1);

#line 452 "ssbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 456 "ssbgst.f"
	    i__4 = *kb - k + 1;
#line 456 "ssbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 457 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 458 "ssbgst.f"
		if (nrt > 0) {
#line 458 "ssbgst.f"
		    slartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2 - m], &work[j2 - m], &ka1);
#line 458 "ssbgst.f"
		}
#line 462 "ssbgst.f"
/* L110: */
#line 462 "ssbgst.f"
	    }

#line 464 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 468 "ssbgst.f"
		i__4 = j1;
#line 468 "ssbgst.f"
		i__2 = ka1;
#line 468 "ssbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__4 : j <= i__4; j += i__2) {
#line 469 "ssbgst.f"
		    i__1 = *n - m;
#line 469 "ssbgst.f"
		    srot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j - m], &work[j 
			    - m]);
#line 471 "ssbgst.f"
/* L120: */
#line 471 "ssbgst.f"
		}
#line 472 "ssbgst.f"
	    }
#line 473 "ssbgst.f"
/* L130: */
#line 473 "ssbgst.f"
	}

#line 475 "ssbgst.f"
	if (update) {
#line 476 "ssbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt,i-kbt+ka+1) outside the */
/*              band and store it in WORK(i-kbt) */

#line 481 "ssbgst.f"
		work[i__ - kbt] = -bb[kb1 - kbt + i__ * bb_dim1] * ra1;
#line 482 "ssbgst.f"
	    }
#line 483 "ssbgst.f"
	}

#line 485 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 486 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 487 "ssbgst.f"
		i__3 = 2, i__2 = k - i0 + 1;
#line 487 "ssbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 488 "ssbgst.f"
	    } else {
/* Computing MAX */
#line 489 "ssbgst.f"
		i__3 = 1, i__2 = k - i0 + 1;
#line 489 "ssbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 490 "ssbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 494 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 495 "ssbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 496 "ssbgst.f"
		if (nrt > 0) {
#line 496 "ssbgst.f"
		    slartv_(&nrt, &ab[l + (j2 - l + 1) * ab_dim1], &inca, &ab[
			    l + 1 + (j2 - l + 1) * ab_dim1], &inca, &work[*n 
			    + j2 - *ka], &work[j2 - *ka], &ka1);
#line 496 "ssbgst.f"
		}
#line 500 "ssbgst.f"
/* L140: */
#line 500 "ssbgst.f"
	    }
#line 501 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 502 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 503 "ssbgst.f"
	    i__3 = j2;
#line 503 "ssbgst.f"
	    i__2 = -ka1;
#line 503 "ssbgst.f"
	    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 504 "ssbgst.f"
		work[j] = work[j - *ka];
#line 505 "ssbgst.f"
		work[*n + j] = work[*n + j - *ka];
#line 506 "ssbgst.f"
/* L150: */
#line 506 "ssbgst.f"
	    }
#line 507 "ssbgst.f"
	    i__2 = j1;
#line 507 "ssbgst.f"
	    i__3 = ka1;
#line 507 "ssbgst.f"
	    for (j = j2; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j) */

#line 512 "ssbgst.f"
		work[j] *= ab[(j + 1) * ab_dim1 + 1];
#line 513 "ssbgst.f"
		ab[(j + 1) * ab_dim1 + 1] = work[*n + j] * ab[(j + 1) * 
			ab_dim1 + 1];
#line 514 "ssbgst.f"
/* L160: */
#line 514 "ssbgst.f"
	    }
#line 515 "ssbgst.f"
	    if (update) {
#line 516 "ssbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 516 "ssbgst.f"
		    work[i__ - k + *ka] = work[i__ - k];
#line 516 "ssbgst.f"
		}
#line 518 "ssbgst.f"
	    }
#line 519 "ssbgst.f"
/* L170: */
#line 519 "ssbgst.f"
	}

#line 521 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 522 "ssbgst.f"
	    i__3 = 1, i__2 = k - i0 + 1;
#line 522 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 523 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 524 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 525 "ssbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 530 "ssbgst.f"
		slargv_(&nr, &ab[j2 * ab_dim1 + 1], &inca, &work[j2], &ka1, &
			work[*n + j2], &ka1);

/*              apply rotations in 2nd set from the right */

#line 535 "ssbgst.f"
		i__3 = *ka - 1;
#line 535 "ssbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 536 "ssbgst.f"
		    slartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &work[*n + j2], 
			    &work[j2], &ka1);
#line 539 "ssbgst.f"
/* L180: */
#line 539 "ssbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 544 "ssbgst.f"
		slar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &work[
			*n + j2], &work[j2], &ka1);

#line 548 "ssbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 552 "ssbgst.f"
	    i__3 = *kb - k + 1;
#line 552 "ssbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 553 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 554 "ssbgst.f"
		if (nrt > 0) {
#line 554 "ssbgst.f"
		    slartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2], &work[j2], &ka1);
#line 554 "ssbgst.f"
		}
#line 558 "ssbgst.f"
/* L190: */
#line 558 "ssbgst.f"
	    }

#line 560 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 564 "ssbgst.f"
		i__3 = j1;
#line 564 "ssbgst.f"
		i__2 = ka1;
#line 564 "ssbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 565 "ssbgst.f"
		    i__4 = *n - m;
#line 565 "ssbgst.f"
		    srot_(&i__4, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j], &work[j]);
#line 567 "ssbgst.f"
/* L200: */
#line 567 "ssbgst.f"
		}
#line 568 "ssbgst.f"
	    }
#line 569 "ssbgst.f"
/* L210: */
#line 569 "ssbgst.f"
	}

#line 571 "ssbgst.f"
	i__2 = *kb - 1;
#line 571 "ssbgst.f"
	for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 572 "ssbgst.f"
	    i__3 = 1, i__4 = k - i0 + 2;
#line 572 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__4) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 576 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 577 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 578 "ssbgst.f"
		if (nrt > 0) {
#line 578 "ssbgst.f"
		    slartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    work[*n + j2 - m], &work[j2 - m], &ka1);
#line 578 "ssbgst.f"
		}
#line 582 "ssbgst.f"
/* L220: */
#line 582 "ssbgst.f"
	    }
#line 583 "ssbgst.f"
/* L230: */
#line 583 "ssbgst.f"
	}

#line 585 "ssbgst.f"
	if (*kb > 1) {
#line 586 "ssbgst.f"
	    i__2 = i__ - *kb + (*ka << 1) + 1;
#line 586 "ssbgst.f"
	    for (j = *n - 1; j >= i__2; --j) {
#line 587 "ssbgst.f"
		work[*n + j - m] = work[*n + j - *ka - m];
#line 588 "ssbgst.f"
		work[j - m] = work[j - *ka - m];
#line 589 "ssbgst.f"
/* L240: */
#line 589 "ssbgst.f"
	    }
#line 590 "ssbgst.f"
	}

#line 592 "ssbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 596 "ssbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 600 "ssbgst.f"
	    bii = bb[i__ * bb_dim1 + 1];
#line 601 "ssbgst.f"
	    i__2 = i1;
#line 601 "ssbgst.f"
	    for (j = i__; j <= i__2; ++j) {
#line 602 "ssbgst.f"
		ab[j - i__ + 1 + i__ * ab_dim1] /= bii;
#line 603 "ssbgst.f"
/* L250: */
#line 603 "ssbgst.f"
	    }
/* Computing MAX */
#line 604 "ssbgst.f"
	    i__2 = 1, i__3 = i__ - *ka;
#line 604 "ssbgst.f"
	    i__4 = i__;
#line 604 "ssbgst.f"
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 605 "ssbgst.f"
		ab[i__ - j + 1 + j * ab_dim1] /= bii;
#line 606 "ssbgst.f"
/* L260: */
#line 606 "ssbgst.f"
	    }
#line 607 "ssbgst.f"
	    i__4 = i__ - 1;
#line 607 "ssbgst.f"
	    for (k = i__ - kbt; k <= i__4; ++k) {
#line 608 "ssbgst.f"
		i__2 = k;
#line 608 "ssbgst.f"
		for (j = i__ - kbt; j <= i__2; ++j) {
#line 609 "ssbgst.f"
		    ab[k - j + 1 + j * ab_dim1] = ab[k - j + 1 + j * ab_dim1] 
			    - bb[i__ - j + 1 + j * bb_dim1] * ab[i__ - k + 1 
			    + k * ab_dim1] - bb[i__ - k + 1 + k * bb_dim1] * 
			    ab[i__ - j + 1 + j * ab_dim1] + ab[i__ * ab_dim1 
			    + 1] * bb[i__ - j + 1 + j * bb_dim1] * bb[i__ - k 
			    + 1 + k * bb_dim1];
#line 614 "ssbgst.f"
/* L270: */
#line 614 "ssbgst.f"
		}
/* Computing MAX */
#line 615 "ssbgst.f"
		i__2 = 1, i__3 = i__ - *ka;
#line 615 "ssbgst.f"
		i__1 = i__ - kbt - 1;
#line 615 "ssbgst.f"
		for (j = max(i__2,i__3); j <= i__1; ++j) {
#line 616 "ssbgst.f"
		    ab[k - j + 1 + j * ab_dim1] -= bb[i__ - k + 1 + k * 
			    bb_dim1] * ab[i__ - j + 1 + j * ab_dim1];
#line 618 "ssbgst.f"
/* L280: */
#line 618 "ssbgst.f"
		}
#line 619 "ssbgst.f"
/* L290: */
#line 619 "ssbgst.f"
	    }
#line 620 "ssbgst.f"
	    i__4 = i1;
#line 620 "ssbgst.f"
	    for (j = i__; j <= i__4; ++j) {
/* Computing MAX */
#line 621 "ssbgst.f"
		i__1 = j - *ka, i__2 = i__ - kbt;
#line 621 "ssbgst.f"
		i__3 = i__ - 1;
#line 621 "ssbgst.f"
		for (k = max(i__1,i__2); k <= i__3; ++k) {
#line 622 "ssbgst.f"
		    ab[j - k + 1 + k * ab_dim1] -= bb[i__ - k + 1 + k * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1];
#line 624 "ssbgst.f"
/* L300: */
#line 624 "ssbgst.f"
		}
#line 625 "ssbgst.f"
/* L310: */
#line 625 "ssbgst.f"
	    }

#line 627 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 631 "ssbgst.f"
		i__4 = *n - m;
#line 631 "ssbgst.f"
		d__1 = 1. / bii;
#line 631 "ssbgst.f"
		sscal_(&i__4, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 632 "ssbgst.f"
		if (kbt > 0) {
#line 632 "ssbgst.f"
		    i__4 = *n - m;
#line 632 "ssbgst.f"
		    i__3 = *ldbb - 1;
#line 632 "ssbgst.f"
		    sger_(&i__4, &kbt, &c_b20, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kbt + 1 + (i__ - kbt) * bb_dim1], &i__3,
			     &x[m + 1 + (i__ - kbt) * x_dim1], ldx);
#line 632 "ssbgst.f"
		}
#line 636 "ssbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 640 "ssbgst.f"
	    ra1 = ab[i1 - i__ + 1 + i__ * ab_dim1];
#line 641 "ssbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 647 "ssbgst.f"
	i__4 = *kb - 1;
#line 647 "ssbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 648 "ssbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 653 "ssbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i-k+ka+1,i) */

#line 657 "ssbgst.f"
		    slartg_(&ab[ka1 - k + i__ * ab_dim1], &ra1, &work[*n + 
			    i__ - k + *ka - m], &work[i__ - k + *ka - m], &ra)
			    ;

/*                 create nonzero element a(i-k+ka+1,i-k) outside the */
/*                 band and store it in WORK(i-k) */

#line 663 "ssbgst.f"
		    t = -bb[k + 1 + (i__ - k) * bb_dim1] * ra1;
#line 664 "ssbgst.f"
		    work[i__ - k] = work[*n + i__ - k + *ka - m] * t - work[
			    i__ - k + *ka - m] * ab[ka1 + (i__ - k) * ab_dim1]
			    ;
#line 666 "ssbgst.f"
		    ab[ka1 + (i__ - k) * ab_dim1] = work[i__ - k + *ka - m] * 
			    t + work[*n + i__ - k + *ka - m] * ab[ka1 + (i__ 
			    - k) * ab_dim1];
#line 668 "ssbgst.f"
		    ra1 = ra;
#line 669 "ssbgst.f"
		}
#line 670 "ssbgst.f"
	    }
/* Computing MAX */
#line 671 "ssbgst.f"
	    i__3 = 1, i__1 = k - i0 + 2;
#line 671 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__1) * ka1;
#line 672 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 673 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 674 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 675 "ssbgst.f"
		i__3 = j2, i__1 = i__ + (*ka << 1) - k + 1;
#line 675 "ssbgst.f"
		j2t = max(i__3,i__1);
#line 676 "ssbgst.f"
	    } else {
#line 677 "ssbgst.f"
		j2t = j2;
#line 678 "ssbgst.f"
	    }
#line 679 "ssbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 680 "ssbgst.f"
	    i__3 = j1;
#line 680 "ssbgst.f"
	    i__1 = ka1;
#line 680 "ssbgst.f"
	    for (j = j2t; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j-m) */

#line 685 "ssbgst.f"
		work[j - m] *= ab[ka1 + (j - *ka + 1) * ab_dim1];
#line 686 "ssbgst.f"
		ab[ka1 + (j - *ka + 1) * ab_dim1] = work[*n + j - m] * ab[ka1 
			+ (j - *ka + 1) * ab_dim1];
#line 687 "ssbgst.f"
/* L320: */
#line 687 "ssbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 692 "ssbgst.f"
	    if (nrt > 0) {
#line 692 "ssbgst.f"
		slargv_(&nrt, &ab[ka1 + (j2t - *ka) * ab_dim1], &inca, &work[
			j2t - m], &ka1, &work[*n + j2t - m], &ka1);
#line 692 "ssbgst.f"
	    }
#line 695 "ssbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 699 "ssbgst.f"
		i__1 = *ka - 1;
#line 699 "ssbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 700 "ssbgst.f"
		    slartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &work[*n + j2 
			    - m], &work[j2 - m], &ka1);
#line 703 "ssbgst.f"
/* L330: */
#line 703 "ssbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 708 "ssbgst.f"
		slar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &work[*n + j2 - m], 
			&work[j2 - m], &ka1);

#line 711 "ssbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 715 "ssbgst.f"
	    i__1 = *kb - k + 1;
#line 715 "ssbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 716 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 717 "ssbgst.f"
		if (nrt > 0) {
#line 717 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2 - m], &work[j2 - m], &ka1);
#line 717 "ssbgst.f"
		}
#line 721 "ssbgst.f"
/* L340: */
#line 721 "ssbgst.f"
	    }

#line 723 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 727 "ssbgst.f"
		i__1 = j1;
#line 727 "ssbgst.f"
		i__3 = ka1;
#line 727 "ssbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 728 "ssbgst.f"
		    i__2 = *n - m;
#line 728 "ssbgst.f"
		    srot_(&i__2, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j - m], &work[j 
			    - m]);
#line 730 "ssbgst.f"
/* L350: */
#line 730 "ssbgst.f"
		}
#line 731 "ssbgst.f"
	    }
#line 732 "ssbgst.f"
/* L360: */
#line 732 "ssbgst.f"
	}

#line 734 "ssbgst.f"
	if (update) {
#line 735 "ssbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt+ka+1,i-kbt) outside the */
/*              band and store it in WORK(i-kbt) */

#line 740 "ssbgst.f"
		work[i__ - kbt] = -bb[kbt + 1 + (i__ - kbt) * bb_dim1] * ra1;
#line 741 "ssbgst.f"
	    }
#line 742 "ssbgst.f"
	}

#line 744 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 745 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 746 "ssbgst.f"
		i__4 = 2, i__3 = k - i0 + 1;
#line 746 "ssbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 747 "ssbgst.f"
	    } else {
/* Computing MAX */
#line 748 "ssbgst.f"
		i__4 = 1, i__3 = k - i0 + 1;
#line 748 "ssbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 749 "ssbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 753 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 754 "ssbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 755 "ssbgst.f"
		if (nrt > 0) {
#line 755 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + (j2 - *ka) * ab_dim1], &
			    inca, &ab[ka1 - l + (j2 - *ka + 1) * ab_dim1], &
			    inca, &work[*n + j2 - *ka], &work[j2 - *ka], &ka1)
			    ;
#line 755 "ssbgst.f"
		}
#line 759 "ssbgst.f"
/* L370: */
#line 759 "ssbgst.f"
	    }
#line 760 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 761 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 762 "ssbgst.f"
	    i__4 = j2;
#line 762 "ssbgst.f"
	    i__3 = -ka1;
#line 762 "ssbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 763 "ssbgst.f"
		work[j] = work[j - *ka];
#line 764 "ssbgst.f"
		work[*n + j] = work[*n + j - *ka];
#line 765 "ssbgst.f"
/* L380: */
#line 765 "ssbgst.f"
	    }
#line 766 "ssbgst.f"
	    i__3 = j1;
#line 766 "ssbgst.f"
	    i__4 = ka1;
#line 766 "ssbgst.f"
	    for (j = j2; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j) */

#line 771 "ssbgst.f"
		work[j] *= ab[ka1 + (j - *ka + 1) * ab_dim1];
#line 772 "ssbgst.f"
		ab[ka1 + (j - *ka + 1) * ab_dim1] = work[*n + j] * ab[ka1 + (
			j - *ka + 1) * ab_dim1];
#line 773 "ssbgst.f"
/* L390: */
#line 773 "ssbgst.f"
	    }
#line 774 "ssbgst.f"
	    if (update) {
#line 775 "ssbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 775 "ssbgst.f"
		    work[i__ - k + *ka] = work[i__ - k];
#line 775 "ssbgst.f"
		}
#line 777 "ssbgst.f"
	    }
#line 778 "ssbgst.f"
/* L400: */
#line 778 "ssbgst.f"
	}

#line 780 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 781 "ssbgst.f"
	    i__4 = 1, i__3 = k - i0 + 1;
#line 781 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 782 "ssbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 783 "ssbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 784 "ssbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 789 "ssbgst.f"
		slargv_(&nr, &ab[ka1 + (j2 - *ka) * ab_dim1], &inca, &work[j2]
			, &ka1, &work[*n + j2], &ka1);

/*              apply rotations in 2nd set from the left */

#line 794 "ssbgst.f"
		i__4 = *ka - 1;
#line 794 "ssbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 795 "ssbgst.f"
		    slartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &work[*n + j2]
			    , &work[j2], &ka1);
#line 798 "ssbgst.f"
/* L410: */
#line 798 "ssbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 803 "ssbgst.f"
		slar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &work[*n + j2], &
			work[j2], &ka1);

#line 806 "ssbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 810 "ssbgst.f"
	    i__4 = *kb - k + 1;
#line 810 "ssbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 811 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 812 "ssbgst.f"
		if (nrt > 0) {
#line 812 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2], &work[j2], &ka1);
#line 812 "ssbgst.f"
		}
#line 816 "ssbgst.f"
/* L420: */
#line 816 "ssbgst.f"
	    }

#line 818 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 822 "ssbgst.f"
		i__4 = j1;
#line 822 "ssbgst.f"
		i__3 = ka1;
#line 822 "ssbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 823 "ssbgst.f"
		    i__1 = *n - m;
#line 823 "ssbgst.f"
		    srot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &work[*n + j], &work[j]);
#line 825 "ssbgst.f"
/* L430: */
#line 825 "ssbgst.f"
		}
#line 826 "ssbgst.f"
	    }
#line 827 "ssbgst.f"
/* L440: */
#line 827 "ssbgst.f"
	}

#line 829 "ssbgst.f"
	i__3 = *kb - 1;
#line 829 "ssbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 830 "ssbgst.f"
	    i__4 = 1, i__1 = k - i0 + 2;
#line 830 "ssbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 834 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 835 "ssbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 836 "ssbgst.f"
		if (nrt > 0) {
#line 836 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &work[*n + 
			    j2 - m], &work[j2 - m], &ka1);
#line 836 "ssbgst.f"
		}
#line 840 "ssbgst.f"
/* L450: */
#line 840 "ssbgst.f"
	    }
#line 841 "ssbgst.f"
/* L460: */
#line 841 "ssbgst.f"
	}

#line 843 "ssbgst.f"
	if (*kb > 1) {
#line 844 "ssbgst.f"
	    i__3 = i__ - *kb + (*ka << 1) + 1;
#line 844 "ssbgst.f"
	    for (j = *n - 1; j >= i__3; --j) {
#line 845 "ssbgst.f"
		work[*n + j - m] = work[*n + j - *ka - m];
#line 846 "ssbgst.f"
		work[j - m] = work[j - *ka - m];
#line 847 "ssbgst.f"
/* L470: */
#line 847 "ssbgst.f"
	    }
#line 848 "ssbgst.f"
	}

#line 850 "ssbgst.f"
    }

#line 852 "ssbgst.f"
    goto L10;

#line 854 "ssbgst.f"
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

#line 872 "ssbgst.f"
    update = TRUE_;
#line 873 "ssbgst.f"
    i__ = 0;
#line 874 "ssbgst.f"
L490:
#line 875 "ssbgst.f"
    if (update) {
#line 876 "ssbgst.f"
	++i__;
/* Computing MIN */
#line 877 "ssbgst.f"
	i__3 = *kb, i__4 = m - i__;
#line 877 "ssbgst.f"
	kbt = min(i__3,i__4);
#line 878 "ssbgst.f"
	i0 = i__ + 1;
/* Computing MAX */
#line 879 "ssbgst.f"
	i__3 = 1, i__4 = i__ - *ka;
#line 879 "ssbgst.f"
	i1 = max(i__3,i__4);
#line 880 "ssbgst.f"
	i2 = i__ + kbt - ka1;
#line 881 "ssbgst.f"
	if (i__ > m) {
#line 882 "ssbgst.f"
	    update = FALSE_;
#line 883 "ssbgst.f"
	    --i__;
#line 884 "ssbgst.f"
	    i0 = m + 1;
#line 885 "ssbgst.f"
	    if (*ka == 0) {
#line 885 "ssbgst.f"
		return 0;
#line 885 "ssbgst.f"
	    }
#line 887 "ssbgst.f"
	    goto L490;
#line 888 "ssbgst.f"
	}
#line 889 "ssbgst.f"
    } else {
#line 890 "ssbgst.f"
	i__ -= *ka;
#line 891 "ssbgst.f"
	if (i__ < 2) {
#line 891 "ssbgst.f"
	    return 0;
#line 891 "ssbgst.f"
	}
#line 893 "ssbgst.f"
    }

#line 895 "ssbgst.f"
    if (i__ < m - kbt) {
#line 896 "ssbgst.f"
	nx = m;
#line 897 "ssbgst.f"
    } else {
#line 898 "ssbgst.f"
	nx = *n;
#line 899 "ssbgst.f"
    }

#line 901 "ssbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 905 "ssbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 909 "ssbgst.f"
	    bii = bb[kb1 + i__ * bb_dim1];
#line 910 "ssbgst.f"
	    i__3 = i__;
#line 910 "ssbgst.f"
	    for (j = i1; j <= i__3; ++j) {
#line 911 "ssbgst.f"
		ab[j - i__ + ka1 + i__ * ab_dim1] /= bii;
#line 912 "ssbgst.f"
/* L500: */
#line 912 "ssbgst.f"
	    }
/* Computing MIN */
#line 913 "ssbgst.f"
	    i__4 = *n, i__1 = i__ + *ka;
#line 913 "ssbgst.f"
	    i__3 = min(i__4,i__1);
#line 913 "ssbgst.f"
	    for (j = i__; j <= i__3; ++j) {
#line 914 "ssbgst.f"
		ab[i__ - j + ka1 + j * ab_dim1] /= bii;
#line 915 "ssbgst.f"
/* L510: */
#line 915 "ssbgst.f"
	    }
#line 916 "ssbgst.f"
	    i__3 = i__ + kbt;
#line 916 "ssbgst.f"
	    for (k = i__ + 1; k <= i__3; ++k) {
#line 917 "ssbgst.f"
		i__4 = i__ + kbt;
#line 917 "ssbgst.f"
		for (j = k; j <= i__4; ++j) {
#line 918 "ssbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] = ab[k - j + ka1 + j * 
			    ab_dim1] - bb[i__ - j + kb1 + j * bb_dim1] * ab[
			    i__ - k + ka1 + k * ab_dim1] - bb[i__ - k + kb1 + 
			    k * bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1] + 
			    ab[ka1 + i__ * ab_dim1] * bb[i__ - j + kb1 + j * 
			    bb_dim1] * bb[i__ - k + kb1 + k * bb_dim1];
#line 923 "ssbgst.f"
/* L520: */
#line 923 "ssbgst.f"
		}
/* Computing MIN */
#line 924 "ssbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 924 "ssbgst.f"
		i__4 = min(i__1,i__2);
#line 924 "ssbgst.f"
		for (j = i__ + kbt + 1; j <= i__4; ++j) {
#line 925 "ssbgst.f"
		    ab[k - j + ka1 + j * ab_dim1] -= bb[i__ - k + kb1 + k * 
			    bb_dim1] * ab[i__ - j + ka1 + j * ab_dim1];
#line 927 "ssbgst.f"
/* L530: */
#line 927 "ssbgst.f"
		}
#line 928 "ssbgst.f"
/* L540: */
#line 928 "ssbgst.f"
	    }
#line 929 "ssbgst.f"
	    i__3 = i__;
#line 929 "ssbgst.f"
	    for (j = i1; j <= i__3; ++j) {
/* Computing MIN */
#line 930 "ssbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 930 "ssbgst.f"
		i__4 = min(i__1,i__2);
#line 930 "ssbgst.f"
		for (k = i__ + 1; k <= i__4; ++k) {
#line 931 "ssbgst.f"
		    ab[j - k + ka1 + k * ab_dim1] -= bb[i__ - k + kb1 + k * 
			    bb_dim1] * ab[j - i__ + ka1 + i__ * ab_dim1];
#line 933 "ssbgst.f"
/* L550: */
#line 933 "ssbgst.f"
		}
#line 934 "ssbgst.f"
/* L560: */
#line 934 "ssbgst.f"
	    }

#line 936 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 940 "ssbgst.f"
		d__1 = 1. / bii;
#line 940 "ssbgst.f"
		sscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 941 "ssbgst.f"
		if (kbt > 0) {
#line 941 "ssbgst.f"
		    i__3 = *ldbb - 1;
#line 941 "ssbgst.f"
		    sger_(&nx, &kbt, &c_b20, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    *kb + (i__ + 1) * bb_dim1], &i__3, &x[(i__ + 1) * 
			    x_dim1 + 1], ldx);
#line 941 "ssbgst.f"
		}
#line 944 "ssbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 948 "ssbgst.f"
	    ra1 = ab[i1 - i__ + ka1 + i__ * ab_dim1];
#line 949 "ssbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 954 "ssbgst.f"
	i__3 = *kb - 1;
#line 954 "ssbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 955 "ssbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 960 "ssbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i+k-ka-1,i) */

#line 964 "ssbgst.f"
		    slartg_(&ab[k + 1 + i__ * ab_dim1], &ra1, &work[*n + i__ 
			    + k - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k-ka-1,i+k) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 970 "ssbgst.f"
		    t = -bb[kb1 - k + (i__ + k) * bb_dim1] * ra1;
#line 971 "ssbgst.f"
		    work[m - *kb + i__ + k] = work[*n + i__ + k - *ka] * t - 
			    work[i__ + k - *ka] * ab[(i__ + k) * ab_dim1 + 1];
#line 973 "ssbgst.f"
		    ab[(i__ + k) * ab_dim1 + 1] = work[i__ + k - *ka] * t + 
			    work[*n + i__ + k - *ka] * ab[(i__ + k) * ab_dim1 
			    + 1];
#line 975 "ssbgst.f"
		    ra1 = ra;
#line 976 "ssbgst.f"
		}
#line 977 "ssbgst.f"
	    }
/* Computing MAX */
#line 978 "ssbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 978 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;
#line 979 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 980 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 981 "ssbgst.f"
	    if (update) {
/* Computing MIN */
#line 982 "ssbgst.f"
		i__4 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 982 "ssbgst.f"
		j2t = min(i__4,i__1);
#line 983 "ssbgst.f"
	    } else {
#line 984 "ssbgst.f"
		j2t = j2;
#line 985 "ssbgst.f"
	    }
#line 986 "ssbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 987 "ssbgst.f"
	    i__4 = j2t;
#line 987 "ssbgst.f"
	    i__1 = ka1;
#line 987 "ssbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__4 : j <= i__4; j += i__1) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(j) */

#line 992 "ssbgst.f"
		work[j] *= ab[(j + *ka - 1) * ab_dim1 + 1];
#line 993 "ssbgst.f"
		ab[(j + *ka - 1) * ab_dim1 + 1] = work[*n + j] * ab[(j + *ka 
			- 1) * ab_dim1 + 1];
#line 994 "ssbgst.f"
/* L570: */
#line 994 "ssbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 999 "ssbgst.f"
	    if (nrt > 0) {
#line 999 "ssbgst.f"
		slargv_(&nrt, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[j1],
			 &ka1, &work[*n + j1], &ka1);
#line 999 "ssbgst.f"
	    }
#line 1002 "ssbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 1006 "ssbgst.f"
		i__1 = *ka - 1;
#line 1006 "ssbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1007 "ssbgst.f"
		    slartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &work[*n 
			    + j1], &work[j1], &ka1);
#line 1010 "ssbgst.f"
/* L580: */
#line 1010 "ssbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1015 "ssbgst.f"
		slar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &work[*n + 
			j1], &work[j1], &ka1);

#line 1019 "ssbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 1023 "ssbgst.f"
	    i__1 = *kb - k + 1;
#line 1023 "ssbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1024 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1025 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1026 "ssbgst.f"
		if (nrt > 0) {
#line 1026 "ssbgst.f"
		    slartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + j1t], &
			    work[j1t], &ka1);
#line 1026 "ssbgst.f"
		}
#line 1030 "ssbgst.f"
/* L590: */
#line 1030 "ssbgst.f"
	    }

#line 1032 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1036 "ssbgst.f"
		i__1 = j2;
#line 1036 "ssbgst.f"
		i__4 = ka1;
#line 1036 "ssbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__1 : j <= i__1; j += i__4) {
#line 1037 "ssbgst.f"
		    srot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + j], &work[j]);
#line 1039 "ssbgst.f"
/* L600: */
#line 1039 "ssbgst.f"
		}
#line 1040 "ssbgst.f"
	    }
#line 1041 "ssbgst.f"
/* L610: */
#line 1041 "ssbgst.f"
	}

#line 1043 "ssbgst.f"
	if (update) {
#line 1044 "ssbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt-ka-1,i+kbt) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1049 "ssbgst.f"
		work[m - *kb + i__ + kbt] = -bb[kb1 - kbt + (i__ + kbt) * 
			bb_dim1] * ra1;
#line 1050 "ssbgst.f"
	    }
#line 1051 "ssbgst.f"
	}

#line 1053 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1054 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 1055 "ssbgst.f"
		i__3 = 2, i__4 = k + i0 - m;
#line 1055 "ssbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1056 "ssbgst.f"
	    } else {
/* Computing MAX */
#line 1057 "ssbgst.f"
		i__3 = 1, i__4 = k + i0 - m;
#line 1057 "ssbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1058 "ssbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 1062 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1063 "ssbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1064 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1065 "ssbgst.f"
		if (nrt > 0) {
#line 1065 "ssbgst.f"
		    slartv_(&nrt, &ab[l + (j1t + *ka) * ab_dim1], &inca, &ab[
			    l + 1 + (j1t + *ka - 1) * ab_dim1], &inca, &work[*
			    n + m - *kb + j1t + *ka], &work[m - *kb + j1t + *
			    ka], &ka1);
#line 1065 "ssbgst.f"
		}
#line 1070 "ssbgst.f"
/* L620: */
#line 1070 "ssbgst.f"
	    }
#line 1071 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1072 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1073 "ssbgst.f"
	    i__3 = j2;
#line 1073 "ssbgst.f"
	    i__4 = ka1;
#line 1073 "ssbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1074 "ssbgst.f"
		work[m - *kb + j] = work[m - *kb + j + *ka];
#line 1075 "ssbgst.f"
		work[*n + m - *kb + j] = work[*n + m - *kb + j + *ka];
#line 1076 "ssbgst.f"
/* L630: */
#line 1076 "ssbgst.f"
	    }
#line 1077 "ssbgst.f"
	    i__4 = j2;
#line 1077 "ssbgst.f"
	    i__3 = ka1;
#line 1077 "ssbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1082 "ssbgst.f"
		work[m - *kb + j] *= ab[(j + *ka - 1) * ab_dim1 + 1];
#line 1083 "ssbgst.f"
		ab[(j + *ka - 1) * ab_dim1 + 1] = work[*n + m - *kb + j] * ab[
			(j + *ka - 1) * ab_dim1 + 1];
#line 1084 "ssbgst.f"
/* L640: */
#line 1084 "ssbgst.f"
	    }
#line 1085 "ssbgst.f"
	    if (update) {
#line 1086 "ssbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1086 "ssbgst.f"
		    work[m - *kb + i__ + k - *ka] = work[m - *kb + i__ + k];
#line 1086 "ssbgst.f"
		}
#line 1088 "ssbgst.f"
	    }
#line 1089 "ssbgst.f"
/* L650: */
#line 1089 "ssbgst.f"
	}

#line 1091 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1092 "ssbgst.f"
	    i__3 = 1, i__4 = k + i0 - m;
#line 1092 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1093 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1094 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1095 "ssbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1100 "ssbgst.f"
		slargv_(&nr, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[m - *
			kb + j1], &ka1, &work[*n + m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the left */

#line 1105 "ssbgst.f"
		i__3 = *ka - 1;
#line 1105 "ssbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 1106 "ssbgst.f"
		    slartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &work[*n 
			    + m - *kb + j1], &work[m - *kb + j1], &ka1);
#line 1109 "ssbgst.f"
/* L660: */
#line 1109 "ssbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1114 "ssbgst.f"
		slar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &work[*n + 
			m - *kb + j1], &work[m - *kb + j1], &ka1);

#line 1118 "ssbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 1122 "ssbgst.f"
	    i__3 = *kb - k + 1;
#line 1122 "ssbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 1123 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1124 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1125 "ssbgst.f"
		if (nrt > 0) {
#line 1125 "ssbgst.f"
		    slartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + m - *kb + 
			    j1t], &work[m - *kb + j1t], &ka1);
#line 1125 "ssbgst.f"
		}
#line 1130 "ssbgst.f"
/* L670: */
#line 1130 "ssbgst.f"
	    }

#line 1132 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1136 "ssbgst.f"
		i__3 = j2;
#line 1136 "ssbgst.f"
		i__4 = ka1;
#line 1136 "ssbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1137 "ssbgst.f"
		    srot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + m - *kb + j], &work[m - *
			    kb + j]);
#line 1139 "ssbgst.f"
/* L680: */
#line 1139 "ssbgst.f"
		}
#line 1140 "ssbgst.f"
	    }
#line 1141 "ssbgst.f"
/* L690: */
#line 1141 "ssbgst.f"
	}

#line 1143 "ssbgst.f"
	i__4 = *kb - 1;
#line 1143 "ssbgst.f"
	for (k = 1; k <= i__4; ++k) {
/* Computing MAX */
#line 1144 "ssbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1144 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 1148 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1149 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1150 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1151 "ssbgst.f"
		if (nrt > 0) {
#line 1151 "ssbgst.f"
		    slartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &work[*n + j1t], &
			    work[j1t], &ka1);
#line 1151 "ssbgst.f"
		}
#line 1155 "ssbgst.f"
/* L700: */
#line 1155 "ssbgst.f"
	    }
#line 1156 "ssbgst.f"
/* L710: */
#line 1156 "ssbgst.f"
	}

#line 1158 "ssbgst.f"
	if (*kb > 1) {
/* Computing MIN */
#line 1159 "ssbgst.f"
	    i__3 = i__ + *kb;
#line 1159 "ssbgst.f"
	    i__4 = min(i__3,m) - (*ka << 1) - 1;
#line 1159 "ssbgst.f"
	    for (j = 2; j <= i__4; ++j) {
#line 1160 "ssbgst.f"
		work[*n + j] = work[*n + j + *ka];
#line 1161 "ssbgst.f"
		work[j] = work[j + *ka];
#line 1162 "ssbgst.f"
/* L720: */
#line 1162 "ssbgst.f"
	    }
#line 1163 "ssbgst.f"
	}

#line 1165 "ssbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 1169 "ssbgst.f"
	if (update) {

/*           Form  inv(S(i))**T * A * inv(S(i)) */

#line 1173 "ssbgst.f"
	    bii = bb[i__ * bb_dim1 + 1];
#line 1174 "ssbgst.f"
	    i__4 = i__;
#line 1174 "ssbgst.f"
	    for (j = i1; j <= i__4; ++j) {
#line 1175 "ssbgst.f"
		ab[i__ - j + 1 + j * ab_dim1] /= bii;
#line 1176 "ssbgst.f"
/* L730: */
#line 1176 "ssbgst.f"
	    }
/* Computing MIN */
#line 1177 "ssbgst.f"
	    i__3 = *n, i__1 = i__ + *ka;
#line 1177 "ssbgst.f"
	    i__4 = min(i__3,i__1);
#line 1177 "ssbgst.f"
	    for (j = i__; j <= i__4; ++j) {
#line 1178 "ssbgst.f"
		ab[j - i__ + 1 + i__ * ab_dim1] /= bii;
#line 1179 "ssbgst.f"
/* L740: */
#line 1179 "ssbgst.f"
	    }
#line 1180 "ssbgst.f"
	    i__4 = i__ + kbt;
#line 1180 "ssbgst.f"
	    for (k = i__ + 1; k <= i__4; ++k) {
#line 1181 "ssbgst.f"
		i__3 = i__ + kbt;
#line 1181 "ssbgst.f"
		for (j = k; j <= i__3; ++j) {
#line 1182 "ssbgst.f"
		    ab[j - k + 1 + k * ab_dim1] = ab[j - k + 1 + k * ab_dim1] 
			    - bb[j - i__ + 1 + i__ * bb_dim1] * ab[k - i__ + 
			    1 + i__ * ab_dim1] - bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1] + ab[
			    i__ * ab_dim1 + 1] * bb[j - i__ + 1 + i__ * 
			    bb_dim1] * bb[k - i__ + 1 + i__ * bb_dim1];
#line 1187 "ssbgst.f"
/* L750: */
#line 1187 "ssbgst.f"
		}
/* Computing MIN */
#line 1188 "ssbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 1188 "ssbgst.f"
		i__3 = min(i__1,i__2);
#line 1188 "ssbgst.f"
		for (j = i__ + kbt + 1; j <= i__3; ++j) {
#line 1189 "ssbgst.f"
		    ab[j - k + 1 + k * ab_dim1] -= bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[j - i__ + 1 + i__ * ab_dim1];
#line 1191 "ssbgst.f"
/* L760: */
#line 1191 "ssbgst.f"
		}
#line 1192 "ssbgst.f"
/* L770: */
#line 1192 "ssbgst.f"
	    }
#line 1193 "ssbgst.f"
	    i__4 = i__;
#line 1193 "ssbgst.f"
	    for (j = i1; j <= i__4; ++j) {
/* Computing MIN */
#line 1194 "ssbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 1194 "ssbgst.f"
		i__3 = min(i__1,i__2);
#line 1194 "ssbgst.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 1195 "ssbgst.f"
		    ab[k - j + 1 + j * ab_dim1] -= bb[k - i__ + 1 + i__ * 
			    bb_dim1] * ab[i__ - j + 1 + j * ab_dim1];
#line 1197 "ssbgst.f"
/* L780: */
#line 1197 "ssbgst.f"
		}
#line 1198 "ssbgst.f"
/* L790: */
#line 1198 "ssbgst.f"
	    }

#line 1200 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 1204 "ssbgst.f"
		d__1 = 1. / bii;
#line 1204 "ssbgst.f"
		sscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 1205 "ssbgst.f"
		if (kbt > 0) {
#line 1205 "ssbgst.f"
		    sger_(&nx, &kbt, &c_b20, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    i__ * bb_dim1 + 2], &c__1, &x[(i__ + 1) * x_dim1 
			    + 1], ldx);
#line 1205 "ssbgst.f"
		}
#line 1208 "ssbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 1212 "ssbgst.f"
	    ra1 = ab[i__ - i1 + 1 + i1 * ab_dim1];
#line 1213 "ssbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 1218 "ssbgst.f"
	i__4 = *kb - 1;
#line 1218 "ssbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 1219 "ssbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 1224 "ssbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i,i+k-ka-1) */

#line 1228 "ssbgst.f"
		    slartg_(&ab[ka1 - k + (i__ + k - *ka) * ab_dim1], &ra1, &
			    work[*n + i__ + k - *ka], &work[i__ + k - *ka], &
			    ra);

/*                 create nonzero element a(i+k,i+k-ka-1) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 1234 "ssbgst.f"
		    t = -bb[k + 1 + i__ * bb_dim1] * ra1;
#line 1235 "ssbgst.f"
		    work[m - *kb + i__ + k] = work[*n + i__ + k - *ka] * t - 
			    work[i__ + k - *ka] * ab[ka1 + (i__ + k - *ka) * 
			    ab_dim1];
#line 1237 "ssbgst.f"
		    ab[ka1 + (i__ + k - *ka) * ab_dim1] = work[i__ + k - *ka] 
			    * t + work[*n + i__ + k - *ka] * ab[ka1 + (i__ + 
			    k - *ka) * ab_dim1];
#line 1239 "ssbgst.f"
		    ra1 = ra;
#line 1240 "ssbgst.f"
		}
#line 1241 "ssbgst.f"
	    }
/* Computing MAX */
#line 1242 "ssbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1242 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;
#line 1243 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1244 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1245 "ssbgst.f"
	    if (update) {
/* Computing MIN */
#line 1246 "ssbgst.f"
		i__3 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1246 "ssbgst.f"
		j2t = min(i__3,i__1);
#line 1247 "ssbgst.f"
	    } else {
#line 1248 "ssbgst.f"
		j2t = j2;
#line 1249 "ssbgst.f"
	    }
#line 1250 "ssbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1251 "ssbgst.f"
	    i__3 = j2t;
#line 1251 "ssbgst.f"
	    i__1 = ka1;
#line 1251 "ssbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(j) */

#line 1256 "ssbgst.f"
		work[j] *= ab[ka1 + (j - 1) * ab_dim1];
#line 1257 "ssbgst.f"
		ab[ka1 + (j - 1) * ab_dim1] = work[*n + j] * ab[ka1 + (j - 1) 
			* ab_dim1];
#line 1258 "ssbgst.f"
/* L800: */
#line 1258 "ssbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1263 "ssbgst.f"
	    if (nrt > 0) {
#line 1263 "ssbgst.f"
		slargv_(&nrt, &ab[ka1 + j1 * ab_dim1], &inca, &work[j1], &ka1,
			 &work[*n + j1], &ka1);
#line 1263 "ssbgst.f"
	    }
#line 1266 "ssbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 1270 "ssbgst.f"
		i__1 = *ka - 1;
#line 1270 "ssbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1271 "ssbgst.f"
		    slartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &work[*n + j1], &
			    work[j1], &ka1);
#line 1273 "ssbgst.f"
/* L810: */
#line 1273 "ssbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1278 "ssbgst.f"
		slar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &work[*n + j1]
			, &work[j1], &ka1);

#line 1282 "ssbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 1286 "ssbgst.f"
	    i__1 = *kb - k + 1;
#line 1286 "ssbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1287 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1288 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1289 "ssbgst.f"
		if (nrt > 0) {
#line 1289 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + j1t], &work[j1t], &ka1);
#line 1289 "ssbgst.f"
		}
#line 1293 "ssbgst.f"
/* L820: */
#line 1293 "ssbgst.f"
	    }

#line 1295 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1299 "ssbgst.f"
		i__1 = j2;
#line 1299 "ssbgst.f"
		i__3 = ka1;
#line 1299 "ssbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 1300 "ssbgst.f"
		    srot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + j], &work[j]);
#line 1302 "ssbgst.f"
/* L830: */
#line 1302 "ssbgst.f"
		}
#line 1303 "ssbgst.f"
	    }
#line 1304 "ssbgst.f"
/* L840: */
#line 1304 "ssbgst.f"
	}

#line 1306 "ssbgst.f"
	if (update) {
#line 1307 "ssbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt,i+kbt-ka-1) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1312 "ssbgst.f"
		work[m - *kb + i__ + kbt] = -bb[kbt + 1 + i__ * bb_dim1] * 
			ra1;
#line 1313 "ssbgst.f"
	    }
#line 1314 "ssbgst.f"
	}

#line 1316 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1317 "ssbgst.f"
	    if (update) {
/* Computing MAX */
#line 1318 "ssbgst.f"
		i__4 = 2, i__3 = k + i0 - m;
#line 1318 "ssbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1319 "ssbgst.f"
	    } else {
/* Computing MAX */
#line 1320 "ssbgst.f"
		i__4 = 1, i__3 = k + i0 - m;
#line 1320 "ssbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1321 "ssbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 1325 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1326 "ssbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1327 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1328 "ssbgst.f"
		if (nrt > 0) {
#line 1328 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + (j1t + l - 1) * ab_dim1], 
			    &inca, &ab[ka1 - l + (j1t + l - 1) * ab_dim1], &
			    inca, &work[*n + m - *kb + j1t + *ka], &work[m - *
			    kb + j1t + *ka], &ka1);
#line 1328 "ssbgst.f"
		}
#line 1333 "ssbgst.f"
/* L850: */
#line 1333 "ssbgst.f"
	    }
#line 1334 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1335 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1336 "ssbgst.f"
	    i__4 = j2;
#line 1336 "ssbgst.f"
	    i__3 = ka1;
#line 1336 "ssbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1337 "ssbgst.f"
		work[m - *kb + j] = work[m - *kb + j + *ka];
#line 1338 "ssbgst.f"
		work[*n + m - *kb + j] = work[*n + m - *kb + j + *ka];
#line 1339 "ssbgst.f"
/* L860: */
#line 1339 "ssbgst.f"
	    }
#line 1340 "ssbgst.f"
	    i__3 = j2;
#line 1340 "ssbgst.f"
	    i__4 = ka1;
#line 1340 "ssbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1345 "ssbgst.f"
		work[m - *kb + j] *= ab[ka1 + (j - 1) * ab_dim1];
#line 1346 "ssbgst.f"
		ab[ka1 + (j - 1) * ab_dim1] = work[*n + m - *kb + j] * ab[ka1 
			+ (j - 1) * ab_dim1];
#line 1347 "ssbgst.f"
/* L870: */
#line 1347 "ssbgst.f"
	    }
#line 1348 "ssbgst.f"
	    if (update) {
#line 1349 "ssbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1349 "ssbgst.f"
		    work[m - *kb + i__ + k - *ka] = work[m - *kb + i__ + k];
#line 1349 "ssbgst.f"
		}
#line 1351 "ssbgst.f"
	    }
#line 1352 "ssbgst.f"
/* L880: */
#line 1352 "ssbgst.f"
	}

#line 1354 "ssbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1355 "ssbgst.f"
	    i__4 = 1, i__3 = k + i0 - m;
#line 1355 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1356 "ssbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1357 "ssbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1358 "ssbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1363 "ssbgst.f"
		slargv_(&nr, &ab[ka1 + j1 * ab_dim1], &inca, &work[m - *kb + 
			j1], &ka1, &work[*n + m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the right */

#line 1368 "ssbgst.f"
		i__4 = *ka - 1;
#line 1368 "ssbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 1369 "ssbgst.f"
		    slartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &work[*n + m - *kb 
			    + j1], &work[m - *kb + j1], &ka1);
#line 1372 "ssbgst.f"
/* L890: */
#line 1372 "ssbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1377 "ssbgst.f"
		slar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &work[*n + m 
			- *kb + j1], &work[m - *kb + j1], &ka1);

#line 1381 "ssbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 1385 "ssbgst.f"
	    i__4 = *kb - k + 1;
#line 1385 "ssbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 1386 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1387 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1388 "ssbgst.f"
		if (nrt > 0) {
#line 1388 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + m - *kb + j1t], &work[m - *kb 
			    + j1t], &ka1);
#line 1388 "ssbgst.f"
		}
#line 1393 "ssbgst.f"
/* L900: */
#line 1393 "ssbgst.f"
	    }

#line 1395 "ssbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1399 "ssbgst.f"
		i__4 = j2;
#line 1399 "ssbgst.f"
		i__3 = ka1;
#line 1399 "ssbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1400 "ssbgst.f"
		    srot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &work[*n + m - *kb + j], &work[m - *
			    kb + j]);
#line 1402 "ssbgst.f"
/* L910: */
#line 1402 "ssbgst.f"
		}
#line 1403 "ssbgst.f"
	    }
#line 1404 "ssbgst.f"
/* L920: */
#line 1404 "ssbgst.f"
	}

#line 1406 "ssbgst.f"
	i__3 = *kb - 1;
#line 1406 "ssbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1407 "ssbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1407 "ssbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 1411 "ssbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1412 "ssbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1413 "ssbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1414 "ssbgst.f"
		if (nrt > 0) {
#line 1414 "ssbgst.f"
		    slartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &work[*n + j1t], &work[j1t], &ka1);
#line 1414 "ssbgst.f"
		}
#line 1418 "ssbgst.f"
/* L930: */
#line 1418 "ssbgst.f"
	    }
#line 1419 "ssbgst.f"
/* L940: */
#line 1419 "ssbgst.f"
	}

#line 1421 "ssbgst.f"
	if (*kb > 1) {
/* Computing MIN */
#line 1422 "ssbgst.f"
	    i__4 = i__ + *kb;
#line 1422 "ssbgst.f"
	    i__3 = min(i__4,m) - (*ka << 1) - 1;
#line 1422 "ssbgst.f"
	    for (j = 2; j <= i__3; ++j) {
#line 1423 "ssbgst.f"
		work[*n + j] = work[*n + j + *ka];
#line 1424 "ssbgst.f"
		work[j] = work[j + *ka];
#line 1425 "ssbgst.f"
/* L950: */
#line 1425 "ssbgst.f"
	    }
#line 1426 "ssbgst.f"
	}

#line 1428 "ssbgst.f"
    }

#line 1430 "ssbgst.f"
    goto L490;

/*     End of SSBGST */

} /* ssbgst_ */

