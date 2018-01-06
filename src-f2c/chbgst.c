#line 1 "chbgst.f"
/* chbgst.f -- translated by f2c (version 20100827).
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

#line 1 "chbgst.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, */
/*                          LDX, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBGST reduces a complex Hermitian-definite banded generalized */
/* > eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y, */
/* > such that C has the same bandwidth as A. */
/* > */
/* > B must have been previously factorized as S**H*S by CPBSTF, using a */
/* > split Cholesky factorization. A is overwritten by C = X**H*A*X, where */
/* > X = S**(-1)*Q and Q is a unitary matrix chosen to preserve the */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first ka+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka). */
/* > */
/* >          On exit, the transformed matrix X**H*A*X, stored in the same */
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
/* >          BB is COMPLEX array, dimension (LDBB,N) */
/* >          The banded factor S from the split Cholesky factorization of */
/* >          B, as returned by CPBSTF, stored in the first kb+1 rows of */
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
/* >          X is COMPLEX array, dimension (LDX,N) */
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
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int chbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work, 
	doublereal *rwork, integer *info, ftnlen vect_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, bb_dim1, bb_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublecomplex t;
    static integer i0, i1, i2, j1, j2;
    static doublecomplex ra;
    static integer nr, nx, ka1, kb1;
    static doublecomplex ra1;
    static integer j1t, j2t;
    static doublereal bii;
    static integer kbt, nrt, inca;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *), 
	    cgerc_(integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    ;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int cgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper, wantx;
    extern /* Subroutine */ int clar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *), clacgv_(integer *, doublecomplex *, 
	    integer *), csscal_(integer *, doublereal *, doublecomplex *, 
	    integer *), claset_(char *, integer *, integer *, doublecomplex *,
	     doublecomplex *, doublecomplex *, integer *, ftnlen), clartg_(
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *), xerbla_(char *, integer *, ftnlen), clargv_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *);
    static logical update;
    extern /* Subroutine */ int clartv_(integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 213 "chbgst.f"
    /* Parameter adjustments */
#line 213 "chbgst.f"
    ab_dim1 = *ldab;
#line 213 "chbgst.f"
    ab_offset = 1 + ab_dim1;
#line 213 "chbgst.f"
    ab -= ab_offset;
#line 213 "chbgst.f"
    bb_dim1 = *ldbb;
#line 213 "chbgst.f"
    bb_offset = 1 + bb_dim1;
#line 213 "chbgst.f"
    bb -= bb_offset;
#line 213 "chbgst.f"
    x_dim1 = *ldx;
#line 213 "chbgst.f"
    x_offset = 1 + x_dim1;
#line 213 "chbgst.f"
    x -= x_offset;
#line 213 "chbgst.f"
    --work;
#line 213 "chbgst.f"
    --rwork;
#line 213 "chbgst.f"

#line 213 "chbgst.f"
    /* Function Body */
#line 213 "chbgst.f"
    wantx = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 214 "chbgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 215 "chbgst.f"
    ka1 = *ka + 1;
#line 216 "chbgst.f"
    kb1 = *kb + 1;
#line 217 "chbgst.f"
    *info = 0;
#line 218 "chbgst.f"
    if (! wantx && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 219 "chbgst.f"
	*info = -1;
#line 220 "chbgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 221 "chbgst.f"
	*info = -2;
#line 222 "chbgst.f"
    } else if (*n < 0) {
#line 223 "chbgst.f"
	*info = -3;
#line 224 "chbgst.f"
    } else if (*ka < 0) {
#line 225 "chbgst.f"
	*info = -4;
#line 226 "chbgst.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 227 "chbgst.f"
	*info = -5;
#line 228 "chbgst.f"
    } else if (*ldab < *ka + 1) {
#line 229 "chbgst.f"
	*info = -7;
#line 230 "chbgst.f"
    } else if (*ldbb < *kb + 1) {
#line 231 "chbgst.f"
	*info = -9;
#line 232 "chbgst.f"
    } else if (*ldx < 1 || wantx && *ldx < max(1,*n)) {
#line 233 "chbgst.f"
	*info = -11;
#line 234 "chbgst.f"
    }
#line 235 "chbgst.f"
    if (*info != 0) {
#line 236 "chbgst.f"
	i__1 = -(*info);
#line 236 "chbgst.f"
	xerbla_("CHBGST", &i__1, (ftnlen)6);
#line 237 "chbgst.f"
	return 0;
#line 238 "chbgst.f"
    }

/*     Quick return if possible */

#line 242 "chbgst.f"
    if (*n == 0) {
#line 242 "chbgst.f"
	return 0;
#line 242 "chbgst.f"
    }

#line 245 "chbgst.f"
    inca = *ldab * ka1;

/*     Initialize X to the unit matrix, if needed */

#line 249 "chbgst.f"
    if (wantx) {
#line 249 "chbgst.f"
	claset_("Full", n, n, &c_b1, &c_b2, &x[x_offset], ldx, (ftnlen)4);
#line 249 "chbgst.f"
    }

/*     Set M to the splitting point m. It must be the same value as is */
/*     used in CPBSTF. The chosen value allows the arrays WORK and RWORK */
/*     to be of dimension (N). */

#line 256 "chbgst.f"
    m = (*n + *kb) / 2;

/*     The routine works in two phases, corresponding to the two halves */
/*     of the split Cholesky factorization of B as S**H*S where */

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
/*     inv(S(i))**H*A*inv(S(i)). This creates a triangular bulge outside */
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

/*     The real cosines and complex sines of the rotations are stored in */
/*     the arrays RWORK and WORK, those of the 1st set in elements */
/*     2:m-kb-1, and those of the 2nd set in elements m-kb+1:n. */

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

#line 317 "chbgst.f"
    update = TRUE_;
#line 318 "chbgst.f"
    i__ = *n + 1;
#line 319 "chbgst.f"
L10:
#line 320 "chbgst.f"
    if (update) {
#line 321 "chbgst.f"
	--i__;
/* Computing MIN */
#line 322 "chbgst.f"
	i__1 = *kb, i__2 = i__ - 1;
#line 322 "chbgst.f"
	kbt = min(i__1,i__2);
#line 323 "chbgst.f"
	i0 = i__ - 1;
/* Computing MIN */
#line 324 "chbgst.f"
	i__1 = *n, i__2 = i__ + *ka;
#line 324 "chbgst.f"
	i1 = min(i__1,i__2);
#line 325 "chbgst.f"
	i2 = i__ - kbt + ka1;
#line 326 "chbgst.f"
	if (i__ < m + 1) {
#line 327 "chbgst.f"
	    update = FALSE_;
#line 328 "chbgst.f"
	    ++i__;
#line 329 "chbgst.f"
	    i0 = m;
#line 330 "chbgst.f"
	    if (*ka == 0) {
#line 330 "chbgst.f"
		goto L480;
#line 330 "chbgst.f"
	    }
#line 332 "chbgst.f"
	    goto L10;
#line 333 "chbgst.f"
	}
#line 334 "chbgst.f"
    } else {
#line 335 "chbgst.f"
	i__ += *ka;
#line 336 "chbgst.f"
	if (i__ > *n - 1) {
#line 336 "chbgst.f"
	    goto L480;
#line 336 "chbgst.f"
	}
#line 338 "chbgst.f"
    }

#line 340 "chbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 344 "chbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 348 "chbgst.f"
	    i__1 = kb1 + i__ * bb_dim1;
#line 348 "chbgst.f"
	    bii = bb[i__1].r;
#line 349 "chbgst.f"
	    i__1 = ka1 + i__ * ab_dim1;
#line 349 "chbgst.f"
	    i__2 = ka1 + i__ * ab_dim1;
#line 349 "chbgst.f"
	    d__1 = ab[i__2].r / bii / bii;
#line 349 "chbgst.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 350 "chbgst.f"
	    i__1 = i1;
#line 350 "chbgst.f"
	    for (j = i__ + 1; j <= i__1; ++j) {
#line 351 "chbgst.f"
		i__2 = i__ - j + ka1 + j * ab_dim1;
#line 351 "chbgst.f"
		i__3 = i__ - j + ka1 + j * ab_dim1;
#line 351 "chbgst.f"
		z__1.r = ab[i__3].r / bii, z__1.i = ab[i__3].i / bii;
#line 351 "chbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 352 "chbgst.f"
/* L20: */
#line 352 "chbgst.f"
	    }
/* Computing MAX */
#line 353 "chbgst.f"
	    i__1 = 1, i__2 = i__ - *ka;
#line 353 "chbgst.f"
	    i__3 = i__ - 1;
#line 353 "chbgst.f"
	    for (j = max(i__1,i__2); j <= i__3; ++j) {
#line 354 "chbgst.f"
		i__1 = j - i__ + ka1 + i__ * ab_dim1;
#line 354 "chbgst.f"
		i__2 = j - i__ + ka1 + i__ * ab_dim1;
#line 354 "chbgst.f"
		z__1.r = ab[i__2].r / bii, z__1.i = ab[i__2].i / bii;
#line 354 "chbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 355 "chbgst.f"
/* L30: */
#line 355 "chbgst.f"
	    }
#line 356 "chbgst.f"
	    i__3 = i__ - 1;
#line 356 "chbgst.f"
	    for (k = i__ - kbt; k <= i__3; ++k) {
#line 357 "chbgst.f"
		i__1 = k;
#line 357 "chbgst.f"
		for (j = i__ - kbt; j <= i__1; ++j) {
#line 358 "chbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 358 "chbgst.f"
		    i__4 = j - k + ka1 + k * ab_dim1;
#line 358 "chbgst.f"
		    i__5 = j - i__ + kb1 + i__ * bb_dim1;
#line 358 "chbgst.f"
		    d_cnjg(&z__5, &ab[k - i__ + ka1 + i__ * ab_dim1]);
#line 358 "chbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 358 "chbgst.f"
		    z__3.r = ab[i__4].r - z__4.r, z__3.i = ab[i__4].i - 
			    z__4.i;
#line 358 "chbgst.f"
		    d_cnjg(&z__7, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 358 "chbgst.f"
		    i__6 = j - i__ + ka1 + i__ * ab_dim1;
#line 358 "chbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 358 "chbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 358 "chbgst.f"
		    i__7 = ka1 + i__ * ab_dim1;
#line 358 "chbgst.f"
		    d__1 = ab[i__7].r;
#line 358 "chbgst.f"
		    i__8 = j - i__ + kb1 + i__ * bb_dim1;
#line 358 "chbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 358 "chbgst.f"
		    d_cnjg(&z__10, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 358 "chbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 358 "chbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 358 "chbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 366 "chbgst.f"
/* L40: */
#line 366 "chbgst.f"
		}
/* Computing MAX */
#line 367 "chbgst.f"
		i__1 = 1, i__2 = i__ - *ka;
#line 367 "chbgst.f"
		i__4 = i__ - kbt - 1;
#line 367 "chbgst.f"
		for (j = max(i__1,i__2); j <= i__4; ++j) {
#line 368 "chbgst.f"
		    i__1 = j - k + ka1 + k * ab_dim1;
#line 368 "chbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 368 "chbgst.f"
		    d_cnjg(&z__3, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 368 "chbgst.f"
		    i__5 = j - i__ + ka1 + i__ * ab_dim1;
#line 368 "chbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 368 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 368 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 371 "chbgst.f"
/* L50: */
#line 371 "chbgst.f"
		}
#line 372 "chbgst.f"
/* L60: */
#line 372 "chbgst.f"
	    }
#line 373 "chbgst.f"
	    i__3 = i1;
#line 373 "chbgst.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MAX */
#line 374 "chbgst.f"
		i__4 = j - *ka, i__1 = i__ - kbt;
#line 374 "chbgst.f"
		i__2 = i__ - 1;
#line 374 "chbgst.f"
		for (k = max(i__4,i__1); k <= i__2; ++k) {
#line 375 "chbgst.f"
		    i__4 = k - j + ka1 + j * ab_dim1;
#line 375 "chbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 375 "chbgst.f"
		    i__5 = k - i__ + kb1 + i__ * bb_dim1;
#line 375 "chbgst.f"
		    i__6 = i__ - j + ka1 + j * ab_dim1;
#line 375 "chbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 375 "chbgst.f"
		    z__1.r = ab[i__1].r - z__2.r, z__1.i = ab[i__1].i - 
			    z__2.i;
#line 375 "chbgst.f"
		    ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 377 "chbgst.f"
/* L70: */
#line 377 "chbgst.f"
		}
#line 378 "chbgst.f"
/* L80: */
#line 378 "chbgst.f"
	    }

#line 380 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 384 "chbgst.f"
		i__3 = *n - m;
#line 384 "chbgst.f"
		d__1 = 1. / bii;
#line 384 "chbgst.f"
		csscal_(&i__3, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 385 "chbgst.f"
		if (kbt > 0) {
#line 385 "chbgst.f"
		    i__3 = *n - m;
#line 385 "chbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 385 "chbgst.f"
		    cgerc_(&i__3, &kbt, &z__1, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kb1 - kbt + i__ * bb_dim1], &c__1, &x[m 
			    + 1 + (i__ - kbt) * x_dim1], ldx);
#line 385 "chbgst.f"
		}
#line 389 "chbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 393 "chbgst.f"
	    i__3 = i__ - i1 + ka1 + i1 * ab_dim1;
#line 393 "chbgst.f"
	    ra1.r = ab[i__3].r, ra1.i = ab[i__3].i;
#line 394 "chbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 400 "chbgst.f"
	i__3 = *kb - 1;
#line 400 "chbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 401 "chbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 406 "chbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i,i-k+ka+1) */

#line 410 "chbgst.f"
		    clartg_(&ab[k + 1 + (i__ - k + *ka) * ab_dim1], &ra1, &
			    rwork[i__ - k + *ka - m], &work[i__ - k + *ka - m]
			    , &ra);

/*                 create nonzero element a(i-k,i-k+ka+1) outside the */
/*                 band and store it in WORK(i-k) */

#line 416 "chbgst.f"
		    i__2 = kb1 - k + i__ * bb_dim1;
#line 416 "chbgst.f"
		    z__2.r = -bb[i__2].r, z__2.i = -bb[i__2].i;
#line 416 "chbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 416 "chbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 417 "chbgst.f"
		    i__2 = i__ - k;
#line 417 "chbgst.f"
		    i__4 = i__ - k + *ka - m;
#line 417 "chbgst.f"
		    z__2.r = rwork[i__4] * t.r, z__2.i = rwork[i__4] * t.i;
#line 417 "chbgst.f"
		    d_cnjg(&z__4, &work[i__ - k + *ka - m]);
#line 417 "chbgst.f"
		    i__1 = (i__ - k + *ka) * ab_dim1 + 1;
#line 417 "chbgst.f"
		    z__3.r = z__4.r * ab[i__1].r - z__4.i * ab[i__1].i, 
			    z__3.i = z__4.r * ab[i__1].i + z__4.i * ab[i__1]
			    .r;
#line 417 "chbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 417 "chbgst.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 420 "chbgst.f"
		    i__2 = (i__ - k + *ka) * ab_dim1 + 1;
#line 420 "chbgst.f"
		    i__4 = i__ - k + *ka - m;
#line 420 "chbgst.f"
		    z__2.r = work[i__4].r * t.r - work[i__4].i * t.i, z__2.i =
			     work[i__4].r * t.i + work[i__4].i * t.r;
#line 420 "chbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 420 "chbgst.f"
		    i__5 = (i__ - k + *ka) * ab_dim1 + 1;
#line 420 "chbgst.f"
		    z__3.r = rwork[i__1] * ab[i__5].r, z__3.i = rwork[i__1] * 
			    ab[i__5].i;
#line 420 "chbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 420 "chbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 422 "chbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 423 "chbgst.f"
		}
#line 424 "chbgst.f"
	    }
/* Computing MAX */
#line 425 "chbgst.f"
	    i__2 = 1, i__4 = k - i0 + 2;
#line 425 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__2,i__4) * ka1;
#line 426 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 427 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 428 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 429 "chbgst.f"
		i__2 = j2, i__4 = i__ + (*ka << 1) - k + 1;
#line 429 "chbgst.f"
		j2t = max(i__2,i__4);
#line 430 "chbgst.f"
	    } else {
#line 431 "chbgst.f"
		j2t = j2;
#line 432 "chbgst.f"
	    }
#line 433 "chbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 434 "chbgst.f"
	    i__2 = j1;
#line 434 "chbgst.f"
	    i__4 = ka1;
#line 434 "chbgst.f"
	    for (j = j2t; i__4 < 0 ? j >= i__2 : j <= i__2; j += i__4) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j-m) */

#line 439 "chbgst.f"
		i__1 = j - m;
#line 439 "chbgst.f"
		i__5 = j - m;
#line 439 "chbgst.f"
		i__6 = (j + 1) * ab_dim1 + 1;
#line 439 "chbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 439 "chbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 440 "chbgst.f"
		i__1 = (j + 1) * ab_dim1 + 1;
#line 440 "chbgst.f"
		i__5 = j - m;
#line 440 "chbgst.f"
		i__6 = (j + 1) * ab_dim1 + 1;
#line 440 "chbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 440 "chbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 441 "chbgst.f"
/* L90: */
#line 441 "chbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 446 "chbgst.f"
	    if (nrt > 0) {
#line 446 "chbgst.f"
		clargv_(&nrt, &ab[j2t * ab_dim1 + 1], &inca, &work[j2t - m], &
			ka1, &rwork[j2t - m], &ka1);
#line 446 "chbgst.f"
	    }
#line 449 "chbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 453 "chbgst.f"
		i__4 = *ka - 1;
#line 453 "chbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 454 "chbgst.f"
		    clartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - m], 
			    &work[j2 - m], &ka1);
#line 457 "chbgst.f"
/* L100: */
#line 457 "chbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 462 "chbgst.f"
		clar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &
			rwork[j2 - m], &work[j2 - m], &ka1);

#line 466 "chbgst.f"
		clacgv_(&nr, &work[j2 - m], &ka1);
#line 467 "chbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 471 "chbgst.f"
	    i__4 = *kb - k + 1;
#line 471 "chbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 472 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 473 "chbgst.f"
		if (nrt > 0) {
#line 473 "chbgst.f"
		    clartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2 - m], &work[j2 - m], &ka1);
#line 473 "chbgst.f"
		}
#line 477 "chbgst.f"
/* L110: */
#line 477 "chbgst.f"
	    }

#line 479 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 483 "chbgst.f"
		i__4 = j1;
#line 483 "chbgst.f"
		i__2 = ka1;
#line 483 "chbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__4 : j <= i__4; j += i__2) {
#line 484 "chbgst.f"
		    i__1 = *n - m;
#line 484 "chbgst.f"
		    d_cnjg(&z__1, &work[j - m]);
#line 484 "chbgst.f"
		    crot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j - m], &z__1);
#line 486 "chbgst.f"
/* L120: */
#line 486 "chbgst.f"
		}
#line 487 "chbgst.f"
	    }
#line 488 "chbgst.f"
/* L130: */
#line 488 "chbgst.f"
	}

#line 490 "chbgst.f"
	if (update) {
#line 491 "chbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt,i-kbt+ka+1) outside the */
/*              band and store it in WORK(i-kbt) */

#line 496 "chbgst.f"
		i__3 = i__ - kbt;
#line 496 "chbgst.f"
		i__2 = kb1 - kbt + i__ * bb_dim1;
#line 496 "chbgst.f"
		z__2.r = -bb[i__2].r, z__2.i = -bb[i__2].i;
#line 496 "chbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 496 "chbgst.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 497 "chbgst.f"
	    }
#line 498 "chbgst.f"
	}

#line 500 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 501 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 502 "chbgst.f"
		i__3 = 2, i__2 = k - i0 + 1;
#line 502 "chbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 503 "chbgst.f"
	    } else {
/* Computing MAX */
#line 504 "chbgst.f"
		i__3 = 1, i__2 = k - i0 + 1;
#line 504 "chbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 505 "chbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 509 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 510 "chbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 511 "chbgst.f"
		if (nrt > 0) {
#line 511 "chbgst.f"
		    clartv_(&nrt, &ab[l + (j2 - l + 1) * ab_dim1], &inca, &ab[
			    l + 1 + (j2 - l + 1) * ab_dim1], &inca, &rwork[j2 
			    - *ka], &work[j2 - *ka], &ka1);
#line 511 "chbgst.f"
		}
#line 515 "chbgst.f"
/* L140: */
#line 515 "chbgst.f"
	    }
#line 516 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 517 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 518 "chbgst.f"
	    i__3 = j2;
#line 518 "chbgst.f"
	    i__2 = -ka1;
#line 518 "chbgst.f"
	    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 519 "chbgst.f"
		i__4 = j;
#line 519 "chbgst.f"
		i__1 = j - *ka;
#line 519 "chbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 520 "chbgst.f"
		rwork[j] = rwork[j - *ka];
#line 521 "chbgst.f"
/* L150: */
#line 521 "chbgst.f"
	    }
#line 522 "chbgst.f"
	    i__2 = j1;
#line 522 "chbgst.f"
	    i__3 = ka1;
#line 522 "chbgst.f"
	    for (j = j2; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j) */

#line 527 "chbgst.f"
		i__4 = j;
#line 527 "chbgst.f"
		i__1 = j;
#line 527 "chbgst.f"
		i__5 = (j + 1) * ab_dim1 + 1;
#line 527 "chbgst.f"
		z__1.r = work[i__1].r * ab[i__5].r - work[i__1].i * ab[i__5]
			.i, z__1.i = work[i__1].r * ab[i__5].i + work[i__1].i 
			* ab[i__5].r;
#line 527 "chbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 528 "chbgst.f"
		i__4 = (j + 1) * ab_dim1 + 1;
#line 528 "chbgst.f"
		i__1 = j;
#line 528 "chbgst.f"
		i__5 = (j + 1) * ab_dim1 + 1;
#line 528 "chbgst.f"
		z__1.r = rwork[i__1] * ab[i__5].r, z__1.i = rwork[i__1] * ab[
			i__5].i;
#line 528 "chbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 529 "chbgst.f"
/* L160: */
#line 529 "chbgst.f"
	    }
#line 530 "chbgst.f"
	    if (update) {
#line 531 "chbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 531 "chbgst.f"
		    i__3 = i__ - k + *ka;
#line 531 "chbgst.f"
		    i__2 = i__ - k;
#line 531 "chbgst.f"
		    work[i__3].r = work[i__2].r, work[i__3].i = work[i__2].i;
#line 531 "chbgst.f"
		}
#line 533 "chbgst.f"
	    }
#line 534 "chbgst.f"
/* L170: */
#line 534 "chbgst.f"
	}

#line 536 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 537 "chbgst.f"
	    i__3 = 1, i__2 = k - i0 + 1;
#line 537 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 538 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 539 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 540 "chbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 545 "chbgst.f"
		clargv_(&nr, &ab[j2 * ab_dim1 + 1], &inca, &work[j2], &ka1, &
			rwork[j2], &ka1);

/*              apply rotations in 2nd set from the right */

#line 550 "chbgst.f"
		i__3 = *ka - 1;
#line 550 "chbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 551 "chbgst.f"
		    clartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2], &
			    work[j2], &ka1);
#line 554 "chbgst.f"
/* L180: */
#line 554 "chbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 559 "chbgst.f"
		clar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &
			rwork[j2], &work[j2], &ka1);

#line 563 "chbgst.f"
		clacgv_(&nr, &work[j2], &ka1);
#line 564 "chbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 568 "chbgst.f"
	    i__3 = *kb - k + 1;
#line 568 "chbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 569 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 570 "chbgst.f"
		if (nrt > 0) {
#line 570 "chbgst.f"
		    clartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2], &work[j2], &ka1);
#line 570 "chbgst.f"
		}
#line 574 "chbgst.f"
/* L190: */
#line 574 "chbgst.f"
	    }

#line 576 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 580 "chbgst.f"
		i__3 = j1;
#line 580 "chbgst.f"
		i__2 = ka1;
#line 580 "chbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 581 "chbgst.f"
		    i__4 = *n - m;
#line 581 "chbgst.f"
		    d_cnjg(&z__1, &work[j]);
#line 581 "chbgst.f"
		    crot_(&i__4, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j], &z__1);
#line 583 "chbgst.f"
/* L200: */
#line 583 "chbgst.f"
		}
#line 584 "chbgst.f"
	    }
#line 585 "chbgst.f"
/* L210: */
#line 585 "chbgst.f"
	}

#line 587 "chbgst.f"
	i__2 = *kb - 1;
#line 587 "chbgst.f"
	for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 588 "chbgst.f"
	    i__3 = 1, i__4 = k - i0 + 2;
#line 588 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__4) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 592 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 593 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 594 "chbgst.f"
		if (nrt > 0) {
#line 594 "chbgst.f"
		    clartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2 - m], &work[j2 - m], &ka1);
#line 594 "chbgst.f"
		}
#line 598 "chbgst.f"
/* L220: */
#line 598 "chbgst.f"
	    }
#line 599 "chbgst.f"
/* L230: */
#line 599 "chbgst.f"
	}

#line 601 "chbgst.f"
	if (*kb > 1) {
#line 602 "chbgst.f"
	    i__2 = j2 + *ka;
#line 602 "chbgst.f"
	    for (j = *n - 1; j >= i__2; --j) {
#line 603 "chbgst.f"
		rwork[j - m] = rwork[j - *ka - m];
#line 604 "chbgst.f"
		i__3 = j - m;
#line 604 "chbgst.f"
		i__4 = j - *ka - m;
#line 604 "chbgst.f"
		work[i__3].r = work[i__4].r, work[i__3].i = work[i__4].i;
#line 605 "chbgst.f"
/* L240: */
#line 605 "chbgst.f"
	    }
#line 606 "chbgst.f"
	}

#line 608 "chbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 612 "chbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 616 "chbgst.f"
	    i__2 = i__ * bb_dim1 + 1;
#line 616 "chbgst.f"
	    bii = bb[i__2].r;
#line 617 "chbgst.f"
	    i__2 = i__ * ab_dim1 + 1;
#line 617 "chbgst.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 617 "chbgst.f"
	    d__1 = ab[i__3].r / bii / bii;
#line 617 "chbgst.f"
	    ab[i__2].r = d__1, ab[i__2].i = 0.;
#line 618 "chbgst.f"
	    i__2 = i1;
#line 618 "chbgst.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 619 "chbgst.f"
		i__3 = j - i__ + 1 + i__ * ab_dim1;
#line 619 "chbgst.f"
		i__4 = j - i__ + 1 + i__ * ab_dim1;
#line 619 "chbgst.f"
		z__1.r = ab[i__4].r / bii, z__1.i = ab[i__4].i / bii;
#line 619 "chbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 620 "chbgst.f"
/* L250: */
#line 620 "chbgst.f"
	    }
/* Computing MAX */
#line 621 "chbgst.f"
	    i__2 = 1, i__3 = i__ - *ka;
#line 621 "chbgst.f"
	    i__4 = i__ - 1;
#line 621 "chbgst.f"
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 622 "chbgst.f"
		i__2 = i__ - j + 1 + j * ab_dim1;
#line 622 "chbgst.f"
		i__3 = i__ - j + 1 + j * ab_dim1;
#line 622 "chbgst.f"
		z__1.r = ab[i__3].r / bii, z__1.i = ab[i__3].i / bii;
#line 622 "chbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 623 "chbgst.f"
/* L260: */
#line 623 "chbgst.f"
	    }
#line 624 "chbgst.f"
	    i__4 = i__ - 1;
#line 624 "chbgst.f"
	    for (k = i__ - kbt; k <= i__4; ++k) {
#line 625 "chbgst.f"
		i__2 = k;
#line 625 "chbgst.f"
		for (j = i__ - kbt; j <= i__2; ++j) {
#line 626 "chbgst.f"
		    i__3 = k - j + 1 + j * ab_dim1;
#line 626 "chbgst.f"
		    i__1 = k - j + 1 + j * ab_dim1;
#line 626 "chbgst.f"
		    i__5 = i__ - j + 1 + j * bb_dim1;
#line 626 "chbgst.f"
		    d_cnjg(&z__5, &ab[i__ - k + 1 + k * ab_dim1]);
#line 626 "chbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 626 "chbgst.f"
		    z__3.r = ab[i__1].r - z__4.r, z__3.i = ab[i__1].i - 
			    z__4.i;
#line 626 "chbgst.f"
		    d_cnjg(&z__7, &bb[i__ - k + 1 + k * bb_dim1]);
#line 626 "chbgst.f"
		    i__6 = i__ - j + 1 + j * ab_dim1;
#line 626 "chbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 626 "chbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 626 "chbgst.f"
		    i__7 = i__ * ab_dim1 + 1;
#line 626 "chbgst.f"
		    d__1 = ab[i__7].r;
#line 626 "chbgst.f"
		    i__8 = i__ - j + 1 + j * bb_dim1;
#line 626 "chbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 626 "chbgst.f"
		    d_cnjg(&z__10, &bb[i__ - k + 1 + k * bb_dim1]);
#line 626 "chbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 626 "chbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 626 "chbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 632 "chbgst.f"
/* L270: */
#line 632 "chbgst.f"
		}
/* Computing MAX */
#line 633 "chbgst.f"
		i__2 = 1, i__3 = i__ - *ka;
#line 633 "chbgst.f"
		i__1 = i__ - kbt - 1;
#line 633 "chbgst.f"
		for (j = max(i__2,i__3); j <= i__1; ++j) {
#line 634 "chbgst.f"
		    i__2 = k - j + 1 + j * ab_dim1;
#line 634 "chbgst.f"
		    i__3 = k - j + 1 + j * ab_dim1;
#line 634 "chbgst.f"
		    d_cnjg(&z__3, &bb[i__ - k + 1 + k * bb_dim1]);
#line 634 "chbgst.f"
		    i__5 = i__ - j + 1 + j * ab_dim1;
#line 634 "chbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 634 "chbgst.f"
		    z__1.r = ab[i__3].r - z__2.r, z__1.i = ab[i__3].i - 
			    z__2.i;
#line 634 "chbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 637 "chbgst.f"
/* L280: */
#line 637 "chbgst.f"
		}
#line 638 "chbgst.f"
/* L290: */
#line 638 "chbgst.f"
	    }
#line 639 "chbgst.f"
	    i__4 = i1;
#line 639 "chbgst.f"
	    for (j = i__; j <= i__4; ++j) {
/* Computing MAX */
#line 640 "chbgst.f"
		i__1 = j - *ka, i__2 = i__ - kbt;
#line 640 "chbgst.f"
		i__3 = i__ - 1;
#line 640 "chbgst.f"
		for (k = max(i__1,i__2); k <= i__3; ++k) {
#line 641 "chbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 641 "chbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 641 "chbgst.f"
		    i__5 = i__ - k + 1 + k * bb_dim1;
#line 641 "chbgst.f"
		    i__6 = j - i__ + 1 + i__ * ab_dim1;
#line 641 "chbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 641 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 641 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 643 "chbgst.f"
/* L300: */
#line 643 "chbgst.f"
		}
#line 644 "chbgst.f"
/* L310: */
#line 644 "chbgst.f"
	    }

#line 646 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 650 "chbgst.f"
		i__4 = *n - m;
#line 650 "chbgst.f"
		d__1 = 1. / bii;
#line 650 "chbgst.f"
		csscal_(&i__4, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 651 "chbgst.f"
		if (kbt > 0) {
#line 651 "chbgst.f"
		    i__4 = *n - m;
#line 651 "chbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 651 "chbgst.f"
		    i__3 = *ldbb - 1;
#line 651 "chbgst.f"
		    cgeru_(&i__4, &kbt, &z__1, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kbt + 1 + (i__ - kbt) * bb_dim1], &i__3,
			     &x[m + 1 + (i__ - kbt) * x_dim1], ldx);
#line 651 "chbgst.f"
		}
#line 655 "chbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 659 "chbgst.f"
	    i__4 = i1 - i__ + 1 + i__ * ab_dim1;
#line 659 "chbgst.f"
	    ra1.r = ab[i__4].r, ra1.i = ab[i__4].i;
#line 660 "chbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 666 "chbgst.f"
	i__4 = *kb - 1;
#line 666 "chbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 667 "chbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 672 "chbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i-k+ka+1,i) */

#line 676 "chbgst.f"
		    clartg_(&ab[ka1 - k + i__ * ab_dim1], &ra1, &rwork[i__ - 
			    k + *ka - m], &work[i__ - k + *ka - m], &ra);

/*                 create nonzero element a(i-k+ka+1,i-k) outside the */
/*                 band and store it in WORK(i-k) */

#line 682 "chbgst.f"
		    i__3 = k + 1 + (i__ - k) * bb_dim1;
#line 682 "chbgst.f"
		    z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 682 "chbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 682 "chbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 683 "chbgst.f"
		    i__3 = i__ - k;
#line 683 "chbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 683 "chbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 683 "chbgst.f"
		    d_cnjg(&z__4, &work[i__ - k + *ka - m]);
#line 683 "chbgst.f"
		    i__2 = ka1 + (i__ - k) * ab_dim1;
#line 683 "chbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 683 "chbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 683 "chbgst.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 685 "chbgst.f"
		    i__3 = ka1 + (i__ - k) * ab_dim1;
#line 685 "chbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 685 "chbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 685 "chbgst.f"
		    i__2 = i__ - k + *ka - m;
#line 685 "chbgst.f"
		    i__5 = ka1 + (i__ - k) * ab_dim1;
#line 685 "chbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 685 "chbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 685 "chbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 687 "chbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 688 "chbgst.f"
		}
#line 689 "chbgst.f"
	    }
/* Computing MAX */
#line 690 "chbgst.f"
	    i__3 = 1, i__1 = k - i0 + 2;
#line 690 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__1) * ka1;
#line 691 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 692 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 693 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 694 "chbgst.f"
		i__3 = j2, i__1 = i__ + (*ka << 1) - k + 1;
#line 694 "chbgst.f"
		j2t = max(i__3,i__1);
#line 695 "chbgst.f"
	    } else {
#line 696 "chbgst.f"
		j2t = j2;
#line 697 "chbgst.f"
	    }
#line 698 "chbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 699 "chbgst.f"
	    i__3 = j1;
#line 699 "chbgst.f"
	    i__1 = ka1;
#line 699 "chbgst.f"
	    for (j = j2t; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j-m) */

#line 704 "chbgst.f"
		i__2 = j - m;
#line 704 "chbgst.f"
		i__5 = j - m;
#line 704 "chbgst.f"
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
#line 704 "chbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 704 "chbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 705 "chbgst.f"
		i__2 = ka1 + (j - *ka + 1) * ab_dim1;
#line 705 "chbgst.f"
		i__5 = j - m;
#line 705 "chbgst.f"
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
#line 705 "chbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 705 "chbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 706 "chbgst.f"
/* L320: */
#line 706 "chbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 711 "chbgst.f"
	    if (nrt > 0) {
#line 711 "chbgst.f"
		clargv_(&nrt, &ab[ka1 + (j2t - *ka) * ab_dim1], &inca, &work[
			j2t - m], &ka1, &rwork[j2t - m], &ka1);
#line 711 "chbgst.f"
	    }
#line 714 "chbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 718 "chbgst.f"
		i__1 = *ka - 1;
#line 718 "chbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 719 "chbgst.f"
		    clartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &rwork[j2 - m]
			    , &work[j2 - m], &ka1);
#line 722 "chbgst.f"
/* L330: */
#line 722 "chbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 727 "chbgst.f"
		clar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &rwork[j2 - m], &
			work[j2 - m], &ka1);

#line 730 "chbgst.f"
		clacgv_(&nr, &work[j2 - m], &ka1);
#line 731 "chbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 735 "chbgst.f"
	    i__1 = *kb - k + 1;
#line 735 "chbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 736 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 737 "chbgst.f"
		if (nrt > 0) {
#line 737 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - 
			    m], &work[j2 - m], &ka1);
#line 737 "chbgst.f"
		}
#line 741 "chbgst.f"
/* L340: */
#line 741 "chbgst.f"
	    }

#line 743 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 747 "chbgst.f"
		i__1 = j1;
#line 747 "chbgst.f"
		i__3 = ka1;
#line 747 "chbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 748 "chbgst.f"
		    i__2 = *n - m;
#line 748 "chbgst.f"
		    crot_(&i__2, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j - m], &work[j - m]
			    );
#line 750 "chbgst.f"
/* L350: */
#line 750 "chbgst.f"
		}
#line 751 "chbgst.f"
	    }
#line 752 "chbgst.f"
/* L360: */
#line 752 "chbgst.f"
	}

#line 754 "chbgst.f"
	if (update) {
#line 755 "chbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt+ka+1,i-kbt) outside the */
/*              band and store it in WORK(i-kbt) */

#line 760 "chbgst.f"
		i__4 = i__ - kbt;
#line 760 "chbgst.f"
		i__3 = kbt + 1 + (i__ - kbt) * bb_dim1;
#line 760 "chbgst.f"
		z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 760 "chbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 760 "chbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 761 "chbgst.f"
	    }
#line 762 "chbgst.f"
	}

#line 764 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 765 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 766 "chbgst.f"
		i__4 = 2, i__3 = k - i0 + 1;
#line 766 "chbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 767 "chbgst.f"
	    } else {
/* Computing MAX */
#line 768 "chbgst.f"
		i__4 = 1, i__3 = k - i0 + 1;
#line 768 "chbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 769 "chbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 773 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 774 "chbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 775 "chbgst.f"
		if (nrt > 0) {
#line 775 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + (j2 - *ka) * ab_dim1], &
			    inca, &ab[ka1 - l + (j2 - *ka + 1) * ab_dim1], &
			    inca, &rwork[j2 - *ka], &work[j2 - *ka], &ka1);
#line 775 "chbgst.f"
		}
#line 779 "chbgst.f"
/* L370: */
#line 779 "chbgst.f"
	    }
#line 780 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 781 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 782 "chbgst.f"
	    i__4 = j2;
#line 782 "chbgst.f"
	    i__3 = -ka1;
#line 782 "chbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 783 "chbgst.f"
		i__1 = j;
#line 783 "chbgst.f"
		i__2 = j - *ka;
#line 783 "chbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 784 "chbgst.f"
		rwork[j] = rwork[j - *ka];
#line 785 "chbgst.f"
/* L380: */
#line 785 "chbgst.f"
	    }
#line 786 "chbgst.f"
	    i__3 = j1;
#line 786 "chbgst.f"
	    i__4 = ka1;
#line 786 "chbgst.f"
	    for (j = j2; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j) */

#line 791 "chbgst.f"
		i__1 = j;
#line 791 "chbgst.f"
		i__2 = j;
#line 791 "chbgst.f"
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
#line 791 "chbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 791 "chbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 792 "chbgst.f"
		i__1 = ka1 + (j - *ka + 1) * ab_dim1;
#line 792 "chbgst.f"
		i__2 = j;
#line 792 "chbgst.f"
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
#line 792 "chbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 792 "chbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 793 "chbgst.f"
/* L390: */
#line 793 "chbgst.f"
	    }
#line 794 "chbgst.f"
	    if (update) {
#line 795 "chbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 795 "chbgst.f"
		    i__4 = i__ - k + *ka;
#line 795 "chbgst.f"
		    i__3 = i__ - k;
#line 795 "chbgst.f"
		    work[i__4].r = work[i__3].r, work[i__4].i = work[i__3].i;
#line 795 "chbgst.f"
		}
#line 797 "chbgst.f"
	    }
#line 798 "chbgst.f"
/* L400: */
#line 798 "chbgst.f"
	}

#line 800 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 801 "chbgst.f"
	    i__4 = 1, i__3 = k - i0 + 1;
#line 801 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 802 "chbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 803 "chbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 804 "chbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 809 "chbgst.f"
		clargv_(&nr, &ab[ka1 + (j2 - *ka) * ab_dim1], &inca, &work[j2]
			, &ka1, &rwork[j2], &ka1);

/*              apply rotations in 2nd set from the left */

#line 814 "chbgst.f"
		i__4 = *ka - 1;
#line 814 "chbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 815 "chbgst.f"
		    clartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &rwork[j2], &
			    work[j2], &ka1);
#line 818 "chbgst.f"
/* L410: */
#line 818 "chbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 823 "chbgst.f"
		clar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &rwork[j2], &work[
			j2], &ka1);

#line 826 "chbgst.f"
		clacgv_(&nr, &work[j2], &ka1);
#line 827 "chbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 831 "chbgst.f"
	    i__4 = *kb - k + 1;
#line 831 "chbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 832 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 833 "chbgst.f"
		if (nrt > 0) {
#line 833 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2], 
			    &work[j2], &ka1);
#line 833 "chbgst.f"
		}
#line 837 "chbgst.f"
/* L420: */
#line 837 "chbgst.f"
	    }

#line 839 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 843 "chbgst.f"
		i__4 = j1;
#line 843 "chbgst.f"
		i__3 = ka1;
#line 843 "chbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 844 "chbgst.f"
		    i__1 = *n - m;
#line 844 "chbgst.f"
		    crot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j], &work[j]);
#line 846 "chbgst.f"
/* L430: */
#line 846 "chbgst.f"
		}
#line 847 "chbgst.f"
	    }
#line 848 "chbgst.f"
/* L440: */
#line 848 "chbgst.f"
	}

#line 850 "chbgst.f"
	i__3 = *kb - 1;
#line 850 "chbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 851 "chbgst.f"
	    i__4 = 1, i__1 = k - i0 + 2;
#line 851 "chbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 855 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 856 "chbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 857 "chbgst.f"
		if (nrt > 0) {
#line 857 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - 
			    m], &work[j2 - m], &ka1);
#line 857 "chbgst.f"
		}
#line 861 "chbgst.f"
/* L450: */
#line 861 "chbgst.f"
	    }
#line 862 "chbgst.f"
/* L460: */
#line 862 "chbgst.f"
	}

#line 864 "chbgst.f"
	if (*kb > 1) {
#line 865 "chbgst.f"
	    i__3 = j2 + *ka;
#line 865 "chbgst.f"
	    for (j = *n - 1; j >= i__3; --j) {
#line 866 "chbgst.f"
		rwork[j - m] = rwork[j - *ka - m];
#line 867 "chbgst.f"
		i__4 = j - m;
#line 867 "chbgst.f"
		i__1 = j - *ka - m;
#line 867 "chbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 868 "chbgst.f"
/* L470: */
#line 868 "chbgst.f"
	    }
#line 869 "chbgst.f"
	}

#line 871 "chbgst.f"
    }

#line 873 "chbgst.f"
    goto L10;

#line 875 "chbgst.f"
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

#line 893 "chbgst.f"
    update = TRUE_;
#line 894 "chbgst.f"
    i__ = 0;
#line 895 "chbgst.f"
L490:
#line 896 "chbgst.f"
    if (update) {
#line 897 "chbgst.f"
	++i__;
/* Computing MIN */
#line 898 "chbgst.f"
	i__3 = *kb, i__4 = m - i__;
#line 898 "chbgst.f"
	kbt = min(i__3,i__4);
#line 899 "chbgst.f"
	i0 = i__ + 1;
/* Computing MAX */
#line 900 "chbgst.f"
	i__3 = 1, i__4 = i__ - *ka;
#line 900 "chbgst.f"
	i1 = max(i__3,i__4);
#line 901 "chbgst.f"
	i2 = i__ + kbt - ka1;
#line 902 "chbgst.f"
	if (i__ > m) {
#line 903 "chbgst.f"
	    update = FALSE_;
#line 904 "chbgst.f"
	    --i__;
#line 905 "chbgst.f"
	    i0 = m + 1;
#line 906 "chbgst.f"
	    if (*ka == 0) {
#line 906 "chbgst.f"
		return 0;
#line 906 "chbgst.f"
	    }
#line 908 "chbgst.f"
	    goto L490;
#line 909 "chbgst.f"
	}
#line 910 "chbgst.f"
    } else {
#line 911 "chbgst.f"
	i__ -= *ka;
#line 912 "chbgst.f"
	if (i__ < 2) {
#line 912 "chbgst.f"
	    return 0;
#line 912 "chbgst.f"
	}
#line 914 "chbgst.f"
    }

#line 916 "chbgst.f"
    if (i__ < m - kbt) {
#line 917 "chbgst.f"
	nx = m;
#line 918 "chbgst.f"
    } else {
#line 919 "chbgst.f"
	nx = *n;
#line 920 "chbgst.f"
    }

#line 922 "chbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 926 "chbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 930 "chbgst.f"
	    i__3 = kb1 + i__ * bb_dim1;
#line 930 "chbgst.f"
	    bii = bb[i__3].r;
#line 931 "chbgst.f"
	    i__3 = ka1 + i__ * ab_dim1;
#line 931 "chbgst.f"
	    i__4 = ka1 + i__ * ab_dim1;
#line 931 "chbgst.f"
	    d__1 = ab[i__4].r / bii / bii;
#line 931 "chbgst.f"
	    ab[i__3].r = d__1, ab[i__3].i = 0.;
#line 932 "chbgst.f"
	    i__3 = i__ - 1;
#line 932 "chbgst.f"
	    for (j = i1; j <= i__3; ++j) {
#line 933 "chbgst.f"
		i__4 = j - i__ + ka1 + i__ * ab_dim1;
#line 933 "chbgst.f"
		i__1 = j - i__ + ka1 + i__ * ab_dim1;
#line 933 "chbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 933 "chbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 934 "chbgst.f"
/* L500: */
#line 934 "chbgst.f"
	    }
/* Computing MIN */
#line 935 "chbgst.f"
	    i__4 = *n, i__1 = i__ + *ka;
#line 935 "chbgst.f"
	    i__3 = min(i__4,i__1);
#line 935 "chbgst.f"
	    for (j = i__ + 1; j <= i__3; ++j) {
#line 936 "chbgst.f"
		i__4 = i__ - j + ka1 + j * ab_dim1;
#line 936 "chbgst.f"
		i__1 = i__ - j + ka1 + j * ab_dim1;
#line 936 "chbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 936 "chbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 937 "chbgst.f"
/* L510: */
#line 937 "chbgst.f"
	    }
#line 938 "chbgst.f"
	    i__3 = i__ + kbt;
#line 938 "chbgst.f"
	    for (k = i__ + 1; k <= i__3; ++k) {
#line 939 "chbgst.f"
		i__4 = i__ + kbt;
#line 939 "chbgst.f"
		for (j = k; j <= i__4; ++j) {
#line 940 "chbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 940 "chbgst.f"
		    i__2 = k - j + ka1 + j * ab_dim1;
#line 940 "chbgst.f"
		    i__5 = i__ - j + kb1 + j * bb_dim1;
#line 940 "chbgst.f"
		    d_cnjg(&z__5, &ab[i__ - k + ka1 + k * ab_dim1]);
#line 940 "chbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 940 "chbgst.f"
		    z__3.r = ab[i__2].r - z__4.r, z__3.i = ab[i__2].i - 
			    z__4.i;
#line 940 "chbgst.f"
		    d_cnjg(&z__7, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 940 "chbgst.f"
		    i__6 = i__ - j + ka1 + j * ab_dim1;
#line 940 "chbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 940 "chbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 940 "chbgst.f"
		    i__7 = ka1 + i__ * ab_dim1;
#line 940 "chbgst.f"
		    d__1 = ab[i__7].r;
#line 940 "chbgst.f"
		    i__8 = i__ - j + kb1 + j * bb_dim1;
#line 940 "chbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 940 "chbgst.f"
		    d_cnjg(&z__10, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 940 "chbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 940 "chbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 940 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 948 "chbgst.f"
/* L520: */
#line 948 "chbgst.f"
		}
/* Computing MIN */
#line 949 "chbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 949 "chbgst.f"
		i__4 = min(i__1,i__2);
#line 949 "chbgst.f"
		for (j = i__ + kbt + 1; j <= i__4; ++j) {
#line 950 "chbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 950 "chbgst.f"
		    i__2 = k - j + ka1 + j * ab_dim1;
#line 950 "chbgst.f"
		    d_cnjg(&z__3, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 950 "chbgst.f"
		    i__5 = i__ - j + ka1 + j * ab_dim1;
#line 950 "chbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 950 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 950 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 953 "chbgst.f"
/* L530: */
#line 953 "chbgst.f"
		}
#line 954 "chbgst.f"
/* L540: */
#line 954 "chbgst.f"
	    }
#line 955 "chbgst.f"
	    i__3 = i__;
#line 955 "chbgst.f"
	    for (j = i1; j <= i__3; ++j) {
/* Computing MIN */
#line 956 "chbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 956 "chbgst.f"
		i__4 = min(i__1,i__2);
#line 956 "chbgst.f"
		for (k = i__ + 1; k <= i__4; ++k) {
#line 957 "chbgst.f"
		    i__1 = j - k + ka1 + k * ab_dim1;
#line 957 "chbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 957 "chbgst.f"
		    i__5 = i__ - k + kb1 + k * bb_dim1;
#line 957 "chbgst.f"
		    i__6 = j - i__ + ka1 + i__ * ab_dim1;
#line 957 "chbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 957 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 957 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 959 "chbgst.f"
/* L550: */
#line 959 "chbgst.f"
		}
#line 960 "chbgst.f"
/* L560: */
#line 960 "chbgst.f"
	    }

#line 962 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 966 "chbgst.f"
		d__1 = 1. / bii;
#line 966 "chbgst.f"
		csscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 967 "chbgst.f"
		if (kbt > 0) {
#line 967 "chbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 967 "chbgst.f"
		    i__3 = *ldbb - 1;
#line 967 "chbgst.f"
		    cgeru_(&nx, &kbt, &z__1, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    *kb + (i__ + 1) * bb_dim1], &i__3, &x[(i__ + 1) * 
			    x_dim1 + 1], ldx);
#line 967 "chbgst.f"
		}
#line 970 "chbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 974 "chbgst.f"
	    i__3 = i1 - i__ + ka1 + i__ * ab_dim1;
#line 974 "chbgst.f"
	    ra1.r = ab[i__3].r, ra1.i = ab[i__3].i;
#line 975 "chbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 980 "chbgst.f"
	i__3 = *kb - 1;
#line 980 "chbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 981 "chbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 986 "chbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i+k-ka-1,i) */

#line 990 "chbgst.f"
		    clartg_(&ab[k + 1 + i__ * ab_dim1], &ra1, &rwork[i__ + k 
			    - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k-ka-1,i+k) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 996 "chbgst.f"
		    i__4 = kb1 - k + (i__ + k) * bb_dim1;
#line 996 "chbgst.f"
		    z__2.r = -bb[i__4].r, z__2.i = -bb[i__4].i;
#line 996 "chbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 996 "chbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 997 "chbgst.f"
		    i__4 = m - *kb + i__ + k;
#line 997 "chbgst.f"
		    i__1 = i__ + k - *ka;
#line 997 "chbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 997 "chbgst.f"
		    d_cnjg(&z__4, &work[i__ + k - *ka]);
#line 997 "chbgst.f"
		    i__2 = (i__ + k) * ab_dim1 + 1;
#line 997 "chbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 997 "chbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 997 "chbgst.f"
		    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 1000 "chbgst.f"
		    i__4 = (i__ + k) * ab_dim1 + 1;
#line 1000 "chbgst.f"
		    i__1 = i__ + k - *ka;
#line 1000 "chbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 1000 "chbgst.f"
		    i__2 = i__ + k - *ka;
#line 1000 "chbgst.f"
		    i__5 = (i__ + k) * ab_dim1 + 1;
#line 1000 "chbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 1000 "chbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 1000 "chbgst.f"
		    ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 1002 "chbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 1003 "chbgst.f"
		}
#line 1004 "chbgst.f"
	    }
/* Computing MAX */
#line 1005 "chbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1005 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;
#line 1006 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1007 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1008 "chbgst.f"
	    if (update) {
/* Computing MIN */
#line 1009 "chbgst.f"
		i__4 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1009 "chbgst.f"
		j2t = min(i__4,i__1);
#line 1010 "chbgst.f"
	    } else {
#line 1011 "chbgst.f"
		j2t = j2;
#line 1012 "chbgst.f"
	    }
#line 1013 "chbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1014 "chbgst.f"
	    i__4 = j2t;
#line 1014 "chbgst.f"
	    i__1 = ka1;
#line 1014 "chbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__4 : j <= i__4; j += i__1) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(j) */

#line 1019 "chbgst.f"
		i__2 = j;
#line 1019 "chbgst.f"
		i__5 = j;
#line 1019 "chbgst.f"
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
#line 1019 "chbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 1019 "chbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 1020 "chbgst.f"
		i__2 = (j + *ka - 1) * ab_dim1 + 1;
#line 1020 "chbgst.f"
		i__5 = j;
#line 1020 "chbgst.f"
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
#line 1020 "chbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 1020 "chbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 1021 "chbgst.f"
/* L570: */
#line 1021 "chbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1026 "chbgst.f"
	    if (nrt > 0) {
#line 1026 "chbgst.f"
		clargv_(&nrt, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[j1],
			 &ka1, &rwork[j1], &ka1);
#line 1026 "chbgst.f"
	    }
#line 1029 "chbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 1033 "chbgst.f"
		i__1 = *ka - 1;
#line 1033 "chbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1034 "chbgst.f"
		    clartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &rwork[
			    j1], &work[j1], &ka1);
#line 1037 "chbgst.f"
/* L580: */
#line 1037 "chbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1042 "chbgst.f"
		clar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &rwork[j1], 
			&work[j1], &ka1);

#line 1046 "chbgst.f"
		clacgv_(&nr, &work[j1], &ka1);
#line 1047 "chbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 1051 "chbgst.f"
	    i__1 = *kb - k + 1;
#line 1051 "chbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1052 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1053 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1054 "chbgst.f"
		if (nrt > 0) {
#line 1054 "chbgst.f"
		    clartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[j1t], &work[
			    j1t], &ka1);
#line 1054 "chbgst.f"
		}
#line 1058 "chbgst.f"
/* L590: */
#line 1058 "chbgst.f"
	    }

#line 1060 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1064 "chbgst.f"
		i__1 = j2;
#line 1064 "chbgst.f"
		i__4 = ka1;
#line 1064 "chbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__1 : j <= i__1; j += i__4) {
#line 1065 "chbgst.f"
		    crot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[j], &work[j]);
#line 1067 "chbgst.f"
/* L600: */
#line 1067 "chbgst.f"
		}
#line 1068 "chbgst.f"
	    }
#line 1069 "chbgst.f"
/* L610: */
#line 1069 "chbgst.f"
	}

#line 1071 "chbgst.f"
	if (update) {
#line 1072 "chbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt-ka-1,i+kbt) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1077 "chbgst.f"
		i__3 = m - *kb + i__ + kbt;
#line 1077 "chbgst.f"
		i__4 = kb1 - kbt + (i__ + kbt) * bb_dim1;
#line 1077 "chbgst.f"
		z__2.r = -bb[i__4].r, z__2.i = -bb[i__4].i;
#line 1077 "chbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 1077 "chbgst.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 1078 "chbgst.f"
	    }
#line 1079 "chbgst.f"
	}

#line 1081 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1082 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 1083 "chbgst.f"
		i__3 = 2, i__4 = k + i0 - m;
#line 1083 "chbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1084 "chbgst.f"
	    } else {
/* Computing MAX */
#line 1085 "chbgst.f"
		i__3 = 1, i__4 = k + i0 - m;
#line 1085 "chbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1086 "chbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 1090 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1091 "chbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1092 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1093 "chbgst.f"
		if (nrt > 0) {
#line 1093 "chbgst.f"
		    clartv_(&nrt, &ab[l + (j1t + *ka) * ab_dim1], &inca, &ab[
			    l + 1 + (j1t + *ka - 1) * ab_dim1], &inca, &rwork[
			    m - *kb + j1t + *ka], &work[m - *kb + j1t + *ka], 
			    &ka1);
#line 1093 "chbgst.f"
		}
#line 1098 "chbgst.f"
/* L620: */
#line 1098 "chbgst.f"
	    }
#line 1099 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1100 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1101 "chbgst.f"
	    i__3 = j2;
#line 1101 "chbgst.f"
	    i__4 = ka1;
#line 1101 "chbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1102 "chbgst.f"
		i__1 = m - *kb + j;
#line 1102 "chbgst.f"
		i__2 = m - *kb + j + *ka;
#line 1102 "chbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 1103 "chbgst.f"
		rwork[m - *kb + j] = rwork[m - *kb + j + *ka];
#line 1104 "chbgst.f"
/* L630: */
#line 1104 "chbgst.f"
	    }
#line 1105 "chbgst.f"
	    i__4 = j2;
#line 1105 "chbgst.f"
	    i__3 = ka1;
#line 1105 "chbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1110 "chbgst.f"
		i__1 = m - *kb + j;
#line 1110 "chbgst.f"
		i__2 = m - *kb + j;
#line 1110 "chbgst.f"
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
#line 1110 "chbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 1110 "chbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 1111 "chbgst.f"
		i__1 = (j + *ka - 1) * ab_dim1 + 1;
#line 1111 "chbgst.f"
		i__2 = m - *kb + j;
#line 1111 "chbgst.f"
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
#line 1111 "chbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 1111 "chbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1112 "chbgst.f"
/* L640: */
#line 1112 "chbgst.f"
	    }
#line 1113 "chbgst.f"
	    if (update) {
#line 1114 "chbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1114 "chbgst.f"
		    i__3 = m - *kb + i__ + k - *ka;
#line 1114 "chbgst.f"
		    i__4 = m - *kb + i__ + k;
#line 1114 "chbgst.f"
		    work[i__3].r = work[i__4].r, work[i__3].i = work[i__4].i;
#line 1114 "chbgst.f"
		}
#line 1116 "chbgst.f"
	    }
#line 1117 "chbgst.f"
/* L650: */
#line 1117 "chbgst.f"
	}

#line 1119 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1120 "chbgst.f"
	    i__3 = 1, i__4 = k + i0 - m;
#line 1120 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1121 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1122 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1123 "chbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1128 "chbgst.f"
		clargv_(&nr, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[m - *
			kb + j1], &ka1, &rwork[m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the left */

#line 1133 "chbgst.f"
		i__3 = *ka - 1;
#line 1133 "chbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 1134 "chbgst.f"
		    clartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &rwork[m 
			    - *kb + j1], &work[m - *kb + j1], &ka1);
#line 1137 "chbgst.f"
/* L660: */
#line 1137 "chbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1142 "chbgst.f"
		clar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &rwork[m - *
			kb + j1], &work[m - *kb + j1], &ka1);

#line 1146 "chbgst.f"
		clacgv_(&nr, &work[m - *kb + j1], &ka1);
#line 1147 "chbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 1151 "chbgst.f"
	    i__3 = *kb - k + 1;
#line 1151 "chbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 1152 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1153 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1154 "chbgst.f"
		if (nrt > 0) {
#line 1154 "chbgst.f"
		    clartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[m - *kb + j1t],
			     &work[m - *kb + j1t], &ka1);
#line 1154 "chbgst.f"
		}
#line 1159 "chbgst.f"
/* L670: */
#line 1159 "chbgst.f"
	    }

#line 1161 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1165 "chbgst.f"
		i__3 = j2;
#line 1165 "chbgst.f"
		i__4 = ka1;
#line 1165 "chbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1166 "chbgst.f"
		    crot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[m - *kb + j], &work[m - *kb + 
			    j]);
#line 1168 "chbgst.f"
/* L680: */
#line 1168 "chbgst.f"
		}
#line 1169 "chbgst.f"
	    }
#line 1170 "chbgst.f"
/* L690: */
#line 1170 "chbgst.f"
	}

#line 1172 "chbgst.f"
	i__4 = *kb - 1;
#line 1172 "chbgst.f"
	for (k = 1; k <= i__4; ++k) {
/* Computing MAX */
#line 1173 "chbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1173 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 1177 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1178 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1179 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1180 "chbgst.f"
		if (nrt > 0) {
#line 1180 "chbgst.f"
		    clartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[j1t], &work[
			    j1t], &ka1);
#line 1180 "chbgst.f"
		}
#line 1184 "chbgst.f"
/* L700: */
#line 1184 "chbgst.f"
	    }
#line 1185 "chbgst.f"
/* L710: */
#line 1185 "chbgst.f"
	}

#line 1187 "chbgst.f"
	if (*kb > 1) {
#line 1188 "chbgst.f"
	    i__4 = i2 - *ka;
#line 1188 "chbgst.f"
	    for (j = 2; j <= i__4; ++j) {
#line 1189 "chbgst.f"
		rwork[j] = rwork[j + *ka];
#line 1190 "chbgst.f"
		i__3 = j;
#line 1190 "chbgst.f"
		i__1 = j + *ka;
#line 1190 "chbgst.f"
		work[i__3].r = work[i__1].r, work[i__3].i = work[i__1].i;
#line 1191 "chbgst.f"
/* L720: */
#line 1191 "chbgst.f"
	    }
#line 1192 "chbgst.f"
	}

#line 1194 "chbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 1198 "chbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 1202 "chbgst.f"
	    i__4 = i__ * bb_dim1 + 1;
#line 1202 "chbgst.f"
	    bii = bb[i__4].r;
#line 1203 "chbgst.f"
	    i__4 = i__ * ab_dim1 + 1;
#line 1203 "chbgst.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 1203 "chbgst.f"
	    d__1 = ab[i__3].r / bii / bii;
#line 1203 "chbgst.f"
	    ab[i__4].r = d__1, ab[i__4].i = 0.;
#line 1204 "chbgst.f"
	    i__4 = i__ - 1;
#line 1204 "chbgst.f"
	    for (j = i1; j <= i__4; ++j) {
#line 1205 "chbgst.f"
		i__3 = i__ - j + 1 + j * ab_dim1;
#line 1205 "chbgst.f"
		i__1 = i__ - j + 1 + j * ab_dim1;
#line 1205 "chbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 1205 "chbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1206 "chbgst.f"
/* L730: */
#line 1206 "chbgst.f"
	    }
/* Computing MIN */
#line 1207 "chbgst.f"
	    i__3 = *n, i__1 = i__ + *ka;
#line 1207 "chbgst.f"
	    i__4 = min(i__3,i__1);
#line 1207 "chbgst.f"
	    for (j = i__ + 1; j <= i__4; ++j) {
#line 1208 "chbgst.f"
		i__3 = j - i__ + 1 + i__ * ab_dim1;
#line 1208 "chbgst.f"
		i__1 = j - i__ + 1 + i__ * ab_dim1;
#line 1208 "chbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 1208 "chbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1209 "chbgst.f"
/* L740: */
#line 1209 "chbgst.f"
	    }
#line 1210 "chbgst.f"
	    i__4 = i__ + kbt;
#line 1210 "chbgst.f"
	    for (k = i__ + 1; k <= i__4; ++k) {
#line 1211 "chbgst.f"
		i__3 = i__ + kbt;
#line 1211 "chbgst.f"
		for (j = k; j <= i__3; ++j) {
#line 1212 "chbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 1212 "chbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 1212 "chbgst.f"
		    i__5 = j - i__ + 1 + i__ * bb_dim1;
#line 1212 "chbgst.f"
		    d_cnjg(&z__5, &ab[k - i__ + 1 + i__ * ab_dim1]);
#line 1212 "chbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 1212 "chbgst.f"
		    z__3.r = ab[i__2].r - z__4.r, z__3.i = ab[i__2].i - 
			    z__4.i;
#line 1212 "chbgst.f"
		    d_cnjg(&z__7, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1212 "chbgst.f"
		    i__6 = j - i__ + 1 + i__ * ab_dim1;
#line 1212 "chbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 1212 "chbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 1212 "chbgst.f"
		    i__7 = i__ * ab_dim1 + 1;
#line 1212 "chbgst.f"
		    d__1 = ab[i__7].r;
#line 1212 "chbgst.f"
		    i__8 = j - i__ + 1 + i__ * bb_dim1;
#line 1212 "chbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 1212 "chbgst.f"
		    d_cnjg(&z__10, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1212 "chbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 1212 "chbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 1212 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1218 "chbgst.f"
/* L750: */
#line 1218 "chbgst.f"
		}
/* Computing MIN */
#line 1219 "chbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 1219 "chbgst.f"
		i__3 = min(i__1,i__2);
#line 1219 "chbgst.f"
		for (j = i__ + kbt + 1; j <= i__3; ++j) {
#line 1220 "chbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 1220 "chbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 1220 "chbgst.f"
		    d_cnjg(&z__3, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1220 "chbgst.f"
		    i__5 = j - i__ + 1 + i__ * ab_dim1;
#line 1220 "chbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 1220 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 1220 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1223 "chbgst.f"
/* L760: */
#line 1223 "chbgst.f"
		}
#line 1224 "chbgst.f"
/* L770: */
#line 1224 "chbgst.f"
	    }
#line 1225 "chbgst.f"
	    i__4 = i__;
#line 1225 "chbgst.f"
	    for (j = i1; j <= i__4; ++j) {
/* Computing MIN */
#line 1226 "chbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 1226 "chbgst.f"
		i__3 = min(i__1,i__2);
#line 1226 "chbgst.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 1227 "chbgst.f"
		    i__1 = k - j + 1 + j * ab_dim1;
#line 1227 "chbgst.f"
		    i__2 = k - j + 1 + j * ab_dim1;
#line 1227 "chbgst.f"
		    i__5 = k - i__ + 1 + i__ * bb_dim1;
#line 1227 "chbgst.f"
		    i__6 = i__ - j + 1 + j * ab_dim1;
#line 1227 "chbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 1227 "chbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 1227 "chbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1229 "chbgst.f"
/* L780: */
#line 1229 "chbgst.f"
		}
#line 1230 "chbgst.f"
/* L790: */
#line 1230 "chbgst.f"
	    }

#line 1232 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 1236 "chbgst.f"
		d__1 = 1. / bii;
#line 1236 "chbgst.f"
		csscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 1237 "chbgst.f"
		if (kbt > 0) {
#line 1237 "chbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 1237 "chbgst.f"
		    cgerc_(&nx, &kbt, &z__1, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    i__ * bb_dim1 + 2], &c__1, &x[(i__ + 1) * x_dim1 
			    + 1], ldx);
#line 1237 "chbgst.f"
		}
#line 1240 "chbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 1244 "chbgst.f"
	    i__4 = i__ - i1 + 1 + i1 * ab_dim1;
#line 1244 "chbgst.f"
	    ra1.r = ab[i__4].r, ra1.i = ab[i__4].i;
#line 1245 "chbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 1250 "chbgst.f"
	i__4 = *kb - 1;
#line 1250 "chbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 1251 "chbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 1256 "chbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i,i+k-ka-1) */

#line 1260 "chbgst.f"
		    clartg_(&ab[ka1 - k + (i__ + k - *ka) * ab_dim1], &ra1, &
			    rwork[i__ + k - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k,i+k-ka-1) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 1266 "chbgst.f"
		    i__3 = k + 1 + i__ * bb_dim1;
#line 1266 "chbgst.f"
		    z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 1266 "chbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 1266 "chbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 1267 "chbgst.f"
		    i__3 = m - *kb + i__ + k;
#line 1267 "chbgst.f"
		    i__1 = i__ + k - *ka;
#line 1267 "chbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 1267 "chbgst.f"
		    d_cnjg(&z__4, &work[i__ + k - *ka]);
#line 1267 "chbgst.f"
		    i__2 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1267 "chbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 1267 "chbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 1267 "chbgst.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 1270 "chbgst.f"
		    i__3 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1270 "chbgst.f"
		    i__1 = i__ + k - *ka;
#line 1270 "chbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 1270 "chbgst.f"
		    i__2 = i__ + k - *ka;
#line 1270 "chbgst.f"
		    i__5 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1270 "chbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 1270 "chbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 1270 "chbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1272 "chbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 1273 "chbgst.f"
		}
#line 1274 "chbgst.f"
	    }
/* Computing MAX */
#line 1275 "chbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1275 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;
#line 1276 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1277 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1278 "chbgst.f"
	    if (update) {
/* Computing MIN */
#line 1279 "chbgst.f"
		i__3 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1279 "chbgst.f"
		j2t = min(i__3,i__1);
#line 1280 "chbgst.f"
	    } else {
#line 1281 "chbgst.f"
		j2t = j2;
#line 1282 "chbgst.f"
	    }
#line 1283 "chbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1284 "chbgst.f"
	    i__3 = j2t;
#line 1284 "chbgst.f"
	    i__1 = ka1;
#line 1284 "chbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(j) */

#line 1289 "chbgst.f"
		i__2 = j;
#line 1289 "chbgst.f"
		i__5 = j;
#line 1289 "chbgst.f"
		i__6 = ka1 + (j - 1) * ab_dim1;
#line 1289 "chbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 1289 "chbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 1290 "chbgst.f"
		i__2 = ka1 + (j - 1) * ab_dim1;
#line 1290 "chbgst.f"
		i__5 = j;
#line 1290 "chbgst.f"
		i__6 = ka1 + (j - 1) * ab_dim1;
#line 1290 "chbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 1290 "chbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 1291 "chbgst.f"
/* L800: */
#line 1291 "chbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1296 "chbgst.f"
	    if (nrt > 0) {
#line 1296 "chbgst.f"
		clargv_(&nrt, &ab[ka1 + j1 * ab_dim1], &inca, &work[j1], &ka1,
			 &rwork[j1], &ka1);
#line 1296 "chbgst.f"
	    }
#line 1299 "chbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 1303 "chbgst.f"
		i__1 = *ka - 1;
#line 1303 "chbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1304 "chbgst.f"
		    clartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &rwork[j1], &work[
			    j1], &ka1);
#line 1306 "chbgst.f"
/* L810: */
#line 1306 "chbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1311 "chbgst.f"
		clar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &rwork[j1], &
			work[j1], &ka1);

#line 1315 "chbgst.f"
		clacgv_(&nr, &work[j1], &ka1);
#line 1316 "chbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 1320 "chbgst.f"
	    i__1 = *kb - k + 1;
#line 1320 "chbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1321 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1322 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1323 "chbgst.f"
		if (nrt > 0) {
#line 1323 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[j1t], &work[j1t], &ka1);
#line 1323 "chbgst.f"
		}
#line 1327 "chbgst.f"
/* L820: */
#line 1327 "chbgst.f"
	    }

#line 1329 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1333 "chbgst.f"
		i__1 = j2;
#line 1333 "chbgst.f"
		i__3 = ka1;
#line 1333 "chbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 1334 "chbgst.f"
		    d_cnjg(&z__1, &work[j]);
#line 1334 "chbgst.f"
		    crot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[j], &z__1);
#line 1336 "chbgst.f"
/* L830: */
#line 1336 "chbgst.f"
		}
#line 1337 "chbgst.f"
	    }
#line 1338 "chbgst.f"
/* L840: */
#line 1338 "chbgst.f"
	}

#line 1340 "chbgst.f"
	if (update) {
#line 1341 "chbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt,i+kbt-ka-1) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1346 "chbgst.f"
		i__4 = m - *kb + i__ + kbt;
#line 1346 "chbgst.f"
		i__3 = kbt + 1 + i__ * bb_dim1;
#line 1346 "chbgst.f"
		z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 1346 "chbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 1346 "chbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 1347 "chbgst.f"
	    }
#line 1348 "chbgst.f"
	}

#line 1350 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1351 "chbgst.f"
	    if (update) {
/* Computing MAX */
#line 1352 "chbgst.f"
		i__4 = 2, i__3 = k + i0 - m;
#line 1352 "chbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1353 "chbgst.f"
	    } else {
/* Computing MAX */
#line 1354 "chbgst.f"
		i__4 = 1, i__3 = k + i0 - m;
#line 1354 "chbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1355 "chbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 1359 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1360 "chbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1361 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1362 "chbgst.f"
		if (nrt > 0) {
#line 1362 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + (j1t + l - 1) * ab_dim1], 
			    &inca, &ab[ka1 - l + (j1t + l - 1) * ab_dim1], &
			    inca, &rwork[m - *kb + j1t + *ka], &work[m - *kb 
			    + j1t + *ka], &ka1);
#line 1362 "chbgst.f"
		}
#line 1367 "chbgst.f"
/* L850: */
#line 1367 "chbgst.f"
	    }
#line 1368 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1369 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1370 "chbgst.f"
	    i__4 = j2;
#line 1370 "chbgst.f"
	    i__3 = ka1;
#line 1370 "chbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1371 "chbgst.f"
		i__1 = m - *kb + j;
#line 1371 "chbgst.f"
		i__2 = m - *kb + j + *ka;
#line 1371 "chbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 1372 "chbgst.f"
		rwork[m - *kb + j] = rwork[m - *kb + j + *ka];
#line 1373 "chbgst.f"
/* L860: */
#line 1373 "chbgst.f"
	    }
#line 1374 "chbgst.f"
	    i__3 = j2;
#line 1374 "chbgst.f"
	    i__4 = ka1;
#line 1374 "chbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1379 "chbgst.f"
		i__1 = m - *kb + j;
#line 1379 "chbgst.f"
		i__2 = m - *kb + j;
#line 1379 "chbgst.f"
		i__5 = ka1 + (j - 1) * ab_dim1;
#line 1379 "chbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 1379 "chbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 1380 "chbgst.f"
		i__1 = ka1 + (j - 1) * ab_dim1;
#line 1380 "chbgst.f"
		i__2 = m - *kb + j;
#line 1380 "chbgst.f"
		i__5 = ka1 + (j - 1) * ab_dim1;
#line 1380 "chbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 1380 "chbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1381 "chbgst.f"
/* L870: */
#line 1381 "chbgst.f"
	    }
#line 1382 "chbgst.f"
	    if (update) {
#line 1383 "chbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1383 "chbgst.f"
		    i__4 = m - *kb + i__ + k - *ka;
#line 1383 "chbgst.f"
		    i__3 = m - *kb + i__ + k;
#line 1383 "chbgst.f"
		    work[i__4].r = work[i__3].r, work[i__4].i = work[i__3].i;
#line 1383 "chbgst.f"
		}
#line 1385 "chbgst.f"
	    }
#line 1386 "chbgst.f"
/* L880: */
#line 1386 "chbgst.f"
	}

#line 1388 "chbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1389 "chbgst.f"
	    i__4 = 1, i__3 = k + i0 - m;
#line 1389 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1390 "chbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1391 "chbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1392 "chbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1397 "chbgst.f"
		clargv_(&nr, &ab[ka1 + j1 * ab_dim1], &inca, &work[m - *kb + 
			j1], &ka1, &rwork[m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the right */

#line 1402 "chbgst.f"
		i__4 = *ka - 1;
#line 1402 "chbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 1403 "chbgst.f"
		    clartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &rwork[m - *kb + j1]
			    , &work[m - *kb + j1], &ka1);
#line 1406 "chbgst.f"
/* L890: */
#line 1406 "chbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1411 "chbgst.f"
		clar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &rwork[m - *
			kb + j1], &work[m - *kb + j1], &ka1);

#line 1415 "chbgst.f"
		clacgv_(&nr, &work[m - *kb + j1], &ka1);
#line 1416 "chbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 1420 "chbgst.f"
	    i__4 = *kb - k + 1;
#line 1420 "chbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 1421 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1422 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1423 "chbgst.f"
		if (nrt > 0) {
#line 1423 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[m - *kb + j1t], &work[m - *kb + 
			    j1t], &ka1);
#line 1423 "chbgst.f"
		}
#line 1428 "chbgst.f"
/* L900: */
#line 1428 "chbgst.f"
	    }

#line 1430 "chbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1434 "chbgst.f"
		i__4 = j2;
#line 1434 "chbgst.f"
		i__3 = ka1;
#line 1434 "chbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1435 "chbgst.f"
		    d_cnjg(&z__1, &work[m - *kb + j]);
#line 1435 "chbgst.f"
		    crot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[m - *kb + j], &z__1);
#line 1437 "chbgst.f"
/* L910: */
#line 1437 "chbgst.f"
		}
#line 1438 "chbgst.f"
	    }
#line 1439 "chbgst.f"
/* L920: */
#line 1439 "chbgst.f"
	}

#line 1441 "chbgst.f"
	i__3 = *kb - 1;
#line 1441 "chbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1442 "chbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1442 "chbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 1446 "chbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1447 "chbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1448 "chbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1449 "chbgst.f"
		if (nrt > 0) {
#line 1449 "chbgst.f"
		    clartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[j1t], &work[j1t], &ka1);
#line 1449 "chbgst.f"
		}
#line 1453 "chbgst.f"
/* L930: */
#line 1453 "chbgst.f"
	    }
#line 1454 "chbgst.f"
/* L940: */
#line 1454 "chbgst.f"
	}

#line 1456 "chbgst.f"
	if (*kb > 1) {
#line 1457 "chbgst.f"
	    i__3 = i2 - *ka;
#line 1457 "chbgst.f"
	    for (j = 2; j <= i__3; ++j) {
#line 1458 "chbgst.f"
		rwork[j] = rwork[j + *ka];
#line 1459 "chbgst.f"
		i__4 = j;
#line 1459 "chbgst.f"
		i__1 = j + *ka;
#line 1459 "chbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 1460 "chbgst.f"
/* L950: */
#line 1460 "chbgst.f"
	    }
#line 1461 "chbgst.f"
	}

#line 1463 "chbgst.f"
    }

#line 1465 "chbgst.f"
    goto L490;

/*     End of CHBGST */

} /* chbgst_ */

