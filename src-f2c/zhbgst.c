#line 1 "zhbgst.f"
/* zhbgst.f -- translated by f2c (version 20100827).
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

#line 1 "zhbgst.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHBGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X, */
/*                          LDX, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBGST reduces a complex Hermitian-definite banded generalized */
/* > eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y, */
/* > such that C has the same bandwidth as A. */
/* > */
/* > B must have been previously factorized as S**H*S by ZPBSTF, using a */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          BB is COMPLEX*16 array, dimension (LDBB,N) */
/* >          The banded factor S from the split Cholesky factorization of */
/* >          B, as returned by ZPBSTF, stored in the first kb+1 rows of */
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
/* >          X is COMPLEX*16 array, dimension (LDX,N) */
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
/* >          WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhbgst_(char *vect, char *uplo, integer *n, integer *ka, 
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
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int zgeru_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical wantx;
    extern /* Subroutine */ int zlar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *), xerbla_(char *, integer *, ftnlen), 
	    zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    static logical update;
    extern /* Subroutine */ int zlacgv_(integer *, doublecomplex *, integer *)
	    , zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zlartg_(
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *), zlargv_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *), zlartv_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, doublecomplex *, integer *);


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

#line 213 "zhbgst.f"
    /* Parameter adjustments */
#line 213 "zhbgst.f"
    ab_dim1 = *ldab;
#line 213 "zhbgst.f"
    ab_offset = 1 + ab_dim1;
#line 213 "zhbgst.f"
    ab -= ab_offset;
#line 213 "zhbgst.f"
    bb_dim1 = *ldbb;
#line 213 "zhbgst.f"
    bb_offset = 1 + bb_dim1;
#line 213 "zhbgst.f"
    bb -= bb_offset;
#line 213 "zhbgst.f"
    x_dim1 = *ldx;
#line 213 "zhbgst.f"
    x_offset = 1 + x_dim1;
#line 213 "zhbgst.f"
    x -= x_offset;
#line 213 "zhbgst.f"
    --work;
#line 213 "zhbgst.f"
    --rwork;
#line 213 "zhbgst.f"

#line 213 "zhbgst.f"
    /* Function Body */
#line 213 "zhbgst.f"
    wantx = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 214 "zhbgst.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 215 "zhbgst.f"
    ka1 = *ka + 1;
#line 216 "zhbgst.f"
    kb1 = *kb + 1;
#line 217 "zhbgst.f"
    *info = 0;
#line 218 "zhbgst.f"
    if (! wantx && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 219 "zhbgst.f"
	*info = -1;
#line 220 "zhbgst.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 221 "zhbgst.f"
	*info = -2;
#line 222 "zhbgst.f"
    } else if (*n < 0) {
#line 223 "zhbgst.f"
	*info = -3;
#line 224 "zhbgst.f"
    } else if (*ka < 0) {
#line 225 "zhbgst.f"
	*info = -4;
#line 226 "zhbgst.f"
    } else if (*kb < 0 || *kb > *ka) {
#line 227 "zhbgst.f"
	*info = -5;
#line 228 "zhbgst.f"
    } else if (*ldab < *ka + 1) {
#line 229 "zhbgst.f"
	*info = -7;
#line 230 "zhbgst.f"
    } else if (*ldbb < *kb + 1) {
#line 231 "zhbgst.f"
	*info = -9;
#line 232 "zhbgst.f"
    } else if (*ldx < 1 || wantx && *ldx < max(1,*n)) {
#line 233 "zhbgst.f"
	*info = -11;
#line 234 "zhbgst.f"
    }
#line 235 "zhbgst.f"
    if (*info != 0) {
#line 236 "zhbgst.f"
	i__1 = -(*info);
#line 236 "zhbgst.f"
	xerbla_("ZHBGST", &i__1, (ftnlen)6);
#line 237 "zhbgst.f"
	return 0;
#line 238 "zhbgst.f"
    }

/*     Quick return if possible */

#line 242 "zhbgst.f"
    if (*n == 0) {
#line 242 "zhbgst.f"
	return 0;
#line 242 "zhbgst.f"
    }

#line 245 "zhbgst.f"
    inca = *ldab * ka1;

/*     Initialize X to the unit matrix, if needed */

#line 249 "zhbgst.f"
    if (wantx) {
#line 249 "zhbgst.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &x[x_offset], ldx, (ftnlen)4);
#line 249 "zhbgst.f"
    }

/*     Set M to the splitting point m. It must be the same value as is */
/*     used in ZPBSTF. The chosen value allows the arrays WORK and RWORK */
/*     to be of dimension (N). */

#line 256 "zhbgst.f"
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

#line 317 "zhbgst.f"
    update = TRUE_;
#line 318 "zhbgst.f"
    i__ = *n + 1;
#line 319 "zhbgst.f"
L10:
#line 320 "zhbgst.f"
    if (update) {
#line 321 "zhbgst.f"
	--i__;
/* Computing MIN */
#line 322 "zhbgst.f"
	i__1 = *kb, i__2 = i__ - 1;
#line 322 "zhbgst.f"
	kbt = min(i__1,i__2);
#line 323 "zhbgst.f"
	i0 = i__ - 1;
/* Computing MIN */
#line 324 "zhbgst.f"
	i__1 = *n, i__2 = i__ + *ka;
#line 324 "zhbgst.f"
	i1 = min(i__1,i__2);
#line 325 "zhbgst.f"
	i2 = i__ - kbt + ka1;
#line 326 "zhbgst.f"
	if (i__ < m + 1) {
#line 327 "zhbgst.f"
	    update = FALSE_;
#line 328 "zhbgst.f"
	    ++i__;
#line 329 "zhbgst.f"
	    i0 = m;
#line 330 "zhbgst.f"
	    if (*ka == 0) {
#line 330 "zhbgst.f"
		goto L480;
#line 330 "zhbgst.f"
	    }
#line 332 "zhbgst.f"
	    goto L10;
#line 333 "zhbgst.f"
	}
#line 334 "zhbgst.f"
    } else {
#line 335 "zhbgst.f"
	i__ += *ka;
#line 336 "zhbgst.f"
	if (i__ > *n - 1) {
#line 336 "zhbgst.f"
	    goto L480;
#line 336 "zhbgst.f"
	}
#line 338 "zhbgst.f"
    }

#line 340 "zhbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 344 "zhbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 348 "zhbgst.f"
	    i__1 = kb1 + i__ * bb_dim1;
#line 348 "zhbgst.f"
	    bii = bb[i__1].r;
#line 349 "zhbgst.f"
	    i__1 = ka1 + i__ * ab_dim1;
#line 349 "zhbgst.f"
	    i__2 = ka1 + i__ * ab_dim1;
#line 349 "zhbgst.f"
	    d__1 = ab[i__2].r / bii / bii;
#line 349 "zhbgst.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 350 "zhbgst.f"
	    i__1 = i1;
#line 350 "zhbgst.f"
	    for (j = i__ + 1; j <= i__1; ++j) {
#line 351 "zhbgst.f"
		i__2 = i__ - j + ka1 + j * ab_dim1;
#line 351 "zhbgst.f"
		i__3 = i__ - j + ka1 + j * ab_dim1;
#line 351 "zhbgst.f"
		z__1.r = ab[i__3].r / bii, z__1.i = ab[i__3].i / bii;
#line 351 "zhbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 352 "zhbgst.f"
/* L20: */
#line 352 "zhbgst.f"
	    }
/* Computing MAX */
#line 353 "zhbgst.f"
	    i__1 = 1, i__2 = i__ - *ka;
#line 353 "zhbgst.f"
	    i__3 = i__ - 1;
#line 353 "zhbgst.f"
	    for (j = max(i__1,i__2); j <= i__3; ++j) {
#line 354 "zhbgst.f"
		i__1 = j - i__ + ka1 + i__ * ab_dim1;
#line 354 "zhbgst.f"
		i__2 = j - i__ + ka1 + i__ * ab_dim1;
#line 354 "zhbgst.f"
		z__1.r = ab[i__2].r / bii, z__1.i = ab[i__2].i / bii;
#line 354 "zhbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 355 "zhbgst.f"
/* L30: */
#line 355 "zhbgst.f"
	    }
#line 356 "zhbgst.f"
	    i__3 = i__ - 1;
#line 356 "zhbgst.f"
	    for (k = i__ - kbt; k <= i__3; ++k) {
#line 357 "zhbgst.f"
		i__1 = k;
#line 357 "zhbgst.f"
		for (j = i__ - kbt; j <= i__1; ++j) {
#line 358 "zhbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 358 "zhbgst.f"
		    i__4 = j - k + ka1 + k * ab_dim1;
#line 358 "zhbgst.f"
		    i__5 = j - i__ + kb1 + i__ * bb_dim1;
#line 358 "zhbgst.f"
		    d_cnjg(&z__5, &ab[k - i__ + ka1 + i__ * ab_dim1]);
#line 358 "zhbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 358 "zhbgst.f"
		    z__3.r = ab[i__4].r - z__4.r, z__3.i = ab[i__4].i - 
			    z__4.i;
#line 358 "zhbgst.f"
		    d_cnjg(&z__7, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 358 "zhbgst.f"
		    i__6 = j - i__ + ka1 + i__ * ab_dim1;
#line 358 "zhbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 358 "zhbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 358 "zhbgst.f"
		    i__7 = ka1 + i__ * ab_dim1;
#line 358 "zhbgst.f"
		    d__1 = ab[i__7].r;
#line 358 "zhbgst.f"
		    i__8 = j - i__ + kb1 + i__ * bb_dim1;
#line 358 "zhbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 358 "zhbgst.f"
		    d_cnjg(&z__10, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 358 "zhbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 358 "zhbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 358 "zhbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 366 "zhbgst.f"
/* L40: */
#line 366 "zhbgst.f"
		}
/* Computing MAX */
#line 367 "zhbgst.f"
		i__1 = 1, i__2 = i__ - *ka;
#line 367 "zhbgst.f"
		i__4 = i__ - kbt - 1;
#line 367 "zhbgst.f"
		for (j = max(i__1,i__2); j <= i__4; ++j) {
#line 368 "zhbgst.f"
		    i__1 = j - k + ka1 + k * ab_dim1;
#line 368 "zhbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 368 "zhbgst.f"
		    d_cnjg(&z__3, &bb[k - i__ + kb1 + i__ * bb_dim1]);
#line 368 "zhbgst.f"
		    i__5 = j - i__ + ka1 + i__ * ab_dim1;
#line 368 "zhbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 368 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 368 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 371 "zhbgst.f"
/* L50: */
#line 371 "zhbgst.f"
		}
#line 372 "zhbgst.f"
/* L60: */
#line 372 "zhbgst.f"
	    }
#line 373 "zhbgst.f"
	    i__3 = i1;
#line 373 "zhbgst.f"
	    for (j = i__; j <= i__3; ++j) {
/* Computing MAX */
#line 374 "zhbgst.f"
		i__4 = j - *ka, i__1 = i__ - kbt;
#line 374 "zhbgst.f"
		i__2 = i__ - 1;
#line 374 "zhbgst.f"
		for (k = max(i__4,i__1); k <= i__2; ++k) {
#line 375 "zhbgst.f"
		    i__4 = k - j + ka1 + j * ab_dim1;
#line 375 "zhbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 375 "zhbgst.f"
		    i__5 = k - i__ + kb1 + i__ * bb_dim1;
#line 375 "zhbgst.f"
		    i__6 = i__ - j + ka1 + j * ab_dim1;
#line 375 "zhbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 375 "zhbgst.f"
		    z__1.r = ab[i__1].r - z__2.r, z__1.i = ab[i__1].i - 
			    z__2.i;
#line 375 "zhbgst.f"
		    ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 377 "zhbgst.f"
/* L70: */
#line 377 "zhbgst.f"
		}
#line 378 "zhbgst.f"
/* L80: */
#line 378 "zhbgst.f"
	    }

#line 380 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 384 "zhbgst.f"
		i__3 = *n - m;
#line 384 "zhbgst.f"
		d__1 = 1. / bii;
#line 384 "zhbgst.f"
		zdscal_(&i__3, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 385 "zhbgst.f"
		if (kbt > 0) {
#line 385 "zhbgst.f"
		    i__3 = *n - m;
#line 385 "zhbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 385 "zhbgst.f"
		    zgerc_(&i__3, &kbt, &z__1, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kb1 - kbt + i__ * bb_dim1], &c__1, &x[m 
			    + 1 + (i__ - kbt) * x_dim1], ldx);
#line 385 "zhbgst.f"
		}
#line 389 "zhbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 393 "zhbgst.f"
	    i__3 = i__ - i1 + ka1 + i1 * ab_dim1;
#line 393 "zhbgst.f"
	    ra1.r = ab[i__3].r, ra1.i = ab[i__3].i;
#line 394 "zhbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 400 "zhbgst.f"
	i__3 = *kb - 1;
#line 400 "zhbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 401 "zhbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 406 "zhbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i,i-k+ka+1) */

#line 410 "zhbgst.f"
		    zlartg_(&ab[k + 1 + (i__ - k + *ka) * ab_dim1], &ra1, &
			    rwork[i__ - k + *ka - m], &work[i__ - k + *ka - m]
			    , &ra);

/*                 create nonzero element a(i-k,i-k+ka+1) outside the */
/*                 band and store it in WORK(i-k) */

#line 416 "zhbgst.f"
		    i__2 = kb1 - k + i__ * bb_dim1;
#line 416 "zhbgst.f"
		    z__2.r = -bb[i__2].r, z__2.i = -bb[i__2].i;
#line 416 "zhbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 416 "zhbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 417 "zhbgst.f"
		    i__2 = i__ - k;
#line 417 "zhbgst.f"
		    i__4 = i__ - k + *ka - m;
#line 417 "zhbgst.f"
		    z__2.r = rwork[i__4] * t.r, z__2.i = rwork[i__4] * t.i;
#line 417 "zhbgst.f"
		    d_cnjg(&z__4, &work[i__ - k + *ka - m]);
#line 417 "zhbgst.f"
		    i__1 = (i__ - k + *ka) * ab_dim1 + 1;
#line 417 "zhbgst.f"
		    z__3.r = z__4.r * ab[i__1].r - z__4.i * ab[i__1].i, 
			    z__3.i = z__4.r * ab[i__1].i + z__4.i * ab[i__1]
			    .r;
#line 417 "zhbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 417 "zhbgst.f"
		    work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 420 "zhbgst.f"
		    i__2 = (i__ - k + *ka) * ab_dim1 + 1;
#line 420 "zhbgst.f"
		    i__4 = i__ - k + *ka - m;
#line 420 "zhbgst.f"
		    z__2.r = work[i__4].r * t.r - work[i__4].i * t.i, z__2.i =
			     work[i__4].r * t.i + work[i__4].i * t.r;
#line 420 "zhbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 420 "zhbgst.f"
		    i__5 = (i__ - k + *ka) * ab_dim1 + 1;
#line 420 "zhbgst.f"
		    z__3.r = rwork[i__1] * ab[i__5].r, z__3.i = rwork[i__1] * 
			    ab[i__5].i;
#line 420 "zhbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 420 "zhbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 422 "zhbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 423 "zhbgst.f"
		}
#line 424 "zhbgst.f"
	    }
/* Computing MAX */
#line 425 "zhbgst.f"
	    i__2 = 1, i__4 = k - i0 + 2;
#line 425 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__2,i__4) * ka1;
#line 426 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 427 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 428 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 429 "zhbgst.f"
		i__2 = j2, i__4 = i__ + (*ka << 1) - k + 1;
#line 429 "zhbgst.f"
		j2t = max(i__2,i__4);
#line 430 "zhbgst.f"
	    } else {
#line 431 "zhbgst.f"
		j2t = j2;
#line 432 "zhbgst.f"
	    }
#line 433 "zhbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 434 "zhbgst.f"
	    i__2 = j1;
#line 434 "zhbgst.f"
	    i__4 = ka1;
#line 434 "zhbgst.f"
	    for (j = j2t; i__4 < 0 ? j >= i__2 : j <= i__2; j += i__4) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j-m) */

#line 439 "zhbgst.f"
		i__1 = j - m;
#line 439 "zhbgst.f"
		i__5 = j - m;
#line 439 "zhbgst.f"
		i__6 = (j + 1) * ab_dim1 + 1;
#line 439 "zhbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 439 "zhbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 440 "zhbgst.f"
		i__1 = (j + 1) * ab_dim1 + 1;
#line 440 "zhbgst.f"
		i__5 = j - m;
#line 440 "zhbgst.f"
		i__6 = (j + 1) * ab_dim1 + 1;
#line 440 "zhbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 440 "zhbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 441 "zhbgst.f"
/* L90: */
#line 441 "zhbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 446 "zhbgst.f"
	    if (nrt > 0) {
#line 446 "zhbgst.f"
		zlargv_(&nrt, &ab[j2t * ab_dim1 + 1], &inca, &work[j2t - m], &
			ka1, &rwork[j2t - m], &ka1);
#line 446 "zhbgst.f"
	    }
#line 449 "zhbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 453 "zhbgst.f"
		i__4 = *ka - 1;
#line 453 "zhbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 454 "zhbgst.f"
		    zlartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - m], 
			    &work[j2 - m], &ka1);
#line 457 "zhbgst.f"
/* L100: */
#line 457 "zhbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 462 "zhbgst.f"
		zlar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &
			rwork[j2 - m], &work[j2 - m], &ka1);

#line 466 "zhbgst.f"
		zlacgv_(&nr, &work[j2 - m], &ka1);
#line 467 "zhbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 471 "zhbgst.f"
	    i__4 = *kb - k + 1;
#line 471 "zhbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 472 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 473 "zhbgst.f"
		if (nrt > 0) {
#line 473 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2 - m], &work[j2 - m], &ka1);
#line 473 "zhbgst.f"
		}
#line 477 "zhbgst.f"
/* L110: */
#line 477 "zhbgst.f"
	    }

#line 479 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 483 "zhbgst.f"
		i__4 = j1;
#line 483 "zhbgst.f"
		i__2 = ka1;
#line 483 "zhbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__4 : j <= i__4; j += i__2) {
#line 484 "zhbgst.f"
		    i__1 = *n - m;
#line 484 "zhbgst.f"
		    d_cnjg(&z__1, &work[j - m]);
#line 484 "zhbgst.f"
		    zrot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j - m], &z__1);
#line 486 "zhbgst.f"
/* L120: */
#line 486 "zhbgst.f"
		}
#line 487 "zhbgst.f"
	    }
#line 488 "zhbgst.f"
/* L130: */
#line 488 "zhbgst.f"
	}

#line 490 "zhbgst.f"
	if (update) {
#line 491 "zhbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt,i-kbt+ka+1) outside the */
/*              band and store it in WORK(i-kbt) */

#line 496 "zhbgst.f"
		i__3 = i__ - kbt;
#line 496 "zhbgst.f"
		i__2 = kb1 - kbt + i__ * bb_dim1;
#line 496 "zhbgst.f"
		z__2.r = -bb[i__2].r, z__2.i = -bb[i__2].i;
#line 496 "zhbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 496 "zhbgst.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 497 "zhbgst.f"
	    }
#line 498 "zhbgst.f"
	}

#line 500 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 501 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 502 "zhbgst.f"
		i__3 = 2, i__2 = k - i0 + 1;
#line 502 "zhbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 503 "zhbgst.f"
	    } else {
/* Computing MAX */
#line 504 "zhbgst.f"
		i__3 = 1, i__2 = k - i0 + 1;
#line 504 "zhbgst.f"
		j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 505 "zhbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 509 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 510 "zhbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 511 "zhbgst.f"
		if (nrt > 0) {
#line 511 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + (j2 - l + 1) * ab_dim1], &inca, &ab[
			    l + 1 + (j2 - l + 1) * ab_dim1], &inca, &rwork[j2 
			    - *ka], &work[j2 - *ka], &ka1);
#line 511 "zhbgst.f"
		}
#line 515 "zhbgst.f"
/* L140: */
#line 515 "zhbgst.f"
	    }
#line 516 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 517 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 518 "zhbgst.f"
	    i__3 = j2;
#line 518 "zhbgst.f"
	    i__2 = -ka1;
#line 518 "zhbgst.f"
	    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 519 "zhbgst.f"
		i__4 = j;
#line 519 "zhbgst.f"
		i__1 = j - *ka;
#line 519 "zhbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 520 "zhbgst.f"
		rwork[j] = rwork[j - *ka];
#line 521 "zhbgst.f"
/* L150: */
#line 521 "zhbgst.f"
	    }
#line 522 "zhbgst.f"
	    i__2 = j1;
#line 522 "zhbgst.f"
	    i__3 = ka1;
#line 522 "zhbgst.f"
	    for (j = j2; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {

/*              create nonzero element a(j-ka,j+1) outside the band */
/*              and store it in WORK(j) */

#line 527 "zhbgst.f"
		i__4 = j;
#line 527 "zhbgst.f"
		i__1 = j;
#line 527 "zhbgst.f"
		i__5 = (j + 1) * ab_dim1 + 1;
#line 527 "zhbgst.f"
		z__1.r = work[i__1].r * ab[i__5].r - work[i__1].i * ab[i__5]
			.i, z__1.i = work[i__1].r * ab[i__5].i + work[i__1].i 
			* ab[i__5].r;
#line 527 "zhbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 528 "zhbgst.f"
		i__4 = (j + 1) * ab_dim1 + 1;
#line 528 "zhbgst.f"
		i__1 = j;
#line 528 "zhbgst.f"
		i__5 = (j + 1) * ab_dim1 + 1;
#line 528 "zhbgst.f"
		z__1.r = rwork[i__1] * ab[i__5].r, z__1.i = rwork[i__1] * ab[
			i__5].i;
#line 528 "zhbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 529 "zhbgst.f"
/* L160: */
#line 529 "zhbgst.f"
	    }
#line 530 "zhbgst.f"
	    if (update) {
#line 531 "zhbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 531 "zhbgst.f"
		    i__3 = i__ - k + *ka;
#line 531 "zhbgst.f"
		    i__2 = i__ - k;
#line 531 "zhbgst.f"
		    work[i__3].r = work[i__2].r, work[i__3].i = work[i__2].i;
#line 531 "zhbgst.f"
		}
#line 533 "zhbgst.f"
	    }
#line 534 "zhbgst.f"
/* L170: */
#line 534 "zhbgst.f"
	}

#line 536 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 537 "zhbgst.f"
	    i__3 = 1, i__2 = k - i0 + 1;
#line 537 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__2) * ka1;
#line 538 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 539 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 540 "zhbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 545 "zhbgst.f"
		zlargv_(&nr, &ab[j2 * ab_dim1 + 1], &inca, &work[j2], &ka1, &
			rwork[j2], &ka1);

/*              apply rotations in 2nd set from the right */

#line 550 "zhbgst.f"
		i__3 = *ka - 1;
#line 550 "zhbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 551 "zhbgst.f"
		    zlartv_(&nr, &ab[ka1 - l + j2 * ab_dim1], &inca, &ab[*ka 
			    - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2], &
			    work[j2], &ka1);
#line 554 "zhbgst.f"
/* L180: */
#line 554 "zhbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 559 "zhbgst.f"
		zlar2v_(&nr, &ab[ka1 + j2 * ab_dim1], &ab[ka1 + (j2 + 1) * 
			ab_dim1], &ab[*ka + (j2 + 1) * ab_dim1], &inca, &
			rwork[j2], &work[j2], &ka1);

#line 563 "zhbgst.f"
		zlacgv_(&nr, &work[j2], &ka1);
#line 564 "zhbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 568 "zhbgst.f"
	    i__3 = *kb - k + 1;
#line 568 "zhbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 569 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 570 "zhbgst.f"
		if (nrt > 0) {
#line 570 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2], &work[j2], &ka1);
#line 570 "zhbgst.f"
		}
#line 574 "zhbgst.f"
/* L190: */
#line 574 "zhbgst.f"
	    }

#line 576 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 580 "zhbgst.f"
		i__3 = j1;
#line 580 "zhbgst.f"
		i__2 = ka1;
#line 580 "zhbgst.f"
		for (j = j2; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
#line 581 "zhbgst.f"
		    i__4 = *n - m;
#line 581 "zhbgst.f"
		    d_cnjg(&z__1, &work[j]);
#line 581 "zhbgst.f"
		    zrot_(&i__4, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j], &z__1);
#line 583 "zhbgst.f"
/* L200: */
#line 583 "zhbgst.f"
		}
#line 584 "zhbgst.f"
	    }
#line 585 "zhbgst.f"
/* L210: */
#line 585 "zhbgst.f"
	}

#line 587 "zhbgst.f"
	i__2 = *kb - 1;
#line 587 "zhbgst.f"
	for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
#line 588 "zhbgst.f"
	    i__3 = 1, i__4 = k - i0 + 2;
#line 588 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__4) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 592 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 593 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 594 "zhbgst.f"
		if (nrt > 0) {
#line 594 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + (j2 + ka1 - l) * ab_dim1], &inca, &
			    ab[l + 1 + (j2 + ka1 - l) * ab_dim1], &inca, &
			    rwork[j2 - m], &work[j2 - m], &ka1);
#line 594 "zhbgst.f"
		}
#line 598 "zhbgst.f"
/* L220: */
#line 598 "zhbgst.f"
	    }
#line 599 "zhbgst.f"
/* L230: */
#line 599 "zhbgst.f"
	}

#line 601 "zhbgst.f"
	if (*kb > 1) {
#line 602 "zhbgst.f"
	    i__2 = j2 + *ka;
#line 602 "zhbgst.f"
	    for (j = *n - 1; j >= i__2; --j) {
#line 603 "zhbgst.f"
		rwork[j - m] = rwork[j - *ka - m];
#line 604 "zhbgst.f"
		i__3 = j - m;
#line 604 "zhbgst.f"
		i__4 = j - *ka - m;
#line 604 "zhbgst.f"
		work[i__3].r = work[i__4].r, work[i__3].i = work[i__4].i;
#line 605 "zhbgst.f"
/* L240: */
#line 605 "zhbgst.f"
	    }
#line 606 "zhbgst.f"
	}

#line 608 "zhbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 612 "zhbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 616 "zhbgst.f"
	    i__2 = i__ * bb_dim1 + 1;
#line 616 "zhbgst.f"
	    bii = bb[i__2].r;
#line 617 "zhbgst.f"
	    i__2 = i__ * ab_dim1 + 1;
#line 617 "zhbgst.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 617 "zhbgst.f"
	    d__1 = ab[i__3].r / bii / bii;
#line 617 "zhbgst.f"
	    ab[i__2].r = d__1, ab[i__2].i = 0.;
#line 618 "zhbgst.f"
	    i__2 = i1;
#line 618 "zhbgst.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 619 "zhbgst.f"
		i__3 = j - i__ + 1 + i__ * ab_dim1;
#line 619 "zhbgst.f"
		i__4 = j - i__ + 1 + i__ * ab_dim1;
#line 619 "zhbgst.f"
		z__1.r = ab[i__4].r / bii, z__1.i = ab[i__4].i / bii;
#line 619 "zhbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 620 "zhbgst.f"
/* L250: */
#line 620 "zhbgst.f"
	    }
/* Computing MAX */
#line 621 "zhbgst.f"
	    i__2 = 1, i__3 = i__ - *ka;
#line 621 "zhbgst.f"
	    i__4 = i__ - 1;
#line 621 "zhbgst.f"
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 622 "zhbgst.f"
		i__2 = i__ - j + 1 + j * ab_dim1;
#line 622 "zhbgst.f"
		i__3 = i__ - j + 1 + j * ab_dim1;
#line 622 "zhbgst.f"
		z__1.r = ab[i__3].r / bii, z__1.i = ab[i__3].i / bii;
#line 622 "zhbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 623 "zhbgst.f"
/* L260: */
#line 623 "zhbgst.f"
	    }
#line 624 "zhbgst.f"
	    i__4 = i__ - 1;
#line 624 "zhbgst.f"
	    for (k = i__ - kbt; k <= i__4; ++k) {
#line 625 "zhbgst.f"
		i__2 = k;
#line 625 "zhbgst.f"
		for (j = i__ - kbt; j <= i__2; ++j) {
#line 626 "zhbgst.f"
		    i__3 = k - j + 1 + j * ab_dim1;
#line 626 "zhbgst.f"
		    i__1 = k - j + 1 + j * ab_dim1;
#line 626 "zhbgst.f"
		    i__5 = i__ - j + 1 + j * bb_dim1;
#line 626 "zhbgst.f"
		    d_cnjg(&z__5, &ab[i__ - k + 1 + k * ab_dim1]);
#line 626 "zhbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 626 "zhbgst.f"
		    z__3.r = ab[i__1].r - z__4.r, z__3.i = ab[i__1].i - 
			    z__4.i;
#line 626 "zhbgst.f"
		    d_cnjg(&z__7, &bb[i__ - k + 1 + k * bb_dim1]);
#line 626 "zhbgst.f"
		    i__6 = i__ - j + 1 + j * ab_dim1;
#line 626 "zhbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 626 "zhbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 626 "zhbgst.f"
		    i__7 = i__ * ab_dim1 + 1;
#line 626 "zhbgst.f"
		    d__1 = ab[i__7].r;
#line 626 "zhbgst.f"
		    i__8 = i__ - j + 1 + j * bb_dim1;
#line 626 "zhbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 626 "zhbgst.f"
		    d_cnjg(&z__10, &bb[i__ - k + 1 + k * bb_dim1]);
#line 626 "zhbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 626 "zhbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 626 "zhbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 632 "zhbgst.f"
/* L270: */
#line 632 "zhbgst.f"
		}
/* Computing MAX */
#line 633 "zhbgst.f"
		i__2 = 1, i__3 = i__ - *ka;
#line 633 "zhbgst.f"
		i__1 = i__ - kbt - 1;
#line 633 "zhbgst.f"
		for (j = max(i__2,i__3); j <= i__1; ++j) {
#line 634 "zhbgst.f"
		    i__2 = k - j + 1 + j * ab_dim1;
#line 634 "zhbgst.f"
		    i__3 = k - j + 1 + j * ab_dim1;
#line 634 "zhbgst.f"
		    d_cnjg(&z__3, &bb[i__ - k + 1 + k * bb_dim1]);
#line 634 "zhbgst.f"
		    i__5 = i__ - j + 1 + j * ab_dim1;
#line 634 "zhbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 634 "zhbgst.f"
		    z__1.r = ab[i__3].r - z__2.r, z__1.i = ab[i__3].i - 
			    z__2.i;
#line 634 "zhbgst.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 637 "zhbgst.f"
/* L280: */
#line 637 "zhbgst.f"
		}
#line 638 "zhbgst.f"
/* L290: */
#line 638 "zhbgst.f"
	    }
#line 639 "zhbgst.f"
	    i__4 = i1;
#line 639 "zhbgst.f"
	    for (j = i__; j <= i__4; ++j) {
/* Computing MAX */
#line 640 "zhbgst.f"
		i__1 = j - *ka, i__2 = i__ - kbt;
#line 640 "zhbgst.f"
		i__3 = i__ - 1;
#line 640 "zhbgst.f"
		for (k = max(i__1,i__2); k <= i__3; ++k) {
#line 641 "zhbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 641 "zhbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 641 "zhbgst.f"
		    i__5 = i__ - k + 1 + k * bb_dim1;
#line 641 "zhbgst.f"
		    i__6 = j - i__ + 1 + i__ * ab_dim1;
#line 641 "zhbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 641 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 641 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 643 "zhbgst.f"
/* L300: */
#line 643 "zhbgst.f"
		}
#line 644 "zhbgst.f"
/* L310: */
#line 644 "zhbgst.f"
	    }

#line 646 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 650 "zhbgst.f"
		i__4 = *n - m;
#line 650 "zhbgst.f"
		d__1 = 1. / bii;
#line 650 "zhbgst.f"
		zdscal_(&i__4, &d__1, &x[m + 1 + i__ * x_dim1], &c__1);
#line 651 "zhbgst.f"
		if (kbt > 0) {
#line 651 "zhbgst.f"
		    i__4 = *n - m;
#line 651 "zhbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 651 "zhbgst.f"
		    i__3 = *ldbb - 1;
#line 651 "zhbgst.f"
		    zgeru_(&i__4, &kbt, &z__1, &x[m + 1 + i__ * x_dim1], &
			    c__1, &bb[kbt + 1 + (i__ - kbt) * bb_dim1], &i__3,
			     &x[m + 1 + (i__ - kbt) * x_dim1], ldx);
#line 651 "zhbgst.f"
		}
#line 655 "zhbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 659 "zhbgst.f"
	    i__4 = i1 - i__ + 1 + i__ * ab_dim1;
#line 659 "zhbgst.f"
	    ra1.r = ab[i__4].r, ra1.i = ab[i__4].i;
#line 660 "zhbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions down toward the bottom of the */
/*        band */

#line 666 "zhbgst.f"
	i__4 = *kb - 1;
#line 666 "zhbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 667 "zhbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 672 "zhbgst.f"
		if (i__ - k + *ka < *n && i__ - k > 1) {

/*                 generate rotation to annihilate a(i-k+ka+1,i) */

#line 676 "zhbgst.f"
		    zlartg_(&ab[ka1 - k + i__ * ab_dim1], &ra1, &rwork[i__ - 
			    k + *ka - m], &work[i__ - k + *ka - m], &ra);

/*                 create nonzero element a(i-k+ka+1,i-k) outside the */
/*                 band and store it in WORK(i-k) */

#line 682 "zhbgst.f"
		    i__3 = k + 1 + (i__ - k) * bb_dim1;
#line 682 "zhbgst.f"
		    z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 682 "zhbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 682 "zhbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 683 "zhbgst.f"
		    i__3 = i__ - k;
#line 683 "zhbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 683 "zhbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 683 "zhbgst.f"
		    d_cnjg(&z__4, &work[i__ - k + *ka - m]);
#line 683 "zhbgst.f"
		    i__2 = ka1 + (i__ - k) * ab_dim1;
#line 683 "zhbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 683 "zhbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 683 "zhbgst.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 686 "zhbgst.f"
		    i__3 = ka1 + (i__ - k) * ab_dim1;
#line 686 "zhbgst.f"
		    i__1 = i__ - k + *ka - m;
#line 686 "zhbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 686 "zhbgst.f"
		    i__2 = i__ - k + *ka - m;
#line 686 "zhbgst.f"
		    i__5 = ka1 + (i__ - k) * ab_dim1;
#line 686 "zhbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 686 "zhbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 686 "zhbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 688 "zhbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 689 "zhbgst.f"
		}
#line 690 "zhbgst.f"
	    }
/* Computing MAX */
#line 691 "zhbgst.f"
	    i__3 = 1, i__1 = k - i0 + 2;
#line 691 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__3,i__1) * ka1;
#line 692 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 693 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 694 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 695 "zhbgst.f"
		i__3 = j2, i__1 = i__ + (*ka << 1) - k + 1;
#line 695 "zhbgst.f"
		j2t = max(i__3,i__1);
#line 696 "zhbgst.f"
	    } else {
#line 697 "zhbgst.f"
		j2t = j2;
#line 698 "zhbgst.f"
	    }
#line 699 "zhbgst.f"
	    nrt = (*n - j2t + *ka) / ka1;
#line 700 "zhbgst.f"
	    i__3 = j1;
#line 700 "zhbgst.f"
	    i__1 = ka1;
#line 700 "zhbgst.f"
	    for (j = j2t; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j-m) */

#line 705 "zhbgst.f"
		i__2 = j - m;
#line 705 "zhbgst.f"
		i__5 = j - m;
#line 705 "zhbgst.f"
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
#line 705 "zhbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 705 "zhbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 706 "zhbgst.f"
		i__2 = ka1 + (j - *ka + 1) * ab_dim1;
#line 706 "zhbgst.f"
		i__5 = j - m;
#line 706 "zhbgst.f"
		i__6 = ka1 + (j - *ka + 1) * ab_dim1;
#line 706 "zhbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 706 "zhbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 707 "zhbgst.f"
/* L320: */
#line 707 "zhbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 712 "zhbgst.f"
	    if (nrt > 0) {
#line 712 "zhbgst.f"
		zlargv_(&nrt, &ab[ka1 + (j2t - *ka) * ab_dim1], &inca, &work[
			j2t - m], &ka1, &rwork[j2t - m], &ka1);
#line 712 "zhbgst.f"
	    }
#line 715 "zhbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 719 "zhbgst.f"
		i__1 = *ka - 1;
#line 719 "zhbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 720 "zhbgst.f"
		    zlartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &rwork[j2 - m]
			    , &work[j2 - m], &ka1);
#line 723 "zhbgst.f"
/* L330: */
#line 723 "zhbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 728 "zhbgst.f"
		zlar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &rwork[j2 - m], &
			work[j2 - m], &ka1);

#line 731 "zhbgst.f"
		zlacgv_(&nr, &work[j2 - m], &ka1);
#line 732 "zhbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 736 "zhbgst.f"
	    i__1 = *kb - k + 1;
#line 736 "zhbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 737 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 738 "zhbgst.f"
		if (nrt > 0) {
#line 738 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - 
			    m], &work[j2 - m], &ka1);
#line 738 "zhbgst.f"
		}
#line 742 "zhbgst.f"
/* L340: */
#line 742 "zhbgst.f"
	    }

#line 744 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 748 "zhbgst.f"
		i__1 = j1;
#line 748 "zhbgst.f"
		i__3 = ka1;
#line 748 "zhbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 749 "zhbgst.f"
		    i__2 = *n - m;
#line 749 "zhbgst.f"
		    zrot_(&i__2, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j - m], &work[j - m]
			    );
#line 751 "zhbgst.f"
/* L350: */
#line 751 "zhbgst.f"
		}
#line 752 "zhbgst.f"
	    }
#line 753 "zhbgst.f"
/* L360: */
#line 753 "zhbgst.f"
	}

#line 755 "zhbgst.f"
	if (update) {
#line 756 "zhbgst.f"
	    if (i2 <= *n && kbt > 0) {

/*              create nonzero element a(i-kbt+ka+1,i-kbt) outside the */
/*              band and store it in WORK(i-kbt) */

#line 761 "zhbgst.f"
		i__4 = i__ - kbt;
#line 761 "zhbgst.f"
		i__3 = kbt + 1 + (i__ - kbt) * bb_dim1;
#line 761 "zhbgst.f"
		z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 761 "zhbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 761 "zhbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 762 "zhbgst.f"
	    }
#line 763 "zhbgst.f"
	}

#line 765 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 766 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 767 "zhbgst.f"
		i__4 = 2, i__3 = k - i0 + 1;
#line 767 "zhbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 768 "zhbgst.f"
	    } else {
/* Computing MAX */
#line 769 "zhbgst.f"
		i__4 = 1, i__3 = k - i0 + 1;
#line 769 "zhbgst.f"
		j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 770 "zhbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 774 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 775 "zhbgst.f"
		nrt = (*n - j2 + *ka + l) / ka1;
#line 776 "zhbgst.f"
		if (nrt > 0) {
#line 776 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + (j2 - *ka) * ab_dim1], &
			    inca, &ab[ka1 - l + (j2 - *ka + 1) * ab_dim1], &
			    inca, &rwork[j2 - *ka], &work[j2 - *ka], &ka1);
#line 776 "zhbgst.f"
		}
#line 780 "zhbgst.f"
/* L370: */
#line 780 "zhbgst.f"
	    }
#line 781 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 782 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 783 "zhbgst.f"
	    i__4 = j2;
#line 783 "zhbgst.f"
	    i__3 = -ka1;
#line 783 "zhbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 784 "zhbgst.f"
		i__1 = j;
#line 784 "zhbgst.f"
		i__2 = j - *ka;
#line 784 "zhbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 785 "zhbgst.f"
		rwork[j] = rwork[j - *ka];
#line 786 "zhbgst.f"
/* L380: */
#line 786 "zhbgst.f"
	    }
#line 787 "zhbgst.f"
	    i__3 = j1;
#line 787 "zhbgst.f"
	    i__4 = ka1;
#line 787 "zhbgst.f"
	    for (j = j2; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+1,j-ka) outside the band */
/*              and store it in WORK(j) */

#line 792 "zhbgst.f"
		i__1 = j;
#line 792 "zhbgst.f"
		i__2 = j;
#line 792 "zhbgst.f"
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
#line 792 "zhbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 792 "zhbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 793 "zhbgst.f"
		i__1 = ka1 + (j - *ka + 1) * ab_dim1;
#line 793 "zhbgst.f"
		i__2 = j;
#line 793 "zhbgst.f"
		i__5 = ka1 + (j - *ka + 1) * ab_dim1;
#line 793 "zhbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 793 "zhbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 794 "zhbgst.f"
/* L390: */
#line 794 "zhbgst.f"
	    }
#line 795 "zhbgst.f"
	    if (update) {
#line 796 "zhbgst.f"
		if (i__ - k < *n - *ka && k <= kbt) {
#line 796 "zhbgst.f"
		    i__4 = i__ - k + *ka;
#line 796 "zhbgst.f"
		    i__3 = i__ - k;
#line 796 "zhbgst.f"
		    work[i__4].r = work[i__3].r, work[i__4].i = work[i__3].i;
#line 796 "zhbgst.f"
		}
#line 798 "zhbgst.f"
	    }
#line 799 "zhbgst.f"
/* L400: */
#line 799 "zhbgst.f"
	}

#line 801 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 802 "zhbgst.f"
	    i__4 = 1, i__3 = k - i0 + 1;
#line 802 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__3) * ka1;
#line 803 "zhbgst.f"
	    nr = (*n - j2 + *ka) / ka1;
#line 804 "zhbgst.f"
	    j1 = j2 + (nr - 1) * ka1;
#line 805 "zhbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 810 "zhbgst.f"
		zlargv_(&nr, &ab[ka1 + (j2 - *ka) * ab_dim1], &inca, &work[j2]
			, &ka1, &rwork[j2], &ka1);

/*              apply rotations in 2nd set from the left */

#line 815 "zhbgst.f"
		i__4 = *ka - 1;
#line 815 "zhbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 816 "zhbgst.f"
		    zlartv_(&nr, &ab[l + 1 + (j2 - l) * ab_dim1], &inca, &ab[
			    l + 2 + (j2 - l) * ab_dim1], &inca, &rwork[j2], &
			    work[j2], &ka1);
#line 819 "zhbgst.f"
/* L410: */
#line 819 "zhbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 824 "zhbgst.f"
		zlar2v_(&nr, &ab[j2 * ab_dim1 + 1], &ab[(j2 + 1) * ab_dim1 + 
			1], &ab[j2 * ab_dim1 + 2], &inca, &rwork[j2], &work[
			j2], &ka1);

#line 827 "zhbgst.f"
		zlacgv_(&nr, &work[j2], &ka1);
#line 828 "zhbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 832 "zhbgst.f"
	    i__4 = *kb - k + 1;
#line 832 "zhbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 833 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 834 "zhbgst.f"
		if (nrt > 0) {
#line 834 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2], 
			    &work[j2], &ka1);
#line 834 "zhbgst.f"
		}
#line 838 "zhbgst.f"
/* L420: */
#line 838 "zhbgst.f"
	    }

#line 840 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 844 "zhbgst.f"
		i__4 = j1;
#line 844 "zhbgst.f"
		i__3 = ka1;
#line 844 "zhbgst.f"
		for (j = j2; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 845 "zhbgst.f"
		    i__1 = *n - m;
#line 845 "zhbgst.f"
		    zrot_(&i__1, &x[m + 1 + j * x_dim1], &c__1, &x[m + 1 + (j 
			    + 1) * x_dim1], &c__1, &rwork[j], &work[j]);
#line 847 "zhbgst.f"
/* L430: */
#line 847 "zhbgst.f"
		}
#line 848 "zhbgst.f"
	    }
#line 849 "zhbgst.f"
/* L440: */
#line 849 "zhbgst.f"
	}

#line 851 "zhbgst.f"
	i__3 = *kb - 1;
#line 851 "zhbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 852 "zhbgst.f"
	    i__4 = 1, i__1 = k - i0 + 2;
#line 852 "zhbgst.f"
	    j2 = i__ - k - 1 + max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 856 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 857 "zhbgst.f"
		nrt = (*n - j2 + l) / ka1;
#line 858 "zhbgst.f"
		if (nrt > 0) {
#line 858 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + j2 * ab_dim1], &inca, &ab[
			    ka1 - l + (j2 + 1) * ab_dim1], &inca, &rwork[j2 - 
			    m], &work[j2 - m], &ka1);
#line 858 "zhbgst.f"
		}
#line 862 "zhbgst.f"
/* L450: */
#line 862 "zhbgst.f"
	    }
#line 863 "zhbgst.f"
/* L460: */
#line 863 "zhbgst.f"
	}

#line 865 "zhbgst.f"
	if (*kb > 1) {
#line 866 "zhbgst.f"
	    i__3 = j2 + *ka;
#line 866 "zhbgst.f"
	    for (j = *n - 1; j >= i__3; --j) {
#line 867 "zhbgst.f"
		rwork[j - m] = rwork[j - *ka - m];
#line 868 "zhbgst.f"
		i__4 = j - m;
#line 868 "zhbgst.f"
		i__1 = j - *ka - m;
#line 868 "zhbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 869 "zhbgst.f"
/* L470: */
#line 869 "zhbgst.f"
	    }
#line 870 "zhbgst.f"
	}

#line 872 "zhbgst.f"
    }

#line 874 "zhbgst.f"
    goto L10;

#line 876 "zhbgst.f"
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

#line 894 "zhbgst.f"
    update = TRUE_;
#line 895 "zhbgst.f"
    i__ = 0;
#line 896 "zhbgst.f"
L490:
#line 897 "zhbgst.f"
    if (update) {
#line 898 "zhbgst.f"
	++i__;
/* Computing MIN */
#line 899 "zhbgst.f"
	i__3 = *kb, i__4 = m - i__;
#line 899 "zhbgst.f"
	kbt = min(i__3,i__4);
#line 900 "zhbgst.f"
	i0 = i__ + 1;
/* Computing MAX */
#line 901 "zhbgst.f"
	i__3 = 1, i__4 = i__ - *ka;
#line 901 "zhbgst.f"
	i1 = max(i__3,i__4);
#line 902 "zhbgst.f"
	i2 = i__ + kbt - ka1;
#line 903 "zhbgst.f"
	if (i__ > m) {
#line 904 "zhbgst.f"
	    update = FALSE_;
#line 905 "zhbgst.f"
	    --i__;
#line 906 "zhbgst.f"
	    i0 = m + 1;
#line 907 "zhbgst.f"
	    if (*ka == 0) {
#line 907 "zhbgst.f"
		return 0;
#line 907 "zhbgst.f"
	    }
#line 909 "zhbgst.f"
	    goto L490;
#line 910 "zhbgst.f"
	}
#line 911 "zhbgst.f"
    } else {
#line 912 "zhbgst.f"
	i__ -= *ka;
#line 913 "zhbgst.f"
	if (i__ < 2) {
#line 913 "zhbgst.f"
	    return 0;
#line 913 "zhbgst.f"
	}
#line 915 "zhbgst.f"
    }

#line 917 "zhbgst.f"
    if (i__ < m - kbt) {
#line 918 "zhbgst.f"
	nx = m;
#line 919 "zhbgst.f"
    } else {
#line 920 "zhbgst.f"
	nx = *n;
#line 921 "zhbgst.f"
    }

#line 923 "zhbgst.f"
    if (upper) {

/*        Transform A, working with the upper triangle */

#line 927 "zhbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 931 "zhbgst.f"
	    i__3 = kb1 + i__ * bb_dim1;
#line 931 "zhbgst.f"
	    bii = bb[i__3].r;
#line 932 "zhbgst.f"
	    i__3 = ka1 + i__ * ab_dim1;
#line 932 "zhbgst.f"
	    i__4 = ka1 + i__ * ab_dim1;
#line 932 "zhbgst.f"
	    d__1 = ab[i__4].r / bii / bii;
#line 932 "zhbgst.f"
	    ab[i__3].r = d__1, ab[i__3].i = 0.;
#line 933 "zhbgst.f"
	    i__3 = i__ - 1;
#line 933 "zhbgst.f"
	    for (j = i1; j <= i__3; ++j) {
#line 934 "zhbgst.f"
		i__4 = j - i__ + ka1 + i__ * ab_dim1;
#line 934 "zhbgst.f"
		i__1 = j - i__ + ka1 + i__ * ab_dim1;
#line 934 "zhbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 934 "zhbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 935 "zhbgst.f"
/* L500: */
#line 935 "zhbgst.f"
	    }
/* Computing MIN */
#line 936 "zhbgst.f"
	    i__4 = *n, i__1 = i__ + *ka;
#line 936 "zhbgst.f"
	    i__3 = min(i__4,i__1);
#line 936 "zhbgst.f"
	    for (j = i__ + 1; j <= i__3; ++j) {
#line 937 "zhbgst.f"
		i__4 = i__ - j + ka1 + j * ab_dim1;
#line 937 "zhbgst.f"
		i__1 = i__ - j + ka1 + j * ab_dim1;
#line 937 "zhbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 937 "zhbgst.f"
		ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 938 "zhbgst.f"
/* L510: */
#line 938 "zhbgst.f"
	    }
#line 939 "zhbgst.f"
	    i__3 = i__ + kbt;
#line 939 "zhbgst.f"
	    for (k = i__ + 1; k <= i__3; ++k) {
#line 940 "zhbgst.f"
		i__4 = i__ + kbt;
#line 940 "zhbgst.f"
		for (j = k; j <= i__4; ++j) {
#line 941 "zhbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 941 "zhbgst.f"
		    i__2 = k - j + ka1 + j * ab_dim1;
#line 941 "zhbgst.f"
		    i__5 = i__ - j + kb1 + j * bb_dim1;
#line 941 "zhbgst.f"
		    d_cnjg(&z__5, &ab[i__ - k + ka1 + k * ab_dim1]);
#line 941 "zhbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 941 "zhbgst.f"
		    z__3.r = ab[i__2].r - z__4.r, z__3.i = ab[i__2].i - 
			    z__4.i;
#line 941 "zhbgst.f"
		    d_cnjg(&z__7, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 941 "zhbgst.f"
		    i__6 = i__ - j + ka1 + j * ab_dim1;
#line 941 "zhbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 941 "zhbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 941 "zhbgst.f"
		    i__7 = ka1 + i__ * ab_dim1;
#line 941 "zhbgst.f"
		    d__1 = ab[i__7].r;
#line 941 "zhbgst.f"
		    i__8 = i__ - j + kb1 + j * bb_dim1;
#line 941 "zhbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 941 "zhbgst.f"
		    d_cnjg(&z__10, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 941 "zhbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 941 "zhbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 941 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 949 "zhbgst.f"
/* L520: */
#line 949 "zhbgst.f"
		}
/* Computing MIN */
#line 950 "zhbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 950 "zhbgst.f"
		i__4 = min(i__1,i__2);
#line 950 "zhbgst.f"
		for (j = i__ + kbt + 1; j <= i__4; ++j) {
#line 951 "zhbgst.f"
		    i__1 = k - j + ka1 + j * ab_dim1;
#line 951 "zhbgst.f"
		    i__2 = k - j + ka1 + j * ab_dim1;
#line 951 "zhbgst.f"
		    d_cnjg(&z__3, &bb[i__ - k + kb1 + k * bb_dim1]);
#line 951 "zhbgst.f"
		    i__5 = i__ - j + ka1 + j * ab_dim1;
#line 951 "zhbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 951 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 951 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 954 "zhbgst.f"
/* L530: */
#line 954 "zhbgst.f"
		}
#line 955 "zhbgst.f"
/* L540: */
#line 955 "zhbgst.f"
	    }
#line 956 "zhbgst.f"
	    i__3 = i__;
#line 956 "zhbgst.f"
	    for (j = i1; j <= i__3; ++j) {
/* Computing MIN */
#line 957 "zhbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 957 "zhbgst.f"
		i__4 = min(i__1,i__2);
#line 957 "zhbgst.f"
		for (k = i__ + 1; k <= i__4; ++k) {
#line 958 "zhbgst.f"
		    i__1 = j - k + ka1 + k * ab_dim1;
#line 958 "zhbgst.f"
		    i__2 = j - k + ka1 + k * ab_dim1;
#line 958 "zhbgst.f"
		    i__5 = i__ - k + kb1 + k * bb_dim1;
#line 958 "zhbgst.f"
		    i__6 = j - i__ + ka1 + i__ * ab_dim1;
#line 958 "zhbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 958 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 958 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 960 "zhbgst.f"
/* L550: */
#line 960 "zhbgst.f"
		}
#line 961 "zhbgst.f"
/* L560: */
#line 961 "zhbgst.f"
	    }

#line 963 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 967 "zhbgst.f"
		d__1 = 1. / bii;
#line 967 "zhbgst.f"
		zdscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 968 "zhbgst.f"
		if (kbt > 0) {
#line 968 "zhbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 968 "zhbgst.f"
		    i__3 = *ldbb - 1;
#line 968 "zhbgst.f"
		    zgeru_(&nx, &kbt, &z__1, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    *kb + (i__ + 1) * bb_dim1], &i__3, &x[(i__ + 1) * 
			    x_dim1 + 1], ldx);
#line 968 "zhbgst.f"
		}
#line 971 "zhbgst.f"
	    }

/*           store a(i1,i) in RA1 for use in next loop over K */

#line 975 "zhbgst.f"
	    i__3 = i1 - i__ + ka1 + i__ * ab_dim1;
#line 975 "zhbgst.f"
	    ra1.r = ab[i__3].r, ra1.i = ab[i__3].i;
#line 976 "zhbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 981 "zhbgst.f"
	i__3 = *kb - 1;
#line 981 "zhbgst.f"
	for (k = 1; k <= i__3; ++k) {
#line 982 "zhbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 987 "zhbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i+k-ka-1,i) */

#line 991 "zhbgst.f"
		    zlartg_(&ab[k + 1 + i__ * ab_dim1], &ra1, &rwork[i__ + k 
			    - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k-ka-1,i+k) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 997 "zhbgst.f"
		    i__4 = kb1 - k + (i__ + k) * bb_dim1;
#line 997 "zhbgst.f"
		    z__2.r = -bb[i__4].r, z__2.i = -bb[i__4].i;
#line 997 "zhbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 997 "zhbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 998 "zhbgst.f"
		    i__4 = m - *kb + i__ + k;
#line 998 "zhbgst.f"
		    i__1 = i__ + k - *ka;
#line 998 "zhbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 998 "zhbgst.f"
		    d_cnjg(&z__4, &work[i__ + k - *ka]);
#line 998 "zhbgst.f"
		    i__2 = (i__ + k) * ab_dim1 + 1;
#line 998 "zhbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 998 "zhbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 998 "zhbgst.f"
		    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 1001 "zhbgst.f"
		    i__4 = (i__ + k) * ab_dim1 + 1;
#line 1001 "zhbgst.f"
		    i__1 = i__ + k - *ka;
#line 1001 "zhbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 1001 "zhbgst.f"
		    i__2 = i__ + k - *ka;
#line 1001 "zhbgst.f"
		    i__5 = (i__ + k) * ab_dim1 + 1;
#line 1001 "zhbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 1001 "zhbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 1001 "zhbgst.f"
		    ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 1003 "zhbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 1004 "zhbgst.f"
		}
#line 1005 "zhbgst.f"
	    }
/* Computing MAX */
#line 1006 "zhbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1006 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;
#line 1007 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1008 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1009 "zhbgst.f"
	    if (update) {
/* Computing MIN */
#line 1010 "zhbgst.f"
		i__4 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1010 "zhbgst.f"
		j2t = min(i__4,i__1);
#line 1011 "zhbgst.f"
	    } else {
#line 1012 "zhbgst.f"
		j2t = j2;
#line 1013 "zhbgst.f"
	    }
#line 1014 "zhbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1015 "zhbgst.f"
	    i__4 = j2t;
#line 1015 "zhbgst.f"
	    i__1 = ka1;
#line 1015 "zhbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__4 : j <= i__4; j += i__1) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(j) */

#line 1020 "zhbgst.f"
		i__2 = j;
#line 1020 "zhbgst.f"
		i__5 = j;
#line 1020 "zhbgst.f"
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
#line 1020 "zhbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 1020 "zhbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 1021 "zhbgst.f"
		i__2 = (j + *ka - 1) * ab_dim1 + 1;
#line 1021 "zhbgst.f"
		i__5 = j;
#line 1021 "zhbgst.f"
		i__6 = (j + *ka - 1) * ab_dim1 + 1;
#line 1021 "zhbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 1021 "zhbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 1022 "zhbgst.f"
/* L570: */
#line 1022 "zhbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1027 "zhbgst.f"
	    if (nrt > 0) {
#line 1027 "zhbgst.f"
		zlargv_(&nrt, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[j1],
			 &ka1, &rwork[j1], &ka1);
#line 1027 "zhbgst.f"
	    }
#line 1030 "zhbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the left */

#line 1034 "zhbgst.f"
		i__1 = *ka - 1;
#line 1034 "zhbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1035 "zhbgst.f"
		    zlartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &rwork[
			    j1], &work[j1], &ka1);
#line 1038 "zhbgst.f"
/* L580: */
#line 1038 "zhbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1043 "zhbgst.f"
		zlar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &rwork[j1], 
			&work[j1], &ka1);

#line 1047 "zhbgst.f"
		zlacgv_(&nr, &work[j1], &ka1);
#line 1048 "zhbgst.f"
	    }

/*           start applying rotations in 1st set from the right */

#line 1052 "zhbgst.f"
	    i__1 = *kb - k + 1;
#line 1052 "zhbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1053 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1054 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1055 "zhbgst.f"
		if (nrt > 0) {
#line 1055 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[j1t], &work[
			    j1t], &ka1);
#line 1055 "zhbgst.f"
		}
#line 1059 "zhbgst.f"
/* L590: */
#line 1059 "zhbgst.f"
	    }

#line 1061 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1065 "zhbgst.f"
		i__1 = j2;
#line 1065 "zhbgst.f"
		i__4 = ka1;
#line 1065 "zhbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__1 : j <= i__1; j += i__4) {
#line 1066 "zhbgst.f"
		    zrot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[j], &work[j]);
#line 1068 "zhbgst.f"
/* L600: */
#line 1068 "zhbgst.f"
		}
#line 1069 "zhbgst.f"
	    }
#line 1070 "zhbgst.f"
/* L610: */
#line 1070 "zhbgst.f"
	}

#line 1072 "zhbgst.f"
	if (update) {
#line 1073 "zhbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt-ka-1,i+kbt) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1078 "zhbgst.f"
		i__3 = m - *kb + i__ + kbt;
#line 1078 "zhbgst.f"
		i__4 = kb1 - kbt + (i__ + kbt) * bb_dim1;
#line 1078 "zhbgst.f"
		z__2.r = -bb[i__4].r, z__2.i = -bb[i__4].i;
#line 1078 "zhbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 1078 "zhbgst.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 1079 "zhbgst.f"
	    }
#line 1080 "zhbgst.f"
	}

#line 1082 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1083 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 1084 "zhbgst.f"
		i__3 = 2, i__4 = k + i0 - m;
#line 1084 "zhbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1085 "zhbgst.f"
	    } else {
/* Computing MAX */
#line 1086 "zhbgst.f"
		i__3 = 1, i__4 = k + i0 - m;
#line 1086 "zhbgst.f"
		j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1087 "zhbgst.f"
	    }

/*           finish applying rotations in 2nd set from the right */

#line 1091 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1092 "zhbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1093 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1094 "zhbgst.f"
		if (nrt > 0) {
#line 1094 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + (j1t + *ka) * ab_dim1], &inca, &ab[
			    l + 1 + (j1t + *ka - 1) * ab_dim1], &inca, &rwork[
			    m - *kb + j1t + *ka], &work[m - *kb + j1t + *ka], 
			    &ka1);
#line 1094 "zhbgst.f"
		}
#line 1099 "zhbgst.f"
/* L620: */
#line 1099 "zhbgst.f"
	    }
#line 1100 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1101 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1102 "zhbgst.f"
	    i__3 = j2;
#line 1102 "zhbgst.f"
	    i__4 = ka1;
#line 1102 "zhbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1103 "zhbgst.f"
		i__1 = m - *kb + j;
#line 1103 "zhbgst.f"
		i__2 = m - *kb + j + *ka;
#line 1103 "zhbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 1104 "zhbgst.f"
		rwork[m - *kb + j] = rwork[m - *kb + j + *ka];
#line 1105 "zhbgst.f"
/* L630: */
#line 1105 "zhbgst.f"
	    }
#line 1106 "zhbgst.f"
	    i__4 = j2;
#line 1106 "zhbgst.f"
	    i__3 = ka1;
#line 1106 "zhbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {

/*              create nonzero element a(j-1,j+ka) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1111 "zhbgst.f"
		i__1 = m - *kb + j;
#line 1111 "zhbgst.f"
		i__2 = m - *kb + j;
#line 1111 "zhbgst.f"
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
#line 1111 "zhbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 1111 "zhbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 1112 "zhbgst.f"
		i__1 = (j + *ka - 1) * ab_dim1 + 1;
#line 1112 "zhbgst.f"
		i__2 = m - *kb + j;
#line 1112 "zhbgst.f"
		i__5 = (j + *ka - 1) * ab_dim1 + 1;
#line 1112 "zhbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 1112 "zhbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1113 "zhbgst.f"
/* L640: */
#line 1113 "zhbgst.f"
	    }
#line 1114 "zhbgst.f"
	    if (update) {
#line 1115 "zhbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1115 "zhbgst.f"
		    i__3 = m - *kb + i__ + k - *ka;
#line 1115 "zhbgst.f"
		    i__4 = m - *kb + i__ + k;
#line 1115 "zhbgst.f"
		    work[i__3].r = work[i__4].r, work[i__3].i = work[i__4].i;
#line 1115 "zhbgst.f"
		}
#line 1117 "zhbgst.f"
	    }
#line 1118 "zhbgst.f"
/* L650: */
#line 1118 "zhbgst.f"
	}

#line 1120 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1121 "zhbgst.f"
	    i__3 = 1, i__4 = k + i0 - m;
#line 1121 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__4) * ka1;
#line 1122 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1123 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1124 "zhbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1129 "zhbgst.f"
		zlargv_(&nr, &ab[(j1 + *ka) * ab_dim1 + 1], &inca, &work[m - *
			kb + j1], &ka1, &rwork[m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the left */

#line 1134 "zhbgst.f"
		i__3 = *ka - 1;
#line 1134 "zhbgst.f"
		for (l = 1; l <= i__3; ++l) {
#line 1135 "zhbgst.f"
		    zlartv_(&nr, &ab[ka1 - l + (j1 + l) * ab_dim1], &inca, &
			    ab[*ka - l + (j1 + l) * ab_dim1], &inca, &rwork[m 
			    - *kb + j1], &work[m - *kb + j1], &ka1);
#line 1138 "zhbgst.f"
/* L660: */
#line 1138 "zhbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1143 "zhbgst.f"
		zlar2v_(&nr, &ab[ka1 + j1 * ab_dim1], &ab[ka1 + (j1 - 1) * 
			ab_dim1], &ab[*ka + j1 * ab_dim1], &inca, &rwork[m - *
			kb + j1], &work[m - *kb + j1], &ka1);

#line 1147 "zhbgst.f"
		zlacgv_(&nr, &work[m - *kb + j1], &ka1);
#line 1148 "zhbgst.f"
	    }

/*           start applying rotations in 2nd set from the right */

#line 1152 "zhbgst.f"
	    i__3 = *kb - k + 1;
#line 1152 "zhbgst.f"
	    for (l = *ka - 1; l >= i__3; --l) {
#line 1153 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1154 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1155 "zhbgst.f"
		if (nrt > 0) {
#line 1155 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[m - *kb + j1t],
			     &work[m - *kb + j1t], &ka1);
#line 1155 "zhbgst.f"
		}
#line 1160 "zhbgst.f"
/* L670: */
#line 1160 "zhbgst.f"
	    }

#line 1162 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1166 "zhbgst.f"
		i__3 = j2;
#line 1166 "zhbgst.f"
		i__4 = ka1;
#line 1166 "zhbgst.f"
		for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
#line 1167 "zhbgst.f"
		    zrot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[m - *kb + j], &work[m - *kb + 
			    j]);
#line 1169 "zhbgst.f"
/* L680: */
#line 1169 "zhbgst.f"
		}
#line 1170 "zhbgst.f"
	    }
#line 1171 "zhbgst.f"
/* L690: */
#line 1171 "zhbgst.f"
	}

#line 1173 "zhbgst.f"
	i__4 = *kb - 1;
#line 1173 "zhbgst.f"
	for (k = 1; k <= i__4; ++k) {
/* Computing MAX */
#line 1174 "zhbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1174 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;

/*           finish applying rotations in 1st set from the right */

#line 1178 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1179 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1180 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1181 "zhbgst.f"
		if (nrt > 0) {
#line 1181 "zhbgst.f"
		    zlartv_(&nrt, &ab[l + j1t * ab_dim1], &inca, &ab[l + 1 + (
			    j1t - 1) * ab_dim1], &inca, &rwork[j1t], &work[
			    j1t], &ka1);
#line 1181 "zhbgst.f"
		}
#line 1185 "zhbgst.f"
/* L700: */
#line 1185 "zhbgst.f"
	    }
#line 1186 "zhbgst.f"
/* L710: */
#line 1186 "zhbgst.f"
	}

#line 1188 "zhbgst.f"
	if (*kb > 1) {
#line 1189 "zhbgst.f"
	    i__4 = i2 - *ka;
#line 1189 "zhbgst.f"
	    for (j = 2; j <= i__4; ++j) {
#line 1190 "zhbgst.f"
		rwork[j] = rwork[j + *ka];
#line 1191 "zhbgst.f"
		i__3 = j;
#line 1191 "zhbgst.f"
		i__1 = j + *ka;
#line 1191 "zhbgst.f"
		work[i__3].r = work[i__1].r, work[i__3].i = work[i__1].i;
#line 1192 "zhbgst.f"
/* L720: */
#line 1192 "zhbgst.f"
	    }
#line 1193 "zhbgst.f"
	}

#line 1195 "zhbgst.f"
    } else {

/*        Transform A, working with the lower triangle */

#line 1199 "zhbgst.f"
	if (update) {

/*           Form  inv(S(i))**H * A * inv(S(i)) */

#line 1203 "zhbgst.f"
	    i__4 = i__ * bb_dim1 + 1;
#line 1203 "zhbgst.f"
	    bii = bb[i__4].r;
#line 1204 "zhbgst.f"
	    i__4 = i__ * ab_dim1 + 1;
#line 1204 "zhbgst.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 1204 "zhbgst.f"
	    d__1 = ab[i__3].r / bii / bii;
#line 1204 "zhbgst.f"
	    ab[i__4].r = d__1, ab[i__4].i = 0.;
#line 1205 "zhbgst.f"
	    i__4 = i__ - 1;
#line 1205 "zhbgst.f"
	    for (j = i1; j <= i__4; ++j) {
#line 1206 "zhbgst.f"
		i__3 = i__ - j + 1 + j * ab_dim1;
#line 1206 "zhbgst.f"
		i__1 = i__ - j + 1 + j * ab_dim1;
#line 1206 "zhbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 1206 "zhbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1207 "zhbgst.f"
/* L730: */
#line 1207 "zhbgst.f"
	    }
/* Computing MIN */
#line 1208 "zhbgst.f"
	    i__3 = *n, i__1 = i__ + *ka;
#line 1208 "zhbgst.f"
	    i__4 = min(i__3,i__1);
#line 1208 "zhbgst.f"
	    for (j = i__ + 1; j <= i__4; ++j) {
#line 1209 "zhbgst.f"
		i__3 = j - i__ + 1 + i__ * ab_dim1;
#line 1209 "zhbgst.f"
		i__1 = j - i__ + 1 + i__ * ab_dim1;
#line 1209 "zhbgst.f"
		z__1.r = ab[i__1].r / bii, z__1.i = ab[i__1].i / bii;
#line 1209 "zhbgst.f"
		ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1210 "zhbgst.f"
/* L740: */
#line 1210 "zhbgst.f"
	    }
#line 1211 "zhbgst.f"
	    i__4 = i__ + kbt;
#line 1211 "zhbgst.f"
	    for (k = i__ + 1; k <= i__4; ++k) {
#line 1212 "zhbgst.f"
		i__3 = i__ + kbt;
#line 1212 "zhbgst.f"
		for (j = k; j <= i__3; ++j) {
#line 1213 "zhbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 1213 "zhbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 1213 "zhbgst.f"
		    i__5 = j - i__ + 1 + i__ * bb_dim1;
#line 1213 "zhbgst.f"
		    d_cnjg(&z__5, &ab[k - i__ + 1 + i__ * ab_dim1]);
#line 1213 "zhbgst.f"
		    z__4.r = bb[i__5].r * z__5.r - bb[i__5].i * z__5.i, 
			    z__4.i = bb[i__5].r * z__5.i + bb[i__5].i * 
			    z__5.r;
#line 1213 "zhbgst.f"
		    z__3.r = ab[i__2].r - z__4.r, z__3.i = ab[i__2].i - 
			    z__4.i;
#line 1213 "zhbgst.f"
		    d_cnjg(&z__7, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1213 "zhbgst.f"
		    i__6 = j - i__ + 1 + i__ * ab_dim1;
#line 1213 "zhbgst.f"
		    z__6.r = z__7.r * ab[i__6].r - z__7.i * ab[i__6].i, 
			    z__6.i = z__7.r * ab[i__6].i + z__7.i * ab[i__6]
			    .r;
#line 1213 "zhbgst.f"
		    z__2.r = z__3.r - z__6.r, z__2.i = z__3.i - z__6.i;
#line 1213 "zhbgst.f"
		    i__7 = i__ * ab_dim1 + 1;
#line 1213 "zhbgst.f"
		    d__1 = ab[i__7].r;
#line 1213 "zhbgst.f"
		    i__8 = j - i__ + 1 + i__ * bb_dim1;
#line 1213 "zhbgst.f"
		    z__9.r = d__1 * bb[i__8].r, z__9.i = d__1 * bb[i__8].i;
#line 1213 "zhbgst.f"
		    d_cnjg(&z__10, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1213 "zhbgst.f"
		    z__8.r = z__9.r * z__10.r - z__9.i * z__10.i, z__8.i = 
			    z__9.r * z__10.i + z__9.i * z__10.r;
#line 1213 "zhbgst.f"
		    z__1.r = z__2.r + z__8.r, z__1.i = z__2.i + z__8.i;
#line 1213 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1219 "zhbgst.f"
/* L750: */
#line 1219 "zhbgst.f"
		}
/* Computing MIN */
#line 1220 "zhbgst.f"
		i__1 = *n, i__2 = i__ + *ka;
#line 1220 "zhbgst.f"
		i__3 = min(i__1,i__2);
#line 1220 "zhbgst.f"
		for (j = i__ + kbt + 1; j <= i__3; ++j) {
#line 1221 "zhbgst.f"
		    i__1 = j - k + 1 + k * ab_dim1;
#line 1221 "zhbgst.f"
		    i__2 = j - k + 1 + k * ab_dim1;
#line 1221 "zhbgst.f"
		    d_cnjg(&z__3, &bb[k - i__ + 1 + i__ * bb_dim1]);
#line 1221 "zhbgst.f"
		    i__5 = j - i__ + 1 + i__ * ab_dim1;
#line 1221 "zhbgst.f"
		    z__2.r = z__3.r * ab[i__5].r - z__3.i * ab[i__5].i, 
			    z__2.i = z__3.r * ab[i__5].i + z__3.i * ab[i__5]
			    .r;
#line 1221 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 1221 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1224 "zhbgst.f"
/* L760: */
#line 1224 "zhbgst.f"
		}
#line 1225 "zhbgst.f"
/* L770: */
#line 1225 "zhbgst.f"
	    }
#line 1226 "zhbgst.f"
	    i__4 = i__;
#line 1226 "zhbgst.f"
	    for (j = i1; j <= i__4; ++j) {
/* Computing MIN */
#line 1227 "zhbgst.f"
		i__1 = j + *ka, i__2 = i__ + kbt;
#line 1227 "zhbgst.f"
		i__3 = min(i__1,i__2);
#line 1227 "zhbgst.f"
		for (k = i__ + 1; k <= i__3; ++k) {
#line 1228 "zhbgst.f"
		    i__1 = k - j + 1 + j * ab_dim1;
#line 1228 "zhbgst.f"
		    i__2 = k - j + 1 + j * ab_dim1;
#line 1228 "zhbgst.f"
		    i__5 = k - i__ + 1 + i__ * bb_dim1;
#line 1228 "zhbgst.f"
		    i__6 = i__ - j + 1 + j * ab_dim1;
#line 1228 "zhbgst.f"
		    z__2.r = bb[i__5].r * ab[i__6].r - bb[i__5].i * ab[i__6]
			    .i, z__2.i = bb[i__5].r * ab[i__6].i + bb[i__5].i 
			    * ab[i__6].r;
#line 1228 "zhbgst.f"
		    z__1.r = ab[i__2].r - z__2.r, z__1.i = ab[i__2].i - 
			    z__2.i;
#line 1228 "zhbgst.f"
		    ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1230 "zhbgst.f"
/* L780: */
#line 1230 "zhbgst.f"
		}
#line 1231 "zhbgst.f"
/* L790: */
#line 1231 "zhbgst.f"
	    }

#line 1233 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by inv(S(i)) */

#line 1237 "zhbgst.f"
		d__1 = 1. / bii;
#line 1237 "zhbgst.f"
		zdscal_(&nx, &d__1, &x[i__ * x_dim1 + 1], &c__1);
#line 1238 "zhbgst.f"
		if (kbt > 0) {
#line 1238 "zhbgst.f"
		    z__1.r = -1., z__1.i = -0.;
#line 1238 "zhbgst.f"
		    zgerc_(&nx, &kbt, &z__1, &x[i__ * x_dim1 + 1], &c__1, &bb[
			    i__ * bb_dim1 + 2], &c__1, &x[(i__ + 1) * x_dim1 
			    + 1], ldx);
#line 1238 "zhbgst.f"
		}
#line 1241 "zhbgst.f"
	    }

/*           store a(i,i1) in RA1 for use in next loop over K */

#line 1245 "zhbgst.f"
	    i__4 = i__ - i1 + 1 + i1 * ab_dim1;
#line 1245 "zhbgst.f"
	    ra1.r = ab[i__4].r, ra1.i = ab[i__4].i;
#line 1246 "zhbgst.f"
	}

/*        Generate and apply vectors of rotations to chase all the */
/*        existing bulges KA positions up toward the top of the band */

#line 1251 "zhbgst.f"
	i__4 = *kb - 1;
#line 1251 "zhbgst.f"
	for (k = 1; k <= i__4; ++k) {
#line 1252 "zhbgst.f"
	    if (update) {

/*              Determine the rotations which would annihilate the bulge */
/*              which has in theory just been created */

#line 1257 "zhbgst.f"
		if (i__ + k - ka1 > 0 && i__ + k < m) {

/*                 generate rotation to annihilate a(i,i+k-ka-1) */

#line 1261 "zhbgst.f"
		    zlartg_(&ab[ka1 - k + (i__ + k - *ka) * ab_dim1], &ra1, &
			    rwork[i__ + k - *ka], &work[i__ + k - *ka], &ra);

/*                 create nonzero element a(i+k,i+k-ka-1) outside the */
/*                 band and store it in WORK(m-kb+i+k) */

#line 1267 "zhbgst.f"
		    i__3 = k + 1 + i__ * bb_dim1;
#line 1267 "zhbgst.f"
		    z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 1267 "zhbgst.f"
		    z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r 
			    * ra1.i + z__2.i * ra1.r;
#line 1267 "zhbgst.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 1268 "zhbgst.f"
		    i__3 = m - *kb + i__ + k;
#line 1268 "zhbgst.f"
		    i__1 = i__ + k - *ka;
#line 1268 "zhbgst.f"
		    z__2.r = rwork[i__1] * t.r, z__2.i = rwork[i__1] * t.i;
#line 1268 "zhbgst.f"
		    d_cnjg(&z__4, &work[i__ + k - *ka]);
#line 1268 "zhbgst.f"
		    i__2 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1268 "zhbgst.f"
		    z__3.r = z__4.r * ab[i__2].r - z__4.i * ab[i__2].i, 
			    z__3.i = z__4.r * ab[i__2].i + z__4.i * ab[i__2]
			    .r;
#line 1268 "zhbgst.f"
		    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 1268 "zhbgst.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 1271 "zhbgst.f"
		    i__3 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1271 "zhbgst.f"
		    i__1 = i__ + k - *ka;
#line 1271 "zhbgst.f"
		    z__2.r = work[i__1].r * t.r - work[i__1].i * t.i, z__2.i =
			     work[i__1].r * t.i + work[i__1].i * t.r;
#line 1271 "zhbgst.f"
		    i__2 = i__ + k - *ka;
#line 1271 "zhbgst.f"
		    i__5 = ka1 + (i__ + k - *ka) * ab_dim1;
#line 1271 "zhbgst.f"
		    z__3.r = rwork[i__2] * ab[i__5].r, z__3.i = rwork[i__2] * 
			    ab[i__5].i;
#line 1271 "zhbgst.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 1271 "zhbgst.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 1273 "zhbgst.f"
		    ra1.r = ra.r, ra1.i = ra.i;
#line 1274 "zhbgst.f"
		}
#line 1275 "zhbgst.f"
	    }
/* Computing MAX */
#line 1276 "zhbgst.f"
	    i__3 = 1, i__1 = k + i0 - m + 1;
#line 1276 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__3,i__1) * ka1;
#line 1277 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1278 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1279 "zhbgst.f"
	    if (update) {
/* Computing MIN */
#line 1280 "zhbgst.f"
		i__3 = j2, i__1 = i__ - (*ka << 1) + k - 1;
#line 1280 "zhbgst.f"
		j2t = min(i__3,i__1);
#line 1281 "zhbgst.f"
	    } else {
#line 1282 "zhbgst.f"
		j2t = j2;
#line 1283 "zhbgst.f"
	    }
#line 1284 "zhbgst.f"
	    nrt = (j2t + *ka - 1) / ka1;
#line 1285 "zhbgst.f"
	    i__3 = j2t;
#line 1285 "zhbgst.f"
	    i__1 = ka1;
#line 1285 "zhbgst.f"
	    for (j = j1; i__1 < 0 ? j >= i__3 : j <= i__3; j += i__1) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(j) */

#line 1290 "zhbgst.f"
		i__2 = j;
#line 1290 "zhbgst.f"
		i__5 = j;
#line 1290 "zhbgst.f"
		i__6 = ka1 + (j - 1) * ab_dim1;
#line 1290 "zhbgst.f"
		z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * ab[i__6]
			.i, z__1.i = work[i__5].r * ab[i__6].i + work[i__5].i 
			* ab[i__6].r;
#line 1290 "zhbgst.f"
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
#line 1291 "zhbgst.f"
		i__2 = ka1 + (j - 1) * ab_dim1;
#line 1291 "zhbgst.f"
		i__5 = j;
#line 1291 "zhbgst.f"
		i__6 = ka1 + (j - 1) * ab_dim1;
#line 1291 "zhbgst.f"
		z__1.r = rwork[i__5] * ab[i__6].r, z__1.i = rwork[i__5] * ab[
			i__6].i;
#line 1291 "zhbgst.f"
		ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 1292 "zhbgst.f"
/* L800: */
#line 1292 "zhbgst.f"
	    }

/*           generate rotations in 1st set to annihilate elements which */
/*           have been created outside the band */

#line 1297 "zhbgst.f"
	    if (nrt > 0) {
#line 1297 "zhbgst.f"
		zlargv_(&nrt, &ab[ka1 + j1 * ab_dim1], &inca, &work[j1], &ka1,
			 &rwork[j1], &ka1);
#line 1297 "zhbgst.f"
	    }
#line 1300 "zhbgst.f"
	    if (nr > 0) {

/*              apply rotations in 1st set from the right */

#line 1304 "zhbgst.f"
		i__1 = *ka - 1;
#line 1304 "zhbgst.f"
		for (l = 1; l <= i__1; ++l) {
#line 1305 "zhbgst.f"
		    zlartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &rwork[j1], &work[
			    j1], &ka1);
#line 1307 "zhbgst.f"
/* L810: */
#line 1307 "zhbgst.f"
		}

/*              apply rotations in 1st set from both sides to diagonal */
/*              blocks */

#line 1312 "zhbgst.f"
		zlar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &rwork[j1], &
			work[j1], &ka1);

#line 1316 "zhbgst.f"
		zlacgv_(&nr, &work[j1], &ka1);
#line 1317 "zhbgst.f"
	    }

/*           start applying rotations in 1st set from the left */

#line 1321 "zhbgst.f"
	    i__1 = *kb - k + 1;
#line 1321 "zhbgst.f"
	    for (l = *ka - 1; l >= i__1; --l) {
#line 1322 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1323 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1324 "zhbgst.f"
		if (nrt > 0) {
#line 1324 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[j1t], &work[j1t], &ka1);
#line 1324 "zhbgst.f"
		}
#line 1328 "zhbgst.f"
/* L820: */
#line 1328 "zhbgst.f"
	    }

#line 1330 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 1st set */

#line 1334 "zhbgst.f"
		i__1 = j2;
#line 1334 "zhbgst.f"
		i__3 = ka1;
#line 1334 "zhbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
#line 1335 "zhbgst.f"
		    d_cnjg(&z__1, &work[j]);
#line 1335 "zhbgst.f"
		    zrot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[j], &z__1);
#line 1337 "zhbgst.f"
/* L830: */
#line 1337 "zhbgst.f"
		}
#line 1338 "zhbgst.f"
	    }
#line 1339 "zhbgst.f"
/* L840: */
#line 1339 "zhbgst.f"
	}

#line 1341 "zhbgst.f"
	if (update) {
#line 1342 "zhbgst.f"
	    if (i2 > 0 && kbt > 0) {

/*              create nonzero element a(i+kbt,i+kbt-ka-1) outside the */
/*              band and store it in WORK(m-kb+i+kbt) */

#line 1347 "zhbgst.f"
		i__4 = m - *kb + i__ + kbt;
#line 1347 "zhbgst.f"
		i__3 = kbt + 1 + i__ * bb_dim1;
#line 1347 "zhbgst.f"
		z__2.r = -bb[i__3].r, z__2.i = -bb[i__3].i;
#line 1347 "zhbgst.f"
		z__1.r = z__2.r * ra1.r - z__2.i * ra1.i, z__1.i = z__2.r * 
			ra1.i + z__2.i * ra1.r;
#line 1347 "zhbgst.f"
		work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 1348 "zhbgst.f"
	    }
#line 1349 "zhbgst.f"
	}

#line 1351 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
#line 1352 "zhbgst.f"
	    if (update) {
/* Computing MAX */
#line 1353 "zhbgst.f"
		i__4 = 2, i__3 = k + i0 - m;
#line 1353 "zhbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1354 "zhbgst.f"
	    } else {
/* Computing MAX */
#line 1355 "zhbgst.f"
		i__4 = 1, i__3 = k + i0 - m;
#line 1355 "zhbgst.f"
		j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1356 "zhbgst.f"
	    }

/*           finish applying rotations in 2nd set from the left */

#line 1360 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1361 "zhbgst.f"
		nrt = (j2 + *ka + l - 1) / ka1;
#line 1362 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1363 "zhbgst.f"
		if (nrt > 0) {
#line 1363 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + (j1t + l - 1) * ab_dim1], 
			    &inca, &ab[ka1 - l + (j1t + l - 1) * ab_dim1], &
			    inca, &rwork[m - *kb + j1t + *ka], &work[m - *kb 
			    + j1t + *ka], &ka1);
#line 1363 "zhbgst.f"
		}
#line 1368 "zhbgst.f"
/* L850: */
#line 1368 "zhbgst.f"
	    }
#line 1369 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1370 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1371 "zhbgst.f"
	    i__4 = j2;
#line 1371 "zhbgst.f"
	    i__3 = ka1;
#line 1371 "zhbgst.f"
	    for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1372 "zhbgst.f"
		i__1 = m - *kb + j;
#line 1372 "zhbgst.f"
		i__2 = m - *kb + j + *ka;
#line 1372 "zhbgst.f"
		work[i__1].r = work[i__2].r, work[i__1].i = work[i__2].i;
#line 1373 "zhbgst.f"
		rwork[m - *kb + j] = rwork[m - *kb + j + *ka];
#line 1374 "zhbgst.f"
/* L860: */
#line 1374 "zhbgst.f"
	    }
#line 1375 "zhbgst.f"
	    i__3 = j2;
#line 1375 "zhbgst.f"
	    i__4 = ka1;
#line 1375 "zhbgst.f"
	    for (j = j1; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {

/*              create nonzero element a(j+ka,j-1) outside the band */
/*              and store it in WORK(m-kb+j) */

#line 1380 "zhbgst.f"
		i__1 = m - *kb + j;
#line 1380 "zhbgst.f"
		i__2 = m - *kb + j;
#line 1380 "zhbgst.f"
		i__5 = ka1 + (j - 1) * ab_dim1;
#line 1380 "zhbgst.f"
		z__1.r = work[i__2].r * ab[i__5].r - work[i__2].i * ab[i__5]
			.i, z__1.i = work[i__2].r * ab[i__5].i + work[i__2].i 
			* ab[i__5].r;
#line 1380 "zhbgst.f"
		work[i__1].r = z__1.r, work[i__1].i = z__1.i;
#line 1381 "zhbgst.f"
		i__1 = ka1 + (j - 1) * ab_dim1;
#line 1381 "zhbgst.f"
		i__2 = m - *kb + j;
#line 1381 "zhbgst.f"
		i__5 = ka1 + (j - 1) * ab_dim1;
#line 1381 "zhbgst.f"
		z__1.r = rwork[i__2] * ab[i__5].r, z__1.i = rwork[i__2] * ab[
			i__5].i;
#line 1381 "zhbgst.f"
		ab[i__1].r = z__1.r, ab[i__1].i = z__1.i;
#line 1382 "zhbgst.f"
/* L870: */
#line 1382 "zhbgst.f"
	    }
#line 1383 "zhbgst.f"
	    if (update) {
#line 1384 "zhbgst.f"
		if (i__ + k > ka1 && k <= kbt) {
#line 1384 "zhbgst.f"
		    i__4 = m - *kb + i__ + k - *ka;
#line 1384 "zhbgst.f"
		    i__3 = m - *kb + i__ + k;
#line 1384 "zhbgst.f"
		    work[i__4].r = work[i__3].r, work[i__4].i = work[i__3].i;
#line 1384 "zhbgst.f"
		}
#line 1386 "zhbgst.f"
	    }
#line 1387 "zhbgst.f"
/* L880: */
#line 1387 "zhbgst.f"
	}

#line 1389 "zhbgst.f"
	for (k = *kb; k >= 1; --k) {
/* Computing MAX */
#line 1390 "zhbgst.f"
	    i__4 = 1, i__3 = k + i0 - m;
#line 1390 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__3) * ka1;
#line 1391 "zhbgst.f"
	    nr = (j2 + *ka - 1) / ka1;
#line 1392 "zhbgst.f"
	    j1 = j2 - (nr - 1) * ka1;
#line 1393 "zhbgst.f"
	    if (nr > 0) {

/*              generate rotations in 2nd set to annihilate elements */
/*              which have been created outside the band */

#line 1398 "zhbgst.f"
		zlargv_(&nr, &ab[ka1 + j1 * ab_dim1], &inca, &work[m - *kb + 
			j1], &ka1, &rwork[m - *kb + j1], &ka1);

/*              apply rotations in 2nd set from the right */

#line 1403 "zhbgst.f"
		i__4 = *ka - 1;
#line 1403 "zhbgst.f"
		for (l = 1; l <= i__4; ++l) {
#line 1404 "zhbgst.f"
		    zlartv_(&nr, &ab[l + 1 + j1 * ab_dim1], &inca, &ab[l + 2 
			    + (j1 - 1) * ab_dim1], &inca, &rwork[m - *kb + j1]
			    , &work[m - *kb + j1], &ka1);
#line 1407 "zhbgst.f"
/* L890: */
#line 1407 "zhbgst.f"
		}

/*              apply rotations in 2nd set from both sides to diagonal */
/*              blocks */

#line 1412 "zhbgst.f"
		zlar2v_(&nr, &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 
			1], &ab[(j1 - 1) * ab_dim1 + 2], &inca, &rwork[m - *
			kb + j1], &work[m - *kb + j1], &ka1);

#line 1416 "zhbgst.f"
		zlacgv_(&nr, &work[m - *kb + j1], &ka1);
#line 1417 "zhbgst.f"
	    }

/*           start applying rotations in 2nd set from the left */

#line 1421 "zhbgst.f"
	    i__4 = *kb - k + 1;
#line 1421 "zhbgst.f"
	    for (l = *ka - 1; l >= i__4; --l) {
#line 1422 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1423 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1424 "zhbgst.f"
		if (nrt > 0) {
#line 1424 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[m - *kb + j1t], &work[m - *kb + 
			    j1t], &ka1);
#line 1424 "zhbgst.f"
		}
#line 1429 "zhbgst.f"
/* L900: */
#line 1429 "zhbgst.f"
	    }

#line 1431 "zhbgst.f"
	    if (wantx) {

/*              post-multiply X by product of rotations in 2nd set */

#line 1435 "zhbgst.f"
		i__4 = j2;
#line 1435 "zhbgst.f"
		i__3 = ka1;
#line 1435 "zhbgst.f"
		for (j = j1; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
#line 1436 "zhbgst.f"
		    d_cnjg(&z__1, &work[m - *kb + j]);
#line 1436 "zhbgst.f"
		    zrot_(&nx, &x[j * x_dim1 + 1], &c__1, &x[(j - 1) * x_dim1 
			    + 1], &c__1, &rwork[m - *kb + j], &z__1);
#line 1438 "zhbgst.f"
/* L910: */
#line 1438 "zhbgst.f"
		}
#line 1439 "zhbgst.f"
	    }
#line 1440 "zhbgst.f"
/* L920: */
#line 1440 "zhbgst.f"
	}

#line 1442 "zhbgst.f"
	i__3 = *kb - 1;
#line 1442 "zhbgst.f"
	for (k = 1; k <= i__3; ++k) {
/* Computing MAX */
#line 1443 "zhbgst.f"
	    i__4 = 1, i__1 = k + i0 - m + 1;
#line 1443 "zhbgst.f"
	    j2 = i__ + k + 1 - max(i__4,i__1) * ka1;

/*           finish applying rotations in 1st set from the left */

#line 1447 "zhbgst.f"
	    for (l = *kb - k; l >= 1; --l) {
#line 1448 "zhbgst.f"
		nrt = (j2 + l - 1) / ka1;
#line 1449 "zhbgst.f"
		j1t = j2 - (nrt - 1) * ka1;
#line 1450 "zhbgst.f"
		if (nrt > 0) {
#line 1450 "zhbgst.f"
		    zlartv_(&nrt, &ab[ka1 - l + 1 + (j1t - ka1 + l) * ab_dim1]
			    , &inca, &ab[ka1 - l + (j1t - ka1 + l) * ab_dim1],
			     &inca, &rwork[j1t], &work[j1t], &ka1);
#line 1450 "zhbgst.f"
		}
#line 1454 "zhbgst.f"
/* L930: */
#line 1454 "zhbgst.f"
	    }
#line 1455 "zhbgst.f"
/* L940: */
#line 1455 "zhbgst.f"
	}

#line 1457 "zhbgst.f"
	if (*kb > 1) {
#line 1458 "zhbgst.f"
	    i__3 = i2 - *ka;
#line 1458 "zhbgst.f"
	    for (j = 2; j <= i__3; ++j) {
#line 1459 "zhbgst.f"
		rwork[j] = rwork[j + *ka];
#line 1460 "zhbgst.f"
		i__4 = j;
#line 1460 "zhbgst.f"
		i__1 = j + *ka;
#line 1460 "zhbgst.f"
		work[i__4].r = work[i__1].r, work[i__4].i = work[i__1].i;
#line 1461 "zhbgst.f"
/* L950: */
#line 1461 "zhbgst.f"
	    }
#line 1462 "zhbgst.f"
	}

#line 1464 "zhbgst.f"
    }

#line 1466 "zhbgst.f"
    goto L490;

/*     End of ZHBGST */

} /* zhbgst_ */

