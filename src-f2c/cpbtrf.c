#line 1 "cpbtrf.f"
/* cpbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "cpbtrf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b21 = -1.;
static doublereal c_b22 = 1.;
static integer c__33 = 33;

/* > \brief \b CPBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBTRF( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBTRF computes the Cholesky factorization of a complex Hermitian */
/* > positive definite band matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**H * U,  if UPLO = 'U', or */
/* >    A = L  * L**H,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**H*U or A = L*L**H of the band */
/* >          matrix A, in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the factorization could not be */
/* >                completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  N = 6, KD = 2, and UPLO = 'U': */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46 */
/* >      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 */
/* >     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 */
/* > */
/* >  Similarly, if UPLO = 'L' the format of A is as follows: */
/* > */
/* >  On entry:                       On exit: */
/* > */
/* >     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66 */
/* >     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   * */
/* >     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    * */
/* > */
/* >  Array elements marked * are not used by the routine. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989 */

/*  ===================================================================== */
/* Subroutine */ int cpbtrf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, i2, i3, ib, nb, ii, jj;
    static doublecomplex work[1056]	/* was [33][32] */;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cherk_(char *, char *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, doublereal *,
	     doublecomplex *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), 
	    cpbtf2_(char *, integer *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen), cpotf2_(char *, integer *, doublecomplex *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 189 "cpbtrf.f"
    /* Parameter adjustments */
#line 189 "cpbtrf.f"
    ab_dim1 = *ldab;
#line 189 "cpbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 189 "cpbtrf.f"
    ab -= ab_offset;
#line 189 "cpbtrf.f"

#line 189 "cpbtrf.f"
    /* Function Body */
#line 189 "cpbtrf.f"
    *info = 0;
#line 190 "cpbtrf.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 192 "cpbtrf.f"
	*info = -1;
#line 193 "cpbtrf.f"
    } else if (*n < 0) {
#line 194 "cpbtrf.f"
	*info = -2;
#line 195 "cpbtrf.f"
    } else if (*kd < 0) {
#line 196 "cpbtrf.f"
	*info = -3;
#line 197 "cpbtrf.f"
    } else if (*ldab < *kd + 1) {
#line 198 "cpbtrf.f"
	*info = -5;
#line 199 "cpbtrf.f"
    }
#line 200 "cpbtrf.f"
    if (*info != 0) {
#line 201 "cpbtrf.f"
	i__1 = -(*info);
#line 201 "cpbtrf.f"
	xerbla_("CPBTRF", &i__1, (ftnlen)6);
#line 202 "cpbtrf.f"
	return 0;
#line 203 "cpbtrf.f"
    }

/*     Quick return if possible */

#line 207 "cpbtrf.f"
    if (*n == 0) {
#line 207 "cpbtrf.f"
	return 0;
#line 207 "cpbtrf.f"
    }

/*     Determine the block size for this environment */

#line 212 "cpbtrf.f"
    nb = ilaenv_(&c__1, "CPBTRF", uplo, n, kd, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     The block size must not exceed the semi-bandwidth KD, and must not */
/*     exceed the limit set by the size of the local array WORK. */

#line 217 "cpbtrf.f"
    nb = min(nb,32);

#line 219 "cpbtrf.f"
    if (nb <= 1 || nb > *kd) {

/*        Use unblocked code */

#line 223 "cpbtrf.f"
	cpbtf2_(uplo, n, kd, &ab[ab_offset], ldab, info, (ftnlen)1);
#line 224 "cpbtrf.f"
    } else {

/*        Use blocked code */

#line 228 "cpbtrf.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Compute the Cholesky factorization of a Hermitian band */
/*           matrix, given the upper triangle of the matrix in band */
/*           storage. */

/*           Zero the upper triangle of the work array. */

#line 236 "cpbtrf.f"
	    i__1 = nb;
#line 236 "cpbtrf.f"
	    for (j = 1; j <= i__1; ++j) {
#line 237 "cpbtrf.f"
		i__2 = j - 1;
#line 237 "cpbtrf.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 238 "cpbtrf.f"
		    i__3 = i__ + j * 33 - 34;
#line 238 "cpbtrf.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 239 "cpbtrf.f"
/* L10: */
#line 239 "cpbtrf.f"
		}
#line 240 "cpbtrf.f"
/* L20: */
#line 240 "cpbtrf.f"
	    }

/*           Process the band matrix one diagonal block at a time. */

#line 244 "cpbtrf.f"
	    i__1 = *n;
#line 244 "cpbtrf.f"
	    i__2 = nb;
#line 244 "cpbtrf.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 245 "cpbtrf.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 245 "cpbtrf.f"
		ib = min(i__3,i__4);

/*              Factorize the diagonal block */

#line 249 "cpbtrf.f"
		i__3 = *ldab - 1;
#line 249 "cpbtrf.f"
		cpotf2_(uplo, &ib, &ab[*kd + 1 + i__ * ab_dim1], &i__3, &ii, (
			ftnlen)1);
#line 250 "cpbtrf.f"
		if (ii != 0) {
#line 251 "cpbtrf.f"
		    *info = i__ + ii - 1;
#line 252 "cpbtrf.f"
		    goto L150;
#line 253 "cpbtrf.f"
		}
#line 254 "cpbtrf.f"
		if (i__ + ib <= *n) {

/*                 Update the relevant part of the trailing submatrix. */
/*                 If A11 denotes the diagonal block which has just been */
/*                 factorized, then we need to update the remaining */
/*                 blocks in the diagram: */

/*                    A11   A12   A13 */
/*                          A22   A23 */
/*                                A33 */

/*                 The numbers of rows and columns in the partitioning */
/*                 are IB, I2, I3 respectively. The blocks A12, A22 and */
/*                 A23 are empty if IB = KD. The upper triangle of A13 */
/*                 lies outside the band. */

/* Computing MIN */
#line 270 "cpbtrf.f"
		    i__3 = *kd - ib, i__4 = *n - i__ - ib + 1;
#line 270 "cpbtrf.f"
		    i2 = min(i__3,i__4);
/* Computing MIN */
#line 271 "cpbtrf.f"
		    i__3 = ib, i__4 = *n - i__ - *kd + 1;
#line 271 "cpbtrf.f"
		    i3 = min(i__3,i__4);

#line 273 "cpbtrf.f"
		    if (i2 > 0) {

/*                    Update A12 */

#line 277 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 277 "cpbtrf.f"
			i__4 = *ldab - 1;
#line 277 "cpbtrf.f"
			ctrsm_("Left", "Upper", "Conjugate transpose", "Non-"\
				"unit", &ib, &i2, &c_b1, &ab[*kd + 1 + i__ * 
				ab_dim1], &i__3, &ab[*kd + 1 - ib + (i__ + ib)
				 * ab_dim1], &i__4, (ftnlen)4, (ftnlen)5, (
				ftnlen)19, (ftnlen)8);

/*                    Update A22 */

#line 284 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 284 "cpbtrf.f"
			i__4 = *ldab - 1;
#line 284 "cpbtrf.f"
			cherk_("Upper", "Conjugate transpose", &i2, &ib, &
				c_b21, &ab[*kd + 1 - ib + (i__ + ib) * 
				ab_dim1], &i__3, &c_b22, &ab[*kd + 1 + (i__ + 
				ib) * ab_dim1], &i__4, (ftnlen)5, (ftnlen)19);
#line 287 "cpbtrf.f"
		    }

#line 289 "cpbtrf.f"
		    if (i3 > 0) {

/*                    Copy the lower triangle of A13 into the work array. */

#line 293 "cpbtrf.f"
			i__3 = i3;
#line 293 "cpbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 294 "cpbtrf.f"
			    i__4 = ib;
#line 294 "cpbtrf.f"
			    for (ii = jj; ii <= i__4; ++ii) {
#line 295 "cpbtrf.f"
				i__5 = ii + jj * 33 - 34;
#line 295 "cpbtrf.f"
				i__6 = ii - jj + 1 + (jj + i__ + *kd - 1) * 
					ab_dim1;
#line 295 "cpbtrf.f"
				work[i__5].r = ab[i__6].r, work[i__5].i = ab[
					i__6].i;
#line 296 "cpbtrf.f"
/* L30: */
#line 296 "cpbtrf.f"
			    }
#line 297 "cpbtrf.f"
/* L40: */
#line 297 "cpbtrf.f"
			}

/*                    Update A13 (in the work array). */

#line 301 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 301 "cpbtrf.f"
			ctrsm_("Left", "Upper", "Conjugate transpose", "Non-"\
				"unit", &ib, &i3, &c_b1, &ab[*kd + 1 + i__ * 
				ab_dim1], &i__3, work, &c__33, (ftnlen)4, (
				ftnlen)5, (ftnlen)19, (ftnlen)8);

/*                    Update A23 */

#line 307 "cpbtrf.f"
			if (i2 > 0) {
#line 307 "cpbtrf.f"
			    z__1.r = -1., z__1.i = -0.;
#line 307 "cpbtrf.f"
			    i__3 = *ldab - 1;
#line 307 "cpbtrf.f"
			    i__4 = *ldab - 1;
#line 307 "cpbtrf.f"
			    cgemm_("Conjugate transpose", "No transpose", &i2,
				     &i3, &ib, &z__1, &ab[*kd + 1 - ib + (i__ 
				    + ib) * ab_dim1], &i__3, work, &c__33, &
				    c_b1, &ab[ib + 1 + (i__ + *kd) * ab_dim1],
				     &i__4, (ftnlen)19, (ftnlen)12);
#line 307 "cpbtrf.f"
			}

/*                    Update A33 */

#line 316 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 316 "cpbtrf.f"
			cherk_("Upper", "Conjugate transpose", &i3, &ib, &
				c_b21, work, &c__33, &c_b22, &ab[*kd + 1 + (
				i__ + *kd) * ab_dim1], &i__3, (ftnlen)5, (
				ftnlen)19);

/*                    Copy the lower triangle of A13 back into place. */

#line 322 "cpbtrf.f"
			i__3 = i3;
#line 322 "cpbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 323 "cpbtrf.f"
			    i__4 = ib;
#line 323 "cpbtrf.f"
			    for (ii = jj; ii <= i__4; ++ii) {
#line 324 "cpbtrf.f"
				i__5 = ii - jj + 1 + (jj + i__ + *kd - 1) * 
					ab_dim1;
#line 324 "cpbtrf.f"
				i__6 = ii + jj * 33 - 34;
#line 324 "cpbtrf.f"
				ab[i__5].r = work[i__6].r, ab[i__5].i = work[
					i__6].i;
#line 325 "cpbtrf.f"
/* L50: */
#line 325 "cpbtrf.f"
			    }
#line 326 "cpbtrf.f"
/* L60: */
#line 326 "cpbtrf.f"
			}
#line 327 "cpbtrf.f"
		    }
#line 328 "cpbtrf.f"
		}
#line 329 "cpbtrf.f"
/* L70: */
#line 329 "cpbtrf.f"
	    }
#line 330 "cpbtrf.f"
	} else {

/*           Compute the Cholesky factorization of a Hermitian band */
/*           matrix, given the lower triangle of the matrix in band */
/*           storage. */

/*           Zero the lower triangle of the work array. */

#line 338 "cpbtrf.f"
	    i__2 = nb;
#line 338 "cpbtrf.f"
	    for (j = 1; j <= i__2; ++j) {
#line 339 "cpbtrf.f"
		i__1 = nb;
#line 339 "cpbtrf.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 340 "cpbtrf.f"
		    i__3 = i__ + j * 33 - 34;
#line 340 "cpbtrf.f"
		    work[i__3].r = 0., work[i__3].i = 0.;
#line 341 "cpbtrf.f"
/* L80: */
#line 341 "cpbtrf.f"
		}
#line 342 "cpbtrf.f"
/* L90: */
#line 342 "cpbtrf.f"
	    }

/*           Process the band matrix one diagonal block at a time. */

#line 346 "cpbtrf.f"
	    i__2 = *n;
#line 346 "cpbtrf.f"
	    i__1 = nb;
#line 346 "cpbtrf.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 347 "cpbtrf.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 347 "cpbtrf.f"
		ib = min(i__3,i__4);

/*              Factorize the diagonal block */

#line 351 "cpbtrf.f"
		i__3 = *ldab - 1;
#line 351 "cpbtrf.f"
		cpotf2_(uplo, &ib, &ab[i__ * ab_dim1 + 1], &i__3, &ii, (
			ftnlen)1);
#line 352 "cpbtrf.f"
		if (ii != 0) {
#line 353 "cpbtrf.f"
		    *info = i__ + ii - 1;
#line 354 "cpbtrf.f"
		    goto L150;
#line 355 "cpbtrf.f"
		}
#line 356 "cpbtrf.f"
		if (i__ + ib <= *n) {

/*                 Update the relevant part of the trailing submatrix. */
/*                 If A11 denotes the diagonal block which has just been */
/*                 factorized, then we need to update the remaining */
/*                 blocks in the diagram: */

/*                    A11 */
/*                    A21   A22 */
/*                    A31   A32   A33 */

/*                 The numbers of rows and columns in the partitioning */
/*                 are IB, I2, I3 respectively. The blocks A21, A22 and */
/*                 A32 are empty if IB = KD. The lower triangle of A31 */
/*                 lies outside the band. */

/* Computing MIN */
#line 372 "cpbtrf.f"
		    i__3 = *kd - ib, i__4 = *n - i__ - ib + 1;
#line 372 "cpbtrf.f"
		    i2 = min(i__3,i__4);
/* Computing MIN */
#line 373 "cpbtrf.f"
		    i__3 = ib, i__4 = *n - i__ - *kd + 1;
#line 373 "cpbtrf.f"
		    i3 = min(i__3,i__4);

#line 375 "cpbtrf.f"
		    if (i2 > 0) {

/*                    Update A21 */

#line 379 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 379 "cpbtrf.f"
			i__4 = *ldab - 1;
#line 379 "cpbtrf.f"
			ctrsm_("Right", "Lower", "Conjugate transpose", "Non"\
				"-unit", &i2, &ib, &c_b1, &ab[i__ * ab_dim1 + 
				1], &i__3, &ab[ib + 1 + i__ * ab_dim1], &i__4,
				 (ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)8);

/*                    Update A22 */

#line 386 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 386 "cpbtrf.f"
			i__4 = *ldab - 1;
#line 386 "cpbtrf.f"
			cherk_("Lower", "No transpose", &i2, &ib, &c_b21, &ab[
				ib + 1 + i__ * ab_dim1], &i__3, &c_b22, &ab[(
				i__ + ib) * ab_dim1 + 1], &i__4, (ftnlen)5, (
				ftnlen)12);
#line 389 "cpbtrf.f"
		    }

#line 391 "cpbtrf.f"
		    if (i3 > 0) {

/*                    Copy the upper triangle of A31 into the work array. */

#line 395 "cpbtrf.f"
			i__3 = ib;
#line 395 "cpbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 396 "cpbtrf.f"
			    i__4 = min(jj,i3);
#line 396 "cpbtrf.f"
			    for (ii = 1; ii <= i__4; ++ii) {
#line 397 "cpbtrf.f"
				i__5 = ii + jj * 33 - 34;
#line 397 "cpbtrf.f"
				i__6 = *kd + 1 - jj + ii + (jj + i__ - 1) * 
					ab_dim1;
#line 397 "cpbtrf.f"
				work[i__5].r = ab[i__6].r, work[i__5].i = ab[
					i__6].i;
#line 398 "cpbtrf.f"
/* L100: */
#line 398 "cpbtrf.f"
			    }
#line 399 "cpbtrf.f"
/* L110: */
#line 399 "cpbtrf.f"
			}

/*                    Update A31 (in the work array). */

#line 403 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 403 "cpbtrf.f"
			ctrsm_("Right", "Lower", "Conjugate transpose", "Non"\
				"-unit", &i3, &ib, &c_b1, &ab[i__ * ab_dim1 + 
				1], &i__3, work, &c__33, (ftnlen)5, (ftnlen)5,
				 (ftnlen)19, (ftnlen)8);

/*                    Update A32 */

#line 410 "cpbtrf.f"
			if (i2 > 0) {
#line 410 "cpbtrf.f"
			    z__1.r = -1., z__1.i = -0.;
#line 410 "cpbtrf.f"
			    i__3 = *ldab - 1;
#line 410 "cpbtrf.f"
			    i__4 = *ldab - 1;
#line 410 "cpbtrf.f"
			    cgemm_("No transpose", "Conjugate transpose", &i3,
				     &i2, &ib, &z__1, work, &c__33, &ab[ib + 
				    1 + i__ * ab_dim1], &i__3, &c_b1, &ab[*kd 
				    + 1 - ib + (i__ + ib) * ab_dim1], &i__4, (
				    ftnlen)12, (ftnlen)19);
#line 410 "cpbtrf.f"
			}

/*                    Update A33 */

#line 419 "cpbtrf.f"
			i__3 = *ldab - 1;
#line 419 "cpbtrf.f"
			cherk_("Lower", "No transpose", &i3, &ib, &c_b21, 
				work, &c__33, &c_b22, &ab[(i__ + *kd) * 
				ab_dim1 + 1], &i__3, (ftnlen)5, (ftnlen)12);

/*                    Copy the upper triangle of A31 back into place. */

#line 425 "cpbtrf.f"
			i__3 = ib;
#line 425 "cpbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 426 "cpbtrf.f"
			    i__4 = min(jj,i3);
#line 426 "cpbtrf.f"
			    for (ii = 1; ii <= i__4; ++ii) {
#line 427 "cpbtrf.f"
				i__5 = *kd + 1 - jj + ii + (jj + i__ - 1) * 
					ab_dim1;
#line 427 "cpbtrf.f"
				i__6 = ii + jj * 33 - 34;
#line 427 "cpbtrf.f"
				ab[i__5].r = work[i__6].r, ab[i__5].i = work[
					i__6].i;
#line 428 "cpbtrf.f"
/* L120: */
#line 428 "cpbtrf.f"
			    }
#line 429 "cpbtrf.f"
/* L130: */
#line 429 "cpbtrf.f"
			}
#line 430 "cpbtrf.f"
		    }
#line 431 "cpbtrf.f"
		}
#line 432 "cpbtrf.f"
/* L140: */
#line 432 "cpbtrf.f"
	    }
#line 433 "cpbtrf.f"
	}
#line 434 "cpbtrf.f"
    }
#line 435 "cpbtrf.f"
    return 0;

#line 437 "cpbtrf.f"
L150:
#line 438 "cpbtrf.f"
    return 0;

/*     End of CPBTRF */

} /* cpbtrf_ */

