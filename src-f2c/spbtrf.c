#line 1 "spbtrf.f"
/* spbtrf.f -- translated by f2c (version 20100827).
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

#line 1 "spbtrf.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 1.;
static doublereal c_b21 = -1.;
static integer c__33 = 33;

/* > \brief \b SPBTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPBTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spbtrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spbtrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spbtrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPBTRF( UPLO, N, KD, AB, LDAB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPBTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite band matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
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
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* >          On exit, if INFO = 0, the triangular factor U or L from the */
/* >          Cholesky factorization A = U**T*U or A = L*L**T of the band */
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

/* > \ingroup realOTHERcomputational */

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
/* Subroutine */ int spbtrf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, i2, i3, ib, nb, ii, jj;
    static doublereal work[1056]	/* was [33][32] */;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     strsm_(char *, char *, char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen), ssyrk_(char *, char *, integer *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), spbtf2_(char *, integer 
	    *, integer *, doublereal *, integer *, integer *, ftnlen), 
	    spotf2_(char *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
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

#line 187 "spbtrf.f"
    /* Parameter adjustments */
#line 187 "spbtrf.f"
    ab_dim1 = *ldab;
#line 187 "spbtrf.f"
    ab_offset = 1 + ab_dim1;
#line 187 "spbtrf.f"
    ab -= ab_offset;
#line 187 "spbtrf.f"

#line 187 "spbtrf.f"
    /* Function Body */
#line 187 "spbtrf.f"
    *info = 0;
#line 188 "spbtrf.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 190 "spbtrf.f"
	*info = -1;
#line 191 "spbtrf.f"
    } else if (*n < 0) {
#line 192 "spbtrf.f"
	*info = -2;
#line 193 "spbtrf.f"
    } else if (*kd < 0) {
#line 194 "spbtrf.f"
	*info = -3;
#line 195 "spbtrf.f"
    } else if (*ldab < *kd + 1) {
#line 196 "spbtrf.f"
	*info = -5;
#line 197 "spbtrf.f"
    }
#line 198 "spbtrf.f"
    if (*info != 0) {
#line 199 "spbtrf.f"
	i__1 = -(*info);
#line 199 "spbtrf.f"
	xerbla_("SPBTRF", &i__1, (ftnlen)6);
#line 200 "spbtrf.f"
	return 0;
#line 201 "spbtrf.f"
    }

/*     Quick return if possible */

#line 205 "spbtrf.f"
    if (*n == 0) {
#line 205 "spbtrf.f"
	return 0;
#line 205 "spbtrf.f"
    }

/*     Determine the block size for this environment */

#line 210 "spbtrf.f"
    nb = ilaenv_(&c__1, "SPBTRF", uplo, n, kd, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     The block size must not exceed the semi-bandwidth KD, and must not */
/*     exceed the limit set by the size of the local array WORK. */

#line 215 "spbtrf.f"
    nb = min(nb,32);

#line 217 "spbtrf.f"
    if (nb <= 1 || nb > *kd) {

/*        Use unblocked code */

#line 221 "spbtrf.f"
	spbtf2_(uplo, n, kd, &ab[ab_offset], ldab, info, (ftnlen)1);
#line 222 "spbtrf.f"
    } else {

/*        Use blocked code */

#line 226 "spbtrf.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*           Compute the Cholesky factorization of a symmetric band */
/*           matrix, given the upper triangle of the matrix in band */
/*           storage. */

/*           Zero the upper triangle of the work array. */

#line 234 "spbtrf.f"
	    i__1 = nb;
#line 234 "spbtrf.f"
	    for (j = 1; j <= i__1; ++j) {
#line 235 "spbtrf.f"
		i__2 = j - 1;
#line 235 "spbtrf.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 236 "spbtrf.f"
		    work[i__ + j * 33 - 34] = 0.;
#line 237 "spbtrf.f"
/* L10: */
#line 237 "spbtrf.f"
		}
#line 238 "spbtrf.f"
/* L20: */
#line 238 "spbtrf.f"
	    }

/*           Process the band matrix one diagonal block at a time. */

#line 242 "spbtrf.f"
	    i__1 = *n;
#line 242 "spbtrf.f"
	    i__2 = nb;
#line 242 "spbtrf.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 243 "spbtrf.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 243 "spbtrf.f"
		ib = min(i__3,i__4);

/*              Factorize the diagonal block */

#line 247 "spbtrf.f"
		i__3 = *ldab - 1;
#line 247 "spbtrf.f"
		spotf2_(uplo, &ib, &ab[*kd + 1 + i__ * ab_dim1], &i__3, &ii, (
			ftnlen)1);
#line 248 "spbtrf.f"
		if (ii != 0) {
#line 249 "spbtrf.f"
		    *info = i__ + ii - 1;
#line 250 "spbtrf.f"
		    goto L150;
#line 251 "spbtrf.f"
		}
#line 252 "spbtrf.f"
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
#line 268 "spbtrf.f"
		    i__3 = *kd - ib, i__4 = *n - i__ - ib + 1;
#line 268 "spbtrf.f"
		    i2 = min(i__3,i__4);
/* Computing MIN */
#line 269 "spbtrf.f"
		    i__3 = ib, i__4 = *n - i__ - *kd + 1;
#line 269 "spbtrf.f"
		    i3 = min(i__3,i__4);

#line 271 "spbtrf.f"
		    if (i2 > 0) {

/*                    Update A12 */

#line 275 "spbtrf.f"
			i__3 = *ldab - 1;
#line 275 "spbtrf.f"
			i__4 = *ldab - 1;
#line 275 "spbtrf.f"
			strsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
				&i2, &c_b18, &ab[*kd + 1 + i__ * ab_dim1], &
				i__3, &ab[*kd + 1 - ib + (i__ + ib) * ab_dim1]
				, &i__4, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
				ftnlen)8);

/*                    Update A22 */

#line 281 "spbtrf.f"
			i__3 = *ldab - 1;
#line 281 "spbtrf.f"
			i__4 = *ldab - 1;
#line 281 "spbtrf.f"
			ssyrk_("Upper", "Transpose", &i2, &ib, &c_b21, &ab[*
				kd + 1 - ib + (i__ + ib) * ab_dim1], &i__3, &
				c_b18, &ab[*kd + 1 + (i__ + ib) * ab_dim1], &
				i__4, (ftnlen)5, (ftnlen)9);
#line 284 "spbtrf.f"
		    }

#line 286 "spbtrf.f"
		    if (i3 > 0) {

/*                    Copy the lower triangle of A13 into the work array. */

#line 290 "spbtrf.f"
			i__3 = i3;
#line 290 "spbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 291 "spbtrf.f"
			    i__4 = ib;
#line 291 "spbtrf.f"
			    for (ii = jj; ii <= i__4; ++ii) {
#line 292 "spbtrf.f"
				work[ii + jj * 33 - 34] = ab[ii - jj + 1 + (
					jj + i__ + *kd - 1) * ab_dim1];
#line 293 "spbtrf.f"
/* L30: */
#line 293 "spbtrf.f"
			    }
#line 294 "spbtrf.f"
/* L40: */
#line 294 "spbtrf.f"
			}

/*                    Update A13 (in the work array). */

#line 298 "spbtrf.f"
			i__3 = *ldab - 1;
#line 298 "spbtrf.f"
			strsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
				&i3, &c_b18, &ab[*kd + 1 + i__ * ab_dim1], &
				i__3, work, &c__33, (ftnlen)4, (ftnlen)5, (
				ftnlen)9, (ftnlen)8);

/*                    Update A23 */

#line 304 "spbtrf.f"
			if (i2 > 0) {
#line 304 "spbtrf.f"
			    i__3 = *ldab - 1;
#line 304 "spbtrf.f"
			    i__4 = *ldab - 1;
#line 304 "spbtrf.f"
			    sgemm_("Transpose", "No Transpose", &i2, &i3, &ib,
				     &c_b21, &ab[*kd + 1 - ib + (i__ + ib) * 
				    ab_dim1], &i__3, work, &c__33, &c_b18, &
				    ab[ib + 1 + (i__ + *kd) * ab_dim1], &i__4,
				     (ftnlen)9, (ftnlen)12);
#line 304 "spbtrf.f"
			}

/*                    Update A33 */

#line 312 "spbtrf.f"
			i__3 = *ldab - 1;
#line 312 "spbtrf.f"
			ssyrk_("Upper", "Transpose", &i3, &ib, &c_b21, work, &
				c__33, &c_b18, &ab[*kd + 1 + (i__ + *kd) * 
				ab_dim1], &i__3, (ftnlen)5, (ftnlen)9);

/*                    Copy the lower triangle of A13 back into place. */

#line 318 "spbtrf.f"
			i__3 = i3;
#line 318 "spbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 319 "spbtrf.f"
			    i__4 = ib;
#line 319 "spbtrf.f"
			    for (ii = jj; ii <= i__4; ++ii) {
#line 320 "spbtrf.f"
				ab[ii - jj + 1 + (jj + i__ + *kd - 1) * 
					ab_dim1] = work[ii + jj * 33 - 34];
#line 321 "spbtrf.f"
/* L50: */
#line 321 "spbtrf.f"
			    }
#line 322 "spbtrf.f"
/* L60: */
#line 322 "spbtrf.f"
			}
#line 323 "spbtrf.f"
		    }
#line 324 "spbtrf.f"
		}
#line 325 "spbtrf.f"
/* L70: */
#line 325 "spbtrf.f"
	    }
#line 326 "spbtrf.f"
	} else {

/*           Compute the Cholesky factorization of a symmetric band */
/*           matrix, given the lower triangle of the matrix in band */
/*           storage. */

/*           Zero the lower triangle of the work array. */

#line 334 "spbtrf.f"
	    i__2 = nb;
#line 334 "spbtrf.f"
	    for (j = 1; j <= i__2; ++j) {
#line 335 "spbtrf.f"
		i__1 = nb;
#line 335 "spbtrf.f"
		for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 336 "spbtrf.f"
		    work[i__ + j * 33 - 34] = 0.;
#line 337 "spbtrf.f"
/* L80: */
#line 337 "spbtrf.f"
		}
#line 338 "spbtrf.f"
/* L90: */
#line 338 "spbtrf.f"
	    }

/*           Process the band matrix one diagonal block at a time. */

#line 342 "spbtrf.f"
	    i__2 = *n;
#line 342 "spbtrf.f"
	    i__1 = nb;
#line 342 "spbtrf.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 343 "spbtrf.f"
		i__3 = nb, i__4 = *n - i__ + 1;
#line 343 "spbtrf.f"
		ib = min(i__3,i__4);

/*              Factorize the diagonal block */

#line 347 "spbtrf.f"
		i__3 = *ldab - 1;
#line 347 "spbtrf.f"
		spotf2_(uplo, &ib, &ab[i__ * ab_dim1 + 1], &i__3, &ii, (
			ftnlen)1);
#line 348 "spbtrf.f"
		if (ii != 0) {
#line 349 "spbtrf.f"
		    *info = i__ + ii - 1;
#line 350 "spbtrf.f"
		    goto L150;
#line 351 "spbtrf.f"
		}
#line 352 "spbtrf.f"
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
#line 368 "spbtrf.f"
		    i__3 = *kd - ib, i__4 = *n - i__ - ib + 1;
#line 368 "spbtrf.f"
		    i2 = min(i__3,i__4);
/* Computing MIN */
#line 369 "spbtrf.f"
		    i__3 = ib, i__4 = *n - i__ - *kd + 1;
#line 369 "spbtrf.f"
		    i3 = min(i__3,i__4);

#line 371 "spbtrf.f"
		    if (i2 > 0) {

/*                    Update A21 */

#line 375 "spbtrf.f"
			i__3 = *ldab - 1;
#line 375 "spbtrf.f"
			i__4 = *ldab - 1;
#line 375 "spbtrf.f"
			strsm_("Right", "Lower", "Transpose", "Non-unit", &i2,
				 &ib, &c_b18, &ab[i__ * ab_dim1 + 1], &i__3, &
				ab[ib + 1 + i__ * ab_dim1], &i__4, (ftnlen)5, 
				(ftnlen)5, (ftnlen)9, (ftnlen)8);

/*                    Update A22 */

#line 381 "spbtrf.f"
			i__3 = *ldab - 1;
#line 381 "spbtrf.f"
			i__4 = *ldab - 1;
#line 381 "spbtrf.f"
			ssyrk_("Lower", "No Transpose", &i2, &ib, &c_b21, &ab[
				ib + 1 + i__ * ab_dim1], &i__3, &c_b18, &ab[(
				i__ + ib) * ab_dim1 + 1], &i__4, (ftnlen)5, (
				ftnlen)12);
#line 384 "spbtrf.f"
		    }

#line 386 "spbtrf.f"
		    if (i3 > 0) {

/*                    Copy the upper triangle of A31 into the work array. */

#line 390 "spbtrf.f"
			i__3 = ib;
#line 390 "spbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 391 "spbtrf.f"
			    i__4 = min(jj,i3);
#line 391 "spbtrf.f"
			    for (ii = 1; ii <= i__4; ++ii) {
#line 392 "spbtrf.f"
				work[ii + jj * 33 - 34] = ab[*kd + 1 - jj + 
					ii + (jj + i__ - 1) * ab_dim1];
#line 393 "spbtrf.f"
/* L100: */
#line 393 "spbtrf.f"
			    }
#line 394 "spbtrf.f"
/* L110: */
#line 394 "spbtrf.f"
			}

/*                    Update A31 (in the work array). */

#line 398 "spbtrf.f"
			i__3 = *ldab - 1;
#line 398 "spbtrf.f"
			strsm_("Right", "Lower", "Transpose", "Non-unit", &i3,
				 &ib, &c_b18, &ab[i__ * ab_dim1 + 1], &i__3, 
				work, &c__33, (ftnlen)5, (ftnlen)5, (ftnlen)9,
				 (ftnlen)8);

/*                    Update A32 */

#line 404 "spbtrf.f"
			if (i2 > 0) {
#line 404 "spbtrf.f"
			    i__3 = *ldab - 1;
#line 404 "spbtrf.f"
			    i__4 = *ldab - 1;
#line 404 "spbtrf.f"
			    sgemm_("No transpose", "Transpose", &i3, &i2, &ib,
				     &c_b21, work, &c__33, &ab[ib + 1 + i__ * 
				    ab_dim1], &i__3, &c_b18, &ab[*kd + 1 - ib 
				    + (i__ + ib) * ab_dim1], &i__4, (ftnlen)
				    12, (ftnlen)9);
#line 404 "spbtrf.f"
			}

/*                    Update A33 */

#line 412 "spbtrf.f"
			i__3 = *ldab - 1;
#line 412 "spbtrf.f"
			ssyrk_("Lower", "No Transpose", &i3, &ib, &c_b21, 
				work, &c__33, &c_b18, &ab[(i__ + *kd) * 
				ab_dim1 + 1], &i__3, (ftnlen)5, (ftnlen)12);

/*                    Copy the upper triangle of A31 back into place. */

#line 418 "spbtrf.f"
			i__3 = ib;
#line 418 "spbtrf.f"
			for (jj = 1; jj <= i__3; ++jj) {
#line 419 "spbtrf.f"
			    i__4 = min(jj,i3);
#line 419 "spbtrf.f"
			    for (ii = 1; ii <= i__4; ++ii) {
#line 420 "spbtrf.f"
				ab[*kd + 1 - jj + ii + (jj + i__ - 1) * 
					ab_dim1] = work[ii + jj * 33 - 34];
#line 421 "spbtrf.f"
/* L120: */
#line 421 "spbtrf.f"
			    }
#line 422 "spbtrf.f"
/* L130: */
#line 422 "spbtrf.f"
			}
#line 423 "spbtrf.f"
		    }
#line 424 "spbtrf.f"
		}
#line 425 "spbtrf.f"
/* L140: */
#line 425 "spbtrf.f"
	    }
#line 426 "spbtrf.f"
	}
#line 427 "spbtrf.f"
    }
#line 428 "spbtrf.f"
    return 0;

#line 430 "spbtrf.f"
L150:
#line 431 "spbtrf.f"
    return 0;

/*     End of SPBTRF */

} /* spbtrf_ */

