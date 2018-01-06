#line 1 "zhbtrd.f"
/* zhbtrd.f -- translated by f2c (version 20100827).
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

#line 1 "zhbtrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHBTRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHBTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbtrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbtrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbtrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KD, LDAB, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBTRD reduces a complex Hermitian band matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q**H * A * Q = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          = 'N':  do not form Q; */
/* >          = 'V':  form Q; */
/* >          = 'U':  update a matrix X, by forming X*Q. */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array.  The */
/* >          j-th column of A is stored in the j-th column of the array AB */
/* >          as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* >          On exit, the diagonal elements of AB are overwritten by the */
/* >          diagonal elements of the tridiagonal matrix T; if KD > 0, the */
/* >          elements on the first superdiagonal (if UPLO = 'U') or the */
/* >          first subdiagonal (if UPLO = 'L') are overwritten by the */
/* >          off-diagonal elements of T; the rest of AB is overwritten by */
/* >          values generated during the reduction. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          On entry, if VECT = 'U', then Q must contain an N-by-N */
/* >          matrix X; if VECT = 'N' or 'V', then Q need not be set. */
/* > */
/* >          On exit: */
/* >          if VECT = 'V', Q contains the N-by-N unitary matrix Q; */
/* >          if VECT = 'U', Q contains the product X*Q; */
/* >          if VECT = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. */
/* >          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by Linda Kaufman, Bell Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublecomplex *q, integer *ldq, doublecomplex *work, integer *info, 
	ftnlen vect_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublecomplex t;
    static integer i2, j1, j2, nq, nr, kd1, ibl, iqb, kdn, jin, nrt, kdm1, 
	    inca, jend, lend, jinc;
    static doublereal abst;
    static integer incx, last;
    static doublecomplex temp;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer j1end, j1inc, iqend;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static logical initq, wantq, upper;
    extern /* Subroutine */ int zlar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *);
    static integer iqaend;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlacgv_(
	    integer *, doublecomplex *, integer *), zlaset_(char *, integer *,
	     integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), zlartg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *), zlargv_(integer *
	    , doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *), zlartv_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 212 "zhbtrd.f"
    /* Parameter adjustments */
#line 212 "zhbtrd.f"
    ab_dim1 = *ldab;
#line 212 "zhbtrd.f"
    ab_offset = 1 + ab_dim1;
#line 212 "zhbtrd.f"
    ab -= ab_offset;
#line 212 "zhbtrd.f"
    --d__;
#line 212 "zhbtrd.f"
    --e;
#line 212 "zhbtrd.f"
    q_dim1 = *ldq;
#line 212 "zhbtrd.f"
    q_offset = 1 + q_dim1;
#line 212 "zhbtrd.f"
    q -= q_offset;
#line 212 "zhbtrd.f"
    --work;
#line 212 "zhbtrd.f"

#line 212 "zhbtrd.f"
    /* Function Body */
#line 212 "zhbtrd.f"
    initq = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 213 "zhbtrd.f"
    wantq = initq || lsame_(vect, "U", (ftnlen)1, (ftnlen)1);
#line 214 "zhbtrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 215 "zhbtrd.f"
    kd1 = *kd + 1;
#line 216 "zhbtrd.f"
    kdm1 = *kd - 1;
#line 217 "zhbtrd.f"
    incx = *ldab - 1;
#line 218 "zhbtrd.f"
    iqend = 1;

#line 220 "zhbtrd.f"
    *info = 0;
#line 221 "zhbtrd.f"
    if (! wantq && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 222 "zhbtrd.f"
	*info = -1;
#line 223 "zhbtrd.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 224 "zhbtrd.f"
	*info = -2;
#line 225 "zhbtrd.f"
    } else if (*n < 0) {
#line 226 "zhbtrd.f"
	*info = -3;
#line 227 "zhbtrd.f"
    } else if (*kd < 0) {
#line 228 "zhbtrd.f"
	*info = -4;
#line 229 "zhbtrd.f"
    } else if (*ldab < kd1) {
#line 230 "zhbtrd.f"
	*info = -6;
#line 231 "zhbtrd.f"
    } else if (*ldq < max(1,*n) && wantq) {
#line 232 "zhbtrd.f"
	*info = -10;
#line 233 "zhbtrd.f"
    }
#line 234 "zhbtrd.f"
    if (*info != 0) {
#line 235 "zhbtrd.f"
	i__1 = -(*info);
#line 235 "zhbtrd.f"
	xerbla_("ZHBTRD", &i__1, (ftnlen)6);
#line 236 "zhbtrd.f"
	return 0;
#line 237 "zhbtrd.f"
    }

/*     Quick return if possible */

#line 241 "zhbtrd.f"
    if (*n == 0) {
#line 241 "zhbtrd.f"
	return 0;
#line 241 "zhbtrd.f"
    }

/*     Initialize Q to the unit matrix, if needed */

#line 246 "zhbtrd.f"
    if (initq) {
#line 246 "zhbtrd.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 246 "zhbtrd.f"
    }

/*     Wherever possible, plane rotations are generated and applied in */
/*     vector operations of length NR over the index set J1:J2:KD1. */

/*     The real cosines and complex sines of the plane rotations are */
/*     stored in the arrays D and WORK. */

#line 255 "zhbtrd.f"
    inca = kd1 * *ldab;
/* Computing MIN */
#line 256 "zhbtrd.f"
    i__1 = *n - 1;
#line 256 "zhbtrd.f"
    kdn = min(i__1,*kd);
#line 257 "zhbtrd.f"
    if (upper) {

#line 259 "zhbtrd.f"
	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working with */
/*           the upper triangle */

#line 264 "zhbtrd.f"
	    nr = 0;
#line 265 "zhbtrd.f"
	    j1 = kdn + 2;
#line 266 "zhbtrd.f"
	    j2 = 1;

#line 268 "zhbtrd.f"
	    i__1 = kd1 + ab_dim1;
#line 268 "zhbtrd.f"
	    i__2 = kd1 + ab_dim1;
#line 268 "zhbtrd.f"
	    d__1 = ab[i__2].r;
#line 268 "zhbtrd.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 269 "zhbtrd.f"
	    i__1 = *n - 2;
#line 269 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th row of matrix to tridiagonal form */

#line 273 "zhbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 274 "zhbtrd.f"
		    j1 += kdn;
#line 275 "zhbtrd.f"
		    j2 += kdn;

#line 277 "zhbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 282 "zhbtrd.f"
			zlargv_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    ZLARTV or ZROT is used */

#line 291 "zhbtrd.f"
			if (nr >= (*kd << 1) - 1) {
#line 292 "zhbtrd.f"
			    i__2 = *kd - 1;
#line 292 "zhbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 293 "zhbtrd.f"
				zlartv_(&nr, &ab[l + 1 + (j1 - 1) * ab_dim1], 
					&inca, &ab[l + j1 * ab_dim1], &inca, &
					d__[j1], &work[j1], &kd1);
#line 296 "zhbtrd.f"
/* L10: */
#line 296 "zhbtrd.f"
			    }

#line 298 "zhbtrd.f"
			} else {
#line 299 "zhbtrd.f"
			    jend = j1 + (nr - 1) * kd1;
#line 300 "zhbtrd.f"
			    i__2 = jend;
#line 300 "zhbtrd.f"
			    i__3 = kd1;
#line 300 "zhbtrd.f"
			    for (jinc = j1; i__3 < 0 ? jinc >= i__2 : jinc <= 
				    i__2; jinc += i__3) {
#line 301 "zhbtrd.f"
				zrot_(&kdm1, &ab[(jinc - 1) * ab_dim1 + 2], &
					c__1, &ab[jinc * ab_dim1 + 1], &c__1, 
					&d__[jinc], &work[jinc]);
#line 304 "zhbtrd.f"
/* L20: */
#line 304 "zhbtrd.f"
			    }
#line 305 "zhbtrd.f"
			}
#line 306 "zhbtrd.f"
		    }


#line 309 "zhbtrd.f"
		    if (k > 2) {
#line 310 "zhbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i,i+k-1) */
/*                       within the band */

#line 315 "zhbtrd.f"
			    zlartg_(&ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1]
				    , &ab[*kd - k + 2 + (i__ + k - 1) * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 318 "zhbtrd.f"
			    i__3 = *kd - k + 3 + (i__ + k - 2) * ab_dim1;
#line 318 "zhbtrd.f"
			    ab[i__3].r = temp.r, ab[i__3].i = temp.i;

/*                       apply rotation from the right */

#line 322 "zhbtrd.f"
			    i__3 = k - 3;
#line 322 "zhbtrd.f"
			    zrot_(&i__3, &ab[*kd - k + 4 + (i__ + k - 2) * 
				    ab_dim1], &c__1, &ab[*kd - k + 3 + (i__ + 
				    k - 1) * ab_dim1], &c__1, &d__[i__ + k - 
				    1], &work[i__ + k - 1]);
#line 325 "zhbtrd.f"
			}
#line 326 "zhbtrd.f"
			++nr;
#line 327 "zhbtrd.f"
			j1 = j1 - kdn - 1;
#line 328 "zhbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 333 "zhbtrd.f"
		    if (nr > 0) {
#line 333 "zhbtrd.f"
			zlar2v_(&nr, &ab[kd1 + (j1 - 1) * ab_dim1], &ab[kd1 + 
				j1 * ab_dim1], &ab[*kd + j1 * ab_dim1], &inca,
				 &d__[j1], &work[j1], &kd1);
#line 333 "zhbtrd.f"
		    }

/*                 apply plane rotations from the left */

#line 340 "zhbtrd.f"
		    if (nr > 0) {
#line 341 "zhbtrd.f"
			zlacgv_(&nr, &work[j1], &kd1);
#line 342 "zhbtrd.f"
			if ((*kd << 1) - 1 < nr) {

/*                    Dependent on the the number of diagonals either */
/*                    ZLARTV or ZROT is used */

#line 347 "zhbtrd.f"
			    i__3 = *kd - 1;
#line 347 "zhbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 348 "zhbtrd.f"
				if (j2 + l > *n) {
#line 349 "zhbtrd.f"
				    nrt = nr - 1;
#line 350 "zhbtrd.f"
				} else {
#line 351 "zhbtrd.f"
				    nrt = nr;
#line 352 "zhbtrd.f"
				}
#line 353 "zhbtrd.f"
				if (nrt > 0) {
#line 353 "zhbtrd.f"
				    zlartv_(&nrt, &ab[*kd - l + (j1 + l) * 
					    ab_dim1], &inca, &ab[*kd - l + 1 
					    + (j1 + l) * ab_dim1], &inca, &
					    d__[j1], &work[j1], &kd1);
#line 353 "zhbtrd.f"
				}
#line 357 "zhbtrd.f"
/* L30: */
#line 357 "zhbtrd.f"
			    }
#line 358 "zhbtrd.f"
			} else {
#line 359 "zhbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 360 "zhbtrd.f"
			    if (j1end >= j1) {
#line 361 "zhbtrd.f"
				i__3 = j1end;
#line 361 "zhbtrd.f"
				i__2 = kd1;
#line 361 "zhbtrd.f"
				for (jin = j1; i__2 < 0 ? jin >= i__3 : jin <=
					 i__3; jin += i__2) {
#line 362 "zhbtrd.f"
				    i__4 = *kd - 1;
#line 362 "zhbtrd.f"
				    zrot_(&i__4, &ab[*kd - 1 + (jin + 1) * 
					    ab_dim1], &incx, &ab[*kd + (jin + 
					    1) * ab_dim1], &incx, &d__[jin], &
					    work[jin]);
#line 365 "zhbtrd.f"
/* L40: */
#line 365 "zhbtrd.f"
				}
#line 366 "zhbtrd.f"
			    }
/* Computing MIN */
#line 367 "zhbtrd.f"
			    i__2 = kdm1, i__3 = *n - j2;
#line 367 "zhbtrd.f"
			    lend = min(i__2,i__3);
#line 368 "zhbtrd.f"
			    last = j1end + kd1;
#line 369 "zhbtrd.f"
			    if (lend > 0) {
#line 369 "zhbtrd.f"
				zrot_(&lend, &ab[*kd - 1 + (last + 1) * 
					ab_dim1], &incx, &ab[*kd + (last + 1) 
					* ab_dim1], &incx, &d__[last], &work[
					last]);
#line 369 "zhbtrd.f"
			    }
#line 373 "zhbtrd.f"
			}
#line 374 "zhbtrd.f"
		    }

#line 376 "zhbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 380 "zhbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 385 "zhbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 386 "zhbtrd.f"
			    i__2 = 0, i__3 = k - 3;
#line 386 "zhbtrd.f"
			    i2 = max(i__2,i__3);
#line 387 "zhbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 388 "zhbtrd.f"
			    if (k == 2) {
#line 388 "zhbtrd.f"
				iqaend += *kd;
#line 388 "zhbtrd.f"
			    }
#line 390 "zhbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 391 "zhbtrd.f"
			    i__2 = j2;
#line 391 "zhbtrd.f"
			    i__3 = kd1;
#line 391 "zhbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 392 "zhbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 393 "zhbtrd.f"
				++i2;
/* Computing MAX */
#line 394 "zhbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 394 "zhbtrd.f"
				iqb = max(i__4,i__5);
#line 395 "zhbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 396 "zhbtrd.f"
				i__4 = iqaend + *kd;
#line 396 "zhbtrd.f"
				iqaend = min(i__4,iqend);
#line 397 "zhbtrd.f"
				d_cnjg(&z__1, &work[j]);
#line 397 "zhbtrd.f"
				zrot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&z__1);
#line 399 "zhbtrd.f"
/* L50: */
#line 399 "zhbtrd.f"
			    }
#line 400 "zhbtrd.f"
			} else {

#line 402 "zhbtrd.f"
			    i__3 = j2;
#line 402 "zhbtrd.f"
			    i__2 = kd1;
#line 402 "zhbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 403 "zhbtrd.f"
				d_cnjg(&z__1, &work[j]);
#line 403 "zhbtrd.f"
				zrot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					z__1);
#line 405 "zhbtrd.f"
/* L60: */
#line 405 "zhbtrd.f"
			    }
#line 406 "zhbtrd.f"
			}

#line 408 "zhbtrd.f"
		    }

#line 410 "zhbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 414 "zhbtrd.f"
			--nr;
#line 415 "zhbtrd.f"
			j2 = j2 - kdn - 1;
#line 416 "zhbtrd.f"
		    }

#line 418 "zhbtrd.f"
		    i__2 = j2;
#line 418 "zhbtrd.f"
		    i__3 = kd1;
#line 418 "zhbtrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) 
			    {

/*                    create nonzero element a(j-1,j+kd) outside the band */
/*                    and store it in WORK */

#line 423 "zhbtrd.f"
			i__4 = j + *kd;
#line 423 "zhbtrd.f"
			i__5 = j;
#line 423 "zhbtrd.f"
			i__6 = (j + *kd) * ab_dim1 + 1;
#line 423 "zhbtrd.f"
			z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * 
				ab[i__6].i, z__1.i = work[i__5].r * ab[i__6]
				.i + work[i__5].i * ab[i__6].r;
#line 423 "zhbtrd.f"
			work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 424 "zhbtrd.f"
			i__4 = (j + *kd) * ab_dim1 + 1;
#line 424 "zhbtrd.f"
			i__5 = j;
#line 424 "zhbtrd.f"
			i__6 = (j + *kd) * ab_dim1 + 1;
#line 424 "zhbtrd.f"
			z__1.r = d__[i__5] * ab[i__6].r, z__1.i = d__[i__5] * 
				ab[i__6].i;
#line 424 "zhbtrd.f"
			ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 425 "zhbtrd.f"
/* L70: */
#line 425 "zhbtrd.f"
		    }
#line 426 "zhbtrd.f"
/* L80: */
#line 426 "zhbtrd.f"
		}
#line 427 "zhbtrd.f"
/* L90: */
#line 427 "zhbtrd.f"
	    }
#line 428 "zhbtrd.f"
	}

#line 430 "zhbtrd.f"
	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

#line 434 "zhbtrd.f"
	    i__1 = *n - 1;
#line 434 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "zhbtrd.f"
		i__3 = *kd + (i__ + 1) * ab_dim1;
#line 435 "zhbtrd.f"
		t.r = ab[i__3].r, t.i = ab[i__3].i;
#line 436 "zhbtrd.f"
		abst = z_abs(&t);
#line 437 "zhbtrd.f"
		i__3 = *kd + (i__ + 1) * ab_dim1;
#line 437 "zhbtrd.f"
		ab[i__3].r = abst, ab[i__3].i = 0.;
#line 438 "zhbtrd.f"
		e[i__] = abst;
#line 439 "zhbtrd.f"
		if (abst != 0.) {
#line 440 "zhbtrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 440 "zhbtrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 441 "zhbtrd.f"
		} else {
#line 442 "zhbtrd.f"
		    t.r = 1., t.i = 0.;
#line 443 "zhbtrd.f"
		}
#line 444 "zhbtrd.f"
		if (i__ < *n - 1) {
#line 444 "zhbtrd.f"
		    i__3 = *kd + (i__ + 2) * ab_dim1;
#line 444 "zhbtrd.f"
		    i__2 = *kd + (i__ + 2) * ab_dim1;
#line 444 "zhbtrd.f"
		    z__1.r = ab[i__2].r * t.r - ab[i__2].i * t.i, z__1.i = ab[
			    i__2].r * t.i + ab[i__2].i * t.r;
#line 444 "zhbtrd.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 444 "zhbtrd.f"
		}
#line 446 "zhbtrd.f"
		if (wantq) {
#line 447 "zhbtrd.f"
		    d_cnjg(&z__1, &t);
#line 447 "zhbtrd.f"
		    zscal_(n, &z__1, &q[(i__ + 1) * q_dim1 + 1], &c__1);
#line 448 "zhbtrd.f"
		}
#line 449 "zhbtrd.f"
/* L100: */
#line 449 "zhbtrd.f"
	    }
#line 450 "zhbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 454 "zhbtrd.f"
	    i__1 = *n - 1;
#line 454 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 455 "zhbtrd.f"
		e[i__] = 0.;
#line 456 "zhbtrd.f"
/* L110: */
#line 456 "zhbtrd.f"
	    }
#line 457 "zhbtrd.f"
	}

/*        copy diagonal elements to D */

#line 461 "zhbtrd.f"
	i__1 = *n;
#line 461 "zhbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 462 "zhbtrd.f"
	    i__3 = i__;
#line 462 "zhbtrd.f"
	    i__2 = kd1 + i__ * ab_dim1;
#line 462 "zhbtrd.f"
	    d__[i__3] = ab[i__2].r;
#line 463 "zhbtrd.f"
/* L120: */
#line 463 "zhbtrd.f"
	}

#line 465 "zhbtrd.f"
    } else {

#line 467 "zhbtrd.f"
	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working with */
/*           the lower triangle */

#line 472 "zhbtrd.f"
	    nr = 0;
#line 473 "zhbtrd.f"
	    j1 = kdn + 2;
#line 474 "zhbtrd.f"
	    j2 = 1;

#line 476 "zhbtrd.f"
	    i__1 = ab_dim1 + 1;
#line 476 "zhbtrd.f"
	    i__3 = ab_dim1 + 1;
#line 476 "zhbtrd.f"
	    d__1 = ab[i__3].r;
#line 476 "zhbtrd.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 477 "zhbtrd.f"
	    i__1 = *n - 2;
#line 477 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th column of matrix to tridiagonal form */

#line 481 "zhbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 482 "zhbtrd.f"
		    j1 += kdn;
#line 483 "zhbtrd.f"
		    j2 += kdn;

#line 485 "zhbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 490 "zhbtrd.f"
			zlargv_(&nr, &ab[kd1 + (j1 - kd1) * ab_dim1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply plane rotations from one side */


/*                    Dependent on the the number of diagonals either */
/*                    ZLARTV or ZROT is used */

#line 499 "zhbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 500 "zhbtrd.f"
			    i__3 = *kd - 1;
#line 500 "zhbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 501 "zhbtrd.f"
				zlartv_(&nr, &ab[kd1 - l + (j1 - kd1 + l) * 
					ab_dim1], &inca, &ab[kd1 - l + 1 + (
					j1 - kd1 + l) * ab_dim1], &inca, &d__[
					j1], &work[j1], &kd1);
#line 504 "zhbtrd.f"
/* L130: */
#line 504 "zhbtrd.f"
			    }
#line 505 "zhbtrd.f"
			} else {
#line 506 "zhbtrd.f"
			    jend = j1 + kd1 * (nr - 1);
#line 507 "zhbtrd.f"
			    i__3 = jend;
#line 507 "zhbtrd.f"
			    i__2 = kd1;
#line 507 "zhbtrd.f"
			    for (jinc = j1; i__2 < 0 ? jinc >= i__3 : jinc <= 
				    i__3; jinc += i__2) {
#line 508 "zhbtrd.f"
				zrot_(&kdm1, &ab[*kd + (jinc - *kd) * ab_dim1]
					, &incx, &ab[kd1 + (jinc - *kd) * 
					ab_dim1], &incx, &d__[jinc], &work[
					jinc]);
#line 511 "zhbtrd.f"
/* L140: */
#line 511 "zhbtrd.f"
			    }
#line 512 "zhbtrd.f"
			}

#line 514 "zhbtrd.f"
		    }

#line 516 "zhbtrd.f"
		    if (k > 2) {
#line 517 "zhbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i+k-1,i) */
/*                       within the band */

#line 522 "zhbtrd.f"
			    zlartg_(&ab[k - 1 + i__ * ab_dim1], &ab[k + i__ * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 524 "zhbtrd.f"
			    i__2 = k - 1 + i__ * ab_dim1;
#line 524 "zhbtrd.f"
			    ab[i__2].r = temp.r, ab[i__2].i = temp.i;

/*                       apply rotation from the left */

#line 528 "zhbtrd.f"
			    i__2 = k - 3;
#line 528 "zhbtrd.f"
			    i__3 = *ldab - 1;
#line 528 "zhbtrd.f"
			    i__4 = *ldab - 1;
#line 528 "zhbtrd.f"
			    zrot_(&i__2, &ab[k - 2 + (i__ + 1) * ab_dim1], &
				    i__3, &ab[k - 1 + (i__ + 1) * ab_dim1], &
				    i__4, &d__[i__ + k - 1], &work[i__ + k - 
				    1]);
#line 531 "zhbtrd.f"
			}
#line 532 "zhbtrd.f"
			++nr;
#line 533 "zhbtrd.f"
			j1 = j1 - kdn - 1;
#line 534 "zhbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 539 "zhbtrd.f"
		    if (nr > 0) {
#line 539 "zhbtrd.f"
			zlar2v_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &ab[j1 * 
				ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 2], &
				inca, &d__[j1], &work[j1], &kd1);
#line 539 "zhbtrd.f"
		    }

/*                 apply plane rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    ZLARTV or ZROT is used */

#line 550 "zhbtrd.f"
		    if (nr > 0) {
#line 551 "zhbtrd.f"
			zlacgv_(&nr, &work[j1], &kd1);
#line 552 "zhbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 553 "zhbtrd.f"
			    i__2 = *kd - 1;
#line 553 "zhbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 554 "zhbtrd.f"
				if (j2 + l > *n) {
#line 555 "zhbtrd.f"
				    nrt = nr - 1;
#line 556 "zhbtrd.f"
				} else {
#line 557 "zhbtrd.f"
				    nrt = nr;
#line 558 "zhbtrd.f"
				}
#line 559 "zhbtrd.f"
				if (nrt > 0) {
#line 559 "zhbtrd.f"
				    zlartv_(&nrt, &ab[l + 2 + (j1 - 1) * 
					    ab_dim1], &inca, &ab[l + 1 + j1 * 
					    ab_dim1], &inca, &d__[j1], &work[
					    j1], &kd1);
#line 559 "zhbtrd.f"
				}
#line 563 "zhbtrd.f"
/* L150: */
#line 563 "zhbtrd.f"
			    }
#line 564 "zhbtrd.f"
			} else {
#line 565 "zhbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 566 "zhbtrd.f"
			    if (j1end >= j1) {
#line 567 "zhbtrd.f"
				i__2 = j1end;
#line 567 "zhbtrd.f"
				i__3 = kd1;
#line 567 "zhbtrd.f"
				for (j1inc = j1; i__3 < 0 ? j1inc >= i__2 : 
					j1inc <= i__2; j1inc += i__3) {
#line 568 "zhbtrd.f"
				    zrot_(&kdm1, &ab[(j1inc - 1) * ab_dim1 + 
					    3], &c__1, &ab[j1inc * ab_dim1 + 
					    2], &c__1, &d__[j1inc], &work[
					    j1inc]);
#line 571 "zhbtrd.f"
/* L160: */
#line 571 "zhbtrd.f"
				}
#line 572 "zhbtrd.f"
			    }
/* Computing MIN */
#line 573 "zhbtrd.f"
			    i__3 = kdm1, i__2 = *n - j2;
#line 573 "zhbtrd.f"
			    lend = min(i__3,i__2);
#line 574 "zhbtrd.f"
			    last = j1end + kd1;
#line 575 "zhbtrd.f"
			    if (lend > 0) {
#line 575 "zhbtrd.f"
				zrot_(&lend, &ab[(last - 1) * ab_dim1 + 3], &
					c__1, &ab[last * ab_dim1 + 2], &c__1, 
					&d__[last], &work[last]);
#line 575 "zhbtrd.f"
			    }
#line 579 "zhbtrd.f"
			}
#line 580 "zhbtrd.f"
		    }



#line 584 "zhbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 588 "zhbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 593 "zhbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 594 "zhbtrd.f"
			    i__3 = 0, i__2 = k - 3;
#line 594 "zhbtrd.f"
			    i2 = max(i__3,i__2);
#line 595 "zhbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 596 "zhbtrd.f"
			    if (k == 2) {
#line 596 "zhbtrd.f"
				iqaend += *kd;
#line 596 "zhbtrd.f"
			    }
#line 598 "zhbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 599 "zhbtrd.f"
			    i__3 = j2;
#line 599 "zhbtrd.f"
			    i__2 = kd1;
#line 599 "zhbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 600 "zhbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 601 "zhbtrd.f"
				++i2;
/* Computing MAX */
#line 602 "zhbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 602 "zhbtrd.f"
				iqb = max(i__4,i__5);
#line 603 "zhbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 604 "zhbtrd.f"
				i__4 = iqaend + *kd;
#line 604 "zhbtrd.f"
				iqaend = min(i__4,iqend);
#line 605 "zhbtrd.f"
				zrot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&work[j]);
#line 607 "zhbtrd.f"
/* L170: */
#line 607 "zhbtrd.f"
			    }
#line 608 "zhbtrd.f"
			} else {

#line 610 "zhbtrd.f"
			    i__2 = j2;
#line 610 "zhbtrd.f"
			    i__3 = kd1;
#line 610 "zhbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 611 "zhbtrd.f"
				zrot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					work[j]);
#line 613 "zhbtrd.f"
/* L180: */
#line 613 "zhbtrd.f"
			    }
#line 614 "zhbtrd.f"
			}
#line 615 "zhbtrd.f"
		    }

#line 617 "zhbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 621 "zhbtrd.f"
			--nr;
#line 622 "zhbtrd.f"
			j2 = j2 - kdn - 1;
#line 623 "zhbtrd.f"
		    }

#line 625 "zhbtrd.f"
		    i__3 = j2;
#line 625 "zhbtrd.f"
		    i__2 = kd1;
#line 625 "zhbtrd.f"
		    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) 
			    {

/*                    create nonzero element a(j+kd,j-1) outside the */
/*                    band and store it in WORK */

#line 630 "zhbtrd.f"
			i__4 = j + *kd;
#line 630 "zhbtrd.f"
			i__5 = j;
#line 630 "zhbtrd.f"
			i__6 = kd1 + j * ab_dim1;
#line 630 "zhbtrd.f"
			z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * 
				ab[i__6].i, z__1.i = work[i__5].r * ab[i__6]
				.i + work[i__5].i * ab[i__6].r;
#line 630 "zhbtrd.f"
			work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 631 "zhbtrd.f"
			i__4 = kd1 + j * ab_dim1;
#line 631 "zhbtrd.f"
			i__5 = j;
#line 631 "zhbtrd.f"
			i__6 = kd1 + j * ab_dim1;
#line 631 "zhbtrd.f"
			z__1.r = d__[i__5] * ab[i__6].r, z__1.i = d__[i__5] * 
				ab[i__6].i;
#line 631 "zhbtrd.f"
			ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 632 "zhbtrd.f"
/* L190: */
#line 632 "zhbtrd.f"
		    }
#line 633 "zhbtrd.f"
/* L200: */
#line 633 "zhbtrd.f"
		}
#line 634 "zhbtrd.f"
/* L210: */
#line 634 "zhbtrd.f"
	    }
#line 635 "zhbtrd.f"
	}

#line 637 "zhbtrd.f"
	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

#line 641 "zhbtrd.f"
	    i__1 = *n - 1;
#line 641 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 642 "zhbtrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 642 "zhbtrd.f"
		t.r = ab[i__2].r, t.i = ab[i__2].i;
#line 643 "zhbtrd.f"
		abst = z_abs(&t);
#line 644 "zhbtrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 644 "zhbtrd.f"
		ab[i__2].r = abst, ab[i__2].i = 0.;
#line 645 "zhbtrd.f"
		e[i__] = abst;
#line 646 "zhbtrd.f"
		if (abst != 0.) {
#line 647 "zhbtrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 647 "zhbtrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 648 "zhbtrd.f"
		} else {
#line 649 "zhbtrd.f"
		    t.r = 1., t.i = 0.;
#line 650 "zhbtrd.f"
		}
#line 651 "zhbtrd.f"
		if (i__ < *n - 1) {
#line 651 "zhbtrd.f"
		    i__2 = (i__ + 1) * ab_dim1 + 2;
#line 651 "zhbtrd.f"
		    i__3 = (i__ + 1) * ab_dim1 + 2;
#line 651 "zhbtrd.f"
		    z__1.r = ab[i__3].r * t.r - ab[i__3].i * t.i, z__1.i = ab[
			    i__3].r * t.i + ab[i__3].i * t.r;
#line 651 "zhbtrd.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 651 "zhbtrd.f"
		}
#line 653 "zhbtrd.f"
		if (wantq) {
#line 654 "zhbtrd.f"
		    zscal_(n, &t, &q[(i__ + 1) * q_dim1 + 1], &c__1);
#line 655 "zhbtrd.f"
		}
#line 656 "zhbtrd.f"
/* L220: */
#line 656 "zhbtrd.f"
	    }
#line 657 "zhbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 661 "zhbtrd.f"
	    i__1 = *n - 1;
#line 661 "zhbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 662 "zhbtrd.f"
		e[i__] = 0.;
#line 663 "zhbtrd.f"
/* L230: */
#line 663 "zhbtrd.f"
	    }
#line 664 "zhbtrd.f"
	}

/*        copy diagonal elements to D */

#line 668 "zhbtrd.f"
	i__1 = *n;
#line 668 "zhbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 669 "zhbtrd.f"
	    i__2 = i__;
#line 669 "zhbtrd.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 669 "zhbtrd.f"
	    d__[i__2] = ab[i__3].r;
#line 670 "zhbtrd.f"
/* L240: */
#line 670 "zhbtrd.f"
	}
#line 671 "zhbtrd.f"
    }

#line 673 "zhbtrd.f"
    return 0;

/*     End of ZHBTRD */

} /* zhbtrd_ */

