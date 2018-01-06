#line 1 "dsbtrd.f"
/* dsbtrd.f -- translated by f2c (version 20100827).
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

#line 1 "dsbtrd.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__1 = 1;

/* > \brief \b DSBTRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSBTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbtrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbtrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbtrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KD, LDAB, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBTRD reduces a real symmetric band matrix A to symmetric */
/* > tridiagonal form T by an orthogonal similarity transformation: */
/* > Q**T * A * Q = T. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          On entry, if VECT = 'U', then Q must contain an N-by-N */
/* >          matrix X; if VECT = 'N' or 'V', then Q need not be set. */
/* > */
/* >          On exit: */
/* >          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q; */
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
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by Linda Kaufman, Bell Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublereal *q, integer *ldq, doublereal *work, integer *info, ftnlen 
	vect_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, j, k, l, i2, j1, j2, nq, nr, kd1, ibl, iqb, kdn, jin, 
	    nrt, kdm1, inca, jend, lend, jinc, incx, last;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer j1end, j1inc, iqend;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical initq, wantq, upper;
    extern /* Subroutine */ int dlar2v_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer iqaend;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen), dlargv_(
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlartv_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 208 "dsbtrd.f"
    /* Parameter adjustments */
#line 208 "dsbtrd.f"
    ab_dim1 = *ldab;
#line 208 "dsbtrd.f"
    ab_offset = 1 + ab_dim1;
#line 208 "dsbtrd.f"
    ab -= ab_offset;
#line 208 "dsbtrd.f"
    --d__;
#line 208 "dsbtrd.f"
    --e;
#line 208 "dsbtrd.f"
    q_dim1 = *ldq;
#line 208 "dsbtrd.f"
    q_offset = 1 + q_dim1;
#line 208 "dsbtrd.f"
    q -= q_offset;
#line 208 "dsbtrd.f"
    --work;
#line 208 "dsbtrd.f"

#line 208 "dsbtrd.f"
    /* Function Body */
#line 208 "dsbtrd.f"
    initq = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 209 "dsbtrd.f"
    wantq = initq || lsame_(vect, "U", (ftnlen)1, (ftnlen)1);
#line 210 "dsbtrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 211 "dsbtrd.f"
    kd1 = *kd + 1;
#line 212 "dsbtrd.f"
    kdm1 = *kd - 1;
#line 213 "dsbtrd.f"
    incx = *ldab - 1;
#line 214 "dsbtrd.f"
    iqend = 1;

#line 216 "dsbtrd.f"
    *info = 0;
#line 217 "dsbtrd.f"
    if (! wantq && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 218 "dsbtrd.f"
	*info = -1;
#line 219 "dsbtrd.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 220 "dsbtrd.f"
	*info = -2;
#line 221 "dsbtrd.f"
    } else if (*n < 0) {
#line 222 "dsbtrd.f"
	*info = -3;
#line 223 "dsbtrd.f"
    } else if (*kd < 0) {
#line 224 "dsbtrd.f"
	*info = -4;
#line 225 "dsbtrd.f"
    } else if (*ldab < kd1) {
#line 226 "dsbtrd.f"
	*info = -6;
#line 227 "dsbtrd.f"
    } else if (*ldq < max(1,*n) && wantq) {
#line 228 "dsbtrd.f"
	*info = -10;
#line 229 "dsbtrd.f"
    }
#line 230 "dsbtrd.f"
    if (*info != 0) {
#line 231 "dsbtrd.f"
	i__1 = -(*info);
#line 231 "dsbtrd.f"
	xerbla_("DSBTRD", &i__1, (ftnlen)6);
#line 232 "dsbtrd.f"
	return 0;
#line 233 "dsbtrd.f"
    }

/*     Quick return if possible */

#line 237 "dsbtrd.f"
    if (*n == 0) {
#line 237 "dsbtrd.f"
	return 0;
#line 237 "dsbtrd.f"
    }

/*     Initialize Q to the unit matrix, if needed */

#line 242 "dsbtrd.f"
    if (initq) {
#line 242 "dsbtrd.f"
	dlaset_("Full", n, n, &c_b9, &c_b10, &q[q_offset], ldq, (ftnlen)4);
#line 242 "dsbtrd.f"
    }

/*     Wherever possible, plane rotations are generated and applied in */
/*     vector operations of length NR over the index set J1:J2:KD1. */

/*     The cosines and sines of the plane rotations are stored in the */
/*     arrays D and WORK. */

#line 251 "dsbtrd.f"
    inca = kd1 * *ldab;
/* Computing MIN */
#line 252 "dsbtrd.f"
    i__1 = *n - 1;
#line 252 "dsbtrd.f"
    kdn = min(i__1,*kd);
#line 253 "dsbtrd.f"
    if (upper) {

#line 255 "dsbtrd.f"
	if (*kd > 1) {

/*           Reduce to tridiagonal form, working with upper triangle */

#line 259 "dsbtrd.f"
	    nr = 0;
#line 260 "dsbtrd.f"
	    j1 = kdn + 2;
#line 261 "dsbtrd.f"
	    j2 = 1;

#line 263 "dsbtrd.f"
	    i__1 = *n - 2;
#line 263 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th row of matrix to tridiagonal form */

#line 267 "dsbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 268 "dsbtrd.f"
		    j1 += kdn;
#line 269 "dsbtrd.f"
		    j2 += kdn;

#line 271 "dsbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 276 "dsbtrd.f"
			dlargv_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    DLARTV or DROT is used */

#line 285 "dsbtrd.f"
			if (nr >= (*kd << 1) - 1) {
#line 286 "dsbtrd.f"
			    i__2 = *kd - 1;
#line 286 "dsbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 287 "dsbtrd.f"
				dlartv_(&nr, &ab[l + 1 + (j1 - 1) * ab_dim1], 
					&inca, &ab[l + j1 * ab_dim1], &inca, &
					d__[j1], &work[j1], &kd1);
#line 290 "dsbtrd.f"
/* L10: */
#line 290 "dsbtrd.f"
			    }

#line 292 "dsbtrd.f"
			} else {
#line 293 "dsbtrd.f"
			    jend = j1 + (nr - 1) * kd1;
#line 294 "dsbtrd.f"
			    i__2 = jend;
#line 294 "dsbtrd.f"
			    i__3 = kd1;
#line 294 "dsbtrd.f"
			    for (jinc = j1; i__3 < 0 ? jinc >= i__2 : jinc <= 
				    i__2; jinc += i__3) {
#line 295 "dsbtrd.f"
				drot_(&kdm1, &ab[(jinc - 1) * ab_dim1 + 2], &
					c__1, &ab[jinc * ab_dim1 + 1], &c__1, 
					&d__[jinc], &work[jinc]);
#line 298 "dsbtrd.f"
/* L20: */
#line 298 "dsbtrd.f"
			    }
#line 299 "dsbtrd.f"
			}
#line 300 "dsbtrd.f"
		    }


#line 303 "dsbtrd.f"
		    if (k > 2) {
#line 304 "dsbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i,i+k-1) */
/*                       within the band */

#line 309 "dsbtrd.f"
			    dlartg_(&ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1]
				    , &ab[*kd - k + 2 + (i__ + k - 1) * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 312 "dsbtrd.f"
			    ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1] = temp;

/*                       apply rotation from the right */

#line 316 "dsbtrd.f"
			    i__3 = k - 3;
#line 316 "dsbtrd.f"
			    drot_(&i__3, &ab[*kd - k + 4 + (i__ + k - 2) * 
				    ab_dim1], &c__1, &ab[*kd - k + 3 + (i__ + 
				    k - 1) * ab_dim1], &c__1, &d__[i__ + k - 
				    1], &work[i__ + k - 1]);
#line 319 "dsbtrd.f"
			}
#line 320 "dsbtrd.f"
			++nr;
#line 321 "dsbtrd.f"
			j1 = j1 - kdn - 1;
#line 322 "dsbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 327 "dsbtrd.f"
		    if (nr > 0) {
#line 327 "dsbtrd.f"
			dlar2v_(&nr, &ab[kd1 + (j1 - 1) * ab_dim1], &ab[kd1 + 
				j1 * ab_dim1], &ab[*kd + j1 * ab_dim1], &inca,
				 &d__[j1], &work[j1], &kd1);
#line 327 "dsbtrd.f"
		    }

/*                 apply plane rotations from the left */

#line 334 "dsbtrd.f"
		    if (nr > 0) {
#line 335 "dsbtrd.f"
			if ((*kd << 1) - 1 < nr) {

/*                    Dependent on the the number of diagonals either */
/*                    DLARTV or DROT is used */

#line 340 "dsbtrd.f"
			    i__3 = *kd - 1;
#line 340 "dsbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 341 "dsbtrd.f"
				if (j2 + l > *n) {
#line 342 "dsbtrd.f"
				    nrt = nr - 1;
#line 343 "dsbtrd.f"
				} else {
#line 344 "dsbtrd.f"
				    nrt = nr;
#line 345 "dsbtrd.f"
				}
#line 346 "dsbtrd.f"
				if (nrt > 0) {
#line 346 "dsbtrd.f"
				    dlartv_(&nrt, &ab[*kd - l + (j1 + l) * 
					    ab_dim1], &inca, &ab[*kd - l + 1 
					    + (j1 + l) * ab_dim1], &inca, &
					    d__[j1], &work[j1], &kd1);
#line 346 "dsbtrd.f"
				}
#line 350 "dsbtrd.f"
/* L30: */
#line 350 "dsbtrd.f"
			    }
#line 351 "dsbtrd.f"
			} else {
#line 352 "dsbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 353 "dsbtrd.f"
			    if (j1end >= j1) {
#line 354 "dsbtrd.f"
				i__3 = j1end;
#line 354 "dsbtrd.f"
				i__2 = kd1;
#line 354 "dsbtrd.f"
				for (jin = j1; i__2 < 0 ? jin >= i__3 : jin <=
					 i__3; jin += i__2) {
#line 355 "dsbtrd.f"
				    i__4 = *kd - 1;
#line 355 "dsbtrd.f"
				    drot_(&i__4, &ab[*kd - 1 + (jin + 1) * 
					    ab_dim1], &incx, &ab[*kd + (jin + 
					    1) * ab_dim1], &incx, &d__[jin], &
					    work[jin]);
#line 358 "dsbtrd.f"
/* L40: */
#line 358 "dsbtrd.f"
				}
#line 359 "dsbtrd.f"
			    }
/* Computing MIN */
#line 360 "dsbtrd.f"
			    i__2 = kdm1, i__3 = *n - j2;
#line 360 "dsbtrd.f"
			    lend = min(i__2,i__3);
#line 361 "dsbtrd.f"
			    last = j1end + kd1;
#line 362 "dsbtrd.f"
			    if (lend > 0) {
#line 362 "dsbtrd.f"
				drot_(&lend, &ab[*kd - 1 + (last + 1) * 
					ab_dim1], &incx, &ab[*kd + (last + 1) 
					* ab_dim1], &incx, &d__[last], &work[
					last]);
#line 362 "dsbtrd.f"
			    }
#line 366 "dsbtrd.f"
			}
#line 367 "dsbtrd.f"
		    }

#line 369 "dsbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 373 "dsbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 378 "dsbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 379 "dsbtrd.f"
			    i__2 = 0, i__3 = k - 3;
#line 379 "dsbtrd.f"
			    i2 = max(i__2,i__3);
#line 380 "dsbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 381 "dsbtrd.f"
			    if (k == 2) {
#line 381 "dsbtrd.f"
				iqaend += *kd;
#line 381 "dsbtrd.f"
			    }
#line 383 "dsbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 384 "dsbtrd.f"
			    i__2 = j2;
#line 384 "dsbtrd.f"
			    i__3 = kd1;
#line 384 "dsbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 385 "dsbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 386 "dsbtrd.f"
				++i2;
/* Computing MAX */
#line 387 "dsbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 387 "dsbtrd.f"
				iqb = max(i__4,i__5);
#line 388 "dsbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 389 "dsbtrd.f"
				i__4 = iqaend + *kd;
#line 389 "dsbtrd.f"
				iqaend = min(i__4,iqend);
#line 390 "dsbtrd.f"
				drot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&work[j]);
#line 392 "dsbtrd.f"
/* L50: */
#line 392 "dsbtrd.f"
			    }
#line 393 "dsbtrd.f"
			} else {

#line 395 "dsbtrd.f"
			    i__3 = j2;
#line 395 "dsbtrd.f"
			    i__2 = kd1;
#line 395 "dsbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 396 "dsbtrd.f"
				drot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					work[j]);
#line 398 "dsbtrd.f"
/* L60: */
#line 398 "dsbtrd.f"
			    }
#line 399 "dsbtrd.f"
			}

#line 401 "dsbtrd.f"
		    }

#line 403 "dsbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 407 "dsbtrd.f"
			--nr;
#line 408 "dsbtrd.f"
			j2 = j2 - kdn - 1;
#line 409 "dsbtrd.f"
		    }

#line 411 "dsbtrd.f"
		    i__2 = j2;
#line 411 "dsbtrd.f"
		    i__3 = kd1;
#line 411 "dsbtrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) 
			    {

/*                    create nonzero element a(j-1,j+kd) outside the band */
/*                    and store it in WORK */

#line 416 "dsbtrd.f"
			work[j + *kd] = work[j] * ab[(j + *kd) * ab_dim1 + 1];
#line 417 "dsbtrd.f"
			ab[(j + *kd) * ab_dim1 + 1] = d__[j] * ab[(j + *kd) * 
				ab_dim1 + 1];
#line 418 "dsbtrd.f"
/* L70: */
#line 418 "dsbtrd.f"
		    }
#line 419 "dsbtrd.f"
/* L80: */
#line 419 "dsbtrd.f"
		}
#line 420 "dsbtrd.f"
/* L90: */
#line 420 "dsbtrd.f"
	    }
#line 421 "dsbtrd.f"
	}

#line 423 "dsbtrd.f"
	if (*kd > 0) {

/*           copy off-diagonal elements to E */

#line 427 "dsbtrd.f"
	    i__1 = *n - 1;
#line 427 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "dsbtrd.f"
		e[i__] = ab[*kd + (i__ + 1) * ab_dim1];
#line 429 "dsbtrd.f"
/* L100: */
#line 429 "dsbtrd.f"
	    }
#line 430 "dsbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 434 "dsbtrd.f"
	    i__1 = *n - 1;
#line 434 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "dsbtrd.f"
		e[i__] = 0.;
#line 436 "dsbtrd.f"
/* L110: */
#line 436 "dsbtrd.f"
	    }
#line 437 "dsbtrd.f"
	}

/*        copy diagonal elements to D */

#line 441 "dsbtrd.f"
	i__1 = *n;
#line 441 "dsbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 442 "dsbtrd.f"
	    d__[i__] = ab[kd1 + i__ * ab_dim1];
#line 443 "dsbtrd.f"
/* L120: */
#line 443 "dsbtrd.f"
	}

#line 445 "dsbtrd.f"
    } else {

#line 447 "dsbtrd.f"
	if (*kd > 1) {

/*           Reduce to tridiagonal form, working with lower triangle */

#line 451 "dsbtrd.f"
	    nr = 0;
#line 452 "dsbtrd.f"
	    j1 = kdn + 2;
#line 453 "dsbtrd.f"
	    j2 = 1;

#line 455 "dsbtrd.f"
	    i__1 = *n - 2;
#line 455 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th column of matrix to tridiagonal form */

#line 459 "dsbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 460 "dsbtrd.f"
		    j1 += kdn;
#line 461 "dsbtrd.f"
		    j2 += kdn;

#line 463 "dsbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 468 "dsbtrd.f"
			dlargv_(&nr, &ab[kd1 + (j1 - kd1) * ab_dim1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply plane rotations from one side */


/*                    Dependent on the the number of diagonals either */
/*                    DLARTV or DROT is used */

#line 477 "dsbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 478 "dsbtrd.f"
			    i__3 = *kd - 1;
#line 478 "dsbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 479 "dsbtrd.f"
				dlartv_(&nr, &ab[kd1 - l + (j1 - kd1 + l) * 
					ab_dim1], &inca, &ab[kd1 - l + 1 + (
					j1 - kd1 + l) * ab_dim1], &inca, &d__[
					j1], &work[j1], &kd1);
#line 482 "dsbtrd.f"
/* L130: */
#line 482 "dsbtrd.f"
			    }
#line 483 "dsbtrd.f"
			} else {
#line 484 "dsbtrd.f"
			    jend = j1 + kd1 * (nr - 1);
#line 485 "dsbtrd.f"
			    i__3 = jend;
#line 485 "dsbtrd.f"
			    i__2 = kd1;
#line 485 "dsbtrd.f"
			    for (jinc = j1; i__2 < 0 ? jinc >= i__3 : jinc <= 
				    i__3; jinc += i__2) {
#line 486 "dsbtrd.f"
				drot_(&kdm1, &ab[*kd + (jinc - *kd) * ab_dim1]
					, &incx, &ab[kd1 + (jinc - *kd) * 
					ab_dim1], &incx, &d__[jinc], &work[
					jinc]);
#line 489 "dsbtrd.f"
/* L140: */
#line 489 "dsbtrd.f"
			    }
#line 490 "dsbtrd.f"
			}

#line 492 "dsbtrd.f"
		    }

#line 494 "dsbtrd.f"
		    if (k > 2) {
#line 495 "dsbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i+k-1,i) */
/*                       within the band */

#line 500 "dsbtrd.f"
			    dlartg_(&ab[k - 1 + i__ * ab_dim1], &ab[k + i__ * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 502 "dsbtrd.f"
			    ab[k - 1 + i__ * ab_dim1] = temp;

/*                       apply rotation from the left */

#line 506 "dsbtrd.f"
			    i__2 = k - 3;
#line 506 "dsbtrd.f"
			    i__3 = *ldab - 1;
#line 506 "dsbtrd.f"
			    i__4 = *ldab - 1;
#line 506 "dsbtrd.f"
			    drot_(&i__2, &ab[k - 2 + (i__ + 1) * ab_dim1], &
				    i__3, &ab[k - 1 + (i__ + 1) * ab_dim1], &
				    i__4, &d__[i__ + k - 1], &work[i__ + k - 
				    1]);
#line 509 "dsbtrd.f"
			}
#line 510 "dsbtrd.f"
			++nr;
#line 511 "dsbtrd.f"
			j1 = j1 - kdn - 1;
#line 512 "dsbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 517 "dsbtrd.f"
		    if (nr > 0) {
#line 517 "dsbtrd.f"
			dlar2v_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &ab[j1 * 
				ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 2], &
				inca, &d__[j1], &work[j1], &kd1);
#line 517 "dsbtrd.f"
		    }

/*                 apply plane rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    DLARTV or DROT is used */

#line 528 "dsbtrd.f"
		    if (nr > 0) {
#line 529 "dsbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 530 "dsbtrd.f"
			    i__2 = *kd - 1;
#line 530 "dsbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 531 "dsbtrd.f"
				if (j2 + l > *n) {
#line 532 "dsbtrd.f"
				    nrt = nr - 1;
#line 533 "dsbtrd.f"
				} else {
#line 534 "dsbtrd.f"
				    nrt = nr;
#line 535 "dsbtrd.f"
				}
#line 536 "dsbtrd.f"
				if (nrt > 0) {
#line 536 "dsbtrd.f"
				    dlartv_(&nrt, &ab[l + 2 + (j1 - 1) * 
					    ab_dim1], &inca, &ab[l + 1 + j1 * 
					    ab_dim1], &inca, &d__[j1], &work[
					    j1], &kd1);
#line 536 "dsbtrd.f"
				}
#line 540 "dsbtrd.f"
/* L150: */
#line 540 "dsbtrd.f"
			    }
#line 541 "dsbtrd.f"
			} else {
#line 542 "dsbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 543 "dsbtrd.f"
			    if (j1end >= j1) {
#line 544 "dsbtrd.f"
				i__2 = j1end;
#line 544 "dsbtrd.f"
				i__3 = kd1;
#line 544 "dsbtrd.f"
				for (j1inc = j1; i__3 < 0 ? j1inc >= i__2 : 
					j1inc <= i__2; j1inc += i__3) {
#line 545 "dsbtrd.f"
				    drot_(&kdm1, &ab[(j1inc - 1) * ab_dim1 + 
					    3], &c__1, &ab[j1inc * ab_dim1 + 
					    2], &c__1, &d__[j1inc], &work[
					    j1inc]);
#line 548 "dsbtrd.f"
/* L160: */
#line 548 "dsbtrd.f"
				}
#line 549 "dsbtrd.f"
			    }
/* Computing MIN */
#line 550 "dsbtrd.f"
			    i__3 = kdm1, i__2 = *n - j2;
#line 550 "dsbtrd.f"
			    lend = min(i__3,i__2);
#line 551 "dsbtrd.f"
			    last = j1end + kd1;
#line 552 "dsbtrd.f"
			    if (lend > 0) {
#line 552 "dsbtrd.f"
				drot_(&lend, &ab[(last - 1) * ab_dim1 + 3], &
					c__1, &ab[last * ab_dim1 + 2], &c__1, 
					&d__[last], &work[last]);
#line 552 "dsbtrd.f"
			    }
#line 556 "dsbtrd.f"
			}
#line 557 "dsbtrd.f"
		    }



#line 561 "dsbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 565 "dsbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 570 "dsbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 571 "dsbtrd.f"
			    i__3 = 0, i__2 = k - 3;
#line 571 "dsbtrd.f"
			    i2 = max(i__3,i__2);
#line 572 "dsbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 573 "dsbtrd.f"
			    if (k == 2) {
#line 573 "dsbtrd.f"
				iqaend += *kd;
#line 573 "dsbtrd.f"
			    }
#line 575 "dsbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 576 "dsbtrd.f"
			    i__3 = j2;
#line 576 "dsbtrd.f"
			    i__2 = kd1;
#line 576 "dsbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 577 "dsbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 578 "dsbtrd.f"
				++i2;
/* Computing MAX */
#line 579 "dsbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 579 "dsbtrd.f"
				iqb = max(i__4,i__5);
#line 580 "dsbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 581 "dsbtrd.f"
				i__4 = iqaend + *kd;
#line 581 "dsbtrd.f"
				iqaend = min(i__4,iqend);
#line 582 "dsbtrd.f"
				drot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&work[j]);
#line 584 "dsbtrd.f"
/* L170: */
#line 584 "dsbtrd.f"
			    }
#line 585 "dsbtrd.f"
			} else {

#line 587 "dsbtrd.f"
			    i__2 = j2;
#line 587 "dsbtrd.f"
			    i__3 = kd1;
#line 587 "dsbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 588 "dsbtrd.f"
				drot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					work[j]);
#line 590 "dsbtrd.f"
/* L180: */
#line 590 "dsbtrd.f"
			    }
#line 591 "dsbtrd.f"
			}
#line 592 "dsbtrd.f"
		    }

#line 594 "dsbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 598 "dsbtrd.f"
			--nr;
#line 599 "dsbtrd.f"
			j2 = j2 - kdn - 1;
#line 600 "dsbtrd.f"
		    }

#line 602 "dsbtrd.f"
		    i__3 = j2;
#line 602 "dsbtrd.f"
		    i__2 = kd1;
#line 602 "dsbtrd.f"
		    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) 
			    {

/*                    create nonzero element a(j+kd,j-1) outside the */
/*                    band and store it in WORK */

#line 607 "dsbtrd.f"
			work[j + *kd] = work[j] * ab[kd1 + j * ab_dim1];
#line 608 "dsbtrd.f"
			ab[kd1 + j * ab_dim1] = d__[j] * ab[kd1 + j * ab_dim1]
				;
#line 609 "dsbtrd.f"
/* L190: */
#line 609 "dsbtrd.f"
		    }
#line 610 "dsbtrd.f"
/* L200: */
#line 610 "dsbtrd.f"
		}
#line 611 "dsbtrd.f"
/* L210: */
#line 611 "dsbtrd.f"
	    }
#line 612 "dsbtrd.f"
	}

#line 614 "dsbtrd.f"
	if (*kd > 0) {

/*           copy off-diagonal elements to E */

#line 618 "dsbtrd.f"
	    i__1 = *n - 1;
#line 618 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 619 "dsbtrd.f"
		e[i__] = ab[i__ * ab_dim1 + 2];
#line 620 "dsbtrd.f"
/* L220: */
#line 620 "dsbtrd.f"
	    }
#line 621 "dsbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 625 "dsbtrd.f"
	    i__1 = *n - 1;
#line 625 "dsbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 626 "dsbtrd.f"
		e[i__] = 0.;
#line 627 "dsbtrd.f"
/* L230: */
#line 627 "dsbtrd.f"
	    }
#line 628 "dsbtrd.f"
	}

/*        copy diagonal elements to D */

#line 632 "dsbtrd.f"
	i__1 = *n;
#line 632 "dsbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 633 "dsbtrd.f"
	    d__[i__] = ab[i__ * ab_dim1 + 1];
#line 634 "dsbtrd.f"
/* L240: */
#line 634 "dsbtrd.f"
	}
#line 635 "dsbtrd.f"
    }

#line 637 "dsbtrd.f"
    return 0;

/*     End of DSBTRD */

} /* dsbtrd_ */

