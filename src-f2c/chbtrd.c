#line 1 "chbtrd.f"
/* chbtrd.f -- translated by f2c (version 20100827).
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

#line 1 "chbtrd.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CHBTRD */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CHBTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbtrd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbtrd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbtrd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, */
/*                          WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO, VECT */
/*       INTEGER            INFO, KD, LDAB, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), E( * ) */
/*       COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBTRD reduces a complex Hermitian band matrix A to real symmetric */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          D is REAL array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ,N) */
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
/* >          WORK is COMPLEX array, dimension (N) */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by Linda Kaufman, Bell Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int chbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
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
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer j1end, j1inc;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer iqend;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical initq, wantq, upper;
    extern /* Subroutine */ int clar2v_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublecomplex *, integer *), clacgv_(integer *, doublecomplex *, 
	    integer *);
    static integer iqaend;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), clartg_(doublecomplex *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *), xerbla_(char *, integer *, 
	    ftnlen), clargv_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *), clartv_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, doublecomplex *, integer *);


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

#line 212 "chbtrd.f"
    /* Parameter adjustments */
#line 212 "chbtrd.f"
    ab_dim1 = *ldab;
#line 212 "chbtrd.f"
    ab_offset = 1 + ab_dim1;
#line 212 "chbtrd.f"
    ab -= ab_offset;
#line 212 "chbtrd.f"
    --d__;
#line 212 "chbtrd.f"
    --e;
#line 212 "chbtrd.f"
    q_dim1 = *ldq;
#line 212 "chbtrd.f"
    q_offset = 1 + q_dim1;
#line 212 "chbtrd.f"
    q -= q_offset;
#line 212 "chbtrd.f"
    --work;
#line 212 "chbtrd.f"

#line 212 "chbtrd.f"
    /* Function Body */
#line 212 "chbtrd.f"
    initq = lsame_(vect, "V", (ftnlen)1, (ftnlen)1);
#line 213 "chbtrd.f"
    wantq = initq || lsame_(vect, "U", (ftnlen)1, (ftnlen)1);
#line 214 "chbtrd.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 215 "chbtrd.f"
    kd1 = *kd + 1;
#line 216 "chbtrd.f"
    kdm1 = *kd - 1;
#line 217 "chbtrd.f"
    incx = *ldab - 1;
#line 218 "chbtrd.f"
    iqend = 1;

#line 220 "chbtrd.f"
    *info = 0;
#line 221 "chbtrd.f"
    if (! wantq && ! lsame_(vect, "N", (ftnlen)1, (ftnlen)1)) {
#line 222 "chbtrd.f"
	*info = -1;
#line 223 "chbtrd.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 224 "chbtrd.f"
	*info = -2;
#line 225 "chbtrd.f"
    } else if (*n < 0) {
#line 226 "chbtrd.f"
	*info = -3;
#line 227 "chbtrd.f"
    } else if (*kd < 0) {
#line 228 "chbtrd.f"
	*info = -4;
#line 229 "chbtrd.f"
    } else if (*ldab < kd1) {
#line 230 "chbtrd.f"
	*info = -6;
#line 231 "chbtrd.f"
    } else if (*ldq < max(1,*n) && wantq) {
#line 232 "chbtrd.f"
	*info = -10;
#line 233 "chbtrd.f"
    }
#line 234 "chbtrd.f"
    if (*info != 0) {
#line 235 "chbtrd.f"
	i__1 = -(*info);
#line 235 "chbtrd.f"
	xerbla_("CHBTRD", &i__1, (ftnlen)6);
#line 236 "chbtrd.f"
	return 0;
#line 237 "chbtrd.f"
    }

/*     Quick return if possible */

#line 241 "chbtrd.f"
    if (*n == 0) {
#line 241 "chbtrd.f"
	return 0;
#line 241 "chbtrd.f"
    }

/*     Initialize Q to the unit matrix, if needed */

#line 246 "chbtrd.f"
    if (initq) {
#line 246 "chbtrd.f"
	claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 246 "chbtrd.f"
    }

/*     Wherever possible, plane rotations are generated and applied in */
/*     vector operations of length NR over the index set J1:J2:KD1. */

/*     The real cosines and complex sines of the plane rotations are */
/*     stored in the arrays D and WORK. */

#line 255 "chbtrd.f"
    inca = kd1 * *ldab;
/* Computing MIN */
#line 256 "chbtrd.f"
    i__1 = *n - 1;
#line 256 "chbtrd.f"
    kdn = min(i__1,*kd);
#line 257 "chbtrd.f"
    if (upper) {

#line 259 "chbtrd.f"
	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working with */
/*           the upper triangle */

#line 264 "chbtrd.f"
	    nr = 0;
#line 265 "chbtrd.f"
	    j1 = kdn + 2;
#line 266 "chbtrd.f"
	    j2 = 1;

#line 268 "chbtrd.f"
	    i__1 = kd1 + ab_dim1;
#line 268 "chbtrd.f"
	    i__2 = kd1 + ab_dim1;
#line 268 "chbtrd.f"
	    d__1 = ab[i__2].r;
#line 268 "chbtrd.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 269 "chbtrd.f"
	    i__1 = *n - 2;
#line 269 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th row of matrix to tridiagonal form */

#line 273 "chbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 274 "chbtrd.f"
		    j1 += kdn;
#line 275 "chbtrd.f"
		    j2 += kdn;

#line 277 "chbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 282 "chbtrd.f"
			clargv_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    CLARTV or CROT is used */

#line 291 "chbtrd.f"
			if (nr >= (*kd << 1) - 1) {
#line 292 "chbtrd.f"
			    i__2 = *kd - 1;
#line 292 "chbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 293 "chbtrd.f"
				clartv_(&nr, &ab[l + 1 + (j1 - 1) * ab_dim1], 
					&inca, &ab[l + j1 * ab_dim1], &inca, &
					d__[j1], &work[j1], &kd1);
#line 296 "chbtrd.f"
/* L10: */
#line 296 "chbtrd.f"
			    }

#line 298 "chbtrd.f"
			} else {
#line 299 "chbtrd.f"
			    jend = j1 + (nr - 1) * kd1;
#line 300 "chbtrd.f"
			    i__2 = jend;
#line 300 "chbtrd.f"
			    i__3 = kd1;
#line 300 "chbtrd.f"
			    for (jinc = j1; i__3 < 0 ? jinc >= i__2 : jinc <= 
				    i__2; jinc += i__3) {
#line 301 "chbtrd.f"
				crot_(&kdm1, &ab[(jinc - 1) * ab_dim1 + 2], &
					c__1, &ab[jinc * ab_dim1 + 1], &c__1, 
					&d__[jinc], &work[jinc]);
#line 304 "chbtrd.f"
/* L20: */
#line 304 "chbtrd.f"
			    }
#line 305 "chbtrd.f"
			}
#line 306 "chbtrd.f"
		    }


#line 309 "chbtrd.f"
		    if (k > 2) {
#line 310 "chbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i,i+k-1) */
/*                       within the band */

#line 315 "chbtrd.f"
			    clartg_(&ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1]
				    , &ab[*kd - k + 2 + (i__ + k - 1) * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 318 "chbtrd.f"
			    i__3 = *kd - k + 3 + (i__ + k - 2) * ab_dim1;
#line 318 "chbtrd.f"
			    ab[i__3].r = temp.r, ab[i__3].i = temp.i;

/*                       apply rotation from the right */

#line 322 "chbtrd.f"
			    i__3 = k - 3;
#line 322 "chbtrd.f"
			    crot_(&i__3, &ab[*kd - k + 4 + (i__ + k - 2) * 
				    ab_dim1], &c__1, &ab[*kd - k + 3 + (i__ + 
				    k - 1) * ab_dim1], &c__1, &d__[i__ + k - 
				    1], &work[i__ + k - 1]);
#line 325 "chbtrd.f"
			}
#line 326 "chbtrd.f"
			++nr;
#line 327 "chbtrd.f"
			j1 = j1 - kdn - 1;
#line 328 "chbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 333 "chbtrd.f"
		    if (nr > 0) {
#line 333 "chbtrd.f"
			clar2v_(&nr, &ab[kd1 + (j1 - 1) * ab_dim1], &ab[kd1 + 
				j1 * ab_dim1], &ab[*kd + j1 * ab_dim1], &inca,
				 &d__[j1], &work[j1], &kd1);
#line 333 "chbtrd.f"
		    }

/*                 apply plane rotations from the left */

#line 340 "chbtrd.f"
		    if (nr > 0) {
#line 341 "chbtrd.f"
			clacgv_(&nr, &work[j1], &kd1);
#line 342 "chbtrd.f"
			if ((*kd << 1) - 1 < nr) {

/*                    Dependent on the the number of diagonals either */
/*                    CLARTV or CROT is used */

#line 347 "chbtrd.f"
			    i__3 = *kd - 1;
#line 347 "chbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 348 "chbtrd.f"
				if (j2 + l > *n) {
#line 349 "chbtrd.f"
				    nrt = nr - 1;
#line 350 "chbtrd.f"
				} else {
#line 351 "chbtrd.f"
				    nrt = nr;
#line 352 "chbtrd.f"
				}
#line 353 "chbtrd.f"
				if (nrt > 0) {
#line 353 "chbtrd.f"
				    clartv_(&nrt, &ab[*kd - l + (j1 + l) * 
					    ab_dim1], &inca, &ab[*kd - l + 1 
					    + (j1 + l) * ab_dim1], &inca, &
					    d__[j1], &work[j1], &kd1);
#line 353 "chbtrd.f"
				}
#line 357 "chbtrd.f"
/* L30: */
#line 357 "chbtrd.f"
			    }
#line 358 "chbtrd.f"
			} else {
#line 359 "chbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 360 "chbtrd.f"
			    if (j1end >= j1) {
#line 361 "chbtrd.f"
				i__3 = j1end;
#line 361 "chbtrd.f"
				i__2 = kd1;
#line 361 "chbtrd.f"
				for (jin = j1; i__2 < 0 ? jin >= i__3 : jin <=
					 i__3; jin += i__2) {
#line 362 "chbtrd.f"
				    i__4 = *kd - 1;
#line 362 "chbtrd.f"
				    crot_(&i__4, &ab[*kd - 1 + (jin + 1) * 
					    ab_dim1], &incx, &ab[*kd + (jin + 
					    1) * ab_dim1], &incx, &d__[jin], &
					    work[jin]);
#line 365 "chbtrd.f"
/* L40: */
#line 365 "chbtrd.f"
				}
#line 366 "chbtrd.f"
			    }
/* Computing MIN */
#line 367 "chbtrd.f"
			    i__2 = kdm1, i__3 = *n - j2;
#line 367 "chbtrd.f"
			    lend = min(i__2,i__3);
#line 368 "chbtrd.f"
			    last = j1end + kd1;
#line 369 "chbtrd.f"
			    if (lend > 0) {
#line 369 "chbtrd.f"
				crot_(&lend, &ab[*kd - 1 + (last + 1) * 
					ab_dim1], &incx, &ab[*kd + (last + 1) 
					* ab_dim1], &incx, &d__[last], &work[
					last]);
#line 369 "chbtrd.f"
			    }
#line 373 "chbtrd.f"
			}
#line 374 "chbtrd.f"
		    }

#line 376 "chbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 380 "chbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 385 "chbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 386 "chbtrd.f"
			    i__2 = 0, i__3 = k - 3;
#line 386 "chbtrd.f"
			    i2 = max(i__2,i__3);
#line 387 "chbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 388 "chbtrd.f"
			    if (k == 2) {
#line 388 "chbtrd.f"
				iqaend += *kd;
#line 388 "chbtrd.f"
			    }
#line 390 "chbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 391 "chbtrd.f"
			    i__2 = j2;
#line 391 "chbtrd.f"
			    i__3 = kd1;
#line 391 "chbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 392 "chbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 393 "chbtrd.f"
				++i2;
/* Computing MAX */
#line 394 "chbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 394 "chbtrd.f"
				iqb = max(i__4,i__5);
#line 395 "chbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 396 "chbtrd.f"
				i__4 = iqaend + *kd;
#line 396 "chbtrd.f"
				iqaend = min(i__4,iqend);
#line 397 "chbtrd.f"
				d_cnjg(&z__1, &work[j]);
#line 397 "chbtrd.f"
				crot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&z__1);
#line 399 "chbtrd.f"
/* L50: */
#line 399 "chbtrd.f"
			    }
#line 400 "chbtrd.f"
			} else {

#line 402 "chbtrd.f"
			    i__3 = j2;
#line 402 "chbtrd.f"
			    i__2 = kd1;
#line 402 "chbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 403 "chbtrd.f"
				d_cnjg(&z__1, &work[j]);
#line 403 "chbtrd.f"
				crot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					z__1);
#line 405 "chbtrd.f"
/* L60: */
#line 405 "chbtrd.f"
			    }
#line 406 "chbtrd.f"
			}

#line 408 "chbtrd.f"
		    }

#line 410 "chbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 414 "chbtrd.f"
			--nr;
#line 415 "chbtrd.f"
			j2 = j2 - kdn - 1;
#line 416 "chbtrd.f"
		    }

#line 418 "chbtrd.f"
		    i__2 = j2;
#line 418 "chbtrd.f"
		    i__3 = kd1;
#line 418 "chbtrd.f"
		    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) 
			    {

/*                    create nonzero element a(j-1,j+kd) outside the band */
/*                    and store it in WORK */

#line 423 "chbtrd.f"
			i__4 = j + *kd;
#line 423 "chbtrd.f"
			i__5 = j;
#line 423 "chbtrd.f"
			i__6 = (j + *kd) * ab_dim1 + 1;
#line 423 "chbtrd.f"
			z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * 
				ab[i__6].i, z__1.i = work[i__5].r * ab[i__6]
				.i + work[i__5].i * ab[i__6].r;
#line 423 "chbtrd.f"
			work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 424 "chbtrd.f"
			i__4 = (j + *kd) * ab_dim1 + 1;
#line 424 "chbtrd.f"
			i__5 = j;
#line 424 "chbtrd.f"
			i__6 = (j + *kd) * ab_dim1 + 1;
#line 424 "chbtrd.f"
			z__1.r = d__[i__5] * ab[i__6].r, z__1.i = d__[i__5] * 
				ab[i__6].i;
#line 424 "chbtrd.f"
			ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 425 "chbtrd.f"
/* L70: */
#line 425 "chbtrd.f"
		    }
#line 426 "chbtrd.f"
/* L80: */
#line 426 "chbtrd.f"
		}
#line 427 "chbtrd.f"
/* L90: */
#line 427 "chbtrd.f"
	    }
#line 428 "chbtrd.f"
	}

#line 430 "chbtrd.f"
	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

#line 434 "chbtrd.f"
	    i__1 = *n - 1;
#line 434 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "chbtrd.f"
		i__3 = *kd + (i__ + 1) * ab_dim1;
#line 435 "chbtrd.f"
		t.r = ab[i__3].r, t.i = ab[i__3].i;
#line 436 "chbtrd.f"
		abst = z_abs(&t);
#line 437 "chbtrd.f"
		i__3 = *kd + (i__ + 1) * ab_dim1;
#line 437 "chbtrd.f"
		ab[i__3].r = abst, ab[i__3].i = 0.;
#line 438 "chbtrd.f"
		e[i__] = abst;
#line 439 "chbtrd.f"
		if (abst != 0.) {
#line 440 "chbtrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 440 "chbtrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 441 "chbtrd.f"
		} else {
#line 442 "chbtrd.f"
		    t.r = 1., t.i = 0.;
#line 443 "chbtrd.f"
		}
#line 444 "chbtrd.f"
		if (i__ < *n - 1) {
#line 444 "chbtrd.f"
		    i__3 = *kd + (i__ + 2) * ab_dim1;
#line 444 "chbtrd.f"
		    i__2 = *kd + (i__ + 2) * ab_dim1;
#line 444 "chbtrd.f"
		    z__1.r = ab[i__2].r * t.r - ab[i__2].i * t.i, z__1.i = ab[
			    i__2].r * t.i + ab[i__2].i * t.r;
#line 444 "chbtrd.f"
		    ab[i__3].r = z__1.r, ab[i__3].i = z__1.i;
#line 444 "chbtrd.f"
		}
#line 446 "chbtrd.f"
		if (wantq) {
#line 447 "chbtrd.f"
		    d_cnjg(&z__1, &t);
#line 447 "chbtrd.f"
		    cscal_(n, &z__1, &q[(i__ + 1) * q_dim1 + 1], &c__1);
#line 448 "chbtrd.f"
		}
#line 449 "chbtrd.f"
/* L100: */
#line 449 "chbtrd.f"
	    }
#line 450 "chbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 454 "chbtrd.f"
	    i__1 = *n - 1;
#line 454 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 455 "chbtrd.f"
		e[i__] = 0.;
#line 456 "chbtrd.f"
/* L110: */
#line 456 "chbtrd.f"
	    }
#line 457 "chbtrd.f"
	}

/*        copy diagonal elements to D */

#line 461 "chbtrd.f"
	i__1 = *n;
#line 461 "chbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 462 "chbtrd.f"
	    i__3 = i__;
#line 462 "chbtrd.f"
	    i__2 = kd1 + i__ * ab_dim1;
#line 462 "chbtrd.f"
	    d__[i__3] = ab[i__2].r;
#line 463 "chbtrd.f"
/* L120: */
#line 463 "chbtrd.f"
	}

#line 465 "chbtrd.f"
    } else {

#line 467 "chbtrd.f"
	if (*kd > 1) {

/*           Reduce to complex Hermitian tridiagonal form, working with */
/*           the lower triangle */

#line 472 "chbtrd.f"
	    nr = 0;
#line 473 "chbtrd.f"
	    j1 = kdn + 2;
#line 474 "chbtrd.f"
	    j2 = 1;

#line 476 "chbtrd.f"
	    i__1 = ab_dim1 + 1;
#line 476 "chbtrd.f"
	    i__3 = ab_dim1 + 1;
#line 476 "chbtrd.f"
	    d__1 = ab[i__3].r;
#line 476 "chbtrd.f"
	    ab[i__1].r = d__1, ab[i__1].i = 0.;
#line 477 "chbtrd.f"
	    i__1 = *n - 2;
#line 477 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Reduce i-th column of matrix to tridiagonal form */

#line 481 "chbtrd.f"
		for (k = kdn + 1; k >= 2; --k) {
#line 482 "chbtrd.f"
		    j1 += kdn;
#line 483 "chbtrd.f"
		    j2 += kdn;

#line 485 "chbtrd.f"
		    if (nr > 0) {

/*                    generate plane rotations to annihilate nonzero */
/*                    elements which have been created outside the band */

#line 490 "chbtrd.f"
			clargv_(&nr, &ab[kd1 + (j1 - kd1) * ab_dim1], &inca, &
				work[j1], &kd1, &d__[j1], &kd1);

/*                    apply plane rotations from one side */


/*                    Dependent on the the number of diagonals either */
/*                    CLARTV or CROT is used */

#line 499 "chbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 500 "chbtrd.f"
			    i__3 = *kd - 1;
#line 500 "chbtrd.f"
			    for (l = 1; l <= i__3; ++l) {
#line 501 "chbtrd.f"
				clartv_(&nr, &ab[kd1 - l + (j1 - kd1 + l) * 
					ab_dim1], &inca, &ab[kd1 - l + 1 + (
					j1 - kd1 + l) * ab_dim1], &inca, &d__[
					j1], &work[j1], &kd1);
#line 504 "chbtrd.f"
/* L130: */
#line 504 "chbtrd.f"
			    }
#line 505 "chbtrd.f"
			} else {
#line 506 "chbtrd.f"
			    jend = j1 + kd1 * (nr - 1);
#line 507 "chbtrd.f"
			    i__3 = jend;
#line 507 "chbtrd.f"
			    i__2 = kd1;
#line 507 "chbtrd.f"
			    for (jinc = j1; i__2 < 0 ? jinc >= i__3 : jinc <= 
				    i__3; jinc += i__2) {
#line 508 "chbtrd.f"
				crot_(&kdm1, &ab[*kd + (jinc - *kd) * ab_dim1]
					, &incx, &ab[kd1 + (jinc - *kd) * 
					ab_dim1], &incx, &d__[jinc], &work[
					jinc]);
#line 511 "chbtrd.f"
/* L140: */
#line 511 "chbtrd.f"
			    }
#line 512 "chbtrd.f"
			}

#line 514 "chbtrd.f"
		    }

#line 516 "chbtrd.f"
		    if (k > 2) {
#line 517 "chbtrd.f"
			if (k <= *n - i__ + 1) {

/*                       generate plane rotation to annihilate a(i+k-1,i) */
/*                       within the band */

#line 522 "chbtrd.f"
			    clartg_(&ab[k - 1 + i__ * ab_dim1], &ab[k + i__ * 
				    ab_dim1], &d__[i__ + k - 1], &work[i__ + 
				    k - 1], &temp);
#line 524 "chbtrd.f"
			    i__2 = k - 1 + i__ * ab_dim1;
#line 524 "chbtrd.f"
			    ab[i__2].r = temp.r, ab[i__2].i = temp.i;

/*                       apply rotation from the left */

#line 528 "chbtrd.f"
			    i__2 = k - 3;
#line 528 "chbtrd.f"
			    i__3 = *ldab - 1;
#line 528 "chbtrd.f"
			    i__4 = *ldab - 1;
#line 528 "chbtrd.f"
			    crot_(&i__2, &ab[k - 2 + (i__ + 1) * ab_dim1], &
				    i__3, &ab[k - 1 + (i__ + 1) * ab_dim1], &
				    i__4, &d__[i__ + k - 1], &work[i__ + k - 
				    1]);
#line 531 "chbtrd.f"
			}
#line 532 "chbtrd.f"
			++nr;
#line 533 "chbtrd.f"
			j1 = j1 - kdn - 1;
#line 534 "chbtrd.f"
		    }

/*                 apply plane rotations from both sides to diagonal */
/*                 blocks */

#line 539 "chbtrd.f"
		    if (nr > 0) {
#line 539 "chbtrd.f"
			clar2v_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &ab[j1 * 
				ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 2], &
				inca, &d__[j1], &work[j1], &kd1);
#line 539 "chbtrd.f"
		    }

/*                 apply plane rotations from the right */


/*                    Dependent on the the number of diagonals either */
/*                    CLARTV or CROT is used */

#line 550 "chbtrd.f"
		    if (nr > 0) {
#line 551 "chbtrd.f"
			clacgv_(&nr, &work[j1], &kd1);
#line 552 "chbtrd.f"
			if (nr > (*kd << 1) - 1) {
#line 553 "chbtrd.f"
			    i__2 = *kd - 1;
#line 553 "chbtrd.f"
			    for (l = 1; l <= i__2; ++l) {
#line 554 "chbtrd.f"
				if (j2 + l > *n) {
#line 555 "chbtrd.f"
				    nrt = nr - 1;
#line 556 "chbtrd.f"
				} else {
#line 557 "chbtrd.f"
				    nrt = nr;
#line 558 "chbtrd.f"
				}
#line 559 "chbtrd.f"
				if (nrt > 0) {
#line 559 "chbtrd.f"
				    clartv_(&nrt, &ab[l + 2 + (j1 - 1) * 
					    ab_dim1], &inca, &ab[l + 1 + j1 * 
					    ab_dim1], &inca, &d__[j1], &work[
					    j1], &kd1);
#line 559 "chbtrd.f"
				}
#line 563 "chbtrd.f"
/* L150: */
#line 563 "chbtrd.f"
			    }
#line 564 "chbtrd.f"
			} else {
#line 565 "chbtrd.f"
			    j1end = j1 + kd1 * (nr - 2);
#line 566 "chbtrd.f"
			    if (j1end >= j1) {
#line 567 "chbtrd.f"
				i__2 = j1end;
#line 567 "chbtrd.f"
				i__3 = kd1;
#line 567 "chbtrd.f"
				for (j1inc = j1; i__3 < 0 ? j1inc >= i__2 : 
					j1inc <= i__2; j1inc += i__3) {
#line 568 "chbtrd.f"
				    crot_(&kdm1, &ab[(j1inc - 1) * ab_dim1 + 
					    3], &c__1, &ab[j1inc * ab_dim1 + 
					    2], &c__1, &d__[j1inc], &work[
					    j1inc]);
#line 571 "chbtrd.f"
/* L160: */
#line 571 "chbtrd.f"
				}
#line 572 "chbtrd.f"
			    }
/* Computing MIN */
#line 573 "chbtrd.f"
			    i__3 = kdm1, i__2 = *n - j2;
#line 573 "chbtrd.f"
			    lend = min(i__3,i__2);
#line 574 "chbtrd.f"
			    last = j1end + kd1;
#line 575 "chbtrd.f"
			    if (lend > 0) {
#line 575 "chbtrd.f"
				crot_(&lend, &ab[(last - 1) * ab_dim1 + 3], &
					c__1, &ab[last * ab_dim1 + 2], &c__1, 
					&d__[last], &work[last]);
#line 575 "chbtrd.f"
			    }
#line 579 "chbtrd.f"
			}
#line 580 "chbtrd.f"
		    }



#line 584 "chbtrd.f"
		    if (wantq) {

/*                    accumulate product of plane rotations in Q */

#line 588 "chbtrd.f"
			if (initq) {

/*                 take advantage of the fact that Q was */
/*                 initially the Identity matrix */

#line 593 "chbtrd.f"
			    iqend = max(iqend,j2);
/* Computing MAX */
#line 594 "chbtrd.f"
			    i__3 = 0, i__2 = k - 3;
#line 594 "chbtrd.f"
			    i2 = max(i__3,i__2);
#line 595 "chbtrd.f"
			    iqaend = i__ * *kd + 1;
#line 596 "chbtrd.f"
			    if (k == 2) {
#line 596 "chbtrd.f"
				iqaend += *kd;
#line 596 "chbtrd.f"
			    }
#line 598 "chbtrd.f"
			    iqaend = min(iqaend,iqend);
#line 599 "chbtrd.f"
			    i__3 = j2;
#line 599 "chbtrd.f"
			    i__2 = kd1;
#line 599 "chbtrd.f"
			    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j 
				    += i__2) {
#line 600 "chbtrd.f"
				ibl = i__ - i2 / kdm1;
#line 601 "chbtrd.f"
				++i2;
/* Computing MAX */
#line 602 "chbtrd.f"
				i__4 = 1, i__5 = j - ibl;
#line 602 "chbtrd.f"
				iqb = max(i__4,i__5);
#line 603 "chbtrd.f"
				nq = iqaend + 1 - iqb;
/* Computing MIN */
#line 604 "chbtrd.f"
				i__4 = iqaend + *kd;
#line 604 "chbtrd.f"
				iqaend = min(i__4,iqend);
#line 605 "chbtrd.f"
				crot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, 
					&q[iqb + j * q_dim1], &c__1, &d__[j], 
					&work[j]);
#line 607 "chbtrd.f"
/* L170: */
#line 607 "chbtrd.f"
			    }
#line 608 "chbtrd.f"
			} else {

#line 610 "chbtrd.f"
			    i__2 = j2;
#line 610 "chbtrd.f"
			    i__3 = kd1;
#line 610 "chbtrd.f"
			    for (j = j1; i__3 < 0 ? j >= i__2 : j <= i__2; j 
				    += i__3) {
#line 611 "chbtrd.f"
				crot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[
					j * q_dim1 + 1], &c__1, &d__[j], &
					work[j]);
#line 613 "chbtrd.f"
/* L180: */
#line 613 "chbtrd.f"
			    }
#line 614 "chbtrd.f"
			}
#line 615 "chbtrd.f"
		    }

#line 617 "chbtrd.f"
		    if (j2 + kdn > *n) {

/*                    adjust J2 to keep within the bounds of the matrix */

#line 621 "chbtrd.f"
			--nr;
#line 622 "chbtrd.f"
			j2 = j2 - kdn - 1;
#line 623 "chbtrd.f"
		    }

#line 625 "chbtrd.f"
		    i__3 = j2;
#line 625 "chbtrd.f"
		    i__2 = kd1;
#line 625 "chbtrd.f"
		    for (j = j1; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) 
			    {

/*                    create nonzero element a(j+kd,j-1) outside the */
/*                    band and store it in WORK */

#line 630 "chbtrd.f"
			i__4 = j + *kd;
#line 630 "chbtrd.f"
			i__5 = j;
#line 630 "chbtrd.f"
			i__6 = kd1 + j * ab_dim1;
#line 630 "chbtrd.f"
			z__1.r = work[i__5].r * ab[i__6].r - work[i__5].i * 
				ab[i__6].i, z__1.i = work[i__5].r * ab[i__6]
				.i + work[i__5].i * ab[i__6].r;
#line 630 "chbtrd.f"
			work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 631 "chbtrd.f"
			i__4 = kd1 + j * ab_dim1;
#line 631 "chbtrd.f"
			i__5 = j;
#line 631 "chbtrd.f"
			i__6 = kd1 + j * ab_dim1;
#line 631 "chbtrd.f"
			z__1.r = d__[i__5] * ab[i__6].r, z__1.i = d__[i__5] * 
				ab[i__6].i;
#line 631 "chbtrd.f"
			ab[i__4].r = z__1.r, ab[i__4].i = z__1.i;
#line 632 "chbtrd.f"
/* L190: */
#line 632 "chbtrd.f"
		    }
#line 633 "chbtrd.f"
/* L200: */
#line 633 "chbtrd.f"
		}
#line 634 "chbtrd.f"
/* L210: */
#line 634 "chbtrd.f"
	    }
#line 635 "chbtrd.f"
	}

#line 637 "chbtrd.f"
	if (*kd > 0) {

/*           make off-diagonal elements real and copy them to E */

#line 641 "chbtrd.f"
	    i__1 = *n - 1;
#line 641 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 642 "chbtrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 642 "chbtrd.f"
		t.r = ab[i__2].r, t.i = ab[i__2].i;
#line 643 "chbtrd.f"
		abst = z_abs(&t);
#line 644 "chbtrd.f"
		i__2 = i__ * ab_dim1 + 2;
#line 644 "chbtrd.f"
		ab[i__2].r = abst, ab[i__2].i = 0.;
#line 645 "chbtrd.f"
		e[i__] = abst;
#line 646 "chbtrd.f"
		if (abst != 0.) {
#line 647 "chbtrd.f"
		    z__1.r = t.r / abst, z__1.i = t.i / abst;
#line 647 "chbtrd.f"
		    t.r = z__1.r, t.i = z__1.i;
#line 648 "chbtrd.f"
		} else {
#line 649 "chbtrd.f"
		    t.r = 1., t.i = 0.;
#line 650 "chbtrd.f"
		}
#line 651 "chbtrd.f"
		if (i__ < *n - 1) {
#line 651 "chbtrd.f"
		    i__2 = (i__ + 1) * ab_dim1 + 2;
#line 651 "chbtrd.f"
		    i__3 = (i__ + 1) * ab_dim1 + 2;
#line 651 "chbtrd.f"
		    z__1.r = ab[i__3].r * t.r - ab[i__3].i * t.i, z__1.i = ab[
			    i__3].r * t.i + ab[i__3].i * t.r;
#line 651 "chbtrd.f"
		    ab[i__2].r = z__1.r, ab[i__2].i = z__1.i;
#line 651 "chbtrd.f"
		}
#line 653 "chbtrd.f"
		if (wantq) {
#line 654 "chbtrd.f"
		    cscal_(n, &t, &q[(i__ + 1) * q_dim1 + 1], &c__1);
#line 655 "chbtrd.f"
		}
#line 656 "chbtrd.f"
/* L220: */
#line 656 "chbtrd.f"
	    }
#line 657 "chbtrd.f"
	} else {

/*           set E to zero if original matrix was diagonal */

#line 661 "chbtrd.f"
	    i__1 = *n - 1;
#line 661 "chbtrd.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 662 "chbtrd.f"
		e[i__] = 0.;
#line 663 "chbtrd.f"
/* L230: */
#line 663 "chbtrd.f"
	    }
#line 664 "chbtrd.f"
	}

/*        copy diagonal elements to D */

#line 668 "chbtrd.f"
	i__1 = *n;
#line 668 "chbtrd.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 669 "chbtrd.f"
	    i__2 = i__;
#line 669 "chbtrd.f"
	    i__3 = i__ * ab_dim1 + 1;
#line 669 "chbtrd.f"
	    d__[i__2] = ab[i__3].r;
#line 670 "chbtrd.f"
/* L240: */
#line 670 "chbtrd.f"
	}
#line 671 "chbtrd.f"
    }

#line 673 "chbtrd.f"
    return 0;

/*     End of CHBTRD */

} /* chbtrd_ */

