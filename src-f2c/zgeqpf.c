#line 1 "zgeqpf.f"
/* zgeqpf.f -- translated by f2c (version 20100827).
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

#line 1 "zgeqpf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZGEQPF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEQPF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeqpf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeqpf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeqpf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine ZGEQP3. */
/* > */
/* > ZGEQPF computes a QR factorization with column pivoting of a */
/* > complex M-by-N matrix A: A*P = Q*R. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. N >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the upper triangle of the array contains the */
/* >          min(M,N)-by-N upper triangular matrix R; the elements */
/* >          below the diagonal, together with the array TAU, */
/* >          represent the unitary matrix Q as a product of */
/* >          min(m,n) elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* >          JPVT is INTEGER array, dimension (N) */
/* >          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* >          to the front of A*P (a leading column); if JPVT(i) = 0, */
/* >          the i-th column of A is a free column. */
/* >          On exit, if JPVT(i) = k, then the i-th column of A*P */
/* >          was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \ingroup complex16GEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of elementary reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(n) */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). */
/* > */
/* >  The matrix P is represented in jpvt as follows: If */
/* >     jpvt(j) = i */
/* >  then the jth column of P is the ith canonical unit vector. */
/* > */
/* >  Partial column norm updating strategy modified by */
/* >    Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* >    University of Zagreb, Croatia. */
/* >  -- April 2011                                                      -- */
/* >  For more details see LAPACK Working Note 176. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgeqpf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, ma, mn;
    static doublecomplex aii;
    static integer pvt;
    static doublereal temp, temp2, tol3z;
    static integer itemp;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), zswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zgeqr2_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern /* Subroutine */ int zunm2r_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *);


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

/*     Test the input arguments */

#line 191 "zgeqpf.f"
    /* Parameter adjustments */
#line 191 "zgeqpf.f"
    a_dim1 = *lda;
#line 191 "zgeqpf.f"
    a_offset = 1 + a_dim1;
#line 191 "zgeqpf.f"
    a -= a_offset;
#line 191 "zgeqpf.f"
    --jpvt;
#line 191 "zgeqpf.f"
    --tau;
#line 191 "zgeqpf.f"
    --work;
#line 191 "zgeqpf.f"
    --rwork;
#line 191 "zgeqpf.f"

#line 191 "zgeqpf.f"
    /* Function Body */
#line 191 "zgeqpf.f"
    *info = 0;
#line 192 "zgeqpf.f"
    if (*m < 0) {
#line 193 "zgeqpf.f"
	*info = -1;
#line 194 "zgeqpf.f"
    } else if (*n < 0) {
#line 195 "zgeqpf.f"
	*info = -2;
#line 196 "zgeqpf.f"
    } else if (*lda < max(1,*m)) {
#line 197 "zgeqpf.f"
	*info = -4;
#line 198 "zgeqpf.f"
    }
#line 199 "zgeqpf.f"
    if (*info != 0) {
#line 200 "zgeqpf.f"
	i__1 = -(*info);
#line 200 "zgeqpf.f"
	xerbla_("ZGEQPF", &i__1, (ftnlen)6);
#line 201 "zgeqpf.f"
	return 0;
#line 202 "zgeqpf.f"
    }

#line 204 "zgeqpf.f"
    mn = min(*m,*n);
#line 205 "zgeqpf.f"
    tol3z = sqrt(dlamch_("Epsilon", (ftnlen)7));

/*     Move initial columns up front */

#line 209 "zgeqpf.f"
    itemp = 1;
#line 210 "zgeqpf.f"
    i__1 = *n;
#line 210 "zgeqpf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "zgeqpf.f"
	if (jpvt[i__] != 0) {
#line 212 "zgeqpf.f"
	    if (i__ != itemp) {
#line 213 "zgeqpf.f"
		zswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1],
			 &c__1);
#line 214 "zgeqpf.f"
		jpvt[i__] = jpvt[itemp];
#line 215 "zgeqpf.f"
		jpvt[itemp] = i__;
#line 216 "zgeqpf.f"
	    } else {
#line 217 "zgeqpf.f"
		jpvt[i__] = i__;
#line 218 "zgeqpf.f"
	    }
#line 219 "zgeqpf.f"
	    ++itemp;
#line 220 "zgeqpf.f"
	} else {
#line 221 "zgeqpf.f"
	    jpvt[i__] = i__;
#line 222 "zgeqpf.f"
	}
#line 223 "zgeqpf.f"
/* L10: */
#line 223 "zgeqpf.f"
    }
#line 224 "zgeqpf.f"
    --itemp;

/*     Compute the QR factorization and update remaining columns */

#line 228 "zgeqpf.f"
    if (itemp > 0) {
#line 229 "zgeqpf.f"
	ma = min(itemp,*m);
#line 230 "zgeqpf.f"
	zgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
#line 231 "zgeqpf.f"
	if (ma < *n) {
#line 232 "zgeqpf.f"
	    i__1 = *n - ma;
#line 232 "zgeqpf.f"
	    zunm2r_("Left", "Conjugate transpose", m, &i__1, &ma, &a[a_offset]
		    , lda, &tau[1], &a[(ma + 1) * a_dim1 + 1], lda, &work[1], 
		    info, (ftnlen)4, (ftnlen)19);
#line 234 "zgeqpf.f"
	}
#line 235 "zgeqpf.f"
    }

#line 237 "zgeqpf.f"
    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of */
/*        work store the exact column norms. */

#line 242 "zgeqpf.f"
	i__1 = *n;
#line 242 "zgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {
#line 243 "zgeqpf.f"
	    i__2 = *m - itemp;
#line 243 "zgeqpf.f"
	    rwork[i__] = dznrm2_(&i__2, &a[itemp + 1 + i__ * a_dim1], &c__1);
#line 244 "zgeqpf.f"
	    rwork[*n + i__] = rwork[i__];
#line 245 "zgeqpf.f"
/* L20: */
#line 245 "zgeqpf.f"
	}

/*        Compute factorization */

#line 249 "zgeqpf.f"
	i__1 = mn;
#line 249 "zgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {

/*           Determine ith pivot column and swap if necessary */

#line 253 "zgeqpf.f"
	    i__2 = *n - i__ + 1;
#line 253 "zgeqpf.f"
	    pvt = i__ - 1 + idamax_(&i__2, &rwork[i__], &c__1);

#line 255 "zgeqpf.f"
	    if (pvt != i__) {
#line 256 "zgeqpf.f"
		zswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
			c__1);
#line 257 "zgeqpf.f"
		itemp = jpvt[pvt];
#line 258 "zgeqpf.f"
		jpvt[pvt] = jpvt[i__];
#line 259 "zgeqpf.f"
		jpvt[i__] = itemp;
#line 260 "zgeqpf.f"
		rwork[pvt] = rwork[i__];
#line 261 "zgeqpf.f"
		rwork[*n + pvt] = rwork[*n + i__];
#line 262 "zgeqpf.f"
	    }

/*           Generate elementary reflector H(i) */

#line 266 "zgeqpf.f"
	    i__2 = i__ + i__ * a_dim1;
#line 266 "zgeqpf.f"
	    aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 267 "zgeqpf.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 267 "zgeqpf.f"
	    i__3 = i__ + 1;
#line 267 "zgeqpf.f"
	    zlarfg_(&i__2, &aii, &a[min(i__3,*m) + i__ * a_dim1], &c__1, &tau[
		    i__]);
#line 269 "zgeqpf.f"
	    i__2 = i__ + i__ * a_dim1;
#line 269 "zgeqpf.f"
	    a[i__2].r = aii.r, a[i__2].i = aii.i;

#line 271 "zgeqpf.f"
	    if (i__ < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

#line 275 "zgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 275 "zgeqpf.f"
		aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 276 "zgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 276 "zgeqpf.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 277 "zgeqpf.f"
		i__2 = *m - i__ + 1;
#line 277 "zgeqpf.f"
		i__3 = *n - i__;
#line 277 "zgeqpf.f"
		d_cnjg(&z__1, &tau[i__]);
#line 277 "zgeqpf.f"
		zlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			z__1, &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (
			ftnlen)4);
#line 279 "zgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 279 "zgeqpf.f"
		a[i__2].r = aii.r, a[i__2].i = aii.i;
#line 280 "zgeqpf.f"
	    }

/*           Update partial column norms */

#line 284 "zgeqpf.f"
	    i__2 = *n;
#line 284 "zgeqpf.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 285 "zgeqpf.f"
		if (rwork[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 290 "zgeqpf.f"
		    temp = z_abs(&a[i__ + j * a_dim1]) / rwork[j];
/* Computing MAX */
#line 291 "zgeqpf.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 291 "zgeqpf.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 292 "zgeqpf.f"
		    d__1 = rwork[j] / rwork[*n + j];
#line 292 "zgeqpf.f"
		    temp2 = temp * (d__1 * d__1);
#line 293 "zgeqpf.f"
		    if (temp2 <= tol3z) {
#line 294 "zgeqpf.f"
			if (*m - i__ > 0) {
#line 295 "zgeqpf.f"
			    i__3 = *m - i__;
#line 295 "zgeqpf.f"
			    rwork[j] = dznrm2_(&i__3, &a[i__ + 1 + j * a_dim1]
				    , &c__1);
#line 296 "zgeqpf.f"
			    rwork[*n + j] = rwork[j];
#line 297 "zgeqpf.f"
			} else {
#line 298 "zgeqpf.f"
			    rwork[j] = 0.;
#line 299 "zgeqpf.f"
			    rwork[*n + j] = 0.;
#line 300 "zgeqpf.f"
			}
#line 301 "zgeqpf.f"
		    } else {
#line 302 "zgeqpf.f"
			rwork[j] *= sqrt(temp);
#line 303 "zgeqpf.f"
		    }
#line 304 "zgeqpf.f"
		}
#line 305 "zgeqpf.f"
/* L30: */
#line 305 "zgeqpf.f"
	    }

#line 307 "zgeqpf.f"
/* L40: */
#line 307 "zgeqpf.f"
	}
#line 308 "zgeqpf.f"
    }
#line 309 "zgeqpf.f"
    return 0;

/*     End of ZGEQPF */

} /* zgeqpf_ */

