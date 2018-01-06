#line 1 "cgeqpf.f"
/* cgeqpf.f -- translated by f2c (version 20100827).
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

#line 1 "cgeqpf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CGEQPF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEQPF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqpf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqpf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqpf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine CGEQP3. */
/* > */
/* > CGEQPF computes a QR factorization with column pivoting of a */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          TAU is COMPLEX array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
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

/* > \ingroup complexGEcomputational */

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
/* Subroutine */ int cgeqpf_(integer *m, integer *n, doublecomplex *a, 
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
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), cswap_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer itemp;
    extern /* Subroutine */ int cgeqr2_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *);
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cunm2r_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), clarfg_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, doublereal *, integer *);


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

#line 191 "cgeqpf.f"
    /* Parameter adjustments */
#line 191 "cgeqpf.f"
    a_dim1 = *lda;
#line 191 "cgeqpf.f"
    a_offset = 1 + a_dim1;
#line 191 "cgeqpf.f"
    a -= a_offset;
#line 191 "cgeqpf.f"
    --jpvt;
#line 191 "cgeqpf.f"
    --tau;
#line 191 "cgeqpf.f"
    --work;
#line 191 "cgeqpf.f"
    --rwork;
#line 191 "cgeqpf.f"

#line 191 "cgeqpf.f"
    /* Function Body */
#line 191 "cgeqpf.f"
    *info = 0;
#line 192 "cgeqpf.f"
    if (*m < 0) {
#line 193 "cgeqpf.f"
	*info = -1;
#line 194 "cgeqpf.f"
    } else if (*n < 0) {
#line 195 "cgeqpf.f"
	*info = -2;
#line 196 "cgeqpf.f"
    } else if (*lda < max(1,*m)) {
#line 197 "cgeqpf.f"
	*info = -4;
#line 198 "cgeqpf.f"
    }
#line 199 "cgeqpf.f"
    if (*info != 0) {
#line 200 "cgeqpf.f"
	i__1 = -(*info);
#line 200 "cgeqpf.f"
	xerbla_("CGEQPF", &i__1, (ftnlen)6);
#line 201 "cgeqpf.f"
	return 0;
#line 202 "cgeqpf.f"
    }

#line 204 "cgeqpf.f"
    mn = min(*m,*n);
#line 205 "cgeqpf.f"
    tol3z = sqrt(slamch_("Epsilon", (ftnlen)7));

/*     Move initial columns up front */

#line 209 "cgeqpf.f"
    itemp = 1;
#line 210 "cgeqpf.f"
    i__1 = *n;
#line 210 "cgeqpf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "cgeqpf.f"
	if (jpvt[i__] != 0) {
#line 212 "cgeqpf.f"
	    if (i__ != itemp) {
#line 213 "cgeqpf.f"
		cswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1],
			 &c__1);
#line 214 "cgeqpf.f"
		jpvt[i__] = jpvt[itemp];
#line 215 "cgeqpf.f"
		jpvt[itemp] = i__;
#line 216 "cgeqpf.f"
	    } else {
#line 217 "cgeqpf.f"
		jpvt[i__] = i__;
#line 218 "cgeqpf.f"
	    }
#line 219 "cgeqpf.f"
	    ++itemp;
#line 220 "cgeqpf.f"
	} else {
#line 221 "cgeqpf.f"
	    jpvt[i__] = i__;
#line 222 "cgeqpf.f"
	}
#line 223 "cgeqpf.f"
/* L10: */
#line 223 "cgeqpf.f"
    }
#line 224 "cgeqpf.f"
    --itemp;

/*     Compute the QR factorization and update remaining columns */

#line 228 "cgeqpf.f"
    if (itemp > 0) {
#line 229 "cgeqpf.f"
	ma = min(itemp,*m);
#line 230 "cgeqpf.f"
	cgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
#line 231 "cgeqpf.f"
	if (ma < *n) {
#line 232 "cgeqpf.f"
	    i__1 = *n - ma;
#line 232 "cgeqpf.f"
	    cunm2r_("Left", "Conjugate transpose", m, &i__1, &ma, &a[a_offset]
		    , lda, &tau[1], &a[(ma + 1) * a_dim1 + 1], lda, &work[1], 
		    info, (ftnlen)4, (ftnlen)19);
#line 234 "cgeqpf.f"
	}
#line 235 "cgeqpf.f"
    }

#line 237 "cgeqpf.f"
    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of */
/*        work store the exact column norms. */

#line 242 "cgeqpf.f"
	i__1 = *n;
#line 242 "cgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {
#line 243 "cgeqpf.f"
	    i__2 = *m - itemp;
#line 243 "cgeqpf.f"
	    rwork[i__] = scnrm2_(&i__2, &a[itemp + 1 + i__ * a_dim1], &c__1);
#line 244 "cgeqpf.f"
	    rwork[*n + i__] = rwork[i__];
#line 245 "cgeqpf.f"
/* L20: */
#line 245 "cgeqpf.f"
	}

/*        Compute factorization */

#line 249 "cgeqpf.f"
	i__1 = mn;
#line 249 "cgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {

/*           Determine ith pivot column and swap if necessary */

#line 253 "cgeqpf.f"
	    i__2 = *n - i__ + 1;
#line 253 "cgeqpf.f"
	    pvt = i__ - 1 + isamax_(&i__2, &rwork[i__], &c__1);

#line 255 "cgeqpf.f"
	    if (pvt != i__) {
#line 256 "cgeqpf.f"
		cswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
			c__1);
#line 257 "cgeqpf.f"
		itemp = jpvt[pvt];
#line 258 "cgeqpf.f"
		jpvt[pvt] = jpvt[i__];
#line 259 "cgeqpf.f"
		jpvt[i__] = itemp;
#line 260 "cgeqpf.f"
		rwork[pvt] = rwork[i__];
#line 261 "cgeqpf.f"
		rwork[*n + pvt] = rwork[*n + i__];
#line 262 "cgeqpf.f"
	    }

/*           Generate elementary reflector H(i) */

#line 266 "cgeqpf.f"
	    i__2 = i__ + i__ * a_dim1;
#line 266 "cgeqpf.f"
	    aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 267 "cgeqpf.f"
	    i__2 = *m - i__ + 1;
/* Computing MIN */
#line 267 "cgeqpf.f"
	    i__3 = i__ + 1;
#line 267 "cgeqpf.f"
	    clarfg_(&i__2, &aii, &a[min(i__3,*m) + i__ * a_dim1], &c__1, &tau[
		    i__]);
#line 269 "cgeqpf.f"
	    i__2 = i__ + i__ * a_dim1;
#line 269 "cgeqpf.f"
	    a[i__2].r = aii.r, a[i__2].i = aii.i;

#line 271 "cgeqpf.f"
	    if (i__ < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

#line 275 "cgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 275 "cgeqpf.f"
		aii.r = a[i__2].r, aii.i = a[i__2].i;
#line 276 "cgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 276 "cgeqpf.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 277 "cgeqpf.f"
		i__2 = *m - i__ + 1;
#line 277 "cgeqpf.f"
		i__3 = *n - i__;
#line 277 "cgeqpf.f"
		d_cnjg(&z__1, &tau[i__]);
#line 277 "cgeqpf.f"
		clarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			z__1, &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (
			ftnlen)4);
#line 279 "cgeqpf.f"
		i__2 = i__ + i__ * a_dim1;
#line 279 "cgeqpf.f"
		a[i__2].r = aii.r, a[i__2].i = aii.i;
#line 280 "cgeqpf.f"
	    }

/*           Update partial column norms */

#line 284 "cgeqpf.f"
	    i__2 = *n;
#line 284 "cgeqpf.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 285 "cgeqpf.f"
		if (rwork[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 290 "cgeqpf.f"
		    temp = z_abs(&a[i__ + j * a_dim1]) / rwork[j];
/* Computing MAX */
#line 291 "cgeqpf.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 291 "cgeqpf.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 292 "cgeqpf.f"
		    d__1 = rwork[j] / rwork[*n + j];
#line 292 "cgeqpf.f"
		    temp2 = temp * (d__1 * d__1);
#line 293 "cgeqpf.f"
		    if (temp2 <= tol3z) {
#line 294 "cgeqpf.f"
			if (*m - i__ > 0) {
#line 295 "cgeqpf.f"
			    i__3 = *m - i__;
#line 295 "cgeqpf.f"
			    rwork[j] = scnrm2_(&i__3, &a[i__ + 1 + j * a_dim1]
				    , &c__1);
#line 296 "cgeqpf.f"
			    rwork[*n + j] = rwork[j];
#line 297 "cgeqpf.f"
			} else {
#line 298 "cgeqpf.f"
			    rwork[j] = 0.;
#line 299 "cgeqpf.f"
			    rwork[*n + j] = 0.;
#line 300 "cgeqpf.f"
			}
#line 301 "cgeqpf.f"
		    } else {
#line 302 "cgeqpf.f"
			rwork[j] *= sqrt(temp);
#line 303 "cgeqpf.f"
		    }
#line 304 "cgeqpf.f"
		}
#line 305 "cgeqpf.f"
/* L30: */
#line 305 "cgeqpf.f"
	    }

#line 307 "cgeqpf.f"
/* L40: */
#line 307 "cgeqpf.f"
	}
#line 308 "cgeqpf.f"
    }
#line 309 "cgeqpf.f"
    return 0;

/*     End of CGEQPF */

} /* cgeqpf_ */

