#line 1 "dgeqpf.f"
/* dgeqpf.f -- translated by f2c (version 20100827).
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

#line 1 "dgeqpf.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DGEQPF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEQPF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqpf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqpf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqpf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGEQP3. */
/* > */
/* > DGEQPF computes a QR factorization with column pivoting of a */
/* > real M-by-N matrix A: A*P = Q*R. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the upper triangle of the array contains the */
/* >          min(M,N)-by-N upper triangular matrix R; the elements */
/* >          below the diagonal, together with the array TAU, */
/* >          represent the orthogonal matrix Q as a product of */
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
/* >          TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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

/* > \ingroup doubleGEcomputational */

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
/* >     H = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
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
/* Subroutine */ int dgeqpf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, ma, mn;
    static doublereal aii;
    static integer pvt;
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2, tol3z;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dgeqr2_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dorm2r_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 183 "dgeqpf.f"
    /* Parameter adjustments */
#line 183 "dgeqpf.f"
    a_dim1 = *lda;
#line 183 "dgeqpf.f"
    a_offset = 1 + a_dim1;
#line 183 "dgeqpf.f"
    a -= a_offset;
#line 183 "dgeqpf.f"
    --jpvt;
#line 183 "dgeqpf.f"
    --tau;
#line 183 "dgeqpf.f"
    --work;
#line 183 "dgeqpf.f"

#line 183 "dgeqpf.f"
    /* Function Body */
#line 183 "dgeqpf.f"
    *info = 0;
#line 184 "dgeqpf.f"
    if (*m < 0) {
#line 185 "dgeqpf.f"
	*info = -1;
#line 186 "dgeqpf.f"
    } else if (*n < 0) {
#line 187 "dgeqpf.f"
	*info = -2;
#line 188 "dgeqpf.f"
    } else if (*lda < max(1,*m)) {
#line 189 "dgeqpf.f"
	*info = -4;
#line 190 "dgeqpf.f"
    }
#line 191 "dgeqpf.f"
    if (*info != 0) {
#line 192 "dgeqpf.f"
	i__1 = -(*info);
#line 192 "dgeqpf.f"
	xerbla_("DGEQPF", &i__1, (ftnlen)6);
#line 193 "dgeqpf.f"
	return 0;
#line 194 "dgeqpf.f"
    }

#line 196 "dgeqpf.f"
    mn = min(*m,*n);
#line 197 "dgeqpf.f"
    tol3z = sqrt(dlamch_("Epsilon", (ftnlen)7));

/*     Move initial columns up front */

#line 201 "dgeqpf.f"
    itemp = 1;
#line 202 "dgeqpf.f"
    i__1 = *n;
#line 202 "dgeqpf.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "dgeqpf.f"
	if (jpvt[i__] != 0) {
#line 204 "dgeqpf.f"
	    if (i__ != itemp) {
#line 205 "dgeqpf.f"
		dswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1],
			 &c__1);
#line 206 "dgeqpf.f"
		jpvt[i__] = jpvt[itemp];
#line 207 "dgeqpf.f"
		jpvt[itemp] = i__;
#line 208 "dgeqpf.f"
	    } else {
#line 209 "dgeqpf.f"
		jpvt[i__] = i__;
#line 210 "dgeqpf.f"
	    }
#line 211 "dgeqpf.f"
	    ++itemp;
#line 212 "dgeqpf.f"
	} else {
#line 213 "dgeqpf.f"
	    jpvt[i__] = i__;
#line 214 "dgeqpf.f"
	}
#line 215 "dgeqpf.f"
/* L10: */
#line 215 "dgeqpf.f"
    }
#line 216 "dgeqpf.f"
    --itemp;

/*     Compute the QR factorization and update remaining columns */

#line 220 "dgeqpf.f"
    if (itemp > 0) {
#line 221 "dgeqpf.f"
	ma = min(itemp,*m);
#line 222 "dgeqpf.f"
	dgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
#line 223 "dgeqpf.f"
	if (ma < *n) {
#line 224 "dgeqpf.f"
	    i__1 = *n - ma;
#line 224 "dgeqpf.f"
	    dorm2r_("Left", "Transpose", m, &i__1, &ma, &a[a_offset], lda, &
		    tau[1], &a[(ma + 1) * a_dim1 + 1], lda, &work[1], info, (
		    ftnlen)4, (ftnlen)9);
#line 226 "dgeqpf.f"
	}
#line 227 "dgeqpf.f"
    }

#line 229 "dgeqpf.f"
    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of */
/*        work store the exact column norms. */

#line 234 "dgeqpf.f"
	i__1 = *n;
#line 234 "dgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {
#line 235 "dgeqpf.f"
	    i__2 = *m - itemp;
#line 235 "dgeqpf.f"
	    work[i__] = dnrm2_(&i__2, &a[itemp + 1 + i__ * a_dim1], &c__1);
#line 236 "dgeqpf.f"
	    work[*n + i__] = work[i__];
#line 237 "dgeqpf.f"
/* L20: */
#line 237 "dgeqpf.f"
	}

/*        Compute factorization */

#line 241 "dgeqpf.f"
	i__1 = mn;
#line 241 "dgeqpf.f"
	for (i__ = itemp + 1; i__ <= i__1; ++i__) {

/*           Determine ith pivot column and swap if necessary */

#line 245 "dgeqpf.f"
	    i__2 = *n - i__ + 1;
#line 245 "dgeqpf.f"
	    pvt = i__ - 1 + idamax_(&i__2, &work[i__], &c__1);

#line 247 "dgeqpf.f"
	    if (pvt != i__) {
#line 248 "dgeqpf.f"
		dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
			c__1);
#line 249 "dgeqpf.f"
		itemp = jpvt[pvt];
#line 250 "dgeqpf.f"
		jpvt[pvt] = jpvt[i__];
#line 251 "dgeqpf.f"
		jpvt[i__] = itemp;
#line 252 "dgeqpf.f"
		work[pvt] = work[i__];
#line 253 "dgeqpf.f"
		work[*n + pvt] = work[*n + i__];
#line 254 "dgeqpf.f"
	    }

/*           Generate elementary reflector H(i) */

#line 258 "dgeqpf.f"
	    if (i__ < *m) {
#line 259 "dgeqpf.f"
		i__2 = *m - i__ + 1;
#line 259 "dgeqpf.f"
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * 
			a_dim1], &c__1, &tau[i__]);
#line 260 "dgeqpf.f"
	    } else {
#line 261 "dgeqpf.f"
		dlarfg_(&c__1, &a[*m + *m * a_dim1], &a[*m + *m * a_dim1], &
			c__1, &tau[*m]);
#line 262 "dgeqpf.f"
	    }

#line 264 "dgeqpf.f"
	    if (i__ < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

#line 268 "dgeqpf.f"
		aii = a[i__ + i__ * a_dim1];
#line 269 "dgeqpf.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 270 "dgeqpf.f"
		i__2 = *m - i__ + 1;
#line 270 "dgeqpf.f"
		i__3 = *n - i__;
#line 270 "dgeqpf.f"
		dlarf_("LEFT", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[(*
			n << 1) + 1], (ftnlen)4);
#line 272 "dgeqpf.f"
		a[i__ + i__ * a_dim1] = aii;
#line 273 "dgeqpf.f"
	    }

/*           Update partial column norms */

#line 277 "dgeqpf.f"
	    i__2 = *n;
#line 277 "dgeqpf.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 278 "dgeqpf.f"
		if (work[j] != 0.) {

/*                 NOTE: The following 4 lines follow from the analysis in */
/*                 Lapack Working Note 176. */

#line 283 "dgeqpf.f"
		    temp = (d__1 = a[i__ + j * a_dim1], abs(d__1)) / work[j];
/* Computing MAX */
#line 284 "dgeqpf.f"
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
#line 284 "dgeqpf.f"
		    temp = max(d__1,d__2);
/* Computing 2nd power */
#line 285 "dgeqpf.f"
		    d__1 = work[j] / work[*n + j];
#line 285 "dgeqpf.f"
		    temp2 = temp * (d__1 * d__1);
#line 286 "dgeqpf.f"
		    if (temp2 <= tol3z) {
#line 287 "dgeqpf.f"
			if (*m - i__ > 0) {
#line 288 "dgeqpf.f"
			    i__3 = *m - i__;
#line 288 "dgeqpf.f"
			    work[j] = dnrm2_(&i__3, &a[i__ + 1 + j * a_dim1], 
				    &c__1);
#line 289 "dgeqpf.f"
			    work[*n + j] = work[j];
#line 290 "dgeqpf.f"
			} else {
#line 291 "dgeqpf.f"
			    work[j] = 0.;
#line 292 "dgeqpf.f"
			    work[*n + j] = 0.;
#line 293 "dgeqpf.f"
			}
#line 294 "dgeqpf.f"
		    } else {
#line 295 "dgeqpf.f"
			work[j] *= sqrt(temp);
#line 296 "dgeqpf.f"
		    }
#line 297 "dgeqpf.f"
		}
#line 298 "dgeqpf.f"
/* L30: */
#line 298 "dgeqpf.f"
	    }

#line 300 "dgeqpf.f"
/* L40: */
#line 300 "dgeqpf.f"
	}
#line 301 "dgeqpf.f"
    }
#line 302 "dgeqpf.f"
    return 0;

/*     End of DGEQPF */

} /* dgeqpf_ */

