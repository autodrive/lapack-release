#line 1 "dlaqp2.f"
/* dlaqp2.f -- translated by f2c (version 20100827).
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

#line 1 "dlaqp2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLAQP2 computes a QR factorization with column pivoting of the matrix block. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAQP2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, */
/*                          WORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, M, N, OFFSET */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            JPVT( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), VN1( * ), VN2( * ), */
/*      $                   WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQP2 computes a QR factorization with column pivoting of */
/* > the block A(OFFSET+1:M,1:N). */
/* > The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized. */
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
/* >          The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* >          OFFSET is INTEGER */
/* >          The number of rows of the matrix A that must be pivoted */
/* >          but no factorized. OFFSET >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is */
/* >          the triangular factor obtained; the elements in block */
/* >          A(OFFSET+1:M,1:N) below the diagonal, together with the */
/* >          array TAU, represent the orthogonal matrix Q as a product of */
/* >          elementary reflectors. Block A(1:OFFSET,1:N) has been */
/* >          accordingly pivoted, but no factorized. */
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
/* > \param[in,out] VN1 */
/* > \verbatim */
/* >          VN1 is DOUBLE PRECISION array, dimension (N) */
/* >          The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* >          VN2 is DOUBLE PRECISION array, dimension (N) */
/* >          The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* >    X. Sun, Computer Science Dept., Duke University, USA */
/* > \n */
/* >  Partial column norm updating strategy modified on April 2011 */
/* >    Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* >    University of Zagreb, Croatia. */

/* > \par References: */
/*  ================ */
/* > */
/* > LAPACK Working Note 176 */

/* > \htmlonly */
/* > <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> */
/* > \endhtmlonly */

/*  ===================================================================== */
/* Subroutine */ int dlaqp2_(integer *m, integer *n, integer *offset, 
	doublereal *a, integer *lda, integer *jpvt, doublereal *tau, 
	doublereal *vn1, doublereal *vn2, doublereal *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, mn;
    static doublereal aii;
    static integer pvt;
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2, tol3z;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static integer offpi, itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 189 "dlaqp2.f"
    /* Parameter adjustments */
#line 189 "dlaqp2.f"
    a_dim1 = *lda;
#line 189 "dlaqp2.f"
    a_offset = 1 + a_dim1;
#line 189 "dlaqp2.f"
    a -= a_offset;
#line 189 "dlaqp2.f"
    --jpvt;
#line 189 "dlaqp2.f"
    --tau;
#line 189 "dlaqp2.f"
    --vn1;
#line 189 "dlaqp2.f"
    --vn2;
#line 189 "dlaqp2.f"
    --work;
#line 189 "dlaqp2.f"

#line 189 "dlaqp2.f"
    /* Function Body */
/* Computing MIN */
#line 189 "dlaqp2.f"
    i__1 = *m - *offset;
#line 189 "dlaqp2.f"
    mn = min(i__1,*n);
#line 190 "dlaqp2.f"
    tol3z = sqrt(dlamch_("Epsilon", (ftnlen)7));

/*     Compute factorization. */

#line 194 "dlaqp2.f"
    i__1 = mn;
#line 194 "dlaqp2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 196 "dlaqp2.f"
	offpi = *offset + i__;

/*        Determine ith pivot column and swap if necessary. */

#line 200 "dlaqp2.f"
	i__2 = *n - i__ + 1;
#line 200 "dlaqp2.f"
	pvt = i__ - 1 + idamax_(&i__2, &vn1[i__], &c__1);

#line 202 "dlaqp2.f"
	if (pvt != i__) {
#line 203 "dlaqp2.f"
	    dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
#line 204 "dlaqp2.f"
	    itemp = jpvt[pvt];
#line 205 "dlaqp2.f"
	    jpvt[pvt] = jpvt[i__];
#line 206 "dlaqp2.f"
	    jpvt[i__] = itemp;
#line 207 "dlaqp2.f"
	    vn1[pvt] = vn1[i__];
#line 208 "dlaqp2.f"
	    vn2[pvt] = vn2[i__];
#line 209 "dlaqp2.f"
	}

/*        Generate elementary reflector H(i). */

#line 213 "dlaqp2.f"
	if (offpi < *m) {
#line 214 "dlaqp2.f"
	    i__2 = *m - offpi + 1;
#line 214 "dlaqp2.f"
	    dlarfg_(&i__2, &a[offpi + i__ * a_dim1], &a[offpi + 1 + i__ * 
		    a_dim1], &c__1, &tau[i__]);
#line 216 "dlaqp2.f"
	} else {
#line 217 "dlaqp2.f"
	    dlarfg_(&c__1, &a[*m + i__ * a_dim1], &a[*m + i__ * a_dim1], &
		    c__1, &tau[i__]);
#line 218 "dlaqp2.f"
	}

#line 220 "dlaqp2.f"
	if (i__ < *n) {

/*           Apply H(i)**T to A(offset+i:m,i+1:n) from the left. */

#line 224 "dlaqp2.f"
	    aii = a[offpi + i__ * a_dim1];
#line 225 "dlaqp2.f"
	    a[offpi + i__ * a_dim1] = 1.;
#line 226 "dlaqp2.f"
	    i__2 = *m - offpi + 1;
#line 226 "dlaqp2.f"
	    i__3 = *n - i__;
#line 226 "dlaqp2.f"
	    dlarf_("Left", &i__2, &i__3, &a[offpi + i__ * a_dim1], &c__1, &
		    tau[i__], &a[offpi + (i__ + 1) * a_dim1], lda, &work[1], (
		    ftnlen)4);
#line 228 "dlaqp2.f"
	    a[offpi + i__ * a_dim1] = aii;
#line 229 "dlaqp2.f"
	}

/*        Update partial column norms. */

#line 233 "dlaqp2.f"
	i__2 = *n;
#line 233 "dlaqp2.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 234 "dlaqp2.f"
	    if (vn1[j] != 0.) {

/*              NOTE: The following 4 lines follow from the analysis in */
/*              Lapack Working Note 176. */

/* Computing 2nd power */
#line 239 "dlaqp2.f"
		d__2 = (d__1 = a[offpi + j * a_dim1], abs(d__1)) / vn1[j];
#line 239 "dlaqp2.f"
		temp = 1. - d__2 * d__2;
#line 240 "dlaqp2.f"
		temp = max(temp,0.);
/* Computing 2nd power */
#line 241 "dlaqp2.f"
		d__1 = vn1[j] / vn2[j];
#line 241 "dlaqp2.f"
		temp2 = temp * (d__1 * d__1);
#line 242 "dlaqp2.f"
		if (temp2 <= tol3z) {
#line 243 "dlaqp2.f"
		    if (offpi < *m) {
#line 244 "dlaqp2.f"
			i__3 = *m - offpi;
#line 244 "dlaqp2.f"
			vn1[j] = dnrm2_(&i__3, &a[offpi + 1 + j * a_dim1], &
				c__1);
#line 245 "dlaqp2.f"
			vn2[j] = vn1[j];
#line 246 "dlaqp2.f"
		    } else {
#line 247 "dlaqp2.f"
			vn1[j] = 0.;
#line 248 "dlaqp2.f"
			vn2[j] = 0.;
#line 249 "dlaqp2.f"
		    }
#line 250 "dlaqp2.f"
		} else {
#line 251 "dlaqp2.f"
		    vn1[j] *= sqrt(temp);
#line 252 "dlaqp2.f"
		}
#line 253 "dlaqp2.f"
	    }
#line 254 "dlaqp2.f"
/* L10: */
#line 254 "dlaqp2.f"
	}

#line 256 "dlaqp2.f"
/* L20: */
#line 256 "dlaqp2.f"
    }

#line 258 "dlaqp2.f"
    return 0;

/*     End of DLAQP2 */

} /* dlaqp2_ */

