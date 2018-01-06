#line 1 "dlaswp.f"
/* dlaswp.f -- translated by f2c (version 20100827).
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

#line 1 "dlaswp.f"
/* > \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASWP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, K1, K2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASWP performs a series of row interchanges on the matrix A. */
/* > One row interchange is initiated for each of rows K1 through K2 of A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the matrix of column dimension N to which the row */
/* >          interchanges will be applied. */
/* >          On exit, the permuted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* > \endverbatim */
/* > */
/* > \param[in] K1 */
/* > \verbatim */
/* >          K1 is INTEGER */
/* >          The first element of IPIV for which a row interchange will */
/* >          be done. */
/* > \endverbatim */
/* > */
/* > \param[in] K2 */
/* > \verbatim */
/* >          K2 is INTEGER */
/* >          (K2-K1+1) is the number of elements of IPIV for which a row */
/* >          interchange will be done. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX)) */
/* >          The vector of pivot indices. Only the elements in positions */
/* >          K1 through K1+(K2-K1)*INCX of IPIV are accessed. */
/* >          IPIV(K) = L implies rows K and L are to be interchanged. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of IPIV.  If IPIV */
/* >          is negative, the pivots are applied in reverse order. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by */
/* >   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaswp_(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    static doublereal temp;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

#line 140 "dlaswp.f"
    /* Parameter adjustments */
#line 140 "dlaswp.f"
    a_dim1 = *lda;
#line 140 "dlaswp.f"
    a_offset = 1 + a_dim1;
#line 140 "dlaswp.f"
    a -= a_offset;
#line 140 "dlaswp.f"
    --ipiv;
#line 140 "dlaswp.f"

#line 140 "dlaswp.f"
    /* Function Body */
#line 140 "dlaswp.f"
    if (*incx > 0) {
#line 141 "dlaswp.f"
	ix0 = *k1;
#line 142 "dlaswp.f"
	i1 = *k1;
#line 143 "dlaswp.f"
	i2 = *k2;
#line 144 "dlaswp.f"
	inc = 1;
#line 145 "dlaswp.f"
    } else if (*incx < 0) {
#line 146 "dlaswp.f"
	ix0 = *k1 + (*k1 - *k2) * *incx;
#line 147 "dlaswp.f"
	i1 = *k2;
#line 148 "dlaswp.f"
	i2 = *k1;
#line 149 "dlaswp.f"
	inc = -1;
#line 150 "dlaswp.f"
    } else {
#line 151 "dlaswp.f"
	return 0;
#line 152 "dlaswp.f"
    }

#line 154 "dlaswp.f"
    n32 = *n / 32 << 5;
#line 155 "dlaswp.f"
    if (n32 != 0) {
#line 156 "dlaswp.f"
	i__1 = n32;
#line 156 "dlaswp.f"
	for (j = 1; j <= i__1; j += 32) {
#line 157 "dlaswp.f"
	    ix = ix0;
#line 158 "dlaswp.f"
	    i__2 = i2;
#line 158 "dlaswp.f"
	    i__3 = inc;
#line 158 "dlaswp.f"
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
#line 159 "dlaswp.f"
		ip = ipiv[ix];
#line 160 "dlaswp.f"
		if (ip != i__) {
#line 161 "dlaswp.f"
		    i__4 = j + 31;
#line 161 "dlaswp.f"
		    for (k = j; k <= i__4; ++k) {
#line 162 "dlaswp.f"
			temp = a[i__ + k * a_dim1];
#line 163 "dlaswp.f"
			a[i__ + k * a_dim1] = a[ip + k * a_dim1];
#line 164 "dlaswp.f"
			a[ip + k * a_dim1] = temp;
#line 165 "dlaswp.f"
/* L10: */
#line 165 "dlaswp.f"
		    }
#line 166 "dlaswp.f"
		}
#line 167 "dlaswp.f"
		ix += *incx;
#line 168 "dlaswp.f"
/* L20: */
#line 168 "dlaswp.f"
	    }
#line 169 "dlaswp.f"
/* L30: */
#line 169 "dlaswp.f"
	}
#line 170 "dlaswp.f"
    }
#line 171 "dlaswp.f"
    if (n32 != *n) {
#line 172 "dlaswp.f"
	++n32;
#line 173 "dlaswp.f"
	ix = ix0;
#line 174 "dlaswp.f"
	i__1 = i2;
#line 174 "dlaswp.f"
	i__3 = inc;
#line 174 "dlaswp.f"
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
#line 175 "dlaswp.f"
	    ip = ipiv[ix];
#line 176 "dlaswp.f"
	    if (ip != i__) {
#line 177 "dlaswp.f"
		i__2 = *n;
#line 177 "dlaswp.f"
		for (k = n32; k <= i__2; ++k) {
#line 178 "dlaswp.f"
		    temp = a[i__ + k * a_dim1];
#line 179 "dlaswp.f"
		    a[i__ + k * a_dim1] = a[ip + k * a_dim1];
#line 180 "dlaswp.f"
		    a[ip + k * a_dim1] = temp;
#line 181 "dlaswp.f"
/* L40: */
#line 181 "dlaswp.f"
		}
#line 182 "dlaswp.f"
	    }
#line 183 "dlaswp.f"
	    ix += *incx;
#line 184 "dlaswp.f"
/* L50: */
#line 184 "dlaswp.f"
	}
#line 185 "dlaswp.f"
    }

#line 187 "dlaswp.f"
    return 0;

/*     End of DLASWP */

} /* dlaswp_ */

