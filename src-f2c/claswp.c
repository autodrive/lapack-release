#line 1 "claswp.f"
/* claswp.f -- translated by f2c (version 20100827).
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

#line 1 "claswp.f"
/* > \brief \b CLASWP performs a series of row interchanges on a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLASWP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claswp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claswp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claswp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, K1, K2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASWP performs a series of row interchanges on the matrix A. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          The last element of IPIV for which a row interchange will */
/* >          be done. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (K2*abs(INCX)) */
/* >          The vector of pivot indices.  Only the elements in positions */
/* >          K1 through K2 of IPIV are accessed. */
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

/* > \date September 2012 */

/* > \ingroup complexOTHERauxiliary */

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
/* Subroutine */ int claswp_(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    static doublecomplex temp;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

#line 140 "claswp.f"
    /* Parameter adjustments */
#line 140 "claswp.f"
    a_dim1 = *lda;
#line 140 "claswp.f"
    a_offset = 1 + a_dim1;
#line 140 "claswp.f"
    a -= a_offset;
#line 140 "claswp.f"
    --ipiv;
#line 140 "claswp.f"

#line 140 "claswp.f"
    /* Function Body */
#line 140 "claswp.f"
    if (*incx > 0) {
#line 141 "claswp.f"
	ix0 = *k1;
#line 142 "claswp.f"
	i1 = *k1;
#line 143 "claswp.f"
	i2 = *k2;
#line 144 "claswp.f"
	inc = 1;
#line 145 "claswp.f"
    } else if (*incx < 0) {
#line 146 "claswp.f"
	ix0 = (1 - *k2) * *incx + 1;
#line 147 "claswp.f"
	i1 = *k2;
#line 148 "claswp.f"
	i2 = *k1;
#line 149 "claswp.f"
	inc = -1;
#line 150 "claswp.f"
    } else {
#line 151 "claswp.f"
	return 0;
#line 152 "claswp.f"
    }

#line 154 "claswp.f"
    n32 = *n / 32 << 5;
#line 155 "claswp.f"
    if (n32 != 0) {
#line 156 "claswp.f"
	i__1 = n32;
#line 156 "claswp.f"
	for (j = 1; j <= i__1; j += 32) {
#line 157 "claswp.f"
	    ix = ix0;
#line 158 "claswp.f"
	    i__2 = i2;
#line 158 "claswp.f"
	    i__3 = inc;
#line 158 "claswp.f"
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
#line 159 "claswp.f"
		ip = ipiv[ix];
#line 160 "claswp.f"
		if (ip != i__) {
#line 161 "claswp.f"
		    i__4 = j + 31;
#line 161 "claswp.f"
		    for (k = j; k <= i__4; ++k) {
#line 162 "claswp.f"
			i__5 = i__ + k * a_dim1;
#line 162 "claswp.f"
			temp.r = a[i__5].r, temp.i = a[i__5].i;
#line 163 "claswp.f"
			i__5 = i__ + k * a_dim1;
#line 163 "claswp.f"
			i__6 = ip + k * a_dim1;
#line 163 "claswp.f"
			a[i__5].r = a[i__6].r, a[i__5].i = a[i__6].i;
#line 164 "claswp.f"
			i__5 = ip + k * a_dim1;
#line 164 "claswp.f"
			a[i__5].r = temp.r, a[i__5].i = temp.i;
#line 165 "claswp.f"
/* L10: */
#line 165 "claswp.f"
		    }
#line 166 "claswp.f"
		}
#line 167 "claswp.f"
		ix += *incx;
#line 168 "claswp.f"
/* L20: */
#line 168 "claswp.f"
	    }
#line 169 "claswp.f"
/* L30: */
#line 169 "claswp.f"
	}
#line 170 "claswp.f"
    }
#line 171 "claswp.f"
    if (n32 != *n) {
#line 172 "claswp.f"
	++n32;
#line 173 "claswp.f"
	ix = ix0;
#line 174 "claswp.f"
	i__1 = i2;
#line 174 "claswp.f"
	i__3 = inc;
#line 174 "claswp.f"
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
#line 175 "claswp.f"
	    ip = ipiv[ix];
#line 176 "claswp.f"
	    if (ip != i__) {
#line 177 "claswp.f"
		i__2 = *n;
#line 177 "claswp.f"
		for (k = n32; k <= i__2; ++k) {
#line 178 "claswp.f"
		    i__4 = i__ + k * a_dim1;
#line 178 "claswp.f"
		    temp.r = a[i__4].r, temp.i = a[i__4].i;
#line 179 "claswp.f"
		    i__4 = i__ + k * a_dim1;
#line 179 "claswp.f"
		    i__5 = ip + k * a_dim1;
#line 179 "claswp.f"
		    a[i__4].r = a[i__5].r, a[i__4].i = a[i__5].i;
#line 180 "claswp.f"
		    i__4 = ip + k * a_dim1;
#line 180 "claswp.f"
		    a[i__4].r = temp.r, a[i__4].i = temp.i;
#line 181 "claswp.f"
/* L40: */
#line 181 "claswp.f"
		}
#line 182 "claswp.f"
	    }
#line 183 "claswp.f"
	    ix += *incx;
#line 184 "claswp.f"
/* L50: */
#line 184 "claswp.f"
	}
#line 185 "claswp.f"
    }

#line 187 "claswp.f"
    return 0;

/*     End of CLASWP */

} /* claswp_ */

