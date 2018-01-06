#line 1 "zlaswp.f"
/* zlaswp.f -- translated by f2c (version 20100827).
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

#line 1 "zlaswp.f"
/* > \brief \b ZLASWP performs a series of row interchanges on a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLASWP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaswp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaswp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaswp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, K1, K2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASWP performs a series of row interchanges on the matrix A. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed. */
/* >          IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be */
/* >          interchanged. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of IPIV. If INCX */
/* >          is negative, the pivots are applied in reverse order. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup complex16OTHERauxiliary */

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
/* Subroutine */ int zlaswp_(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    static doublecomplex temp;


/*  -- LAPACK auxiliary routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows */
/*     K1 through K2. */

#line 142 "zlaswp.f"
    /* Parameter adjustments */
#line 142 "zlaswp.f"
    a_dim1 = *lda;
#line 142 "zlaswp.f"
    a_offset = 1 + a_dim1;
#line 142 "zlaswp.f"
    a -= a_offset;
#line 142 "zlaswp.f"
    --ipiv;
#line 142 "zlaswp.f"

#line 142 "zlaswp.f"
    /* Function Body */
#line 142 "zlaswp.f"
    if (*incx > 0) {
#line 143 "zlaswp.f"
	ix0 = *k1;
#line 144 "zlaswp.f"
	i1 = *k1;
#line 145 "zlaswp.f"
	i2 = *k2;
#line 146 "zlaswp.f"
	inc = 1;
#line 147 "zlaswp.f"
    } else if (*incx < 0) {
#line 148 "zlaswp.f"
	ix0 = *k1 + (*k1 - *k2) * *incx;
#line 149 "zlaswp.f"
	i1 = *k2;
#line 150 "zlaswp.f"
	i2 = *k1;
#line 151 "zlaswp.f"
	inc = -1;
#line 152 "zlaswp.f"
    } else {
#line 153 "zlaswp.f"
	return 0;
#line 154 "zlaswp.f"
    }

#line 156 "zlaswp.f"
    n32 = *n / 32 << 5;
#line 157 "zlaswp.f"
    if (n32 != 0) {
#line 158 "zlaswp.f"
	i__1 = n32;
#line 158 "zlaswp.f"
	for (j = 1; j <= i__1; j += 32) {
#line 159 "zlaswp.f"
	    ix = ix0;
#line 160 "zlaswp.f"
	    i__2 = i2;
#line 160 "zlaswp.f"
	    i__3 = inc;
#line 160 "zlaswp.f"
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
#line 161 "zlaswp.f"
		ip = ipiv[ix];
#line 162 "zlaswp.f"
		if (ip != i__) {
#line 163 "zlaswp.f"
		    i__4 = j + 31;
#line 163 "zlaswp.f"
		    for (k = j; k <= i__4; ++k) {
#line 164 "zlaswp.f"
			i__5 = i__ + k * a_dim1;
#line 164 "zlaswp.f"
			temp.r = a[i__5].r, temp.i = a[i__5].i;
#line 165 "zlaswp.f"
			i__5 = i__ + k * a_dim1;
#line 165 "zlaswp.f"
			i__6 = ip + k * a_dim1;
#line 165 "zlaswp.f"
			a[i__5].r = a[i__6].r, a[i__5].i = a[i__6].i;
#line 166 "zlaswp.f"
			i__5 = ip + k * a_dim1;
#line 166 "zlaswp.f"
			a[i__5].r = temp.r, a[i__5].i = temp.i;
#line 167 "zlaswp.f"
/* L10: */
#line 167 "zlaswp.f"
		    }
#line 168 "zlaswp.f"
		}
#line 169 "zlaswp.f"
		ix += *incx;
#line 170 "zlaswp.f"
/* L20: */
#line 170 "zlaswp.f"
	    }
#line 171 "zlaswp.f"
/* L30: */
#line 171 "zlaswp.f"
	}
#line 172 "zlaswp.f"
    }
#line 173 "zlaswp.f"
    if (n32 != *n) {
#line 174 "zlaswp.f"
	++n32;
#line 175 "zlaswp.f"
	ix = ix0;
#line 176 "zlaswp.f"
	i__1 = i2;
#line 176 "zlaswp.f"
	i__3 = inc;
#line 176 "zlaswp.f"
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
#line 177 "zlaswp.f"
	    ip = ipiv[ix];
#line 178 "zlaswp.f"
	    if (ip != i__) {
#line 179 "zlaswp.f"
		i__2 = *n;
#line 179 "zlaswp.f"
		for (k = n32; k <= i__2; ++k) {
#line 180 "zlaswp.f"
		    i__4 = i__ + k * a_dim1;
#line 180 "zlaswp.f"
		    temp.r = a[i__4].r, temp.i = a[i__4].i;
#line 181 "zlaswp.f"
		    i__4 = i__ + k * a_dim1;
#line 181 "zlaswp.f"
		    i__5 = ip + k * a_dim1;
#line 181 "zlaswp.f"
		    a[i__4].r = a[i__5].r, a[i__4].i = a[i__5].i;
#line 182 "zlaswp.f"
		    i__4 = ip + k * a_dim1;
#line 182 "zlaswp.f"
		    a[i__4].r = temp.r, a[i__4].i = temp.i;
#line 183 "zlaswp.f"
/* L40: */
#line 183 "zlaswp.f"
		}
#line 184 "zlaswp.f"
	    }
#line 185 "zlaswp.f"
	    ix += *incx;
#line 186 "zlaswp.f"
/* L50: */
#line 186 "zlaswp.f"
	}
#line 187 "zlaswp.f"
    }

#line 189 "zlaswp.f"
    return 0;

/*     End of ZLASWP */

} /* zlaswp_ */

