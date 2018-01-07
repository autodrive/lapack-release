#line 1 "sasum.f"
/* sasum.f -- translated by f2c (version 20100827).
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

#line 1 "sasum.f"
/* > \brief \b SASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SASUM(N,SX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SASUM takes the sum of the absolute values. */
/* >    uses unrolled loops for increment equal to one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal sasum_(integer *n, doublereal *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, m, mp1, nincx;
    static doublereal stemp;


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 96 "sasum.f"
    /* Parameter adjustments */
#line 96 "sasum.f"
    --sx;
#line 96 "sasum.f"

#line 96 "sasum.f"
    /* Function Body */
#line 96 "sasum.f"
    ret_val = 0.;
#line 97 "sasum.f"
    stemp = 0.;
#line 98 "sasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 98 "sasum.f"
	return ret_val;
#line 98 "sasum.f"
    }
#line 99 "sasum.f"
    if (*incx == 1) {
/*        code for increment equal to 1 */


/*        clean-up loop */

#line 105 "sasum.f"
	m = *n % 6;
#line 106 "sasum.f"
	if (m != 0) {
#line 107 "sasum.f"
	    i__1 = m;
#line 107 "sasum.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 108 "sasum.f"
		stemp += (d__1 = sx[i__], abs(d__1));
#line 109 "sasum.f"
	    }
#line 110 "sasum.f"
	    if (*n < 6) {
#line 111 "sasum.f"
		ret_val = stemp;
#line 112 "sasum.f"
		return ret_val;
#line 113 "sasum.f"
	    }
#line 114 "sasum.f"
	}
#line 115 "sasum.f"
	mp1 = m + 1;
#line 116 "sasum.f"
	i__1 = *n;
#line 116 "sasum.f"
	for (i__ = mp1; i__ <= i__1; i__ += 6) {
#line 117 "sasum.f"
	    stemp = stemp + (d__1 = sx[i__], abs(d__1)) + (d__2 = sx[i__ + 1],
		     abs(d__2)) + (d__3 = sx[i__ + 2], abs(d__3)) + (d__4 = 
		    sx[i__ + 3], abs(d__4)) + (d__5 = sx[i__ + 4], abs(d__5)) 
		    + (d__6 = sx[i__ + 5], abs(d__6));
#line 120 "sasum.f"
	}
#line 121 "sasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 125 "sasum.f"
	nincx = *n * *incx;
#line 126 "sasum.f"
	i__1 = nincx;
#line 126 "sasum.f"
	i__2 = *incx;
#line 126 "sasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 127 "sasum.f"
	    stemp += (d__1 = sx[i__], abs(d__1));
#line 128 "sasum.f"
	}
#line 129 "sasum.f"
    }
#line 130 "sasum.f"
    ret_val = stemp;
#line 131 "sasum.f"
    return ret_val;
} /* sasum_ */

