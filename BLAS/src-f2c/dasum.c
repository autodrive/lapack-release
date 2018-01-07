#line 1 "dasum.f"
/* dasum.f -- translated by f2c (version 20100827).
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

#line 1 "dasum.f"
/* > \brief \b DASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DASUM takes the sum of the absolute values. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup double_blas_level1 */

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
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, m, mp1;
    static doublereal dtemp;
    static integer nincx;


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 75 "dasum.f"
    /* Parameter adjustments */
#line 75 "dasum.f"
    --dx;
#line 75 "dasum.f"

#line 75 "dasum.f"
    /* Function Body */
#line 75 "dasum.f"
    ret_val = 0.;
#line 76 "dasum.f"
    dtemp = 0.;
#line 77 "dasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 77 "dasum.f"
	return ret_val;
#line 77 "dasum.f"
    }
#line 78 "dasum.f"
    if (*incx == 1) {
/*        code for increment equal to 1 */


/*        clean-up loop */

#line 84 "dasum.f"
	m = *n % 6;
#line 85 "dasum.f"
	if (m != 0) {
#line 86 "dasum.f"
	    i__1 = m;
#line 86 "dasum.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 87 "dasum.f"
		dtemp += (d__1 = dx[i__], abs(d__1));
#line 88 "dasum.f"
	    }
#line 89 "dasum.f"
	    if (*n < 6) {
#line 90 "dasum.f"
		ret_val = dtemp;
#line 91 "dasum.f"
		return ret_val;
#line 92 "dasum.f"
	    }
#line 93 "dasum.f"
	}
#line 94 "dasum.f"
	mp1 = m + 1;
#line 95 "dasum.f"
	i__1 = *n;
#line 95 "dasum.f"
	for (i__ = mp1; i__ <= i__1; i__ += 6) {
#line 96 "dasum.f"
	    dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1],
		     abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = 
		    dx[i__ + 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) 
		    + (d__6 = dx[i__ + 5], abs(d__6));
#line 99 "dasum.f"
	}
#line 100 "dasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 104 "dasum.f"
	nincx = *n * *incx;
#line 105 "dasum.f"
	i__1 = nincx;
#line 105 "dasum.f"
	i__2 = *incx;
#line 105 "dasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 106 "dasum.f"
	    dtemp += (d__1 = dx[i__], abs(d__1));
#line 107 "dasum.f"
	}
#line 108 "dasum.f"
    }
#line 109 "dasum.f"
    ret_val = dtemp;
#line 110 "dasum.f"
    return ret_val;
} /* dasum_ */

