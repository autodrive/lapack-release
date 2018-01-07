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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
#line 76 "sasum.f"
    /* Parameter adjustments */
#line 76 "sasum.f"
    --sx;
#line 76 "sasum.f"

#line 76 "sasum.f"
    /* Function Body */
#line 76 "sasum.f"
    ret_val = 0.;
#line 77 "sasum.f"
    stemp = 0.;
#line 78 "sasum.f"
    if (*n <= 0 || *incx <= 0) {
#line 78 "sasum.f"
	return ret_val;
#line 78 "sasum.f"
    }
#line 79 "sasum.f"
    if (*incx == 1) {
/*        code for increment equal to 1 */


/*        clean-up loop */

#line 85 "sasum.f"
	m = *n % 6;
#line 86 "sasum.f"
	if (m != 0) {
#line 87 "sasum.f"
	    i__1 = m;
#line 87 "sasum.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 88 "sasum.f"
		stemp += (d__1 = sx[i__], abs(d__1));
#line 89 "sasum.f"
	    }
#line 90 "sasum.f"
	    if (*n < 6) {
#line 91 "sasum.f"
		ret_val = stemp;
#line 92 "sasum.f"
		return ret_val;
#line 93 "sasum.f"
	    }
#line 94 "sasum.f"
	}
#line 95 "sasum.f"
	mp1 = m + 1;
#line 96 "sasum.f"
	i__1 = *n;
#line 96 "sasum.f"
	for (i__ = mp1; i__ <= i__1; i__ += 6) {
#line 97 "sasum.f"
	    stemp = stemp + (d__1 = sx[i__], abs(d__1)) + (d__2 = sx[i__ + 1],
		     abs(d__2)) + (d__3 = sx[i__ + 2], abs(d__3)) + (d__4 = 
		    sx[i__ + 3], abs(d__4)) + (d__5 = sx[i__ + 4], abs(d__5)) 
		    + (d__6 = sx[i__ + 5], abs(d__6));
#line 100 "sasum.f"
	}
#line 101 "sasum.f"
    } else {

/*        code for increment not equal to 1 */

#line 105 "sasum.f"
	nincx = *n * *incx;
#line 106 "sasum.f"
	i__1 = nincx;
#line 106 "sasum.f"
	i__2 = *incx;
#line 106 "sasum.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 107 "sasum.f"
	    stemp += (d__1 = sx[i__], abs(d__1));
#line 108 "sasum.f"
	}
#line 109 "sasum.f"
    }
#line 110 "sasum.f"
    ret_val = stemp;
#line 111 "sasum.f"
    return ret_val;
} /* sasum_ */

