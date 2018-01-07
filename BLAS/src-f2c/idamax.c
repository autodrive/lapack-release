#line 1 "idamax.f"
/* idamax.f -- translated by f2c (version 20100827).
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

#line 1 "idamax.f"
/* > \brief \b IDAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IDAMAX(N,DX,INCX) */

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
/* >    IDAMAX finds the index of the first element having maximum absolute value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup aux_blas */

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
integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__;


/*  -- Reference BLAS level1 routine (version 3.6.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
#line 75 "idamax.f"
    /* Parameter adjustments */
#line 75 "idamax.f"
    --dx;
#line 75 "idamax.f"

#line 75 "idamax.f"
    /* Function Body */
#line 75 "idamax.f"
    ret_val = 0;
#line 76 "idamax.f"
    if (*n < 1 || *incx <= 0) {
#line 76 "idamax.f"
	return ret_val;
#line 76 "idamax.f"
    }
#line 77 "idamax.f"
    ret_val = 1;
#line 78 "idamax.f"
    if (*n == 1) {
#line 78 "idamax.f"
	return ret_val;
#line 78 "idamax.f"
    }
#line 79 "idamax.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 83 "idamax.f"
	dmax__ = abs(dx[1]);
#line 84 "idamax.f"
	i__1 = *n;
#line 84 "idamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 85 "idamax.f"
	    if ((d__1 = dx[i__], abs(d__1)) > dmax__) {
#line 86 "idamax.f"
		ret_val = i__;
#line 87 "idamax.f"
		dmax__ = (d__1 = dx[i__], abs(d__1));
#line 88 "idamax.f"
	    }
#line 89 "idamax.f"
	}
#line 90 "idamax.f"
    } else {

/*        code for increment not equal to 1 */

#line 94 "idamax.f"
	ix = 1;
#line 95 "idamax.f"
	dmax__ = abs(dx[1]);
#line 96 "idamax.f"
	ix += *incx;
#line 97 "idamax.f"
	i__1 = *n;
#line 97 "idamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 98 "idamax.f"
	    if ((d__1 = dx[ix], abs(d__1)) > dmax__) {
#line 99 "idamax.f"
		ret_val = i__;
#line 100 "idamax.f"
		dmax__ = (d__1 = dx[ix], abs(d__1));
#line 101 "idamax.f"
	    }
#line 102 "idamax.f"
	    ix += *incx;
#line 103 "idamax.f"
	}
#line 104 "idamax.f"
    }
#line 105 "idamax.f"
    return ret_val;
} /* idamax_ */

