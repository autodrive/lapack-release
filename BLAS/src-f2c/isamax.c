#line 1 "isamax.f"
/* isamax.f -- translated by f2c (version 20100827).
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

#line 1 "isamax.f"
/* > \brief \b ISAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ISAMAX(N,SX,INCX) */

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
/* >    ISAMAX finds the index of element having max. absolute value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

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
integer isamax_(integer *n, doublereal *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal smax;


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
#line 75 "isamax.f"
    /* Parameter adjustments */
#line 75 "isamax.f"
    --sx;
#line 75 "isamax.f"

#line 75 "isamax.f"
    /* Function Body */
#line 75 "isamax.f"
    ret_val = 0;
#line 76 "isamax.f"
    if (*n < 1 || *incx <= 0) {
#line 76 "isamax.f"
	return ret_val;
#line 76 "isamax.f"
    }
#line 77 "isamax.f"
    ret_val = 1;
#line 78 "isamax.f"
    if (*n == 1) {
#line 78 "isamax.f"
	return ret_val;
#line 78 "isamax.f"
    }
#line 79 "isamax.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 83 "isamax.f"
	smax = abs(sx[1]);
#line 84 "isamax.f"
	i__1 = *n;
#line 84 "isamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 85 "isamax.f"
	    if ((d__1 = sx[i__], abs(d__1)) > smax) {
#line 86 "isamax.f"
		ret_val = i__;
#line 87 "isamax.f"
		smax = (d__1 = sx[i__], abs(d__1));
#line 88 "isamax.f"
	    }
#line 89 "isamax.f"
	}
#line 90 "isamax.f"
    } else {

/*        code for increment not equal to 1 */

#line 94 "isamax.f"
	ix = 1;
#line 95 "isamax.f"
	smax = abs(sx[1]);
#line 96 "isamax.f"
	ix += *incx;
#line 97 "isamax.f"
	i__1 = *n;
#line 97 "isamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 98 "isamax.f"
	    if ((d__1 = sx[ix], abs(d__1)) > smax) {
#line 99 "isamax.f"
		ret_val = i__;
#line 100 "isamax.f"
		smax = (d__1 = sx[ix], abs(d__1));
#line 101 "isamax.f"
	    }
#line 102 "isamax.f"
	    ix += *incx;
#line 103 "isamax.f"
	}
#line 104 "isamax.f"
    }
#line 105 "isamax.f"
    return ret_val;
} /* isamax_ */

