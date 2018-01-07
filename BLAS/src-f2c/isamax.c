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
/* >    ISAMAX finds the index of the first element having maximum absolute value. */
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
#line 95 "isamax.f"
    /* Parameter adjustments */
#line 95 "isamax.f"
    --sx;
#line 95 "isamax.f"

#line 95 "isamax.f"
    /* Function Body */
#line 95 "isamax.f"
    ret_val = 0;
#line 96 "isamax.f"
    if (*n < 1 || *incx <= 0) {
#line 96 "isamax.f"
	return ret_val;
#line 96 "isamax.f"
    }
#line 97 "isamax.f"
    ret_val = 1;
#line 98 "isamax.f"
    if (*n == 1) {
#line 98 "isamax.f"
	return ret_val;
#line 98 "isamax.f"
    }
#line 99 "isamax.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 103 "isamax.f"
	smax = abs(sx[1]);
#line 104 "isamax.f"
	i__1 = *n;
#line 104 "isamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 105 "isamax.f"
	    if ((d__1 = sx[i__], abs(d__1)) > smax) {
#line 106 "isamax.f"
		ret_val = i__;
#line 107 "isamax.f"
		smax = (d__1 = sx[i__], abs(d__1));
#line 108 "isamax.f"
	    }
#line 109 "isamax.f"
	}
#line 110 "isamax.f"
    } else {

/*        code for increment not equal to 1 */

#line 114 "isamax.f"
	ix = 1;
#line 115 "isamax.f"
	smax = abs(sx[1]);
#line 116 "isamax.f"
	ix += *incx;
#line 117 "isamax.f"
	i__1 = *n;
#line 117 "isamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 118 "isamax.f"
	    if ((d__1 = sx[ix], abs(d__1)) > smax) {
#line 119 "isamax.f"
		ret_val = i__;
#line 120 "isamax.f"
		smax = (d__1 = sx[ix], abs(d__1));
#line 121 "isamax.f"
	    }
#line 122 "isamax.f"
	    ix += *incx;
#line 123 "isamax.f"
	}
#line 124 "isamax.f"
    }
#line 125 "isamax.f"
    return ret_val;
} /* isamax_ */

