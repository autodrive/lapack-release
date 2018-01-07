#line 1 "icamax.f"
/* icamax.f -- translated by f2c (version 20100827).
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

#line 1 "icamax.f"
/* > \brief \b ICAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ICAMAX(N,CX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX CX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)| */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

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
integer icamax_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal smax;
    extern doublereal scabs1_(doublecomplex *);


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
#line 76 "icamax.f"
    /* Parameter adjustments */
#line 76 "icamax.f"
    --cx;
#line 76 "icamax.f"

#line 76 "icamax.f"
    /* Function Body */
#line 76 "icamax.f"
    ret_val = 0;
#line 77 "icamax.f"
    if (*n < 1 || *incx <= 0) {
#line 77 "icamax.f"
	return ret_val;
#line 77 "icamax.f"
    }
#line 78 "icamax.f"
    ret_val = 1;
#line 79 "icamax.f"
    if (*n == 1) {
#line 79 "icamax.f"
	return ret_val;
#line 79 "icamax.f"
    }
#line 80 "icamax.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 84 "icamax.f"
	smax = scabs1_(&cx[1]);
#line 85 "icamax.f"
	i__1 = *n;
#line 85 "icamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 86 "icamax.f"
	    if (scabs1_(&cx[i__]) > smax) {
#line 87 "icamax.f"
		ret_val = i__;
#line 88 "icamax.f"
		smax = scabs1_(&cx[i__]);
#line 89 "icamax.f"
	    }
#line 90 "icamax.f"
	}
#line 91 "icamax.f"
    } else {

/*        code for increment not equal to 1 */

#line 95 "icamax.f"
	ix = 1;
#line 96 "icamax.f"
	smax = scabs1_(&cx[1]);
#line 97 "icamax.f"
	ix += *incx;
#line 98 "icamax.f"
	i__1 = *n;
#line 98 "icamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 99 "icamax.f"
	    if (scabs1_(&cx[ix]) > smax) {
#line 100 "icamax.f"
		ret_val = i__;
#line 101 "icamax.f"
		smax = scabs1_(&cx[ix]);
#line 102 "icamax.f"
	    }
#line 103 "icamax.f"
	    ix += *incx;
#line 104 "icamax.f"
	}
#line 105 "icamax.f"
    }
#line 106 "icamax.f"
    return ret_val;
} /* icamax_ */

