#line 1 "izamax.f"
/* izamax.f -- translated by f2c (version 20100827).
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

#line 1 "izamax.f"
/* > \brief \b IZAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IZAMAX(N,ZX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 ZX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)| */
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
/* >     jack dongarra, 1/15/85. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
integer izamax_(integer *n, doublecomplex *zx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__;
    extern doublereal dcabs1_(doublecomplex *);


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
#line 76 "izamax.f"
    /* Parameter adjustments */
#line 76 "izamax.f"
    --zx;
#line 76 "izamax.f"

#line 76 "izamax.f"
    /* Function Body */
#line 76 "izamax.f"
    ret_val = 0;
#line 77 "izamax.f"
    if (*n < 1 || *incx <= 0) {
#line 77 "izamax.f"
	return ret_val;
#line 77 "izamax.f"
    }
#line 78 "izamax.f"
    ret_val = 1;
#line 79 "izamax.f"
    if (*n == 1) {
#line 79 "izamax.f"
	return ret_val;
#line 79 "izamax.f"
    }
#line 80 "izamax.f"
    if (*incx == 1) {

/*        code for increment equal to 1 */

#line 84 "izamax.f"
	dmax__ = dcabs1_(&zx[1]);
#line 85 "izamax.f"
	i__1 = *n;
#line 85 "izamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 86 "izamax.f"
	    if (dcabs1_(&zx[i__]) > dmax__) {
#line 87 "izamax.f"
		ret_val = i__;
#line 88 "izamax.f"
		dmax__ = dcabs1_(&zx[i__]);
#line 89 "izamax.f"
	    }
#line 90 "izamax.f"
	}
#line 91 "izamax.f"
    } else {

/*        code for increment not equal to 1 */

#line 95 "izamax.f"
	ix = 1;
#line 96 "izamax.f"
	dmax__ = dcabs1_(&zx[1]);
#line 97 "izamax.f"
	ix += *incx;
#line 98 "izamax.f"
	i__1 = *n;
#line 98 "izamax.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 99 "izamax.f"
	    if (dcabs1_(&zx[ix]) > dmax__) {
#line 100 "izamax.f"
		ret_val = i__;
#line 101 "izamax.f"
		dmax__ = dcabs1_(&zx[ix]);
#line 102 "izamax.f"
	    }
#line 103 "izamax.f"
	    ix += *incx;
#line 104 "izamax.f"
	}
#line 105 "izamax.f"
    }
#line 106 "izamax.f"
    return ret_val;
} /* izamax_ */

