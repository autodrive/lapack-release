#line 1 "srot.f"
/* srot.f -- translated by f2c (version 20100827).
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

#line 1 "srot.f"
/* > \brief \b SROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S) */

/*       .. Scalar Arguments .. */
/*       REAL C,S */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    applies a plane rotation. */
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
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int srot_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
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
#line 73 "srot.f"
    /* Parameter adjustments */
#line 73 "srot.f"
    --sy;
#line 73 "srot.f"
    --sx;
#line 73 "srot.f"

#line 73 "srot.f"
    /* Function Body */
#line 73 "srot.f"
    if (*n <= 0) {
#line 73 "srot.f"
	return 0;
#line 73 "srot.f"
    }
#line 74 "srot.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */

#line 78 "srot.f"
	i__1 = *n;
#line 78 "srot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 79 "srot.f"
	    stemp = *c__ * sx[i__] + *s * sy[i__];
#line 80 "srot.f"
	    sy[i__] = *c__ * sy[i__] - *s * sx[i__];
#line 81 "srot.f"
	    sx[i__] = stemp;
#line 82 "srot.f"
	}
#line 83 "srot.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 88 "srot.f"
	ix = 1;
#line 89 "srot.f"
	iy = 1;
#line 90 "srot.f"
	if (*incx < 0) {
#line 90 "srot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 90 "srot.f"
	}
#line 91 "srot.f"
	if (*incy < 0) {
#line 91 "srot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 91 "srot.f"
	}
#line 92 "srot.f"
	i__1 = *n;
#line 92 "srot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 93 "srot.f"
	    stemp = *c__ * sx[ix] + *s * sy[iy];
#line 94 "srot.f"
	    sy[iy] = *c__ * sy[iy] - *s * sx[ix];
#line 95 "srot.f"
	    sx[ix] = stemp;
#line 96 "srot.f"
	    ix += *incx;
#line 97 "srot.f"
	    iy += *incy;
#line 98 "srot.f"
	}
#line 99 "srot.f"
    }
#line 100 "srot.f"
    return 0;
} /* srot_ */

