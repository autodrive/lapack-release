#line 1 "drot.f"
/* drot.f -- translated by f2c (version 20100827).
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

#line 1 "drot.f"
/* > \brief \b DROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION C,S */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*),DY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DROT applies a plane rotation. */
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
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
    static doublereal dtemp;


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
#line 73 "drot.f"
    /* Parameter adjustments */
#line 73 "drot.f"
    --dy;
#line 73 "drot.f"
    --dx;
#line 73 "drot.f"

#line 73 "drot.f"
    /* Function Body */
#line 73 "drot.f"
    if (*n <= 0) {
#line 73 "drot.f"
	return 0;
#line 73 "drot.f"
    }
#line 74 "drot.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */

#line 78 "drot.f"
	i__1 = *n;
#line 78 "drot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 79 "drot.f"
	    dtemp = *c__ * dx[i__] + *s * dy[i__];
#line 80 "drot.f"
	    dy[i__] = *c__ * dy[i__] - *s * dx[i__];
#line 81 "drot.f"
	    dx[i__] = dtemp;
#line 82 "drot.f"
	}
#line 83 "drot.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 88 "drot.f"
	ix = 1;
#line 89 "drot.f"
	iy = 1;
#line 90 "drot.f"
	if (*incx < 0) {
#line 90 "drot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 90 "drot.f"
	}
#line 91 "drot.f"
	if (*incy < 0) {
#line 91 "drot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 91 "drot.f"
	}
#line 92 "drot.f"
	i__1 = *n;
#line 92 "drot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 93 "drot.f"
	    dtemp = *c__ * dx[ix] + *s * dy[iy];
#line 94 "drot.f"
	    dy[iy] = *c__ * dy[iy] - *s * dx[ix];
#line 95 "drot.f"
	    dx[ix] = dtemp;
#line 96 "drot.f"
	    ix += *incx;
#line 97 "drot.f"
	    iy += *incy;
#line 98 "drot.f"
	}
#line 99 "drot.f"
    }
#line 100 "drot.f"
    return 0;
} /* drot_ */

