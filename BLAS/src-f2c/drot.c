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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] DX */
/* > \verbatim */
/* >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of DX */
/* > \endverbatim */
/* > */
/* > \param[in,out] DY */
/* > \verbatim */
/* >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of DY */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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
#line 114 "drot.f"
    /* Parameter adjustments */
#line 114 "drot.f"
    --dy;
#line 114 "drot.f"
    --dx;
#line 114 "drot.f"

#line 114 "drot.f"
    /* Function Body */
#line 114 "drot.f"
    if (*n <= 0) {
#line 114 "drot.f"
	return 0;
#line 114 "drot.f"
    }
#line 115 "drot.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */

#line 119 "drot.f"
	i__1 = *n;
#line 119 "drot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 120 "drot.f"
	    dtemp = *c__ * dx[i__] + *s * dy[i__];
#line 121 "drot.f"
	    dy[i__] = *c__ * dy[i__] - *s * dx[i__];
#line 122 "drot.f"
	    dx[i__] = dtemp;
#line 123 "drot.f"
	}
#line 124 "drot.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 129 "drot.f"
	ix = 1;
#line 130 "drot.f"
	iy = 1;
#line 131 "drot.f"
	if (*incx < 0) {
#line 131 "drot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 131 "drot.f"
	}
#line 132 "drot.f"
	if (*incy < 0) {
#line 132 "drot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 132 "drot.f"
	}
#line 133 "drot.f"
	i__1 = *n;
#line 133 "drot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 134 "drot.f"
	    dtemp = *c__ * dx[ix] + *s * dy[iy];
#line 135 "drot.f"
	    dy[iy] = *c__ * dy[iy] - *s * dx[ix];
#line 136 "drot.f"
	    dx[ix] = dtemp;
#line 137 "drot.f"
	    ix += *incx;
#line 138 "drot.f"
	    iy += *incy;
#line 139 "drot.f"
	}
#line 140 "drot.f"
    }
#line 141 "drot.f"
    return 0;
} /* drot_ */

