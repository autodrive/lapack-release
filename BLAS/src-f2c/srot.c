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

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
/* > \endverbatim */
/* > */
/* > \param[in,out] SY */
/* > \verbatim */
/* >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of SY */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

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
#line 114 "srot.f"
    /* Parameter adjustments */
#line 114 "srot.f"
    --sy;
#line 114 "srot.f"
    --sx;
#line 114 "srot.f"

#line 114 "srot.f"
    /* Function Body */
#line 114 "srot.f"
    if (*n <= 0) {
#line 114 "srot.f"
	return 0;
#line 114 "srot.f"
    }
#line 115 "srot.f"
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */

#line 119 "srot.f"
	i__1 = *n;
#line 119 "srot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 120 "srot.f"
	    stemp = *c__ * sx[i__] + *s * sy[i__];
#line 121 "srot.f"
	    sy[i__] = *c__ * sy[i__] - *s * sx[i__];
#line 122 "srot.f"
	    sx[i__] = stemp;
#line 123 "srot.f"
	}
#line 124 "srot.f"
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

#line 129 "srot.f"
	ix = 1;
#line 130 "srot.f"
	iy = 1;
#line 131 "srot.f"
	if (*incx < 0) {
#line 131 "srot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 131 "srot.f"
	}
#line 132 "srot.f"
	if (*incy < 0) {
#line 132 "srot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 132 "srot.f"
	}
#line 133 "srot.f"
	i__1 = *n;
#line 133 "srot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 134 "srot.f"
	    stemp = *c__ * sx[ix] + *s * sy[iy];
#line 135 "srot.f"
	    sy[iy] = *c__ * sy[iy] - *s * sx[ix];
#line 136 "srot.f"
	    sx[ix] = stemp;
#line 137 "srot.f"
	    ix += *incx;
#line 138 "srot.f"
	    iy += *incy;
#line 139 "srot.f"
	}
#line 140 "srot.f"
    }
#line 141 "srot.f"
    return 0;
} /* srot_ */

