#line 1 "zdrot.f"
/* zdrot.f -- translated by f2c (version 20100827).
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

#line 1 "zdrot.f"
/* > \brief \b ZDROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       DOUBLE PRECISION   C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         CX( * ), CY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Applies a plane rotation, where the cos and sin (c and s) are real */
/* > and the vectors cx and cy are complex. */
/* > jack dongarra, linpack, 3/11/78. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the vectors cx and cy. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* >          CX is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array CX must contain the n */
/* >           element vector cx. On exit, CX is overwritten by the updated */
/* >           vector cx. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           CX. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* >          CY is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array CY must contain the n */
/* >           element vector cy. On exit, CY is overwritten by the updated */
/* >           vector cy. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           CY. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >           On entry, C specifies the cosine, cos. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >           On entry, S specifies the sine, sin. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int zdrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*  -- Reference BLAS level1 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 122 "zdrot.f"
    /* Parameter adjustments */
#line 122 "zdrot.f"
    --cy;
#line 122 "zdrot.f"
    --cx;
#line 122 "zdrot.f"

#line 122 "zdrot.f"
    /* Function Body */
#line 122 "zdrot.f"
    if (*n <= 0) {
#line 122 "zdrot.f"
	return 0;
#line 122 "zdrot.f"
    }
#line 124 "zdrot.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 128 "zdrot.f"
	i__1 = *n;
#line 128 "zdrot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 129 "zdrot.f"
	    i__2 = i__;
#line 129 "zdrot.f"
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 129 "zdrot.f"
	    i__3 = i__;
#line 129 "zdrot.f"
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
#line 129 "zdrot.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 129 "zdrot.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 130 "zdrot.f"
	    i__2 = i__;
#line 130 "zdrot.f"
	    i__3 = i__;
#line 130 "zdrot.f"
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 130 "zdrot.f"
	    i__4 = i__;
#line 130 "zdrot.f"
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
#line 130 "zdrot.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 130 "zdrot.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 131 "zdrot.f"
	    i__2 = i__;
#line 131 "zdrot.f"
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 132 "zdrot.f"
	}
#line 133 "zdrot.f"
    } else {

/*        code for unequal increments or equal increments not equal */
/*          to 1 */

#line 138 "zdrot.f"
	ix = 1;
#line 139 "zdrot.f"
	iy = 1;
#line 140 "zdrot.f"
	if (*incx < 0) {
#line 140 "zdrot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 140 "zdrot.f"
	}
#line 142 "zdrot.f"
	if (*incy < 0) {
#line 142 "zdrot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 142 "zdrot.f"
	}
#line 144 "zdrot.f"
	i__1 = *n;
#line 144 "zdrot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 145 "zdrot.f"
	    i__2 = ix;
#line 145 "zdrot.f"
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 145 "zdrot.f"
	    i__3 = iy;
#line 145 "zdrot.f"
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
#line 145 "zdrot.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 145 "zdrot.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 146 "zdrot.f"
	    i__2 = iy;
#line 146 "zdrot.f"
	    i__3 = iy;
#line 146 "zdrot.f"
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 146 "zdrot.f"
	    i__4 = ix;
#line 146 "zdrot.f"
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
#line 146 "zdrot.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 146 "zdrot.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 147 "zdrot.f"
	    i__2 = ix;
#line 147 "zdrot.f"
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 148 "zdrot.f"
	    ix += *incx;
#line 149 "zdrot.f"
	    iy += *incy;
#line 150 "zdrot.f"
	}
#line 151 "zdrot.f"
    }
#line 152 "zdrot.f"
    return 0;
} /* zdrot_ */

