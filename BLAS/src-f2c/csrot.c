#line 1 "csrot.f"
/* csrot.f -- translated by f2c (version 20100827).
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

#line 1 "csrot.f"
/* > \brief \b CSROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER           INCX, INCY, N */
/*       REAL              C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX           CX( * ), CY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSROT applies a plane rotation, where the cos and sin (c and s) are real */
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
/* >          CX is COMPLEX array, dimension at least */
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
/* >          CY is COMPLEX array, dimension at least */
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
/* >          C is REAL */
/* >           On entry, C specifies the cosine, cos. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is REAL */
/* >           On entry, S specifies the sine, sin. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int csrot_(integer *n, doublecomplex *cx, integer *incx, 
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

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 122 "csrot.f"
    /* Parameter adjustments */
#line 122 "csrot.f"
    --cy;
#line 122 "csrot.f"
    --cx;
#line 122 "csrot.f"

#line 122 "csrot.f"
    /* Function Body */
#line 122 "csrot.f"
    if (*n <= 0) {
#line 122 "csrot.f"
	return 0;
#line 122 "csrot.f"
    }
#line 124 "csrot.f"
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

#line 128 "csrot.f"
	i__1 = *n;
#line 128 "csrot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 129 "csrot.f"
	    i__2 = i__;
#line 129 "csrot.f"
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 129 "csrot.f"
	    i__3 = i__;
#line 129 "csrot.f"
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
#line 129 "csrot.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 129 "csrot.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 130 "csrot.f"
	    i__2 = i__;
#line 130 "csrot.f"
	    i__3 = i__;
#line 130 "csrot.f"
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 130 "csrot.f"
	    i__4 = i__;
#line 130 "csrot.f"
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
#line 130 "csrot.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 130 "csrot.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 131 "csrot.f"
	    i__2 = i__;
#line 131 "csrot.f"
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 132 "csrot.f"
	}
#line 133 "csrot.f"
    } else {

/*        code for unequal increments or equal increments not equal */
/*          to 1 */

#line 138 "csrot.f"
	ix = 1;
#line 139 "csrot.f"
	iy = 1;
#line 140 "csrot.f"
	if (*incx < 0) {
#line 140 "csrot.f"
	    ix = (-(*n) + 1) * *incx + 1;
#line 140 "csrot.f"
	}
#line 142 "csrot.f"
	if (*incy < 0) {
#line 142 "csrot.f"
	    iy = (-(*n) + 1) * *incy + 1;
#line 142 "csrot.f"
	}
#line 144 "csrot.f"
	i__1 = *n;
#line 144 "csrot.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 145 "csrot.f"
	    i__2 = ix;
#line 145 "csrot.f"
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 145 "csrot.f"
	    i__3 = iy;
#line 145 "csrot.f"
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
#line 145 "csrot.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 145 "csrot.f"
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 146 "csrot.f"
	    i__2 = iy;
#line 146 "csrot.f"
	    i__3 = iy;
#line 146 "csrot.f"
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 146 "csrot.f"
	    i__4 = ix;
#line 146 "csrot.f"
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
#line 146 "csrot.f"
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 146 "csrot.f"
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 147 "csrot.f"
	    i__2 = ix;
#line 147 "csrot.f"
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 148 "csrot.f"
	    ix += *incx;
#line 149 "csrot.f"
	    iy += *incy;
#line 150 "csrot.f"
	}
#line 151 "csrot.f"
    }
#line 152 "csrot.f"
    return 0;
} /* csrot_ */

