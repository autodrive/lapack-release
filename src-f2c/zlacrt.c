#line 1 "zlacrt.f"
/* zlacrt.f -- translated by f2c (version 20100827).
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

#line 1 "zlacrt.f"
/* > \brief \b ZLACRT performs a linear transformation of a pair of complex vectors. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLACRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       COMPLEX*16         C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         CX( * ), CY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACRT performs the operation */
/* > */
/* >    (  c  s )( x )  ==> ( x ) */
/* >    ( -s  c )( y )      ( y ) */
/* > */
/* > where c and s are complex and the vectors x and y are complex. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of elements in the vectors CX and CY. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* >          CX is COMPLEX*16 array, dimension (N) */
/* >          On input, the vector x. */
/* >          On output, CX is overwritten with c*x + s*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of CX.  INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* >          CY is COMPLEX*16 array, dimension (N) */
/* >          On input, the vector y. */
/* >          On output, CY is overwritten with -s*x + c*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between successive values of CY.  INCY <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is COMPLEX*16 */
/* >          C and S define the matrix */
/* >             [  C   S  ]. */
/* >             [ -S   C  ] */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlacrt_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *c__, doublecomplex *
	s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 129 "zlacrt.f"
    /* Parameter adjustments */
#line 129 "zlacrt.f"
    --cy;
#line 129 "zlacrt.f"
    --cx;
#line 129 "zlacrt.f"

#line 129 "zlacrt.f"
    /* Function Body */
#line 129 "zlacrt.f"
    if (*n <= 0) {
#line 129 "zlacrt.f"
	return 0;
#line 129 "zlacrt.f"
    }
#line 131 "zlacrt.f"
    if (*incx == 1 && *incy == 1) {
#line 131 "zlacrt.f"
	goto L20;
#line 131 "zlacrt.f"
    }

/*     Code for unequal increments or equal increments not equal to 1 */

#line 136 "zlacrt.f"
    ix = 1;
#line 137 "zlacrt.f"
    iy = 1;
#line 138 "zlacrt.f"
    if (*incx < 0) {
#line 138 "zlacrt.f"
	ix = (-(*n) + 1) * *incx + 1;
#line 138 "zlacrt.f"
    }
#line 140 "zlacrt.f"
    if (*incy < 0) {
#line 140 "zlacrt.f"
	iy = (-(*n) + 1) * *incy + 1;
#line 140 "zlacrt.f"
    }
#line 142 "zlacrt.f"
    i__1 = *n;
#line 142 "zlacrt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 143 "zlacrt.f"
	i__2 = ix;
#line 143 "zlacrt.f"
	z__2.r = c__->r * cx[i__2].r - c__->i * cx[i__2].i, z__2.i = c__->r * 
		cx[i__2].i + c__->i * cx[i__2].r;
#line 143 "zlacrt.f"
	i__3 = iy;
#line 143 "zlacrt.f"
	z__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i, z__3.i = s->r * cy[
		i__3].i + s->i * cy[i__3].r;
#line 143 "zlacrt.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 143 "zlacrt.f"
	ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 144 "zlacrt.f"
	i__2 = iy;
#line 144 "zlacrt.f"
	i__3 = iy;
#line 144 "zlacrt.f"
	z__2.r = c__->r * cy[i__3].r - c__->i * cy[i__3].i, z__2.i = c__->r * 
		cy[i__3].i + c__->i * cy[i__3].r;
#line 144 "zlacrt.f"
	i__4 = ix;
#line 144 "zlacrt.f"
	z__3.r = s->r * cx[i__4].r - s->i * cx[i__4].i, z__3.i = s->r * cx[
		i__4].i + s->i * cx[i__4].r;
#line 144 "zlacrt.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 144 "zlacrt.f"
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 145 "zlacrt.f"
	i__2 = ix;
#line 145 "zlacrt.f"
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 146 "zlacrt.f"
	ix += *incx;
#line 147 "zlacrt.f"
	iy += *incy;
#line 148 "zlacrt.f"
/* L10: */
#line 148 "zlacrt.f"
    }
#line 149 "zlacrt.f"
    return 0;

/*     Code for both increments equal to 1 */

#line 153 "zlacrt.f"
L20:
#line 154 "zlacrt.f"
    i__1 = *n;
#line 154 "zlacrt.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 155 "zlacrt.f"
	i__2 = i__;
#line 155 "zlacrt.f"
	z__2.r = c__->r * cx[i__2].r - c__->i * cx[i__2].i, z__2.i = c__->r * 
		cx[i__2].i + c__->i * cx[i__2].r;
#line 155 "zlacrt.f"
	i__3 = i__;
#line 155 "zlacrt.f"
	z__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i, z__3.i = s->r * cy[
		i__3].i + s->i * cy[i__3].r;
#line 155 "zlacrt.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 155 "zlacrt.f"
	ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 156 "zlacrt.f"
	i__2 = i__;
#line 156 "zlacrt.f"
	i__3 = i__;
#line 156 "zlacrt.f"
	z__2.r = c__->r * cy[i__3].r - c__->i * cy[i__3].i, z__2.i = c__->r * 
		cy[i__3].i + c__->i * cy[i__3].r;
#line 156 "zlacrt.f"
	i__4 = i__;
#line 156 "zlacrt.f"
	z__3.r = s->r * cx[i__4].r - s->i * cx[i__4].i, z__3.i = s->r * cx[
		i__4].i + s->i * cx[i__4].r;
#line 156 "zlacrt.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 156 "zlacrt.f"
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 157 "zlacrt.f"
	i__2 = i__;
#line 157 "zlacrt.f"
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
#line 158 "zlacrt.f"
/* L30: */
#line 158 "zlacrt.f"
    }
#line 159 "zlacrt.f"
    return 0;
} /* zlacrt_ */

