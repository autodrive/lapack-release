#line 1 "zrot.f"
/* zrot.f -- translated by f2c (version 20100827).
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

#line 1 "zrot.f"
/* > \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZROT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f"
> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f"
> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f"
> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       DOUBLE PRECISION   C */
/*       COMPLEX*16         S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         CX( * ), CY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZROT   applies a plane rotation, where the cos (C) is real and the */
/* > sin (S) is complex, and the vectors CX and CY are complex. */
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
/* >          On input, the vector X. */
/* >          On output, CX is overwritten with C*X + S*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of CY.  INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* >          CY is COMPLEX*16 array, dimension (N) */
/* >          On input, the vector Y. */
/* >          On output, CY is overwritten with -CONJG(S)*X + C*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >          The increment between successive values of CY.  INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is COMPLEX*16 */
/* >          C and S define a rotation */
/* >             [  C          S  ] */
/* >             [ -conjg(S)   C  ] */
/* >          where C*C + S*CONJG(S) = 1.0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex stemp;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 131 "zrot.f"
    /* Parameter adjustments */
#line 131 "zrot.f"
    --cy;
#line 131 "zrot.f"
    --cx;
#line 131 "zrot.f"

#line 131 "zrot.f"
    /* Function Body */
#line 131 "zrot.f"
    if (*n <= 0) {
#line 131 "zrot.f"
	return 0;
#line 131 "zrot.f"
    }
#line 133 "zrot.f"
    if (*incx == 1 && *incy == 1) {
#line 133 "zrot.f"
	goto L20;
#line 133 "zrot.f"
    }

/*     Code for unequal increments or equal increments not equal to 1 */

#line 138 "zrot.f"
    ix = 1;
#line 139 "zrot.f"
    iy = 1;
#line 140 "zrot.f"
    if (*incx < 0) {
#line 140 "zrot.f"
	ix = (-(*n) + 1) * *incx + 1;
#line 140 "zrot.f"
    }
#line 142 "zrot.f"
    if (*incy < 0) {
#line 142 "zrot.f"
	iy = (-(*n) + 1) * *incy + 1;
#line 142 "zrot.f"
    }
#line 144 "zrot.f"
    i__1 = *n;
#line 144 "zrot.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 145 "zrot.f"
	i__2 = ix;
#line 145 "zrot.f"
	z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 145 "zrot.f"
	i__3 = iy;
#line 145 "zrot.f"
	z__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i, z__3.i = s->r * cy[
		i__3].i + s->i * cy[i__3].r;
#line 145 "zrot.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 145 "zrot.f"
	stemp.r = z__1.r, stemp.i = z__1.i;
#line 146 "zrot.f"
	i__2 = iy;
#line 146 "zrot.f"
	i__3 = iy;
#line 146 "zrot.f"
	z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 146 "zrot.f"
	d_cnjg(&z__4, s);
#line 146 "zrot.f"
	i__4 = ix;
#line 146 "zrot.f"
	z__3.r = z__4.r * cx[i__4].r - z__4.i * cx[i__4].i, z__3.i = z__4.r * 
		cx[i__4].i + z__4.i * cx[i__4].r;
#line 146 "zrot.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 146 "zrot.f"
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 147 "zrot.f"
	i__2 = ix;
#line 147 "zrot.f"
	cx[i__2].r = stemp.r, cx[i__2].i = stemp.i;
#line 148 "zrot.f"
	ix += *incx;
#line 149 "zrot.f"
	iy += *incy;
#line 150 "zrot.f"
/* L10: */
#line 150 "zrot.f"
    }
#line 151 "zrot.f"
    return 0;

/*     Code for both increments equal to 1 */

#line 155 "zrot.f"
L20:
#line 156 "zrot.f"
    i__1 = *n;
#line 156 "zrot.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 157 "zrot.f"
	i__2 = i__;
#line 157 "zrot.f"
	z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
#line 157 "zrot.f"
	i__3 = i__;
#line 157 "zrot.f"
	z__3.r = s->r * cy[i__3].r - s->i * cy[i__3].i, z__3.i = s->r * cy[
		i__3].i + s->i * cy[i__3].r;
#line 157 "zrot.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 157 "zrot.f"
	stemp.r = z__1.r, stemp.i = z__1.i;
#line 158 "zrot.f"
	i__2 = i__;
#line 158 "zrot.f"
	i__3 = i__;
#line 158 "zrot.f"
	z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
#line 158 "zrot.f"
	d_cnjg(&z__4, s);
#line 158 "zrot.f"
	i__4 = i__;
#line 158 "zrot.f"
	z__3.r = z__4.r * cx[i__4].r - z__4.i * cx[i__4].i, z__3.i = z__4.r * 
		cx[i__4].i + z__4.i * cx[i__4].r;
#line 158 "zrot.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 158 "zrot.f"
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
#line 159 "zrot.f"
	i__2 = i__;
#line 159 "zrot.f"
	cx[i__2].r = stemp.r, cx[i__2].i = stemp.i;
#line 160 "zrot.f"
/* L30: */
#line 160 "zrot.f"
    }
#line 161 "zrot.f"
    return 0;
} /* zrot_ */

