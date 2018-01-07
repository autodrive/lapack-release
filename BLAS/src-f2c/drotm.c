#line 1 "drotm.f"
/* drotm.f -- translated by f2c (version 20100827).
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

#line 1 "drotm.f"
/* > \brief \b DROTM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DPARAM(5),DX(*),DY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */
/* > */
/* >    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN */
/* >    (DY**T) */
/* > */
/* >    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/* >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY. */
/* >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */
/* > */
/* >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */
/* >    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM. */
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
/* > \param[in] DPARAM */
/* > \verbatim */
/* >          DPARAM is DOUBLE PRECISION array, dimension (5) */
/* >     DPARAM(1)=DFLAG */
/* >     DPARAM(2)=DH11 */
/* >     DPARAM(3)=DH21 */
/* >     DPARAM(4)=DH12 */
/* >     DPARAM(5)=DH22 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal w, z__;
    static integer kx, ky;
    static doublereal dh11, dh12, dh21, dh22, dflag;
    static integer nsteps;


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
/*     .. Data statements .. */
#line 118 "drotm.f"
    /* Parameter adjustments */
#line 118 "drotm.f"
    --dparam;
#line 118 "drotm.f"
    --dy;
#line 118 "drotm.f"
    --dx;
#line 118 "drotm.f"

#line 118 "drotm.f"
    /* Function Body */
/*     .. */

#line 121 "drotm.f"
    dflag = dparam[1];
#line 122 "drotm.f"
    if (*n <= 0 || dflag + two == zero) {
#line 122 "drotm.f"
	return 0;
#line 122 "drotm.f"
    }
#line 123 "drotm.f"
    if (*incx == *incy && *incx > 0) {

#line 125 "drotm.f"
	nsteps = *n * *incx;
#line 126 "drotm.f"
	if (dflag < zero) {
#line 127 "drotm.f"
	    dh11 = dparam[2];
#line 128 "drotm.f"
	    dh12 = dparam[4];
#line 129 "drotm.f"
	    dh21 = dparam[3];
#line 130 "drotm.f"
	    dh22 = dparam[5];
#line 131 "drotm.f"
	    i__1 = nsteps;
#line 131 "drotm.f"
	    i__2 = *incx;
#line 131 "drotm.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 132 "drotm.f"
		w = dx[i__];
#line 133 "drotm.f"
		z__ = dy[i__];
#line 134 "drotm.f"
		dx[i__] = w * dh11 + z__ * dh12;
#line 135 "drotm.f"
		dy[i__] = w * dh21 + z__ * dh22;
#line 136 "drotm.f"
	    }
#line 137 "drotm.f"
	} else if (dflag == zero) {
#line 138 "drotm.f"
	    dh12 = dparam[4];
#line 139 "drotm.f"
	    dh21 = dparam[3];
#line 140 "drotm.f"
	    i__2 = nsteps;
#line 140 "drotm.f"
	    i__1 = *incx;
#line 140 "drotm.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 141 "drotm.f"
		w = dx[i__];
#line 142 "drotm.f"
		z__ = dy[i__];
#line 143 "drotm.f"
		dx[i__] = w + z__ * dh12;
#line 144 "drotm.f"
		dy[i__] = w * dh21 + z__;
#line 145 "drotm.f"
	    }
#line 146 "drotm.f"
	} else {
#line 147 "drotm.f"
	    dh11 = dparam[2];
#line 148 "drotm.f"
	    dh22 = dparam[5];
#line 149 "drotm.f"
	    i__1 = nsteps;
#line 149 "drotm.f"
	    i__2 = *incx;
#line 149 "drotm.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 150 "drotm.f"
		w = dx[i__];
#line 151 "drotm.f"
		z__ = dy[i__];
#line 152 "drotm.f"
		dx[i__] = w * dh11 + z__;
#line 153 "drotm.f"
		dy[i__] = -w + dh22 * z__;
#line 154 "drotm.f"
	    }
#line 155 "drotm.f"
	}
#line 156 "drotm.f"
    } else {
#line 157 "drotm.f"
	kx = 1;
#line 158 "drotm.f"
	ky = 1;
#line 159 "drotm.f"
	if (*incx < 0) {
#line 159 "drotm.f"
	    kx = (1 - *n) * *incx + 1;
#line 159 "drotm.f"
	}
#line 160 "drotm.f"
	if (*incy < 0) {
#line 160 "drotm.f"
	    ky = (1 - *n) * *incy + 1;
#line 160 "drotm.f"
	}

#line 162 "drotm.f"
	if (dflag < zero) {
#line 163 "drotm.f"
	    dh11 = dparam[2];
#line 164 "drotm.f"
	    dh12 = dparam[4];
#line 165 "drotm.f"
	    dh21 = dparam[3];
#line 166 "drotm.f"
	    dh22 = dparam[5];
#line 167 "drotm.f"
	    i__2 = *n;
#line 167 "drotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 168 "drotm.f"
		w = dx[kx];
#line 169 "drotm.f"
		z__ = dy[ky];
#line 170 "drotm.f"
		dx[kx] = w * dh11 + z__ * dh12;
#line 171 "drotm.f"
		dy[ky] = w * dh21 + z__ * dh22;
#line 172 "drotm.f"
		kx += *incx;
#line 173 "drotm.f"
		ky += *incy;
#line 174 "drotm.f"
	    }
#line 175 "drotm.f"
	} else if (dflag == zero) {
#line 176 "drotm.f"
	    dh12 = dparam[4];
#line 177 "drotm.f"
	    dh21 = dparam[3];
#line 178 "drotm.f"
	    i__2 = *n;
#line 178 "drotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 179 "drotm.f"
		w = dx[kx];
#line 180 "drotm.f"
		z__ = dy[ky];
#line 181 "drotm.f"
		dx[kx] = w + z__ * dh12;
#line 182 "drotm.f"
		dy[ky] = w * dh21 + z__;
#line 183 "drotm.f"
		kx += *incx;
#line 184 "drotm.f"
		ky += *incy;
#line 185 "drotm.f"
	    }
#line 186 "drotm.f"
	} else {
#line 187 "drotm.f"
	    dh11 = dparam[2];
#line 188 "drotm.f"
	    dh22 = dparam[5];
#line 189 "drotm.f"
	    i__2 = *n;
#line 189 "drotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 190 "drotm.f"
		w = dx[kx];
#line 191 "drotm.f"
		z__ = dy[ky];
#line 192 "drotm.f"
		dx[kx] = w * dh11 + z__;
#line 193 "drotm.f"
		dy[ky] = -w + dh22 * z__;
#line 194 "drotm.f"
		kx += *incx;
#line 195 "drotm.f"
		ky += *incy;
#line 196 "drotm.f"
	    }
#line 197 "drotm.f"
	}
#line 198 "drotm.f"
    }
#line 199 "drotm.f"
    return 0;
} /* drotm_ */

