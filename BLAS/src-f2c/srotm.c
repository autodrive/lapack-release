#line 1 "srotm.f"
/* srotm.f -- translated by f2c (version 20100827).
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

#line 1 "srotm.f"
/* > \brief \b SROTM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SPARAM(5),SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */
/* > */
/* >    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN */
/* >    (SX**T) */
/* > */
/* >    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/* >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY. */
/* >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */
/* > */
/* >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/* >    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM. */
/* > */
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
/* > \param[in] SPARAM */
/* > \verbatim */
/* >          SPARAM is REAL array, dimension (5) */
/* >     SPARAM(1)=SFLAG */
/* >     SPARAM(2)=SH11 */
/* >     SPARAM(3)=SH21 */
/* >     SPARAM(4)=SH12 */
/* >     SPARAM(5)=SH22 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup single_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int srotm_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy, doublereal *sparam)
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
    static doublereal sh11, sh12, sh21, sh22, sflag;
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
#line 119 "srotm.f"
    /* Parameter adjustments */
#line 119 "srotm.f"
    --sparam;
#line 119 "srotm.f"
    --sy;
#line 119 "srotm.f"
    --sx;
#line 119 "srotm.f"

#line 119 "srotm.f"
    /* Function Body */
/*     .. */

#line 122 "srotm.f"
    sflag = sparam[1];
#line 123 "srotm.f"
    if (*n <= 0 || sflag + two == zero) {
#line 123 "srotm.f"
	return 0;
#line 123 "srotm.f"
    }
#line 124 "srotm.f"
    if (*incx == *incy && *incx > 0) {

#line 126 "srotm.f"
	nsteps = *n * *incx;
#line 127 "srotm.f"
	if (sflag < zero) {
#line 128 "srotm.f"
	    sh11 = sparam[2];
#line 129 "srotm.f"
	    sh12 = sparam[4];
#line 130 "srotm.f"
	    sh21 = sparam[3];
#line 131 "srotm.f"
	    sh22 = sparam[5];
#line 132 "srotm.f"
	    i__1 = nsteps;
#line 132 "srotm.f"
	    i__2 = *incx;
#line 132 "srotm.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 133 "srotm.f"
		w = sx[i__];
#line 134 "srotm.f"
		z__ = sy[i__];
#line 135 "srotm.f"
		sx[i__] = w * sh11 + z__ * sh12;
#line 136 "srotm.f"
		sy[i__] = w * sh21 + z__ * sh22;
#line 137 "srotm.f"
	    }
#line 138 "srotm.f"
	} else if (sflag == zero) {
#line 139 "srotm.f"
	    sh12 = sparam[4];
#line 140 "srotm.f"
	    sh21 = sparam[3];
#line 141 "srotm.f"
	    i__2 = nsteps;
#line 141 "srotm.f"
	    i__1 = *incx;
#line 141 "srotm.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 142 "srotm.f"
		w = sx[i__];
#line 143 "srotm.f"
		z__ = sy[i__];
#line 144 "srotm.f"
		sx[i__] = w + z__ * sh12;
#line 145 "srotm.f"
		sy[i__] = w * sh21 + z__;
#line 146 "srotm.f"
	    }
#line 147 "srotm.f"
	} else {
#line 148 "srotm.f"
	    sh11 = sparam[2];
#line 149 "srotm.f"
	    sh22 = sparam[5];
#line 150 "srotm.f"
	    i__1 = nsteps;
#line 150 "srotm.f"
	    i__2 = *incx;
#line 150 "srotm.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 151 "srotm.f"
		w = sx[i__];
#line 152 "srotm.f"
		z__ = sy[i__];
#line 153 "srotm.f"
		sx[i__] = w * sh11 + z__;
#line 154 "srotm.f"
		sy[i__] = -w + sh22 * z__;
#line 155 "srotm.f"
	    }
#line 156 "srotm.f"
	}
#line 157 "srotm.f"
    } else {
#line 158 "srotm.f"
	kx = 1;
#line 159 "srotm.f"
	ky = 1;
#line 160 "srotm.f"
	if (*incx < 0) {
#line 160 "srotm.f"
	    kx = (1 - *n) * *incx + 1;
#line 160 "srotm.f"
	}
#line 161 "srotm.f"
	if (*incy < 0) {
#line 161 "srotm.f"
	    ky = (1 - *n) * *incy + 1;
#line 161 "srotm.f"
	}

#line 163 "srotm.f"
	if (sflag < zero) {
#line 164 "srotm.f"
	    sh11 = sparam[2];
#line 165 "srotm.f"
	    sh12 = sparam[4];
#line 166 "srotm.f"
	    sh21 = sparam[3];
#line 167 "srotm.f"
	    sh22 = sparam[5];
#line 168 "srotm.f"
	    i__2 = *n;
#line 168 "srotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "srotm.f"
		w = sx[kx];
#line 170 "srotm.f"
		z__ = sy[ky];
#line 171 "srotm.f"
		sx[kx] = w * sh11 + z__ * sh12;
#line 172 "srotm.f"
		sy[ky] = w * sh21 + z__ * sh22;
#line 173 "srotm.f"
		kx += *incx;
#line 174 "srotm.f"
		ky += *incy;
#line 175 "srotm.f"
	    }
#line 176 "srotm.f"
	} else if (sflag == zero) {
#line 177 "srotm.f"
	    sh12 = sparam[4];
#line 178 "srotm.f"
	    sh21 = sparam[3];
#line 179 "srotm.f"
	    i__2 = *n;
#line 179 "srotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 180 "srotm.f"
		w = sx[kx];
#line 181 "srotm.f"
		z__ = sy[ky];
#line 182 "srotm.f"
		sx[kx] = w + z__ * sh12;
#line 183 "srotm.f"
		sy[ky] = w * sh21 + z__;
#line 184 "srotm.f"
		kx += *incx;
#line 185 "srotm.f"
		ky += *incy;
#line 186 "srotm.f"
	    }
#line 187 "srotm.f"
	} else {
#line 188 "srotm.f"
	    sh11 = sparam[2];
#line 189 "srotm.f"
	    sh22 = sparam[5];
#line 190 "srotm.f"
	    i__2 = *n;
#line 190 "srotm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 191 "srotm.f"
		w = sx[kx];
#line 192 "srotm.f"
		z__ = sy[ky];
#line 193 "srotm.f"
		sx[kx] = w * sh11 + z__;
#line 194 "srotm.f"
		sy[ky] = -w + sh22 * z__;
#line 195 "srotm.f"
		kx += *incx;
#line 196 "srotm.f"
		ky += *incy;
#line 197 "srotm.f"
	    }
#line 198 "srotm.f"
	}
#line 199 "srotm.f"
    }
#line 200 "srotm.f"
    return 0;
} /* srotm_ */

