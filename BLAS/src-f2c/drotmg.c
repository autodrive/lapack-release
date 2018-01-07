#line 1 "drotmg.f"
/* drotmg.f -- translated by f2c (version 20100827).
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

#line 1 "drotmg.f"
/* > \brief \b DROTMG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION DD1,DD2,DX1,DY1 */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DPARAM(5) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/* >    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T. */
/* >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */
/* > */
/* >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */
/* >    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22 */
/* >    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE */
/* >    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.) */
/* > */
/* >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE */
/* >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE */
/* >    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in,out] DD1 */
/* > \verbatim */
/* >          DD1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DD2 */
/* > \verbatim */
/* >          DD2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DX1 */
/* > \verbatim */
/* >          DX1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] DY1 */
/* > \verbatim */
/* >          DY1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DPARAM */
/* > \verbatim */
/* >          DPARAM is DOUBLE PRECISION array, dimension 5 */
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

/* > \date November 2011 */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int drotmg_(doublereal *dd1, doublereal *dd2, doublereal *
	dx1, doublereal *dy1, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal gam = 4096.;
    static doublereal gamsq = 16777216.;
    static doublereal rgamsq = 5.9604645e-8;

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal du, dp1, dp2, dq1, dq2, dh11, dh12, dh21, dh22, dflag, 
	    dtemp;


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */

#line 116 "drotmg.f"
    /* Parameter adjustments */
#line 116 "drotmg.f"
    --dparam;
#line 116 "drotmg.f"

#line 116 "drotmg.f"
    /* Function Body */
/*     .. */
#line 120 "drotmg.f"
    if (*dd1 < zero) {
/*        GO ZERO-H-D-AND-DX1.. */
#line 122 "drotmg.f"
	dflag = -one;
#line 123 "drotmg.f"
	dh11 = zero;
#line 124 "drotmg.f"
	dh12 = zero;
#line 125 "drotmg.f"
	dh21 = zero;
#line 126 "drotmg.f"
	dh22 = zero;

#line 128 "drotmg.f"
	*dd1 = zero;
#line 129 "drotmg.f"
	*dd2 = zero;
#line 130 "drotmg.f"
	*dx1 = zero;
#line 131 "drotmg.f"
    } else {
/*        CASE-DD1-NONNEGATIVE */
#line 133 "drotmg.f"
	dp2 = *dd2 * *dy1;
#line 134 "drotmg.f"
	if (dp2 == zero) {
#line 135 "drotmg.f"
	    dflag = -two;
#line 136 "drotmg.f"
	    dparam[1] = dflag;
#line 137 "drotmg.f"
	    return 0;
#line 138 "drotmg.f"
	}
/*        REGULAR-CASE.. */
#line 140 "drotmg.f"
	dp1 = *dd1 * *dx1;
#line 141 "drotmg.f"
	dq2 = dp2 * *dy1;
#line 142 "drotmg.f"
	dq1 = dp1 * *dx1;

#line 144 "drotmg.f"
	if (abs(dq1) > abs(dq2)) {
#line 145 "drotmg.f"
	    dh21 = -(*dy1) / *dx1;
#line 146 "drotmg.f"
	    dh12 = dp2 / dp1;

#line 148 "drotmg.f"
	    du = one - dh12 * dh21;

#line 150 "drotmg.f"
	    if (du > zero) {
#line 151 "drotmg.f"
		dflag = zero;
#line 152 "drotmg.f"
		*dd1 /= du;
#line 153 "drotmg.f"
		*dd2 /= du;
#line 154 "drotmg.f"
		*dx1 *= du;
#line 155 "drotmg.f"
	    }
#line 156 "drotmg.f"
	} else {
#line 158 "drotmg.f"
	    if (dq2 < zero) {
/*              GO ZERO-H-D-AND-DX1.. */
#line 160 "drotmg.f"
		dflag = -one;
#line 161 "drotmg.f"
		dh11 = zero;
#line 162 "drotmg.f"
		dh12 = zero;
#line 163 "drotmg.f"
		dh21 = zero;
#line 164 "drotmg.f"
		dh22 = zero;

#line 166 "drotmg.f"
		*dd1 = zero;
#line 167 "drotmg.f"
		*dd2 = zero;
#line 168 "drotmg.f"
		*dx1 = zero;
#line 169 "drotmg.f"
	    } else {
#line 170 "drotmg.f"
		dflag = one;
#line 171 "drotmg.f"
		dh11 = dp1 / dp2;
#line 172 "drotmg.f"
		dh22 = *dx1 / *dy1;
#line 173 "drotmg.f"
		du = one + dh11 * dh22;
#line 174 "drotmg.f"
		dtemp = *dd2 / du;
#line 175 "drotmg.f"
		*dd2 = *dd1 / du;
#line 176 "drotmg.f"
		*dd1 = dtemp;
#line 177 "drotmg.f"
		*dx1 = *dy1 * du;
#line 178 "drotmg.f"
	    }
#line 179 "drotmg.f"
	}
/*     PROCEDURE..SCALE-CHECK */
#line 182 "drotmg.f"
	if (*dd1 != zero) {
#line 183 "drotmg.f"
	    while(*dd1 <= rgamsq || *dd1 >= gamsq) {
#line 184 "drotmg.f"
		if (dflag == zero) {
#line 185 "drotmg.f"
		    dh11 = one;
#line 186 "drotmg.f"
		    dh22 = one;
#line 187 "drotmg.f"
		    dflag = -one;
#line 188 "drotmg.f"
		} else {
#line 189 "drotmg.f"
		    dh21 = -one;
#line 190 "drotmg.f"
		    dh12 = one;
#line 191 "drotmg.f"
		    dflag = -one;
#line 192 "drotmg.f"
		}
#line 193 "drotmg.f"
		if (*dd1 <= rgamsq) {
/* Computing 2nd power */
#line 194 "drotmg.f"
		    d__1 = gam;
#line 194 "drotmg.f"
		    *dd1 *= d__1 * d__1;
#line 195 "drotmg.f"
		    *dx1 /= gam;
#line 196 "drotmg.f"
		    dh11 /= gam;
#line 197 "drotmg.f"
		    dh12 /= gam;
#line 198 "drotmg.f"
		} else {
/* Computing 2nd power */
#line 199 "drotmg.f"
		    d__1 = gam;
#line 199 "drotmg.f"
		    *dd1 /= d__1 * d__1;
#line 200 "drotmg.f"
		    *dx1 *= gam;
#line 201 "drotmg.f"
		    dh11 *= gam;
#line 202 "drotmg.f"
		    dh12 *= gam;
#line 203 "drotmg.f"
		}
#line 204 "drotmg.f"
	    }
#line 205 "drotmg.f"
	}
#line 207 "drotmg.f"
	if (*dd2 != zero) {
#line 208 "drotmg.f"
	    while(abs(*dd2) <= rgamsq || abs(*dd2) >= gamsq) {
#line 209 "drotmg.f"
		if (dflag == zero) {
#line 210 "drotmg.f"
		    dh11 = one;
#line 211 "drotmg.f"
		    dh22 = one;
#line 212 "drotmg.f"
		    dflag = -one;
#line 213 "drotmg.f"
		} else {
#line 214 "drotmg.f"
		    dh21 = -one;
#line 215 "drotmg.f"
		    dh12 = one;
#line 216 "drotmg.f"
		    dflag = -one;
#line 217 "drotmg.f"
		}
#line 218 "drotmg.f"
		if (abs(*dd2) <= rgamsq) {
/* Computing 2nd power */
#line 219 "drotmg.f"
		    d__1 = gam;
#line 219 "drotmg.f"
		    *dd2 *= d__1 * d__1;
#line 220 "drotmg.f"
		    dh21 /= gam;
#line 221 "drotmg.f"
		    dh22 /= gam;
#line 222 "drotmg.f"
		} else {
/* Computing 2nd power */
#line 223 "drotmg.f"
		    d__1 = gam;
#line 223 "drotmg.f"
		    *dd2 /= d__1 * d__1;
#line 224 "drotmg.f"
		    dh21 *= gam;
#line 225 "drotmg.f"
		    dh22 *= gam;
#line 226 "drotmg.f"
		}
#line 227 "drotmg.f"
	    }
#line 228 "drotmg.f"
	}
#line 230 "drotmg.f"
    }
#line 232 "drotmg.f"
    if (dflag < zero) {
#line 233 "drotmg.f"
	dparam[2] = dh11;
#line 234 "drotmg.f"
	dparam[3] = dh21;
#line 235 "drotmg.f"
	dparam[4] = dh12;
#line 236 "drotmg.f"
	dparam[5] = dh22;
#line 237 "drotmg.f"
    } else if (dflag == zero) {
#line 238 "drotmg.f"
	dparam[3] = dh21;
#line 239 "drotmg.f"
	dparam[4] = dh12;
#line 240 "drotmg.f"
    } else {
#line 241 "drotmg.f"
	dparam[2] = dh11;
#line 242 "drotmg.f"
	dparam[5] = dh22;
#line 243 "drotmg.f"
    }
#line 245 "drotmg.f"
    dparam[1] = dflag;
#line 246 "drotmg.f"
    return 0;
} /* drotmg_ */

