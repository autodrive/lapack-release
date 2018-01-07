#line 1 "srotmg.f"
/* srotmg.f -- translated by f2c (version 20100827).
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

#line 1 "srotmg.f"
/* > \brief \b SROTMG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM) */

/*       .. Scalar Arguments .. */
/*       REAL SD1,SD2,SX1,SY1 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SPARAM(5) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/* >    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*>    SY2)**T. */
/* >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */
/* > */
/* >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/* >    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22 */
/* >    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE */
/* >    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.) */
/* > */
/* >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE */
/* >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE */
/* >    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in,out] SD1 */
/* > \verbatim */
/* >          SD1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] SD2 */
/* > \verbatim */
/* >          SD2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX1 */
/* > \verbatim */
/* >          SX1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SY1 */
/* > \verbatim */
/* >          SY1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] SPARAM */
/* > \verbatim */
/* >          SPARAM is REAL array, dimension 5 */
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

/* > \date November 2011 */

/* > \ingroup single_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int srotmg_(doublereal *sd1, doublereal *sd2, doublereal *
	sx1, doublereal *sy1, doublereal *sparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal gam = 4096.;
    static doublereal gamsq = 16777200.;
    static doublereal rgamsq = 5.96046e-8;

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal su, sp1, sp2, sq1, sq2, sh11, sh12, sh21, sh22, sflag, 
	    stemp;


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

#line 116 "srotmg.f"
    /* Parameter adjustments */
#line 116 "srotmg.f"
    --sparam;
#line 116 "srotmg.f"

#line 116 "srotmg.f"
    /* Function Body */
/*     .. */
#line 120 "srotmg.f"
    if (*sd1 < zero) {
/*        GO ZERO-H-D-AND-SX1.. */
#line 122 "srotmg.f"
	sflag = -one;
#line 123 "srotmg.f"
	sh11 = zero;
#line 124 "srotmg.f"
	sh12 = zero;
#line 125 "srotmg.f"
	sh21 = zero;
#line 126 "srotmg.f"
	sh22 = zero;

#line 128 "srotmg.f"
	*sd1 = zero;
#line 129 "srotmg.f"
	*sd2 = zero;
#line 130 "srotmg.f"
	*sx1 = zero;
#line 131 "srotmg.f"
    } else {
/*        CASE-SD1-NONNEGATIVE */
#line 133 "srotmg.f"
	sp2 = *sd2 * *sy1;
#line 134 "srotmg.f"
	if (sp2 == zero) {
#line 135 "srotmg.f"
	    sflag = -two;
#line 136 "srotmg.f"
	    sparam[1] = sflag;
#line 137 "srotmg.f"
	    return 0;
#line 138 "srotmg.f"
	}
/*        REGULAR-CASE.. */
#line 140 "srotmg.f"
	sp1 = *sd1 * *sx1;
#line 141 "srotmg.f"
	sq2 = sp2 * *sy1;
#line 142 "srotmg.f"
	sq1 = sp1 * *sx1;

#line 144 "srotmg.f"
	if (abs(sq1) > abs(sq2)) {
#line 145 "srotmg.f"
	    sh21 = -(*sy1) / *sx1;
#line 146 "srotmg.f"
	    sh12 = sp2 / sp1;

#line 148 "srotmg.f"
	    su = one - sh12 * sh21;

#line 150 "srotmg.f"
	    if (su > zero) {
#line 151 "srotmg.f"
		sflag = zero;
#line 152 "srotmg.f"
		*sd1 /= su;
#line 153 "srotmg.f"
		*sd2 /= su;
#line 154 "srotmg.f"
		*sx1 *= su;
#line 155 "srotmg.f"
	    }
#line 156 "srotmg.f"
	} else {
#line 158 "srotmg.f"
	    if (sq2 < zero) {
/*              GO ZERO-H-D-AND-SX1.. */
#line 160 "srotmg.f"
		sflag = -one;
#line 161 "srotmg.f"
		sh11 = zero;
#line 162 "srotmg.f"
		sh12 = zero;
#line 163 "srotmg.f"
		sh21 = zero;
#line 164 "srotmg.f"
		sh22 = zero;

#line 166 "srotmg.f"
		*sd1 = zero;
#line 167 "srotmg.f"
		*sd2 = zero;
#line 168 "srotmg.f"
		*sx1 = zero;
#line 169 "srotmg.f"
	    } else {
#line 170 "srotmg.f"
		sflag = one;
#line 171 "srotmg.f"
		sh11 = sp1 / sp2;
#line 172 "srotmg.f"
		sh22 = *sx1 / *sy1;
#line 173 "srotmg.f"
		su = one + sh11 * sh22;
#line 174 "srotmg.f"
		stemp = *sd2 / su;
#line 175 "srotmg.f"
		*sd2 = *sd1 / su;
#line 176 "srotmg.f"
		*sd1 = stemp;
#line 177 "srotmg.f"
		*sx1 = *sy1 * su;
#line 178 "srotmg.f"
	    }
#line 179 "srotmg.f"
	}
/*     PROCESURE..SCALE-CHECK */
#line 182 "srotmg.f"
	if (*sd1 != zero) {
#line 183 "srotmg.f"
	    while(*sd1 <= rgamsq || *sd1 >= gamsq) {
#line 184 "srotmg.f"
		if (sflag == zero) {
#line 185 "srotmg.f"
		    sh11 = one;
#line 186 "srotmg.f"
		    sh22 = one;
#line 187 "srotmg.f"
		    sflag = -one;
#line 188 "srotmg.f"
		} else {
#line 189 "srotmg.f"
		    sh21 = -one;
#line 190 "srotmg.f"
		    sh12 = one;
#line 191 "srotmg.f"
		    sflag = -one;
#line 192 "srotmg.f"
		}
#line 193 "srotmg.f"
		if (*sd1 <= rgamsq) {
/* Computing 2nd power */
#line 194 "srotmg.f"
		    d__1 = gam;
#line 194 "srotmg.f"
		    *sd1 *= d__1 * d__1;
#line 195 "srotmg.f"
		    *sx1 /= gam;
#line 196 "srotmg.f"
		    sh11 /= gam;
#line 197 "srotmg.f"
		    sh12 /= gam;
#line 198 "srotmg.f"
		} else {
/* Computing 2nd power */
#line 199 "srotmg.f"
		    d__1 = gam;
#line 199 "srotmg.f"
		    *sd1 /= d__1 * d__1;
#line 200 "srotmg.f"
		    *sx1 *= gam;
#line 201 "srotmg.f"
		    sh11 *= gam;
#line 202 "srotmg.f"
		    sh12 *= gam;
#line 203 "srotmg.f"
		}
#line 204 "srotmg.f"
	    }
#line 205 "srotmg.f"
	}
#line 207 "srotmg.f"
	if (*sd2 != zero) {
#line 208 "srotmg.f"
	    while(abs(*sd2) <= rgamsq || abs(*sd2) >= gamsq) {
#line 209 "srotmg.f"
		if (sflag == zero) {
#line 210 "srotmg.f"
		    sh11 = one;
#line 211 "srotmg.f"
		    sh22 = one;
#line 212 "srotmg.f"
		    sflag = -one;
#line 213 "srotmg.f"
		} else {
#line 214 "srotmg.f"
		    sh21 = -one;
#line 215 "srotmg.f"
		    sh12 = one;
#line 216 "srotmg.f"
		    sflag = -one;
#line 217 "srotmg.f"
		}
#line 218 "srotmg.f"
		if (abs(*sd2) <= rgamsq) {
/* Computing 2nd power */
#line 219 "srotmg.f"
		    d__1 = gam;
#line 219 "srotmg.f"
		    *sd2 *= d__1 * d__1;
#line 220 "srotmg.f"
		    sh21 /= gam;
#line 221 "srotmg.f"
		    sh22 /= gam;
#line 222 "srotmg.f"
		} else {
/* Computing 2nd power */
#line 223 "srotmg.f"
		    d__1 = gam;
#line 223 "srotmg.f"
		    *sd2 /= d__1 * d__1;
#line 224 "srotmg.f"
		    sh21 *= gam;
#line 225 "srotmg.f"
		    sh22 *= gam;
#line 226 "srotmg.f"
		}
#line 227 "srotmg.f"
	    }
#line 228 "srotmg.f"
	}
#line 230 "srotmg.f"
    }
#line 232 "srotmg.f"
    if (sflag < zero) {
#line 233 "srotmg.f"
	sparam[2] = sh11;
#line 234 "srotmg.f"
	sparam[3] = sh21;
#line 235 "srotmg.f"
	sparam[4] = sh12;
#line 236 "srotmg.f"
	sparam[5] = sh22;
#line 237 "srotmg.f"
    } else if (sflag == zero) {
#line 238 "srotmg.f"
	sparam[3] = sh21;
#line 239 "srotmg.f"
	sparam[4] = sh12;
#line 240 "srotmg.f"
    } else {
#line 241 "srotmg.f"
	sparam[2] = sh11;
#line 242 "srotmg.f"
	sparam[5] = sh22;
#line 243 "srotmg.f"
    }
#line 245 "srotmg.f"
    sparam[1] = sflag;
#line 246 "srotmg.f"
    return 0;
} /* srotmg_ */

