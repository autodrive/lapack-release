#line 1 "slartg.f"
/* slartg.f -- translated by f2c (version 20100827).
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

#line 1 "slartg.f"
/* > \brief \b SLARTG generates a plane rotation with real cosine and real sine. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARTG( F, G, CS, SN, R ) */

/*       .. Scalar Arguments .. */
/*       REAL               CS, F, G, R, SN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARTG generate a plane rotation so that */
/* > */
/* >    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1. */
/* >    [ -SN  CS  ]     [ G ]     [ 0 ] */
/* > */
/* > This is a slower, more accurate version of the BLAS1 routine SROTG, */
/* > with the following other differences: */
/* >    F and G are unchanged on return. */
/* >    If G=0, then CS=1 and SN=0. */
/* >    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any */
/* >       floating point operations (saves work in SBDSQR when */
/* >       there are zeros on the diagonal). */
/* > */
/* > If F exceeds G in magnitude, CS will be positive. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is REAL */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is REAL */
/* >          The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is REAL */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is REAL */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is REAL */
/* >          The nonzero component of the rotated vector. */
/* > */
/* >  This version has a few statements commented out for thread safety */
/* >  (machine parameters are computed on each entry). 10 feb 03, SJH. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slartg_(doublereal *f, doublereal *g, doublereal *cs, 
	doublereal *sn, doublereal *r__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal f1, g1, eps, scale;
    static integer count;
    static doublereal safmn2, safmx2;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     LOGICAL            FIRST */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2 */
/*     .. */
/*     .. Data statements .. */
/*     DATA               FIRST / .TRUE. / */
/*     .. */
/*     .. Executable Statements .. */

/*     IF( FIRST ) THEN */
#line 140 "slartg.f"
    safmin = slamch_("S", (ftnlen)1);
#line 141 "slartg.f"
    eps = slamch_("E", (ftnlen)1);
#line 142 "slartg.f"
    d__1 = slamch_("B", (ftnlen)1);
#line 142 "slartg.f"
    i__1 = (integer) (log(safmin / eps) / log(slamch_("B", (ftnlen)1)) / 2.);
#line 142 "slartg.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 144 "slartg.f"
    safmx2 = 1. / safmn2;
/*        FIRST = .FALSE. */
/*     END IF */
#line 147 "slartg.f"
    if (*g == 0.) {
#line 148 "slartg.f"
	*cs = 1.;
#line 149 "slartg.f"
	*sn = 0.;
#line 150 "slartg.f"
	*r__ = *f;
#line 151 "slartg.f"
    } else if (*f == 0.) {
#line 152 "slartg.f"
	*cs = 0.;
#line 153 "slartg.f"
	*sn = 1.;
#line 154 "slartg.f"
	*r__ = *g;
#line 155 "slartg.f"
    } else {
#line 156 "slartg.f"
	f1 = *f;
#line 157 "slartg.f"
	g1 = *g;
/* Computing MAX */
#line 158 "slartg.f"
	d__1 = abs(f1), d__2 = abs(g1);
#line 158 "slartg.f"
	scale = max(d__1,d__2);
#line 159 "slartg.f"
	if (scale >= safmx2) {
#line 160 "slartg.f"
	    count = 0;
#line 161 "slartg.f"
L10:
#line 162 "slartg.f"
	    ++count;
#line 163 "slartg.f"
	    f1 *= safmn2;
#line 164 "slartg.f"
	    g1 *= safmn2;
/* Computing MAX */
#line 165 "slartg.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 165 "slartg.f"
	    scale = max(d__1,d__2);
#line 166 "slartg.f"
	    if (scale >= safmx2) {
#line 166 "slartg.f"
		goto L10;
#line 166 "slartg.f"
	    }
/* Computing 2nd power */
#line 168 "slartg.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 168 "slartg.f"
	    d__2 = g1;
#line 168 "slartg.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 169 "slartg.f"
	    *cs = f1 / *r__;
#line 170 "slartg.f"
	    *sn = g1 / *r__;
#line 171 "slartg.f"
	    i__1 = count;
#line 171 "slartg.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 172 "slartg.f"
		*r__ *= safmx2;
#line 173 "slartg.f"
/* L20: */
#line 173 "slartg.f"
	    }
#line 174 "slartg.f"
	} else if (scale <= safmn2) {
#line 175 "slartg.f"
	    count = 0;
#line 176 "slartg.f"
L30:
#line 177 "slartg.f"
	    ++count;
#line 178 "slartg.f"
	    f1 *= safmx2;
#line 179 "slartg.f"
	    g1 *= safmx2;
/* Computing MAX */
#line 180 "slartg.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 180 "slartg.f"
	    scale = max(d__1,d__2);
#line 181 "slartg.f"
	    if (scale <= safmn2) {
#line 181 "slartg.f"
		goto L30;
#line 181 "slartg.f"
	    }
/* Computing 2nd power */
#line 183 "slartg.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 183 "slartg.f"
	    d__2 = g1;
#line 183 "slartg.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 184 "slartg.f"
	    *cs = f1 / *r__;
#line 185 "slartg.f"
	    *sn = g1 / *r__;
#line 186 "slartg.f"
	    i__1 = count;
#line 186 "slartg.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 187 "slartg.f"
		*r__ *= safmn2;
#line 188 "slartg.f"
/* L40: */
#line 188 "slartg.f"
	    }
#line 189 "slartg.f"
	} else {
/* Computing 2nd power */
#line 190 "slartg.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 190 "slartg.f"
	    d__2 = g1;
#line 190 "slartg.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 191 "slartg.f"
	    *cs = f1 / *r__;
#line 192 "slartg.f"
	    *sn = g1 / *r__;
#line 193 "slartg.f"
	}
#line 194 "slartg.f"
	if (abs(*f) > abs(*g) && *cs < 0.) {
#line 195 "slartg.f"
	    *cs = -(*cs);
#line 196 "slartg.f"
	    *sn = -(*sn);
#line 197 "slartg.f"
	    *r__ = -(*r__);
#line 198 "slartg.f"
	}
#line 199 "slartg.f"
    }
#line 200 "slartg.f"
    return 0;

/*     End of SLARTG */

} /* slartg_ */

