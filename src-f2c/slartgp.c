#line 1 "slartgp.f"
/* slartgp.f -- translated by f2c (version 20100827).
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

#line 1 "slartgp.f"
/* Table of constant values */

static doublereal c_b6 = 1.;

/* > \brief \b SLARTGP generates a plane rotation so that the diagonal is nonnegative. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARTGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARTGP( F, G, CS, SN, R ) */

/*       .. Scalar Arguments .. */
/*       REAL               CS, F, G, R, SN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARTGP generates a plane rotation so that */
/* > */
/* >    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1. */
/* >    [ -SN  CS  ]     [ G ]     [ 0 ] */
/* > */
/* > This is a slower, more accurate version of the Level 1 BLAS routine SROTG, */
/* > with the following other differences: */
/* >    F and G are unchanged on return. */
/* >    If G=0, then CS=(+/-)1 and SN=0. */
/* >    If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1. */
/* > */
/* > The sign is chosen so that R >= 0. */
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

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slartgp_(doublereal *f, doublereal *g, doublereal *cs, 
	doublereal *sn, doublereal *r__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), d_sign(
	    doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal f1, g1, eps, scale;
    static integer count;
    static doublereal safmn2, safmx2;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
#line 138 "slartgp.f"
    safmin = slamch_("S", (ftnlen)1);
#line 139 "slartgp.f"
    eps = slamch_("E", (ftnlen)1);
#line 140 "slartgp.f"
    d__1 = slamch_("B", (ftnlen)1);
#line 140 "slartgp.f"
    i__1 = (integer) (log(safmin / eps) / log(slamch_("B", (ftnlen)1)) / 2.);
#line 140 "slartgp.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 142 "slartgp.f"
    safmx2 = 1. / safmn2;
/*        FIRST = .FALSE. */
/*     END IF */
#line 145 "slartgp.f"
    if (*g == 0.) {
#line 146 "slartgp.f"
	*cs = d_sign(&c_b6, f);
#line 147 "slartgp.f"
	*sn = 0.;
#line 148 "slartgp.f"
	*r__ = abs(*f);
#line 149 "slartgp.f"
    } else if (*f == 0.) {
#line 150 "slartgp.f"
	*cs = 0.;
#line 151 "slartgp.f"
	*sn = d_sign(&c_b6, g);
#line 152 "slartgp.f"
	*r__ = abs(*g);
#line 153 "slartgp.f"
    } else {
#line 154 "slartgp.f"
	f1 = *f;
#line 155 "slartgp.f"
	g1 = *g;
/* Computing MAX */
#line 156 "slartgp.f"
	d__1 = abs(f1), d__2 = abs(g1);
#line 156 "slartgp.f"
	scale = max(d__1,d__2);
#line 157 "slartgp.f"
	if (scale >= safmx2) {
#line 158 "slartgp.f"
	    count = 0;
#line 159 "slartgp.f"
L10:
#line 160 "slartgp.f"
	    ++count;
#line 161 "slartgp.f"
	    f1 *= safmn2;
#line 162 "slartgp.f"
	    g1 *= safmn2;
/* Computing MAX */
#line 163 "slartgp.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 163 "slartgp.f"
	    scale = max(d__1,d__2);
#line 164 "slartgp.f"
	    if (scale >= safmx2) {
#line 164 "slartgp.f"
		goto L10;
#line 164 "slartgp.f"
	    }
/* Computing 2nd power */
#line 166 "slartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 166 "slartgp.f"
	    d__2 = g1;
#line 166 "slartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 167 "slartgp.f"
	    *cs = f1 / *r__;
#line 168 "slartgp.f"
	    *sn = g1 / *r__;
#line 169 "slartgp.f"
	    i__1 = count;
#line 169 "slartgp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 170 "slartgp.f"
		*r__ *= safmx2;
#line 171 "slartgp.f"
/* L20: */
#line 171 "slartgp.f"
	    }
#line 172 "slartgp.f"
	} else if (scale <= safmn2) {
#line 173 "slartgp.f"
	    count = 0;
#line 174 "slartgp.f"
L30:
#line 175 "slartgp.f"
	    ++count;
#line 176 "slartgp.f"
	    f1 *= safmx2;
#line 177 "slartgp.f"
	    g1 *= safmx2;
/* Computing MAX */
#line 178 "slartgp.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 178 "slartgp.f"
	    scale = max(d__1,d__2);
#line 179 "slartgp.f"
	    if (scale <= safmn2) {
#line 179 "slartgp.f"
		goto L30;
#line 179 "slartgp.f"
	    }
/* Computing 2nd power */
#line 181 "slartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 181 "slartgp.f"
	    d__2 = g1;
#line 181 "slartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 182 "slartgp.f"
	    *cs = f1 / *r__;
#line 183 "slartgp.f"
	    *sn = g1 / *r__;
#line 184 "slartgp.f"
	    i__1 = count;
#line 184 "slartgp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "slartgp.f"
		*r__ *= safmn2;
#line 186 "slartgp.f"
/* L40: */
#line 186 "slartgp.f"
	    }
#line 187 "slartgp.f"
	} else {
/* Computing 2nd power */
#line 188 "slartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 188 "slartgp.f"
	    d__2 = g1;
#line 188 "slartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 189 "slartgp.f"
	    *cs = f1 / *r__;
#line 190 "slartgp.f"
	    *sn = g1 / *r__;
#line 191 "slartgp.f"
	}
#line 192 "slartgp.f"
	if (*r__ < 0.) {
#line 193 "slartgp.f"
	    *cs = -(*cs);
#line 194 "slartgp.f"
	    *sn = -(*sn);
#line 195 "slartgp.f"
	    *r__ = -(*r__);
#line 196 "slartgp.f"
	}
#line 197 "slartgp.f"
    }
#line 198 "slartgp.f"
    return 0;

/*     End of SLARTG */

} /* slartgp_ */

