#line 1 "dlartgp.f"
/* dlartgp.f -- translated by f2c (version 20100827).
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

#line 1 "dlartgp.f"
/* Table of constant values */

static doublereal c_b6 = 1.;

/* > \brief \b DLARTGP generates a plane rotation so that the diagonal is nonnegative. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARTGP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartgp
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartgp
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartgp
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARTGP( F, G, CS, SN, R ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   CS, F, G, R, SN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTGP generates a plane rotation so that */
/* > */
/* >    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1. */
/* >    [ -SN  CS  ]     [ G ]     [ 0 ] */
/* > */
/* > This is a slower, more accurate version of the Level 1 BLAS routine DROTG, */
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
/* >          F is DOUBLE PRECISION */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is DOUBLE PRECISION */
/* >          The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] CS */
/* > \verbatim */
/* >          CS is DOUBLE PRECISION */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] SN */
/* > \verbatim */
/* >          SN is DOUBLE PRECISION */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is DOUBLE PRECISION */
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
/* Subroutine */ int dlartgp_(doublereal *f, doublereal *g, doublereal *cs, 
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
    extern doublereal dlamch_(char *, ftnlen);
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
#line 138 "dlartgp.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 139 "dlartgp.f"
    eps = dlamch_("E", (ftnlen)1);
#line 140 "dlartgp.f"
    d__1 = dlamch_("B", (ftnlen)1);
#line 140 "dlartgp.f"
    i__1 = (integer) (log(safmin / eps) / log(dlamch_("B", (ftnlen)1)) / 2.);
#line 140 "dlartgp.f"
    safmn2 = pow_di(&d__1, &i__1);
#line 142 "dlartgp.f"
    safmx2 = 1. / safmn2;
/*        FIRST = .FALSE. */
/*     END IF */
#line 145 "dlartgp.f"
    if (*g == 0.) {
#line 146 "dlartgp.f"
	*cs = d_sign(&c_b6, f);
#line 147 "dlartgp.f"
	*sn = 0.;
#line 148 "dlartgp.f"
	*r__ = abs(*f);
#line 149 "dlartgp.f"
    } else if (*f == 0.) {
#line 150 "dlartgp.f"
	*cs = 0.;
#line 151 "dlartgp.f"
	*sn = d_sign(&c_b6, g);
#line 152 "dlartgp.f"
	*r__ = abs(*g);
#line 153 "dlartgp.f"
    } else {
#line 154 "dlartgp.f"
	f1 = *f;
#line 155 "dlartgp.f"
	g1 = *g;
/* Computing MAX */
#line 156 "dlartgp.f"
	d__1 = abs(f1), d__2 = abs(g1);
#line 156 "dlartgp.f"
	scale = max(d__1,d__2);
#line 157 "dlartgp.f"
	if (scale >= safmx2) {
#line 158 "dlartgp.f"
	    count = 0;
#line 159 "dlartgp.f"
L10:
#line 160 "dlartgp.f"
	    ++count;
#line 161 "dlartgp.f"
	    f1 *= safmn2;
#line 162 "dlartgp.f"
	    g1 *= safmn2;
/* Computing MAX */
#line 163 "dlartgp.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 163 "dlartgp.f"
	    scale = max(d__1,d__2);
#line 164 "dlartgp.f"
	    if (scale >= safmx2) {
#line 164 "dlartgp.f"
		goto L10;
#line 164 "dlartgp.f"
	    }
/* Computing 2nd power */
#line 166 "dlartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 166 "dlartgp.f"
	    d__2 = g1;
#line 166 "dlartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 167 "dlartgp.f"
	    *cs = f1 / *r__;
#line 168 "dlartgp.f"
	    *sn = g1 / *r__;
#line 169 "dlartgp.f"
	    i__1 = count;
#line 169 "dlartgp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 170 "dlartgp.f"
		*r__ *= safmx2;
#line 171 "dlartgp.f"
/* L20: */
#line 171 "dlartgp.f"
	    }
#line 172 "dlartgp.f"
	} else if (scale <= safmn2) {
#line 173 "dlartgp.f"
	    count = 0;
#line 174 "dlartgp.f"
L30:
#line 175 "dlartgp.f"
	    ++count;
#line 176 "dlartgp.f"
	    f1 *= safmx2;
#line 177 "dlartgp.f"
	    g1 *= safmx2;
/* Computing MAX */
#line 178 "dlartgp.f"
	    d__1 = abs(f1), d__2 = abs(g1);
#line 178 "dlartgp.f"
	    scale = max(d__1,d__2);
#line 179 "dlartgp.f"
	    if (scale <= safmn2) {
#line 179 "dlartgp.f"
		goto L30;
#line 179 "dlartgp.f"
	    }
/* Computing 2nd power */
#line 181 "dlartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 181 "dlartgp.f"
	    d__2 = g1;
#line 181 "dlartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 182 "dlartgp.f"
	    *cs = f1 / *r__;
#line 183 "dlartgp.f"
	    *sn = g1 / *r__;
#line 184 "dlartgp.f"
	    i__1 = count;
#line 184 "dlartgp.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "dlartgp.f"
		*r__ *= safmn2;
#line 186 "dlartgp.f"
/* L40: */
#line 186 "dlartgp.f"
	    }
#line 187 "dlartgp.f"
	} else {
/* Computing 2nd power */
#line 188 "dlartgp.f"
	    d__1 = f1;
/* Computing 2nd power */
#line 188 "dlartgp.f"
	    d__2 = g1;
#line 188 "dlartgp.f"
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
#line 189 "dlartgp.f"
	    *cs = f1 / *r__;
#line 190 "dlartgp.f"
	    *sn = g1 / *r__;
#line 191 "dlartgp.f"
	}
#line 192 "dlartgp.f"
	if (*r__ < 0.) {
#line 193 "dlartgp.f"
	    *cs = -(*cs);
#line 194 "dlartgp.f"
	    *sn = -(*sn);
#line 195 "dlartgp.f"
	    *r__ = -(*r__);
#line 196 "dlartgp.f"
	}
#line 197 "dlartgp.f"
    }
#line 198 "dlartgp.f"
    return 0;

/*     End of DLARTGP */

} /* dlartgp_ */

