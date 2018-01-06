#line 1 "slasq6.f"
/* slasq6.f -- translated by f2c (version 20100827).
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

#line 1 "slasq6.f"
/* > \brief \b SLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASQ6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq6.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq6.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq6.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2 ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I0, N0, PP */
/*       REAL               DMIN, DMIN1, DMIN2, DN, DNM1, DNM2 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASQ6 computes one dqd (shift equal to zero) transform in */
/* > ping-pong form, with protection against underflow and overflow. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I0 */
/* > \verbatim */
/* >          I0 is INTEGER */
/* >        First index. */
/* > \endverbatim */
/* > */
/* > \param[in] N0 */
/* > \verbatim */
/* >          N0 is INTEGER */
/* >        Last index. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension ( 4*N ) */
/* >        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid */
/* >        an extra argument. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >        PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is REAL */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is REAL */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is REAL */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DN */
/* > \verbatim */
/* >          DN is REAL */
/* >        d(N0), the last value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DNM1 */
/* > \verbatim */
/* >          DNM1 is REAL */
/* >        d(N0-1). */
/* > \endverbatim */
/* > */
/* > \param[out] DNM2 */
/* > \verbatim */
/* >          DNM2 is REAL */
/* >        d(N0-2). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int slasq6_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
	 doublereal *dn, doublereal *dnm1, doublereal *dnm2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer j4, j4p2;
    static doublereal emin, temp;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameter .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Function .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 154 "slasq6.f"
    /* Parameter adjustments */
#line 154 "slasq6.f"
    --z__;
#line 154 "slasq6.f"

#line 154 "slasq6.f"
    /* Function Body */
#line 154 "slasq6.f"
    if (*n0 - *i0 - 1 <= 0) {
#line 154 "slasq6.f"
	return 0;
#line 154 "slasq6.f"
    }

#line 157 "slasq6.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 158 "slasq6.f"
    j4 = (*i0 << 2) + *pp - 3;
#line 159 "slasq6.f"
    emin = z__[j4 + 4];
#line 160 "slasq6.f"
    d__ = z__[j4];
#line 161 "slasq6.f"
    *dmin__ = d__;

#line 163 "slasq6.f"
    if (*pp == 0) {
#line 164 "slasq6.f"
	i__1 = *n0 - 3 << 2;
#line 164 "slasq6.f"
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 165 "slasq6.f"
	    z__[j4 - 2] = d__ + z__[j4 - 1];
#line 166 "slasq6.f"
	    if (z__[j4 - 2] == 0.) {
#line 167 "slasq6.f"
		z__[j4] = 0.;
#line 168 "slasq6.f"
		d__ = z__[j4 + 1];
#line 169 "slasq6.f"
		*dmin__ = d__;
#line 170 "slasq6.f"
		emin = 0.;
#line 171 "slasq6.f"
	    } else if (safmin * z__[j4 + 1] < z__[j4 - 2] && safmin * z__[j4 
		    - 2] < z__[j4 + 1]) {
#line 173 "slasq6.f"
		temp = z__[j4 + 1] / z__[j4 - 2];
#line 174 "slasq6.f"
		z__[j4] = z__[j4 - 1] * temp;
#line 175 "slasq6.f"
		d__ *= temp;
#line 176 "slasq6.f"
	    } else {
#line 177 "slasq6.f"
		z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
#line 178 "slasq6.f"
		d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
#line 179 "slasq6.f"
	    }
#line 180 "slasq6.f"
	    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 181 "slasq6.f"
	    d__1 = emin, d__2 = z__[j4];
#line 181 "slasq6.f"
	    emin = min(d__1,d__2);
#line 182 "slasq6.f"
/* L10: */
#line 182 "slasq6.f"
	}
#line 183 "slasq6.f"
    } else {
#line 184 "slasq6.f"
	i__1 = *n0 - 3 << 2;
#line 184 "slasq6.f"
	for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
#line 185 "slasq6.f"
	    z__[j4 - 3] = d__ + z__[j4];
#line 186 "slasq6.f"
	    if (z__[j4 - 3] == 0.) {
#line 187 "slasq6.f"
		z__[j4 - 1] = 0.;
#line 188 "slasq6.f"
		d__ = z__[j4 + 2];
#line 189 "slasq6.f"
		*dmin__ = d__;
#line 190 "slasq6.f"
		emin = 0.;
#line 191 "slasq6.f"
	    } else if (safmin * z__[j4 + 2] < z__[j4 - 3] && safmin * z__[j4 
		    - 3] < z__[j4 + 2]) {
#line 193 "slasq6.f"
		temp = z__[j4 + 2] / z__[j4 - 3];
#line 194 "slasq6.f"
		z__[j4 - 1] = z__[j4] * temp;
#line 195 "slasq6.f"
		d__ *= temp;
#line 196 "slasq6.f"
	    } else {
#line 197 "slasq6.f"
		z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
#line 198 "slasq6.f"
		d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
#line 199 "slasq6.f"
	    }
#line 200 "slasq6.f"
	    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
#line 201 "slasq6.f"
	    d__1 = emin, d__2 = z__[j4 - 1];
#line 201 "slasq6.f"
	    emin = min(d__1,d__2);
#line 202 "slasq6.f"
/* L20: */
#line 202 "slasq6.f"
	}
#line 203 "slasq6.f"
    }

/*     Unroll last two steps. */

#line 207 "slasq6.f"
    *dnm2 = d__;
#line 208 "slasq6.f"
    *dmin2 = *dmin__;
#line 209 "slasq6.f"
    j4 = (*n0 - 2 << 2) - *pp;
#line 210 "slasq6.f"
    j4p2 = j4 + (*pp << 1) - 1;
#line 211 "slasq6.f"
    z__[j4 - 2] = *dnm2 + z__[j4p2];
#line 212 "slasq6.f"
    if (z__[j4 - 2] == 0.) {
#line 213 "slasq6.f"
	z__[j4] = 0.;
#line 214 "slasq6.f"
	*dnm1 = z__[j4p2 + 2];
#line 215 "slasq6.f"
	*dmin__ = *dnm1;
#line 216 "slasq6.f"
	emin = 0.;
#line 217 "slasq6.f"
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
#line 219 "slasq6.f"
	temp = z__[j4p2 + 2] / z__[j4 - 2];
#line 220 "slasq6.f"
	z__[j4] = z__[j4p2] * temp;
#line 221 "slasq6.f"
	*dnm1 = *dnm2 * temp;
#line 222 "slasq6.f"
    } else {
#line 223 "slasq6.f"
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 224 "slasq6.f"
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
#line 225 "slasq6.f"
    }
#line 226 "slasq6.f"
    *dmin__ = min(*dmin__,*dnm1);

#line 228 "slasq6.f"
    *dmin1 = *dmin__;
#line 229 "slasq6.f"
    j4 += 4;
#line 230 "slasq6.f"
    j4p2 = j4 + (*pp << 1) - 1;
#line 231 "slasq6.f"
    z__[j4 - 2] = *dnm1 + z__[j4p2];
#line 232 "slasq6.f"
    if (z__[j4 - 2] == 0.) {
#line 233 "slasq6.f"
	z__[j4] = 0.;
#line 234 "slasq6.f"
	*dn = z__[j4p2 + 2];
#line 235 "slasq6.f"
	*dmin__ = *dn;
#line 236 "slasq6.f"
	emin = 0.;
#line 237 "slasq6.f"
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
#line 239 "slasq6.f"
	temp = z__[j4p2 + 2] / z__[j4 - 2];
#line 240 "slasq6.f"
	z__[j4] = z__[j4p2] * temp;
#line 241 "slasq6.f"
	*dn = *dnm1 * temp;
#line 242 "slasq6.f"
    } else {
#line 243 "slasq6.f"
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
#line 244 "slasq6.f"
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
#line 245 "slasq6.f"
    }
#line 246 "slasq6.f"
    *dmin__ = min(*dmin__,*dn);

#line 248 "slasq6.f"
    z__[j4 + 2] = *dn;
#line 249 "slasq6.f"
    z__[(*n0 << 2) - *pp] = emin;
#line 250 "slasq6.f"
    return 0;

/*     End of SLASQ6 */

} /* slasq6_ */

