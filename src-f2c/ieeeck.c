#line 1 "ieeeck.f"
/* ieeeck.f -- translated by f2c (version 20100827).
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

#line 1 "ieeeck.f"
/* > \brief \b IEEECK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IEEECK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ISPEC */
/*       REAL               ONE, ZERO */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > IEEECK is called from the ILAENV to verify that Infinity and */
/* > possibly NaN arithmetic is safe (i.e. will not trap). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >          Specifies whether to test just for inifinity arithmetic */
/* >          or whether to test for infinity and NaN arithmetic. */
/* >          = 0: Verify infinity arithmetic only. */
/* >          = 1: Verify infinity and NaN arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[in] ZERO */
/* > \verbatim */
/* >          ZERO is REAL */
/* >          Must contain the value 0.0 */
/* >          This is passed to prevent the compiler from optimizing */
/* >          away this code. */
/* > \endverbatim */
/* > */
/* > \param[in] ONE */
/* > \verbatim */
/* >          ONE is REAL */
/* >          Must contain the value 1.0 */
/* >          This is passed to prevent the compiler from optimizing */
/* >          away this code. */
/* > */
/* >  RETURN VALUE:  INTEGER */
/* >          = 0:  Arithmetic failed to produce the correct answers */
/* >          = 1:  Arithmetic produced the correct answers */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
integer ieeeck_(integer *ispec, doublereal *zero, doublereal *one)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static doublereal nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, 
	    negzro, newzro;


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
#line 102 "ieeeck.f"
    ret_val = 1;

#line 104 "ieeeck.f"
    posinf = *one / *zero;
#line 105 "ieeeck.f"
    if (posinf <= *one) {
#line 106 "ieeeck.f"
	ret_val = 0;
#line 107 "ieeeck.f"
	return ret_val;
#line 108 "ieeeck.f"
    }

#line 110 "ieeeck.f"
    neginf = -(*one) / *zero;
#line 111 "ieeeck.f"
    if (neginf >= *zero) {
#line 112 "ieeeck.f"
	ret_val = 0;
#line 113 "ieeeck.f"
	return ret_val;
#line 114 "ieeeck.f"
    }

#line 116 "ieeeck.f"
    negzro = *one / (neginf + *one);
#line 117 "ieeeck.f"
    if (negzro != *zero) {
#line 118 "ieeeck.f"
	ret_val = 0;
#line 119 "ieeeck.f"
	return ret_val;
#line 120 "ieeeck.f"
    }

#line 122 "ieeeck.f"
    neginf = *one / negzro;
#line 123 "ieeeck.f"
    if (neginf >= *zero) {
#line 124 "ieeeck.f"
	ret_val = 0;
#line 125 "ieeeck.f"
	return ret_val;
#line 126 "ieeeck.f"
    }

#line 128 "ieeeck.f"
    newzro = negzro + *zero;
#line 129 "ieeeck.f"
    if (newzro != *zero) {
#line 130 "ieeeck.f"
	ret_val = 0;
#line 131 "ieeeck.f"
	return ret_val;
#line 132 "ieeeck.f"
    }

#line 134 "ieeeck.f"
    posinf = *one / newzro;
#line 135 "ieeeck.f"
    if (posinf <= *one) {
#line 136 "ieeeck.f"
	ret_val = 0;
#line 137 "ieeeck.f"
	return ret_val;
#line 138 "ieeeck.f"
    }

#line 140 "ieeeck.f"
    neginf *= posinf;
#line 141 "ieeeck.f"
    if (neginf >= *zero) {
#line 142 "ieeeck.f"
	ret_val = 0;
#line 143 "ieeeck.f"
	return ret_val;
#line 144 "ieeeck.f"
    }

#line 146 "ieeeck.f"
    posinf *= posinf;
#line 147 "ieeeck.f"
    if (posinf <= *one) {
#line 148 "ieeeck.f"
	ret_val = 0;
#line 149 "ieeeck.f"
	return ret_val;
#line 150 "ieeeck.f"
    }




/*     Return if we were only asked to check infinity arithmetic */

#line 157 "ieeeck.f"
    if (*ispec == 0) {
#line 157 "ieeeck.f"
	return ret_val;
#line 157 "ieeeck.f"
    }

#line 160 "ieeeck.f"
    nan1 = posinf + neginf;

#line 162 "ieeeck.f"
    nan2 = posinf / neginf;

#line 164 "ieeeck.f"
    nan3 = posinf / posinf;

#line 166 "ieeeck.f"
    nan4 = posinf * *zero;

#line 168 "ieeeck.f"
    nan5 = neginf * negzro;

#line 170 "ieeeck.f"
    nan6 = nan5 * *zero;

#line 172 "ieeeck.f"
    if (nan1 == nan1) {
#line 173 "ieeeck.f"
	ret_val = 0;
#line 174 "ieeeck.f"
	return ret_val;
#line 175 "ieeeck.f"
    }

#line 177 "ieeeck.f"
    if (nan2 == nan2) {
#line 178 "ieeeck.f"
	ret_val = 0;
#line 179 "ieeeck.f"
	return ret_val;
#line 180 "ieeeck.f"
    }

#line 182 "ieeeck.f"
    if (nan3 == nan3) {
#line 183 "ieeeck.f"
	ret_val = 0;
#line 184 "ieeeck.f"
	return ret_val;
#line 185 "ieeeck.f"
    }

#line 187 "ieeeck.f"
    if (nan4 == nan4) {
#line 188 "ieeeck.f"
	ret_val = 0;
#line 189 "ieeeck.f"
	return ret_val;
#line 190 "ieeeck.f"
    }

#line 192 "ieeeck.f"
    if (nan5 == nan5) {
#line 193 "ieeeck.f"
	ret_val = 0;
#line 194 "ieeeck.f"
	return ret_val;
#line 195 "ieeeck.f"
    }

#line 197 "ieeeck.f"
    if (nan6 == nan6) {
#line 198 "ieeeck.f"
	ret_val = 0;
#line 199 "ieeeck.f"
	return ret_val;
#line 200 "ieeeck.f"
    }

#line 202 "ieeeck.f"
    return ret_val;
} /* ieeeck_ */

