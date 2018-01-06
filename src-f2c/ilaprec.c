#line 1 "ilaprec.f"
/* ilaprec.f -- translated by f2c (version 20100827).
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

#line 1 "ilaprec.f"
/* > \brief \b ILAPREC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILAPREC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaprec
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaprec
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaprec
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILAPREC( PREC ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          PREC */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine translated from a character string specifying an */
/* > intermediate precision to the relevant BLAST-specified integer */
/* > constant. */
/* > */
/* > ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a */
/* > character indicating a supported intermediate precision.  Otherwise */
/* > ILAPREC returns the constant value corresponding to PREC. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */


/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
integer ilaprec_(char *prec, ftnlen prec_len)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */
#line 83 "ilaprec.f"
    if (lsame_(prec, "S", (ftnlen)1, (ftnlen)1)) {
#line 84 "ilaprec.f"
	ret_val = 211;
#line 85 "ilaprec.f"
    } else if (lsame_(prec, "D", (ftnlen)1, (ftnlen)1)) {
#line 86 "ilaprec.f"
	ret_val = 212;
#line 87 "ilaprec.f"
    } else if (lsame_(prec, "I", (ftnlen)1, (ftnlen)1)) {
#line 88 "ilaprec.f"
	ret_val = 213;
#line 89 "ilaprec.f"
    } else if (lsame_(prec, "X", (ftnlen)1, (ftnlen)1) || lsame_(prec, "E", (
	    ftnlen)1, (ftnlen)1)) {
#line 90 "ilaprec.f"
	ret_val = 214;
#line 91 "ilaprec.f"
    } else {
#line 92 "ilaprec.f"
	ret_val = -1;
#line 93 "ilaprec.f"
    }
#line 94 "ilaprec.f"
    return ret_val;

/*     End of ILAPREC */

} /* ilaprec_ */

