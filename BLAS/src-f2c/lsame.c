#line 1 "lsame.f"
/* lsame.f -- translated by f2c (version 20100827).
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

#line 1 "lsame.f"
/* > \brief \b LSAME */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       LOGICAL FUNCTION LSAME(CA,CB) */

/*       .. Scalar Arguments .. */
/*       CHARACTER CA,CB */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/* > case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] CA */
/* > \verbatim */
/* >          CA is CHARACTER*1 */
/* > \endverbatim */
/* > */
/* > \param[in] CB */
/* > \verbatim */
/* >          CB is CHARACTER*1 */
/* >          CA and CB specify the single characters to be compared. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup aux_blas */

/*  ===================================================================== */
logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer inta, intb, zcode;


/*  -- Reference BLAS level1 routine (version 3.1) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */

/*     Test if the characters are equal */

#line 76 "lsame.f"
    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
#line 77 "lsame.f"
    if (ret_val) {
#line 77 "lsame.f"
	return ret_val;
#line 77 "lsame.f"
    }

/*     Now test for equivalence if both characters are alphabetic. */

#line 81 "lsame.f"
    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

#line 88 "lsame.f"
    inta = *(unsigned char *)ca;
#line 89 "lsame.f"
    intb = *(unsigned char *)cb;

#line 91 "lsame.f"
    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

#line 96 "lsame.f"
	if (inta >= 97 && inta <= 122) {
#line 96 "lsame.f"
	    inta += -32;
#line 96 "lsame.f"
	}
#line 97 "lsame.f"
	if (intb >= 97 && intb <= 122) {
#line 97 "lsame.f"
	    intb += -32;
#line 97 "lsame.f"
	}

#line 99 "lsame.f"
    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

#line 104 "lsame.f"
	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
#line 104 "lsame.f"
	    inta += 64;
#line 104 "lsame.f"
	}
#line 107 "lsame.f"
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
#line 107 "lsame.f"
	    intb += 64;
#line 107 "lsame.f"
	}

#line 111 "lsame.f"
    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

#line 116 "lsame.f"
	if (inta >= 225 && inta <= 250) {
#line 116 "lsame.f"
	    inta += -32;
#line 116 "lsame.f"
	}
#line 117 "lsame.f"
	if (intb >= 225 && intb <= 250) {
#line 117 "lsame.f"
	    intb += -32;
#line 117 "lsame.f"
	}
#line 118 "lsame.f"
    }
#line 119 "lsame.f"
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

#line 125 "lsame.f"
    return ret_val;
} /* lsame_ */

