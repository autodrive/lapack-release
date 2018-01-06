#line 1 "slamrg.f"
/* slamrg.f -- translated by f2c (version 20100827).
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

#line 1 "slamrg.f"
/* > \brief \b SLAMRG creates a permutation list to merge the entries of two independently sorted sets into a 
single set sorted in ascending order. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAMRG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slamrg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slamrg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slamrg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAMRG( N1, N2, A, STRD1, STRD2, INDEX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            N1, N2, STRD1, STRD2 */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            INDEX( * ) */
/*       REAL               A( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAMRG will create a permutation list which will merge the elements */
/* > of A (which is composed of two independently sorted sets) into a */
/* > single set which is sorted in ascending order. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* >          N2 is INTEGER */
/* >         These arguments contain the respective lengths of the two */
/* >         sorted lists to be merged. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (N1+N2) */
/* >         The first N1 elements of A contain a list of numbers which */
/* >         are sorted in either ascending or descending order.  Likewise */
/* >         for the final N2 elements. */
/* > \endverbatim */
/* > */
/* > \param[in] STRD1 */
/* > \verbatim */
/* >          STRD1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] STRD2 */
/* > \verbatim */
/* >          STRD2 is INTEGER */
/* >         These are the strides to be taken through the array A. */
/* >         Allowable strides are 1 and -1.  They indicate whether a */
/* >         subset of A is sorted in ascending (STRDx = 1) or descending */
/* >         (STRDx = -1) order. */
/* > \endverbatim */
/* > */
/* > \param[out] INDEX */
/* > \verbatim */
/* >          INDEX is INTEGER array, dimension (N1+N2) */
/* >         On exit this array will contain a permutation such that */
/* >         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be */
/* >         sorted in ascending order. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int slamrg_(integer *n1, integer *n2, doublereal *a, integer 
	*strd1, integer *strd2, integer *index)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ind1, ind2, n1sv, n2sv;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 122 "slamrg.f"
    /* Parameter adjustments */
#line 122 "slamrg.f"
    --index;
#line 122 "slamrg.f"
    --a;
#line 122 "slamrg.f"

#line 122 "slamrg.f"
    /* Function Body */
#line 122 "slamrg.f"
    n1sv = *n1;
#line 123 "slamrg.f"
    n2sv = *n2;
#line 124 "slamrg.f"
    if (*strd1 > 0) {
#line 125 "slamrg.f"
	ind1 = 1;
#line 126 "slamrg.f"
    } else {
#line 127 "slamrg.f"
	ind1 = *n1;
#line 128 "slamrg.f"
    }
#line 129 "slamrg.f"
    if (*strd2 > 0) {
#line 130 "slamrg.f"
	ind2 = *n1 + 1;
#line 131 "slamrg.f"
    } else {
#line 132 "slamrg.f"
	ind2 = *n1 + *n2;
#line 133 "slamrg.f"
    }
#line 134 "slamrg.f"
    i__ = 1;
/*     while ( (N1SV > 0) & (N2SV > 0) ) */
#line 136 "slamrg.f"
L10:
#line 137 "slamrg.f"
    if (n1sv > 0 && n2sv > 0) {
#line 138 "slamrg.f"
	if (a[ind1] <= a[ind2]) {
#line 139 "slamrg.f"
	    index[i__] = ind1;
#line 140 "slamrg.f"
	    ++i__;
#line 141 "slamrg.f"
	    ind1 += *strd1;
#line 142 "slamrg.f"
	    --n1sv;
#line 143 "slamrg.f"
	} else {
#line 144 "slamrg.f"
	    index[i__] = ind2;
#line 145 "slamrg.f"
	    ++i__;
#line 146 "slamrg.f"
	    ind2 += *strd2;
#line 147 "slamrg.f"
	    --n2sv;
#line 148 "slamrg.f"
	}
#line 149 "slamrg.f"
	goto L10;
#line 150 "slamrg.f"
    }
/*     end while */
#line 152 "slamrg.f"
    if (n1sv == 0) {
#line 153 "slamrg.f"
	i__1 = n2sv;
#line 153 "slamrg.f"
	for (n1sv = 1; n1sv <= i__1; ++n1sv) {
#line 154 "slamrg.f"
	    index[i__] = ind2;
#line 155 "slamrg.f"
	    ++i__;
#line 156 "slamrg.f"
	    ind2 += *strd2;
#line 157 "slamrg.f"
/* L20: */
#line 157 "slamrg.f"
	}
#line 158 "slamrg.f"
    } else {
/*     N2SV .EQ. 0 */
#line 160 "slamrg.f"
	i__1 = n1sv;
#line 160 "slamrg.f"
	for (n2sv = 1; n2sv <= i__1; ++n2sv) {
#line 161 "slamrg.f"
	    index[i__] = ind1;
#line 162 "slamrg.f"
	    ++i__;
#line 163 "slamrg.f"
	    ind1 += *strd1;
#line 164 "slamrg.f"
/* L30: */
#line 164 "slamrg.f"
	}
#line 165 "slamrg.f"
    }

#line 167 "slamrg.f"
    return 0;

/*     End of SLAMRG */

} /* slamrg_ */

