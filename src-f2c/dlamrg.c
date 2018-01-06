#line 1 "dlamrg.f"
/* dlamrg.f -- translated by f2c (version 20100827).
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

#line 1 "dlamrg.f"
/* > \brief \b DLAMRG creates a permutation list to merge the entries of two independently sorted sets into a 
single set sorted in ascending order. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAMRG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlamrg.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlamrg.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlamrg.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            DTRD1, DTRD2, N1, N2 */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            INDEX( * ) */
/*       DOUBLE PRECISION   A( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAMRG will create a permutation list which will merge the elements */
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
/* >          A is DOUBLE PRECISION array, dimension (N1+N2) */
/* >         The first N1 elements of A contain a list of numbers which */
/* >         are sorted in either ascending or descending order.  Likewise */
/* >         for the final N2 elements. */
/* > \endverbatim */
/* > */
/* > \param[in] DTRD1 */
/* > \verbatim */
/* >          DTRD1 is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] DTRD2 */
/* > \verbatim */
/* >          DTRD2 is INTEGER */
/* >         These are the strides to be taken through the array A. */
/* >         Allowable strides are 1 and -1.  They indicate whether a */
/* >         subset of A is sorted in ascending (DTRDx = 1) or descending */
/* >         (DTRDx = -1) order. */
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
/* Subroutine */ int dlamrg_(integer *n1, integer *n2, doublereal *a, integer 
	*dtrd1, integer *dtrd2, integer *index)
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

#line 122 "dlamrg.f"
    /* Parameter adjustments */
#line 122 "dlamrg.f"
    --index;
#line 122 "dlamrg.f"
    --a;
#line 122 "dlamrg.f"

#line 122 "dlamrg.f"
    /* Function Body */
#line 122 "dlamrg.f"
    n1sv = *n1;
#line 123 "dlamrg.f"
    n2sv = *n2;
#line 124 "dlamrg.f"
    if (*dtrd1 > 0) {
#line 125 "dlamrg.f"
	ind1 = 1;
#line 126 "dlamrg.f"
    } else {
#line 127 "dlamrg.f"
	ind1 = *n1;
#line 128 "dlamrg.f"
    }
#line 129 "dlamrg.f"
    if (*dtrd2 > 0) {
#line 130 "dlamrg.f"
	ind2 = *n1 + 1;
#line 131 "dlamrg.f"
    } else {
#line 132 "dlamrg.f"
	ind2 = *n1 + *n2;
#line 133 "dlamrg.f"
    }
#line 134 "dlamrg.f"
    i__ = 1;
/*     while ( (N1SV > 0) & (N2SV > 0) ) */
#line 136 "dlamrg.f"
L10:
#line 137 "dlamrg.f"
    if (n1sv > 0 && n2sv > 0) {
#line 138 "dlamrg.f"
	if (a[ind1] <= a[ind2]) {
#line 139 "dlamrg.f"
	    index[i__] = ind1;
#line 140 "dlamrg.f"
	    ++i__;
#line 141 "dlamrg.f"
	    ind1 += *dtrd1;
#line 142 "dlamrg.f"
	    --n1sv;
#line 143 "dlamrg.f"
	} else {
#line 144 "dlamrg.f"
	    index[i__] = ind2;
#line 145 "dlamrg.f"
	    ++i__;
#line 146 "dlamrg.f"
	    ind2 += *dtrd2;
#line 147 "dlamrg.f"
	    --n2sv;
#line 148 "dlamrg.f"
	}
#line 149 "dlamrg.f"
	goto L10;
#line 150 "dlamrg.f"
    }
/*     end while */
#line 152 "dlamrg.f"
    if (n1sv == 0) {
#line 153 "dlamrg.f"
	i__1 = n2sv;
#line 153 "dlamrg.f"
	for (n1sv = 1; n1sv <= i__1; ++n1sv) {
#line 154 "dlamrg.f"
	    index[i__] = ind2;
#line 155 "dlamrg.f"
	    ++i__;
#line 156 "dlamrg.f"
	    ind2 += *dtrd2;
#line 157 "dlamrg.f"
/* L20: */
#line 157 "dlamrg.f"
	}
#line 158 "dlamrg.f"
    } else {
/*     N2SV .EQ. 0 */
#line 160 "dlamrg.f"
	i__1 = n1sv;
#line 160 "dlamrg.f"
	for (n2sv = 1; n2sv <= i__1; ++n2sv) {
#line 161 "dlamrg.f"
	    index[i__] = ind1;
#line 162 "dlamrg.f"
	    ++i__;
#line 163 "dlamrg.f"
	    ind1 += *dtrd1;
#line 164 "dlamrg.f"
/* L30: */
#line 164 "dlamrg.f"
	}
#line 165 "dlamrg.f"
    }

#line 167 "dlamrg.f"
    return 0;

/*     End of DLAMRG */

} /* dlamrg_ */

