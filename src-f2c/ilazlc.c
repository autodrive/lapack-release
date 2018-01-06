#line 1 "ilazlc.f"
/* ilazlc.f -- translated by f2c (version 20100827).
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

#line 1 "ilazlc.f"
/* > \brief \b ILAZLC scans a matrix for its last non-zero column. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILAZLC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILAZLC( M, N, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            M, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILAZLC scans A for its last non-zero column. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The m by n matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
integer ilazlc_(integer *m, integer *n, doublecomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1, i__2;

    /* Local variables */
    static integer i__;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick test for the common case where one corner is non-zero. */
#line 105 "ilazlc.f"
    /* Parameter adjustments */
#line 105 "ilazlc.f"
    a_dim1 = *lda;
#line 105 "ilazlc.f"
    a_offset = 1 + a_dim1;
#line 105 "ilazlc.f"
    a -= a_offset;
#line 105 "ilazlc.f"

#line 105 "ilazlc.f"
    /* Function Body */
#line 105 "ilazlc.f"
    if (*n == 0) {
#line 106 "ilazlc.f"
	ret_val = *n;
#line 107 "ilazlc.f"
    } else /* if(complicated condition) */ {
#line 107 "ilazlc.f"
	i__1 = *n * a_dim1 + 1;
#line 107 "ilazlc.f"
	i__2 = *m + *n * a_dim1;
#line 107 "ilazlc.f"
	if (a[i__1].r != 0. || a[i__1].i != 0. || (a[i__2].r != 0. || a[i__2]
		.i != 0.)) {
#line 108 "ilazlc.f"
	    ret_val = *n;
#line 109 "ilazlc.f"
	} else {
/*     Now scan each column from the end, returning with the first non-zero. */
#line 111 "ilazlc.f"
	    for (ret_val = *n; ret_val >= 1; --ret_val) {
#line 112 "ilazlc.f"
		i__1 = *m;
#line 112 "ilazlc.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 113 "ilazlc.f"
		    i__2 = i__ + ret_val * a_dim1;
#line 113 "ilazlc.f"
		    if (a[i__2].r != 0. || a[i__2].i != 0.) {
#line 113 "ilazlc.f"
			return ret_val;
#line 113 "ilazlc.f"
		    }
#line 114 "ilazlc.f"
		}
#line 115 "ilazlc.f"
	    }
#line 116 "ilazlc.f"
	}
#line 116 "ilazlc.f"
    }
#line 117 "ilazlc.f"
    return ret_val;
} /* ilazlc_ */

