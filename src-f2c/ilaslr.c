#line 1 "ilaslr.f"
/* ilaslr.f -- translated by f2c (version 20100827).
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

#line 1 "ilaslr.f"
/* > \brief \b ILASLR scans a matrix for its last non-zero row. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILASLR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILASLR( M, N, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            M, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILASLR scans A for its last non-zero row. */
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
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
integer ilaslr_(integer *m, integer *n, doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1;

    /* Local variables */
    static integer i__, j;


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
#line 105 "ilaslr.f"
    /* Parameter adjustments */
#line 105 "ilaslr.f"
    a_dim1 = *lda;
#line 105 "ilaslr.f"
    a_offset = 1 + a_dim1;
#line 105 "ilaslr.f"
    a -= a_offset;
#line 105 "ilaslr.f"

#line 105 "ilaslr.f"
    /* Function Body */
#line 105 "ilaslr.f"
    if (*m == 0) {
#line 106 "ilaslr.f"
	ret_val = *m;
#line 107 "ilaslr.f"
    } else if (a[*m + a_dim1] != 0. || a[*m + *n * a_dim1] != 0.) {
#line 108 "ilaslr.f"
	ret_val = *m;
#line 109 "ilaslr.f"
    } else {
/*     Scan up each column tracking the last zero row seen. */
#line 111 "ilaslr.f"
	ret_val = 0;
#line 112 "ilaslr.f"
	i__1 = *n;
#line 112 "ilaslr.f"
	for (j = 1; j <= i__1; ++j) {
#line 113 "ilaslr.f"
	    i__ = *m;
#line 114 "ilaslr.f"
	    while(a[max(i__,1) + j * a_dim1] == 0. && i__ >= 1) {
#line 115 "ilaslr.f"
		--i__;
#line 116 "ilaslr.f"
	    }
#line 117 "ilaslr.f"
	    ret_val = max(ret_val,i__);
#line 118 "ilaslr.f"
	}
#line 119 "ilaslr.f"
    }
#line 120 "ilaslr.f"
    return ret_val;
} /* ilaslr_ */

