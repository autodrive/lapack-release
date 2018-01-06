#line 1 "ilaclr.f"
/* ilaclr.f -- translated by f2c (version 20100827).
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

#line 1 "ilaclr.f"
/* > \brief \b ILACLR scans a matrix for its last non-zero row. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ILACLR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION ILACLR( M, N, A, LDA ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            M, N, LDA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ILACLR scans A for its last non-zero row. */
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
/* >          A is array, dimension (LDA,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
integer ilaclr_(integer *m, integer *n, doublecomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1, i__2;

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
#line 105 "ilaclr.f"
    /* Parameter adjustments */
#line 105 "ilaclr.f"
    a_dim1 = *lda;
#line 105 "ilaclr.f"
    a_offset = 1 + a_dim1;
#line 105 "ilaclr.f"
    a -= a_offset;
#line 105 "ilaclr.f"

#line 105 "ilaclr.f"
    /* Function Body */
#line 105 "ilaclr.f"
    if (*m == 0) {
#line 106 "ilaclr.f"
	ret_val = *m;
#line 107 "ilaclr.f"
    } else /* if(complicated condition) */ {
#line 107 "ilaclr.f"
	i__1 = *m + a_dim1;
#line 107 "ilaclr.f"
	i__2 = *m + *n * a_dim1;
#line 107 "ilaclr.f"
	if (a[i__1].r != 0. || a[i__1].i != 0. || (a[i__2].r != 0. || a[i__2]
		.i != 0.)) {
#line 108 "ilaclr.f"
	    ret_val = *m;
#line 109 "ilaclr.f"
	} else {
/*     Scan up each column tracking the last zero row seen. */
#line 111 "ilaclr.f"
	    ret_val = 0;
#line 112 "ilaclr.f"
	    i__1 = *n;
#line 112 "ilaclr.f"
	    for (j = 1; j <= i__1; ++j) {
#line 113 "ilaclr.f"
		i__ = *m;
#line 114 "ilaclr.f"
		for(;;) { /* while(complicated condition) */
#line 114 "ilaclr.f"
		    i__2 = max(i__,1) + j * a_dim1;
#line 114 "ilaclr.f"
		    if (!(a[i__2].r == 0. && a[i__2].i == 0. && i__ >= 1))
#line 114 "ilaclr.f"
		    	break;
#line 115 "ilaclr.f"
		    --i__;
#line 116 "ilaclr.f"
		}
#line 117 "ilaclr.f"
		ret_val = max(ret_val,i__);
#line 118 "ilaclr.f"
	    }
#line 119 "ilaclr.f"
	}
#line 119 "ilaclr.f"
    }
#line 120 "ilaclr.f"
    return ret_val;
} /* ilaclr_ */

