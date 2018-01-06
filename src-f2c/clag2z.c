#line 1 "clag2z.f"
/* clag2z.f -- translated by f2c (version 20100827).
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

#line 1 "clag2z.f"
/* > \brief \b CLAG2Z converts a complex single precision matrix to a complex double precision matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLAG2Z + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clag2z.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clag2z.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clag2z.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAG2Z( M, N, SA, LDSA, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDSA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            SA( LDSA, * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAG2Z converts a COMPLEX matrix, SA, to a COMPLEX*16 matrix, A. */
/* > */
/* > Note that while it is possible to overflow while converting */
/* > from double to single, it is not possible to overflow when */
/* > converting from single to double. */
/* > */
/* > This is an auxiliary routine so there is no argument checking. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of lines of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* >          SA is COMPLEX array, dimension (LDSA,N) */
/* >          On entry, the M-by-N coefficient matrix SA. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSA */
/* > \verbatim */
/* >          LDSA is INTEGER */
/* >          The leading dimension of the array SA.  LDSA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On exit, the M-by-N coefficient matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clag2z_(integer *m, integer *n, doublecomplex *sa, 
	integer *ldsa, doublecomplex *a, integer *lda, integer *info)
{
    /* System generated locals */
    integer sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

#line 126 "clag2z.f"
    /* Parameter adjustments */
#line 126 "clag2z.f"
    sa_dim1 = *ldsa;
#line 126 "clag2z.f"
    sa_offset = 1 + sa_dim1;
#line 126 "clag2z.f"
    sa -= sa_offset;
#line 126 "clag2z.f"
    a_dim1 = *lda;
#line 126 "clag2z.f"
    a_offset = 1 + a_dim1;
#line 126 "clag2z.f"
    a -= a_offset;
#line 126 "clag2z.f"

#line 126 "clag2z.f"
    /* Function Body */
#line 126 "clag2z.f"
    *info = 0;
#line 127 "clag2z.f"
    i__1 = *n;
#line 127 "clag2z.f"
    for (j = 1; j <= i__1; ++j) {
#line 128 "clag2z.f"
	i__2 = *m;
#line 128 "clag2z.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 129 "clag2z.f"
	    i__3 = i__ + j * a_dim1;
#line 129 "clag2z.f"
	    i__4 = i__ + j * sa_dim1;
#line 129 "clag2z.f"
	    a[i__3].r = sa[i__4].r, a[i__3].i = sa[i__4].i;
#line 130 "clag2z.f"
/* L10: */
#line 130 "clag2z.f"
	}
#line 131 "clag2z.f"
/* L20: */
#line 131 "clag2z.f"
    }
#line 132 "clag2z.f"
    return 0;

/*     End of CLAG2Z */

} /* clag2z_ */

