#line 1 "dlag2s.f"
/* dlag2s.f -- translated by f2c (version 20100827).
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

#line 1 "dlag2s.f"
/* > \brief \b DLAG2S converts a double precision matrix to a single precision matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAG2S + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlag2s.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlag2s.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlag2s.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDSA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               SA( LDSA, * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAG2S converts a DOUBLE PRECISION matrix, SA, to a SINGLE */
/* > PRECISION matrix, A. */
/* > */
/* > RMAX is the overflow for the SINGLE PRECISION arithmetic */
/* > DLAG2S checks that all the entries of A are between -RMAX and */
/* > RMAX. If not the convertion is aborted and a flag is raised. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N coefficient matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] SA */
/* > \verbatim */
/* >          SA is REAL array, dimension (LDSA,N) */
/* >          On exit, if INFO=0, the M-by-N coefficient matrix SA; if */
/* >          INFO>0, the content of SA is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSA */
/* > \verbatim */
/* >          LDSA is INTEGER */
/* >          The leading dimension of the array SA.  LDSA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          = 1:  an entry of the matrix A is greater than the SINGLE */
/* >                PRECISION overflow threshold, in this case, the content */
/* >                of SA in exit is unspecified. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlag2s_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *sa, integer *ldsa, integer *info)
{
    /* System generated locals */
    integer sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal rmax;
    extern doublereal slamch_(char *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 136 "dlag2s.f"
    /* Parameter adjustments */
#line 136 "dlag2s.f"
    a_dim1 = *lda;
#line 136 "dlag2s.f"
    a_offset = 1 + a_dim1;
#line 136 "dlag2s.f"
    a -= a_offset;
#line 136 "dlag2s.f"
    sa_dim1 = *ldsa;
#line 136 "dlag2s.f"
    sa_offset = 1 + sa_dim1;
#line 136 "dlag2s.f"
    sa -= sa_offset;
#line 136 "dlag2s.f"

#line 136 "dlag2s.f"
    /* Function Body */
#line 136 "dlag2s.f"
    rmax = slamch_("O", (ftnlen)1);
#line 137 "dlag2s.f"
    i__1 = *n;
#line 137 "dlag2s.f"
    for (j = 1; j <= i__1; ++j) {
#line 138 "dlag2s.f"
	i__2 = *m;
#line 138 "dlag2s.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 139 "dlag2s.f"
	    if (a[i__ + j * a_dim1] < -rmax || a[i__ + j * a_dim1] > rmax) {
#line 140 "dlag2s.f"
		*info = 1;
#line 141 "dlag2s.f"
		goto L30;
#line 142 "dlag2s.f"
	    }
#line 143 "dlag2s.f"
	    sa[i__ + j * sa_dim1] = a[i__ + j * a_dim1];
#line 144 "dlag2s.f"
/* L10: */
#line 144 "dlag2s.f"
	}
#line 145 "dlag2s.f"
/* L20: */
#line 145 "dlag2s.f"
    }
#line 146 "dlag2s.f"
    *info = 0;
#line 147 "dlag2s.f"
L30:
#line 148 "dlag2s.f"
    return 0;

/*     End of DLAG2S */

} /* dlag2s_ */

