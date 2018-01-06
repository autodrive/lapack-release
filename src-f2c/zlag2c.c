#line 1 "zlag2c.f"
/* zlag2c.f -- translated by f2c (version 20100827).
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

#line 1 "zlag2c.f"
/* > \brief \b ZLAG2C converts a complex double precision matrix to a complex single precision matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAG2C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlag2c.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlag2c.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlag2c.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAG2C( M, N, A, LDA, SA, LDSA, INFO ) */

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
/* > ZLAG2C converts a COMPLEX*16 matrix, SA, to a COMPLEX matrix, A. */
/* > */
/* > RMAX is the overflow for the SINGLE PRECISION arithmetic */
/* > ZLAG2C checks that all the entries of A are between -RMAX and */
/* > RMAX. If not the conversion is aborted and a flag is raised. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          SA is COMPLEX array, dimension (LDSA,N) */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlag2c_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *sa, integer *ldsa, integer *info)
{
    /* System generated locals */
    integer sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal rmax;
    extern doublereal slamch_(char *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 138 "zlag2c.f"
    /* Parameter adjustments */
#line 138 "zlag2c.f"
    a_dim1 = *lda;
#line 138 "zlag2c.f"
    a_offset = 1 + a_dim1;
#line 138 "zlag2c.f"
    a -= a_offset;
#line 138 "zlag2c.f"
    sa_dim1 = *ldsa;
#line 138 "zlag2c.f"
    sa_offset = 1 + sa_dim1;
#line 138 "zlag2c.f"
    sa -= sa_offset;
#line 138 "zlag2c.f"

#line 138 "zlag2c.f"
    /* Function Body */
#line 138 "zlag2c.f"
    rmax = slamch_("O", (ftnlen)1);
#line 139 "zlag2c.f"
    i__1 = *n;
#line 139 "zlag2c.f"
    for (j = 1; j <= i__1; ++j) {
#line 140 "zlag2c.f"
	i__2 = *m;
#line 140 "zlag2c.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 141 "zlag2c.f"
	    i__3 = i__ + j * a_dim1;
#line 141 "zlag2c.f"
	    i__4 = i__ + j * a_dim1;
#line 141 "zlag2c.f"
	    if (a[i__3].r < -rmax || a[i__4].r > rmax || d_imag(&a[i__ + j * 
		    a_dim1]) < -rmax || d_imag(&a[i__ + j * a_dim1]) > rmax) {
#line 145 "zlag2c.f"
		*info = 1;
#line 146 "zlag2c.f"
		goto L30;
#line 147 "zlag2c.f"
	    }
#line 148 "zlag2c.f"
	    i__3 = i__ + j * sa_dim1;
#line 148 "zlag2c.f"
	    i__4 = i__ + j * a_dim1;
#line 148 "zlag2c.f"
	    sa[i__3].r = a[i__4].r, sa[i__3].i = a[i__4].i;
#line 149 "zlag2c.f"
/* L10: */
#line 149 "zlag2c.f"
	}
#line 150 "zlag2c.f"
/* L20: */
#line 150 "zlag2c.f"
    }
#line 151 "zlag2c.f"
    *info = 0;
#line 152 "zlag2c.f"
L30:
#line 153 "zlag2c.f"
    return 0;

/*     End of ZLAG2C */

} /* zlag2c_ */

