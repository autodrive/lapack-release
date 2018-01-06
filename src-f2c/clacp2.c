#line 1 "clacp2.f"
/* clacp2.f -- translated by f2c (version 20100827).
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

#line 1 "clacp2.f"
/* > \brief \b CLACP2 copies all or part of a real two-dimensional array to a complex array. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLACP2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacp2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacp2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacp2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLACP2( UPLO, M, N, A, LDA, B, LDB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            LDA, LDB, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ) */
/*       COMPLEX            B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACP2 copies all or part of a real two-dimensional matrix A to a */
/* > complex matrix B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies the part of the matrix A to be copied to B. */
/* >          = 'U':      Upper triangular part */
/* >          = 'L':      Lower triangular part */
/* >          Otherwise:  All of the matrix A */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          The m by n matrix A.  If UPLO = 'U', only the upper trapezium */
/* >          is accessed; if UPLO = 'L', only the lower trapezium is */
/* >          accessed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On exit, B = A in the locations specified by UPLO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,M). */
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
/* Subroutine */ int clacp2_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublecomplex *b, integer *ldb, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 135 "clacp2.f"
    /* Parameter adjustments */
#line 135 "clacp2.f"
    a_dim1 = *lda;
#line 135 "clacp2.f"
    a_offset = 1 + a_dim1;
#line 135 "clacp2.f"
    a -= a_offset;
#line 135 "clacp2.f"
    b_dim1 = *ldb;
#line 135 "clacp2.f"
    b_offset = 1 + b_dim1;
#line 135 "clacp2.f"
    b -= b_offset;
#line 135 "clacp2.f"

#line 135 "clacp2.f"
    /* Function Body */
#line 135 "clacp2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 136 "clacp2.f"
	i__1 = *n;
#line 136 "clacp2.f"
	for (j = 1; j <= i__1; ++j) {
#line 137 "clacp2.f"
	    i__2 = min(j,*m);
#line 137 "clacp2.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 138 "clacp2.f"
		i__3 = i__ + j * b_dim1;
#line 138 "clacp2.f"
		i__4 = i__ + j * a_dim1;
#line 138 "clacp2.f"
		b[i__3].r = a[i__4], b[i__3].i = 0.;
#line 139 "clacp2.f"
/* L10: */
#line 139 "clacp2.f"
	    }
#line 140 "clacp2.f"
/* L20: */
#line 140 "clacp2.f"
	}

#line 142 "clacp2.f"
    } else if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 143 "clacp2.f"
	i__1 = *n;
#line 143 "clacp2.f"
	for (j = 1; j <= i__1; ++j) {
#line 144 "clacp2.f"
	    i__2 = *m;
#line 144 "clacp2.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 145 "clacp2.f"
		i__3 = i__ + j * b_dim1;
#line 145 "clacp2.f"
		i__4 = i__ + j * a_dim1;
#line 145 "clacp2.f"
		b[i__3].r = a[i__4], b[i__3].i = 0.;
#line 146 "clacp2.f"
/* L30: */
#line 146 "clacp2.f"
	    }
#line 147 "clacp2.f"
/* L40: */
#line 147 "clacp2.f"
	}

#line 149 "clacp2.f"
    } else {
#line 150 "clacp2.f"
	i__1 = *n;
#line 150 "clacp2.f"
	for (j = 1; j <= i__1; ++j) {
#line 151 "clacp2.f"
	    i__2 = *m;
#line 151 "clacp2.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 152 "clacp2.f"
		i__3 = i__ + j * b_dim1;
#line 152 "clacp2.f"
		i__4 = i__ + j * a_dim1;
#line 152 "clacp2.f"
		b[i__3].r = a[i__4], b[i__3].i = 0.;
#line 153 "clacp2.f"
/* L50: */
#line 153 "clacp2.f"
	    }
#line 154 "clacp2.f"
/* L60: */
#line 154 "clacp2.f"
	}
#line 155 "clacp2.f"
    }

#line 157 "clacp2.f"
    return 0;

/*     End of CLACP2 */

} /* clacp2_ */

