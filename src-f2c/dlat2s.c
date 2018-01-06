#line 1 "dlat2s.f"
/* dlat2s.f -- translated by f2c (version 20100827).
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

#line 1 "dlat2s.f"
/* > \brief \b DLAT2S converts a double-precision triangular matrix to a single-precision triangular matrix. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAT2S + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlat2s.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlat2s.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlat2s.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAT2S( UPLO, N, A, LDA, SA, LDSA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDSA, N */
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
/* > DLAT2S converts a DOUBLE PRECISION triangular matrix, SA, to a SINGLE */
/* > PRECISION triangular matrix, A. */
/* > */
/* > RMAX is the overflow for the SINGLE PRECISION arithmetic */
/* > DLAS2S checks that all the entries of A are between -RMAX and */
/* > RMAX. If not the convertion is aborted and a flag is raised. */
/* > */
/* > This is an auxiliary routine so there is no argument checking. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows and columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the N-by-N triangular coefficient matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SA */
/* > \verbatim */
/* >          SA is REAL array, dimension (LDSA,N) */
/* >          Only the UPLO part of SA is referenced.  On exit, if INFO=0, */
/* >          the N-by-N coefficient matrix SA; if INFO>0, the content of */
/* >          the UPLO part of SA is unspecified. */
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
/* >                of the UPLO part of SA in exit is unspecified. */
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
/* Subroutine */ int dlat2s_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *sa, integer *ldsa, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal rmax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
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

#line 142 "dlat2s.f"
    /* Parameter adjustments */
#line 142 "dlat2s.f"
    a_dim1 = *lda;
#line 142 "dlat2s.f"
    a_offset = 1 + a_dim1;
#line 142 "dlat2s.f"
    a -= a_offset;
#line 142 "dlat2s.f"
    sa_dim1 = *ldsa;
#line 142 "dlat2s.f"
    sa_offset = 1 + sa_dim1;
#line 142 "dlat2s.f"
    sa -= sa_offset;
#line 142 "dlat2s.f"

#line 142 "dlat2s.f"
    /* Function Body */
#line 142 "dlat2s.f"
    rmax = slamch_("O", (ftnlen)1);
#line 143 "dlat2s.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 144 "dlat2s.f"
    if (upper) {
#line 145 "dlat2s.f"
	i__1 = *n;
#line 145 "dlat2s.f"
	for (j = 1; j <= i__1; ++j) {
#line 146 "dlat2s.f"
	    i__2 = j;
#line 146 "dlat2s.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 147 "dlat2s.f"
		if (a[i__ + j * a_dim1] < -rmax || a[i__ + j * a_dim1] > rmax)
			 {
#line 149 "dlat2s.f"
		    *info = 1;
#line 150 "dlat2s.f"
		    goto L50;
#line 151 "dlat2s.f"
		}
#line 152 "dlat2s.f"
		sa[i__ + j * sa_dim1] = a[i__ + j * a_dim1];
#line 153 "dlat2s.f"
/* L10: */
#line 153 "dlat2s.f"
	    }
#line 154 "dlat2s.f"
/* L20: */
#line 154 "dlat2s.f"
	}
#line 155 "dlat2s.f"
    } else {
#line 156 "dlat2s.f"
	i__1 = *n;
#line 156 "dlat2s.f"
	for (j = 1; j <= i__1; ++j) {
#line 157 "dlat2s.f"
	    i__2 = *n;
#line 157 "dlat2s.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 158 "dlat2s.f"
		if (a[i__ + j * a_dim1] < -rmax || a[i__ + j * a_dim1] > rmax)
			 {
#line 160 "dlat2s.f"
		    *info = 1;
#line 161 "dlat2s.f"
		    goto L50;
#line 162 "dlat2s.f"
		}
#line 163 "dlat2s.f"
		sa[i__ + j * sa_dim1] = a[i__ + j * a_dim1];
#line 164 "dlat2s.f"
/* L30: */
#line 164 "dlat2s.f"
	    }
#line 165 "dlat2s.f"
/* L40: */
#line 165 "dlat2s.f"
	}
#line 166 "dlat2s.f"
    }
#line 167 "dlat2s.f"
L50:

#line 169 "dlat2s.f"
    return 0;

/*     End of DLAT2S */

} /* dlat2s_ */

