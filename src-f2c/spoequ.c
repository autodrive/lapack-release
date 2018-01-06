#line 1 "spoequ.f"
/* spoequ.f -- translated by f2c (version 20100827).
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

#line 1 "spoequ.f"
/* > \brief \b SPOEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spoequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spoequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spoequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOEQU( N, A, LDA, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOEQU computes row and column scalings intended to equilibrate a */
/* > symmetric positive definite matrix A and reduce its condition number */
/* > (with respect to the two-norm).  S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          The N-by-N symmetric positive definite matrix whose scaling */
/* >          factors are to be computed.  Only the diagonal elements of A */
/* >          are referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is REAL */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is REAL */
/* >          Absolute value of largest matrix element.  If AMAX is very */
/* >          close to overflow or very close to underflow, the matrix */
/* >          should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the i-th diagonal element is nonpositive. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int spoequ_(integer *n, doublereal *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal smin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 148 "spoequ.f"
    /* Parameter adjustments */
#line 148 "spoequ.f"
    a_dim1 = *lda;
#line 148 "spoequ.f"
    a_offset = 1 + a_dim1;
#line 148 "spoequ.f"
    a -= a_offset;
#line 148 "spoequ.f"
    --s;
#line 148 "spoequ.f"

#line 148 "spoequ.f"
    /* Function Body */
#line 148 "spoequ.f"
    *info = 0;
#line 149 "spoequ.f"
    if (*n < 0) {
#line 150 "spoequ.f"
	*info = -1;
#line 151 "spoequ.f"
    } else if (*lda < max(1,*n)) {
#line 152 "spoequ.f"
	*info = -3;
#line 153 "spoequ.f"
    }
#line 154 "spoequ.f"
    if (*info != 0) {
#line 155 "spoequ.f"
	i__1 = -(*info);
#line 155 "spoequ.f"
	xerbla_("SPOEQU", &i__1, (ftnlen)6);
#line 156 "spoequ.f"
	return 0;
#line 157 "spoequ.f"
    }

/*     Quick return if possible */

#line 161 "spoequ.f"
    if (*n == 0) {
#line 162 "spoequ.f"
	*scond = 1.;
#line 163 "spoequ.f"
	*amax = 0.;
#line 164 "spoequ.f"
	return 0;
#line 165 "spoequ.f"
    }

/*     Find the minimum and maximum diagonal elements. */

#line 169 "spoequ.f"
    s[1] = a[a_dim1 + 1];
#line 170 "spoequ.f"
    smin = s[1];
#line 171 "spoequ.f"
    *amax = s[1];
#line 172 "spoequ.f"
    i__1 = *n;
#line 172 "spoequ.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 173 "spoequ.f"
	s[i__] = a[i__ + i__ * a_dim1];
/* Computing MIN */
#line 174 "spoequ.f"
	d__1 = smin, d__2 = s[i__];
#line 174 "spoequ.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 175 "spoequ.f"
	d__1 = *amax, d__2 = s[i__];
#line 175 "spoequ.f"
	*amax = max(d__1,d__2);
#line 176 "spoequ.f"
/* L10: */
#line 176 "spoequ.f"
    }

#line 178 "spoequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 182 "spoequ.f"
	i__1 = *n;
#line 182 "spoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 183 "spoequ.f"
	    if (s[i__] <= 0.) {
#line 184 "spoequ.f"
		*info = i__;
#line 185 "spoequ.f"
		return 0;
#line 186 "spoequ.f"
	    }
#line 187 "spoequ.f"
/* L20: */
#line 187 "spoequ.f"
	}
#line 188 "spoequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 193 "spoequ.f"
	i__1 = *n;
#line 193 "spoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "spoequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 195 "spoequ.f"
/* L30: */
#line 195 "spoequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 199 "spoequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 200 "spoequ.f"
    }
#line 201 "spoequ.f"
    return 0;

/*     End of SPOEQU */

} /* spoequ_ */

